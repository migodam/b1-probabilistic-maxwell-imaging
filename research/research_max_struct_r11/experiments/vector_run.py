"""R11 no-truth inversion runner, requires explicit root-frozen authorization."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import time
ENTRY=time.perf_counter()
from pathlib import Path
import sys,json,argparse,traceback
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.experiments.vector_run import initialize as original_initialize,configure,dump,stamp,save_physical

def initialize(config,output,side):
    if not config['r11'].get('retain_allocator_cache',False):return original_initialize(config,output,side)
    import torch
    from research.research_max_struct_r10.src.vector_grid import VectorGrid
    from research.research_max_struct_r10.src.vector_resources import CooperativeGate
    from research.research_max_struct_r11.src.vector_cuda_dense import CachedDenseResidentVectorModel
    if config['execution_backend']!='cuda-dense' or side!=12:raise ValueError('Thin allocator change only for same12³ retained-dense CUDA model')
    gate=CooperativeGate(output/'GPU_BOUNDARIES.jsonl',memory_limit_bytes=config['execution']['memory_limit_bytes'])
    gate.wait('before_any_model_or_CUDA');cpu=VectorGrid(side,config['physics']['width'],config['acquisition'])
    model=CachedDenseResidentVectorModel(cpu,gate);model.progress_file=output/'RHS_PROGRESS.jsonl'
    dump(output/'SHARED_SETUP.json',dict(setup=model.setup,CPU_counts=dict(cpu.counts),CUDA_counts=dict(model.counts),
        nested_in='PROCESS_COST.json',additive=False,all_source_vector_fields=True,source_shape=list(model.inc.shape),receiver_shape=list(model.S.shape)))
    print(json.dumps(dict(event='SHARED_MODEL_READY',side=side,vector_dimension=model.M,allocated_bytes=torch.cuda.memory_allocated(),reserved_bytes=torch.cuda.memory_reserved(),allocator_cache_retained=True)),flush=True)
    return model

def invert(config, output, queue, maximum, offset=0):
    import numpy as np
    import torch
    from research.research_max_struct_r10.src.vector_cuda_model import ExactReceiverGeometry,NoiseReceiverSubspace
    from research.research_max_struct_r11.src.vector_cie import auxiliary_solve
    from research.research_max_struct_r11.src.vector_denominator import VectorLowRankDenominator
    from research.research_max_struct_r10.src.vector_physical import PhysicalObjective,physical_corrector
    from research.research_max_struct_r10.src.vector_energy import current_receiver_energy,affine_representation_energy,scalar_material_energy
    tasks=queue['tasks'][offset:offset+maximum]
    if not tasks or offset+maximum>queue['authorized_task_count']:raise ValueError('Queue subset exceeds frozen authorization')
    global ACTIVE_MODEL
    model=initialize(config,output,12);ACTIVE_MODEL=model;geometry=ExactReceiverGeometry(model)
    basis_wall,basis_cpu=time.perf_counter(),time.process_time()
    basis_path=REPO/config['r11']['basis_path']
    with np.load(basis_path,allow_pickle=False) as z:
        material_operator=VectorLowRankDenominator(model.array(z['right_vectors']),model.array(z['images']),config['r11']['projection'],model.counts)
    if model.device.type=='cuda':torch.cuda.synchronize()
    basis_wall=time.perf_counter()-basis_wall;basis_cpu=time.process_time()-basis_cpu
    dump(output/'DENOMINATOR_SETUP.json',dict(basis_path=str(basis_path),rank=material_operator.rank,M=model.M,N=model.N,projection=material_operator.projection,allstrong_certified=False,GT_access=False,
        CPU_seconds=basis_cpu,wall_seconds=basis_wall,scope='Basis load/upload/GalerkinB cold assembly; nested in process shared setup; source-basis construction on Mac separately charged'))
    dump(output/'RECEIVER_EXACT_SVD.json',geometry.info)
    dump(output/'SHARED_BACKGROUND_COST.json',dict(CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-ENTRY,
        basis_CPU_seconds=basis_cpu,basis_wall_seconds=basis_wall,
        common_CPU_excluding_basis_seconds=time.process_time()-basis_cpu,
        common_wall_excluding_basis_seconds=time.perf_counter()-ENTRY-basis_wall,
        gate_helper_CPU_seconds=getattr(model.gate,'helper_CPU_seconds',0.),
        shared_pause_wall_seconds=getattr(model.gate,'paused_wall_seconds',0.),
        scope='Whole cold process prefix: imports/configuration/sharedG/E/S/FFT/receiverSVD/basis upload and IO; nested in whole process, not additive to it',
        standalone='Ordinary CIE methods use common setup excluding basis; matrix methods add basis; publicY DPCS cold setup is case-specific and recorded in receiver info'))
    settings=dict(mf=config['inverse']['native_MF'],beta_material=config['inverse']['beta_material'],
        beta_current_dual=config['inverse']['beta_current_dual'],max_iterations_per_stage=500,
        minimum_iterations_per_stage=20,coefficient_delta_RMS_tolerance=.001,
        initial_chi=.25,material_upper=8.,material_mode='ordinary')
    current_id=None;terminals=[];case_views=[]
    for task_index,task in enumerate(tasks):
        boundary_wall,boundary_cpu=time.perf_counter(),time.process_time()
        boundary_helper=getattr(model.gate,'helper_CPU_seconds',0.);boundary_pause=getattr(model.gate,'paused_wall_seconds',0.)
        model.gate.wait('before_task_'+task['task_id'],torch.cuda.memory_reserved() if model.device.type=='cuda' else 0)
        boundary_wall=time.perf_counter()-boundary_wall;boundary_cpu=time.process_time()-boundary_cpu
        if task['public_id']!=current_id:
            current_id=task['public_id'];public=REPO/config['r11']['public_root']/current_id
            with np.load(public/'observations.npz',allow_pickle=False) as z:
                if z.files!=['Y']:raise ValueError('Only public Y is permitted in inverse NPZ')
                Y=z['Y'].copy()
            noise=json.loads((public/'PUBLIC.json').read_text())
            if Y.shape!=(180,60) or not np.isfinite(Y).all():raise ValueError('Invalid all-source vector observations')
            receiver=NoiseReceiverSubspace(model,geometry,Y,noise)
            case_views.append(dict(public_id=current_id,receiver=receiver.info,
                input_evidence=[stamp(public/'observations.npz'),stamp(public/'PUBLIC.json')]))
        # Share immutable G/FFT/E/S/SVD/DPCS only. Each arm starts a fresh
        # material-dependent state cache and independent optimizer history.
        if hasattr(model,'dense_state'):model.dense_state=None
        objective=PhysicalObjective(model,Y,upper=8.)
        dest=ROOT/'delegated/vector/inversions'/task['task_id'];dest.mkdir(parents=True,exist_ok=False)
        wall,cpu=time.perf_counter()-boundary_wall,time.process_time()-boundary_cpu;first_event=len(objective.events);before=dict(model.counts)
        first_batch=len(getattr(model,'batch_records',[]));helper_cpu_start=boundary_helper
        original_evaluate=objective.evaluate
        def paid_evaluate(*args,**kwargs):
            record=original_evaluate(*args,**kwargs)
            event=objective.events[-1]
            helper_cpu=getattr(model.gate,'helper_CPU_seconds',0.)-helper_cpu_start
            event.update(completion_arm_CPU_seconds=time.process_time()-cpu+helper_cpu,
                completion_arm_process_CPU_seconds=time.process_time()-cpu,
                completion_arm_GPU_query_helper_CPU_seconds=helper_cpu,
                completion_arm_wall_seconds=time.perf_counter()-wall,
                completion_arm_physical_LU=int(model.counts['physical_LU']-before.get('physical_LU',0)),
                sigma_whitened_squared_norm=float(np.sum(abs(record['prediction']-Y)**2)/float(noise['sigma'])**2))
            return record
        objective.evaluate=paid_evaluate
        device_start=device_end=None;pause_start=boundary_pause
        if model.device.type=='cuda':
            device_start=torch.cuda.Event(enable_timing=True);device_end=torch.cuda.Event(enable_timing=True);device_start.record()
        trace=[]
        def emit(row):
            trace.append(row)
            with (dest/'TRACE.jsonl').open('a') as handle:handle.write(json.dumps(row)+'\n')
            if row.get('kind')!='joint' or row.get('iteration',0)%20==0:print(json.dumps(dict(task_id=task['task_id'],**row)),flush=True)
        initial=objective.evaluate(np.full(model.N,.25),task['task_id']+'_uniform_initial',False)
        save_physical(dest,'physical_initial',initial)
        def callback(record,m,J,Q,gamma,basis):
            emit(dict(kind='joint',**record))
        def endpoint(label,m,state):
            arrays=dict(material=m.cpu().numpy(),J=state['J'].cpu().numpy(),Q=state['Q'].cpu().numpy(),
                gamma=np.empty((0,60),complex) if state['gamma'] is None else state['gamma'].cpu().numpy())
            np.savez_compressed(dest/(label+'_joint.npz'),**arrays)
            energy=dict(current=current_receiver_energy(state['J'],receiver,model),
                material=scalar_material_energy(arrays['material'],model.N,model.cpu.h))
            if state['gamma'] is not None:
                energy['representation']=affine_representation_energy(state['J'],state['gamma'],state['basis'],receiver,model)
                paid=objective.evaluate(arrays['material'],task['task_id']+'_'+label+'_physical',False)
                save_physical(dest,label+'_physical',paid)
                energy['physical_Phi']=paid['Phi']
            dump(dest/(label+'_ENERGY.json'),energy)
        settings['material_mode']=task['material_mode']
        if task['material_mode']=='ordinary':
            from research.research_max_struct_r10.src.vector_cie import auxiliary_solve as original_auxiliary_solve
            result=original_auxiliary_solve(model,receiver,model.array(Y),task['arm'],settings,callback,endpoint)
        else:
            result=auxiliary_solve(model,receiver,model.array(Y),'GlobalDual',settings,callback,endpoint,material_operator)
        last=result['material'].cpu().numpy()
        best,corrector=physical_corrector(objective,last,config['inverse']['corrector'],task['task_id'],
            lambda label,record:save_physical(dest,label,record),emit)
        save_physical(dest,'prediction',best)
        # All arrays here are already-paid CPU states; no forward/GT read.
        paid_states=list(objective.records.values())
        np.savez_compressed(dest/'PAID_PHYSICAL_MATERIALS.npz',material=np.stack([r['m'] for r in paid_states]))
        dump(dest/'PAID_PHYSICAL_STATE_INDEX.json',dict(states=[dict(row=i,state_id=r['id'],Phi=r['Phi']) for i,r in enumerate(paid_states)],
            events=objective.events,selection='No GT; exact paid physical materials and completion CPU/wall/LU coordinates',
            timing_scope='Arm includes initialization, auxiliary stages and all physical evaluations through each completed event; background setup separately charged'))
        if hasattr(model,'batch_records'):
            dump(dest/'PHYSICAL_BATCH_PROFILE.json',dict(records=model.batch_records[first_batch:],
                scope='Already measured R10 LU and sameLU multiRHS batches; nested within physical objectives; not additive to arm CPU/wall'))
        if config.get('cost_prefixes'):
            prefix_rows=[];materials=[];by_id={r['id']:r for r in paid_states}
            axes=dict(CPU_seconds='completion_arm_CPU_seconds',wall_seconds='completion_arm_wall_seconds',physical_LU='completion_arm_physical_LU')
            for axis,field in axes.items():
                for budget in config['cost_prefixes'][axis]:
                    eligible=[e for e in objective.events if e[field]<=budget]
                    if eligible:
                        event=min(eligible,key=lambda e:e['Phi']);r=by_id[event['state_id']];materials.append(r['m'])
                        prefix_rows.append(dict(axis=axis,budget=budget,status='PAID_PREFIX_SELECTED',material_row=len(materials)-1,
                            state_id=r['id'],Phi=r['Phi'],completed_cost=event[field],paid_event=event['index']))
                    else:prefix_rows.append(dict(axis=axis,budget=budget,status='NO_ELIGIBLE_PAID_STATE',material_row=None))
            np.savez_compressed(dest/'COST_PREFIX_MATERIALS.npz',material=np.stack(materials) if materials else np.empty((0,model.N)))
            dump(dest/'COST_PREFIX_SEAL.json',dict(status='PUBLIC_COST_PHI_SELECTION_SEALED_WITHOUT_GT',prefixes=prefix_rows,truth_access=False,
                selection='Min physical Phi among completed paid states under frozen CPU/wall/LU budgets; no interpolation/extra solve',CPU_includes_GPU_query_helper=True,
                shared_background_cost='Separate unique process setup and stand-alone cold views; arm coordinates do not duplicate shared setup'))
        device_ms=None
        if device_start is not None:
            device_end.record();device_end.synchronize();device_ms=device_start.elapsed_time(device_end)
        counts={k:int(v-before.get(k,0)) for k,v in model.counts.items() if v-before.get(k,0)}
        receipt=dict(status='SEALED_PUBLIC_PREDICTION_NO_SCORE',task=task,stages=result['stages'],physical_corrector=corrector,
            physical_Phi=best['Phi'],PG_inf=best['projected_gradient_inf'],counts=counts,
            CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
            GPU_query_helper_CPU_seconds=getattr(model.gate,'helper_CPU_seconds',0.)-helper_cpu_start,
            CPU_including_GPU_query_helper_seconds=time.process_time()-cpu+getattr(model.gate,'helper_CPU_seconds',0.)-helper_cpu_start,
            before_task_gate_wall_seconds=boundary_wall,shared_background_setup_excluded_from_arm=True,
            CUDA_event_stream_interval_ms=device_ms,
            CUDA_event_scope='Default-stream timestamp interval across this task; includes host launch/IO/idle and cooperative pause gaps; NOT active-kernel time or SM occupancy',
            cooperative_pause_wall_seconds=getattr(model.gate,'paused_wall_seconds',0.)-pause_start,
            GPU_peak_allocated_bytes_to_date=torch.cuda.max_memory_allocated() if model.device.type=='cuda' else 0,
            GPU_peak_reserved_bytes_to_date=torch.cuda.max_memory_reserved() if model.device.type=='cuda' else 0,
            additive=False,included_in='PROCESS_COST.json',paid_events=objective.events[first_event:],
            current_metric='ordinary_CIE',current_aggregation='sourcewise' if task['arm'] in ('sourceCIE','sourceDual') else 'global',
            current_beta=settings['beta_material'] if task['arm']=='sourceCIE' else settings['beta_current_dual'],
            material_beta=settings['beta_material'],material_metric=task['material_mode'],
            comparison_scope='Full/lagged match GlobalDual current; sourceDual controls sourcewise PR at dual beta; sourceCIE uses sourcewise PR and material beta',
            backend=model.setup['backend'],truth_access=False,
            auxiliary_state='Saved independent J/Q/gamma and real material at every native MF stage; not replaced with physical J',
            cross_arm_cache_scope='Immutable operators/SVD/publicY DPCS only; material-dependent LU/J/gradient and optimizer histories reset',
            physical_selection='Lowest public Phi among paid physical corrector states, no GT',
            representation_domain='J/Q/S-right modes are complex3N; material is realN, never directly put in dyadic modes')
        dump(dest/'RECEIPT.json',receipt)
        dump(dest/'SEALED.json',dict(status='PREDICTION_SEALED_BEFORE_ROOT_OFFLINE_SCORING',files=[stamp(p) for p in sorted(dest.iterdir()) if p.is_file()],truth_access=False))
        terminals.append(dict(task_id=task['task_id'],status=receipt['status']))
        dump(output/'TERMINALS.json',dict(terminals=terminals,requested=len(tasks),cases=case_views,
            all_requested_terminal=len(terminals)==len(tasks),truth_scoring='NOT_RUN_ROOT_ONLY'))
        print(json.dumps(dict(event='TASK_TERMINAL',task_id=task['task_id'],Phi=best['Phi'],CPU_seconds=receipt['CPU_seconds'],counts=counts)),flush=True)
    return model


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--config',required=True);parser.add_argument('--queue',required=True)
    parser.add_argument('--output',required=True);parser.add_argument('--max-tasks',type=int,required=True);parser.add_argument('--offset',type=int,default=0)
    args=parser.parse_args();config=json.loads((ROOT/args.config).read_text());queue=json.loads((ROOT/args.queue).read_text())
    if config.get('r11',{}).get('execution_authorized') is not True:raise RuntimeError('Root numerical configuration not approved')
    output=ROOT/args.output;output.mkdir(parents=True,exist_ok=False);model=None;failure=None;status='FAILED'
    try:
        import torch
        torch.set_num_threads(1);torch.set_num_interop_threads(1);process_config,job=configure()
        dump(output/'FROZEN_CONFIG.json',config);dump(output/'PROCESS_CONFIGURATION.json',process_config)
        if config['r11'].get('science_contract_path'):
            science_path=REPO/config['r11']['science_contract_path'];science=json.loads(science_path.read_text())
            if science['status']!='ROOT_FROZEN_BEFORE_NEW_GENERATION' or not science['execution_authorized']:raise RuntimeError('Missing frozen science contract')
            dump(output/'FROZEN_SCIENCE_CONTRACT.json',science);dump(output/'SCIENCE_CONTRACT_SOURCE.json',stamp(science_path))
        dump(output/'SOURCE_EVIDENCE.json',[stamp(p) for p in sorted((ROOT/'src').glob('vector_*.py'))]+[stamp(Path(__file__))])
        model=invert(config,output,queue,args.max_tasks,args.offset);status='COMPLETE'
    except BaseException as exc:
        failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());print(json.dumps(failure),flush=True)
    finally:
        import torch
        if model is None:model=globals().get('ACTIVE_MODEL')
        dump(output/'PROCESS_COST.json',dict(status=status,failure=failure,wall_seconds=time.perf_counter()-ENTRY,
            counts=dict(model.counts) if model else None,truth_access=False,whole_OS_CPU_seconds=None,
            CPU_seconds=time.process_time(),CPU_scope='Whole current Python process, including interpreter imports/setup/all arms/exports; time.process_time Windows kernel+user',
            GPU_peak_allocated_bytes=torch.cuda.max_memory_allocated() if torch.cuda.is_initialized() else 0,
            GPU_peak_reserved_bytes=torch.cuda.max_memory_reserved() if torch.cuda.is_initialized() else 0,
            gate_helper_CPU_seconds=model.gate.helper_CPU_seconds if model else None,
            gate_query_count=model.gate.query_count if model else None,additive=False))
    return 0 if status=='COMPLETE' else 1

if __name__=='__main__':raise SystemExit(main())
