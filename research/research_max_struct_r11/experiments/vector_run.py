"""R11 no-truth inversion runner, requires explicit root-frozen authorization."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import time
ENTRY=time.perf_counter()
from pathlib import Path
import sys,json,argparse,traceback,gc,subprocess
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.experiments.vector_run import initialize as original_initialize,configure,dump,stamp,save_physical

def initialize(config,output,side):
    if not config['r11'].get('retain_allocator_cache',False):
        model=original_initialize(config,output,side)
        if model.device.type=='cpu' and not hasattr(model,'batch_records'):model.batch_records=[]
        return model
    import torch
    from research.research_max_struct_r10.src.vector_grid import VectorGrid
    from research.research_max_struct_r10.src.vector_resources import CooperativeGate
    from research.research_max_struct_r11.src.vector_cuda_dense import CachedDenseResidentVectorModel
    if config['execution_backend']!='cuda-dense' or side!=12:raise ValueError('Thin allocator change only for same12³ retained-dense CUDA model')
    gate=CooperativeGate(output/'GPU_BOUNDARIES.jsonl',memory_limit_bytes=config['execution']['memory_limit_bytes'],
        minimum_query_interval_seconds=float(config['r11'].get('idle_snapshot_cache_seconds',1.)))
    gate.wait('before_any_model_or_CUDA');cpu=VectorGrid(side,config['physics']['width'],config['acquisition'])
    model=CachedDenseResidentVectorModel(cpu,gate);model.progress_file=output/'RHS_PROGRESS.jsonl'
    dump(output/'SHARED_SETUP.json',dict(setup=model.setup,CPU_counts=dict(cpu.counts),CUDA_counts=dict(model.counts),
        nested_in='PROCESS_COST.json',additive=False,all_source_vector_fields=True,source_shape=list(model.inc.shape),receiver_shape=list(model.S.shape)))
    print(json.dumps(dict(event='SHARED_MODEL_READY',side=side,vector_dimension=model.M,allocated_bytes=torch.cuda.memory_allocated(),reserved_bytes=torch.cuda.memory_reserved(),allocator_cache_retained=True)),flush=True)
    return model

def current_rss():
    if sys.platform=='win32':
        import psutil
        return int(psutil.Process().memory_info().rss),0.
    import resource
    process=subprocess.Popen(['/bin/ps','-o','rss=','-p',str(os.getpid())],stdout=subprocess.PIPE,text=True)
    value=process.stdout.read();_,status,usage=os.wait4(process.pid,0);process.returncode=os.waitstatus_to_exitcode(status)
    return int(value.strip())*1024,usage.ru_utime+usage.ru_stime

def invert(config, output, queue, maximum, offset=0):
    global ACTIVE_MODEL,ORIGINAL_REQUESTED_TASKS,ACTUAL_DISPATCHED_TASKS,ACTIVE_TERMINALS
    import numpy as np
    import torch
    from research.research_max_struct_r10.src.vector_cuda_model import ExactReceiverGeometry,NoiseReceiverSubspace
    from research.research_max_struct_r11.src.vector_cie import auxiliary_solve
    from research.research_max_struct_r11.src.vector_denominator import VectorLowRankDenominator
    from research.research_max_struct_r10.src.vector_physical import PhysicalObjective,physical_corrector
    from research.research_max_struct_r10.src.vector_energy import current_receiver_energy,affine_representation_energy,scalar_material_energy
    tasks=queue['tasks'][offset:offset+maximum]
    if not tasks or offset+maximum>queue['authorized_task_count']:raise ValueError('Queue subset exceeds frozen authorization')
    original_requested=len(tasks);host_role=config['r11'].get('host_dispatch_role')
    if host_role:
        dispatch_path=REPO/config['r11']['host_dispatch_path'];dispatch=json.loads(dispatch_path.read_text())
        mapping={t['task_id']:t['host'] for t in dispatch['tasks']}
        if len(mapping)!=50 or len(dispatch['tasks'])!=50:raise ValueError('Host dispatch must cover exactly50 unique frozen tasks')
        tasks=[t for t in tasks if mapping[t['task_id']]==host_role]
        dump(output/'HOST_DISPATCH.json',dict(source=stamp(dispatch_path),host=host_role,original_requested=original_requested,
            actual_dispatched=len(tasks),tasks=tasks,global_frozen_tasks=50,scientific_budget_unchanged=True))
    ORIGINAL_REQUESTED_TASKS=original_requested;ACTUAL_DISPATCHED_TASKS=len(tasks)
    if not tasks:
        ACTIVE_TERMINALS=[]
        dump(output/'ENGINEERING_DISPATCH_NOOP.json',dict(status='ZERO_FW_HOST_TRANSFER_NOOP',original_requested=original_requested,
            actual_scientific_tasks=0,new_physical_LU=0,new_Maxwell_RHS=0,host=host_role,
            reason='All requested object tasks moved to other host by pre-GT root engineering dispatch; not scientific terminal'))
        dump(output/'TERMINALS.json',dict(terminals=[],requested=0,original_requested=original_requested,actual_dispatched=0,
            all_requested_terminal=True,scientific_terminals=0,engineering_noop=True,truth_scoring='NOT_RUN'))
        print(json.dumps(dict(event='HOST_DISPATCH_ZERO_FW_NOOP',original_requested=original_requested,actual_dispatched=0,scientific_terminals=0)),flush=True)
        return None
    model=initialize(config,output,12);ACTIVE_MODEL=model
    engineering=None
    if config['r11'].get('public_linalg_comparison_states'):
        from research.research_max_struct_r11.src.vector_backend_selection import compare_public_linalg_libraries
        engineering=compare_public_linalg_libraries(model,REPO/config['r11']['public_linalg_comparison_states'],output/'PUBLIC_LINALG_COMPARISON')
        model.setup['confirmation_selected_linalg_library']=engineering['selected_library']
    elif config['r11'].get('reuse_public_linalg_selection') and model.device.type=='cuda':
        selection_path=REPO/config['r11']['reuse_public_linalg_selection'];selection=json.loads(selection_path.read_text())
        if selection['status']!='EXISTING_LIBRARY_SELECTED_BEFORE_ANY_CONFIRMATION_INVERSE':raise RuntimeError('Missing completed public6LU library qualification')
        torch.backends.cuda.preferred_linalg_library(selection['selected_library'])
        model.setup['confirmation_selected_linalg_library']=selection['selected_library']
        dump(output/'REUSED_PUBLIC_LINALG_SELECTION.json',dict(source=stamp(selection_path),receipt=selection,new_probe_LU=0,new_GT_access=False))
    geometry=ExactReceiverGeometry(model)
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
        engineering_CPU_seconds=engineering['CPU_seconds'] if engineering else 0.,engineering_wall_seconds=engineering['wall_seconds'] if engineering else 0.,
        common_CPU_excluding_basis_and_engineering_seconds=time.process_time()-basis_cpu-(engineering['CPU_seconds'] if engineering else 0.),
        common_wall_excluding_basis_and_engineering_seconds=time.perf_counter()-ENTRY-basis_wall-(engineering['wall_seconds'] if engineering else 0.),
        gate_helper_CPU_seconds=getattr(model.gate,'helper_CPU_seconds',0.),
        shared_pause_wall_seconds=getattr(model.gate,'paused_wall_seconds',0.),
        scope='Whole cold process prefix: imports/configuration/sharedG/E/S/FFT/receiverSVD/basis upload and IO; nested in whole process, not additive to it',
        standalone='Ordinary CIE methods use common setup excluding basis; matrix methods add basis; publicY DPCS cold setup is case-specific and recorded in receiver info'))
    settings=dict(mf=config['inverse']['native_MF'],beta_material=config['inverse']['beta_material'],
        beta_current_dual=config['inverse']['beta_current_dual'],max_iterations_per_stage=500,
        minimum_iterations_per_stage=20,coefficient_delta_RMS_tolerance=.001,
        initial_chi=.25,material_upper=8.,material_mode='ordinary')
    current_id=None;terminals=[];ACTIVE_TERMINALS=terminals;case_views=[]
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
        material_operator.latest_factor=None
        if model.device.type=='cpu' and config['r11'].get('spool_historical_currents'):
            from research.research_max_struct_r11.src.vector_spooled_physical import SpooledPhysicalObjective
            capacity=int(config['inverse']['corrector']['strict_maxfun'])+len(settings['mf'])+1
            objective=SpooledPhysicalObjective(model,Y,upper=8.,capacity=capacity,
                directory=ROOT/'delegated/vector/temporary_current_spool',owner=task['task_id'],verify_first=task_index==0)
        else:objective=PhysicalObjective(model,Y,upper=8.)
        dest=ROOT/config['r11'].get('inversion_output_root','delegated/vector/inversions')/task['task_id'];dest.mkdir(parents=True,exist_ok=False)
        wall,cpu=time.perf_counter()-boundary_wall,time.process_time()-boundary_cpu;first_event=len(objective.events);before=dict(model.counts)
        first_batch=len(getattr(model,'batch_records',[]));helper_cpu_start=boundary_helper
        original_evaluate=objective.evaluate
        def paid_evaluate(*args,**kwargs):
            record=original_evaluate(*args,**kwargs)
            event=objective.events[-1]
            if model.device.type=='cpu':
                if event['counts_delta'].get('physical_LU',0):model.batch_records.append(dict(model.last_dense_batch))
                if event['counts_delta'].get('physical_adjoint_RHS',0):model.batch_records.append(dict(model.last_dense_adjoint_batch))
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
        if hasattr(objective,'verify_initial_archive'):objective.verify_initial_archive(dest/'physical_initial.npz',initial)
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
                shared_background_cost='Separate unique process setup and stand-alone cold views; arm coordinates do not duplicate shared setup',
                infrastructure_attempt_metadata=task.get('infrastructure_attempt_metadata'),
                cost_scope=task.get('infrastructure_attempt_metadata',{}).get('cost_prefix_scope','NORMAL_FIRST_ATTEMPT_END_TO_END_ARM')))
        device_ms=None
        if device_start is not None:
            device_end.record();device_end.synchronize();device_ms=device_start.elapsed_time(device_end)
        counts={k:int(v-before.get(k,0)) for k,v in model.counts.items() if v-before.get(k,0)}
        receipt=dict(status='SEALED_PUBLIC_PREDICTION_NO_SCORE',task=task,infrastructure_attempt_metadata=task.get('infrastructure_attempt_metadata'),stages=result['stages'],physical_corrector=corrector,
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
            torch_linalg_library=str(torch.backends.cuda.preferred_linalg_library()) if model.device.type=='cuda' else 'SciPy CPU retained LU; CUDA library not applicable',
            allocator_policy=model.setup.get('R11_allocator_policy','Inherited R10 allocator/CPU policy'),host_platform=sys.platform,host_dispatch_role=host_role,
            auxiliary_state='Saved independent J/Q/gamma and real material at every native MF stage; not replaced with physical J',
            cross_arm_cache_scope='Immutable operators/SVD/publicY DPCS only; material-dependent LU/J/gradient and optimizer histories reset',
            physical_selection='Lowest public Phi among paid physical corrector states, no GT',
            representation_domain='J/Q/S-right modes are complex3N; material is realN, never directly put in dyadic modes')
        dump(dest/'RECEIPT.json',receipt)
        dump(dest/'SEALED.json',dict(status='PREDICTION_SEALED_BEFORE_ROOT_OFFLINE_SCORING',files=[stamp(p) for p in sorted(dest.iterdir()) if p.is_file()],truth_access=False))
        terminals.append(dict(task_id=task['task_id'],status=receipt['status']))
        dump(output/'TERMINALS.json',dict(terminals=terminals,requested=len(tasks),original_requested=original_requested,actual_dispatched=len(tasks),host_dispatch_role=host_role,cases=case_views,
            all_requested_terminal=len(terminals)==len(tasks),truth_scoring='NOT_RUN_ROOT_ONLY'))
        print(json.dumps(dict(event='TASK_TERMINAL',task_id=task['task_id'],Phi=best['Phi'],CPU_seconds=receipt['CPU_seconds'],counts=counts)),flush=True)
        if config['r11'].get('release_completed_arm_storage'):
            rss_before,helper1=current_rss()
            spool_report=objective.finish_spool() if hasattr(objective,'finish_spool') else None
            paid_states=None;by_id=None;objective=None;original_evaluate=None;paid_evaluate=None;best=None;initial=None
            model.dense_state=None;material_operator.latest_factor=None;gc.collect()
            rss_after,helper2=current_rss()
            if spool_report is not None:dump(dest/'MATERIAL_HISTORY_SPOOL.json',spool_report)
            dump(dest/'POST_ARCHIVE_STORAGE_RELEASE.json',dict(RSS_before_bytes=rss_before,RSS_after_bytes=rss_after,
                RSS_monitor_helper_CPU_seconds=helper1+helper2,all_archived_arrays_preserved=True,no_single_arm_state_or_optimizer_changed=True,
                scope='After scientific prediction/prefix seal; archive cleanup and monitoring nested in whole process'))
        stop_path=config['r11'].get('stop_after_arm_request_path')
        if stop_path and (ROOT/stop_path).exists():
            globals()['ENGINEERING_STOP_AFTER_ARM']=True
            dump(output/'CONTROLLED_STOP_AFTER_ARM.json',dict(status='NATURAL_SEALED_ARM_BOUNDARY_ENGINEERING_STOP',last_task_id=task['task_id'],
                completed_terminal_tasks=len(terminals),actual_dispatched_tasks=len(tasks),counts=dict(model.counts),
                no_active_optimization_interrupted=True,all_current_arm_arrays_archived=True,request_path=stop_path,truth_access=False))
            break
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
        if config['r11'].get('retry_own_visibility_timeout'):
            from research.research_max_struct_r11.src.vector_resource_visibility import install
            dump(output/'RESOURCE_VISIBILITY_ENGINEERING_POLICY.json',dict(**install(output/'VISIBILITY_TIMEOUT_RETRIES.jsonl'),
                successful_idle_snapshot_cache_seconds=float(config['r11'].get('idle_snapshot_cache_seconds',1.)),
                maximum_new_external_process_detection_delay_seconds=float(config['r11'].get('idle_snapshot_cache_seconds',1.)),
                external_or_expired_snapshot_requires_valid_query_before_CUDA=True))
        if config['r11'].get('science_contract_path'):
            science_path=REPO/config['r11']['science_contract_path'];science=json.loads(science_path.read_text())
            if science['status']!='ROOT_FROZEN_BEFORE_NEW_GENERATION' or not science['execution_authorized']:raise RuntimeError('Missing frozen science contract')
            dump(output/'FROZEN_SCIENCE_CONTRACT.json',science);dump(output/'SCIENCE_CONTRACT_SOURCE.json',stamp(science_path))
        dump(output/'SOURCE_EVIDENCE.json',[stamp(p) for p in sorted((ROOT/'src').glob('vector_*.py'))]+[stamp(Path(__file__))])
        model=invert(config,output,queue,args.max_tasks,args.offset)
        status=('COMPLETE_ENGINEERING_STOP_AFTER_ARM' if globals().get('ENGINEERING_STOP_AFTER_ARM') else
            'COMPLETE_ENGINEERING_DISPATCH_NOOP' if model is None and globals().get('ACTUAL_DISPATCHED_TASKS')==0 else 'COMPLETE')
    except BaseException as exc:
        failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());print(json.dumps(failure),flush=True)
    finally:
        import torch
        if model is None:model=globals().get('ACTIVE_MODEL')
        dump(output/'PROCESS_COST.json',dict(status=status,failure=failure,wall_seconds=time.perf_counter()-ENTRY,
            original_requested_tasks=globals().get('ORIGINAL_REQUESTED_TASKS',args.max_tasks),
            actual_dispatched_tasks=globals().get('ACTUAL_DISPATCHED_TASKS'),completed_terminal_tasks=len(globals().get('ACTIVE_TERMINALS',[])),
            host_dispatch_role=config['r11'].get('host_dispatch_role'),global_frozen_tasks=50 if config['r11'].get('host_dispatch_role') else None,
            counts=dict(model.counts) if model else None,truth_access=False,whole_OS_CPU_seconds=None,
            CPU_seconds=time.process_time(),CPU_scope='Whole current Python process, including interpreter imports/setup/all arms/exports; time.process_time Windows kernel+user',
            GPU_peak_allocated_bytes=torch.cuda.max_memory_allocated() if torch.cuda.is_initialized() else 0,
            GPU_peak_reserved_bytes=torch.cuda.max_memory_reserved() if torch.cuda.is_initialized() else 0,
            gate_helper_CPU_seconds=model.gate.helper_CPU_seconds if model else None,
            gate_query_count=model.gate.query_count if model else None,
            gate_cached_boundary_checks=model.gate.cached_boundary_checks if model else None,
            gate_successful_idle_cache_seconds=model.gate.minimum_query_interval_seconds if model else None,additive=False))
    return 0 if status.startswith('COMPLETE') else 1

if __name__=='__main__':raise SystemExit(main())
