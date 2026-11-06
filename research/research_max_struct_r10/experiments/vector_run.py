"""Finite R10 vector probe/generation/inversion; each mode is a new process.

Inverse mode reads public observation/physics files only; material generation
is a separate process and lazy import. No scorer is called by this runner.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[name]='1'
import time
ENTRY=time.perf_counter()
from pathlib import Path
import sys,json,argparse,traceback,datetime,importlib.util
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))


def dump(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');tmp.replace(path)


def stamp(path):
    path=Path(path);s=path.stat()
    return dict(path=str(path.relative_to(REPO)),bytes=s.st_size,mtime_ns=s.st_mtime_ns)


def configure():
    if os.name!='nt':
        return dict(process=dict(pid=os.getpid(),platform=sys.platform,BLAS_threads=1),job=dict(active=False,scope='Local CPU fallback only')),None
    from research.research_max_struct_r10.src.vector_resources import configure_own_windows_process
    priority=configure_own_windows_process()
    if os.environ.get('R10_PARENT_JOB_ACTIVE')=='1':return dict(process=priority,job=dict(active=True,inherited_connected_parent=True)),None
    path=REPO/'research/research_max_struct_r9/delegated/3d_extension/windows_resolved.py'
    spec=importlib.util.spec_from_file_location('own_job_helper',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return dict(process=priority,job=module.job_object()),module


def initialize(config, output, side):
    import torch
    from research.research_max_struct_r10.src.vector_grid import VectorGrid
    from research.research_max_struct_r10.src.vector_cuda_model import ResidentVectorModel
    from research.research_max_struct_r10.src.vector_resources import CooperativeGate
    gate=CooperativeGate(output/'GPU_BOUNDARIES.jsonl',memory_limit_bytes=config['execution']['memory_limit_bytes'])
    if config.get('execution_backend','cuda').startswith('cuda'):gate.wait('before_any_model_or_CUDA')
    cpu=VectorGrid(side,config['physics']['width'],config['acquisition'])
    if config.get('execution_backend','cuda') in ('cpu','cpu-dense'):
        from research.research_max_struct_r10.src.vector_cpu_model import CpuVectorModel,CpuDenseVectorModel
        model=CpuDenseVectorModel(cpu) if config['execution_backend']=='cpu-dense' and side==12 else CpuVectorModel(cpu)
    elif config['execution_backend']=='cuda-dense':
        from research.research_max_struct_r10.src.vector_cuda_dense import DenseResidentVectorModel
        model=DenseResidentVectorModel(cpu,gate)
    else:model=ResidentVectorModel(cpu,gate)
    model.progress_file=output/'RHS_PROGRESS.jsonl'
    dump(output/'SHARED_SETUP.json',dict(setup=model.setup,CPU_counts=dict(cpu.counts),CUDA_counts=dict(model.counts),
        nested_in='PROCESS_COST.json',additive=False,all_source_vector_fields=True,source_shape=list(model.inc.shape),receiver_shape=list(model.S.shape)))
    print(json.dumps(dict(event='SHARED_MODEL_READY',side=side,vector_dimension=model.M,
        allocated_bytes=torch.cuda.memory_allocated() if model.device.type=='cuda' else 0,reserved_bytes=torch.cuda.memory_reserved() if model.device.type=='cuda' else 0)),flush=True)
    return model


def probe(config, output):
    import torch
    model=initialize(config,output,12)
    # Necessary working-set observation only: no material, current solve, GT or identity replay.
    receipt=dict(status='MODEL_READY_NO_SCIENTIFIC_SOLVE',
        CPU_setup=model.setup,allocated_bytes=torch.cuda.memory_allocated() if model.device.type=='cuda' else 0,reserved_bytes=torch.cuda.memory_reserved() if model.device.type=='cuda' else 0,
        estimated_additional_PR_and_Krylov_bytes=350000000,limit_bytes=1500000000,
        full_source_shape=list(model.inc.shape),full_receiver_shape=list(model.S.shape),
        uncertainty='CUDA context/WDDM driver memory may be unavailable perprocess; runtime boundary gate checks tensor reserved bytes',
        new_primal_RHS=0,new_adjoint_RHS=0,private_data_generated=False)
    dump(output/'RESOURCE_RECEIPT.json',receipt);print(json.dumps(dict(event='RESOURCE_RECEIPT',**receipt)),flush=True)
    return model


def generate(config, output, design_file, ids, sides=(12,18)):
    import numpy as np
    import torch
    # This module is never imported by inverse mode.
    from research.research_max_struct_r10.src.vector_materials import material_from_design
    designs=json.loads(design_file.read_text());chosen={c['design_id']:c for c in designs['cases'] if c['design_id'] in ids}
    if set(chosen)!=set(ids):raise ValueError('Unrecognized frozen generation ID')
    from collections import Counter
    from types import SimpleNamespace
    models={};generation=[];public=[];allcounts=Counter();allrecords=[];helper_cpu=0.;helper_queries=0;helper_pause=0.;helper_cached=0
    for side in sides:
        model=initialize(config,output/('grid'+str(side)),side);models[side]=model
        start_event=end_event=None
        if model.device.type=='cuda':
            start_event=torch.cuda.Event(enable_timing=True);end_event=torch.cuda.Event(enable_timing=True);start_event.record()
        for ordinal,case_id in enumerate(ids):
            private=ROOT/'sealed_truth/vector'/case_id;private.mkdir(parents=True,exist_ok=True)
            x=material_from_design(chosen[case_id],model.cpu.points,model.cpu.h)
            chi=torch.as_tensor(x,device=model.device,dtype=torch.float64)
            J=model.solve(chi,torch.repeat_interleave(chi,3)[:,None]*model.inc,label='generation_'+case_id+'_grid'+str(side))
            Q=model.inc+model.green(J);clean=model.S@J
            private_file=private/('grid'+str(side)+'.npz')
            if private_file.exists():raise RuntimeError('Writeonce private generation')
            np.savez_compressed(private_file,material=x,J=J.cpu().numpy(),Q=Q.cpu().numpy(),Y_clean=clean.cpu().numpy())
            sigma=.03*float(torch.linalg.vector_norm(clean).item())/np.sqrt(clean.numel())
            seed=202611070000+int(case_id[1:])*100+side
            rng=np.random.default_rng(seed);noise=sigma/np.sqrt(2)*(rng.standard_normal((180,60))+1j*rng.standard_normal((180,60)))
            folder=ROOT/'data_public/vector'/(case_id+('_matched12' if side==12 else '_cross18to12'))
            folder.mkdir(parents=True,exist_ok=False)
            np.savez_compressed(folder/'observations.npz',Y=clean.cpu().numpy()+noise)
            meta=dict(model='iid_circular_complex',sigma=sigma,expected_relative_Frobenius=.03,
                convention='Eabs(n)^2=sigma^2',whitening=False,seed=seed,realization_renormalized=False,
                public_id=folder.name,source_positions_resolution_independent=True,
                inverse_side=12,generation_side=side,truth_access_required_for_inverse=False)
            dump(folder/'PUBLIC.json',meta);public.append(stamp(folder/'PUBLIC.json'))
            generation.append(dict(case_id=case_id,side=side,private=stamp(private_file),public=stamp(folder/'observations.npz'),
                counts=dict(model.counts),private_fields_not_inverse_inputs=True))
            print(json.dumps(dict(event='PUBLIC_OBSERVATION_GENERATED',public_id=folder.name,side=side,
                sigma=sigma,allocated_bytes=torch.cuda.memory_allocated() if model.device.type=='cuda' else 0,reserved_bytes=torch.cuda.memory_reserved() if model.device.type=='cuda' else 0)),flush=True)
        device_ms=None
        if start_event is not None:
            end_event.record();end_event.synchronize();device_ms=start_event.elapsed_time(end_event)
        dump(output/('GRID'+str(side)+'_PROCESS_VIEW.json'),dict(counts=dict(model.counts),solve_records=model.solve_records,
            CUDA_event_stream_interval_ms=device_ms,CUDA_event_scope='Default stream elapsed including host/IO/pause gaps; not active kernel time',
            additive=False,included_in='PROCESS_COST.json'))
        allcounts.update(model.counts);allrecords.extend(model.solve_records)
        helper_cpu+=model.gate.helper_CPU_seconds;helper_queries+=model.gate.query_count
        helper_pause+=getattr(model.gate,'paused_wall_seconds',0.);helper_cached+=getattr(model.gate,'cached_boundary_checks',0)
        # Grid18 reuses no grid12 tensors; release before next constructor.
        del models[side];del model,J,Q,clean,chi
        if torch.cuda.is_initialized():torch.cuda.empty_cache()
    dump(output/'GENERATION_SEAL.json',dict(status='REQUESTED_OBSERVATIONS_GENERATED_WITH_PRIVATE_FIELDS_SEALED',
        entries=generation,public_metadata=public,inverse_GT_access=False,selection_by_score=False))
    return SimpleNamespace(counts=allcounts,solve_records=allrecords,
        gate=SimpleNamespace(helper_CPU_seconds=helper_cpu,query_count=helper_queries,unmeasured_helper_invocations=0,
            paused_wall_seconds=helper_pause,cached_boundary_checks=helper_cached))


def save_physical(folder,label,record):
    import numpy as np
    payload=dict(material=record['m'],prediction=record['prediction'],J=record['J'],Q=record['Q'])
    if record.get('gradient') is not None:payload['gradient']=record['gradient']
    np.savez_compressed(folder/(label+'.npz'),**payload)
    dump(folder/(label+'.json'),dict(state_id=record['id'],Phi=record['Phi'],PG_inf=record.get('projected_gradient_inf'),
        first_owner=record['first_owner'],uses=record['uses'].copy()))


def invert(config, output, queue, maximum, offset=0):
    import numpy as np
    import torch
    from research.research_max_struct_r10.src.vector_cuda_model import ExactReceiverGeometry,NoiseReceiverSubspace
    from research.research_max_struct_r10.src.vector_cie import auxiliary_solve
    from research.research_max_struct_r10.src.vector_physical import PhysicalObjective,physical_corrector
    from research.research_max_struct_r10.src.vector_energy import current_receiver_energy,affine_representation_energy,scalar_material_energy
    tasks=queue['tasks'][offset:offset+maximum]
    if not tasks or offset+maximum>queue['authorized_task_count']:raise ValueError('Queue subset exceeds frozen authorization')
    model=initialize(config,output,12);geometry=ExactReceiverGeometry(model)
    dump(output/'RECEIVER_EXACT_SVD.json',geometry.info)
    settings=dict(mf=config['inverse']['native_MF'],beta_material=config['inverse']['beta_material'],
        beta_current_dual=config['inverse']['beta_current_dual'],max_iterations_per_stage=500,
        minimum_iterations_per_stage=20,coefficient_delta_RMS_tolerance=.001,
        initial_chi=.25,material_upper=8.)
    current_id=None;terminals=[];case_views=[]
    for task_index,task in enumerate(tasks):
        model.gate.wait('before_task_'+task['task_id'],torch.cuda.memory_reserved() if model.device.type=='cuda' else 0)
        if task['public_id']!=current_id:
            current_id=task['public_id'];public=ROOT/'data_public/vector'/current_id
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
        dest=ROOT/'RESULTS/vector'/task['task_id'];dest.mkdir(parents=True,exist_ok=False)
        wall,cpu=time.perf_counter(),time.process_time();first_event=len(objective.events);before=dict(model.counts)
        device_start=device_end=None;pause_start=getattr(model.gate,'paused_wall_seconds',0.)
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
        result=auxiliary_solve(model,receiver,model.array(Y),task['arm'],settings,callback,endpoint)
        last=result['material'].cpu().numpy()
        best,corrector=physical_corrector(objective,last,config['inverse']['corrector'],task['task_id'],
            lambda label,record:save_physical(dest,label,record),emit)
        save_physical(dest,'prediction',best)
        device_ms=None
        if device_start is not None:
            device_end.record();device_end.synchronize();device_ms=device_start.elapsed_time(device_end)
        counts={k:int(v-before.get(k,0)) for k,v in model.counts.items() if v-before.get(k,0)}
        receipt=dict(status='SEALED_PUBLIC_PREDICTION_NO_SCORE',task=task,stages=result['stages'],physical_corrector=corrector,
            physical_Phi=best['Phi'],PG_inf=best['projected_gradient_inf'],counts=counts,
            CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
            CUDA_event_stream_interval_ms=device_ms,
            CUDA_event_scope='Default-stream timestamp interval across this task; includes host launch/IO/idle and cooperative pause gaps; NOT active-kernel time or SM occupancy',
            cooperative_pause_wall_seconds=getattr(model.gate,'paused_wall_seconds',0.)-pause_start,
            GPU_peak_allocated_bytes_to_date=torch.cuda.max_memory_allocated() if model.device.type=='cuda' else 0,
            GPU_peak_reserved_bytes_to_date=torch.cuda.max_memory_reserved() if model.device.type=='cuda' else 0,
            additive=False,included_in='PROCESS_COST.json',paid_events=objective.events[first_event:],
            current_metric='ordinary_CIE',current_aggregation='sourcewise' if task['arm']=='sourceCIE' else 'global',
            current_beta=settings['beta_material'] if task['arm']=='sourceCIE' else settings['beta_current_dual'],
            material_beta=settings['beta_material'],material_metric='Aplus' if task['arm']=='AplusGlobal' else 'ordinary_CIE',
            comparison_scope='sourceCIE differs in current beta AND sourcewise PR; GlobalDual matches AplusGlobal current exactly',
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
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['probe','generate','invert'])
    parser.add_argument('--backend',choices=['cuda','cuda-dense','cpu','cpu-dense'],default='cuda')
    parser.add_argument('--config',default='configs/vector_actual_contract.json');parser.add_argument('--output',required=True)
    parser.add_argument('--generation-sides',type=int,nargs='+',default=[12,18])
    parser.add_argument('--design',default='configs/vector_material_design.json');parser.add_argument('--ids',nargs='*',default=['v02','v09'])
    parser.add_argument('--task-offset',type=int,default=0)
    parser.add_argument('--queue',default='configs/vector_first_queue.json');parser.add_argument('--max-tasks',type=int,default=3)
    args=parser.parse_args();output=ROOT/args.output;output.mkdir(parents=True,exist_ok=False)
    model=None;status='FAILED';failure=None;process=None
    try:
        process_config,job_module=configure()
        import torch
        if os.name=='nt':
            import psutil
            process=psutil.Process()
        else:
            import resource
            from types import SimpleNamespace
            process=SimpleNamespace(cpu_times=lambda:SimpleNamespace(user=resource.getrusage(resource.RUSAGE_SELF).ru_utime,system=resource.getrusage(resource.RUSAGE_SELF).ru_stime),memory_info=lambda:SimpleNamespace(peak_wset=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        torch.set_num_threads(1);torch.set_num_interop_threads(1)
        config_path=ROOT/args.config;config=json.loads(config_path.read_text());config['execution_backend']=args.backend
        dump(output/'FROZEN_CONFIG.json',config);dump(output/'PROCESS_CONFIGURATION.json',process_config)
        dump(output/'SOURCE_EVIDENCE.json',[stamp(p) for p in sorted((ROOT/'src').glob('vector_*.py'))]+[stamp(Path(__file__)),stamp(config_path)])
        if args.mode=='probe':model=probe(config,output)
        elif args.mode=='generate':model=generate(config,output,ROOT/args.design,args.ids,args.generation_sides)
        else:model=invert(config,output,json.loads((ROOT/args.queue).read_text()),args.max_tasks,args.task_offset)
        status='COMPLETE'
    except BaseException as exc:
        failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());print(json.dumps(failure),flush=True)
    finally:
        if process:
            t=process.cpu_times();cpu=t.user+t.system;peak=getattr(process.memory_info(),'peak_wset',None)
        else:cpu=peak=None
        cost=dict(status=status,mode=args.mode,platform='XINAN' if os.name=='nt' else 'Mac',CPU_seconds=cpu,wall_seconds=time.perf_counter()-ENTRY,
            peak_RSS_bytes=peak,additive=True,included_in=None,scope='Wholeprocess imports/sharedsetup/gates/all requested work/IO; nested views not additive',
            failure=failure,counts=dict(model.counts) if model else None,
            peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated() if 'torch' in locals() and torch.cuda.is_initialized() else None,
            peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved() if 'torch' in locals() and torch.cuda.is_initialized() else None,
            GT_access=args.mode=='generate',GT_access_scope='Private generation only; inverse mode has no design/GT file read',
            actual_solve_records=model.solve_records if model else None,
            GPU_gate_helper_CPU_seconds=model.gate.helper_CPU_seconds if model else None,
            GPU_gate_helper_invocations=model.gate.query_count if model else None,
            GPU_gate_cached_boundary_checks=getattr(model.gate,'cached_boundary_checks',0) if model else None,
            GPU_gate_paused_wall_seconds=getattr(model.gate,'paused_wall_seconds',0.) if model else None,
            GPU_gate_helper_CPU_scope='Separate nvidia-smi children excluded from main psutil CPU; add once, not all snapshots again')
        dump(output/'PROCESS_COST.json',cost)
    return 0 if status=='COMPLETE' else 1


if __name__=='__main__':raise SystemExit(main())
