"""Execute exactly one public R6 case/arm task; no batch or GT access.

CLI: --manifest FILE --task-id ID [--config FILE] [--output-dir DIR]
Alternative --public-list FILE accepts the same public task-list JSON schema.
The parent freezes/creates the manifest. This module has no execution on import.
"""
from pathlib import Path
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[_key]='1'
import argparse
from datetime import datetime, timezone
import json
import sys
import time
import traceback

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
if str(REPO) not in sys.path:sys.path.insert(0,str(REPO))


def dump(path, value):
    import numpy as np
    def convert(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False,default=convert)+'\n')


def resolve(path, base=REPO):
    path=Path(path)
    return path if path.is_absolute() else base/path


def load_public(files):
    import numpy as np
    public={}
    forbidden={'truth','alpha','mask','yclean','y_clean','seed','family','gt','ground_truth',
               'm_true','m_truth','true_material'}
    for filename in files:
        path=resolve(filename)
        if any(part.lower() in ('data_private','private','truth','gt','sealed_gt') for part in path.parts):
            raise ValueError('Private path is not a public input: '+str(path))
        with np.load(path,allow_pickle=False) as data:
            if forbidden & {key.lower() for key in data.files}:
                raise ValueError('Unexpected private payload key in public file: '+str(path))
            for key in data.files:
                if key in public and not np.array_equal(public[key],data[key]):
                    raise ValueError('Conflicting public value: '+key)
                public[key]=data[key]
    if 'E0' not in public and 'E' in public:public['E0']=public['E']
    return public


def choose_task(manifest, task_id):
    tasks=manifest if isinstance(manifest,list) else manifest.get('tasks',[])
    if isinstance(tasks,dict):
        if task_id not in tasks:raise ValueError('Task ID outside public manifest')
        return dict(tasks[task_id],task_id=task_id)
    found=[t for t in tasks if t.get('task_id')==task_id]
    if len(found)!=1:raise ValueError('Task ID must occur exactly once in public manifest')
    return found[0]


def execute(args):
    wall,cpu=time.perf_counter(),time.process_time()
    manifest_path=resolve(args.manifest or args.public_list)
    manifest=json.loads(manifest_path.read_text())
    task=choose_task(manifest,args.task_id)
    config_ref=args.config or task.get('config') or (manifest.get('config') if isinstance(manifest,dict) else None)
    if not config_ref:raise ValueError('An explicit frozen config is required')
    config_path=resolve(config_ref)
    cfg=json.loads(config_path.read_text())
    arm=task['arm'];case_id=task['case_id']
    files=task.get('public_files')
    if not files or not isinstance(files,list):raise ValueError('Task must declare public_files list')
    out=resolve(args.output_dir or task.get('output_dir') or
                str(ROOT/'results'/cfg['run_id']/args.task_id))
    out.mkdir(parents=True,exist_ok=False)
    receipt=dict(status='RUNNING',task_id=args.task_id,case_id=case_id,arm=arm,
        started_utc=datetime.now(timezone.utc).isoformat(),executable=sys.executable,
        manifest_path=str(manifest_path),config_path=str(config_path),task=task,configuration=cfg,
        truth_access='NONE',scoring='NOT_RUN',seed=None,
        seed_reason='No random draws; deterministic solve of public observations',
        failed_attempts=[],device='CPU',threads_per_process=1,
        unknown_costs=None)
    dump(out/'FROZEN_CONFIG.json',cfg);dump(out/'PUBLIC_TASK.json',task)
    dump(out/'RECEIPT.json',receipt)
    model=observer=None
    iteration_observations=[]
    try:
        import numpy as np
        from threadpoolctl import threadpool_limits,threadpool_info
        from research.research_fbg_r3_experiments.src.physics import MaxwellModel
        from research.research_max_struct_r6.src.benchmark_solver import (
            ARMS,SOURCE_PATHS,arm_config,run_auxiliary,PaidStateObserver,refine_observed,difference)
        if arm not in ARMS:raise ValueError('Unknown fixed R6 arm')
        with threadpool_limits(limits=1):
            receipt['threads']=threadpool_info()
            setup_wall,setup_cpu=time.perf_counter(),time.process_time()
            public=load_public(files)
            for key,value in cfg.get('expected_public',{}).items():
                if not np.array_equal(np.asarray(public[key]),np.asarray(value)):
                    raise ValueError('Frozen public mismatch: '+key)
            model=MaxwellModel(public);Y=np.asarray(public['Y'],complex)
            if model.E.shape!=(576,12) or model.S.shape!=(48,576) or Y.shape!=(48,12):
                raise ValueError('R6 fixed benchmark requires 576 real cells, 12 sources, 48 receivers')
            if tuple(model.shape)!=(24,24):raise ValueError('R6 fixed grid requires 24x24')
            m0=np.full(model.n,cfg['common']['host'],float)
            numerical=arm_config(arm,cfg['common'],Y,model.E,cfg.get('continuation_stages'))
            np.savez_compressed(out/'SOURCE_METRICS.npz',
                data_weights=numerical['data_weights'],state_weights=numerical['state_weights'])
            dump(out/'NUMERICAL_SCHEDULE.json',{key:value for key,value in numerical.items()
                                              if key not in ('data_weights','state_weights')})
            receipt.update(setup_wall_seconds=time.perf_counter()-setup_wall,
                setup_cpu_seconds=time.process_time()-setup_cpu,numerical_sources=SOURCE_PATHS,
                public_files=[str(resolve(f)) for f in files],initialization=
                'm=shared .1 host; J=DPCS pseudoinverse current; APCS A=0 in original Fourier coordinates')
            def snapshot(name,row):
                arrays={key:row[key] for key in ('m','j','Q','Y_prediction','gradient') if key in row}
                np.savez_compressed(out/(name+'.npz'),**arrays)
                dump(out/(name+'.json'),{key:value for key,value in row.items() if key not in arrays})
            observer=PaidStateObserver(model,Y,snapshot,
                cfg.get('corrector_evaluation_checkpoints',[100,500,1000]))
            aux_wall,aux_cpu=time.perf_counter(),time.process_time()
            before=dict(model.cost)
            def iteration(row):
                iteration_observations.append(dict(stage=row['stage'],iteration=row['iteration'],
                    cpu_seconds=time.process_time()-aux_cpu,wall_seconds=time.perf_counter()-aux_wall,
                    physical_model_cost=dict(model.cost),operation_counts=row.get('operation_counts')))
            def stage(row):
                index=row['stage']
                np.savez_compressed(out/f'stage_{index:02d}_joint.npz',
                    m=row['m'],J=row['J'],Q=row['Q'],A=row['coefficients'])
                dump(out/f'stage_{index:02d}_joint.json',
                    dict({key:value for key,value in row.items() if key not in ('m','J','Q','coefficients')},
                         cpu_seconds=time.process_time()-aux_cpu,
                         wall_seconds=time.perf_counter()-aux_wall,
                         physical_model_cost=dict(model.cost)))
            numerical.update(regularizer=lambda m:model.reg.evaluate(m)[:2],
                physical_sigma=float(np.linalg.norm(Y)),iteration_callback=iteration,stage_callback=stage,
                receipt_path=str(out/'AUXILIARY_EXECUTION.json'))
            auxiliary=run_auxiliary(observer,Y,m0,arm,numerical)
            receipt.update(auxiliary_wall_seconds=time.perf_counter()-aux_wall,
                auxiliary_cpu_seconds=time.process_time()-aux_cpu,
                auxiliary_physical_model_cost=difference(dict(model.cost),before),
                auxiliary_execution_receipt=auxiliary['execution_receipt'],auxiliary_status=auxiliary['status'])
            history=auxiliary['history']
            times={(r['stage'],r['iteration']):r for r in iteration_observations}
            for row in history:
                if 'iteration' in row:row['observed_cost']=times[(row['stage'],row['iteration'])]
            dump(out/'HISTORY.json',history);dump(out/'AUXILIARY_COST.json',auxiliary['cost'])
            if len(observer.validation_records)!=len(auxiliary['candidates']):
                raise RuntimeError('Paid physical validations do not align with returned candidates')
            for a,b in zip(observer.validation_records,auxiliary['candidates']):
                if not np.array_equal(a['m'],b['m']) or not np.isclose(a['objective'],b['objective'],rtol=1e-13,atol=1e-15):
                    raise RuntimeError('Recorded candidate differs from original physical validation')
            selected=min(auxiliary['candidates'],key=lambda row:row['objective'])
            receipt['auxiliary_candidate_scores']=[{key:value for key,value in row.items() if key!='m'}
                                                   for row in auxiliary['candidates']]
            receipt['auxiliary_selected_stage']=selected['stage']
            snapshot('selected_before_corrector',min(observer.validation_records,key=lambda row:row['objective']))
            np.savez_compressed(out/'AUXILIARY_CANDIDATES.npz',
                m=np.array([row['m'] for row in auxiliary['candidates']]),
                original_Phi=np.array([row['objective'] for row in auxiliary['candidates']]))
            maximum=cfg['corrector_maxiter']
            if isinstance(maximum,bool) or not isinstance(maximum,int) or maximum<=0:
                raise ValueError('Original physical corrector allowance must be an explicit positive integer')
            refined,corrector=refine_observed(observer,Y,selected['m'],maximum=maximum,
                accepted_checkpoints=cfg.get('corrector_accepted_checkpoints',[100,500,1000]))
            dump(out/'CORRECTOR.json',corrector)
            receipt['corrector']={key:value for key,value in corrector.items()
                                  if key not in ('history','accepted_history','evaluation_history')}
            # The observer retains the original-Phi incumbent from charged calls.
            # No finalization forward is needed to export its exact prediction.
            best=observer.best
            if best is None:raise RuntimeError('No paid physical incumbent was retained')
            snapshot('prediction',best)
            receipt.update(status='SEALED_PUBLIC_PREDICTION',physical_objective=best['objective'],
                physical_relative_residual=best['physical_relative_residual'],
                linear_relative_residual=best['state_relative_residual'],
                auxiliary_iterations=sum('iteration' in row for row in history),
                corrector_improved_original_Phi=bool(best['objective']<selected['objective']),
                finalization_physical_model_cost={},
                selection='Minimum original Phi of host + all true stage endpoints + original R4 corrector trajectory incumbent',
                snapshot_policy='Initial and stage physical states; stage joint m/J/Q/A; configured paid corrector points; final paid incumbent',
                continuation='Explicit appended mf12 stages in same invocation; original Fourier transport, no J reset')
    except Exception as exc:
        receipt.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc())
        receipt['failed_attempts'].append(dict(error=repr(exc),charged=True))
        if hasattr(exc,'cie_execution_receipt'):receipt['failed_auxiliary_execution']=exc.cie_execution_receipt
        if iteration_observations:dump(out/'PARTIAL_ITERATION_COST.json',iteration_observations)
    finally:
        receipt.update(finished_utc=datetime.now(timezone.utc).isoformat(),
            wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
            physical_model_cost=dict(model.cost) if model is not None else None,
            physical_operation_wall=dict(model.timings) if model is not None else None,
            physical_operation_cpu=dict(model._physics.cpu_timings) if model is not None else None,
            cost_scope='Process imports/public load/setup, complete configured auxiliary schedule, original host/stage physical validation, original R4 corrector and observation/serialization. Nested costs overlap; do not add twice.')
        dump(out/'RECEIPT.json',receipt)
    return receipt,out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    source=parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--manifest');source.add_argument('--public-list')
    parser.add_argument('--task-id',required=True)
    parser.add_argument('--config');parser.add_argument('--output-dir')
    args=parser.parse_args()
    receipt,out=execute(args)
    print(json.dumps(dict(task_id=args.task_id,status=receipt['status'],output_dir=str(out),
        wall_seconds=receipt['wall_seconds'],cpu_seconds=receipt['cpu_seconds'])),flush=True)
    return 1 if receipt['status']=='FAILED' else 0


if __name__=='__main__':raise SystemExit(main())
