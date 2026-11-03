"""Execute exactly one frozen auxiliary-only R5 replay; never opens truth."""
from pathlib import Path
import os
os.environ['OMP_NUM_THREADS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
os.environ['VECLIB_MAXIMUM_THREADS']='1'
import argparse
import json
import sys
import time
import traceback
from datetime import datetime, timezone

WALL, CPU = time.perf_counter(), time.process_time()
ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_max_struct_r5.src.trajectory_replay import (
    TrajectoryRecorder, compare_reference, dump)

parser=argparse.ArgumentParser()
parser.add_argument('--run-id',required=True)
args=parser.parse_args()
manifest=json.loads((ROOT/'configs/trajectory_replay_frozen.json').read_text())
cfg=next(row for row in manifest['runs'] if row['run_id']==args.run_id)
out=ROOT/'RESULTS/trajectory_replay'/cfg['run_id']
out.mkdir(parents=True,exist_ok=False)
dump(out/'FROZEN_CONFIG.json',cfg)
receipt=dict(status='RUNNING', run_id=cfg['run_id'], configuration=cfg,
    started_utc=datetime.now(timezone.utc).isoformat(), executable=sys.executable,
    seed=None, seed_reason='Historical inverse solve is deterministic; no random draws; public observations reused unchanged.',
    truth_access='NONE', corrector='NOT_RUN', parameter_tuning=False,
    failed_attempts=[], threads=1, device='CPU', callbacks=None)
dump(out/'RECEIPT.json',receipt)
model=recorder=None
try:
    import numpy as np
    from threadpoolctl import threadpool_limits, threadpool_info
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from research.research_max_struct_r2.src.cie_global import run_cie
    from research.research_max_struct_r4.src.cie_onehop import run_onehop_cie
    with threadpool_limits(limits=1):
        receipt['threadpools']=threadpool_info()
        tick,ct=time.perf_counter(),time.process_time()
        public={}
        for name in cfg['public_files']:
            with np.load(REPO/name,allow_pickle=False) as data:
                if any(key.lower() in {'truth','alpha','mask','yclean','seed','family'} for key in data.files):
                    raise ValueError('Unexpected private payload key')
                public.update({key:data[key] for key in data.files})
        model=MaxwellModel(public)
        Y=public['Y'];m0=np.full(model.n,cfg['numerical']['host'])
        recorder=TrajectoryRecorder(out,model,m0)
        receipt['load_setup_wall_seconds']=time.perf_counter()-tick
        receipt['load_setup_cpu_seconds']=time.process_time()-ct
        numerical={key:value for key,value in cfg['numerical'].items() if key != 'host'}
        numerical.update(method='cie',regularizer=lambda m:model.reg.evaluate(m)[:2],
            physical_sigma=float(np.linalg.norm(Y)),iteration_callback=recorder.iteration,
            stage_callback=recorder.stage,receipt_path=str(out/'AUXILIARY_EXECUTION.json'))
        tick,ct=time.perf_counter(),time.process_time()
        auxiliary=(run_cie if cfg['method']=='global' else run_onehop_cie)(model,Y,m0,numerical)
        receipt['auxiliary_wall_seconds']=time.perf_counter()-tick
        receipt['auxiliary_cpu_seconds']=time.process_time()-ct
        dump(out/'HISTORY.json',auxiliary['history'])
        dump(out/'AUXILIARY_COST.json',auxiliary['cost'])
        dump(out/'REFERENCE_COMPARISON.json',compare_reference(REPO,cfg,auxiliary))
        np.savez_compressed(out/'AUXILIARY_CANDIDATES.npz',
            m=np.array([c['m'] for c in auxiliary['candidates']]),
            objective=np.array([c['objective'] for c in auxiliary['candidates']]),
            physical_relative_residual=np.array([c['physical_relative_residual'] for c in auxiliary['candidates']]))
        receipt.update(status='COMPLETED_AUX_ONLY', auxiliary_status=auxiliary['status'],
            iterations=sum('iteration' in row for row in auxiliary['history']),
            auxiliary_execution_receipt=auxiliary['execution_receipt'],
            original_physics_validations=[{key:value for key,value in row.items() if key!='m'} for row in auxiliary['candidates']],
            reference_comparison='REFERENCE_COMPARISON.json')
except Exception as exc:
    receipt.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc())
    receipt['failed_attempts'].append(dict(error=repr(exc),charged=True))
    if hasattr(exc,'cie_execution_receipt'):
        receipt['failed_auxiliary_execution']=exc.cie_execution_receipt
finally:
    receipt.update(finished_utc=datetime.now(timezone.utc).isoformat(),
        wall_seconds=time.perf_counter()-WALL,cpu_seconds=time.process_time()-CPU,
        physical_model_cost=dict(model.cost) if model else {},
        physical_operation_wall=dict(model.timings) if model else {},
        physical_operation_cpu=dict(model._physics.cpu_timings) if model else {},
        callbacks=recorder.receipt() if recorder else None,
        cost_scope='Process imports/load/setup/unchanged aux optimizer+its original initial and stage physical validations/callbacks/output serialization. Final receipt write excluded. No corrector, scoring, truth reads, random draws, or formula tests.')
    dump(out/'RECEIPT.json',receipt)
print(json.dumps(dict(run_id=args.run_id,status=receipt['status'],wall_seconds=receipt['wall_seconds'],cpu_seconds=receipt['cpu_seconds'])),flush=True)
if receipt['status']=='FAILED':
    raise SystemExit(1)
