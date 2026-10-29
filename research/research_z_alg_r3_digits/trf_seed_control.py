#!/usr/bin/env python3
"""Post-GT matched public coarse-seed TRF control; original coarse cost is reference only."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,json,time,shutil
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from threadpoolctl import threadpool_limits
from src.solver import load_public,BudgetStop,StageStop,json_safe
from trf_diagnostic import TRFRun
ROOT=Path(__file__).resolve().parent
CONFIG=ROOT/'audit/FULLJ_DIAGNOSTIC_MATCHED_SEED_CONFIG.json'
SNAP=ROOT/'audit/FULLJ_DIAGNOSTIC_MATCHED_SEED_SOURCE_SNAPSHOT'

def write(path,obj):
 p=Path(path)
 if p.exists():raise FileExistsError(str(p))
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(json_safe(obj),indent=2,allow_nan=False)+'\n')

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--freeze',action='store_true');ap.add_argument('--public');ap.add_argument('--operators');ap.add_argument('--output');a=ap.parse_args()
 if a.freeze:
  if SNAP.exists():raise FileExistsError(str(SNAP))
  SNAP.mkdir(parents=True);rows=[]
  for name in ('trf_seed_control.py','trf_diagnostic.py','audit/FULLJ_DIAGNOSTIC_MATCHED_SEED_CONFIG.json','src/solver.py','src/physics.py','src/multiscale.py','src/__init__.py'):
   p=SNAP/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,p);rows.append(dict(path=name,bytes=p.stat().st_size))
  write(SNAP/'FREEZE.json',dict(frozen_utc=datetime.now(timezone.utc).isoformat(),files=rows,status='POST_GT_MATCHED_PUBLIC_SEED_NOT_BLIND',truth_access='NONE',hashes='NOT_COMPUTED'));return
 manifest=json.loads((SNAP/'FREEZE.json').read_text())
 for row in manifest['files']:
  if (ROOT/row['path']).read_bytes()!=(SNAP/row['path']).read_bytes():raise RuntimeError('Frozen source changed')
 contract=json.loads(CONFIG.read_text());originalpath=ROOT/contract['original_path'];original=json.loads(originalpath.read_text())
 st=next(s for s in original['stages'] if s['stage']==24);seed=np.asarray(st['initial_material'],float);coarsewall=float(st['started_seconds']);remaining=max(0.,90-coarsewall)
 begun,cpu=time.perf_counter(),time.process_time();cfg=dict(original['config']);cfg.update(wall_seconds=remaining,common_coarse_native=False,max_gn_iterations=1000000000)
 with threadpool_limits(limits=1):
  run=TRFRun(load_public(a.public,a.operators),cfg);run.started,run.cpu_started=begun,cpu;run.active_stage=24;run.stage_deadline=remaining
  status='OPTIMIZER_FINISHED';error=None
  try:run.gn(seed,np.eye(576),remaining,origin='MATCHED_PUBLIC_COARSE_SEED_TRF')
  except (BudgetStop,StageStop):status='REMAINING_BUDGET_STOP';run.status='BUDGET_STOP';run.stop_reason='ALLOCATED_FINE_REMAINING'
  except BaseException as exc:status='FAILED_PAID';error=repr(exc);run.failures[type(exc).__name__]+=1
  run.events.append(dict(kind='MULTISCALE_STAGE',stage=24,parameter_dimension=576,fixed_physics_dimension=576,
    started_seconds=0.,wall_seconds=run.elapsed(),status=status,initial_material=seed,
    counts=dict(run.model.counts),has_evaluated_candidate=bool(run.archive),optimizer='matchedseedSciPyTRF'))
  result=run.finish()
 result.update(method='trf_matched_public_seed',development_status='POST_GT_MATCHED_PUBLIC_COARSE_SEED_NOT_BLIND',error=error,
   original_seed_input_path=str(originalpath),original_seed_input_bytes=originalpath.stat().st_size,
   seed_provenance='EXACT stored original stage24 initial_material; no truth or label',
   original_coarse_wall_reference_only=coarsewall,original_coarse_counts_reference_only=st['coarse_seed_counts'],
   allocated_fine_budget_seconds=remaining,allocated_total_wall_seconds=coarsewall+result['wall_seconds'],
   cost_aggregation='ROOT ACCOUNTING IS NEW FINE TRF ONLY; original coarse reference must not be added again to total study costs',
   source_freeze_utc=manifest['frozen_utc'])
 write(a.output,result);print(json.dumps(dict(output=a.output,status=result['status'],wall_seconds=result['wall_seconds'],allocated_total_wall_seconds=result['allocated_total_wall_seconds'],residual=(result.get('top1') or {}).get('relative_residual'),accounting=result['accounting'])),flush=True)
if __name__=='__main__':main()
