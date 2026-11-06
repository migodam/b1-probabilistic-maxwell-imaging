"""Seal80 paid predictions and declared portfolios before any confirmation GT read."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,time,sys
import numpy as np
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];OUT=ROOT/'delegated/diagnostics/bothblocks_confirmation'

def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def main():
 wall,cpu=time.perf_counter(),time.process_time();q=json.loads((OUT/'QUEUE.json').read_text());tasks=json.loads((ROOT/'configs/twod_bothblocks16_all_tasks.json').read_text())['tasks'];rule=json.loads((OUT/'PORTFOLIO_RULE_BEFORE_DATA_SCORE.json').read_text())
 if q['status'] not in ('COMPLETED','COMPLETED_WITH_ERRORS') or q.get('pending') or len(q['jobs'])!=80:raise RuntimeError('All80 terminal required; no truth opened')
 if len({j['task_id'] for j in q['jobs']})!=80 or any(j['status']=='RUNNING' for j in q['jobs']):raise RuntimeError('Ambiguous terminal outcomes')
 jobs={j['task_id']:j for j in q['jobs']};entries=[];lookup={}
 for t in tasks:
  if t['task_id'] not in jobs:raise RuntimeError('Manifest/queue mismatch')
  d=REPO/t['output_dir'];rp=d/'RECEIPT.json';r=json.loads(rp.read_text()) if rp.exists() else {}
  if r.get('status')=='RUNNING':raise RuntimeError('Running task remains')
  pred=d/'R11_CAP_PREDICTION.npz';meta=d/'R11_CAP_PREDICTION.json'
  available=bool(pred.exists() and meta.exists() and r.get('cap_prediction_status')=='SEALED_PAID_INCUMBENT')
  row=dict(task=t,receipt_status=r.get('status','MISSING_RECEIPT'),terminal_job=jobs[t['task_id']],cap_prediction_available=available,prediction_path=str(pred),prediction_meta=str(meta),
   cap_original_Phi=r.get('cap_physical_objective'),truth_access_in_solver=r.get('truth_access'),prediction_bytes=pred.stat().st_size if pred.exists() else None,prediction_mtime_ns=pred.stat().st_mtime_ns if pred.exists() else None)
  if r.get('truth_access') not in ('NONE',None):raise RuntimeError('Unexpected online truth access')
  entries.append(row);lookup[t['case_id'],t['method']]=row
 dump(OUT/'ALL80_TERMINAL_PREDICTION_SEAL.json',dict(status='ALL80_TERMINAL_BEFORE_ANY_CONFIRMATION_GT_SCORING',entries=entries,available_predictions=sum(e['cap_prediction_available'] for e in entries),truth_access='NONE_SEAL_PROCESS'))
 ports=[];mats=[]
 for cid in sorted({t['case_id'] for t in tasks}):
  for portfolio in rule['portfolios']:
   methods=portfolio['always_run'];parts=[lookup[cid,m] for m in methods]
   for view,cap in [('completed_480_each',480),('prefix240_each',240)]:
    candidates=[]
    for e in parts:
     t=e['task'];d=REPO/t['output_dir']
     if cap==480:
      if not e['cap_prediction_available']:continue
      with np.load(d/'R11_CAP_PREDICTION.npz') as z:m=z['m'].copy()
      candidates.append(dict(method=t['method'],Phi=e['cap_original_Phi'],m=m,source=str(d/'R11_CAP_PREDICTION.npz'),paid_CPU_seconds=None))
     else:
      cp=json.loads((d/'R11_CPU_CHECKPOINTS.json').read_text()) if (d/'R11_CPU_CHECKPOINTS.json').exists() else []
      state=next((p for p in cp if p['CPU_cap_seconds']==cap and 'material_index' in p),None)
      if state is None:continue
      with np.load(d/'R11_CHECKPOINT_MATERIALS.npz') as z:m=z['m'][state['material_index']].copy()
      candidates.append(dict(method=t['method'],Phi=state['original_Phi'],m=m,source=str(d/'R11_CHECKPOINT_MATERIALS.npz'),paid_CPU_seconds=state['state_paid_CPU_seconds']))
    actual=[e['terminal_job']['child_process_CPU_seconds'] for e in parts];natural=all(x<=cap for x in actual)
    row=dict(case_id=cid,portfolio=portfolio['name'],view=view,component_methods=methods,component_allowed_CPU_caps=[cap,cap],total_allowed_CPU_cap=2*cap,
      component_actual_complete_child_CPU_seconds=actual,complete_deployment_CPU_seconds=sum(actual) if cap==480 or natural else None,
      cost_status='EXACT_ACTUAL_SUM' if cap==480 or natural else 'SAVED_PREFIX_ALLOWED_CAP; EXACT_TRUNCATED_DEPLOYMENT_SPEND_NOT_OBSERVED',
      wall_status='NO_SINGLE_DEPLOYMENT_WALL_MEASURED; SHARED_QUEUE_WALL_NOT_USED',candidate_outcomes=[dict(method=e['task']['method'],receipt_status=e['receipt_status'],cap_prediction_available=e['cap_prediction_available']) for e in parts])
    if not candidates:row.update(status='NO_ELIGIBLE_PAID_STATE',material_index=None)
    else:
     # Input ordering CIE first plus strict Phi improvement implements declared tie policy.
     chosen=candidates[0]
     for c in candidates[1:]:
      if c['Phi']<chosen['Phi']:chosen=c
     row.update(status='SEALED_PUBLIC_PHI_SELECTION',selected_method=chosen['method'],original_Phi=chosen['Phi'],selected_source=chosen['source'],material_index=len(mats),eligible_candidates=len(candidates),selected_paid_CPU_seconds=chosen['paid_CPU_seconds']);mats.append(chosen['m'])
    ports.append(row)
 np.savez_compressed(OUT/'PORTFOLIO_SEALED_MATERIALS.npz',m=np.asarray(mats))
 dump(OUT/'PORTFOLIO_SELECTION_BEFORE_GT.json',dict(status='SEALED_BEFORE_ANY_CONFIRMATION_GT_READ',rule_file='PORTFOLIO_RULE_BEFORE_DATA_SCORE.json',rows=ports,selection_GT_access='NONE',extra_solver_runs=0))
 dump(OUT/'PREDICTION_SEAL_RECEIPT.json',dict(status='COMPLETE',terminal_tasks=80,predictions=sum(e['cap_prediction_available'] for e in entries),portfolio_predictions=len(mats),physical_forward=0,physical_adjoint=0,truth_access='NONE',CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall))
 print(json.dumps(dict(status='SEALED_BEFORE_GT',predictions=sum(e['cap_prediction_available'] for e in entries),portfolios=len(mats))))
if __name__=='__main__':main()
