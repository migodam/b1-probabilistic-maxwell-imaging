"""All64-terminal public prediction/prefix/actual-end-cost seal; no GT reads."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,time,itertools,shutil
import numpy as np
HERE=Path(__file__).absolute().parent;ROOT=HERE.parents[2];REPO=ROOT.parents[1]
def read(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def resolve(p):
 p=Path(str(p).replace('\\','/'));return p if p.is_absolute() else REPO/p
def arrays(p):
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def method(t):return t.get('material_readout',{}).get('label',t['arm'])
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--queue',type=Path,action='append',required=True);ap.add_argument('--output',type=Path,default=HERE/'public_seal');a=ap.parse_args();wall=time.perf_counter()
 tasks=read(HERE/'TASKS_V2.json')['tasks'];rule=read(HERE/'SCORING_CONFIG_V2.json');cfg=read(HERE/'SOLVER_CONFIG_V2.json');jobs={}
 if len(tasks)!=64 or len({t['task_id'] for t in tasks})!=64 or len({t['case_id'] for t in tasks})!=8:raise ValueError('Exactly64 original tasks/8 materials required')
 expected={(c,g,m) for c in {t['case_id'] for t in tasks} for g in rule['grid_conditions'] for m in rule['methods']}
 if {(t['case_id'],t['grid_condition'],method(t)) for t in tasks}!=expected:raise ValueError('Unbalanced paired64 contract')
 for path in a.queue:
  q=read(path)
  if q.get('status') not in ('COMPLETED','COMPLETED_WITH_ERRORS','STOPPED_DISPATCH_CHILDREN_FINISHED') or q.get('running'):raise RuntimeError('All dispatch phases must be terminal before public seal')
  if q['status']!='STOPPED_DISPATCH_CHILDREN_FINISHED' and q.get('pending'):raise RuntimeError('Pending tasks')
  for j in q['jobs']:
   if j['task_id'] in jobs or j['status'] not in ('FINISHED','FAILED','LAUNCH_FAILED'):raise ValueError('Duplicate/nonterminal original identity; attempts require separate accounting')
   jobs[j['task_id']]=dict(j,source_queue_receipt=str(path))
 if set(jobs)!={t['task_id'] for t in tasks}:raise RuntimeError('All64 terminal outcomes required; never seal a subset')
 a.output.mkdir(parents=True,exist_ok=False);dump(a.output/'PUBLIC_GATE.json',dict(status='ALL64_TERMINAL_BEFORE_ANY_SCORER_PRIVATE_READ',truth_access='NONE',tasks=64,independent_materials=8,created_utc=datetime.now(timezone.utc).isoformat()))
 mats=[];ypreds=[];rows=[];pairs=[];sources=[]
 def select(t,data,B,view,pair_id=None,side=None):
  j=jobs[t['task_id']];valid,receipt,trace,pred=data;chosen=None
  indices=np.flatnonzero((trace['CPU_seconds']<=B)&(trace['CPU_seconds']<=480)) if valid and trace is not None and B is not None else np.array([],int)
  if indices.size:
   k=int(indices[np.argmin(trace['original_Phi'][indices])]);chosen=dict(material_index=len(mats),paid_CPU_seconds=float(trace['CPU_seconds'][k]),original_Phi=float(trace['original_Phi'][k]),paid_LU=int(trace['LU'][k]),paid_complex_RHS=int(trace['complex_RHS'][k]));mats.append(trace['m'][k].copy())
  r=dict(case_id=t['case_id'],grid_condition=t['grid_condition'],method=method(t),task_id=t['task_id'],view=view,pair_id=pair_id,side=side,budget_B=B,selected=chosen,task_status=receipt.get('status',j['status']),actual_complete_child_CPU_seconds=j.get('child_process_CPU_seconds'),actual_physical_cost=receipt.get('physical_model_cost'),source_queue_receipt=j['source_queue_receipt'],public_files=t['public_files'])
  if view=='FINAL_CAP_ELIGIBLE' and valid:
   if chosen is None or chosen['original_Phi']!=receipt['cap_physical_objective'] or not np.array_equal(mats[chosen['material_index']],pred['m']):raise ValueError('Final cap export differs from public minPhi trajectory')
   r['Y_prediction_index']=len(ypreds);ypreds.append(pred['Y_prediction'].copy())
  rows.append(r);return r
 by={(t['case_id'],t['grid_condition'],method(t)):t for t in tasks}
 for cid in sorted({t['case_id'] for t in tasks}):
  loaded={}
  for condition in rule['grid_conditions']:
   for m in rule['methods']:
    t=by[cid,condition,m];d=resolve(t['output_dir']);j=jobs[t['task_id']];r=read(d/'RECEIPT.json') if (d/'RECEIPT.json').exists() else {}
    if r.get('status') in ('RUNNING','PREPARING','PREPARED'):raise RuntimeError('Nonterminal child receipt')
    valid=j['status']=='FINISHED' and r.get('status') in ('SEALED_PUBLIC_PREDICTION','CPU_BUDGET_REACHED') and r.get('cap_prediction_status')=='SEALED_PAID_INCUMBENT'
    trace=pred=None
    if valid:
     if (r.get('task_id'),r.get('case_id'),r.get('arm'))!=(t['task_id'],cid,t['arm']):raise ValueError('Paid identity mismatch')
     trace=arrays(d/'R11_INCUMBENT_TRAJECTORY.npz');pred=arrays(d/'R11_CAP_PREDICTION.npz')
     if trace['m'].ndim!=2 or trace['m'].shape[1]!=576 or not np.isfinite(trace['m']).all() or not np.isfinite(trace['CPU_seconds']).all() or not np.isfinite(trace['original_Phi']).all():raise ValueError('Malformed public trajectory')
    loaded[condition,m]=(valid,r,trace,pred);sources.append(dict(task=t,job=j,receipt_status=r.get('status'),valid_cap=valid))
    select(t,loaded[condition,m],480,'FINAL_CAP_ELIGIBLE')
    for cap in rule['CPU_checkpoints']:select(t,loaded[condition,m],cap,'CPU_'+str(cap))
  specifications=[]
  for condition in rule['grid_conditions']:
   for l,r in rule['actual_spend']['within_condition_method_pairs']:specifications.append(((condition,l),(condition,r),'METHOD_PAIR_'+condition+'_'+l+'__'+r))
  for m in rule['methods']:specifications.append((('matched24',m),('cross48to24',m),'GRID_PAIR_'+m))
  for left,right,pid in specifications:
   lt,rt=by[cid,*left],by[cid,*right];lT=jobs[lt['task_id']].get('child_process_CPU_seconds');rT=jobs[rt['task_id']].get('child_process_CPU_seconds');B=min(lT,rT) if lT is not None and rT is not None else None
   lr=select(lt,loaded[left],B,'ACTUAL_CPU_PAIRED',pid,'left');rr=select(rt,loaded[right],B,'ACTUAL_CPU_PAIRED',pid,'right');pairs.append(dict(case_id=cid,pair_id=pid,common_B=B,left_selection_row=len(rows)-2,right_selection_row=len(rows)-1,both_available=lr['selected'] is not None and rr['selected'] is not None,exact_equal_execution_spend_claimed=False,scope='Both full original children charged once; public trajectory readout at commonB, not equal actual consumedCPU'))
 np.savez_compressed(a.output/'SELECTED_MATERIALS.npz',m=np.stack(mats) if mats else np.empty((0,576)),Y_prediction=np.stack(ypreds) if ypreds else np.empty((0,48,12),complex))
 dump(a.output/'SELECTIONS.json',dict(rows=rows,pairs=pairs,sources=sources,truth_access='NONE'))
 for f in ('TASKS_V2.json','SOLVER_CONFIG_V2.json','SCORING_CONFIG_V2.json'):shutil.copy2(HERE/f,a.output/f)
 receipt=dict(status='ALL64_PUBLIC_PAIRED_GRID_CHOICES_SEALED_BEFORE_GT',finished_utc=datetime.now(timezone.utc).isoformat(),tasks=64,independent_materials=8,grid_conditions=2,methods=4,selections=len(rows),pairs=len(pairs),truth_access='NONE',new_solver=0,new_FW=0,CPU_seconds_near_write=time.process_time(),wall_seconds=time.perf_counter()-wall,actual_total_child_CPU_seconds=sum(j['child_process_CPU_seconds'] for j in jobs.values()) if all(j.get('child_process_CPU_seconds') is not None for j in jobs.values()) else None,cost_scope='Each original child charged once; pair/deployment views nonadditive; future external sealer cost replaces nested self timer')
 dump(a.output/'RECEIPT.json',receipt);print(json.dumps(receipt))
if __name__=='__main__':main()
