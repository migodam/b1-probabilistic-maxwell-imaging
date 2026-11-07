"""Separate all288-terminal root-authorized scoring. Import performs no GT read."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import argparse,csv,json,time,itertools
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
METHODS=('cie_single','cie_dual','global_dual','fullG19','freezeA','contractedT')
def read(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def resolve(p):
 p=Path(p);return p if p.is_absolute() else REPO/p
def arrays(p):
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def method(t):return t.get('material_readout',{}).get('label',t['arm'])
def writecsv(p,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with Path(p).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def metrics(m,truth,high):
 m=np.asarray(m);truth=np.asarray(truth)
 if m.shape!=(576,) or truth.shape!=(576,) or np.iscomplexobj(m) or not np.isfinite(m).all():raise ValueError('576 finite real cells required')
 amplitude=float(np.max(truth)-.1);d=m-truth;rmse=float(np.linalg.norm(d)/24);den=np.count_nonzero((m>1)|(truth>1));iou=float(np.count_nonzero((m>1)&(truth>1))/den) if den else 1.
 legal=bool(np.min(m)>=.02-1e-12 and np.max(m)<=20+1e-12)
 row=dict(RMSE=rmse,NRMSE=rmse/amplitude,fixed_support_IoU=iou,material_legal=legal,material_success=bool(legal and rmse/amplitude<=.15 and iou>=.7),material_min=float(np.min(m)),material_max=float(np.max(m)))
 for label,value in [('low',2.5),('high',float(high)),('background',.1)]:
  mask=np.isclose(truth,value,rtol=0,atol=1e-12);row['true_'+label+'_cells']=int(mask.sum());row['RMSE_true_'+label+'_region']=float(np.sqrt(np.mean(d[mask]**2))) if mask.any() else None
 return row

def gate(a):
 if a.authorize_postseal!='ROOT_EXPLICIT_POSTSEAL_AUTHORIZATION' or not a.authorization_note.strip():raise PermissionError('Root postseal authorization required; no GT opened')
 frozen=read(a.freeze)
 if frozen['status']!='FROZEN_BEFORE_GENERATION_AND_CONFIRMATION_GT_ACCESS':raise ValueError('Scoring freeze required')
 if (frozen['primary_support']['threshold'],frozen['primary_success']['nrmse_threshold'],frozen['primary_success']['fixed_support_iou_threshold'])!=(1.,.15,.7):raise ValueError('Frozen metric differs')
 tasks=[];jobs={}
 for p in a.manifest:tasks.extend(read(p)['tasks'])
 for p in a.queue:
  q=read(p)
  if q['status'] not in ('COMPLETED','COMPLETED_WITH_ERRORS','STOPPED_DISPATCH_CHILDREN_FINISHED') or q.get('running'):raise RuntimeError('Every dispatch must be terminal before GT access')
  if q['status']!='STOPPED_DISPATCH_CHILDREN_FINISHED' and (q.get('pending') or q.get('not_started_task_ids')):raise RuntimeError('Unfinished active queue')
  # A deliberately stopped dispatcher contributes only its terminal jobs.
  # Aggregate unique-job equality with the original288-task manifest below
  # proves all its formerly unstarted tasks also ended in a successor queue.
  for j in q['jobs']:
   if j['task_id'] in jobs or j['status'] not in ('FINISHED','FAILED','LAUNCH_FAILED'):raise RuntimeError('Duplicate/nonterminal job')
   jobs[j['task_id']]=j
 ids=[t['task_id'] for t in tasks];cases=sorted({t['case_id'] for t in tasks})
 if len(tasks)!=288 or len(set(ids))!=288 or len(cases)!=48 or set(jobs)!=set(ids):raise ValueError('Require exactly288 terminal outcomes/48 cases')
 entries=[]
 for case in cases:
  if sorted(method(t) for t in tasks if t['case_id']==case)!=sorted(METHODS):raise ValueError('Six balanced methods required')
 for t in tasks:
  d=resolve(t['output_dir']);p=d/'RECEIPT.json';r=read(p) if p.exists() else {};j=jobs[t['task_id']]
  if r.get('status') in ('RUNNING','PREPARING','PREPARED'):raise RuntimeError('Nonterminal prediction receipt')
  if j['status']=='FINISHED' and r.get('status') not in ('SEALED_PUBLIC_PREDICTION','CPU_BUDGET_REACHED'):raise RuntimeError('False finished child')
  sealed=r.get('cap_prediction_status')=='SEALED_PAID_INCUMBENT'
  if sealed:
   if (r.get('task_id'),r.get('case_id'),r.get('arm'))!=(t['task_id'],t['case_id'],t['arm']):raise ValueError('Receipt identity mismatch before GT')
   public_task=read(d/'PUBLIC_TASK.json')
   if (public_task.get('task_id'),public_task.get('case_id'),method(public_task))!=(t['task_id'],t['case_id'],method(t)):raise ValueError('Task method mismatch before GT')
   for name in ('R11_CAP_PREDICTION.npz','R11_INCUMBENT_TRAJECTORY.npz','R11_CPU_CHECKPOINTS.json'):
    if not (d/name).is_file():raise RuntimeError('Missing seal output '+str(d/name))
  entries.append(dict(task=t,job=j,receipt=r,directory=d,sealed=sealed))
 return frozen,entries,cases

def paired(rows,replicates=10000,seed=202611071108):
 rng=np.random.default_rng(seed);out=[]
 views=sorted({r['view'] for r in rows})
 for view in views:
  lookup={(r['case_id'],r['method']):r for r in rows if r['view']==view}
  for left,right in itertools.combinations(METHODS,2):
   pairs=[(c,lookup[(c,left)],lookup[(c,right)]) for c in sorted({c for c,m in lookup}) if (c,left) in lookup and (c,right) in lookup]
   available=[(c,l,r) for c,l,r in pairs if l.get('NRMSE') is not None and r.get('NRMSE') is not None]
   base=dict(view=view,left=left,right=right,sign='left minus right',objects=len(pairs),paired_finite=len(available),missing_pairs=len(pairs)-len(available),left_only_success=sum(bool(l['material_success']) and not bool(r['material_success']) for c,l,r in pairs),right_only_success=sum(bool(r['material_success']) and not bool(l['material_success']) for c,l,r in pairs))
   if available:
    diffs=np.array([l['NRMSE']-r['NRMSE'] for c,l,r in available]);regimes=np.array([l['regime'] for c,l,r in available]);boot=[]
    for regime in ('high','mixed'):
     values=diffs[regimes==regime]
     if values.size:boot.append(values[rng.integers(0,values.size,size=(replicates,values.size))])
    sample=np.concatenate(boot,axis=1);means=sample.mean(axis=1);medians=np.median(sample,axis=1)
    base.update(mean_NRMSE_difference=float(diffs.mean()),median_NRMSE_difference=float(np.median(diffs)),mean95low=float(np.quantile(means,.025)),mean95high=float(np.quantile(means,.975)),median95low=float(np.quantile(medians,.025)),median95high=float(np.quantile(medians,.975)),bootstrap_replicates=replicates,bootstrap_seed=seed,stratified_by='high/mixed independent objects; only finite paired scores, missingness reported')
   out.append(base)
 return out

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--manifest',type=Path,action='append',required=True);ap.add_argument('--queue',type=Path,action='append',required=True);ap.add_argument('--freeze',type=Path,default=ROOT/'configs/twod_independent48_scoring.json');ap.add_argument('--truth-dir',type=Path,default=ROOT/'sealed_truth/independent48');ap.add_argument('--output',type=Path,default=ROOT/'delegated/twod/independent48/scored');ap.add_argument('--authorize-postseal',required=True);ap.add_argument('--authorization-note',required=True);a=ap.parse_args();wall=time.perf_counter();frozen,entries,cases=gate(a)
 a.output.mkdir(parents=True,exist_ok=False)
 dump(a.output/'PUBLIC_GATE_BEFORE_TRUTH.json',dict(status='PASS_ALL288_TERMINAL',truth_opened=False,cases=cases,tasks=288,authorization_note=a.authorization_note))
 # First private read occurs only after all public gate checks and audit write.
 meta={r['case_id']:r for r in read(a.truth_dir/'generator_metadata.json')['records']};truths={c:arrays(a.truth_dir/(c+'.npz'))['m'] for c in cases}
 rows=[];trajectories={};frontier_path=a.output/'ACTUAL_COST_FRONTIER.csv';frontier_file=frontier_path.open('w',newline='');frontier_writer=None
 try:
  for e in entries:
   t,j,r,d=e['task'],e['job'],e['receipt'],e['directory'];c=t['case_id'];private=meta[c]
   base=dict(case_id=c,method=method(t),family=private['family'],regime=private['regime'],high_chi=private['high_chi'],task_status=r.get('status',j['status']),actual_child_CPU_seconds=j.get('child_process_CPU_seconds'),final_physical_LU=(r.get('physical_model_cost') or {}).get('LU'),final_complex_RHS=(r.get('physical_model_cost') or {}).get('actual_complex_RHS'))
   trace=arrays(d/'R11_INCUMBENT_TRAJECTORY.npz') if e['sealed'] else None;trajectories[(c,method(t))]=(e,base,trace)
   final=dict(base,view='FINAL_CAP_ELIGIBLE',material_success=False,available=False)
   if e['sealed']:final.update(available=True,original_Phi=r['cap_physical_objective'],**metrics(arrays(d/'R11_CAP_PREDICTION.npz')['m'],truths[c],private['high_chi']))
   rows.append(final)
   for cap in frozen['CPU_checkpoints_seconds']:
    row=dict(base,view='CPU_'+str(cap),CPU_cap_seconds=cap,material_success=False,available=False)
    eligible=np.flatnonzero(trace['CPU_seconds']<=cap) if trace is not None else np.array([],int)
    if eligible.size:
     k=int(eligible[np.argmin(trace['original_Phi'][eligible])]);row.update(available=True,paid_CPU_seconds=float(trace['CPU_seconds'][k]),paid_LU=int(trace['LU'][k]),original_Phi=float(trace['original_Phi'][k]),**metrics(trace['m'][k],truths[c],private['high_chi']))
    rows.append(row)
   if trace is not None:
    for k in range(len(trace['m'])):
     fr=dict(base,paid_CPU_seconds=float(trace['CPU_seconds'][k]),paid_LU=int(trace['LU'][k]),paid_complex_RHS=int(trace['complex_RHS'][k]),original_Phi=float(trace['original_Phi'][k]),**metrics(trace['m'][k],truths[c],private['high_chi']))
     fr['frontier_point']='PAID_INCUMBENT';fr['eligible_at480']=bool(fr['paid_CPU_seconds']<=480)
     if frontier_writer is None:frontier_writer=csv.DictWriter(frontier_file,fieldnames=list(fr));frontier_writer.writeheader()
     frontier_writer.writerow(fr)
    # A late atomic state stays excluded from CPU checkpoints, but its actual
    # paid cost/result is visible on the actual-cost frontier, without GT choice.
    actual_path=d/'R11_ACTUAL_PAID_INCUMBENT.json'
    if actual_path.exists():
     actual=read(actual_path);paid=read(d/'R11_PAID_PHYSICAL_STATES.json');matching=[q for q in paid if q['original_Phi']==actual['objective']];state_cost=min(matching,key=lambda q:q['CPU_seconds']) if matching else None
     if state_cost is not None and state_cost['CPU_seconds']>480:
      actual_m=arrays(d/'R11_ACTUAL_PAID_INCUMBENT.npz')['m'];pc=state_cost['physical_model_cost']
      fr=dict(base,paid_CPU_seconds=state_cost['CPU_seconds'],paid_LU=pc.get('LU'),paid_complex_RHS=pc.get('actual_complex_RHS'),original_Phi=actual['objective'],**metrics(actual_m,truths[c],private['high_chi']),frontier_point='ATOMIC_OVERRUN_ACTUAL_ONLY',eligible_at480=False);frontier_writer.writerow(fr)
    # Show all actual terminal expenditure even when last improvement was earlier.
    if final.get('available'):
     fr=dict(base,paid_CPU_seconds=base['actual_child_CPU_seconds'],paid_LU=base['final_physical_LU'],paid_complex_RHS=base['final_complex_RHS'],original_Phi=final['original_Phi'],**{k:final[k] for k in metrics(arrays(d/'R11_CAP_PREDICTION.npz')['m'],truths[c],private['high_chi'])},frontier_point='TERMINAL_SPEND_CAP_ELIGIBLE_RESULT',eligible_at480=True);frontier_writer.writerow(fr)
  for c in cases:
   group=[trajectories[(c,m)] for m in METHODS];valid=[int(b['final_physical_LU']) for e,b,tr in group if b.get('final_physical_LU') is not None];common=min(valid) if len(valid)==6 else None
   for e,b,tr in group:
    row=dict(b,view='MATCHED_PHYSICAL_LU',common_LU_cap=common,material_success=False,available=False)
    eligible=np.flatnonzero(tr['LU']<=common) if tr is not None and common is not None else np.array([],int)
    if eligible.size:
     k=int(eligible[np.argmin(tr['original_Phi'][eligible])]);row.update(available=True,paid_CPU_seconds=float(tr['CPU_seconds'][k]),paid_LU=int(tr['LU'][k]),original_Phi=float(tr['original_Phi'][k]),**metrics(tr['m'][k],truths[c],meta[c]['high_chi']))
    rows.append(row)
 finally:frontier_file.close()
 writecsv(a.output/'SCORES_AND_FIXED_BUDGETS.csv',rows)
 summary=[]
 for view in sorted({r['view'] for r in rows}):
  for m in METHODS:
   for key,value in [('overall','ALL')]+[(key,value) for key in ('regime','family','high_chi') for value in sorted({r[key] for r in rows},key=str)]:
    selected=[r for r in rows if r['view']==view and r['method']==m and (key=='overall' or r[key]==value)];finite=[r['NRMSE'] for r in selected if r.get('NRMSE') is not None]
    summary.append(dict(view=view,method=m,group=key,value=value,objects=len(selected),finite_scores=len(finite),missing=len(selected)-len(finite),material_success=sum(bool(r['material_success']) for r in selected),mean_NRMSE=float(np.mean(finite)) if finite else None,median_NRMSE=float(np.median(finite)) if finite else None))
 writecsv(a.output/'AGGREGATES.csv',summary);writecsv(a.output/'PAIRED_STRATIFIED_BOOTSTRAP.csv',paired(rows));dump(a.output/'RECEIPT.json',dict(status='SCORED_POST_ALL288_TERMINAL',truth_access='ROOT_EXPLICIT_POSTSEAL_ONLY',GT_selected_candidate=False,GT_selected_iterate=False,physics_runs=0,CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall,objects=48,tasks=288,frontier_scope='All cap-eligible paid incumbents, explicit atomic-overrun actual-only point if present, and terminal actual-spend endpoint. No overshoot used for CPU checkpoint prediction.',unknown_costs=['External whole scorer process CPU after final write not recorded by self timer']))
 print(json.dumps(dict(status='SCORED_POST_ALL288_TERMINAL',output=str(a.output))))
if __name__=='__main__':main()
