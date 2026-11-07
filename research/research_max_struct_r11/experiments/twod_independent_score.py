"""Separate all288-terminal root-authorized scoring. Import performs no GT read."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import argparse,csv,json,time,itertools
SCORER_MODULE_CPU_START=time.process_time()
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.stats import chi2
from twod_recovery import active_task,active_job,queue_receipt,attempt_cost_info,affected_objects
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
METHODS=('cie_single','cie_dual','global_dual','fullG19','freezeA','contractedT')
def read(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def resolve(p):
 p=Path(str(p).replace('\\','/'));return p if p.is_absolute() else REPO/p
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


def public_observation(task):
 """Read only public Y/sigma, without loading G or calling any physics."""
 observed={}
 forbidden={'truth','gt','ground_truth','m_true','m_truth','true_material','y_clean','yclean','family','seed','mask','alpha'}
 for filename in task['public_files']:
  path=resolve(filename)
  if {'data_private','private','truth','gt','sealed_gt','sealed_truth'} & {part.lower() for part in path.parts}:raise ValueError('Private path in public residual input')
  with np.load(path,allow_pickle=False) as z:
   if forbidden & {key.lower() for key in z.files}:raise ValueError('Private payload in public residual input')
   for key in ('Y','sigma'):
    if key in z.files:
     value=z[key].copy()
     if key in observed and not np.array_equal(observed[key],value):raise ValueError('Conflicting public noise/observations')
     observed[key]=value
 if set(observed)!={'Y','sigma'}:raise ValueError('Public Y and scalar sigma required')
 sigma=np.asarray(observed['sigma'])
 if sigma.size!=1 or np.iscomplexobj(sigma) or not np.isfinite(sigma).all() or float(sigma.reshape(-1)[0])<=0:raise ValueError('Positive finite scalar raw CN sigma required')
 return np.asarray(observed['Y'],complex),float(sigma.reshape(-1)[0])

def data_metrics(prediction,Y,sigma):
 prediction=np.asarray(prediction,complex)
 if prediction.shape!=Y.shape or not np.isfinite(prediction).all():raise ValueError('Finite sealed prediction with public observation shape required')
 norm=float(np.linalg.norm(prediction-Y));ynorm=float(np.linalg.norm(Y))
 if ynorm<=0:raise ValueError('Nonzero public data norm required')
 radius=float(np.sqrt(chi2.ppf(.99,2*Y.size)/2))
 return dict(data_relative_residual=norm/ynorm,data_residual_Frobenius_norm=norm,
     data_residual_over_sigma=norm/sigma,CN99_radius_over_sigma=radius,
     CN99_pass=bool(norm/sigma<=radius),complex_observations=int(Y.size),public_sigma=sigma,
     data_residual_source='Exact sealed paid Y_prediction; no forward/adjoint; E|CN|^2=sigma^2')

def gate(a):
 if a.authorize_postseal!='ROOT_EXPLICIT_POSTSEAL_AUTHORIZATION' or not a.authorization_note.strip():raise PermissionError('Root postseal authorization required; no GT opened')
 frozen=read(a.freeze)
 if frozen['status']!='FROZEN_BEFORE_GENERATION_AND_CONFIRMATION_GT_ACCESS':raise ValueError('Scoring freeze required')
 if (frozen['primary_support']['threshold'],frozen['primary_success']['nrmse_threshold'],frozen['primary_success']['fixed_support_iou_threshold'])!=(1.,.15,.7):raise ValueError('Frozen metric differs')
 tasks=[];jobs={}
 for p in a.manifest:tasks.extend(active_task(t) for t in read(p)['tasks'])
 for p in a.queue:
  p=queue_receipt(p);q=read(p)
  if q['status'] not in ('COMPLETED','COMPLETED_WITH_ERRORS','STOPPED_DISPATCH_CHILDREN_FINISHED','INTERRUPTED_CONTROLLER_RETAINED_TERMINAL_PREFIX') or q.get('running'):raise RuntimeError('Every dispatch must be terminal before GT access')
  if q['status'] not in ('STOPPED_DISPATCH_CHILDREN_FINISHED','INTERRUPTED_CONTROLLER_RETAINED_TERMINAL_PREFIX') and (q.get('pending') or q.get('not_started_task_ids')):raise RuntimeError('Unfinished active queue')
  # A deliberately stopped dispatcher contributes only its terminal jobs.
  # Aggregate unique-job equality with the original288-task manifest below
  # proves all its formerly unstarted tasks also ended in a successor queue.
  for j in q['jobs']:
   if not active_job(j):continue
   if j['task_id'] in jobs or j['status'] not in ('FINISHED','FAILED','LAUNCH_FAILED'):raise RuntimeError('Duplicate/nonterminal job')
   jobs[j['task_id']]=dict(j,source_queue_receipt=str(p))
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
 secondary=read(a.secondary_public_seal/'RECEIPT.json')
 if secondary.get('status')!='ALL96_PUBLIC_DEPLOYMENTS_SEALED_BEFORE48_GT' or secondary.get('deployments')!=96 or secondary.get('truth_access')!='NONE':raise RuntimeError('Root-required public secondary96 seal must precede any48GT read')
 choices=read(a.secondary_public_seal/'DEPLOYMENT_SELECTIONS.json')['rows']
 if len(choices)!=96 or sorted({r['case_id'] for r in choices})!=cases:raise RuntimeError('Secondary seal identity mismatch before GT')
 global_appendix=read(a.global_appendix_seal/'RECEIPT.json')
 if global_appendix.get('status')!='ADDITIONAL48_GLOBAL_PUBLIC_CHOICES_SEALED_BEFORE48_GT' or global_appendix.get('additional_choices')!=48 or not global_appendix.get('original96_unchanged') or global_appendix.get('truth_access')!='NONE':raise RuntimeError('Append-only48 Global choices must seal before48GT')
 global_choices=read(a.global_appendix_seal/'DEPLOYMENT_SELECTIONS.json')['rows']
 if len(global_choices)!=48 or sorted({r['case_id'] for r in global_choices})!=cases or any(r['method']!='CIE_then_GlobalDual' for r in global_choices):raise RuntimeError('Global appendix choice identity mismatch')
 fair=read(a.actual_cpu_fairness_seal/'RECEIPT.json')
 if fair.get('status')!='ALL384_PUBLIC_ACTUAL_CPU_PAIRS_SEALED_BEFORE48_GT' or fair.get('pairs')!=384 or fair.get('truth_access')!='NONE':raise RuntimeError('Root-required384 actual-CPU pair choices must public-seal before48GT')
 pairs=read(a.actual_cpu_fairness_seal/'PAIRED_SELECTIONS.json')['rows']
 expected={(c,l,r) for c in cases for l in ('cie_single','cie_dual') for r in ('fullG19','freezeA','contractedT','global_dual')}
 if len(pairs)!=384 or {(p['case_id'],p['baseline'],p['comparator']) for p in pairs}!=expected:raise RuntimeError('Actual-CPU public selection identity mismatch')
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


def score_actual_cpu_pairs(directory,truths,metadata):
 """Selections are already public-sealed; no GT-based budget or state choice."""
 pairs=read(directory/'PAIRED_SELECTIONS.json')['rows'];mats=arrays(directory/'PAIRED_MATERIALS.npz')['m'];rows=[]
 for pair in pairs:
  case=pair['case_id']
  for side in ('left','right'):
   member=pair[side];selected=member['selected'];meta=metadata[case]
   row=dict(engineering_affected_object=pair.get('engineering_affected_object',False),exact_end_to_end_cost_qualified=pair.get('exact_end_to_end_cost_qualified',True),end_to_end_same_CPU_status=pair.get('end_to_end_same_CPU_status'),case_id=case,method=member['method'],view='ACTUAL_CPU_PAIR_'+pair['pair_id'],pair_id=pair['pair_id'],side=side,family=meta['family'],regime=meta['regime'],high_chi=meta['high_chi'],common_actual_end_CPU_B=pair['common_actual_end_CPU_B'],actual_full_child_end_CPU_seconds=member['actual_full_child_end_CPU_seconds'],nominal_unconsumed_CPU_allowance=member['nominal_unconsumed_CPU_allowance'],actual_CPU_overrun=member['actual_CPU_overrun'],actual_last_recorded_paid_physical_CPU_seconds=member['actual_last_recorded_paid_physical_CPU_seconds'],last_cap_eligible_paid_physical_CPU_seconds_at_B=member['last_cap_eligible_paid_physical_CPU_seconds_at_B'],B_minus_last_cap_eligible_physical_record=member['B_minus_last_cap_eligible_physical_record'],task_status=member['task_status'],source_queue_receipt=member['source_queue_receipt'],cpu_logical=member['cpu_logical'],material_success=False,available=selected is not None,data_metrics_at_B='NOT_RETAINED_AT_EVERY_INCUMBENT; no new forward')
   if selected is not None:
    index=int(selected['material_index']);paid=selected['paid_CPU_seconds'];B=pair['common_actual_end_CPU_B']
    if B is None or paid>B:raise RuntimeError('Late fairness state forbidden')
    row.update(attempt_cost_info(member['task_id'],member['actual_full_child_end_CPU_seconds']))
    row.update(original_Phi=selected['original_Phi'],paid_CPU_seconds=paid,paid_LU=selected['physical_LU'],paid_complex_RHS=selected['physical_complex_RHS'],B_minus_selected_paid_CPU_seconds=selected['B_minus_selected_paid_CPU_seconds'],**metrics(mats[index],truths[case],meta['high_chi']))
   rows.append(row)
 return rows

def main():
 scoring_body_CPU_start=time.process_time();scoring_started_utc=datetime.now(timezone.utc).isoformat()
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--manifest',type=Path,action='append',required=True);ap.add_argument('--queue',type=Path,action='append',required=True);ap.add_argument('--freeze',type=Path,default=ROOT/'configs/twod_independent48_scoring.json');ap.add_argument('--truth-dir',type=Path,default=ROOT/'sealed_truth/independent48');ap.add_argument('--output',type=Path,default=ROOT/'delegated/twod/independent48/scored');ap.add_argument('--authorize-postseal',required=True);ap.add_argument('--authorization-note',required=True);ap.add_argument('--secondary-public-seal',type=Path,default=ROOT/'delegated/twod/independent48/secondary/public_seal');ap.add_argument('--actual-cpu-fairness-seal',type=Path,default=ROOT/'delegated/twod/independent48/actual_cpu_fairness/public_seal');ap.add_argument('--global-appendix-seal',type=Path,default=ROOT/'delegated/twod/independent48/secondary/global_appendix/public_seal');a=ap.parse_args();wall=time.perf_counter();frozen,entries,cases=gate(a)
 a.output.mkdir(parents=True,exist_ok=False)
 dump(a.output/'PUBLIC_GATE_BEFORE_TRUTH.json',dict(status='PASS_ALL288_TERMINAL',created_utc=datetime.now(timezone.utc).isoformat(),truth_opened=False,cases=cases,tasks=288,authorization_note=a.authorization_note))
 # First private read occurs only after all public gate checks and audit write.
 meta={r['case_id']:r for r in read(a.truth_dir/'generator_metadata.json')['records']};truths={c:arrays(a.truth_dir/(c+'.npz'))['m'] for c in cases}
 rows=[];trajectories={};public_cache={};frontier_path=a.output/'ACTUAL_COST_FRONTIER.csv';frontier_file=frontier_path.open('w',newline='');frontier_writer=None
 try:
  for e in entries:
   t,j,r,d=e['task'],e['job'],e['receipt'],e['directory'];c=t['case_id'];private=meta[c]
   base=dict(case_id=c,method=method(t),source_queue_receipt=j.get('source_queue_receipt'),cpu_logical=j.get('cpu'),corrector_stop=r.get('corrector',{}).get('stop'),corrector_iterations=r.get('corrector',{}).get('iterations'),auxiliary_iterations=r.get('auxiliary_iterations'),family=private['family'],regime=private['regime'],high_chi=private['high_chi'],task_status=r.get('status',j['status']),actual_child_CPU_seconds=j.get('child_process_CPU_seconds'),final_physical_LU=(r.get('physical_model_cost') or {}).get('LU'),final_complex_RHS=(r.get('physical_model_cost') or {}).get('actual_complex_RHS'))
   base.update(attempt_cost_info(t['task_id'],j.get('child_process_CPU_seconds')));base['engineering_affected_object']=c in affected_objects()
   trace=arrays(d/'R11_INCUMBENT_TRAJECTORY.npz') if e['sealed'] else None;trajectories[(c,method(t))]=(e,base,trace)
   final=dict(base,view='FINAL_CAP_ELIGIBLE',material_success=False,available=False)
   if e['sealed']:
    sealed=arrays(d/'R11_CAP_PREDICTION.npz')
    if c not in public_cache:public_cache[c]=public_observation(t)
    Y,sigma=public_cache[c]
    final.update(available=True,original_Phi=r['cap_physical_objective'],**metrics(sealed['m'],truths[c],private['high_chi']),**data_metrics(sealed['Y_prediction'],Y,sigma))
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
     if frontier_writer is None:frontier_writer=csv.DictWriter(frontier_file,fieldnames=list(dict.fromkeys(list(fr)+[key for entry in entries for key in attempt_cost_info(entry['task']['task_id'],entry['job'].get('child_process_CPU_seconds'))])));frontier_writer.writeheader()
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
 fairness_rows=score_actual_cpu_pairs(a.actual_cpu_fairness_seal,truths,meta)
 writecsv(a.output/'ACTUAL_CPU_FAIRNESS_SCORES.csv',fairness_rows)
 fairness_paired=[r for r in paired(fairness_rows) if r['objects']>0]
 writecsv(a.output/'ACTUAL_CPU_FAIRNESS_PAIRED_BOOTSTRAP.csv',fairness_paired)
 affected=set(affected_objects())
 writecsv(a.output/'PAIRED_SENSITIVITY_EXCLUDING_INTERRUPTED_OBJECTS.csv',paired([r for r in rows if r['case_id'] not in affected]))
 writecsv(a.output/'ACTUAL_CPU_FAIRNESS_PAIRED_SENSITIVITY_EXCLUDING_INTERRUPTED_OBJECTS.csv',[p for p in paired([r for r in fairness_rows if r['case_id'] not in affected]) if p['objects']>0])
 dump(a.output/'ENGINEERING_SENSITIVITY_CONTRACT.json',dict(primary_all_objects_retained=48,excluded_objects_secondary_only=sorted(affected),sensitivity_objects=48-len(affected),reason='Root-authorized control-abort sensitivity, registered before48GT, not outcome-based deletion',exact_end_to_end_cost_unknown_for_affected=True))
 summary=[]
 for view in sorted({r['view'] for r in rows}):
  for m in METHODS:
   for key,value in [('overall','ALL')]+[(key,value) for key in ('regime','family','high_chi') for value in sorted({r[key] for r in rows},key=str)]:
    selected=[r for r in rows if r['view']==view and r['method']==m and (key=='overall' or r[key]==value)];finite=[r['NRMSE'] for r in selected if r.get('NRMSE') is not None]
    summary.append(dict(view=view,method=m,group=key,value=value,objects=len(selected),finite_scores=len(finite),missing=len(selected)-len(finite),material_success=sum(bool(r['material_success']) for r in selected),mean_NRMSE=float(np.mean(finite)) if finite else None,median_NRMSE=float(np.median(finite)) if finite else None))
 writecsv(a.output/'AGGREGATES.csv',summary);writecsv(a.output/'PAIRED_STRATIFIED_BOOTSTRAP.csv',paired(rows));dump(a.output/'RECEIPT.json',dict(status='SCORED_POST_ALL288_TERMINAL',scoring_started_utc=scoring_started_utc,scoring_finished_utc=datetime.now(timezone.utc).isoformat(),truth_access='ROOT_EXPLICIT_POSTSEAL_ONLY',GT_selected_candidate=False,GT_selected_iterate=False,physics_runs=0,CPU_seconds=time.process_time(),CPU_clock='Absolute whole scorer process_time since child creation; includes all imports and scoring through receipt preparation',module_entry_CPU_seconds=SCORER_MODULE_CPU_START,body_entry_CPU_seconds=scoring_body_CPU_start,body_elapsed_CPU_seconds=time.process_time()-scoring_body_CPU_start,wall_seconds=time.perf_counter()-wall,objects=48,tasks=288,secondary_conditional_deployments=144,original96_deployment_protocol_unchanged=True,secondary_actual_CPU_pairs=384,secondary_actual_CPU_rule='CommonB=min two actual full child end CPU, min originalPhi withinB; originalT_end and state granularity gaps retained; not exact equal billed CPU',frontier_scope='All cap-eligible paid incumbents, explicit atomic-overrun actual-only point if present, and terminal actual-spend endpoint. No overshoot used for CPU checkpoint prediction.',unknown_costs=['External whole scorer process CPU after final write not recorded by self timer']))
 print(json.dumps(dict(status='SCORED_POST_ALL288_TERMINAL',output=str(a.output))))
if __name__=='__main__':main()
