#!/usr/bin/env python3
"""Portable replay of saved, already-EXPOSED R4 material-risk actions.

Only Python+NumPy. Explicit evidence root; no original absolute paths, source
identity gates, model imports, inference, Maxwell, model selection or new tests.
The output directory must be new. This is offline arithmetic reproduction,
not reproduction of acquisition/training or posterior qualification.
"""
from pathlib import Path
import argparse,csv,json,os,re,resource,sys,time
import numpy as np
STAGE='p3l_r4_conditional_transport';BOOT_SEED=20261025099;DRAWS=4000;LEVEL=.9833333333333333
REFS=('frozen_r2','strong_spectral','fusion_diag','fusion_permuted')

def load_json(path):
 with Path(path).open() as f:return json.load(f)
def write_new(path,value):
 with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
def locate(root):
 # Fixed layouts only, never a recursive filesystem search.
 candidates=[root/'research'/STAGE,root/'workspace/research'/STAGE]
 if root.name==STAGE:candidates.append(root)
 found=[p for p in candidates if (p/'results/fresh/locked/risk45/scores/RESULTS.json').is_file()]
 if len(found)!=1:raise ValueError('Specify repository/workspace root, extraction root containing workspace/, or the R4 stage directory; exactly one matching fixed layout required')
 return found[0]
def safe_id(value):
 if not isinstance(value,str) or not re.fullmatch(r'[A-Za-z0-9_]+',value):raise ValueError('unsafe method or parent ID')
 return value
def bootstrap(candidate,reference):
 rng=np.random.default_rng(BOOT_SEED);absolute=[];relative=[]
 for begin in range(0,DRAWS,100):
  index=rng.integers(0,len(candidate),size=(min(100,DRAWS-begin),len(candidate)))
  cm=candidate[index].mean(1);rm=reference[index].mean(1);absolute.extend(rm-cm);relative.extend((rm-cm)/rm)
 tail=(1-LEVEL)/2
 return dict(n_parent=len(candidate),candidate_mean=float(candidate.mean()),reference_mean=float(reference.mean()),improvement_absolute=float((reference-candidate).mean()),improvement_relative=float((reference.mean()-candidate.mean())/reference.mean()),absolute_CI=np.quantile(absolute,[tail,1-tail]).tolist(),relative_CI=np.quantile(relative,[tail,1-tail]).tolist())
def write_csv(path,rows):
 with Path(path).open('x',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
def run(args,stage):
 used=[];methods_table=[];family_table=[];comparisons=[];numerical=[];failures=[];source_reads=0
 def use(path,role):
  nonlocal source_reads
  source_reads+=1;used.append(dict(path=path.relative_to(stage).as_posix(),role=role,bytes=path.stat().st_size))
 for cohort,prefix,n in [('locked','R4L1',512),('confirmation','R4L2',256)]:
  base=stage/'results/fresh'/cohort;boundpath=base/'BOUND_RISK.json';scorepath=base/'risk45/scores/RESULTS.json';use(boundpath,'compact frozen method/seed/parent table');use(scorepath,'compact original closed numerical results')
  bound=load_json(boundpath);closed=load_json(scorepath);methods=bound['methods'];parents=[f'{prefix}_{i:05d}' for i in range(n)]
  if bound['parents']!=parents or closed['status']!='OFFLINE_SCORE_COMPLETE' or closed['n_parent']!=n or len(methods)!=45 or set(closed['methods'])!={m['id'] for m in methods}:raise ValueError('fixed45 methods/512or256 parent roster differs from closed tables')
  if bound['risk_rule']['primary_family']!='fusion_spline' or bound['risk_rule']['adapt_after_L1'] is not False:raise ValueError('unexpected frozen primary rule; replay never chooses another winner')
  truth=[]
  for parent in parents:
   safe_id(parent);path=stage/'inputs/offline_only'/(parent+'.npz');use(path,'saved exposed material truth; only chi key')
   with np.load(path,allow_pickle=False) as z:chi=z['chi']
   if chi.shape!=(576,) or not np.iscomplexobj(chi) or not np.isfinite(chi).all():raise ValueError('finite full576 complex chi required')
   truth.append(chi)
  truth=np.asarray(truth,dtype=np.complex128);loss={}
  for m in methods:
   mid=safe_id(m['id']);path=base/'risk45/online'/(mid+'.npz');use(path,'saved canonical risk export; parent/raw/clipped/fallback keys only');loss[mid]={}
   with np.load(path,allow_pickle=False) as z:
    if z['parent'].astype(str).tolist()!=parents:raise ValueError('saved action parent order mismatch')
    fallback=z['fallback_anchor'].astype(bool)
    if fallback.shape!=(n,) or int(fallback.sum())!=closed['methods'][mid]['fallback_count']:raise ValueError('fallback denominator differs from closed result')
    if fallback.any():failures.append(dict(cohort=cohort,method=mid,parents=np.asarray(parents)[fallback].tolist(),policy='original anchor fallback retained; no parent dropping'))
    for action in ['raw','clipped']:
     pred=z['prediction_'+action]
     if pred.shape!=(n,576) or not np.iscomplexobj(pred) or not np.isfinite(pred).all():raise ValueError('invalid stored action; cannot drop parents')
     error=np.asarray(pred,dtype=np.complex128)-truth;square=error.real*error.real+error.imag*error.imag
     metrics=dict(full_mse=square.mean(1),full_mae=np.sqrt(square).mean(1));loss[mid][action]=metrics
     for metric,values in metrics.items():
      original=closed['methods'][mid][action][metric]['mean'];difference=abs(float(values.mean())-original);numerical.append(dict(cohort=cohort,method=mid,action=action,metric=metric,max_absolute_difference=difference))
      if difference>1e-12:raise ValueError('saved-array mean disagrees with compact closed result')
     methods_table.append(dict(cohort=cohort,method=mid,family=m['family'],seed=m['seed'],action=action,parents=n,fallback_count=int(fallback.sum()),full_mse=float(metrics['full_mse'].mean()),full_mae=float(metrics['full_mae'].mean())))
  family_loss={}
  for family in dict.fromkeys(m['family'] for m in methods):
   selected=[m for m in methods if m['family']==family]
   if len({m['seed'] for m in selected})!=len(selected):raise ValueError('duplicate seed in family; no independent-seed inflation')
   family_loss[family]={}
   for action in ['raw','clipped']:
    metrics={k:np.mean([loss[m['id']][action][k] for m in selected],axis=0) for k in ['full_mse','full_mae']};family_loss[family][action]=metrics
    family_table.append(dict(cohort=cohort,family=family,action=action,parents=n,seeds=len(selected),full_mse=float(metrics['full_mse'].mean()),full_mae=float(metrics['full_mae'].mean()),aggregation='mean seed losses per parent, not an ensemble image'))
  # Exactly the existing fixed primary/control comparisons, no new family or
  # pooled hypothesis. Extra controls keep their original descriptive scope.
  for action in ['raw','clipped']:
   for reference in REFS:
    previous=[c for c in closed['comparisons'] if c['candidate_family']=='fusion_spline' and c['reference_family']==reference and c['action']==action]
    frozen=[c for c in bound['comparisons'] if c['candidate_family']=='fusion_spline' and c['reference_family']==reference and c['action']==action]
    if len(previous)!=1 or len(frozen)!=1:raise ValueError('fixed comparison missing from compact source tables')
    observed=bootstrap(family_loss['fusion_spline'][action]['full_mse'],family_loss[reference][action]['full_mse'])
    difference=max(float(np.max(abs(np.asarray(observed[k])-np.asarray(previous[0][k])))) for k in ['candidate_mean','reference_mean','improvement_absolute','improvement_relative','absolute_CI','relative_CI'])
    if difference>1e-12:raise ValueError('independent parent-bootstrap result differs from closed table')
    comparisons.append(dict(cohort=cohort,candidate_family='fusion_spline',reference_family=reference,action=action,**observed,comparison_max_abs_difference=difference,draws=DRAWS,seed=BOOT_SEED,CI_level=LEVEL,scope='original primary risk comparison' if reference=='frozen_r2' and action=='clipped' else 'existing supporting descriptive interval; no added familywise claim'))
  np.savez_compressed(args.out/(cohort+'_family_parent_losses.npz'),parent=np.asarray(parents),**{family+'__'+action+'__'+metric:values for family,actions in family_loss.items() for action,metrics in actions.items() for metric,values in metrics.items()})
 write_csv(args.out/'METHODS.csv',methods_table);write_csv(args.out/'FAMILIES.csv',family_table)
 write_csv(args.out/'COMPARISONS.csv',[{**{k:v for k,v in row.items() if k not in ['absolute_CI','relative_CI']},'absolute_CI_low':row['absolute_CI'][0],'absolute_CI_high':row['absolute_CI'][1],'relative_CI_low':row['relative_CI'][0],'relative_CI_high':row['relative_CI'][1]} for row in comparisons])
 write_new(args.out/'RESULTS.json',dict(status='SAVED_ARRAY_ARITHMETIC_MATCHES_COMPACT_CLOSED_RESULTS',exposed_offline_only=True,parents=dict(locked=512,confirmation=256),methods_per_cohort=45,family_table=family_table,comparisons=comparisons,failures=failures,max_mean_absolute_difference=max(v['max_absolute_difference'] for v in numerical),max_comparison_absolute_difference=max(v['comparison_max_abs_difference'] for v in comparisons),numerical_reconciliation=numerical,source_files=used,source_files_read=source_reads,relative_inputs_only=True,original_absolute_paths_dereferenced=False,model_loads=0,model_inference=0,physics_RHS=0,new_hypotheses=0,interpretation='recomputed scores of already saved actions only; no new posterior/inverse-solver/architecture qualification'))
 return dict(source_files_read=source_reads,failures=len(failures),max_mean_absolute_difference=max(v['max_absolute_difference'] for v in numerical),max_comparison_absolute_difference=max(v['comparison_max_abs_difference'] for v in comparisons))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--offline-exposed-data',action='store_true',help='acknowledge reading already-exposed saved truth for offline arithmetic replay; no fresh blind test')
 args=p.parse_args()
 if not args.offline_exposed_data:p.error('explicit --offline-exposed-data acknowledgment is required')
 if args.out.exists():raise FileExistsError('NEW output directory required; preserve prior replay')
 root=args.evidence_root.expanduser().absolute();stage=locate(root);args.out=args.out.expanduser().absolute();args.out.mkdir(parents=True,exist_ok=False);wall=time.perf_counter();cpu=time.process_time();status='FAILED';error=None;facts={}
 try:facts=run(args,stage);status='EXPOSED_OFFLINE_REPLAY_COMPLETE'
 except BaseException as exc:error=repr(exc);raise
 finally:
  write_new(args.out/'RECEIPT.json',dict(status=status,error=error,evidence_root_argument=str(root),located_stage=str(stage),output=str(args.out),offline_exposed_data_acknowledged=True,wall_seconds=time.perf_counter()-wall,CPU_seconds=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),RHS=0,new_training=0,new_inference=0,new_hypotheses=0,Book_used=False,research_ledger_modified=False,old_absolute_source_checks_used=False,hash_checks=False,**facts))
 print(json.dumps(dict(status=status,out=str(args.out),**facts)))
if __name__=='__main__':main()
