"""Seal all8 absorptive DEV predictions before any scoring GT read; fixed exposure only."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time,argparse,csv
from collections import Counter
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];OUT=ROOT/'delegated/diagnostics/absorptive_dev';sys.path.insert(0,str(REPO));R10=ROOT.parents[2]/'r10_worktree/research/research_max_struct_r10'
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else v.tolist())+'\n')
def csvwrite(p,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)
def seal():
 q=json.loads((OUT/'QUEUE.json').read_text());tasks=json.loads((OUT/'TASKS.json').read_text())['tasks'];jobs={j['task_id']:j for j in q['jobs']}
 if q['status']!='COMPLETED' or len(jobs)!=8 or q.get('pending'):raise RuntimeError('All8 must finish; preserve failures and report to root before GT')
 entries=[]
 for t in tasks:
  d=REPO/t['output_dir'];r=json.loads((d/'RECEIPT.json').read_text());p=d/'R11_CAP_PREDICTION.npz'
  if r.get('cap_prediction_status')!='SEALED_PAID_INCUMBENT' or r.get('truth_access')!='NONE' or r['status'] not in ('SEALED_PUBLIC_PREDICTION','CPU_BUDGET_REACHED'):raise RuntimeError('Unsealed or failed task; no GT read')
  entries.append(dict(task=t,job=jobs[t['task_id']],receipt_status=r['status'],prediction=str(p),prediction_bytes=p.stat().st_size,prediction_mtime_ns=p.stat().st_mtime_ns,original_Phi=r['cap_physical_objective'],truth_access_in_solver=r['truth_access']))
 dump(OUT/'ALL8_PREDICTIONS_SEALED_BEFORE_GT.json',dict(status='ALL8_SEALED_BEFORE_ANY_ABSORPTIVE_DEV_GT_SCORING',entries=entries,classification='EXPOSED_DEV_NOT_CONFIRMATION',selection='Existing same originalPhi eligible paid incumbent, no GT selection',new_forward=0,GT_access='NONE_THIS_PHASE'));print('All8 predictions sealed before scoring')
def score():
 wall,cpu=time.perf_counter(),time.process_time();seal=json.loads((OUT/'ALL8_PREDICTIONS_SEALED_BEFORE_GT.json').read_text())
 if seal['status']!='ALL8_SEALED_BEFORE_ANY_ABSORPTIVE_DEV_GT_SCORING' or len(seal['entries'])!=8:raise RuntimeError('All8 public seals required')
 from research.research_max_struct_r10.experiments.score_twod_comparison import material_metrics
 meta={r['case_id']:r for r in json.loads((R10/'configs/twod_generator.json').read_text())['cases']};meta.update({r['case_id']:r for r in json.loads((ROOT/'configs/twod_bothblocks16_generator.json').read_text())['cases']});truth={};rows=[];physical=Counter();aux=Counter();cps=[]
 for cid in sorted({e['task']['case_id'] for e in seal['entries']}):
  path=R10/'sealed_truth/twod'/f'{cid}.npz' if cid.startswith('r10_') else ROOT/'sealed_truth/bothblocks16'/f'{cid}.npz'
  with np.load(path) as z:truth[cid]=z['m'].copy()
 for e in seal['entries']:
  t=e['task'];cid=t['case_id'];d=REPO/t['output_dir'];r=json.loads((d/'RECEIPT.json').read_text());pub={}
  for p in t['public_files']:
   if Path(p).name=='operators.npz':continue
   with np.load(REPO/p) as z:pub.update(dict(z))
  with np.load(e['prediction']) as z:m=z['m'];yp=z['Y_prediction']
  gt=truth[cid];metrics=material_metrics(m,gt,float(gt.max()-.1),meta[cid]['high_chi']);cn=float(np.linalg.norm(yp-pub['Y'])/float(pub['sigma']));threshold=float(np.sqrt(chi2.ppf(.99,2*pub['Y'].size)/2));pc=r.get('physical_model_cost',{});ac=(r.get('auxiliary_execution_receipt') or {}).get('operation_counts',{});physical.update(pc);aux.update({k:v for k,v in ac.items() if not any(tag in k for tag in ('seconds','microseconds','bytes'))})
  row=dict(case_id=cid,method=t['material_readout']['label'],rank=19 if 'Green' in t['material_readout']['label'] else 0,classification='FIXED_EXPOSED_DEV_NOT_CONFIRMATION',status=r['status'],amplitude_NRMSE=metrics['amplitude_NRMSE'],fixed_support_IoU=metrics['fixed_support_IoU'],material_success=metrics['material_success'],material_legal=metrics['material_legal'],original_Phi=e['original_Phi'],whitened_residual=cn,CN99_threshold=threshold,data_CN99_pass=bool(cn<=threshold),external_complete_child_CPU_seconds=e['job']['external_complete_child_CPU_seconds'],near_exit_whole_child_CPU_seconds=r.get('whole_child_CPU_seconds'),original_aux_selected_stage=r.get('auxiliary_selected_stage'),corrector_stop=(r.get('corrector') or {}).get('stop'),corrector_iterations=(r.get('corrector') or {}).get('iterations'),physical_cost=pc,auxiliary_count_sums=ac,GT_access='Only after all8 absorptive predictions sealed; exposed old/publicly scored GT');rows.append(row)
  points=json.loads((d/'R11_CPU_CHECKPOINTS.json').read_text())
  with np.load(d/'R11_CHECKPOINT_MATERIALS.npz') as z:materials=z['m']
  for p in points:
   x=dict(case_id=cid,method=row['method'],**p)
   if 'material_index' in p:x.update(**material_metrics(materials[p['material_index']],gt,float(gt.max()-.1),meta[cid]['high_chi']))
   cps.append(x)
 dump(OUT/'ALL8_POSTSEAL_SCORES.json',rows);csvwrite(OUT/'ALL8_POSTSEAL_SCORES.csv',rows);dump(OUT/'CPU_CHECKPOINT_SCORES.json',cps)
 # References are already paid/scored endpoints, never new inverse inputs or extra costs.
 old=json.loads((R10/'delegated/twod_score/scored/SCORES_FINAL.json').read_text());realboth=json.loads((ROOT/'delegated/diagnostics/BOTHBLOCKS_DEV_SCORES.json').read_text());sixteen=json.loads((ROOT/'delegated/diagnostics/bothblocks_confirmation/SINGLE_ARM_SCORES.json').read_text());comparison=[]
 for cid in sorted(truth):
  ref=([dict(r,reference_label=r['method']) for r in old if r['case_id']==cid and r['method'] in ('cie_single','cie_dual','global_dual','matrix_green19')]+[dict(r,reference_label='real_GreenBothBlocks') for r in realboth if r['case_id']==cid]) if cid.startswith('r10_') else [dict(r,reference_label=r['method']) for r in sixteen if r['case_id']==cid]
  for x in ref:comparison.append(dict(case_id=cid,method=x['reference_label'],role='EXISTING_PAID_REFERENCE_NONADDITIVE',amplitude_NRMSE=x['amplitude_NRMSE'],fixed_support_IoU=x['fixed_support_IoU'],material_success=x['material_success'],original_Phi=x['original_Phi'],cost_scope='Existing R10/realBoth/16 reference; never add its cost again or compare cross-host rawCPU'))
 comparison+=[dict(case_id=r['case_id'],method=r['method'],role='NEW_FIXED_ABSORPTIVE_DEV',amplitude_NRMSE=r['amplitude_NRMSE'],fixed_support_IoU=r['fixed_support_IoU'],material_success=r['material_success'],original_Phi=r['original_Phi'],cost_scope='Actual new child measured separately') for r in rows];dump(OUT/'COMPARISON_WITH_EXISTING_REFERENCES.json',comparison);csvwrite(OUT/'COMPARISON_WITH_EXISTING_REFERENCES.csv',comparison)
 summary=[]
 for method in ('absorptiveGreenBothBlocks','absorptiveScalarBothBlocks'):
  rr=[r for r in rows if r['method']==method];summary.append(dict(method=method,objects=4,material_PASS=sum(r['material_success'] for r in rr),CN99_PASS=sum(r['data_CN99_pass'] for r in rr),mean_amplitude_NRMSE=float(np.mean([r['amplitude_NRMSE'] for r in rr])),actual_complete_child_CPU_seconds_sum=sum(r['external_complete_child_CPU_seconds'] for r in rr),scope='Fixed four exposed DEV objects; scalar control is not mature original CIE; no novelty/generalization claim'))
 dump(OUT/'SUMMARY.json',summary);q=json.loads((OUT/'QUEUE.json').read_text());wire=json.loads((OUT/'WIRING_EXTERNAL_COST.json').read_text());dump(OUT/'COST_LEDGER.json',dict(status='COMPLETE',classification='FIXED_EXPOSED_DEV',tasks=8,failed_tasks=0,single_local_lane=True,GPU_used=False,new_observation_generation=0,exact_inversion_child_CPU_seconds=q['external_child_CPU_seconds_sum'],controller_own_CPU_seconds=q['controller_own_CPU_seconds'],wiring_external_CPU_seconds=wire['CPU_seconds'],wiring_physical_FW=0,physical_cost=dict(physical),auxiliary_count_sums=dict(aux),reference_costs_nonadditive=True,scoring_new_forward=0,scoring_new_adjoint=0,scoring_CPU_seconds=time.process_time()-cpu,scoring_wall_seconds=time.perf_counter()-wall,wall_scope='Campaign shared-machine elapsed, not guaranteed throughput or active-GPU time'))
 print(json.dumps(summary))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--phase',choices=('seal','score'),required=True);a=p.parse_args();seal() if a.phase=='seal' else score()
