"""All80-terminal postseal scoring, paired object bootstrap and declared portfolios."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,csv,time
from collections import Counter,defaultdict
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];OUT=ROOT/'delegated/diagnostics/bothblocks_confirmation'
sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.experiments.score_twod_comparison import material_metrics
METHODS=('cie_single','cie_dual','global_dual','fullG19','fullGreenBothBlocks')
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else v.tolist())+'\n')
def csvwrite(p,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
  for r in rows:w.writerow({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()})
def main():
 wall,cpu=time.perf_counter(),time.process_time();seal=json.loads((OUT/'ALL80_TERMINAL_PREDICTION_SEAL.json').read_text());portseal=json.loads((OUT/'PORTFOLIO_SELECTION_BEFORE_GT.json').read_text());q=json.loads((OUT/'QUEUE.json').read_text());sc=json.loads((ROOT/'configs/twod_bothblocks16_scoring.json').read_text())
 if seal['status']!='ALL80_TERMINAL_BEFORE_ANY_CONFIRMATION_GT_SCORING' or len(seal['entries'])!=80 or portseal['status']!='SEALED_BEFORE_ANY_CONFIRMATION_GT_READ' or q['status'] not in ('COMPLETED','COMPLETED_WITH_ERRORS'):raise RuntimeError('All80 prediction/portfolio seals required; no truth opened')
 if sc['bootstrap_seed']!=202611071117 or sc['bootstrap_replicates']!=10000:raise RuntimeError('Frozen bootstrap differs')
 cfg=json.loads((ROOT/'configs/twod_bothblocks16_generator.json').read_text());meta={r['case_id']:r for r in cfg['cases']};truths={};pubs={};rows=[];prefix=[];physical=Counter();aux=Counter();armcpu=defaultdict(float)
 # First actual private read occurs only after both all80 and portfolio seals above.
 for cid in meta:
  with np.load(ROOT/'sealed_truth/bothblocks16'/f'{cid}.npz') as z:truths[cid]=z['m'].copy()
  with np.load(ROOT/'data_public/bothblocks16'/f'{cid}.npz') as z:pubs[cid]=dict(Y=z['Y'],sigma=float(z['sigma']))
 for e in seal['entries']:
  t=e['task'];cid=t['case_id'];method=t['method'];d=REPO/t['output_dir'];rec=json.loads((d/'RECEIPT.json').read_text()) if (d/'RECEIPT.json').exists() else {};gt=truths[cid];amp=float(gt.max()-.1);pc=rec.get('physical_model_cost',{});ac=(rec.get('auxiliary_execution_receipt') or {}).get('operation_counts',{})
  physical.update(pc);aux.update({k:v for k,v in ac.items() if not any(tag in k for tag in ('seconds','microseconds','bytes'))});armcpu[method]+=e['terminal_job']['child_process_CPU_seconds']
  row=dict(case_id=cid,method=method,family=meta[cid]['family'],regime=meta[cid]['regime'],high_chi=meta[cid]['high_chi'],status=rec.get('status','MISSING_RECEIPT'),prediction_available=e['cap_prediction_available'],material_success=False,
   original_Phi=e['cap_original_Phi'],external_complete_child_CPU_seconds=e['terminal_job']['child_process_CPU_seconds'],near_end_whole_child_CPU_seconds=rec.get('whole_child_CPU_seconds'),CPU_cap=480,CPU_overrun_seconds=rec.get('CPU_overrun_seconds'),
   selected_auxiliary_stage=rec.get('auxiliary_selected_stage'),corrector_stop=(rec.get('corrector') or {}).get('stop'),corrector_iterations=(rec.get('corrector') or {}).get('iterations'),physical_cost=pc,auxiliary_cost=ac,truth_access='POST_ALL80_TERMINAL_AND_PORTFOLIO_SEAL')
  if e['cap_prediction_available']:
   with np.load(e['prediction_path']) as z:m=z['m'].copy();pred=z['Y_prediction'].copy()
   metrics=material_metrics(m,gt,amp,meta[cid]['high_chi']);Y,sigma=pubs[cid]['Y'],pubs[cid]['sigma'];cn=float(np.linalg.norm(pred-Y)/sigma);bound=float(np.sqrt(chi2.ppf(.99,2*Y.size)/2))
   row.update(**metrics,whitened_residual=cn,CN99_threshold=bound,data_CN99_pass=bool(cn<=bound))
  rows.append(row)
  cps=json.loads((d/'R11_CPU_CHECKPOINTS.json').read_text()) if (d/'R11_CPU_CHECKPOINTS.json').exists() else []
  cpmap={x['CPU_cap_seconds']:x for x in cps}
  if (d/'R11_CHECKPOINT_MATERIALS.npz').exists():
   with np.load(d/'R11_CHECKPOINT_MATERIALS.npz') as z:ms=z['m']
  else:ms=None
  for cap in sc['CPU_prefixes']:
   x=cpmap.get(cap,{});p=dict(case_id=cid,method=method,family=meta[cid]['family'],regime=meta[cid]['regime'],CPU_allowed_cap=cap,endpoint_status=rec.get('status','MISSING_RECEIPT'),material_success=False,
    actual_natural_end_child_CPU_seconds=row['external_complete_child_CPU_seconds'] if row['external_complete_child_CPU_seconds']<=cap else None,
    actual_prefix_spend_status='EXACT_NATURAL_COMPLETION' if row['external_complete_child_CPU_seconds']<=cap else 'SAVED_PAID_PREFIX; INTERRUPTED_EXECUTION_SPEND_NOT_OBSERVED')
   if 'material_index' in x and ms is not None:p.update(status='AVAILABLE',original_Phi=x['original_Phi'],paid_state_CPU_seconds=x['state_paid_CPU_seconds'],phase=x['phase'],**material_metrics(ms[x['material_index']],gt,amp,meta[cid]['high_chi']))
   else:p.update(status='NO_ELIGIBLE_PAID_STATE')
   prefix.append(p)
 dump(OUT/'SINGLE_ARM_SCORES.json',rows);csvwrite(OUT/'SINGLE_ARM_SCORES.csv',rows);dump(OUT/'CPU_PREFIX_SCORES.json',prefix);csvwrite(OUT/'CPU_PREFIX_SCORES.csv',prefix)
 ports=[]
 with np.load(OUT/'PORTFOLIO_SEALED_MATERIALS.npz') as z:pm=z['m']
 for p in portseal['rows']:
  r=dict(p,regime=meta[p['case_id']]['regime'],family=meta[p['case_id']]['family'],material_success=False,truth_access='POST_PUBLIC_PHI_SELECTION_SEAL')
  if p.get('material_index') is not None:
   gt=truths[p['case_id']];r.update(**material_metrics(pm[p['material_index']],gt,float(gt.max()-.1),meta[p['case_id']]['high_chi']))
   if p['view']=='completed_480_each':
    chosen=next(x for x in rows if x['case_id']==p['case_id'] and x['method']==p['selected_method']);r.update(data_CN99_pass=chosen.get('data_CN99_pass'),whitened_residual=chosen.get('whitened_residual'))
  ports.append(r)
 dump(OUT/'PORTFOLIO_SCORES.json',ports);csvwrite(OUT/'PORTFOLIO_SCORES.csv',ports)
 summary=[]
 for method in METHODS:
  for regime in ('all','high','mixed'):
   rs=[x for x in rows if x['method']==method and (regime=='all' or x['regime']==regime)];valid=[x for x in rs if x.get('amplitude_NRMSE') is not None]
   summary.append(dict(method=method,regime=regime,objects=len(rs),predictions=len(valid),material_passes=sum(x['material_success'] for x in rs),CN99_passes=sum(x.get('data_CN99_pass',False) for x in rs),
    mean_amplitude_NRMSE=float(np.mean([x['amplitude_NRMSE'] for x in valid])) if valid else None,median_amplitude_NRMSE=float(np.median([x['amplitude_NRMSE'] for x in valid])) if valid else None,
    median_fixed_IoU=float(np.median([x['fixed_support_IoU'] for x in valid])) if valid else None,actual_child_CPU_seconds_sum=sum(x['external_complete_child_CPU_seconds'] for x in rs)))
 dump(OUT/'SUMMARY_BY_ARM_REGIME.json',summary);csvwrite(OUT/'SUMMARY_BY_ARM_REGIME.csv',summary)
 # Same resampled object indices for all paired quantities; stratify eight high/eight mixed.
 cases=sorted(meta);idx={c:i for i,c in enumerate(cases)};rng=np.random.default_rng(sc['bootstrap_seed']);strata=[[idx[c] for c in cases if meta[c]['regime']==reg] for reg in ('high','mixed')]
 boot=np.concatenate([rng.choice(s,size=(sc['bootstrap_replicates'],len(s)),replace=True) for s in strata],axis=1)
 def paired(a,b,label):
  aa={x['case_id']:x for x in a};bb={x['case_id']:x for x in b};valid=[c for c in cases if aa[c].get('amplitude_NRMSE') is not None and bb[c].get('amplitude_NRMSE') is not None]
  result=dict(comparison=label,objects=16,complete_metric_pairs=len(valid),bootstrap_unit='object; high/mixed stratified',bootstrap_replicates=10000,bootstrap_seed=sc['bootstrap_seed'])
  if len(valid)==16:
   delta=np.array([aa[c]['amplitude_NRMSE']-bb[c]['amplitude_NRMSE'] for c in cases]);bs=delta[boot].mean(axis=1)
   result.update(mean_amplitude_NRMSE_difference=float(delta.mean()),mean_difference_CI95=np.quantile(bs,[.025,.975]).tolist(),median_amplitude_NRMSE_difference=float(np.median(delta)),NRMSE_wins=int(np.count_nonzero(delta<0)),NRMSE_ties=int(np.count_nonzero(delta==0)),NRMSE_losses=int(np.count_nonzero(delta>0)))
  sd=np.array([int(aa[c].get('material_success',False))-int(bb[c].get('material_success',False)) for c in cases]);bs=sd[boot].mean(axis=1)
  result.update(material_success_rate_difference=float(sd.mean()),success_rate_difference_CI95=np.quantile(bs,[.025,.975]).tolist(),material_gate_improvements=int(np.count_nonzero(sd>0)),material_gate_regressions=int(np.count_nonzero(sd<0)))
  return result
 comparisons=[]
 for baseline in ('cie_single','cie_dual','global_dual','fullG19'):
  comparisons.append(paired([x for x in rows if x['method']=='fullGreenBothBlocks'],[x for x in rows if x['method']==baseline],'fullGreenBothBlocks minus '+baseline))
 for view in ('completed_480_each','prefix240_each'):
  comparisons.append(paired([x for x in ports if x['portfolio']=='CIE_plus_BothBlocks' and x['view']==view],[x for x in ports if x['portfolio']=='CIE_plus_Dual' and x['view']==view],'portfolio CIE+Both minus CIE+Dual '+view))
  cie=[x for x in prefix if x['method']=='cie_single' and x['CPU_allowed_cap']==480]
  comparisons.append(paired([x for x in ports if x['portfolio']=='CIE_plus_BothBlocks' and x['view']==view],cie,'portfolio CIE+Both '+view+' minus CIE480'))
 dump(OUT/'PAIRED_UNCERTAINTY.json',comparisons);csvwrite(OUT/'PAIRED_UNCERTAINTY.csv',comparisons)
 generation=json.loads((OUT/'GENERATION_RECEIPT.json').read_text());check=json.loads((ROOT/'delegated/diagnostics/BOTHBLOCKS_WIRING_CHECK.json').read_text())
 dump(OUT/'COST_LEDGER.json',dict(status='COMPLETE',inversions=80,failed_jobs=sum(j['status']!='FINISHED' for j in q['jobs']),exact_inversion_child_CPU_seconds=q['exact_inversion_child_CPU_seconds_sum'],exact_generation_child_CPU_seconds=q['generation_job']['external_child_CPU_seconds'],controller_CPU_seconds=q['controller_CPU_seconds'],campaign_shared_wall_seconds=q['wall_seconds'],per_method_child_CPU_seconds=dict(armcpu),physical_inversion_cost=dict(physical),auxiliary_inversion_cost_count_sums=dict(aux),auxiliary_scope='Per-run full counters retained in SINGLE_ARM_SCORES; aggregated counts exclude timers and memory-byte peaks; nested timers are not added to child CPU',generation_physical_cost=generation['physical_cost'],
  scoring_CPU_seconds=time.process_time()-cpu,scoring_wall_seconds=time.perf_counter()-wall,scoring_new_FW=0,scoring_new_adjoint=0,portfolio_extra_solver_runs=0,portfolio_cost_scope='Each deployment sums both paid candidates. Campaign80 contains shared CIE once; never add overlapping portfolio costs to campaign total.',hardware='Same Mac two lanes BLAS1 for all5; raw CPU never pooled with Windows48',bootstrap_replicates=10000,GT_scope='First scoring read only after all80 terminal and portfolio seal'))
 dump(OUT/'SCORING_RECEIPT.json',dict(status='COMPLETE',single_rows=len(rows),prefix_rows=len(prefix),portfolio_rows=len(ports),comparisons=len(comparisons),physical_forward=0,physical_adjoint=0,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall))
 print(json.dumps(dict(status='COMPLETE',summary=[x for x in summary if x['regime']=='all'])))
if __name__=='__main__':main()
