"""Public-only32pair seal then separate postseal GT scoring of budget controls."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time,argparse,csv
from collections import Counter
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];BASE=ROOT/'delegated/diagnostics/bothblocks_confirmation';OUT=BASE/'posthoc_budget_control';sys.path.insert(0,str(REPO))
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else v.tolist())+'\n')
def csvwrite(p,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)
def seal():
 wall,cpu=time.perf_counter(),time.process_time();m=json.loads((OUT/'FROZEN_MANIFEST.json').read_text());q=json.loads((OUT/'QUEUE.json').read_text());jobs={j['pair_id']:j for j in q['jobs']}
 if q['status'] not in ('COMPLETED','COMPLETED_WITH_ERRORS') or q['pending'] or len(jobs)!=sum(r['trigger'] for r in m['pairs']):raise RuntimeError('All triggers must terminate before any GT read')
 rows=[];ms=[];ys=[]
 for r in m['pairs']:
  x=dict(r);job=jobs.get(r['pair_id']);extra=0.
  if r['trigger']:
   d=Path(r['output_dir']);receipt=json.loads((d/'RECEIPT.json').read_text());extra=job['new_child_actual_CPU_seconds'];source=receipt['selected_prediction'];phi=receipt['selected_original_Phi'];selected_kind=receipt['selected_kind'];eligible=receipt['selected_cumulative_eligible_CPU_seconds']
   if selected_kind=='NEW_PAID_ELIGIBLE_LOWER_PHI' and not (eligible<=r['target_G19_actual_CPU_seconds'] and phi<r['original_Phi']):raise RuntimeError('Illegal late/nonimproving candidate; never backfill')
   x.update(continuation_status=receipt['status'],continuation_corrector_stop=json.loads((d/'CORRECTOR.json').read_text())['stop'] if (d/'CORRECTOR.json').exists() else 'CPU_INTERRUPT_NO_RETURN',new_physical_cost=receipt['new_physical_cost'])
  else:source=r['original_prediction'];phi=r['original_Phi'];selected_kind='ORIGINAL_RETAINED_NO_TRIGGER';eligible=r['original_actual_child_CPU_seconds'];x.update(continuation_status='NOT_RUN',new_physical_cost={})
  total=r['original_actual_child_CPU_seconds']+extra
  x.update(selected_prediction_source=source,selected_original_Phi=phi,selected_kind=selected_kind,selected_cumulative_eligible_CPU_seconds=eligible,new_child_actual_CPU_seconds=extra,total_actual_CPU_seconds=total,unused_CPU_seconds=max(0.,r['target_G19_actual_CPU_seconds']-total),charged_new_atomic_or_exit_overrun_CPU_seconds=max(0.,extra-r['remaining_new_child_CPU_allowance_seconds']) if r['trigger'] else 0.,original_CPU_already_above_reference_seconds=max(0.,r['original_actual_child_CPU_seconds']-r['target_G19_actual_CPU_seconds']),material_index=len(ms),selection_GT_access='NONE_THIS_SELECTION',cost_scope='Actual old+new CPU charged; no idle-fill, no assertion of equal actual spend')
  with np.load(source) as z:ms.append(z['m'].copy());ys.append(z['Y_prediction'].copy())
  rows.append(x)
 np.savez_compressed(OUT/'ALL32_SELECTED_PREDICTIONS.npz',m=np.asarray(ms),Y_prediction=np.asarray(ys))
 dump(OUT/'ALL32_SELECTIONS_BEFORE_GT_SCORE.json',dict(status='ALL_TRIGGERED_TERMINAL_AND_ALL32_PUBLIC_PHI_SELECTIONS_SEALED',classification='POSTHOC_BUDGET_CONTROL',rows=rows,original80_and48_unchanged=True,GT_access='NONE_THIS_PHASE',new_forward=0,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall));print(json.dumps(dict(status='SEALED_BEFORE_GT',pairs=32,triggered=len(jobs))))
def score():
 wall,cpu=time.perf_counter(),time.process_time();sel=json.loads((OUT/'ALL32_SELECTIONS_BEFORE_GT_SCORE.json').read_text())
 if sel['status']!='ALL_TRIGGERED_TERMINAL_AND_ALL32_PUBLIC_PHI_SELECTIONS_SEALED':raise RuntimeError('All32 seal required')
 from research.research_max_struct_r10.experiments.score_twod_comparison import material_metrics
 original=json.loads((BASE/'SINGLE_ARM_SCORES.json').read_text());ix={(x['case_id'],x['method']):x for x in original};truth={};public={};rows=[];newpc=Counter()
 for cid in sorted({r['case_id'] for r in sel['rows']}):
  with np.load(ROOT/'sealed_truth/bothblocks16'/f'{cid}.npz') as z:truth[cid]=z['m'].copy()
  with np.load(ROOT/'data_public/bothblocks16'/f'{cid}.npz') as z:public[cid]=dict(Y=z['Y'].copy(),sigma=float(z['sigma']))
 with np.load(OUT/'ALL32_SELECTED_PREDICTIONS.npz') as z:ms=z['m'];ys=z['Y_prediction']
 for r in sel['rows']:
  cid,method=r['case_id'],r['method'];old=ix[cid,method];g=ix[cid,'fullG19'];gt=truth[cid];met=material_metrics(ms[r['material_index']],gt,float(gt.max()-.1),old['high_chi']);cn=float(np.linalg.norm(ys[r['material_index']]-public[cid]['Y'])/public[cid]['sigma']);bound=float(np.sqrt(chi2.ppf(.99,2*public[cid]['Y'].size)/2));newpc.update(r['new_physical_cost']);rows.append(dict(r,family=old['family'],regime=old['regime'],original_material_PASS=old['material_success'],original_amplitude_NRMSE=old['amplitude_NRMSE'],original_fixed_IoU=old['fixed_support_IoU'],controlled_material_PASS=met['material_success'],controlled_amplitude_NRMSE=met['amplitude_NRMSE'],controlled_fixed_IoU=met['fixed_support_IoU'],controlled_CN99_pass=bool(cn<=bound),controlled_whitened_residual=cn,CN99_threshold=bound,delta_NRMSE_control_minus_original=met['amplitude_NRMSE']-old['amplitude_NRMSE'],G19_material_PASS=g['material_success'],G19_amplitude_NRMSE=g['amplitude_NRMSE'],G19_fixed_IoU=g['fixed_support_IoU'],delta_NRMSE_control_minus_G19=met['amplitude_NRMSE']-g['amplitude_NRMSE'],GT_access='Only after all triggers and all32 selections sealed'))
 dump(OUT/'ALL32_POSTSEAL_SCORES.json',rows);csvwrite(OUT/'ALL32_POSTSEAL_SCORES.csv',rows);summary=[]
 for method in ('cie_single','cie_dual'):
  rr=[r for r in rows if r['method']==method];summary.append(dict(method=method,objects=16,triggered=sum(r['trigger'] for r in rr),original_material_PASS=sum(r['original_material_PASS'] for r in rr),controlled_material_PASS=sum(r['controlled_material_PASS'] for r in rr),controlled_CN99_pass=sum(r['controlled_CN99_pass'] for r in rr),mean_original_NRMSE=float(np.mean([r['original_amplitude_NRMSE'] for r in rr])),mean_controlled_NRMSE=float(np.mean([r['controlled_amplitude_NRMSE'] for r in rr])),mean_G19_NRMSE=float(np.mean([r['G19_amplitude_NRMSE'] for r in rr])),material_gate_improvements=sum(r['controlled_material_PASS'] and not r['original_material_PASS'] for r in rr),material_gate_regressions=sum(r['original_material_PASS'] and not r['controlled_material_PASS'] for r in rr),new_actual_child_CPU_seconds=sum(r['new_child_actual_CPU_seconds'] for r in rr),total_control_actual_CPU_seconds=sum(r['total_actual_CPU_seconds'] for r in rr),G19_reference_actual_CPU_seconds_sum=sum(r['target_G19_actual_CPU_seconds'] for r in rr),unused_CPU_allowance_seconds_sum=sum(r['unused_CPU_seconds'] for r in rr),charged_new_overrun_seconds_sum=sum(r['charged_new_atomic_or_exit_overrun_CPU_seconds'] for r in rr),not_equal_actual_spend=True))
 dump(OUT/'SUMMARY_BY_METHOD.json',summary);csvwrite(OUT/'SUMMARY_BY_METHOD.csv',summary);q=json.loads((OUT/'QUEUE.json').read_text());dump(OUT/'COST_LEDGER.json',dict(status='COMPLETE_POSTHOC_BUDGET_CONTROL',pairs=32,actual_continuation_children=len(q['jobs']),failed_children=sum(j['status']=='FAILED' for j in q['jobs']),new_exact_child_CPU_seconds=q['new_child_actual_CPU_seconds_sum'],controller_CPU_seconds=q['controller_CPU_seconds'],new_physical_cost=dict(newpc),original80_cost_added_again=False,original80_and48_unchanged=True,selection_and_scoring_new_forward=0,scoring_CPU_seconds=time.process_time()-cpu,scoring_wall_seconds=time.perf_counter()-wall,cost_scope='Only7 new continuation children are new campaign spend; original32 and G19 references are deployment views of already paid80, not additive'))
 print(json.dumps(summary))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--phase',choices=('seal','score'),required=True);a=p.parse_args();seal() if a.phase=='seal' else score()
