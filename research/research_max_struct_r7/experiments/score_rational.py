"""Post-primary-seal offline scoring of disclosed rational-action anchors.
No GT online; no new physical calls. R6 mapping is opened only after explicit
ALL_PRIMARY_SEALED authorization and verification of all180 primary outputs.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,time
import numpy as np
try:from .score_primary import load,write,csvwrite,resolve,stamp,authorize,verify_seal,material_metrics,noise_metrics,public_observation
except ImportError:from score_primary import load,write,csvwrite,resolve,stamp,authorize,verify_seal,material_metrics,noise_metrics,public_observation
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
CASES=('obs_40ac5b79e8a74ac3','obs_82ac794e228243bb','obs_42cb233e4f224bb3')
LABELS=('anchor_continuation','rank2_finite','rank8_finite','local_exact_origin')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--primary-seal',required=True);p.add_argument('--authorization',required=True);p.add_argument('--output',required=True);p.add_argument('--rational-root',default=str(ROOT/'RESULTS/rational'));a=p.parse_args()
 authorize(a.authorization);verify_seal(a.primary_seal);start,cpu=time.perf_counter(),time.process_time();out=Path(a.output);rational=Path(a.rational_root)
 if out.exists():raise FileExistsError('Write-once rational scoring output exists')
 # Verify every requested public endpoint exists before any R6 truth/mapping open.
 sealed=[]
 for case in CASES:
  root=rational/case;receipt=load(root/'RECEIPT.json')
  if receipt['status']!='SEALED_NO_GT':raise RuntimeError('Incomplete rational action case '+case)
  for label in LABELS:
   folder=root/label
   for name in ['prediction.npz','BUDGET_STATES.npz','RECEIPT.json','CORRECTOR.json']:
    if not(folder/name).is_file():raise RuntimeError('Missing sealed rational endpoint '+str(folder/name))
  sealed.extend(stamp(x) for x in sorted(root.rglob('*')) if x.is_file())
 oldroot=REPO/'research/research_max_struct_r6';cases={r['case_id']:r for r in load(oldroot/'data_public/cases.json')['records']}
 historical=load(oldroot/'sealed_truth/scoring/historical_mapping_PRIVATE.json');fresh=load(oldroot/'sealed_truth/confirmation/mapping_PRIVATE.json')['records'];mapping={r['case_id']:r for r in historical+fresh}
 tasks=load(oldroot/'configs/primary_final_manifest.json')['tasks'];rows=[];pictures={};truths={};costs=[]
 for case in CASES:
  meta=mapping[case];c=cases[case];task=next(t for t in tasks if t['case_id']==case and t['arm']=='cie_dual');public=public_observation(task)
  with np.load(resolve(meta['truth_path']),allow_pickle=False) as z:truth=z['m'].ravel().copy();support=z['mask'].ravel().astype(bool)
  truths[case]=truth
  with np.load(resolve(task['output_dir'])/'prediction.npz',allow_pickle=False) as z:anchor=z['m'].ravel().copy()
  anchor_metrics=material_metrics(anchor,truth,support,c['amplitude']);pictures[(case,'R6_sealed_anchor')]=anchor
  base=load(rational/case/'RECEIPT.json');costs.append(dict(case_id=case,cpu_seconds=base['cpu_seconds'],wall_seconds=base['wall_seconds'],physical_counts=base['physical_counts'],small_objective_queries=base['small_objective_queries'],setup_cpu_seconds=base['setup_cpu_seconds'],setup_wall_seconds=base['setup_wall_seconds'],status=base['status']))
  for label in LABELS:
   folder=rational/case/label;receipt=load(folder/'RECEIPT.json')
   with np.load(folder/'prediction.npz',allow_pickle=False) as z:final={k:z[k].copy() for k in z.files}
   with np.load(folder/'BUDGET_STATES.npz',allow_pickle=False) as z:endpoint=z['m'][0].copy()
   pre=material_metrics(endpoint,truth,support,c['amplitude']);post=material_metrics(final['m'],truth,support,c['amplitude']);noise=noise_metrics(final['Y_prediction'],public,True)
   source=receipt.get('source_action') or {};row=dict(case_id=case,label=label,classification='EXPOSED_R6_FAILED_ANCHOR_DEVELOPMENT_NOT_NEW_CONFIRMATION',shape=meta.get('shape',meta.get('family',case)),anchor_metrics=anchor_metrics,finite_endpoint_material_metrics=pre,corrected_material_metrics=post,**noise,anchor_original_Phi=base['anchor_original_Phi'],finite_endpoint_original_Phi=receipt['finite_endpoint_original_Phi'],corrected_original_Phi=receipt['corrected_original_Phi'],material_NRMSE_improvement_over_anchor=anchor_metrics['material_nrmse']-post['material_nrmse'],physical_cost=receipt['cost'],stop=receipt['stop'],action_kind=source.get('kind'),state_rank=source.get('state_rank'),path_max_above_origin=source.get('path_max_above_origin'),finite_action_material_precision='initial saved budget incumbent float32; final prediction float64',same_old_joint_material_CN99=bool(post['material_image_success'] and noise['noise_CN99_pass']))
   rows.append(row);pictures[(case,label)]=final['m'].ravel()
 out.mkdir(parents=True);write(out/'RATIONAL_PUBLIC_SEAL.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),files=sealed,truth_access_at_seal='NONE',hash_checks='NONE'));write(out/'RATIONAL_SCORES.json',dict(rows=rows,anchor_cases=3,independent_confirmation_objects=0,selection='Top endpoints chosen and fully corrected with original publicPhi before independent scorer accessedGT; no scorer-selected substitute endpoint',physical_scoring_calls=0));csvwrite(out/'RATIONAL_SCORES.csv',rows);write(out/'RATIONAL_COSTS.json',dict(rows=costs,scope='One whole sequential search/corrector process percase. Endpoint budgetCPU snapshots include earlier search and earlier endpoints; not separate comparable optimizationCPU. Full costs include anchor setup/small queries/validation and correction.',extra_scoring_LU=0,extra_scoring_RHS=0))
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axes=plt.subplots(3,6,figsize=(16,8),layout='constrained')
 for i,case in enumerate(CASES):
  amp=cases[case]['amplitude'];axes[i,0].imshow(truths[case].reshape(24,24),vmin=.02,vmax=.1+amp,cmap='viridis');axes[i,0].set_title(case+'\nR6 GT',fontsize=8)
  for j,label in enumerate(('R6_sealed_anchor',)+LABELS,1):
   m=pictures[(case,label)];axes[i,j].imshow(m.reshape(24,24),vmin=.02,vmax=.1+amp,cmap='viridis')
   if label=='R6_sealed_anchor':nrmse=next(r['anchor_metrics']['material_nrmse'] for r in rows if r['case_id']==case)
   else:nrmse=next(r['corrected_material_metrics']['material_nrmse'] for r in rows if r['case_id']==case and r['label']==label)
   axes[i,j].set_title(label+'\nNRMSE %.3f'%nrmse,fontsize=8)
  for ax in axes[i]:ax.set_xticks([]);ax.set_yticks([])
 fig.suptitle('Disclosed rational-action development anchors: all four corrected endpoints; not new independent holdout')
 fig.savefig(ROOT/'FIGURES/primary/rational_disclosed_endpoints.png',dpi=160);plt.close(fig)
 write(out/'SCORE_RECEIPT.json',dict(status='COMPLETE',created_utc=datetime.now(timezone.utc).isoformat(),primary_seal=a.primary_seal,authorization='ALL_PRIMARY_SEALED',truth_access='R6 mapping/scorer-only after complete180primary seal',rows=len(rows),new_physical_calls=0,new_physical_LU=0,new_physical_RHS=0,wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu))
 print(json.dumps(dict(status='RATIONAL_SCORED',rows=len(rows),out=str(out))))
if __name__=='__main__':main()
