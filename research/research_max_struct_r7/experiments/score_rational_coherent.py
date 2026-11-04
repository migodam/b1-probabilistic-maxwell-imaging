"""Separate postseal scoring of3 coherent rank32 endpoints; no physics calls."""
from pathlib import Path
from datetime import datetime,timezone
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(ROOT/'experiments'))
from score_primary import load,write,csvwrite,resolve,stamp,unchanged,authorize,verify_seal,material_metrics,noise_metrics,public_observation,fixed_prior
CASES=('obs_40ac5b79e8a74ac3','obs_82ac794e228243bb','obs_42cb233e4f224bb3')
def main():
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--authorization',required=True);p.add_argument('--primary-seal',required=True);p.add_argument('--output',default=str(ROOT/'RESULTS/scoring/rational_coherent'));a=p.parse_args()
 authorize(a.authorization);verify_seal(a.primary_seal);tick,cpu=time.perf_counter(),time.process_time();out=Path(a.output)
 if out.exists():raise FileExistsError('Preserve existing coherent scores')
 old12=[stamp(f) for f in sorted((ROOT/'RESULTS/scoring/rational').rglob('*')) if f.is_file()];evidence=[];receipts={}
 for case in CASES:
  folder=ROOT/'RESULTS/rational_coherent'/case;r=load(folder/'RECEIPT.json')
  if r['status']!='SEALED_NO_GT' or r['physical_dimension']!=576 or r['selected']['state_rank']!=32:raise RuntimeError('Incomplete/not rank32 coherent endpoint '+case)
  for name in ['prediction.npz','BUDGET_STATES.npz','CORRECTOR.json']:
   if not(folder/name).is_file():raise RuntimeError('Missing coherent endpoint '+case+'/'+name)
  receipts[case]=r;evidence.extend(stamp(f) for f in sorted(folder.iterdir()) if f.is_file())
 out.mkdir(parents=True);write(out/'PUBLIC_SEAL.json',dict(status='THREE_COHERENT_PUBLIC_ENDPOINTS_SEALED_BEFORE_SCORING',created_utc=datetime.now(timezone.utc).isoformat(),primary_seal=a.primary_seal,files=evidence,original12_scoring_files=old12,truth_access_at_this_seal='NONE',hash_checks='NONE'))
 old=REPO/'research/research_max_struct_r6';mapping={r['case_id']:r for r in load(old/'sealed_truth/scoring/historical_mapping_PRIVATE.json')+load(old/'sealed_truth/confirmation/mapping_PRIVATE.json')['records']};cases={r['case_id']:r for r in load(old/'data_public/cases.json')['records']};tasks=load(old/'configs/primary_final_manifest.json')['tasks'];rows=[];pictures={}
 for case in CASES:
  folder=ROOT/'RESULTS/rational_coherent'/case;r=receipts[case];meta=mapping[case];task=next(t for t in tasks if t['case_id']==case and t['arm']=='cie_dual');public=public_observation(task)
  with np.load(resolve(meta['truth_path']),allow_pickle=False) as z:truth=z['m'].ravel().copy();support=z['mask'].ravel().astype(bool)
  with np.load(folder/'prediction.npz',allow_pickle=False) as z:pred={k:z[k].copy() for k in z.files}
  with np.load(folder/'BUDGET_STATES.npz',allow_pickle=False) as z:finite=z['m'][0].copy()
  with np.load(resolve(task['output_dir'])/'prediction.npz',allow_pickle=False) as z:anchor=z['m'].ravel().copy()
  m=pred['m'].ravel();met=material_metrics(m,truth,support,cases[case]['amplitude']);noise=noise_metrics(pred['Y_prediction'],public,True);phi=.5*noise['relative_residual']**2+fixed_prior(m);corr=load(folder/'CORRECTOR.json')
  row=dict(case_id=case,endpoint='rank32_coherent_corrected',classification='DISCLOSED_OLD_FAILED_ANCHOR_NOT_NEW_CONFIRMATION',state_rank=32,amplitude_dimension=r['selected']['amplitude_dimension'],material_cells=int(m.size),**met,**noise,anchor_material_metrics=material_metrics(anchor,truth,support,cases[case]['amplitude']),finite_endpoint_material_metrics=material_metrics(finite,truth,support,cases[case]['amplitude']),anchor_original_Phi=r['original_Phi'],finite_endpoint_original_Phi=r['best_finite_Phi'],corrected_Phi_receipt=r['corrected_Phi'],corrected_Phi_recomputed_from_saved_Y=phi,Phi_receipt_difference=phi-r['corrected_Phi'],residual_receipt_difference=noise['relative_residual']-r['physical_relative_residual'],full_physical_dimension=r['physical_dimension'],experiment_physical_counts=r['physical_counts'],experiment_cpu_seconds=r['cpu_seconds'],experiment_wall_seconds=r['wall_seconds'],corrector_stop=corr.get('stop'),corrector_success=corr.get('success'),new_scoring_physical_calls=0)
  rows.append(row);pictures[case]=(truth,anchor,m)
 write(out/'SCORES.json',dict(rows=rows,endpoints=3,material_passes=sum(r['material_image_success'] for r in rows),independent_new_objects=0,scope='Only3 coherent rank32 endpoints. Original12 separate endpoints include continued-anchor and local-search controls; combined15 must not be called15 new rational proposals.',selection='Score fixed original-Phi-chosen coherent endpoint after existing parent ALL_PRIMARY_SEALED authorization; no GT candidate selection',new_physical_calls=0));csvwrite(out/'SCORES.csv',rows)
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axes=plt.subplots(3,3,figsize=(9,9),layout='constrained')
 for i,case in enumerate(CASES):
  truth,anchor,m=pictures[case];vmax=max(float(np.max(x)) for x in [truth,anchor,m])
  for j,(label,x) in enumerate(zip(['GT','Old sealed Dual anchor','Rank32 coherent corrected'],[truth,anchor,m])):
   axes[i,j].imshow(x.reshape(24,24),vmin=.02,vmax=vmax,cmap='viridis');axes[i,j].set_xticks([]);axes[i,j].set_yticks([]);axes[i,j].set_title(label+'\n'+case if j==0 else label+'\nNRMSE %.3f'%material_metrics(x,truth,truth>.5,cases[case]['amplitude'])['material_nrmse'],fontsize=9)
 fig.suptitle('Three separately scored coherent rank32 endpoints\nFull576-variable maps; disclosed old objects, no new independent confirmation',fontsize=11);fig.savefig(out/'FULL_MATERIAL_MAPS.png',dpi=160);plt.close(fig)
 if not all(unchanged(s) for s in evidence+old12):raise RuntimeError('Public endpoints or original12 scoring files changed during scoring')
 if max(abs(r['Phi_receipt_difference']) for r in rows)>1e-12:raise RuntimeError('Saved-Y Phi mismatch')
 write(out/'RECEIPT.json',dict(status='COMPLETE',created_utc=datetime.now(timezone.utc).isoformat(),authorization='Parent ALL_PRIMARY_SEALED remains authorized; new explicit coherent3 scoring task',rows=3,material_passes=sum(r['material_image_success'] for r in rows),coherent_predictions_and_Phi_unchanged=True,original12_seal_and_results_unchanged=True,new_physical_calls=0,new_physical_LU=0,new_physical_RHS=0,wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-cpu,hash_checks='NONE'))
 print(json.dumps([dict(case_id=r['case_id'],NRMSE=r['material_nrmse'],IoU=r['otsu_iou'],material_pass=r['material_image_success'],Phi_gap=r['Phi_receipt_difference']) for r in rows]))
if __name__=='__main__':main()
