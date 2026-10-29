"""Postprediction two-case physical contrast sensitivity diagnostic only."""
import json,time,shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from score_predictions import metrics,independent
ROOT=Path(__file__).resolve().parents[1]
def main():
 ids=sorted(p.stem for p in (ROOT/'public_data/contrast_control').glob('*.npz'))
 if len(ids)!=2:raise RuntimeError('Exactly two diagnostic observations required')
 dest=ROOT/'experiments/frozen_predictions/contrast_control'
 if dest.exists():raise RuntimeError('Refusing overwrite contrast prediction freeze')
 files=[]
 for ident in ids:
  p=ROOT/'results/contrast_control'/f'{ident}.json';r=json.loads(p.read_text());q=dest/f'{ident}.json';q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q);files.append(dict(observation_id=ident,source=str(p.relative_to(ROOT)),snapshot=str(q.relative_to(ROOT)),bytes=q.stat().st_size))
 freeze=dict(frozen_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),predictions=2,files=files,scope='POSTGT selected two-case sensitivity; not a new blind main test',solver_truth_access='anonymous newY only',new_truth_scored_after_prediction_snapshot=True)
 (ROOT/'audit/EXPERIMENT_PREDICTIONS_FREEZE_CONTRAST_CONTROL.json').write_text(json.dumps(freeze,indent=2)+'\n')
 mapping=json.loads((ROOT/'private_truth/contrast_control/mapping_PRIVATE.json').read_text());lookup={r['observation_id']:r for r in mapping['records']};ops=dict(np.load(ROOT/'public_data/operators.npz'));rows=[];plot=[];validation=[];public=ROOT/'truth_after_freeze/contrast_control';public.mkdir(parents=True,exist_ok=True)
 for ident in ids:
  meta=lookup[ident];truth=np.load(ROOT/'private_truth/contrast_control'/f'{ident}.npz');t=truth['m'];mask=truth['mask'].astype(bool);r=json.loads((dest/f'{ident}.json').read_text());m=np.array(r['m']);Y=np.load(ROOT/'public_data/contrast_control'/f'{ident}.npz')['Y'];met=metrics(m,t,mask);v,_=independent(m,ops,Y);validation.append(v)
  support=m>=1.65;union=np.logical_or(support,mask).sum();row=dict(observation_id=ident,source_observation_id=meta['source_observation_id'],digit=meta['digit'],status=r['status'],RMSE=met['RMSE'],MAE=met['MAE'],foreground_RMS=met['foreground_RMS'],background_RMS=met['background_RMS'],contrast_normalized_RMSE=met['RMSE']/1.3,lowcontrast_support_IoU_at1_65=float(np.logical_and(support,mask).sum()/max(union,1)),main_fixed_threshold5_IoU=met['IoU'],main_IoU5_applicable=False,solver_relative_residual=r['top1']['relative_residual'],all_unknowns_released=r['all_unknowns_released'],wall_seconds=r['wall_seconds'],cpu_seconds=r['cpu_seconds'],accounting=r['accounting'],**v);rows.append(row)
  np.savez_compressed(public/f'{ident}.npz',m=t,mask=mask);plot.append((meta,m,t))
 out=dict(scope='POSTGT physical-regime sensitivity on2publicfailure-selected morphologies; not main recovery/not a new blind benchmark',formula='m_low=1+1.3*(mGT-min)/(max-min), allpixels',source_and_config_frozen_before_new_Y=True,rows=rows,metrics_note='RMSE/MAE primary; low-support IoU threshold1.65 is diagnostic only and cannot replace main predeclared threshold5. Perfect lowGT has allm<5, so main IoU5 is inapplicable',independent_validation=dict(LU=2,explicit_source_RHS=24,zgecon_calls=2,opaque_estimator_internal_RHS='NOT_EXPOSED_NOT_ZERO',seconds=sum(v['independent_seconds'] for v in validation)))
 (ROOT/'results/contrast_control_scores.json').write_text(json.dumps(out,indent=2)+'\n');(public/'mapping_seedless.json').write_text(json.dumps(dict(records=mapping['records'],scope=out['scope']),indent=2)+'\n')
 fig,axes=plt.subplots(2,4,figsize=(10,5),squeeze=False,layout='constrained')
 for i,(meta,m,t) in enumerate(sorted(plot,key=lambda x:x[0]['digit'])):
  high=np.load(ROOT/'truth_after_freeze/holdout'/f'{meta["source_observation_id"]}.npz')['m'];old=json.loads((ROOT/'experiments/frozen_predictions/holdout/fullj_native'/f'{meta["source_observation_id"]}.json').read_text());row=next(r for r in rows if r['observation_id']==meta['observation_id'])
  for j,(a,title) in enumerate([(high,'Original GT'),(np.array(old['m']),'Original custom GN'),(t,'Rescaled GT'),(m,'Same custom GN')]):
   lo,hi=(.5,20) if j<2 else (1,2.3);im=axes[i,j].imshow(a.reshape(24,24),vmin=lo,vmax=hi,cmap='viridis',interpolation='nearest');axes[i,j].set_xticks([]);axes[i,j].set_yticks([])
   if i==0:axes[i,j].set_title(title,fontsize=10)
   if j==0:axes[i,j].set_ylabel(f'Digit {meta["digit"]}')
   if j==3:axes[i,j].set_xlabel(f'RMSE {row["RMSE"]:.3f}; relative range {row["contrast_normalized_RMSE"]:.3f}',fontsize=8)
   fig.colorbar(im,ax=axes[i,j],shrink=.7)
 fig.suptitle('POSTGT sensitivity only: original [0.5,20] colors versus rescaled [1,2.3] colors',fontsize=11);fig.savefig(ROOT/'figures/contrast_regime_sensitivity.png',dpi=170);plt.close(fig)
 print(json.dumps(dict(scope=out['scope'],rows=[{k:r[k] for k in ['digit','RMSE','MAE','contrast_normalized_RMSE','lowcontrast_support_IoU_at1_65','solver_relative_residual','all_unknowns_released','wall_seconds']} for r in rows],validation=out['independent_validation'])))
if __name__=='__main__':main()
