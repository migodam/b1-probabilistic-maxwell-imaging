"""Pure post-score aggregation; never reads private truth or modifies frozen files."""
import json,time
from collections import Counter
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
METHODS=['fullj_native','fullj_reciprocal','globaldual']
LABELS=['Full-J native','Full-J reciprocal','Global dual + Full-J']
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def main():
 freeze=json.loads((ROOT/'audit/HOLDOUT_FREEZE.json').read_text());checks=[]
 for f in freeze['files']:
  live=(ROOT/f['path']).read_bytes();saved=(ROOT/f['snapshot']).read_bytes();checks.append(dict(path=f['path'],current_bytes=len(live),snapshot_bytes=len(saved),byte_equal=live==saved))
 dump(ROOT/'audit/EXPERIMENT_SOURCE_EQUALITY_HOLDOUT.json',dict(checked_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),method='full file bytes, no hashes',files=checks,all_16_unchanged=len(checks)==16 and all(c['byte_equal'] for c in checks)))
 if not all(c['byte_equal'] for c in checks):raise RuntimeError('Frozen source/config/protocol drift')
 scores=json.loads((ROOT/'results/holdout_scores.json').read_text());lookup={(r['observation_id'],r['method']):r for r in scores['rows']};rows=[]
 iou_reference=[]
 for ident in sorted(set(r['observation_id'] for r in scores['rows'])):
  truth=np.load(ROOT/'truth_after_freeze/holdout'/f'{ident}.npz');material=truth['m'];mask=truth['mask'].astype(bool);physical_mask=material>=5.
  iou_reference.append(dict(observation_id=ident,digit=next(r['digit'] for r in scores['rows'] if r['observation_id']==ident),GT_material_threshold5_self_IoU=float(np.logical_and(physical_mask,mask).sum()/np.logical_or(physical_mask,mask).sum()),truth_alpha_mask_pixels=int(mask.sum()),GT_material_threshold5_pixels=int(physical_mask.sum())))
 dump(ROOT/'results/holdout_iou_reference.json',dict(auxiliary_only=True,computed_after_all_prediction_freeze=True,definition='IoU of exact GT material>=5 against stored truth alpha>=0.5 mask',interpretation='reference only; not a global upper bound, model score, acceptance gate, or alteration of predeclared metrics',rows=iou_reference,mean_GT_self_IoU=float(np.mean([r['GT_material_threshold5_self_IoU'] for r in iou_reference]))))
 for p in (ROOT/'experiments/frozen_predictions/holdout').glob('*/*.json'):
  r=json.loads(p.read_text());s=lookup[(p.stem,p.parent.name)];a=r.get('accounting',{});c=a.get('counts',{});stages=r.get('stages',[])
  rows.append(dict(observation_id=p.stem,digit=s['digit'],method=p.parent.name,status=r['status'],RMSE=s['RMSE'],IoU=s['IoU'],SSIM=s['SSIM'],wall_seconds=r.get('wall_seconds',0),cpu_seconds=r.get('cpu_seconds',0),full_factorizations=a.get('full_factorizations',0),physical_complex_RHS=a.get('full_complex_rhs',0),source_RHS=a.get('source_rhs',0),source_adjoint_RHS=a.get('adjoint_rhs',0),FullJ_receiver_RHS=a.get('fullj_receiver_transpose_rhs',0),smallmatrix_Gram_RHS=c.get('Gram_RHS',0),smallmatrix_Gram_factorizations=c.get('Gram_factorizations',0),smallgram_gradients=c.get('smallgram_gradients',0),dual_accepted_iterations=c.get('dual_optimizer_accepted_iterations',0),fullj_builds=a.get('fullj_builds',0),all_unknowns_released=r.get('all_unknowns_released',False),final576_derivatives=s.get('final576_fullj_builds',0)+s.get('final576_smallgram_gradients',0),coarse_seconds=sum(z.get('wall_seconds',0) for z in stages if z.get('stage') in [8,12]),final576_seconds=sum(z.get('wall_seconds',0) for z in stages if z.get('stage')==24),archive_RMS1_distinct=s.get('archive_distinct_greedy_RMS1'),archive_pairwise_material_RMS_max=s.get('archive_pairwise_material_RMS_max')))
 group={}
 totals=['wall_seconds','cpu_seconds','full_factorizations','physical_complex_RHS','source_RHS','source_adjoint_RHS','FullJ_receiver_RHS','smallmatrix_Gram_RHS','smallmatrix_Gram_factorizations','smallgram_gradients','dual_accepted_iterations','fullj_builds','coarse_seconds','final576_seconds']
 for method in METHODS:
  mr=[r for r in rows if r['method']==method];valid=[r for r in mr if r['RMSE'] is not None]
  group[method]=dict(count=len(mr),finite_endpoint_count=len(valid),all_unknowns_released_count=sum(r['all_unknowns_released'] for r in mr),positive_final576_derivative_count=sum(r['final576_derivatives']>0 for r in mr),status_counts=dict(Counter(r['status'] for r in mr)),totals={key:sum(r[key] for r in mr) for key in totals},quality={key:dict(mean=float(np.mean([r[key] for r in valid])),median=float(np.median([r[key] for r in valid])),min=float(min(r[key] for r in valid)),max=float(max(r[key] for r in valid))) if valid else None for key in ['RMSE','IoU','SSIM']},archive_RMS1_distinct_counts=dict(Counter(str(r['archive_RMS1_distinct']) for r in mr)))
 batch=json.loads((ROOT/'audit/EXPERIMENT_BATCH_HOLDOUT.json').read_text());validation=json.loads((ROOT/'audit/EXPERIMENT_VALIDATION_HOLDOUT.json').read_text());generation=json.loads((ROOT/'audit/EXPERIMENT_GENERATION_HOLDOUT.json').read_text())
 dump(ROOT/'results/holdout_delivery_aggregation.json',dict(split='holdout',source_files_unchanged=True,cases=9,methods=3,expected_predictions=27,actual_predictions=len(rows),by_method=group,rows=rows,batch_wall_seconds=batch['total_batch_wall_seconds'],generation=generation,validation=validation,accounting_notes=['Physical PDE RHS separate from 12x12 Gram-system RHS','Condition estimates perform separately timed opaque LAPACK work beyond explicit source RHS; internal estimator RHS not exposed','Archive K are optimizer iterates, not branches or independent samples','No accepted-image gate inferred from numerical completion']))
 fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
 for method,label in zip(METHODS,LABELS):
  mr=[r for r in rows if r['method']==method and r['RMSE'] is not None]
  for ax,key in zip(axes,['cpu_seconds','physical_complex_RHS']):
   ax.scatter([r[key] for r in mr],[r['RMSE'] for r in mr],label=label,alpha=.85,s=28)
   for r in mr:ax.annotate(str(r['digit']),(r[key],r['RMSE']),xytext=(3,3),textcoords='offset points',fontsize=7)
 for ax in axes:ax.set_ylabel('Top1 material RMSE');ax.grid(alpha=.2)
 axes[0].set_xlabel('Paid optimizer CPU seconds per image');axes[1].set_xlabel('Paid physical PDE complex RHS per image');axes[1].legend(fontsize=8)
 fig.suptitle('Nine frozen holdout targets: image error and paid work; labels are digits');fig.savefig(ROOT/'figures/holdout_cost_vs_error.png',dpi=170);plt.close(fig)
 fig,axes=plt.subplots(1,3,figsize=(10,3.5),layout='constrained')
 for ax,key,title in zip(axes,['cpu_seconds','physical_complex_RHS','smallmatrix_Gram_RHS'],['Optimizer CPU seconds','Physical PDE complex RHS','12x12 Gram system RHS']):
  ax.bar(np.arange(3),[group[m]['totals'][key] for m in METHODS],color=['C0','C1','C2']);ax.set_xticks(np.arange(3),['Native','Reciprocal','Dual + Full-J'],rotation=20,fontsize=8);ax.set_title(title,fontsize=10);ax.grid(axis='y',alpha=.2)
 fig.suptitle('Main holdout totals across nine images; setup/validation reported separately');fig.savefig(ROOT/'figures/holdout_method_costs.png',dpi=170);plt.close(fig)
 print(json.dumps(dict(predictions=len(rows),source_unchanged=True,by_method=group,batch_wall_seconds=batch['total_batch_wall_seconds'])))
if __name__=='__main__':main()
