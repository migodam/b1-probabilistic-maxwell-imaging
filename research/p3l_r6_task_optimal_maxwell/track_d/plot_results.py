"""Post-score Track D figures. Preparation only until root explicitly authorizes.

No physics, model loader, trainer or evaluator is imported. D truths are opened
only for fixed/display cases after GT_UNLOCK + SCORING_COMPLETE + root CLI flag.
"""
from __future__ import annotations
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[_key]='1'
from pathlib import Path
import argparse,csv,json,sys,time,shutil
import numpy as np
TRACK=Path(__file__).resolve().parent;OUT=TRACK.parent;sys.path.insert(0,str(OUT))
from track_d.common import freeze,identity,write,GROUPS,SEEDS
LABELS={'full_deterministic':'Full Det','selected_deterministic':'Low KKT-q','low_actual_deterministic':'Low actual','full_kkt_deterministic':'Full KKT-q','full_distill':'Full Distill','full_strong':'Full Strong'}
COLORS={'full_deterministic':'#334155','selected_deterministic':'#007f73','low_actual_deterministic':'#cf7b19','full_kkt_deterministic':'#6871bd','full_distill':'#c0517d','full_strong':'#8a949e'}
CONTRAST_LABELS={'risk_NI_1pct_same_frozen_weights':'Risk NI (1%)','KKT_vs_same_low_actual_gradient_risk':'KKT vs same-low actual','zero_shot_risk_vs_full_Distill_SECONDARY':'vs Full Distill (secondary)','whole_system_wall_20pct':'Wall saving (20%)'}
PHASES=('host','transition','pure')

def csvwrite(path,rows):
 if not rows:raise ValueError('nonempty source rows required')
 with Path(path).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def load_completed(root_authorized):
 if not root_authorized:raise RuntimeError('NOT_RUN: await explicit root authorization, then pass --root-authorized')
 for p in (TRACK/'GT_UNLOCK.json',TRACK/'scores/SCORING_COMPLETE.json',TRACK/'ALL_512_PREDICTIONS_SEALED.json'):
  if not p.exists():raise RuntimeError('NOT_RUN: root GT unlock and complete main score required before plotting')
 d=freeze();complete=json.loads((TRACK/'scores/SCORING_COMPLETE.json').read_text())
 if complete['status']!='COMPLETE' or complete['global_family_size']!=8:raise ValueError('complete frozen global8 score required')
 globalseal=TRACK/'ALL_512_PREDICTIONS_SEALED.json';unlock=json.loads((TRACK/'GT_UNLOCK.json').read_text());globalrecord=json.loads(globalseal.read_text())
 if unlock.get('status')!='ROOT_UNLOCKED_AFTER_SEAL' or unlock.get('global_seal_identity')!=identity(globalseal):raise ValueError('root unlock must bind exact global seal')
 if globalrecord['freeze_identity']!=identity(TRACK/'D_FREEZE.json'):raise ValueError('freeze changed since global seal')
 all_data={};sources=[]
 for group in GROUPS:
  path=TRACK/'scores'/f'{group}_RESULTS.json';res=json.loads(path.read_text())
  if res['statistics']!=d['statistics'] or res['n']!=256:raise ValueError('result statistics/parent count')
  lossfile=TRACK/'scores'/f'{group}_PARENT_LOSSES.npz'
  with np.load(lossfile,allow_pickle=False) as z:
   ids=z['parent'].astype(str).tolist();loss={m:z[m].copy() for m in d['estimator_arms']};counts=z['phase_counts'].copy()
  sealfile=TRACK/'predictions'/group/'PREDICTIONS_SEALED.json';seal=json.loads(sealfile.read_text())
  if ids!=seal['parents'] or len(ids)!=256 or len(set(ids))!=256:raise ValueError('parent order/independence')
  records={};pred={};seconds={}
  for rec in seal['records']:
   artifact=rec['identity']
   if identity(artifact['path'])!=artifact:raise ValueError('sealed prediction changed')
   with np.load(artifact['path'],allow_pickle=False) as z:
    if z['parent'].astype(str).tolist()!=ids:raise ValueError('prediction order')
    pred.setdefault(rec['method'],{})[rec['seed']]=z['prediction'].copy()
   records[(rec['method'],rec['seed'])]=rec
  for method in d['estimator_arms']:
   if set(pred[method])!=set(SEEDS) or loss[method].shape!=(256,3):raise ValueError('three frozen seed contract')
   p=np.stack([pred[method][s] for s in SEEDS],1)
   if not np.isfinite(p).all():raise ValueError('nonfinite sealed predictions')
   pred[method]=p
   seconds[method]=np.mean([np.array(records[(method,s)]['per_parent_cold_anchor_seconds'])+records[(method,s)]['per_parent_warm_learning_seconds'] for s in SEEDS],axis=0)
  contexts={};reports={};context_identities=[]
  for item in seal['context_artifacts']:
   for key in ('context_identity','report_identity'):
    if identity(item[key]['path'])!=item[key]:raise ValueError('sealed context/report changed')
   reports[(item['parent'],item['context'])]=json.loads(Path(item['report_identity']['path']).read_text())
   contexts[(item['parent'],item['context'])]=Path(item['context_identity']['path']);context_identities.append(item)
  all_data[group]=dict(result=res,ids=ids,loss=loss,counts=counts,pred=pred,seconds=seconds,contexts=contexts,reports=reports)
  sources.extend([identity(path),identity(lossfile),identity(sealfile),identity(TRACK/'scores'/f'{group}_PER_PARENT_SEED_METRICS.csv'),identity(TRACK/'scores'/f'{group}_BOOTSTRAP_DRAWS.npz')])
 return d,all_data,sources

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root-authorized',action='store_true');parser.add_argument('--output',type=Path,default=TRACK/'figures');args=parser.parse_args()
 d,data,sources=load_completed(args.root_authorized)
 destination=args.output.resolve()
 if not destination.is_relative_to(TRACK.resolve()):raise ValueError('figures must stay under track_d')
 if destination.exists():raise FileExistsError('figure output exists; no silent overwrite/re-selection')
 destination.mkdir(parents=True);source_dir=destination/'sources';source_dir.mkdir()
 for group in GROUPS:
  for suffix in ('PARENT_LOSSES.npz','PER_PARENT_SEED_METRICS.csv','RESULTS.json'):
   shutil.copy2(TRACK/'scores'/f'{group}_{suffix}',source_dir/f'{group}_{suffix}')
 # Plot imports only after all authorization/completion guards have passed.
 os.environ.setdefault('MPLCONFIGDIR',str(destination/'mpl_cache'))
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white','pdf.fonttype':42})
 methods=list(d['estimator_arms']);full=d['primary_full'];selected=d['primary_selected'];actual=d['actualgradient_reference'];distill=d['distill_reference'];manifest=[];started=time.time()
 def finish(fig,name,caption,source_names,selection='PREDECLARED_OR_ALL_PARENTS'):
  fig.savefig(destination/f'{name}.png',dpi=180,bbox_inches='tight');fig.savefig(destination/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
  manifest.append(dict(name=name,caption=caption,selection=selection,png=identity(destination/f'{name}.png'),pdf=identity(destination/f'{name}.pdf'),sources=[identity(source_dir/x) for x in source_names],QA='RENDERED_NOT_VISUALLY_REVIEWED'))
 # 1. Risk/cost, both conditions shown separately; all six arms included.
 pareto=[];fig,axes=plt.subplots(1,2,figsize=(12,4.5))
 for ax,(group,a) in zip(axes,data.items()):
  for j,method in enumerate(methods):
   risk=float(a['loss'][method].mean());wall=float(a['seconds'][method].mean());ax.scatter(wall,risk,color=COLORS[method],s=52,marker='*' if method==selected else 'o',zorder=3,label=LABELS[method])
   pareto.append(dict(group=group,method=method,mean_seed_loss_MSE=risk,mean_anchor_plus_learning_seconds=wall))
  pts=[r for r in pareto if r['group']==group];front=sorted([r for r in pts if not any(t['mean_seed_loss_MSE']<=r['mean_seed_loss_MSE'] and t['mean_anchor_plus_learning_seconds']<=r['mean_anchor_plus_learning_seconds'] and (t['mean_seed_loss_MSE']<r['mean_seed_loss_MSE'] or t['mean_anchor_plus_learning_seconds']<r['mean_anchor_plus_learning_seconds']) for t in pts)],key=lambda r:r['mean_anchor_plus_learning_seconds'])
  ax.plot([r['mean_anchor_plus_learning_seconds'] for r in front],[r['mean_seed_loss_MSE'] for r in front],color='#8c979f',ls='--',lw=.8);ax.set(title=f'{group}: noise x{GROUPS[group]:g}',xlabel='Actual anchor + warm learning time (s / parent)',ylabel='Primary raw MSE: mean of three seed losses');ax.grid(alpha=.15)
 fig.suptitle('Original execution-window risk / cost — frozen candidate, descriptive Pareto front');fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=6,frameon=False);fig.tight_layout(rect=[0,.10,1,.93]);csvwrite(source_dir/'RISK_COST.csv',pareto);finish(fig,'01_RISK_COST_PARETO','Each condition separate; no noise pooling, no new candidate selection. One measured batch64 inference per seed.',['RISK_COST.csv'])
 # 2. All eight frozen contrasts; estimates and one-sided upper bounds only.
 contrast_rows=[];fig,axes=plt.subplots(2,2,figsize=(12,7),gridspec_kw={'height_ratios':[3,1]})
 risk_names=['risk_NI_1pct_same_frozen_weights','KKT_vs_same_low_actual_gradient_risk','zero_shot_risk_vs_full_Distill_SECONDARY']
 for col,(group,a) in enumerate(data.items()):
  for row,names in [(0,risk_names),(1,['whole_system_wall_20pct'])]:
   ax=axes[row,col]
   for i,name in enumerate(names):
    b=a['result']['comparisons'][name];v,u=b['estimate'],b['one_sided_upper'];color='#007f73' if b['gate'] else '#ad5255'
    ax.plot([v,u],[i,i],color=color,lw=2);ax.scatter([v],[i],color=color,s=34);ax.scatter([u],[i],facecolors='none',edgecolors=color,marker='s',s=38)
    contrast_rows.append(dict(group=group,contrast=name,estimate=v,one_sided_upper=u,confidence=b['one_sided_confidence'],gate=b['gate']))
   ax.axvline(0,color='#555',lw=.9);ax.set_yticks(range(len(names)),[CONTRAST_LABELS[x] for x in names]);ax.invert_yaxis();ax.set_xlabel('MSE difference (negative favors selected)' if row==0 else 'Seconds difference: selected - .8 full');ax.grid(axis='x',alpha=.15)
   if row==0:ax.set_title(f'{group}: noise x{GROUPS[group]:g}')
 fig.suptitle('Global eight-comparison family — one-sided 99.375% upper bounds\nCircle: estimate; open square: upper bound; gate requires upper < 0');fig.tight_layout();csvwrite(source_dir/'GLOBAL8_CONTRASTS.csv',contrast_rows);finish(fig,'02_GLOBAL8_CONTRASTS','All four contrasts per condition. Secondary Distill contrast cannot select a new champion.',['GLOBAL8_CONTRASTS.csv'])
 # Shared map gallery helper. Truth is accessed ONLY here after score/unlock guards.
 cfg=json.loads((OUT/'inputs/Research_P3LR5T/decoder_sources/ONLINE_TARGET.json').read_text());host=complex(*cfg['host_chi']);limits=[(host.real,host.real+cfg['inclusion_increment_real_interval'][1]),(host.imag,host.imag+cfg['inclusion_increment_imag_interval'][1])]
 gallery_methods=['GT','Full physics',full,selected,actual,distill]
 def gallery(group,indices,name,posthoc=False):
  a=data[group];ids=[a['ids'][i] for i in indices];maps=[]
  for i,parent in zip(indices,ids):
   with np.load(OUT/'raw/offline_truth'/f'{parent}.npz',allow_pickle=False) as z:
    if str(z['parent'])!=parent:raise ValueError('display truth parent order')
    truth=z['chi'].copy()
   with np.load(a['contexts'][(parent,'full')],allow_pickle=False) as z:physical=z['chi0'].copy()
   maps.append(np.stack([truth,physical,*[a['pred'][method][i].mean(0) for method in (full,selected,actual,distill)]]))
  maps=np.array(maps);fig,axes=plt.subplots(2*len(indices),6,figsize=(15,4.5*len(indices)),squeeze=False)
  for ri,parent in enumerate(ids):
   for component in (0,1):
    for col,method in enumerate(gallery_methods):
     ax=axes[2*ri+component,col];image=(maps[ri,col].real if component==0 else maps[ri,col].imag).reshape(24,24);im=ax.imshow(image,cmap='viridis',vmin=limits[component][0],vmax=limits[component][1],origin='lower');ax.set_xticks([]);ax.set_yticks([])
     if ri==0 and component==0:ax.set_title(method if method in ('GT','Full physics') else LABELS[method],fontsize=10)
     if col==0:ax.set_ylabel(f'{parent}\n'+('Real' if component==0 else 'Imag'),fontsize=8)
     if col==5:fig.colorbar(im,ax=ax,fraction=.046,pad=.02)
  qualifier='OFFLINE_ERROR_SELECTED post-hoc cases' if posthoc else 'Fixed parent indices 0000 / 0001'
  fig.suptitle(f'{group}: {qualifier}\nRaw three-seed prediction ensembles for display; primary risk averages SEED LOSSES. Fixed legal-truth color scales; raw overshoots may saturate.');fig.tight_layout(rect=[0,0,1,.94])
  source=f'{name}_DISPLAY_ARRAYS.npz';np.savez_compressed(source_dir/source,parent=np.array(ids),method=np.array(gallery_methods),material_maps=maps,component_limits=np.array(limits));finish(fig,name,qualifier+'; display ensembles are not the primary risk estimator.',[source],'OFFLINE_ERROR_SELECTED' if posthoc else 'FIXED_PARENT_0000_0001')
 for group in data:gallery(group,[0,1],f'03_FIXED_{group}')
 # 4. Phase risks use full-anchor classification, never GT-defined phase.
 phase_rows=[];fig,axes=plt.subplots(1,2,figsize=(13,4.6))
 for ax,(group,a) in zip(axes,data.items()):
  x=np.arange(3);width=.12
  for j,method in enumerate(methods):
   values=[]
   for phase in PHASES:
    r=a['result']['secondary_seed_mean_metrics'][method]['raw']['anchor_phase'][phase];values.append(np.nan if r['cell_weighted_mse'] is None else r['cell_weighted_mse']);phase_rows.append(dict(group=group,method=method,phase=phase,cell_count=r['cell_count'],cell_weighted_mse=r['cell_weighted_mse'],full_grid_mse_contribution=r['full_grid_mse_contribution']))
   ax.bar(x+(j-(len(methods)-1)/2)*width,values,width,color=COLORS[method],label=LABELS[method])
  counts=a['result']['phase_ROI_counts_host_transition_pure'];ax.set_xticks(x,[f'{p}\n{n:,} cells' for p,n in zip(PHASES,counts)]);ax.set(title=group,ylabel='Raw MSE per phase cell');ax.grid(axis='y',alpha=.15)
 fig.legend(*axes[1].get_legend_handles_labels(),loc='lower center',ncol=6,frameon=False);fig.suptitle('Risk by FULL-anchor phase — seed-mean losses, fixed ROI, no truth phase labels');fig.tight_layout(rect=[0,.10,1,.93]);csvwrite(source_dir/'ANCHOR_PHASE_RISK.csv',phase_rows);finish(fig,'04_ANCHOR_PHASE_RISK','Phase count-weighted risk; full-grid contributions retained in source table.',['ANCHOR_PHASE_RISK.csv'])
 # 5. Residual, risk and cost scatter; finite-best nonconvergence is retained.
 scatter=[];fig,axes=plt.subplots(2,3,figsize=(14,8))
 for row,(group,a) in enumerate(data.items()):
  risk=a['loss'][selected].mean(1);delta=risk-a['loss'][full].mean(1);phi={k:[] for k in ('full','low')};low_success=[]
  for i,parent in enumerate(a['ids']):
   for kind in ('full','low'):
    with np.load(a['contexts'][(parent,kind)],allow_pickle=False) as z:phi[kind].append(float(z['phi']))
   low_success.append(bool(a['reports'][(parent,'low')]['original_probit_report']['optimizer'].get('success',False)))
   scatter.append(dict(group=group,parent=parent,full_phi=phi['full'][-1],low_phi=phi['low'][-1],selected_seed_mean_MSE=float(risk[i]),selected_minus_full_MSE=float(delta[i]),selected_wall_seconds=float(a['seconds'][selected][i]),low_scipy_success=low_success[-1]))
  color=np.where(low_success,'#007f73','#ba5956')
  for col,kind in enumerate(('full','low')):axes[row,col].scatter(phi[kind],risk,c=color,s=14,alpha=.55);axes[row,col].set(xlabel=f'{kind.title()} anchor Phi',ylabel='Selected raw seed-mean MSE',title=group)
  axes[row,2].scatter(a['seconds'][selected],delta,c=color,s=14,alpha=.55);axes[row,2].axhline(0,lw=.8,color='gray');axes[row,2].set(xlabel='Selected anchor + learning time (s)',ylabel='Selected - Full raw MSE',title=group)
 fig.suptitle('Per-parent diagnostics — green: low Probit SciPy success; red: finite-best retained\nResidual is not material risk or posterior qualification; no points excluded.');fig.tight_layout(rect=[0,0,1,.94]);csvwrite(source_dir/'RESIDUAL_RISK_COST.csv',scatter);finish(fig,'05_RESIDUAL_RISK_COST','All independent parents; descriptive associations only.',['RESIDUAL_RISK_COST.csv'])
 # 6. Finite seed disagreement decomposition, not population training variance.
 disagreement=[];fig,axes=plt.subplots(1,2,figsize=(13,4.8))
 harmcases=[]
 for ax,(group,a) in zip(axes,data.items()):
  ensemble_means=[];dis_means=[]
  for method in methods:
   p=a['pred'][method];dis=np.mean(abs(p-p.mean(1,keepdims=True))**2,axis=(1,2));seed_loss=a['loss'][method].mean(1);ensloss=seed_loss-dis
   if np.min(ensloss)<-1e-10:raise ValueError('MSE disagreement identity violated')
   ensemble_means.append(float(ensloss.mean()));dis_means.append(float(dis.mean()))
   for i,parent in enumerate(a['ids']):disagreement.append(dict(group=group,parent=parent,method=method,seed_mean_MSE=float(seed_loss[i]),ensemble_MSE_from_identity=float(ensloss[i]),finite_three_seed_disagreement=float(dis[i])))
  xx=np.arange(len(methods));ax.bar(xx,ensemble_means,color='#b8c6d5',label='Ensemble MSE');ax.bar(xx,dis_means,bottom=ensemble_means,color='#526d8d',label='Fixed-seed disagreement');ax.set_xticks(xx,[LABELS[m] for m in methods],rotation=30,ha='right');ax.set(title=group,ylabel='MSE: seed-mean = ensemble + disagreement')
  delta=a['loss'][selected].mean(1)-a['loss'][full].mean(1);worst=np.argsort(-delta,kind='stable')[:2].tolist()
  for i in worst:harmcases.append(dict(group=group,parent=a['ids'][i],parent_index=i,selection='OFFLINE_ERROR_SELECTED largest selected-minus-full raw seed-mean MSE; may be negative',selected_minus_full_MSE=float(delta[i])))
  gallery(group,worst,f'07_POSTHOC_WORST_DELTA_{group}',posthoc=True)
 fig.suptitle('Three frozen checkpoints: exact finite-ensemble decomposition\nNot population training variance; ensemble is a separate descriptive action.');fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False);fig.tight_layout(rect=[0,.10,1,.93]);csvwrite(source_dir/'SEED_DISAGREEMENT.csv',disagreement);csvwrite(source_dir/'POSTHOC_CASE_SELECTION.csv',harmcases);finish(fig,'06_SEED_DISAGREEMENT','All parent/method fixed-seed disagreements; no seed sample is an independent parent.',['SEED_DISAGREEMENT.csv','POSTHOC_CASE_SELECTION.csv'])
 shutil.copy2(Path(__file__),source_dir/'plot_results.py')
 write(destination/'FIGURE_MANIFEST.json',dict(status='RENDERED_REQUIRES_VIEW_IMAGE_QA',created_utc=time.time(),root_explicit_authorization_flag=True,main_score_complete=True,GT_unlock_identity=identity(TRACK/'GT_UNLOCK.json'),global_prediction_seal_identity=identity(TRACK/'ALL_512_PREDICTIONS_SEALED.json'),source_script=identity(source_dir/'plot_results.py'),primary_score_sources=sources,copied_primary_sources=[identity(p) for p in source_dir.glob('*') if p.name.endswith(('PARENT_LOSSES.npz','PER_PARENT_SEED_METRICS.csv','RESULTS.json'))],figures=manifest,new_physics_RHS=0,new_training=0,no_model_tuning=True,truth_access_scope='fixed0000/0001 and two largest selected-minus-full risk-delta display cases per condition only',elapsed_seconds=time.time()-started,labels='English; Chinese captions/report supplied externally'))
 print(json.dumps(dict(status='RENDERED_REQUIRES_VIEW_IMAGE_QA',figures=len(manifest),output=str(destination),new_physics_RHS=0)),flush=True)
if __name__=='__main__':main()
