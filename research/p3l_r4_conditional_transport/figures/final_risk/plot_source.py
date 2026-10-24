#!/usr/bin/env python3
"""Standalone figures from closed canonical actions and nowEXPOSED truth only."""
from pathlib import Path
import json,sys,os,time,csv,shutil,resource
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p3r4.evaluation import checkpoint_identity,write_new_json
S=ROOT/'research/p3l_r4_conditional_transport';F=S/'results/fresh';O=S/'figures/final_risk';AUD=F/'independent_risk_statistical_audit_v1'
LEVEL=.9833333333333333;SEED=20261025099
FAMILIES=[('physics_probit64','Probit\nanchor'),('frozen_r2','Frozen R2\nCNN'),('old512_cnn','Reoptimized CNN\n512 contexts'),('strong_spectral','CNN\n1024 contexts'),('fusion_spline','CNN + Spline\nfeatures'),('fusion_diag','CNN + Diag\nfeatures'),('fusion_permuted','CNN + permuted\nSpline features'),('fusion_field','CNN + field\nmoment features')]
COLORS=['#6f7881','#353a40','#7a9dc4','#356a9b','#d47a24','#368777','#ad7a9b','#7597a8']
GALLERY=[('physics_probit64','Probit anchor\nonline context'),('frozen_r2','Frozen R2 CNN\n5842 parameters'),('old512_cnn_s2026102501','Reoptimized CNN\n512 contexts'),('strong_spectral_s2026102501','CNN · 1024 contexts\nsame 5842 parameters'),('fusion_spline_s2026102501','CNN + Spline features\n6418 parameters'),('fusion_diag_s2026102501','CNN + Diag features\n6418 parameters')]
def read(p):return json.loads(Path(p).read_text())
def check(v):
 if checkpoint_identity(v['path'])!=v:raise ValueError('closed source changed')
def csvwrite(path,rows):
 with path.open('x',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def save(fig,name):
 fig.savefig(O/(name+'.png'),dpi=180,facecolor='white');fig.savefig(O/(name+'.pdf'),facecolor='white');plt.close(fig)
def interval(v,ids):return np.quantile(v[ids].mean(1),[(1-LEVEL)/2,1-(1-LEVEL)/2])
def main():
 wall,cpu=time.monotonic(),time.process_time();audit=read(AUD/'AUDIT.json');receipt=read(S/'IV/jobs/independent_heldout_risk_audit_v1/CLOSED.json')
 if receipt['status']!='COMPLETE' or audit['status']!='INDEPENDENT_READ_ONLY_STATISTICAL_AUDIT_COMPLETE' or audit['fallback_failure_count']!=0:raise ValueError('independent closed audit required before plot claims')
 if O.exists():raise FileExistsError('new plot output only')
 O.mkdir(parents=True);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight','pdf.fonttype':42})
 sources=[checkpoint_identity(AUD/'AUDIT.json'),checkpoint_identity(S/'IV/jobs/independent_heldout_risk_audit_v1/CLOSED.json')];data={};truths={};parentlists={};gates={};source_rows=[];barrows=[];phase_rows=[];galleryrows=[]
 for group,n in [('locked',512),('confirmation',256)]:
  base=F/group/'risk45';gate=read(base/'online/ONLINE_COMPLETE.json');parents=gate['parents'];parentlists[group]=parents;gates[group]=gate
  if len(parents)!=n:raise ValueError('fixed cohort roster')
  for identity in gate['sealed_artifact_identities']:check(identity)
  with np.load(AUD/(group+'_seed_mean_parent_losses.npz')) as z:data[group]={k:z[k].copy() for k in z.files if k!='parent'}
  sources.append(checkpoint_identity(AUD/(group+'_seed_mean_parent_losses.npz')));sources.append(checkpoint_identity(base/'online/ONLINE_COMPLETE.json'))
  truth=[]
  for parent in parents:
   path=S/'inputs/offline_only'/(parent+'.npz');sources.append(checkpoint_identity(path))
   with np.load(path,allow_pickle=False) as z:truth.append(z['chi'])
  truths[group]=np.array(truth)
  ids=np.random.default_rng(SEED).integers(0,n,size=(4000,n));data[group]['bootstrap_ids']=ids
  for family,label in FAMILIES:
   v=data[group][family+'__clipped__full_mse'];ci=interval(v,ids);row=dict(cohort=group,family=family,action='clipped',n_parent=n,seed_loss_aggregation='mean perparent losses; no ensemble image',mean_full_complex_mse=float(v.mean()),pointwise_mean_CI_low=float(ci[0]),pointwise_mean_CI_high=float(ci[1]),CI_level=LEVEL,bootstrap_draws=4000,bootstrap_seed=SEED);barrows.append(row)
  for parent_index,parent in enumerate(parents):
   r2=data[group]['frozen_r2__clipped__full_mse'][parent_index];spline=data[group]['fusion_spline__clipped__full_mse'][parent_index]
   source_rows.append(dict(cohort=group,parent=parent,R2_clipped_full_mse=float(r2),Spline_seed_loss_mean_clipped_full_mse=float(spline),delta_R2_minus_Spline=float(r2-spline)))
 # Predefined four parents, fixed single seed; never an averaged image action.
 base=F/'locked/risk45';gallery_actions={}
 for mid,label in GALLERY:
  path=base/'online'/(mid+'.npz');sources.append(checkpoint_identity(path))
  with np.load(path,allow_pickle=False) as z:gallery_actions[mid]=z['prediction_clipped'][:4].copy()
 for i in range(4):
  parent=parentlists['locked'][i];truth=truths['locked'][i];images=[truth]+[gallery_actions[mid][i] for mid,_ in GALLERY];labels=['Ground truth\nnow EXPOSED']+[label for _,label in GALLERY]
  fig=plt.figure(figsize=(17.6,5.8));gs=fig.add_gridspec(2,8,width_ratios=[1]*7+[.055],left=.025,right=.985,bottom=.14,top=.83,hspace=.25,wspace=.13)
  for row,(component,lo,hi,tag) in enumerate([('real',.1,1.9,r'Re $\chi$'),('imag',.04,.22,r'Im $\chi$')]):
   for j,(image,label) in enumerate(zip(images,labels)):
    ax=fig.add_subplot(gs[row,j]);im=ax.imshow(getattr(image,component).reshape(24,24),origin='lower',interpolation='nearest',cmap='viridis',vmin=lo,vmax=hi);ax.set_xticks([]);ax.set_yticks([])
    if row==0:ax.set_title(label,fontsize=10,pad=7)
    if j==0:ax.set_ylabel(tag,fontsize=11)
    if row==1 and j>0:
     mse=float(np.mean(abs(image-truth)**2));ax.text(.5,-.10,f'full complex MSE {mse:.5f}',ha='center',va='top',transform=ax.transAxes,fontsize=8)
   cb=fig.colorbar(im,cax=fig.add_subplot(gs[row,7]));cb.set_ticks(np.linspace(lo,hi,5));cb.ax.tick_params(labelsize=8)
  fig.suptitle(f'{parent}: fixed first-four illustration | host susceptibility = 0.10 + 0.04i',fontsize=14,y=.96)
  fig.text(.5,.045,r'Clipped material risk actions $\widehat{\chi}$, not posterior means. Spline/Diag/CNN visualizations use fixed seed 2026102501; R2 uses seed 33.',ha='center',fontsize=10)
  fig.text(.5,.012,'512 → 1024 context expansion retrains the same CNN architecture; it is not a new architecture. Color scales fixed identically for all parents/methods.',ha='center',fontsize=9,color='#444444');save(fig,'gallery_'+parent)
  for j,(image,label) in enumerate(zip(images,labels)):galleryrows.append(dict(parent=parent,column=j,label=label.replace('\n',' '),action='truth' if j==0 else 'clipped fixedseed action',full_complex_mse=float(np.mean(abs(image-truth)**2)),Re_scale_min=.1,Re_scale_max=1.9,Im_scale_min=.04,Im_scale_max=.22,selection='first4 predefined L1, no performance selection'))
 # Both cohort bars, mean loss perparent over seeds, no image ensemble.
 fig,axes=plt.subplots(1,2,figsize=(17,5.7),sharey=True)
 for ax,group,title in zip(axes,['locked','confirmation'],['L1 · 512 parents','Independent confirmation · 256 parents']):
  rows=[r for r in barrows if r['cohort']==group];mean=np.array([r['mean_full_complex_mse'] for r in rows])*1000;lo=np.array([r['pointwise_mean_CI_low'] for r in rows])*1000;hi=np.array([r['pointwise_mean_CI_high'] for r in rows])*1000
  ax.bar(np.arange(len(rows)),mean,color=COLORS,yerr=np.stack((mean-lo,hi-mean)),capsize=3,error_kw=dict(lw=1));ax.set_xticks(np.arange(len(rows)),[label for _,label in FAMILIES],fontsize=8,rotation=25,ha='right');ax.set_title(title);ax.set_ylim(bottom=0);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
  for x,y in enumerate(mean):ax.text(x,y+.10,f'{y:.3f}',ha='center',va='bottom',fontsize=8)
 axes[0].set_ylabel(r'Full 576-cell complex MSE ($\times 10^{-3}$)');fig.suptitle('Frozen clipped risk actions · seed losses averaged within each parent',fontsize=14);fig.subplots_adjust(bottom=.29,top=.85,wspace=.08);fig.text(.5,.025,'Whiskers: pointwise 98.333% bootstrap intervals for parent-mean risk, 4000 draws. Extra controls do not add a familywise claim.',ha='center',fontsize=10);save(fig,'two_cohort_risk_bars')
 # Paired parent scatter/hist include every parent and all gain/loss tails.
 allr=np.concatenate([data[g]['frozen_r2__clipped__full_mse'] for g in ['locked','confirmation']]);alls=np.concatenate([data[g]['fusion_spline__clipped__full_mse'] for g in ['locked','confirmation']]);lims=[min(allr.min(),alls.min())*.8,max(allr.max(),alls.max())*1.2];deltas=(allr-alls)*1000;edges=np.linspace(deltas.min(),deltas.max(),41)
 fig,axes=plt.subplots(2,2,figsize=(12,9))
 for col,group,title in [(0,'locked','L1 · 512'),(1,'confirmation','Confirmation · 256')]:
  r2=data[group]['frozen_r2__clipped__full_mse'];spl=data[group]['fusion_spline__clipped__full_mse'];ax=axes[0,col];ax.scatter(r2,spl,s=15,alpha=.45,color=COLORS[4],edgecolors='none');ax.plot(lims,lims,'--',lw=1,color='#555555');ax.set(xscale='log',yscale='log',xlim=lims,ylim=lims,xlabel='Frozen R2 parent loss',ylabel='Spline family mean seed loss',title=title);ax.set_aspect('equal','box');ax.grid(alpha=.13)
  comp=next(c for c in audit['comparisons'] if c['panel']==group and c['reference_family']=='frozen_r2' and c['action']=='clipped');ax.text(.04,.96,f"Mean risk improvement {100*comp['improvement_relative']:.2f}%\n98.333% CI [{100*comp['relative_CI'][0]:.2f}, {100*comp['relative_CI'][1]:.2f}]%",transform=ax.transAxes,va='top',fontsize=10)
  ax=axes[1,col];diff=(r2-spl)*1000;ax.hist(diff,bins=edges,color=COLORS[4],alpha=.8);ax.axvline(0,color='#444444',lw=1);ax.axvline(diff.mean(),color=COLORS[1],lw=1.5,ls='--');ax.set(xlabel=r'$\ell_{R2}-\overline{\ell}_{Spline\ seeds}$ ($\times10^{-3}$); positive = gain',ylabel='Parents',title=f'All paired parent differences · {title}');ax.grid(axis='y',alpha=.12)
 fig.suptitle('Paired saved-action risk comparisons · no seed-ensemble image',fontsize=14);fig.tight_layout(rect=[0,.03,1,.95]);fig.text(.5,.012,'All 768 parents shown across separate frozen cohorts; no outlier trimming and no performance-selected gallery. Pooled views are descriptive only.',ha='center',fontsize=9);save(fig,'paired_R2_to_Spline')
 # Anchor-defined transition vs ALLothercells; full576 denominator ensures
 # contributions add exactly to total gain, even when transition mask empty.
 for group in ['locked','confirmation']:
  base=F/group/'risk45';g=gates[group];cfg=g['geometry']['config'];roi=np.array(g['geometry']['roi'],int);A=np.array(g['geometry']['A']);n=len(parentlists[group]);gt=truths[group]
  path=base/'online/online_anchor.npz';sources.append(checkpoint_identity(path))
  with np.load(path,allow_pickle=False) as z:q0=z['q0']
  anchor_s=(cfg['field_mean']+q0[:,:64]@A[roi].T)/cfg['phase_transition_halfwidth'];mask=np.zeros((n,576),bool);mask[:,roi]=(anchor_s>-1)&(anchor_s<1)
  def family_cell_loss(family):
   methods=[m for m in g['freeze']['methods'] if m['family']==family];values=[]
   for m in methods:
    p=base/'online'/(m['id']+'.npz');sources.append(checkpoint_identity(p))
    with np.load(p,allow_pickle=False) as z:pred=z['prediction_clipped']
    values.append(abs(pred-gt)**2)
   return np.mean(values,axis=0)
  spl=family_cell_loss('fusion_spline')
  for ref in ['frozen_r2','strong_spectral']:
   gain=family_cell_loss(ref)-spl;transition=(gain*mask).sum(1)/576;other=(gain*(~mask)).sum(1)/576;total=gain.mean(1);expected=data[group][ref+'__clipped__full_mse']-data[group]['fusion_spline__clipped__full_mse']
   if not np.allclose(transition+other,total,atol=1e-15,rtol=0) or not np.allclose(total,expected,atol=1e-15,rtol=0):raise ValueError('phase decomposition not equal to full-grid seed-loss gain')
   for i,parent in enumerate(parentlists[group]):phase_rows.append(dict(cohort=group,parent=parent,reference=ref,candidate='fusion_spline',transition_cells=int(mask[i].sum()),other_cells=int((~mask[i]).sum()),transition_gain_full576=float(transition[i]),other_gain_full576=float(other[i]),total_gain_full576=float(total[i]),mask_source='online anchor s; -1<s<1; allother576cells complement',seed_aggregation='percell loss mean, no averaged image'))
 phase_summary=[];fig,axes=plt.subplots(1,2,figsize=(12,5),sharey=True)
 for ax,group,title in zip(axes,['locked','confirmation'],['L1 · 512','Confirmation · 256']):
  for xpos,ref,label in [(0,'frozen_r2','Frozen R2 → Spline'),(1,'strong_spectral','Strong CNN → Spline')]:
   rows=[r for r in phase_rows if r['cohort']==group and r['reference']==ref];tr=np.array([r['transition_gain_full576'] for r in rows]);ot=np.array([r['other_gain_full576'] for r in rows]);total=tr+ot;ci=interval(total,data[group]['bootstrap_ids']);a,b=tr.mean()*1000,ot.mean()*1000
   ax.bar(xpos-.15,a,width=.3,color='#cc7728',label='Anchor-transition cells' if xpos==0 else None);ax.bar(xpos+.15,b,width=.3,color='#658d9e',label='All other cells' if xpos==0 else None);ax.plot([xpos-.35,xpos+.35],[total.mean()*1000]*2,color='#333333',lw=1.5);ax.text(xpos,max(a,b,total.mean()*1000)+.025,f'total {total.mean()*1000:.3f}',ha='center',fontsize=9)
   phase_summary.append(dict(cohort=group,reference=ref,candidate='fusion_spline',mean_transition_gain_full576=float(tr.mean()),mean_other_gain_full576=float(ot.mean()),mean_total_gain_full576=float(total.mean()),total_gain_CI_low=float(ci[0]),total_gain_CI_high=float(ci[1]),mean_transition_cell_count=float(np.mean([r['transition_cells'] for r in rows])),identity_max_abs_error=float(np.max(abs(total-(data[group][ref+'__clipped__full_mse']-data[group]['fusion_spline__clipped__full_mse']))))))
  ax.set_xticks([0,1],['Frozen R2 → Spline','Strong CNN → Spline']);ax.set_title(title);ax.axhline(0,color='#444444',lw=.8);ax.grid(axis='y',alpha=.12);ax.set_axisbelow(True)
 axes[0].set_ylabel(r'Mean MSE gain contribution ($\times10^{-3}$; denominator576)');axes[0].legend(frameon=False,fontsize=9);fig.suptitle('Anchor-phase gain decomposition · observed risk, not a mechanism proof',fontsize=13);fig.subplots_adjust(bottom=.23,top=.84,wspace=.10);fig.text(.5,.065,'Transition mask: online anchor −1 < s < 1; complement includes every other cell. Contributions sum to each parent’s full-grid risk gain.',ha='center',fontsize=9);fig.text(.5,.025,'Saved clipped actions; seed losses averaged at each cell. No truth-defined mask selection, no new inference, and no posterior claim.',ha='center',fontsize=9);save(fig,'anchor_transition_gain_decomposition')
 csvwrite(O/'plot_source_metrics.csv',source_rows);csvwrite(O/'risk_bar_source_metrics.csv',barrows);csvwrite(O/'predefined_gallery_source_metrics.csv',galleryrows);csvwrite(O/'phase_gain_per_parent.csv',phase_rows);csvwrite(O/'phase_gain_summary.csv',phase_summary)
 for v in sources:check(v)
 shutil.copyfile(__file__,O/'plot_source.py')
 write_new_json(O/'plot_source_metrics.json',dict(status='PLOTS_FROM_CLOSED_ARRAYS_ONLY',heldout_label_status='EXPOSED_AFTER_CLOSED_SCORING',failures=0,source_identities=sources,predefined_gallery_parents=parentlists['locked'][:4],gallery_fixed_seed=2026102501,R2_seed=33,action='clipped',seed_statistics='perparent/percell seed loss means; galleries fixed single seed; never average image action',host_chi=[.1,.04],fixed_color_scales=dict(Re=[.1,1.9],Im=[.04,.22]),risk_labels='material risk action chi-hat; no posterior mean/sample/density qualification',architecture='512to1024 retraining same5842parameterCNN; fusion6418 usesadditional4providerchannels',phase_summary=phase_summary,bootstrap=dict(draws=4000,seed=SEED,CI_level=LEVEL,extra_plot_intervals='pointwise descriptive; no newfamilywiseclaim'),model_inference=0,new_physics_RHS=0,new_fits=0,posthoc_gallery=False,wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,peak_RSS_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024)))
 print(json.dumps(dict(status='PLOTS_FROM_CLOSED_ARRAYS_ONLY',path=str(O),png_files=len(list(O.glob('*.png'))),wall_seconds=time.monotonic()-wall)))
if __name__=='__main__':main()
