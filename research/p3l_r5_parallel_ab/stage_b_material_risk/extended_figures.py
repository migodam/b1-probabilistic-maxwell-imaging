import sys,json,time
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];OLD=ROOT.parent/'b1-p3l-r4-conditional-transport';sys.path.insert(0,str(OLD/'src'))
import numpy as np,torch,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p3r4.risk import full_material,roi_components,RiskActionPredictor
from p3r2.residual_risk import read_anchor,online_features
R=json.loads((OUT/'RESULTS.json').read_text());d=dict(np.load(OUT/'prepared_validation.npz'));g=torch.load(OUT/'geometry.pt',map_location='cpu',weights_only=False)['geometry'];torch.set_num_threads(1);seeds=[2026102601,2026102602,2026102603];oldseeds=[2026102501,2026102502,2026102503]
labels=['GT (existing VAL)','Physics anchor','Strong CNN ensemble','K32 ensemble','Det. matched ensemble','Distill128 ensemble'];material=[d['truth'],full_material(d['anchor'],g),np.mean([full_material(d[f'strong_s{s}'],g) for s in oldseeds],0),np.mean([full_material(d[f'K32_s{s}'],g) for s in oldseeds],0),np.mean([np.load(OUT/'experiments'/f'deterministic_matched_s{s}'/'validation.npz')['prediction'] for s in seeds],0),np.mean([np.load(OUT/'learned_context'/f'distilled_w128_n32768_s{s}'/'validation.npz')['prediction'] for s in seeds],0)];protocol=dict(indices=[0,1,2,3],parents=d['parent'][:4].tolist(),selection='fixedVALindices0..3; no loss-basedselection; notnewlockedtruth',prediction='3seedensembleimages; separatefrommeanseedlossprimarymetric',raw_predictions=True);(OUT/'GALLERY_PROTOCOL.json').write_text(json.dumps(protocol,indent=2)+'\n')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'figure.dpi':150})
for comp,attr,vmin,vmax,cmap in [('Re','real',.1,1.9,'viridis'),('Im','imag',.04,.22,'magma')]:
 fig,axs=plt.subplots(4,6,figsize=(15,10.5),layout='constrained')
 for i in range(4):
  for j,(name,arr) in enumerate(zip(labels,material)):
   im=axs[i,j].imshow(getattr(arr[i],attr).reshape(24,24),cmap=cmap,vmin=vmin,vmax=vmax,origin='lower',interpolation='nearest');axs[i,j].set_xticks([]);axs[i,j].set_yticks([])
   if i==0:axs[i,j].set_title(name,fontsize=9)
   if j==0:axs[i,j].set_ylabel(str(d['parent'][i]),fontsize=9)
 fig.colorbar(im,ax=axs,shrink=.8,label=f'{comp}(chi), shared original material scale',extend='both');fig.suptitle(f'Fixed existing VAL indices0–3 · {comp} susceptibility · raw ensemble predictions, no new truth');fig.savefig(OUT/f'MATERIAL_GALLERY_{comp.upper()}.png',bbox_inches='tight');plt.close(fig)
# Regionper-cellrisk panel, using original y-onlyanchor phase.
fams=['strong_CNN','original_K32_frozen','deterministic_matched','distilled_w128_n32768'];names=['StrongCNN','K32','Det.matched','Distill128'];fig,axs=plt.subplots(1,3,figsize=(13,4.3),layout='constrained')
for ax,p in zip(axs,R['phase']):
 vals=[p['risks'][f] for f in fams];ax.bar(range(4),vals,color=['#828b95','#e8a75b','#3599a9','#185f82']);ax.set_xticks(range(4),names,rotation=22,ha='right');ax.set_title(f"{p['region']} · {p['ROI_cell_fraction']*100:.1f}% ROI");ax.set_ylabel('Complex MSE per selected cell');ax.grid(axis='y',alpha=.2)
fig.suptitle('Per-cell material risk in original anchor-defined phase regions');fig.savefig(OUT/'PHASE_PER_CELL_RISK.png',bbox_inches='tight');plt.close(fig)
# Frozen singular spectra andallrankriskcurves, ridge=1 frozen acrosscontrols.
arms=['raw_y','decoder_only','physical_residual','same_spectrum_random'];curve={};records=[]
for arm in arms:
 rows=[json.loads((OUT/'low_rank_correction'/f'nested_{arm}_s{s}'/'RESULT.json').read_text()) for s in seeds];spec=np.mean([r['selected']['train_singular_values'] for r in rows],0);val=[];oof=[]
 for k in [2,4,8,16]:
  rr=[next(c for c in r['candidates'] if c['rank']==k and c['ridge']==1.) for r in rows];val.append(float(np.mean([r['validation_full_MSE'] for r in rr])));oof.append(float(np.mean([r['nested_corrected_OOF_MSE'] for r in rr])))
 curve[arm]=dict(singular_values=spec.tolist(),ranks=[2,4,8,16],validation_MSE=val,nested_OOF_MSE=oof,nested_baseline_OOF_MSE=float(np.mean([r['selected']['nested_baseline_OOF_MSE'] for r in rows])))
# Cost-onlyshape-equivalent rankkernels: repeatstoredrank2factors tolargerrank;
# riskvaluescomeonlyfromactualsavedfitresults, neverthese syntheticweights.
S=OLD/'research/p3l_r4_conditional_transport';a=read_anchor(S/'inputs/strong_context_1024_256/anchor_validation.npz');anc={k:a[k][:64] for k in ['observation','residual','q0','chi0','J']};y=d['y'][:64];baseline=RiskActionPredictor(S/'models/risk/strongrisk_spectral_s2026102501/BEST.pt')
for arm in arms:
 fit=dict(np.load(OUT/'low_rank_correction'/f'nested_{arm}_s2026102601'/'BEST.npz'))
 for k in [2,4,8,16]:
  F=np.tile(fit['input_factor'],(1,(k+1)//2))[:,:k];G=np.tile(fit['output_factor'],((k+1)//2,1))[:k]
  def call():
   base=baseline.predict(y,anc,batch_size=64)
   if arm in ['raw_y','same_spectrum_random']:raw=y.astype(float)
   elif arm=='decoder_only':raw=np.concatenate((roi_components(anc['chi0'],g['roi']).reshape(64,-1),anc['q0']),1).astype(float)
   else:raw=np.concatenate((online_features(anc,g)['image'][:,2:6].reshape(64,-1),anc['q0']),1).astype(float)
   if arm=='same_spectrum_random':X=(((raw-fit['raw_mean'])/fit['raw_scale']@fit['random_projection'])-fit['projection_mean'])@fit['spectrum_factor']
   else:X=(raw-fit['raw_mean'])/fit['raw_scale']@fit['projection']
   corr=((X-fit['input_mean'])@F@G+fit['output_mean'])*24;return full_material(roi_components(base,g['roi'])+corr.reshape(64,2,-1),g)
  call();times=[]
  for rep in range(15):start=time.monotonic();call();times.append(time.monotonic()-start)
  records.append(dict(arm=arm,rank=k,batch=64,raw_seconds=times,median_seconds_per_parent=float(np.median(times)/64),scope='completeCPUbaseline+realfeaturetransform+shapeequivalentdense rankk correction+outputassembly; storedrank2factorsrepeatedforcostonly, trainedrankkriskcomesfromsavedactualfit',new_fits=0,physics_RHS=0))
(OUT/'RANK_COST_CURVES.json').write_text(json.dumps(dict(status='COMPLETE',records=records,new_fits=0,physics_RHS=0),indent=2)+'\n');(OUT/'RISK_SPECTRUM_CURVES.json').write_text(json.dumps(dict(status='COMPLETE',ridge=1.,curves=curve,scope='TRAINcrosscovsingularspectrum andactualfrozenrank2/4/8/16nestedOOF/VALrisk',physics_RHS=0),indent=2)+'\n')
fig,axs=plt.subplots(2,2,figsize=(12.5,9),layout='constrained');colors=['#488faa','#c2913c','#b84c53','#626ca5']
for arm,color in zip(arms,colors):
 c=curve[arm];ss=np.array(c['singular_values']);axs[0,0].semilogy(np.arange(1,33),ss[:32]**2,label=arm,color=color);axs[0,1].plot(c['ranks'],c['nested_OOF_MSE'],'o-',label=arm,color=color);axs[1,0].plot(c['ranks'],c['validation_MSE'],'o-',label=arm,color=color);times=[next(r['median_seconds_per_parent']*1e6 for r in records if r['arm']==arm and r['rank']==k) for k in c['ranks']];axs[1,1].plot(times,c['validation_MSE'],'o-',label=arm,color=color)
 for x,y,k in zip(times,c['validation_MSE'],c['ranks']):axs[1,1].annotate(str(k),(x,y),fontsize=7)
axs[0,0].set_title('Predictable residual spectrum (TRAIN crosscov)');axs[0,0].set_xlabel('Singular index');axs[0,0].set_ylabel('Squared singular value');axs[0,0].legend(fontsize=8)
axs[0,1].axhline(curve['raw_y']['nested_baseline_OOF_MSE'],color='gray',ls='--',label='Baseline outerOOF');axs[0,1].set_title('Honest nested OOF material risk');axs[0,1].set_xlabel('Correction rank, ridge=1');axs[0,1].set_ylabel('Complex full-grid MSE');axs[0,1].legend(fontsize=8)
axs[1,0].axhline(next(r['seed_mean_parent_MSE'] for r in R['summary'] if r['family']=='strong_CNN'),color='gray',ls='--');axs[1,0].set_title('Existing VAL256 risk, not locked confirmation');axs[1,0].set_xlabel('Correction rank, ridge=1');axs[1,0].set_ylabel('Complex full-grid MSE');axs[1,1].set_title('Risk–cost curve (rankkernel shape benchmark)');axs[1,1].set_xlabel('Complete CPU µs / parent, batch64');axs[1,1].set_ylabel('Existing VAL MSE')
for ax in axs.flat:ax.grid(alpha=.2)
fig.suptitle('Frozen RRR rank-grid: training spectrum does not predict held-out gain',fontsize=13);fig.savefig(OUT/'RRR_SPECTRUM_RANK_RISK_COST.png',bbox_inches='tight');fig.savefig(OUT/'RRR_SPECTRUM_RANK_RISK_COST.svg',bbox_inches='tight');plt.close(fig)
print('fixedVALgallery, phasepanel and RRRspectrum/rank/risk/cost figures written')
