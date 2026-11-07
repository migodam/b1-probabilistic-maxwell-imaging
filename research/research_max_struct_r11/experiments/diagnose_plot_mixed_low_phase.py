"""All8 mixed low-phase visible display plus paired low-region RMSE/mean, zero solver."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,csv,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).absolute().parents[1];BASE=ROOT/'delegated/diagnostics/bothblocks_confirmation';FIG=ROOT/'FIGURES/bothblocks_confirmation';AUDIT=BASE/'finite_noise_prior_audit/CANDIDATE_DATA_PRIOR_ALL80.csv'
METHODS=('cie_single','cie_dual','global_dual','fullG19','fullGreenBothBlocks');LABELS=('CIE','Dual','Global','G19 material','BothBlocks');COLORS=plt.get_cmap('tab10').colors

def main():
 wall=time.perf_counter()
 if not (BASE/'SCORING_RECEIPT.json').exists():raise RuntimeError('Only already authorized and sealed16 results may be plotted')
 cfg=json.loads((ROOT/'configs/twod_bothblocks16_generator.json').read_text());cases=[c for c in cfg['cases'] if c['regime']=='mixed'];seal=json.loads((BASE/'ALL80_TERMINAL_PREDICTION_SEAL.json').read_text());src={(e['task']['case_id'],e['task']['method']):e['prediction_path'] for e in seal['entries']}
 with AUDIT.open() as f:allrows=list(csv.DictReader(f))
 rows=[r for r in allrows if r['regime']=='mixed'];ix={(r['case_id'],r['method']):r for r in rows}
 if len(cases)!=8 or len(rows)!=40 or any((c['case_id'],m) not in ix for c in cases for m in METHODS):raise RuntimeError('Must include all8 mixed x5, with no success/GT-score filtering')
 originals={p.name:dict(bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in FIG.glob('*') if p.name in ('MATERIALS_HIGH.png','MATERIALS_HIGH.pdf','MATERIALS_MIXED.png','MATERIALS_MIXED.pdf')}
 cmap=plt.get_cmap('viridis').copy();cmap.set_over('#d95f02')
 fig,axes=plt.subplots(8,6,figsize=(17,23),layout='constrained')
 for ri,c in enumerate(cases):
  cid=c['case_id']
  with np.load(ROOT/'sealed_truth/bothblocks16'/f'{cid}.npz') as z:gt=z['m'].copy()
  materials=[gt]
  for method in METHODS:
   with np.load(src[cid,method]) as z:materials.append(z['m'].copy())
  for col,m in enumerate(materials):
   ax=axes[ri,col];im=ax.imshow(m.reshape(24,24),origin='lower',vmin=0,vmax=3.5,cmap=cmap,interpolation='nearest');ax.set_xticks([]);ax.set_yticks([])
   if col==0:ax.set_title('GT | low chi 2.5',fontsize=11);ax.set_ylabel(cid[-4:]+' '+c['family'].replace('_',' ')+'\ntrue high chi '+str(c['high_chi']),fontsize=10)
   else:
    row=ix[cid,METHODS[col-1]];ax.set_title(f"{LABELS[col-1]}\nlow RMSE {float(row['low_region_RMSE']):.3f} | mean {float(row['low_region_mean_estimate']):.3f}",fontsize=9.5,fontweight='bold' if col==5 else 'normal')
 fig.colorbar(im,ax=axes,shrink=.6,extend='max',label='Low-phase display chi 0–3.5 only; orange = all values >3.5 saturated')
 fig.suptitle('All8 mixed objects | low-phase visibility scale 0–3.5\nHigh phase intentionally saturated in orange; this is NOT a full-contrast material display',fontsize=14)
 fig.supxlabel('All5 frozen methods and all8 mixed objects shown. Original quantitative 0–20 contrast plates remain unchanged. No solver/FW; no48GT.',fontsize=10)
 files=[]
 for ext in ('png','pdf'):
  p=FIG/f'MIXED_LOW_PHASE_SCALE_0_3p5.{ext}';fig.savefig(p,dpi=150);files.append(p.name)
 plt.close(fig)
 fig,axes=plt.subplots(1,2,figsize=(17,6.5),layout='constrained');x=np.arange(8);offsets=np.linspace(-.2,.2,5);labels=[c['case_id'][-4:]+'\n'+c['family'].replace('_',' ') for c in cases]
 for j,(key,title,ylabel) in enumerate([('low_region_RMSE','Spatial error on truth-low cells','Low-region RMSE (absolute chi units)'),('low_region_mean_estimate','Regional mean does not measure spatial fidelity','Mean estimated chi on truth-low cells')]):
  ax=axes[j]
  for ci,c in enumerate(cases):
   vals=[float(ix[c['case_id'],m][key]) for m in METHODS];ax.vlines(ci,min(vals),max(vals),color='#bbbbbb',lw=1,zorder=1)
  for mi,method in enumerate(METHODS):
   vals=[float(ix[c['case_id'],method][key]) for c in cases];ax.scatter(x+offsets[mi],vals,color=COLORS[mi],s=42,marker=('o','s','^','D','P')[mi],label=LABELS[mi],zorder=3)
  ax.set_xticks(x,labels,fontsize=9);ax.set_ylabel(ylabel);ax.set_title(title,fontsize=12);ax.grid(axis='y',alpha=.2)
  if j==1:ax.axhline(2.5,color='black',ls='--',lw=1,label='GT low chi = 2.5')
 handles,legendlabels=axes[1].get_legend_handles_labels();fig.legend(handles,legendlabels,loc='outside lower center',ncol=6,fontsize=10)
 fig.suptitle('All8 mixed objects, paired by object | cached finite-noise/prior audit fields\nEach five-marker group uses the same truth-low mask (GT chi=2.5); no gate or success-based selection',fontsize=13)
 for ext in ('png','pdf'):
  p=FIG/f'MIXED_LOW_REGION_RMSE_AND_MEAN.{ext}';fig.savefig(p,dpi=160);files.append(p.name)
 plt.close(fig)
 metrics=[]
 for mi,method in enumerate(METHODS):
  v=[ix[c['case_id'],method] for c in cases];metrics.append(dict(method=method,objects=8,median_low_region_RMSE=float(np.median([float(r['low_region_RMSE']) for r in v])),median_low_region_mean_estimate=float(np.median([float(r['low_region_mean_estimate']) for r in v])),truth_low_chi=2.5))
 (FIG/'MIXED_LOW_PHASE_METRICS.json').write_text(json.dumps(metrics,indent=2)+'\n')
 unchanged=all((FIG/name).stat().st_size==v['bytes'] and (FIG/name).stat().st_mtime_ns==v['mtime_ns'] for name,v in originals.items())
 receipt=dict(status='COMPLETE',purpose='Offline low-phase visibility and spatial-error inspection, not new contrast inference or algorithm',objects=8,methods=5,all_mixed_included=True,success_filter=False,display_range=[0,3.5],high_phase_display='All material estimates and GT >3.5 use one orange saturation color; true high12/18 not quantitatively shown',full_contrast_0_20_originals_unchanged=unchanged,original_full_contrast_files=originals,source_cached_region_table=str(AUDIT.relative_to(ROOT)),source_materials='Previously sealed80 predictions, only8 mixed GT reused',files=files,prior_or_threshold_changes=0,new_solver=0,new_forward=0,new_adjoint=0,GT48_reads=0,GT_scope='Already root-authorized16GT, only mixed8 for this plot',near_exit_child_CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall)
 (FIG/'MIXED_LOW_PHASE_PLOT_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
