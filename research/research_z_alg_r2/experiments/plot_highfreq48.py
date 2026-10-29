"""Private-data postfreeze plotting only; no solver calls or parameter tuning."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize
P=Path(__file__).resolve().parents[1]
def main():
 if not (P/'experiments/HIGHFREQ48_FREEZE_RECEIPT.json').exists():raise RuntimeError('Need prediction freeze before truth plot')
 c=json.loads((P/'experiments/HIGHFREQ48_CONTRACT.json').read_text());scores=json.loads((P/'results/highfreq48_scores.json').read_text());truth=np.load(P/'private/highfreq48/sealed_targets.npz');centers=np.array(c['geometry']['centers']);radius=c['geometry']['radius'];matnorm=Normalize(.5,4);errnorm=Normalize(-4,4);families=c['private_case_families'];arms=[('exactSVD','Raw exact SVD'),('FullJ_midpoint','Full-J midpoint'),('SVD_FullJ','SVD + Full-J')];fig,axs=plt.subplots(6,7,figsize=(21,14),layout='constrained')
 def draw(ax,values,error=False):
  coll=PatchCollection([Circle(x,radius) for x in centers],cmap='RdBu_r' if error else 'viridis',norm=errnorm if error else matnorm,edgecolors='.3',linewidths=.35);coll.set_array(values);ax.add_collection(coll);ax.set(xlim=(-.39,.39),ylim=(-.29,.29),aspect='equal');ax.set_xticks([-.3,0,.3]);ax.set_yticks([-.2,0,.2]);ax.tick_params(labelsize=6);ax.grid(alpha=.13)
 for i in range(3):
  case=f'HF48_{i:02d}';gt=truth[case]
  for mi,mode in enumerate(['clean','relative_complex_1e-3']):
   row=2*i+mi;draw(axs[row,0],gt);axs[row,0].set_title(f'Truth: {families[i]}\n'+('Clean' if mode=='clean' else 'Same image, noise 0.001'),fontsize=9);axs[row,0].set_ylabel('y (model length)',fontsize=6)
   for j,(arm,label) in enumerate(arms):
    r=json.loads((P/f'results/highfreq48_{case}_{mode}_{arm}.json').read_text());m=np.array(r['m']);score=next(x for x in scores['rows'] if x['case']==case and x['noise_mode']==mode and x['method']==arm);col=1+2*j;draw(axs[row,col],m);draw(axs[row,col+1],m-gt,True);legal=r['validation']['legal'];axs[row,col].set_title(label+(' (outside box)' if not legal else ''),fontsize=9);axs[row,col+1].set_title(f"Error m−truth\nRMS {score['material_RMS']:.2e}; data {r['validation']['response_relative']:.1e}",fontsize=8)
   for ax in axs[row]:ax.set_xlabel('x (model length)',fontsize=6)
 fig.colorbar(plt.cm.ScalarMappable(norm=matnorm,cmap='viridis'),ax=axs[:,[0,1,3,5]].ravel().tolist(),shrink=.65,pad=.01,extend='both',label='real material, common fixed [0.5,4] scale')
 fig.colorbar(plt.cm.ScalarMappable(norm=errnorm,cmap='RdBu_r'),ax=axs[:,[2,4,6]].ravel().tolist(),shrink=.65,pad=.01,extend='both',label='signed error, common scale ±4')
 fig.suptitle('NEW high-frequency 48-pulse model: all 48 cells unknown; k = 20π; receiver condition 3.72×10⁴\nFrozen before GT score. Raw current SVD is accurate without noise and unstable at fixed 0.001 noise; bounded Full-J repairs noisy data.',fontsize=13)
 for ext in ['png','svg','pdf']:fig.savefig(P/f'figures/highfreq48_imaging.{ext}',dpi=160)
 plt.close(fig)
 # Separate cost/quality view prevents combining failed noisy raw SVD with successes.
 fig,axs=plt.subplots(1,3,figsize=(14,4),layout='constrained');colors=['#4477aa','#cc6677','#228833'];methods=['FullJ_midpoint','exactSVD','SVD_FullJ'];labels=['Full-J midpoint','Raw SVD','SVD+Full-J'];tot=scores['method_totals']
 for ax,mode in zip(axs[:2],['clean','relative_complex_1e-3']):
  for j,method in enumerate(methods):
   rows=[r for r in scores['rows'] if r['method']==method and r['noise_mode']==mode];ax.scatter([r['incremental_wall_seconds']*1000 for r in rows],[r['material_RMS'] for r in rows],color=colors[j],label=labels[j],s=42)
  ax.set(xscale='log',yscale='log',xlabel='incremental solve + validation (ms)',ylabel='material RMS',title='Clean' if mode=='clean' else 'Fixed relative noise 0.001');ax.grid(alpha=.2);ax.legend(fontsize=8)
 x=np.arange(3);setup=[tot[m]['shared_setup_wall_seconds'] for m in methods];incremental=[tot[m]['incremental_wall_seconds'] for m in methods];axs[2].bar(x,np.array(incremental)*1000,color=colors,label='6 incremental runs');axs[2].bar(x,np.array(setup)*1000,bottom=np.array(incremental)*1000,color='.65',label='one shared setup');axs[2].set(xticks=x,xticklabels=labels,ylabel='wall time ms',title='Setup once + six fixed inputs');axs[2].tick_params(axis='x',rotation=15);axs[2].legend(fontsize=8);fig.suptitle('Measured single tiny-system batch; raw SVD noisy quality FAIL; timings do not establish universal scaling')
 for ext in ['png','svg','pdf']:fig.savefig(P/f'figures/highfreq48_cost_quality.{ext}',dpi=160)
 plt.close(fig);print('Highfreq48 imaging and cost/quality figures written; no solver rerun.')
if __name__=='__main__':main()
