"""Post-freeze visualization process; intentionally reads truth only for plotting.
Never import from solver processes. No solver is run or tuned here.
"""
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
 if not (P/'experiments/FREEZE_RECEIPT.json').exists() or not (P/'results/blind_scores.json').exists():raise RuntimeError('Plot only after prediction freeze and scoring')
 truth=np.load(P/'private/sealed_targets.npz');c=json.loads((P/'experiments/contract.json').read_text());scores=json.loads((P/'results/blind_scores.json').read_text());score={(r['case'],r['method']):r for r in scores['rows']}
 cases=[f'INFO12_{i:02d}' for i in range(3)]+[f'E1_LOW6_{i:02d}' for i in range(3)]
 fig,axs=plt.subplots(6,7,figsize=(20,15),layout='constrained');methods=[('receiver_svd','Receiver SVD'),('FullJ_adjoint_TRF','Exact Full-J TRF'),('frozen_J_TRF','Frozen-J TRF')];material_norm=Normalize(.5,4);error_norm=Normalize(-1.5,1.5)
 def draw(ax,centers,radius,values,error=False):
  coll=PatchCollection([Circle(x,radius) for x in centers],cmap='RdBu_r' if error else 'viridis',norm=error_norm if error else material_norm,edgecolors='.25',linewidths=.5);coll.set_array(values);ax.add_collection(coll);ax.set(xlim=(-.45,.5),ylim=(-.38,.4),aspect='equal');ax.set_xticks([-.3,0,.3]);ax.set_yticks([-.3,0,.3]);ax.tick_params(labelsize=6);ax.grid(alpha=.12);return coll
 for row,case in enumerate(cases):
  preset=case.rsplit('_',1)[0];g=c['presets'][preset];centers=np.array(g['centers']);gt=truth[case];draw(axs[row,0],centers,.04,gt);axs[row,0].set_title(f'{case}\nHidden truth',fontsize=10);axs[row,0].set_ylabel('y / wavelength',fontsize=7)
  for idx,(method,label) in enumerate(methods):
   col=1+2*idx;f=P/'results'/f'blind_{case}_{method}.json';r=json.loads(f.read_text())
   if 'm' not in r:
    for ax in axs[row,col:col+2]:ax.axis('off');ax.text(.5,.5,'Rank deficient\nNo candidate',ha='center',va='center',transform=ax.transAxes,fontsize=9,color='.4')
    axs[row,col].set_title(label,fontsize=10);continue
   m=np.array(r['m']);draw(axs[row,col],centers,.04,m);draw(axs[row,col+1],centers,.04,m-gt,True);axs[row,col].set_title(label,fontsize=10);rms=score[case,method]['material_RMS'];resp=r['validation']['response_relative'];axs[row,col+1].set_title(f'Error m−truth\nRMS {rms:.2e}; data {resp:.1e}',fontsize=8)
  for ax in axs[row]:ax.set_xlabel('x / wavelength',fontsize=7)
 material_map=plt.cm.ScalarMappable(norm=material_norm,cmap='viridis');error_map=plt.cm.ScalarMappable(norm=error_norm,cmap='RdBu_r');fig.colorbar(material_map,ax=axs[:,[0,1,3,5]].ravel().tolist(),shrink=.7,pad=.015,label='real material contrast (common scale)');fig.colorbar(error_map,ax=axs[:,[2,4,6]].ravel().tolist(),shrink=.7,pad=.015,label='signed material error (common scale)')
 fig.suptitle('Frozen blind predictions scored after freeze | INFO12: exact identifiable imaging; E1: response fits do not recover hidden material\nNo noise; independent state/response validation; numerical-error panels use the same ±1.5 scale as underobserved E1',fontsize=13)
 (P/'figures').mkdir(exist_ok=True)
 for ext in ['png','svg','pdf']:fig.savefig(P/f'figures/blind_imaging.{ext}',dpi=160)
 plt.close(fig)
 print('Wrote frozen-score GT/reconstruction/error figure in PNG, SVG, PDF; no solver rerun.')
if __name__=='__main__':main()
