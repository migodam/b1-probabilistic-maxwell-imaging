"""Postscore visualization ONLY: read private GT after prediction freeze.
No forward evaluations, solver calls, root reordering, selection, or tuning.
All independently qualified quotient materials appear in original solver order.
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
def read(path):return json.loads(path.read_text())
def main():
 if not (P/'results/general6_PREDICTION_FREEZE_RECEIPT.json').exists():raise RuntimeError('Need frozen predictions')
 if not (P/'results/general6_scores.json').exists():raise RuntimeError('GT visualization allowed only after scoring')
 contract=read(P/'experiments/GENERAL6_CONTRACT.json');scores=read(P/'results/general6_scores.json');summary=read(P/'results/general6_summary.json');truth=np.load(P/'private/general6/sealed_targets.npz');g=contract['geometry'];centers=np.asarray(g['centers']);norm=Normalize(.5,30);cmap=plt.colormaps['viridis'];methods=['general_quotient20','ReducedFullJ64'];pred={};qualified={};scoremap={}
 for i in range(3):
  case=f'GENERAL6_{i:02d}'
  for method in methods:
   pred[case,method]=read(P/f'results/general6_{case}_{method}.json');qualified[case,method]=[x for x in pred[case,method]['independently_validated'] if x['status']=='QUALIFIED'];scoremap[case,method]=next(x for x in scores['rows'] if x['case']==case and x['method']==method)
 if [len(qualified[f'GENERAL6_{i:02d}','general_quotient20']) for i in range(3)]!=[2,1,1]:raise RuntimeError('Expected frozen2/1/1 roots;do not silently select')
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.labelsize':9,'pdf.fonttype':42,'svg.fonttype':'none'})
 fig,axs=plt.subplots(3,4,figsize=(16,11.4),layout='constrained');mapaxes=[]
 def draw(ax,m,title,detail):
  m=np.asarray(m);coll=PatchCollection([Circle(x,g['radius']) for x in centers],cmap=cmap,norm=norm,edgecolors='#333333',linewidths=.65);coll.set_array(m);ax.add_collection(coll)
  for (x,y),v in zip(centers,m):
   color=cmap(norm(v));lum=.2126*color[0]+.7152*color[1]+.0722*color[2];ax.text(x,y,f'{v:.1f}',ha='center',va='center',color='white' if lum<.5 else '#111111',fontsize=7.5,fontweight='normal')
  ax.set(xlim=(-.43,.50),ylim=(-.36,.40),aspect='equal',xlabel='x (model length)',ylabel='y (model length)');ax.set_xticks([-.3,0,.3]);ax.set_yticks([-.3,0,.3]);ax.tick_params(labelsize=8);ax.grid(alpha=.15);ax.set_title(title+'\n'+detail,pad=7);mapaxes.append(ax)
 def drawroot(ax,case,method,index,title):
  root=qualified[case,method][index];rms=scoremap[case,method]['all_qualified_root_RMS_to_hidden_truth'][index];res=root['validation']['response_relative'];draw(ax,root['m'],title,f'data residual {res:.1e}  |  GT RMS {rms:.2e}')
 # Case00: truth, EVERY quotient root in returned order, and the baseline root.
 case='GENERAL6_00';draw(axs[0,0],truth[case],'(a) Target 00: hidden truth','Postscore display; six independent unknowns');drawroot(axs[0,1],case,'general_quotient20',0,'(b) Quotient root 1 — returned order');drawroot(axs[0,2],case,'general_quotient20',1,'(c) Quotient root 2 — returned order');drawroot(axs[0,3],case,'ReducedFullJ64',0,'(d) Reduced Full-J: one returned root')
 for row in [1,2]:
  case=f'GENERAL6_{row:02d}';draw(axs[row,0],truth[case],f'Target 0{row}: hidden truth','Postscore display; same frozen geometry');drawroot(axs[row,1],case,'general_quotient20',0,'Quotient: one returned physical root');drawroot(axs[row,2],case,'ReducedFullJ64',0,'Reduced Full-J: one returned root')
 # Acquisition geometry and frozen costs are secondary to the material maps.
 ax=axs[1,3]
 for i,(x,y) in enumerate(centers):ax.add_patch(Circle((x,y),g['radius'],fc='#4477aa',ec='#222222',lw=.5));ax.text(x+.055,y,str(i+1),fontsize=8)
 src=np.asarray(g['sources']);rec=np.asarray(g['receivers']);ax.scatter(*src.T,marker='*',s=120,c='#228833',label='3 labeled sources');ax.scatter(*rec.T,marker='s',s=55,c='#ee7733',label='1 complex receiver')
 for i,(x,y) in enumerate(src):ax.annotate(f's{i+1}',(x,y),xytext=(5,4),textcoords='offset points',fontsize=8)
 ax.set(xlim=(-.85,.85),ylim=(-.85,1.0),aspect='equal',xlabel='x (model length)',ylabel='y (model length)',title='Frozen acquisition: 6 unknowns, 6 real data\nk = 2π; disk radius 0.04; box [0.5, 30]');ax.grid(alpha=.15);ax.annotate('R',tuple(rec[0]),xytext=(6,-2),textcoords='offset points',fontsize=8)
 ax=axs[2,3];entries={r['method']:r for r in summary['methods']};names=['Quotient 20','Reduced Full-J 64'];cost=[entries[m]['total_wall_including_logging'] for m in methods];bars=ax.bar(names,cost,color=['#4477aa','#aa6677'],width=.55)
 for b,v in zip(bars,cost):ax.text(b.get_x()+b.get_width()/2,v+.025,f'{v:.5f} s',ha='center',fontsize=10)
 ax.set(ylim=(0,1.83),ylabel='Measured wall seconds',title='Frozen three-target batch, setup included');ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True);ax.text(.03,.95,'Returned root counts\nQuotient: 2 / 1 / 1\nFull-J 64: 1 / 1 / 1',transform=ax.transAxes,ha='left',va='top',fontsize=9)
 d=float(np.linalg.norm(np.asarray(qualified['GENERAL6_00','general_quotient20'][0]['m'])-np.asarray(qualified['GENERAL6_00','general_quotient20'][1]['m']))/np.sqrt(6))
 fig.suptitle('Blind six-cell material recovery: an additional, distinct real material fits the same observations\nAll returned physical roots shown; labels preserve solver order. GT used only after the frozen prediction set was scored.',fontsize=15)
 cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=mapaxes,orientation='horizontal',shrink=.75,pad=.025,aspect=60);cb.set_label(f'Real material, common fixed [0.5, 30] scale  |  Target 00 root-to-root RMS = {d:.6f}');cb.set_ticks([.5,5,10,15,20,25,30])
 for ext in ['png','pdf','svg']:fig.savefig(P/f'figures/general6_blind_roots.{ext}',dpi=200)
 plt.close(fig);print('GENERAL6 all physical roots plotted in solver order. Private GT read only in postscore visualization; no physics/solver rerun.')
if __name__=='__main__':main()
