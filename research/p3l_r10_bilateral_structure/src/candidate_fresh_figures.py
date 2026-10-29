"""Post-seal, post-score scientific figures; no physics."""
from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parent));import candidate_fresh as cf
P=cf.P;OUT=cf.OUT;pv=cf.pv

def run():
 stats=json.loads((OUT/'STATISTICS.json').read_text());sealed=json.loads((OUT/'PREDICTIONS_SEALED.json').read_text());rows=stats['rows'];arms=cf.ALLARMS;colors=plt.get_cmap('tab10').colors;labels=['Full66 RML','Shape target','Bilateral target','Fixed Q8','Dynamic Q8','Random66 ray'];fig,ax=plt.subplots(2,2,figsize=(13,8))
 for i,(arm,label) in enumerate(zip(arms,labels)):
  rs=[r for r in rows if r['arm']==arm and r['K']==8];xx=np.arange(10)+(i-2.5)*.04
  ax[0,0].plot(xx,[r['oracle_RMS'] for r in rs],'o-',lw=.8,ms=3,color=colors[i],label=label);ax[0,1].plot(xx,[r['uniform_action_MSE'] for r in rs],'o-',lw=.8,ms=3,color=colors[i]);ax[1,0].plot(xx,[r['qualified_rate'] for r in rs],'o-',lw=.8,ms=3,color=colors[i]);ax[1,1].plot(xx,[r['pair_diameter'] for r in rs],'o-',lw=.8,ms=3,color=colors[i])
 ax[0,0].axhline(.075,ls='--',color='black',lw=1);ax[0,0].set_ylabel('Oracle material RMS (offline)');ax[0,1].set_ylabel('Uniform finite-cloud action MSE');ax[1,0].set_ylabel('Qualified proposal fraction');ax[1,1].set_ylabel('Qualified cloud pair diameter (RMS)')
 for a in ax.ravel():a.set_xlabel('Independent fresh parent index');a.set_xticks(range(10));a.grid(alpha=.25)
 ax[0,0].legend(fontsize=7,ncol=2);fig.suptitle('Fresh10, K8: bounded mechanism check; methods have different cost and radius constraints',fontsize=11);fig.tight_layout();fig.savefig(P/'figures/candidate_fresh10_parent_metrics.png',dpi=170);plt.close(fig)
 fig,ax=plt.subplots(1,2,figsize=(11,4))
 for i,(arm,label) in enumerate(zip(arms,labels)):
  series=[]
  for K in (1,2,4,8):
   rs=[r for r in rows if r['arm']==arm and r['K']==K];series.append((np.mean([r['RHS256'] for r in rs]),np.mean([r['oracle_MSE'] for r in rs if r['oracle_MSE'] is not None]),np.mean([r['uniform_action_MSE'] for r in rs])))
  a=np.array(series);ax[0].plot(a[:,0],a[:,1],'o-',color=colors[i],label=label);ax[1].plot(a[:,0],a[:,2],'o-',color=colors[i])
 for a in ax:a.set_xlabel('Mean actual proposal RHS256 per parent (excludes shared anchor)');a.grid(alpha=.25)
 ax[0].set_ylabel('Mean offline oracle MSE');ax[1].set_ylabel('Mean uniform action MSE');ax[0].legend(fontsize=7);fig.suptitle('Prefix K1/2/4/8: descriptive cost curves, no equal-budget interpolation',fontsize=11);fig.tight_layout();fig.savefig(P/'figures/candidate_fresh10_cost_prefix.png',dpi=170);plt.close(fig)
 # Preselected fresh parent000, actions determined without truth; color scale fixed by prior physical range.
 parent='R10BI_FRESH_000';truth=pv.load_npz(cf.PRIVATE/'truth'/f'{parent}.npz')['material'];anchor=pv.load_npz(OUT/'parents'/parent/'anchor/ANCHOR.npz')['material'];fields=[truth,anchor]
 for arm in ('FULL66_RML','BILATERAL_TARGET_RML_BLOCK'):
  X=[anchor] if sealed['anchors'][parent]['Phi']<=282.8737921622988 else []
  for r in sealed['rows']:
   if r['parent']==parent and r['arm']==arm and r.get('qualified') and r.get('path'):
    m=pv.load_npz(r['path'])['material']
    if not X or min(np.linalg.norm(m-b)/24 for b in X)>1e-9:X.append(m)
  fields.append(np.mean(X,axis=0) if X else anchor)
 fig,ax=plt.subplots(1,4,figsize=(12,3));titles=['Simulated truth','Shared anchor','Full66 uniform action','Bilateral uniform action']
 for a,m,title in zip(ax,fields,titles):im=a.imshow(m.real.reshape(24,24),origin='lower',vmin=.1,vmax=1.9,cmap='viridis');a.set_title(title,fontsize=9);a.set_xticks([]);a.set_yticks([])
 fig.colorbar(im,ax=ax.ravel().tolist(),shrink=.75,label='Real susceptibility');fig.suptitle('Fixed fresh parent000; actions use qualified clouds without truth selection',fontsize=10);fig.savefig(P/'figures/candidate_fresh10_material_actions.png',dpi=170,bbox_inches='tight');plt.close(fig)
 print('FIGURES_COMPLETE',flush=True)
if __name__=='__main__':run()
