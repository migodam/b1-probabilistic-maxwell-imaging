"""Predetermined parent000 gallery of all eight sealed proposal records; no GT IO or physics."""
from pathlib import Path
import json,time,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[1];O=P/'results/candidates_fresh10';parent='R10BI_FRESH_000';ARMS=['FULL66_RML','BILATERAL_TARGET_RML_BLOCK','FIXED_Q8_RAY'];LABELS=['Full66 RML','Bilateral target','Fixed Q8'];sigma=.00119211798514519;gate=282.8737921622988

def load(path):
 with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def identity(path):
 p=Path(path);return dict(path=str(p),bytes=p.stat().st_size,blake2b_256=hashlib.blake2b(p.read_bytes(),digest_size=32).hexdigest())

def run():
 sealpath=O/'PREDICTIONS_SEALED.json';seal=json.loads(sealpath.read_text());assert len(seal['rows'])==480 and seal['GT_used'] is False
 ap=O/'parents'/parent/'anchor/ANCHOR.npz';anchor=load(ap);y=load(O/'inputs/online'/f'{parent}.npz')['observation'];aphi=float(np.sum(abs(anchor['prediction']-y)**2)/sigma**2);fig,axes=plt.subplots(3,9,figsize=(19,7));fig.subplots_adjust(left=.075,right=.99,top=.87,bottom=.14,wspace=.08,hspace=.34);metadata=[];im=None
 for ri,(arm,label) in enumerate(zip(ARMS,LABELS)):
  records=sorted([r for r in seal['rows'] if r['parent']==parent and r['arm']==arm],key=lambda x:x['k']);assert [r['k'] for r in records]==list(range(8));Y=axes[ri,0].get_position();fig.text(.008,(Y.y0+Y.y1)/2,label,fontsize=10,ha='left',va='center',rotation=90)
  for ci in range(9):
   ax=axes[ri,ci];ax.set_xticks([]);ax.set_yticks([]);path=ap if ci==0 else records[ci-1].get('path');record=None if ci==0 else records[ci-1];available=bool(path and Path(path).exists());qualified=bool(aphi<=gate) if ci==0 else bool(record.get('qualified'));actual=None
   if available:
    a=anchor if ci==0 else load(path);actual=float(np.sum(abs(a['prediction']-y)**2)/sigma**2);im=ax.imshow(a['material'].real.reshape(24,24),origin='lower',cmap='viridis',vmin=.1,vmax=1.9,interpolation='nearest');state='qualified' if qualified else 'REJECTED';title=('shared anchor' if ci==0 else f'k={ci-1}')+' | '+state+'\n'+f'Phi={actual:.3f}'
   else:
    ax.set_facecolor('#eeeeee');ax.text(.5,.5,'UNAVAILABLE\nradius bracket\nwithin latent cap',ha='center',va='center',transform=ax.transAxes,fontsize=8,color='#8b2020');title=f'k={ci-1} | unavailable\nPhi=not evaluated'
   color='#237443' if qualified else '#a62e2e';ax.set_title(title,fontsize=8,color=color,pad=3)
   for sp in ax.spines.values():sp.set_color(color);sp.set_linewidth(1.3)
   metadata.append(dict(parent=parent,arm=arm,column=ci,proposal_index=None if ci==0 else ci-1,kind='ANCHOR' if ci==0 else 'SEALED_PROPOSAL_ENDPOINT',qualified=qualified,available=available,Phi_actual_original_observation=actual,Phi_sealed=None if ci==0 else record.get('Phi'),status='SHARED_ANCHOR' if ci==0 else record['status'],source=None if not available else identity(path),selection='predetermined parent000; all records in original index order, no truth or oracle selection'))
 cax=fig.add_axes([.31,.065,.48,.022]);cb=fig.colorbar(im,cax=cax,orientation='horizontal');cb.set_label('Real susceptibility (one common physical scale for all panels)',fontsize=9);cb.ax.tick_params(labelsize=8);fig.suptitle('Finite material possibilities: fixed fresh parent000, shared anchor + all K8 records',fontsize=13,y=.975);fig.text(.5,.932,'Green: qualified. Red: rejected or unavailable. True-observation gate Phi <=282.874; old Fixed Q8 also retains radius/progress constraints.',ha='center',fontsize=9)
 out=P/'figures/candidate_fresh10_cloud_gallery.png';fig.savefig(out,dpi=180);fig.savefig(out.with_suffix('.pdf'));plt.close(fig);d=dict(utc=time.time(),source_script=identity(Path(__file__)),parent_selection='predetermined index000, no GT selection',GT_access=False,new_RHS=0,new_LU=0,color_limits=[.1,1.9],shape=[3,9],arms=ARMS,prediction_seal=identity(sealpath),panels=metadata,rejected_and_unavailable_retained=True,outputs=[str(out),str(out.with_suffix('.pdf'))]);(P/'figures/candidate_fresh10_cloud_gallery.metadata.json').write_text(json.dumps(d,indent=2)+'\n');print('GALLERY_COMPLETE',out)
if __name__=='__main__':run()
