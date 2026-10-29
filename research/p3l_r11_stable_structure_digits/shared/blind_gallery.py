from core import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 key={};rng=np.random.default_rng(11223307)
 for ds in ['N','P']:
  order=rng.permutation(9);key[ds]=[f'{ds}_{i:03d}' for i in order]
  for method in ['FULL256','AGENT_A','AGENT_B','AGENT_C']:
   fig,axs=plt.subplots(3,3,figsize=(9,10));fig.suptitle(f'{ds} — {method} — ANONYMIZED RECONSTRUCTIONS',fontsize=14)
   for j,(i,ax) in enumerate(zip(order,axs.ravel())):
    p=P/'results'/ds/method/f'{ds}_{i:03d}'/'PREDICTION.npz';a=np.load(p)['occupancy'].reshape(24,24)
    ax.imshow(a,cmap='gray_r',vmin=0,vmax=1,interpolation='nearest');ax.set_title(chr(65+j));ax.axis('off')
   fig.tight_layout();fig.savefig(P/'figures'/f'BLIND_{ds}_{method}.png',dpi=120);plt.close(fig)
 write(P/'audit/BLIND_KEY_SEALED.json',dict(order=key,labels_hidden=True,GT_files_read=False,created_utc=time.time(),note='Openonlyafterblindreadrecord'))
 write(P/'audit/ALL_PREDICTIONS_SEALED.json',dict(utc=time.time(),expected=72,files=[dict(path=str(p.relative_to(P)),blake2b=blake(p)) for p in (P/'results').glob('*/*/*/PREDICTION.npz')],GT_scoring_done=False))
if __name__=='__main__':main()
