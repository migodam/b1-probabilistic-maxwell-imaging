"""Plot the exposed public-state witness; no optimizer, truth or new forward."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import json,time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    cpu,wall=time.process_time(),time.perf_counter()
    root=Path(__file__).resolve().parents[1]
    witness=root/'delegated/full_current/same_shadow_witness'
    receipt=json.loads((witness/'WITNESS.json').read_text())
    with np.load(witness/'STATES.npz') as z:states={k:z[k].copy() for k in z.files}
    paths=[root/c['run'] for c in receipt['cases']]
    with np.load(paths[0]/'EVENT_GEOMETRY.npz') as z:
        thresholds=z['thresholds'].copy()
    first=np.sort(thresholds[np.isfinite(thresholds)])[:2]
    final=[]
    for p in paths:
        with np.load(p/'primary.npz') as z:final.append(z['m'].copy())
    fig,axes=plt.subplots(2,3,figsize=(14.2,9.1))
    ax=axes[0,0]
    for t in first:ax.axvline(t,color='gray',ls=':',lw=1)
    for c,color in zip(receipt['cases'],['tab:blue','tab:orange']):
        ax.scatter([c['target_t']],[0],color=color,s=70)
    ax.annotate('0.1% event overshoot',(receipt['cases'][0]['target_t'],0),xytext=(-30,52),textcoords='offset points',fontsize=9,arrowprops={'arrowstyle':'-'})
    ax.annotate('first-interval midpoint',(receipt['cases'][1]['target_t'],0),xytext=(-22,-54),textcoords='offset points',fontsize=9,arrowprops={'arrowstyle':'-'})
    ax.set_xlim(first[0]-.001,first[1]+.001);ax.set_ylim(-.3,.3);ax.set_yticks([]);ax.set_xlabel('normal-ray parameter t')
    ax.set_title('Both cross the same first event',fontsize=10)
    for j,key in enumerate(['fractional_m','midpoint_m'],1):
        im=axes[0,j].imshow(states[key].reshape(24,24),origin='lower',vmin=0,vmax=20,cmap='viridis')
        axes[0,j].set_title(('Event overshoot' if j==1 else 'Midpoint')+' material shadow',fontsize=10)
    difference=np.sqrt(np.sum(abs(states['fractional_X']-states['midpoint_X'])**2+abs(states['fractional_Q']-states['midpoint_Q'])**2,axis=1))
    dm=axes[1,0].imshow(difference.reshape(24,24),origin='lower',cmap='magma')
    axes[1,0].set_title('Cellwise full-state difference\nrelative Frobenius distance = 0.0552%',fontsize=10)
    fig.colorbar(dm,ax=axes[1,0],fraction=.04,pad=.03,label='scaled-state norm per cell')
    for j,(m,c) in enumerate(zip(final,receipt['cases']),1):
        axes[1,j].imshow(m.reshape(24,24),origin='lower',vmin=0,vmax=20,cmap='viridis')
        axes[1,j].set_title('Future paid primary\nresidual / delta = '+f"{c['future_primary']['residual_over_delta']:.3f}",fontsize=10)
    for ax in [axes[0,1],axes[0,2],*axes[1]]:ax.set_xticks([]);ax.set_yticks([])
    fig.suptitle('Exposed B__r12_confirm24_0010: material shadow alone does not determine the future',fontsize=13,y=.985)
    fig.text(.04,.025,'Max material-shadow difference = 9.64e-14. Same reference arrays; no GT used here.\nDifferent subsequent full-state paths; this is not a material nonuniqueness or chaos theorem.',fontsize=10)
    fig.subplots_adjust(left=.04,right=.91,top=.90,bottom=.14,wspace=.50,hspace=.48)
    cax=fig.add_axes([.94,.22,.012,.57])
    fig.colorbar(im,cax=cax,label='real contrast chi; same [0,20] scale')
    out=root/'analysis/mechanism_public/FIGURES';out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/'same_shadow_different_future.png',dpi=180);fig.savefig(out/'same_shadow_different_future.pdf');plt.close(fig)
    (out.parent/'PLOT_RECEIPT.json').write_text(json.dumps(dict(status='EXPOSED_PUBLIC_DIAGNOSTIC',truth_access='NONE',new_Maxwell_LU=0,source=str(witness.relative_to(root)),CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall),indent=2)+'\n')
    (out.parent/'REPORT.md').write_text('# Same material shadow, different future\n\nExposed public-state diagnostic. Figure uses already-paid physical candidates selected by original Phi; it does not use material GT. Exact reference arrays, material-shadow difference9.64e-14 and full-state difference0.0552% are retained in the witness source. No claim of physical material nonuniqueness, chaos or a global basin theorem.\n')
if __name__=='__main__':main()
