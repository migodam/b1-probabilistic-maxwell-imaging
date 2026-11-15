"""Public-operator retained energies; no truth, forward, optimizer or SVD."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import csv,json,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    root=Path(__file__).resolve().parents[1];out=root/'analysis/spectral_selection_v1';out.mkdir(parents=True,exist_ok=False)
    c=time.process_time();w=time.perf_counter();done=set();rows=[];plots=[]
    for p in sorted((root/'audit/rank_all24/complex').glob('EXECUTION_*.json')):
        meta=json.loads(p.read_text());group=meta['operator_group']
        if group in done:continue
        done.add(group);case=Path(meta['public_path']).stem
        with np.load(root/'data_public/independent24'/(case+'.npz'),allow_pickle=False) as z:G=z['G'];S=z['S']
        with np.load(p.parent/group/'SPACES.npz',allow_pickle=False) as z:a={k:z[k].copy() for k in z.files}
        norms={'G':float(np.vdot(G,G).real),'S':float(np.vdot(S,S).real)}
        choices={'G19':a['G_right_frame'][:,:19],'S19':a['S_right_frame'][:,:19],**{k:a[k] for k in ('U11','U10','U01','U00')}}
        for name,U in choices.items():
            GU=G@U;SU=S@U;g=float(np.vdot(GU,GU).real);s=float(np.vdot(SU,SU).real)
            gal=U.conj().T@GU if name in ('G19','S19','U11') else None
            rows.append(dict(operator_group=group,representative_case=case,space=name,rank=U.shape[1],
                G_input_energy=g,G_total_energy=norms['G'],G_input_energy_fraction=g/norms['G'],
                S_input_energy=s,S_total_energy=norms['S'],S_input_energy_fraction=s/norms['S'],
                G_Galerkin_energy=None if gal is None else float(np.vdot(gal,gal).real),
                G_Galerkin_energy_fraction=None if gal is None else float(np.vdot(gal,gal).real)/norms['G'],
                interpretation='Background current operator Frobenius energy; not material-Jacobian information, not a basin predictor',truth_access='NONE'))
        plots.append((group,a))
    with (out/'RETAINED_ENERGIES.csv').open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    fig,axes=plt.subplots(len(plots),2,figsize=(13,8),squeeze=False)
    for i,(group,a) in enumerate(plots):
        for key,color in [('G','C0'),('S','C1')]:
            sigma=a[key+'_singular_values'];x=np.arange(1,len(sigma)+1)
            axes[i,0].semilogy(x,sigma/sigma[0],label=key,color=color)
            axes[i,1].plot(x,np.cumsum(sigma**2)/np.sum(sigma**2),label=key,color=color)
        axes[i,0].axhline(.1,color='gray',ls=':');axes[i,0].set_ylim(1e-8,1.2)
        for ax in axes[i]:ax.axvline(19,color='black',ls='--',lw=.8);ax.set_xlim(1,80);ax.grid(alpha=.2);ax.legend();ax.set_xlabel('Independent descending singular index of each operator')
        axes[i,0].set_title(group+': singular values / own maximum');axes[i,1].set_title(group+': cumulative squared singular values')
    fig.suptitle('Public current operators: matching singular indices does not match directions')
    fig.tight_layout();figs=root/'FIGURES/geometry';figs.mkdir(parents=True,exist_ok=True)
    for suffix in ('png','pdf'):fig.savefig(figs/('background_spectra_and_energy.'+suffix),dpi=160,bbox_inches='tight')
    plt.close(fig)
    receipt=dict(status='COMPLETE',operator_groups=len(plots),rows=len(rows),truth_access='NONE',physical_LU=0,physical_complex_RHS=0,new_SVD=0,
        CPU_seconds=time.process_time()-c,whole_process_CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-w,
        scope='Actual cached-frame energy accounting, no repeated identity qualification or algorithm tuning')
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt))

if __name__=='__main__':main()
