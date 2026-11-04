"""Original scientific figures from completed receipts; no GT or physics calls."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,time
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]

def main():
    start,cpu=time.perf_counter(),time.process_time();out=ROOT/'FIGURES/mechanisms';out.mkdir(parents=True,exist_ok=True)
    ids=['obs_40ac5b79e8a74ac3','obs_82ac794e228243bb','obs_42cb233e4f224bb3']
    labels=['Slit shell: reported gradient stop','Thin ring: inherited budget stop','42cb: inherited budget stop']
    fig,axes=plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
    summary=[]
    for ax,ident,label in zip(axes,ids,labels):
        small=json.loads((ROOT/'RESULTS/rational'/ident/'RECEIPT.json').read_text())
        large=json.loads((ROOT/'RESULTS/rational_coherent'/ident/'RECEIPT.json').read_text())
        for entry in small['validations']:
            c=json.loads((ROOT/'RESULTS/rational'/ident/entry['label']/'CORRECTOR.json').read_text())
            h=np.array(c['history']);ax.plot(np.arange(len(h)),np.minimum.accumulate(h[:,0]),label=entry['label'].replace('_',' '),linewidth=1.1)
        c=json.loads((ROOT/'RESULTS/rational_coherent'/ident/'CORRECTOR.json').read_text());h=np.array(c['history'])
        ax.plot(np.arange(len(h)),np.minimum.accumulate(h[:,0]),label='rank32 finite',linewidth=1.2,linestyle='--')
        ax.set_title(label,fontsize=11);ax.set_xlabel('Paid physical corrector evaluations');ax.set_ylabel('Original Phi incumbent (public selection)');ax.grid(alpha=.2)
        summary.append(dict(case_id=ident,anchor=small['anchor_original_Phi'],rank32_best_finite=large['best_finite_Phi'],rank32_after_corrector=large['corrected_Phi'],lowest_free_Hessian=large['lowest_free_physical_Hessian_eigenvalues']))
    axes[1].legend(fontsize=8);fig.suptitle('Exact finite-action search vs continuing the same disclosed physical anchor (equal corrector caps; costs differ)')
    fig.savefig(out/'rational_vs_continuation.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for k,(ax,name) in enumerate(zip(axes,['rational','rational_coherent'])):
        rows=json.loads((ROOT/'RESULTS'/name/ids[0]/'SEARCH_RECORDS.json').read_text())
        chosen=min(rows,key=lambda r:r['best_Phi']);index=chosen['action_id']
        if name=='rational':
            with np.load(ROOT/'RESULTS'/name/ids[0]/'OBJECTIVE_MAPS.npz') as z:x=z[f'x_{index}'];v=z[f'Phi_{index}']
        else:
            with np.load(ROOT/'RESULTS'/name/ids[0]/f'map_{index:02d}.npz') as z:x=z['x'];v=z['Phi']
        tri=ax.tricontourf(x[:,0],x[:,1],np.log10(np.maximum(v,1e-16)),levels=24,cmap='viridis');fig.colorbar(tri,ax=ax,label='log10(original Phi)')
        ax.scatter(0,0,c='red',s=25,label='original anchor');ax.scatter(*chosen['best_x'],c='white',marker='x',s=40,label='selected finite endpoint')
        ax.set_title(f"Slit shell: selected rank{chosen['state_rank']} map");ax.set_xlabel('Legal amplitude 1');ax.set_ylabel('Legal amplitude 2');ax.legend(fontsize=8)
    fig.suptitle('Search maps show the tested finite family, not a certificate over all legal materials')
    fig.savefig(out/'finite_action_objective_maps.png',dpi=180);plt.close(fig)
    # Geometry figure is a schematic of the new conditional formula, not a physical experiment.
    b=.28-.38j;m=np.linspace(.02,20,500);z=m/(1+b*m);target=1.7+1.1j
    candidates=[0,-1,np.argmin(abs(z-target))];fig,ax=plt.subplots(figsize=(5.5,4.8),layout='constrained');ax.plot(z.real,z.imag,label='legal real-material Mobius arc')
    ax.scatter([target.real],[target.imag],marker='x',c='black',label='unconstrained state LS point')
    ax.scatter(z[candidates].real,z[candidates].imag,c=['grey','grey','red'],label='endpoints / closest arc point');ax.set_aspect('equal');ax.grid(alpha=.2);ax.set_xlabel('Re z');ax.set_ylabel('Im z');ax.legend(fontsize=8);ax.set_title('Coordinate geometry only; prior-MM acceptance is additional')
    fig.savefig(out/'mobius_coordinate_schematic.png',dpi=180);plt.close(fig)
    (ROOT/'audit/MECHANISM_FIGURE_RECEIPT.json').write_text(json.dumps(dict(wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,new_physical_calls=0,GT_access='NONE',summary=summary),indent=2)+'\n')

if __name__=='__main__':main()
