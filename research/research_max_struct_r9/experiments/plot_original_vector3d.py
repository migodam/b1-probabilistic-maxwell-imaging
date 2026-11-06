"""Offline fixed views of sealed original-scale predictions; zero new physics."""
import os
os.environ['MPLBACKEND']='Agg'
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
import json,time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
BASE=Path(__file__).resolve().parents[1];OUT=BASE/'RESULTS'/'original_vector3d';FIG=BASE/'FIGURES'
ARMS=['sourceCIE','sourceDual','GlobalDual','AplusGlobal','LSIEsource']
LABELS=['Source CIE','Source Dual','Global Dual','R9 vector A+','Source LSIE']

def main():
    wall,cpu=time.perf_counter(),time.process_time()
    # This guard precedes opening any truth. Scoring is separately authorized offline.
    scores=json.loads((OUT/'OFFLINE_SCORES.json').read_text())['scores']
    for arm in ARMS:
        if not (OUT/'inversions'/arm/'SEALED.json').exists():raise RuntimeError('All five predictions must be sealed first')
    truth=np.load(OUT/'sealed_truth'/'truth.npz')['inverse_material'];n=round(truth.size**(1/3))
    mats=[truth]+[np.load(OUT/'inversions'/a/'stage2.npz')['material'] for a in ARMS]
    summaries=[None]+[next(s for s in scores if s['arm']==a and s['name']=='stage2') for a in ARMS]
    coordinates=(np.arange(n)-(n-1)/2)*(3/n);slices=[int(np.argmin(abs(coordinates-z))) for z in [-.7,.7]]
    FIG.mkdir(exist_ok=True);plt.rcParams.update({'font.size':10,'axes.titlesize':11})
    fig,axes=plt.subplots(3,6,figsize=(17.5,8),layout='constrained',sharex=True,sharey=True)
    for col,(name,m,score) in enumerate(zip(['Truth']+LABELS,mats,summaries)):
        cube=m.reshape(n,n,n);images=[cube[:,:,i].T for i in slices]+[cube.max(axis=2).T]
        for row,image in enumerate(images):
            im=axes[row,col].imshow(image,origin='lower',extent=[-1.5,1.5,-1.5,1.5],cmap='viridis',vmin=0,vmax=4,interpolation='nearest')
            axes[row,col].set_xticks([-1,0,1]);axes[row,col].set_yticks([-1,0,1])
            if row==2:axes[row,col].set_xlabel('x / λ')
        title=name
        if score:title+='\nErr %.3f; physical RMS %.3f'%(score['paper_Err'],score['relative_source_RMS'])
        axes[0,col].set_title(title)
    for row,label in enumerate(['z / λ = %.2f'%coordinates[slices[0]],'z / λ = %.2f'%coordinates[slices[1]],'Maximum along z']):
        axes[row,0].set_ylabel(label+'\ny / λ')
    fig.colorbar(im,ax=axes,shrink=.8,label='Real contrast χ (common box 0–4)')
    fig.suptitle('PUBLIC DEVELOPMENT · original 3λ four cubes · 60³ data / 30³ inverse · independent CM+RR DDA',fontsize=13)
    for ext in ['png','pdf']:fig.savefig(FIG/('original_vector3d_materials.'+ext),dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    costs=[]
    for number,(arm,label) in enumerate(zip(ARMS,LABELS)):
        dest=OUT/'inversions'/arm
        hist=[json.loads(line) for line in (dest/'history.jsonl').read_text().splitlines()]
        ep=json.loads((dest/'ENDPOINTS.json').read_text());cost=json.loads((OUT/('invert_'+arm+'_COST.json')).read_text());costs.append(cost)
        color=plt.rcParams['axes.prop_cycle'].by_key()['color'][number]
        axes[0,0].plot([e['cumulative_CPU_seconds']/60 for e in ep],[e['relative_source_RMS'] for e in ep],'-o',label=label,color=color)
        axes[0,1].semilogy([h['cumulative_CPU_seconds']/60 for h in hist],[h['auxiliary_sum'] for h in hist],label=label,color=color)
        good=[h for h in hist if h['coefficient_delta_RMS'] is not None and h['coefficient_delta_RMS']>0]
        axes[1,0].semilogy([h['cumulative_CPU_seconds']/60 for h in good],[h['coefficient_delta_RMS'] for h in good],label=label,color=color)
    axes[0,0].set_title('Full physical data RMS: stage endpoints only');axes[0,0].set_ylabel('Mean-source relative RMS')
    axes[0,1].set_title('Auxiliary objective (arm/stage-specific β)');axes[0,1].set_ylabel('Relative data + current sum')
    axes[1,0].set_title('Literal paper Eq16 coefficient change');axes[1,0].set_ylabel('RMS relative γ change')
    axes[1,0].axhline(.001,color='black',linestyle=':',label='Paper stop 0.001')
    for ax in [axes[0,0],axes[0,1],axes[1,0]]:
        ax.set_xlabel('Actual cumulative CPU minutes');ax.grid(alpha=.2)
    x=np.arange(5);axes[1,1].bar(x-.18,[c['CPU_seconds']/60 for c in costs],width=.36,label='CPU')
    axes[1,1].bar(x+.18,[c['wall_seconds']/60 for c in costs],width=.36,label='Wall (contention dependent)')
    axes[1,1].set_xticks(x,LABELS,rotation=20);axes[1,1].set_ylabel('Minutes');axes[1,1].set_title('Inverse setup + auxiliary + endpoint validation')
    axes[0,0].legend(fontsize=9);axes[1,0].legend(fontsize=8);axes[1,1].legend(fontsize=9)
    fig.suptitle('PUBLIC DEVELOPMENT · no physical corrector · all stage endpoints retained · no GT iterate selection',fontsize=12)
    for ext in ['png','pdf']:fig.savefig(FIG/('original_vector3d_convergence.'+ext),dpi=180)
    plt.close(fig)
    (OUT/'PLOT_COST.json').write_text(json.dumps(dict(scope='OFFLINE_SEALED_PREDICTIONS; truth only for displayed reference and precomputed score labels',
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,physics_evaluations=0,
        views='fixed two z slices closest to ±0.7λ and max-z projection; common χ scale0..4',
        inputs='sealed final-stage materials, offline scores, paid histories and endpoint/cost receipts'),indent=2)+'\n')

if __name__=='__main__':main()
