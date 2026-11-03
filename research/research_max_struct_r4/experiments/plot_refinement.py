"""Post-seal refinement comparison. No physical computations or selection."""
from pathlib import Path
import json,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
def main():
    tick=time.perf_counter()
    dev=json.loads((ROOT/'results/DEVELOPMENT_SCORES.json').read_text())['rows']
    conf=json.loads((ROOT/'results/CONFIRMATION_SCORES.json').read_text())['rows']
    add=json.loads((ROOT/'results/DEVELOPMENT_POST_CONFIRMATION_ADDENDA_V1.json').read_text())['rows']
    cases=['obs_b20d5a55677d4eca','obs_42cb233e4f224bb3']
    fig,axes=plt.subplots(2,3,figsize=(9,6.7),layout='constrained')
    for i,case in enumerate(cases):
        truthroot=ROOT.parent/'research_max_struct_r2/sealed_truth' if i==0 else ROOT/'sealed_truth/confirmation'
        original=ROOT/'results/onehop_absorption/prediction.npz' if i==0 else ROOT/'results/confirmation_remote/predictions'/case/'onehop.npz'
        original_score=next(r for r in dev if r['arm']=='onehop_absorption') if i==0 else next(r for r in conf if r['observation_id']==case and r['method']=='onehop')
        avg_score=next(r for r in add if r['observation_id']==case and r['method'].startswith('partition_average'))
        paths=[truthroot/(case+'.npz'),original,ROOT/'results/partition_average'/case/'prediction.npz']
        for j,path in enumerate(paths):
            m=np.load(path,allow_pickle=False)['m'];ax=axes[i,j]
            im=ax.imshow(m.reshape(24,24),vmin=0,vmax=20,origin='lower',cmap='viridis',interpolation='nearest')
            if j==0:label=('Original development' if i==0 else 'Post-confirmation development')+'\nTruth'
            else:
                score=original_score if j==1 else avg_score
                label=('Fixed partition' if j==1 else 'Four-shift average')+'\nNRMSE %.3f | residual %.4f'%(score['amplitude_normalized_nrmse'],score['relative_residual'])
            ax.set_title(label,fontsize=10);ax.set_xticks([]);ax.set_yticks([])
    fig.colorbar(im,ax=axes,label='Physical real contrast',shrink=.8)
    fig.suptitle('Partition-origin repair: all two frozen development attempts retained')
    for ext in ('png','pdf'):fig.savefig(ROOT/'figures'/('partition_refinement.'+ext),dpi=160)
    plt.close(fig)
    (ROOT/'audit/REFINEMENT_PLOT_RECEIPT.json').write_text(json.dumps(dict(status='RENDERED',truth_access='POST_SEAL_EXPOSED_DEVELOPMENT',new_physical_LU=0,new_complex_RHS=0,wall_seconds=time.perf_counter()-tick,selection='none; both prescribed cases shown'),indent=2))
if __name__=='__main__':main()
