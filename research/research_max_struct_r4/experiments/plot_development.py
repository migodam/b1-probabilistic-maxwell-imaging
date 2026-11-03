"""Post-seal figures only. Uses disclosed development truth; no physics."""
from pathlib import Path
import json,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
def main():
    tick=time.perf_counter();out=ROOT/'figures';out.mkdir(exist_ok=True)
    rows=json.loads((ROOT/'results/DEVELOPMENT_SCORES.json').read_text())['rows']
    scores={row['arm']:row for row in rows}
    truth=np.load(ROOT.parent/'research_max_struct_r2/sealed_truth/obs_b20d5a55677d4eca.npz')['m']
    specs=[('Truth',None),('Native CIE / 8000','native_extended'),('Global CIE / 1500','global_common_corrector'),
        ('Self only','self_absorption'),('A: exact block','A_block_absorption'),('A: one exchange','onehop_absorption'),
        ('B: finite port profile','B_lambda1_port'),('C: constrained finite ray','box_cooperative')]
    fig,axes=plt.subplots(2,4,figsize=(12,6.7),layout='constrained')
    for ax,(label,key) in zip(axes.flat,specs):
        m=truth if key is None else np.load(ROOT/scores[key]['prediction'])['m']
        im=ax.imshow(m.reshape(24,24),origin='lower',vmin=0,vmax=20,cmap='viridis',interpolation='nearest')
        title=label if key is None else label+'\nNRMSE %.3f | IoU %.3f'%(scores[key]['amplitude_normalized_nrmse'],scores[key]['otsu_iou'])
        ax.set_title(title,fontsize=10);ax.set_xticks([]);ax.set_yticks([])
    fig.colorbar(im,ax=axes,label='Physical real contrast',shrink=.8)
    fig.suptitle('One disclosed development object — full physical endpoint, identical color scale',fontsize=13)
    fig.savefig(out/'development_materials.png',dpi=170);fig.savefig(out/'development_materials.pdf');plt.close(fig)
    data=json.loads((ROOT/'audit/DEVELOPMENT_STAGE_DIAGNOSTIC.json').read_text())['rows']
    fig,axes=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
    for arm,label in [('onehop_absorption','A: one exchange'),('block_absorption','A: exact block'),('self_absorption','Self only'),('global_common_corrector','Global CIE')]:
        r=[x for x in data if x['arm']==arm and x['label']!='historical_R2_corrector_reference']
        axes[0].plot(range(5),[x['amplitude_normalized_nrmse'] for x in r],'o-',label=label)
        axes[1].semilogy(range(5),[x['relative_physical_residual'] for x in r],'o-',label=label)
    for ax in axes:
        ax.set_xticks(range(5),['Host','mf3','mf6','mf10','Corrector']);ax.grid(alpha=.25)
    axes[0].set_ylabel('Material RMSE / 12');axes[1].set_ylabel('Original Maxwell relative residual')
    axes[0].legend(fontsize=9);fig.suptitle('OFFLINE_DIAGNOSTIC: material recovery and residual can move in opposite directions')
    fig.savefig(out/'development_paths.png',dpi=170);fig.savefig(out/'development_paths.pdf');plt.close(fig)
    (ROOT/'audit/DEVELOPMENT_PLOT_RECEIPT.json').write_text(json.dumps(dict(status='RENDERED',truth_access='POST_SEAL_DISCLOSED_DEVELOPMENT',new_physical_LU=0,wall_seconds=time.perf_counter()-tick,figures=['development_materials','development_paths']),indent=2))
if __name__=='__main__':main()
