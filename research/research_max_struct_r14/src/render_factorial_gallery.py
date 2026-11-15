"""Legible post-seal gallery with dedicated row colorbar columns; no selection."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import argparse,json,time,sys
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--history-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();r=a.root.resolve();h=a.history_root.resolve();out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    sys.path[:0]=[str(r/'src'),str(h/'src')]
    import score_dualmetric as base
    _,tasks,_,_=base.authorize(r,r/'audit/PREDICTION_SEAL_AUTHORIZED.json')
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import score_results as helper
    c=time.process_time();w=time.perf_counter()
    scores=json.loads((r/'analysis/exposed48_v2/ALL_VIEWS.json').read_text());scores={(x['case_id'],x['solver']):x for x in scores if x['view']=='primary'}
    cases=sorted({x['case_id'] for x in tasks});lookup={(x['case_id'],x['pair']):x for x in tasks};artifacts=[]
    with PdfPages(out/'all48_readable.pdf') as pdf:
        for page,start in enumerate(range(0,48,6),1):
            fig=plt.figure(figsize=(13.8,16.4),layout='constrained');grid=fig.add_gridspec(6,6,width_ratios=[1,1,1,1,1,.065],wspace=.025,hspace=.08)
            for i,cid in enumerate(cases[start:start+6]):
                truth=helper.load_truth(helper.truth_spec(h,dict(case_id=cid)))['m'];arrays=[truth]
                for pair in base.PAIRS:
                    task=lookup[(cid,pair)]
                    with np.load(r/'RESULTS'/task['campaign']/task['task_id']/'primary.npz',allow_pickle=False) as z:arrays.append(z['m'].copy())
                vmax=max(float(x.max()) for x in arrays)
                for j,m in enumerate(arrays):
                    ax=fig.add_subplot(grid[i,j]);image=ax.imshow(m.reshape(24,24),vmin=.02,vmax=vmax,cmap='viridis',interpolation='nearest')
                    ax.set_xticks([]);ax.set_yticks([])
                    if j==0:ax.set_ylabel(cid,fontsize=9)
                    if i==0:ax.set_title('Offline GT' if j==0 else 'pair'+base.PAIRS[j-1],fontsize=12)
                    if j:
                        s=scores[(cid,'pair'+base.PAIRS[j-1])];ax.set_xlabel(f"NRMSE {s['NRMSE']:.3f} · IoU {s['IoU']:.3f}",fontsize=9)
                cb=fig.colorbar(image,cax=fig.add_subplot(grid[i,5]));cb.ax.tick_params(labelsize=8)
            fig.suptitle('Exposed48: same observations, four metric pairs, public min-Φ primary\nχ = εr − 1; one shared material color scale per row; no interpolation',fontsize=14)
            name=f'materials_{page:02d}.png';fig.savefig(out/name,dpi=160);pdf.savefig(fig);plt.close(fig);artifacts.append(name)
    receipt=dict(status='READABLE_GALLERY_RENDERED',cases=48,selection_changed=False,
        source='analysis/exposed48_v2; original plotting outputs retained',change='Dedicated colorbar column fixes overlap; same original arrays/metrics/public selections',
        CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,new_forward=0,new_RHS=0,artifacts=artifacts+['all48_readable.pdf'])
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))


if __name__=='__main__':main()
