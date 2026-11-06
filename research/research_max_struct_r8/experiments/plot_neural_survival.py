"""Post-seal visualization of every neural initializer and paid endpoints.

No solver, training, GT-based selection or new metric definition.
"""
from pathlib import Path
import json,time
CPU=time.process_time();WALL=time.perf_counter()
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
load=lambda p:json.loads(p.read_text())

def main():
    original=load(ROOT/'RESULTS/warm_scoring/RECEIPT.json')
    anchor=load(ROOT/'delegated/warm_anchor_analysis/CONFIRMATION_RECEIPT.json')
    if original['status']!='SCORED_ALL_FROZEN_PREDICTIONS' or anchor['status']!='SEPARATE126_ANCHOR_OFFLINE_SCORED':
        raise RuntimeError('Both complete offline scorings required before visualizing test truth')
    rows=load(ROOT/'RESULTS/warm_scoring/SCORES.json')['rows'];by={(x['case_id'],x['arm'],x['initialization']):x for x in rows}
    tasks={t['task_id']:t for t in load(ROOT/'configs/warm_all_manifest.json')['tasks']}
    out=ROOT/'RESULTS/neural_survival_visualization';out.mkdir(exist_ok=False)
    figdir=ROOT/'FIGURES';fixed={f'r8_test_{i:04d}' for i in (0,3,24,27,36,39)}
    arms=[('physical_lbfgs','uniform'),('physical_lbfgs','nn'),('a_plus','nn'),('a_plus_sourcewise','nn')]
    with PdfPages(figdir/'neural_initial_to_final_all42.pdf') as pdf:
        for case in sorted({x['case_id'] for x in rows}):
            with np.load(ROOT/f'sealed_truth/warm_test/{case}.npz',allow_pickle=False) as z:truth=z['m'].reshape(24,24).copy()
            nn=by[case,'physical_lbfgs','nn'];d=REPO/tasks[nn['task_id']]['output_dir']
            with np.load(d/'initial.npz',allow_pickle=False) as z:init=z['m'].reshape(24,24).copy()
            images=[truth,init];titles=['GT (offline)',f"NN initialization\nNRMSE={nn['initial_material_nrmse']:.4f}"]
            for arm,initial in arms:
                rr=by[case,arm,initial];p=REPO/tasks[rr['task_id']]['output_dir']
                with np.load(p/'prediction.npz',allow_pickle=False) as z:images.append(z['m'].reshape(24,24).copy())
                titles.append(f"{arm} / {initial}\nNRMSE={rr['material_nrmse']:.4f}")
            vmin=min(.02,*(float(x.min()) for x in images));vmax=max(float(x.max()) for x in images)
            fig,axes=plt.subplots(2,3,figsize=(12,8.5))
            for ax,m,title in zip(axes.ravel(),images,titles):
                im=ax.imshow(m,vmin=vmin,vmax=vmax,cmap='viridis');ax.axis('off');ax.set_title(title,fontsize=10)
            fig.suptitle(f"{nn['cohort']} / {case}\nFrozen initializer and final paid-Phi candidates; full shared contrast range",fontsize=12)
            fig.subplots_adjust(top=.83,left=.025,right=.91,bottom=.04,hspace=.27,wspace=.08)
            cax=fig.add_axes([.94,.13,.014,.60]);fig.colorbar(im,cax=cax,label='Material contrast')
            pdf.savefig(fig)
            if case in fixed:fig.savefig(figdir/f'neural_initial_to_final_fixed_{case}.png',dpi=160)
            plt.close(fig)
    (out/'RECEIPT.json').write_text(json.dumps(dict(status='ALL42_OFFLINE_VISUALIZED',objects=42,
        fixed_PNG_cases=sorted(fixed),selection='Six indices fixed before scoring; complete all42 PDF; no GT selection',
        truth_access='AUTHORIZED_POST_BOTH_COMPLETE_SCORES',new_physical_calls=0,LU=0,complex_RHS=0,
        cpu_seconds=time.process_time()-CPU,wall_seconds=time.perf_counter()-WALL),indent=2)+'\n')

if __name__=='__main__':main()
