"""Postscore rendering only: stable colors and distinct NOT_RUN placeholders.
The sealed scorer/physical code and all predictions remain unchanged.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import time
CPU,WALL=time.process_time(),time.perf_counter()
from pathlib import Path
import json
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'RESULTS/vector_extension';FIG=ROOT/'FIGURES/vector_extension'
rows=json.loads((OUT/'SCORES.json').read_text())['rows'];arms=json.loads((ROOT/'configs/vector_extension_phase2.json').read_text())['arms'];cases=list(dict.fromkeys(r['case_id'] for r in rows))
colors=dict(zip(arms,['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray']))
for view in ('central_slice','maximum_projection'):
    fig,axes=plt.subplots(len(cases),9,figsize=(20,9.2),squeeze=False,gridspec_kw=dict(wspace=.04,hspace=.05))
    for i,case in enumerate(cases):
        with np.load(ROOT/'sealed_truth/vector_extension'/case/'truth.npz',allow_pickle=False) as z:truth=z['chi'].reshape((8,)*3)
        for j,arm in enumerate(['truth']+arms):
            ax=axes[i,j];ax.set_xticks([]);ax.set_yticks([])
            if i==0:ax.set_title(arm.replace('_','\n'),fontsize=8)
            if j==0:ax.set_ylabel(case.replace('_',' '),fontsize=8)
            if arm=='truth':m=truth
            elif not any(r['case_id']==case and r['arm']==arm for r in rows):
                ax.set_facecolor('#eeeeee');ax.text(.5,.5,'NOT_RUN\n(control contract)',ha='center',va='center',fontsize=8);continue
            elif (OUT/case/arm/'prediction.npz').exists():
                with np.load(OUT/case/arm/'prediction.npz',allow_pickle=False) as z:m=z['m'].reshape((8,)*3)
            else:ax.text(.5,.5,'FAILED',ha='center',va='center');continue
            image=m[:,:,4] if view=='central_slice' else np.max(m,axis=2)
            ax.imshow(image.T,origin='lower',cmap='magma',vmin=0,vmax=8)
    fig.subplots_adjust(left=.035,right=.95,bottom=.02,top=.88)
    fig.suptitle('Strict real vector DDA: '+view.replace('_',' ')+'; chi=epsilon-1, all512unknowns\nScaled finite-model comparison; continuum accuracy NOT certified',y=.98)
    cax=fig.add_axes([.965,.2,.008,.6]);fig.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(0,8),cmap='magma'),cax=cax)
    fig.savefig(FIG/(view+'.png'),dpi=160);plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(13,9))
for ax,case in zip(axes.ravel(),cases):
    for row in [r for r in rows if r['case_id']==case]:
        path=OUT/'offline_curves'/f'{case}__{row["arm"]}.json'
        if not path.exists():continue
        c=json.loads(path.read_text());ax.plot([v['CPU_seconds'] for v in c],[v['chi_relative_L2_NRMSE'] for v in c],label=row['arm'],color=colors[row['arm']])
    ax.set_title(case+'\n'+next(r['platform'] for r in rows if r['case_id']==case));ax.set_xlabel('Online CPU s (setup/import/export separately accounted)');ax.set_ylabel('Chi relative L2 NRMSE');ax.grid(alpha=.2);ax.legend(fontsize=6)
fig.suptitle('Offline-scored paid original-Phi incumbents; no GT endpoint selection');fig.tight_layout();fig.savefig(FIG/'paid_cost_curves.png',dpi=160);plt.close(fig)
(OUT/'RENDER_COST.json').write_text(json.dumps(dict(CPU_seconds=time.process_time()-CPU,wall_seconds=time.perf_counter()-WALL,
    scope='Postscore visual clarity only, no physics/rescoring/algorithm change; sealed scorer source unchanged',changes=['Consistent per-arm colors across hardware/case panels','NOT_RUN distinct from FAILED','Tighter material image layout']),indent=2)+'\n')
