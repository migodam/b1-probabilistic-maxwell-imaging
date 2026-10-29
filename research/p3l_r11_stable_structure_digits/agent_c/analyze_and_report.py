"""Zero-solve postprocessing: transparent figures and Chinese bound report."""
from __future__ import annotations
import csv,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

D=Path(__file__).resolve().parent;R=D/'path_results'
sys.path.insert(0,str(D.parent/'shared'));import core


def read(name):return list(csv.DictReader((R/name).open()))
def arr(rows,key):return np.array([float(r[key]) for r in rows])


def main():
    rows=read('NODE_RESULTS.csv');attempts=read('ATTEMPT_RESULTS.csv')
    summary=json.loads((R/'SUMMARY.json').read_text());receipt=json.loads((R/'RECEIPT.json').read_text())
    figs=D/'figures';figs.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'savefig.dpi':180})
    colors={'actual':'#222222','global':'#c16b2b','local':'#25799b'}
    ranks=[32,64,96,128];rank_rows=[]
    for rank in ranks:
        for mode in ('global','local'):
            aa=[r for r in attempts if r['mode']==mode and int(r['rank'])==rank]
            rank_rows.append(dict(mode=mode,rank=rank,n=len(aa),
                actual_min=float(min(arr(aa,'actual_white_error'))),
                actual_median=float(np.median(arr(aa,'actual_white_error'))),
                actual_max=float(max(arr(aa,'actual_white_error'))),
                white_bound_median=float(np.median(arr(aa,'white_bound'))),
                potential_bound_median=float(np.median(arr(aa,'potential_bound'))),
                white_actual_pass=sum(float(r['actual_white_error'])<=.01 for r in aa),
                gate_pass=sum(r['accepted']=='True' for r in aa)))
    with (D/'RANK_SUMMARY.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rank_rows[0]));w.writeheader();w.writerows(rank_rows)
    fig,axs=plt.subplots(1,3,figsize=(15,4.6),layout='constrained')
    ax=axs[0]
    for mode,key,label in [('global','actual_white_error','Actual port error'),('global','white_bound','Global bound'),('local','white_bound','Local bound')]:
        median=[];lo=[];hi=[]
        for rank in ranks:
            q=np.quantile(arr([r for r in attempts if r['mode']==mode and int(r['rank'])==rank],key),[.25,.5,.75])
            lo.append(q[0]);median.append(q[1]);hi.append(q[2])
        c=colors['actual'] if key=='actual_white_error' else colors[mode]
        ax.plot(ranks,median,'o-',label=label,color=c);ax.fill_between(ranks,lo,hi,color=c,alpha=.12)
    ax.axhline(.01,color='#8b2635',ls='--',label='White-error gate = 0.01')
    ax.set(yscale='log',xlabel='Fixed Green-space rank',ylabel='Whitened response error / bound',
        title='A  Fixed propagation space remains inaccurate',xticks=ranks)
    ax.legend(fontsize=8,loc='lower left');ax.grid(axis='y',alpha=.15)
    ax=axs[1]
    aa=[r for r in attempts if r['mode']=='global' and int(r['rank'])==128]
    for mode in ('global','local'):
        bb=[r for r in attempts if r['mode']==mode and int(r['rank'])==128]
        ax.scatter(arr(bb,'actual_white_error'),arr(bb,'white_bound'),s=17,alpha=.55,c=colors[mode],label=f'{mode.capitalize()} bound')
    lim=[.05,1e3];ax.plot(lim,lim,'--',c='#888888',lw=1,label='Bound = actual')
    ax.set(xscale='log',yscale='log',xlabel='Actual rank-128 white error',ylabel='Bound on the same error',
        title='B  Conservative, but not the only failure')
    ax.legend(fontsize=8);ax.grid(alpha=.15)
    ax=axs[2]
    shape_points=[]
    for pi in range(8):
        parent=f'R10BI_FRESH_{pi:03d}'
        for ti in range(3):
            ss=[r for r in aa if r['parent']==parent and int(r['t_index'])==ti]
            shape_points.append(dict(parent=parent,t_index=ti,new_support_cells=int(ss[0]['new_support_cells']),
                mean_rank128_actual_white_error=float(np.mean(arr(ss,'actual_white_error')))))
    for pi in range(8):
        ss=[r for r in shape_points if r['parent']==f'R10BI_FRESH_{pi:03d}']
        ax.scatter([r['new_support_cells'] for r in ss],[r['mean_rank128_actual_white_error'] for r in ss],
            s=35,label=f'Parent {pi}')
    ax.axhline(.01,color='#8b2635',ls='--');ax.set(yscale='log',xlabel='Newly active ROI cells relative to anchor',
        ylabel='Shape mean actual rank-128 white error',title='C  Failures also occur without support birth')
    ax.text(.04,.05,'24 shapes, grouped within 8 exposed parents',transform=ax.transAxes,fontsize=8)
    ax.grid(axis='y',alpha=.15)
    fig.suptitle('R11 global/local dissipative control on the same R10 paths\n192 amplitude nodes per arm; full residual gates use the actual observed data',fontsize=12)
    for ext in ('png','pdf'):fig.savefig(figs/f'GLOBAL_LOCAL_PATHS.{ext}')
    plt.close(fig)
    costs=[]
    for mode in ('global','local'):
        rr=[r for r in rows if r['mode']==mode];aa=[r for r in attempts if r['mode']==mode]
        costs.append(dict(mode=mode,basis_setup_seconds=summary['modes'][mode]['basis_setup_seconds'],
            active_bound_fixed_setup_seconds=summary['modes'][mode]['bound']['setup_seconds'],
            projection_seconds=float(sum(arr(rr,'projection_seconds'))),
            approximate_solve_seconds=float(sum(arr(aa,'solve_seconds'))),
            bound_seconds=float(sum(arr(aa,'bound_seconds'))),
            exact_fallback_seconds=float(sum(arr(rr,'exact_seconds'))),
            node_wall_seconds=float(sum(arr(rr,'wall_seconds'))),
            conservative_RHS=int(sum(arr(rr,'conservative_RHS')))))
    fig,axs=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
    labels=['Fixed setup','Shape projection','Reduced solves','Residual gates','Exact fallback','Other node work']
    cols=['#577590','#43aa8b','#90be6d','#f9c74f','#f9844a','#b7b7b7']
    for i,c in enumerate(costs):
        setup=c['basis_setup_seconds']+c['active_bound_fixed_setup_seconds']
        known=c['projection_seconds']+c['approximate_solve_seconds']+c['bound_seconds']+c['exact_fallback_seconds']
        vals=[setup,c['projection_seconds'],c['approximate_solve_seconds'],c['bound_seconds'],c['exact_fallback_seconds'],max(0,c['node_wall_seconds']-known)]
        bottom=0.
        for label,col,v in zip(labels,cols,vals):
            axs[0].bar(i,v,bottom=bottom,color=col,label=label if i==0 else None);bottom+=v
        axs[1].bar(i,c['conservative_RHS'],color=colors[c['mode']])
    axs[0].set(xticks=[0,1],xticklabels=['Global','Local'],ylabel='Measured solver component time (s)',title='Cold solver components, including all fallbacks')
    axs[0].legend(fontsize=8,loc='upper left',bbox_to_anchor=(0,1))
    axs[1].set(xticks=[0,1],xticklabels=['Global','Local'],ylabel='New conservative complex RHS',title='Both arms: 0 accepted / 192 nodes')
    axs[1].bar_label(axs[1].containers[0],padding=3);axs[1].bar_label(axs[1].containers[1],padding=3)
    axs[1].set_ylim(0,9000)
    fig.suptitle('Timing is descriptive: one controlled window, without an isolated statistical speedup study',fontsize=11)
    for ext in ('png','pdf'):fig.savefig(figs/f'GLOBAL_LOCAL_COSTS.{ext}')
    plt.close(fig)
    core.write(D/'ANALYSIS.json',dict(status='COMPLETE_ZERO_SOLVE_POSTPROCESSING',rank_summary=rank_rows,
        shapes=shape_points,costs=costs,total_path_wall_seconds=summary['wall_seconds'],
        total_solver_components_seconds=sum(summary['modes'][m]['full_cold_task_seconds'] for m in ('global','local')),
        remaining_observation_io_metric_audit_and_instrumentation_seconds=summary['wall_seconds']-sum(summary['modes'][m]['full_cold_task_seconds'] for m in ('global','local')),
        instrumentation='Captured existing approximate fields/residuals; extra unused bound construction belongs to study overhead, not a standalone solver',
        no_new_physics=True,GT_access=False,total_C_conservative_RHS=receipt['RHS_conservative_total']+80,
        API_verification_extra_RHS=80,RHS576=0))


if __name__=='__main__':main()
