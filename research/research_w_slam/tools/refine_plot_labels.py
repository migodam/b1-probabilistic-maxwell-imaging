#!/usr/bin/env python3
"""Presentation-only repair of two saved pilot plots. No solver or experiments.

Original fixed parent IDs and descriptive two-noise posterior averages stay
unchanged. Existing CSVs and every other figure/document are left untouched.
Run through the accounting wrapper; array rendering still costs CPU.
"""
from pathlib import Path
from collections import defaultdict
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

ROOT=Path(__file__).resolve().parents[1]
PARENTS=(0,3,6,9)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--figures',type=Path)
    args=parser.parse_args();root=args.root.resolve();out=(args.figures or root/'figures/pilot_actual').resolve()
    targets=[out/(name+suffix) for name in ['E2_fixed_beliefs','E3_error_vs_CPU'] for suffix in ['.png','.pdf']]
    if not out.is_dir() or any(not path.exists() for path in targets):
        raise FileNotFoundError('This label refinement replaces only the four existing target image files')
    inference=root/'results/e2e3/inference';gridpath=root/'results/e2e3/finite_bank/online/FINITE_GRID.npz'
    with np.load(inference/'POSTERIORS.npz',allow_pickle=False) as z:
        if not bool(z['complete']) or not z['completed_mask'].all():raise RuntimeError('Complete saved reference required')
        refs=z['weights'].copy();methods=z['methods'].astype(str).tolist()
    with np.load(gridpath,allow_pickle=False) as z:support=z['support'].copy().reshape(-1,4);pose=z['pose'].copy()
    if refs.shape!=(12,2,3,2950) or support.shape!=(2950,4) or pose.shape!=(2950,2):
        raise ValueError('Unexpected finite reference/grid schema')
    logs=json.loads((inference/'METHOD_RESULTS.json').read_text())
    setup_path=root/'results/e2e3/finite_bank/online/JOB_COST_FINAL.json'
    setup=json.loads(setup_path.read_text())
    if setup.get('complete') is not True:raise RuntimeError('Complete terminal setup cost receipt required')
    setup_cpu=float(setup['measured_process_cpu_seconds'])
    if not np.isfinite(setup_cpu):raise ValueError('Nonfinite setup CPU')
    if len(methods)!=11 or len(logs)!=12*2*3*len(methods):raise ValueError('Expected full eleven-method panel')
    # Canonical standalone cached attribution: sum the three updates within each
    # noise case, average repeats within parent, then average twelve parents.
    by=defaultdict(list)
    for row in logs:by[(row['parent'],row['repeat'],row['method'])].append(row)
    points=[]
    for index,method in enumerate(methods,1):
        parent_cpu=[];parent_regret=[];parent_cdf=[]
        for parent in range(12):
            cpu=[];regret=[];cdf=[]
            for repeat in (0,1):
                rr=sorted(by[(parent,repeat,method)],key=lambda r:r['time'])
                if [r['time'] for r in rr]!=[1,2,3]:raise ValueError('Incomplete/duplicate method prefix rows')
                values=[float(r['standalone_cpu_seconds']) for r in rr]
                if not np.isfinite(values).all():raise ValueError('Canonical standalone CPU is missing/nonfinite')
                cpu.append(sum(values))
                final=rr[-1]['posterior_metrics']
                regret.append(float(final['action_regret']));cdf.append(float(max(final['coordinate_cdf_sup'][-2:])))
            parent_cpu.append(float(np.mean(cpu)));parent_regret.append(float(np.mean(regret)));parent_cdf.append(float(np.mean(cdf)))
        points.append(dict(ID=index,method=method,cached_CPU_ms=1000*float(np.mean(parent_cpu)),
            mean_action_regret=float(np.mean(parent_regret)),mean_target_CDF_sup=float(np.mean(parent_cdf)),
            parent_count=12,noise_cases=24))
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.bbox':'tight'})
    # E2 keeps the original statistic exactly: mean of two saved posterior
    # distributions, never inference conditioned on an averaged observation.
    norm=Normalize(0.,1.);xv=np.unique(pose[:,0]);yv=np.unique(pose[:,1])
    pose_ids=np.searchsorted(yv,pose[:,1])*len(xv)+np.searchsorted(xv,pose[:,0])
    e2values=[]
    fig,axes=plt.subplots(4,6,figsize=(16,12),squeeze=False)
    for ip,parent in enumerate(PARENTS):
        for prefix in range(3):
            w=refs[parent,:,prefix].mean(axis=0)
            if not np.isfinite(w).all() or not np.isclose(w.sum(),1.,rtol=0,atol=1e-10):
                raise ValueError('Saved descriptive posterior average is not normalized')
            occupied=w@support
            pmass=np.bincount(pose_ids,weights=w,minlength=len(xv)*len(yv)).reshape(len(yv),len(xv))
            a,b=axes[ip,2*prefix:2*prefix+2]
            a.imshow(occupied.reshape(2,2),origin='lower',extent=(-.5,1.5,-.5,1.5),norm=norm,cmap='viridis',interpolation='nearest')
            for row,col in np.ndindex((2,2)):
                value=occupied[2*row+col]
                a.text(col,row,f'z{2*row+col}\n{value:.3f}',ha='center',va='center',fontsize=7,color='white' if value<.5 else 'black')
            a.set_xticks([0,1]);a.set_yticks([0,1]);a.set_xlabel('Public region column');a.set_ylabel('Public region row')
            a.set_title(f'P{parent}, prefix{prefix+1}: occupancy',fontsize=9)
            xx,yy=np.meshgrid(xv,yv)
            b.scatter(xx.ravel(),yy.ravel(),c=pmass.ravel(),s=80,norm=norm,cmap='magma',edgecolors='#777777',linewidths=.3)
            b.set_xlim(-1.1,1.1);b.set_ylim(1.5,2.7);b.set_xticks([-1,0,1]);b.set_yticks([1.6,2.1,2.6])
            b.set_xlabel(r'$p_x/\lambda_0$');b.set_ylabel(r'$p_y/\lambda_0$')
            b.set_title(f'P{parent}, prefix{prefix+1}: pose atom mass',fontsize=9)
            e2values.append(dict(parent=parent,prefix=prefix+1,occupancy=occupied.tolist(),pose_atom_mass=pmass.tolist()))
    fig.subplots_adjust(left=.055,right=.88,bottom=.065,top=.925,wspace=.65,hspace=.65)
    ca=fig.add_axes([.915,.565,.012,.32]);cb=fig.add_axes([.915,.145,.012,.32])
    fig.colorbar(ScalarMappable(norm=norm,cmap='viridis'),cax=ca,ticks=[0,.25,.5,.75,1],label='Occupancy probability (common0–1)')
    fig.colorbar(ScalarMappable(norm=norm,cmap='magma'),cax=cb,ticks=[0,.25,.5,.75,1],label='Finite pose atom mass (common0–1)')
    fig.suptitle('Fixed parents0/3/6/9 · descriptive two-noise posterior average\nOccupancy axes are public2×2 region indices, not physical wall shape; pose colors are atom masses',fontsize=12,y=.985)
    for suffix in ['.png','.pdf']:fig.savefig(out/('E2_fixed_beliefs'+suffix),dpi=180)
    plt.close(fig)
    # E3 numbered IDs have dedicated, evenly spaced label lanes; no long method
    # names or overlapping labels are placed at near-zero data coordinates.
    fig=plt.figure(figsize=(15,9.5))
    gs=fig.add_gridspec(2,2,height_ratios=[1.,1.05],left=.065,right=.92,bottom=.035,top=.875,hspace=.42,wspace=.5)
    axes=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1])]
    palette=plt.get_cmap('tab20')
    for ax,key,title in zip(axes,['mean_action_regret','mean_target_CDF_sup'],['Final action regret under finite reference','Final target CDF sup error']):
        ordered=sorted(points,key=lambda r:(r[key],r['cached_CPU_ms'],r['ID']))
        lanes={row['ID']:.05+i*.09 for i,row in enumerate(ordered)}
        for row in points:
            color=palette((row['ID']-1)%20)
            ax.scatter(row['cached_CPU_ms'],row[key],s=34,color=color,zorder=3)
            ax.annotate(str(row['ID']),xy=(row['cached_CPU_ms'],row[key]),xycoords='data',
                xytext=(1.035,lanes[row['ID']]),textcoords='axes fraction',ha='left',va='center',fontsize=8,
                color=color,bbox=dict(boxstyle='round,pad=.16',facecolor='white',edgecolor=color,linewidth=.6),
                arrowprops=dict(arrowstyle='-',color=color,linewidth=.55,alpha=.75),annotation_clip=False)
        ax.set_title(title);ax.set_xlabel('Attributed cached CPU per3-prefix case (ms)');ax.set_ylabel(title)
        ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,3));ax.margins(x=.09,y=.12);ax.grid(alpha=.18)
    legend=fig.add_subplot(gs[1,:]);legend.axis('off')
    cells=[[str(r['ID']),r['method'],f"{r['cached_CPU_ms']:.4g}",f"{r['mean_action_regret']:.4g}",f"{r['mean_target_CDF_sup']:.4g}"] for r in points]
    table=legend.table(cellText=cells,colLabels=['ID','All11 methods (source order)','Cached CPU (ms)','Mean regret','Mean target CDF sup'],
                       colWidths=[.055,.445,.16,.16,.18],cellLoc='left',bbox=[0,0,1,1])
    table.auto_set_font_size(False);table.set_fontsize(8)
    for i,row in enumerate(points,1):
        color=palette((row['ID']-1)%20);table[i,0].get_text().set_color(color);table[i,0].get_text().set_weight('bold')
    fig.suptitle(f'Shared physical bank:11,800 RHS · 2,950 FEM factors · {setup_cpu:,.1f} CPU s (excluded from cached x-axis)\nCanonical standalone attribution includes likelihood/reference dependence; CPU is not latency and no online acceleration is claimed',fontsize=11,y=.98)
    for suffix in ['.png','.pdf']:fig.savefig(out/('E3_error_vs_CPU'+suffix),dpi=180)
    plt.close(fig)
    note=dict(scope='PRESENTATION_ONLY_TWO_EXISTING_PLOTS',changed_files=[p.name for p in targets],
        existing_CSV_policy='All existing CSV files preserved; none opened for writing',other_figures_or_documents_modified=False,
        source_arrays_changed=False,new_PDE=False,new_rhs=0,new_factorizations=0,new_gate=False,
        E2=dict(fixed_parent_ids=list(PARENTS),prefixes=[1,2,3],statistic='refs[parent,:,prefix].mean(axis=0), same as original generator',
            interpretation='descriptive average of two posterior distributions; not posterior conditioned on averaged data',
            occupancy_coordinates='2x2 public region row/column indices, not physical shape',
            pose_coordinates='p_x/lambda0,p_y/lambda0',occupancy_color_limits=[0,1],pose_atom_mass_color_limits=[0,1],
            plotted_saved_values=e2values),
        E3=dict(source='METHOD_RESULTS.standalone_cpu_seconds',statistic='sum3updates pernoise case, average2repeats perparent, average12parents',
            x_units='milliseconds of attributed cached CPU',method_ID_table=points,
            common_setup_rhs=11800,common_setup_FEM_factors=2950,common_setup_CPU_seconds=setup_cpu,
            common_setup_CPU_source=str(setup_path.relative_to(root)),setup_CPU_excluded_from_cached_axis=True,
            latency_or_speed_advantage_claim=False),
        provenance=['results/e2e3/inference/POSTERIORS.npz','results/e2e3/finite_bank/online/FINITE_GRID.npz',
                    'results/e2e3/inference/METHOD_RESULTS.json','results/e2e3/finite_bank/online/JOB_COST_FINAL.json'])
    (out/'LABEL_REFINEMENT.json').write_text(json.dumps(note,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(figures_directory=str(out),changed=[p.name for p in targets],CSV_modified=False,new_PDE=False)))


if __name__=='__main__':main()
