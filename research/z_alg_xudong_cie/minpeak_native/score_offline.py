"""Join sealed primal certificates with the offline material and clean-port references."""

import argparse
import csv
import json

import numpy as np

from search import HERE, STUDY, clocks, elapsed, full_ring_inputs


def scores(estimate, truth, y_fit, y_clean):
    inside=truth>0
    half_estimate=estimate>=.5*estimate.max() if estimate.max()>0 else np.zeros(len(estimate),bool)
    half_truth=truth>=.5*truth.max()
    return dict(material_relative_l2=float(np.linalg.norm(estimate-truth)/np.linalg.norm(truth)),
        inside_relative_l2=float(np.linalg.norm(estimate[inside]-truth[inside])/np.linalg.norm(truth[inside])),
        background_rms=float(np.sqrt(np.mean(estimate[~inside]**2))),
        half_peak_iou=float(np.count_nonzero(half_estimate&half_truth)/np.count_nonzero(half_estimate|half_truth)),
        relative_clean_port_error=float(np.linalg.norm(y_fit-y_clean)/np.linalg.norm(y_clean)))


def read_json(path):
    return json.loads(path.read_text())


def render(rows, fields, output):
    import matplotlib.pyplot as plt

    objects=sorted(fields)
    fig,axes=plt.subplots(len(objects),3,figsize=(9,2.75*len(objects)),layout='constrained')
    fig.suptitle('Native 96² minimum-peak bracket | common contrast scale',fontsize=14)
    peak_scale=max(float(x.max()) for val in fields.values() for x in val.values())
    image=None
    for j,case in enumerate(objects):
        for k,arm in enumerate(('truth','1','2')):
            ax=axes[j,k]
            ax.set_xticks([])
            ax.set_yticks([])
            if arm not in fields[case]:
                ax.text(.5,.5,'No saved numerical candidate',ha='center',va='center')
                continue
            a=fields[case][arm]
            image=ax.imshow(a.reshape(96,96),vmin=0,vmax=peak_scale,cmap='viridis',interpolation='nearest')
            if arm=='truth':
                ax.set_title(f'{case} | offline truth {a.max():.3f}')
            else:
                r=next(r for r in rows if r['object']==case and r['noise_pct']==float(arm))
                status='verified' if r['upper_peak_online'] is not None else 'outside or unverified'
                ax.set_title(f'{arm}% | {status}\npeak {a.max():.3f}, residual {r["selected_numeric_residual"]:.3g}')
    fig.colorbar(image,ax=axes[:,2],shrink=.60,label='Contrast')
    fig.savefig(output,dpi=145)
    plt.close(fig)


def render_brackets(rows, output):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig,axes=plt.subplots(1,2,figsize=(10,4.3),sharex=True,sharey=True,layout='constrained')
    for noise,ax in zip((1.,2.),axes,strict=True):
        subset=[r for r in rows if r['noise_pct']==noise]
        for j,row in enumerate(subset):
            lo=row['lower_peak_all_inherited']
            up=row['best_known_peak_upper']
            ax.hlines(j,lo,up,color='#778592',linewidth=2,zorder=1)
            ax.scatter(lo,j,facecolors='white',edgecolors='#243b52',s=42,zorder=2)
            if row['upper_peak_online'] is not None:
                ax.scatter(row['upper_peak_online'],j,marker='o',color='#087b74',s=37,zorder=3)
            ax.scatter(row['upper_peak_offline_true_fallback'],j,marker='d',facecolors='none',
                       edgecolors='#bd6445',s=45,zorder=4)
        ax.set(title=f'{noise:g}% noise',xlabel='Maximum material contrast',xlim=(0,4.15),
               yticks=range(len(subset)),yticklabels=[r['object'] for r in subset])
        ax.grid(axis='x',alpha=.18)
    axes[0].invert_yaxis()
    marks=[Line2D([],[],marker='o',mfc='white',mec='#243b52',linestyle='none',label='Strict lower endpoint'),
           Line2D([],[],marker='o',color='#087b74',linestyle='none',label='Verified online candidate'),
           Line2D([],[],marker='d',mfc='none',mec='#bd6445',linestyle='none',label='Generating material (offline)')]
    fig.suptitle('Minimum achievable contrast peak on the original native data ball')
    fig.legend(handles=marks,loc='outside lower center',ncols=3,frameon=False)
    fig.savefig(output,dpi=165)
    plt.close(fig)


def main(args):
    start=clocks()
    _,_,_,observations=full_ring_inputs()
    lower=read_json(HERE/'lower_bounds.json')['observations']
    archive=STUDY/'synthetic/evaluation'
    offline=read_json(STUDY/'submission/curve/observation_audit.json')['observations']
    old_radius={r['id']:r for r in offline}
    by_id={}
    for path in args.verified:
        for case,rows in read_json(path)['cases'].items():
            for row in rows:
                by_id.setdefault(case,{})[row['source']]=row
    fields={}
    summaries=[]
    for obs in observations:
        if obs.noise_pct not in (1.,2.):
            continue
        object_data=fields.setdefault(obs.object_id,{})
        with np.load(archive/(obs.object_id+'.npz'),allow_pickle=False) as file:
            truth=file['m_fine']
            object_data['truth']=truth
        with np.load(archive/(obs.object_id+'__pw360_10x20__n02.npz'),allow_pickle=False) as file:
            y_clean=file['Y_fine_clean']
        certified=by_id.get(obs.id,{})
        candidates=[]
        for folder in args.inputs:
            path=folder/(obs.id+'.json')
            if path.exists():
                r=read_json(path)
                r['source']=folder.name
                r['file']=str(path.with_suffix('.npz').resolve().relative_to(HERE))
                candidates.append(r)
        passes=[(r,certified[r['source']]) for r in candidates if
                certified.get(r['source'],{}).get('status')=='CERTIFIED_FEASIBLE']
        selected=min(passes,key=lambda w:w[1]['peak']) if passes else None
        decision=selected[0] if selected else min(candidates,key=lambda r:r['selected_numeric_residual'],default=None)
        fallback=old_radius[obs.id]
        if fallback['true_distance_upper']>obs.padded_radius or fallback['admissibility_margin_lower']<=0:
            raise RuntimeError('Original offline true-material fallback does not have a verified ball inclusion.')
        true_peak=float(truth.max())
        metrics={}
        if decision is not None:
            with np.load(HERE/decision['file'],allow_pickle=False) as p:
                m,y=p['m_selected'],p['Y_selected_numeric']
            object_data[str(int(obs.noise_pct))]=m
            metrics=scores(m,truth,y,y_clean)
        bound=lower[obs.id]
        upper=selected[1]['peak'] if selected else None
        best_certified_cap=bound['strongest_including_original_paper1']['excluded_cap']
        if upper is not None and upper <= best_certified_cap:
            raise RuntimeError('Independent upper witness conflicts with the certified exclusion.')
        row=dict(observation=obs.id,object=obs.object_id,noise_pct=obs.noise_pct,radius=obs.padded_radius,
                 lower_peak_fixed_plus_gradient64=bound['strongest_fixed_plus_gradient64']['excluded_cap'],
                 lower_peak_all_inherited=best_certified_cap,upper_peak_online=upper,
                 upper_peak_offline_true_fallback=true_peak,
                 offline_true_port_residual_upper=fallback['true_distance_upper'],
                 best_known_peak_upper=min(true_peak,upper) if upper is not None else true_peak,
                 exhibited_peak_reduction_over_truth=max(0.,true_peak-upper) if upper is not None else None,
                 unknown_online_upper=upper is None,
                 online_branch=decision['source'] if decision else None,
                 online_material_file=decision['file'] if decision else None,
                 selected_numeric_peak=decision['selected_peak'] if decision else None,
                 selected_numeric_residual=decision['selected_numeric_residual'] if decision else None,
                 selected_strict_residual_upper=selected[1]['distance_upper'] if selected else None,
                 selected_positive_cells=decision['selected_positive_cells'] if decision else None,
                 **metrics)
        summaries.append(row)
    args.output.write_text(json.dumps({'scope':'Offline labels opened after saved online candidates and strict feasibility decisions. The truth fallback is a separate preverified feasible material, never an algorithm output.',
                                       'rows':summaries,**elapsed(start)},indent=2,allow_nan=False)+'\n')
    headers=list(dict.fromkeys(k for row in summaries for k in row))
    with args.output.with_suffix('.csv').open('w',newline='') as file:
        writer=csv.DictWriter(file,fieldnames=headers)
        writer.writeheader()
        writer.writerows(summaries)
    if args.figure is not None:
        render(summaries,fields,args.figure)
    if args.bracket_figure is not None:
        render_brackets(summaries,args.bracket_figure)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verified',nargs='+',type=type(HERE),required=True)
    parser.add_argument('--inputs',nargs='+',type=type(HERE),default=[HERE/'runs',HERE/'gn',HERE/'continuation'])
    parser.add_argument('--output',type=type(HERE),default=HERE/'brackets.json')
    parser.add_argument('--figure',type=type(HERE))
    parser.add_argument('--bracket-figure',type=type(HERE))
    main(parser.parse_args())
