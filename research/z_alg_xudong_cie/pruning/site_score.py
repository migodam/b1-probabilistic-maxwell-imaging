"""Decode protected compact-site locations only after all sixteen online grids are sealed."""

import argparse
from collections import Counter
import json

import numpy as np

from native_box import HERE,clock,elapsed


def draw(rows,output):
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    objects=sorted({r['object_id'] for r in rows})
    views=(('PW360_10x20',.02),('PW360_10x20',.1),('PW180_10x20',.02),('PW180_10x20',.1))
    fig,axes=plt.subplots(4,4,figsize=(9.5,9.8),layout='constrained')
    for i,identity in enumerate(objects):
        for j,(ap,noise) in enumerate(views):
            ax=axes[i,j]
            row=next(r for r in rows if (r['object_id'],r['aperture'],r['noise_relative'])==(identity,ap,noise))
            shade=np.zeros((4,4))
            for y,x in row['survivors']:
                shade[y,x]=1
            ax.imshow(shade,vmin=0,vmax=1,cmap=ListedColormap(['#d9dfe5','#147a6b']),interpolation='nearest')
            y,x=row['true_tile']
            ax.add_patch(Rectangle((x-.43,y-.43),.86,.86,fill=False,edgecolor='#142537',linewidth=2.5))
            z,w=row['born_best_tile']
            ax.plot(w,z,'*',color='#f0ad32',markersize=10,markeredgecolor='#734809',markeredgewidth=.6)
            ax.set_xticks(np.arange(-.5,4,1.),minor=True)
            ax.set_yticks(np.arange(-.5,4,1.),minor=True)
            ax.tick_params(which='both',bottom=False,left=False,labelbottom=False,labelleft=False)
            ax.grid(which='minor',linewidth=.8,color='white')
            if i==0:
                ax.set_title(f'{ap[2:5]}° | {100*noise:g}%')
            if j==0:
                ax.set_ylabel(identity,fontsize=9)
    legends=[Line2D([],[],marker='s',color='#147a6b',linestyle='none',label='Retained tile'),
             Line2D([],[],marker='s',mfc='none',mec='#142537',linestyle='none',label='Generating tile (offline)'),
             Line2D([],[],marker='*',color='#f0ad32',linestyle='none',label='First-Born top rank')]
    fig.suptitle('All 16 native site hypotheses for four independent compact materials')
    fig.legend(handles=legends,loc='outside lower center',ncols=3,frameon=False)
    fig.savefig(output,dpi=155)
    plt.close(fig)


def classify(online,labels,truth_admissible):
    rows=[]
    for identity,pred in online['predictions'].items():
        oid=pred['object_id']
        truth=labels[oid]
        survivors=pred['survivors']
        receive=pred['receive_only_survivors']
        source=pred['illumination_only_survivors']
        born=pred['born_ranking']
        born_rank=next(i for i,row in enumerate(born,1) if row['tile']==truth)
        found=truth in survivors
        rows.append(dict(observation_id=identity,object_id=oid,aperture=pred['aperture'],
                         noise_relative=pred['noise_relative'],true_tile=truth,
                         radius=pred['radius'],truth_in_ball=truth_admissible,
                         survivors=survivors,true_retained=found,
                         exact_site_certified=found and len(survivors)==1,incorrect_pruned=16-len(survivors)-int(not found),
                         receive_only_survivors=receive,receive_only_true_retained=truth in receive,
                         receive_only_exact_site=receive==[truth],illumination_only_survivors=source,
                         illumination_only_true_retained=truth in source,illumination_only_exact_site=source==[truth],
                         born_best_tile=born[0]['tile'],born_true_rank=born_rank))
    return rows


def aggregate(rows,key,value):
    selected=[r for r in rows if key is None or r[key]==value]
    return dict(views=len(selected),independent_materials=len({r['object_id'] for r in selected}),
                truth_in_ball=sum(r['truth_in_ball'] for r in selected),
                true_retained=sum(r['true_retained'] for r in selected),
                false_prunes=sum(r['truth_in_ball'] and not r['true_retained'] for r in selected),
                exact_site_certified=sum(r['exact_site_certified'] for r in selected),
                wrong_tiles_strictly_pruned=sum(r['incorrect_pruned'] for r in selected),
                all_wrong_tiles=len(selected)*15,
                retained_counts=[len(r['survivors']) for r in selected],
                receive_only_exact_sites=sum(r['receive_only_exact_site'] for r in selected),
                receive_only_wrong_tiles_retained=sum(len(r['receive_only_survivors'])-r['receive_only_true_retained'] for r in selected),
                illumination_only_exact_sites=sum(r['illumination_only_exact_site'] for r in selected),
                illumination_only_wrong_tiles_retained=sum(len(r['illumination_only_survivors'])-r['illumination_only_true_retained'] for r in selected),
                born_top_one=sum(r['born_true_rank']==1 for r in selected),
                born_average_rank=float(np.mean([r['born_true_rank'] for r in selected])))


def run(args):
    start=clock()
    online=json.loads(args.online.read_text())
    if online['state']!='SEALED_PUBLIC_SITE_OUTPUTS' or len(online['predictions'])!=16:
        raise ValueError('Opening true sites requires a full online seal.')
    if any(len(pred['tile_scores'])!=16 for pred in online['predictions'].values()):
        raise ValueError('All 16 competing boxes are required for every observation.')
    labels=json.loads(args.labels.read_text())
    if len(labels)!=4 or any(len(site)!=2 or any(not isinstance(i,int) or not 0<=i<4 for i in site)
                             for site in labels.values()):
        raise ValueError('Expected only four anonymous object-to-tile labels.')
    acquisition=json.loads(args.acquisition.read_text())
    truth_admissible=(acquisition['public_observations']==acquisition['truth_inside_padded_ball']==16 and
                      acquisition['minimum_interval_margin']>0)
    if not truth_admissible:
        raise ValueError('Require independent admissibility of the full compact cohort.')
    rows=classify(online,labels,truth_admissible)
    primary=born=ablations=0.
    workflows=Counter()
    for result in online['predictions'].values():
        born+=result['born_comparator_cost']['cpu_seconds']
        for branch in result['tile_scores']:
            check=branch['primary']
            workflows[(check['rigor'],check['decision'])]+=1
            primary+=check['check_time']['cpu_seconds']
            ablations+=sum(c['cpu_seconds'] for c in branch['one_sided'])
    cost=dict(strict_all_tiles_cpu=primary,born_all_tiles_scoring_cpu=born,ablation_two_arms_cpu=ablations,
              geometry_and_born_setup=online['aperture_preparations'],online_with_artifact_io=online['online_cost_to_date'])
    results=dict(scope='Matched declared one-tile inclusion class with unchanged physical data and noise balls.',
                 online_sealed_at=online['sealed_at'],overall=aggregate(rows,None,None),
                 acquisition_receipt=args.acquisition.name,
                 by_aperture={k:aggregate(rows,'aperture',k) for k in sorted({r['aperture'] for r in rows})},
                 by_noise={v:aggregate(rows,'noise_relative',v) for v in (.02,.1)},
                 workflows=[dict(method=method,decision=decision,count=n) for (method,decision),n in workflows.items()],
                 costs=cost,case_rows=rows,**elapsed(start))
    args.output.write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
    if args.figure is not None:
        draw(rows,args.figure)
    print(results['overall'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--online',type=type(HERE),default=HERE/'site_online.json')
    parser.add_argument('--labels',type=type(HERE),required=True)
    parser.add_argument('--acquisition',type=type(HERE),default=HERE/'site_control/ACQUISITION_RECEIPT.json')
    parser.add_argument('--output',type=type(HERE),default=HERE/'site_scored.json')
    parser.add_argument('--figure',type=type(HERE))
    run(parser.parse_args())
