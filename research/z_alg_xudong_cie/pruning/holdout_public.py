"""Post-freeze controls using only public B3 data and already recorded branch boxes."""

import argparse
import json
import sys
from time import sleep

from threadpoolctl import threadpool_limits

from ablation import compare
from branches import as_branch,central_safety_control
from fast_support import TileProjection
from native_box import HERE,NativeBox,clock,elapsed

sys.path.insert(0,str(HERE.parent))
from holdout.data import (APERTURES,DEFAULT_ROOT,get_native_case,get_native_model,manifest_path,
                          planned_cases,write_json)


def independent_screen(row,nsource):
    profiles=row.get('screen',[])
    return bool(profiles) and all(c.get('status')=='compact_rank_budget' or
                                 c.get('nu',-1)>=nsource for c in profiles)


def broad(root,geometries):
    branch=central_safety_control()
    cache={}
    rows={}
    t=clock()
    for meta in planned_cases():
        case,public=get_native_case(meta['case'],root)
        geom=geometries[case['array']]
        rad=float(public['delta'])
        start=clock()
        if case['array'] in cache:
            result=dict(branch=branch.name,decision='RETAIN',residual_lower=0.,radius=rad,
                        reuse_reason='every frozen profile exceeded the numerical negative-rank budget or compact rank budget',
                        reused_from=cache[case['array']],**elapsed(start))
        else:
            result=geom.score(public['Y'],rad,branch,max_proofs=1,max_rank=850)
            if independent_screen(result,geom.ns):
                cache[case['array']]=case['case']
        result['tile_caps']=branch.caps
        rows[case['case']]=result
    return dict(scope='Central broad geometry declared before any truth access. There are no offline coverage labels in this output.',
                cases=rows,**elapsed(t))


def ablate(root,predictions,geometries,*,follow=False,checkpoint=None):
    apps={label:TileProjection(g) for label,g in geometries.items()}
    rows={}
    pending={row['case'] for row in planned_cases()}
    start=clock()
    report=dict(scope='Three separately computed port restrictions on the exact same previously nominated online boxes.')
    while pending:
        completed=[]
        for identity in sorted(pending):
            receipt=predictions/(identity+'__pruned.json')
            if not receipt.is_file():
                continue
            saved=json.loads(receipt.read_text())
            entries=[c for c in saved['certificates'] if c.get('branch') in ('ranked_support','shifted_support')]
            if len(entries)!=2 or any('tile_caps' not in e for e in entries):
                continue
            case,public=get_native_case(identity,root)
            app=apps[case['array']]
            rows[identity]=[compare(app,public['Y'],float(public['delta']),as_branch(c['branch'],c['tile_caps']))
                            for c in entries]
            completed.append(identity)
        pending.difference_update(completed)
        report.update(missing_online_boxes=sorted(pending),cases=rows,**elapsed(start))
        if checkpoint and completed:
            checkpoint(report)
            print('public ablation',len(rows),'complete;',len(pending),'pending',flush=True)
        if not follow or not pending:
            break
        sleep(10.)
    return report


def main(args):
    if not manifest_path(args.root).is_file():
        raise FileNotFoundError('Await the frozen public B3 observation directory.')
    start=clock()
    geometries={}
    for label in APERTURES:
        g=get_native_model(label,args.root)
        geometries[label]=NativeBox(g['spectrum'],g['S'],g['E'])
    geometry_setup=elapsed(start)
    payload={}
    if args.mode in ('broad','both'):
        payload['broad']=broad(args.root,geometries)
    if args.mode in ('ablation','both'):
        def checkpoint(report):
            write_json(args.output,dict(public_observations_only=True,geometry_setup=geometry_setup,
                                        **{**payload,'ablation':report},total_cost=elapsed(start)))
        payload['ablation']=ablate(args.root,args.predictions,geometries,follow=args.follow,checkpoint=checkpoint)
    write_json(args.output,dict(public_observations_only=True,geometry_setup=geometry_setup,
                               **payload,total_cost=elapsed(start)))
    print({k:(len(v['cases']),len(v.get('missing_online_boxes',[]))) for k,v in payload.items()},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=type(HERE),default=DEFAULT_ROOT)
    parser.add_argument('--predictions',type=type(HERE),default=HERE.parent/'benchmark/predictions')
    parser.add_argument('--output',type=type(HERE),default=HERE/'holdout_public_controls.json')
    parser.add_argument('--mode',choices=('broad','ablation','both'),default='both')
    parser.add_argument('--follow',action='store_true',help='Wait for the remaining public online boxes.')
    with threadpool_limits(1):
        main(parser.parse_args())
