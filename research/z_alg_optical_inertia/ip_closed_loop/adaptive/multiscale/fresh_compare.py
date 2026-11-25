"""Matched fresh multiscale and exhaustive evaluations of the same weak data."""
import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from centre_tree import (BACKGROUND, HIGH, LOW, MAPS, OUT, PANEL, SCALES, WEIGHTS,
                          Geometry, SectorGeometry, clock, elapsed, geom_path,
                          local_certificate, margin_distances, read, response_reduction,
                          run_object, save_json)
from certified_map import local_models

FRESH=OUT/'fresh'
ROUTES=('adaptive','exhaustive')


def add_cost(cost,category,tick):
    for unit,value in elapsed(tick).items():
        cost[f'{category}_{unit}']+=value


def fresh_singletons(g,sector,shape,needs_leaf,Y,r):
    tick=clock()
    n=g.n
    costs=defaultdict(float)
    phase=clock()
    # This is the same bank used by the archived same-contract primary sector.
    registry=read(MAPS/'results'/f'sector_{g.side}_{shape}_primary.json')
    weights=tuple(registry['weights'])
    scales=tuple(registry['scales'])
    geometries=sorted(local_models(g,shape),key=lambda t:len(t[0]),reverse=True)
    # Parent branches choose which fresh leaf proof to execute. The exhaustive
    # comparator sees every centre. Neither mask comes from the archive.
    masks=dict(adaptive=needs_leaf,exhaustive=np.ones(n,dtype=bool))
    margin={key:np.full(n,-np.inf) for key in ROUTES}
    negative={key:np.full(n,-1,dtype=np.int16) for key in ROUTES}
    cert_count={key:0 for key in ROUTES}
    cert_use={key:0 for key in ROUTES}
    visits={key:0 for key in ROUTES}
    built=valid=0
    adaptive_state_wall=adaptive_state_cpu=0.
    add_cost(costs,'clip_partition',phase)
    for class_id,(ids,cells,profile) in enumerate(geometries):
        for j,w in enumerate(weights):
            q={key:np.flatnonzero(masks[key][ids]&(margin[key][ids]<=0)) for key in ROUTES}
            if not len(q['exhaustive']) and not len(q['adaptive']):
                break
            setup=clock()
            # Keep one state block in memory. Both routes consume this exact
            # same checked factor before it is released.
            verified=local_certificate(sector,cells[0],profile,LOW,HIGH,w)
            dt=elapsed(setup)
            costs['one_common_state_bank_wall_s']+=dt['wall_s']
            costs['one_common_state_bank_cpu_s']+=dt['cpu_s']
            built+=1
            for key in ROUTES:
                cert_use[key]+=int(len(q[key])>0)
            if len(q['adaptive']):
                adaptive_state_wall+=dt['wall_s']
                adaptive_state_cpu+=dt['cpu_s']
            if verified is None:
                continue
            valid+=1
            # Rotate order so warm data cache order is paired by type or weight.
            order=ROUTES if (class_id+j)%2==0 else ROUTES[::-1]
            for key in order:
                if len(q[key]):
                    cert_count[key]+=1
                for k in range(0,len(q[key]),128):
                    start=clock()
                    selected=q[key][k:k+128]
                    reduced=response_reduction(sector,cells[selected],verified)
                    values,rank,_=margin_distances(reduced,Y,r,scales)
                    ii=ids[selected]
                    choose=values>margin[key][ii]
                    margin[key][ii[choose]]=values[choose]
                    negative[key][ii[choose]]=rank[choose]
                    visits[key]+=len(selected)
                    add_cost(costs,f'{key}_fresh_leaf',start)
        if g.side==96 and class_id%64==0:
            print('fresh',g.side,shape,'clip',class_id,'/',len(geometries),
                  'fresh exclusions',[(k,int((margin[k]>0).sum())) for k in ROUTES],
                  'state attempts',built,flush=True)
    return dict(arrays=dict(adaptive_margin=margin['adaptive'],exhaustive_margin=margin['exhaustive'],
                            adaptive_inertia=negative['adaptive'],exhaustive_inertia=negative['exhaustive']),
                receipt=dict(singleton_weights=weights,singleton_scales=scales,
                             clip_types=len(geometries),common_state_bank_attempts=built,
                             common_state_bank_verified=valid,
                             route_state_banks_needed=cert_use,route_verified_states_consumed=cert_count,
                             route_centre_weight_evaluations=visits,costs=dict(costs),
                             adaptive_alone_chargeable_state_wall_s=adaptive_state_wall,
                             adaptive_alone_chargeable_state_cpu_s=adaptive_state_cpu,
                             **elapsed(tick)))


def main():
    a=argparse.ArgumentParser()
    a.add_argument('side',type=int,choices=(32,96))
    a.add_argument('--only',choices=tuple(x for p in PANEL.values() for x in p))
    v=a.parse_args()
    if v.only and v.only not in PANEL[v.side]:
        a.error('Object side does not match the selected grid side.')
    started=clock()
    FRESH.mkdir(exist_ok=True)
    g=Geometry(geom_path(v.side))
    sector=SectorGeometry(g,BACKGROUND)
    global_setup=elapsed(started)
    roster=read(MAPS/'public'/'acquisition.json')['observations']
    rows={}
    objects=(v.only,) if v.only else PANEL[v.side]
    for oid,meta in roster.items():
        if meta['side']!=v.side or meta['scene']!='weak' or oid.split('__')[1] not in objects:
            continue
        result=run_object(g,sector,oid,meta,leaf_replay=fresh_singletons,destination=FRESH)
        rows[oid]=result
    if len(rows)!=len(objects):
        raise ArithmeticError('One or more requested data cases were not present.')
    tag=v.only or 'all'
    save_json(FRESH/f'panel_{v.side}_{tag}.json',dict(side=v.side,global_geometry_setup=global_setup,
              parent_weights=WEIGHTS[v.side],parent_scales=SCALES,rows=rows,**elapsed(started)))


if __name__=='__main__':
    with threadpool_limits(limits=1,user_api='blas'):main()
