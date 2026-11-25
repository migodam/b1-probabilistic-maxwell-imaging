"""Branch-and-bound the entire raster with verified physical parent-box tests."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np
from threadpoolctl import threadpool_limits

OUT=Path(__file__).resolve().parent
MAPS=OUT.parents[1]/'maps'
UPDATES=OUT.parent/'updated_witnesses'
sys.path.insert(0,str(UPDATES))
sys.path.insert(0,str(MAPS))
from common import FREEZE, Geometry, clock, elapsed, geom_path, pattern, save_json
from sector_map import SectorGeometry, local_certificate, margin_distances, response_reduction
from recheck import recheck_using_inverse

PANEL={32:('M32-D01','M32-R01','M32-P01','M32-A01'),
       96:('M96-R01','M96-P01')}
BACKGROUND=FREEZE['model']['background_upper']
LOW,HIGH=FREEZE['material_classes']['primary_core_range']
WEIGHTS={32:(.125,.5,1.,2.,8.),96:(.125,.03125,.5)}
SCALES=(.5,1.,.25,2.,.125,4.,.0625,.03125,.015625)
GROUP_SHAPES32={(8,8),(4,4),(2,2),(2,1),(1,2)}
GROUP_SHAPES96={(6,6),(3,3),(2,2),(2,1),(1,2)}


def read(path):
    return json.loads(path.read_text())


def tile_centres(tile,side):
    x,y,w,h=tile
    return (np.arange(y,y+h)[:,None]*side+np.arange(x,x+w)[None,:]).ravel()


def split(tile):
    x,y,w,h=tile
    if w==h==1:
        return []
    if w==2 and h==2:
        return [(x,y,2,1),(x,y+1,2,1)]
    widths=(w,) if w==1 else (w//2,w-w//2)
    heights=(h,) if h==1 else (h//2,h-h//2)
    return [(int(xx),int(yy),ww,hh)
            for xx,ww in zip(x+np.cumsum((0,*widths[:-1])),widths,strict=True)
            for yy,hh in zip(y+np.cumsum((0,*heights[:-1])),heights,strict=True)]


def child_box_hull(offsets,profile,w,h):
    """Coordinatewise envelope of the union, including exterior [0,b] for each child."""
    slots={}
    offsets_and_profile=[(int(a),int(b),float(p)) for (a,b),p in zip(offsets,profile,strict=True)]
    for sy in range(h):
        for sx in range(w):
            for dx,dy,p in offsets_and_profile:
                key=(sx+dx,sy+dy)
                if key in slots:
                    count,low,high=slots[key]
                    slots[key]=(count+1,min(low,LOW*p),max(high,HIGH*p))
                else:
                    slots[key]=(1,LOW*p,HIGH*p)
    keys=sorted(slots,key=lambda xy:(xy[1],xy[0]))
    bounds=np.array([(low if count==w*h else 0., high if count==w*h else max(BACKGROUND,high))
                     for count,low,high in (slots[x] for x in keys)])
    return np.array(keys,dtype=int),bounds[:,0],bounds[:,1]


def inset_group(tiles,g,offset):
    coordinates=np.array([[x,y] for x,y,_,_ in tiles])[:,None,:]+offset[None,:,:]
    linear=coordinates[...,1]*g.side+coordinates[...,0]
    if np.any(coordinates<0) or np.any(coordinates>=g.side):
        raise ArithmeticError('Clipped parent was put in an interior translation class.')
    return linear


def process_level(eligible,g,sector,footprint,Y,r,library,costs):
    answer={}
    for key,tiles in eligible.items():
        w,h=key
        if key not in library:
            offset,low,high=child_box_hull(*footprint,w,h)
            canonical=inset_group(tiles[:1],g,offset)[0]
            proofs=[]
            for weight in WEIGHTS[g.side]:
                t=clock()
                verified=local_certificate(sector,canonical,np.ones(len(canonical)),low,high,weight)
                dt=elapsed(t)
                costs['parent_local_wall_s']+=dt['wall_s']
                costs['parent_local_cpu_s']+=dt['cpu_s']
                proofs.append((weight,verified))
            library[key]=(offset,low,high,proofs)
        offset,_,_,proofs=library[key]
        cells=inset_group(tiles,g,offset)
        best=np.full(len(tiles),-np.inf)
        index=np.full(len(tiles),-1)
        selected_scale=np.zeros(len(tiles))
        selected_weight=np.zeros(len(tiles))
        for weight,verified in proofs:
            pending=np.flatnonzero(best<=0)
            if not len(pending):
                break
            if verified is None:
                continue
            for q in np.array_split(pending,max(1,(len(pending)+47)//48)):
                t=clock()
                batch=response_reduction(sector,cells[q],verified)
                margin,k,scale=margin_distances(batch,Y,r,SCALES)
                dt=elapsed(t)
                costs['parent_online_wall_s']+=dt['wall_s']
                costs['parent_online_cpu_s']+=dt['cpu_s']
                costs['parent_node_weight_evaluations']+=len(q)
                improved=margin>best[q]
                row=q[improved]
                best[row]=margin[improved]
                index[row]=k[improved]
                selected_scale[row]=scale[improved]
                selected_weight[row]=weight
        for j,tile in enumerate(tiles):
            answer[tile]=(best[j],int(index[j]),selected_scale[j],selected_weight[j],
                          any(v is not None for _,v in proofs))
    return answer


def run_object(g,sector,oid,metadata,*,leaf_replay=None,destination=OUT):
    started=clock()
    shape=metadata['shape']
    with np.load(MAPS/'public'/f'{oid}.npz',allow_pickle=False) as p:
        Y,r=p['Y'],float(p['r'])
    offsets,profile=pattern(g.side,shape)
    lower,upper=offsets.min(axis=0),offsets.max(axis=0)
    sizes=GROUP_SHAPES32 if g.side==32 else GROUP_SHAPES96
    library={}
    costs=defaultdict(float)
    wave=[(0,0,g.side,g.side)]
    depth=0
    tested=0
    pruned=np.zeros(g.n,dtype=bool)
    fallback=np.zeros(g.n,dtype=bool)
    records=[]
    stages=defaultdict(lambda: defaultdict(int))
    for t in sizes:
        stages[t]
    while wave:
        later=[]
        eligible=defaultdict(list)
        for node in wave:
            x,y,w,h=node
            if w==h==1:
                fallback[x+y*g.side]=True
                records.append((*node,depth,0,-np.inf,-1,0.,0.))
                continue
            interior=(np.all(np.r_[x,y]+lower>=0) and
                      np.all(np.r_[x+w-1,y+h-1]+upper<g.side))
            if (w,h) in sizes and interior:
                eligible[w,h].append(node)
                continue
            if (w,h) in sizes and not interior:
                stages[w,h]['edge_subdivided']+=1
            records.append((*node,depth,1 if not interior else 2,-np.inf,-1,0.,0.))
            later.extend(split(node))
        computed=process_level(eligible,g,sector,(offsets,profile),Y,r,library,costs)
        tested+=len(computed)
        for node,(margin,negative,scale,weight,has_state) in computed.items():
            w,h=node[2:]
            row=stages[w,h]
            row['interior_attempted']+=1
            row['verified_state_count']+=int(has_state)
            row['finite_margin_count']+=int(np.isfinite(margin))
            if margin>0:
                ids=tile_centres(node,g.side)
                if np.any(pruned[ids]):
                    raise ArithmeticError('A centre was already covered by a pruned parent.')
                pruned[ids]=True
                row['strict_parents']+=1
                row['pruned_centres']+=len(ids)
                row['strict_negative_inertia_parents']+=int(negative>0)
            else:
                later.extend(split(node))
            records.append((*node,depth,3 if margin>0 else 4 if has_state else 5,margin,negative,scale,weight))
        if eligible:
            print(oid.split('__')[1], 'depth',depth,'groups',len(computed),'branch-pruned',int(pruned.sum()),flush=True)
        wave=later
        depth+=1
    if not np.all(pruned^fallback):
        raise ArithmeticError('Recursive parent tests and singleton queue did not partition the field.')

    parent_time=elapsed(started)
    fresh=None
    if leaf_replay is not None:
        fresh=leaf_replay(g,sector,shape,fallback,Y,r)
    tick=clock()
    archive_path=MAPS/'results'/f'{oid}__b{BACKGROUND:g}__primary_sector.npz'
    cache=np.load(archive_path,allow_pickle=False)
    individual=cache['margin_sector']>0
    closed=pruned.copy()
    if fresh is None:
        closed[fallback]=individual[fallback]
    else:
        closed[fallback]=fresh['arrays']['adaptive_margin'][fallback]>0
    lookup_time=elapsed(tick)
    one=read(MAPS/'results'/f'witnesses_{g.side}_{shape}.json')['witnesses'][oid][str(BACKGROUND)]['witness']
    with np.load(MAPS/'results'/one['file'],allow_pickle=False) as obj:
        cx,cy=map(int,obj['centre'])
        if r!=float(one['noise_radius']):
            raise ArithmeticError('Old witness radius changed.')
        audit=recheck_using_inverse(sector,obj['material'],g.B,obj['current'],Y,r,one['full_field_inverse_upper'])
    if not audit['verified'] or closed[cy*g.side+cx]:
        raise ArithmeticError('Adaptive tree does not retain verified original physical witness.')
    chronology=read(MAPS/'results'/f'sector_{g.side}_{shape}_primary.json')
    arrays=dict(records=np.array(records),parent_pruned=pruned,leaf_fallback=fallback,
                excluded=closed,archived_excluded=individual)
    for (w,h),(offset,l,u,proofs) in library.items():
        arrays[f'tile_{w}x{h}__offset']=offset
        arrays[f'tile_{w}x{h}__lower']=l
        arrays[f'tile_{w}x{h}__upper']=u
        arrays[f'tile_{w}x{h}__state_available']=np.array([q is not None for _,q in proofs])
    if fresh is not None:
        arrays.update(fresh['arrays'])
    np.savez_compressed(destination/f'{oid}.npz',**arrays)
    receipt=dict(oid=oid,shape=shape,side=g.side,total_centres=g.n,
                 tested_parent_nodes=tested,strict_parent_certificates=sum(x['strict_parents'] for x in stages.values()),
                 pruned_by_positive_parent=int(np.sum(pruned)),
                 fallback_singletons=int(np.sum(fallback)),
                 archived_exhaustive_exclusions=int(np.sum(individual)),
                 final_exclusions=int(np.sum(closed)),
                 new_exclusions_beyond_archived_singletons=int(np.sum(closed&~individual)),
                 archive_exclusions_lost=int(np.sum(individual&~closed)),
                 final_survivors=int(np.sum(~closed)),
                 fixed_translation_types=len(library),
                 singleton_archive=str(archive_path.relative_to(OUT.parents[1])),
                 archived_singleton_reference_read=lookup_time,
                 uses_archived_singleton_for_decisions=fresh is None,
                 parent_stage_total=parent_time,
                 archived_exhaustive_reference=dict(wall_s=chronology['elapsed']['wall_s'],
                     cpu_s=chronology['elapsed']['cpu_s'],observation_count=len(chronology['summary']),
                     costs_are_batched_and_not_comparable= True), counts_by_tile={
                     f'{w}x{h}':dict(v) for (w,h),v in sorted(stages.items(),reverse=True)},
                 new_parent_costs=dict(costs),
                 old_physical_witness=dict(file=one['file'],centre=[cx,cy],retained=True,strict_check=audit),
                 **elapsed(started))
    if fresh is not None:
        exhaustive=fresh['arrays']['exhaustive_margin']>0
        same=np.array_equal(exhaustive,closed)
        receipt['fresh_singleton_replay']=dict(**fresh['receipt'],
            fresh_exhaustive_exclusions=int(exhaustive.sum()),
            same_as_fresh_exhaustive=bool(same),
            adaptive_additions_vs_fresh_exhaustive=int((closed&~exhaustive).sum()),
            adaptive_losses_vs_fresh_exhaustive=int((exhaustive&~closed).sum()),
            fresh_exhaustive_parity_with_archive=bool(np.array_equal(exhaustive,individual)))
    print(oid.split('__')[1],receipt['strict_parent_certificates'],receipt['pruned_by_positive_parent'],
          'of',g.n,'archive/new excluded',receipt['archived_exhaustive_exclusions'],receipt['final_exclusions'],flush=True)
    save_json(destination/f'{oid}.json',receipt)
    return receipt


def main():
    a=argparse.ArgumentParser()
    a.add_argument('side',type=int,choices=(32,96))
    v=a.parse_args()
    start=clock()
    g=Geometry(geom_path(v.side))
    s=SectorGeometry(g,BACKGROUND)
    setup=elapsed(start)
    roster=read(MAPS/'public'/'acquisition.json')['observations']
    rows={}
    for oid,info in roster.items():
        if info['side']==v.side and info['scene']=='weak' and oid.split('__')[1] in PANEL[v.side]:
            rows[oid]=run_object(g,s,oid,info)
    if len(rows)!=len(PANEL[v.side]):
        raise ArithmeticError('Archived full-field observation panel incomplete.')
    save_json(OUT/f'panel_{v.side}.json',dict(side=v.side,global_setup=setup,rows=rows,
                  weights=WEIGHTS[v.side],scales=SCALES,prior_background=BACKGROUND,
                  setup_and_total=elapsed(start),
                  proof='Convex coordinatewise parent material box includes all independently bounded child classes; use inherited one-sided signed sector tests before inspecting cached strict singleton proofs. Whole boundary coverage with clipped groups split, eligible interior translates share the same one-sided local certificate.',
                  record_columns=('x','y','width','height','depth','status','margin','negative_index','scale','weight'),
                  status=dict(leaf=0,edge_split=1,large_or_unselected_split=2,rejected_parent=3,
                              unresolved_parent=4,state_verification_unavailable=5)))


if __name__=='__main__':
    with threadpool_limits(limits=1,user_api='blas'):main()
