"""Bound new-measurement material diameters from old observations.

Products of the exact discrete VIE transfer numerators retain their cancellation
before material intervals are substituted. Only the acquisition/scoring programs
ever access the generating materials.
"""
import argparse
from dataclasses import dataclass
from functools import lru_cache
import heapq
import json
import math
from time import perf_counter, process_time

import numpy as np
from flint import acb, arb
from numba import njit

from data14 import ACTIONS, BANDS, ROOT, load_cases
from material_controls import GRID, HI, K, LO, RESULTS, SIGNS, quantize, restrict, wave_table

THREE = 3**K
POWER = 1 << 18
BASIS = np.sum(SIGNS * (3**np.arange(K))[None, :], axis=1)
WEIGHT_PRODUCT = np.array([2**sum((n//(3**j))%3 == 1 for j in range(K)) for n in range(THREE)])
TAUS = (.1, .5, 1.)


@dataclass
class RatioTable:
    data: np.ndarray
    exp: int
    degree: int


@lru_cache(maxsize=512)
def ratio_table(geometry, d1, d0):
    a = wave_table(geometry, d1)
    if d0 is None:
        q, num = zip(*a)
        degree = 1
    else:
        b = wave_table(geometry, d0)
        q = [acb(0) for _ in range(THREE)]
        num = [[acb(0) for _ in a[0][1]] for _ in range(THREE)]
        for i, (qa, na) in enumerate(a):
            for j, (qb, nb) in enumerate(b):
                h = BASIS[i] + BASIS[j]
                q[h] += qa * qb
                for t, (x, y) in enumerate(zip(na, nb)):
                    num[h][t] += x * qb - y * qa
        for h, weight in enumerate(WEIGHT_PRODUCT.tolist()):
            q[h] /= weight
            num[h] = [v/weight for v in num[h]]
        degree = 2
    den, eq = quantize([[v.real, v.imag] for v in q])
    ctrl, en = quantize([[v.real for v in row]+[v.imag for v in row] for row in num])
    return RatioTable(np.column_stack([den, ctrl]), en-eq, degree)


@njit(cache=True)
def mix(a, b, t):
    return a + (b-a)*t//GRID


@njit(cache=True)
def restrict_two(table, lo, hi):
    out = table.copy()
    stride = 1
    for axis in range(K):
        l, u = lo[axis], hi[axis]
        for v in range(THREE):
            if (v//stride) % 3:
                continue
            for j in range(table.shape[1]):
                a, b, c = out[v,j], out[v+stride,j], out[v+2*stride,j]
                al, bl = mix(a,b,l), mix(b,c,l)
                au, bu = mix(a,b,u), mix(b,c,u)
                out[v,j] = mix(al,bl,l)
                out[v+stride,j] = mix(al,bl,u)
                out[v+2*stride,j] = mix(au,bu,u)
        stride *= 3
    return out


@njit(cache=True)
def ratio_enclosure(table, exponent, degree, boxes):
    d = (table.shape[1]-2)//2
    centers = np.zeros((len(boxes),d),np.complex128)
    radii = np.full((len(boxes),d),np.inf)
    err = 1 + K*degree
    for i in range(len(boxes)):
        lo, hi = boxes[i,0], boxes[i,1]
        if degree==1:
            z = restrict(table,lo,hi)
        else:
            z = restrict_two(table,lo,hi)
        q = complex(z[:,0].mean(), z[:,1].mean())
        if abs(q)==0:
            continue
        # An integer projection furnishes a positive lower bound for |Q|.
        dx,dy = int(round(POWER*q.real/abs(q))),int(round(POWER*q.imag/abs(q)))
        projmin = np.min(dx*z[:,0]+dy*z[:,1]) - err*(abs(dx)+abs(dy))
        # An integer L1 norm is a verified upper bound for the direction's L2 norm.
        # A tighter integer norm comes from sqrt plus two integer units.
        dnorm = int(np.sqrt(float(dx*dx+dy*dy)))+2
        if projmin<=0:
            continue
        for t in range(d):
            a = complex(z[:,2+t].mean(),z[:,2+d+t].mean())/q
            if max(abs(a.real),abs(a.imag))>32:
                continue
            cr,ci = int(round(POWER*a.real)),int(round(POWER*a.imag))
            rr = POWER*z[:,2+t]-cr*z[:,0]+ci*z[:,1]
            ii = POWER*z[:,2+d+t]-ci*z[:,0]-cr*z[:,1]
            rounding = err*(POWER+abs(cr)+abs(ci))
            mr = max(abs(np.min(rr)),abs(np.max(rr))) + rounding
            mi = max(abs(np.min(ii)),abs(np.max(ii))) + rounding
            # Remainders lie in this rectangle. The slack exceeds binary64 sqrt,
            # multiply, and divide roundoff; the final radius is rounded upward.
            rad = math.sqrt(float(mr)**2+float(mi)**2)*dnorm/(float(projmin)*POWER)
            rad = np.nextafter((rad*(1.+2.e-14)+1.e-300)*2.**exponent,np.inf)
            centers[i,t] = complex(cr,ci)*2.**(exponent-18)
            radii[i,t] = rad
    return centers, radii


def group_cover(boxes, count):
    """Coarsen a union for pair comparisons; every original box stays inside a group."""
    groups=[]
    def add(ids):
        lo=boxes[ids,0].min(axis=0);hi=boxes[ids,1].max(axis=0)
        weight=np.max(hi-lo) if len(ids)>1 else 0
        heapq.heappush(groups,(-int(weight),int(ids[0]),ids,lo,hi))
    add(np.arange(len(boxes)))
    while len(groups)<count and groups and groups[0][0]<0:
        _,_,ids,lo,hi=heapq.heappop(groups)
        mid=boxes[ids,0]+boxes[ids,1]
        axis=np.argmax(np.ptp(mid,axis=0))
        order=ids[np.argsort(mid[:,axis],kind="stable")]
        cut=len(ids)//2
        add(order[:cut]);add(order[cut:])
    return np.array([[g[3],g[4]] for g in groups])


def next_down(x):
    return np.nextafter(x,-np.inf)


def disk_gap(center,radius,i,j):
    """Rigorously lowered per-channel gap; handles large, almost equal centers."""
    gap=np.abs(center[i]-center[j])
    uncertainty=radius[i]+radius[j]+4.e-14*(np.abs(center[i])+np.abs(center[j])+1.e-300)
    return np.maximum(0.,next_down(gap-uncertainty))


def physical_disk(value,gamma):
    center,rad=value
    return center, np.nextafter((1.+3e-14)*((1.+gamma)*rad+gamma*np.abs(center))+1.e-300,np.inf)


def score_case(case, counts=64, source="seeded", bands=BANDS):
    tic,cpu=perf_counter(),process_time()
    if source=="seeded":
        prefix=RESULTS/"material_covers_seeded"
        inner=np.load(prefix/(case.case_id+".npz"))["boxes"]
        oldwork=json.loads((prefix/(case.case_id+".json")).read_text())
    else:
        with np.load(RESULTS/"pivot_old_regions.npz") as src:
            origin=src[case.case_id+"/"+source]
        x=(origin-LO)*(GRID/(HI-LO))
        inner=np.stack([np.maximum(0,np.floor(x[:,0])-1),np.minimum(GRID,np.ceil(x[:,1])+1)],axis=1).astype(np.int64)
        oldwork={}
    groups=group_cover(inner,counts)
    i,j=np.triu_indices(len(groups))
    dist_axes=np.maximum(groups[i,1]-groups[j,0],groups[j,1]-groups[i,0])
    upper=np.nextafter((HI-LO)/GRID*np.sqrt(np.sum(dist_axes**2,axis=1)/K)*(1.+1.e-14),np.inf)
    no_new=float(upper.max(initial=0))
    geometry=json.dumps(case.geometry,sort_keys=True)
    gamma=min(case.gain_radius,case.hard_gain_radius or 1.)
    eps=case.additive_radius
    # g F0 = D-e. Both alternatives use the same measured D.
    # The old ellipsoid implies ||e||<=eps. R12 also has the original product calibration model,
    # but the primary comparison retains the larger additive outer bound.
    records={}
    evals=0
    cache={}
    def disks(d1,d0):
        nonlocal evals
        key=(d1,d0)
        if key not in cache:
            table=ratio_table(geometry,d1,d0)
            vv=ratio_enclosure(table.data,table.exp,table.degree,groups)
            cache[key]=physical_disk(vv,gamma)
            evals+=len(groups)
        return cache[key]
    def vnorm(x):
        return np.maximum(0.,next_down(np.sqrt(np.sum(x*x,axis=1))*(1.-1.e-13)-1.e-300))
    def compare(name, lower, plain):
        action=ACTIONS[name]
        r=float(np.nextafter(case.new_noise_radius(action)*(1.+1.e-13),np.inf))
        possible=lower<=2*r
        standard=plain<=2*r
        k=int(np.argmax(np.where(possible,upper,-1)))
        return dict(D_upper=float(upper[k]),absolute_only_D_upper=float(np.max(upper[standard],initial=0)),
            unresolved_pairs=int(possible.sum()),excluded_pairs=int((~possible).sum()),
            absolute_excluded=int((~standard).sum()),extra_contrast_exclusions=int(np.sum(standard&~possible)),
            max_separation_lower=float(lower.max(initial=0.)),new_noise_radius=r,
            worst_groups=[int(i[k]),int(j[k])],half_band=action.max_offset,frequency_blocks=action.exposures)
    repeat=compare("repeat",np.zeros_like(upper),np.zeros_like(upper))
    records["repeat"]=repeat
    for b in bands:
        one=[]
        for delta,name in ((-b,f"one_minus_{b:g}"),(b,f"one_plus_{b:g}")):
            dc,dr=disks(delta,None)
            direct=vnorm(disk_gap(dc,dr,i,j))
            ac,ar=disks(delta,0.)
            inc=disk_gap(ac,ar,i,j)
            # In Hilbert norm the difference of the two old additive errors is <=2eps.
            lower=np.maximum(direct,np.maximum(0.,next_down(vnorm(inc)-2*eps)))
            records[name]=compare(name,lower,direct)
            one.append((direct,lower,dc,dr,ac,ar))
        dc,dr=disks(b,-b)
        contrast=vnorm(disk_gap(dc,dr,i,j))
        # Orthogonal sum/difference of the two equal-energy exposures.
        ac=(one[0][4]+one[1][4])/2
        ar=(one[0][5]+one[1][5])/2 + 5.e-14*(abs(one[0][4])+abs(one[1][4])+1.e-300)
        avg=np.maximum(0.,next_down(vnorm(disk_gap(ac,ar,i,j))-2*eps))
        sep=np.sqrt(avg**2+.25*contrast**2)*(1.-1.e-13)
        raw=np.sqrt(.5*(one[0][0]**2+one[1][0]**2))*(1.-1.e-13)
        paired=np.sqrt(.5*(one[0][1]**2+one[1][1]**2))*(1.-1.e-13)
        records[f"split_{b:g}"]=compare(f"split_{b:g}",np.maximum(sep,paired),raw)
    def choose(tau):
        if no_new<=tau:
            return {"action":"old_only","fallback":False,"D_upper":no_new}
        eligible=[n for n,v in records.items() if v["D_upper"]<=tau]
        if not eligible:
            return {"action":"split_0.02","fallback":True,"D_upper":records["split_0.02"]["D_upper"]}
        n=min(eligible,key=lambda v:(ACTIONS[v].exposures,ACTIONS[v].max_offset,records[v]["D_upper"],v))
        return {"action":n,"fallback":False,"D_upper":records[n]["D_upper"]}
    return dict(case_id=case.case_id,cohort=case.cohort,method=source,old_work=oldwork,num_old_boxes=len(inner),
                group_count=len(groups),pair_count=len(i),groups=groups.tolist(),
                zero_new_D_upper=no_new,actions=records,selected={str(t):choose(t) for t in TAUS},
                group_response_ranges=evals,wall_s=perf_counter()-tic,cpu_s=process_time()-cpu)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--case",action="append")
    p.add_argument("--groups",type=int,default=64)
    p.add_argument("--source",default="seeded")
    p.add_argument("--bands",nargs="+",type=float)
    p.add_argument("--tag",default="")
    args=p.parse_args()
    path=RESULTS/("frequency_covers"+args.tag)
    path.mkdir(exist_ok=True)
    for case in load_cases():
        if args.case and case.case_id not in args.case:
            continue
        row=score_case(case,args.groups,args.source,args.bands or BANDS)
        (path/(case.case_id+".json")).write_text(json.dumps(row,allow_nan=False)+"\n")
        print(case.case_id,"old",round(row["zero_new_D_upper"],3),"new",
              round(row["actions"]["split_0.02"]["D_upper"],3),row["selected"]["0.5"]["action"],
              round(row["wall_s"],2),flush=True)

if __name__=="__main__":
    main()
