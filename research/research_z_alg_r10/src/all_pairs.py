"""Combine all two-coordinate enclosures; score exclusions with whole-box distances."""
import argparse
from itertools import combinations
import json
import time
import numpy as np
from flint import arb
from bounded_certificate import ROOT, feasible_witness
from pivot_certificate import cell_tables, cover, tensor, intervals

BOX=(0.5,30.)
DEPTH=8


def distance(pixels,pair,z,step):
    if not len(pixels):
        raise ArithmeticError("Empty certified set: investigate data/model incompatibility")
    lower=BOX[0]+pixels*step
    upper=lower+step
    point=z[list(pair)]
    # Dyadic box edges are exact. Round the final binary subtraction down by one ulp.
    delta=np.maximum(np.maximum(lower-point,point-upper),0)
    return max(0.,float(np.nextafter(np.min(np.max(delta,axis=1)),0.)))


def analyze(d,name,rule):
    tic=time.process_time()
    key=f'{name}/{rule}'
    B,H,R,host,y=(d[key+'/B'],d[name+'/H'],d[name+'/R'],d[key+'/y_host'],d[key+'/y'])
    truth=d[name+'/truth']
    witnesses=[truth,*[z for z in d[name+'/clean_witnesses'] if np.linalg.norm(z-truth)>1e-4]]
    pivot=set(d[name+'/pivots'].tolist())
    table=cell_tables(H,R,B,host,y,BOX)
    k=len(R)
    polys={triple:tensor(table,list(triple)) for triple in combinations(range(k),3)}
    resolution=2**DEPTH;step=(BOX[1]-BOX[0])/resolution
    marginal=np.ones((k,resolution),dtype=bool)
    gap=np.zeros(len(witnesses));first=np.zeros(len(witnesses))
    binary_masks=[];area=[];checked=0
    for pair in combinations(range(k),2):
        pix,visits=cover(table,list(pair),DEPTH,polys)
        checked+=visits
        mask=np.zeros((resolution,resolution),dtype=bool)
        mask[tuple(pix.T)]=True
        binary_masks.append(np.packbits(mask))
        area.append(len(pix)/resolution**2)
        for j,i in enumerate(pair):
            marginal[i]&=np.any(mask,axis=1-j)
        dist=np.array([distance(pix,pair,z,step) for z in witnesses])
        gap=np.maximum(gap,dist)
        if set(pair)==pivot:
            first=dist
    projected=[intervals(np.flatnonzero(v),BOX,DEPTH) for v in marginal]
    hull_widths=[s[-1][1]-s[0][0] for s in projected]
    with_all=[feasible_witness(z,H,R,B,host,y) for z in witnesses]
    with_old=[feasible_witness(z,H,R,B[:,:-1],host[:-1],y[:-1]) for z in witnesses]
    rec=dict(candidates=d[name+'/central_roots'].tolist(), old_witness_count=len(witnesses),
             old_witnesses=with_old, new_witnesses=with_all,
             witness_vectors=[z.tolist() for z in witnesses],one_chart_linf_radii=first.tolist(),
             all_charts_linf_radii=gap.tolist(), coordinate_sets=projected,mean_pair_area=float(np.mean(area)),
             min_pair_area=min(area),max_pair_area=max(area),
             rms_diameter_bound=float(np.nextafter(float((sum((arb(w)**2 for w in hull_widths),arb(0))/k).sqrt().upper()),np.inf)),
             inspected_pivot_boxes=checked,cpu_seconds=time.process_time()-tic)
    return rec,np.asarray(binary_masks)


def run(cohort):
    d=np.load(ROOT/'results'/f'{cohort}_inputs.npz')
    results={};masks={}
    for name in d['case_ids'].tolist():
        previous=[]
        for rule in ('active','random','fisher'):
            key=name+'/'+rule
            identical=next((p for p in previous if np.array_equal(d[p+'/B'],d[key+'/B']) and np.array_equal(d[p+'/y'],d[key+'/y'])),None)
            if identical:
                results[key]=dict(results[identical],same_data_as=identical,cpu_seconds=0.0)
                masks[key]=masks[identical]
            else:
                results[key],masks[key]=analyze(d,name,rule)
                previous.append(key)
        print(name,'all-chart active certified radii',np.round(results[name+'/active']['all_charts_linf_radii'],3),flush=True)
    out=dict(precision_bits=128,box=BOX,depth=DEPTH,pairs=[list(v) for v in combinations(range(6),2)],
             noise_scales=[.0002,.01],records=results)
    (ROOT/'results'/f'{cohort}_all_pairs.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    np.savez_compressed(ROOT/'results'/f'{cohort}_all_pairs_masks.npz',**masks)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('cohort',choices=['natural','confirm'],default='natural',nargs='?')
    run(p.parse_args().cohort)

