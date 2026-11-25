"""Strict upper bounds for every sparse draw from two old-compatible materials."""
import argparse
from itertools import combinations
import json

import numpy as np
from flint import acb, arb, ctx
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from beam_common import HERE, atomic_write, clocks, elapsed, geometry, herm, save
from interval_certify import ExactBox, ProofGrid, ball, hi, midpoint, norm_upper, positive_lower, tr
from native_paired import green_kernel_roundoff, pad_all
from native_residual import ResidualBoundBox
from paired_replay import dimblock, matrix_values


SUPPORT_SIZES=(1,2,3)


def spectral_upper(C,W,ratio,energy,noise):
    cand=herm(midpoint(C))
    gram=herm(midpoint(W))
    eig=float(eigh(cand,gram,check_finite=False,subset_by_index=(len(cand)-1,len(cand)-1),eigvals_only=True)[0])
    if eig<=0.:
        eig=1.e-20
    for relative in (1.000001,1.00001,1.0001,1.001,1.01,1.1):
        bound=arb(eig)*arb(relative)
        if positive_lower(acb(bound)*W-C) is not None:
            v=(bound*energy).sqrt()/noise
            return hi(v),bool(v<ratio)
    return None,False


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('side',type=int,choices=(32,96))
    args=parser.parse_args()
    side=args.side
    src=HERE.parent/'pair'/f'power_{side}'
    record=json.loads((src/'proposal.json').read_text())
    bank=np.load(src/'public_classes.npz')
    old=np.load(src/'old_readings.npz')
    ctx.prec=128
    ctx.threads=2
    start=clocks()
    geom=geometry(side)
    grid=ProofGrid(geom)
    b=ball(geom.B)
    power=tr(b)*b
    power=acb(1/power[0,0].real)*power
    eps=None if geom.spectrum is None else green_kernel_roundoff(grid)
    rows={}
    storage={}
    for i,r in record['rows'].items():
        t=clocks()
        key=i+'__'
        cells=bank[key+'cells']
        numeric,exact=grid.restricted(cells)
        ctr=bank[key+'centers']
        width=bank[key+'widths']
        models=[ExactBox(exact,numeric,c,width) if eps is None
                else ResidualBoundBox(exact,numeric,c,width,eps*len(cells)) for c in ctr]
        y=[pad_all(x.Y,getattr(x,'y_error',arb(0))) for x in models]
        v0=ball(np.eye(b.ncols())[:,[0]])
        old_y=ball(old[key+'Y'])
        ro=arb(r['old_radius'])
        verified_old=[hi(norm_upper(v*v0-old_y)/ro) for v in y]
        diff=y[1]-y[0]
        powerdiff=tr(diff)*diff
        ene=arb(r['fluence'])
        noise=arb(r['new_radius'])
        counts={}
        for k in SUPPORT_SIZES:
            start_k=clocks()
            all_values=[]
            for sub in combinations(range(1,b.ncols()),k):
                eig,over=spectral_upper(dimblock(powerdiff,sub,sub),dimblock(power,sub,sub),2,ene,noise)
                all_values.append(dict(indices=list(sub),maximum_center_gap_over_noise_upper=eig,
                                       common_reading_implied_for_any_coefficients=bool(over and max(verified_old)<1.)))
            numbers=[q['maximum_center_gap_over_noise_upper'] for q in all_values]
            finite=[x for x in numbers if x is not None]
            counts[str(k)]=dict(all_combinations=len(all_values),all_verified=len(finite)==len(numbers),
                                 number_guaranteed_ambiguous=sum(q['common_reading_implied_for_any_coefficients'] for q in all_values),
                                 max_upper_over_all_supports=max(finite) if len(finite)==len(numbers) else None,
                                 all_supports_guaranteed_ambiguous=all(q['common_reading_implied_for_any_coefficients'] for q in all_values),
                                 rows=all_values,**elapsed(start_k))
        rows[i]=dict(old_center_residual_over_radius_upper=verified_old,fluence=r['fluence'],caps=counts,**elapsed(t))
        storage[key+'difference_midpoint'],storage[key+'difference_radii']=matrix_values(diff)
        print(i,[(k,v['number_guaranteed_ambiguous'],v['all_combinations'],v['max_upper_over_all_supports']) for k,v in counts.items()],flush=True)
    name=HERE/f'sparse_physical_power_{side}'
    atomic_write(name.with_suffix('.json'),dict(side=side,
                 scope='Same old-compatible whole materials; all exact support-specific original response difference norm upper bounds at the registered common fluence. Strict values below two certify a common legal reading for every complex coefficient vector on that support.',
                 rows=rows,**elapsed(start)))
    save(name.with_suffix('.npz'),**storage)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
