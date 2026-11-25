"""Original-operator counterexamples for every permitted plain plane wave."""
import argparse
import json

import numpy as np
from flint import acb, arb, ctx
from threadpoolctl import threadpool_limits

from beam_common import HERE, atomic_write, clocks, elapsed, geometry, save
from interval_certify import (ExactBox, ProofGrid, ball, hi, midpoint, norm_upper, tr)
from native_paired import green_kernel_roundoff
from native_residual import ResidualBoundBox


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('study')
    args=parser.parse_args()
    side=int(args.study.rsplit('_',1)[-1])
    source=HERE.parent/'pair'/args.study
    prop=json.loads((source/'proposal.json').read_text())
    bank=np.load(source/'public_classes.npz')
    old=np.load(source/'old_readings.npz')
    start=clocks()
    ctx.prec=128
    ctx.threads=3
    geo=geometry(side)
    grid=ProofGrid(geo)
    pu=tr(ball(geo.B))*ball(geo.B)
    pu=acb(1/pu[0,0].real)*pu
    eps=None if geo.spectrum is None else green_kernel_roundoff(grid)
    rows={}
    stored={}
    for i,record in prop['rows'].items():
        p=i+'__'
        t=clocks()
        cells=bank[p+'cells']
        numeric,exact=grid.restricted(cells)
        ctr=bank[p+'centers']
        widths=bank[p+'widths']
        models=[ExactBox(exact,numeric,c,widths) if eps is None
                else ResidualBoundBox(exact,numeric,c,widths,eps*len(cells)) for c in ctr]
        yc=[u.Y for u in models]
        uerrs=[getattr(u,'y_error',arb(0)) for u in models]
        at=old[p+'Y']
        rr=record['new_radius']
        ro=record['old_radius']
        e=record['fluence']
        angles=[]
        stored_mid=[]
        witness_old=[]
        for j in range(geo.B.shape[1]):
            if j==0:
                arm=ball(np.eye(geo.B.shape[1],dtype=complex)[:,[j]])
                witness_old=[hi((norm_upper(y*arm-ball(at))+err)/arb(ro))
                             for y,err in zip(yc,uerrs,strict=True)]
                continue
            amp=(arb(e)/pu[j,j].real).sqrt()
            arm=acb(amp)*ball(np.eye(geo.B.shape[1],dtype=complex)[:,[j]])
            centers=[x*arm for x in yc]
            guess=ball(midpoint(acb(arb(1)/2)*(centers[0]+centers[1])))
            candidate=[hi((norm_upper(v-guess)+err*amp)/arb(rr)) for v,err in zip(centers,uerrs,strict=True)]
            angles.append(dict(index=j,angle_degrees=j*360/32,verified_residual_over_new_radius_upper=candidate,
                               midpoint_is_legal_for_both=bool(max(candidate+witness_old)<1)))
            stored_mid.append(midpoint(guess))
        rows[i]=dict(old_residual_over_radius_upper=witness_old,
                     all_permitted_individual_angles_have_a_common_legal_reading=all(
                          a['midpoint_is_legal_for_both'] for a in angles),
                     maximum_midpoint_residual_over_new_radius_upper=max(max(a['verified_residual_over_new_radius_upper']) for a in angles),
                     fluence=e,new_radius=rr,angles=angles,**elapsed(t))
        stored[i+'__common_midpoint_by_angle']=np.column_stack(stored_mid)
        print(args.study,i,sum(a['midpoint_is_legal_for_both'] for a in angles),len(angles),flush=True)
    name=HERE/f'angle_bound_{args.study}'
    atomic_write(name.with_suffix('.json'),dict(side=side,
                 scope='All 15 permitted new single plane waves at each fixed common fluence; uses two original-operator, old-compatible material center witnesses.',
                 rows=rows,**elapsed(start)))
    save(name.with_suffix('.npz'),**stored)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
