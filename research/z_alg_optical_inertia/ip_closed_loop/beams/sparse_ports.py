"""Restrict verified prospective quadratic forms to small coherent subsets."""
import argparse
from itertools import combinations
import json
from math import comb

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import eigh, cho_factor, cho_solve
from threadpoolctl import threadpool_limits

from beam_common import HERE, atomic_write, clocks, elapsed, geometry, herm, save
from interval_certify import ball, endpoints, tr
from paired_replay import projected_separation


SIZES=(1,2,3,4,15)


def restore(midpoint, radii):
    return acb_mat(midpoint.shape[0],midpoint.shape[1],
                   [acb(arb(float(a.real),float(r)),arb(float(a.imag),float(r)))
                    for a,r in zip(midpoint.flat,radii.flat,strict=True)])


def reduced_matrix(forms,dual,m):
    p=forms[0].shape[0]-m
    I=np.eye(p)
    a,b=forms
    h0,h1=dual
    receiving=np.block([[I+h0*a[m:,m:],-I],[-I,I+h1*b[m:,m:]]])
    coupling=np.vstack((h0*a[m:,:m],h1*b[m:,:m]))
    pure=h0*a[:m,:m]+h1*b[:m,:m]
    inverse=cho_factor(herm(receiving),check_finite=False)
    return herm(pure-coupling.conj().T@cho_solve(inverse,coupling,check_finite=False))


def best_support(H,gram,available,k):
    highest=-np.inf
    support=None
    beam=None
    for js in combinations(available,k):
        mm=H[np.ix_(js,js)]
        pp=gram[np.ix_(js,js)]
        values,vectors=eigh(mm,pp,check_finite=False,subset_by_index=(k-1,k-1))
        if values[0]>highest:
            highest=float(values[0])
            support=js
            beam=vectors[:,0]
    a=np.zeros(len(H),complex)
    a[list(support)]=beam
    return highest,support,a


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('side',type=int,choices=(32,96))
    args=parser.parse_args()
    side=args.side
    root=HERE.parent/'pair'/f'power_{side}'
    p=json.loads((root/'proposal.json').read_text())
    array=np.load(root/'public_classes.npz')
    prev=json.loads((HERE/f'paired_replay_power_{side}.json').read_text())
    all_angle=json.loads((HERE/f'angle_bound_power_{side}.json').read_text())
    certs=np.load(HERE/f'paired_replay_power_{side}.npz')
    ctx.prec=128
    ctx.threads=1
    start=clocks()
    geom=geometry(side)
    b=ball(geom.B)
    power=tr(b)*b
    power=acb(1/power[0,0].real)*power
    m=b.ncols()
    gram=array['power']
    rows={}
    storage={}
    for i,r in p['rows'].items():
        t=clocks()
        cases={}
        one_angle=all_angle['rows'][i]
        for mode in ('indefinite','positive'):
            mt=r['cases'][mode]
            if not all(x['old_saturation'] for x in prev['rows'][i]['outcomes'][mode]['forms']):
                continue
            dual=mt['dual_weights']
            mid=[certs[f'{i}__{mode}__{j}__q_midpoint'] for j in range(2)]
            rad=[certs[f'{i}__{mode}__{j}__q_entry_radii'] for j in range(2)]
            h=reduced_matrix(mid,dual,m)
            intervals=[restore(a,b) for a,b in zip(mid,rad,strict=True)]
            sigma=r['new_radius']
            results={}
            for k in SIZES:
                cand=clocks()
                value,indices,direction=best_support(h,gram,tuple(range(1,m)),k)
                a=ball(direction[:,None])
                energy=arb(r['fluence'])
                scale=(energy/(tr(a)*power*a)[0,0].real).sqrt()
                action=acb(scale)*a
                proof=projected_separation(intervals,action,dual,sigma,energy)
                if k==1 and proof['strict_separation'] and one_angle['all_permitted_individual_angles_have_a_common_legal_reading']:
                    raise ArithmeticError('Single-port certificate contradicts a legal common observation.')
                results[str(k)]=dict(indices=list(indices),angles_degrees=[j*360/32 for j in indices],
                                     supports_searched=comb(m-1,k),spectral_proposal=value,
                                     exact_common_fluence=r['fluence'],
                                     exact_scalar_normalization_interval=endpoints(scale),
                                     result=proof,**elapsed(cand))
                storage[f'{i}__{mode}__k{k}__unscaled'] = direction
            cases[mode]=dict(verified_sufficient_port_count=next((k for k in SIZES if results[str(k)]['result']['strict_separation']),None),
                             supports=results)
        rows[i]=dict(at_least_two_ports_required=one_angle['all_permitted_individual_angles_have_a_common_legal_reading'],
                     exact_allowed_fluence=r['fluence'],cases=cases,**elapsed(t))
        print(side,i,[(n,v['verified_sufficient_port_count'],
                       [round(x['result'].get('lower_distance_over_radius',0.),5) for x in v['supports'].values()])
                      for n,v in cases.items()],flush=True)
    name=HERE/f'sparse_ports_power_{side}'
    atomic_write(name.with_suffix('.json'),
                 dict(side=side,method='same frozen conditioned pair and distance multipliers; exhaustively propose a generalized Rayleigh maximum for each allowed support of sizes 1,2,3,4,15; independently verify one candidate for each size by restored certified Q enclosures',
                      scope='worst-case guarantee for fixed energy budget; unverified supports are unknown',rows=rows,**elapsed(start)))
    save(name.with_suffix('.npz'),**storage)


if __name__=='__main__':
    with threadpool_limits(limits=1):main()
