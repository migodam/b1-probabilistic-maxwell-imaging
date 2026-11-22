"""Replay chosen profiles on all 9216 cells and independently certify the distances."""

import argparse
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from flint import arb, ctx
from threadpoolctl import threadpool_limits

from curve_inputs import CAPS, HERE, full_ring_inputs
from envelope import FourierEnvelope
from native import Geometry, Native, Weights, low
from search import Projected, local_arms, times, save_json


def from_trial(cap, scale, r):
    real, optical = tuple(r['weights']),tuple(r['optical'])
    if len(set(real))==len(set(optical))==1:
        scale=real[0]
    return Weights(cap,scale,real,optical)


def rank(r):
    return (r['d'],r.get('first_retained_eigenvalue',0.))


def candidate_pool(row, arm):
    u=row['cap']
    global_proposals=[r for r in row['arms'].values() if r is not None]
    choice=[]
    if arm == 'global':
        if u==.5:
            # Inherited shared controls maximize efficiency and guarantee a stable
            # exact comparison to the original radius. Smaller controls are fallbacks.
            for mu in (.0025,.001,.00035):
                choice.append((Weights.global_(u,mu,0.),0,'inherited_real'))
        else:
            for r in sorted(global_proposals+row['best_global_profiles'],key=rank,reverse=True):
                choice.append((Weights.global_(u,r['mu'],r['t']),r['nu'],
                               'global_real' if r['t']==0. else 'global_optical'))
    if arm.startswith('local') or arm=='matched':
        trials=[]
        for data in row['local'].values():
            if data['status']!='evaluated':
                continue
            for z in data['trials']:
                if not z['same_inertia']:
                    continue
                keep={'local':z['arm'].endswith('blocks'),
                      'local_real':z['arm']=='real_blocks',
                      'local_optical':z['arm']=='optical_blocks',
                      'local_joint':z['arm']=='joint_blocks',
                      'matched':z['arm'].endswith('scalar')}[arm]
                if keep:
                    trials.append((data['baseline']['mu'],z))
        for scale,r in sorted(trials,key=lambda data:rank(data[1]),reverse=True):
            choice.append((from_trial(u,scale,r),r['nu'],r['arm']))
    return list(dict.fromkeys(choice))


def evaluate(geometry,envelope,w,nu_expected,y,profiles):
    if w not in profiles:
        start,cpu=perf_counter(),process_time()
        counters=geometry.cost.copy()
        native=Native(geometry,w)
        modes=native.modes(nu_expected)
        modes['native_mode_wall_seconds'] = modes.pop('wall_seconds')
        modes['native_mode_cpu_seconds'] = modes.pop('cpu_seconds')
        nu=native.nu
        if nu>=geometry.s:
            neg=native.ensure_negative(geometry.s)
            record=dict(status='illumination_budget_exhausted_verified',
                        n_negative_at_least=geometry.s, negative_block_upper=neg)
        else:
            proof=envelope.verify(w,nu)
            neg=native.ensure_negative(nu)
            record=dict(status='native_gap_and_inertia_verified',exact_nu=nu,
                        negative_block_upper=neg,**proof)
        record=dict(**w.as_dict(),full_dimension=geometry.n,expected_projection_nu=nu_expected,
                    **modes,**record,profile_operations=dict(geometry.cost-counters),
                    **{f'profile_{k}':v for k,v in times(start,cpu).items()})
        profiles[w]=native,record
    native,record=profiles[w]
    if record['status']=='illumination_budget_exhausted_verified':
        return dict(distance_lower=0., status=record['status'])
    start,cpu=perf_counter(),process_time()
    counters=geometry.cost.copy()
    nu=record['exact_nu']
    xx,v,vals,diag=native.numerical_distance(y,nu)
    result=native.ensure_source(y,nu,record['positive_gap_lower'],xx,v,vals)
    return dict(**diag,**result,exact_nu=nu,status='native_distance_verified',
                operations=dict(geometry.cost-counters),
                **{f'case_{k}':v for k,v in times(start,cpu).items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cap',type=float,choices=CAPS,required=True)
    parser.add_argument('--arm',choices=('global','local','local_real','local_optical','local_joint','matched'),required=True)
    parser.add_argument('--object',type=str,action='append')
    parser.add_argument('--noise',type=float,action='append')
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--search',type=Path,default=HERE/'search.json')
    parser.add_argument('--stable-recovery',action='store_true')
    args=parser.parse_args()
    ctx.prec,ctx.threads=112,2
    start,cpu=perf_counter(),process_time()
    spec,s,b,observations=full_ring_inputs()
    geometry=Geometry(spec,s,b)
    envelope=FourierEnvelope(geometry)
    source=json.loads(args.search.read_text())
    relevant=[r for r in source['rows'] if r['cap']==args.cap
              and (not args.object or r['object_id'] in args.object)
              and (not args.noise or r['noise_pct'] in args.noise)]
    # The extra local crossing is rechecked first, then all positive global rows.
    relevant.sort(key=lambda r: (0 if r['object_id']=='optfresh_04' and r['noise_pct']==1 else 1,
                                 r['object_id'],r['noise_pct']))
    y={r.id:r.y for r in observations}
    obs={r.id:r for r in observations}
    projected=Projected(spec,s,b,source['q']) if args.stable_recovery and args.arm!='global' else None
    profiles={}
    irrecoverable_profiles={}
    receipts=[]
    for row in relevant:
        candidates=candidate_pool(row,args.arm)
        if args.stable_recovery:
            if args.arm=='global':
                candidates=[(Weights.global_(args.cap, .0001, .1), 4, 'stable_global_recovery')]
            else:
                stable=local_arms(projected,args.cap,dict(mu=.0001,t=.1),obs[row['id']])
                row['local']={'stable_global':stable}
                candidates=candidate_pool(row,args.arm)
        if not candidates:
            receipts.append(dict(id=row['id'],cap=args.cap,noise_pct=row['noise_pct'],
                                 radius=row['radius'],status='no_surviving_proposal',distance_lower=0.))
            continue
        failures=[]
        outcome=None
        for w,nu,method in candidates[:3]:
            if w in irrecoverable_profiles:
                failures.append(dict(error=irrecoverable_profiles[w],reused_failure=True,**w.as_dict()))
                continue
            try:
                result=evaluate(geometry,envelope,w,nu,y[row['id']],profiles)
            except (ArithmeticError,ValueError,RuntimeError) as e:
                if w not in profiles:
                    irrecoverable_profiles[w]=str(e)
                failures.append(dict(error=str(e),**w.as_dict()))
                continue
            distance=result['distance_lower']
            margin=low(arb(distance)-arb(row['radius']))
            result=dict(**result,margin_lower=margin,strictly_excludes=margin>0.,method=method,
                        profile_key=w.as_dict())
            if outcome is None or distance>outcome['distance_lower']:
                outcome=result
            if margin>0. or result.get('distance',0.)<=row['radius']:
                break
        if outcome is None:
            outcome=dict(status='verification_unresolved',distance_lower=0.,strictly_excludes=False)
        receipts.append(dict(id=row['id'],object_id=row['object_id'],cap=args.cap,noise_pct=row['noise_pct'],
                             radius=row['radius'],**outcome,unverified_profiles=failures,
                             elapsed_through_case=times(start,cpu)))
        save_json(args.out,dict(scope='Exact represented archived 96² Green multiplier and B,S,Y; native interval certificate',
                                arm=args.arm,cap=args.cap,search=args.search.name,
                                stable_recovery=args.stable_recovery,
                                precision_bits=ctx.prec,profiles=[d for v,d in profiles.values()],
                                rows=receipts,total_operations=dict(geometry.cost),**times(start,cpu)))
        print(args.arm, row['id'],outcome.get('exact_nu'),outcome.get('distance'),
              outcome['distance_lower'],row['radius'],outcome['status'],
              outcome.get('method'),len(failures),flush=True)


if __name__=='__main__':
    with threadpool_limits(1):
        main()
