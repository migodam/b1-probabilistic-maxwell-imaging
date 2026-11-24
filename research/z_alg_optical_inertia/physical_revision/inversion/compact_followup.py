"""Replay exact discrepancy and exposed-object cap continuation."""

import json

import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

from compact import (DenseModel, HIGH, LOW, NATIVE, OUT, STAGES, enclosure, get_observations,
                     norm, ball, safe_float)
from fit import native, search_stage
from gate import clocks, elapsed, write

MENU = (.5, .6, .75, 1., 1.25, 2.)
CALLS_PER_CAP = 16


def audit(geometry, obs, candidate):
    p, rem, receipt = enclosure(geometry, candidate['material'])
    d = norm(p-ball(obs.y))
    return dict(residual_lower=max(0.,safe_float(d-rem,upper=False)),
                residual_upper=safe_float(d+rem,upper=True), **receipt)


def continuation(geometry, obs, skip):
    start = clocks()
    x = None
    result = None
    rows = []
    for cap in MENU[1 if skip else 0:]:
        regular, candidate, row = search_stage(geometry,obs,x,9,CALLS_PER_CAP,cap)
        x = regular['x']
        proof = audit(geometry,obs,candidate)
        rows.append(dict(**row,physical_residual=proof))
        result = candidate
        if proof['residual_upper']<=obs.padded_radius:
            break
    return result, dict(stages=rows, found=rows[-1]['physical_residual']['residual_upper']<=obs.padded_radius,
                        tested_cap=rows[-1]['cap'], peak=result['peak'],
                        calls=sum(row['actual_calls'] for row in rows), **elapsed(start))


def exact_policy(geometry, obs):
    start = clocks()
    initial, minres, first = search_stage(geometry,obs,None,*STAGES[0],LOW)
    proof = audit(geometry,obs,minres)
    take = proof['residual_lower']>obs.padded_radius
    stages = [first]
    for side, calls in STAGES[1:]:
        initial, minres, stage = search_stage(geometry,obs,initial['x'],side,calls,HIGH if take else LOW)
        stages.append(stage)
    return initial, dict(stages=stages,releases=take,first_stage_audit=proof,
                         final_audit=audit(geometry,obs,initial),
                         calls=sum(r['actual_calls'] for r in stages), **elapsed(start))


def equal_calls(geometry, obs, cut):
    if cut>STAGES[0][1]:
        raise ArithmeticError('Matched-call checkpoint lies beyond stage 1.')
    best, _, row = search_stage(geometry,obs,None,STAGES[0][0],cut,HIGH)
    return best, dict(stage=row, final_audit=audit(geometry,obs,best), calls=row['actual_calls'])


def assess(geometry, observations, folder):
    with np.load(OUT/'offline_truth/reference.npz',allow_pickle=False) as f:
        truth, clean, held_clean = f['material'], f['clean'], f['held']
        held_geometry = geometry[0], f['held_sensor'], f['held_incident']
    scores = {}
    for obs in observations:
        with np.load(folder/(obs.id+'.npz'),allow_pickle=False) as data:
            for name in (key.removesuffix('_material') for key in data if key.endswith('_material')):
                estimate = data[name+'_material']
                p, r, a = enclosure(geometry,estimate)
                h, rh, ah = enclosure(held_geometry,estimate)
                delta = estimate-truth
                inside = truth>0
                dp = norm(p-ball(obs.y))
                scores.setdefault(obs.id,{})[name] = dict(
                    full_relative_l2=float(np.linalg.norm(delta)/np.linalg.norm(truth)),
                    inside_relative_l2=float(np.linalg.norm(delta[inside])/np.linalg.norm(truth[inside])),
                    background_rms=float(np.sqrt(np.mean(delta[~inside]**2))), peak=float(estimate.max()),
                    clean_train_relative_l2=float(norm(p-ball(clean))/norm(ball(clean))),
                    withheld_relative_l2=float(norm(h-ball(held_clean))/norm(ball(held_clean))),
                    residual_over_radius_lower=max(0.,safe_float((dp-r)/obs.padded_radius,upper=False)),
                    residual_over_radius_upper=safe_float((dp+r)/obs.padded_radius,upper=True),
                    train_enclosure=a, held_enclosure=ah)
    write(folder/'offline_scores.json',scores)


def main():
    ctx.prec,ctx.threads = 112,2
    geometry, observations = get_observations()
    native.SIDE,native.N,native.Model = NATIVE,NATIVE**2,DenseModel
    folder = OUT/'exposed_followup'
    folder.mkdir(exist_ok=True)
    gates = json.loads((OUT/'gate.json').read_text())
    for obs,result in zip(observations,gates,strict=True):
        initial_outcome = json.loads((OUT/'fits'/(obs.id+'.json')).read_text())
        cut = initial_outcome['reports']['fixed_low']['charge_full_objective_calls']
        branches = {'exact_discrepancy':exact_policy(geometry,obs),
                    'high_from_start_matched_calls':equal_calls(geometry,obs,cut),
                    'class_guided_continuation':continuation(geometry,obs,result['strict_exclusion']),
                    'ordinary_continuation':continuation(geometry,obs,False)}
        np.savez_compressed(folder/(obs.id+'.npz'),Y=obs.y,
                            **{name+'_material':candidate['material'] for name,(candidate,row) in branches.items()},
                            **{name+'_ports':candidate['ports'] for name,(candidate,row) in branches.items()})
        write(folder/(obs.id+'.json'),dict(
            id=obs.id,radius=obs.padded_radius,branches={name:row for name,(candidate,row) in branches.items()}))
    assess(geometry,observations,folder)


if __name__=='__main__':
    with threadpool_limits(1):
        main()
