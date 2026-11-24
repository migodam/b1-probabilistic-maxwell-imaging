"""Replay the pre-inversion strict native contrast gate on the fixed full ring."""

import csv
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
CURVE = STUDY/'submission/curve'
sys.path.insert(0, str(CURVE))
from curve_inputs import Observation, full_ring_inputs

REFERENCE = STUDY/'validation/srap/online.py'
spec = spec_from_file_location('srap_cached_verifier', REFERENCE)
cached = module_from_spec(spec)
spec.loader.exec_module(cached)


def clocks():
    return perf_counter(), process_time()


def elapsed(start):
    return dict(wall_seconds=perf_counter()-start[0], cpu_seconds=process_time()-start[1])


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def original():
    s, p, b, paired = full_ring_inputs()
    observations = [case for case in paired if case.noise_pct in (1., 2.)]
    feas = cached.read(STUDY/'results/theory/fine_feasibility.json')
    public = STUDY/'synthetic/public'
    for row in feas['observations']:
        oid = row['id']
        if row['geometry'] != 'PW360_10x20':
            continue
        with np.load(public/(oid+'.npz'), allow_pickle=False) as raw:
            y = raw['Y']
        if row['noise_fraction'] == .02:
            prior = next(case for case in observations if case.id == oid[:-3]+'pair2pct')
            if not np.array_equal(y, prior.y) or row['inflated_radius'] != prior.padded_radius:
                raise ArithmeticError('Paired 2% inputs differ from original 2% inputs.')
        elif row['noise_fraction'] == .10:
            observations.append(Observation(oid, oid.split('__')[0], 10.,
                                            row['sample_radius_binary64'], row['inflated_radius'], y))
    return (s, p, b), observations


def inherited_distances(observations):
    arms = ('strict_global_distance_lower', 'strict_local_distance_lower',
            'fixed_0_distance_lower', 'fixed_2_distance_lower')
    keyed = {case.id:case for case in observations}
    with (CURVE/'noise_curve.csv').open() as raw:
        for row in csv.DictReader(raw):
            obs = keyed.get(row['id'])
            cap = float(row['cap'])
            if obs is None or cap not in (.5, .75):
                continue
            if obs.padded_radius != float(row['radius']):
                raise ArithmeticError('Verification radius differs from the input.')
            for name in arms:
                if row[name]:
                    yield obs.id, cap, float(row[name]), name


def input_rows(observations, cap):
    for obs in observations:
        yield obs.y, dict(id=obs.id, object_id=obs.object_id, cap=cap, noise_pct=obs.noise_pct,
                          padded_radius=obs.padded_radius, sample_radius=obs.sample_radius)


def replay(geometry, observations, output, archived=True):
    ctx.prec, ctx.threads = 112, 2
    start = clocks()
    geometries, profiles, best = {}, [], {}
    r_to_use = {obs.id:obs.padded_radius for obs in observations}
    for profile in cached.choices():
        t = clocks()
        sub = profile['block_side']
        if sub not in geometries:
            timer = clocks()
            geometries[sub] = (cached.Geometry(*geometry, block_side=sub), elapsed(timer))
        g, loading = geometries[sub]
        pre = g.cost.copy()
        p = cached.saved_operator(g, profile)
        nu = profile['nu']
        ngap, neg_record = None, None
        offline = clocks()
        if nu:
            p.modes(nu)
            ng, ngap = p.ensure_negative(nu, return_gap=True)
            neg_record = dict(negative_congruence_upper=ng, negative_gap_lower=ngap, **elapsed(offline))
        cache_start = clocks()
        m, defect, cache_work = cached.exact_common_response(p, nu)
        resp_time = elapsed(cache_start)
        rows = []
        for y, meta in input_rows(observations, profile['cap']):
            online_start = clocks()
            result = cached.bound_source(y, p, nu, profile['positive_gap_lower'], ngap, m, defect, meta)
            result['source_verification'] = elapsed(online_start)
            rows.append(result)
            key = meta['id'], meta['cap']
            d = result['joint_distance_lower']
            if key not in best or d > best[key]['strict_distance_lower']:
                best[key] = dict(strict_distance_lower=d, source=profile['label'])
        profiles.append(dict(**profile, geometry_preparation=loading, negative_reverification=neg_record,
                             shared_stationary_work=cache_work, shared_stationary_total=resp_time,
                             offline_operations=dict(g.cost-pre), rows=rows, total=elapsed(t)))
        print(profile['label'], sum(row['joint_status']=='excluded' for row in rows),
              'of', len(rows), 'cache s', round(resp_time['wall_seconds'], 3), flush=True)
    if archived:
        for oid, cap, d, source in inherited_distances(observations):
            key = oid, cap
            if d > best[key]['strict_distance_lower']:
                best[key] = dict(strict_distance_lower=d, source='v2/'+source)
    decisions = []
    for (oid, cap), row in sorted(best.items()):
        radius = r_to_use[oid]
        d = row['strict_distance_lower']
        decisions.append(dict(id=oid, cap=cap, radius=radius, **row, strict_ratio=d/radius,
                              strict_exclusion=d>radius, requested_window_impossible=d>2.5*radius))
    write(output, dict(scope='Stored native 96² lossless scalar TM with specified complex Frobenius ball',
                       reused_proof_source=str(REFERENCE.relative_to(STUDY)), bits=ctx.prec, flint_threads=ctx.threads,
                       profiles=profiles, decisions=decisions, **elapsed(start)))
    return decisions


if __name__ == '__main__':
    with threadpool_limits(1):
        geometry, observations = original()
        replay(geometry, observations, HERE/'original_gate.json')
