"""Reconcile archived native lower bounds on exactly the original padded observation balls."""

import csv
import json

from search import HERE, STUDY, full_ring_inputs

PILOT=STUDY.parent/'z_alg_neural_dual_pilot/results/cases.csv'
CURVE=STUDY/'submission/curve/noise_curve.csv'


def strong(rows):
    if not rows:
        return None
    return max(rows,key=lambda row:(row['excluded_cap'],row['proven_distance_lower']))


def main():
    _,_,_,observations=full_ring_inputs()
    with PILOT.open() as fh:
        inherited=list(csv.DictReader(fh))
    with CURVE.open() as fh:
        additional=list(csv.DictReader(fh))
    results={}
    for obs in observations:
        if obs.noise_pct not in (1.,2.):
            continue
        candidates=[]
        all_attempts=[]
        for row in inherited:
            if (row['cohort']!='Fresh6' or row['oid']!=obs.object_id or float(row['noise'])!=obs.noise_pct or
                    (row['arm']!='fixed_portfolio' and (row['arm'],row['side'])!=('gradient1','8'))):
                continue
            if float(row['radius'])!=obs.padded_radius:
                raise ArithmeticError('Stored pilot certificate uses a different radius.')
            low=float(row['lower']) if row['lower'] else None
            attempts=dict(excluded_cap=float(row['u']),proven_distance_lower=low,pilot_arm=row['arm'],
                          archived_id=row['id'],radius=obs.padded_radius,status=row['decision'],
                          archived_proposal_wall=float(row['proposal_wall']),archived_proof_wall=float(row['proof_wall']),
                          archived_charged_wall=float(row['charged_wall']),archived_charged_cpu=float(row['charged_cpu']))
            all_attempts.append(attempts)
            if low is not None and low>obs.padded_radius:
                candidates.append(attempts)
        original=[]
        for row in additional:
            if row['id']!=obs.id:
                continue
            if float(row['radius'])!=obs.padded_radius:
                raise ArithmeticError('Stored Paper 1 certificate uses a different radius.')
            arms={'strict_global_distance_lower':'original_paper1_global',
                  'strict_local_distance_lower':'original_paper1_local_2x2',
                  'fixed_0_distance_lower':'original_paper1_fixed_4x4_profile_0',
                  'fixed_2_distance_lower':'original_paper1_fixed_4x4_profile_2'}
            for key,arm in arms.items():
                value=float(row[key]) if row[key] else 0.
                if value>obs.padded_radius:
                    original.append(dict(excluded_cap=float(row['cap']),proven_distance_lower=value,
                                         paper1_arm=arm,archived_id=obs.id,radius=obs.padded_radius))
        results[obs.id]=dict(object=obs.object_id,noise_pct=obs.noise_pct,radius=obs.padded_radius,
                                  strongest_fixed_plus_gradient64=strong(candidates),
                                  strongest_including_original_paper1=strong(candidates+original),
                                  native_attempts=all_attempts,original_native_passing_rows=original)
    (HERE/'lower_bounds.json').write_text(json.dumps(dict(observations=results,
           fixed_plus_gradient_receipt='research/z_alg_neural_dual_pilot/results/cases.csv',
           paper1_receipt='research/z_alg_optical_inertia/submission/curve/noise_curve.csv'),
           indent=2,allow_nan=False)+'\n')
    for key,row in results.items():
        print(key,row['strongest_fixed_plus_gradient64']['excluded_cap'],
              row['strongest_including_original_paper1']['excluded_cap'])


if __name__=='__main__':
    main()
