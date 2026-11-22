#!/usr/bin/env python3
"""Report matched parameters and verify independent noise-only exclusions."""
from scan import CAP, OUT, load, read, save

import math
from flint import arb
from core.certificate import certify, safe_float
from core.inertia import Prepared, clock, timing


def matched_parameters(result):
    observations = result['observations_by_frequency']
    matched = {}
    for name, branch in result['joint_parameter_results'].items():
        differences = []
        for point in branch:
            controls = [next(row for row in frequency if row['mu'] == point['mu'] and row['ratio'] == point['ratio'])
                        for frequency in observations]
            reps = zip(*(row['results'] for row in controls), point['results'], strict=True)
            differences.extend(j['lower_bound_numeric'] - math.hypot(a['lower_bound_numeric'], b['lower_bound_numeric'])
                               for a,b,j in reps)
        matched[name] = {'point_observation_pairs':len(differences), 'minimum_gain':min(differences),
                         'maximum_gain':max(differences), 'strict_numerical_gains':sum(d>1e-10 for d in differences),
                         'numerical_ties':sum(abs(d)<=1e-10 for d in differences)}
    return matched


def independent_proofs(result):
    started = clock()
    ids, frequencies = load()
    selected = [case for case,v in result['observation_scores'].items() if v['independent_excludes_sample_noise_fro_numeric']]
    verified = {}
    for observation in selected:
        k = ids.index(observation)
        source = result['observation_scores'][observation]
        certs = []
        radii = []
        for freq,use in zip(frequencies,source['independent_terms'],strict=True):
            problem = Prepared(freq.G,freq.S,freq.B,freq.Y[k])
            rad = freq.obs[k]['radius_noise_fro']
            radii.append(rad)
            report = certify(problem,CAP,use['mu'],use['mu']*use['ratio'],radii={'sample':rad},threads=1)
            certs.append(report)
        is_verified = all(report['status']=='Arb_verified' for report in certs)
        joint_lower = sum((arb(report['certified_frobenius_residual_lower'])**2 for report in certs),arb(0)) if is_verified else arb(0)
        joint_radius = sum((arb(rad)**2 for rad in radii),arb(0))
        verified[observation] = dict(independent_components=certs,
            status='Arb_verified' if is_verified else 'Arb_incomplete',
            joint_lower=safe_float(joint_lower.sqrt(),upper=False),
            sample_radius=safe_float(joint_radius.sqrt(),upper=True),
            strict_sample_ball_exclusion=bool(joint_lower>joint_radius))
        print(observation, verified[observation]['status'], verified[observation]['strict_sample_ball_exclusion'],flush=True)
        save('independent_audit.json', {'cases':verified, **timing(started)})


if __name__ == '__main__':
    numerical = read(OUT/'numerical.json')
    save('matched_parameters.json', matched_parameters(numerical))
    independent_proofs(numerical)
