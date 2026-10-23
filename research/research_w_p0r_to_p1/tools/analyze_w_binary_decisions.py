"""Exact two-fixed-mean Gaussian decision control, using saved endpoints only.

This is a separately assigned two-atom prior, NOT the continuous P1 family,
the unknown-wall posterior, or an online reconstruction experiment. No PDE.
"""
import json
from pathlib import Path
import time
import numpy as np
from scipy.special import expit, ndtr

ROOT = Path(__file__).resolve().parents[1]


def binary_risk(distance, separation, order):
    nodes, weights = np.polynomial.hermite.hermgauss(order)
    # Under world 1, log(p(y|1)/p(y|0)) = d^2/2 + d Z, Z~N(0,1).
    log_odds = distance * np.sqrt(2) * nodes + distance**2 / 2
    variance = expit(log_odds) * expit(-log_odds)
    return float(separation**2 * (weights @ variance) / np.sqrt(np.pi))


def main():
    start = time.process_time()
    source = ROOT / 'results/B2/validation/joint_seed2_narrow_fullj/joint8_PAIR.npz'
    with np.load(source) as a:
        distance = float(np.linalg.norm((a['fine1'] - a['fine0']) / .01))
    separation = .1
    orders = [41, 81, 161]
    risks = [binary_risk(distance, separation, n) for n in orders]
    error = float(ndtr(-distance / 2))
    prior_risk = separation**2 / 4
    point_risk = separation**2 * error
    result = dict(
        status='EXACT_BINARY_GAUSSIAN_CONTROL_NUMERICAL_INTEGRATION',
        source=str(source.relative_to(ROOT)),
        model='Two fixed full-Maxwell means with known common Gaussian noise, gain=1; explicitly assigned masses .5/.5.',
        independent_of_continuous_P1_reference=True,
        mean_distance_white=distance, target_position_separation=separation,
        optimal_classification_error=error,
        posterior_mean_squared_risk_by_GH_order=dict(zip(map(str, orders), risks)),
        last_order_absolute_change=abs(risks[-1]-risks[-2]),
        zero_observation_prior_mean_squared_risk=prior_risk,
        hard_endpoint_MAP_squared_risk=point_risk,
        posterior_mean_relative_gain_vs_prior=1-risks[-1]/prior_risk,
        posterior_mean_relative_gain_vs_endpoint_MAP=1-risks[-1]/point_risk,
        nearest_endpoint_two_point_risk_lower_bound=separation**2*error/4,
        interpretation='A mean can beat a hard selected world by loss-aware averaging alone. Such a gain is not evidence of progressive scheduling or a new Maxwell algorithm.',
        exclusions=['Not continuous-P1 prior risk', 'Not measured held-out performance',
                    'Not unknown-wall posterior mass', 'Not a new physical endpoint'],
        new_rhs=0, prior_source_physics_cost='Already paid and retained in B2.',
        pure_analysis_cpu_seconds=time.process_time()-start)
    out = ROOT / 'results/B2/binary_probability_diagnostic'
    out.mkdir(parents=True, exist_ok=True)
    (out/'RESULTS.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result, allow_nan=False))


if __name__ == '__main__':
    main()
