"""Matched physical-call comparison with a full-material, task-constrained GN step."""

import argparse
from collections import Counter
import json
from pathlib import Path
import time

from run import AMPLITUDES, OUT, audited_forward, ball_status, box_margin, dump_json, layout, norm_interval, strongest_near_null
from archives import AUS, RING, austria_observation, r10_observation
from replay import SOURCE

import numpy as np
from scipy.linalg import svdvals
from scipy.optimize import brentq

SIGNED_TARGET = {RING: .75, AUS: .175}
GN_DAMPING = 1e-4
GN_ITERATIONS = 256


def project_box_and_task(trial, weights, target):
    relevant = weights > 0
    a = weights[relevant]
    x = trial[relevant]
    if not abs(target) < np.sum(a):
        raise ArithmeticError('The requested task increment exceeds the symmetric material box')
    low = np.min((x - 1) / a) - 1
    high = np.max((x + 1) / a) + 1
    shift = brentq(lambda mu: a @ np.clip(x - mu * a, -1, 1) - target,
                   low, high, xtol=1e-13)
    result = np.clip(trial - shift * weights, -1, 1)
    if abs(weights @ result - target) > 1e-11:
        raise ArithmeticError('Task equality could not be solved to tolerance')
    return result


def gauss_newton(J, residual, margin, q, target):
    """One damped GN model, with two signed box-and-task least-squares solves."""
    active = margin > 0
    D = margin[active]
    scaled = J[:, active] * D
    sigma = float(svdvals(scaled, check_finite=False)[0])
    scaled /= sigma
    rhs = residual / sigma
    H = scaled.T @ scaled
    H.flat[::len(H)+1] += GN_DAMPING
    linear = scaled.T @ rhs
    step = 1 / (1 + GN_DAMPING)
    weights = q[active] * D
    directions, notes = {}, {}
    for sign in (-1, 1):
        signed_target = sign * target
        z = project_box_and_task(np.zeros_like(D), weights, signed_target)
        extrapolated, momentum = z, 1.
        for _ in range(GN_ITERATIONS):
            updated = project_box_and_task(extrapolated - step * (H @ extrapolated + linear), weights, signed_target)
            next_momentum = (1 + np.sqrt(1 + 4 * momentum**2)) / 2
            extrapolated = updated + (momentum - 1) / next_momentum * (updated - z)
            z, momentum = updated, next_momentum
        delta = np.zeros_like(margin)
        delta[active] = D * z
        directions[sign] = delta
        stationarity = z - project_box_and_task(z - step * (H @ z + linear), weights, signed_target)
        notes[str(sign)] = {
            'projected_gradient_step_error': float(np.linalg.norm(stationarity)),
            'linearized_data_residual_norm': float(np.linalg.norm(residual + J @ delta)),
            'task_equality_abs_error': float(abs(q @ delta - signed_target)),
            'material_rms_step': float(np.linalg.norm(delta) / np.sqrt(len(margin))),
        }
    return directions, {'sigma_max_box_scaled_J': sigma, 'ridge_on_scaled_coordinates': GN_DAMPING,
                        'projected_least_squares_iterations_each_sign': GN_ITERATIONS, 'signed_solves': notes}


def material_and_field(obs, m, q, name, arrays):
    lo, hi = obs.bounds
    if np.min(m) < lo - 1e-13 or np.max(m) > hi + 1e-13:
        raise ArithmeticError(f'Illegal material at {name}')
    m = np.clip(m, lo, hi)
    state, _, audit = audited_forward(obs.model, m, name)
    y = state['y']
    arrays[f'{name}__material'] = m
    arrays[f'{name}__prediction'] = y
    archived_distance = norm_interval(y, obs.Y, audit['field_error_bound'])
    return state, {
        'array_prefix': name,
        'regional_chi': float(q @ m),
        'archived_distance_interval': archived_distance,
        'relative_observation_residual': float(np.linalg.norm(y - obs.Y) / np.linalg.norm(obs.Y)),
        'numerics': audit,
    }


def run_arm(obs, center, q, method, arrays):
    started = time.perf_counter(), time.process_time()
    before = Counter(obs.model.counts)
    m0 = obs.endpoints[center].m
    margin = box_margin(m0, obs.bounds)
    base, center_record = material_and_field(obs, m0, q, f'{method}__initial', arrays)
    J = obs.model.jacobian(base)
    target = SIGNED_TARGET[obs.id]
    proposal_started = time.perf_counter(), time.process_time()
    if method == 'weak_projection':
        h, _, details = strongest_near_null(J, margin, q)
        gain = target / (q @ h)
        if not 0 < gain <= 1:
            raise ArithmeticError('Frozen signed target exceeds the projected direction range')
        plus = h * gain
        proposals = {-1: -plus, 1: plus}
        details['scale_of_projected_box_extremum'] = float(gain)
    else:
        proposals, details = gauss_newton(J, obs.model.pack(base['y'] - obs.Y), margin, q, target)
    details['total_proposal_wall_seconds'] = time.perf_counter() - proposal_started[0]
    details['total_proposal_CPU_seconds'] = time.process_time() - proposal_started[1]
    candidates = [{'signed_task_target': 0., 'backtracking_fraction': 0., **center_record}]
    for sign, h in proposals.items():
        for alpha in AMPLITUDES:
            label = f'{method}__{"plus" if sign > 0 else "minus"}__{int(alpha*100):03d}'
            _, candidate = material_and_field(obs, m0 + alpha * h, q, label, arrays)
            candidates.append({'signed_task_target': sign * target * alpha, 'backtracking_fraction': alpha, **candidate})
    cost = dict(obs.model.counts - before)
    if cost['LU'] != 7 or cost['FullJ'] != 1:
        raise RuntimeError('The matched physical-call budget was exceeded')
    return {'proposal': details, 'cost': cost, 'wall_seconds': time.perf_counter() - started[0],
            'cpu_seconds': time.process_time() - started[1], 'candidates': candidates}


def feasible_span(candidates, radius):
    admitted = [row for row in candidates if ball_status(row['archived_distance_interval'], radius) == 'inside']
    if not admitted:
        return {'certified_inside': 0, 'feasible_span': None}
    minimum = min(admitted, key=lambda row: row['regional_chi'])
    maximum = max(admitted, key=lambda row: row['regional_chi'])
    span = maximum['regional_chi'] - minimum['regional_chi']
    return {'certified_inside': len(admitted), 'feasible_span': span, 'two_point_error_lower_bound': span / 2,
            'minimum': minimum, 'maximum': maximum}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args()
    output = args.output
    original_seal = json.loads((output / 'ENDPOINT_SEAL.json').read_text(encoding='utf-8'))
    arms = ('weak_projection', 'task_constrained_gauss_newton')
    frozen = {
        'context': 'Comparison added after viewing initial study summaries. Neither arm can read a material label.',
        'original_online_seal': original_seal['cases'],
        'signed_task_target': SIGNED_TARGET,
        'common_pipeline': 'Each arm evaluates the same historical center, builds one full Jacobian, and makes '
                           'three full-wave calls per sign. Both use the same symmetric per-pixel material box. '
                           'All seven fields per arm get separate inverse-residual audits.',
        'backtracking_fractions': list(AMPLITUDES), 'matrix_GN_iterations_each_sign': GN_ITERATIONS,
        'result': {},
    }
    for make, center in ((r10_observation, 'r12_g19'), (austria_observation, 'cie_dual')):
        obs = make()
        _, q = layout(obs)
        arrays = {'Y_archive': obs.Y, 'task_vector': q}
        outcomes = {method: run_arm(obs, center, q, method, arrays) for method in arms}
        if outcomes[arms[0]]['cost'] != outcomes[arms[1]]['cost']:
            raise RuntimeError('The actual physical call counts differ between methods')
        path = output / f'{obs.id}__comparison_arrays.npz'
        np.savez_compressed(path, **arrays)
        frozen['result'][obs.id] = {'initial': center, 'archive_norm': float(np.linalg.norm(obs.Y)),
                                   'radius': obs.radius, 'arrays': path.name, 'methods': outcomes}
    dump_json(output / 'GN_COMPARISON_ENDPOINT_SEAL.json', frozen)

    offline_radii = json.loads((SOURCE / 'OFFLINE_SCORE.json').read_text(encoding='utf-8'))
    for case, record in frozen['result'].items():
        old_audit = json.loads((output / f'{case}.json').read_text(encoding='utf-8'))
        old_candidates = [{'array_prefix': name, 'regional_chi': hist['regional_chi'],
                           'archived_distance_interval': hist['distance_to_archive_interval'],
                           'source': hist['source']['file']}
                          for name, hist in old_audit['historical_baselines'].items()]
        radii = {'historical_archive_ball': record['radius']}
        if case == AUS:
            radii['original_R8_noise_ball'] = offline_radii[case]['original_r8_noise_absolute_radius_reconstructed_offline']
        comparands = {method: row['candidates'] for method, row in record['methods'].items()}
        comparands['historical_saved_endpoints'] = old_candidates
        record['spans'] = {
            name: {'radius': radius, 'relative_radius': radius / record['archive_norm'],
                   'spans_by_method': {method: feasible_span(candidates, radius)
                                       for method, candidates in comparands.items()}}
            for name, radius in radii.items()
        }
        for name, row in record['spans'].items():
            print(case, name, {method: summary['feasible_span']
                               for method, summary in row['spans_by_method'].items()}, flush=True)
    dump_json(output / 'GN_BUDGET_COMPARISON.json', frozen)


if __name__ == '__main__':
    main()
