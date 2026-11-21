"""Full-wave material secants on saved 24²/32² Maxwell acquisitions."""

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import time

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[name] = '1'

import numpy as np
from scipy.linalg import lu_solve, svd

from archives import AUS, RING, SRAP_CASES, WORKSPACE, observations, srap_observation

OUT = Path(__file__).resolve().parent / 'results'
AMPLITUDES = (0.25, 0.5, 1.0)
RELATIVE_CUTOFF = 1e-3


def layout(obs):
    center = (obs.positions.min(axis=0) + obs.positions.max(axis=0)) / 2
    mask = np.max(np.abs(obs.positions - center), axis=1) <= obs.width / 6
    return mask, mask.astype(float) / np.sum(mask)


def dump_json(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def box_margin(m, bounds):
    lower, upper = bounds
    if np.min(m) < lower or np.max(m) > upper:
        raise ValueError('Archived center falls outside its material box')
    return np.minimum(m - lower, upper - m)


def strongest_near_null(J, margin, task):
    wall, cpu = time.perf_counter(), time.process_time()
    scaled = J * margin
    _, singular, Vh = svd(scaled, full_matrices=False, check_finite=False, lapack_driver='gesdd')
    cutoff = RELATIVE_CUTOFF * singular[0]
    rank = int(np.sum(singular > cutoff))
    target = margin * task
    strong = Vh[:rank]
    weak = target - strong.T @ (strong @ target)
    max_abs = np.max(np.abs(weak))
    direction = None if max_abs == 0 else margin * (weak / max_abs)
    return direction, singular, {
        'material_unknowns': int(len(margin)),
        'real_observation_rows': int(J.shape[0]),
        'relative_singular_cutoff': RELATIVE_CUTOFF,
        'absolute_singular_cutoff': float(cutoff),
        'retained_strong_directions': rank,
        'weak_dimension_at_cutoff': int(len(margin) - rank),
        'positive_symmetric_margin_pixels': int(np.count_nonzero(margin > 0)),
        'roi_pixels_with_positive_symmetric_margin': int(np.count_nonzero(target)),
        'weak_component_norm_over_box_scaled_task_norm': (
            float(np.linalg.norm(weak) / np.linalg.norm(target)) if np.any(target) else None
        ),
        'wall_seconds': time.perf_counter() - wall,
        'cpu_seconds': time.process_time() - cpu,
    }


def pair_directions(targeted, margin, mask, task, amplitude):
    roi = margin * mask
    yield 'projected_task', amplitude * targeted, 'same nominal fraction of symmetric box allowance'
    yield 'direct_roi', amplitude * roi, 'same nominal fraction of symmetric box allowance'
    match = (task @ targeted) * amplitude / (task @ roi)
    if match <= 1:
        yield 'direct_roi_equal_task', match * roi, 'same regional mean gap'


def audited_forward(model, m, label):
    """Enclose the field error for the represented discrete operators in double precision."""
    state = model.forward(m, label=label)
    factor = model.factor
    with model.operation('INVERSE_RESIDUAL_AUDIT', rhs=model.n, inverse_check_RHS=model.n):
        A = model.identity - m[:, None] * model.G
        b = m[:, None] * model.E0
        B = lu_solve(factor, model.identity, check_finite=False)
        measured_inverse_defect = np.linalg.norm(A @ B - model.identity)
        measured_state_defect = np.linalg.norm(A @ state['j'] - b)

        # Frobenius norm dominates the spectral norm. The 32(n+1)eps factor
        # exceeds the complex dot-product error and covers norm evaluation.
        u = np.finfo(np.float64).eps
        dot = 32 * (model.n + 1) * u
        A_error = 32 * u * (np.sqrt(model.n) + np.linalg.norm(m[:, None] * model.G))
        Bnorm, Anorm, jnorm = (np.linalg.norm(value) for value in (B, A, state['j']))
        inverse_defect = 2 * measured_inverse_defect + (dot * Anorm + A_error) * Bnorm
        if not np.isfinite(inverse_defect) or inverse_defect >= 1:
            raise ArithmeticError(f'Unverified discrete inverse for {label}: {inverse_defect:.3g}')
        inverse_bound = 2 * Bnorm / (1 - inverse_defect)
        state_defect = 2 * measured_state_defect + (dot * Anorm + A_error) * jnorm + 32 * u * np.linalg.norm(b)
        Snorm = np.linalg.norm(model.S)
        error = Snorm * (inverse_bound * state_defect + dot * jnorm)

    return state, factor, {
        'field_error_bound': float(error),
        'state_relative_residual': state['state_relative_residual'],
        'inverse_defect_bound': float(inverse_defect),
        'inverse_2norm_bound': float(inverse_bound),
    }


def interval(value, error):
    rounding = 128 * np.finfo(np.float64).eps * abs(value)
    return [float(max(0, value - error - rounding)), float(value + error + rounding)]


def norm_interval(a, b, a_error, b_error=0.):
    # Includes subtraction rounding when fields nearly cancel.
    amount = np.linalg.norm(a - b)
    subtraction_error = 8 * np.finfo(np.float64).eps * (np.linalg.norm(a) + np.linalg.norm(b))
    return interval(amount, a_error + b_error + subtraction_error)


def ball_status(distance, radius):
    if distance[1] <= radius:
        return 'inside'
    if distance[0] > radius:
        return 'outside'
    return 'uncertain'


def qualification(plus, minus, radius, task_gap):
    p, m = (ball_status(distance, radius) for distance in (plus, minus))
    return {
        'plus': p, 'minus': m, 'both_inside': p == m == 'inside',
        'max_task_error_lower_bound': task_gap / 2 if p == m == 'inside' else None,
    }


def thresholds_for_center(obs, state):
    # The 10% synthetic center-ball is an independent thought experiment.
    # Use the actual observed-data budget for the common-center checks.
    return {
        'archived_observation_sphere': obs.radius,
        'same_budget_center_sphere': obs.radius,
        'nominal_10pct_center_sphere': float(.1 * np.linalg.norm(state['y'])),
    }


def evaluate_pair(obs, m, state, state_bound, J, h, label, task, thresholds, arrays):
    lo, hi = obs.bounds
    # Clip at a rounding-only boundary crossing; any substantive crossing is an error.
    plus_material, minus_material = m + h, m - h
    for material in (plus_material, minus_material):
        if np.min(material) < lo - 1e-13 or np.max(material) > hi + 1e-13:
            raise ValueError(f'Illegal endpoint: {label}')
    plus_material = np.clip(plus_material, lo, hi)
    minus_material = np.clip(minus_material, lo, hi)
    plus, factor, pcheck = audited_forward(obs.model, plus_material, label + '_plus')
    minus, _, mcheck = audited_forward(obs.model, minus_material, label + '_minus')
    pe, me = pcheck['field_error_bound'], mcheck['field_error_bound']
    ce = state_bound['field_error_bound']
    py, my, cy = plus['y'], minus['y'], state['y']

    with obs.model.operation('EXACT_SECANT_IDENTITY_AUDIT', rhs=obs.model.ns, secant_check_RHS=obs.model.ns):
        secant_rhs = (plus_material - minus_material)[:, None] * minus['Q']
        identity_field = obs.model.S @ lu_solve(factor, secant_rhs, check_finite=False)
    actual_secant = py - my
    linear = obs.model.pack(actual_secant)
    linearized = J @ (plus_material - minus_material)
    gap = float(abs(task @ (plus_material - minus_material)))
    history_plus = norm_interval(py, obs.Y, pe)
    history_minus = norm_interval(my, obs.Y, me)
    center_plus = norm_interval(py, cy, pe, ce)
    center_minus = norm_interval(my, cy, me, ce)
    separation = norm_interval(py, my, pe, me)
    for k, v in [('m_plus', plus_material), ('m_minus', minus_material),
                 ('Y_plus', py), ('Y_minus', my)]:
        arrays[f'{label}__{k}'] = v
    return {
        'array_prefix': label,
        'regional_chi_plus': float(task @ plus_material),
        'regional_chi_minus': float(task @ minus_material),
        'task_gap': gap,
        'material_rms_between_ends': float(np.linalg.norm(plus_material - minus_material) / np.sqrt(len(m))),
        'material_linf_between_ends': float(np.max(np.abs(plus_material - minus_material))),
        'plus_min_max': [float(np.min(plus_material)), float(np.max(plus_material))],
        'minus_min_max': [float(np.min(minus_material)), float(np.max(minus_material))],
        'field_separation_interval': separation,
        'plus_distance_to_center_interval': center_plus,
        'minus_distance_to_center_interval': center_minus,
        'even_midpoint_drift_interval': norm_interval((py + my) / 2, cy, (pe + me) / 2, ce),
        'plus_distance_to_archive_interval': history_plus,
        'minus_distance_to_archive_interval': history_minus,
        'linearized_field_separation': float(np.linalg.norm(linearized)),
        'nonlinear_discrepancy_from_linearized': float(np.linalg.norm(linear - linearized)),
        'exact_cross_end_secant_formula_error_norm': float(np.linalg.norm(actual_secant - identity_field)),
        'plus_numerics': pcheck,
        'minus_numerics': mcheck,
        'archive_ball': qualification(history_plus, history_minus, thresholds['archived_observation_sphere'], gap),
        'same_budget_center_ball': qualification(center_plus, center_minus, thresholds['same_budget_center_sphere'], gap),
        'nominal_10pct_center_ball': qualification(center_plus, center_minus, thresholds['nominal_10pct_center_sphere'], gap),
        'midpoint_same_budget_overlap_proven': bool(separation[1] <= 2 * thresholds['same_budget_center_sphere']),
        'midpoint_same_budget_max_task_error_lower_bound': (
            gap / 2 if separation[1] <= 2 * thresholds['same_budget_center_sphere'] else None
        ),
    }


def common_material_prior(m, side):
    grid = m.reshape(side, side)
    difference = np.r_[np.diff(grid, axis=0).ravel(), np.diff(grid, axis=1).ravel()] / 19.5
    return 1e-4 * float(np.mean(difference**2) + np.mean(np.sqrt(difference**2 + .01**2) - .01))


def audited_history(obs, q, arrays):
    rows = {}
    for arm, endpoint in obs.endpoints.items():
        state, _, numerics = audited_forward(obs.model, endpoint.m, arm)
        y = state['y']
        gap = norm_interval(y, obs.Y, numerics['field_error_bound'])
        arrays[f'history__{arm}__m'] = endpoint.m
        arrays[f'history__{arm}__Y'] = y
        rows[arm] = {
            'regional_chi': float(q @ endpoint.m),
            'relative_residual': float(np.linalg.norm(y - obs.Y) / np.linalg.norm(obs.Y)),
            'archive_ball': ball_status(gap, obs.radius),
            'distance_to_archive_interval': gap,
            'common_phi': float(.5 * np.linalg.norm(y - obs.Y)**2 / np.linalg.norm(obs.Y)**2 +
                                common_material_prior(endpoint.m, obs.side)),
            'lower_active': int(np.sum(endpoint.m <= obs.bounds[0] + 1e-12)),
            'upper_active': int(np.sum(endpoint.m >= obs.bounds[1] - 1e-12)),
            'replay_saved_field_rel_error': None if endpoint.stored_prediction is None else float(
                np.linalg.norm(y - endpoint.stored_prediction) /
                max(np.linalg.norm(endpoint.stored_prediction), 1e-300)
            ),
            'numerics': numerics,
            'source': endpoint.provenance,
        }
    archived_pairs = []
    eligible = [arm for arm, row in rows.items() if row['archive_ball'] == 'inside']
    for first_index, first in enumerate(eligible):
        for second in eligible[first_index + 1:]:
            separation = abs(rows[first]['regional_chi'] - rows[second]['regional_chi'])
            archived_pairs.append({
                'arms': [first, second], 'task_gap': separation,
                'max_task_error_lower_bound': separation / 2,
            })
    return rows, sorted(archived_pairs, key=lambda row: row['task_gap'], reverse=True)


def disjoint_array_audit(obs, arrays, result):
    """Replay a disjoint aperture of the same already saved synthetic material."""
    half_case, full_case = SRAP_CASES[:2]
    if obs.id != half_case:
        return
    full = srap_observation(full_case)
    assert np.array_equal(obs.positions, full.positions)
    half_meta, full_meta = (record.provenance['acquisition'] for record in (obs, full))
    half_tx, full_tx = (np.array(meta['tx_plane_wave_propagation_deg']) for meta in (half_meta, full_meta))
    half_rx, full_rx = (np.array(meta['rx_m']) for meta in (half_meta, full_meta))
    tx = np.min(np.abs(half_tx[:, None] - full_tx[None, :]), axis=0) > 1e-8
    rx = np.min(np.linalg.norm(half_rx[:, None] - full_rx[None, :], axis=2), axis=0) > 1e-8
    assert np.any(tx) and np.any(rx)
    model = full.model
    model.E0 = model.E0[:, tx]
    model.S = model.S[rx]
    model.ns, model.nr = np.count_nonzero(tx), np.count_nonzero(rx)
    y = full.Y[np.ix_(rx, tx)]
    audited = []
    for row in result['pairs']:
        if row['direction'] != 'projected_task':
            continue
        key = row['array_prefix']
        predictions = []
        for suffix in ('m_plus', 'm_minus'):
            state, _, audit = audited_forward(model, arrays[f'{key}__{suffix}'], f'{key}__heldout_{suffix}')
            predictions.append(state['y'])
            arrays[f'{key}__heldout_{suffix[2:]}_Y'] = state['y']
        audited.append({
            'pair': key,
            'plus_relative_residual': float(np.linalg.norm(predictions[0] - y) / np.linalg.norm(y)),
            'minus_relative_residual': float(np.linalg.norm(predictions[1] - y) / np.linalg.norm(y)),
            'pair_gap_over_heldout_data_norm': float(np.linalg.norm(predictions[0] - predictions[1]) / np.linalg.norm(y)),
        })
    arrays['heldout_Y'] = y
    result['same_material_disjoint_aperture'] = {
        'observation': full_case,
        'transmitter_indices': np.flatnonzero(tx).tolist(),
        'receiver_indices': np.flatnonzero(rx).tolist(),
        'fields': audited,
        'cost': dict(model.counts),
        'operator_cpu_seconds': dict(model.cpu_timings),
    }


def case_run(obs):
    t0, c0 = time.perf_counter(), time.process_time()
    mask, q = layout(obs)
    arrays = {'Y_archive': obs.Y, 'positions': obs.positions, 'task_mask': mask, 'task_vector': q}
    history, archived_pairs = audited_history(obs, q, arrays)
    cost_history = dict(obs.model.counts)
    result = {
        'case': obs.id, 'side': obs.side, 'number_of_real_material_cells': obs.side**2,
        'ports': list(obs.Y.shape), 'domain_width': obs.width, 'material_box': list(obs.bounds),
        'acquisition': obs.provenance, 'archive_norm': float(np.linalg.norm(obs.Y)),
        'archive_radius': float(obs.radius), 'archive_relative_radius': float(obs.radius / np.linalg.norm(obs.Y)),
        'roi_pixel_count': int(np.sum(mask)), 'historical_baselines': history,
        'historical_feasible_pairs': archived_pairs, 'centers': {}, 'pairs': [],
        'cost_historical_endpoint_replays': cost_history,
        'historical_endpoint_replay_cpu_seconds': time.process_time() - c0,
        'historical_endpoint_replay_wall_seconds': time.perf_counter() - t0,
    }
    for arm in obs.centers:
        w, c = time.perf_counter(), time.process_time()
        prior_cost = Counter(obs.model.counts)
        m = obs.endpoints[arm].m
        margin = box_margin(m, obs.bounds)
        state, _, check = audited_forward(obs.model, m, f'{arm}_center')
        J = obs.model.jacobian(state)
        direction, singular, direction_summary = strongest_near_null(J, margin, q)
        thresholds = thresholds_for_center(obs, state)
        arrays[f'{arm}__singular_values'] = singular
        arrays[f'{arm}__symmetric_margin'] = margin
        center_gap = norm_interval(state['y'], obs.Y, check['field_error_bound'])
        result['centers'][arm] = {
            'historical_endpoint': arm, 'numerics': check, 'thresholds': thresholds,
            'scaled_J_svd': direction_summary, 'archive_distance_interval': center_gap,
            'archive_radius_slack_lower_bound': float(obs.radius - center_gap[1]),
            'archive_ball': ball_status(center_gap, obs.radius),
        }
        if direction is None:
            result['centers'][arm].update({
                'stop': 'regional task has no component in the symmetric weak search space',
                'cost': dict(obs.model.counts - prior_cost), 'cpu_seconds': time.process_time() - c,
                'wall_seconds': time.perf_counter() - w,
            })
            continue
        arrays[f'{arm}__direction'] = direction

        center_predict = obs.model.pack(state['y'])
        h = 1e-4 * direction
        state_p = obs.model.forward(m + h, label=f'{arm}_Jcheck_plus')
        state_m = obs.model.forward(m - h, label=f'{arm}_Jcheck_minus')
        difference = obs.model.pack(state_p['y']) - obs.model.pack(state_m['y'])
        tangent = J @ (2 * h)
        result['centers'][arm]['jacobian_central_difference_error_over_center_norm'] = float(
            np.linalg.norm(difference - tangent) / max(np.linalg.norm(center_predict), 1e-300)
        )
        for amplitude in AMPLITUDES:
            for method, h, matched_quantity in pair_directions(direction, margin, mask, q, amplitude):
                prefix = f'{arm}__{method}__{int(100*amplitude):03d}'
                row = evaluate_pair(obs, m, state, check, J, h, prefix, q, thresholds, arrays)
                result['pairs'].append({
                    'center': arm, 'direction': method, 'amplitude': amplitude,
                    'matching_rule': matched_quantity, **row,
                })
        result['centers'][arm].update({
            'cost': dict(obs.model.counts - prior_cost), 'cpu_seconds': time.process_time() - c,
            'wall_seconds': time.perf_counter() - w,
        })

    disjoint_array_audit(obs, arrays, result)
    result['operation_cost'] = dict(obs.model.counts)
    result['operator_wall_seconds'] = dict(obs.model.timings)
    result['operator_cpu_seconds'] = dict(obs.model.cpu_timings)
    result['case_compute_wall_seconds'] = time.perf_counter() - t0
    result['case_compute_cpu_seconds'] = time.process_time() - c0
    return result, arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    t0, c0 = time.perf_counter(), time.process_time()
    names = []
    for obs in observations():
        result, arrays = case_run(obs)
        temp = args.output / f'{obs.id}.tmp.npz'
        np.savez_compressed(temp, **arrays)
        temp.replace(args.output / f'{obs.id}.npz')
        dump_json(args.output / f'{obs.id}.json', result)
        names.append(obs.id)
        print(f'{obs.id}: {len(result["pairs"])} paired end-point rows, {result["case_compute_cpu_seconds"]:.2f} CPU s', flush=True)
    expected = [RING, AUS, *SRAP_CASES]
    if names != expected:
        raise RuntimeError('Incomplete observation cohort')
    dump_json(args.output / 'ENDPOINT_SEAL.json', {
        'cases': names, 'truth_read': False,
        'physics': 'original matching worktree full nonlinear scalar VIE, complex dense solves',
        'roi': 'every cell center inside the public central square of domain-width/3',
        'directions': ['box-scaled relative truncated-SVD projection', 'direct ROI', 'equal task direct ROI if feasible'],
        'alpha': list(AMPLITUDES), 'relative_singular_cutoff': RELATIVE_CUTOFF,
        'execution_wall_seconds_including_serialization': time.perf_counter() - t0,
        'execution_cpu_seconds_including_serialization': time.process_time() - c0,
        'workspace': str(WORKSPACE),
    })


if __name__ == '__main__':
    main()
