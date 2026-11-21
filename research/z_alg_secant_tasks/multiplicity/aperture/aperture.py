"""Finite-aperture contrasts: moving components and moving one receiver."""

from collections import Counter
from itertools import combinations, product
from pathlib import Path
import sys
from time import perf_counter, process_time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'continuum'))
from cylinders import BASE, K, RADIUS, SPACING, Cylinders
from separated import encoded, l2, newton, write_json

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from threadpoolctl import threadpool_limits


ORDERS = (8, 12)
COUNTS = (2, 3)
DIRECTIONS = {'x': None, '45deg': (2**-.5, 2**-.5), 'y': (0., 1.)}
RATIOS = (.5, 1., 2., 4.)
OFFSETS = np.array([0., .01, .02, .04, .10, .20, -.01, -.02, -.04, -.10, -.20])
NOISE = .10
REFLECTED_OFFSETS = np.r_[0, np.arange(6, 11), np.arange(1, 6)]


def summarize_layout(roots, model):
    if len(roots) < 2:
        return None
    ports = np.array([port for _, _, port in roots])
    cross = ports[:, ~model.local]
    center = np.mean(ports, axis=0)
    center[model.local] = model.target
    pairwise = [{
        'branches': [name_a, name_b],
        'material_L2': l2(first-second),
        'full_port_Frobenius': l2(a-b),
        'local_port_L2': l2((a-b)[model.local]),
        'cross_port_L2': l2((a-b)[~model.local]),
    } for (name_a, first, a), (name_b, second, b) in combinations(roots, 2)]
    distances = [row['full_port_Frobenius'] for row in pairwise]
    r_lower = max(distances)/2
    r_upper = max(l2(matrix-center) for matrix in ports)
    r_near = min(distances)/2
    reference = l2(center)
    delta = NOISE*reference
    return {
        'full_reference_matrix': encoded(center), 'full_reference_Frobenius': reference,
        'fixed_repeated_local_target_L2': l2(model.target),
        'cross_reference_L2': l2(center[~model.local]),
        'cross_L2_by_branch': {root[0]: l2(cross[idx]) for idx, root in enumerate(roots)},
        'radius_of_covering_all_roots_bounds': [r_lower, r_upper],
        'nearest_pair_covering_radius': r_near,
        'covering_radius_fraction_full_reference_bounds': [r_lower/reference, r_upper/reference],
        'nearest_pair_radius_fraction_full_reference': r_near/reference,
        'radius_10_percent_same_full_reference': delta,
        'all_roots_in_fixed_reference_ball': r_upper <= delta,
        'all_roots_cannot_share_any_ball': r_lower > delta,
        'nearest_pair_has_intersecting_balls': r_near <= delta,
        'pairwise': pairwise,
    }


def layout_trials(order, saved):
    one = Cylinders(1, 0, order)
    one.target = one.solve(BASE)[0][one.local]
    partner, _, _, attempt = newton(one, BASE[::-1])
    seeds = {'A': BASE}
    if attempt['status'] == 'numerical_root' and l2(BASE-partner) > 1e-5:
        seeds['B'] = partner
    cases = []
    for count, (label, direction) in product(COUNTS, DIRECTIONS.items()):
        model = Cylinders(count, SPACING, order, direction=direction)
        model.target = np.repeat(one.target, count)
        key = f'layout_M{order}_L{count}_{label}'
        saved[f'{key}__centers'] = model.centers
        saved[f'{key}__sources'] = model.sources
        saved[f'{key}__receivers'] = model.receivers
        roots, attempts = [], {}
        for branch in map(''.join, product(seeds, repeat=count)):
            initial = np.concatenate([seeds[tag] for tag in branch])
            material, matrix, _, details = newton(model, initial)
            status = details['status'] == 'numerical_root' and all(l2(material-r[1]) > 1e-5 for r in roots)
            details['accepted_as_distinct'] = status
            attempts[branch] = details
            if status:
                roots.append((branch, material, matrix))
                saved[f'{key}_{branch}__material'] = material
                saved[f'{key}_{branch}__ports'] = matrix
        cases.append({'direction': label, 'components': count, 'key': key,
                      'cross_indices_receiver_and_source': np.argwhere(~model.local).tolist(),
                      'root_count': len(roots), 'expected_root_count': 2**count,
                      'attempts': attempts, 'cost': dict(model.cost), 'readings': summarize_layout(roots, model)})
    for case in cases:
        x_layout = next(c for c in cases if c['components'] == case['components'] and c['direction'] == 'x')
        anchor = x_layout['readings']['full_reference_Frobenius'] if x_layout['readings'] else None
        row = case['readings']
        if row is not None and anchor is not None:
            row['shared_x_layout_reference_Frobenius'] = anchor
            row['covering_radius_fraction_shared_x_reference_bounds'] = [v/anchor for v in row['radius_of_covering_all_roots_bounds']]
            row['nearest_pair_radius_fraction_shared_x_reference'] = row['nearest_pair_covering_radius']/anchor
    return {'order': order, 'component_target_port': encoded(one.target.reshape(1, 1)),
            'component_partner_attempt': attempt, 'single_component_cost': dict(one.cost), 'cases': cases}


def mirror_ports(model, material, receiver_map, source_amplitude):
    tau, _ = model.transmission(material, differentiate=False)
    state = model.identity-tau[:, None]*model.G
    incident = tau*model.B[:, 0]*source_amplitude
    coefficients = lu_solve(lu_factor(state, check_finite=False), incident, check_finite=False)
    model.cost.update({'state_LU': 1, 'source_RHS': 1})
    return receiver_map@coefficients


def symmetric_receiver_pair(readings, axis_reference):
    result = []
    for h, positive, negative in zip(OFFSETS[1:6], readings[:, 1:6].T, readings[:, 6:].T, strict=True):
        z = np.column_stack((positive, negative))
        even = (z+z[:, ::-1])/2
        odd = (z-z[:, ::-1])/2
        even_norm = np.linalg.norm(even, axis=1)
        odd_norm = np.linalg.norm(odd, axis=1)
        mixed_norm = np.minimum(even_norm, odd_norm)
        delta = NOISE*l2(z[0])
        axis_delta = NOISE*(2**.5)*axis_reference
        result.append({
            'positive_offset': float(h), 'z_by_branch_positive_negative': encoded(z),
            'even_by_branch': encoded(even), 'odd_by_branch': encoded(odd),
            'even_L2_by_branch': even_norm.tolist(), 'odd_L2_by_branch': odd_norm.tolist(),
            'clean_reference_noise_radius': delta,
            'fixed_repeated_axis_reference_noise_radius': axis_delta,
            'odd_over_clean_radius_for_calibrated_gain': float(odd_norm[0]/delta),
            'min_even_odd_over_clean_radius_for_unknown_gain': float(mixed_norm[0]/delta),
            'calibrated_odd_margin_above_radius': float(odd_norm[0]-delta),
            'unknown_gain_min_even_odd_margin_above_radius': float(mixed_norm[0]-delta),
            'calibrated_odd_over_fixed_axis_radius': float(odd_norm[0]/axis_delta),
            'unknown_gain_min_even_odd_over_fixed_axis_radius': float(mixed_norm[0]/axis_delta),
            'calibrated_test_rules_out_joint_feasibility_of_mirror_pair': bool(odd_norm[0] > delta),
            'unknown_gain_test_rules_out_joint_feasibility_of_mirror_pair': bool(mixed_norm[0] > delta),
        })
    return result


def mirror_trials(order, saved):
    cases = []
    for count, ratio in product(COUNTS, RATIOS):
        model = Cylinders(count, SPACING, order, rx_x=0., k=K*ratio, direction=DIRECTIONS['y'])
        points = np.c_[OFFSETS, np.full(len(OFFSETS), .95)]
        new_map = model.outgoing(points).reshape(len(OFFSETS), -1)
        receiver_map = np.vstack([model.S[1:], new_map])
        scale = ratio**.5
        key = f'mirror_M{order}_L{count}_k{ratio:g}'
        saved[f'{key}__centers'] = model.centers
        saved[f'{key}__fixed_source'] = model.sources[0]
        saved[f'{key}__old_receivers'] = model.receivers[1:]
        saved[f'{key}__new_receiver_cases'] = points
        pairs = []
        for tags in map(''.join, product('AB', repeat=count-1)):
            label = 'A'+tags
            reverse = label.translate(str.maketrans('AB', 'BA'))
            material = np.concatenate([BASE if tag == 'A' else BASE[::-1] for tag in label])
            reflected = material.reshape(-1, 2)[:, ::-1].ravel()
            ports = np.array([mirror_ports(model, mat, receiver_map, scale) for mat in (material, reflected)])
            old, new = ports[:, :count-1], ports[:, count-1:]
            old_gap = l2(old[0]-old[1])
            gap = np.abs(new[0]-new[1])
            axis_mean = np.mean(new[:, 0])
            marginal_ref = float(abs(axis_mean))
            joint_ref = l2(np.r_[np.mean(old, axis=0), axis_mean])
            joint_gap = np.hypot(old_gap, gap)
            odd = (new[0]-new[1])[1:6] - (new[0]-new[1])[6:]
            derivative = odd/(2*OFFSETS[1:6])
            avg_gap = (gap[1:6]+gap[6:])/2
            slope = float(np.polyfit(np.log(OFFSETS[1:4]), np.log(avg_gap[:3]), 1)[0])
            pairs.append({
                'branches': [label, reverse], 'material_L2': l2(material-reflected),
                'fixed_old_axis_gap_L2': old_gap, 'axis_marginal_reference_magnitude': marginal_ref,
                'axis_joint_reference_L2': joint_ref, 'extra_port_gap': gap.tolist(),
                'marginal_half_gap_over_fixed_10_percent_radius': (gap/(2*NOISE*marginal_ref)).tolist(),
                'joint_half_gap_over_fixed_10_percent_radius': (joint_gap/(2*NOISE*joint_ref)).tolist(),
                'joint_pair_can_be_excluded_at_fixed_10_percent_noise': (joint_gap > 2*NOISE*joint_ref).tolist(),
                'reflected_opposite_offset_L2_field_mismatch': l2(new[0]-new[1, REFLECTED_OFFSETS]),
                'first_three_offset_log_gap_slope': slope,
                'odd_secant_derivative_complex_by_positive_offset': encoded(derivative.reshape(1, -1))[0],
                'symmetric_receiver_pair': symmetric_receiver_pair(new, marginal_ref),
            })
            saved[f'{key}_{label}__material_pair'] = np.stack([material, reflected])
            saved[f'{key}_{label}__old_and_new_port_cases'] = ports
        cases.append({'components': count, 'frequency_ratio': ratio, 'source_amplitude_relative_to_k0': scale,
                      'key': key, 'pairs': pairs, 'cost': dict(model.cost)})
    for case in cases:
        anchor = next(c for c in cases if c['components'] == case['components'] and c['frequency_ratio'] == 1.)
        for row, at_k0 in zip(case['pairs'], anchor['pairs'], strict=True):
            old = row['fixed_old_axis_gap_L2']
            gap = np.hypot(old, row['extra_port_gap'])
            row['joint_half_gap_over_fixed_k0_10_percent_radius'] = (gap/(2*NOISE*at_k0['axis_joint_reference_L2'])).tolist()
    return {'order': order, 'cases': cases}


def truncation(layouts, mirrors, saved):
    result = {'layouts': [], 'mirrors': []}
    for coarse, fine in zip(layouts[0]['cases'], layouts[1]['cases'], strict=True):
        p, q = coarse['key'], fine['key']
        shared = [tag for tag in fine['attempts'] if f'{q}_{tag}__material' in saved and f'{p}_{tag}__material' in saved]
        result['layouts'].append({
            'key': q, 'common_roots': len(shared),
            'max_material_L2_change': max((l2(saved[f'{p}_{tag}__material']-saved[f'{q}_{tag}__material']) for tag in shared), default=None),
            'max_full_port_change': max((l2(saved[f'{p}_{tag}__ports']-saved[f'{q}_{tag}__ports']) for tag in shared), default=None),
        })
    for coarse, fine in zip(mirrors[0]['cases'], mirrors[1]['cases'], strict=True):
        p, q = coarse['key'], fine['key']
        keys = [f"{row['branches'][0]}__old_and_new_port_cases" for row in fine['pairs']]
        result['mirrors'].append({'key': q, 'max_full_port_change': max(l2(saved[f'{p}_{tag}']-saved[f'{q}_{tag}']) for tag in keys)})
    return result


def total_cost(rows):
    total = Counter()
    for row in rows:
        total.update(row.get('single_component_cost', {}))
        for case in row['cases']:
            total.update(case['cost'])
    return dict(total)


def main():
    saved = {}
    begin, cpu = perf_counter(), process_time()
    layouts = [layout_trials(order, saved) for order in ORDERS]
    timing = {'layout_wall_seconds': perf_counter()-begin, 'layout_cpu_seconds': process_time()-cpu}
    begin, cpu = perf_counter(), process_time()
    mirrors = [mirror_trials(order, saved) for order in ORDERS]
    timing.update(mirror_wall_seconds=perf_counter()-begin, mirror_cpu_seconds=process_time()-cpu)
    result = {
        'scope': 'Numerical Fourier–Bessel angular truncations, scattered fields, no rigorous tail bounds',
        'units': 'reference vacuum wavelength, lengths held fixed across frequency',
        'disk_radius': RADIUS, 'disk_centers_relative': [[-.18, 0], [.18, 0]],
        'component_spacing': SPACING, 'frequency_base_k': K, 'orders': list(ORDERS),
        'layout_receiver_offset_x': .013, 'aperture_receiver_offset_x_grid': OFFSETS.tolist(),
        'noise_fraction': NOISE, 'layout_trials': layouts, 'mirror_trials': mirrors,
        'angular_order_comparison': truncation(layouts, mirrors, saved),
        'layout_cost': total_cost(layouts), 'mirror_cost': total_cost(mirrors), 'timing': timing,
    }
    np.savez_compressed(HERE/'materials_and_ports.npz', **saved)
    write_json(HERE/'results.json', result)


if __name__ == '__main__':
    with threadpool_limits(limits=1, user_api='blas'):
        main()
