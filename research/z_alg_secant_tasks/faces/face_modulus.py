"""Test exact opposing faces of the frozen material specification boxes."""
import json
import os
from pathlib import Path
import sys
from collections import Counter
from time import perf_counter, process_time

for var in ('OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[var] = '1'

import numpy as np
from flint import arb, acb, acb_mat, ctx
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import lsq_linear

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bounds'))
sys.path.insert(0, str(ROOT / 'array'))
from box_modulus import ball_array, from_mat, eye
from array_secant import pack, split_observation

ctx.prec = 128
ARMS = ('base', 'angle', 'frequency', 'rich')
WIDTHS = (.05, .25, .75)
BOTH_SIGNS = (-1, 1)


def interval(value):
    return [float(np.nextafter(float(value.lower()), -np.inf)),
            float(np.nextafter(float(value.upper()), np.inf))]


def c_norm(vector):
    return sum(abs(v) ** 2 for v in vector).sqrt()


def acb_channels(responses, arm):
    base_rx, base_tx = [0, 1, 0], [0, 0, 1]
    first, second = responses
    base = first[base_rx, base_tx]
    energy = 1 / arb(2).sqrt()
    if arm == 'base':
        return base
    if arm == 'angle':
        return energy * np.r_[base, first[base_rx, [2, 2, 3]]]
    if arm == 'frequency':
        return energy * np.r_[base, (arb(6) / 5).sqrt() * second[base_rx, base_tx]]
    return energy * first.T.ravel()


def classification(fields, center, radius):
    distances = [c_norm(field - center) for field in fields]
    separation = c_norm(fields[1] - fields[0])
    # min(δ-d_-, δ-d_+) encloses the same margin without comparing overlapping balls.
    margins = [radius - distance for distance in distances]
    margin_low = min(interval(entry)[0] for entry in margins)
    margin_high = min(interval(entry)[1] for entry in margins)
    return {
        'distances_to_center': [interval(d) for d in distances],
        'end_to_end_distance': interval(separation),
        'center_ball_both': bool(all(value > 0 for value in margins)),
        'some_common_midpoint_ball': bool(separation < 2 * radius),
        'no_common_ball_for_these_two': bool(separation > 2 * radius),
        'remaining_radius_about_F0': [margin_low, margin_high],
        'remaining_radius_about_midpoint': interval(radius - separation / 2),
    }


class FrozenVIE:
    def __init__(self, archive, n, cost):
        self.n, self.cost = n, cost
        stem = f'disks{n:02d}__'
        host = archive[stem + 'known_host_chi']
        self.host = host
        self.host_exact = np.array([arb(float(v)) for v in host], dtype=object)
        self.bands = [[archive[f'{stem}{band}_{key}'] for key in ('G', 'S', 'B')] for band in ('f0', 'f1')]
        self.exact_bands = [[ball_array(matrix) for matrix in band] for band in self.bands]
        self.identity = np.eye(n + len(host))
        self.identity_exact = eye(len(self.identity))

    def certify(self, exact_material, phase):
        m = np.r_[exact_material, self.host_exact]
        response = []
        for G, S, B in self.exact_bands:
            A = acb_mat((self.identity_exact - m[:, None] * G).tolist())
            b = acb_mat((m[:, None] * B).tolist())
            j = from_mat(A.solve(b))
            self.cost[f'{phase}_Arb_solve'] += 1
            self.cost[f'{phase}_Arb_RHS'] += B.shape[1]
            response.append(S @ j)
        return response

    def fit_state(self, material):
        m = np.r_[material, self.host]
        responses, jacobians = [], []
        for G, S, B in self.bands:
            a = self.identity - m[:, None] * G
            factor = lu_factor(a, check_finite=False)
            j = lu_solve(factor, m[:, None] * B, check_finite=False)
            t = lu_solve(factor, S.T, trans=1, check_finite=False).T[:, :self.n]
            q = (B + G @ j)[:self.n]
            responses.append(S @ j)
            jacobians.append(np.einsum('ri,is->rsi', t, q))
            self.cost['optional_GN_float_LU'] += 1
            self.cost['optional_GN_float_source_RHS'] += B.shape[1]
            self.cost['optional_GN_float_receiver_RHS'] += S.shape[0]
        return split_observation(*responses), split_observation(*jacobians)


def inside_binary(value, center, lower, upper):
    # Move each representable end inward once before rounding the trial material.
    return np.clip(value, np.nextafter(lower, center), np.nextafter(upper, center))


def exact_face(outside, exact_center, exact_radius, mask, sign):
    material = np.array([arb(float(v)) for v in outside], dtype=object)
    material[mask] = exact_center[mask] + sign * exact_radius[mask]
    for value, center, radius in zip(material, exact_center, exact_radius, strict=True):
        if value < 0.5 or value > 30 or not (value >= center - radius and value <= center + radius):
            raise ArithmeticError('A reconstructed face coefficient violates the exact material box')
    return material


def bounded_fit(action, target, lower, upper):
    answer = lsq_linear(action, target, bounds=(lower, upper), method='bvls',
                        tol=1e-12, max_iter=200, verbose=0)
    return answer.x, {'status': int(answer.status), 'iterations': int(answer.nit),
                      'linear_residual': float(np.linalg.norm(answer.fun)),
                      'optimality': float(answer.optimality)}


def run_case(case, radius_lookup, archive, rows, materials, cost):
    name, n = case['id'], case['n']
    model = FrozenVIE(archive, n, cost)
    center = np.array(case['material_center'], dtype=float)
    mask = np.array([i in case['roi_indices'] for i in range(n)], dtype=bool)
    if not np.array_equal(center, archive[f'{name}__material_center']):
        raise ArithmeticError('The frozen center does not match its stored material array')
    if not np.array_equal(mask, archive[f'{name}__roi_mask']):
        raise ArithmeticError('The frozen ROI does not match its stored index mask')
    exact_center = np.array([arb(float(v)) for v in center], dtype=object)
    raw_center = model.certify(exact_center, 'center')
    f0 = {arm: acb_channels(raw_center, arm) for arm in ARMS}
    center_discrepancy = {
        arm: float(np.max(np.abs(np.array([complex(v) for v in f0[arm]]) - archive[f'{name}__{arm}__center'])))
        for arm in ARMS
    }
    if max(center_discrepancy.values()) >= 1e-12:
        raise ArithmeticError('The channel ordering or energy scaling no longer matches the saved center field')
    materials[f'{name}__center'] = center
    materials[f'{name}__ROI'] = mask
    for width in WIDTHS:
        exact_radius = np.array([min(arb(float(width)), c - arb(.5), arb(30) - c)
                                 for c in exact_center], dtype=object)
        r = np.array([float(v) for v in exact_radius])
        lo, hi = center - r, center + r
        task_prior_diameter = 2 * sum(exact_radius[mask]) / mask.sum().item()
        direct_material = [exact_face(center, exact_center, exact_radius, mask, sign) for sign in BOTH_SIGNS]
        direct_fields = [model.certify(m, 'roi_only') for m in direct_material]
        for arm in ARMS:
            t, c = perf_counter(), process_time()
            rad = arb(radius_lookup[name, arm])
            observed_center = f0[arm]
            prefix = f'{name}__{width}__{arm}'
            control = classification([acb_channels(field, arm) for field in direct_fields], observed_center, rad)
            j = archive[f'{name}__{arm}__J']
            displacement, linear_diagnostics = bounded_fit(j[:, ~mask], -(j[:, mask] @ r[mask]), -r[~mask], r[~mask])
            cost['outside_only_linear_fits'] += 1
            compensated = []
            for sign in BOTH_SIGNS:
                trial = center.copy()
                trial[~mask] = inside_binary(center[~mask] + sign * displacement, center[~mask], lo[~mask], hi[~mask])
                compensated.append(trial)
            cheq = [exact_face(m, exact_center, exact_radius, mask, sign)
                    for m, sign in zip(compensated, BOTH_SIGNS, strict=True)]
            fields = [acb_channels(model.certify(m, 'linear_compensated'), arm) for m in cheq]
            linear = classification(fields, observed_center, rad)
            additional = None
            selected = linear
            if not linear['center_ball_both']:
                # Frozen fallback is one GN correction per side. Outside coefficients
                # alone move, using the original box and reaching the same ROI face.
                refined_fields, gn_diagnostics = [], []
                for m, sign in zip(compensated, BOTH_SIGNS, strict=True):
                    point = m.copy()
                    point[mask] = center[mask] + sign * r[mask]
                    predicted, jacobians = model.fit_state(point)
                    center_float = np.array([complex(v) for v in observed_center])
                    free = point[~mask]
                    update, details = bounded_fit(pack(jacobians[arm])[:, ~mask],
                                                  pack(center_float - predicted[arm]), lo[~mask] - free, hi[~mask] - free)
                    candidate = point.copy()
                    candidate[~mask] = inside_binary(free + update, center[~mask], lo[~mask], hi[~mask])
                    materials[f'{prefix}__GN_outside_{sign:+}'] = candidate[~mask]
                    arb_material = exact_face(candidate, exact_center, exact_radius, mask, sign)
                    refined_fields.append(acb_channels(model.certify(arb_material, 'additional_GN'), arm))
                    gn_diagnostics.append(details)
                    cost['optional_GN_outside_fits'] += 1
                additional = classification(refined_fields, observed_center, rad)
                additional['linearized_solves'] = gn_diagnostics
                if additional['remaining_radius_about_F0'][0] > linear['remaining_radius_about_F0'][0]:
                    selected = additional
            for point, sign in zip(compensated, BOTH_SIGNS, strict=True):
                materials[f'{prefix}__LS_outside_{sign:+}'] = point[~mask]
            hypotheses = [('roi_only', control), ('outside_linear', linear), ('outside_GN', additional)]
            best_center = max(((method, row) for method, row in hypotheses if row is not None),
                              key=lambda entry: entry[1]['remaining_radius_about_F0'][0])
            best_midpoint = max(((method, row) for method, row in hypotheses if row is not None),
                                key=lambda entry: entry[1]['remaining_radius_about_midpoint'][0])
            rows.append({
                'case': name, 'width': width, 'arm': arm, 'box_roi_task_prior_width': interval(task_prior_diameter),
                'delta': float(rad), 'roi_only': control, 'linear': linear, 'GN': additional,
                'center_replay_max_abs_discrepancy_from_saved_binary64': center_discrepancy[arm],
                'linearized_fit': linear_diagnostics,
                'best_center_ball': best_center[0], 'best_midpoint_ball': best_midpoint[0],
                'center_ball_both': best_center[1]['center_ball_both'],
                'some_common_midpoint_ball': best_midpoint[1]['some_common_midpoint_ball'],
                'roi_only_no_common_but_compensated_center_witness': bool(
                    control['no_common_ball_for_these_two'] and selected['center_ball_both']),
                'best_remaining_radius_about_F0': best_center[1]['remaining_radius_about_F0'],
                'best_remaining_radius_about_midpoint': best_midpoint[1]['remaining_radius_about_midpoint'],
                'row_wall_s': perf_counter() - t, 'row_CPU_s': process_time() - c,
            })


def main():
    start, cpu = perf_counter(), process_time()
    data = json.loads((ROOT / 'array/results/summary.json').read_text())
    radius = {}
    for row in data['pairs']:
        key = row['case'], row['arm']
        if key in radius and radius[key] != row['radius']:
            raise ArithmeticError('Nonconstant archived radius for a frozen case and acquisition')
        radius[key] = row['radius']
    rows, parameters, cost = [], {}, Counter()
    with np.load(ROOT / 'array/results/physical_arrays.npz', allow_pickle=False) as archive:
        for case in data['cases']:
            run_case(case, radius, archive, rows, parameters, cost)
            print(case['id'], 'completed rows:', len(rows), flush=True)
    out = ROOT / 'faces'
    assert len(rows) == 72
    np.savez_compressed(out / 'face_parameters.npz', **parameters)
    common = lambda group: sum(group[k]['center_ball_both'] for k in ('roi_only', 'linear')) > 0
    summary = {
        'center_ball_both': sum(row['center_ball_both'] for row in rows),
        'midpoint_only': sum(row['some_common_midpoint_ball'] and not row['center_ball_both'] for row in rows),
        'no_current_witness': sum(not row['some_common_midpoint_ball'] for row in rows),
        'roi_only_center_ball_both': sum(row['roi_only']['center_ball_both'] for row in rows),
        'at_most_one_linear_fit_center_ball_both': sum(common(row) for row in rows),
        'roi_only_no_common_but_compensated_center_witness': sum(row['roi_only_no_common_but_compensated_center_witness'] for row in rows),
    }
    result = {'metadata': {
        'precision_bits': ctx.prec,
        'scope': 'frozen binary64 G/S/B/host/centers/radii, Arb forward VIE, exact energy scaling',
        'target': 'mean of frozen left ROI; every ROI coefficient fixed at its exact opposing spec-box face',
        'control': 'full-wave same envelope with all outside pixels held exactly at the center',
        'linear': 'one symmetric bounded least-squares solve for outside coefficients on the saved real J',
        'optional_GN': 'only if the linear pair misses central ball, one outside-only fixed-task GN per sign',
        'outside_parameters': 'stored absolute binary64 material coefficients; reconstruct each ROI slot in Arb as c±exact clipped r',
        'wall_s': perf_counter() - start, 'CPU_s': process_time() - cpu,
        'counts': dict(cost),
    }, 'summary': summary, 'rows': rows}
    (out / 'results.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(summary)


if __name__ == '__main__':
    main()
