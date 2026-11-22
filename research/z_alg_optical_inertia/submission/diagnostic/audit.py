#!/usr/bin/env python3
"""Independently enclose stored GN endpoints and score their materials offline."""
import argparse
from pathlib import Path
import sys

from run import CASES, HIGH, HERE, NATIVE, RING, STUDY, clocks, elapsed, load_case, read_json, write_json

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from flint import acb, acb_mat, arb, ctx

sys.path.insert(0, str(STUDY))
from core.certificate import ball, frobenius_squared
from synthetic.nonvacuity import ExactFourierGeometry, StableActiveModel, hi, lo

R10 = STUDY.parents[1] / 'research/delegated/max_struct_r1/r10_worktree/research/research_max_struct_r10'
RING_PROOF = HIGH / 'sector_inertia/certificate_u5.json'
NATIVE_PROOF = STUDY / 'results/theory/fine_feasibility.json'


def norm(x):
    return frobenius_squared(x).sqrt()


def extrema(m):
    exact_material = [arb(float(value)) for value in m]
    return dict(minimum_lo=lo(min(exact_material)), maximum_hi=hi(max(exact_material)),
                minimum_nonnegative=bool(min(exact_material) >= 0))


def distances(corrected, error, target, online):
    distance = norm(ball(target) - corrected)
    online_error = norm(ball(online) - corrected).upper() + error
    return dict(full_wave_distance_lo=max(0., lo(distance.lower() - error)),
                full_wave_distance_hi=hi(distance.upper() + error),
                online_forward_error_hi=hi(online_error))


def dense_audit(model, m, y, online):
    """Use an interval Neumann defect on the matched 24² active system."""
    started = clocks()
    active = np.flatnonzero(m)
    values = m[active]
    G = model.G[np.ix_(active, active)]
    B = model.E0[active]
    S = model.S[:, active]
    n = len(active)
    A = -values[:, None] * G
    A.flat[::n+1] += 1.
    lu = lu_factor(A)
    X = lu_solve(lu, values[:, None] * B)
    W = lu_solve(lu, S.T, trans=1).T
    R = lu_solve(lu, np.eye(n, dtype=complex))
    M = [acb(float(value)) for value in values]
    G, B, S, X, W, R = map(ball, (G, B, S, X, W, R))
    I = acb_mat(n, n)
    for i in range(n):
        I[i, i] = 1
    A = I - acb_mat(n, n, [M[i] * G[i, j] for i in range(n) for j in range(n)])
    rhs = acb_mat(n, B.ncols(), [M[i] * B[i, j] for i in range(n) for j in range(B.ncols())])
    inverse_defect = norm(I - R*A).upper()
    if not inverse_defect < 1:
        raise ArithmeticError('Arb inverse defect did not establish invertibility.')
    inv_bound = (norm(R).upper()/(1 - inverse_defect)).upper()
    primal = rhs - A*X
    dual = S - W*A
    err = (norm(dual)*inv_bound*norm(primal)).upper()
    corrected = S*X + W*primal
    return dict(active_cells=n, state_inverse_defect_hi=hi(inverse_defect), state_inverse_norm_hi=hi(inv_bound),
                corrected_port_remainder_hi=hi(err), inverse_rhs=n, source_rhs=B.ncols(), receiver_rhs=S.nrows(),
                **distances(corrected, err, y, online), **elapsed(started))


def fourier_audit(geometry, m, y, online):
    """Reuse the independently rounded active-cell exact-Fourier replay."""
    started = clocks()
    model = StableActiveModel(geometry, m)
    inv_bound, inverse = model.inverse_bound()
    corrected, error, audit = model.enclose_field(inv_bound, 128)
    return dict(active_cells=model.n, inverse=inverse, state_and_dual=audit,
                inverse_rhs=model.n, source_rhs=geometry.B.shape[1], receiver_rhs=geometry.S.shape[0],
                **distances(corrected, error, y, online), **elapsed(started))


def reference_and_metrics(case, m, ring_reference):
    if case == RING:
        with np.load(ring_reference, allow_pickle=False) as raw:
            truth = raw['m']
        # R10 has a low host and several high-contrast ring components.
        inside = truth > .1
    else:
        truth_path = STUDY / 'synthetic/evaluation/optfresh_01.npz'
        with np.load(truth_path, allow_pickle=False) as raw:
            truth = raw['m_fine']
        inside = truth > 0.
    error = m - truth
    return truth, dict(reference_source=('authorized original R10 synthetic reference' if case == RING
                                         else str(truth_path.relative_to(STUDY.parents[1]))),
                       truth_peak=float(truth.max()),
                       global_relative_l2=float(np.linalg.norm(error) / np.linalg.norm(truth)),
                       interior_relative_l2=float(np.linalg.norm(error[inside]) / np.linalg.norm(truth[inside])),
                       background_rmse=float(np.sqrt(np.mean(error[~inside]**2))),
                       interior_pixels=int(np.count_nonzero(inside)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=CASES)
    parser.add_argument('--r10-reference', type=Path,
                        help='Optional original synthetic reference, available only in an authorized workspace.')
    args = parser.parse_args()
    ctx.prec, ctx.threads = 128, 1

    case_output = HERE / 'interval_and_offline.json'
    report = read_json(case_output) if case_output.exists() else {}
    for case in CASES:
        if args.case and case != args.case:
            continue
        start = clocks()
        online = read_json(HERE / (case + '_online.json'))
        turn = online['accepted_indices'][-1]
        with np.load(HERE / (case + '_path.npz'), allow_pickle=False) as sealed:
            m, y, y_online = (sealed[key] for key in ('materials', 'Y', 'predictions'))
        model, measurement = load_case(case)
        if not np.array_equal(y, measurement):
            raise ValueError('Sealed observation differs from the full-wave model input.')

        if case == RING:
            proof = read_json(RING_PROOF)
            radius = proof['radii']['frozen_observation_ball']
            cap, floor = proof['u'], proof['certified_frobenius_residual_lower']
            with np.load(R10 / 'data_public/twod' / (case + '.npz'), allow_pickle=False) as archived:
                if not np.array_equal(y, archived['Y']):
                    raise ValueError('Original R10 noise draw differs from the shared public archive.')
            records = {'terminal': dense_audit(model, m[turn], y, y_online[turn])}
            proof_file = RING_PROOF
        else:
            proof = next(row for row in read_json(NATIVE_PROOF)['checked_certificate_margins']
                         if row['id'] == case and row['tested_cap'] == .5)
            radius, cap, floor = (proof[key] for key in ('inflated_radius', 'tested_cap', 'distance_lower'))
            geometry = ExactFourierGeometry(STUDY / 'synthetic/public/fine_geometry.npz')
            # Enclose precisely the aperture used in this observation.
            receivers, sources = geometry.inds['PW180_10x20']
            geometry.S = geometry.S[receivers]
            geometry.B = geometry.B[:, sources]
            records = {'terminal': fourier_audit(geometry, m[turn], y, y_online[turn])}
            print(case, 'terminal', records['terminal'], flush=True)
            # Include the chronological first through third updates, across the
            # transition. The primary result is always the frozen terminal iterate.
            for k in (1, 2, 3):
                idx = online['accepted_indices'][k]
                records[f'iterate_{k}'] = fourier_audit(geometry, m[idx], y, y_online[idx])
                records[f'iterate_{k}'].update(extrema(m[idx]))
                print(case, f'iterate_{k}', records[f'iterate_{k}'], flush=True)
            proof_file = NATIVE_PROOF
        records['terminal'].update(extrema(m[turn]))
        for row in records.values():
            row.update(in_cap=row['minimum_nonnegative'] and row['maximum_hi'] <= cap,
                       point_rejected=row['full_wave_distance_lo'] > radius)
        audit_cost = elapsed(start)
        offline = clocks()
        if case == RING and args.r10_reference is None:
            metrics = dict(report.get(case, {}).get('offline_reference', {}))
            metrics['score_origin'] = ('archived scalar scores of the frozen reconstruction'
                                       if metrics else 'reference scores unavailable')
        else:
            reference, metrics = reference_and_metrics(case, m[turn], args.r10_reference)
            metrics['score_origin'] = 'current offline scoring'
            if case == RING:
                state = model.forward(reference)
                metrics['reference_interval'] = dense_audit(model, reference, y, state['y'])
            else:
                feasibility = next(row for row in read_json(NATIVE_PROOF)['observations'] if row['id'] == case)
                metrics['reference_interval'] = dict(full_wave_distance_lo=feasibility['true_forward_to_Y_lower'],
                                                     full_wave_distance_hi=feasibility['true_forward_to_Y_upper'])
                np.savez_compressed(HERE / (case + '_offline_reference.npz'), m=reference)
            metrics['reference_inside_ball'] = metrics['reference_interval']['full_wave_distance_hi'] <= radius
        report[case] = dict(proof=str(proof_file.relative_to(STUDY.parent)), cap=cap, radius=radius,
                            cap_wide_residual_floor=floor, endpoint_idx=turn, intervals=records,
                            audit_setup_and_work=audit_cost, offline_reference=metrics, offline_cost=elapsed(offline))
        write_json(case_output, report)
        print(case, 'completed', report[case], flush=True)


if __name__ == '__main__':
    main()
