"""Separate the shifted interval effect from support and check a nonzero fixed host."""
import argparse
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.linalg import solve
from threadpoolctl import threadpool_limits

from run import CONTROL_ID, HERE, NATIVE_IDS, ball_parameters, box, clocks, coordinates, load,     read, scan_box, timing, write
from verify import (ProofGrid, ball, diagonal_times, frob, hi, offline_witness, prove,
                    selection)
from run import CONTROL
from flint import arb, ctx


def matched_box(obs, params):
    return box(obs.geometry, params['center'], params['radius_fraction'], params['peak'], params['width'])


def shadow_weight(obs, ids, low, high):
    G, _, B = obs.geometry.restricted(ids)
    m, _ = ball_parameters(low, high)
    field = solve(np.eye(len(ids)) - G * m[None, :], B, check_finite=False)
    power = np.sum(abs(field)**2, axis=1)
    return np.clip(np.median(power) / np.maximum(power, 1.e-30), .25, 4.)


def candidates(obs, scan, arm):
    n = obs.geometry.side
    if arm == 'broad':
        _, candidate, _ = next(selection(scan, 'high'))
        ids, lower, upper = matched_box(obs, candidate)
        params = {k: candidate[k] for k in ('center', 'peak', 'width', 'radius_fraction', 'placement')}
        return [(params, ids, np.zeros(len(lower)), upper, None)]
    if arm == 'disk':
        candidate = next(b for b in scan['boxes'] if b['peak'] is None)
        shapes = []
        for width in (1., 1.5, 2., 2.5, 3.):
            params = dict(center=candidate['center'], radius_fraction=candidate['radius_fraction'],
                          placement=candidate['placement'], width=width, peak=4.)
            ids, lower, upper = matched_box(obs, params)
            shapes.append((params, ids, lower, upper, None))
        return shapes
    if arm == 'host':
        candidate = next(b for b in scan['boxes'] if b['placement'] == 'born_peak')
        center = np.array(candidate['center'])
        xy = coordinates(n)
        shapes = []
        for name, sign in [('host_left', -1), ('host_right', 1)]:
            host_center = center + np.array([.42*sign, 0.])
            host = np.flatnonzero(np.sum((xy - host_center)**2, axis=1) <= .105**2 + 1.e-15)
            for width in (1., 1.5, 2., 2.5):
                params = dict(center=center.tolist(), radius_fraction=.245, placement='born_peak',
                              width=width, peak=4., host=name, host_center=host_center.tolist(),
                              host_radius=.105, host_fixed_contrast=3.25, host_cells=len(host))
                p, lower, upper = matched_box(obs, params)
                if np.intersect1d(p, host).size or not len(host):
                    raise ArithmeticError('The fixed host must be nonempty and disjoint from the slot.')
                ids = np.r_[p, host]
                lower, upper = np.r_[lower, np.full(len(host), 3.25)], np.r_[upper, np.full(len(host), 3.25)]
                shapes.append((dict(params, metric='identity'), ids, lower, upper, None))
                weight = shadow_weight(obs, ids, lower, upper)
                shapes.append((dict(params, metric='inverse_center_illumination'), ids, lower, upper, weight))
        return shapes
    raise ValueError(arm)


def choose(rows, arm):
    if arm == 'broad':
        yield 'broad', 0, max(rows[0]['by_mu'], key=lambda r: r['distance_over_radius_numeric'] or 0.)
        return
    acceptable = [(i, r) for i, b in enumerate(rows) for r in b['by_mu']
                  if b['norm_T_numeric'] > 1 and 0 < r['nu'] and (r['distance_over_radius_numeric'] or 0.) > 1]
    if not acceptable:
        return
    for metric in ('identity', 'inverse_center_illumination') if arm == 'host' else ('disk',):
        subset = [(i, r) for i, r in acceptable if rows[i]['parameters'].get('metric', 'disk') == metric]
        if not subset:
            continue
        i, r = max(subset, key=lambda v: (rows[v[0]]['parameters']['width'], v[1]['distance_over_radius_numeric']))
        yield 'rank_tail_' + metric, i, r


def host_check(obs, grid, ids, lower, upper, center_port):
    """Replay the exact fixed-host elimination and compare with the original solve."""
    start = clocks()
    local = lower < upper
    k = int(np.count_nonzero(local))
    G_np, _, B_np = obs.geometry.restricted(ids)
    m, _ = ball_parameters(lower, upper)
    A = np.eye(len(ids)) - m[:, None] * G_np
    J = solve(A, m[:, None] * B_np, check_finite=False)
    _, exact = grid.restricted(ids)
    G, S, B = exact
    from flint import acb_mat
    def sub(a, r0, r1, c0, c1):
        return acb_mat(r1-r0, c1-c0, [a[i, j] for i in range(r0, r1) for j in range(c0, c1)])
    n, s = B.nrows(), B.ncols()
    p = S.nrows()
    Gpp, Ghh, Gph, Ghp = (sub(G, *v) for v in ((0, k, 0, k), (k,n,k,n),(0,k,k,n),(k,n,0,k)))
    Sp, Sh = sub(S, 0, p, 0, k), sub(S, 0, p, k, n)
    Bp, Bh = sub(B, 0, k, 0, s), sub(B, k, n, 0, s)
    m_h = [arb(float(v)) for v in lower[k:]]
    m_p = [(arb(float(a)) + arb(float(b)))/2 for a, b in zip(lower[:k], upper[:k], strict=True)]
    Ah = ball(np.eye(n-k)) - diagonal_times(m_h, Ghh)
    Jh0 = Ah.solve(diagonal_times(m_h, Bh))
    Jhp = Ah.solve(diagonal_times(m_h, Ghp))
    G_eff = Gpp + Gph * Jhp
    B_eff = Bp + Gph * Jh0
    S_eff = Sp + Sh * Jhp
    y_host = Sh * Jh0
    A_eff = ball(np.eye(k)) - diagonal_times(m_p, G_eff)
    Jp = A_eff.solve(diagonal_times(m_p, B_eff))
    port = y_host + S_eff * Jp
    comparison = frob(port - ball(center_port))
    return dict(exact_host_invertibility=True, host_cells=n-k, slot_cells=k,
                exact_host_rhs=s+k, exact_slot_source_rhs=s,
                port_vs_lifted_center_forward_upper=hi(comparison),
                original_system_center_residual_numeric=float(np.linalg.norm(A @ J - m[:, None] * B_np)),
                **timing(start))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('arm', choices=('broad', 'disk', 'host'))
    ap.add_argument('--out', type=Path, default=HERE / 'results')
    ap.add_argument('--scan-only', action='store_true')
    ap.add_argument('--certify-from-scan', action='store_true')
    ap.add_argument('--case', action='append')
    args = ap.parse_args()
    ctx.prec, ctx.threads = 128, 2
    start = clocks()
    path = args.out / ('followup_' + args.arm)
    scan = read(args.out / 'scan.json')['cases']
    if args.arm == 'disk':
        wanted = [CONTROL_ID]
    elif args.arm == 'host':
        wanted = list(NATIVE_IDS)
    else:
        wanted = [v for v in scan if v != CONTROL_ID]
    obs_set = load(args.out, args.case or wanted)
    receipts = read(path.with_suffix('.json'))['cases'] if args.certify_from_scan else {}
    for name, obs in obs_set.items():
        if args.certify_from_scan:
            continue
        started = clocks()
        plans = candidates(obs, scan[name], args.arm)
        count, rows = Counter(), []
        for params, ids, lower, upper, weight in plans:
            block = scan_box(obs, ids, lower, upper, count, W=weight)
            rows.append(dict(parameters=params, **block))
            strongest = max((r['distance_over_radius_numeric'] or 0. for r in block['by_mu'] if r['nu'] > 0),
                            default=0.)
            print(args.arm, name, params, 'non-contraction', block['norm_T_numeric'],
                  'indef', strongest, flush=True)
        receipts[name] = dict(proposals=rows, operations=dict(count), scan_cost=timing(started))
        write(path.with_suffix('.json'), dict(cases=receipts, scan_only=True, **timing(start)))
    if args.scan_only:
        return
    grids, saved = {}, {}
    for name, obs in obs_set.items():
        started = clocks()
        plans = candidates(obs, scan[name], args.arm)
        if obs.geometry.side not in grids:
            grids[obs.geometry.side] = ProofGrid(obs.geometry)
            receipts[name]['strict_geometry_cost'] = timing(started)
        grid = grids[obs.geometry.side]
        checks = {}
        for label, i, row in choose(receipts[name]['proposals'], args.arm):
            params, ids, lower, upper, weight = plans[i]
            attempt = clocks()
            try:
                result, payload = prove(obs, grid, ids, lower, upper, row['mu'],
                                        weight=weight, try_dark=args.arm == 'disk')
                result.update(status='outward_verified', parameters=params)
                if args.arm == 'host':
                    result['fixed_host_elimination'] = host_check(
                        obs, grid, ids, lower, upper, payload['center_ports'])
                if args.arm == 'disk':
                    broad = dict(active=ids, lower=np.zeros(len(ids)), upper=np.full(len(ids), .75))
                    result['same_support_low_physical_witness'] = offline_witness(
                        obs, grid, broad, CONTROL / (CONTROL_ID + '.npz'))
                saved.update({name + '__' + label + '__' + k: v for k, v in payload.items()})
            except (ArithmeticError, ValueError, ZeroDivisionError, np.linalg.LinAlgError) as exc:
                result = dict(status='unresolved', parameters=params, reason=str(exc), **timing(attempt))
            checks[label] = result
            print(args.arm, name, label, result['status'], result.get('nu'), result.get('distance_over_radius_lower'),
                  result.get('reason', ''), flush=True)
        receipts[name]['strict_checks'] = checks
        receipts[name]['strict_cost'] = timing(started)
        write(path.with_suffix('.json'), dict(cases=receipts, scan_only=False, **timing(start)))
        np.savez_compressed(path.parent / (path.name + '_witnesses.npz'), **saved)


if __name__ == '__main__':
    with threadpool_limits(1):
        main()
