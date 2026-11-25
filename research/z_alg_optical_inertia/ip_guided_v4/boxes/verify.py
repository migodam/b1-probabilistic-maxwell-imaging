"""Enclose selected translated-sector bounds and centers on the original grids."""
import argparse
from collections import Counter
from pathlib import Path

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import cholesky, eigh, solve
from scipy.sparse.linalg import svds
from threadpoolctl import threadpool_limits

from run import (CONTROL, CONTROL_ID, HERE, SMALL_ID, arrays, ball_parameters, box, clocks,
                 herm, load, read, star, timing, write)


def ball(x):
    return acb_mat(x.shape[0], x.shape[1], [acb(float(v.real), float(v.imag)) for v in x.flat])


def tr(x):
    return x.conjugate().transpose()


def frob2(a):
    return sum((z.real * z.real + z.imag * z.imag for z in a.entries()), arb(0))


def frob(a):
    return frob2(a).sqrt()


def hi(x):
    return float(np.nextafter(float(x.upper()), np.inf))


def lo(x):
    return float(np.nextafter(float(x.lower()), -np.inf))


def midpoint(A):
    return np.reshape([complex(v) for v in A.entries()], (A.nrows(), A.ncols()))


def gersh(a, *, upper, n=None):
    n = a.nrows() if n is None else n
    sign = 1 if upper else -1
    endpoints = []
    for i in range(n):
        off = sum((abs(a[i, j]) for j in range(n) if i != j), arb(0))
        edge = a[i, i].real + sign * off
        endpoints.append(edge.upper() if upper else edge.lower())
    return max(endpoints) if upper else min(endpoints)


def diagonal_times(diag, b):
    n, m = b.nrows(), b.ncols()
    return acb_mat(n, m, [diag[i] * b[i, j] for i in range(n) for j in range(m)])


class ProofGrid:
    def __init__(self, geom):
        self.geom = geom
        self.kernel = None
        if geom.spectrum is not None:
            # The binary64 multiplier is enclosed by a separable Arb inverse DFT.
            spectrum = geom.spectrum
            p = spectrum.shape[0]
            horizontal = [acb.dft([acb(float(z.real), float(z.imag)) for z in row], inverse=True)
                          for row in spectrum]
            vertical = [acb.dft([row[j] for row in horizontal], inverse=True) for j in range(p)]
            self.kernel = [vertical[j][i] for i in range(p) for j in range(p)]

    def restricted(self, ids):
        G, S, B = self.geom.restricted(ids)
        if self.kernel is None:
            exact_G = ball(G)
        else:
            iy, ix = self.geom.index(ids)
            flat = iy * self.geom.spectrum.shape[0] + ix
            exact_G = acb_mat(len(ids), len(ids), [self.kernel[int(i)] for i in flat.flat])
        return (G, S, B), (exact_G, ball(S), ball(B))

    def full_scatter(self, current):
        if self.geom.G is not None:
            return self.geom.G @ current
        n = self.geom.side
        p = self.geom.spectrum.shape[0]
        padded = np.zeros((p, p, current.shape[1]), dtype=complex)
        padded[:n, :n] = current.reshape(n, n, -1)
        transformed = np.fft.fft2(padded, axes=(0, 1))
        result = np.fft.ifft2(self.geom.spectrum[:, :, None] * transformed, axes=(0, 1))
        return result[:n, :n].reshape(n*n, -1)


def stationary_proposal(numeric, Y, m, r, mu, weight):
    G, S, B = numeric
    A, RG, C, P, D, S2, SY, YY = arrays(G, S, B, Y, m, r, weight)
    K = herm(S2 + mu * C)
    L = SY + mu * P
    d, V = eigh(K, check_finite=False, driver='evd')
    nu = int(np.count_nonzero(d < 0))
    if nu >= B.shape[1]:
        raise ArithmeticError('Negative directions exhaust the source dimension.')
    gap = float(.45 * np.min(abs(d)))
    N = V[:, :nu]
    Wc = N * np.sqrt(-2*d[:nu])
    F = cholesky(K + Wc @ star(Wc) - gap * np.eye(len(K)), lower=True, check_finite=False)
    X = V @ ((star(V) @ L) / d[:, None])
    h, Hv = eigh(herm(YY + mu * D - star(L) @ X), check_finite=False)
    return A, RG, C, nu, gap, N, Wc, F, X, h[::-1], Hv[:, ::-1]


def state_check(A, RG, B, exact_A, exact_RG, m_exact, numeric, exact, grid, ids, y, r_obs):
    """Compute E_bg and independent original-grid state and receiver checks."""
    started = clocks()
    _, sensor, _ = numeric
    _, S, exact_B = exact
    n = len(A)
    m = np.array([float(v) for v in m_exact])
    T = solve(A.T, RG.T, check_finite=False).T
    _, sig, vh = svds(T, k=1, which='LM', tol=1.e-11, rng=117)
    probe = solve(A, vh[0].conj(), check_finite=False)
    p = ball(probe[:, None])
    norm_lower = frob(exact_RG * p).lower() / frob(exact_A * p).upper()
    # A verified square solve establishes invertibility on the exact finite model.
    exact_current = exact_A.solve(diagonal_times(m_exact, exact_B))
    center_distance = frob(S * exact_current - ball(y))
    J = np.zeros((grid.geom.side**2, B.shape[1]), dtype=complex)
    J[ids] = midpoint(exact_current)
    field = grid.geom.B + grid.full_scatter(J)
    m_full = np.zeros(len(J))
    m_full[ids] = m
    full_defect = J - m_full[:, None] * field
    port = grid.geom.S @ J
    local_field = B + numeric[0] @ J[ids]
    # Enclose the source-distance radius for the direct Neumann comparison where possible.
    contraction_certified = False
    kappa_upper = None
    triangle = None
    if sig[0] < 1:
        target = float(sig[0]) + 1.e-5 * (1. - float(sig[0]))
        Kc0 = herm(target**2 * (star(A) @ A) - star(RG) @ RG)
        eigen_floor = float(eigh(Kc0, subset_by_index=(0, 0), check_finite=False)[0][0])
        # Check a *downward shifted* positive factor directly.
        shift = eigen_floor * .45
        shifted = ball(cholesky(Kc0 - shift * np.eye(n), lower=True, check_finite=False))
        i_ball = ball(np.eye(n))
        erre = frob(Kc0_to_exact(exact_A, exact_RG, target) - acb(shift) * i_ball
                    - shifted * tr(shifted)).upper()
        if arb(shift) > erre:
            contraction_certified = True
            kappa_upper = target
            # Use exact A^-1 action on an explicit receiver adjoint and center
            # source.  These are finite-state solves on the exact original model.
            sa = tr(tr(exact_A).solve(tr(S)))
            # Outside this branch the uniform triangle-Neumann radius is undefined.
            # The caller adds R*E_bg and evaluates the actual lower bound.
            triangle = (sa, exact_B)
    return dict(T_norm_numeric=float(sig[0]), T_norm_lower=lo(norm_lower),
                center_distance_lower=lo(center_distance), center_distance_upper=hi(center_distance),
                center_distance_over_radius_lower=lo(center_distance / arb(r_obs)),
                center_distance_over_radius_upper=hi(center_distance / arb(r_obs)),
                exact_center_state_solve=True, contraction_proved=contraction_certified,
                T_norm_upper=kappa_upper,
                full_state_defect_numeric=float(np.linalg.norm(full_defect)),
                restricted_vs_full_field_difference_numeric=float(np.linalg.norm(local_field - field[ids])),
                full_numeric_center_forward_fft_columns=B.shape[1] if grid.kernel is not None else 0,
                **timing(started)), J, field, port, exact_current, triangle


def Kc0_to_exact(A, RG, target):
    return acb(arb(float(target))**2) * (tr(A) * A) - tr(RG) * RG


def dark_obstruction(C, S, exact_C, exact_S):
    """Try a strict receiver-null witness to exclude all positive-definite shifts."""
    started = clocks()
    p, n = S.shape
    if p >= n:
        return dict(status='no_guaranteed_receiver_kernel', **timing(started))
    _, _, vh = np.linalg.svd(S, full_matrices=True)
    Z = star(vh[p:])
    dark = herm(star(Z) @ (C @ Z))
    d, candidate = eigh(dark, subset_by_index=(0, 0), check_finite=False)
    result = dict(dark_sector_value_numeric=float(d[0]))
    if d[0] >= 0:
        return {**result, 'status': 'no_negative_trial', **timing(started)}
    prior = ctx.prec
    ctx.prec = max(256, prior)
    try:
        v = ball(Z @ candidate)
        St = tr(exact_S)
        exact_v = v - St * (exact_S * St).solve(exact_S * v)
        value = (tr(exact_v) * exact_C * exact_v)[0, 0].real
        norm = frob(exact_v)
        result.update(dark_sector_value_upper=hi(value), dark_vector_norm_lower=lo(norm),
                      every_mu_positive_indefinite_as_lower_direction=bool(value < 0 and norm > 0),
                      status='negative_verified' if value < 0 and norm > 0 else 'unresolved')
    except (ZeroDivisionError, ValueError) as exc:
        result.update(status='unresolved', error=str(exc))
    finally:
        ctx.prec = prior
    return {**result, **timing(started)}


def spectral_moment(A, RG, exact_A, exact_RG):
    """Rigorous homogeneous-endpoint radius from tr(T^k), if the moment suffices."""
    started = clocks()
    n = len(A)
    powers = (4, 8, 16, 32)
    T = solve(A.T, RG.T, check_finite=False).T
    spectral = float(np.max(np.abs(np.linalg.eigvals(T))))
    if spectral <= 1.01:
        return dict(rho_numeric=spectral, status='strict_divergence_not_proposed', **timing(started))
    best, chosen = 0., None
    t = T
    for k in (2, *powers):
        t = t @ t
        candidate = float((abs(np.trace(t)) / n)**(1./k))
        if candidate > best:
            best, chosen = candidate, k
        if candidate > 1.02:
            break
    if best <= 1.:
        return dict(rho_numeric=spectral, largest_numeric_moment=best,
                    status='moment_numerically_inconclusive', **timing(started))
    Tb = exact_A.transpose().solve(exact_RG.transpose()).transpose()
    k = 1
    while k < chosen:
        Tb = Tb * Tb
        k *= 2
    trace = sum((Tb[i, i] for i in range(n)), acb(0))
    divided_lower = abs(trace).lower() / n
    strict = bool(divided_lower > 1)
    floor = (divided_lower.log() / chosen).exp() if divided_lower > 0 else arb(0)
    return dict(rho_numeric=spectral, moment_power=chosen, rho_lower=lo(floor),
                status='rho_above_one' if strict else 'moment_unresolved',
                full_Arb_products=int(np.log2(chosen)), exact_Arb_T_RHS=n, **timing(started))


def prove(obs, grid, ids, lower, upper, mu, *, weight=None, try_spectral=False, try_dark=False):
    start = clocks()
    n = len(ids)
    if weight is None:
        weight = np.ones(n)
    numeric, exact = grid.restricted(ids)
    G_np, S_np, B_np = numeric
    G, S, B = exact
    m, r = ball_parameters(lower, upper)
    started = clocks()
    A_np, RG_np, C_np, nu, shift, N, Wc, Fc, X, h, Hv = stationary_proposal(
        numeric, obs.Y, m, r, mu, weight)
    proposal = timing(started)
    started = clocks()
    m_exact = [(arb(float(a)) + arb(float(b))) / 2 for a, b in zip(lower, upper, strict=True)]
    r_exact = [(arb(float(b)) - arb(float(a))) / 2 for a, b in zip(lower, upper, strict=True)]
    weight_exact = list(map(lambda x: arb(float(x)), weight))
    Y, I = ball(obs.Y), ball(np.eye(n))
    A = I - diagonal_times(m_exact, G)
    F = diagonal_times(m_exact, B)
    RG, RB = (diagonal_times(r_exact, x) for x in (G, B))
    C = tr(A) * diagonal_times(weight_exact, A) - tr(RG) * diagonal_times(weight_exact, RG)
    P = tr(A) * diagonal_times(weight_exact, F) + tr(RG) * diagonal_times(weight_exact, RB)
    D = tr(F) * diagonal_times(weight_exact, F) - tr(RB) * diagonal_times(weight_exact, RB)
    mu_b = acb(float(mu))
    K, L, YY = tr(S) * S + mu_b * C, tr(S) * Y + mu_b * P, tr(Y) * Y + mu_b * D
    Fc, X, Hv = tuple(map(ball, (Fc, X, Hv)))
    if nu:
        N, Wc = ball(N), ball(Wc)
        K_plus = K + Wc * tr(Wc)
        negative = gersh(tr(N) * (K * N + acb(shift) * N), upper=True)
        if not negative < 0:
            raise ArithmeticError('Arb negative subspace unresolved.')
    else:
        K_plus = K
        negative = None
    cholesky_residual = frob(K_plus - acb(shift) * I - Fc * tr(Fc)).upper()
    alpha = (arb(shift) - cholesky_residual).lower()
    if not alpha > 0:
        raise ArithmeticError('Arb shifted Cholesky remainder could not enclose inertia.')
    inertia = timing(started)
    started = clocks()
    KX = K * X
    defect = frob2(KX - L).upper()
    correction = (defect / alpha).upper()
    H0 = YY - tr(L) * X - tr(X) * L + tr(X) * KX
    projected_H, gram = tr(Hv) * H0 * Hv, tr(Hv) * Hv
    terms = []
    squared = arb(0)
    for j in range(nu + 1, B.ncols() + 1):
        numerator = gersh(projected_H, upper=False, n=j)
        den = gersh(gram, upper=True, n=j)
        term = max(arb(0), (numerator / den - correction).lower()) if numerator > 0 and den > 0 else arb(0)
        squared += term
        terms.append(max(0., lo(term)))
    lower_bound = max(arb(0), squared.lower()).sqrt().lower()
    endpoint = timing(started)
    center, J, E, port, J_exact, triangle = state_check(
        A_np, RG_np, B_np, A, RG, m_exact, numeric, exact, grid, ids, obs.Y, obs.radius)
    if triangle is not None:
        sa, exact_B = triangle
        e_center = exact_B + G * J_exact
        forcing = diagonal_times(r_exact, e_center)
        radius_upper = (frob(sa).upper() * frob(forcing).upper()
                        / (1 - arb(center['T_norm_upper']))).upper()
        exact_distance = frob(S * J_exact - Y)
        neumann_distance = max(arb(0), (exact_distance.lower() - radius_upper).lower())
        center.update(neumann_residual_over_radius_lower=max(0., lo(neumann_distance / arb(obs.radius))),
                      neumann_port_variation_upper=hi(radius_upper))
    dark = dark_obstruction(C_np, S_np, C, S) if try_dark else None
    moment = spectral_moment(A_np, RG_np, A, RG) if try_spectral else None
    record = dict(nu=nu, mu=mu, cells=n, source_dimension=B.ncols(),
                  peak_lower=float(lower.max()), peak_upper=float(upper.max()),
                  verified_K_abs_gap_lower=lo(alpha),
                  shifted_Cholesky_defect_upper=hi(cholesky_residual),
                  restricted_negative_shifted_edge_upper=None if negative is None else hi(negative),
                  stationary_defect_upper=hi(defect.sqrt()), H_scalar_correction_upper=hi(correction),
                  H_numeric_eigenvalues_descending=h.tolist(),
                  H_retained_eigenvalue_lowers=terms,
                  distance_lower=max(0., lo(lower_bound)),
                  distance_over_radius_lower=max(0., lo(lower_bound / arb(obs.radius))),
                  strictly_excludes=bool(squared > arb(obs.radius)**2),
                  neumann_2norm_criterion_fails=bool(center['T_norm_lower'] > 1),
                  center=center, exact_dark_check=dark, homogeneous_endpoint_moment=moment,
                  cost=dict(numeric_proposal=proposal, arb_state_inertia=inertia, arb_endpoint=endpoint,
                            arb_full_state_gram_products=3, exact_center_source_rhs=B.ncols(),
                            **timing(start)))
    return record, dict(active=ids, lower=lower, upper=upper, weight=weight,
                        center_material_active=m, center_current_full=J, center_field_full=E, center_ports=port)


def selection(case, which):
    def eligible(b, positive=None):
        v = [x for x in b['by_mu'] if x['distance_over_radius_numeric'] is not None]
        if positive is not None:
            v = [x for x in v if (x['nu'] == 0) == positive]
        return max(v, default=None, key=lambda x: x['distance_over_radius_numeric'])

    boxes = case['boxes']
    if which == 'control':
        box0 = next(b for b in boxes if b['peak'] is None)
        yield 'control', box0, eligible(box0)
        return
    candidate = [b for b in boxes if b['peak'] == 4 and b['norm_T_numeric'] > 1
                 and eligible(b, False) is not None
                 and eligible(b, False)['distance_over_radius_numeric'] > 1]
    if not candidate:
        raise ArithmeticError('No eligible indefinite high-peak box in the numerical proposals.')
    largest = max(candidate, key=lambda b: (b['width'], b['radius_fraction'],
                                          eligible(b, False)['distance_over_radius_numeric']))
    yield 'rank_tail', largest, eligible(largest, False)
    same = eligible(largest, True)
    if same is not None:
        yield 'same_box_positive', largest, same
    smaller = [b for b in boxes if b['placement'] == largest['placement'] and
               b['radius_fraction'] == largest['radius_fraction'] and b['peak'] == largest['peak']
               and b['norm_T_numeric'] < 1 and eligible(b, True) is not None]
    best_contractive = max(smaller, key=lambda b: b['width'])
    yield 'last_contractive', best_contractive, eligible(best_contractive, True)


def offline_witness(obs, grid, payload, data_path):
    start = clocks()
    with np.load(data_path, allow_pickle=False) as data:
        m = data['material'].ravel()
    physical_support = np.flatnonzero(m)
    _, (G, S, B) = grid.restricted(physical_support)
    ms = [arb(float(m[i])) for i in physical_support]
    A = ball(np.eye(len(ms))) - diagonal_times(ms, G)
    J = A.solve(diagonal_times(ms, B))
    residual = frob(S * J - ball(obs.Y))
    operator_error = bool(residual.upper() <= arb(obs.radius))
    inside_box = None
    if payload:
        active = payload['active']
        inside_box = bool(np.all(m[active] >= payload['lower']) and np.all(m[active] <= payload['upper'])
                          and len(np.setdiff1d(physical_support, active)) == 0)
    return dict(peak=float(np.max(m)), cells=len(physical_support), agrees_with_own_data_ball=operator_error,
                residual_over_radius_upper=hi(residual.upper() / arb(obs.radius)), belongs_to_box=inside_box,
                exact_material_single_solve=True, source_rhs=obs.Y.shape[1], **timing(start))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=HERE / 'results')
    ap.add_argument('--scan', default='scan.json')
    ap.add_argument('--name', default='verified')
    ap.add_argument('--case', action='append')
    ap.add_argument('--spectral', action='store_true')
    ap.add_argument('--dark', action='store_true')
    ap.add_argument('--precision', type=int, default=128)
    ap.add_argument('--threads', type=int, default=2)
    args = ap.parse_args()
    ctx.prec, ctx.threads = args.precision, args.threads
    scan = read(args.out / args.scan)['cases']
    start = clocks()
    geos, rows, saved = {}, {}, {}
    for obs_id, obs in load(args.out, args.case).items():
        case = scan[obs_id]
        begun = clocks()
        if obs.geometry.side not in geos:
            geos[obs.geometry.side] = ProofGrid(obs.geometry)
            setup = timing(begun)
        else:
            setup = None
        grid = geos[obs.geometry.side]
        outputs, attempted, control_payload = {}, [], None
        for label, b, row in selection(case, 'control' if obs_id == CONTROL_ID else 'high'):
            params = {k: b[k] for k in ('placement', 'radius_fraction', 'peak', 'width')}
            params['placement_center'] = b['center']
            if b['peak'] is None:
                ids, lower, upper = box(obs.geometry, b['center'], b['radius_fraction'], .375, .375)
                lower[:], upper[:] = 0., .75
            else:
                ids, lower, upper = box(obs.geometry, b['center'], b['radius_fraction'], b['peak'], b['width'])
            a = clocks()
            try:
                result, payload = prove(obs, grid, ids, lower, upper, row['mu'],
                                        try_spectral=args.spectral and label == 'rank_tail',
                                        try_dark=args.dark and label == 'rank_tail')
                result.update(params, status='outward_verified', factor=row['factor'])
                saved.update({obs_id + '__' + label + '__' + k: v for k, v in payload.items()})
                if label == 'control':
                    control_payload = payload
            except (ArithmeticError, ValueError, ZeroDivisionError, np.linalg.LinAlgError) as exc:
                result = dict(params, status='unresolved', reason=str(exc), **timing(a))
            attempted.append(label)
            outputs[label] = result
            print(obs_id, label, result['status'], result.get('nu'), result.get('distance_over_radius_lower'),
                  result.get('center', {}).get('T_norm_lower'), result.get('reason', ''), flush=True)
        rows[obs_id] = dict(side=obs.geometry.side, radius=obs.radius, setup=setup, tests=outputs,
                            attempts=attempted, **timing(begun))
        if obs_id in (SMALL_ID, CONTROL_ID):
            data = ((args.out / (SMALL_ID + '__offline.npz')) if obs_id == SMALL_ID
                    else (CONTROL / (CONTROL_ID + '.npz')))
            rows[obs_id]['offline_generating_witness'] = offline_witness(obs, grid, control_payload, data)
        write(args.out / (args.name + '.json'),
              dict(cases=rows, arb_precision=args.precision, arb_threads=args.threads, **timing(start)))
        np.savez_compressed(args.out / (args.name + '_witnesses.npz'), **saved)


if __name__ == '__main__':
    with threadpool_limits(1):
        main()
