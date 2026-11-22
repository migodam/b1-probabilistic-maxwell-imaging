"""Outward-rounded Arb proof on the exact saved binary64 finite model."""
import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import cholesky, eigh

try:
    from .inertia import adjoint, clock, diagonals, hermitian, timing
except ImportError:
    from inertia import adjoint, clock, diagonals, hermitian, timing


def ball(array):
    n, m = array.shape
    return acb_mat(n, m, [acb(float(x.real), float(x.imag)) for x in array.flat])


def star(a):
    return a.conjugate().transpose()


def frobenius_squared(a):
    return sum((v.real * v.real + v.imag * v.imag for v in a.entries()), arb(0))


def gershgorin(a, *, upper):
    sign = 1 if upper else -1
    endpoints = []
    for i in range(a.nrows()):
        off = sum((abs(a[i, j]) for j in range(a.ncols()) if i != j), arb(0))
        edge = a[i, i].real + sign * off
        endpoints.append(edge.upper() if upper else edge.lower())
    return max(endpoints) if upper else min(endpoints)


def safe_float(a, *, upper):
    value = a.upper() if upper else a.lower()
    return float(np.nextafter(float(value), np.inf if upper else -np.inf))


def _arb_inputs(problem, u, mu, lam):
    n = problem.n
    u, mu, lam = diagonals(n, u, mu, lam)
    G, S, B, Y = map(ball, (problem.G, problem.S, problem.B, problem.Y))
    m = [acb(float(x)) for x in mu]
    z = [acb(arb(float(a)) * arb(float(b)), arb(float(c))) / 2
         for a, b, c in zip(u, mu, lam, strict=True)]
    if all(np.all(a == a[0]) for a in (u, mu, lam)):
        ZG, ZB = G * z[0], B * z[0]
    else:
        ZG = acb_mat(n, n, [z[i] * G[i, j] for i in range(n) for j in range(n)])
        ZB = acb_mat(n, problem.s, [z[i] * B[i, j] for i in range(n) for j in range(problem.s)])
    K = star(S) * S - ZG - star(ZG)
    for i, value in enumerate(m):
        K[i, i] += value
    return K, star(S) * Y + ZB, star(Y) * Y


def certify(problem, u, mu, lam=0., radii=None, precision=96, threads=2):
    """Verify exact finite-model inertia and a lower residual bound.

    Any benchmark confidence/noise assumption is a property of the supplied radii.
    Numeric proposals are disposable and enter only through verified residual identities.
    """
    started = clock()
    radii = {} if radii is None else {k: float(v) for k, v in radii.items()}
    if not all(np.isfinite(value) and value >= 0. for value in radii.values()):
        raise ValueError('Radii must be nonnegative and finite.')
    n, s = problem.n, problem.s
    try:
        start = clock()
        K, L = problem.matrices(u, mu, lam)
        d, V = eigh(K, check_finite=False, driver='evd')
        nu = int(np.count_nonzero(d < 0.))
        sigma = float(np.min(np.abs(d)) * .5)
        if sigma <= 0. or not np.isfinite(sigma):
            raise ArithmeticError('Unable to propose a strictly positive gap.')
        N = V[:, :nu]
        W = N * np.sqrt(-2. * d[:nu])
        pos = K + W @ adjoint(W)
        pos.flat[::n + 1] -= sigma
        F = cholesky(pos, lower=True, check_finite=False)
        if nu < s:
            X = V @ ((adjoint(V) @ L) / d[:, None])
            hd, hv = eigh(hermitian(problem.Y2 - adjoint(L) @ X), check_finite=False)
            hd, hv = hd[::-1], hv[:, ::-1]
        proposal = timing(start)

        old_precision, old_threads = ctx.prec, ctx.threads
        ctx.prec, ctx.threads = precision, threads
        try:
            start = clock()
            K, L, Y2 = _arb_inputs(problem, u, mu, lam)
            I = acb_mat(n, n)
            for i in range(n):
                I[i, i] = 1
            F = ball(F)
            residual = K - acb(sigma) * I - F * star(F)
            if nu:
                W = ball(W)
                residual += W * star(W)
            positive_residual = frobenius_squared(residual).sqrt().upper()
            alpha = (arb(sigma) - positive_residual).lower()
            if not alpha > 0.:
                raise ArithmeticError('Shifted Cholesky residual did not verify the positive gap.')
            negative_upper = None
            if nu:
                N = ball(N)
                negative_upper = gershgorin(star(N) * (K * N + acb(sigma) * N), upper=True)
                if not negative_upper < 0.:
                    raise ArithmeticError('Projected negative block did not verify the negative gap.')
            inertia_time = timing(start)

            start = clock()
            correction = arb(0)
            defect = arb(0)
            trace_floor = arb(0)
            terms = []
            individual_floor = []
            coherent_ceiling = arb(0)
            individual_ceiling = arb(0)
            if nu < s:
                X = ball(X)
                KX = K * X
                defect = frobenius_squared(KX - L).upper()
                correction = (defect / alpha).upper()
                H0 = Y2 - star(L) * X - star(X) * L + star(X) * KX
                if nu == 0:
                    individual_floor = [max(arb(0), (H0[j, j].real - correction).lower())
                                        for j in range(s)]
                    individual_ceiling = max(arb(0), *(H0[j, j].real.upper() for j in range(s)))
                    V = ball(hv)
                    gram_floor = gershgorin(star(V) * V, upper=False)
                    if gram_floor > 0.:
                        top_ceiling = gershgorin(star(V) * H0 * V, upper=True)
                        coherent_ceiling = max(arb(0), top_ceiling) / gram_floor
                    else:
                        coherent_ceiling = max(arb(0), gershgorin(H0, upper=True))
                for j in range(nu + 1, s + 1):
                    if hd[j - 1] <= 0.:
                        break
                    Vj = ball(hv[:, :j])
                    numerator_floor = gershgorin(star(Vj) * H0 * Vj, upper=False)
                    denominator_ceiling = gershgorin(star(Vj) * Vj, upper=True)
                    eigen_floor = arb(0)
                    if numerator_floor > 0 and denominator_ceiling > 0:
                        eigen_floor = max(arb(0), (numerator_floor / denominator_ceiling - correction).lower())
                    trace_floor += eigen_floor
                    terms.append({'index_descending_one_based': j,
                                  'lower': max(0., safe_float(eigen_floor, upper=False))})
            trace_floor = max(arb(0), trace_floor.lower())
            lower = max(0., safe_float(trace_floor.sqrt(), upper=False))
            return {
                'status': 'Arb_verified', 'precision_bits': precision, 'Arb_threads': threads,
                'verified_negative_eigenvalues': nu,
                'binary64_proposed_shift': sigma,
                'positive_shifted_cholesky_defect_Frobenius_upper': safe_float(positive_residual, upper=True),
                'negative_subspace_shifted_gershgorin_upper': (
                    None if negative_upper is None else safe_float(negative_upper, upper=True)),
                'certified_K_absolute_spectral_gap_lower': safe_float(alpha, upper=False),
                'KX_minus_L_Frobenius_upper': safe_float(defect.sqrt(), upper=True),
                'H_operator_correction_upper': safe_float(correction, upper=True),
                'retained_H_eigenvalue_lower_bounds': terms,
                'certified_frobenius_residual_lower': lower,
                'strongest_coherent_source_lower': (max(0., safe_float(arb(terms[0]['lower']).sqrt(), upper=False))
                                                    if nu == 0 and terms else 0.),
                'strongest_original_source_lower': (
                    max(0., safe_float(max(individual_floor).sqrt(), upper=False)) if individual_floor else 0.),
                'original_source_lower': [max(0., safe_float(x.sqrt(), upper=False)) for x in individual_floor],
                'strongest_coherent_source_upper_same_multipliers': (
                    safe_float(coherent_ceiling.sqrt(), upper=True) if coherent_ceiling > 0. else 0.),
                'strongest_original_source_upper_same_multipliers': (
                    safe_float(individual_ceiling.sqrt(), upper=True) if individual_ceiling > 0. else 0.),
                'radii': radii,
                'strictly_excludes': {name: bool(trace_floor > arb(radius)**2) for name, radius in radii.items()},
                'margin_lower': {name: safe_float(arb(lower) - arb(radius), upper=False)
                                 for name, radius in radii.items()},
                'cost': {'double_proposal': proposal, 'Arb_inertia': inertia_time,
                         'Arb_H': timing(start), 'full_dimension_cubic_Arb_products': 1,
                         'negative_rank_update_Arb_products': int(nu > 0),
                         'double_K_eigendecompositions': 1, 'double_shifted_cholesky': 1,
                         'double_H_eigendecompositions': int(nu < s), **timing(started)},
            }
        finally:
            ctx.prec, ctx.threads = old_precision, old_threads
    except (ArithmeticError, ValueError, np.linalg.LinAlgError) as error:
        return {'status': 'Arb_unresolved', 'reason': str(error), 'cost': timing(started)}
