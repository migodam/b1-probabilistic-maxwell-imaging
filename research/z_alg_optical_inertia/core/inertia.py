"""Numerical proposals for sector and real-material exclusion certificates."""
from collections import Counter
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import eigh, eigvalsh
from scipy.linalg.lapack import get_lapack_funcs


def adjoint(matrix):
    return matrix.conj().T


def hermitian(matrix):
    return (matrix + adjoint(matrix)) * .5


def clock():
    return perf_counter(), process_time()


def timing(start):
    return dict(zip(('wall_seconds', 'cpu_seconds'),
                    (float(t - t0) for t, t0 in zip(clock(), start))))


def diagonals(n, u, mu, lam):
    """Validate diagonal sector caps, nonnegative sector weights and real equality weights."""
    values = []
    for arg in (u, mu, lam):
        array = np.asarray(arg)
        if not np.isrealobj(array) or array.shape not in ((), (n,)):
            raise ValueError('Each parameter must be a real scalar or a real state-length vector.')
        array = np.broadcast_to(array.astype(float), (n,))
        if not np.all(np.isfinite(array)):
            raise ValueError('Caps and weights must be finite.')
        values.append(array)
    if np.any(values[0] < 0.) or np.any(values[1] < 0.):
        raise ValueError('Sector caps and weights must be nonnegative.')
    return values


class Prepared:
    """Fixed complex observations, with the same reception operator for every source."""

    def __init__(self, G, S, B, Y):
        self.G, self.S, self.B, self.Y = [np.asarray(x, dtype=complex) for x in (G, S, B, Y)]
        self.n, self.s = self.B.shape
        self.r = self.S.shape[0]
        if self.G.shape != (self.n, self.n) or self.S.shape != (self.r, self.n):
            raise ValueError('G and S dimensions must agree with B.')
        if self.Y.shape != (self.r, self.s):
            raise ValueError('Use only common measured rows, and preserve B/Y source order.')
        if min(self.n, self.s, self.r) == 0 or not all(
                np.isfinite(x).all() for x in (self.G, self.S, self.B, self.Y)):
            raise ValueError('Inputs must be nonempty finite matrices.')
        self.Gh = hermitian(self.G)
        self.Ga = (self.G - adjoint(self.G)) / (2j)
        self.S2 = hermitian(adjoint(self.S) @ self.S)
        self.SY = adjoint(self.S) @ self.Y
        self.Y2 = hermitian(adjoint(self.Y) @ self.Y)
        self.operations = Counter(receiver_grams=1, receiver_observation_products=1, observation_grams=1)

    def matrices(self, u, mu, lam=0.):
        u, mu, lam = diagonals(self.n, u, mu, lam)
        z = .5 * (mu * u + 1j * lam)
        K = self.S2 - 2. * hermitian(z[:, None] * self.G)
        K.flat[::self.n + 1] += mu
        return K, self.SY + z[:, None] * self.B

    def evaluate_ldl(self, u, mu, lam=0., *, return_solution=False):
        """Use LDL proposals for spatial search and fall back near unresolved pivots."""
        started = clock()
        K, L = self.matrices(u, mu, lam)
        factor, work_query, solve = get_lapack_funcs(('hetrf', 'hetrf_lwork', 'hetrs'), (K,))
        work, _ = work_query(self.n, lower=True)
        lu, pivots, info = factor(K, lower=True, lwork=int(work.real))
        self.operations['dense_K_LDL_factorizations'] += 1
        if info != 0:
            return self.evaluate(u, mu, lam, return_solution=return_solution)
        d = lu.diagonal().real.copy()
        if np.any(pivots < 0):
            i = 0
            while i < self.n:
                if pivots[i] < 0:
                    diagonal_mean = (d[i] + d[i + 1]) * .5
                    deviation = np.hypot((d[i] - d[i + 1]) * .5, abs(lu[i + 1, i]))
                    d[i], d[i + 1] = diagonal_mean - deviation, diagonal_mean + deviation
                    i += 1
                i += 1
        if np.min(np.abs(d)) <= spectral_tolerance(d):
            return self.evaluate(u, mu, lam, return_solution=return_solution)
        nu = int(np.count_nonzero(d < 0.))
        if nu >= self.s:
            return {**result_row(nu, self.s), 'method': 'LDL_numeric', **timing(started)}
        X, info = solve(lu, pivots, L, lower=True)
        self.operations['LDL_source_rhs'] += self.s
        defect = float(np.linalg.norm(K @ X - L))
        if info != 0 or defect > 1.e-7 * max(np.linalg.norm(L), 1.e-300):
            return self.evaluate(u, mu, lam, return_solution=return_solution)
        self.operations['H_eigendecompositions'] += 1
        row = {**result_row(nu, self.s, self.Y2 - adjoint(L) @ X), 'method': 'LDL_numeric',
               'KX_minus_L_Frobenius': defect, **timing(started)}
        if return_solution:
            row['solution_current'] = X
            row['linear_term'] = L
        return row

    def evaluate(self, u, mu, lam=0., *, return_solution=False):
        """Independent dense numerical proposal. Certificate must rebuild the matrices."""
        started = clock()
        K, L = self.matrices(u, mu, lam)
        d, V = eigh(K, check_finite=False, driver='evd')
        self.operations['dense_K_eigendecompositions'] += 1
        gap = float(np.min(np.abs(d)))
        tolerance = spectral_tolerance(d)
        row = {'K_min_eigenvalue': float(d[0]), 'K_max_eigenvalue': float(d[-1]),
               'K_min_abs_eigenvalue': gap, 'K_singularity_tolerance': tolerance}
        if gap <= tolerance:
            return {**row, 'nu': None, 'status': 'numeric_inertia_unresolved',
                    'lower_bound_numeric': None, **timing(started)}
        nu = int(np.count_nonzero(d < 0.))
        if nu >= self.s:
            return {**result_row(nu, self.s), **row, **timing(started)}
        X = V @ ((adjoint(V) @ L) / d[:, None])
        self.operations['source_rhs'] += self.s
        self.operations['H_eigendecompositions'] += 1
        row.update(result_row(nu, self.s, self.Y2 - adjoint(L) @ X),
                   KX_minus_L_Frobenius=float(np.linalg.norm(K @ X - L)), **timing(started))
        if return_solution:
            row['solution_current'] = X
            row['linear_term'] = L
        return row


def spectral_tolerance(values):
    return float(max(np.finfo(float).tiny, 128. * len(values) * np.finfo(float).eps
                     * np.max(np.abs(values))))


def result_row(nu, s, H=None):
    if H is None:
        eig = None
        bound = 0.
        status = 'negative_directions_exhaust_sources_numeric'
    else:
        eig = eigvalsh(hermitian(H), check_finite=False)[::-1]
        bound = float(np.sqrt(np.maximum(eig[nu:], 0.).sum()))
        status = 'positive_definite_numeric' if nu == 0 else 'indefinite_with_source_directions_numeric'
    return {'nu': nu, 'available_source_directions': max(0, s - nu), 'status': status,
            'H_eigenvalues_descending': None if eig is None else eig.tolist(),
            'lower_bound_numeric': bound}


class SectorBasis:
    """Reuse the eigendecomposition of A in K = mu*A + S* S across positive mu."""

    def __init__(self, G, u, ratio, weight=1., equality_weight=1.):
        self.u, self.ratio, self.weight, self.equality_weight = u, ratio, weight, equality_weight
        n = G.shape[0]
        upper, w, e = diagonals(n, u, weight, equality_weight)
        self.z = .5 * (w * upper + 1j * ratio * e)
        A = -2. * hermitian(self.z[:, None] * G)
        A.flat[::n + 1] += w
        if not np.any(A.imag):
            A = np.array(A.real, order='F')
        self.d, self.V = eigh(A, check_finite=False, driver='evd', overwrite_a=True)
        self.n_negative = int(np.count_nonzero(self.d < 0.))
        self.min_gap = float(np.min(np.abs(self.d)))
        self.usable = self.min_gap > spectral_tolerance(self.d)

    def scan(self, problem, mus):
        """Diagonal solve plus receiver-space Schur update, with a dense fallback."""
        costs = problem.operations
        costs['A_basis_uses'] += 1
        if not self.usable:
            return [self._dense(problem, mu) for mu in mus]
        Vstar = adjoint(self.V)
        CV = Vstar @ adjoint(problem.S)
        DY = CV @ problem.Y
        ZB = Vstar @ (self.z[:, None] * problem.B)
        DA_CV = CV / self.d[:, None]
        receiver_update = hermitian(adjoint(CV) @ DA_CV)
        dr, Vr = eigh(receiver_update, check_finite=False)
        costs['receiver_update_eigendecompositions'] += 1
        rows = []
        for mu in mus:
            row_start = clock()
            if mu <= 0.:
                raise ValueError('The reusable sector basis requires strictly positive mu.')
            eigen_update = 1. + dr / mu
            if np.min(np.abs(eigen_update)) <= max(1.e-10, spectral_tolerance(eigen_update)):
                rows.append(self._dense(problem, mu))
                continue
            nu = self.n_negative - int(np.count_nonzero(eigen_update < 0.))
            if not 0 <= nu <= problem.n:
                rows.append(self._dense(problem, mu))
                continue
            if nu >= problem.s:
                rows.append({**result_row(nu, problem.s), 'method': 'receiver_update_numeric',
                             **timing(row_start)})
                continue
            Lprime = DY + mu * ZB
            Xprime0 = Lprime / (mu * self.d[:, None])
            rec_rhs = adjoint(CV) @ Xprime0
            rec_sol = Vr @ ((adjoint(Vr) @ rec_rhs) / eigen_update[:, None])
            Xprime = Xprime0 - (DA_CV @ rec_sol) / mu
            # Assess the residual in the diagonal basis before accepting the proposal.
            defect = mu * self.d[:, None] * Xprime + CV @ (adjoint(CV) @ Xprime) - Lprime
            rel = float(np.linalg.norm(defect) / max(np.linalg.norm(Lprime), 1.e-300))
            costs['receiver_update_source_rhs'] += problem.s
            if rel > 1.e-7:
                rows.append(self._dense(problem, mu))
                continue
            H = problem.Y2 - adjoint(Lprime) @ Xprime
            costs['H_eigendecompositions'] += 1
            rows.append({**result_row(nu, problem.s, H), 'method': 'receiver_update_numeric',
                         'relative_transformed_solve_defect': rel, **timing(row_start)})
        return rows

    def _dense(self, problem, mu):
        return {**problem.evaluate(self.u, mu * self.weight, mu * self.ratio * self.equality_weight),
                'method': 'dense_fallback_numeric'}
