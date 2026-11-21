"""Exclude upper material caps using the scalar VIE sector and old fields."""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from replay import OBSERVATION_RADII, SOURCE, load_operators, read_json

import numpy as np
from scipy.linalg import LinAlgError, cho_solve, cholesky, eigh, eigvalsh, svdvals
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits


CAPS = (.1, .5, 1., 2., 3., 4.)
WEIGHTS = (1e-4, 1e-3, 1e-2, 1e-1, 1., 10.)
CHECKS = {'r10_twod_0001': (2., .01), 'austria_eps4_400MHz': (.1, 1.)}
PROJECTED_CHECKS = {'r10_twod_0001': (4., .0001), 'austria_eps4_400MHz': (.5, .1),
                    's2d_007__pw360_10x20': (.5, .01)}


def adjoint(x):
    return x.conj().T


def norm_squared(x):
    return float(np.vdot(x, x).real)


def trust_amplitude(eigenvalues, components, ellipsoid_radius_squared):
    """Upper-bound a quadratic over the Euclidean ball via its scalar dual."""
    weights = eigenvalues*np.sum(abs(components)**2, axis=1)
    largest = eigenvalues.max()
    lower = largest + 1e-10*max(largest, 1.)
    upper = lower + max(largest, 1.)

    def derivative(t):
        return ellipsoid_radius_squared-np.sum(weights/(t-eigenvalues)**2)

    multiplier = lower
    if derivative(lower) < 0.:
        while derivative(upper) < 0.:
            upper = largest + 2*(upper-largest)
        multiplier = brentq(derivative, lower, upper)
    cap_squared = (norm_squared(components) + multiplier*ellipsoid_radius_squared
                   + np.sum(weights/(multiplier-eigenvalues)))
    return float(np.sqrt(cap_squared)), float(multiplier)


class Sector:
    def __init__(self, case, observation, physical, radii, spectral_cache):
        self.case, self.Y, self.radii = case, observation, radii
        self.G, self.S, self.B = (physical[field] for field in ('G', 'S', 'E0'))
        self.I = np.eye(self.G.shape[0])
        self.y_squared = norm_squared(self.Y)
        self.Y2 = adjoint(self.Y)@self.Y
        self.H = (self.G+adjoint(self.G))/2
        self.S2 = adjoint(self.S)@self.S
        self.SY = adjoint(self.S)@self.Y
        self.Snorm = float(svdvals(self.S)[0])
        self.Bnorm = float(np.linalg.norm(self.B))
        self.cost = Counter(receiver_svd=1)
        dim = self.G.shape[0]
        if dim not in spectral_cache:
            spectral_cache[dim] = (self.G, float(svdvals(self.G)[0]), *eigh(self.H))
            self.cost.update(green_svd=1, green_hermitian_eigendecomposition=1)
        stored_g, self.Gnorm, self.h, v = spectral_cache[dim]
        if not np.array_equal(self.G, stored_g):
            raise ValueError('The same state dimension is expected to carry the same frozen Green matrix.')
        self.b_mode, self.s_mode = adjoint(v)@self.B, self.S@v

    def hit(self, distance):
        return {key: bool(distance > radius) for key, radius in self.radii.items()}

    def magnitude_controls(self, u):
        y_norm = np.sqrt(self.y_squared)
        contraction = 1-u*self.Gnorm
        direct = {'small_gain_margin': contraction, 'status': 'noncontractive'}
        if contraction > 0.:
            cap = self.Snorm*u*self.Bnorm/contraction
            d = max(0., y_norm-cap)
            direct.update(status='usable', receiver_amplitude_cap=cap,
                          distance_lower_bound=d, excludes=self.hit(d))

        spectrum = 1-u*self.h
        sector = {'C_lambda_min': float(spectrum.min()), 'status': 'indefinite_no_finite_amplitude_bound'}
        if spectrum.min() > 0.:
            center = (u/2)*self.b_mode/spectrum[:, None]
            energy_radius_squared = float(np.sum(spectrum[:, None]*abs(center)**2))
            max_current, lam_j = trust_amplitude(1/spectrum, center, energy_radius_squared)
            covariance = (self.s_mode/spectrum)@adjoint(self.s_mode)
            eigenvalues, vectors = eigh(covariance)
            self.cost['receiver_sector_eigendecomposition'] += 1
            readout_center = adjoint(vectors)@(self.s_mode@center)
            max_receiver, lam_y = trust_amplitude(np.maximum(eigenvalues, 0.), readout_center, energy_radius_squared)
            d_current = max(0., y_norm-self.Snorm*max_current)
            d_receiver = max(0., y_norm-max_receiver)
            sector.update(status='usable', current_amplitude_cap=max_current,
                          current_trust_multiplier=lam_j, receiver_amplitude_cap=max_receiver,
                          receiver_trust_multiplier=lam_y, ellipsoid_radius_squared=energy_radius_squared,
                          distance_lower_using_current_amplitude=d_current,
                          distance_lower_using_receiver_amplitude=d_receiver,
                          current_amplitude_excludes=self.hit(d_current), receiver_amplitude_excludes=self.hit(d_receiver))
        return {'traditional_neumann': direct, 'direction_free_sector_amplitude': sector}

    def quadratic(self, u, mu):
        return self.S2+mu*(self.I-u*self.H), self.SY+mu*u*self.B/2

    def numerical(self, u, mu, gram_matrices):
        begin, cpu = perf_counter(), process_time()
        K, L = self.quadratic(u, mu)
        row = {'u': u, 'mu': mu}
        self.cost['sector_cholesky_attempt'] += 1
        try:
            factored = cholesky(K, lower=True, check_finite=False)
        except LinAlgError:
            row.update(status='no_positive_definite_cholesky')
        else:
            self.cost.update(sector_cholesky_success=1, sector_forward_RHS=L.shape[1])
            candidate = cho_solve((factored, True), L, check_finite=False)
            gram = self.Y2-adjoint(L)@candidate
            gram = (gram+adjoint(gram))/2
            spectrum = eigvalsh(gram, check_finite=False)
            value = float(gram.trace().real)
            bound = np.sqrt(max(0., value))
            positive_squared = float(np.maximum(spectrum, 0.).sum())
            positive_bound = np.sqrt(positive_squared)
            coherent_squared = float(max(0., spectrum[-1]))
            selected_columns_squared = float(np.maximum(gram.diagonal().real, 0.).sum())
            key = f'{self.case}__u_{u:g}__mu_{mu:g}'
            gram_matrices[key] = gram
            self.cost['illumination_gram_eigenvalues'] += 1
            row.update(status='numerical_only', objective_bound_squared_unclipped=value,
                       distance_lower_bound=bound, excludes=self.hit(bound),
                       trace_bound_squared_clipped=max(0., value),
                       positive_gram_bound_squared=positive_squared, positive_gram_distance_lower_bound=positive_bound,
                       best_single_coherent_squared=coherent_squared,
                       best_single_coherent_distance_lower_bound=np.sqrt(coherent_squared),
                       best_single_coherent_excludes=self.hit(np.sqrt(coherent_squared)),
                       positive_original_columns_squared=selected_columns_squared,
                       positive_original_columns_distance_lower_bound=np.sqrt(selected_columns_squared),
                       positive_original_columns_excludes=self.hit(np.sqrt(selected_columns_squared)),
                       positive_gram_excludes=self.hit(positive_bound), gram_array_key=key,
                       gram_eigenvalues=spectrum.tolist(),
                       residual_of_KX_minus_L_Frobenius=float(np.linalg.norm(K@candidate-L)))
        row.update(wall_seconds=perf_counter()-begin, cpu_seconds=process_time()-cpu)
        return row


def ball_matrix(array, acb, acb_mat):
    rows, columns = array.shape
    return acb_mat(rows, columns, [acb(float(z.real), float(z.imag)) for z in array.flat])


def arb_norm_squared(matrix, arb):
    return sum((z.real*z.real+z.imag*z.imag for z in matrix.entries()), arb(0))


def audit(sector, u, mu, project_positive_gram=False):
    """With Arb, show K>=alpha I and correct the quadratic for the solve defect."""
    from flint import acb, acb_mat, arb, ctx
    ctx.prec, ctx.threads = 96, 1
    begin, cpu = perf_counter(), process_time()
    nominal_u, nominal_mu = u, mu
    K, L = sector.quadratic(u, mu)
    n, source_rhs = L.shape
    observable_minimum = eigvalsh(K, subset_by_index=(0, 0), check_finite=False)[0]
    shift = float(observable_minimum/2)
    if shift <= 0.:
        raise ArithmeticError('Could not select a positive Cholesky shift for the requested audit.')
    factor = cholesky(K-shift*sector.I, lower=True, check_finite=False)
    X = cho_solve((cholesky(K, lower=True, check_finite=False), True), L, check_finite=False)
    positive_modes = None
    if project_positive_gram:
        gram = sector.Y2-adjoint(L)@X
        values, basis = eigh((gram+adjoint(gram))/2, check_finite=False)
        positive_modes = basis[:, values > 0.]
    proposal_t = {'wall_seconds': perf_counter()-begin, 'cpu_seconds': process_time()-cpu}
    ball_start, ball_cpu = perf_counter(), process_time()
    matrix = lambda z: ball_matrix(z, acb, acb_mat)
    S, G, B, Y, X, C = map(matrix, (sector.S, sector.G, sector.B, sector.Y, X, factor))
    I = matrix(sector.I)
    a, b, c = acb(arb(str(u))), acb(arb(str(mu))), acb(shift)
    H = (G+G.conjugate().transpose())/2
    St = S.conjugate().transpose()
    K = St*S + b*(I-a*H)
    L = St*Y + (a*b/2)*B
    difference = K-c*I - C*C.conjugate().transpose()
    frobenius_error = arb_norm_squared(difference, arb).sqrt().upper()
    alpha = (arb(shift)-frobenius_error).lower()
    if not alpha > 0.:
        raise ArithmeticError('Arb cannot prove positive definiteness from the shifted Cholesky residual.')
    KX = K*X
    rhs_squared = arb_norm_squared(KX-L, arb).upper()
    Xt = X.conjugate().transpose()
    raw = arb_norm_squared(Y, arb)+(Xt*(KX-2*L)).trace().real
    lower = (raw-rhs_squared/alpha).lower()

    def bounds(value):
        return [float(np.nextafter(float(value.lower()), -np.inf)),
                float(np.nextafter(float(value.upper()), np.inf))]

    delta_comparison = {key: bool(lower > arb(radius)**2) for key, radius in sector.radii.items()}
    result = {
        'u': nominal_u, 'mu': nominal_mu, 'status': 'Arb_verified', 'precision_bits': ctx.prec,
        'nominal_eigenvalue_used_only_to_propose_shift': float(observable_minimum),
        'exact_dyadic_Cholesky_shift': shift, 'K_minus_shift_I_minus_CCstar_Frobenius_upper': bounds(frobenius_error)[1],
        'K_lambda_min_strict_lower': bounds(alpha)[0], 'KX_minus_L_Frobenius_upper': bounds(rhs_squared.sqrt())[1],
        'quadratic_value_at_X_interval': bounds(raw), 'solve_defect_squared_over_alpha_upper': bounds(rhs_squared/alpha)[1],
        'distance_lower_bound_squared': max(0., bounds(lower)[0]),
        'distance_lower_bound': max(0., bounds(max(arb(0), lower).sqrt())[0]),
        'strictly_excludes': delta_comparison,
        'numeric_proposal': {**proposal_t, 'eig_min': 1, 'cholesky': 2, 'source_RHS': source_rhs,
                             'illumination_gram_eigh': int(project_positive_gram)},
        'Arb': {'wall_seconds': perf_counter()-ball_start, 'cpu_seconds': process_time()-ball_cpu,
                'matrix_dimension': n, 'source_RHS_for_residual': source_rhs,
                'full_square_ball_matrix_products': 1},
    }
    if positive_modes is not None:
        def project(basis):
            subspace_begin, subspace_cpu = perf_counter(), process_time()
            rank = basis.shape[1]
            V = matrix(basis)
            Vt = V.conjugate().transpose()
            gram_defect = arb_norm_squared(Vt*V-matrix(np.eye(rank)), arb).sqrt().upper()
            beta = (1+gram_defect).upper()
            XV, YV = X*V, Y*V
            value = arb_norm_squared(YV, arb)+(XV.conjugate().transpose()*((KX-2*L)*V)).trace().real
            defect = arb_norm_squared((KX-L)*V, arb).upper()
            bound = ((value-defect/alpha)/beta).lower()
            return {
                'proposal_rank': rank, 'gram_defect_Frobenius_upper': bounds(gram_defect)[1],
                'projector_beta_upper': bounds(beta)[1], 'quadratic_value_at_XV_interval': bounds(value),
                'solve_defect_squared_over_alpha_upper': bounds(defect/alpha)[1],
                'distance_lower_bound_squared': max(0., bounds(bound)[0]),
                'distance_lower_bound': max(0., bounds(max(arb(0), bound).sqrt())[0]),
                'strictly_excludes': {key: bool(bound > arb(radius)**2) for key, radius in sector.radii.items()},
                'wall_seconds': perf_counter()-subspace_begin, 'cpu_seconds': process_time()-subspace_cpu,
            }

        result['positive_gram_projection'] = project(positive_modes)
        result['single_coherent_projection'] = project(positive_modes[:, -1:])
    result.update(wall_seconds=perf_counter()-begin, cpu_seconds=process_time()-cpu)
    return result


def scan_case(case, radii, spectral_cache, numeric_only, gram_matrices):
    begin, cpu = perf_counter(), process_time()
    with np.load(SOURCE/f'{case}.npz', allow_pickle=False) as data:
        Y = data['Y_archive']
    sector = Sector(case, Y, load_operators(case), radii, spectral_cache)
    preparations = {'wall_seconds': perf_counter()-begin, 'cpu_seconds': process_time()-cpu}
    grid = []
    for u in CAPS:
        candidates = [sector.numerical(u, mu, gram_matrices) for mu in WEIGHTS]
        usable = [row for row in candidates if row['status'] == 'numerical_only']
        best = max(usable, key=lambda row: row['objective_bound_squared_unclipped']) if usable else None
        best_positive = max(usable, key=lambda row: row['positive_gram_bound_squared']) if usable else None
        best_single = max(usable, key=lambda row: row['best_single_coherent_squared']) if usable else None
        best_columns = max(usable, key=lambda row: row['positive_original_columns_squared']) if usable else None
        grid.append({'u': u, 'bounds_by_mu': candidates, 'best_grid_result': best, 'best_positive_gram_result': best_positive,
                     'best_single_coherent_result': best_single, 'best_positive_original_columns_result': best_columns,
                     'controls': sector.magnitude_controls(u)})
    result = {
        'state_dimension': sector.B.shape[0], 'receiver_source_shape': list(Y.shape),
        'norm_Y': float(np.linalg.norm(Y)), 'archived_radii': sector.radii,
        'operator_norms_numerical': {'G_2': sector.Gnorm, 'S_2': sector.Snorm, 'B_Frobenius': sector.Bnorm,
                                     'G_H_lambda_max': float(sector.h[-1])},
        'prepare_timing': preparations, 'grid': grid, 'numeric_cost': dict(sector.cost),
        'numeric_wall_seconds': perf_counter()-begin, 'numeric_cpu_seconds': process_time()-cpu,
    }
    if case in CHECKS and not numeric_only:
        result['independent_certificate'] = audit(sector, *CHECKS[case])
    if case in PROJECTED_CHECKS and not numeric_only:
        result['positive_gram_projection_certificate'] = audit(sector, *PROJECTED_CHECKS[case], project_positive_gram=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE/'results.json')
    parser.add_argument('--numerics-only', action='store_true', help='skip the Arb certificates')
    args = parser.parse_args()
    begin, cpu = perf_counter(), process_time()
    radii_by_case = read_json(OBSERVATION_RADII)
    spectral_cache = {}
    gram_matrices = {}
    results = {}
    with threadpool_limits(limits=1, user_api='blas'):
        for case, radii in radii_by_case.items():
            results[case] = scan_case(case, radii, spectral_cache, args.numerics_only, gram_matrices)
            case_result = results[case]
            audits = sum(key in case_result for key in ('independent_certificate', 'positive_gram_projection_certificate'))
            print(case, f'{audits} audited grid entries', flush=True)
    array_name = 'gram_matrices.npz' if args.output == HERE/'results.json' else f'{args.output.stem}_gram_matrices.npz'
    result = {
        'scope': 'Frozen binary64 scalar TM operators and archived observation radii; material caps used alone, no labels',
        'caps': CAPS, 'multipliers': WEIGHTS, 'certified_parameters': CHECKS if not args.numerics_only else {},
        'projected_certified_parameters': PROJECTED_CHECKS if not args.numerics_only else {},
        'cases': results, 'gram_arrays': array_name,
        'wall_seconds_before_write': perf_counter()-begin, 'cpu_seconds_before_write': process_time()-cpu,
    }
    np.savez_compressed(args.output.with_name(array_name), **gram_matrices)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
