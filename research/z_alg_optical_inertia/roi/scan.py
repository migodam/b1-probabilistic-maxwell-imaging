#!/usr/bin/env python3
"""ROI proposal sweep on strictly public observation arrays."""
import os

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'

import argparse
from collections import Counter
from pathlib import Path
import sys
import json

import numpy as np
from scipy.linalg import eigh, ldl, solve_triangular

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
RESEARCH = STUDY.parent
WORKSPACE = RESEARCH.parent
OUTPUT = STUDY / 'results/roi'
HIGH = RESEARCH / 'z_alg_secant_tasks/highdim'
LA = RESEARCH / 'limited_aperture_routed_dr_20261120'
SRAP = RESEARCH / 'srap_cie_array_study_20261121'
R8 = RESEARCH / 'delegated/max_struct_r1/r8_worktree'
sys.path.insert(0, str(STUDY))
from core.inertia import Prepared, adjoint, clock, hermitian, result_row, spectral_tolerance, timing

CANONICAL = ('austria_eps4_400MHz', 's2d_007__pw180_10x20',
             's2d_007__pw360_10x20', 's2d_010__pw180_10x20')
LEGACY = ('la180_eccentric_shell_eps4p5', 'la180_core_shell_eps5p0', 'la180_double_rings_eps4p0')
CAPS = (0., 0.1, 0.5, 1., 2.)
RATIOS = (0., 0.5, 2., 8., 32., 128.)
ALPHAS = tuple(float(2.**i) for i in (-10, -8, -6, -4, -2, 0, 2, 4))
OPTICAL_CONTRAST_WEIGHTS = ((1., 2.), (4., 1.), (4., 2.))
REGIONS = np.arange(1024).reshape(4, 8, 4, 8).transpose(0, 2, 1, 3).reshape(16, 64)
MASK_OTHER = [np.setdiff1d(np.arange(1024), inside) for inside in REGIONS]


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def write(path, obj):
    temporary = path.with_suffix('.partial')
    temporary.write_text(json.dumps(obj, indent=1, ensure_ascii=False) + '\n', encoding='utf8')
    temporary.replace(path)


def observations():
    radii = read(HIGH / 'observation_radii.json')
    with np.load(HIGH / 'public_vie_inputs.npz', allow_pickle=False) as operators:
        G = operators['G_1']
        for case in CANONICAL:
            path = HIGH / f'results/{case}.npz'
            with np.load(path, allow_pickle=False) as data:
                Y = data['Y_archive']
            r = {}
            if 'original_R8_noise_ball' in radii[case]:
                r['sample_noise_ball'] = radii[case]['original_R8_noise_ball']
            r['historical_observation_ball'] = radii[case]['frozen_observation_ball']
            yield case, Prepared(G, operators[f'{case}__S'], operators[f'{case}__E0'], Y), r, str(path.relative_to(WORKSPACE))
    sys.path.insert(0, str(R8))
    from research.research_max_struct_r8.src.paper_physics import acquisition
    G = acquisition(32, 400e6)['G']
    for case in LEGACY:
        path = LA / f'.runs/public/{case}.npz'
        with np.load(path, allow_pickle=False) as data:
            yield case, Prepared(G, data['S'], data['E'], data['Y']), {
                'sample_noise_ball': float(data['delta'])}, str(path.relative_to(WORKSPACE))


def shared_geometry(observation_list):
    groups = []
    for row in observation_list:
        _, problem, _, _ = row
        for group in groups:
            ref = group[0][1]
            if all(np.array_equal(a, b) for a, b in ((problem.G, ref.G), (problem.S, ref.S), (problem.B, ref.B))):
                group.append(row)
                break
        else:
            groups.append([row])
    return groups


def invert_hermitian(K, update_limit):
    """LDL and explicit inverse for reused proposal bases; compact block pivots."""
    factor, D, perm = ldl(K, lower=True, hermitian=True, check_finite=False)
    n = K.shape[0]
    sub = np.flatnonzero(np.diag(D, k=1))
    whole = np.setdiff1d(np.arange(n), np.r_[sub, sub+1])
    diagonal = D.diagonal()[whole].real
    tol = spectral_tolerance(np.r_[diagonal, D[sub, sub+1]])
    if np.any(np.abs(diagonal) <= tol):
        return None
    nnegative = int(np.count_nonzero(diagonal < 0))
    for k in sub:
        block = D[k:k+2, k:k+2]
        eigen = np.linalg.eigvalsh(block)
        if np.min(np.abs(eigen)) <= tol:
            return None
        nnegative += int(np.count_nonzero(eigen < 0))
    if nnegative >= update_limit:
        return None, nnegative
    tri_inv = solve_triangular(factor[perm], np.eye(n), lower=True, unit_diagonal=True, check_finite=False)
    inverse_times_factor = tri_inv.copy()
    inverse_times_factor[whole] /= diagonal[:, None]
    for k in sub:
        inverse_times_factor[k:k+2] = np.linalg.solve(D[k:k+2, k:k+2], tri_inv[k:k+2])
    index = np.argsort(perm)
    inverse = hermitian(adjoint(tri_inv) @ inverse_times_factor)
    return inverse[np.ix_(index, index)], nnegative


def settings(cap, U, rho):
    yield 'real_global' if rho == 0 else 'optical_global', 1., 1.
    if cap == 0:
        yield 'exact_zero_current', None, None
    choices = [(4., 1.)] if rho == 0 else list(OPTICAL_CONTRAST_WEIGHTS)
    if cap > 0:
        choices += [(U / cap, 1.)]
    for weight, equality in dict.fromkeys(choices):
        yield 'real_two_zone' if rho == 0 else 'optical_two_zone', weight, equality


def pivot_solve(M, rhs):
    values, vectors = eigh(hermitian(M), check_finite=False, driver='evd')
    if np.min(np.abs(values)) <= spectral_tolerance(values):
        return None
    nnegative = int(np.count_nonzero(values < 0.))
    return vectors @ ((adjoint(vectors) @ rhs) / values[:, None]), nnegative


def scores(H, nu, sources):
    if nu >= sources:
        return result_row(nu, sources)
    return result_row(nu, sources, H)


def candidate(cap, U, mu, rho, w, q, inside, family, base_inverse, inverse_Gstar, problem,
              base_Y, base_L, base_X, nu0, small_blocks):
    sources = problem.s
    # Illumination blocks are stacked horizontally so one LDL update covers all observed objects.
    C, C12, C22 = small_blocks
    XR = base_X[inside]
    if family == 'exact_zero_current':
        solution = pivot_solve(C, XR)
        if solution is None:
            return None
        correction, removed = solution
        nu = nu0 - removed
        if nu < 0:
            return None
        X = base_X - base_inverse[:, inside] @ correction
        L = base_L
    else:
        delta_z = 0.5 * mu * (w * cap - U + 1j * rho * (q - 1.))
        diagonal = mu * (w - 1.)
        L = base_L.copy()
        L[inside] += delta_z * base_Y
        X0 = base_X + base_inverse[:, inside] @ (delta_z * base_Y)
        if delta_z == 0.:
            if diagonal == 0.:
                nu, X = nu0, X0
            else:
                correction_matrix = C.copy()
                correction_matrix.flat[::len(inside)+1] += 1./diagonal
                solution = pivot_solve(correction_matrix, X0[inside])
                if solution is None:
                    return None
                correction, nn = solution
                nu = nu0 + int(diagonal < 0.) - nn
                X = X0 - base_inverse[:, inside] @ correction
        else:
            G_R = problem.G[inside]
            N = len(inside)
            extent = 2*N
            Dinv = np.zeros((extent, extent), complex)
            # Invert each two-by-two coupling block [[d, -z], [-z*, 0]].
            k = np.arange(N)
            Dinv[k, N+k] = -1./delta_z.conjugate()
            Dinv[N+k, k] = -1./delta_z
            Dinv[N+k, N+k] = -diagonal / abs(delta_z)**2
            Dinv[:N,:N] += C
            Dinv[:N,N:] += C12
            Dinv[N:,:N] += adjoint(C12)
            Dinv[N:,N:] += C22
            solution = pivot_solve(Dinv, np.concatenate((X0[inside], G_R @ X0)))
            if solution is None:
                return None
            correction, nn = solution
            nu = nu0 + N - nn
            X = X0 - base_inverse[:, inside] @ correction[:N] - inverse_Gstar[:, inside] @ correction[N:]
        if not 0 <= nu <= problem.n:
            return None
    lower = []
    for xn, ln, y2 in zip(np.split(X, X.shape[1]//sources, axis=1),
                          np.split(L, L.shape[1]//sources, axis=1), base_Y.grams, strict=True):
        lower.append(scores(y2-adjoint(ln)@xn, nu, sources))
    return lower


class IlluminationStack:
    def __init__(self, geometry, count):
        self.geometry = geometry
        self.count = count

    def __rmul__(self, coefficient):
        return np.tile(coefficient * self.geometry, (1, self.count))


def scan(group, U):
    started = clock()
    ref = group[0][1]
    count = len(group)
    results = {case: {'source': path, 'radii': radii, 'Y_norm': float(np.linalg.norm(problem.Y)),
             'tested_caps': list(CAPS), 'regions': [{'row': i//4, 'col': i%4,
             'trials': {}} for i in range(16)]} for case, problem, radii, path in group}
    operations = Counter()
    mu_scale = float(np.vdot(ref.S, ref.S).real / ref.n)
    base_inertias = []
    for rho in RATIOS:
        for alpha in ALPHAS:
            mu = mu_scale * alpha
            K, Lref = ref.matrices(U, mu, mu*rho)
            base_solve = invert_hermitian(K, ref.s + len(REGIONS[0]))
            operations['base_LDL'] += 1
            if base_solve is None:
                operations['base_pivot_unresolved'] += 1
                continue
            Kinv, nu0 = base_solve
            # Each region update has at most 64 positive eigenvalues, even if rank 128.
            # Removing 64 coordinates admits the same lower bound by interlacing.
            base_inertias.append(dict(mu=mu, ratio=rho, nu_base=nu0,
                                      nu_lower_all_regions=max(0, nu0-len(REGIONS[0]))))
            if Kinv is None:
                operations['base_rank_rejected'] += 1
                continue
            operations['base_inverse'] += 1
            KGstar = Kinv @ adjoint(ref.G)
            all_L = np.concatenate([Lref + p.SY-ref.SY for _, p, _, _ in group], axis=1)
            all_X = Kinv @ all_L
            defect = np.linalg.norm(K @ all_X-all_L)/max(np.linalg.norm(all_L), np.finfo(float).tiny)
            if not defect < 1e-6:
                operations['base_solve_residual_rejected'] += 1
                continue
            for index, inside in enumerate(REGIONS):
                illuminations = IlluminationStack(ref.B[inside], count)
                illuminations.grams = tuple(p.Y2 for _, p, _, _ in group)
                indices = np.ix_(inside, inside)
                small_blocks = Kinv[indices], KGstar[indices], ref.G[inside] @ KGstar[:, inside]
                for cap in CAPS:
                    for arm, w, q in settings(cap, U, rho):
                        out = candidate(cap, U, mu, rho, w, q, inside, arm, Kinv, KGstar,
                                        ref, illuminations, all_L, all_X, nu0, small_blocks)
                        operations['rank_update_trials'] += 1
                        if out is None:
                            operations['rank_update_unresolved'] += 1
                            continue
                        for (case, _, _, _), row in zip(group, out, strict=True):
                            branch = results[case]['regions'][index]['trials'].setdefault(str(cap), {})
                            current = branch.get(arm)
                            order = (row['lower_bound_numeric'], -row['nu'])
                            if current is None or order > (current['lower_bound_numeric'], -current['nu']):
                                branch[arm] = {**row, 'mu_background': mu, 'lambda_over_mu_background': rho,
                                               'mu_region_over_mu_background': w,
                                               'lambda_region_over_lambda_background': q,
                                               'base_relative_solve_defect': float(defect)}
        print(f'U={U:g}, cases={list(results)}, lambda/mu={rho:g}', flush=True)
    for result in results.values():
        for region in result['regions']:
            for cap, trials in region['trials'].items():
                if trials:
                    region.setdefault('winners', {})[cap] = max(trials, key=lambda x: (
                        trials[x]['lower_bound_numeric'], -trials[x]['nu']))
    return {'case_results': results, 'mu_scale': mu_scale, 'base_inertias': base_inertias,
            'operations': dict(operations), **timing(started)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--U', nargs='+', type=float, default=[20., 4.])
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--cases', nargs='+')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    groups = shared_geometry([row for row in observations() if not args.cases or row[0] in args.cases])
    for U in args.U:
        payload = dict(U=U, tested_caps=CAPS, ratios=RATIOS, mu_scale_multipliers=ALPHAS,
                       optical_two_zone=OPTICAL_CONTRAST_WEIGHTS, geometries=[])
        path = args.output / f'scan_U{U:g}.json'
        for group in groups:
            payload['geometries'].append(scan(group, U))
            write(path, payload)
            print(f'Wrote {path.name}', flush=True)


if __name__ == '__main__':
    main()
