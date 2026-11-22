#!/usr/bin/env python3
"""Freeze inputs and screen the same MAX-STRUCT geometry with a shared factor per profile."""
import argparse
from collections import Counter
from itertools import product
import json
import os
from pathlib import Path
import sys

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[name] = '1'

import numpy as np
from scipy.linalg import eigh, eigvalsh
from scipy.linalg.lapack import get_lapack_funcs
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
ROOT = STUDY.parents[1]
ARCHIVE = ROOT/'research/delegated/max_struct_r1/r14_worktree/research/research_max_struct_r14'
OLD_OPERATOR = ROOT/'research/z_alg_secant_tasks/highdim/public_vie_inputs.npz'
sys.path.insert(0, str(STUDY))
from core.inertia import adjoint, clock, hermitian, spectral_tolerance, timing

R10_REFERENCE = (5., .001, 0.)


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def write(path, data):
    pending = path.with_suffix('.partial')
    pending.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    pending.replace(path)


def export_public(freeze, archive_dir, output):
    manifest = read(archive_dir/'MANIFEST.json')
    info = {r['case_id']: r for r in manifest['case_files']}
    arrays = ('G', 'S', 'E0', 'centers', 'sources', 'receivers', 'grid_shape', 'k', 'radius', 'cell_radius')
    first = None
    ys, meta = [], []
    for case in freeze['cases']:
        original = info[case]
        with np.load(archive_dir/f'{case}.npz', allow_pickle=False) as f:
            geo = {key: f[key] for key in arrays}
            # The published dimensions and the persisted NPZ both participate.
            (
                side, fine, noise, delta, sigma, lb, ub
            ) = (float(f[key]) for key in ('inverse_side', 'forward_side', 'noise_relative_rms', 'delta',
                                          'sigma', 'lower', 'upper'))
            y = f['Y']
        if first is None:
            first = geo
        if (
            not all(np.array_equal(geo[key], first[key]) for key in arrays)
            or (side, fine, noise, y.shape) != (24, 24, .01, (48, 12))
            or delta != original['delta_noise_only'] or sigma != original['sigma']
        ):
            raise ValueError(f'Observation disagrees with the frozen common geometry or public noise manifest: {case}')
        ys.append(y)
        meta.append(dict(case=case, relative_noise_rms=noise, sigma=sigma, delta=delta,
                         lower=lb, upper=ub, observation_version='24to24-12tx-48rx',
                         source=original['input_file'], norm_Y=float(np.linalg.norm(y))))
    r10_key_matches = {}
    if OLD_OPERATOR.exists():
        with np.load(OLD_OPERATOR, allow_pickle=False) as old:
            r10_key = old['r10_twod_0001__G_key'].item()
            r10_key_matches = {key: bool(np.array_equal(first[key], old[name])) for key, name in (
                ('G', r10_key), ('S', 'r10_twod_0001__S'), ('E0', 'r10_twod_0001__E0'))}
    np.savez_compressed(output, **first, Y=np.stack(ys), cases=np.array(freeze['cases']),
                        delta=np.array([row['delta'] for row in meta]),
                        sigma=np.array([row['sigma'] for row in meta]))
    report = dict(source_archive=str(archive_dir.relative_to(ROOT)) if archive_dir.is_relative_to(ROOT) else str(archive_dir),
                  noise_model=manifest['noise_model'], delta_definition=manifest['delta_definition'],
                  seven_case_geometry_exact_shared=True, matches_archived_r10_arrays_exactly=r10_key_matches,
                  observations=meta)
    write(HERE/'public_provenance.json', report)


def load_public(path, contract):
    with np.load(path, allow_pickle=False) as f:
        if f['cases'].tolist() != contract['cases']:
            raise ValueError('Public case list does not match the freeze.')
        return f['G'], f['S'], f['E0'], f['Y'], f['delta']


def decode_ldl(lu, pivots):
    diag = lu.diagonal().real.copy()
    negative_pivots = np.flatnonzero(pivots < 0)
    if negative_pivots.size:
        paired_runs = np.split(negative_pivots, np.flatnonzero(np.diff(negative_pivots) > 1) + 1)
        for pair in paired_runs:
            for i in pair[::2]:
                middle = (diag[i] + diag[i+1])*.5
                offset = np.hypot((diag[i] - diag[i+1])*.5, abs(lu[i+1,i]))
                diag[i:i+2] = (middle-offset, middle+offset)
    return diag


def profile_result(K, L, shape, solver, total, cost):
    m, s = shape
    factor, work, solve = solver
    lu, pivots, code = factor(K, lower=True, lwork=work, overwrite_a=False)
    cost['geometry_LDL_factorizations'] += 1
    if not code:
        pivvals = decode_ldl(lu, pivots)
        if np.min(abs(pivvals)) > spectral_tolerance(pivvals):
            nu = int((pivvals < 0).sum())
            if nu >= s:
                return nu, 'numeric_source_exhaustion', []
            X, err = solve(lu, pivots, L, lower=True)
            cost['LDL_RHS'] += L.shape[1]
            defect = float(np.linalg.norm(K@X-L) / np.linalg.norm(L))
            if not err and defect < 1.e-8:
                return nu, 'numeric_shared_LDL', spectra(L, X, total, m, s, cost)
    vals, sub = eigh(K, check_finite=False, driver='evd')
    cost['fallback_eigendecompositions'] += 1
    if np.min(abs(vals)) <= spectral_tolerance(vals):
        return None, 'numeric_near_singular', []
    nu = int((vals < 0).sum())
    if nu >= s:
        return nu, 'numeric_source_exhaustion_fallback', []
    X = sub @ ((sub.conj().T @ L) / vals[:,None])
    cost['eigen_RHS'] += L.shape[1]
    return nu, 'numeric_eigen_fallback', spectra(L, X, total, m, s, cost)


def spectra(L, X, Y2, m, s, cost):
    case_matrices = []
    for ix in range(m):
        part = slice(ix*s, (ix+1)*s)
        H = hermitian(Y2[ix] - L[:,part].conj().T @ X[:,part])
        ds = eigvalsh(H, check_finite=False)[::-1]
        case_matrices.append(dict(eig=ds.tolist(), diagonal=H.diagonal().real.tolist()))
        cost['source_Grams_and_eigenspectra'] += 1
    return case_matrices


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--public', type=Path, default=HERE/'public_inputs.npz')
    parser.add_argument('--export-from', type=Path, help='Read original public R14 independent24 NPZ files.')
    args = parser.parse_args()
    freeze = read(HERE/'frozen.json')
    started = clock()
    if args.export_from:
        export_public(freeze, args.export_from, args.public)
    G, S, B, Y, _ = load_public(args.public, freeze)
    m, _, s = Y.shape
    Y2 = [hermitian(adjoint(y) @ y) for y in Y]
    SY = np.column_stack([adjoint(S) @ y for y in Y])
    S2 = hermitian(adjoint(S) @ S)
    H = hermitian(G)
    O = (G-adjoint(G))/(2j)
    I = np.eye(len(G))
    ldlt, query, solve = get_lapack_funcs(('hetrf','hetrf_lwork','hetrs'), (G,))
    work, _ = query(len(G), lower=True)
    solver = ldlt, int(work.real), solve
    costs = Counter(receiver_Gram=1, observation_Grams=m, receiver_observation_products=m)
    setup = timing(started)
    values = []
    for u in freeze['caps']:
        for t, mu in product(freeze['lambda_over_mu'], freeze['mu']):
            K = S2 + mu*(I-u*H+t*O)
            lam = mu*t
            zB = (mu*u+1j*lam)*.5*B
            L = SY + np.tile(zB, m)
            nu, status, case_matrices = profile_result(K, L, (m,s), solver, Y2, costs)
            values.append(dict(profile=len(values), cap=u, mu=mu, lambda_over_mu=t, lambda_=lam, nu=nu,
                               status=status, historical_r10_profile=(u,mu,t)==R10_REFERENCE,
                               case_matrices=case_matrices))
        print('done cap', u, 'after', round(timing(started)['wall_seconds'],2), 's', flush=True)
    write(HERE/'numerical_grid.json', dict(cases=freeze['cases'], values=values, setup=setup,
                                           operations=dict(costs), **timing(started)))


if __name__ == '__main__':
    with threadpool_limits(1):
        main()
