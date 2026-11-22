"""Numerically resolve the dark-subspace threshold and the saved Green convention."""
import argparse
from pathlib import Path

import numpy as np
from scipy.linalg import cholesky, eigh, eigvalsh, qr, svdvals
from scipy.special import j0, jv
from scipy.spatial.distance import cdist
from threadpoolctl import threadpool_limits

from inertia import adjoint, clock, hermitian, timing
from scan import RESULTS, load_cases, save_json

RATIOS = (-100., -10., -1., 0., .1, .3, 1., 3., 10., 30., 100., 300., 1000., 10_000.)


def threshold(problem, ratios):
    start = clock()
    s = svdvals(problem.S, check_finite=False)
    q, _ = qr(adjoint(problem.S), mode='full', check_finite=False)
    dark = q[:, problem.r:]
    H = hermitian(adjoint(dark) @ (problem.Gh @ dark))
    R = hermitian(adjoint(dark) @ (problem.Ga @ dark))
    d = H.shape[0]
    eye = np.eye(d)
    r_eig = eigvalsh(R, check_finite=False)
    rows = []
    for t in ratios:
        began = clock()
        P = eye + t * R
        try:
            cholesky(P, check_finite=False)
        except np.linalg.LinAlgError:
            rows.append({'lambda_over_mu': t, 'u_star_numeric': None,
                         'status': 'I_plus_t_ImH_on_nullspace_not_PD_numeric', **timing(began)})
            continue
        biggest = float(eigh(H, P, subset_by_index=(d - 1, d - 1), eigvals_only=True,
                             check_finite=False, driver='gvx')[0])
        rows.append({'lambda_over_mu': t, 'u_star_numeric': None if biggest <= 0. else 1. / biggest,
                     'generalized_max_eigenvalue': biggest, 'status': 'numeric', **timing(began)})
    zero = next(row['u_star_numeric'] for row in rows if row['lambda_over_mu'] == 0.)
    for row in rows:
        if row['u_star_numeric'] is not None and zero is not None:
            row['gain_over_zero_numeric'] = row['u_star_numeric'] - zero
    return {'dark_dimension_nominal': d, 'receiver_smallest_singular_value_numeric': float(s[-1]),
            'receiver_largest_singular_value_numeric': float(s[0]),
            'S_times_dark_basis_Frobenius_numeric': float(np.linalg.norm(problem.S @ dark)),
            'dark_ImH_min_eigenvalue_numeric': float(r_eig[0]),
            'dark_ImH_max_eigenvalue_numeric': float(r_eig[-1]),
            'rows': rows, 'generalized_eigenproblems': sum(row['status'] == 'numeric' for row in rows),
            'receiver_singular_value_decompositions': 1, 'receiver_full_QR_factorizations': 1,
            'dark_radiation_eigenproblems': 1, **timing(start)}


def convention(problem, positions, meta):
    """Compare H1 disk radiation with original positive imaginary coefficient, unweighted."""
    acq = meta.get('acquisition', {})
    k = float(acq.get('k', 2. * np.pi * 400.e6 / 3.e8))
    side = round(np.sqrt(problem.n))
    spacing = np.ptp(positions, axis=0).max() / (side - 1)
    radius = float(acq.get('disk_radius', spacing / np.sqrt(np.pi)))
    coefficient = .5 * np.pi * k * radius * jv(1, k * radius)
    analytic = coefficient * j0(k * cdist(positions, positions))
    defect = np.linalg.norm(problem.Ga - analytic)
    diagonal = np.diag(problem.Ga).real
    return {'k': k, 'disk_radius': radius, 'j0_coefficient': coefficient,
            'actual_Im_H_diagonal_min_max': [float(diagonal.min()), float(diagonal.max())],
            'relative_Frobenius_difference_vs_H1_disk_radiation_formula': float(
                defect / np.linalg.norm(analytic)),
            'relative_Frobenius_difference_vs_conjugate_H2_disk_radiation_formula': float(
                np.linalg.norm(problem.Ga + analytic) / np.linalg.norm(analytic)),
            'operator_coordinate_order': 'as saved',
            'mass_normalization': 'constant pulse radius within the grid; no additional cell area multiplier'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=RESULTS / 'dark_thresholds.json')
    parser.add_argument('--cohort', choices=('archive', 'synthetic', 'all'), default='all')
    args = parser.parse_args()
    started = clock()
    output = {'numeric_only': True,
              'definition': 'first loss of positivity of I-u Re_H(G)+lambda_over_mu Im_H(G) restricted to ker(S)',
              'cases': {}, 'geometries': {}}
    with threadpool_limits(limits=1, user_api='blas'):
        for group in load_cases(args.cohort).values():
            previously_run = []
            for case, (p, _, pos, meta) in group.items():
                previous = next((key for S, key in previously_run if np.array_equal(S, p.S)), None)
                if previous is None:
                    output['geometries'][case] = threshold(p, RATIOS)
                    output['geometries'][case]['convention'] = convention(p, pos, meta)
                    previously_run.append((p.S, case))
                    row = output['geometries'][case]
                    concise = [(r['lambda_over_mu'], r['u_star_numeric']) for r in row['rows']
                               if r['lambda_over_mu'] in (0., 10., 100., 1000., 10000.)]
                    print(case, concise, flush=True)
                output['cases'][case] = previous or case
                output.update(timing(started))
                save_json(args.output, output)


if __name__ == '__main__':
    main()
