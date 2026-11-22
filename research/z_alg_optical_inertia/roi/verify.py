#!/usr/bin/env python3
"""Audit local winners and common negative subspaces from saved binary64 arrays."""
from scan import MASK_OTHER, OUTPUT, REGIONS, observations, read, write

import argparse
import numpy as np
from scipy.linalg import cholesky, eigh
from flint import ctx

from core.certificate import _arb_inputs, ball, certify, frobenius_squared, gershgorin, safe_float, star
from core.inertia import Prepared, adjoint, clock, hermitian, timing


def negative_gap(M, proposed):
    """A verified positive lower bound for -M on the tested coordinates."""
    lower = gershgorin(-M, upper=False)
    if lower > 0.:
        return float(safe_float(lower, upper=False))
    values = np.linalg.eigvalsh(hermitian(-proposed))
    shift = values[0] / 2.
    if shift <= 0.:
        raise ArithmeticError('Endpoint restriction is not numerically negative.')
    F = ball(cholesky(-hermitian(proposed) - shift * np.eye(len(proposed)), lower=True))
    residual = -M - F * star(F)
    for j in range(len(proposed)):
        residual[j,j] -= shift
    lower = shift - frobenius_squared(residual).sqrt().upper()
    if lower <= 0.:
        raise ArithmeticError('Endpoint restricted Cholesky failed the enclosure.')
    return safe_float(lower, upper=False)


def dome(problem, U, mu_min, rho_min, rho_max, requested):
    """Prove at least k negative background directions over a continuous parameter box.

    The same projected subspace is checked at the two rho endpoints at minimum mu.
    Dividing by positive mu, increasing mu subtracts a PSD multiple of S*S.
    Any local diagonal weight/cap update shares the 960-square background principal submatrix.
    """
    started = clock()
    spectral_endpoint, _ = problem.matrices(U, mu_min, mu_min*rho_max)
    d, V = eigh(spectral_endpoint, subset_by_index=(0, requested-1), check_finite=False)
    if d[-1] >= 0.:
        return {'status': 'numeric_subspace_insufficient', 'proposed_final_negative_eigenvalue': float(d[-1]),
                'requested_common_negative_directions': requested, 'cost': timing(started)}
    old_bits, old_threads = ctx.prec, ctx.threads
    ctx.prec, ctx.threads = 96, 1
    try:
        N = ball(V)
        gaps = {}
        for rho in (rho_min, rho_max):
            K, _, _ = _arb_inputs(problem, U, mu_min, mu_min*rho)
            K_double, _ = problem.matrices(U, mu_min, mu_min*rho)
            gap = negative_gap(star(N) * (K * N), adjoint(V) @ K_double @ V)
            gaps[str(rho)] = gap
        return {'status': 'Arb_verified_continuous_rank_obstruction',
                'U': U, 'minimum_mu': mu_min, 'ratio_interval': [rho_min, rho_max],
                'verified_common_negative_directions': requested,
                'verified_negative_directions_each_region_lower': requested-64,
                'sources': problem.s, 'restricted_negative_endpoint_gaps': gaps,
                'observational_input_required': False, 'arb_threads': 1, 'cost': timing(started)}
    except ArithmeticError as e:
        return {'status': 'Arb_unresolved', 'reason': str(e), 'cost': timing(started)}
    finally:
        ctx.prec, ctx.threads = old_bits, old_threads



def background_proofs(problem, U, mu_min, rho_min, rho_max):
    """Check room for negative directions directly in each background principal block."""
    started = clock()
    old_bits, old_threads = ctx.prec, ctx.threads
    ctx.prec, ctx.threads = 96, 1
    requested = problem.s + 4
    reports = []
    try:
        exact = [_arb_inputs(problem, U, mu_min, mu_min*rho)[0] for rho in (rho_min, rho_max)]
        numerical = [problem.matrices(U, mu_min, mu_min*rho)[0] for rho in (rho_min, rho_max)]
        for index, outside in enumerate(MASK_OTHER):
            proposal = numerical[1][np.ix_(outside, outside)]
            eig, restricted = eigh(proposal, subset_by_index=(0,requested-1), check_finite=False)
            row = {'region': index, 'row': index // 4, 'column': index % 4}
            if eig[-1] >= 0.:
                row.update(status='numeric_subspace_insufficient', projected_max_eigenvalue=float(eig[-1]))
            else:
                V = np.zeros((problem.n, requested), dtype=complex)
                V[outside] = restricted
                N = ball(V)
                try:
                    gaps = [negative_gap(star(N)*(K*N), adjoint(V) @ K_float @ V)
                            for K, K_float in zip(exact, numerical, strict=True)]
                    row.update(status='Arb_verified', negative_directions_lower=requested,
                               verified_restricted_negative_endpoint_gaps=gaps)
                except ArithmeticError as e:
                    row.update(status='Arb_unresolved', reason=str(e))
            reports.append(row)
    finally:
        ctx.prec, ctx.threads = old_bits, old_threads
    complete = all(row['status']=='Arb_verified' for row in reports)
    return {'status': 'Arb_verified_continuous_rank_obstruction' if complete else 'Arb_incomplete_background_proofs',
            'U': U, 'minimum_mu': mu_min, 'ratio_interval':[rho_min, rho_max],
            'verified_negative_directions_each_region_lower':requested if complete else None,
            'sources':problem.s, 'observational_input_required':False, 'arb_threads':1,
            'region_proofs':reports, 'cost':timing(started)}


def parameters(problem, U, cap, arm, record, index):
    mu_bg, rho = record['mu_background'], record['lambda_over_mu_background']
    if arm == 'exact_zero_current':
        outside = MASK_OTHER[index]
        subproblem = Prepared(problem.G[np.ix_(outside, outside)], problem.S[:, outside],
                              problem.B[outside], problem.Y)
        return subproblem, U, mu_bg, mu_bg*rho
    inside = REGIONS[index]
    u = np.full(problem.n, U)
    mu = np.full(problem.n, mu_bg)
    lam = np.full(problem.n, mu_bg*rho)
    u[inside] = cap
    mu[inside] *= record['mu_region_over_mu_background']
    lam[inside] *= record['lambda_region_over_lambda_background']
    return problem, u, mu, lam


def verify(U, folder):
    scan = read(folder / f'scan_U{U:g}.json')
    jobs = {case: (p, radii) for case, p, radii, _ in observations()}
    result = {'U': U, 'shared_geometry_audits': [], 'selected_region_certificates': []}
    for group in scan['geometries']:
        cases = list(group['case_results'])
        problem, _ = jobs[cases[0]]
        scores = group['base_inertias']
        min_mu = min(s['mu'] for s in scores)
        min_nu = min(s['nu_base'] for s in scores)
        if min_nu >= problem.s + 64:
            requested = min(min_nu-4, problem.s+64+96)
            proof = dome(problem, U, min_mu, min(scan['ratios']), max(scan['ratios']), requested)
        else:
            proof = background_proofs(problem, U, min_mu, min(scan['ratios']), max(scan['ratios']))
        result['shared_geometry_audits'].append({'cases': cases, 'minimum_numeric_base_nu': min_nu, **proof})
        print(cases, proof['status'], proof.get('verified_negative_directions_each_region_lower'), flush=True)
        write(folder / f'verified_U{U:g}.json', result)
        for case, doc in group['case_results'].items():
            problem, radii = jobs[case]
            for index, region in enumerate(doc['regions']):
                trials = []
                for cap, arms in region['trials'].items():
                    arm = region['winners'][cap]
                    record = arms[arm]
                    if record['lower_bound_numeric'] > min(radii.values()):
                        trials.append((float(cap), arm, record))
                for cap, arm, record in sorted(trials, reverse=True):
                    work, u, mu, lam = parameters(problem, U, cap, arm, record, index)
                    verified = certify(work, u, mu, lam, radii=radii, threads=1)
                    result['selected_region_certificates'].append(dict(case=case, region=index,
                        cap=cap, arm=arm, proposal=record, verification=verified))
                    if verified.get('status') == 'Arb_verified' and any(verified['strictly_excludes'].values()):
                        break
    write(folder / f'verified_U{U:g}.json', result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--U', nargs='+', type=float, default=[20., 4.])
    parser.add_argument('--output', type=type(OUTPUT), default=OUTPUT)
    args = parser.parse_args()
    for U in args.U:
        verify(U, args.output)


if __name__ == '__main__':
    main()
