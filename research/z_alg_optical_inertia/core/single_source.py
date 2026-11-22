"""Test all scalar multipliers against an exact dark mode and a public-data current."""
import argparse
from pathlib import Path

import numpy as np
from flint import acb, arb, ctx
from scipy.linalg import eigh, eigvalsh, qr, svd
from threadpoolctl import threadpool_limits

from certificate import ball, frobenius_squared, gershgorin, safe_float, star
from dark_witness import negative_witness
from inertia import adjoint, clock, hermitian, timing
from scan import RESULTS, load_cases, save_json

CAP = .75
CASES = ('optfresh_01__pw360_10x20__n02', 'optfresh_05__pw360_10x20__n02')


def numeric_blocks(p, J, u):
    E = p.B + p.G @ J
    JE = adjoint(J) @ E
    Qre = hermitian(adjoint(J) @ J - u * JE)
    Qim = (JE - adjoint(JE)) / (2j)
    return Qre, Qim


def fit_proposals(p, radius, u, max_ratio):
    started = clock()
    U, d, Vstar = svd(p.S, full_matrices=False, check_finite=False)
    V = adjoint(Vstar)
    Q, _ = qr(adjoint(p.S), mode='full', check_finite=False)
    N = Q[:, p.r:]
    dark_eig, basis = eigh(hermitian(adjoint(N) @ (p.Ga @ N)), driver='evd', check_finite=False)
    basis = N @ basis
    coord_B = adjoint(basis) @ (-.5j * p.B)
    coord_G = adjoint(basis) @ p.Ga
    Uy = adjoint(U) @ p.Y
    lower_best, per_source_best = None, [None] * p.s
    records = []
    for truncation in range(p.r, 0, -1):
        J0 = V[:, :truncation] @ (Uy[:truncation] / d[:truncation, None])
        defect = float(np.linalg.norm(p.S @ J0 - p.Y))
        if defect >= .98 * radius:
            continue
        rhs = coord_G @ J0 + coord_B
        for beta in np.logspace(-13., 1., 43):
            J = J0 - basis @ (rhs / (dark_eig[:, None] + beta))
            Qr, Qi = numeric_blocks(p, J, u)
            T = -Qi
            R = max_ratio * T - Qr
            eig_T, eig_R = eigvalsh(T, check_finite=False), eigvalsh(R, check_finite=False)
            ts, rs = T.diagonal().real, R.diagonal().real
            defect = float(np.linalg.norm(p.S @ J - p.Y))
            score = min(eig_T[0] / max(abs(eig_T[-1]), 1.), eig_R[0] / max(abs(eig_R[-1]), 1.),
                        (radius - defect) / radius)
            row = {'truncation': truncation, 'beta': float(beta), 'fit': defect, 'score': float(score),
                   'min_T': float(eig_T[0]), 'min_ratio_T_minus_Qre': float(eig_R[0]),
                   'both_positive_diagonal_count': int(np.count_nonzero((ts > 0.) & (rs > 0.)))}
            records.append(row)
            if lower_best is None or score > lower_best[0]['score']:
                lower_best = (row, J)
            for s in range(p.s):
                value = min(float(ts[s] / max(abs(eig_T[-1]), 1.)),
                            float(rs[s] / max(abs(eig_R[-1]), 1.)), (radius - defect) / radius)
                old = per_source_best[s]
                if old is None or value > old[0]:
                    per_source_best[s] = (value, J[:, s:s+1], row)
    return lower_best, per_source_best, records, timing(started)


def verify_positive(A, Anp):
    _, modes = eigh(Anp, check_finite=False)
    V = ball(modes)
    return gershgorin(star(V) * A * V, upper=False)


def verify_fit(p, Jnp, u, t, radius, columns=None):
    if columns is None:
        B_np, Y_np = p.B, p.Y
    else:
        B_np, Y_np = p.B[:, columns], p.Y[:, columns]
    G, S, B, Y, J = map(ball, (p.G, p.S, B_np, Y_np, Jnp))
    JE = star(J) * (B + G * J)
    Re = star(J) * J - arb(u) * (JE + star(JE)) / 2
    T = (star(JE) - JE) / acb(0., 2.)
    R = arb(t) * T - Re
    E_np = B_np + p.G @ Jnp
    JE_np = adjoint(Jnp) @ E_np
    T_np = -(JE_np - adjoint(JE_np)) / (2j)
    R_np = t * T_np - hermitian(adjoint(Jnp) @ Jnp - u * JE_np)
    T_low, R_low = verify_positive(T, T_np), verify_positive(R, R_np)
    data_upper = frobenius_squared(S * J - Y).sqrt().upper()
    return {'strict_fit': bool(data_upper < arb(radius)), 'strict_T_positive_definite': bool(T_low > 0.),
            'strict_t_T_minus_Qre_positive_definite': bool(R_low > 0.),
            'T_congruence_lower': safe_float(T_low, upper=False),
            't_T_minus_Qre_congruence_lower': safe_float(R_low, upper=False),
            'Frobenius_data_error_upper': safe_float(data_upper, upper=True), 'noise_radius': radius}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=RESULTS / 'single_source_obstruction.json')
    parser.add_argument('--case', action='append', help='select full-ring 2% observations for the fixed u=.75 test')
    args = parser.parse_args()
    cases = args.case or CASES
    if any('__pw360_10x20__n02' not in case for case in cases):
        parser.error('The dark obstruction in this experiment is for the full-ring 2% release.')
    started = clock()
    data = next(iter(load_cases('synthetic').values()))
    results = {'cap': CAP, 'geometry': '32x32 PW360 10x20 exact frozen operators', 'cases': {}}
    t = 10_000.
    frozen = {}
    with threadpool_limits(limits=1, user_api='blas'):
        ctx.prec, ctx.threads = 256, 2
        witness, frozen['dark_proposal'] = negative_witness(data[CASES[0]][0], CAP, t, return_proposal=True)
        results['dark_witness_through_ratio'] = witness
        results['ratio'] = t
        for case in cases:
            p, radii, *_ = data[case]
            radius = radii['noise_only_fro']
            best, by_source, records, elapsed = fit_proposals(p, radius, CAP, t)
            result = {'search': records, 'numeric_cost': elapsed}
            if best is not None:
                row, X = best
                result['best_joint_numeric'] = row
                if row['score'] > 0.:
                    frozen[case] = X
                    result['joint_certificate'] = verify_fit(p, X, CAP, t, radius)
            if 'joint_certificate' not in result:
                result['individual_certificates'] = {}
                for index, proposal in enumerate(by_source):
                    if proposal is not None and proposal[0] > 0.:
                        value, X, row = proposal
                        result['individual_certificates'][str(index)] = {
                            'proposal': row, **verify_fit(p, X, CAP, t, radius, [index])}
            results['cases'][case] = result
            print(case, result.get('best_joint_numeric'), result.get('joint_certificate'),
                  'individual certificates', len(result.get('individual_certificates', {})), flush=True)
        arrays_path = args.output.with_suffix('.npz')
        np.savez_compressed(arrays_path, **frozen)
        results['proposal_arrays'] = arrays_path.name
        results.update(timing(started))
        save_json(args.output, results)


if __name__ == '__main__':
    main()
