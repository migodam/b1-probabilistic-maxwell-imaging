"""Compare a single fixed-inertia 4×4 optical gradient and matched scalar updates."""
import argparse
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from inertia import adjoint, clock, hermitian, timing
from scan import RESULTS, block_ids, compact, key, load_cases, save_json
from summarize import load_numeric

STEPS = (.06, .18, .45)


def compare(problem, order, u, mu, lam):
    d = problem.evaluate_ldl(u, mu, lam, return_solution=True)
    X, L = d.pop('solution_current', None), d.pop('linear_term', None)
    baseline = {**compact(d), 'u': u, 'mu': mu, 'lambda': lam}
    nu = baseline['nu']
    if X is None or key(baseline)[0] <= 0.:
        return {'baseline': baseline, 'status': 'no_positive_retained_H_mode', 'rows': []}
    H = hermitian(problem.Y2 - adjoint(L) @ X)
    eigenvalues, modes = eigh(H, check_finite=False)
    eigenvalues, modes = eigenvalues[::-1], modes[:, ::-1]
    active = (np.arange(problem.s) >= nu) & (eigenvalues > 0.)
    V = modes[:, active]
    E = problem.B + problem.G @ X
    cell_gradient = np.sum(((X @ V).conj() * (E @ V)).imag, axis=1)
    gradients = np.bincount(order, weights=cell_gradient, minlength=16)
    scalar_gradient = float(gradients.sum())
    scale = max(abs(lam), mu)
    rows = []
    for multiplier in STEPS:
        magnitude = multiplier * scale
        candidates = (
            ('same_count_scalar', float(lam + np.sign(scalar_gradient) * magnitude)),
            ('gradient_blocks', lam + np.sign(gradients) * magnitude),
        )
        for method, trial in candidates:
            space = method == 'gradient_blocks'
            row = problem.evaluate_ldl(u, mu, trial[order] if space else trial)
            rows.append({**compact(row), 'method': method, 'step': magnitude,
                         'lambda': lam if space else trial,
                         'block_lambda': trial.tolist() if space else [], 'same_inertia': row['nu'] == nu})
    diff_cutoff = float(eigenvalues[nu - 1] - eigenvalues[nu]) if nu else None
    return {'baseline': baseline, 'status': 'tested', 'retained_positive_count': int(np.count_nonzero(active)),
            'tail_cutoff_eigen_gap_numeric': diff_cutoff, 'scalar_gradient': scalar_gradient,
            'block_gradients': gradients.tolist(), 'rows': rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=RESULTS / 'gradient_blocks.json')
    parser.add_argument('--case', action='append')
    args = parser.parse_args()
    scan, _ = load_numeric()
    started = clock()
    output = {'step_factors': STEPS, 'matched_scalar_count': len(STEPS), 'cases': {}}
    with threadpool_limits(limits=1, user_api='blas'):
        for geometry in load_cases('all').values():
            for case, (p, radii, positions, _) in geometry.items():
                if args.case is not None and case not in args.case:
                    continue
                caps = (5., 5.25) if case.startswith('r10') else (.5, .75, 1.)
                order = block_ids(positions)
                start = clock()
                result = {'radii': radii, 'caps': {}}
                for u in caps:
                    base = scan[case]['caps'][str(u)]['scalar_envelope']['global_signed']
                    row = compare(p, order, u, base['mu'], base['lambda'])
                    result['caps'][str(u)] = row
                    if row['status'] == 'tested':
                        b = row['baseline']
                        arms = {a: max((r for r in row['rows'] if r['method'] == a and r['same_inertia']),
                                       key=key, default=None)
                                for a in ('same_count_scalar', 'gradient_blocks')}
                        strongest = max((b, arms['same_count_scalar']), key=key)
                        space = arms['gradient_blocks']
                        print(case, u, b['nu'], key(b)[0], key(strongest)[0], key(space)[0],
                              'new passes', [r for r, delta in radii.items() if key(strongest)[0] <= delta < key(space)[0]],
                              flush=True)
                result['operations'] = dict(p.operations)
                result.update(timing(start))
                output['cases'][case] = result
                output.update(timing(started))
                save_json(args.output, output)


if __name__ == '__main__':
    main()
