"""Verify selected block candidates, and profile a finer local scalar control."""
import argparse
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from inertia import SectorBasis, clock, timing
from scan import RESULTS, block_ids, compact, key, load_cases, read_json, save_json

CASES = ('optfresh_03__pw180_10x20__n10', 'optfresh_04__pw360_10x20__n02')
CAPS = (.5, .75)


def selected(grad, case, u, source):
    row = grad['cases'][case]['caps'][str(u)]
    subset = [r for r in row['rows'] if r['method'] == source and r['same_inertia']]
    return row['baseline'], max(subset, key=key, default=None)


def tighter_control():
    start = clock()
    output = {'scope': 'observed-only local scalar control around the two numerical block crossings', 'cases': {}}
    problems = next(iter(load_cases('synthetic').values()))
    grad = read_json(RESULTS / 'gradient_blocks.json')
    for (case, u), rg in zip(zip(CASES, CAPS, strict=True), ((.32, .48, 65), (-.08, .2, 57)), strict=True):
        p, radii, _, _ = problems[case]
        base, block = selected(grad, case, u, 'gradient_blocks')
        ratios = np.linspace(*rg)
        mus = base['mu'] * np.geomspace(.4, 2.5, 41)
        rows = []
        for t in ratios:
            sb = SectorBasis(p.G, u, t)
            for mu, result in zip(mus, sb.scan(p, mus), strict=True):
                rows.append({**compact(result), 'u': u, 'mu': mu, 'lambda': mu * t, 'lambda_over_mu': t})
        winner = max(rows, key=key)
        output['cases'][case] = {'cap': u, 'radii': radii, 'numerical_grid': rows,
                                 'block': block, 'best_scalar': winner, 'basis_decompositions': len(ratios),
                                 'operations': dict(p.operations)}
        print(case, 'refined scalar', winner, 'block', key(block)[0], 'radius', radii, flush=True)
    output.update(timing(start))
    save_json(RESULTS / 'gradient_scalar_control.json', output)


def polish():
    from scipy.optimize import minimize

    started = clock()
    output = {'scope': 'deterministic two-parameter scalar optimization seeded by the local control grid',
              'cases': {}}
    data = next(iter(load_cases('synthetic').values()))
    first = read_json(RESULTS / 'gradient_scalar_control.json')
    for case, record in first['cases'].items():
        p, radii, *_ = data[case]
        cap, base = record['cap'], record['best_scalar']
        start = np.array([np.log(base['mu']), base['lambda_over_mu']])
        evaluations = []

        def objective(x):
            mu, t = float(np.exp(x[0])), float(x[1])
            row = p.evaluate_ldl(cap, mu, mu * t)
            result = {**compact(row), 'u': cap, 'mu': mu, 'lambda': mu * t, 'lambda_over_mu': t}
            evaluations.append(result)
            return -key(result)[0]

        solution = minimize(objective, start, method='Nelder-Mead',
                            bounds=((start[0] - 1., start[0] + 1.), (start[1] - .1, start[1] + .1)),
                            options={'maxfev': 60, 'xatol': 1.e-7, 'fatol': 1.e-9,
                                     'initial_simplex': np.array([start, start + [.04, 0.], start + [0., .003]])})
        best = max(evaluations, key=key)
        output['cases'][case] = {'best_scalar': best, 'evaluations': evaluations, 'radii': radii,
                                 'optimizer_stop': str(solution.message), 'operations': dict(p.operations)}
        print(case, best['lower_bound_numeric'], best['nu'], radii, solution.message, flush=True)
    output.update(timing(started))
    save_json(RESULTS / 'gradient_scalar_polish.json', output)


def verify():
    from certificate import certify
    from verify import comparison_radii

    start = clock()
    output = {'scope': 'block candidates preselected solely by numerical observed objective and noise-only radius',
              'cases': {}}
    grad = read_json(RESULTS / 'gradient_blocks.json')
    data = {case: row for group in load_cases('all').values() for case, row in group.items()}
    cases = dict(zip(CASES, CAPS, strict=True))
    cases['r10_twod_0001'] = 5.
    for case, u in cases.items():
        base, block = selected(grad, case, u, 'gradient_blocks')
        _, scalar = selected(grad, case, u, 'same_count_scalar')
        p, radii, pos, _ = data[case]
        radii, scope = comparison_radii(case, radii)
        lam = np.asarray(block['block_lambda'])[block_ids(pos)]
        control = max((base, scalar), key=key)
        block_proof = certify(p, u, base['mu'], lam, radii=radii)
        ctrl_proof = certify(p, u, base['mu'], control['lambda'], radii=radii)
        output['cases'][case] = {'cap': u, 'block_proposal': {'u': u, 'mu': base['mu'], **block},
                                 'block_receipt': block_proof, 'scalar_proposal': {'u': u, 'mu': base['mu'], **control},
                                 'scalar_receipt': ctrl_proof, 'radii': radii, 'proxy_scope': scope}
        output.update(timing(start))
        save_json(RESULTS / 'gradient_certificates.json', output)
        print(case, u, 'block', block_proof.get('certified_frobenius_residual_lower'),
              block_proof.get('strictly_excludes'), 'scalar', ctrl_proof.get('certified_frobenius_residual_lower'),
              ctrl_proof.get('strictly_excludes'), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('refine', 'polish', 'verify'))
    args = parser.parse_args()
    with threadpool_limits(limits=1, user_api='blas'):
        if args.mode == 'refine':
            tighter_control()
        elif args.mode == 'polish':
            polish()
        else:
            verify()


if __name__ == '__main__':
    main()
