"""Focused scalar controls on the archived and fresh ring apertures."""
import argparse
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from inertia import SectorBasis, clock, timing
from scan import RESULTS, compact, envelope, load_cases, save_json

CAP = .5
RATIOS = (-3., -2., -1., -.5, -.2, -.1, -.05, 0., .05, .1, .15, .2, .25, .3, .35,
          .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, .9, .95, 1., 1.05, 1.1,
          1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2., 2.5, 3., 4., 5., 7., 10.)
MUS = tuple(map(float, np.logspace(-5., 0., 31)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=RESULTS / 'focused_scalar.json')
    parser.add_argument('--caps', type=float, nargs='+', default=(CAP,))
    parser.add_argument('--array', choices=('half', 'all'), default='half')
    args = parser.parse_args()
    started = clock()
    output = {'caps': args.caps, 'array': args.array, 'ratio_grid': RATIOS, 'mu_grid': MUS, 'cases': {}}
    with threadpool_limits(limits=1, user_api='blas'):
        groups = load_cases('all')
        archive_group = groups['archived/G_1']
        synth_group = next(group for name, group in groups.items() if name.startswith('cross_grid/'))
        selected = {key: data for key, data in (archive_group | synth_group).items()
                    if args.array == 'all' or 'pw180' in key}
        operator_classes = []
        for name, (p, radii, _, meta) in selected.items():
            group = next((members for g, members in operator_classes if np.array_equal(g, p.G)), None)
            if group is None:
                group = {}
                operator_classes.append((p.G, group))
            group[name] = p
            output['cases'][name] = {'radii': radii, 'cohort': meta['cohort'],
                                     'caps': {str(u): [] for u in args.caps}}
        for G, group in operator_classes:
            for u in args.caps:
                for t in RATIOS:
                    basis = SectorBasis(G, u, t)
                    for name, problem in group.items():
                        for mu, row in zip(MUS, basis.scan(problem, MUS), strict=True):
                            output['cases'][name]['caps'][str(u)].append(compact({
                                **row, 'u': u, 'mu': mu, 'lambda': mu * t, 'lambda_over_mu': t}))
                for name in group:
                    x = output['cases'][name]
                    x.setdefault('envelope', {})[str(u)] = envelope(x['caps'][str(u)])
                    best = x['envelope'][str(u)]['global_signed']
                    print(name, u, best['lower_bound_numeric'], best['nu'], best['mu'], best['lambda_over_mu'],
                          flush=True)
                output['shared_A_basis_eigendecompositions'] = len(RATIOS) * len(operator_classes) * len(args.caps)
                output.update(timing(started))
                save_json(args.output, output)
            for name, problem in group.items():
                output['cases'][name]['operations'] = dict(problem.operations)
        output.update(timing(started))
        save_json(args.output, output)


if __name__ == '__main__':
    main()
