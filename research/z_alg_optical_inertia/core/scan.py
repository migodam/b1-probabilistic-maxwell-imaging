"""Scan scalar and 4-by-4 block multipliers on released observations only."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigvalsh
from threadpoolctl import threadpool_limits

from inertia import Prepared, SectorBasis, clock, timing

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ARCHIVE = ROOT.parent / 'z_alg_secant_tasks' / 'highdim'
SYNTHETIC = ROOT / 'synthetic' / 'public'
RESULTS = ROOT / 'results' / 'core'
CAPS = (.1, .25, .5, .75, 1., 1.25, 1.5, 2., 2.5, 3., 3.5, 4.,
        4.25, 4.5, 4.75, 5., 5.25, 5.5, 6.)
ARCHIVE_CAPS = {'r10_twod_0001': CAPS, 'others': CAPS[:12]}
SPATIAL_CAPS = (.25, .5, .75, 1., 4., 5.)
MUS = (1.e-4, 1.e-3, 1.e-2, 1.e-1, 1., 10.)
T_GRID = (0., .03, -.03, .1, -.1, .3, -.3, 1., -1., 3., -3., 10., -10.,
          30., -30., 100., -100., 300., -300., 1000., -1000., 10_000.)


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save_json(path, data):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(data, allow_nan=False, separators=(',', ':')) + '\n', encoding='utf-8')
    temp.replace(path)


def load_cases(cohort):
    cases = defaultdict(dict)
    if cohort in ('archive', 'all'):
        sys.path.insert(0, str(ARCHIVE))
        from replay import OBSERVATION_RADII, OPERATORS, SOURCE, load_operators

        with np.load(OPERATORS, allow_pickle=False) as operators:
            for case, radii in read_json(OBSERVATION_RADII).items():
                case_path = SOURCE / f'{case}.npz'
                doc = read_json(case_path.with_suffix('.json'))
                geom = 'archived/' + operators[f'{case}__G_key'].item()
                model = load_operators(case)
                with np.load(case_path, allow_pickle=False) as ob:
                    problem = Prepared(model['G'], model['S'], model['E0'], ob['Y_archive'])
                    cases[geom][case] = (problem, radii, ob['positions'], {
                        'cohort': 'exposed_archive', 'positions_file': str(case_path.relative_to(ROOT.parent)),
                        'case': case, 'acquisition': doc['acquisition'], 'archived_material_box': doc['material_box'],
                        'tested_caps': ARCHIVE_CAPS.get(case, ARCHIVE_CAPS['others']),
                    })
    if cohort in ('synthetic', 'all'):
        index = read_json(SYNTHETIC / 'manifest.json')
        with np.load(SYNTHETIC / index['green_file'], allow_pickle=False) as data:
            G = data['G']
            positions = data['positions'] if 'positions' in data else data['centers']
        geometry = {}
        for name, row in index['geometries'].items():
            with np.load(SYNTHETIC / row['file'], allow_pickle=False) as source:
                geometry[name] = source['S'], source['B']
        for row in index['observations']:
            S, B = geometry[row['geometry']]
            with np.load(SYNTHETIC / row['file'], allow_pickle=False) as source:
                problem = Prepared(G, S, B, source['Y'])
            cases['cross_grid/' + index['green_file']][row['id']] = (
                problem, {'noise_only_fro': row['radius_noise_fro']}, positions,
                {'cohort': 'fresh_cross_grid', **row,
                 'tested_caps': tuple(u for u in CAPS if u <= row['upper_contrast'])},
            )
    return cases


def block_ids(positions):
    start = positions.min(axis=0)
    frac = (positions - start) / np.ptp(positions, axis=0)
    grid = np.minimum((frac * 4.).astype(int), 3)
    return 4 * grid[:, 1] + grid[:, 0]


def key(row):
    if row is None:
        return -1., -np.inf
    lower = row['lower_bound_numeric']
    return -1. if lower is None else lower, -np.inf if row['nu'] is None else -row['nu']


def best(rows, predicate):
    return max((row for row in rows if predicate(row)), key=key, default=None)


def compact(row, *, h=False):
    drop = {'wall_seconds', 'cpu_seconds'}
    if not h:
        drop.add('H_eigenvalues_descending')
    return {k: v for k, v in row.items() if k not in drop}


def envelope(rows):
    return {name: best(rows, predicate) for name, predicate in (
        ('sector_PD', lambda r: r['lambda'] == 0. and r['nu'] == 0),
        ('sector', lambda r: r['lambda'] == 0.),
        ('global_plus', lambda r: r['lambda'] >= 0.),
        ('global_signed', lambda r: True),
        ('strictly_positive_lambda', lambda r: r['lambda'] > 0.),
        ('strictly_negative_lambda', lambda r: r['lambda'] < 0.),
    )}


def scalar_scan(geometry, caps, ratios, mus):
    initial = next(iter(geometry.values()))[0]
    start = clock()
    diagnostics = eigvalsh(initial.Ga.real if not np.any(initial.Ga.imag) else initial.Ga)
    shared = {'green_symmetry_defect_max': float(np.max(np.abs(initial.G - initial.G.T))),
              'ImH_min_eigenvalue_numeric': float(diagnostics[0]),
              'ImH_max_eigenvalue_numeric': float(diagnostics[-1]),
              'ImH_negative_count_raw_numeric': int(np.count_nonzero(diagnostics < 0.)),
              'ImH_count_larger_1e-12_times_max_numeric': int(np.count_nonzero(
                  diagnostics > diagnostics[-1] * 1.e-12)),
              'geometry_spectrum_cost': timing(start)}
    results = {case: {} for case in geometry}
    basis_count, basis_wall, basis_cpu = 0, 0., 0.
    for u in caps:
        for ratio in ratios:
            start = clock()
            basis = SectorBasis(initial.G, u, ratio)
            basis_count += 1
            spent = timing(start)
            basis_wall += spent['wall_seconds']
            basis_cpu += spent['cpu_seconds']
            for case, (p, _, _, meta) in geometry.items():
                if u not in meta['tested_caps']:
                    continue
                rows = results[case].setdefault(u, [])
                for mu, value in zip(mus, basis.scan(p, mus), strict=True):
                    rows.append({**value, 'u': u, 'mu': float(mu), 'lambda': float(mu * ratio),
                                 'lambda_over_mu': ratio})
    shared['A_basis_eigendecompositions'] = basis_count
    shared['A_basis_cost'] = {'wall_seconds': basis_wall, 'cpu_seconds': basis_cpu}
    return results, shared


def spatial_proposals(problem, positions, u, global_rows, *, coordinate):
    """Coordinate searches use the observed dual objective alone."""
    started = clock()
    available = [global_rows['global_signed'], global_rows['global_plus'], global_rows['sector']]
    order = block_ids(positions)
    uniform = np.ones(16)
    seed = max(available, key=key)
    ref_mu, ref_lam = seed['mu'], seed['lambda']
    mu_values, lam_values = ref_mu * uniform, ref_lam * uniform
    incumbent = dict(seed, block_mu=mu_values.tolist(), block_lambda=lam_values.tolist())
    evaluations = []

    def evaluate(mu, lam, label):
        nonlocal mu_values, lam_values, incumbent
        proposal = problem.evaluate_ldl(u, mu[order], lam[order])
        row = {**proposal, 'u': u, 'mu': ref_mu, 'lambda': ref_lam,
               'block_mu': mu.tolist(), 'block_lambda': lam.tolist(), 'proposal': label}
        evaluations.append(compact(row))
        if key(row) > key(incumbent):
            mu_values, lam_values, incumbent = mu, lam, row

    # Four balanced geometry-only directions give a bounded affordable first pass.
    ids = np.arange(16)
    directions = (ids % 4 >= 2, ids // 4 >= 2, np.isin(ids, [5, 6, 9, 10]),
                  ((ids % 4) // 2 + (ids // 8)) == 1)
    for mask in directions:
        sign = 2. * mask - 1.
        for direction in (-1., 1.):
            evaluate(mu_values * np.exp(.25 * direction * sign), lam_values, 'balanced_mu')
            if coordinate:
                evaluate(mu_values, lam_values + .25 * max(abs(ref_lam), ref_mu) * direction * sign,
                         'balanced_lambda')
    if coordinate:
        for block in range(16):
            for factor, direction in ((.8, -1.), (1.25, 1.)):
                mu, lam = mu_values.copy(), lam_values.copy()
                mu[block] *= factor
                lam[block] += .25 * max(abs(ref_lam), ref_mu) * direction
                evaluate(mu, lam, 'block_mu_and_lambda')
    return {'best': compact(incumbent, h=True), 'evaluations': evaluations, 'coordinate_sweep': coordinate,
            'numeric_only': True, **timing(started)}


def spatial_export(full):
    """Retain only the spatial additions when using an existing scalar scan."""
    cases = {}
    for case, row in full['cases'].items():
        caps = {u: {'spatial': cap['spatial']} for u, cap in row['caps'].items() if 'spatial' in cap}
        if caps:
            cases[case] = {'caps': caps, 'spatial_operations': row.get('spatial_operations', {})}
    return {'numeric_only': True, 'spatial_cost': full['spatial_cost'], 'cases': cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cohort', choices=('archive', 'synthetic', 'all'), default='all')
    parser.add_argument('--case', action='append', help='only run one released observation; repeat as needed')
    parser.add_argument('--output', type=Path, default=RESULTS / 'scan.json')
    parser.add_argument('--quick', action='store_true', help='quick diagnostic subset of ratios and caps')
    parser.add_argument('--spatial', action='store_true', help='add bounded block-coordinate search')
    parser.add_argument('--scalar-input', type=Path, help='reuse existing scalar rows for the spatial stage')
    args = parser.parse_args()
    start = clock()
    output = read_json(args.scalar_input) if args.scalar_input else {
        'numeric_only': True, 'cohort': args.cohort, 'quick': args.quick,
        'cases': {}, 'geometries': {}, 'mu_grid': MUS,
        'ratios': (0., -10., -1., 1., 10., 100.) if args.quick else T_GRID,
    }
    with threadpool_limits(limits=1, user_api='blas'):
        cases = load_cases(args.cohort)
        for geom, group in cases.items():
            if args.case:
                group = {case: data for case, data in group.items() if case in args.case}
            if not group:
                continue
            caps = sorted(set().union(*(data[3]['tested_caps'] for data in group.values())))
            if args.quick:
                caps = [u for u in caps if u in (.1, .5, 1., 2., 4., 5., 5.5)]
            if not args.scalar_input:
                data, shared = scalar_scan(group, caps, output['ratios'], MUS)
                output['geometries'][geom] = shared
                for case, (p, radii, _, meta) in group.items():
                    record = {'metadata': meta, 'geometry_of_G': geom, 'state_dimension': p.n,
                              'receiver_source_shape': list(p.Y.shape), 'radii': radii,
                              'observed_Frobenius': float(np.linalg.norm(p.Y)),
                              'scalar_operations': dict(p.operations), 'caps': {}}
                    for u, rows in data[case].items():
                        record['caps'][str(u)] = {'scalar_envelope': envelope(rows),
                                                   'scalar_rows': [compact(row) for row in rows]}
                    output['cases'][case] = record
                    highest = {arm: max((float(u) for u, d in record['caps'].items()
                                         if d['scalar_envelope'][arm] is not None
                                         and key(d['scalar_envelope'][arm])[0] > min(radii.values())), default=None)
                               for arm in ('sector', 'global_plus', 'global_signed')}
                    print(case, 'numeric exclusion caps', highest, flush=True)
            if args.spatial:
                for case, (p, radii, positions, _) in group.items():
                    record = output['cases'][case]
                    initial = p.operations.copy()
                    active_caps = set(SPATIAL_CAPS)
                    # Extra proposals at the first grid point above the best scalar exclusion.
                    cap_values = sorted(float(u) for u in record['caps'])
                    qualifying = [u for u in cap_values if key(record['caps'][str(u)][
                        'scalar_envelope']['global_signed'])[0] > min(radii.values())]
                    coordinate_cap = min(cap_values)
                    if qualifying:
                        top = max(qualifying)
                        coordinate_cap = next((u for u in cap_values if u > top), top)
                        active_caps.update((top, coordinate_cap))
                    for k, row in record['caps'].items():
                        cap = float(k)
                        if cap in active_caps:
                            row['spatial'] = spatial_proposals(p, positions, cap, row['scalar_envelope'],
                                                               coordinate=(cap == coordinate_cap and not args.quick))
                    record['spatial_operations'] = dict(p.operations - initial)
                    print(case, 'spatial search complete', flush=True)
                    output['spatial_cost'] = timing(start)
                    save_json(args.output, spatial_export(output) if args.scalar_input else output)
            output['scalar_cost' if not args.scalar_input else 'spatial_cost'] = timing(start)
            save_json(args.output, spatial_export(output) if args.scalar_input else output)


if __name__ == '__main__':
    main()
