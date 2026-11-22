"""Freeze observed dual proposals before opening cross-grid proxy error data, then certify."""
import argparse
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from certificate import certify
from inertia import clock, timing
from scan import RESULTS, SYNTHETIC, block_ids, envelope, key, load_cases, read_json, save_json

ARM_NAMES = ('sector_PD', 'sector', 'global_plus', 'global_signed', 'spatial')


def rows_by_arm(record, name):
    for cap in record['caps'].values():
        row = cap.get('spatial', {}).get('best') if name == 'spatial' else cap['scalar_envelope'][name]
        if row is not None and row['lower_bound_numeric'] is not None:
            yield row


def proposals(record):
    frozen = {}

    def add(row, tag):
        kp = (row['u'], row['mu'], row['lambda'], tuple(row.get('block_mu', [])),
              tuple(row.get('block_lambda', [])))
        if 'block_mu' in row and len(set(row['block_mu'])) == len(set(row['block_lambda'])) == 1:
            kp = (row['u'], row['block_mu'][0], row['block_lambda'][0], (), ())
        if kp not in frozen:
            u, m, t, ms, ts = kp
            frozen[kp] = {'u': u, 'mu': m, 'lambda': t, 'nu_numeric': row['nu'],
                          'numeric_lower': row['lower_bound_numeric'],
                          'block_mu': ms, 'block_lambda': ts, 'tags': []}
        if tag not in frozen[kp]['tags']:
            frozen[kp]['tags'].append(tag)

    baselines = {}
    for arm in ARM_NAMES:
        baselines[arm] = {}
        for name, radius in record['radii'].items():
            feasible = [row for row in rows_by_arm(record, arm) if row['lower_bound_numeric'] > radius]
            if not feasible:
                continue
            chosen = max(feasible, key=lambda x: (x['u'], x['lower_bound_numeric']))
            baselines[arm][name] = chosen
            # Record the sector reference and the strongest scalar/spatial endpoints. Other
            # arms need independent verification only when they add an exclusion cap.
            sector = baselines.get('sector', {}).get(name)
            sector_u = sector['u'] if sector is not None else -1.
            global_best_u = max((baselines.get(a, {}).get(name, {}).get('u', -1.)
                                 for a in ('global_plus', 'global_signed')), default=-1.)
            if arm == 'sector' or (arm in ('global_plus', 'global_signed') and chosen['u'] > sector_u):
                add(chosen, f'largest_discovered_{arm}_exclusion__{name}')
            if arm == 'sector_PD' and record['metadata']['cohort'] == 'exposed_archive':
                add(chosen, f'largest_discovered_{arm}_exclusion__{name}')
            if arm == 'spatial' and chosen['u'] > global_best_u:
                add(chosen, f'largest_discovered_{arm}_exclusion__{name}')

    case = record['metadata'].get('id', record['metadata'].get('case'))
    presentations = {
        'r10_twod_0001': (5., .001, 3.),
        'austria_eps4_400MHz': (.75, .01, 30.),
        's2d_007__pw180_10x20': (.5, .01, 1.),
        's2d_010__pw180_10x20': (.5, .0001, 1.),
        's2d_007__pw360_10x20': (.75, .001, 10.),
        'optfresh_01__pw180_10x20__n02': (.5, .01, 1.),
        'optfresh_01__pw360_10x20__n02': (.75, .001, 10.),
    }
    if case in ('austria_eps4_400MHz', 'optfresh_01__pw360_10x20__n02',
                'optfresh_05__pw360_10x20__n02'):
        comparison = record['caps']['0.75']['scalar_envelope']
        add(comparison['sector'], 'focused_sector_margin_reference')
        add(comparison['global_signed'], 'focused_signed_margin_reference')
    if case in presentations:
        u, mu, ratio = presentations[case]
        cap = record['caps'][str(u)]
        for row in cap['scalar_rows']:
            if row['mu'] == mu and row['lambda_over_mu'] in (0., ratio):
                add(row, 'fixed_mu_inertia_pair')
        if cap.get('spatial'):
            s = cap['spatial']['best']
            g = cap['scalar_envelope']['global_signed']
            if key(s)[0] > key(g)[0]:
                add(s, 'spatial_residual_margin_reference')
                add(g, 'scalar_residual_margin_reference')
    return list(frozen.values()), {arm: {r: v['u'] for r, v in active.items()}
                                   for arm, active in baselines.items()}


def comparison_radii(case, radii):
    if 'noise_only_fro' not in radii:
        return radii, {}
    eval_file = SYNTHETIC.parent / 'evaluation' / f'{case}.npz'
    with np.load(eval_file, allow_pickle=False) as offline:
        discrepancy = float(offline['discrepancy_fro'])
        direct = float(offline['proxy_error_fro'])
    return {**radii, 'noise_plus_this_proxy_discrepancy': float(np.nextafter(
        discrepancy + radii['noise_only_fro'], np.inf)), 'exact_this_proxy_error_upper': direct}, {
        'scope': 'the saved candidate proxy only; no uniform discretization envelope',
        'coarse_proxy_discrepancy_upper': discrepancy, 'coarse_proxy_error_upper': direct,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=RESULTS / 'scalar.json')
    parser.add_argument('--spatial', type=Path, action='append', help='merge corresponding frozen spatial proposals')
    parser.add_argument('--focused', type=Path, action='append', help='merge frozen focused scalar controls')
    parser.add_argument('--cohort', choices=('archive', 'synthetic', 'all'), default='all')
    parser.add_argument('--case', action='append', help='only prove selected released observations')
    parser.add_argument('--output', type=Path, default=RESULTS / 'certificates.json')
    parser.add_argument('--prior', type=Path, action='append', help='reuse prior verified receipts for identical proposals')
    parser.add_argument('--precision', type=int, default=96)
    parser.add_argument('--threads', type=int, default=2)
    args = parser.parse_args()
    scan = read_json(args.input)
    for path in args.spatial or ():
        for case, row in read_json(path)['cases'].items():
            for cap, results in row['caps'].items():
                if 'spatial' in results and case in scan['cases']:
                    scan['cases'][case]['caps'][cap]['spatial'] = results['spatial']
    for path in args.focused or ():
        for case, record in read_json(path)['cases'].items():
            if case not in scan['cases']:
                continue
            cap_rows = {str(record['rows'][0]['u']): record['rows']} if 'rows' in record else record['caps']
            for cap, rows in cap_rows.items():
                if not rows:
                    continue
                entry = scan['cases'][case]['caps'][cap]
                entry['scalar_rows'] += rows
                entry['scalar_envelope'] = envelope(entry['scalar_rows'])
    prior = [read_json(path) for path in args.prior or ()]
    started = clock()
    output = {'cohort': args.cohort, 'selection_requires_observed_arrays_only': True,
              'proposal_freeze': {}, 'cases': {}}
    with threadpool_limits(limits=1, user_api='blas'):
        cases = {case: data for group in load_cases(args.cohort).values() for case, data in group.items()}
        for case in cases:
            if case not in scan['cases'] or (args.case and case not in args.case):
                continue
            selected, caps = proposals(scan['cases'][case])
            output['proposal_freeze'][case] = {'proposals': selected, 'numeric_discovered_caps': caps}
        save_json(args.output, output)
        for case in cases:
            if case not in output['proposal_freeze']:
                continue
            p, radii, positions, _ = cases[case]
            case_proposals = output['proposal_freeze'][case]['proposals']
            radii, scope = comparison_radii(case, radii)
            receipts = []
            result = {'radii': radii, 'proxy_evaluation': scope, 'receipts': receipts}
            output['cases'][case] = result
            for idx, row in enumerate(case_proposals):
                u, mu, lam = row['u'], row['mu'], row['lambda']
                order = block_ids(positions)
                if row['block_mu']:
                    mu, lam = np.asarray(row['block_mu'])[order], np.asarray(row['block_lambda'])[order]
                proof = None
                if prior:
                    coefficients = ('u', 'mu', 'lambda', 'block_mu', 'block_lambda')
                    to_compare = [row[x] if not isinstance(row[x], tuple) else list(row[x]) for x in coefficients]
                    for report in prior:
                        if case not in report['cases'] or report['cases'][case]['radii'] != radii:
                            continue
                        for old in report['cases'][case]['receipts']:
                            old_proposal = report['proposal_freeze'][case]['proposals'][old['proposal_index']]
                            if [old_proposal[x] for x in coefficients] == to_compare and old['status'] == 'Arb_verified':
                                proof = {**old, 'reused_verified_receipt': True}
                                break
                        if proof is not None:
                            break
                if proof is None:
                    proof = certify(p, u, mu, lam, radii=radii, precision=args.precision, threads=args.threads)
                if proof['status'] == 'Arb_unresolved':
                    retry = certify(p, u, mu, lam, radii=radii, precision=max(160, args.precision),
                                    threads=args.threads)
                    proof = {**retry, 'first_attempt': proof}
                receipts.append({**proof, 'proposal_index': idx})
                print(case, u, row['tags'], proof.get('verified_negative_eigenvalues'),
                      proof.get('certified_frobenius_residual_lower'), proof.get('strictly_excludes'), flush=True)
                output.update(timing(started))
                save_json(args.output, output)


if __name__ == '__main__':
    main()
