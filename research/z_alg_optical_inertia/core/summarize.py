"""Create concise tables and figures from the numerical scans and Arb receipts."""
import argparse
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from scan import RESULTS, ROOT, envelope, key, read_json, save_json

FIGURES = ROOT / 'figures'
ARMS = ('sector_PD', 'sector', 'global_plus', 'global_signed', 'spatial')
PLOT_ARMS = ('sector', 'global_plus', 'global_signed', 'spatial')
COLORS = ('#193c66', '#148878', '#8f55a3', '#cf693c')
LABELS = ('Sector only', 'Sector + optical (λ≥0)', 'Sector + real equality (signed λ)', '4×4 blocks')


def report_positions(paths):
    out = Counter()
    for path in paths:
        if path.is_file():
            out[str(path.name)] += 1
    return dict(out)


def load_numeric():
    raw = read_json(RESULTS / 'scalar.json')
    costs = {'scalar': raw['scalar_cost'], 'scalar_geometry': raw['geometries']}
    for name in ('spatial_archive', 'spatial_synthetic'):
        path = RESULTS / f'{name}.json'
        if not path.is_file():
            continue
        extension = read_json(path)
        costs[name] = extension['spatial_cost']
        for case, new in extension['cases'].items():
            for u, cap in new['caps'].items():
                if 'spatial' in cap:
                    raw['cases'][case]['caps'][u]['spatial'] = cap['spatial']
            if 'spatial_operations' in new:
                raw['cases'][case]['spatial_operations'] = new['spatial_operations']
    for name in ('focused_scalar', 'focused_upper'):
        path = RESULTS / f'{name}.json'
        if not path.is_file():
            continue
        extension = read_json(path)
        costs[name] = {k: extension[k] for k in ('wall_seconds', 'cpu_seconds', 'shared_A_basis_eigendecompositions')}
        for case, new in extension['cases'].items():
            rows = {'0.5': new['rows']} if 'rows' in new else new['caps']
            for u, points in rows.items():
                if not points:
                    continue
                cap = raw['cases'][case]['caps'][u]
                cap['scalar_rows'] += points
                cap['scalar_envelope'] = envelope(cap['scalar_rows'])
            raw['cases'][case].setdefault('focused_operations', {})[name] = new.get('operations', {})
    return raw['cases'], costs


def merge_gradient(cases, costs):
    numeric = RESULTS / 'gradient_blocks.json'
    local = RESULTS / 'gradient_scalar_control.json'
    if numeric.is_file():
        data = read_json(numeric)
        costs['gradient_blocks'] = {k: data[k] for k in ('wall_seconds', 'cpu_seconds')}
        for case, experiment in data['cases'].items():
            for u, changes in experiment['caps'].items():
                rows = changes['rows']
                if not rows:
                    continue
                interior = cases[case]['caps'][u]
                reference = changes['baseline']
                mu = reference['mu']
                scalars = [{**r, 'u': float(u), 'mu': mu, 'lambda_over_mu': r['lambda'] / mu}
                           for r in rows if r['method'] == 'same_count_scalar' and r['same_inertia']]
                interior['scalar_rows'] += scalars
                interior['scalar_envelope'] = envelope(interior['scalar_rows'])
                candidates = [{**r, 'u': float(u), 'mu': mu, 'block_mu': [mu] * 16}
                              for r in rows if r['method'] == 'gradient_blocks' and r['same_inertia']]
                if candidates:
                    interior['gradient_best'] = max(candidates, key=key)
                    interior['gradient_attempts'] = [r for r in rows if r['method'] == 'gradient_blocks']
    if local.is_file():
        data = read_json(local)
        costs['gradient_scalar_control'] = {k: data[k] for k in ('wall_seconds', 'cpu_seconds')}
        for case, row in data['cases'].items():
            interior = cases[case]['caps'][str(row['cap'])]
            interior['scalar_rows'] += row['numerical_grid']
            interior['scalar_envelope'] = envelope(interior['scalar_rows'])

    polished = RESULTS / 'gradient_scalar_polish.json'
    if polished.is_file():
        data = read_json(polished)
        costs['gradient_scalar_polish'] = {k: data[k] for k in ('wall_seconds', 'cpu_seconds')}
        for case, row in data['cases'].items():
            cap = str(row['best_scalar']['u'])
            interior = cases[case]['caps'][cap]
            interior['scalar_rows'] += row['evaluations']
            interior['scalar_envelope'] = envelope(interior['scalar_rows'])


def gradient_evidence():
    file = RESULTS / 'gradient_certificates.json'
    if not file.is_file():
        return None, {}
    source = read_json(file)
    output = {'cases': {}, 'proposal_freeze': {}}
    for case, row in source['cases'].items():
        values, receipts = [], []
        for arm in ('block', 'scalar'):
            trial = row[f'{arm}_proposal']
            is_block = arm == 'block'
            values.append({'u': row['cap'], 'mu': trial['mu'], 'lambda': trial['lambda'],
                           'block_lambda': trial.get('block_lambda', []) if is_block else [],
                           'block_mu': [trial['mu']] * 16 if is_block else [],
                           'tags': ['fixed_inertia_gradient_control']})
            receipts.append({**row[f'{arm}_receipt'], 'proposal_index': len(receipts)})
        output['proposal_freeze'][case] = {'proposals': values}
        output['cases'][case] = {'receipts': receipts}
    return output, {k: source[k] for k in ('wall_seconds', 'cpu_seconds')}


def arms_at(cap):
    current = cap['scalar_envelope'].copy()
    if 'spatial' in cap or 'gradient_best' in cap:
        # The spatial family also contains scalar points. Keep the block search receipt separately.
        row = cap.get('spatial', {}).get('best')
        current['spatial'] = max((row, cap.get('gradient_best'), current['global_signed']), key=key)
    return current


def number_rows(rows):
    sets = {(r['mu'], r['lambda']) for r in rows}
    return len(rows), len(sets)


def certs_for(case, evidence):
    unique = {}
    for doc in evidence:
        if case not in doc['cases']:
            continue
        for row in doc['cases'][case]['receipts']:
            proposal = doc['proposal_freeze'][case]['proposals'][row['proposal_index']]
            kp = (proposal['u'], proposal['mu'], proposal['lambda'], tuple(proposal['block_mu']),
                  tuple(proposal['block_lambda']))
            previous = unique.get(kp)
            if previous is None or row['status'] == 'Arb_verified':
                unique[kp] = {'proposal': proposal, **row}
    return list(unique.values())


def arm_mask(p, nu):
    if p['block_mu']:
        return ('spatial',)
    arms = ['global_signed', 'spatial']
    if p['lambda'] >= 0.:
        arms.append('global_plus')
    if p['lambda'] == 0.:
        arms.append('sector')
        if nu == 0:
            arms.append('sector_PD')
    return arms


def summarize_case(case, record, evidence):
    checks = certs_for(case, evidence)
    numeric = {arm: {} for arm in ARMS}
    strict = {arm: {} for arm in ARMS}
    points, unique = 0, 0
    block_rows, changed = 0, []
    for u, cap in record['caps'].items():
        u = float(u)
        current = arms_at(cap)
        count, unique_count = number_rows(cap['scalar_rows'])
        points += count
        unique += unique_count
        for arm, row in current.items():
            if arm in numeric and row is not None and row['lower_bound_numeric'] is not None:
                for r, value in record['radii'].items():
                    if row['lower_bound_numeric'] > value:
                        prior = numeric[arm].get(r, {})
                        if u > prior.get('u', -1.):
                            numeric[arm][r] = {'u': u, 'nu': row['nu'], 'lower': row['lower_bound_numeric'],
                                               'mu': row['mu'], 'lambda': row['lambda']}
        block_rows += len(cap.get('spatial', {}).get('evaluations', [])) + len(cap.get('gradient_attempts', []))
        if key(current.get('spatial'))[0] > key(current['global_signed'])[0]:
            changed.append({'u': u, 'scalar': key(current['global_signed'])[0],
                            'block': key(current['spatial'])[0]})
    for proof in checks:
        if proof['status'] != 'Arb_verified':
            continue
        p, nu = proof['proposal'], proof['verified_negative_eigenvalues']
        for arm in arm_mask(p, nu):
            for r, excluded in proof['strictly_excludes'].items():
                if excluded:
                    best = strict[arm].get(r, {})
                    item = {'u': p['u'], 'nu': nu, 'lower': proof['certified_frobenius_residual_lower'],
                            'mu': p['mu'], 'lambda': p['lambda'], 'margin_lower': proof['margin_lower'][r],
                            'coherent_single_lower_same_multipliers': proof.get('strongest_coherent_source_lower'),
                            'original_single_lower_same_multipliers': proof.get('strongest_original_source_lower'),
                            'coherent_single_upper_same_multipliers': proof.get('strongest_coherent_source_upper_same_multipliers'),
                            'original_single_upper_same_multipliers': proof.get('strongest_original_source_upper_same_multipliers'),
                            'radius': proof['radii'][r], 'block_mu': p['block_mu'], 'block_lambda': p['block_lambda']}
                    if (item['u'], item['lower']) > (best.get('u', -1.), best.get('lower', -1.)):
                        strict[arm][r] = item
    return {'numeric_discovered_caps': numeric, 'strict_discovered_caps': strict,
            'radii': record['radii'], 'numeric_rows': points, 'unique_scalar_points': unique,
            'spatial_trials': block_rows, 'block_improvements_over_sampled_scalar': changed,
            'strict_verification_receipts': len(checks), 'strict_unresolved_receipts':
            sum(c['status'] != 'Arb_verified' for c in checks),
            'operations': {k: record[k] for k in ('scalar_operations', 'spatial_operations', 'focused_operations')
                           if k in record}}


def plot_envelope(ax, record, title, xlim=None):
    caps = [(float(u), arms_at(cap)) for u, cap in record['caps'].items()]
    caps.sort(key=lambda c: c[0])
    radius = next(iter(record['radii'].values()))
    for arm, color, label, linestyle in zip(PLOT_ARMS, COLORS, LABELS, ('-', '--', ':', '-.'), strict=True):
        points = [(u, d[arm]['lower_bound_numeric'] / radius) for u, d in caps if arm in d and
                  d[arm] is not None and d[arm]['lower_bound_numeric'] is not None]
        us, values = zip(*points)
        ax.plot(us, values, marker='.', color=color, lw=1.45, ms=3, label=label, linestyle=linestyle)
    ax.axhline(1., color='#9c2a22', lw=.9, alpha=.8, label='Noise radius')
    ax.set_title(title, fontsize=9, loc='left', fontweight='normal')
    ax.set_yscale('symlog', linthresh=.15)
    ax.set_ylim(bottom=0.)
    if xlim:
        ax.set_xlim(xlim)
    ax.grid(True, alpha=.2)


def plot_envelopes(cases):
    archived = {k: v for k, v in cases.items() if v['metadata']['cohort'] == 'exposed_archive'}
    names = ('R10 ring | 48×12', 'Austria | 40×20', 'SRAP shell | half ring',
             'SRAP shell | full ring', 'SRAP double ring | half ring')
    f, axes = plt.subplots(2, 3, figsize=(11, 6.6), constrained_layout=True)
    for i, ((_, case), title) in enumerate(zip(archived.items(), names, strict=True)):
        plot_envelope(axes.flat[i], case, title, (0., 6.15) if i == 0 else (0., 1.55))
    handles, labels = axes.flat[0].get_legend_handles_labels()
    axes.flat[-1].axis('off')
    axes.flat[-1].legend(handles, labels, loc='center', frameon=False)
    f.supxlabel('Hypothesized upper contrast u')
    f.supylabel('Numerical lower residual / archive noise radius')
    f.suptitle('Frozen archives · multiplier envelopes', fontweight='bold')
    f.savefig(FIGURES / 'core_archive_envelopes.png', dpi=185)
    plt.close(f)

    f, axes = plt.subplots(2, 3, figsize=(11, 6.6), constrained_layout=True)
    for i in range(6):
        obs = f'optfresh_{i+1:02d}__pw180_10x20__n02'
        plot_envelope(axes.flat[i], cases[obs], f'Fresh object {i+1} | half ring | 2% noise', (0., 1.15))
    handles, labels = axes.flat[0].get_legend_handles_labels()
    f.legend(handles, labels, loc='outside lower center', ncol=3, frameon=False, fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel('Hypothesized upper contrast u')
    f.supylabel('Numerical lower residual / noise-only radius')
    f.suptitle('Six independent cross-grid materials · half-ring numerical scan', fontweight='bold')
    f.savefig(FIGURES / 'core_fresh_envelopes.png', dpi=185)
    plt.close(f)


def plot_inertia(cases):
    selectors = [('r10_twod_0001', 5., .001, 'R10 ring'),
                 ('austria_eps4_400MHz', .75, .01, 'Austria'),
                 ('s2d_007__pw180_10x20', .5, .01, 'SRAP shell · half ring'),
                 ('s2d_007__pw360_10x20', .75, .001, 'SRAP shell · full ring')]
    f, axes = plt.subplots(2, 2, figsize=(10.5, 6.1), constrained_layout=True)
    for ax, (case, u, mu, name) in zip(axes.flat, selectors, strict=True):
        record = cases[case]
        cap = record['caps'][str(u)]
        rows = {x['lambda_over_mu']: x for x in cap['scalar_rows'] if x['mu'] == mu and
                x['nu'] is not None and x['lambda_over_mu'] >= 0.}
        ratios = sorted(rows)
        nu = [rows[t]['nu'] for t in ratios]
        radius = next(iter(record['radii'].values()))
        bound = [rows[t]['lower_bound_numeric'] / radius for t in ratios]
        h, = ax.step(ratios, nu, where='mid', c=COLORS[0], lw=1.8, label='Negative directions ν')
        h2 = ax.axhline(record['receiver_source_shape'][1], ls=':', c=COLORS[0], lw=.8, label='Source budget s')
        ax2 = ax.twinx()
        k, = ax2.plot(ratios, bound, '.-', c=COLORS[-1], lw=1.4, label='Lower residual / noise radius')
        k2 = ax2.axhline(1., ls=':', c=COLORS[-1], lw=.8, label='Noise boundary')
        ax.set_xscale('symlog', linthresh=.1)
        ax.set_xlim(-.005, 1.1e4)
        ax.set_title(f'{name} | u={u:g}, μ={mu:g}', fontsize=9, loc='left')
        ax.grid(True, alpha=.2)
        ax.set_ylim(bottom=0.)
        ax2.set_ylim(bottom=0.)
    f.legend((h, h2, k, k2), (h.get_label(), h2.get_label(), k.get_label(), k2.get_label()),
             loc='outside lower center', ncol=4, frameon=False, fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel('Positive optical ratio λ/μ')
    f.suptitle('Negative directions and residual margins', fontweight='bold')
    f.savefig(FIGURES / 'core_inertia_curves.png', dpi=185)
    plt.close(f)


def plot_dark():
    data = read_json(RESULTS / 'dark_thresholds.json')
    f, axes = plt.subplots(1, 2, figsize=(9.6, 3.95), constrained_layout=True,
                           gridspec_kw={'width_ratios': (1., 2.4)})
    items = (('r10_twod_0001', 'R10 ring', axes[0]), ('austria_eps4_400MHz', 'Austria', axes[1]),
             ('s2d_007__pw360_10x20', '20-site full ring', axes[1]),
             ('s2d_007__pw180_10x20', '20-site half ring', axes[1]))
    for i, (case, label, ax) in enumerate(items):
        points = [(r['lambda_over_mu'], r['u_star_numeric']) for r in data['geometries'][case]['rows']
                  if r['lambda_over_mu'] >= 0. and r['u_star_numeric'] is not None]
        t, y = zip(*points)
        ax.plot(t, y, '.-', c=COLORS[i], label=label)
        ax.set_xscale('symlog', linthresh=.1)
        ax.set_xlim(-.005, 1.1e4)
        ax.grid(True, alpha=.22)
    axes[0].set_ylim(4.10, 4.13)
    axes[1].set_ylim(.3, .78)
    for ax in axes:
        ax.legend(loc='lower right', frameon=False, fontsize=8)
    f.supxlabel('Positive optical ratio λ/μ')
    f.supylabel('Numerical dark threshold u*')
    f.suptitle('Receiver geometry sets the room for scalar optical braking', fontweight='bold')
    f.savefig(FIGURES / 'core_dark_thresholds.png', dpi=185)
    plt.close(f)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificates', type=Path, action='append')
    args = parser.parse_args()
    cases, costs = load_numeric()
    paths = args.certificates or [RESULTS / name for name in ('certificates.json', 'certificates_focused_upper.json')]
    evidence = []
    for path in paths:
        if path.is_file():
            block = read_json(path)
            evidence.append(block)
            costs[path.stem] = {k: block[k] for k in ('wall_seconds', 'cpu_seconds')}
    for name in ('dark_thresholds', 'single_source_obstruction', 'single_source_obstruction_block_case'):
        path = RESULTS / f'{name}.json'
        if path.is_file():
            block = read_json(path)
            costs[name] = {k: block[k] for k in ('wall_seconds', 'cpu_seconds')}
    merge_gradient(cases, costs)
    grade, grade_cost = gradient_evidence()
    if grade:
        evidence.append(grade)
        costs['gradient_certificates'] = grade_cost
        paths.append(RESULTS / 'gradient_certificates.json')
    result = {'input_receipts': report_positions(paths), 'costs': costs, 'cases': {}}
    for case, record in cases.items():
        result['cases'][case] = summarize_case(case, record, evidence)
    save_json(RESULTS / 'summary.json', result)
    plot_envelopes(cases)
    plot_inertia(cases)
    plot_dark()


if __name__ == '__main__':
    main()
