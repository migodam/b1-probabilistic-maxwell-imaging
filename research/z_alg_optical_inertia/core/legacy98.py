"""Pair the historical 180 sector points with their frozen positive-optical controls."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
OLD = RESEARCH / 'z_alg_secant_tasks/highdim/sector_inertia/results.json'
NEW = HERE.parent / 'results/core/scalar.json'
OUT = HERE.parent / 'results/core/legacy98.json'


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def summarize(rows, radii):
    exhausted = [r for r in rows if r['old_nu_exhausted_sources']]
    return {'old_grid_points': len(rows), 'old_exhausted_points': len(exhausted),
            'recovered_below_source_count': sum(r['recovered_below_source_count'] for r in exhausted),
            'positive_lower_bound': sum(r['positive_lower_bound'] for r in exhausted),
            'exceeds_any_old_radius': sum(any(r['exceeds_old_radius'].values()) for r in exhausted),
            'exceeds_by_radius': {name: sum(r['exceeds_old_radius'][name] for r in exhausted)
                                  for name in radii},
            'positive_ratio_points_for_exhausted': sum(len(r['positive_lambda_trials']) for r in exhausted),
            'minimum_nu_reduction_when_recovered': min((r['old_nu'] - r['best_nu_trial']['nu']
                                                       for r in exhausted if r['recovered_below_source_count']),
                                                      default=None),
            'maximum_nu_reduction_when_recovered': max((r['old_nu'] - r['best_nu_trial']['nu']
                                                       for r in exhausted if r['recovered_below_source_count']),
                                                      default=None),
            'best_residual_lower_over_exhausted': max(
                (r['best_lower_bound_trial']['lower_bound_numeric'] for r in exhausted), default=None),
            'old_new_zero_lambda_nu_disagreements': sum(r['old_nu'] != r['zero_lambda_nu'] for r in rows),
            'old_new_zero_lambda_max_absolute_bound_difference': max(
                (abs(r['old_lower_bound'] - r['zero_lambda_lower_bound']) for r in rows), default=None)}


def trial(row):
    return {key: row[key] for key in ('lambda', 'lambda_over_mu', 'nu', 'lower_bound_numeric')}


def match(case, prior, current, ratios):
    if prior['receiver_source_shape'] != current['receiver_source_shape'] or prior['radii'] != current['radii']:
        raise ValueError(f'{case}: mismatched source geometry or historical radii.')
    src_count = prior['receiver_source_shape'][1]
    radii = prior['radii']
    paired = []
    for cap in prior['grid']:
        u = cap['u']
        for old in cap['rows']:
            mu = old['mu']
            now = [r for r in current['caps'][str(u)]['scalar_rows'] if r['mu'] == mu]
            zero, = [r for r in now if r['lambda'] == 0.]
            pos = [trial(r) for r in now if r['lambda'] > 0.]
            if {r['lambda_over_mu'] for r in pos} != set(ratios) or len(pos) != len(ratios):
                raise ValueError(f'{case}: incomplete positive-lambda coverage at u={u} and mu={mu}.')
            if any(r['nu'] is None or r['lower_bound_numeric'] is None for r in [zero, *pos]):
                raise ArithmeticError(f'{case}: unresolved numerical point at u={u} and mu={mu}.')
            best_nu = min(pos, key=lambda r: (r['nu'], -r['lower_bound_numeric'], r['lambda']))
            best_bound = max(pos, key=lambda r: (r['lower_bound_numeric'], -r['nu'], -r['lambda']))
            bound = best_bound['lower_bound_numeric']
            paired.append({'u': u, 'mu': mu, 'old_nu': old['nu'], 'zero_lambda_nu': zero['nu'],
                           'old_lower_bound': old['lower_bound_numeric'],
                           'zero_lambda_lower_bound': zero['lower_bound_numeric'],
                           'old_nu_exhausted_sources': old['nu'] >= src_count,
                           'recovered_below_source_count': best_nu['nu'] < src_count,
                           'positive_lower_bound': bound > 0.,
                           'exceeds_old_radius': {name: bound > radius for name, radius in radii.items()},
                           'best_nu_trial': best_nu, 'best_lower_bound_trial': best_bound,
                           'positive_lambda_trials': pos})
    return {'source_count': src_count, 'old_radii': radii, **summarize(paired, radii),
            'by_cap': {str(cap['u']): summarize([r for r in paired if r['u'] == cap['u']], radii)
                       for cap in prior['grid']}, 'paired_points': paired}


def main():
    old, new = read_json(OLD), read_json(NEW)
    ratios = [t for t in new['ratios'] if t > 0.]
    cases = {case: match(case, record, new['cases'][case], ratios)
             for case, record in old['cases'].items()}
    aggregate = {key: sum(row[key] for row in cases.values()) for key in (
        'old_grid_points', 'old_exhausted_points', 'recovered_below_source_count', 'positive_lower_bound',
        'exceeds_any_old_radius', 'positive_ratio_points_for_exhausted', 'old_new_zero_lambda_nu_disagreements')}
    if aggregate['old_grid_points'] != 180 or aggregate['old_exhausted_points'] != 98:
        raise ArithmeticError('The historical 180-point / 98-exhausted reference changed.')
    if aggregate['old_new_zero_lambda_nu_disagreements']:
        raise ArithmeticError('The new zero-multiplier numerical baseline disagrees with historical inertia.')
    result = {'scope': 'finite binary64 primary scalar grid; same view, cap and mu; no interval claim',
              'unit': 'historical (view,u,mu) point, counted once if any of the sampled positive lambdas qualifies',
              'old_input': str(OLD.relative_to(RESEARCH)), 'new_input': str(NEW.relative_to(RESEARCH)),
              'positive_lambda_over_mu_ratios': ratios, 'historical_caps': old['caps'],
              'historical_mu': old['multipliers'], 'aggregate': aggregate, 'cases': cases}
    OUT.write_text(json.dumps(result, separators=(',', ':'), allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(aggregate, ensure_ascii=False))
    for case, row in cases.items():
        print(case, *(row[k] for k in ('old_exhausted_points', 'recovered_below_source_count',
                                     'positive_lower_bound', 'exceeds_any_old_radius',
                                     'best_residual_lower_over_exhausted')))


if __name__ == '__main__':
    main()
