"""Read the archived generating materials only after the old-only regions are saved."""
from collections import Counter
import json

import numpy as np

from data14 import R10, R12, ROOT, load_cases
from pivot_controls import profile_witness
from pivot_experiment import METHODS, make_case


def truth_data():
    layout = {x['label']: x for x in json.loads((R12 / 'results/spectral_new_layouts.json').read_text())}
    cases = {}
    for i, row in enumerate(json.loads((R12 / 'results/spectral_online.json').read_text())):
        cases[f'spectral-{i:02d}'] = (np.asarray(layout[row['label']]['pair'][row['truth_side']]),
                                    len(row['design']['central_roots']))
    for file in ('natural_inputs.npz', 'confirm_inputs.npz'):
        with np.load(R10 / 'results' / file) as data:
            for name in data['case_ids']:
                cases[str(name)] = (data[f'{name}/truth'], None)
    return cases


def score_region(boxes, m):
    if not len(boxes):
        return dict(gt_retained=False, gt_in_hull=False, gt_distance_to_outer_union_linf=None)
    lo, hi = boxes[:, 0], boxes[:, 1]
    covered = np.all((lo <= m) & (m <= hi), axis=1)
    return dict(gt_retained=bool(np.any(covered)),
                gt_in_hull=bool(np.all(lo.min(axis=0) <= m) and np.all(m <= hi.max(axis=0))),
                gt_distance_to_outer_union_linf=float(np.min(np.maximum(np.maximum(lo - m, m - hi), 0.).max(axis=1))))


def main():
    summary = json.loads((ROOT / 'results/pivot_old.json').read_text())
    scored = {}
    with np.load(ROOT / 'results/pivot_old_regions.npz') as archive:
        known = truth_data()
        for public in load_cases():
            name = public.case_id
            m, central_roots = known[name]
            score = profile_witness(make_case(public), m)
            method = {a: score_region(archive[f'{name}/{a}'], m) for a in METHODS}
            scored[name] = dict(cohort=summary['cases'][name]['cohort'],
                                joint_contract_score=str(score), generated_feasible=bool(score < 1),
                                generated_infeasible=bool(score > 1), legacy_central_root_count=central_roots,
                                methods=method)
    counts = {}
    for group in ('spectral', 'natural', 'aperture', 'spectral_zero_central_roots'):
        if group == 'spectral_zero_central_roots':
            cases = [row for row in scored.values() if row['legacy_central_root_count'] == 0]
        else:
            cases = [row for row in scored.values() if row['cohort'] == group]
        counts[group] = dict(cases=len(cases), feasible=sum(x['generated_feasible'] for x in cases),
                             infeasible=sum(x['generated_infeasible'] for x in cases),
                             legacy_roots_distribution=dict(Counter(x['legacy_central_root_count'] for x in cases))
                             if group.startswith('spectral') else None,
                             coverage={a: dict(retained=sum(x['methods'][a]['gt_retained'] for x in cases),
                                 feasible_but_excluded=sum(x['generated_feasible'] and not x['methods'][a]['gt_retained']
                                                           for x in cases)) for a in METHODS})
    (ROOT / 'results/pivot_score.json').write_text(json.dumps(dict(cohorts=counts, cases=scored), indent=2,
        allow_nan=False)+'\n')
    print(json.dumps(counts, indent=2))


if __name__ == '__main__':
    main()
