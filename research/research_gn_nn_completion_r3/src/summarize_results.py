"""Summarize already-sealed offline scores; never train or select candidates."""
from pathlib import Path
import argparse
import csv
import json
import re
from collections import defaultdict
import numpy as np


def mean(values):
    values = [x for x in values if x is not None]
    return float(np.mean(values)) if values else None


def metric(row, candidate, field):
    obj = row.get(candidate)
    return obj.get(field) if obj else None


def summarize(scores, by_amplitude=False):
    grouped = defaultdict(list)
    names = ('run', 'level', 'distribution', 'noise') + (('amplitude',) if by_amplitude else ())
    for row in scores['rows']:
        key = (row['run'], row.get('level', 'low'), row.get('distribution', 'ID'), row['noise'])
        grouped[key + ((row['amplitude'],) if by_amplitude else ())].append(row)
    summaries = []
    for key, rows in grouped.items():
        result = dict(zip(names, key))
        selected = [r for r in rows if r.get('top1')]
        result.update(n_observations=len(rows), n_parents=len({r['parent_id'] for r in rows}),
                      cost_scope=rows[0].get('cost_scope', 'INCREMENTAL_WARM_MODEL_POST_GN'),
                      n_selected=len(selected), n_joint_success=sum(r.get('joint_success', False) for r in rows),
                      n_within_incremental_budget=sum(r.get('accounted_within_budget', False) for r in rows),
                      n_not_run=sum(r.get('status') == 'NOT_RUN' for r in rows),
                      physical_qualification_rate=len(selected)/len(rows),
                      joint_success_rate=sum(r.get('joint_success', False) for r in rows)/len(rows),
                      anchor_NRMSE=mean([metric(r, 'anchor', 'NRMSE') for r in rows]),
                      wall_seconds_mean=mean([r.get('wall_seconds') for r in rows]),
                      total_LU=sum(r.get('costcounts', {}).get('LU', 0) for r in rows),
                      total_complex_RHS=sum(r.get('costcounts', {}).get('actual_complex_RHS', 0) for r in rows))
        for name in ('top1', 'old_score_top1', 'qualified_bestK_OFFLINE', 'all_bestK_OFFLINE'):
            result[name+'_NRMSE'] = mean([metric(r, name, 'NRMSE') for r in rows])
            result[name+'_IoU'] = mean([metric(r, name, 'IoU_Otsu') for r in rows])
            result[name+'_physical_rate'] = sum(bool(r.get(name) and r[name].get('physical_data_feasible')) for r in rows)/len(rows)
            result[name+'_joint_rate'] = sum(bool(r.get(name) and r[name].get('joint_success')) for r in rows)/len(rows)
        raw = [next((c for c in r.get('candidates', []) if c['label'] == 'raw_1'), None) for r in rows]
        result['raw_full_NRMSE'] = mean([c['NRMSE'] for c in raw if c])
        result['raw_full_physical_rate'] = sum(bool(c and c['physical_data_feasible']) for c in raw)/len(rows)
        result['top1_paired_anchor_NRMSE'] = mean([r['anchor']['NRMSE'] for r in selected])
        result['top1_absolute_NRMSE_gain'] = mean([r['anchor']['NRMSE']-r['top1']['NRMSE'] for r in selected])
        base = result['top1_paired_anchor_NRMSE']
        result['top1_relative_NRMSE_gain'] = result['top1_absolute_NRMSE_gain']/base if base else None
        summaries.append(result)
    return summaries


def paired_comparisons(scores):
    """Mother-object clustered paired bootstrap; noise is not a new sample."""
    lookup = {(r['run'], r['obs_id']): r for r in scores['rows']}
    runs = list(dict.fromkeys(r['run'] for r in scores['rows']))
    out = []
    for run in runs:
        tag = re.fullmatch(r'(.*)_(v_image|v_physics|j_full|j_null_current|j_null)_seed(\d+)', run)
        peer_runs = []
        if tag:
            for other in runs:
                peer = re.fullmatch(r'(.*)_(v_image|v_physics|j_full|j_null_current|j_null)_seed(\d+)', other)
                if peer and peer.group(1) == tag.group(1) and peer.group(3) == tag.group(3):
                    peer_runs.append(other)
        for comparator in ['anchor', 'continue_trf', 'weakreg45_anchor'] + peer_runs:
            if run == comparator:
                continue
            for level in ('low', 'stress'):
                parents = defaultdict(list)
                joint_parents = defaultdict(list)
                planned = 0
                for row in scores['rows']:
                    if row['run'] != run or row.get('level', 'low') != level:
                        continue
                    planned += 1
                    top = row.get('top1')
                    other = row.get('anchor') if comparator == 'anchor' else lookup.get((comparator, row['obs_id']), {}).get('top1')
                    other_row = lookup.get((comparator, row['obs_id']))
                    if comparator == 'anchor':
                        anchor_candidate = next((c for c in row.get('candidates',[]) if c['kind']=='anchor'),None)
                        if anchor_candidate is not None:
                            joint_parents[row['parent_id']].append(int(row.get('joint_success',False))-int(anchor_candidate.get('joint_success',False)))
                    elif other_row is not None and row.get('status')!='NOT_RUN' and other_row.get('status')!='NOT_RUN':
                        joint_parents[row['parent_id']].append(int(row.get('joint_success',False))-int(other_row.get('joint_success',False)))
                    if top and other:
                        parents[row['parent_id']].append(other['NRMSE']-top['NRMSE'])
                values = np.array([np.mean(v) for v in parents.values()])
                if not len(values) and not joint_parents:
                    continue
                rng = np.random.default_rng(20261107)
                boot = np.mean(values[rng.integers(0, len(values), size=(2000, len(values)))], axis=1)if len(values)else None
                joint_values=np.array([np.mean(v)for v in joint_parents.values()])
                joint_boot=np.mean(joint_values[rng.integers(0,len(joint_values),size=(2000,len(joint_values)))],axis=1)if len(joint_values)else np.array([np.nan])
                out.append(dict(run=run, comparator=comparator, level=level,
                                planned_observations=planned,
                                paired_observations=sum(map(len, parents.values())),
                                independent_parents=len(values),
                                mean_parent_NRMSE_improvement=float(values.mean())if len(values)else None,
                                paired_parent_bootstrap95_low=float(np.quantile(boot, .025))if len(values)else None,
                                paired_parent_bootstrap95_high=float(np.quantile(boot, .975))if len(values)else None,
                                joint_success_difference=float(joint_values.mean())if len(joint_values)else None,
                                joint_bootstrap95_low=float(np.quantile(joint_boot,.025))if len(joint_values)else None,
                                joint_bootstrap95_high=float(np.quantile(joint_boot,.975))if len(joint_values)else None,
                                joint_independent_parents=len(joint_values),
                                population='available paired top1; selection failures reported separately, never silently imputed as recovery',
                                sign='positive means the named run has lower NRMSE'))
    return out


def write_csv(path, rows):
    if not rows:
        path.write_text('')
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scores', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    scores = json.loads(args.scores.read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    summaries = summarize(scores)
    amplitude_summary = summarize(scores, by_amplitude=True)
    comparisons = paired_comparisons(scores)
    flat = []
    for row in scores['rows']:
        item = {k: row.get(k) for k in ('obs_id', 'parent_id', 'run', 'split', 'noise', 'level', 'distribution', 'family', 'digit', 'amplitude', 'status', 'selection_status', 'accounting_status', 'cost_scope', 'wall_seconds', 'joint_success', 'accounted_within_budget')}
        for name in ('anchor', 'top1', 'old_score_top1', 'qualified_bestK_OFFLINE', 'all_bestK_OFFLINE'):
            for field in ('NRMSE', 'RMSE', 'IoU_Otsu', 'relative_residual'):
                item[name+'_'+field] = metric(row, name, field)
        item.update(LU=row.get('costcounts', {}).get('LU', 0),
                    complex_RHS=row.get('costcounts', {}).get('actual_complex_RHS', 0))
        for k, value in (row.get('current_metrics') or {}).items():
            if isinstance(value, (int, float)):
                item[k] = value
        flat.append(item)
    write_csv(args.output/'RESULTS_BY_CASE.csv', flat)
    write_csv(args.output/'COHORT_SUMMARY.csv', summaries)
    write_csv(args.output/'AMPLITUDE_SUMMARY.csv', amplitude_summary)
    write_csv(args.output/'PAIRED_PARENT_COMPARISONS.csv', comparisons)
    result = dict(source_scores=str(args.scores.resolve()),
                  state='MEASURED_OFFLINE_AFTER_SEAL', cohort_summary=summaries, amplitude_summary=amplitude_summary,
                  paired_parent_comparisons=comparisons,
                  oracle_rows=scores.get('oracle_rows', []),
                  interpretation='Top1 is frozen online selection. Best@K and current oracle use GT offline. Available-top1 means must be read with selection and joint success rates.')
    (args.output/'RESULTS_SUMMARY.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(output=str(args.output), cohorts=len(summaries), paired_comparisons=len(comparisons), rows=len(flat))))


if __name__ == '__main__':
    main()
