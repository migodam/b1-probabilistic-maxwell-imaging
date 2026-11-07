"""Apply the declared VAL-only rule to completed, sealed model evaluations."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import time

STUDY = Path(__file__).resolve().parents[1]


def summaries(path):
    data = json.loads(Path(path).read_text())
    if any(r.get('split') != 'val' for r in data['rows']):
        raise ValueError('Selection consumes VAL only')
    grouped = defaultdict(list)
    for row in data['rows']:
        grouped[row['run']].append(row)
    out = []
    for run, rows in grouped.items():
        sealed = [r for r in rows if r.get('transport_status') == 'SEALED']
        parents = defaultdict(list)
        qualified = 0
        for row in sealed:
            top = row.get('top1')
            qualified += int(top is not None)
            parents[row['parent_id']].append((top or row['anchor'])['NRMSE'])
        complete = len(rows) == 156 and len(sealed) == 156 and len(parents) == 78 and all(len(x) == 2 for x in parents.values())
        mean = sum(sum(x)/len(x) for x in parents.values())/len(parents) if parents else None
        out.append(dict(run=run, complete=complete, planned=len(rows), completed=len(sealed),
                        parents=len(parents), qualified=qualified, physical_rate=qualified/156,
                        failure_filled_mean_NRMSE=mean,
                        selection_only_failures_filled_with_common_GN=sum(r.get('top1') is None for r in sealed)))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores', type=Path, required=True)
    p.add_argument('--mode', choices=['alpha', 'constrained', 'loss'], required=True)
    p.add_argument('--runs', nargs='+', required=True)
    p.add_argument('--reference-run')
    p.add_argument('--reference-scores',type=Path,nargs='*',default=[])
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    rows = summaries(a.scores)
    existing={r['run']for r in rows}
    for path in a.reference_scores:
        for row in summaries(path):
            if row['run']not in existing:
                rows.append(row);existing.add(row['run'])
    lookup = {r['run']: r for r in rows}
    candidates = [lookup[x] for x in a.runs if x in lookup and lookup[x]['complete']]
    reference = lookup.get(a.reference_run)
    if a.mode in ('alpha', 'loss'):
        if reference is None or not reference['complete']:
            raise ValueError('Complete VAL reference required')
        candidates = [r for r in candidates if r['physical_rate'] >= reference['physical_rate']-.05-1e-12]
    if a.mode == 'loss':
        candidates = [r for r in candidates if r['failure_filled_mean_NRMSE'] <= .95*reference['failure_filled_mean_NRMSE']]
    def key(r):
        alpha = 1. if r['run'] == a.reference_run else float(r['run'].rsplit('_posthoc_a', 1)[1]) if a.mode == 'alpha' else 0.
        return (round(r['failure_filled_mean_NRMSE'], 12), -alpha if a.mode == 'alpha' else -r['physical_rate'], r['run'])
    selected = min(candidates, key=key) if candidates else None
    result = dict(status='SELECTED' if selected else 'NO_ELIGIBLE_COMPLETED_ARM', mode=a.mode,
                  declared_contract='audit/VAL_SELECTION_CONTRACT.json', time_unix=time.time(),
                  reference_score_files=[str(path)for path in a.reference_scores],
                  truth_scope='VAL_ONLY', runs=a.runs, reference=reference,
                  summaries=[r for r in rows if r['run'] in a.runs or r['run'] == a.reference_run], selected=selected)
    if a.mode == 'alpha' and selected:
        result['selected_alpha'] = 1. if selected['run'] == a.reference_run else float(selected['run'].rsplit('_posthoc_a', 1)[1])
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
