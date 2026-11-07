"""Frozen VAL subset: unchanged R4/R3 public corrector, then separate GT score.

Default --plan reads metadata only. --execute requires all five final prediction
bundles. No training, GT-based candidate selection, or new physics backend.
"""
from __future__ import annotations

import os
for _key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_key] = '1'

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import csv
import json
import multiprocessing
import re
from pathlib import Path
import sys
import time
import numpy as np

STUDY = Path(__file__).resolve().parents[1]
REPO = STUDY.parents[1]
R4 = STUDY.parent/'research_gn_nn_preservation_r4'
sys.path.insert(0, str(REPO))
ARMS = ('pixel_closed_coupled', 'haar_closed_coupled', 'db2_closed_coupled',
        'haar_closed_no_fullfield', 'haar_free_coupled')
PLAN_NAME = 'CORRECTOR_GROUP_SELECTION.json'
EXTRA_ARM = 'haar_closed_noiseball_coupled'


def tagged(name, tag=''):
    """Stable filename suffix; defaults preserve all primary paths."""
    path = Path(name)
    return f'{path.stem}.{tag}{path.suffix}' if tag else name


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    from research.research_gn_search_space_r2.src.search import safe
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(safe(value), indent=2, allow_nan=False)+'\n')


def read_npz(path):
    with np.load(path, allow_pickle=False) as bundle:
        return {key: np.array(bundle[key]) for key in bundle.files}


def inventory(path):
    path = Path(path).resolve()
    stat = path.stat()
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def records_from(path):
    value = read_json(path)
    records = value if isinstance(value, list) else value['records']
    # Reject accidental test inclusion rather than silently scoring it.
    if any(str(r.get('split', '')).lower() not in ('val', 'validation') for r in records):
        raise ValueError('R5 VAL manifest must contain VAL records only')
    if len({r['obs_id'] for r in records}) != len(records):
        raise ValueError('Duplicate observation IDs')
    return records


def category(record):
    family = record.get('family')
    if family == 'digits' and str(record.get('digit')) in ('2', '8'):
        return 'digit'+str(record['digit'])
    if family in ('islands', 'primitive_composites'):
        return family
    return None


def prepare_plan(study, r4, *, result_root=None, audit_root=None, anchor_root=None, phase="strict", arms=ARMS, tag=""):
    """Predeclare solely from labels, never GT arrays or predicted quality."""
    result_root = Path(result_root) if result_root is not None else study/'results'
    audit_root = Path(audit_root) if audit_root is not None else study/'audit'
    anchor_root = Path(anchor_root) if anchor_root is not None else r4/'anchors'
    manifest = result_root/'VAL_MANIFEST.json'
    if manifest.exists():
        records = records_from(manifest)
        origin = 'R5 VAL metadata manifest'
    else:
        manifest = r4/'data/manifests/evaluation.json'
        records = [r for r in read_json(manifest)['records'] if r['split'] == 'val']
        origin = 'Explicit metadata-only R4 VAL fallback; R5 manifest not yet present'
    selected = [r for r in records if category(r)]
    groups = ('digit2', 'digit8', 'islands', 'primitive_composites')
    counts = {g: dict(observations=sum(category(r) == g for r in selected),
                     mothers=len({r['parent_id'] for r in selected if category(r) == g})) for g in groups}
    missing = [g for g in groups if not counts[g]['mothers']]
    if missing:
        raise ValueError('Requested families unavailable; no replacement allowed: '+str(missing))
    parent_labels = {}
    for r in records:
        parent_labels.setdefault(r['parent_id'], set()).add(category(r))
    if any(len(labels) != 1 for labels in parent_labels.values()):
        raise ValueError('Mother metadata has inconsistent category labels')
    chosen_parents = {r['parent_id'] for r in selected}
    if {r['obs_id'] for r in records if r['parent_id'] in chosen_parents} != {r['obs_id'] for r in selected}:
        raise ValueError('Selection did not retain all observations of each chosen mother')
    allowed = ('obs_id', 'parent_id', 'split', 'noise', 'public_path', 'family', 'digit', 'group')
    plan = dict(schema='r5-corrector-predeclared-val-v1', arms=list(arms),
                metadata_source=str(manifest.resolve()), metadata_origin=origin,
                selection_rule='All VAL mothers labelled digits 2/8, islands, primitive_composites; all observations per mother; no quality or GT selection',
                all_val_observations=len(records), counts_by_category=counts,
                observations=len(selected), mothers=len(chosen_parents),
                nominal_48_observations_24_mothers_match=len(selected) == 48 and len(chosen_parents) == 24,
                nominal_count_note='Actual matching metadata determines count; no truncation or silent family replacement',
                family_observation_inventory=dict(Counter(str(r.get('family')) for r in records)),
                records=[{k: r[k] for k in allowed if k in r} for r in selected],
                public_corrector='research.research_gn_nn_preservation_r4.src.evaluate.evaluate_arm -> unchanged R3 evaluator',
                evaluator_arm='v_physics', raw_scales=[.5, 1.],
                incremental_budget=dict(wall_seconds=10., max_LU=100, max_complex_RHS=6000),
                each_TRF=dict(max_nfev=8, seconds=3.),
                online_ranking='Existing evaluator: qualified nearest full unclipped proposal, else qualified anchor, else failure',
                anchor_accounting='Original necessary GN solve and feature receipt separately added to deployment; not free cached GN',
                truth_arrays_read=False)
    if phase == 'noiseaware':
        plan.update(phase='noiseaware', data_scope='Known DEVELOPMENT VAL only; original156observations/78mothers; no TEST',
                    anchor_root=str(anchor_root.resolve()), result_root=str(result_root.resolve()),
                    phase_II_change='Only common GN stopping changed to every-scale noise-aware discrepancy; unchanged architecture, lambda construction, rho=1, loss, depths and same12epoch maximum training budget; all five arms new weights from scratch',
                    method_tuning='No after-method hyperparameter tuning or depth selection from GT')
    if tag:
        plan.update(tag=tag, extension_scope='Only explicitly listed arms; primary five predictions/plans/seals/CSVs are not rerun or overwritten')
        if EXTRA_ARM in arms:
            plan['additional_mechanism_contract'] = 'Parent-declared Haar noise-ball data objective tolerance only; same lambda_field=100, all other network/loss/depth/training caps fixed, new weights from scratch. This additional arm is separate from primary5.'
    path = audit_root/tagged(PLAN_NAME, tag)
    if path.exists():
        frozen = read_json(path)
        for key in (('arms', 'records', 'incremental_budget', 'raw_scales', 'each_TRF', 'anchor_root', 'phase_II_change') if phase == 'noiseaware' else ('arms', 'records', 'incremental_budget', 'raw_scales', 'each_TRF')):
            if frozen[key] != plan[key]:
                raise ValueError('Existing predeclared subset/protocol differs: '+key)
        return frozen
    write_json(path, plan)
    return plan


def _allocated_seconds(bundle, row, total):
    if 'inference_seconds_amortized' not in bundle:
        raise ValueError('Missing actual amortized NN inference cost; refusing free preparation')
    values = np.asarray(bundle['inference_seconds_amortized'], float).reshape(-1)
    if len(values) not in (1, total):
        raise ValueError('Inference cost must be allocated per-case or scalar per-case allocation')
    seconds = float(values[0] if len(values) == 1 else values[row])
    if not np.isfinite(seconds) or seconds < 0:
        raise ValueError('Invalid inference preparation seconds')
    return seconds


def _online_worker(task):
    """Public-only boundary: no truth path, material/support truth or GT import."""
    from research.research_gn_nn_preservation_r4.src.evaluate import evaluate_arm, save_result
    study, r4 = Path(task['study']), Path(task['r4'])
    anchor_root = Path(task['anchor_root'])
    output_dir = Path(task['output_dir'])
    record = task['record']
    obs_id, arm = record['obs_id'], task['arm']
    directory = output_dir/arm
    if (directory/f'{obs_id}.json').exists() or (directory/f'{obs_id}.npz').exists():
        raise FileExistsError('Refusing to overwrite existing candidate evidence: '+str(directory/obs_id))
    started = time.perf_counter()
    receipt, arrays = None, None
    try:
        operators = read_npz(r4/'data/operators.npz')
        public_path = Path(record['public_path'])
        observation = read_npz(public_path if public_path.is_absolute() else r4/'data'/public_path)
        public = dict(operators, Y=observation['Y'])
        anchor = read_npz(anchor_root/f'{obs_id}.npz')
        if 'm0' not in anchor or 'j0' not in anchor:
            raise ValueError('Original GN anchor m0/j0 missing')
        prep = dict(wall_seconds=task['inference_seconds'], LU=0, actual_complex_RHS=0,
                    complete=True, scope='actual batched NN inference allocated per case',
                    inference_allocation='inference_seconds_amortized in prediction bundle; scalar means per-case allocation, not whole-batch wall')
        receipt, arrays = evaluate_arm(public, anchor, {'m_proposal': np.asarray(task['material'])},
                                      arm='v_physics', noise='clean' if float(record['noise']) == 0 else 'noisy',
                                      preparation_cost=prep)
        receipt.update(arm=arm, evaluator_arm='v_physics', obs_id=obs_id,
                       parent_id=record['parent_id'], candidate_selection_uses_GT=False,
                       source_attribution='Unchanged R4/R3 evaluate_arm; full real-material TRF, no R5 solver backend')
        # Necessary anchor construction remains an explicit additional invoice.
        base = read_json(anchor_root/f'{obs_id}.json')
        feature = base.get('feature_counts', {})
        counts = base['counts']
        anchor_wall = float(base['solve_wall_seconds'])+float(base.get('feature_wall_seconds', 0.))
        anchor_lu = int(counts.get('LU', 0))+int(feature.get('LU', 0))
        anchor_rhs = int(counts.get('actual_complex_RHS', 0))+int(feature.get('actual_complex_RHS', 0))
        receipt['necessary_anchor_deployment'] = dict(wall_seconds=anchor_wall, LU=anchor_lu,
            actual_complex_RHS=anchor_rhs, source=str(anchor_root/f'{obs_id}.json'),
            scope='original GN solve plus NN feature preparation; historical reused artifact is not free deployment')
        receipt['deployment_total'] = dict(wall_seconds=receipt['wall_seconds']+anchor_wall,
            LU=receipt['costcounts'].get('LU', 0)+anchor_lu,
            actual_complex_RHS=receipt['costcounts'].get('actual_complex_RHS', 0)+anchor_rhs)
    except Exception as exc:
        if receipt is None:
            receipt = dict(arm=arm, obs_id=obs_id, parent_id=record['parent_id'], status='FAILED',
                           stop_reason=f'{type(exc).__name__}: {exc}', top1=None, candidates=[], attempts=[],
                           trace=[], truth_access='NONE_PUBLIC_INPUT_ONLY', accounting_status='HOLD_FAILED_PREPARATION_OR_EVALUATION',
                           wall_seconds=time.perf_counter()-started, costcounts={})
            arrays = dict(candidate_materials=np.empty((0, 576)))
        else:
            # Preserve all already-paid solver evidence if deployment accounting fails.
            receipt.update(accounting_status='HOLD_DEPLOYMENT_ACCOUNTING',
                           accounting_error=f'{type(exc).__name__}: {exc}')
    receipt['worker_total_seconds_including_io_and_setup'] = time.perf_counter()-started
    save_result(directory, obs_id, receipt, arrays)
    return dict(arm=arm, obs_id=obs_id, status=receipt['status'])


def execute_online(study, r4, plan, workers, *, result_root=None, audit_root=None, anchor_root=None, tag=""):
    """Freeze all five public prediction inputs and complete online ranking first."""
    result_root = Path(result_root) if result_root is not None else study/'results'
    audit_root = Path(audit_root) if audit_root is not None else study/'audit'
    anchor_root = Path(anchor_root) if anchor_root is not None else r4/'anchors'
    output_dir = result_root/(f'CORRECTED_ONLINE.{tag}' if tag else 'CORRECTED_ONLINE')
    manifest = result_root/'VAL_MANIFEST.json'
    records = records_from(manifest)
    if {r['obs_id'] for r in records} != {r['obs_id'] for r in read_json(r4/'data/manifests/evaluation.json')['records'] if r['split'] == 'val'}:
        raise ValueError('R5 manifest must represent the unchanged entire original VAL cohort')
    tasks, input_inventory = [], [inventory(manifest)]
    for arm in plan['arms']:
        path = result_root/'VAL_PREDICTIONS'/f'{arm}.npz'
        bundle = read_npz(path)
        ids = [str(x.decode() if isinstance(x, bytes) else x) for x in bundle['obs_id']]
        if len(set(ids)) != len(ids) or set(ids) != {r['obs_id'] for r in records}:
            raise ValueError('Each arm requires the identical complete VAL cohort: '+arm)
        if bundle['m_final'].shape != (len(ids), 576) or not np.isfinite(bundle['m_final']).all():
            raise ValueError('Finite raw m_final shape(N,576) required: '+arm)
        by_id = {obs: i for i, obs in enumerate(ids)}
        input_inventory.append(inventory(path))
        for record in plan['records']:
            i = by_id[record['obs_id']]
            tasks.append(dict(study=str(study), r4=str(r4), arm=arm, record=record,
                              anchor_root=str(anchor_root.resolve()), output_dir=str(output_dir.resolve()),
                              material=bundle['m_final'][i], inference_seconds=_allocated_seconds(bundle, i, len(ids))))
    write_json(audit_root/tagged('CORRECTOR_INPUT_INVENTORY.json', tag), dict(inputs=input_inventory,
        predictions_frozen_before_GT=True, integrity_method='size/mtime; no array hashes'))
    if workers == 1:
        statuses = [_online_worker(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn')) as pool:
            statuses = list(pool.map(_online_worker, tasks, chunksize=1))
    outputs = []
    for task in tasks:
        prefix = output_dir/task['arm']/task['record']['obs_id']
        outputs.extend([inventory(prefix.with_suffix('.json')), inventory(prefix.with_suffix('.npz'))])
    seal = dict(status='ALL_PUBLIC_CANDIDATES_AND_ONLINE_RANKING_SAVED',
                created_unix=time.time(), selection_plan=inventory(audit_root/tagged(PLAN_NAME, tag)),
                inputs=input_inventory, outputs=outputs, statuses=statuses, truth_arrays_read=False)
    write_json(audit_root/tagged('CORRECTOR_ONLINE_SEALED.json', tag), seal)
    return seal


def csv_write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value) if isinstance(value, (dict, list)) else value for key, value in row.items()})


def score_val_only(study, r4, plan, seal, *, result_root=None, tag=""):
    """Offline DEV GT: import scorer and open truth only after global online seal."""
    result_root = Path(result_root) if result_root is not None else study/'results'
    if seal['status'] != 'ALL_PUBLIC_CANDIDATES_AND_ONLINE_RANKING_SAVED':
        raise ValueError('Online completion seal required before truth access')
    for item in [seal['selection_plan']]+seal['inputs']+seal['outputs']:
        if inventory(item['path']) != item:
            raise ValueError('Sealed input/output changed before scoring: '+item['path'])
    from research.research_gn_nn_preservation_r4.src.score import material_metrics
    records = {r['obs_id']: r for r in records_from(result_root/'VAL_MANIFEST.json')}
    rows = []
    truth_cache = {}
    for arm in plan['arms']:
        for public_record in plan['records']:
            obs_id = public_record['obs_id']
            record = records[obs_id]
            if str(record['split']).lower() not in ('val', 'validation'):
                raise ValueError('Offline scorer refuses non-VAL truth')
            prefix = result_root/(f'CORRECTED_ONLINE.{tag}' if tag else 'CORRECTED_ONLINE')/arm/obs_id
            receipt, arrays = read_json(prefix.with_suffix('.json')), read_npz(prefix.with_suffix('.npz'))
            truth_path = Path(record['truth_path'])
            if not truth_path.is_absolute():
                truth_path = r4/'data'/truth_path
            if truth_path not in truth_cache:
                with np.load(truth_path, allow_pickle=False) as bundle:
                    truth_cache[truth_path] = {key: np.array(bundle[key]) for key in ('m', 'amplitude', 'alpha')}
            truth = truth_cache[truth_path]
            scored = []
            for candidate, material in zip(receipt['candidates'], arrays['candidate_materials']):
                metrics = material_metrics(material, truth['m'], truth['amplitude'], truth['alpha'])
                scored.append(dict(candidate, **metrics))
            selected = receipt.get('top1')
            top = None if selected is None else next(c for c in scored if c['index'] == selected['index'])
            qualified = [c for c in scored if c['physical_data_feasible']]
            all_best = min(scored, key=lambda c: c['NRMSE']) if scored else None
            qual_best = min(qualified, key=lambda c: c['NRMSE']) if qualified else None
            cost = receipt.get('costcounts', {})
            row = dict(arm=arm, obs_id=obs_id, parent_id=record['parent_id'], noise=record['noise'],
                family=record.get('family'), digit=record.get('digit'), status=receipt['status'],
                stop_reason=receipt.get('stop_reason'), accounting_status=receipt['accounting_status'],
                NRMSE=None if top is None else top['NRMSE'], RMSE=None if top is None else top['RMSE'],
                IoU=None if top is None else top['IoU_Otsu'],
                data_residual=None if top is None else top['relative_residual'],
                qualified=bool(top is not None and top['physical_data_feasible']),
                joint=bool(top is not None and top['physical_data_feasible'] and top['NRMSE'] <= .15 and top['IoU_Otsu'] >= .7 and receipt['accounting_status'] == 'COMPLETE'),
                top1index=None if top is None else top['index'],
                bestKall=all_best, bestKqualified=qual_best,
                bestKall_NRMSE=None if all_best is None else all_best['NRMSE'],
                bestKqualified_NRMSE=None if qual_best is None else qual_best['NRMSE'],
                saved_candidate_count=len(scored), wall=receipt['wall_seconds'],
                LU=cost.get('LU'), RHS=cost.get('actual_complex_RHS'),
                preparation=receipt.get('preparation_cost'), necessary_anchor_deployment=receipt.get('necessary_anchor_deployment'),
                deployment_total=receipt.get('deployment_total'),
                worker_total_seconds_including_io_and_setup=receipt.get('worker_total_seconds_including_io_and_setup'))
            row['accounted_within_incremental_budget'] = bool(
                receipt['accounting_status'] == 'COMPLETE' and row['wall'] <= 10.
                and row['LU'] is not None and row['LU'] <= 100
                and row['RHS'] is not None and row['RHS'] <= 6000)
            row['budgeted_joint'] = row['joint'] and row['accounted_within_incremental_budget']
            row['wall_overshoot_seconds'] = receipt.get('wall_overshoot_seconds')
            deployment = receipt.get('deployment_total', {})
            row.update(deployment_wall=deployment.get('wall_seconds'), deployment_LU=deployment.get('LU'),
                       deployment_RHS=deployment.get('actual_complex_RHS'))
            rows.append(row)
            write_json(prefix.with_suffix('.offline.json'), dict(scope='VAL_KNOWN_DEV_OFFLINE_ONLY',
                online_ranking_unchanged=True, candidates=scored, summary=row,
                bestK_scope='Offline minimum NRMSE over actual saved finite candidates, not an online selector'))
    csv_write(result_root/tagged('CORRECTED_BY_CASE.csv', tag), rows)
    summary = []
    for arm in plan['arms']:
        arm_rows = [r for r in rows if r['arm'] == arm]
        valid = [r for r in arm_rows if r['NRMSE'] is not None]
        item = dict(arm=arm, cases=len(arm_rows), selected=len(valid),
                    qualified=sum(r['qualified'] for r in arm_rows), joint=sum(r['joint'] for r in arm_rows),
                    budgeted_joint=sum(r['budgeted_joint'] for r in arm_rows),
                    status_counts=dict(Counter(r['status'] for r in arm_rows)),
                    offline_bestK_warning='Finite offline oracle diagnostics; never used for online selection')
        for key in ('NRMSE', 'RMSE', 'IoU', 'data_residual', 'bestKall_NRMSE', 'bestKqualified_NRMSE', 'wall', 'LU', 'RHS', 'deployment_wall', 'deployment_LU', 'deployment_RHS'):
            values = [r[key] for r in arm_rows if r.get(key) is not None]
            item[key+'_mean'] = float(np.mean(values)) if values else None
        summary.append(item)
    csv_write(result_root/tagged('CORRECTED_SUMMARY.csv', tag), summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--plan', action='store_true', help='Metadata only; default')
    mode.add_argument('--execute', action='store_true', help='Run correction only after all predictions are ready, then offline VAL score')
    parser.add_argument('--workers', type=int, choices=(1, 2), default=2)
    parser.add_argument('--study', type=Path, default=STUDY)
    parser.add_argument('--r4', type=Path, default=R4, help='Original physics/data/truth root; unchanged for Phase II')
    parser.add_argument('--phase', choices=('strict', 'noiseaware'), default='strict')
    parser.add_argument('--anchor-root', type=Path, help='Explicit public GN anchor directory; passed to spawned workers')
    parser.add_argument('--arms', nargs='+', choices=(*ARMS, EXTRA_ARM), default=list(ARMS),
                        help='Explicit finite method list; nonprimary lists require a unique --tag')
    parser.add_argument('--tag', default='', help='Artifact suffix, e.g. noiseball; isolates plans/seals/candidates/CSVs')
    args = parser.parse_args()
    if len(set(args.arms)) != len(args.arms):
        parser.error('--arms cannot contain duplicates')
    if args.tag and re.fullmatch(r'[A-Za-z0-9_-]+', args.tag) is None:
        parser.error('--tag accepts letters, digits, underscore and hyphen only')
    if tuple(args.arms) != ARMS and not args.tag:
        parser.error('A nonprimary --arms list requires a unique --tag to preserve primary artifacts')
    study, r4 = args.study.resolve(), args.r4.resolve()
    result_root = study/'results' if args.phase == 'strict' else study/'results/noiseaware'
    audit_root = study/'audit' if args.phase == 'strict' else study/'audit/noiseaware'
    anchor_root = args.anchor_root.resolve() if args.anchor_root else (r4/'anchors' if args.phase == 'strict' else study/'inputs/noiseaware_anchors')
    plan = prepare_plan(study, r4, result_root=result_root, audit_root=audit_root,
                        anchor_root=anchor_root, phase=args.phase, arms=args.arms, tag=args.tag)
    print(json.dumps(dict(mode='EXECUTE' if args.execute else 'PLAN_ONLY', observations=plan['observations'],
                         mothers=plan['mothers'], groups=plan['counts_by_category'],
                         nominal_count_match=plan['nominal_48_observations_24_mothers_match'])), flush=True)
    if args.execute:
        seal = execute_online(study, r4, plan, args.workers, result_root=result_root,
                              audit_root=audit_root, anchor_root=anchor_root, tag=args.tag)
        print(json.dumps(score_val_only(study, r4, plan, seal, result_root=result_root, tag=args.tag)), flush=True)


if __name__ == '__main__':
    main()
