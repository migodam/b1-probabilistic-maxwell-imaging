"""Parent-authorized offline R14 scoring/plots; never optimizes or runs physics.

--make-seal inventories public prediction files only and always leaves
 authorized_scoring=false. Only the parent may approve the resulting seal.
"""
from __future__ import annotations
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
import csv
import html
import importlib.util
import json
import math
from pathlib import Path
import sys
import time

CAMPAIGNS = ('dualmetric_pilot12_v1', 'dualmetric_remaining36_v1')
PAIRS = ('00', '01', '10', '11')
VIEWS = ('initial', 'raw', 'corrected', 'primary', 'secondary')
PREFIXES = (15, 30, 60, 90, 120)
OLD_METHODS = ('cie', 'dual_cie', 'g19', 'dr', 'a_plus_mm_r13', 'a_plus_arc_r7')
CASES = tuple(f'r13_fresh_{i:04d}' for i in range(48))
PILOT = (2, 5, 7, 16, 17, 24, 29, 37, 38, 42, 44, 45)
METRICS = ('NRMSE', 'RMSE', 'relativeL2', 'IoU', 'Phi', 'D', 'R',
           'residual_over_delta', 'CPU_seconds')


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write_json(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False, ensure_ascii=False)
        stream.write('\n')


def table(out, name, rows):
    write_json(out / (name + '.json'), rows)
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with (out / (name + '.csv')).open('x', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v
                             for k, v in row.items()})


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def as_number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def as_bool(value):
    if isinstance(value, bool):
        return value
    if str(value).lower() in ('true', '1'):
        return True
    if str(value).lower() in ('false', '0'):
        return False
    return None


def stamp(path):
    path = Path(path).resolve()
    return dict(path=str(path), exists=path.is_file(),
                bytes=path.stat().st_size if path.is_file() else None)


def checked_path(root, path):
    result = Path(path)
    result = (result if result.is_absolute() else root / result).resolve()
    if not result.is_relative_to(root):
        raise ValueError('Prediction/manifest path outside requested R14 root: ' + str(result))
    if any(p.lower() in ('sealed_truth', 'truth', 'private', 'data_private') for p in result.parts):
        raise ValueError('Private path is forbidden in a prediction seal')
    return result


def validate_manifests(root, entries):
    if {entry['campaign'] for entry in entries} != set(CAMPAIGNS) or len(entries) != 2:
        raise ValueError('Exactly the pilot12 and remaining36 campaigns are required')
    tasks, protocols = [], []
    for entry in entries:
        path = checked_path(root, entry['manifest_path'])
        cfg = read(path)
        if cfg != entry['manifest_snapshot'] or cfg['campaign'] != entry['campaign']:
            raise ValueError('Frozen campaign manifest changed')
        protocols.append(cfg['protocol'])
        wanted = {CASES[i] for i in PILOT}
        if entry['campaign'] == CAMPAIGNS[1]:
            wanted = set(CASES) - wanted
        actual = {(t['case_id'], str(t['pair'])) for t in cfg['tasks']}
        expected = {(cid, pair) for cid in wanted for pair in PAIRS}
        if actual != expected or len(cfg['tasks']) != len(expected):
            raise ValueError('Frozen case/pair partition differs from the registered 12+36 set')
        for task in cfg['tasks']:
            if task['campaign'] != cfg['campaign'] or task['task_id'] != task['case_id'] + '__pair' + str(task['pair']):
                raise ValueError('Unexpected task identity')
            checked_path(root, task['input_file'])
        tasks.extend(deepcopy(cfg['tasks']))
    if protocols[0] != protocols[1]:
        raise ValueError('Pilot and remainder use different numerical protocols')
    if len(tasks) != 192 or len({t['task_id'] for t in tasks}) != 192:
        raise ValueError('Exactly 192 unique predictions are required; no implicit partial-score exception')
    protocol = protocols[0]
    if float(protocol['CPU_cap']) != 120 or float(protocol['proposal_fraction']) != .6:
        raise ValueError('Unexpected frozen CPU budget')
    if tuple(protocol.get('material_box', ())) != (.02, 20.):
        raise ValueError('Unexpected physical material box')
    return tasks, protocol


def make_seal(root, seal_path, source_note=''):
    root = root.resolve()
    entries = []
    for campaign in CAMPAIGNS:
        manifest = root / 'configs' / (campaign + '.json')
        entries.append(dict(campaign=campaign, manifest_path=str(manifest),
                            manifest_snapshot=read(manifest), manifest_file=stamp(manifest)))
    tasks, protocol = validate_manifests(root, entries)
    by_campaign = {entry['campaign']: entry for entry in entries}
    for entry in entries:
        entry['tasks'] = []
        queue = root / 'queues' / entry['campaign'] / 'RECEIPT.json'
        queue_value = read(queue) if queue.is_file() else None
        entry['queue_receipt'] = stamp(queue)
        entry['queue_status'] = queue_value.get('status') if queue_value else None
    for task in tasks:
        folder = root / 'RESULTS' / task['campaign'] / task['task_id']
        receipt = folder / 'RECEIPT.json'
        receipt_value, receipt_error = None, None
        if receipt.is_file():
            try:
                receipt_value = read(receipt)
            except Exception as exc:
                receipt_error = repr(exc)
        item = dict(task_id=task['task_id'], case_id=task['case_id'], pair=str(task['pair']),
                    receipt_path=str(receipt), receipt_file=stamp(receipt),
                    status=receipt_value.get('status') if receipt_value else None,
                    receipt_error=receipt_error,
                    predictions={view: dict(**stamp(folder / (view + '.npz')),
                        metadata=stamp(folder / (view + '.json'))) for view in VIEWS},
                    history_path=str(folder / 'HISTORY.json'), history_file=stamp(folder / 'HISTORY.json'),
                    candidate_pool_path=str(folder / 'CANDIDATE_POOL.npz'),
                    candidate_pool_file=stamp(folder / 'CANDIDATE_POOL.npz'),
                    native_history_file=stamp(folder / 'NATIVE_HISTORY.json'))
        by_campaign[task['campaign']]['tasks'].append(item)
    value = dict(schema_version='r14_prediction_seal_v1', root=str(root),
        created_at_utc=datetime.now(timezone.utc).isoformat(), authorized_scoring=False,
        expected_tasks=192, campaigns=entries, source_note=source_note,
        authorization_note='Parent alone sets authorized_scoring=true after reviewing all192 task records; failures/missing predictions remain in denominator.',
        truth_access='NONE', provenance_scope='Path/existence/byte sizes only; no SHA or mtime checks',
        fixed_predictions=list(VIEWS), scope='EXISTING_EXPOSED48_NOT_FRESH_BLIND')
    if seal_path.resolve().is_relative_to(root / 'RESULTS') or seal_path.resolve().is_relative_to(root / 'queues'):
        raise ValueError('Seal belongs outside prediction/queue directories')
    seal_path.parent.mkdir(parents=True, exist_ok=True)
    write_json(seal_path, value)
    return value


def validate_stamp(root, item):
    path = checked_path(root, item['path'])
    exists = path.is_file()
    if exists != item['exists']:
        raise RuntimeError('Sealed file existence changed: ' + str(path))
    if exists and path.stat().st_size != item['bytes']:
        raise RuntimeError('Sealed file byte count changed: ' + str(path))
    return path


def authorize(root, seal_path):
    """Stdlib-only: all authorization and public identity checks precede GT imports."""
    seal = read(seal_path)
    if seal.get('schema_version') != 'r14_prediction_seal_v1':
        raise ValueError('Unsupported seal schema')
    if seal.get('authorized_scoring') is not True:
        raise PermissionError('Parent has not approved this prediction seal: authorized_scoring=true required')
    if Path(seal['root']).resolve() != root or seal.get('expected_tasks') != 192:
        raise ValueError('Seal root/task count mismatch')
    if set(seal.get('fixed_predictions', [])) != set(VIEWS):
        raise ValueError('All five prediction views must be explicitly sealed, including missing views')
    tasks, protocol = validate_manifests(root, seal['campaigns'])
    account = {}
    for entry in seal['campaigns']:
        if validate_stamp(root, entry['manifest_file']) != checked_path(root,entry['manifest_path']):
            raise ValueError('Manifest snapshot file path mismatch')
        queue_path = validate_stamp(root, entry['queue_receipt'])
        if queue_path != root/'queues'/entry['campaign']/'RECEIPT.json':
            raise ValueError('Queue receipt path mismatch')
        if queue_path.is_file():
            q = read(queue_path)
            if q.get('status') != entry['queue_status'] or q.get('status') == 'RUNNING' or q.get('running') or q.get('pending_tasks', q.get('pending', 0)):
                raise ValueError('A campaign queue is not terminal or changed after sealing')
        for row in entry['tasks']:
            ident = row['task_id']
            if ident in account:
                raise ValueError('Duplicate sealed task')
            folder = root / 'RESULTS' / entry['campaign'] / ident
            if Path(row['receipt_path']).resolve() != folder / 'RECEIPT.json':
                raise ValueError('Receipt identity path mismatch')
            receipt_path = validate_stamp(root, row['receipt_file'])
            if receipt_path != folder/'RECEIPT.json':
                raise ValueError('Receipt snapshot path mismatch')
            if receipt_path.is_file() and row.get('receipt_error') is None:
                receipt = read(receipt_path)
                if receipt.get('status') != row['status'] or receipt.get('task', {}).get('task_id') != ident:
                    raise ValueError('Receipt task/status differs from seal')
                if receipt.get('truth_access') != 'NONE':
                    raise ValueError('Unexpected online truth declaration')
            if set(row['predictions']) != set(VIEWS):
                raise ValueError('Incomplete fixed-prediction accounting')
            for view, item in row['predictions'].items():
                if validate_stamp(root, item) != folder / (view + '.npz'):
                    raise ValueError('Fixed prediction path mismatch')
                if validate_stamp(root, item['metadata']) != folder / (view + '.json'):
                    raise ValueError('Fixed prediction metadata mismatch')
            for field, filename in [('history_file', 'HISTORY.json'), ('candidate_pool_file', 'CANDIDATE_POOL.npz'), ('native_history_file', 'NATIVE_HISTORY.json')]:
                if validate_stamp(root, row[field]) != folder / filename:
                    raise ValueError('Task context path mismatch')
            if Path(row['history_path']).resolve() != folder / 'HISTORY.json' or Path(row['candidate_pool_path']).resolve() != folder / 'CANDIDATE_POOL.npz':
                raise ValueError('Task context alias mismatch')
            account[ident] = row
    if set(account) != {t['task_id'] for t in tasks}:
        raise ValueError('Seal must account for every expected task, including failures')
    for task in tasks:
        row = account[task['task_id']]
        if row['case_id'] != task['case_id'] or row['pair'] != str(task['pair']):
            raise ValueError('Sealed pair/case identity differs from manifest')
    return seal, tasks, protocol, account


def eligible(row, cutoff=120.):
    return (row.get('eligible') is True and row.get('legal') is True
            and finite(row.get('Phi')) and finite(row.get('CPU_seconds'))
            and row['CPU_seconds'] <= cutoff and isinstance(row.get('candidate'), int))


def pick(history, cutoff):
    rows = [row for row in history if eligible(row, cutoff)]
    return min(rows, key=lambda r: (r['Phi'], r['candidate'])) if rows else None


def raw_initial_choice(item):
    """GT-free incumbent repair for reporting; archived arrays are never changed."""
    raw=item['predictions']['raw'];initial=item['predictions']['initial']
    raw_meta=read(raw['metadata']['path']) if raw['metadata']['exists'] else None
    initial_meta=read(initial['metadata']['path']) if initial['metadata']['exists'] else None
    initial_ok=(initial['exists'] and initial_meta is not None and eligible(initial_meta))
    raw_ok=(raw['exists'] and raw_meta is not None and eligible(raw_meta))
    use_initial=initial_ok and (not raw_ok or initial_meta['Phi']<raw_meta['Phi'])
    return ('initial' if use_initial else 'raw'),dict(
        raw_original_source=raw['path'],raw_original_candidate=raw_meta.get('candidate') if raw_meta else None,
        raw_original_Phi=raw_meta.get('Phi') if raw_meta else None,
        initial_Phi=initial_meta.get('Phi') if initial_meta else None,
        raw_initial_fallback=bool(use_initial),raw_audit_comparable=initial_ok and raw_ok)


def strict48_bracket(history,receipt):
    termination=(receipt or {}).get('termination_reason',(receipt or {}).get('termination',{}))
    stop=termination.get('proposal') if isinstance(termination,dict) else None
    proposal=[r['CPU_seconds'] for r in history if r.get('phase')=='proposal' and finite(r.get('CPU_seconds'))]
    corrector=[r['CPU_seconds'] for r in history if r.get('phase')=='corrector' and finite(r.get('CPU_seconds'))]
    lower=upper=None;reason='UNAVAILABLE_PROPOSAL_STOP_OR_START_BRACKET'
    direct=((receipt or {}).get('task',{}).get('algorithm',(receipt or {}).get('task',{}).get('method'))=='direct'
            or any(r.get('phase')=='direct' for r in history))
    if direct:return dict(strict48_lower_deadline=None,strict48_upper_deadline=None,strict48_deadline_uncertainty_seconds=None,
        strict48_reason='NOT_APPLICABLE_DIRECT_WHOLE120_PHYSICAL_OPTIMIZATION',proposal_termination=stop,
        strict48_scope='Direct has no proposal/corrector hierarchy;ordinaryCPU/physicalcountprefixes remain valid')
    if isinstance(stop,str) and 'CPU' in stop:
        lower=upper=120.;reason='PROPOSAL_CPU_BUDGET_REACHED_72_PLUS48_CAPPED120'
    elif stop=='SCHEDULE_COMPLETE':
        lower=min(120.,max(proposal)+48.) if proposal else None
        upper=min(120.,min(corrector)+48.) if corrector else None
        reason='EARLY_SCHEDULE_END_START_BETWEEN_LAST_PROPOSAL_AND_FIRST_CORRECTOR_AUDIT'
    return dict(strict48_lower_deadline=lower,strict48_upper_deadline=upper,
        strict48_deadline_uncertainty_seconds=upper-lower if lower is not None and upper is not None else None,
        strict48_reason=reason,proposal_termination=stop,
        strict48_scope='Paid-pool cutoff robustness only; not an actual stopped rerun; all spent CPU charged')


def cumulative_physics(cost):
    """Original physical solve counts only; missing counters are unknown,never0."""
    if not isinstance(cost,dict):return dict(LU=None,complex_RHS=None,LU_source=None,RHS_source=None)
    lu_key='physical_LU' if 'physical_LU' in cost else 'LU' if 'LU' in cost else None
    lu=cost.get(lu_key) if lu_key else None
    if 'actual_complex_RHS' in cost:
        rhs=cost['actual_complex_RHS'];rhs_source='actual_complex_RHS aggregate'
    elif 'source_RHS' in cost and 'adjoint_RHS' in cost:
        rhs=cost['source_RHS']+cost['adjoint_RHS'];rhs_source='source_RHS+adjoint_RHS'
    else:rhs=None;rhs_source=None
    return dict(LU=lu if finite(lu) else None,complex_RHS=rhs if finite(rhs) else None,
                LU_source=lu_key,RHS_source=rhs_source)


def equal_physics_plans(caches,cases,expected_methods):
    plans=[]
    for cid in cases:
        group=[c for c in caches.values() if c['task']['case_id']==cid]
        finals=[cumulative_physics((c['receipt'] or {}).get('physical_cost')) for c in group]
        known=(len(group)==expected_methods and all(v['LU'] is not None and v['complex_RHS'] is not None for v in finals))
        lu=min(v['LU'] for v in finals) if known else None
        rhs=min(v['complex_RHS'] for v in finals) if known else None
        for c,final in zip(group,finals):
            usable=[]
            for row in c['history']:
                counts=cumulative_physics(row.get('physics_cost'))
                if known and eligible(row) and counts['LU'] is not None and counts['complex_RHS'] is not None and counts['LU']<=lu and counts['complex_RHS']<=rhs:usable.append(row)
            selected=min(usable,key=lambda r:(r['Phi'],r['candidate'])) if usable else None
            t=c['task'];plans.append(dict(task_id=t['task_id'],case_id=cid,
                pair=str(t.get('pair',t['method'])),solver=t.get('solver_label',t['method']),campaign=t['campaign'],
                view='common_physics_pool',selection_kind='pool',metadata=selected,
                candidate=selected['candidate'] if selected else None,
                source_file=c['item']['candidate_pool_path'],source_exists=c['item']['candidate_pool_file']['exists'],
                common_physical_LU=lu,common_physical_complex_RHS=rhs,
                actual_final_physics_counts=final,physics_counter_limits_available=known,
                normalization_scope='Componentwise minima of actual final paid original-Maxwell LU/RHS across all same-case methods; candidate completed and eligible before both limits. AuxFFT,4cell actions,SVD and CPU separately charged,not normalized/free; not a rerun'))
    return plans


def public_plans(root,tasks,account,issues,cases=None,expected_methods=4):
    """Public selection before truth; reusable for independently sealed cohorts."""
    cases=CASES if cases is None else tuple(cases)
    plans,caches=[],{}
    for task in tasks:
        ident=task['task_id'];item=account[ident];history=[];receipt=None
        try:
            if item['history_file']['exists']:
                history=read(item['history_path'])
                if not isinstance(history,list) or any(r.get('candidate')!=i for i,r in enumerate(history)):
                    raise ValueError('History/pool index schema is not consecutive')
        except Exception as exc:
            issues.append(dict(task_id=ident,issue='INVALID_PUBLIC_HISTORY',error=repr(exc)));history=[]
        try:
            if item['receipt_file']['exists'] and item.get('receipt_error') is None:receipt=read(item['receipt_path'])
        except Exception as exc:issues.append(dict(task_id=ident,issue='INVALID_PUBLIC_RECEIPT',error=repr(exc)))
        solver=task.get('solver_label',task.get('method','pair'+str(task.get('pair',''))))
        pair=str(task.get('pair',solver));base=dict(task_id=ident,case_id=task['case_id'],pair=pair,solver=solver,campaign=task['campaign'])
        caches[ident]=dict(task=task,history=history,receipt=receipt,item=item)
        try:raw_source,raw_audit=raw_initial_choice(item)
        except Exception as exc:
            raw_source='raw';raw_audit=dict(raw_initial_fallback=False,raw_audit_comparable=False,raw_audit_error=repr(exc))
            issues.append(dict(task_id=ident,issue='RAW_INITIAL_METADATA_ERROR',error=repr(exc)))
        for view in VIEWS:
            source_view,fallback=view,None
            if view=='raw':source_view=raw_source
            if view=='secondary' and not item['predictions'][view]['exists']:
                source_view,fallback='primary','NO_SAVED_CN99_SECONDARY_USE_PRIMARY'
            source=item['predictions'][source_view];meta=None
            try:
                if source['metadata']['exists']:meta=read(source['metadata']['path'])
            except Exception as exc:issues.append(dict(task_id=ident,view=view,issue='INVALID_VIEW_METADATA',error=repr(exc)))
            valid=[h for h in history if eligible(h)]
            if view=='raw':expected=pick([h for h in valid if h.get('phase') in ('initial','proposal')],120.)
            elif source_view=='primary':expected=pick(history,120.)
            elif source_view=='corrected':expected=pick(history,120.) if task.get('method')=='direct' else pick([h for h in valid if h.get('phase')=='corrector'],120.)
            elif source_view=='secondary':
                feasible=[h for h in valid if h.get('noise_feasible') is True and finite(h.get('R'))]
                expected=min(feasible,key=lambda h:(h['R'],h['Phi'],h['candidate'])) if feasible else None
            else:expected=next((h for h in history if h.get('phase')=='initial'),None)
            expected_index=expected.get('candidate') if expected else None
            consistency=meta is not None and meta.get('candidate')==expected_index
            if source['exists'] and not consistency:
                issues.append(dict(task_id=ident,view=view,issue='ARCHIVED_SELECTOR_METADATA_MISMATCH',
                    archived_candidate=meta.get('candidate') if meta else None,expected_public_candidate=expected_index))
            plans.append(dict(base,view=view,selection_kind='saved_view',expected_public_candidate=expected_index,
                archived_selector_consistent=consistency,source_view=source_view,secondary_fallback=fallback,
                source_file=source['path'],source_exists=source['exists'],metadata=meta,
                candidate=meta.get('candidate') if meta else None,**(raw_audit if view=='raw' else {})))
        for cutoff in PREFIXES:
            selected=pick(history,cutoff)
            plans.append(dict(base,view='paid_prefix',selection_kind='pool',CPU_prefix=cutoff,
                cutoff_CPU_seconds=float(cutoff),metadata=selected,candidate=selected['candidate'] if selected else None,
                source_file=item['candidate_pool_path'],source_exists=item['candidate_pool_file']['exists']))
        bracket=strict48_bracket(history,receipt)
        for label in ('lower','upper'):
            cutoff=bracket['strict48_'+label+'_deadline'];selected=pick(history,cutoff) if cutoff is not None else None
            plans.append(dict(base,view='strict48_'+label,selection_kind='pool',cutoff_CPU_seconds=cutoff,
                metadata=selected,candidate=selected['candidate'] if selected else None,
                source_file=item['candidate_pool_path'],source_exists=item['candidate_pool_file']['exists'],**bracket))
    for cid in cases:
        group=[c for c in caches.values() if c['task']['case_id']==cid]
        ends=[(c['receipt'] or {}).get('CPU_seconds',(c['receipt'] or {}).get('whole_process_CPU_seconds')) for c in group]
        common=min([120.]+ends) if len(group)==expected_methods and all(finite(t) for t in ends) else None
        for cutoff in PREFIXES:
            limit=min(float(cutoff),common) if common is not None else None
            for c in group:
                task=c['task'];selected=pick(c['history'],limit) if limit is not None else None
                plans.append(dict(task_id=task['task_id'],case_id=cid,pair=str(task.get('pair',task['method'])),
                    solver=task.get('solver_label',task['method']),campaign=task['campaign'],view='common_CPU_prefix',selection_kind='pool',
                    CPU_prefix=cutoff,cutoff_CPU_seconds=limit,common_recorded_process_CPU=common,
                    metadata=selected,candidate=selected['candidate'] if selected else None,
                    source_file=c['item']['candidate_pool_path'],source_exists=c['item']['candidate_pool_file']['exists']))
    plans.extend(equal_physics_plans(caches,cases,expected_methods))
    return plans,caches


def historical_rows(history_root, issues):
    path = history_root / 'supplement_aplus48/analysis/stack_comparison/SCOREBOARD_COMPACT.csv'
    with path.open(newline='', encoding='utf-8') as stream:
        original = list(csv.DictReader(stream))
    lookup = {}
    for row in original:
        if row.get('initialization') != 'uniform' or row.get('view') != 'primary' or row.get('solver') not in OLD_METHODS:
            continue
        key = row['case_id'], row['solver']
        if key in lookup:
            raise ValueError('Historical uniform primary rows are not unique: ' + str(key))
        lookup[key] = row
    output=[]
    for cid in CASES:
        for solver in OLD_METHODS:
            original = lookup.get((cid, solver)); row = dict(original or {})
            row.update(case_id=cid, solver=solver, initialization='uniform', view='primary',
                       provenance='REUSED_HISTORICAL_R13_UNIFORM48_PRIMARY', source_table=str(path),
                       budget_contract='Historical C2 strict<=48corrector after<=72proposal;not identical to new quartet absolute120',causal_comparison=False)
            for metric in METRICS + ('external_child_CPU_seconds','physical_LU_at_candidate','physical_complex_RHS_at_candidate'):
                row[metric] = as_number(row.get(metric))
            row['available'] = as_bool(row.get('available')) is True
            for key in ('noise_feasible', 'noise_certificate_eligible'):
                row[key] = as_bool(row.get(key))
            if original is None:
                row['availability_reason']='MISSING_HISTORICAL_ROW'
                issues.append(dict(case_id=cid, solver=solver, issue='MISSING_HISTORICAL_ROW'))
            output.append(row)
    return output, path


def accuracy_available(row):
    return row.get('available') is True and finite(row.get('NRMSE')) and finite(row.get('IoU'))


def recovered(row):
    return (row.get('available') is True and finite(row.get('NRMSE')) and finite(row.get('IoU'))
            and row['NRMSE'] <= .15 and row['IoU'] >= .70)


def summaries(rows):
    groups={}
    for row in rows:
        key=(row.get('solver'),row.get('view'),row.get('CPU_prefix'))
        groups.setdefault(key,[]).append(row)
    result=[]
    for (solver,view,prefix), group in groups.items():
        item=dict(solver=solver,view=view,CPU_prefix=prefix,expected_cases=48,
            row_count=len(group),available_cases=sum(r.get('available') is True for r in group),
            accuracy_available_cases=sum(accuracy_available(r) for r in group),
            failed_or_missing_run_cases=[r['case_id'] for r in group if r.get('run_status') in ('FAILED','MISSING',None,'')],
            recovered_posthoc_count=sum(recovered(r) for r in group),
            recovery_denominator=48,recovery_rule='POSTHOC_DESCRIPTIVE_NRMSE<=.15_AND_IoU>=.70_NOT_SELECTION',
            missing_failed_cases=[r['case_id'] for r in group if not r.get('available')],
            noise_feasible_count=sum(r.get('noise_feasible') is True for r in group),
            noise_certificate_eligible_count=sum(r.get('noise_certificate_eligible') is True for r in group))
        for metric in METRICS:
            values=sorted(r[metric] for r in group if finite(r.get(metric)))
            item[metric+'_n']=len(values)
            item[metric+'_mean']=sum(values)/len(values) if values else None
            item[metric+'_median']=(values[(len(values)-1)//2]+values[len(values)//2])/2 if values else None
        result.append(item)
    return result


def score_plans(plans, caches, truths, material_metrics, physical_counts, np, issues):
    pools={}; pictures={}; output=[]
    for ident,c in caches.items():
        if c['item']['candidate_pool_file']['exists']:
            try:
                with np.load(c['item']['candidate_pool_path'],allow_pickle=False) as z:
                    pools[ident]=z['m'].copy()
                if len(pools[ident]) != len(c['history']):
                    raise ValueError('Material pool count differs from paid history')
            except Exception as exc:
                pools.pop(ident,None);issues.append(dict(task_id=ident,issue='INVALID_MATERIAL_POOL',error=repr(exc)))
    for plan in plans:
        c=caches[plan['task_id']]; task=c['task']; meta=plan.get('metadata') or {}; m=None; reason=None
        try:
            if not plan['source_exists']:
                reason='MISSING_FIXED_PREDICTION_OR_POOL'
            elif plan['selection_kind']=='saved_view':
                with np.load(plan['source_file'],allow_pickle=False) as z:m=z['m'].copy()
                candidate=plan['candidate']; pool=pools.get(plan['task_id'])
                if isinstance(candidate,int) and pool is not None and 0<=candidate<len(pool):
                    if not np.array_equal(m,pool[candidate]):
                        issues.append(dict(task_id=plan['task_id'],view=plan['view'],issue='SAVED_VIEW_M_DIFFERS_FROM_POOL'))
                else:
                    issues.append(dict(task_id=plan['task_id'],view=plan['view'],issue='SAVED_VIEW_CANDIDATE_CONTEXT_UNAVAILABLE'))
            elif plan['candidate'] is not None and plan['task_id'] in pools:
                m=pools[plan['task_id']][plan['candidate']].copy()
            else:reason='NO_ELIGIBLE_COMPLETED_PAID_CANDIDATE'
        except Exception as exc:
            m=None;reason='PREDICTION_READ_ERROR';issues.append(dict(task_id=plan['task_id'],view=plan['view'],issue=reason,error=repr(exc)))
        metrics=material_metrics(m,truths.get(plan['case_id']),.02,20.)
        row=dict(case_id=plan['case_id'],solver=plan.get('solver','pair'+str(plan['pair'])),pair=plan['pair'],view=plan['view'],
            campaign=plan['campaign'],task_id=plan['task_id'],run_status=c['item']['status'],
            initialization='uniform',provenance='NEW_R14_FIXED_PREDICTION',
            source_file=plan['source_file'],source_view=plan.get('source_view'),
            secondary_fallback=plan.get('secondary_fallback'),selection_candidate=plan['candidate'],
            expected_public_candidate=plan.get('expected_public_candidate'),
            archived_selector_consistent=plan.get('archived_selector_consistent'),
            CPU_prefix=plan.get('CPU_prefix'),cutoff_CPU_seconds=plan.get('cutoff_CPU_seconds'),
            common_recorded_process_CPU=plan.get('common_recorded_process_CPU'),
            **metrics)
        for key in ('raw_original_source','raw_original_candidate','raw_original_Phi','initial_Phi','raw_initial_fallback','raw_audit_comparable','strict48_lower_deadline','strict48_upper_deadline','strict48_deadline_uncertainty_seconds','strict48_reason','strict48_scope','proposal_termination','common_physical_LU','common_physical_complex_RHS','actual_final_physics_counts','physics_counter_limits_available','normalization_scope'):
            if key in plan:row[key]=plan[key]
        if reason:row['availability_reason']=reason
        for key in ('CPU_seconds','wall_seconds','Phi','D','R','relative_residual','residual_over_delta',
                    'noise_feasible','legal','eligible','phase','deadline','lower_active','upper_active'):
            row[key]=meta.get(key)
        row['physical_cost_at_candidate']=physical_counts(meta.get('physics_cost',{}))
        contract=task.get('observation_contract',{})
        row['model_mismatch_status']=contract.get('model_mismatch_status')
        row['noise_certificate_eligible']=contract.get('model_mismatch_status')=='MATCHED_DISCRETE_MODEL'
        row['noise_qualification']='Original raw CN99 noise-only ball; cross-grid mismatch not covered'
        row['truth_available']=truths.get(plan['case_id']) is not None
        row['accuracy_available']=accuracy_available(row)
        row['recovered_posthoc']=recovered(row)
        output.append(row)
        if plan['view']=='primary' and m is not None:pictures[(plan['case_id'],plan['pair'])]=m
    return output,pictures


def comparisons(new_rows,old_rows):
    primary=[r for r in new_rows if r['view']=='primary']
    old={(r['case_id'],r['solver']):r for r in old_rows}
    new={(r['case_id'],r['solver']):r for r in primary}
    diffs=[]; sets=[]
    references=list(OLD_METHODS)+['pair'+p for p in PAIRS]
    for pair in PAIRS:
        solver='pair'+pair
        for ref in references:
            if ref==solver:continue
            gains=[];losses=[];both=[];neither=[];unavailable=[]
            for cid in CASES:
                a=new[(cid,solver)];b=(old if ref in OLD_METHODS else new)[(cid,ref)]
                paired=accuracy_available(a) and accuracy_available(b)
                row=dict(case_id=cid,solver=solver,reference=ref,paired_available=paired,
                         candidate_recovered=recovered(a),reference_recovered=recovered(b),
                         reference_provenance=b.get('provenance'))
                for metric in METRICS:
                    av,bv=a.get(metric),b.get(metric)
                    row['candidate_'+metric]=av;row['reference_'+metric]=bv
                    row['delta_'+metric]=av-bv if paired and finite(av) and finite(bv) else None
                diffs.append(row)
                if not paired:unavailable.append(cid)
                elif recovered(a) and not recovered(b):gains.append(cid)
                elif recovered(b) and not recovered(a):losses.append(cid)
                elif recovered(a):both.append(cid)
                else:neither.append(cid)
            sets.append(dict(solver=solver,reference=ref,expected_denominator=48,
                candidate_recovery_count=sum(recovered(new[(cid,solver)]) for cid in CASES),
                reference_recovery_count=sum(recovered((old if ref in OLD_METHODS else new)[(cid,ref)]) for cid in CASES),
                gain_cases=gains,loss_cases=losses,both_cases=both,neither_cases=neither,
                unavailable_paired_cases=unavailable,paired_available=48-len(unavailable),
                rule='POSTHOC_DESCRIPTIVE_ONLY; unavailable pairs excluded from gain/loss, retained in48 denominator'))
    return diffs,sets


def corrector_effects(rows,caches,np,pictures):
    lookup={(r['case_id'],r['solver'],r['view']):r for r in rows if r['view'] in ('raw','corrected','primary')}
    effects=[]
    for cid in CASES:
        for pair in PAIRS:
            solver='pair'+pair;a=lookup[(cid,solver,'raw')];b=lookup[(cid,solver,'corrected')]
            row=dict(case_id=cid,solver=solver,raw_available=a.get('available'),corrected_available=b.get('available'),
                raw_candidate=a.get('selection_candidate'),corrected_candidate=b.get('selection_candidate'),
                corrected_view_missing=not b.get('available'),material_L2_change=None,material_RMS_change=None)
            for key in METRICS:
                av,bv=a.get(key),b.get(key);row['raw_'+key]=av;row['corrected_'+key]=bv
                row['delta_'+key]=bv-av if finite(av) and finite(bv) else None
            if a.get('available') and b.get('available'):
                try:
                    with np.load(a['source_file'],allow_pickle=False) as za, np.load(b['source_file'],allow_pickle=False) as zb:
                        delta=zb['m']-za['m'];row['material_L2_change']=float(np.linalg.norm(delta));row['material_RMS_change']=float(np.sqrt(np.mean(delta**2)))
                except Exception as exc:row['material_change_error']=repr(exc)
            effects.append(row)
    return effects


def cost_rows(root,tasks,caches,physical_counts,issues):
    external={}
    for campaign in CAMPAIGNS:
        path=root/'queues'/campaign/'RECEIPT.json'
        if path.is_file():
            for job in read(path).get('jobs',[]):external[job['task_id']]=job
    result=[];history_rows=[];native_rows=[]
    for task in tasks:
        ident=task['task_id'];c=caches[ident];receipt=c['receipt'] or {};job=external.get(ident,{})
        value=(job.get('external') or {}).get('child_process_CPU_seconds',job.get('child_process_CPU_seconds'))
        physical=physical_counts(receipt.get('physical_cost',{}))
        result.append(dict(task_id=ident,case_id=task['case_id'],solver='pair'+str(task['pair']),campaign=task['campaign'],
            run_status=c['item']['status'],external_child_CPU_seconds=value,
            external_CPU_scope=(job.get('external') or {}).get('CPU_scope'),
            internal_final_process_CPU_seconds=receipt.get('CPU_seconds'),
            internal_wall_seconds=receipt.get('wall_seconds'),external_wall_seconds=job.get('wall_seconds'),
            selected_logical_CPU=job.get('cpu'),physical_cost=physical,auxiliary_cost=receipt.get('auxiliary_cost'),
            peak_sampled_RSS_bytes=receipt.get('peak_sampled_RSS_bytes'),
            RSS_scope=receipt.get('RSS_scope'),termination_reason=receipt.get('termination_reason'),
            cost_scope='External full child lifetime; internal timers/counters nested, not summed twice'))
        for h in c['history']:history_rows.append(dict(task_id=ident,case_id=task['case_id'],pair=task['pair'],**h))
        item=c['item']['native_history_file']
        if item['exists']:
            try:
                for row in read(item['path']):native_rows.append({**row,'task_id':ident,'case_id':task['case_id'],'pair':task['pair']})
            except Exception as exc:issues.append(dict(task_id=ident,issue='NATIVE_HISTORY_UNREADABLE',error=repr(exc)))
    return result,history_rows,native_rows


def plots(out,rows,old,pictures,truths,effects,costs,np):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    figures=out/'FIGURES';figures.mkdir()
    paths=[]
    primary={(r['case_id'],r['pair']):r for r in rows if r['view']=='primary'}
    pdf=figures/'all48_materials.pdf'
    with PdfPages(pdf) as pages:
        for offset in range(0,48,6):
            fig,axes=plt.subplots(6,5,figsize=(13,15),squeeze=False)
            for i,cid in enumerate(CASES[offset:offset+6]):
                truth=truths.get(cid);gt=truth['m'] if truth else None
                arrays=[gt]+[pictures.get((cid,p)) for p in PAIRS]
                vmax=max([float(np.max(x)) for x in arrays if x is not None and x.size]+[.1]);im=None
                for j,m in enumerate(arrays):
                    ax=axes[i,j];ax.set_xticks([]);ax.set_yticks([])
                    if m is not None and m.size==576:im=ax.imshow(m.reshape(24,24),vmin=.02,vmax=vmax,cmap='viridis',origin='upper')
                    else:ax.text(.5,.5,'MISSING',ha='center',va='center',transform=ax.transAxes)
                    if j==0:ax.set_ylabel(cid,fontsize=8)
                    if i==0:ax.set_title('GT' if j==0 else 'pair'+PAIRS[j-1])
                    if j and primary[(cid,PAIRS[j-1])].get('available'):
                        r=primary[(cid,PAIRS[j-1])]
                        nr=r.get('NRMSE');io=r.get('IoU')
                        ax.set_xlabel(('NRMSE '+format(nr,'.3f') if finite(nr) else 'NRMSE unavailable')+' | '+('IoU '+format(io,'.3f') if finite(io) else 'IoU unavailable'),fontsize=8)
                if im is not None:fig.colorbar(im,ax=axes[i,:].tolist(),fraction=.012,pad=.01)
            fig.suptitle('Existing exposed48 — frozen primary selections; row-shared material scale')
            fig.subplots_adjust(top=.95,bottom=.04,left=.08,right=.92,hspace=.24,wspace=.08)
            page=figures/f'materials_{offset//6+1:02d}.png';fig.savefig(page,dpi=150);pages.savefig(fig);plt.close(fig);paths.append(str(page.relative_to(out)))
    paths.append(str(pdf.relative_to(out)))
    methods=['pair'+p for p in PAIRS]+list(OLD_METHODS)
    lookup={(r['case_id'],r['solver']):r for r in rows if r['view']=='primary'}
    lookup.update({(r['case_id'],r['solver']):r for r in old})
    membership=np.array([[1 if recovered(lookup[(cid,s)]) else 0 if accuracy_available(lookup[(cid,s)]) else -1 for s in methods] for cid in CASES])
    fig,ax=plt.subplots(figsize=(12,13));im=ax.imshow(membership,aspect='auto',vmin=-1,vmax=1,cmap='RdYlGn')
    ax.set_xticks(range(len(methods)),methods,rotation=45,ha='right');ax.set_yticks(range(48),CASES,fontsize=6)
    ax.set_title('Posthoc recovery membership: −1 unavailable,0 fails criterion,1 passes; six historical references reused')
    fig.tight_layout();path=figures/'recovery_membership.png';fig.savefig(path,dpi=160);plt.close(fig);paths.append(str(path.relative_to(out)))
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for axis,metric in zip(axes.flat,['NRMSE','IoU','Phi','D']):
        for p in PAIRS:
            points=[r for r in summaries(rows) if r['solver']=='pair'+p and r['view']=='paid_prefix']
            points.sort(key=lambda x:x['CPU_prefix'])
            axis.plot([r['CPU_prefix'] for r in points],[r.get(metric+'_mean') if r.get(metric+'_mean') is not None else np.nan for r in points],marker='o',label='pair'+p)
        axis.set_xlabel('Native whole-process CPU cutoff (seconds)');axis.set_ylabel(metric+' mean (available cases)');axis.grid(alpha=.2)
    axes[0,0].legend();fig.suptitle('Paid-prefix curves; candidate chosen solely by original physical Phi')
    fig.tight_layout();path=figures/'paid_prefix_curves.png';fig.savefig(path,dpi=160);plt.close(fig);paths.append(str(path.relative_to(out)))
    fig,axes=plt.subplots(1,4,figsize=(16,4))
    for ax,key in zip(axes,['D','R','Phi','material_RMS_change']):
        for pi,pair in enumerate(PAIRS):
            field=key if key=='material_RMS_change' else 'delta_'+key
            values=[e[field] for e in effects if e['solver']=='pair'+pair and finite(e.get(field))]
            if values:ax.scatter([pi]*len(values),values,s=10,alpha=.5)
        ax.axhline(0,color='k',lw=.8);ax.set_xticks(range(4),['pair'+p for p in PAIRS]);ax.set_ylabel('material RMS movement' if key=='material_RMS_change' else 'corrected − raw '+key)
    fig.suptitle('Corrector effects on original physical objective terms');fig.tight_layout();path=figures/'corrector_effects.png';fig.savefig(path,dpi=160);plt.close(fig);paths.append(str(path.relative_to(out)))
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    for ax,metric in zip(axes,['NRMSE','IoU']):
        for pi,pair in enumerate(PAIRS):
            values=[]
            for cid in CASES:
                a=lookup[(cid,'pair'+pair)];b=lookup[(cid,'a_plus_arc_r7')]
                if a.get('available') and b.get('available') and finite(a.get(metric)) and finite(b.get(metric)):values.append(a[metric]-b[metric])
            if values:ax.scatter([pi]*len(values),values,s=10,alpha=.6)
        ax.axhline(0,color='k',lw=.8);ax.set_xticks(range(4),['pair'+p for p in PAIRS]);ax.set_ylabel(metric+' new − reused A+-Arc')
    fig.suptitle('Paired continuous metrics; no recovery threshold in selection');fig.tight_layout();path=figures/'paired_metric_differences.png';fig.savefig(path,dpi=160);plt.close(fig);paths.append(str(path.relative_to(out)))
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    for ax,metric in zip(axes,['NRMSE','IoU']):
        for pair in PAIRS:
            points=[r for r in summaries(rows) if r['solver']=='pair'+pair and r['view']=='common_CPU_prefix']
            points.sort(key=lambda r:r['CPU_prefix'])
            ax.plot([r['CPU_prefix'] for r in points],[r.get(metric+'_mean') if r.get(metric+'_mean') is not None else np.nan for r in points],marker='o',label='pair'+pair)
        ax.set_xlabel('Requested prefix, clipped per case to common native CPU endpoint (s)')
        ax.set_ylabel(metric+' mean (available cases)');ax.grid(alpha=.2)
    axes[0].legend();fig.suptitle('Four-way common recorded process-CPU comparisons; see per-case actual cutoffs')
    fig.tight_layout();path=figures/'common_CPU_prefix_curves.png';fig.savefig(path,dpi=160);plt.close(fig);paths.append(str(path.relative_to(out)))
    fig,ax=plt.subplots(figsize=(8,4))
    for pi,pair in enumerate(PAIRS):
        values=[r['external_child_CPU_seconds'] for r in costs if r['solver']=='pair'+pair and finite(r.get('external_child_CPU_seconds'))]
        if values:ax.scatter([pi]*len(values),values,s=10,alpha=.5)
    ax.set_xticks(range(4),['pair'+p for p in PAIRS]);ax.set_ylabel('External whole-child lifetime CPU seconds')
    ax.set_title('Actual external costs, including imports/export; missing measurements omitted explicitly')
    fig.tight_layout();path=figures/'external_CPU_costs.png';fig.savefig(path,dpi=160);plt.close(fig);paths.append(str(path.relative_to(out)))
    return paths


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--history-root',type=Path)
    parser.add_argument('--seal',type=Path)
    parser.add_argument('--out',type=Path)
    parser.add_argument('--make-seal',action='store_true')
    parser.add_argument('--source-note',default='')
    args=parser.parse_args();root=args.root.resolve();seal_path=(args.seal or root/'audit/PREDICTION_SEAL.json').resolve()
    if args.make_seal:
        value=make_seal(root,seal_path,args.source_note)
        print(json.dumps(dict(status='SEALED_NOT_AUTHORIZED',seal=str(seal_path),expected_tasks=value['expected_tasks'])))
        return 0
    if args.history_root is None or args.out is None:parser.error('--history-root and --out are required for offline scoring')
    history=args.history_root.resolve();out=args.out.resolve();start_cpu=time.process_time();start_wall=time.perf_counter()
    seal,tasks,protocol,account=authorize(root,seal_path)
    if out.exists():raise FileExistsError('New analysis output directory required')
    if out==root or out==history or out.is_relative_to(history) or any(out.is_relative_to(root/name) for name in ('RESULTS','queues','sealed_truth','data_public','src','configs','audit')):
        raise ValueError('Analysis output must not overlap historical reports or new predictions/data/source/seals')
    out.mkdir(parents=True,exist_ok=False)
    issues=[]
    plans,caches=public_plans(root,tasks,account,issues)
    # This output is fixed before even importing the truth helper, and records every GT-free selector.
    write_json(out/'PUBLIC_SELECTIONS_BEFORE_TRUTH.json',plans)
    raw_public=[r for r in plans if r['view']=='raw']
    write_json(out/'PUBLIC_RAW_INITIAL_COUNTS_BEFORE_TRUTH.json',dict(expected_tasks=192,audited_rows=len(raw_public),
        comparable=sum(r.get('raw_audit_comparable') is True for r in raw_public),
        initial_fallbacks=sum(r.get('raw_initial_fallback') is True for r in raw_public),
        unavailable_or_uncomparable=sum(r.get('raw_audit_comparable') is not True for r in raw_public),truth_access='NONE'))
    old,old_path=historical_rows(history,issues)
    import numpy as np
    helper_path=history/'src/score_results.py'
    sys.path.insert(0,str(history/'src'))
    spec=importlib.util.spec_from_file_location('r14_reused_offline_metrics',helper_path)
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    truths={};truth_sources=[]
    task_by_case={t['case_id']:t for t in tasks}
    for cid in CASES:
        # Do not use R14 data paths to resolve historical scoring truth.
        historical_task=dict(case_id=cid)
        spec_truth=helper.truth_spec(history,historical_task)
        truth_sources.append(dict(case_id=cid,**spec_truth))
        try:truths[cid]=helper.load_truth(spec_truth)
        except Exception as exc:truths[cid]=None;issues.append(dict(case_id=cid,issue='TRUTH_UNAVAILABLE_AFTER_AUTHORIZATION',error=repr(exc)))
    rows,pictures=score_plans(plans,caches,truths,helper.material_metrics,helper.physical_counts,np,issues)
    diffs,sets=comparisons(rows,old)
    effects=corrector_effects(rows,caches,np,pictures)
    costs,histories,native=cost_rows(root,tasks,caches,helper.physical_counts,issues)
    summary=summaries(rows+old)
    raw_audit=[{k:v for k,v in r.items() if k in ('task_id','case_id','solver','raw_original_source','raw_original_candidate','raw_original_Phi','initial_Phi','raw_initial_fallback','raw_audit_comparable','source_view','selection_candidate')} for r in rows if r['view']=='raw']
    table(out,'RAW_INITIAL_AUDIT',raw_audit)
    write_json(out/'RAW_INITIAL_AUDIT_COUNTS.json',dict(expected_tasks=192,audited_rows=len(raw_audit),comparable=sum(r.get('raw_audit_comparable') is True for r in raw_audit),initial_fallbacks=sum(r.get('raw_initial_fallback') is True for r in raw_audit),unavailable_or_uncomparable=sum(r.get('raw_audit_comparable') is not True for r in raw_audit)))
    for name,data in [('ALL_VIEWS',[r for r in rows if r['view'] in VIEWS]),
        ('PAID_PREFIXES',[r for r in rows if r['view']=='paid_prefix']),
        ('COMMON_CPU_PREFIXES',[r for r in rows if r['view']=='common_CPU_prefix']),
        ('STRICT48_PAID_POOL_BRACKETS',[r for r in rows if r['view'] in ('strict48_lower','strict48_upper')]),
        ('COMMON_PHYSICS_PAID_POOL',[r for r in rows if r['view']=='common_physics_pool']),
        ('HISTORICAL_PRIMARY',old),('PAIRED_METRIC_DIFFERENCES',diffs),
        ('RECOVERY_GAIN_LOSS',sets),('CORRECTOR_EFFECTS',effects),('COST_LEDGER',costs),
        ('PAID_CANDIDATE_HISTORY',histories),('NATIVE_STAGE_HISTORY',native),
        ('SUMMARY',summary),('ISSUES',issues),('TRUTH_SOURCES',truth_sources)]:table(out,name,data)
    artifacts=plots(out,rows,old,pictures,truths,effects,costs,np)
    measured=[r['external_child_CPU_seconds'] for r in costs]
    receipt=dict(status='AUTHORIZED_OFFLINE_SCORING_COMPLETE',schema='r14_offline_score_v1',
        root=str(root),history_root=str(history),seal=str(seal_path),seal_schema=seal['schema_version'],
        expected_material_objects=48,expected_new_tasks=192,new_view_rows=960,new_prefix_rows=960,
        actual_new_view_rows=sum(r['view'] in VIEWS for r in rows),actual_new_prefix_rows=sum(r['view']=='paid_prefix' for r in rows),
        historical_primary_rows=len(old),historical_source=stamp(old_path),historical_filter=dict(initialization='uniform',view='primary',solvers=list(OLD_METHODS)),
        protocol_snapshot=protocol,task_status_counts=dict(Counter(c['item']['status'] or 'MISSING' for c in caches.values())),
        issues=len(issues),new_physical_LU=0,new_RHS=0,new_optimization=0,
        scope='EXISTING_EXPOSED48_NOT_FRESH_BLIND',scientific_review='PARENT_REQUIRED',
        threshold_rule='NRMSE<=.15 and IoU>=.70 only posthoc descriptive membership; not online selection/stopping',
        selection='Archived saved views; prefixes min eligible legal completed original Phi from public histories fixed before truth',
        secondary='Missing saved secondary explicitly falls back to saved primary; corrected missing remains unavailable',
        common_CPU='Min native receipt final CPU_seconds across four methods, clipped to120; prefix applied on the same native CPU clock; external child CPU retained separately, no exact per-candidate external clock claimed',
        model='576 unknown real chi; original scalar-TM finite disk-pulse VIE, actual2x2 Green blocks, original edge smooth/TV prior',
        equal_physics='Componentwise same-case min of actual final original Maxwell LU/RHS across four methods;eligible paid candidates within both limits only;auxiliary FFT/block actions/SVD/CPU separate,not normalized or free;missing counters null',
        budget_deviation='Main R14 quartet actual absolute120 total with unused proposal reallocated; historical C2 references had48corrector limit. Strict48 pool lower/upper cutoffs are robustness views, not actual stopped reruns; full CPU remains charged',
        raw_reporting='Raw min archived eligible initial/proposal by public Phi; original raw NPZ untouched and original path/index/Phi retained',
        stage_scope='M1 beta stages1/.2/.2/.2; matched metrics define fixed-stage loss only; stage switches do not define one global potential',
        noise='Original raw circular-complex CN99 ball; cross-grid model mismatch not noise qualified',
        external_child_CPU_total=sum(measured) if all(finite(x) for x in measured) else None,
        measured_external_child_CPU_subtotal=sum(x for x in measured if finite(x)),
        missing_external_child_CPU_tasks=[r['task_id'] for r in costs if not finite(r['external_child_CPU_seconds'])],
        CPU_seconds=time.process_time()-start_cpu,wall_seconds=time.perf_counter()-start_wall,
        scoring_cost_scope='New offline metrics/plots only; solver external costs separate; no optimizer/forward',
        artifacts=artifacts,provenance='File existence/byte sizes, no SHA or mtime validation; same-size replacements are not detectable')
    write_json(out/'RECEIPT.json',receipt)
    body=['<!doctype html><meta charset="utf-8"><title>R14 dual-metric offline comparison</title>',
        '<style>body{font:16px system-ui;max-width:1400px;margin:auto;padding:30px}img{max-width:100%}table{border-collapse:collapse}td,th{padding:7px;border:1px solid #ccc}</style>',
        '<h1>R14 双度量：已曝光48对象离线比较</h1>',
        '<p>四组合的新封存预测与六组历史 uniform48 primary 参考。科学结论由主线程复核。所有缺失/失败仍保留48分母；历史数据并非新盲测。</p>',
        '<p>预算偏差：新四臂采用总CPU120、可重分配proposal空余时间；旧C2参考corrector最多48秒。原48口径只作已支付候选池的保守上下界截取，不是实际重跑；全部费用保留。</p>',
        '<p>NRMSE≤.15 且 IoU≥.70 仅为事后描述性统计；从未用于候选选择。连续误差、原始物理目标、噪声球资格及费用分别列出。</p>',
        '<table><tr><th>Method</th><th>Available/48</th><th>Posthoc recovered/48</th><th>Mean NRMSE</th><th>Mean IoU</th></tr>']
    for row in summary:
        if row['view']=='primary':body.append('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in [row['solver'],row['available_cases'],row['recovered_posthoc_count'],row.get('NRMSE_mean'),row.get('IoU_mean')])+'</tr>')
    body.append('</table><p>Missing corrected views stay missing; absent secondary uses archived primary with an explicit fallback flag. Common CPU comparisons use native process-clock cutoffs; external child lifetime CPU remains a separate cost.</p>')
    for name in ('SUMMARY','ALL_VIEWS','PAID_PREFIXES','COMMON_CPU_PREFIXES','HISTORICAL_PRIMARY','PAIRED_METRIC_DIFFERENCES','RECOVERY_GAIN_LOSS','CORRECTOR_EFFECTS','COST_LEDGER','PAID_CANDIDATE_HISTORY','NATIVE_STAGE_HISTORY','ISSUES'):
        body.append(f'<a href="{name}.csv">{name} CSV</a> · <a href="{name}.json">JSON</a><br>')
    for artifact in artifacts:
        body.append(f'<p><a href="{html.escape(artifact)}">{html.escape(artifact)}</a></p>')
        if artifact.endswith('.png'):body.append(f'<img src="{html.escape(artifact)}">')
    with (out/'report.html').open('x',encoding='utf-8') as stream:stream.write('\n'.join(body))
    print(json.dumps(dict(status=receipt['status'],out=str(out),issues=len(issues),new_physical_LU=0)))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
