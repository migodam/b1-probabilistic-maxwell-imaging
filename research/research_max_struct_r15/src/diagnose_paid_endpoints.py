"""Parent-authorized no-GT validation of 44 fixed, already-selected primary endpoints.

No optimization, candidate generation, material Jacobian, selection or scoring.
The saved gradients and fields are validation artifacts, never online inputs.
"""
from __future__ import annotations

import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
import argparse
from collections import Counter
from datetime import datetime, timezone
import csv
import json
import math
from pathlib import Path
import shutil
import sys
import time
import traceback

sys.dont_write_bytecode = True
CAMPAIGN = 'twofold_exposed24_v1'
SEAL_SCHEMA = 'r15_twofold_exposed24_prediction_seal_v1'
CASES = tuple(f'r14_new_{i:04d}' for i in (9, 20, 21, 23))
METHODS = ('g19', 'g19_s19', 'twofold', 'twofold_hop', 'twofold_dr', 'twofold_ap',
           'twofold_full_dr', 'twofold_reset', 'cie', 'dual_cie', 'dr')
VIEWS = ('initial', 'raw', 'corrected', 'primary', 'secondary')
SCOPE = 'POSTSEAL_PAID_PRIMARY_VALIDATION_NO_GT_NO_SELECTION_FEEDBACK'


def checked(root, value):
    path = Path(value)
    path = (path if path.is_absolute() else root / path).resolve()
    if not path.is_relative_to(root) or any(x.lower() in
        ('sealed_truth', 'truth', 'data_private', 'private', 'private_operational', 'credentials', 'tokens') for x in path.parts):
        raise ValueError('Unsafe/private diagnostic path')
    return path


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def stamp(path):
    return dict(path=str(path.resolve()), exists=path.is_file(), bytes=path.stat().st_size if path.is_file() else None)


def validate_stamp(root, item):
    path = checked(root, item['path'])
    if stamp(path) != item:
        raise ValueError('Sealed original existence/byte size changed: ' + str(path))
    return path


def authorize(root, path):
    """Match the actual analyze_campaign.py prediction seal; no scoring import."""
    seal = read(path)
    if seal.get('schema_version') != SEAL_SCHEMA or seal.get('authorized') is not True:
        raise PermissionError('Require the NEW parent-authorized R15 main264 prediction seal')
    if (Path(seal['root']).resolve() != root or seal.get('campaign') != CAMPAIGN or
        seal.get('expected_tasks') != 264 or seal.get('expected_cases') != 24 or
        set(seal.get('fixed_predictions', [])) != set(VIEWS)):
        raise ValueError('Prediction seal root/campaign/matrix mismatch')
    cfg = read(validate_stamp(root, seal['manifest_file']))
    tasks = cfg['tasks']
    wanted = {(f'r14_new_{i:04d}', method) for i in range(24) for method in METHODS}
    if (cfg != seal['manifest_snapshot'] or cfg.get('campaign') != CAMPAIGN or len(tasks) != 264 or
        {(t['case_id'], t['method']) for t in tasks} != wanted or
        len({t['task_id'] for t in tasks}) != 264 or set(cfg.get('methods', [])) != set(METHODS)):
        raise ValueError('Original frozen24 x11 manifest changed')
    originals = {t['task_id']: t for t in tasks}
    jobs = {}
    for entry in seal['queue_receipts']:
        queue_path = validate_stamp(root, entry['file'])
        if not queue_path.is_relative_to(root / 'queues' / CAMPAIGN):
            raise ValueError('Seal queue outside current main campaign')
        queue = read(queue_path)
        if (queue != entry['snapshot'] or queue.get('status') not in ('COMPLETED', 'COMPLETED_WITH_ERRORS') or
            queue.get('running') or queue.get('pending_tasks', 0) or queue.get('pending_cases', 0)):
            raise ValueError('Main queue not terminal or changed after seal')
        for job in queue.get('jobs', []):
            if job['task_id'] not in originals:
                raise ValueError('Unknown terminal task')
            jobs[job['task_id']] = job
    if set(jobs) != set(originals):
        raise ValueError('All264 terminal jobs, including failures, required')
    account = {}
    for item in seal['tasks']:
        ident = item['task_id']
        if ident in account or ident not in originals:
            raise ValueError('Unknown/duplicate sealed task')
        task = originals[ident]
        if ident != task['case_id'] + '__' + task['method'] or task['campaign'] != CAMPAIGN:
            raise ValueError('Changed task identity')
        folder = root / 'RESULTS' / CAMPAIGN / ident
        actual = [stamp(p) for p in sorted(folder.rglob('*')) if p.is_file()] if folder.is_dir() else []
        if actual != item['all_result_files']:
            raise ValueError('Original result inventory changed since seal: ' + ident)
        validate_stamp(root, item['receipt_file'])
        for view in VIEWS:
            entry = item['predictions'][view]
            prediction = validate_stamp(root, {k: entry[k] for k in ('path', 'exists', 'bytes')})
            if prediction != folder / (view + '.npz'):
                raise ValueError('Unexpected sealed prediction path')
            validate_stamp(root, entry['metadata'])
        if checked(root, task['input_file']) != root / 'data_public/independent24' / (task['case_id'] + '.npz'):
            raise ValueError('Unexpected public input path')
        account[ident] = item
    if set(account) != set(originals):
        raise ValueError('All264 sealed task rows required; missing/failed may not be omitted')
    return seal, cfg, account, jobs


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if hasattr(value, 'item'):
        return clean(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def dump(path, value):
    path.write_text(json.dumps(clean(value), indent=2, allow_nan=False, ensure_ascii=False) + '\n', encoding='utf-8')


def table(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(clean(v), ensure_ascii=False) if isinstance(v, (dict, list)) else clean(v) for k, v in row.items()})


def delta(after, before):
    return {k: after.get(k, 0) - before.get(k, 0) for k in set(after) | set(before)
            if after.get(k, 0) != before.get(k, 0)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'legacy-repo', 'seal'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    cpu, wall = time.process_time(), time.perf_counter()
    root, legacy = args.root.resolve(), args.legacy_repo.resolve()
    seal_path = checked(root, args.seal)
    out = checked(root, args.out or root / 'audit/paid_endpoint_diagnostic_v1')
    if not out.is_relative_to(root / 'audit') or out.exists():
        raise FileExistsError('Choose a NEW audit output directory; original predictions remain sealed')
    seal, cfg, account, jobs = authorize(root, seal_path)
    selected = [dict(case_id=case, method=method, task_id=case + '__' + method) for case in CASES for method in METHODS]
    out.mkdir(parents=True, exist_ok=False)
    (out / 'states').mkdir()
    shutil.copyfile(seal_path, out / 'PARENT_SEAL_AT_ACTIONS.json')
    dump(out / 'FIXED_PLAN.json', dict(scope=SCOPE, expected=44, cases=CASES, methods=METHODS, items=selected,
        seal=str(seal_path), authorized=True, authorized_scoring_flag=seal.get('authorized_scoring'),
        truth_access='NONE', selection='Use existing sealed primary.npz m only; no new candidate selection'))
    rows, events, failures, case_costs = [], [], [], []
    total = Counter()
    active_case_id = None
    try:
        # Same API/import routes as the actual R15 runners; no alternate backend.
        sys.path[:0] = [str(legacy), str(legacy / 'vendor'), str(root.parents[1])]
        import numpy as np
        from research.research_fbg_r3_experiments.src.physics import MaxwellModel
        physics_source = str(Path(sys.modules[MaxwellModel.__module__].__file__).resolve())
        for case in CASES:
            active_case_id = case
            case_cpu, case_wall = time.process_time(), time.perf_counter()
            model, public, setup_error = None, None, None
            try:
                with np.load(checked(root, 'data_public/independent24/' + case + '.npz'), allow_pickle=False) as z:
                    public = {key: z[key].copy() for key in z.files}
                model = MaxwellModel(public)
                Y = np.asarray(public['Y'])
                sigma = float(public['sigma']) if 'sigma' in public else None
                noise_delta = float(public['delta']) if 'delta' in public else None
                lo, hi = float(public.get('lower', .02)), float(public.get('upper', 20.))
                if model.n != 576 or public.get('components', 1) != 1 or lo >= hi:
                    raise ValueError('Expected original576 scalar-material model and box')
                if not np.isfinite(Y).all() or float(np.vdot(Y, Y).real) <= 0:
                    raise ValueError('Invalid public observation normalization')
            except Exception as exc:
                setup_error = repr(exc)
                failures.append(dict(case_id=case, method='CASE_SETUP', error=setup_error, traceback=traceback.format_exc()))
            for item in [t for t in selected if t['case_id'] == case]:
                ident, method = item['task_id'], item['method']
                start_cpu, start_wall = time.process_time(), time.perf_counter()
                before = Counter(model.cost) if model is not None else Counter()
                before_timings = dict(model.timings) if model is not None else {}
                before_cpu_timings = dict(model._physics.cpu_timings) if model is not None else {}
                row = dict(item, status='RUNNING', scope=SCOPE, diagnostic_available=False,
                    D=None, R=None, Phi=None, g_data_L2=None, g_prior_L2=None, g_total_L2=None,
                    projected_gradient_L2=None, projected_gradient_Linf=None, lower_active=None, upper_active=None,
                    sigma_complex_RMS=None, delta=None, residual_over_delta=None, residual_over_sigma_sqrt_ports=None,
                    original_termination_reason=None, original_whole_process_CPU_seconds=None,
                    original_queue_child_CPU_seconds=jobs[ident].get('external', {}).get('child_process_CPU_seconds'),
                    original_queue_status=jobs[ident].get('status'), unknown_fields_are_null=True)
                try:
                    if setup_error:
                        raise RuntimeError('Case setup failed: ' + setup_error)
                    record = account[ident]
                    rp = validate_stamp(root, record['receipt_file'])
                    original = read(rp) if rp.is_file() else {}
                    row.update(original_child_status=original.get('status'),
                        original_termination_reason=original.get('termination_reason', original.get('termination')),
                        original_whole_process_CPU_seconds=original.get('whole_process_CPU_seconds'),
                        original_internal_CPU_seconds=original.get('CPU_seconds'), original_wall_seconds=original.get('wall_seconds'),
                        original_physical_cost=original.get('physical_cost', original.get('physics_cost')),
                        original_primary_candidate=original.get('primary_candidate'))
                    if original and (original.get('truth_access') != 'NONE' or original.get('task', {}).get('task_id') != ident):
                        raise ValueError('Original receipt task/truth declaration mismatch')
                    entry = record['predictions']['primary']
                    prediction = validate_stamp(root, {k: entry[k] for k in ('path', 'exists', 'bytes')})
                    row['primary_path'] = str(prediction)
                    if not prediction.is_file():
                        row.update(status='MISSING_PRIMARY', missing_reason='Sealed primary does not exist; no fallback generated')
                        continue
                    with np.load(prediction, allow_pickle=False) as z:
                        raw_m = z['m']
                        if np.iscomplexobj(raw_m):
                            raise ValueError('Primary m is complex; no implicit real projection')
                        m = np.asarray(raw_m, dtype=float).reshape(model.n).copy()
                    if not np.isfinite(m).all():
                        raise ValueError('Nonfinite primary material')
                    prior_contract = original.get('prior') or {}
                    for key in ('smooth', 'tv', 'epsilon'):
                        if prior_contract.get(key) is not None and float(prior_contract[key]) != float(getattr(model.reg, key)):
                            raise ValueError('Inherited physical prior differs from original receipt: ' + key)
                    forward_before = Counter(model.cost)
                    state = model.state(m)
                    forward_cost = delta(model.cost, forward_before)
                    cache_before = Counter(model.cost)
                    repeat = model.state(m.copy())
                    cache_cost = delta(model.cost, cache_before)
                    cache_same = repeat is state and repeat['lu'] is state['lu'] and not cache_cost
                    if not cache_same:
                        raise RuntimeError('Original state/LU cache was not reused; charged failure retained')
                    gradient_before = Counter(model.cost)
                    phi, g_total = model.data_value_gradient(m, Y)
                    gradient_cost = delta(model.cost, gradient_before)
                    if gradient_cost.get('LU', 0) or gradient_cost.get('FULL_FORWARD', 0):
                        raise RuntimeError('data_value_gradient unexpectedly refactorized cached state; cost retained')
                    R, g_prior, _ = model.reg.evaluate(m)
                    g_total = np.asarray(g_total, dtype=float)
                    g_prior = np.asarray(g_prior, dtype=float)
                    g_data = g_total - g_prior
                    residual = state['Y'] - Y
                    rr2, yy2 = float(np.vdot(residual, residual).real), float(np.vdot(Y, Y).real)
                    D = .5 * rr2 / yy2
                    projected_gradient = m - np.clip(m - g_total, lo, hi)
                    lower_mask, upper_mask = m <= lo + 1e-8, m >= hi - 1e-8
                    values = dict(m=m, J=state['j'], Q=state['Q'], Y_pred=state['Y'], Y_observed=Y,
                        residual=residual, g_data=g_data, g_prior=g_prior, g_total=g_total,
                        box_projected_gradient=projected_gradient, lower_active_mask=lower_mask,
                        upper_active_mask=upper_mask, box=np.array([lo, hi]))
                    if not all(np.isfinite(v).all() for v in values.values()):
                        raise FloatingPointError('Nonfinite diagnostic arrays')
                    np.savez_compressed(out / 'states' / (ident + '.npz'), **values)
                    row.update(status='COMPLETE', diagnostic_available=True, D=D, R=float(R), Phi=float(phi),
                        Phi_minus_D_minus_R=float(phi-D-R), physical_relative_residual=float(np.sqrt(rr2/yy2)),
                        physical_model_relative_residual=state.get('state_relative_residual'),
                        residual_raw_L2=float(np.sqrt(rr2)), observation_energy=yy2, sources=model.ns, receivers=model.nr,
                        material_dof=model.n, material_min=float(m.min()), material_max=float(m.max()), box_lower=lo, box_upper=hi,
                        box_active_tolerance=1e-8, lower_active=int(lower_mask.sum()), upper_active=int(upper_mask.sum()),
                        box_violations=int(np.sum((m < lo) | (m > hi))),
                        sigma_complex_RMS=sigma, delta=noise_delta,
                        sigma_over_delta=float(sigma/noise_delta) if sigma is not None and noise_delta is not None and noise_delta > 0 else None,
                        residual_over_delta=float(np.sqrt(rr2)/noise_delta) if noise_delta is not None and noise_delta > 0 else None,
                        residual_over_sigma_sqrt_ports=float(np.sqrt(rr2)/(sigma*np.sqrt(Y.size))) if sigma is not None and sigma > 0 else None,
                        delta_over_sigma_sqrt_ports=float(noise_delta/(sigma*np.sqrt(Y.size))) if sigma is not None and sigma > 0 and noise_delta is not None else None,
                        noise_qualification=str(public['sigma_role']) if 'sigma_role' in public else original.get('noise_qualification'),
                        noise_ratio_null_reason=None if sigma is not None and sigma > 0 and noise_delta is not None and noise_delta > 0 else 'Public sigma/delta missing or nonpositive; no threshold fabricated',
                        projected_gradient_L2=float(np.linalg.norm(projected_gradient)),
                        projected_gradient_Linf=float(np.max(np.abs(projected_gradient))),
                        projected_gradient_step_scale=1., forward_cost=forward_cost, gradient_cost=gradient_cost,
                        inherited_prior=dict(name=getattr(model.reg, 'name', None), smooth=getattr(model.reg, 'smooth', None),
                                             tv=getattr(model.reg, 'tv', None), epsilon=getattr(model.reg, 'epsilon', None)),
                        retained_LU_reused=cache_same, repeated_state_new_cost=cache_cost,
                        saved_arrays=str((out / 'states' / (ident + '.npz')).relative_to(root)))
                    for name, gradient in [('g_data', g_data), ('g_prior', g_prior), ('g_total', g_total)]:
                        row[name + '_L2'] = float(np.linalg.norm(gradient))
                        row[name + '_Linf'] = float(np.max(np.abs(gradient)))
                except Exception as exc:
                    row.update(status='FAILED', error=repr(exc))
                    failures.append(dict(item, error=repr(exc), traceback=traceback.format_exc()))
                finally:
                    row.update(CPU_seconds=time.process_time()-start_cpu, wall_seconds=time.perf_counter()-start_wall,
                        physical_cost=delta(model.cost, before) if model is not None else {},
                        physical_wall_timings=delta(model.timings, before_timings) if model is not None else {},
                        physical_CPU_timings=delta(model._physics.cpu_timings, before_cpu_timings) if model is not None else {})
                    rows.append(row)
                    events.append(dict(item, status=row['status'], CPU_seconds=row['CPU_seconds'], wall_seconds=row['wall_seconds'], physical_cost=row['physical_cost']))
                    dump(out / 'PROGRESS.json', dict(expected=44, terminal=len(rows), completed=sum(r['status']=='COMPLETE' for r in rows), last=ident))
            counts = dict(model.cost) if model is not None else {}
            total.update(counts)
            case_costs.append(dict(case_id=case, setup_error=setup_error, physical_cost=counts,
                CPU_seconds=time.process_time()-case_cpu, wall_seconds=time.perf_counter()-case_wall))
            active_case_id = None
        table(out / 'ENDPOINTS.csv', rows)
        dump(out / 'ENDPOINTS.json', rows)
        dump(out / 'FAILURES.json', failures)
        table(out / 'COST_EVENTS.csv', events)
        dump(out / 'CASE_COSTS.json', case_costs)
        receipt = dict(status='COMPLETE_WITH_RETAINED_ISSUES' if failures or any(r['status']!='COMPLETE' for r in rows) else 'COMPLETE',
            scope=SCOPE, expected=44, terminal=len(rows), completed=sum(r['status']=='COMPLETE' for r in rows),
            missing=sum(r['status']=='MISSING_PRIMARY' for r in rows), failed=sum(r['status']=='FAILED' for r in rows),
            denominator=44, seal=str(seal_path), physics_source=physics_source, physical_cost=dict(total),
            CPU_seconds=time.process_time()-cpu, whole_process_CPU_at_receipt_seconds=time.process_time(),
            wall_seconds=time.perf_counter()-wall, case_costs=case_costs, BLAS_threads=1,
            truth_access='NONE', full_material_Jacobians=0, optimization_calls=0, candidate_generation=0,
            primary_selection_changed=False, validation_feedback_to_online_selection=False,
            cost_scope='Physical cost sums each case model once including setup/failures. Endpoint costs are nested deltas; ns adjoint RHS reuse retained forward LU. CPU case/endpoint/whole process are nested and must not be added.')
        dump(out / 'RECEIPT.json', receipt)
        lines = ['# Paid primary endpoint validation (no GT)', '',
            f'Fixed denominator: 4 cases ×11 main methods =44. Complete {receipt["completed"]}; missing {receipt["missing"]}; failed {receipt["failed"]}. Missing/failed endpoints remain in ENDPOINTS.csv and are not replaced.', '',
            'Only the new parent-authorized main264 prediction seal is accepted. Each endpoint uses the existing primary.npz real m. No candidate, optimizer, finite step experiment, GT score, or selection feedback is produced.', '',
            'Original MaxwellModel.state(m) supplies J/Q/Y_pred and a retained LU. One repeated identical-m state call checks object/LU identity and zero added physical counters; data_value_gradient then reuses this LU for one ns-RHS adjoint. All actual counters, including a failed operation or unexpected extra factorization, remain charged. No full576 material Jacobian is requested.', '',
            'D=.5||Y_pred−Y||_F²/||Y||_F²; R is the inherited smooth/TV prior; Phi and total gradient come from original data_value_gradient. g_prior comes from the unchanged model.reg.evaluate; g_data=g_total−g_prior. Unit-step box gradient mapping is m−clip(m−g_total,lower,upper); it is a diagnostic with declared scale1, not a new optimizer or candidate. Box activity uses inherited1e-8 reporting tolerance.', '',
            'Raw sigma, delta and residual ratios use the original public fields. Missing/nonpositive denominators yield null, never an invented noise threshold. Unknown termination/cost fields stay null. Original endpoint selection CPU and new validation CPU are separate; OS child queue cost is retained as provenance and not refunded.', '',
            'A small gradient or projected gradient does not establish a strict/local minimum, global optimality, information nullspace, identifiability or recovery. Scientific interpretation belongs to the parent.', '',
            'Arrays: states/<case>__<method>.npz contains m/J/Q/Y_pred/Y_observed/residual, g_data/g_prior/g_total, box gradient mapping and active masks; no truth arrays. Reports: ENDPOINTS.csv/json, COST_EVENTS.csv, CASE_COSTS.json, FAILURES.json, RECEIPT.json.', '',
            f'Actual validation CPU at receipt: {receipt["CPU_seconds"]:.9f}s; whole-process CPU at receipt {receipt["whole_process_CPU_at_receipt_seconds"]:.9f}s; wall {receipt["wall_seconds"]:.9f}s. Physical counters: {json.dumps(dict(total), sort_keys=True)}', '']
        (out / 'REPORT.md').write_text('\n'.join(lines), encoding='utf-8')
        print(json.dumps({k: receipt[k] for k in ('status', 'expected', 'completed', 'missing', 'failed', 'physical_cost', 'CPU_seconds', 'wall_seconds')}))
    except BaseException as exc:
        dump(out / 'FAILED_RECEIPT.json', dict(status='FAILED_RETAIN_PARTIAL_VALIDATION', error=repr(exc), rows=rows,
            case_costs=case_costs, completed_case_physical_cost=dict(total),
            active_case_id=active_case_id,
            active_case_physical_cost=dict(model.cost) if active_case_id is not None and 'model' in locals() and model is not None else {},
            CPU_seconds=time.process_time()-cpu, whole_process_CPU_at_error_seconds=time.process_time(),
            wall_seconds=time.perf_counter()-wall, truth_access='NONE', optimization=0))
        raise


if __name__ == '__main__':
    main()
