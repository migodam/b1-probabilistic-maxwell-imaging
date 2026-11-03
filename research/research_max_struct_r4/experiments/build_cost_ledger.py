"""Read receipts only; no NumPy, solver, truth, hashes, or scientific execution.

Run rows are additive only once per recorded invocation. Component rows are
nonadditive views of those runs. A receipt cannot prove unrecorded retries absent.
"""
from pathlib import Path
import argparse
import csv
import json
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
COST_KEYS = ('physical_model_cost', 'physical_cost', 'physical_and_dense_cost')
SMALL_KEYS = ('small4x4_factorizations', 'small4x4_complex_RHS',
              'small4x4_matrix_apply_complex_RHS', 'scalar_complex_divisions',
              'onehop_block_matrix_complex_RHS')


def read(path):
    return json.loads(path.read_text())


def make_row(path, d, phase='mixed', charged=True):
    cost = next((d[k] for k in COST_KEYS if isinstance(d.get(k), dict)), None)
    if cost is None and any(k in d for k in ('new_physical_LU', 'physical_LU', 'LU')):
        cost = dict(LU=d.get('new_physical_LU', d.get('physical_LU', d.get('LU'))),
                    actual_complex_RHS=d.get('new_complex_RHS', d.get('physical_RHS', d.get('actual_complex_RHS'))))
    cost = cost or {}
    aux = d.get('auxiliary_operation_counts', d.get('auxiliary_cost', {}))
    aux = aux if isinstance(aux, dict) else {}
    wiring = d.get('gradient_check', {})
    wiring_ops = wiring.get('operation_counts', {}) if isinstance(wiring, dict) else {}
    locality_lu = cost.get('locality_LU', 0)
    locality_rhs = cost.get('locality_dense_RHS', 0)
    total_lu, total_rhs = cost.get('LU'), cost.get('actual_complex_RHS')
    row = dict(run_id=d.get('run_id', str(path.relative_to(ROOT))), receipt=str(path.relative_to(ROOT)),
        status=d.get('status', 'RECORDED'), phase=phase, charge_recorded_invocation=charged,
        physical_LU=(total_lu-locality_lu if total_lu is not None else None),
        model_complex_RHS_excluding_locality=(total_rhs-locality_rhs if total_rhs is not None else None),
        counter_LU_including_locality=total_lu, counter_actual_complex_RHS=total_rhs,
        nonlocal_locality_LU=locality_lu, nonlocal_locality_RHS=locality_rhs,
        wall_seconds=d.get('wall_seconds'), cpu_seconds=d.get('cpu_seconds'),
        amount=None, currency=None, hash_verified=False,
        timer_scope=d.get('timer_scope', d.get('cost_scope', 'receipt timer scope not separately documented')),
        physical_factor_dimension_sum=cost.get('full_factorization_dimension'),
        full_counts=cost,
        failure_details={k:d[k] for k in ('failed_attempts','failures','error','diagnostic_failed_attempts','remote_attempts') if k in d},
        notes=[])
    for k in SMALL_KEYS:
        # Embedded gradient check is independent of the recorded auxiliary counters.
        row[k] = aux.get(k, 0) + wiring_ops.get(k, 0)
    if 'saved_R2_global' in d or 'historical_global_physical_cost' in d:
        row['notes'].append('Historical R2 reference excluded from new-work charge.')
    if 'shared_geometry' in d:
        row['notes'].append(d['shared_geometry'])
    if 'shared_reference_receipt' in d:
        row['notes'].append('Shared port reference charged once in its independent setup row.')
    if path.parent.name == 'current_invariants':
        row['notes'].append('Saved-current diagnostic; explicitly zero new physical LU/RHS.')
        row['physical_LU'] = d.get('new_physical_LU')
        row['model_complex_RHS_excluding_locality'] = d.get('new_complex_RHS')
    return row


def build(output):
    rows, components = [], []
    # Only immediate development run receipts. Never traverse confirmation or truth.
    paths = sorted((ROOT/'results').glob('*/RECEIPT.json'))
    paths = [p for p in paths if not p.parent.name.startswith('confirmation')]
    paths += [ROOT/'results/port_profile_setup/REFERENCE.json',
              ROOT/'results/port_profile_setup/GRADIENT_WIRING.json',
              ROOT/'results/tmatrix_control/STARTUP_FAILURE.json']
    paths += sorted((ROOT/'results').glob('*/POSTPROCESS_RECEIPT.json'))
    paths += [ROOT/'audit'/name for name in (
        'COMPUTE_PROBE_STARTUP_FAILURE.json', 'COMPUTE_RESOURCE_PROBE.json',
        'DEVELOPMENT_COST_AUDIT.json', 'DEVELOPMENT_PLOT_RECEIPT.json',
        'DEVELOPMENT_STAGE_DIAGNOSTIC.json')]
    for path in paths:
        if not path.is_file():
            continue
        d = read(path)
        phase = ('validation' if 'GRADIENT' in path.name or path.parent.name in ('curvature','current_invariants')
                 else 'setup' if path.name == 'REFERENCE.json' or 'STARTUP' in path.name or 'PROBE' in path.name
                 else 'scoring' if path.parent.name == 'audit' or path.name == 'POSTPROCESS_RECEIPT.json'
                 else 'mixed')
        row = make_row(path, d, phase)
        rows.append(row)
        # All component views are NONADDITIVE, even when their counters overlap.
        for key, label in [('startup','setup'), ('bound_setup_operation_counts','setup'),
                           ('auxiliary_operation_counts','auxiliary'), ('auxiliary_cost','auxiliary'),
                           ('completion_cost','auxiliary'), ('corrector_cost','corrector'),
                           ('ray_cost','auxiliary'), ('gradient_check','validation')]:
            if isinstance(d.get(key), dict):
                components.append(dict(parent_receipt=row['receipt'], phase=label,
                    source_field=key, additive=False, values=d[key]))
        corrector = path.parent/'CORRECTOR.json'
        if path.name == 'RECEIPT.json' and corrector.is_file():
            c = read(corrector)
            components.append(dict(parent_receipt=row['receipt'], phase='corrector',
                source_field=str(corrector.relative_to(ROOT)), additive=False,
                values={k:v for k,v in c.items() if k in ('cost','wall_seconds','iterations','evaluations','stop')}))
    numerical = ('physical_LU','model_complex_RHS_excluding_locality','nonlocal_locality_LU',
                 'nonlocal_locality_RHS', *SMALL_KEYS)
    totals = {k:sum(row[k] for row in rows if isinstance(row.get(k), (int,float))) for k in numerical}
    missing = {k:[r['receipt'] for r in rows if r.get(k) is None] for k in ('physical_LU','model_complex_RHS_excluding_locality','wall_seconds','cpu_seconds')}
    payload = dict(status='DEVELOPMENT_RECEIPT_INVENTORY_CONFIRMATION_PENDING',
        created_utc=datetime.now(timezone.utc).isoformat(), rows=rows, components=components,
        recorded_known_subtotals=totals, missing_fields=missing,
        confirmation=dict(status='PENDING_PARENT_REBUILD', included=False),
        amount=None, currency=None, hash_verified=False, git_base='ordinary provenance only; not a content verification',
        time_policy='Per-invocation wall/CPU retained; no sum reported as parallel elapsed.',
        completeness='Lower-bound receipt inventory, not total project expenditure. Overwritten/unrecorded retries, generation, transport and historical operator assembly not inferred.')
    output.mkdir(parents=True, exist_ok=True)
    (output/'COST_LEDGER.json').write_text(json.dumps(payload,indent=2)+'\n')
    columns = [k for k in rows[0] if k not in ('full_counts','failure_details','notes')]
    with (output/'COST_LEDGER.csv').open('w', newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    print(json.dumps(dict(status=payload['status'], recorded_rows=len(rows), output=str(output), known_subtotals=totals)))


def build_full(output):
    """Extend saved development snapshot; never rewrite its historical files."""
    base_path = ROOT/'audit/cost_ledger/COST_LEDGER.json'
    base = read(base_path)
    rows = list(base['rows'])
    components = list(base['components'])
    seen = {r['receipt'] for r in rows}
    def append(path, d=None, phase='mixed'):
        key = str(path.relative_to(ROOT))
        if key in seen or not path.is_file():
            return
        d = read(path) if d is None else d
        row = make_row(path, d, phase)
        for k in ('caller_wall_seconds','caller_cpu_seconds','solve_wall_s','total_worker_s'):
            if k in d:
                row[k] = d[k]
        rows.append(row); seen.add(key)
        for count_key in ('auxiliary_operation_counts','bound_setup_counts','bound_setup_operation_counts','gradient_check'):
            if isinstance(d.get(count_key),dict):
                components.append(dict(parent_receipt=key,phase='auxiliary' if count_key=='auxiliary_operation_counts' else 'setup/validation',additive=False,source_field=count_key,values=d[count_key]))
        if isinstance(d.get('solver_import_charge'),dict):
            components.append(dict(parent_receipt=key,phase='setup',additive=False,
                source_field='solver_import_charge',values=d['solver_import_charge']))
        for stage in ('setup','initial','auxiliary','ranking','corrector','finalization'):
            count_key = stage+'_physical_model_cost'
            if count_key in d:
                components.append(dict(parent_receipt=key,phase=stage,additive=False,
                    source_field=count_key,values=d[count_key],
                    wall_seconds=d.get(stage+'_wall_seconds'),cpu_seconds=d.get(stage+'_cpu_seconds')))
    # Parent-approved public-result receipts only; never traverse sealed_truth.
    remote_dir = ROOT/'results/confirmation_remote/results/confirmation'
    for path in sorted(remote_dir.glob('obs_*.json')):
        outer = read(path)
        d = dict(outer.get('result', {}))
        d['run_id'] = outer.get('task_id', path.stem)
        d['worker_status'] = outer.get('status')
        d['solve_wall_s'] = outer.get('solve_wall_s')
        d['total_worker_s'] = outer.get('total_worker_s')
        append(path, d)
    for rel in ('audit/confirmation/GENERATION_RECEIPT.json',
                'grid32_diagnostic/audit/generation/GENERATION_RECEIPT.json'):
        path=ROOT/rel
        if path.is_file():
            d=read(path); d['physical_model_cost']=d.get('physical_counts', {})
            append(path,d,'setup')
    for path in sorted((ROOT/'grid32_diagnostic/results').glob('*/RECEIPT.json')):
        append(path)
    # Any development top-level arm created since the preserved snapshot (e.g. pi).
    for path in sorted((ROOT/'results').glob('*/RECEIPT.json')):
        if not path.parent.name.startswith('confirmation'):
            append(path)
    for pattern in ('partition_average/*/RECEIPT.json','failed_endpoint_continuation/*/RECEIPT.json','cross_platform_remote/development_windows/RECEIPT.json'):
        for path in sorted((ROOT/'results').glob(pattern)):
            append(path)
    score_paths=sorted((ROOT/'results').glob('*SCORES*.json')) + [ROOT/'results/DEVELOPMENT_POST_CONFIRMATION_ADDENDA_V1.json']
    for path in score_paths:
        # Only recorded cost fields, not material scores, are exported.
        d=read(path)
        keep={k:v for k,v in d.items() if k in ('status','wall_seconds','cpu_seconds','timer_scope','new_physical_LU','new_complex_RHS','hash_checks')}
        if 'scoring_new_physical_LU' in d:
            keep['new_physical_LU']=d['scoring_new_physical_LU']
        if 'scoring_new_complex_RHS' in d:
            keep['new_complex_RHS']=d['scoring_new_complex_RHS']
        append(path,keep,'scoring')
    for rel in ('audit/ENVIRONMENT.json','audit/confirmation/SCORING_PREPARATION_V1.json',
                'audit/VALIDATION_PLOT_RECEIPT.json','audit/REFINEMENT_PLOT_RECEIPT.json',
                'audit/confirmation/CONFIRMATION_COMPARISON_AUDIT.json',
                'audit/confirmation/FAILURE_STOP_DIAGNOSTIC.json'):
        append(ROOT/rel,phase='validation')
    transport=[]
    for name in ('PRELAUNCH_ACTIVITY.json','UPLOAD_RECEIPT.json','REMOTE_PREPARATION_RECEIPT.json',
                 'REMOTE_RUN_RECEIPT.json','REMOTE_RESULT_ARCHIVE_RECEIPT.json','DOWNLOAD_RECEIPT.json',
                 'WINDOWS_ZIP_PATH_NORMALIZATION.json','WINDOWS_VALIDATION_ZIP_PATH_NORMALIZATION.json',
                 'WINDOWS_VALIDATION_RUN_STARTUP_FAILURE.json',
                 'WINDOWS_VALIDATION_RUN.json','WINDOWS_VALIDATION_TRANSPORT.json'):
        path=ROOT/'audit/confirmation'/name
        if not path.is_file():continue
        d=read(path)
        entry=dict(receipt=str(path.relative_to(ROOT)),phase='transport/orchestration',
                   wall_seconds=d.get('wall_seconds'),returncode=d.get('returncode'),
                   amount=None,currency=None,add_to_physical_totals=False)
        if name.startswith('WINDOWS_VALIDATION_'):
            execution=d.get('execution',{})
            entry.update(status=('FAILED_IMPORT_PAID_ATTEMPT' if 'STARTUP_FAILURE' in name else 'TRANSPORT_OR_EXECUTION_RECORDED'),
                         upload=d.get('upload'),execution=execution,details=d,
                         physical_LU=(0 if 'STARTUP_FAILURE' in name else None),
                         complex_RHS=(0 if 'STARTUP_FAILURE' in name else None))
        elif name.endswith('ZIP_PATH_NORMALIZATION.json'):
            entry.update(reason=d.get('reason'),entries=len(d.get('entries',[])),
                         original_zip_preserved=d.get('original_zip_preserved'),physical_LU=d.get('physical_LU'))
        else:
            entry.update(stderr=d.get('stderr'),scope=d.get('scope'),
                         actual_cuda_workers=d.get('actual_cuda_workers'))
        transport.append(entry)
    numerical=('physical_LU','model_complex_RHS_excluding_locality','nonlocal_locality_LU',
               'nonlocal_locality_RHS',*SMALL_KEYS)
    totals={k:sum(r[k] for r in rows if isinstance(r.get(k),(int,float))) for k in numerical}
    dimensions={}
    for r in rows:
        lu=r.get('physical_LU'); dims=r.get('physical_factor_dimension_sum')
        dimension = str(int(dims/lu)) if lu and dims is not None and dims/lu==int(dims/lu) else 'unknown_or_no_LU'
        bucket=dimensions.setdefault(dimension,dict(physical_LU=0,model_complex_RHS_excluding_locality=0))
        for k in bucket:
            if isinstance(r.get(k),(int,float)):bucket[k]+=r[k]
    run=read(ROOT/'audit/confirmation/REMOTE_RUN_RECEIPT.json')
    payload=dict(status='FULL_V2_RECORDED_SCIENTIFIC_WORK_COMPLETE',created_utc=datetime.now(timezone.utc).isoformat(),
        development_base=str(base_path.relative_to(ROOT)),rows=rows,components=components,
        recorded_known_subtotals=totals,physical_work_by_dimension=dimensions,
        missing_fields={k:[r['receipt'] for r in rows if r.get(k) is None] for k in ('physical_LU','model_complex_RHS_excluding_locality','wall_seconds','cpu_seconds')},
        confirmation=dict(downloaded_task_receipts=len(list(remote_dir.glob('obs_*.json'))),
            remote_run_elapsed_seconds=run['wall_seconds'],actual_cuda_workers=run.get('actual_cuda_workers'),
            timing_scope=run.get('scope'),worker_wall_sum_is_parallel_elapsed=False),
        transport_and_orchestration=transport,
        cross_platform_windows=dict(status='COMPLETED_RECORDED',receipt='results/cross_platform_remote/development_windows/RECEIPT.json',physical_LU=1032,complex_RHS=24684,solver_wall_seconds=376.1966163000325,orchestration_timing_source='audit/confirmation/WINDOWS_VALIDATION_RUN.json',separate_from_confirmation_elapsed=True),
        amount=None,currency=None,llm_cost=None,hash_verified=False,
        total_project_parallel_elapsed_seconds=None,
        time_policy='Remote run elapsed includes queue startup/all tasks/shutdown. Transport and worker wall are separate views and are never added to it. No full-project concurrency timer exists.',
        limits='Known receipt subtotal only. Failed/refused attempts retained where recorded; missing costs not assumed zero. No confirmation GT or private NPZ read. Ledger rebuild runs no science.')
    output.mkdir(parents=True,exist_ok=True)
    (output/'COST_LEDGER_FULL_V2.json').write_text(json.dumps(payload,indent=2)+'\n')
    columns=sorted({k for row in rows for k in row if k not in ('full_counts','failure_details','notes')})
    with (output/'COST_LEDGER_FULL_V2.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    print(json.dumps(dict(status=payload['status'],rows=len(rows),subtotals=totals,by_dimension=dimensions,output=str(output))))


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT/'audit/cost_ledger')
    parser.add_argument('--full',action='store_true',help='Extend saved development snapshot; output separate FULL_V2 final inventory.')
    args=parser.parse_args()
    (build_full if args.full else build)(args.output.resolve())
