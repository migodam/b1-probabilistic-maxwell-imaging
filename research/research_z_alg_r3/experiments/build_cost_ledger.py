"""Aggregate every original paid receipt, avoiding frozen snapshot duplication."""
import datetime
import csv
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    runs=[]
    for path in sorted((ROOT/'results').glob('*.json')):
        receipt=json.loads(path.read_text())
        if not isinstance(receipt,dict) or 'accounting' not in receipt or 'candidates' not in receipt:
            continue
        accounting=receipt['accounting']; counts=accounting.get('counts',{})
        fits=receipt.get('fits',[])
        split='confirm_v7' if path.name.startswith('confirm_') else ('holdout' if path.name.startswith('holdout') else 'dev')
        runs.append(dict(file=str(path.relative_to(ROOT)),split=split,
                         method=receipt['method'],revision=receipt.get('solver_revision'),
                         status=receipt['status'],stop_reason=receipt.get('stop_reason'),
                         wall_seconds=receipt['wall_seconds'],requested_wall_seconds=receipt['config'].get('wall_seconds'),
                         full_dimension=receipt.get('material_dimension'),
                         full_complex_rhs=accounting['full_complex_rhs'],
                         reduced_complex_rhs=accounting['reduced_complex_rhs'],
                         full_factorizations=accounting['full_factorizations'],
                         reduced_factorizations=accounting['reduced_factorizations'],
                         algebraic_complex_rhs=counts.get('algebraic_complex_rhs',0),
                         algebraic_QR_factorizations=counts.get('algebraic_QR',0),
                         quotient_linear_solves=counts.get('quotient_linear_solves',0),
                         quotient_recovery_least_squares=counts.get('quotient_recovery_least_squares',0),
                         candidate_count=len(receipt['candidates']),
                         full_least_squares_fit_invocations=len(fits),
                         full_fit_status_counts=dict(Counter(str(f.get('status')) for f in fits)),
                         full_fit_nfev_sum=sum(f.get('nfev',0) or 0 for f in fits),
                         full_fit_njev_sum=sum(f.get('njev',0) or 0 for f in fits),
                         operation_counts=counts,operation_seconds=accounting.get('operation_seconds',{}),
                         rejections=receipt.get('rejections',{})))
    generation=[]
    for path in sorted((ROOT/'audit').glob('*GENERATION_RECEIPT.json')):
        r=json.loads(path.read_text())
        generation.append(dict(file=str(path.relative_to(ROOT)),split=r['split'],case_count=len(r['cases']),
                               full_complex_rhs=r['generation_counts'].get('actual_complex_RHS',0),
                               full_factorizations=r['generation_counts'].get('LU',0),
                               counts=r['generation_counts'],operation_seconds=r['generation_timings'],
                               wall_seconds=sum(r['generation_timings'].values()),
                               measured_operation_seconds=sum(r['generation_timings'].values())))
    validation=[]
    scoring=[]
    for path in sorted((ROOT/'results').glob('*.json')):
        r=json.loads(path.read_text())
        if isinstance(r,dict) and 'validated_utc' in r:
            validation.append(dict(file=str(path.relative_to(ROOT)),full_complex_rhs=r['full_complex_rhs'],
                                   full_factorizations=r['full_factorizations'],wall_seconds=r['wall_seconds']))
        if isinstance(r,dict) and r.get('scoring_only'):
            scoring.append(dict(file=str(path.relative_to(ROOT)),full_complex_rhs=r['new_physics_rhs'],
                                wall_seconds=r.get('scoring_wall_seconds'),
                                legacy_wall_not_recorded='scoring_wall_seconds' not in r))
    def totals(rows):
        fields=('full_complex_rhs','reduced_complex_rhs','algebraic_complex_rhs','full_factorizations','reduced_factorizations','algebraic_QR_factorizations','quotient_linear_solves','quotient_recovery_least_squares','wall_seconds')
        return {field:sum(r.get(field,0) or 0 for r in rows) for field in fields}
    failed=[dict(stage='DEV generation attempt 0',status='FAILED_BEFORE_GENERATION',
                 full_complex_rhs=0,full_factorizations=0,wall_seconds=None,
                 evidence='audit/HISTORY_AUDIT_CLOSURE.json',
                 reason='Relative contract path mismatch before material generation; time not retained')]
    process_wall=[]
    for filename,split in [('HOLDOUT_EXECUTION_RECEIPT.json','holdout_v6'),('CONFIRM_V7_EXECUTION_RECEIPT.json','confirm_v7')]:
        execution_path=ROOT/'audit'/filename
        if not execution_path.exists():continue
        execution=json.loads(execution_path.read_text())
        process_wall.append(dict(split=split,subprocess_wall_seconds=sum(x['process_wall_seconds'] for x in execution['rows']),
                                 batch_wall_seconds=execution.get('batch_wall_seconds'),
                                 boundary='Each complete subprocess includes interpreter/imports, numerical solve and result file IO; overlaps numerical solver wall, never add both',
                                 rows=execution['rows']))
    result=dict(schema='R3_ALL_PAID_COST_LEDGER_V1',
                built_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                units='Actual complex RHS columns; all full systems dimension 12; reduced RHS separate, not full equivalent',
                solver_wall_boundary='Starts before cached public operator import; includes full fits, quotient setup/enumeration/filtering, archive evaluation and corrections; excludes Python process/import and file IO',
                setup_boundary='Green setup + receiver SVD charged once in generation; cached operators shared by all methods',
                missing_costs=['Python startup/import/file IO not individually timed','Initial failed generation elapsed time not retained','Legacy DEV scoring elapsed times not retained'],
                failed_attempts=failed,generation=generation,solver_runs=runs,
                independent_validation=validation,scoring=scoring,
                end_to_end_process_wall=process_wall,
                totals=dict(dev_solver=totals([r for r in runs if r['split']=='dev']),
                            holdout_solver=totals([r for r in runs if r['split']=='holdout']),
                            confirm_v7_solver=totals([r for r in runs if r['split']=='confirm_v7']),
                            all_solver=totals(runs),generation=totals(generation),
                            independent_validation=totals(validation),scoring=totals(scoring)),
                duplication_boundary='Each original results/*.json solver receipt counted once; immutable frozen copies are not additional runs',
                no_gt_read=True,no_historical_R2_reruns=True)
    result['grand_actual_complex_rhs']=result['totals']['all_solver']['full_complex_rhs']+result['totals']['generation']['full_complex_rhs']+result['totals']['independent_validation']['full_complex_rhs']
    result['grand_reduced_complex_rhs']=result['totals']['all_solver']['reduced_complex_rhs']
    result['grand_algebraic_complex_LS_rhs']=result['totals']['all_solver']['algebraic_complex_rhs']
    (ROOT/'COST_LEDGER.json').write_text(json.dumps(result,indent=2)+'\n')
    process_lookup={row['output']:row['process_wall_seconds'] for batch in process_wall for row in batch['rows']}
    flat=[]
    fields=['category','split','method','file','status','full_dimension','full_complex_rhs','reduced_complex_rhs','algebraic_complex_rhs','full_factorizations','reduced_factorizations','algebraic_QR_factorizations','quotient_linear_solves','quotient_recovery_least_squares','wall_seconds','end_to_end_process_wall_seconds','candidate_count']
    for category,items in [('solver',runs),('generation',generation),('independent_validation',validation),('scoring',scoring),('failed_attempt',failed)]:
        for item in items:
            row={field:item.get(field) for field in fields};row['category']=category
            row['end_to_end_process_wall_seconds']=process_lookup.get(item.get('file'))
            flat.append(row)
    with (ROOT/'COST_LEDGER.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(flat)
    lines=['# R3 cost ledger','','All RHS below are actual complex columns for 12-dimensional full systems. Reduced and algebraic least-squares RHS are reported separately, never converted to full RHS equivalents. Every original DEV and holdout receipt is charged once; frozen copies are not new runs.','','| Stage | Full RHS | Reduced RHS | Algebraic LS RHS | Full factors | Reduced factors | Inner QR | solver/validation wall s |','|---|---:|---:|---:|---:|---:|---:|---:|']
    for stage,values in result['totals'].items():
        lines.append(f"| {stage} | {values['full_complex_rhs']} | {values['reduced_complex_rhs']} | {values['algebraic_complex_rhs']} | {values['full_factorizations']} | {values['reduced_factorizations']} | {values['algebraic_QR_factorizations']} | {values['wall_seconds']:.6f} |")
    for process in process_wall:lines.extend(['',f"{process['split']} end-to-end subprocess time: {process['subprocess_wall_seconds']:.6f} s including interpreter/imports and result IO. This contains the numerical wall above and is not added to it."])
    lines.extend(['','Generation operation time is separately retained in JSON. Shared geometry setup is charged there. Solver time includes cached operator import and all numerical work; Python startup/import and file IO are excluded. The failed initial generation consumed zero RHS; its elapsed time and legacy DEV scoring times were not retained.','','The JSON provides per-run method/revision, termination status, all operation counts/times, fit invocations, candidate counts and rejection counts. No R2 experiment is rerun or charged as new work.'])
    lines.extend(['',f"Grand paid physical full RHS including all DEV/V6/V7 solver work, target generation and independent validation: {result['grand_actual_complex_rhs']}. Separately, reduced RHS: {result['grand_reduced_complex_rhs']}; algebraic complex least-squares RHS: {result['grand_algebraic_complex_LS_rhs']}. The stage subtotals and all-solver total are overlapping views and must not be added together."])
    (ROOT/'COST_LEDGER.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(result['totals']))


if __name__=='__main__':main()
