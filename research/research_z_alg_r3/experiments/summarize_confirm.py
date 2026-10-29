"""Summarize the final frozen fresh V7 trial without inventing a GT gate."""
import csv
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def main():
    payload=json.loads((ROOT/'results'/'confirm_v7_scores.json').read_text())
    valid=json.loads((ROOT/'results'/'confirm_v7_validation.json').read_text())
    replay=json.loads((ROOT/'results'/'confirm_v7_public_replay_scores.json').read_text())
    if payload['rows']!=replay['rows']:raise RuntimeError('Public m-only replay differs from frozen scoring')
    rows=[];predictions={};phase_checks=[]
    for scored in payload['rows']:
        result=json.loads((ROOT/scored['file']).read_text())
        predictions[(scored['case_id'],result['method'])]=result
        phases=result['phase_receipts']
        cost=result['accounting'];countsum={}
        for phase in phases.values():
            for key,value in phase['accounting']['counts'].items():countsum[key]=countsum.get(key,0)+value
        exact=countsum==cost['counts'] and all(sum(p['accounting'][key] for p in phases.values())==cost[key] for key in ('full_complex_rhs','full_factorizations','reduced_complex_rhs','reduced_factorizations'))
        phase_checks.append(dict(file=scored['file'],phase_counts_sum_exact=exact))
        guard=phases['guard'];continuation=phases.get('continuation')
        rows.append(dict(case_id=scored['case_id'],method=result['method'],status=result['status'],stop_reason=result['stop_reason'],
                         top1_RMS=scored['top1_RMS'],top1_relative_residual=scored['top1_relative_residual'],
                         data_fit=bool(scored['top1_relative_residual']<=1e-6),strict_search_goal=bool(scored['top1_relative_residual']<=1e-10),
                         oracle_best_RMS=scored['oracle_best_RMS'],candidate_count=scored['candidate_count'],
                         data_accepted_candidate_count=scored['data_accepted_count'],
                         outer_wall_seconds=result['wall_seconds'],reported_overrun_seconds=result['outer_budget']['reported_overrun_seconds'],
                         top1_discovery_seconds=result['top1']['discovery_seconds'],full_complex_rhs=cost['full_complex_rhs'],
                         inner_QR=cost['counts'].get('algebraic_QR',0),inner_complex_LS_RHS=cost['counts'].get('algebraic_complex_rhs',0),
                         guard_relative_residual=guard['top1']['relative_residual'],guard_core_wall_seconds=guard['wall_seconds'],
                         guard_full_complex_rhs=guard['accounting']['full_complex_rhs'],
                         continuation_core_wall_seconds=continuation['wall_seconds'] if continuation else 0,
                         continuation_full_complex_rhs=continuation['accounting']['full_complex_rhs'] if continuation else 0))
    summaries=[]
    for method in sorted(set(r['method'] for r in rows)):
        subset=[r for r in rows if r['method']==method]
        summaries.append(dict(method=method,cases=len(subset),top1_data_fit_cases=sum(r['data_fit'] for r in subset),
                              strict_search_goal_cases=sum(r['strict_search_goal'] for r in subset),
                              mean_top1_RMS=float(np.mean([r['top1_RMS'] for r in subset])),
                              median_top1_RMS=float(np.median([r['top1_RMS'] for r in subset])),
                              mean_oracle_best_RMS=float(np.mean([r['oracle_best_RMS'] for r in subset])),
                              mean_outer_wall_seconds=float(np.mean([r['outer_wall_seconds'] for r in subset])),
                              maximum_outer_overrun_seconds=max(r['reported_overrun_seconds'] for r in subset),
                              full_complex_rhs=sum(r['full_complex_rhs'] for r in subset),
                              inner_QR=sum(r['inner_QR'] for r in subset),inner_complex_LS_RHS=sum(r['inner_complex_LS_RHS'] for r in subset)))
    common_guard=[]
    for case in sorted(set(r['case_id'] for r in rows)):
        a=predictions[(case,'guarded_globaldual')]['phase_receipts']['guard']
        b=predictions[(case,'guarded_reciprocal')]['phase_receipts']['guard']
        common_guard.append(dict(case_id=case,identical_numerical_guard=bool(a['top1']['m']==b['top1']['m'] and a['top1']['relative_residual']==b['top1']['relative_residual'] and a['accounting']['counts']==b['accounting']['counts']),
                                 guard_top1_relative_residual=a['top1']['relative_residual'],guard_stopped_on_goal=a['top1']['relative_residual']<=1e-10))
    result=dict(schema='R3_FROZEN_CONFIRM_V7_SUMMARY',contract='experiments/CONFIRM_V7_CONTRACT.json',prediction_freeze='audit/CONFIRM_V7_FINAL_FREEZE.json',
                rows=rows,method_summaries=summaries,phase_accounting_checks=phase_checks,common_guard_checks=common_guard,
                independent_validation=dict(candidates=sum(r['validated_count'] for r in valid['rows']),all_pass=all(r['all_validated_pass'] for r in valid['rows']),full_complex_rhs=valid['full_complex_rhs'],
                                            maximum_stored_residual_disagreement=max(c.get('stored_residual_abs_disagreement',0) for r in valid['rows'] for c in r['rows'])),
                GT_material_RMS_success_threshold=None,
                threshold_boundary='Data-fit threshold1e-6 was declared before new-confirm generation. Strict optimizer stopping goal1e-10 is separate. No GT RMS success gate was predeclared.',
                public_replay=dict(receipt='results/confirm_v7_public_replay_scores.json',all_rows_identical=True,new_physics_rhs=0,material_only_truth_root='truth_after_freeze',target_seeds_disclosed=False),
                scientific_disposition='One additional fresh case at the declared data-fit threshold (3/4 versus2/4). Four finite targets and mixed earlier V6 evidence; independent-paper novelty/support remains HOLD.',
                expansion='NOT_RUN: no N24, NN, GPU, Full66 or further physical experiments after science freeze')
    (ROOT/'results'/'CONFIRM_V7_SUMMARY.json').write_text(json.dumps(result,indent=2)+'\n')
    with (ROOT/'results'/'CONFIRM_V7_METRICS.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    lines=['# Final fresh V7 confirm result','','The separately frozen guarded comparison fits one additional fresh target at the declared all-data residual threshold: GuardedGlobalDual 3/4 versus GuardedReciprocal 2/4. The sample is four targets. The independent-paper scientific disposition remains HOLD.','','The material RMS values below are descriptive. No material-RMS success gate was predeclared; none is introduced retrospectively. A data fit means relative residual ≤1e-6. The stricter optimizer stopping goal is 1e-10; both methods meet that goal on 2/4 cases.','','| Case | GuardedGlobalDual top1 RMS | GuardedGlobalDual residual | GuardedReciprocal top1 RMS | GuardedReciprocal residual |','|---|---:|---:|---:|---:|']
    for case in sorted(set(r['case_id'] for r in rows)):
        a=next(r for r in rows if r['case_id']==case and r['method']=='guarded_globaldual');b=next(r for r in rows if r['case_id']==case and r['method']=='guarded_reciprocal')
        lines.append(f"| {case} | {a['top1_RMS']:.6g} | {a['top1_relative_residual']:.6g} | {b['top1_RMS']:.6g} | {b['top1_relative_residual']:.6g} |")
    lines.extend(['','Case 00 supplies the additional data fit. Both methods fail on 01. Both fit 02 and 03; GuardedGlobalDual case 02 reaches residual 3.75e-10, which passes the declared data-fit criterion but misses its stricter search stop goal. None of these four fresh cases stops at the common native guard; all invoke a continuation. The paired guard materials, residuals and operation counts are identical, and every root cost ledger equals the exact sum of its phases.','','| Method | Data-fit cases | Strict search-goal cases | Mean top1 RMS | Median top1 RMS | Mean outer s | Full RHS | Inner QR | Inner LS RHS |','|---|---:|---:|---:|---:|---:|---:|---:|---:|'])
    for s in summaries:lines.append(f"| {s['method']} | {s['top1_data_fit_cases']}/4 | {s['strict_search_goal_cases']}/4 | {s['mean_top1_RMS']:.6g} | {s['median_top1_RMS']:.6g} | {s['mean_outer_wall_seconds']:.6g} | {s['full_complex_rhs']} | {s['inner_QR']} | {s['inner_complex_LS_RHS']} |")
    lines.extend(['','The total outer budget is 3 s including guard/setup/continuation/wrapper accounting. Actual terminal/bookkeeping overrun is reported: maximum 0.0128075 s. All 968 returned candidates pass independent full-VIE validation, charged 2,904 full RHS. Online top1 remains residual-selected, while best@K is labeled as offline oracle coverage.','','V6 predictions/scores remain unchanged and serve as development evidence for this V7 revision. The four NEWCONFIRM12 targets were generated after the V7 source/config freeze. Their final eight predictions froze before validation or GT scoring. Target material arrays were disclosed only after scoring in `truth_after_freeze/`; random target-generation seeds stay local. Public disclosed-truth scoring reproduces every row exactly and is explicitly not a new blind trial.','','See `CONFIRM_V7_METRICS.csv`, `CONFIRM_V7_SUMMARY.json`, `../COST_LEDGER.json`, and `../EXPERIMENT_TECHNICAL_NOTES.md` for the exact phase costs, original receipts and replay boundary. No physical experiments continue after `audit/EXPERIMENT_SCIENCE_FREEZE.json`.'])
    (ROOT/'results'/'CONFIRM_V7_SUMMARY.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(methods=summaries,phase_costs_exact=all(x['phase_counts_sum_exact'] for x in phase_checks),common_guard_identical=all(x['identical_numerical_guard'] for x in common_guard))))


if __name__=='__main__':main()
