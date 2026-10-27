"""Read saved G0 evidence and record its numerical scope, not scientific hype."""
from pathlib import Path
import csv,json
import numpy as np
from bgga.accounting import BudgetBook,_default
ROOT=Path(__file__).resolve().parents[1]
def main():
    qual=json.loads((ROOT/'results/g0_qualification/RESULTS.json').read_text())
    ops=json.loads((ROOT/'results/g0_operators/RESULTS.json').read_text())
    stress=json.loads((ROOT/'results/g0_stress/RESULTS.json').read_text())
    assert len(qual)==12 and len(ops)==35 and len(stress)==12
    exact=[r for r in ops if r['ordinary_exact_same_discrete']]
    numerical=all(r['numerical_qualification'] for r in qual) and all(r['strict_algebra_status']=='PASS' for r in exact)
    receipts=[json.loads(p.read_text()) for p in (ROOT/'results/jobs').glob('*/RECEIPT.json')]
    assert all(json.loads((ROOT/'results/g0_stress/COMPLETE.json').read_text())[k]==v for k,v in [('fixtures',12),('unique_physical_parents',9)])
    table=[]
    for r in stress:
        d=r['diagnostics']['pose_intrinsic_instrument']['data_only']
        p=r['diagnostics']['pose_only']['data_only']
        row=dict(case=r['case'],family=r['family'],physical_parent=r['physical_parent'],
            empirical_operator_error=r['empirical_operator_error'],rank_status=d['status'],
            wall_rank=d['wall_rank_audit']['rank'],nuisance_rank=d['nuisance_rank_audit']['rank'],
            retained_min=min(d['generalized_retention_eigenvalues']) if d['generalized_retention_eigenvalues'] else None,
            retained_max=max(d['generalized_retention_eigenvalues']) if d['generalized_retention_eigenvalues'] else None,
            pose_only_rank_status=p['status'],
            pose_only_retained_min=min(p['generalized_retention_eigenvalues']) if p['generalized_retention_eigenvalues'] else None,
            pose_only_retained_max=max(p['generalized_retention_eigenvalues']) if p['generalized_retention_eigenvalues'] else None,
            absolute_info_min=min(d['wall_absolute_eigenvalues']),schur_info_min=min(d['wall_data_only_schur_eigenvalues']),
            finite_wall_variance_fraction=r['finite_posterior']['wall_variance_fraction_of_prior'],
            joint_pair_distance=r['finite_pair']['current_norm'],same_fixed_j0_task_distance=r['finite_pair']['task_norm_value'],
            future25_max=r.get('future25',{}).get('max'),
            future25_independent=r.get('future25',{}).get('independent_argmax_distance'),
            future25_empirical_margin=r.get('future25',{}).get('empirical_margin_lower'))
        table.append(row)
    with (ROOT/'results/G0_STRESS_TABLE.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(table[0]));wr.writeheader();wr.writerows(table)
    timings={m:[r['solve_timing']['measured_kernel_cpu_seconds'] for r in ops if r['method']==m]
        for m in ['full','woodbury','schur_from_paid_ports','krylov','goal_rom']}
    result=dict(status='QUALIFIED_DISCRETE_FULLWAVE_FOR_TESTED_CONTRACT' if numerical else 'FAIL_NUMERICAL',
        qualification_count=12,qualification_whitened_min=min(r['fixed_vs_interface_fitted_white'] for r in qual),
        qualification_whitened_max=max(r['fixed_vs_interface_fitted_white'] for r in qual),
        exact_rows=len(exact),maximum_exact_field_relative=max(r['field_difference_to_full']['relative'] for r in exact),
        maximum_exact_response_relative=max(r['response_difference_to_full']['relative'] for r in exact),
        maximum_exact_backward=max(r['full_equation_backward']['backward'] for r in exact),
        operator_query_mean_cpu={m:float(np.mean(v)) for m,v in timings.items()},
        timings_scope='measured query kernels, common process; not independently measured cold deployment',
        goal_ROM_heldout_whitened_errors=[r['full_response_difference_whitened_l2'] for r in ops if r['method']=='goal_rom' and r['rom_partition']=='HELD_OUT'],
        stress_table=table,stress_fixtures=12,independent_physical_worlds=9,
        C2_C3_identity='matched measurement contracts on the same three physical parents',
        gaussian_scope='proper tangent Gaussian reference only, distinct from boxed finite15 prior',
        finite15_scope='conditional known-support mechanism fixtures; no unknown-support A0 claim',
        qualified_continuum=False,uniform_prior_forward_bound=False,
        new_BGGA_cold_advantage=False,
        G1_allowed='finite-prior logic/accuracy replay only; cannot turn replay counts into fresh PDE savings',
        G2_status='NOT_RUN_NO_NET_SAVING_MECHANISM',
        physical_usage=BudgetBook().snapshot('G0'),
        missing_metadata=[dict(case='C1_0',fields=['original exact duality scalar','three-trial corrector history'],reason='pre-serialization software failure; raw arrays reused and charges retained')])
    (ROOT/'results/G0_DECISION.json').write_text(json.dumps(result,default=_default,ensure_ascii=False,indent=2)+'\n')
    gate=json.loads((ROOT/'GATE_DECISION.json').read_text())
    gate.update(overall='G0_NUMERICAL_QUALIFIED_G1_LOGIC_NEXT',G0=result['status'],
        EXACT_DISCRETE_ASSIMILATION='QUALIFIED_ON_TESTED_FIXED_SPACE',
        TARGET_INVARIANT_PROPAGATION_CACHE='QUALIFIED_DISCRETE_TESTED_CONFIGURATIONS; not uniform continuum',
        ANGLE_ONLY_SAFETY='inspect reported spectra/finite15 uncertainty; no safety certificate from retention alone',
        COLD_COST='NO_GO_CURRENT_FULL_PORT_IMPLEMENTATION; no independently measured deployment advantage',
        resource_gate='human grant recorded; stage caps remain independent')
    (ROOT/'GATE_DECISION.json').write_text(json.dumps(gate,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],query_CPU=result['operator_query_mean_cpu'],
        usage=result['physical_usage'],stress_table=table),default=_default,ensure_ascii=False))
if __name__=='__main__':main()
