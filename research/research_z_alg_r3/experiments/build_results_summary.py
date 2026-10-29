"""Aggregate final V6/V7 descriptive evidence and all-paid costs."""
import datetime
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    v6=json.loads((ROOT/'results'/'HOLDOUT_V6_SUMMARY.json').read_text())
    v7=json.loads((ROOT/'results'/'CONFIRM_V7_SUMMARY.json').read_text())
    costs=json.loads((ROOT/'COST_LEDGER.json').read_text())
    freeze=json.loads((ROOT/'audit'/'EXPERIMENT_SCIENCE_FREEZE.json').read_text())
    result=dict(schema='R3_FINAL_RESULTS_SUMMARY_V1',built_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                scientific_disposition='HOLD independent-paper novelty/support; finite mechanism and one additional fresh-case result only',
                science_freeze=freeze,
                model=dict(type='Finite nonlinear scalar-TM disk-pulse VIE',material_unknowns=12,sources=3,receivers=3,complex_observations=9,
                           public_material_box=[.5,30.],all_12_materials_unknown=True,geometry_known=True,noise='none',scope_excludes=['continuous-PDE certification','3D Maxwell','Full66','global root completeness','global optimum','NN/GPU evidence','scale-study evidence']),
                metric_contract=dict(online_selection='Minimum original all-observation relative residual over paid candidate archive',
                                     predeclared_data_fit_relative_residual=1e-6,strict_optimizer_stop_goal=1e-10,
                                     material_error='RMS over all12 material values; scored only after prediction freeze',
                                     predeclared_GT_material_RMS_success_gate=None,
                                     offline_oracle='best@K/minimum material RMS among returned candidates; uses GT after freeze, never an online selector'),
                DEV=dict(targets=2,role='Exploratory algorithm development, not confirmatory evidence',
                         paid_costs=costs['totals']['dev_solver'],
                         successful_pivot_receipt='results/dev_pivot_v6_scores.json',
                         private_material_values_or_target_seeds_supplied_to_solver=False),
                V6=dict(targets=4,source_freeze='audit/HOLDOUT_SOURCE_FREEZE.json',prediction_freeze='audit/HOLDOUT_V6_FINAL_FREEZE.json',
                        role='Initially fresh V6 holdouts; exposed outcomes used as development for V7, never counted as fresh V7 evidence',
                        method_summaries=v6['method_summaries'],metrics_csv='results/HOLDOUT_V6_METRICS.csv',
                        detailed_summary='results/HOLDOUT_V6_SUMMARY.json',
                        negative_evidence='GlobalDual and reciprocal FullJ each fit2/4; GlobalDual loses case01 that every other method solves; no dominance',
                        independent_validation=v6['independent_all_candidate_validation'],
                        paid_costs=costs['totals']['holdout_solver'],
                        figures=['figures/holdout_materials.png','figures/holdout_cost_coverage.png']),
                V7=dict(targets=4,source_freeze='audit/CONFIRM_V7_SOURCE_FREEZE.json',prediction_freeze='audit/CONFIRM_V7_FINAL_FREEZE.json',
                        role='Four new NEWCONFIRM12 targets generated after the V7 source/config freeze',
                        method_summaries=v7['method_summaries'],metrics_csv='results/CONFIRM_V7_METRICS.csv',
                        detailed_summary='results/CONFIRM_V7_SUMMARY.json',
                        primary_finite_result='GuardedGlobalDual fits one additional fresh case at declared residual1e-6:3/4 versus2/4; both meet strict stop1e-10 on2/4',
                        negative_evidence='Both methods fail case01; no predeclared GT accuracy gate or statistical/global/general-performance conclusion; all4 invoke continuation',
                        paired_guard_identical=all(x['identical_numerical_guard'] for x in v7['common_guard_checks']),
                        phase_accounting_exact=all(x['phase_counts_sum_exact'] for x in v7['phase_accounting_checks']),
                        independent_validation=v7['independent_validation'],
                        public_replay=v7['public_replay'],paid_costs=costs['totals']['confirm_v7_solver'],
                        figures=['figures/confirm_v7_materials.png','figures/confirm_v7_cost_coverage.png']),
                all_paid_costs=dict(receipts='COST_LEDGER.json',csv='COST_LEDGER.csv',stage_totals=costs['totals'],
                                    grand_full_physics_complex_rhs=costs['grand_actual_complex_rhs'],
                                    grand_reduced_complex_rhs=costs['grand_reduced_complex_rhs'],
                                    grand_algebraic_complex_LS_rhs=costs['grand_algebraic_complex_LS_rhs'],
                                    end_to_end_process_wall=costs['end_to_end_process_wall'],
                                    failed_attempts=costs['failed_attempts'],missing_costs=costs['missing_costs'],
                                    interpretation='Actual complex columns for12-dimensional full systems; reduced and inner LS RHS separate; process wall overlaps numerical wall and is not added'),
                disclosure=dict(material_only_truth_root='truth_after_freeze',material_files=10,
                                private_truth_retained_local=True,target_generation_seeds_disclosed=False,
                                public_optimizer_seeds='Public deterministic solver config, distinct from hidden target-generation seeds',
                                replay_is_new_blind_evidence=False),
                post_freeze_physical_experiments='NOT_RUN',
                technical_notes='EXPERIMENT_TECHNICAL_NOTES.md')
    (ROOT/'RESULTS_SUMMARY.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output='RESULTS_SUMMARY.json',grand_full_rhs=costs['grand_actual_complex_rhs'],disposition='HOLD')))


if __name__=='__main__':main()
