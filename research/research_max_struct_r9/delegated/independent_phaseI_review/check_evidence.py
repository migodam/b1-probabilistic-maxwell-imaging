"""Independent completed-PhaseI metadata audit; no arrays/GT/physics imported."""
from pathlib import Path
import json,csv,statistics,datetime,time
CPU=time.process_time();ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
readset=[]
def load(name):
    path=ROOT/name;readset.append(name);return json.loads(path.read_text())
scores=load('RESULTS/confirmation_scoring/SCORES.json')['rows']
summary=load('RESULTS/summary/CONFIRMATION_SUMMARY.json')
diagnostic=load('audit/CONFIRMATION_FAILURE_DIAGNOSTIC.json')
weak=load('RESULTS/weak_handoff_scoring/SCORES.json')['rows'][0]
weak_receipt=load('RESULTS/development/obs_40ac5b79e8a74ac3__weak_profile_handoff/RECEIPT.json')
profile=load('RESULTS/development/obs_40ac5b79e8a74ac3__profile_homotopy/RECEIPT.json')
budget=load('RESULTS/development/obs_40ac5b79e8a74ac3__profile_homotopy/BUDGET.json')
ledger=load('audit/ALL_COSTS.json')
comparisons=[];matrix=[]
for method in ('cie_single','cie_dual','global_dual','a_plus'):
    for rule in ('fixed_001','noise_ball_99'):
        rows=[r for r in scores if r['method']==method and r['dpcs_rule']==rule]
        assert len(rows)==6 and len({r['case_id'] for r in rows})==6
        success=sorted(r['case_id'] for r in rows if r['material_image_success'])
        old=next(r for r in summary if r['method']==method and r['rule']==rule)
        recomputed=dict(method=method,rule=rule,n=6,successes=len(success),success_cases=success,
            mean_NRMSE=statistics.mean(r['material_nrmse'] for r in rows),mean_CPU_seconds=statistics.mean(r['cpu_seconds'] for r in rows),
            image_and_noise_success_cases=sorted(r['case_id'] for r in rows if r['same_grid_material_noise_success']))
        assert len(success)==old['successes'] and abs(recomputed['mean_NRMSE']-old['mean_NRMSE'])<1e-12 and abs(recomputed['mean_CPU_seconds']-old['mean_CPU_seconds'])<1e-12
        comparisons.append(recomputed)
        matrix.extend({k:r[k] for k in ('task_id','case_id','method','dpcs_rule','material_nrmse','otsu_iou','material_image_success','noise_CN99_pass','same_grid_material_noise_success','cpu_seconds','original_Phi')} for r in rows)
sets={r['method']:set(r['success_cases']) for r in comparisons if r['rule']=='noise_ball_99'}
unique=sorted(sets['a_plus']-set.union(*(sets[m] for m in ('cie_single','cie_dual','global_dual'))))
assert unique==[] and all(len(r['success_cases'])==(2 if r['rule']=='fixed_001' else 4) for r in comparisons)
endpoint=next(r for r in budget if r['label']=='profile_0_0.01_endpoint')
prefix=profile['pre_solver_CPU_seconds']+endpoint['CPU_seconds'];new=weak_receipt['process_CPU_seconds'];standalone=new+prefix
assert abs(standalone-weak['standalone_CPU_seconds'])<1e-12
actual_weak=next(r for r in ledger['inverse_receipts'] if r['task_id']==weak['task_id'])
assert abs(actual_weak['CPU_seconds']-new)<1e-12
recipe=next(r for r in ledger['standalone_weak_recipes_NOT_ADDITIONAL_ACTUAL_COST'] if r['task_id']==weak['task_id'])
witness=[]
for case in diagnostic['records']:
    witness.append(dict(case_id=case['case_id'],GT_feasible_Phi=case['GT_feasible_Phi'],GT_data_term=case['GT_data_term'],GT_prior=case['GT_prior'],
        arithmetic_difference=case['GT_feasible_Phi']-case['GT_data_term']-case['GT_prior'],GT_legal_box=case['GT_legal_box'],same_discretization=case['same_discretization'],source=case['stored_clean_field_is_feasible_prediction']))
    assert abs(witness[-1]['arithmetic_difference'])<1e-15 and case['GT_legal_box'] and case['same_discretization']
example=next(r for r in scores if r['task_id']=='r9_confirmation_0000__cie_single__fixed_001')
result=dict(status='COMPLETED_PHASEI_METADATA_REVIEW',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    fatal_errors_found=[],wording_issues=['THEORY section1 G_off should explicitly denote local disjoint-block off-diagonal operator, not entire physical offdiagonal G.',
        'Standalone19.63CPU is a measured-prefix attribution recipe, not a separately executed cold end-to-end19.63CPU process.'],
    noise_comparisons=comparisons,Aplus_unique_success_cases=unique,
    Aplus_vs_noise_single_CIE_CPU_percent=(next(r['mean_CPU_seconds'] for r in comparisons if r['method']=='a_plus' and r['rule']=='noise_ball_99')/next(r['mean_CPU_seconds'] for r in comparisons if r['method']=='cie_single' and r['rule']=='noise_ball_99')-1)*100,
    weak_handoff=dict(material_nrmse=weak['material_nrmse'],IoU=weak['otsu_iou'],new_execution_CPU_seconds=new,
        prior_profile_paid_startup_CPU_seconds=profile['pre_solver_CPU_seconds'],prior_profile_paid_firststage_CPU_seconds=endpoint['CPU_seconds'],
        reused_prefix_CPU_seconds=prefix,standalone_method_CPU_seconds=standalone,original_full_profile_CPU_seconds=profile['process_CPU_seconds'],
        actual_ledger_weak_charge=actual_weak['CPU_seconds'],standalone_recipe=recipe,new_execution_LU=weak['LU'],standalone_LU=weak['standalone_LU'],
        current_Cholesky_added_to_standalone=weak['standalone_current_Cholesky'],prefix_not_added_twice=True),
    feasible_witness=witness,common_failures=diagnostic['common_failures'],witness_qualification=diagnostic['qualification'],
    denominator_example=dict(task=example['task_id'],RMSE=example['RMSE'],declared_amplitude=example['amplitude'],amplitude_normalized_NRMSE=example['material_nrmse'],
        relative_L2_material_error=example['relative_l2_material_error'],scope='Distinct denominators, as reportsection2 correctly states; no threshold interchange'),
    maximum_saved_Phi_receipt_difference=max(abs(r.get('Phi_receipt_difference',0)) for r in scores),
    scope=dict(completed_confirmation_rows=66,actual_independent_material_objects=6,ordinary_plus_Aplus_rows_rechecked=48,
        no_new_forward=True,no_new_scoring=True,no_private_GT_arrays_read=True,no_inflight_classic_original_GT_read=True,no_3D_review=True,raw_files_modified=False),
    metadata_readset=readset,review_CPU_seconds=time.process_time()-CPU)
(OUT/'EVIDENCE.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
with (OUT/'PAIRED_SUCCESS_MATRIX.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(matrix[0]));writer.writeheader();writer.writerows(matrix)
print(json.dumps(dict(status=result['status'],fatal_errors=0,wording_issues=2,Aplus_unique_success_cases=unique,weak_standalone_CPU_seconds=standalone),ensure_ascii=False))
