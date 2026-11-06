"""POSTHOC saved-endpoint deployment diagnostic; selector uses original Phi only."""
from pathlib import Path
import json,csv,time
from collections import Counter
import numpy as np
P=Path(__file__).resolve().parent
def csvwrite(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def main():
    start=time.perf_counter();cpu=time.process_time();rr=json.loads((P/'SCORES_FINAL.json').read_text());ix={(r['case_id'],r['method']):r for r in rr};cases=sorted({r['case_id'] for r in rr});rows=[]
    for case in cases:
        a,b=ix[case,'cie_single'],ix[case,'matrix_green19']
        # Selection intentionally touches only public original-Phi fields.
        chosen=a if a['original_Phi']<=b['original_Phi'] else b
        rows.append(dict(case_id=case,regime=a['regime'],status='POSTHOC_DIAGNOSTIC_NOT_CONFIRMATION_ALGORITHM',
            selection='minimum final saved original Phi; exact ties choose CIE',selected_method=chosen['method'],
            CIE_original_Phi=a['original_Phi'],G19_original_Phi=b['original_Phi'],selected_original_Phi=chosen['original_Phi'],
            amplitude_NRMSE=chosen['amplitude_NRMSE'],fixed_support_IoU=chosen['fixed_support_IoU'],material_success=chosen['material_success'],
            legacy_Otsu_IoU=chosen['legacy_Otsu_IoU'],Otsu_full_support_IoU=chosen['Otsu_full_support_IoU'],
            CIE_amplitude_NRMSE=a['amplitude_NRMSE'],CIE_fixed_support_IoU=a['fixed_support_IoU'],CIE_material_success=a['material_success'],
            CIE_existing_success_lost=bool(a['material_success'] and not chosen['material_success']),
            new_material_success_over_CIE=bool(not a['material_success'] and chosen['material_success']),
            all_two_arms_full_child_CPU_seconds=a['external_full_child_CPU_seconds']+b['external_full_child_CPU_seconds'],
            all_two_arms_physical_LU=a['physical_LU']+b['physical_LU'],all_two_arms_physical_complex_RHS=a['physical_complex_RHS']+b['physical_complex_RHS'],
            G19_small_LU=b['denominator_small_LU'],G19_small_LU_dimension=b['denominator_small_LU_dimension'],
            G19_small_triangular_complex_RHS=b['denominator_small_triangular_complex_RHS']))
    summary=[]
    for regime in ('all','low','high','mixed'):
        z=[r for r in rows if regime=='all' or r['regime']==regime]
        cie=[ix[r['case_id'],'cie_single'] for r in z]
        summary.append(dict(regime=regime,cases=len(z),selected_methods=dict(Counter(r['selected_method'] for r in z)),
            mean_amplitude_NRMSE=float(np.mean([r['amplitude_NRMSE'] for r in z])),median_amplitude_NRMSE=float(np.median([r['amplitude_NRMSE'] for r in z])),
            mean_fixed_support_IoU=float(np.mean([r['fixed_support_IoU'] for r in z])),successes=sum(r['material_success'] for r in z),
            CIE_successes=sum(r['material_success'] for r in cie),existing_CIE_success_lost=[r['case_id'] for r in z if r['CIE_existing_success_lost']],
            new_successes=[r['case_id'] for r in z if r['new_material_success_over_CIE']],
            all_two_arms_full_child_CPU_seconds=sum(r['all_two_arms_full_child_CPU_seconds'] for r in z),
            child_CPU_ratio_to_all_CIE=sum(r['all_two_arms_full_child_CPU_seconds'] for r in z)/sum(r['external_full_child_CPU_seconds'] for r in cie),
            all_two_arms_physical_LU=sum(r['all_two_arms_physical_LU'] for r in z),all_two_arms_physical_complex_RHS=sum(r['all_two_arms_physical_complex_RHS'] for r in z)))
    csvwrite(P/'POSTHOC_PORTFOLIO_OBJECTS.csv',rows)
    (P/'POSTHOC_PORTFOLIO_SUMMARY.json').write_text(json.dumps(dict(status='POSTHOC_DIAGNOSTIC_NOT_NEW_CONFIRMATION_ALGORITHM',
        root_authorization='Explicit posthoc CIE-single+G19 minimum original Phi diagnostic after fresh scoring',rows=summary,
        selection_uses_GT=False,GT_metrics_reused_from_frozen_scoring=True,new_GT_file_reads=0,new_forward_calls=0,new_inverse_runs=0,
        cost_scope='Every object pays both completed child lifetimes, physical LU/RHS and G19 small systems, including unselected arm. Shared study generation/transport/controller costs remain separate, not silently prorated.',
        scope='Compensatory deployment diagnostic created after seeing outcomes; no mathematical novelty or new confirmation claim. No thresholds/budgets/candidates changed.',
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start),indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
