"""POSTHOC finite-noise feasible pairs from existing GT currents; no Maxwell solve."""
from pathlib import Path
import json,csv,time
import numpy as np
from scipy.stats import chi2
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
def main():
    start=time.perf_counter();cpu=time.process_time();scores=json.loads((P/'SCORES_FINAL.json').read_text());cie=[r for r in scores if r['method']=='cie_single'];rows=[]
    with np.load(ROOT/'data_public/twod/operators.npz',allow_pickle=False) as z:S=z['S'].copy()
    for r in sorted(cie,key=lambda r:r['case_id']):
        case=r['case_id']
        with np.load(ROOT/'data_public/twod'/(case+'.npz'),allow_pickle=False) as z:Y=z['Y'].copy();sigma=float(z['sigma'])
        with np.load(ROOT/'delegated/twod_energy/gt_current_cache'/(case+'.npz'),allow_pickle=False) as z:gtm=z['m_gt'].copy();jgt=z['j_gt'].copy()
        with np.load(ROOT/'RESULTS/twod'/r['task_id']/'prediction.npz',allow_pickle=False) as z:yp=z['Y_prediction'].copy();m=z['m'].copy()
        ygt=S@jgt;radius=float(np.sqrt(chi2.ppf(.99,2*Y.size)/2));gt_res=float(np.linalg.norm(ygt-Y)/sigma);est_res=float(np.linalg.norm(yp-Y)/sigma)
        glegal=bool(np.isfinite(gtm).all() and not np.iscomplexobj(gtm) and np.min(gtm)>=.02-1e-12 and np.max(gtm)<=20+1e-12)
        both=bool(gt_res<=radius and est_res<=radius and glegal and r['material_legal'])
        rows.append(dict(case_id=case,regime=r['regime'],family=r['family'],amplitude_NRMSE=r['amplitude_NRMSE'],fixed_support_IoU=r['fixed_support_IoU'],
            material_success=r['material_success'],GT_legal_box=glegal,CIE_legal_box=r['material_legal'],sigma=sigma,complex_observations=Y.size,
            single_case_CN99_radius=radius,GT_whitened_residual=gt_res,CIE_whitened_residual=est_res,GT_CN99_pass=bool(gt_res<=radius),CIE_CN99_pass=bool(est_res<=radius),
            both_legal_materials_in_actual_y_CN99_ball=both,material_failed_pair_witness=bool(both and not r['material_success']),
            material_difference_squared_energy=float(np.linalg.norm(m-gtm)**2),pair_data_separation_in_sigma=float(np.linalg.norm(yp-ygt)/sigma),
            original_Phi_CIE=r['original_Phi'],soft_prior_semantics='Shared soft smooth/TV penalty is an objective preference, not a hard admissibility condition; equal/noise-compatible data do not imply equal prior or Phi'))
    with (P/'CIE_GT_CN99_FEASIBLE_PAIRS.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary=dict(status='POSTHOC_OFFLINE_NOISE_COMPATIBILITY_DIAGNOSTIC',cases=24,
        CIE_CN99_failures=[r['case_id'] for r in rows if not r['CIE_CN99_pass']],GT_CN99_failures=[r['case_id'] for r in rows if not r['GT_CN99_pass']],
        CIE_material_successes=sum(r['material_success'] for r in rows),CIE_CN99_pass_count=sum(r['CIE_CN99_pass'] for r in rows),
        CIE_data_pass_material_fail_count=sum(r['CIE_CN99_pass'] and not r['material_success'] for r in rows),
        material_failed_pair_witnesses=[r['case_id'] for r in rows if r['material_failed_pair_witness']],
        scope='Each case is a separate99% circular-complex-noise ball around actual observedY. No24-object joint99% event is claimed. Feasible finite-noise pairs are not exact multiple solutions or information-theoretic lower bounds.',
        prior_scope='Both satisfy real box; shared soft prior can still prefer one. No hard prior-budget constraint or equality of posterior/regularized objective inferred.',
        GT_source='Previously authorized offline GT-current cache; new sealed_truth reads0',new_Maxwell_forward_calls=0,new_Maxwell_LU=0,new_inverse_runs=0,
        existing_current_GT_receiver_operator_actions=24,existing_current_GT_receiver_complex_columns=288,
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start)
    (P/'CIE_GT_CN99_FEASIBLE_PAIRS_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
