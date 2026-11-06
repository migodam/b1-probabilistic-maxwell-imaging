"""Offline only after all nine seals; stored predictions, zero new physics."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import json,csv,time
from pathlib import Path
import numpy as np
from scipy import linalg as la
BASE=Path(__file__).resolve().parents[1];OUT=BASE/'RESULTS'/'original_vector3d'

def main():
    wall=time.perf_counter();protocol=json.loads((BASE/'configs/original_cm_nine.json').read_text())
    seal=json.loads((OUT/'SEALED_ALL_NINE.json').read_text())
    if len(seal['terminals'])!=9:raise RuntimeError('All9 endpoints must be sealed before any material GT read')
    for c in protocol['cohorts']:
        for arm in c['arms']:
            if not (OUT/'cohorts'/c['name']/'inversions'/arm/'SEALED.json').exists():raise RuntimeError('Missing seal')
    # First material-GT read occurs only after all guards above.
    truth=np.load(OUT/'sealed_truth/truth.npz')['inverse_material'];background=truth<=1e-10
    gt_prediction=np.load(OUT/'sealed_truth/coarse_truth_response.npz')['prediction'];rows=[];all_stages=[];checks=[];witness={}
    cache={e['id']:e for e in json.loads((OUT/'NINE_PHYSICAL_CACHE.json').read_text())['entries']}
    for c in protocol['cohorts']:
        name=c['name'];Y=np.load(OUT/'cohorts'/name/'public/observations.npz')['Y']
        gtphi=float(np.mean(np.sum(abs(gt_prediction-Y)**2,axis=0)/np.sum(abs(Y)**2,axis=0)))
        witness[name]=dict(Phi_at_same_known_feasible_material=gtphi,prior=0,
            scope='Stored paid CM30 truth-response; no new forward or iterate/parameter selection; not a minimum or uniqueness claim')
        for arm in c['arms']:
            d=OUT/'cohorts'/name/'inversions'/arm;summary=json.loads((d/'SUMMARY.json').read_text());cost=json.loads((d/'COST.json').read_text())
            final_history=json.loads((d/'history.jsonl').read_text().splitlines()[-1])
            endpoints=json.loads((d/'ENDPOINTS.json').read_text())
            reused=[cache[e['physical_cache_id']] for e in endpoints if e['physical_cache_hit']]
            refillCPU=sum(e['evaluation_CPU_seconds'] for e in reused);refillWall=sum(e['evaluation_wall_seconds'] for e in reused)
            for ep in endpoints:
                m=np.load(d/(ep['name']+'.npz'))['material'];pred=np.load(d/(ep['name']+'_prediction.npz'))['prediction']
                phi=float(np.mean(np.sum(abs(pred-Y)**2,axis=0)/np.sum(abs(Y)**2,axis=0)))
                checks.append(dict(cohort=name,arm=arm,endpoint=ep['name'],stored_prediction_Phi_difference=phi-ep['data_loss_mean']))
                row=dict(cohort=name,data_id=c['data_id'],arm=arm,endpoint=ep['name'],primary=ep['name']=='stage2',
                    paper_Err=float(np.sqrt(np.mean(abs(m-truth)**2/abs(1+truth)**2))),contrast_NRMSE=float(la.norm(m-truth)/la.norm(truth)),
                    Phi=phi,physical_data_RMS=float(np.sqrt(phi)),Phi_GT_feasible=gtphi,Phi_minus_GT_witness=phi-gtphi,
                    peak_contrast=float(m.max()),minimum_contrast=float(m.min()),lower_bound_fraction=float(np.mean(m<=1e-10)),
                    upper_bound_fraction=float(np.mean(m>=4-1e-10)),integrated_contrast_lambda3=float(m.sum()*.001),
                    effective_volume_lambda3=float(m.sum()**2/max(np.sum(m*m),1e-300)*.001),
                    background_mean=float(m[background].mean()),background_median=float(np.median(m[background])),
                    background_p95=float(np.quantile(m[background],.95)),background_max=float(m[background].max()),
                    max_linear_state_residual=ep['max_relative_linear_residual'],physical_cache_hit=ep['physical_cache_hit'],
                    primary_selection='Finalstage2 fixed before GT; stage1 not selected by score')
                all_stages.append(row)
                if row['primary']:
                    row.update(stage1_iterations=summary['stages'][0]['iterations'],stage2_iterations=summary['stages'][1]['iterations'],
                        stage1_stop=summary['stages'][0]['stop'],stage2_stop=summary['stages'][1]['stop'],
                        final_auxiliary_data=final_history['auxiliary_data'],final_auxiliary_state=final_history['auxiliary_state'],
                        final_auxiliary_sum=final_history['auxiliary_sum'],final_Eq16_delta_RMS=final_history['coefficient_delta_RMS'],
                        final_Eq16_delta_squared_mean=final_history['coefficient_delta_squared_mean'],
                        raw_CPU_seconds=cost['CPU_seconds'],setup_amortized_CPU_seconds=cost['setup_amortized_CPU_seconds'],
                        raw_wall_seconds=cost['wall_seconds'],setup_amortized_wall_seconds=cost['setup_amortized_wall_seconds'],
                        reused_prediction_CPU_attribution=refillCPU,reused_prediction_wall_attribution=refillWall,
                        reconstructed_method_CPU_estimate=cost['setup_amortized_CPU_seconds']+refillCPU,
                        reconstructed_method_wall_estimate=cost['setup_amortized_wall_seconds']+refillWall,
                        actual_physical_RHS=cost['counts'].get('primal_RHS',0),physical_cache_hits=cost['counts'].get('physical_prediction_cache_hits',0),
                        cost_scope='Nonadditive per-arm view; reconstructed adds paid first-evaluation cost for cache hits with same declared setup shares; attribution estimate, not independent cold-start timing or standalone speed ranking')
                    rows.append(row)
    qa={'stored_Phi_checks':checks,'maximum_absolute_stored_Phi_difference':max(abs(x['stored_prediction_Phi_difference']) for x in checks),
        'sourcewise_monotonicity':[],'same_initial_sourceCIE_Dual_control':None}
    for c in protocol['cohorts']:
        for arm in ['sourceCIE','LSIEsource']:
            d=OUT/'cohorts'/c['name']/'inversions'/arm;h=[json.loads(x) for x in (d/'history.jsonl').read_text().splitlines()]
            for stage in [0,1]:
                v=np.array([x['auxiliary_sum'] for x in h if x['stage']==stage]);dv=np.diff(v)
                qa['sourcewise_monotonicity'].append(dict(cohort=c['name'],arm=arm,stage=stage,
                    maximum_signed_change=float(dv.max()) if len(dv) else None,max_absolute_increase=float(max(0,dv.max())) if len(dv) else 0.))
    cross=OUT/'cohorts/cross_grid/inversions'
    a=np.load(cross/'sourceCIE/stage1.npz')['material'];b=np.load(cross/'sourceDual/stage1.npz')['material']
    qa['same_initial_sourceCIE_Dual_control']=dict(max_absolute_material_difference=float(np.max(abs(a-b))),
        relative_material_difference=float(la.norm(a-b)/max(la.norm(a),1e-300)),scope='Same first-stage method by frozen formula; no new solve')
    old=OUT/'engineering_attempt1_strided4/inversions/sourceCIE/stage1.npz'
    if old.exists():
        b=np.load(old)['material'];qa['packed_vs_old_strided_paid_stage1']=dict(max_absolute_material_difference=float(np.max(abs(a-b))),
            relative_material_difference=float(la.norm(a-b)/max(la.norm(a),1e-300)),scope='Existing paid endpoint comparison, no GT choice or rerun')
    result=dict(status='OFFLINE_SCORED_AFTER_ALL_NINE_SEALED',primary_rows=rows,all_paid_stage_rows=all_stages,
        witnesses=witness,QA=qa,physics_evaluations=0,scope=protocol['scope'],
        initial_air_reference=[{k:r[k] for k in ['cohort','paper_Err','contrast_NRMSE','Phi']} for r in all_stages if r['endpoint']=='air' and r['arm']=='sourceCIE'],
        dimensions=dict(real_material_unknowns=27000,complex_observations=180*60,real_data_component_upper_bound=21600,
            scope='Count is an upper information-dimension bound, not global uniqueness; box-boundary GT is not proven nonunique by count alone; do not pool with small3D512unknowns'),
        Phi_definition='mean_source sum_receivers|prediction-Y|^2 / sum_receivers|Y|^2; no1/2; R=0; cannot numerically compare to PhaseI Phi',
        data_preparation_access='Sealed container was opened only for Y_clean/noise/old coarse prediction; no material fields read; inverse receives only public Y/operator/box/grid',
        limitation='Two separate data cohorts; matched finite-grid diagnostic is not same-data algorithm gain or continuum/paper reproduction')
    (OUT/'NINE_OFFLINE_SCORES.json').write_text(json.dumps(result,indent=2)+'\n')
    with (OUT/'NINE_PRIMARY_SCORES.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    (OUT/'NINE_SCORING_COST.json').write_text(json.dumps(dict(CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall,
        platform='Mac',physics_evaluations=0,scope='Offline material scoring from stored predictions after all9seals; CPU includes imports',additive=True,included_in=None),indent=2)+'\n')

if __name__=='__main__':main()
