"""Compact release summary from sealed scores and saved energy tables only."""
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime,timezone
import csv,json,statistics
S=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text())
def mean(v):return statistics.fmean(v)if v else None
def main():
    score=read(S/'results/TEST_SCORES.json');final=read(S/'audit/FINAL_RUNS.json')
    groups=defaultdict(list)
    for r in score['rows']:
        for group in ('ALL',r['group']):groups[r['run'],r['noise'],group].append(r)
    table=[]
    for (run,noise,group),rows in groups.items():
        available=[r for r in rows if r.get('top1')is not None]
        def values(key):return [r[key]for r in rows if r.get(key)is not None]
        def candidate_values(key,metric='NRMSE'):return [r[key][metric]for r in rows if r.get(key)is not None]
        raw=[next((c for c in r.get('candidates',[])if c['label']=='raw_1'),None)for r in rows]
        table.append(dict(run=run,noise=noise,group=group,planned_observations=len(rows),independent_mothers=len({r['parent_id']for r in rows}),
            saved=sum(r['transport_status']=='SEALED'for r in rows),qualified=len(available),not_run=sum(r['transport_status']=='NOT_RUN'for r in rows),
            top1_NRMSE_mean=mean(candidate_values('top1')),top1_NRMSE_median=statistics.median(candidate_values('top1'))if available else None,
            top1_RMSE_mean=mean(candidate_values('top1','RMSE')),top1_IoU_mean=mean(candidate_values('top1','IoU_Otsu')),
            top1_full_VIE_relative_residual_mean=mean(candidate_values('top1','relative_residual')),
            raw_1_NRMSE_mean=mean([r['NRMSE']for r in raw if r is not None]),
            joint_success=sum(bool(r['joint_success'])for r in rows),joint_success_total55=sum(bool(r.get('joint_success_total55'))for r in rows),
            old_data_TV_top1_NRMSE_mean=mean(candidate_values('old_score_top1')),
            qualified_bestK_NRMSE_mean=mean(candidate_values('qualified_bestK_OFFLINE')),
            all_bestK_NRMSE_mean=mean(candidate_values('all_bestK_OFFLINE')),
            actual_saved_K_histogram=dict(Counter(len(r.get('candidates',[]))for r in rows)),
            online_total_wall_seconds_mean=mean(values('online_total_wall_seconds')),
            online_total_LU_mean=mean(values('online_total_LU')),online_total_RHS_mean=mean(values('online_total_RHS')),
            method_status_counts=dict(Counter(r['status']for r in rows)),
            CPU_seconds='NOT_RECORDED_FOR_MAIN_SOLVERS'))
    energy=defaultdict(lambda:dict(cases=0,before=0.,after=0.,delta=0.))
    bands=('P_in_strong19','P_in_weak557','S_strong9','S_weak39','S_null528','D509','B67','G_field','S_field')
    with(S/'results/test_energy/grouped_absolute_energy.csv').open()as f:
        for r in csv.DictReader(f):
            if r['quantity']!='Maxwell_current'or r['stage']!='anchor_to_final_top1'or r['projection']not in bands:continue
            d=energy[r['run'],float(r['noise']),r['projection']];d['cases']+=int(r['cases'])
            for key in ('before','after','delta'):d[key]+=float(r[key+'_sum'])
    eout=[dict(run=run,noise=noise,projection=band,cases=d['cases'],**{k+'_mean':d[k]/d['cases']for k in ('before','after','delta')})for(run,noise,band),d in energy.items()]
    calibration=score['GT_shape_metric_calibration']['parents']
    result=dict(schema='gn-nn-r4-compact-results-v1',status='ALL_26_RUNS_216_OBSERVATIONS_SCORED_AFTER_SEAL',created_utc=datetime.now(timezone.utc).isoformat(),
        physical_model=dict(type='nonlinear2DscalarTMdiskpulseVIE',grid=[24,24],unknown_real_materials=576,sources=12,receivers=48,k='2pi',material_box=[.02,20.],training_amplitudes=[.1,.2,.4]),
        data=dict(TRAIN_mothers=504,VAL_mothers=78,TEST_mothers=108,TEST_observations=216,groups=dict(new_instance=60,unseen_class=18,unseen_style=30),held_out_classes=['digit4','digit7','curvedbar'],test_images_in_training=False),
        runs=final['runs'],NN_weight_sets=len(final['nn_runs']),planned_run_observations=26*216,actual_scored_rows=len(score['rows']),
        transport_counts=dict(Counter(r['transport_status']for r in score['rows'])),method_status_counts=dict(Counter(r['status']for r in score['rows'])),
        results=table,true_current_absolute_energy=eout,
        GT_self_Otsu_shape_gate=dict(passing=sum(r['GT_self_shape_gate']for r in calibration),parents=len(calibration),interpretation='metric calibration only, not a strict attainable ceiling; frozen gate unchanged'),
        comparisons='paired independent mothers; noise variants and model seeds are related repetitions; no test-based model selection',
        selection='physical qualification then proximity to full NN proposal; offline bestK is method-specific actual finiteK, no posterior or root enumeration',
        cost_scope='measured original hybridM4/i9/MPS pipeline; shared strictGN included forNN, unusedNNfeatures excluded from ordinary deployment; training/setup separate; LU/RHS omit fixedmatmul/NN/CIEoperations billed in ledger',
        interpretation=dict(algorithm_paper='HOLD',reason='Current-space constraints/feedback do not consistently beat Voxel or the efficient noise-aware ordinary baseline; generic current/material unrolling has direct SOM-Net/CSI precedent. Mechanism evidence and restricted positive cases are retained.'),
        artifacts=dict(scores='TEST_SCORES.json',per_case='TEST_SCORES.csv',cost='COST_LEDGER.csv',energy='test_energy/',independent_review='../delegated/review/INDEPENDENT_TEST_FINDINGS.md',main_images='../figures/final/main/',atlas='../figures/all_test_cases/GALLERY.html',recurrent='../RECURRENT_FINDINGS_ZH.md'),
        NOT_RUN=['distinctinteriorSoftDark(alpha1selected,reusesJfull)','feedback2training','trueA17Schurfeedbackpreconditioner','recurrentseed17/27','coldstartNNwithoutper-targetGN','newR4highcontrastconfirmation'])
    p=S/'results/RESULTS_SUMMARY.json';p.write_text(json.dumps(result,indent=2)+'\n');print(p)
if __name__=='__main__':main()
