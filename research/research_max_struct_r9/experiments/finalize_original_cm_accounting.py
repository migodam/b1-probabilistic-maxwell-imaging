"""Read-only receipt aggregation. No GT/arrays, solves, or GPU events as CPU."""
import json,csv
from pathlib import Path
BASE=Path(__file__).resolve().parents[1];OUT=BASE/'RESULTS'/'original_vector3d'

def main():
    if not (OUT/'SEALED_ALL_NINE.json').exists() or not (OUT/'run_nine_CM_PROCESS_COST.json').exists():raise RuntimeError('Nine seals and final process receipt required')
    if json.loads((OUT/'run_nine_CM_PROCESS_COST.json').read_text())['status']!='COMPLETE':raise RuntimeError('Nine process not complete')
    rows=[]
    historical=OUT/'ACCOUNTING_BEFORE_PHYSICS_DIAG.json'
    if historical.exists():
        for r in json.loads(historical.read_text())['rows']:
            r=dict(r);r.update(additive=True,included_in=None)
            detail=json.loads((OUT/r['artifact']).read_text())
            r['counts']=detail.get('counts',{})
            if 'CPU_counts' in detail:r['counts']={'CPU':detail['CPU_counts'],'CUDA':detail['CUDA_counts'],'completed_primal_RHS':2}
            if 'completed_auxiliary_steps' in detail:r['counts']={k:detail[k] for k in ['completed_auxiliary_steps','completed_physical_RHS','recorded_physical_Arnoldi_steps','actual_full_counter_snapshot_available']}
            rows.append(r)
    else:
        # A fresh reproduction need not repeat historical aborted engineering.
        for f,label in [('CUDA_WHOLE_RHS_BENCHMARK.json','Complete source59 CPU/GPU qualification'),('generate_CUDA_COST.json','CM60 data generation'),('qualify-discretization_CUDA_COST.json','Paid CM30 reference')]:
            if not (OUT/f).exists():continue
            a=json.loads((OUT/f).read_text());rows.append(dict(artifact=f,operation=label,platform='XINAN',CPU_seconds=a['CPU_seconds'],wall_seconds=a['wall_seconds'],terminal_ledger_status=a.get('status','COMPLETE'),scope='Disjoint recorded process; see timer scope',additive=True,included_in=None,counts=a.get('counts',{})))
    for f in ['igt_source59/prepare_PROCESS_COST.json','igt_source59/solve_PROCESS_COST.json','run_nine_CM_PROCESS_COST.json','NINE_SCORING_COST.json','NINE_PLOT_COST_ATTEMPT1_CPU_AXIS.json','NINE_PLOT_COST.json']:
        if not (OUT/f).exists():continue
        a=json.loads((OUT/f).read_text());process_counts=a.get('counts',{})
        if f.startswith('igt_source59/'):
            child='SOLVE_COST.json' if '/solve_' in f else 'PREPARE_COST.json'
            case_counts={str(n):json.loads((OUT/'igt_source59'/str(n)/child).read_text()).get('counts',{}) for n in [30,60]}
            process_counts={'by_grid_nonadditive_count_breakdown':case_counts,'completed_primal_RHS':2 if child=='SOLVE_COST.json' else 0,
                'CPU_original_matrix_verification_Green_RHS':2 if child=='SOLVE_COST.json' else 0}
        machine='XINAN' if f.startswith('igt_source59/') or f=='run_nine_CM_PROCESS_COST.json' else a.get('platform','Mac')
        rows.append(dict(artifact=f,operation=f,platform=machine,
            CPU_seconds=a['CPU_seconds'],wall_seconds=a['wall_seconds'],terminal_ledger_status=a.get('status','COMPLETE'),scope=a.get('scope','Whole process; child timers included'),
            additive=True,included_in=None,counts=process_counts))
    for f in sorted((OUT/'resumed_generation_attempts').glob('attempt*/generate_CUDA_COST.json')):
        a=json.loads(f.read_text());rows.append(dict(artifact=str(f.relative_to(OUT)),operation='Preserved generation attempt before checkpoint continuation',
            platform='XINAN',CPU_seconds=a['CPU_seconds'],wall_seconds=a['wall_seconds'],terminal_ledger_status=a.get('status','RECORDED_INTERRUPTED_ATTEMPT'),
            scope='Disjoint prior generation process; qualified state may be reused but original work remains charged',additive=True,included_in=None,counts=a.get('counts',{})))
    views=[]
    files=[OUT/'NINE_SHARED_OPERATOR_COST.json',OUT/'NINE_DATA_PREPARATION_COST.json']+list((OUT/'cohorts').glob('*/RECEIVER_Y_SETUP_COST.json'))+list((OUT/'cohorts').glob('*/inversions/*/COST.json'))
    for f in files:
        if not f.exists():continue
        a=json.loads(f.read_text());views.append(dict(artifact=str(f.relative_to(OUT)),platform='XINAN',CPU_seconds=a.get('CPU_seconds'),wall_seconds=a.get('wall_seconds'),
            additive=False,included_in='run_nine_CM_PROCESS_COST.json',scope='Measured nonadditive interval/setup view; physical cache order applies',counts=a.get('counts',{})))
    scored=OUT/'NINE_OFFLINE_SCORES.json';attribution=[]
    if scored.exists():
        for a in json.loads(scored.read_text())['primary_rows']:
            attribution.append({k:a[k] for k in ['cohort','arm','raw_CPU_seconds','setup_amortized_CPU_seconds','raw_wall_seconds','setup_amortized_wall_seconds',
                'reused_prediction_CPU_attribution','reused_prediction_wall_attribution','reconstructed_method_CPU_estimate','reconstructed_method_wall_estimate','cost_scope']})
    result=dict(status='FINAL_RECORDED_ACCOUNTING_AFTER_NINE',rows=rows,nonadditive_views=views,
        CPU_seconds_disjoint_sum=sum(r['CPU_seconds'] for r in rows),platform_CPU_seconds={p:sum(r['CPU_seconds'] for r in rows if r['platform']==p) for p in sorted(set(r['platform'] for r in rows))},
        wall_scope='Per-process wall values overlap and old paused wall is included; no sum-as-makespan claim',
        physical_cache_attribution_estimates=attribution,attribution_additive=False,
        attribution_scope='Paid first-evaluation costs restored to cache-hit methods with same declared setup shares; not new independent cold-start timings; never add estimates to process total',
        unknowns=json.loads(historical.read_text()).get('unknowns',[]) if historical.exists() else [],
        reproduction_scope='Fresh output directory reproduces final protocol; historical aborted engineering is preserved but need not be rerun')
    (OUT/'ACCOUNTING_FINAL.json').write_text(json.dumps(result,indent=2)+'\n')
    with (OUT/'ACCOUNTING_FINAL_PROCESS_ROWS.csv').open('w',newline='') as f:
        fields=['artifact','platform','CPU_seconds','wall_seconds','terminal_ledger_status','additive','included_in','scope'];w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
    print(json.dumps({'status':result['status'],'additive_rows':len(rows),'CPU_seconds':result['CPU_seconds_disjoint_sum'],'platforms':result['platform_CPU_seconds']}))

if __name__=='__main__':main()
