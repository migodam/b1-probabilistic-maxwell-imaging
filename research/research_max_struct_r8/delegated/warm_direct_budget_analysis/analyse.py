"""Post-result descriptive extension; scored data only, no physics or GT access.

The paired_budget_endpoints implementation is copied from the existing
experiments/analyse_warm_results.py without changing selection or accounting.
Only its explicit contrast list and output names differ. Original scripts and
sealed results are unchanged. Run using the project's existing Python.
"""
from pathlib import Path
from collections import defaultdict, Counter
import csv,json,time
START_CPU=time.process_time()
START_WALL=time.perf_counter()
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
METRICS=('material_nrmse','otsu_iou','original_Phi','relative_residual',
         'warm_deployment_cpu_seconds','solver_cpu_seconds','LU','complex_RHS')
CONTRASTS=[('physical_lbfgs / nn',arm+' / '+init)
           for init in ('uniform','nn')
           for arm in ('cie_single','cie_dual','global_dual','a_plus','a_plus_sourcewise')]

def finite(x):
    return x is not None and np.isfinite(x)

def save(name,value):
    (OUT/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def table(name,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        w.writeheader();w.writerows(rows)

def paired_budget_endpoints(rows, budgets):
    """Per-object common allowance min(final A,final B); sparse paid incumbents.

    The allowance is identical. Actual checkpoint spend can differ, so unused
    allowance and the next checkpoint gap are explicit. A stopped branch never
    receives invented iterations. Comparisons are stratified by actual hardware.
    """
    meta={r['task_id']:r for r in rows}
    by={(r['case_id'],r['method']):r for r in rows}
    tracks=defaultdict(list)
    for b in budgets:
        if b['task_id'] in meta:tracks[b['task_id']].append(b)
    for bs in tracks.values():bs.sort(key=lambda b:b['index'])
    contrasts=CONTRASTS
    axes={'CPU_seconds':'warm_deployment_cpu_seconds',
          'solver_CPU_seconds':'solver_cpu_seconds','LU':'LU','complex_RHS':'complex_RHS'}
    def coordinate(b,axis):
        if axis not in ('CPU_seconds','solver_CPU_seconds'):return b.get(axis)
        r=meta[b['task_id']]
        if r.get('deployment_preprocessing_hardware')!=r['hardware'] or not finite(r.get('deployment_prefix_cpu_seconds')):
            return None
        return b['CPU_seconds']-(r['deployment_prefix_cpu_seconds'] if axis=='solver_CPU_seconds' else 0)
    def choose(rr,axis,allowance):
        bs=tracks[rr['task_id']]
        if not bs or any(not finite(coordinate(b,axis)) for b in bs):return None
        paid=[b for b in bs if coordinate(b,axis)<=allowance]
        if not paid:return None
        last=max(paid,key=lambda b:b['index'])
        incumbent=min(paid,key=lambda b:(b['original_Phi'],-b['index']))
        nextcosts=[coordinate(b,axis) for b in bs if coordinate(b,axis)>allowance]
        return dict(incumbent=incumbent,last_paid=coordinate(last,axis),
                    unused=max(0.,allowance-coordinate(last,axis)),
                    next_checkpoint_cost=min(nextcosts) if nextcosts else None,
                    gap=(min(nextcosts)-coordinate(last,axis)) if nextcosts else None,
                    age=coordinate(last,axis)-coordinate(incumbent,axis),
                    at_terminal_checkpoint=last['index']==max(b['index'] for b in bs))
    result=[]
    for method,reference in contrasts:
        for case in sorted({r['case_id'] for r in rows}):
            a,b=by.get((case,method)),by.get((case,reference))
            if not a or not b:continue
            for axis,finalkey in axes.items():
                row=dict(case_id=case,cohort=a['cohort'],hardware=a['hardware'],method=method,
                         reference=reference,axis=axis,method_status=a['status'],reference_status=b['status'],
                         selection='Original-Phi incumbent among checkpoints <= common allowance; latest index tie-break; never GT',
                         interpretation='Same at-most allowance, not necessarily equal actual expenditure')
                if a['hardware']!=b['hardware']:
                    row.update(status='MIXED_HARDWARE_PAIR_OMITTED',reference_hardware=b['hardware']);result.append(row);continue
                if not finite(a.get(finalkey)) or not finite(b.get(finalkey)):
                    row.update(status='MISSING_FINAL_COST');result.append(row);continue
                allowance=min(a[finalkey],b[finalkey]);row['common_allowance']=allowance
                aa,bb=choose(a,axis,allowance),choose(b,axis,allowance)
                if aa is None or bb is None:
                    row.update(status='MISSING_PAID_CHECKPOINT_AT_ALLOWANCE');result.append(row);continue
                row['status']='PAIRED_PAID_INCUMBENTS'
                for prefix,state in [('method',aa),('reference',bb)]:
                    row[prefix+'_last_paid_checkpoint_cost']=state['last_paid']
                    row[prefix+'_unused_allowance']=state['unused']
                    row[prefix+'_next_checkpoint_cost']=state['next_checkpoint_cost']
                    row[prefix+'_checkpoint_bracketing_gap']=state['gap']
                    row[prefix+'_incumbent_age_in_paid_cost']=state['age']
                    row[prefix+'_terminal_checkpoint_reached']=state['at_terminal_checkpoint']
                    row[prefix+'_incumbent_index']=state['incumbent']['index']
                    for metric in ('original_Phi','material_nrmse','otsu_iou'):
                        row[prefix+'_'+metric]=state['incumbent'].get(metric)
                for metric in ('original_Phi','material_nrmse','otsu_iou'):
                    va,vb=row.get('method_'+metric),row.get('reference_'+metric)
                    row['difference_'+metric]=va-vb if finite(va) and finite(vb) else None
                result.append(row)
    table('object_matched_at_most_budget_endpoints.csv',result)
    groups=defaultdict(list)
    for r in result:groups[(r['cohort'],r['hardware'],r['method'],r['reference'],r['axis'])].append(r)
    effects=[];rng=np.random.default_rng(2026110608)
    for (cohort,hardware,method,reference,axis),rr in sorted(groups.items()):
        for metric in ('original_Phi','material_nrmse','otsu_iou'):
            valid=[r for r in rr if finite(r.get('difference_'+metric))]
            delta=np.array([r['difference_'+metric] for r in valid])
            ci=[None,None]
            if len(delta):
                boot=np.mean(delta[rng.integers(0,len(delta),(10000,len(delta)))],axis=1)
                ci=np.quantile(boot,[.025,.975]).tolist()
            effects.append(dict(cohort=cohort,hardware=hardware,method=method,reference=reference,axis=axis,
                metric=metric,objects=len(rr),paired_objects=len(valid),missing_pairs=len(rr)-len(valid),
                mean_paired_difference=float(np.mean(delta)) if len(delta) else None,
                percentile95_low=ci[0],percentile95_high=ci[1],
                mean_method_unused_allowance=float(np.mean([r['method_unused_allowance'] for r in valid])) if valid else None,
                mean_reference_unused_allowance=float(np.mean([r['reference_unused_allowance'] for r in valid])) if valid else None,
                sampling_unit='independent object; 10000 descriptive paired resamples',
                budget_rule='Per-object min(final actual cost of the two complete methods); same allowance, sparse checkpoints, no interpolation'))
    table('object_matched_budget_paired_bootstrap.csv',effects)
    return result,effects


def final_pairs(rows):
    by={(r['case_id'],r['method']):r for r in rows}
    result=[]
    for method,reference in CONTRASTS:
        for case in sorted({r['case_id'] for r in rows}):
            a,b=by[(case,method)],by[(case,reference)]
            row=dict(case_id=case,cohort=a['cohort'],hardware=a['hardware'],method=method,
                     reference=reference,method_status=a['status'],reference_status=b['status'],
                     method_material_success=a['material_image_success'],reference_material_success=b['material_image_success'])
            if a['hardware']!=b['hardware']:
                row['status']='MIXED_HARDWARE_PAIR_OMITTED';result.append(row);continue
            row['status']='SAME_HARDWARE_FINAL_PAIR'
            for metric in METRICS:
                av,bv=a.get(metric),b.get(metric)
                row['method_'+metric]=av;row['reference_'+metric]=bv
                row['difference_'+metric]=av-bv if finite(av) and finite(bv) else None
            cost=a['warm_deployment_cpu_seconds'];ref=b['warm_deployment_cpu_seconds']
            row['method_over_reference_deployment_CPU_ratio']=cost/ref if finite(cost) and finite(ref) and ref>0 else None
            result.append(row)
    table('final_object_pairs.csv',result)
    groups=defaultdict(list)
    for r in result:groups[(r['cohort'],r['hardware'],r['method'],r['reference'])].append(r)
    effects=[];rng=np.random.default_rng(2026110619)
    for (cohort,hardware,method,reference),rr in sorted(groups.items()):
        for metric in METRICS:
            valid=[r for r in rr if finite(r.get('difference_'+metric))]
            d=np.array([r['difference_'+metric] for r in valid]);ci=[None,None]
            if len(d):
                boot=np.mean(d[rng.integers(0,len(d),(10000,len(d)))],axis=1)
                ci=np.quantile(boot,[.025,.975]).tolist()
            effects.append(dict(cohort=cohort,hardware=hardware,method=method,reference=reference,
                metric=metric,objects=len(rr),paired_objects=len(valid),missing_pairs=len(rr)-len(valid),
                method_mean=float(np.mean([r['method_'+metric] for r in valid])) if valid else None,
                reference_mean=float(np.mean([r['reference_'+metric] for r in valid])) if valid else None,
                mean_paired_difference=float(np.mean(d)) if len(d) else None,
                percentile95_low=ci[0],percentile95_high=ci[1],
                method_material_successes=sum(r['method_material_success'] for r in rr),
                reference_material_successes=sum(r['reference_material_success'] for r in rr),
                sampling_unit='independent object; 10000 descriptive paired resamples'))
    table('final_paired_effects.csv',effects)
    return result,effects

def main():
    score=json.loads((ROOT/'RESULTS/warm_scoring/SCORES.json').read_text())
    budget=json.loads((ROOT/'RESULTS/warm_scoring/BUDGET_SCORES.json').read_text())
    rows=[r for r in score['rows'] if r['primary']]
    assert len(rows)==630 and len({r['case_id'] for r in rows})==42
    cohorts={c:len({r['case_id'] for r in rows if r['cohort']==c}) for c in {r['cohort'] for r in rows}}
    assert sorted(cohorts.values())==[6,12,24]
    final,fe=final_pairs(rows)
    endpoints,ee=paired_budget_endpoints(rows,budget['rows'])
    assert len(final)==420 and len(endpoints)==1680
    status=dict(Counter(r['status'] for r in endpoints))
    report=['# Post-result direct-physics descriptive cost contrasts','',
        'These ten contrasts were added after seeing the primary results. They are descriptive post-result analyses, not new algorithms, confirmation tests, or preregistered significance tests. No training, solver, physics, interpolation, GT selection, or raw truth access was performed. Existing scored records alone were read.','',
        'All 42 primary objects are retained: ID 24, held-out family 12, cross-grid 6. Cohorts and actual hardware remain separate. Positive differences mean physical L-BFGS / NN is larger; lower error, Phi and cost are favorable, while higher IoU is favorable.','',
        'CPU is measured same-hardware deployment CPU, including frozen preprocessing and inference; solver-only CPU is additionally reported. Training and cold loading are excluded from this deployment comparison. LU and complex RHS are separate axes. Final CPU is full paid run cost, including failed or ineffective work.','',
        'For each object and axis the common allowance is min(final actual cost A, final actual cost B). Each endpoint is the minimum original-Phi paid checkpoint at or below that allowance, with latest-index tie break. No GT chooses a checkpoint. Actual checkpoint costs can differ, so unused allowance, next checkpoint costs, bracketing gaps, and incumbent ages are reported. A stopped branch receives no invented work.','',
        'Intervals are 95% percentile intervals from 10,000 independent-object paired bootstrap resamples within cohort/hardware; they are descriptive, unadjusted for multiple comparisons, and may be unstable with small hardware strata. Missing pairs remain explicit.','',
        'A quality gain against a uniform direct method that failed or stopped quickly cannot establish acceleration: the NN physical solver may spend more CPU to attain its better final result. Use both final costs and common-allowance endpoints. No cost-to-GT-target claim is inferred.','',
        'Endpoint status counts: `'+json.dumps(status,sort_keys=True)+'`.','',
        '| Cohort / hardware | Reference | n | final CPU A / B (s) | final NRMSE delta [95%] | common CPU NRMSE delta [95%] |',
        '|---|---|---:|---:|---:|---:|']
    def fmt(x): return 'NA' if x is None else f'{x:.5g}'
    def effect(r):return fmt(r['mean_paired_difference'])+' ['+fmt(r['percentile95_low'])+', '+fmt(r['percentile95_high'])+']'
    for r in fe:
        if r['metric']!='material_nrmse':continue
        key=(r['cohort'],r['hardware'],r['reference'])
        cost=next(x for x in fe if (x['cohort'],x['hardware'],x['reference'])==key and x['metric']=='warm_deployment_cpu_seconds')
        end=next(x for x in ee if (x['cohort'],x['hardware'],x['reference'])==key and x['axis']=='CPU_seconds' and x['metric']=='material_nrmse')
        report.append(f"| {r['cohort']} / {r['hardware']} | {r['reference']} | {r['paired_objects']} | {fmt(cost['method_mean'])} / {fmt(cost['reference_mean'])} | {effect(r)} | {effect(end)} |")
    report += ['','Machine-readable tables retain original Phi, NRMSE, IoU, final residual and final costs where available; the copied endpoint logic reports Phi, NRMSE and IoU. `final_object_pairs.csv` includes individual cost ratios. No cross-hardware absolute CPU pooling is used.','',
        'Artifacts: `final_object_pairs.csv`, `final_paired_effects.csv`, `object_matched_at_most_budget_endpoints.csv`, `object_matched_budget_paired_bootstrap.csv`, `ANALYSIS.json`, `COST.json`.']
    (OUT/'REPORT.md').write_text('\n'.join(report)+'\n')
    save('ANALYSIS.json',dict(status='DESCRIPTIVE_POST_RESULT_ANALYSIS_COMPLETE',contrasts=CONTRASTS,cohort_objects=cohorts,
        total_primary_objects=42,final_pairs=len(final),endpoint_pairs=len(endpoints),endpoint_status=status,
        final_effects=fe,at_most_budget_effects=ee,bootstrap_resamples=10000,
        new_physical_calls=0,new_LU=0,new_complex_RHS=0,raw_truth_access=False,
        selection='Original-Phi incumbent, never GT; no interpolation',
        provenance='Added after seeing primary results; not preregistered confirmation or significance'))
    save('COST.json',dict(cpu_seconds=time.process_time()-START_CPU,wall_seconds=time.perf_counter()-START_WALL,
        new_LU=0,new_complex_RHS=0,new_physical_calls=0,training_calls=0,
        scope='This descriptive extension process, including imports and table generation; no child process'))
    print(json.dumps(dict(status='complete',final_pairs=len(final),endpoint_pairs=len(endpoints),endpoint_status=status,output=str(OUT))))

if __name__=='__main__':main()
