#!/usr/bin/env python3
"""Read R7 JSON receipts only. No array/truth reads, solver calls, or hashing."""
import argparse
import collections
import datetime
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def collect(root):
    records, warnings, consumed = [], [], set()
    def load(p):
        return json.loads((root / p).read_text())
    def add(p, study, d, counts=None, dimension=None, timing=True, note=''):
        key = (str(p), study)
        if key in consumed:
            raise ValueError('duplicate accounting entry: ' + str(key))
        consumed.add(key)
        c = counts if counts is not None else d.get('physical_model_cost', d.get('physical_counts', d.get('counts', {})))
        # This field accumulates dimension once per factorization; it is NOT n.
        counts_clean = {k:v for k,v in c.items() if isinstance(v,(int,float)) and not isinstance(v,bool)}
        wall = d.get('wall_seconds', d.get('wall_s')) if timing else None
        cpu = d.get('cpu_seconds', d.get('cpu_s')) if timing else None
        rss = d.get('peak_process_rss_bytes', d.get('process_peak_rss_bytes', d.get('process_peak_RSS_bytes')))
        records.append(dict(source=str(p), study=study, matrix_dimension=dimension, counts=counts_clean,
            summed_factorization_dimensions=c.get('full_factorization_dimension'),
            task_wall_seconds=wall, cpu_seconds=cpu, peak_process_rss_bytes=rss,
            missing_cost_reason={k:'not measured in this receipt' for k,v in [('wall_seconds',wall),('cpu_seconds',cpu),('peak_process_rss_bytes',rss)] if v is None},
            status=d.get('status'), failed_attempts=d.get('failed_attempts', c.get('failed_calls', d.get('failures', []))), note=note))
    for pattern,study in [('RESULTS/primary/*/RECEIPT.json','primary'),('RESULTS/targeted/*/RECEIPT.json','targeted_confirmation'),('RESULTS/initialization/*/RECEIPT.json','initialization'),('RESULTS/regression/*/RECEIPT.json','historical_digit_regression'),('delegated/aplus/trajectories/*/RECEIPT.json','aplus_development')]:
        for p in sorted(root.glob(pattern)):
            rel=p.relative_to(root);d=load(rel)
            if 'physical_model_cost' not in d:
                warnings.append({'source':str(rel),'reason':'unfinished receipt; no final cumulative physical_model_cost; excluded from numeric totals'})
                continue
            add(rel,study,d,dimension=576,note='Top-level process totals only; nested auxiliary/corrector costs are overlapping. Actual n=24^2 from frozen inversion-grid contract.')
            records[-1]['auxiliary_operations']=d.get('auxiliary_execution_receipt',{}).get('operation_counts',{})
    for study in ['rational','rational_coherent']:
        for p in sorted(root.glob('RESULTS/'+study+'/*/RECEIPT.json')):
            rel=p.relative_to(root);d=load(rel);add(rel,study,d,dimension=d.get('physical_dimension'),note='Cumulative study receipt includes all serial comparator correctors; nested endpoint receipts excluded.')
            records[-1]['small_objective_queries']=d.get('small_objective_queries')
            records[-1]['small_system_dimensions']=[2,8] if study=='rational' else [32]
            if study=='rational':
                actions=p.parent/'SEARCH_RECORDS.json'
                # Per-action query counts may be absent; do not divide a mixed total.
                records[-1]['small_queries_by_dimension']=None
                records[-1]['small_query_reason']='Top-level receipt gives a combined rank2/rank8 query total; no invented allocation.'
                if actions.exists():
                    queries=collections.Counter()
                    for action in json.loads(actions.read_text()):queries[str(action['state_rank'])]+=action['queries']
                    records[-1]['small_queries_by_dimension']=dict(queries)
                    records[-1]['small_query_reason']='Exact per-action queries from SEARCH_RECORDS; path sampling or other operations are not inferred.'
    for pattern,study,dim in [('RESULTS/vector_3d/*/*_receipt.json','vector_initial',375),('RESULTS/vector_continuation/grid216_sphere/*_receipt.json','vector_216_initial',648),('RESULTS/vector_continuation/*/*/RECEIPT.json','vector_continuation',None)]:
        for p in sorted(root.glob(pattern)):
            rel=p.relative_to(root);d=load(rel);n=dim or (648 if '216' in str(rel) else 375)
            add(rel,study,d,dimension=n,timing=False,note='Only new counts; previous_counts excluded. Timing is represented once by whole-process COMPLETE receipt.')
    for rel in ['RESULTS/vector_3d/SIMULATION_COSTS.json','RESULTS/vector_continuation/grid216_sphere/SIMULATION_COSTS.json']:
        if not (root/rel).exists(): continue
        d=load(rel)
        for i,row in enumerate(d if isinstance(d,list) else [d]):
            add(rel+'#'+str(i),'vector_data_generation',row,dimension=648 if '216' in rel else 375,timing=False)
    for rel,key in [('RESULTS/vector_3d/QUALIFICATION.json','new_qualification_counts'),('RESULTS/vector_3d/NONROOT_DERIVATIVES.json','counts')]:
        if (root/rel).exists():
            d=load(rel);add(rel,'vector_derivative_qualification',d,d[key],375,timing='NONROOT' in rel,note='Inherited R1 Mie receipt is reused and excluded from new R7 costs.')
    if (root/'RESULTS/vector_3d/QUALIFICATION.json').exists():
        add('BENCHMARK_3D.md#CM-coordinate-solve','vector_independent_CM_check',{},dict(physical_LU=1,physical_source_RHS=6),375,False,'Separate CM-coordinate solve explicitly reported in BENCHMARK_3D; included in initial process timing, absent from adapter counters.')
    for rel,key,study in [('RESULTS/vector_continuation/MIE_ACTUAL_GRIDS.json','counts','vector_actual_grid_Mie'),('RESULTS/vector_fractional_validation/RECEIPT.json','physical_counts','vector_fractional_Mie')]:
        if (root/rel).exists():
            d=load(rel)
            for i,row in enumerate(d['rows']):
                add(rel+'#rows/'+str(i),study,row,row[key],row['setup']['vector_dimension'],False)
    for rel in ['RESULTS/vector_3d/COMPLETE.json','RESULTS/vector_continuation/COMPLETE.json','RESULTS/vector_fractional_validation/RECEIPT.json']:
        if (root/rel).exists():add(rel,'vector_process_timing',load(rel),{},note='Whole-process timing only; nested per-run/setup timing excluded from totals.')
    for p in sorted(root.glob('RESULTS/fresnel_*/RECEIPT.json')):
        rel=p.relative_to(root);d=load(rel)
        c=collections.Counter(d.get('physical_cost_training',{}));c.update(d.get('physical_cost_heldout_prediction',{}))
        add(rel,'fresnel_inversion',d,dict(c),576,note='Training + heldout physical counters are disjoint; nested auxiliary counters already included in training.')
        records[-1]['auxiliary_operations']=d.get('auxiliary',{}).get('operation_counts',{})
    rel='audit/GENERATION_RECEIPT.json'
    if (root/rel).exists():
        d=load(rel)
        for i,row in enumerate(d['per_case']):add(rel+'#per_case/'+str(i),'generation30',row,dict(LU=row['LU'],complex_RHS=row['complex_RHS']),row['dimension'],False)
        add(rel,'generation_process_timing',d,{},note='Includes setup/generation/serialization; per-case solve timers excluded from total.')
    rel='audit/targeted/GENERATION_RECEIPT.json'
    if (root/rel).exists():
        d=load(rel)
        for i,row in enumerate(d['per_case']):add(rel+'#per_case/'+str(i),'targeted_generation20',row,dict(LU=row['LU'],complex_RHS=row['complex_RHS']),row['dimension'],False)
        add(rel,'targeted_generation_process_timing',d,{},note='Includes all20 generation forwards/setup/serialization; per-case timers excluded from total.')
    rel='delegated/aplus/WIRING_CHECKS.json'
    if (root/rel).exists():
        d=load(rel);add(rel,'aplus_wiring_audit',d,{k:d[k] for k in ['physical_LU','physical_RHS'] if k in d},576,note='Nonphysics wiring work; recorded zero physical solves do not imply zero CPU expenditure.')
    queues=[]
    for p in sorted(root.glob('audit/queue*/RECEIPT.json')):
        d=json.loads(p.read_text());queues.append(dict(source=str(p.relative_to(root)),status=d.get('status'),started_utc=d.get('started_utc'),finished_utc=d.get('finished_utc'),wall_seconds=d.get('wall_seconds',d.get('batch_wall_seconds')),tasks=len(d.get('tasks',[])),failed_tasks=[t for t in d.get('tasks',[]) if t.get('returncode') not in [0,None]]))
    groups={}
    for r in records:
        g=groups.setdefault(r['study'],dict(records=0,counts_by_dimension={},cpu_seconds_observed_sum=0.,task_wall_seconds_observed_sum=0.,missing_cpu_records=0,missing_wall_records=0))
        g['records']+=1
        for field,target,missing in [('cpu_seconds','cpu_seconds_observed_sum','missing_cpu_records'),('task_wall_seconds','task_wall_seconds_observed_sum','missing_wall_records')]:
            if r[field] is None:g[missing]+=1
            else:g[target]+=r[field]
        b=g['counts_by_dimension'].setdefault(str(r['matrix_dimension']),{})
        for k,v in r['counts'].items():b[k]=b.get(k,0)+v
    for g in groups.values():
        for total,missing in [('cpu_seconds_observed_sum','missing_cpu_records'),('task_wall_seconds_observed_sum','missing_wall_records')]:
            if g[missing]==g['records']:g[total]=None
    return dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope='NEW_R7_ONLY_SNAPSHOT',expected_completed_receipts=dict(primary=180,targeted_confirmation=80,initialization=18,historical_digit_regression=9,aplus_development=20,rational=3,rational_coherent=3,fresnel_inversion=6),observed_complete_receipts={s:sum(r['study']==s for r in records) for s in ['primary','targeted_confirmation','initialization','historical_digit_regression','aplus_development','rational','rational_coherent','fresnel_inversion']},groups=groups,records=records,queue_intervals=queues,warnings=warnings,
        global_elapsed_seconds=None,global_elapsed_reason='Concurrent queues overlap; neither task-wall sums nor queue elapsed sums are global elapsed.',energy_joules=None,money=None,unavailable_reason='Energy and money were not metered. Missing per-process RSS/times are null, never zero.',additional_unmetered_costs=[dict(source='RESULTS/fresnel_auxiliary_selection_audit.json',description='Two additional host-objective audit evaluations; LU/RHS/time unrecorded, not zero.'),dict(source='delegated/public_literature/FRESNEL_FRESH_DOWNLOAD.json',description='Retrieval/download costs and prior failed access attempts; complete CPU/wall/RSS not metered.'),dict(description='Literature reading, source/audit scripting, report/figure/packaging work not comprehensively metered.')])


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    data=collect(args.root.resolve());args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data['observed_complete_receipts']))
if __name__=='__main__':main()
