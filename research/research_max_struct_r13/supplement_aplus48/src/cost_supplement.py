"""Ledger new campaigns exactly once, from public terminal receipts; no physics."""
from pathlib import Path
import json,csv,time,collections
C=time.process_time();T=time.perf_counter();S=Path(__file__).resolve().parents[1]
O=S/'analysis/costs';O.mkdir(exist_ok=True)
campaigns=[(S,'aplus48_extension32'),(S/'dynamics_run','aplus48_dynamics8'),(S/'stack_run','aplus_stack112')]
rows=[];summaries=[];operators=[]
for root,campaign in campaigns:
    q=json.loads((root/'queues'/campaign/'RECEIPT.json').read_text());jobs={x['task_id']:x for x in q['jobs']};counter=collections.defaultdict(list)
    before=len(rows)
    for p in sorted((root/'RESULTS'/campaign).glob('*/RECEIPT.json')):
        d=json.loads(p.read_text());tid=d['task']['task_id'];j=jobs[tid];pc=d.get('physical_cost') or{}
        rhs=pc.get('actual_complex_RHS')
        if rhs is None:
            parts=[pc[k] for k in ('source_RHS','adjoint_RHS','tangent_RHS','receiver_transpose_RHS') if k in pc];rhs=sum(parts) if parts else None
        row=dict(campaign=campaign,task_id=tid,case_id=d['task']['case_id'],method=d['task']['method'],role='instrumented_diagnostic' if campaign=='aplus48_dynamics8' else'new_unchanged_completion' if campaign=='aplus48_extension32' else'fixed_version_or_NN_control',status=d['status'],external_process_CPU_seconds=j['external'].get('child_process_CPU_seconds'),internal_CPU_seconds_reference=d.get('whole_process_CPU_seconds'),job_wall_seconds=j.get('wall_seconds'),physical_LU=pc.get('physical_LU',pc.get('LU')),physical_complex_RHS=rhs,material_DOF=d.get('material_DOF'),endpoint_RSS_bytes=d.get('endpoint_RSS_bytes'),RSS_scope='endpoint, not a sampled peak',receipt=str(p),cost_policy='External process CPU once; internal/history/auxiliary timers overlap and are never added')
        rows.append(row)
        aux=d.get('auxiliary_receipt') or{}
        for group,cost in [('physical',pc),('fixed_projector',d.get('projector_cost') or{}),('auxiliary_operations',aux.get('operation_counts') or{})]:
            for key,value in cost.items():
                if isinstance(value,(int,float)) and not isinstance(value,bool):counter[group,key].append(value)
    rr=rows[before:]
    summaries.append(dict(campaign=campaign,expected_tasks=len(jobs),recorded_tasks=len(rr),failed_tasks=sum(x['status']=='FAILED' for x in rr),external_process_CPU_seconds=sum(x['external_process_CPU_seconds'] for x in rr),CPU_unknown=sum(x['external_process_CPU_seconds'] is None for x in rr),physical_LU=sum(x['physical_LU'] for x in rr if x['physical_LU'] is not None),LU_unknown=sum(x['physical_LU'] is None for x in rr),physical_complex_RHS=sum(x['physical_complex_RHS'] for x in rr if x['physical_complex_RHS'] is not None),RHS_unknown=sum(x['physical_complex_RHS'] is None for x in rr),queue_wall_seconds=q['wall_seconds'],queue_controller_CPU_seconds=q['controller_CPU_seconds'],scope='New campaign only; original re-used experiments not charged again'))
    operators.extend(dict(campaign=campaign,group=g,field=k,known_sum=sum(values),receipts_reporting_field=len(values),meaning='Disjoint fields only within declared scope; RHS role fields overlap aggregate and must not be added; missing/non-applicable not imputed0') for (g,k),values in sorted(counter.items()))
def table(name,rs):
    (O/(name+'.json')).write_text(json.dumps(rs,indent=2)+'\n')
    keys=list(dict.fromkeys(k for r in rs for k in r))
    with (O/(name+'.csv')).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys,lineterminator='\n');w.writeheader();w.writerows(rs)
table('NEW_RUN_COSTS',rows);table('NEW_CAMPAIGN_TOTALS',summaries);table('DECLARED_OPERATOR_COUNTERS',operators)
(O/'RECEIPT.json').write_text(json.dumps(dict(status='PUBLIC_RECEIPT_COST_ACCOUNTING_COMPLETE',new_tasks=len(rows),new_objects=0,new_physics=0,new_scoring=0,CPU_seconds=time.process_time()-C,wall_seconds=time.perf_counter()-T,policy='All original re-use excluded from new consumption; offline diagnostics separately measured'),indent=2)+'\n')
print(json.dumps(summaries))
