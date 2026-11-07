"""Conditional final4 cost metadata only; prefix charged once in production."""
from pathlib import Path
import json,time
ROOT=Path(__file__).resolve().parents[1]


def main():
    cpu=time.process_time();folder=ROOT/'delegated/twod/conditional_confirmation'
    d=json.loads((folder/'RECEIPT.json').read_text());rows=[];views=[];unknown=[]
    def add(ident,source,field,platform,scope):
        p=ROOT/source
        if not p.exists():unknown.append(dict(scope=scope,source=source,reason='Receipt absent'));return
        z=json.loads(p.read_text());v=z.get(field)
        rows.append(dict(id=ident,artifact=source,CPU_seconds=v,platform=platform,scope=scope,status=z.get('status','RECORDED'),additive=True))
        if v is None:unknown.append(dict(source=source,reason='CPU unmeasured, not zero'))
    add('new_four_generation','delegated/twod/conditional_confirmation/GENERATION_RECEIPT.json','CPU_seconds','Mac','Independent new4 data generation once, not another charge for primary24')
    add('inverse_package','delegated/twod/conditional_transport/PACKAGE_RECEIPT.json','CPU_seconds','Mac','Conditional inverse-only bundle')
    for field,label in [('local_controller_CPU_seconds','transport_controller'),('local_SSH_SCP_CPU_seconds','transport_SSH_SCP')]:
        add(label,'delegated/twod/conditional_transport/TRANSPORT_RECEIPT.json',field,'Mac','Preparation transfer only')
    add('remote_plan','delegated/twod/CONDITIONAL_PLAN_RECEIPT.json','controller_CPU_seconds','XINAN','Own resource/manifest preflight')
    for tag in ['CONDITIONAL','CONDITIONAL_REPAIR']:
        for field,label in [('local_controller_CPU_seconds','controller'),('local_SSH_CPU_seconds','SSH')]:
            add(tag+'_'+label,'delegated/twod/'+tag+'_ACTIVE.json',field,'Mac','Own foreground connection; remote science charged separately')
    rows.append(dict(id='pipeline_controller_all_segments',artifact='delegated/twod/conditional_confirmation/RECEIPT.json',CPU_seconds=d['pipeline_controller_CPU_seconds'],platform='XINAN',scope='Sum disjoint original and controlled startup-repair controller process timers; queues/children excluded',status=d['status'],additive=True))
    for phase in d['queues']:
        rows.append(dict(id=phase['name']+'_queue_process',artifact=phase['receipt'],CPU_seconds=phase['external_queue_process_cost']['child_process_CPU_seconds'],platform='XINAN',scope='Queue process full GetProcessTimes; its worker children excluded',status=phase['status'],additive=True))
        for j in phase['jobs']:
            c=j.get('physical_model_cost') or {}
            startup='All R10 outputs must stay inside' in j.get('error_log_tail','') and not c
            rows.append(dict(id=phase['name']+':'+j['task_id'],artifact=phase['receipt'],CPU_seconds=j.get('child_process_CPU_seconds'),platform='XINAN',scope='Whole worker attempt; shared prefix charged once, including preserved zero-physics failed startup',status=j['status'],additive=True,LU=c.get('LU',0 if startup else None),complex_RHS=c.get('actual_complex_RHS',0 if startup else None),engineering_startup=startup))
    views.extend(d.get('deployments',[]))
    for field,label in [('local_controller_CPU_seconds','fetch_controller'),('local_SSH_SCP_CPU_seconds','fetch_SSH_SCP')]:
        add(label,'delegated/twod/fetch_conditional/FETCH_RECEIPT.json',field,'Mac','Terminal prediction transfer/unpack')
    add('remote_archive','delegated/twod/fetch_conditional/ARCHIVE_RECEIPT.json','CPU_seconds','XINAN','Terminal public prediction/metadata archive, zero physics')
    unknown.append(dict(scope='Status queries, resume-code upload, uninstrumented initial helpers',reason='Small engineering CPU not fully measured, never set to zero'))
    result=dict(status='FINAL' if d['status']=='ALL_ACTUAL_CONDITIONAL_PREDICTIONS_SEALED' else 'PARTIAL',cohort='conditional_confirmation',
        rows=rows,deployment_views=views,nested_views=views,unknown_costs=unknown,
        CPU_seconds_known=sum(r['CPU_seconds'] for r in rows if r['CPU_seconds'] is not None),
        actual_inverse_runs=d.get('actual_inverse_runs'),actual_worker_process_attempts=d.get('actual_worker_process_attempts'),
        prefix_runs=d.get('prefix_runs'),triggered_cases=d.get('triggered_cases'),fallback_runs=d.get('fallback_runs'),
        engineering_startup_failures=d.get('engineering_startup_failures'),truth_access='NONE',new_physical_calls=0,
        policy='Production unique phase/worker attempts once; same-task failed startup+valid first solve retained separately; deployment methods each contain prefix+own fallback but overlap and must not be summed',
        analysis_CPU_seconds=time.process_time()-cpu)
    (ROOT/'delegated/twod/CONDITIONAL_COST_LEDGER.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('status','CPU_seconds_known','actual_inverse_runs','prefix_runs','triggered_cases','fallback_runs')}))


if __name__=='__main__':main()
