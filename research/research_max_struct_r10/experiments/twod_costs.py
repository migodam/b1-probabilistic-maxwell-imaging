"""Compile paid JSON receipts only; no arrays, truth, scoring or physical calls."""
from pathlib import Path
import argparse,json,time
ROOT=Path(__file__).resolve().parents[1]


def read(path):return json.loads(Path(path).read_text())


def main():
    cpu,wall=time.process_time(),time.perf_counter()
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohort',choices=['baseline','denominator','right'],default='baseline')
    ap.add_argument('--queue-receipt',type=Path)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    tag={'baseline':'REMOTE','denominator':'DENOMINATOR','right':'RIGHT'}[args.cohort]
    engineering={'baseline':'transport','denominator':'denominator_transport','right':'right_transport'}[args.cohort]
    args.queue_receipt=args.queue_receipt or ROOT/('delegated/twod/'+tag+'_QUEUE_RECEIPT.json')
    args.output=args.output or ROOT/('delegated/twod/'+{'baseline':'COST_LEDGER.json','denominator':'DENOMINATOR_COST_LEDGER.json','right':'RIGHT_COST_LEDGER.json'}[args.cohort])
    queue=read(args.queue_receipt);rows=[];views=[];unknown=[]
    def add(ident,path,field,platform,scope):
        p=Path(path)
        if not p.exists():unknown.append(dict(scope=scope,reason='receipt not fetched or not written'));return
        d=read(p)
        if field not in d:unknown.append(dict(scope=scope,reason='CPU field absent'));return
        rows.append(dict(id=ident,artifact=str(p.relative_to(ROOT)),CPU_seconds=d[field],
            platform=platform,additive=True,status=d.get('status','RECORDED'),scope=scope))
    if args.cohort=='baseline':
        add('generation',ROOT/'delegated/twod/GENERATION_RECEIPT.json','CPU_seconds','Mac','24 formal material draws/forwards/noise/serialization; reused by denominator cohort without additional generation')
    add('inverse_package',ROOT/('delegated/twod/'+engineering+'/PACKAGE_RECEIPT.json'),'CPU_seconds','Mac','Inverse-only package and existing threadpoolctl copy')
    transport=ROOT/('delegated/twod/'+engineering+'/TRANSPORT_RECEIPT.json')
    for field,label in [('local_controller_CPU_seconds','transport_controller'),('local_SSH_SCP_CPU_seconds','transport_SSH_SCP')]:
        add(label,transport,field,'Mac','Prelaunch transport; local process component')
    add('remote_plan',ROOT/('delegated/twod/'+tag+'_PLAN_RECEIPT.json'),'controller_CPU_seconds','XINAN','Resource/manifest/fresh-output preflight; no inversion')
    if args.cohort=='right':
        add('busy_core_plan_preserved',ROOT/'delegated/twod/right_transport/BUSY_P0P2P4_PLAN_RECEIPT.json','controller_CPU_seconds','XINAN','First resource plan found P0/2/4 busy; no inverse launched, preserved separately')
    for field,label in [('local_controller_CPU_seconds','live_SSH_controller'),('local_SSH_CPU_seconds','live_SSH')]:
        active=ROOT/('delegated/twod/'+tag+'_ACTIVE.json')
        if active.exists() and field in read(active):add(label,active,field,'Mac','Foreground SSH lifetime; no remote compute double count')
    rows.append(dict(id='remote_queue_controller',artifact=str(args.queue_receipt.relative_to(ROOT)),
        CPU_seconds=queue['controller_CPU_seconds'],platform='XINAN',additive=True,
        status=queue['status'],scope='Own queue process only; child CPU excluded'))
    for j in queue['jobs']:
        value=j.get('child_process_CPU_seconds')
        if value is None:unknown.append(dict(scope=j['task_id'],reason='External child CPU unavailable'));continue
        rows.append(dict(id=j['task_id'],artifact=str(args.queue_receipt.relative_to(ROOT)),CPU_seconds=value,
            platform='XINAN',additive=True,status=j['status'],scope='Full independent child lifetime GetProcessTimes',
            LU=(j.get('physical_model_cost') or {}).get('LU'),
            complex_RHS=(j.get('physical_model_cost') or {}).get('actual_complex_RHS')))
        views.append(dict(task_id=j['task_id'],arm=j['arm'],algorithm_label=j.get('algorithm_label',j['arm']),wall_seconds=j['wall_seconds'],
            internal_CPU_seconds=j.get('child_internal_CPU_seconds'),additive=False,
            scope='Nested internal timing; per-task wall is not queue makespan',
            auxiliary_iterations=j.get('auxiliary_iterations'),corrector=j.get('corrector')))
    status_path=ROOT/('delegated/twod/'+('STATUS_FETCH_COSTS.jsonl' if tag=='REMOTE' else tag+'_STATUS_FETCH_COSTS.jsonl'))
    if status_path.exists():
        fetch=[json.loads(line) for line in status_path.read_text().splitlines() if line]
        rows.append(dict(id='status_fetch_controller',artifact=str(status_path.relative_to(ROOT)),
            CPU_seconds=sum(x['local_controller_CPU_seconds'] for x in fetch),platform='Mac',additive=True,
            status='METADATA_ONLY',scope='Metadata fetch local controllers; no physics'))
        rows.append(dict(id='status_fetch_SCP_known',artifact=str(status_path.relative_to(ROOT)),
            CPU_seconds=sum(x.get('child_transfer_CPU_seconds') or 0 for x in fetch),
            platform='Mac',additive=True,status='METADATA_ONLY',
            scope='Only explicitly measured later SCP children; earlier missing CPU stays unknown'))
        if any(x.get('child_transfer_CPU_seconds') is None for x in fetch):
            unknown.append(dict(scope='metadata fetch SCP children',reason='Earlier fetch CPU unmeasured; do not treat as zero'))
    fetched=ROOT/('delegated/twod/fetch_'+args.cohort)
    if (fetched/'FETCH_RECEIPT.json').exists():
        for field,label in [('local_controller_CPU_seconds','fetch_local_controller'),('local_SSH_SCP_CPU_seconds','fetch_SSH_SCP')]:
            add(label,fetched/'FETCH_RECEIPT.json',field,'Mac','Terminal public-prediction transfer/unpack; no scoring')
        add('remote_archive',fetched/'ARCHIVE_RECEIPT.json','CPU_seconds','XINAN','Terminal prediction/metadata archiving; no physical calls')
    unknown.append(dict(scope='remote unpack/initial SSH resource queries',reason='CPU not measured; engineering only, no physics'))
    if queue.get('running'):
        unknown.append(dict(scope='running child CPU',reason='Not included until process terminal; queue snapshot total is a paid lower bound'))
    result=dict(status='FINAL' if queue['status']=='COMPLETED' else 'PARTIAL_PAID_LOWER_BOUND',
        cohort=args.cohort,requested=queue['tasks_requested'],terminal=len(queue['jobs']),failed=sum(j['status']!='FINISHED' for j in queue['jobs']),
        rows=rows,nested_views=views,unknown_costs=unknown,CPU_seconds_known=sum(r['CPU_seconds'] for r in rows),
        queue_makespan_seconds=queue.get('wall_seconds'),parallel_wall_not_additive=True,
        truth_access='NONE',new_physical_calls=0,monetary_cost=None,
        analysis_CPU_seconds=time.process_time()-cpu,analysis_wall_seconds=time.perf_counter()-wall,
        analysis_cost_scope='This compiler process is separate metadata work; not included in sum of prior paid rows')
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('status','terminal','failed','CPU_seconds_known','queue_makespan_seconds')}))


if __name__=='__main__':main()
