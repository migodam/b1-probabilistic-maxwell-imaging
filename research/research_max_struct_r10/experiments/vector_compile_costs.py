"""R10 vector campaign metadata-only accounting; no NPZ/GT/solver calls.

wait4/GetProcessTimes whole children supersede nested main/perarm timers.
Active counters are cumulative lower bounds and replaced by their final row.
Device events and norm/Fourier diagnostics are never added as host CPU.
"""
import time
T=time.perf_counter();C=time.process_time()
from pathlib import Path
import json,subprocess,re,argparse
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'RESULTS/vector';D=ROOT/'delegated/vector'


def load(path):
    path=Path(path)
    if 'sealed_truth' in path.parts or path.suffix not in ('.json','.jsonl'):
        raise ValueError('Only explicit nontruth metadata allowed')
    if not path.exists():return None
    text=path.read_text()
    return [json.loads(x) for x in text.splitlines() if x.strip()] if path.suffix=='.jsonl' else json.loads(text)


def cpu_time(text):
    values=text.strip().split(':')
    try:return sum(float(x)*60**i for i,x in enumerate(reversed(values)))
    except ValueError:return None


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--final',action='store_true');args=parser.parse_args()
    rows=[];views=[];gaps=[];running=[]
    def charge(identity,platform,cpu,wall,source,scope,counts=None):
        rows.append(dict(identity=identity,platform=platform,CPU_seconds=cpu,wall_seconds=wall,
            source=str(Path(source).relative_to(ROOT)),scope=scope,counts=counts,additive=True,
            wall_scope='Own interval, not summed as campaign makespan',operation=identity,
            included_in=None,terminal_ledger_status='MEASURED_KNOWN_LOWER_BOUND' if cpu is not None else 'CPU_UNKNOWN'))
    # Completed local scientific children: exact whole process authority.
    for queue in ('cpu_first_queue','mac_extension_queue','mac_cross_queue'):
        folder=D/queue;receipts=load(folder/'PROCESS_RECEIPTS.json') or []
        completed={row['label'] for row in receipts}
        for row in receipts:charge(queue+'-'+row['label'],'Mac',row['CPU_seconds'],row['wall_seconds'],folder/'PROCESS_RECEIPTS.json',row['scope'])
        parent=load(folder/'CONTROLLER_COST.json')
        if parent:
            charge(queue+'-controller','Mac',parent['CPU_seconds'],parent['wall_seconds'],folder/'CONTROLLER_COST.json','Own controller only; children already charged independently')
            if parent['status']=='FAILED' and load(folder/'RESOURCE_HANDOFF_GPU_WINDOW.json'):
                views.append(dict(source=str((folder/'RESOURCE_HANDOFF_GPU_WINDOW.json').relative_to(ROOT)),status='CANCELLED_WAIT_ONLY_RESOURCE_HANDOFF_NOT_ALGORITHM_FAILURE'))
        active=load(folder/'ACTIVE.json')
        if active and active['label'] not in completed:
            p=subprocess.run(['/bin/ps','-p',str(active['pid']),'-o','time='],capture_output=True,text=True)
            value=cpu_time(p.stdout)
            charge(queue+'-'+active['label']+'-active','Mac',value,None,folder/'ACTIVE.json','Active process ps CPU lower bound; replaces no completed row and must be superseded at completion')
            running.append(queue+'-'+active['label'])
    folder=D/'v04_mac_primary';p=folder/'PROCESS_COST.json';d=load(p)
    if d:
        charge('v04-Mac-new-primary','Mac',d['CPU_seconds'],d['wall_seconds'],p,d['scope'])
        charge('v04-Mac-controller','Mac',d['controller_CPU_seconds'],None,p,'Own wait4 launcher only')
    elif load(folder/'ACTIVE.json'):
        active=load(folder/'ACTIVE.json');snapshot=subprocess.run(['/bin/ps','-p',str(active['pid']),'-o','time='],capture_output=True,text=True)
        charge('v04-Mac-new-primary-active','Mac',cpu_time(snapshot.stdout),None,folder/'ACTIVE.json','New whole primary active CPU lower bound, replaces no completed row')
        running.append('v04-Mac-new-primary')
    # One necessary completed CPU backend benchmark, use wait4 not its inner CPU.
    p=D/'CPU_PROBE_ATTEMPT2_PROCESS.json';d=load(p)
    if d:
        charge('CPU-backend-probe','Mac',d['whole_child_CPU_seconds'],d['whole_child_wall_seconds'],p,'Whole public uniform VG child; internal setup/forward/adjoint not re-added',dict(LU=1,primal_RHS=60,adjoint_RHS=60,dimension=5184))
        charge('CPU-probe-controller','Mac',d['local_controller_CPU_seconds'],None,p,'Own benchmark launcher only')
    p=OUT/'cpu_backend_probe_failed_attempt1/RECEIPT.json';d=load(p)
    if d:charge('CPU-probe-failed-receipt','Mac',None,None,p,'Paid uniform forward completed before receipt-assembly TypeError; CPU unavailable, not zero',d['counts_inferred_from_reached_traceback'])
    for name in ('operator_spectral','representation_rank'):
        p=OUT/name/'RECEIPT.json';d=load(p)
        if d:charge(name,'Mac',d['CPU_seconds'],d['wall_seconds'],p,'Own background-only analysis process; no material GT or physical solves',dict(physical_RHS=0))
    # Local connected transport processes exclude all remote scientific CPU.
    for folder in sorted(D.glob('*transport*'))+([D/'remote_first'] if (D/'remote_first').exists() else []):
        for filename in ('TRANSPORT_COST.json','TRANSPORT_CONTROLLER_COST.json'):
            p=folder/filename;d=load(p)
            if d:
                charge(folder.name+'-parent','Mac',d['local_controller_CPU_seconds'],d['wall_seconds'],p,'Own local transport parent, remote science excluded')
                charge(folder.name+'-SSH-SCP','Mac',d['local_transport_child_CPU_seconds'],None,p,'Own local SSH/SCP children, remote science excluded')
    # Historical stopped pre-CUDA probe; not a scientific algorithm failure.
    p=D/'remote_first/STOPPED_RESOURCE_HANDOFF_COST.json';d=load(p)
    if d:
        charge('GPU-wait-stopped-main','XINAN',d['CPU_seconds'],d['wall_seconds'],p,'Exact retained-handle CPU after stopping own wait-only probe; zero CUDA context/physical RHS')
        charge('GPU-wait-query-children','XINAN',d['GPU_gate_helper_CPU_seconds'],None,p,'Completed own nvidia-smi child CPU aggregate; not inside main CPU')
        charge('GPU-wait-stop-helper','XINAN',d['own_stop_helper_CPU_seconds'],d['own_stop_helper_wall_seconds'],p,'Own resource-handoff helper only, clock-resolution lower bound')
        gaps.append(dict(source=str(p.relative_to(ROOT)),CPU_seconds=None,reason='Possible in-flight query tail and unmeasured manual status/upload helpers remain unknown'))
    # Read-only Windows producer snapshot. Completed rows replace active counters.
    p=D/'REMOTE_STATUS_SNAPSHOT.json';remote=load(p)
    if remote:
        done={(x['source'].replace('\\','/').split('/')[2],x['label']) for x in remote['completed_process_rows']}
        for i,row in enumerate(remote['completed_process_rows']):
            charge('Windows-child-'+str(i),'XINAN',row['CPU_seconds'],row['wall_seconds'],p,'Whole Windows child GetProcessTimes: '+row['source']+' / '+row['label'])
        for i,row in enumerate(remote['controller_rows']):
            charge('Windows-controller-'+str(i),'XINAN',row['CPU_seconds'],row['wall_seconds'],p,'Own Windows finite parent only: '+row['source'])
            helper=row.get('GPU_gate_helper_CPU_seconds')
            if helper is not None:charge('Windows-controller-query-'+str(i),'XINAN',helper,None,p,'Own parent query subprocesses only')
        for i,row in enumerate(remote['active_process_snapshots']):
            key=(row['source'].replace('\\','/').split('/')[2],row['label'])
            if key in done:continue
            charge('Windows-active-'+str(i),'XINAN',row['CPU_seconds'],row['wall_seconds'],p,'Cumulative active main CPU lower bound; later final whole child replaces this row')
            running.append(row['label'])
        for i,row in enumerate(remote['process_cost_headers']):
            # Whole scientific child main CPU already represented above. Only
            # its separately measured query child aggregate is another charge.
            helper=row.get('GPU_gate_helper_CPU_seconds')
            if helper is not None:charge('Windows-query-children-'+str(i),'XINAN',helper,None,p,'Own completed scientific process query children, main CPU already charged')
        views.extend(remote['endpoints'])
    p=D/'STATUS_HELPERS.jsonl';d=load(p)
    if d:
        for i,row in enumerate(d):charge('Windows-status-helper-'+str(i),'XINAN',row['CPU_seconds'],row['wall_seconds'],p,row['scope'])
    # Stopped FFT worker had no final process header. Its completed query rows
    # are disjoint child CPU, while STOP GetProcessTimes duplicates the queue row.
    p=ROOT/'RESULTS/vector_engineering_attempts/v08_CUDA_FFT_GMRES_attempt1/case_v08_gpu_process/GPU_BOUNDARIES.jsonl'
    d=load(p)
    if d:
        measured=[x['query_helper_CPU_seconds'] for x in d if x.get('query_helper_CPU_seconds') is not None]
        charge('old-FFT-query-children','XINAN',sum(measured),None,p,'Completed own nvidia-smi children; stopped main CPU already charged by final queue row')
        gaps.append(dict(source=str(p.relative_to(ROOT)),CPU_seconds=None,reason='Any interrupted unrecorded query child tail is unknown; recorded rows only'))
    for name,helperkey in (('gpu_backend_probe','gate_helper_CPU_seconds'),('gpu_dense_backend_probe','query_child_CPU_seconds')):
        p=OUT/name/'RECEIPT.json';d=load(p)
        if d:
            charge(name,'XINAN',d['CPU_seconds'],d['wall_seconds'],p,'Own necessary public m0 backend assessment main process; all timing subintervals nested',d.get('counts'))
            charge(name+'-query-children','XINAN',d.get(helperkey),None,p,'Own query subprocess CPU, separate from main process')
    for p in sorted((D/'remote_repository_costs').glob('*.json')):
        d=load(p)
        charge('remote-repository-helper-'+p.stem,'XINAN',d.get('CPU_seconds'),d.get('wall_seconds'),p,'Own extraction/byte packaging helper; measured zero means clock-resolution lower bound')
    for name,key,wallkey in (('FINAL_HARVEST_COST.json','CPU_seconds','wall_seconds'),('SLOW_QUEUE_MIGRATION.json','helper_CPU_seconds','helper_wall_seconds')):
        p=D/name;d=load(p)
        if d:charge(name,'XINAN',d[key],d[wallkey],p,'Own byte harvest or stop/archive helper only; stopped child CPU counted in queue, not repeated here')
    for tag in ('v04_Windows_CPU_dense_attempt1','v08_CUDA_FFT_GMRES_attempt1'):
        p=ROOT/'RESULTS/vector_engineering_attempts'/tag/'STOP_AND_MIGRATION_RECEIPT.json';d=load(p)
        if d:views.append(dict(source=str(p.relative_to(ROOT)),CPU_seconds=d['CPU_seconds'],additive=False,included_in='Windows-child final queue row',scope='Same terminated child GetProcessTimes, not a second charge'))
    # Per-arm scientific receipts are views of already paid whole processes.
    local_sealed=[]
    protocol=load(ROOT/'configs/vector_all18.json')
    for task in protocol['tasks']:
        folder=OUT/task['task_id'];receipt=load(folder/'RECEIPT.json')
        if receipt:
            views.append(dict(task_id=task['task_id'],CPU_seconds=receipt['CPU_seconds'],wall_seconds=receipt['wall_seconds'],counts=receipt['counts'],additive=False,
                scope='Nested perarm view inside whole scientific worker; not another charge'))
        if (folder/'SEALED.json').exists():local_sealed.append(task['task_id'])
    allsealed=sorted(set(local_sealed+(remote['sealed'] if remote else [])))
    history_path=OUT/'COST_COMPILER_INVOCATIONS.json';history=load(history_path) or []
    history.append(dict(CPU_seconds=time.process_time()-C,wall_seconds=time.perf_counter()-T,platform='Mac'))
    history_path.write_text(json.dumps(history,indent=2)+'\n')
    for i,row in enumerate(history):charge('vector-compiler-'+str(i),'Mac',row['CPU_seconds'],row['wall_seconds'],history_path,'Own metadata compiler invocation only')
    totals={platform:sum(x['CPU_seconds'] for x in rows if x['platform']==platform and x['CPU_seconds'] is not None) for platform in ('Mac','XINAN')}
    gaps.append(dict(CPU_seconds=None,reason='Manual source editing, ordinary inspection, initial unmetered fetch/stop commands, failed compiler invocation, and public endpoint comparison analysis CPU were not individually measured; not zero. GPU events are not host CPU or active device time.'))
    if args.final and (running or len(allsealed)!=18):raise RuntimeError('Final accounting requires all18 terminal with no active primary process')
    result=dict(status='FINAL_KNOWN_LOWER_BOUND' if args.final else ('IN_PROGRESS' if running or len(allsealed)<18 else 'ALL18_SEALED_ACCOUNTING_KNOWN_LOWER_BOUND'),rows=rows,
        by_platform_CPU_lower_bounds=totals,running=running,sealed=allsealed,expected=18,
        nonadditive_views=views,unknowns=gaps,GT_read=False,scoring='NOT_RUN_ROOT_AUTHORIZATION_REQUIRED',
        counter_scope='Never pool LU/RHS of field dimension5184 vs17496 or scalar2D; auxiliary/current matvecs distinct from physical RHS',
        CPU_total_is_known_lower_bound=True,wall_not_summed_as_makespan=True)
    result['finalized_for_parent_compiler']=args.final
    result['scope']='Own R10 vector18 primary, generation, probes, stopped/replaced attempts, background analysis, measured transport/helpers; root scoring/figures are separate'
    result['GPU_event_scope']='Perarm stream timestamp intervals include host/IO/pause gaps; never added to CPU or asserted active kernel time'
    (OUT/'ACCOUNTING_CURRENT.json').write_text(json.dumps(result,indent=2)+'\n')
    if args.final:(OUT/'ACCOUNTING_FINAL.json').write_text(json.dumps(result,indent=2)+'\n')
    lines=['# Vector R10 live cost ledger','',f"Status: {result['status']}; sealed {len(allsealed)}/18; GT scoring NOT_RUN.",'',
        f"Known CPU lower bounds: Mac {totals['Mac']:.6f}s; XINAN {totals['XINAN']:.6f}s. Active cumulative counters are replaced, never added to final same-process CPU.",'',
        'Whole wait4/GetProcessTimes children are additive. Setup/arm/forward/adjoint/device-event views are nested. CUDA events include host/pause gaps and are not active kernel time. Different vector/scalar dimensions are not reduced to one speed metric.','',
        'Historical GPU waiting and engineering failures remain charged/unknown. Publication and old R9 costs excluded. Raw source ledgers are unchanged.','',
        'Running: '+', '.join(running),'','Unmeasured helper/failed-receipt CPU remains null, not zero.']
    (D/'COST_LEDGER_CURRENT.md').write_text('\n'.join(lines)+'\n')
    if args.final:(D/'COST_LEDGER_FINAL.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(status=result['status'],sealed=len(allsealed),by_platform=totals,running=running)))


if __name__=='__main__':main()
