"""Finite predeclared TEST host partition; public inputs only, never scoring."""
from pathlib import Path
from collections import deque,Counter
from concurrent.futures import ThreadPoolExecutor,wait,FIRST_COMPLETED
from datetime import datetime,timezone
import argparse,json,os,platform,re,signal,subprocess,sys,threading,time
STUDY=Path(__file__).resolve().parents[1]
DEADLINE_TEXT='2026-11-07T10:44:27Z'
DEADLINE=datetime(2026,11,7,10,44,27,tzinfo=timezone.utc).timestamp()
INDEPENDENT={'strict_gn45_anchor':'','weakreg45_anchor':'_weakreg1e-5','discrepancy45_anchor':'_discrepancy','discrepancy_allscales45_anchor':'_discrepancy_allscales'}
SAFE=re.compile(r'^[A-Za-z0-9_-]+$')

def read_json(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def write_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');tmp.replace(path)
def complete_pair(base):return base.with_suffix('.npz').is_file() and base.with_suffix('.json').is_file()
def pair_host(base):return read_json(base.with_suffix('.json')).get('execution',{}).get('execution_partition') if complete_pair(base) else None
def online_records():
    records=read_json(STUDY/'data/manifests/online.json')['records']
    if any(set(r)!={'obs_id','parent_id','split','noise','public_path'} for r in records):raise ValueError('Unexpected public manifest fields')
    return records

def load_partition():
    allocation=read_json(STUDY/'audit/TEST_RESOURCE_PARTITION.json')
    records=online_records();test=[r for r in records if r['split']=='test'];parents=list(dict.fromkeys(r['parent_id'] for r in test))
    if len(parents)!=108 or len(test)!=216:raise ValueError('Require all108 TEST mothers /216 observations')
    indices={r['obs_id']:i for i,r in enumerate(records) if r['split']=='test'}
    if allocation['global_online_indices']!=indices:raise ValueError('Global online index map differs from public manifest')
    for host in ('local','remote'):
        wanted=[p for i,p in enumerate(parents) if (i%3==0)==(host=='local')]
        obs=[r['obs_id'] for r in test if r['parent_id'] in set(wanted)]
        if allocation[host]['parent_ids']!=wanted or allocation[host]['obs_ids']!=obs:raise ValueError('Predeclared deterministic host allocation differs')
        if any(sum(r['parent_id']==p for r in test)!=2 for p in wanted):raise ValueError('Mother clean/noisy pairing incomplete')
    return allocation,test

def load_contract():
    freeze=read_json(STUDY/'audit/METHOD_FREEZE.json');final=read_json(STUDY/'audit/FINAL_RUNS.json')
    plans=final['plans'];runs=final['runs'];nn=final['nn_runs']
    if not isinstance(plans,list) or not isinstance(runs,list) or not isinstance(nn,list):raise ValueError('Frozen plans/runs/nn_runs must be lists')
    if runs!=[p['run'] for p in plans] or len(set(runs))!=len(runs):raise ValueError('FINAL_RUNS plans/run order or uniqueness mismatch')
    if any(not isinstance(r,str) or not SAFE.fullmatch(r) for r in runs):raise ValueError('Unsafe run name')
    if any(p['kind'] not in ('NN','continue_trf','independent_baseline','cie') for p in plans):raise ValueError('Unknown frozen method kind')
    if nn!=[p['run'] for p in plans if p['kind']=='NN']:raise ValueError('nn_runs differs from NN plans')
    if freeze['test_runs']!=runs or freeze['nn_runs']!=nn:raise ValueError('METHOD_FREEZE and FINAL_RUNS disagree')
    if freeze['test_observations']!=216 or freeze['test_mothers']!=108:raise ValueError('Frozen TEST denominator differs')
    if not isinstance(freeze['baselines'],list) or any(r not in runs or r in nn for r in freeze['baselines']):raise ValueError('Invalid frozen baseline names')
    for p in plans:
        if p['kind']=='continue_trf' and p['run']!='continue_trf':raise ValueError('Unknown continuation entry')
        if p['kind']=='independent_baseline' and p['run'] not in INDEPENDENT:raise ValueError('Unknown independent baseline')
        if p['kind']=='cie' and (p['method'] not in ('cie_single','cie_dual') or p['run']!=p['method']+'_cold55'):raise ValueError('Unknown CIE plan')
    allocation,test=load_partition()
    return freeze,final,allocation,test

def require_host(host,workers):
    limit=2 if host=='local' else 4
    if not 1<=workers<=limit:raise ValueError('At most2 local /4 remote BLAS1 workers')
    expected='darwin' if host=='local' else 'win32'
    if sys.platform!=expected:raise ValueError(f'{host} allocation must execute on its declared host platform{expected}')

def baseline744_complete():
    records=[r for r in online_records() if r['split'] in ('val','test')]
    return sum(complete_pair(STUDY/'anchors'/(r['obs_id']+suffix)) for r in records for suffix in ('_discrepancy','_weakreg1e-5'))==744

def stop_own(children):
    for child in children:
        if child.poll() is not None:continue
        try:
            if sys.platform=='win32':subprocess.run(['taskkill','/PID',str(child.pid),'/T','/F'],capture_output=True,timeout=8)
            else:os.killpg(child.pid,signal.SIGTERM)
        except (ProcessLookupError,subprocess.TimeoutExpired):pass
    end=time.monotonic()+2
    while any(c.poll() is None for c in children) and time.monotonic()<end:time.sleep(.05)
    for child in children:
        if child.poll() is None:
            try:
                if sys.platform=='win32':child.kill()
                else:os.killpg(child.pid,signal.SIGKILL)
            except ProcessLookupError:pass

def run_fifo(tasks,workers,receipt_path,*,label):
    """Finite ready-observation FIFO. Refill on every completion, no batch barrier."""
    pending=deque(tasks);results=[];active={};lock=threading.Lock();began=time.perf_counter();expired=False
    env=os.environ.copy()
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):env[name]='1'
    def one(task):
        tick=time.perf_counter();row={k:v for k,v in task.items() if k!='command'}
        try:
            path=Path(task['log']);path.parent.mkdir(parents=True,exist_ok=True)
            with path.open('a') as log:
                with lock:
                    if time.time()>=DEADLINE:return dict(row,status='NOT_RUN',reason='ABSOLUTE_DEADLINE_BEFORE_LAUNCH',new_physics_counts={})
                    options={'creationflags':subprocess.CREATE_NEW_PROCESS_GROUP} if sys.platform=='win32' else {'start_new_session':True}
                    child=subprocess.Popen(task['command'],stdout=log,stderr=subprocess.STDOUT,env=env,**options);active[child.pid]=child
                rc=child.wait()
                with lock:active.pop(child.pid,None)
            return dict(row,status='CHILD_COMPLETED' if rc==0 else 'CHILD_FAILED_OR_INTERRUPTED',pid=child.pid,returncode=rc,subprocess_wall_seconds=time.perf_counter()-tick,inflight_cost_status='UNKNOWN_NOT_ZERO' if rc else 'READ_CANONICAL_CASE_RECEIPTS')
        except Exception as exc:
            if 'child' in locals():
                stop_own([child])
                with lock:active.pop(child.pid,None)
            return dict(row,status='LAUNCH_OR_WAIT_ERROR',error=repr(exc),subprocess_wall_seconds=time.perf_counter()-tick,inflight_cost_status='UNKNOWN_NOT_ZERO')
    def snapshot(status):
        write_json(receipt_path,dict(status=status,label=label,controller_pid=os.getpid(),host_system=platform.system(),machine=platform.machine(),workers=workers,BLAS_threads=1,absolute_deadline_utc=DEADLINE_TEXT,planned=len(tasks),completed_child_jobs=len(results),remaining_unlaunched=len(pending),queue_wall_seconds=time.perf_counter()-began,results=results,truth_access='NONE',score_executed=False))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures={}
        for _ in range(min(workers,len(pending))):
            task=pending.popleft();futures[pool.submit(one,task)]=task
        snapshot('FINITE_FIFO_RUNNING')
        while futures:
            remaining=DEADLINE-time.time()
            if remaining<=0:
                expired=True
                with lock:children=list(active.values())
                stop_own(children)
                for f in list(futures):
                    results.append(f.result());futures.pop(f)
                break
            done,_=wait(futures,timeout=min(1.,remaining),return_when=FIRST_COMPLETED)
            for future in done:
                futures.pop(future);results.append(future.result())
                if pending and time.time()<DEADLINE:
                    task=pending.popleft();futures[pool.submit(one,task)]=task
                snapshot('FINITE_FIFO_RUNNING')
        for task in pending:results.append(dict(obs_id=task.get('obs_id'),run=task.get('run'),status='NOT_RUN',reason='ABSOLUTE_DEADLINE_UNLAUNCHED',new_physics_counts={}))
    snapshot('ABSOLUTE_DEADLINE_OWN_CHILDREN_STOPPED' if expired else 'FINITE_FIFO_CLOSED')
    return results

def canonical_costs(runs,obs_ids,host=None):
    rows=[];totals=Counter();seconds=0.
    for obs in obs_ids:
        for run in runs:
            base=STUDY/'evaluation'/run/obs
            if not complete_pair(base):rows.append(dict(run=run,obs_id=obs,status='NOT_RUN',complete_and_assigned=False));continue
            receipt=read_json(base.with_suffix('.json'));counts=receipt.get('costcounts',{})
            for k,v in counts.items():
                if isinstance(v,(int,float)):totals[k]+=v
            wall=float(receipt.get('wall_seconds',0.));seconds+=wall
            assigned=host is None or pair_host(base)==host
            rows.append(dict(run=run,obs_id=obs,status=receipt.get('status') if assigned else 'HOST_MISMATCH_HOLD',original_case_status=receipt.get('status'),complete_and_assigned=assigned,costcounts=counts,charged_case_wall_seconds=wall,execution=receipt.get('execution',{})))
    return dict(rows=rows,canonical_completed=sum(r['complete_and_assigned'] for r in rows),planned=len(runs)*len(obs_ids),costcounts=dict(totals),sum_charged_case_wall_seconds=seconds,scope='one canonical saved case per method and obs; some charged totals include preparation, do not add it again; queue/setup/transport overhead is separate')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--execution-partition',required=True,choices=['local','remote']);parser.add_argument('--workers',type=int);parser.add_argument('--execute',action='store_true');args=parser.parse_args()
    started=time.perf_counter();freeze,final,allocation,test=load_contract();host=args.execution_partition;workers=args.workers or (2 if host=='local' else 4);require_host(host,workers)
    obs_ids=allocation[host]['obs_ids'];plans=final['plans'];cie=[p for p in plans if p['kind']=='cie'];post=[p['run'] for p in plans if p['kind']!='cie']
    stamp=time.time_ns();directory=STUDY/'delegated/data/test_queues';directory.mkdir(parents=True,exist_ok=True)
    receipt_path=directory/f'PARTITION_{host}_{stamp}.json'
    plan=dict(status='NOT_RUN',partition=host,workers=workers,BLAS_threads=1,absolute_deadline_utc=DEADLINE_TEXT,mothers=len(allocation[host]['parent_ids']),observations=len(obs_ids),global_mothers=108,global_observations=216,runs=final['runs'],obs_ids=obs_ids,setup_wall_seconds=time.perf_counter()-started,method_order='NN/continue/independent order rotated by full online-manifest index; CIE precedes postprocessing',truth_access='NONE',score_executed=False,resource_rule='baseline744 complete and own numerical P cores free before remote launch; local previous2CPU work ended',scientific_budgets='unchanged existing per-method budgets; fixed absolute deadline may leave NOT_RUN')
    if not args.execute:print(json.dumps(plan,indent=2));return
    if time.time()>=DEADLINE:write_json(receipt_path,dict(plan,status='NOT_RUN_ABSOLUTE_DEADLINE'));return
    wrong_host=[dict(obs_id=obs,run=run,reported_host=pair_host(STUDY/'evaluation'/run/obs)) for obs in obs_ids for run in final['runs'] if complete_pair(STUDY/'evaluation'/run/obs) and pair_host(STUDY/'evaluation'/run/obs)!=host]
    if wrong_host:write_json(receipt_path,dict(plan,status='HOLD_EXISTING_HOST_PROVENANCE_PRESERVED',wrong_host=wrong_host));return
    if not baseline744_complete():raise RuntimeError('Wait until all744 preexisting traditional anchor baseline pairs complete')
    if host=='local':
        prep=read_json(STUDY/'audit/LOCAL_PREPARATION_PROGRESS.json')
        if prep.get('status') not in ('LOCAL_PHYSICS_COMPLETE','LOCAL_COMPLETE_UPLOADED','LOCAL_COMPLETE_UPLOAD_INCOMPLETE') or prep.get('local_process_pid') is not None:raise RuntimeError('M4 preparation still occupies its2CPU allocation')
        val_path=STUDY/'audit/VAL_EVALUATION_CONTROLLER.json'
        if val_path.exists() and read_json(val_path).get('status')=='RUNNING_TWO_VAL_PARTITIONS':raise RuntimeError('M4 VAL correction still occupies its2CPU allocation')
    write_json(receipt_path,dict(plan,status='EXECUTION_REQUESTED'))
    phase_results=[]
    for p in cie:
        command=[sys.executable,str(STUDY/'experiments/run_cie_reference.py'),'--split','test','--scope','all108','--methods',p['method'],'--workers',str(workers),'--execution-partition',host,'--execute','--test-methods-frozen']
        task=dict(run=p['run'],command=command,log=str(directory/f'CIE_{host}_{p["run"]}_{stamp}.log'))
        # One queue-controller child owns at most workers CIE numeric children.
        result=run_fifo([task],1,directory/f'CIE_CONTROLLER_{host}_{p["run"]}_{stamp}.json',label='CIE_FIRST_'+host);phase_results+=result
        if any(r['status'] not in ('CHILD_COMPLETED',) for r in result):write_json(receipt_path,dict(plan,status='CIE_PHASE_FAILED_OR_DEADLINE',phases=phase_results,costs=canonical_costs(final['runs'],obs_ids,host)));return
    missing=[]
    for obs in obs_ids:
        for p in plans:
            if p['kind']=='NN' and not complete_pair(STUDY/'predictions'/p['run']/obs):missing.append(dict(obs_id=obs,run=p['run'],input='PREDICTION_PAIR'))
            if p['kind']!='cie' and not complete_pair(STUDY/'anchors'/obs):missing.append(dict(obs_id=obs,run=p['run'],input='COMMON_ANCHOR_PAIR'))
            if p['kind']=='independent_baseline' and not complete_pair(STUDY/'anchors'/(obs+INDEPENDENT[p['run']])):missing.append(dict(obs_id=obs,run=p['run'],input='INDEPENDENT_ANCHOR_PAIR'))
    if missing:write_json(receipt_path,dict(plan,status='HOLD_MISSING_PUBLIC_INPUTS_NO_CASE_SELECTION',missing=missing,phases=phase_results));return
    tasks=[]
    if post:
        for obs in obs_ids:
            tasks.append(dict(obs_id=obs,command=[sys.executable,str(STUDY/'experiments/run_evaluation.py'),'--runs',*post,'--splits','test','--obs-id',obs,'--execution-partition',host],log=str(directory/f'POST_{host}_{obs}_{stamp}.log')))
    before=canonical_costs(post,obs_ids,host);queue_start=time.perf_counter()
    results=run_fifo(tasks,workers,directory/f'POST_FIFO_{host}_{stamp}.json',label='TEST_POST_'+host)
    costs=canonical_costs(final['runs'],obs_ids,host)
    status='PARTITION_COMPLETE' if costs['canonical_completed']==costs['planned'] else 'PARTITION_INCOMPLETE_NOT_RUN_OR_FAILED'
    write_json(receipt_path,dict(plan,status=status,phases=phase_results,post_jobs=results,post_queue_wall_seconds=time.perf_counter()-queue_start,total_controller_wall_seconds=time.perf_counter()-started,preexisting_post_cases=before['canonical_completed'],costs=costs,speed_claim='none; hardware, setup, queue, transport and actual case costs reported separately',interrupted_cost_policy='preserve every exception and UNKNOWN_NOT_ZERO; no automatic physics retries'))
    print(receipt_path,flush=True)
if __name__=='__main__':main()
