"""Finite public-only VAL queue, gated on completion of this M4 preparation.

Only explicit request entries run; two BLAS1 partitions; no truth or score API.
Absolute deadline is fixed and only this controller's own children may be stopped.
"""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import json,os,platform,re,signal,subprocess,sys,time

STUDY=Path(__file__).resolve().parents[2]
REQUEST=STUDY/'audit/VAL_EVALUATION_REQUEST.json'
RECEIPT=STUDY/'audit/VAL_EVALUATION_CONTROLLER.json'
PREPARATION=STUDY/'audit/LOCAL_PREPARATION_PROGRESS.json'
DEADLINE_TEXT='2026-11-07T10:44:27Z'
DEADLINE=datetime(2026,11,7,10,44,27,tzinfo=timezone.utc).timestamp()
POLL_SECONDS=15
TERMINAL_PREPARATION={'LOCAL_PHYSICS_COMPLETE','LOCAL_COMPLETE_UPLOADED','LOCAL_COMPLETE_UPLOAD_INCOMPLETE'}
RUN_RE=re.compile(r'^[A-Za-z0-9_-]+$')
INITIAL_RUNS=['v_physics_seed7_cascade2','j_full_seed7_cascade2','j_null_seed7_cascade2','g_weak_seed7_cascade2','g_dark_seed7_cascade2']

def write_json(path,value):
    tmp=path.with_suffix('.controller.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    tmp.replace(path)

def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))

def public_val_ids():
    records=read_json(STUDY/'data/manifests/online.json')['records']
    allowed={'obs_id','parent_id','split','noise','public_path'}
    if any(set(r)!=allowed for r in records):raise ValueError('Unexpected online manifest fields')
    ids=[r['obs_id'] for r in records if r['split']=='val']
    if len(ids)!=156 or len(set(ids))!=156:raise ValueError('Require all156 public VAL observations')
    return ids

def request_runs():
    request=read_json(REQUEST)
    runs=request.get('runs')
    if not isinstance(runs,list) or any(not isinstance(r,str) or not RUN_RE.fullmatch(r) for r in runs):
        raise ValueError('Request runs must be an explicit list of safe run names')
    if len(set(runs))!=len(runs):raise ValueError('Duplicate run in request')
    if not isinstance(request.get('closed',False),bool):raise ValueError('closed must be boolean')
    return request,runs

def preparation_finished():
    try:p=read_json(PREPARATION)
    except (OSError,ValueError):return False,None
    ready=p.get('status') in TERMINAL_PREPARATION and p.get('local_process_pid') is None
    return ready,p

def pair_count(directory,ids):
    return sum((directory/f'{obs}.npz').is_file() and (directory/f'{obs}.json').is_file() for obs in ids)

def prediction_ready(run,ids):
    directory=STUDY/'predictions'/run
    n=pair_count(directory,ids)
    if n!=156:return False,n
    # The producer writes JSON last; require valid complete metadata before launch.
    try:
        for obs in ids:read_json(directory/f'{obs}.json')
    except (OSError,ValueError):return False,n
    return True,n

def canonical_evaluation_cost(run,ids):
    directory=STUDY/'evaluation'/run
    costs=Counter();statuses=Counter();accounting=Counter();seconds=0.;completed=[]
    for obs in ids:
        p=directory/f'{obs}.json'
        if not p.is_file() or not p.with_suffix('.npz').is_file():continue
        r=read_json(p);completed.append(obs);statuses[r.get('status','UNKNOWN')]+=1
        accounting[r.get('accounting_status','UNKNOWN')]+=1
        for k,v in r.get('costcounts',{}).items():
            if isinstance(v,(int,float)):costs[k]+=v
        seconds+=float(r.get('wall_seconds',0.))
    return dict(completed=len(completed),expected=156,costcounts=dict(costs),sum_case_charged_wall_seconds=seconds,status_histogram=dict(statuses),accounting_status_histogram=dict(accounting),accounting_scope='canonical completed VAL case receipts once; existing runner charged totals may include proposal preparation; do not add that preparation again',truth_access='NONE',score_executed=False)

def stop_own_children(children):
    # A fixed absolute deadline must release only CPU processes started here.
    for child in children:
        if child.poll() is None:
            try:os.killpg(child.pid,signal.SIGTERM)
            except ProcessLookupError:pass
    end=time.monotonic()+3
    while any(c.poll() is None for c in children) and time.monotonic()<end:time.sleep(.1)
    for child in children:
        if child.poll() is None:
            try:os.killpg(child.pid,signal.SIGKILL)
            except ProcessLookupError:pass
    for child in children:
        try:child.wait(timeout=1)
        except subprocess.TimeoutExpired:pass

def main():
    ids=public_val_ids();started=time.time();completed={};history=[];active=None
    if RECEIPT.exists():
        previous=read_json(RECEIPT)
        history=list(previous.get('history',[]))
    state=dict(schema='gn-nn-r4-finite-val-public-controller-v1',pid=os.getpid(),started_unix=started,absolute_deadline_utc=DEADLINE_TEXT,max_numeric_children=2,BLAS_threads_per_child=1,hardware=platform.platform(),truth_access='NONE',score_executed=False,request_path=str(REQUEST),completed_runs=completed,history=history)
    def publish(status,**more):
        state.update(status=status,updated_unix=time.time(),active_run=active,**more)
        write_json(RECEIPT,state)
        print(json.dumps(dict(status=status,active_run=active,completed_runs=list(completed),**more)),flush=True)
    while time.time()<DEADLINE:
        try:request,runs=request_runs()
        except (OSError,ValueError) as exc:
            publish('WAIT_VALID_REQUEST',request_error=str(exc));time.sleep(min(POLL_SECONDS,max(0.,DEADLINE-time.time())));continue
        # Evaluate canonical completed receipts once per requested run, even on resume.
        for run in runs:
            if run not in completed and pair_count(STUDY/'evaluation'/run,ids)==156:
                completed[run]=canonical_evaluation_cost(run,ids)
        if request.get('closed',False) and all(r in completed for r in runs):
            publish('CLOSED_ALL_REQUESTED_RUNS_COMPLETED',requested_runs=runs);return 0
        ready,prep=preparation_finished()
        if not ready:
            publish('WAIT_M4_PREPARATION_END',requested_runs=runs,local_preparation_status=prep.get('status') if prep else None)
            time.sleep(min(POLL_SECONDS,max(0.,DEADLINE-time.time())));continue
        selected=None;prediction_counts={}
        for run in runs:
            if run in completed:continue
            ready,n=prediction_ready(run,ids);prediction_counts[run]=n
            if ready:selected=run;break
        if selected is None:
            publish('WAIT_EXPLICIT_RUN_PREDICTIONS',requested_runs=runs,prediction_pair_counts=prediction_counts)
            time.sleep(min(POLL_SECONDS,max(0.,DEADLINE-time.time())));continue
        # Re-read the owner-controlled request immediately before launching.
        try:latest,latest_runs=request_runs()
        except (OSError,ValueError):continue
        if selected not in latest_runs:continue
        active=selected;before=canonical_evaluation_cost(active,ids);children=[];handles=[];commands=[]
        env=os.environ.copy()
        for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):env[name]='1'
        tick=time.monotonic();launch_unix=time.time()
        for part in (0,1):
            if time.time()>=DEADLINE:break
            command=[sys.executable,str(STUDY/'experiments/run_evaluation.py'),'--runs',active,'--splits','val','--partition-count','2','--partition-index',str(part)]
            log_path=STUDY/f'delegated/data/val_evaluation_{active}_P{part}_{time.time_ns()}.log'
            handle=log_path.open('w');handles.append(handle)
            child=subprocess.Popen(command,stdout=handle,stderr=subprocess.STDOUT,env=env,start_new_session=True)
            children.append(child);commands.append(dict(partition_index=part,pid=child.pid,command=command,log=str(log_path)))
        publish('RUNNING_TWO_VAL_PARTITIONS',children=commands,preexisting_evaluations=before['completed'])
        expired=False
        while any(c.poll() is None for c in children):
            if time.time()>=DEADLINE:
                expired=True;stop_own_children(children);break
            time.sleep(min(5.,max(0.,DEADLINE-time.time())))
        for handle in handles:handle.close()
        after=canonical_evaluation_cost(active,ids)
        entry=dict(run=active,launch_unix=launch_unix,controller_batch_wall_seconds=time.monotonic()-tick,children=commands,returncodes=[c.returncode for c in children],preexisting_evaluations=before['completed'],new_completed_evaluations=after['completed']-before['completed'],canonical_cost=after,deadline_stop=expired,inflight_at_stop_cost='UNKNOWN_NOT_ZERO' if expired else 'NONE_KNOWN')
        history.append(entry)
        if expired:
            publish('ABSOLUTE_DEADLINE_STOPPED_OWN_CHILDREN',last_run=entry);return 0
        if len(children)!=2 or any(c.returncode for c in children) or after['completed']!=156:
            publish('HOLD_INCOMPLETE_OR_FAILED_RUN_NO_AUTOMATIC_RETRY',last_run=entry);return 1
        completed[active]=after
        publish('RUN_COMPLETED_156_AND_ACCOUNTED',last_run=entry)
        active=None
    publish('ABSOLUTE_DEADLINE_EXIT_NO_ACTIVE_CHILDREN');return 0

if __name__=='__main__':
    raise SystemExit(main())
