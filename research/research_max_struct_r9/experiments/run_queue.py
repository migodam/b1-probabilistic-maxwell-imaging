"""Finite R9 queue: isolated child processes, immutable manifests, no scoring.
A signal stops new dispatch only; existing children finish without termination.
"""
import argparse
from collections import OrderedDict
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import traceback

ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
DEFAULT_PYTHON="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python"
THREAD_KEYS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')

def utc():return datetime.now(timezone.utc).isoformat()
def resolve(p):
    p=Path(p).expanduser();return (p if p.is_absolute() else REPO/p).resolve()
def dump(p,value):
    p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def stamp(p):
    p=Path(p);s=p.stat();return dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def task_order(rows):
    # First four baseline jobs use distinct case/method pairs; nothing is tuned
    # from outcomes. Then round-robin the remaining public method groups.
    first=[];cases=set();methods=set()
    for row in rows:
        t=row['task'];method=t.get('arm',t.get('algorithm'))
        if row['kind']=='baseline' and t['case_id'] not in cases and method not in methods:
            first.append(row);cases.add(t['case_id']);methods.add(method)
            if len(first)==4:break
    selected={r['task']['task_id'] for r in first};groups=OrderedDict()
    for row in rows:
        if row['task']['task_id'] not in selected:
            key=(row['kind'],row['task'].get('arm',row['task'].get('algorithm')),row['task'].get('dpcs_rule'))
            groups.setdefault(key,[]).append(row)
    order=first[:]
    while any(groups.values()):
        for group in groups.values():
            if group:order.append(group.pop(0))
    return order

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifests',nargs='+',default=[str(ROOT/'configs/baseline_development.json'),str(ROOT/'configs/landscape_development.json'),str(ROOT/'configs/baseline_historical_fixed.json')])
    ap.add_argument('--python',default=DEFAULT_PYTHON)
    ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--r8-repo',default=None,help='Explicit inherited R8 checkout; forwarded only to baseline children')
    ap.add_argument('--receipt-dir',default=str(ROOT/'delegated/development_queue'))
    args=ap.parse_args()
    if not 1<=args.workers<=4:raise ValueError('Authorized local concurrency is1–4')
    # Preserve the venv symlink entry: resolve() would silently select the
    # base interpreter and discard its environment.
    python=Path(args.python).expanduser()
    if not python.is_absolute():python=REPO/python
    python=python.absolute()
    if not python.is_file():raise FileNotFoundError(python)
    out=resolve(args.receipt_dir)
    if not out.is_relative_to(ROOT):raise ValueError('Queue artifacts must stay in R9')
    out.mkdir(parents=True,exist_ok=False);(out/'logs').mkdir()
    rows=[];sources=[]
    for i,name in enumerate(args.manifests):
        source=resolve(name);data=json.loads(source.read_text());copy=out/f'manifest_{i:02d}.json';dump(copy,data);sources.append(stamp(source))
        for task in data['tasks']:
            if 'arm' in task:kind,script='baseline','run_baseline.py'
            elif 'algorithm' in task:kind,script='landscape','run_landscape.py'
            else:raise ValueError('Unknown task contract: '+task['task_id'])
            output=resolve(task['output_dir'])
            if not output.is_relative_to(ROOT):raise ValueError('Task output outside R9')
            command=[str(python),str(ROOT/'experiments'/script),'--manifest',str(copy),'--task-id',task['task_id']]
            if kind=='baseline' and args.r8_repo is not None:
                command+=['--r8-repo',args.r8_repo]
            rows.append(dict(task=task,kind=kind,manifest=str(copy),output_dir=str(output),command=command))
    ids=[r['task']['task_id'] for r in rows]
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate task IDs across manifests')
    pending=task_order(rows);start=time.perf_counter();stopping=[False];running={};finished=[]
    plan=dict(status='FROZEN_FINITE_QUEUE',created_utc=utc(),tasks=len(rows),workers=args.workers,
        executable=str(python),manifests=sources,ordered_tasks=pending,threads_per_child=1,
        truth_access='NONE',scoring='NOT_RUN',automatic_retry=False,hash_checks='NONE')
    dump(out/'PLAN.json',plan)
    receipt=dict(status='RUNNING',started_utc=utc(),tasks=len(rows),workers=args.workers,jobs=finished,
        truth_access='NONE',scoring='NOT_RUN',resource_policy='Only own new children; no other processes inspected or interrupted')
    events=(out/'EVENTS.jsonl').open('x')
    def emit(event):events.write(json.dumps(event,allow_nan=False)+'\n');events.flush();print(json.dumps(event),flush=True)
    def save():
        receipt.update(wall_seconds=time.perf_counter()-start,pending_tasks=len(pending),
            running_tasks=[dict(task_id=v['row']['task']['task_id'],pid=p.pid,started_utc=v['utc'],wall_seconds=time.perf_counter()-v['start']) for p,v in running.items()],jobs=finished)
        dump(out/'RECEIPT.json',receipt)
    def stop_dispatch(signum,frame):
        stopping[0]=True;emit(dict(event='STOP_NEW_DISPATCH_ONLY',utc=utc(),signal=signum,children_terminated=False))
    signal.signal(signal.SIGINT,stop_dispatch);signal.signal(signal.SIGTERM,stop_dispatch)
    env=os.environ.copy();env.update({k:'1' for k in THREAD_KEYS});env['PYTHONUNBUFFERED']='1'
    while pending or running:
        while pending and len(running)<args.workers and not stopping[0]:
            row=pending.pop(0);ident=row['task']['task_id'];log=out/'logs'/(ident+'.log');t=time.perf_counter();when=utc()
            info=dict(task_id=ident,case_id=row['task']['case_id'],kind=row['kind'],command=row['command'],
                output_dir=row['output_dir'],log=str(log),started_utc=when)
            if Path(row['output_dir']).exists():
                info.update(status='NOT_LAUNCHED_EXISTING_OUTPUT',returncode=None,ended_utc=utc(),wall_seconds=0.,error='Existing scientific output preserved; no overwrite or automatic retry')
                finished.append(info);emit(dict(event='ERROR',**info));continue
            fh=log.open('x')
            try:
                p=subprocess.Popen(row['command'],cwd=REPO,env=env,stdout=fh,stderr=subprocess.STDOUT)
                running[p]=dict(row=row,fh=fh,start=t,utc=when,info=info)
                emit(dict(event='START',pid=p.pid,**info))
            except Exception as exc:
                fh.write(traceback.format_exc());fh.close();info.update(status='LAUNCH_FAILED',returncode=None,ended_utc=utc(),wall_seconds=time.perf_counter()-t,error=repr(exc));finished.append(info);emit(dict(event='ERROR',**info))
            save()
        for p,job in list(running.items()):
            code=p.poll()
            if code is None:continue
            job['fh'].close();info=job['info'];info.update(returncode=code,ended_utc=utc(),wall_seconds=time.perf_counter()-job['start'])
            rp=Path(info['output_dir'])/'RECEIPT.json'
            if not rp.exists():rp=Path(info['output_dir']+'__adapter_inputs')/'ADAPTER_PREPARATION.json'
            try:
                rec=json.loads(rp.read_text());info.update(child_receipt=str(rp),child_status=rec.get('status'),child_error=rec.get('error'))
            except Exception as exc:info.update(child_receipt_error=repr(exc))
            bad=code!=0 or str(info.get('child_status','')).startswith('FAILED') or 'child_receipt_error' in info
            info['status']='FAILED' if bad else 'FINISHED'
            if bad:info['error_log_tail']=Path(info['log']).read_text(errors='replace')[-8000:]
            finished.append(info);del running[p];emit(dict(event='END',**info));save()
        if stopping[0] and not running:break
        save()
        if running:time.sleep(.5)
    receipt.update(status='STOPPED_DISPATCH_CHILDREN_FINISHED' if pending else ('COMPLETED_WITH_ERRORS' if any(r['status']!='FINISHED' for r in finished) else 'COMPLETED'),
        finished_utc=utc(),not_started_task_ids=[r['task']['task_id'] for r in pending],failed_tasks=sum(r['status']!='FINISHED' for r in finished))
    save();events.close();return 1 if receipt['status']!='COMPLETED' else 0

if __name__=='__main__':raise SystemExit(main())
