"""Finite Apple/POSIX queue with per-child wait4 accounting, at most four workers.

No psutil, affinity claim, PID matching, worker retries, truth access or daemon.
Popen.poll/wait are NEVER used: wait4 owns child reaping and exact rusage.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,json,platform,subprocess,sys,time,traceback,signal
from pathlib import Path


def reap(process,flags=os.WNOHANG):
    """Return this child's external exit accounting, or None while it runs."""
    pid,status,usage=os.wait4(process.pid,flags)
    if pid==0:return None
    code=os.waitstatus_to_exitcode(status)
    process.returncode=code # Prevent Popen destructor from attempting a second reap.
    return dict(returncode=code,external=dict(child_process_CPU_seconds=usage.ru_utime+usage.ru_stime,user_CPU_seconds=usage.ru_utime,system_CPU_seconds=usage.ru_stime,peak_RSS_bytes=int(usage.ru_maxrss)*(1 if platform.system()=='Darwin' else 1024),CPU_scope='os.wait4 per-child user+system rusage, including interpreter startup/exit and serialization',memory_scope='wait4 per-child peak RSS; Darwin bytes, Linux KiB converted to bytes'))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--legacy-repo',type=Path);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    if not hasattr(os,'wait4'):raise RuntimeError('POSIX wait4 required; this is the Apple/POSIX controller.')
    if not 1<=a.workers<=4:raise ValueError('Worker count must be between1and4')
    root=a.root.resolve();manifest=a.manifest.resolve();legacy=(a.legacy_repo or root.parents[1]).resolve();data=json.loads(manifest.read_text());campaign=data['campaign']
    if data.get('parent_globaldual_patch_ready') is not True:raise RuntimeError('Parent has not frozen the physical GlobalDual runner patch; no dispatch.')
    if not data.get('source_freeze') or str(data['source_freeze']).startswith('PENDING'):raise RuntimeError('Parent source_freeze must be concrete before dispatch.')
    if 'physical_globaldual' not in(root/'src/run_experiment.py').read_text():raise RuntimeError('Required parent physical GlobalDual runner entry is absent.')
    tasks=list(data['tasks'])
    if len({t['task_id'] for t in tasks})!=len(tasks):raise ValueError('Duplicate task IDs')
    for t in tasks:
        if t['campaign']!=campaign:raise ValueError('Campaign mismatch')
        if(root/'RESULTS'/campaign/t['task_id']).exists():raise FileExistsError('Task directory already exists; no overwrite: '+t['task_id'])
        if not(root/t['input_file']).is_file():raise FileNotFoundError(root/t['input_file'])
    queue=root/'queues'/campaign;queue.mkdir(parents=True,exist_ok=False);(queue/'logs').mkdir()
    snapshots=[]
    for f in [manifest]+sorted((root/'src').glob('*.py'))+sorted((root/'runtime').glob('*.py')):
        s=f.stat();snapshots.append(dict(path=str(f),bytes=s.st_size,mtime_ns=s.st_mtime_ns))
    (queue/'SOURCE_SNAPSHOT.json').write_text(json.dumps(dict(source_freeze=data['source_freeze'],files=snapshots,no_SHA_checks=True),indent=2)+'\n')
    env=os.environ.copy();env['PYTHONPATH']=os.pathsep.join([str(legacy),str(legacy/'vendor')]);pending=tasks.copy();running={};jobs=[];events=(queue/'EVENTS.jsonl').open('x');start=time.perf_counter();cpu0=time.process_time();sequence=0;signature=None;fatal=None
    def event(value):events.write(json.dumps(value)+'\n');events.flush()
    def save(status):
        nonlocal sequence,signature
        current=(status,len(pending),len(jobs),tuple(sorted(running)))
        if signature==current:return
        signature=current;sequence+=1
        record=dict(status=status,campaign=campaign,pending=len(pending),pending_task_ids=[t['task_id'] for t in pending],running=[dict(task_id=x['task']['task_id'],pid=pid) for pid,x in running.items()],jobs=jobs,wall_seconds=time.perf_counter()-start,controller_CPU_seconds=time.process_time()-cpu0,finished_child_CPU_seconds=sum(j['external']['child_process_CPU_seconds'] for j in jobs if j.get('external',{}).get('child_process_CPU_seconds') is not None),workers=a.workers,BLAS_threads=1,affinity='NOT_SET_NO_CORE_TYPE_OR_AFFINITY_CLAIM',host='Apple/POSIX local; separate timing from remote i9',python=sys.executable,platform=platform.platform(),source_freeze=data['source_freeze'],truth_access='NONE',error=fatal)
        target=queue/('RECEIPT.json' if status!='RUNNING' else f'STATUS_{sequence:06d}.json')
        with target.open('x') as stream:json.dump(record,stream,indent=2);stream.write('\n')
    def complete(pid,result,interrupted=False):
        run=running.pop(pid);run['log'].close();path=root/'RESULTS'/campaign/run['task']['task_id']/'RECEIPT.json';child=None;receipt_error=None
        if path.is_file():
            try:child=json.loads(path.read_text())
            except Exception as ex:receipt_error=repr(ex)
        row=dict(task_id=run['task']['task_id'],case_id=run['task']['case_id'],method=run['task']['method'],pid=pid,command=run['command'],wall_seconds=time.perf_counter()-run['start'],child_status=child.get('status') if child else None,receipt_error=receipt_error,interrupted_by_controller_cleanup=interrupted,**result)
        jobs.append(row);event(dict(event='END',**row));print(json.dumps(row),flush=True)
    try:
        while pending or running:
            while pending and len(running)<a.workers:
                t=pending.pop(0);command=[sys.executable,str(root/'src/run_experiment.py'),'--root',str(root),'--legacy-repo',str(legacy),'--manifest',str(manifest),'--task-id',t['task_id']];log=(queue/'logs'/(t['task_id']+'.log')).open('x')
                try:process=subprocess.Popen(command,env=env,cwd=root,stdout=log,stderr=subprocess.STDOUT)
                except Exception as ex:
                    log.close();row=dict(task_id=t['task_id'],case_id=t['case_id'],method=t['method'],command=command,returncode=None,child_status='LAUNCH_FAILED',error=repr(ex),external=dict(child_process_CPU_seconds=None,CPU_scope='No child process launched; controller failure cost retained separately'));jobs.append(row);event(dict(event='LAUNCH_FAILED',**row));continue
                running[process.pid]=dict(process=process,task=t,command=command,log=log,start=time.perf_counter());event(dict(event='START',task_id=t['task_id'],pid=process.pid,command=command))
            for pid in list(running):
                result=reap(running[pid]['process'])
                if result is not None:complete(pid,result)
            save('RUNNING')
            if running:time.sleep(.2)
    except BaseException:
        fatal=traceback.format_exc()
        # Clean up ONLY children launched by this controller, never other jobs.
        for x in running.values():
            try:os.kill(x['process'].pid,signal.SIGTERM)
            except ProcessLookupError:pass
        for pid in list(running):complete(pid,reap(running[pid]['process'],flags=0),interrupted=True)
        save('ABORTED_WITH_EVIDENCE');raise
    else:save('COMPLETED' if all(j.get('returncode')==0 for j in jobs) else 'COMPLETED_WITH_ERRORS')
    finally:events.close()
if __name__=='__main__':main()
