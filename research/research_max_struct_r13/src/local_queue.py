"""Finite sequential Apple-CPU control with whole-child resource accounting."""
from pathlib import Path
import argparse,json,subprocess,sys,time,resource


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);a=p.parse_args();r=a.root.resolve()
    d=json.loads(a.manifest.read_text());q=r/'queues'/d['campaign'];q.mkdir(parents=True,exist_ok=False)
    jobs=[];wall=time.perf_counter()
    for task in d['tasks']:
        before=resource.getrusage(resource.RUSAGE_CHILDREN);started=time.perf_counter()
        command=[sys.executable,str(r/'src/run_experiment.py'),'--root',str(r),'--legacy-repo',str(r.parents[1]),'--manifest',str(a.manifest.resolve()),'--task-id',task['task_id']]
        with (q/(task['task_id']+'.log')).open('x') as stream:code=subprocess.call(command,stdout=stream,stderr=subprocess.STDOUT)
        after=resource.getrusage(resource.RUSAGE_CHILDREN);cpu=(after.ru_utime-before.ru_utime)+(after.ru_stime-before.ru_stime)
        path=r/'RESULTS'/task['campaign']/task['task_id']/'RECEIPT.json';child=json.loads(path.read_text()) if path.exists() else None
        row=dict(task_id=task['task_id'],case_id=task['case_id'],method=task['method'],returncode=code,child_status=child.get('status') if child else None,
            wall_seconds=time.perf_counter()-started,external=dict(child_process_CPU_seconds=cpu,CPU_scope='RUSAGE_CHILDREN delta for one sequential child; includes interpreter startup/exit/export'))
        jobs.append(row);print(json.dumps(row),flush=True)
    receipt=dict(status='COMPLETED' if all(j['returncode']==0 for j in jobs) else 'COMPLETED_WITH_ERRORS',pending=0,running=[],jobs=jobs,
        finished_child_CPU_seconds=sum(j['external']['child_process_CPU_seconds'] for j in jobs),wall_seconds=time.perf_counter()-wall,
        controller_CPU_seconds=time.process_time(),host='Apple local CPU; not same-host timing against remote i9 development',truth_access='NONE')
    (q/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':main()
