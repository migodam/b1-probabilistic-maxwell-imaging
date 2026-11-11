"""Finite four-P-core queue. No retries, other-job interference or GT reads."""
import os,sys,time,json,argparse,subprocess,importlib.util
from pathlib import Path
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k]='1'


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--legacy-repo',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--cores',default='0,2,4,6');a=p.parse_args()
    import psutil
    root=a.root.resolve();legacy=a.legacy_repo.resolve();manifest=a.manifest.resolve()
    data=json.loads(manifest.read_text());name=data['campaign'];q=root/'queues'/name;q.mkdir(parents=True,exist_ok=False)
    (q/'logs').mkdir();cores=list(map(int,a.cores.split(',')))
    spec=importlib.util.spec_from_file_location('paid_windows_cpu',legacy/'research/research_max_struct_r11/experiments/twod_queue.py')
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    pending=list(data['tasks']);running={};jobs=[];start=time.perf_counter()
    env=os.environ.copy();env['PYTHONPATH']=os.pathsep.join([str(legacy),str(legacy/'vendor')])
    probe=dict(selected_cores=cores,CPU_percent=psutil.cpu_percent(interval=1,percpu=True),
               other_python=[dict(pid=p.pid,affinity=p.cpu_affinity()) for p in psutil.process_iter()
                             if 'python' in p.name().lower() and p.pid!=os.getpid()])
    (q/'RESOURCE_PROBE.json').write_text(json.dumps(probe,indent=2))
    events=(q/'EVENTS.jsonl').open('x')
    sequence=0;last_signature=None
    def save(status):
        nonlocal sequence,last_signature
        signature=(status,len(pending),len(jobs),tuple(p.pid for p in running))
        if signature==last_signature:return
        last_signature=signature
        record=dict(status=status,campaign=name,pending=len(pending),running=[dict(task_id=x['task']['task_id'],pid=p.pid,cpu=x['cpu']) for p,x in running.items()],
            jobs=jobs,wall_seconds=time.perf_counter()-start,controller_CPU_seconds=time.process_time(),
            finished_child_CPU_seconds=sum(x['external'].get('child_process_CPU_seconds') or 0 for x in jobs))
        # Immutable snapshots avoid Windows rename/delete sharing violations
        # when another read-only observer has opened the previous JSON.
        sequence+=1
        path=q/('RECEIPT.json' if status!='RUNNING' else f'STATUS_{sequence:06d}.json')
        with path.open('x') as stream:json.dump(record,stream,indent=2)
    while pending or running:
        for core in cores:
            if not pending:break
            if core in [x['cpu'] for x in running.values()]:continue
            task=pending.pop(0)
            argv=[str(root/'src/run_experiment.py'),'--root',str(root),'--legacy-repo',str(legacy),
                '--manifest',str(manifest),'--task-id',task['task_id'],'--cpu',str(core)]
            boot='import sys,runpy,psutil;psutil.Process().cpu_affinity(['+str(core)+']);psutil.Process().nice(psutil.NORMAL_PRIORITY_CLASS);sys.argv='+repr(argv)+';runpy.run_path('+repr(argv[0])+',run_name="__main__")'
            log=(q/'logs'/(task['task_id']+'.log')).open('x')
            process=subprocess.Popen([sys.executable,'-c',boot],env=env,cwd=root,stdout=log,stderr=subprocess.STDOUT)
            running[process]=dict(task=task,cpu=core,log=log,start=time.perf_counter())
            events.write(json.dumps(dict(event='START',task_id=task['task_id'],pid=process.pid,cpu=core))+'\n');events.flush()
        for process,run in list(running.items()):
            code=process.poll()
            if code is None:continue
            run['log'].close();rpath=root/'RESULTS'/run['task']['campaign']/run['task']['task_id']/'RECEIPT.json'
            child=json.loads(rpath.read_text()) if rpath.exists() else None
            row=dict(task_id=run['task']['task_id'],case_id=run['task']['case_id'],method=run['task']['method'],cpu=run['cpu'],
                returncode=code,external=old.external_cpu(process),wall_seconds=time.perf_counter()-run['start'],
                child_status=child.get('status') if child else None)
            jobs.append(row);del running[process]
            events.write(json.dumps(dict(event='END',**row))+'\n');events.flush();print(json.dumps(row),flush=True)
        save('RUNNING')
        if running:time.sleep(1)
    save('COMPLETED' if all(x['returncode']==0 for x in jobs) else 'COMPLETED_WITH_ERRORS');events.close()


if __name__=='__main__':main()
