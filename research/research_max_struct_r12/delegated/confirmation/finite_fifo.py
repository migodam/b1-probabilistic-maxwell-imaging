"""Finite eight-P-core queue, R11 external process accounting, public inputs.

No automatic retries, GT, scheduler service, or changes to unrelated processes.
Sequential manifests finish all active children before starting the next stage.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import argparse,importlib.util,json,sys,time,subprocess,traceback
ROOT=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest',type=Path,nargs='+',required=True)
    ap.add_argument('--legacy-repo',type=Path,required=True)
    ap.add_argument('--receipt-dir',type=Path,required=True)
    ap.add_argument('--python',default=sys.executable)
    ap.add_argument('--cpus',type=int,nargs='+',default=[0,2,4,6,8,10,12,14])
    ap.add_argument('--coordinated-queues',type=Path,nargs='*',default=[])
    ap.add_argument('--execute',action='store_true');a=ap.parse_args()
    if len(set(a.cpus))!=len(a.cpus) or not 1<=len(a.cpus)<=8 or not set(a.cpus)<=set(range(0,16,2)):
        raise ValueError('Only root-allocated independent P-core IDs0,2,...14')
    legacy=a.legacy_repo.resolve()
    source=legacy/'research/research_max_struct_r11/experiments/twod_queue.py'
    spec=importlib.util.spec_from_file_location('r11_existing_queue_accounting',source)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    dump,utc,external_cpu=old.dump,old.utc,old.external_cpu
    out=a.receipt_dir.resolve()
    if not out.is_relative_to(ROOT/'queues'):raise ValueError('Queue receipts must be confirmation-owned')
    for predecessor in a.coordinated_queues:
        if not predecessor.resolve().is_relative_to(ROOT/'queues'):raise ValueError('Only own finite predecessor queue leases may be coordinated')
    out.mkdir(parents=True,exist_ok=False);(out/'logs').mkdir()
    allids=set();stages=[]
    for path in a.manifest:
        manifest=json.loads(path.read_text());tasks=manifest['tasks']
        for task in tasks:
            if task['task_id'] in allids:raise ValueError('Duplicate task ID')
            allids.add(task['task_id'])
            result=Path(task['output_dir']).resolve()
            if not result.is_relative_to(ROOT/'RESULTS'):raise ValueError('Output outside R12 owned tree')
            if result.exists():raise FileExistsError('Existing paid output never overwritten: '+str(result))
        stages.append(dict(stage=manifest['stage'],manifest=str(path.resolve()),tasks=tasks))
    for stage in stages:dump(out/('FROZEN_'+stage['stage'].upper()+'_MANIFEST.json'),stage)
    import psutil
    probe=dict(utc=utc(),CPU_percent_by_logical=psutil.cpu_percent(interval=1,percpu=True),
        python_processes=[dict(pid=p.pid,name=p.name(),command=p.cmdline(),affinity=p.cpu_affinity())
            for p in psutil.process_iter() if 'python' in p.name().lower() and p.pid!=os.getpid()],
        selected_cores=a.cpus,allocation='Root explicitly released previous basin allocation and authorized this eight-core finite queue',
        scope='Read-only observation; bind/nice only own new children; no unrelated process edits')
    dump(out/'RESOURCE_PROBE.json',probe)
    env=os.environ.copy();env.update(R12_LEGACY_REPO=str(legacy),PYTHONUNBUFFERED='1')
    env['PYTHONPATH']=os.pathsep.join([str(ROOT.parent),str(legacy),str(legacy/'vendor')])
    receipt=dict(status='PREPARED_NOT_LAUNCHED',started_utc=utc(),jobs=[],running=[],stages=[s['stage'] for s in stages],
        requested=len(allids),cpus=a.cpus,truth_access='NONE',scoring='NOT_RUN',automatic_retry=False,
        coordinated_predecessor_receipts=[str(p.resolve()) for p in a.coordinated_queues],
        cost_policy='Existing R11 external GetProcessTimes full child + own controller CPU; never add nested child timers')
    dump(out/'RECEIPT.json',receipt)
    if not a.execute:
        print(json.dumps(dict(status=receipt['status'],tasks=len(allids),receipt=str(out/'RECEIPT.json'))));return
    events=(out/'EVENTS.jsonl').open('x');jobs=[];running={};start=time.perf_counter()
    def emit(row):
        events.write(json.dumps(row)+'\n');events.flush();print(json.dumps(row),flush=True)
    def save(stage,pending):
        receipt.update(current_stage=stage,pending=len(pending),jobs=jobs,
            running=[dict(task_id=j['task']['task_id'],pid=p.pid,cpu=j['cpu']) for p,j in running.items()],
            wall_seconds=time.perf_counter()-start,controller_CPU_seconds=time.process_time(),
            completed_child_CPU_seconds=sum(j.get('child_process_CPU_seconds') or 0 for j in jobs))
        dump(out/'RECEIPT.json',receipt)
    receipt['status']='RUNNING'
    def predecessor_running():
        return [j for p in a.coordinated_queues for j in json.loads(p.read_text()).get('running',[])]
    for stageindex,stage in enumerate(stages):
        pending=list(stage['tasks']);frozen=out/('FROZEN_'+stage['stage'].upper()+'_MANIFEST.json')
        while pending or running or (stageindex==0 and predecessor_running()):
            stopped=(out/'STOP_DISPATCH').exists()
            occupied={j['cpu'] for j in running.values()}|{j['cpu'] for j in predecessor_running()}
            for core in a.cpus:
                if not pending or stopped:break
                if core in occupied:continue
                task=pending.pop(0);ident=task['task_id'];runner=Path(task['runner'])
                if task['executor_kind']=='baseline':
                    argv=[str(runner),'--manifest',str(frozen),'--config',task['config'],'--task-id',ident,'--vendor-path',str(legacy/'vendor'),'--cpu',str(core)]
                else:
                    argv=[str(runner),'--public-dir',task['public_dir'],'--case',task['case_id'],'--method',task['method'],'--output',task['output_dir'],'--cpu',str(core)]
                    if task['method']=='nn_direct':argv+=['--weights',task['weights']]
                # Affinity is set before importing NumPy/SciPy/Torch or the scientific runner.
                bootstrap=('import os,sys,runpy,psutil;psutil.Process().cpu_affinity(['+str(core)+']);'
                    'sys.argv='+repr(argv)+';runpy.run_path('+repr(str(runner))+',run_name="__main__")')
                command=[str(Path(a.python).absolute()),'-c',bootstrap]
                log=out/'logs'/(ident+'.log');handle=log.open('x')
                info=dict(task_id=ident,case_id=task['case_id'],method=task['algorithm_label'],cpu=core,
                    command=command,output_dir=task['output_dir'],log=str(log),stage=stage['stage'],started_utc=utc())
                try:
                    process=subprocess.Popen(command,cwd=ROOT.parent,env=env,stdout=handle,stderr=subprocess.STDOUT)
                    running[process]=dict(task=task,cpu=core,handle=handle,info=info,start=time.perf_counter())
                    emit(dict(event='START',pid=process.pid,**info))
                except Exception as exc:
                    handle.write(traceback.format_exc());handle.close();info.update(status='LAUNCH_FAILED',error=repr(exc),child_process_CPU_seconds=0.)
                    dump(out/'TERMINAL_JOBS'/(ident+'.json'),info);jobs.append(info);emit(dict(event='END',**info))
                save(stage['stage'],pending)
            for process,job in list(running.items()):
                code=process.poll()
                if code is None:continue
                info=job['info'];info.update(external_cpu(process));job['handle'].close()
                info.update(returncode=code,finished_utc=utc(),wall_seconds=time.perf_counter()-job['start'])
                path=Path(info['output_dir'])/'RECEIPT.json'
                if path.exists():
                    try:
                        child=json.loads(path.read_text());info.update(child_status=child.get('status'),
                            child_internal_CPU_seconds=child.get('whole_child_CPU_seconds',child.get('cpu_seconds')),
                            physical_model_cost=child.get('physical_model_cost',child.get('physics_cost')),
                            objective=child.get('cap_physical_objective',child.get('objective')))
                    except Exception as exc:info['receipt_read_error']=repr(exc)
                info['status']='FINISHED' if code==0 and info.get('child_status') in ('SEALED_PUBLIC_PREDICTION','CPU_BUDGET_REACHED','COMPLETE') else 'FAILED'
                if info['status']=='FAILED':
                    info['error_log_tail']=Path(info['log']).read_text(errors='replace')[-6000:]
                    # An actual failed child always gets a terminal controller receipt,
                    # even if failure happened before its own scientific receipt existed.
                dump(out/'TERMINAL_JOBS'/(info['task_id']+'.json'),info);jobs.append(info)
                del running[process];emit(dict(event='END',**info));save(stage['stage'],pending)
            if stopped and not running:break
            if running or predecessor_running():time.sleep(.5)
        if pending:
            receipt['not_started_task_ids']=[r['task_id'] for r in pending]+[r['task_id'] for s in stages[stageindex+1:] for r in s['tasks']]
            receipt['status']='STOPPED_DISPATCH_OWN_CHILDREN_FINISHED';save(stage['stage'],pending);events.close();return 1
        emit(dict(event='STAGE_FINISHED',stage=stage['stage'],jobs=len(jobs)))
    receipt.update(status='COMPLETED_WITH_ERRORS' if any(j['status']!='FINISHED' for j in jobs) else 'COMPLETED',finished_utc=utc())
    save(stages[-1]['stage'],[]);events.close();return 1 if receipt['status']!='COMPLETED' else 0

if __name__=='__main__':raise SystemExit(main())
