"""Finite CPU-only R11 queue; no GT, scoring, daemon or automatic retries."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import argparse,ctypes,json,subprocess,sys,time,traceback
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]


def utc():return datetime.now(timezone.utc).isoformat()


def dump(path,value):
    p=Path(path);tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');tmp.replace(p)


def external_cpu(process):
    if os.name!='nt':return dict(child_process_CPU_seconds=None,CPU_scope='UNKNOWN_EXTERNAL_NONWINDOWS')
    from ctypes import wintypes
    create,exit,kernel,user=[wintypes.FILETIME() for _ in range(4)]
    api=ctypes.WinDLL('kernel32',use_last_error=True).GetProcessTimes
    api.argtypes=[wintypes.HANDLE]+[ctypes.POINTER(wintypes.FILETIME)]*4
    api.restype=wintypes.BOOL
    if not api(int(process._handle),ctypes.byref(create),ctypes.byref(exit),ctypes.byref(kernel),ctypes.byref(user)):
        raise ctypes.WinError(ctypes.get_last_error())
    def seconds(v):return ((v.dwHighDateTime<<32)|v.dwLowDateTime)/1e7
    return dict(child_process_CPU_seconds=seconds(user)+seconds(kernel),
        user_CPU_seconds=seconds(user),kernel_CPU_seconds=seconds(kernel),
        CPU_scope='Windows GetProcessTimes full child lifetime; includes imports and final serialization')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest',type=Path,default=ROOT/'configs/twod_dev_tasks.json')
    ap.add_argument('--config',type=Path,default=ROOT/'configs/twod_dev_solver.json')
    ap.add_argument('--receipt-dir',type=Path,default=ROOT/'delegated/twod/queue')
    ap.add_argument('--python',default=sys.executable)
    ap.add_argument('--cpus',type=int,nargs='+',default=[0,2,4,6])
    ap.add_argument('--vendor-path',type=Path)
    ap.add_argument('--execute',action='store_true')
    args=ap.parse_args()
    if len(set(args.cpus))!=len(args.cpus) or len(args.cpus)>4 or not set(args.cpus).issubset({0,2,4,6,8,10}):raise ValueError('At most four root-authorized independent P-core logical IDs; no HT/E/vector cores')
    # Preserve a venv symlink entry rather than resolve() to its base interpreter.
    python=Path(args.python).expanduser().absolute()
    out=args.receipt_dir.absolute()
    if not out.is_relative_to(ROOT/'delegated/twod'):raise ValueError('Queue output outside owned tree')
    out.mkdir(parents=True,exist_ok=False);(out/'logs').mkdir()
    import psutil
    proc=psutil.Process()
    if os.name=='nt':proc.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    load=psutil.cpu_percent(interval=1,percpu=True)
    probe=dict(utc=utc(),CPU_percent_by_logical=load,
        selected_CPU_percent={str(i):load[i] for i in args.cpus},
        other_python=[dict(pid=p.pid,name=p.name(),affinity=p.cpu_affinity(),priority=p.nice())
            for p in psutil.process_iter() if 'python' in p.name().lower() and p.pid!=os.getpid()],
        priority='BELOW_NORMAL own queue/workers on Windows',
        contention='Shared host; other video workloads may compete; selected cores are not exclusive',GPU=False)
    dump(out/'RESOURCE_PROBE.json',probe)
    cfg=json.loads(args.config.read_text());manifest=json.loads(args.manifest.read_text())
    if cfg['status'] not in ('FROZEN_FOR_INVERSION','FROZEN_BEFORE_GENERATION_AND_CONFIRMATION_GT_ACCESS'):raise ValueError('Inverse contract not frozen')
    dump(out/'FROZEN_MANIFEST.json',manifest);dump(out/'FROZEN_SOLVER.json',cfg)
    tasks=manifest['tasks'];ids=[r['task_id'] for r in tasks]
    if len(set(ids))!=len(ids):raise ValueError('Duplicate task IDs')
    commands=[]
    for row in tasks:
        result=REPO/row['output_dir']
        diagnostic=(row.get('experiment_scope')=='EXPOSED_POSTHOC_STRUCTURAL_DIAGNOSTIC' and
                    result.is_relative_to(ROOT/'RESULTS/right_projection_diagnostic'))
        conditional=(row.get('experiment_scope')=='FRESH_CONDITIONAL_CONFIRMATION' and
                     result.is_relative_to(ROOT/'RESULTS/conditional_confirmation'))
        if not result.is_relative_to(ROOT/'RESULTS'):raise ValueError('Task output outside owned tree')
        if result.exists() or (result.parent/(result.name+'__adapter_inputs')).exists():
            raise FileExistsError('Existing scientific output is never overwritten: '+str(result))
        command=[str(python),str(ROOT/'experiments/twod_run.py'),
            '--manifest',str(out/'FROZEN_MANIFEST.json'),'--task-id',row['task_id'],
            '--config',str(args.config.absolute())]
        if args.vendor_path:command+=['--vendor-path',str(args.vendor_path.absolute())]
        commands.append(dict(task=row,command=command))
    receipt=dict(status='PREPARED_NOT_LAUNCHED',started_utc=utc(),tasks_requested=len(tasks),
        jobs=[],running=[],cpus=args.cpus,threads_per_worker=1,GPU=False,
        truth_access='NONE',scoring='NOT_RUN',automatic_retry=False,
        resource_probe='RESOURCE_PROBE.json',platform=sys.platform,
        cost_policy='Sum disjoint full child GetProcessTimes plus own controller CPU; child internal timers are nested views, never added twice')
    receipt['controller_CPU_seconds']=time.process_time()
    dump(out/'PLAN.json',dict(tasks=commands,config=cfg,resource=probe))
    dump(out/'RECEIPT.json',receipt)
    print(json.dumps(dict(event='RESOURCE_AND_LAUNCH_PLAN_SAVED',receipt=str(out/'RECEIPT.json'),
        execute=args.execute,tasks=len(tasks),resource=probe)),flush=True)
    if not args.execute:return
    # Foreground SSH must remain connected on this Windows host. STOP_DISPATCH
    # is a task-local handoff: active children finish, no other process touched.
    start,cpu=time.perf_counter(),time.process_time();pending=commands[:];running={};jobs=[]
    env=os.environ.copy();env['PYTHONUNBUFFERED']='1';receipt['status']='RUNNING'
    events=(out/'EVENTS.jsonl').open('x')
    def emit(value):events.write(json.dumps(value)+'\n');events.flush();print(json.dumps(value),flush=True)
    def save():
        receipt.update(jobs=jobs,pending=len(pending),
            running=[dict(task_id=j['row']['task']['task_id'],pid=p.pid,cpu=j['cpu']) for p,j in running.items()],
            wall_seconds=time.perf_counter()-start,controller_CPU_seconds=time.process_time(),
            completed_child_CPU_seconds=sum(j.get('child_process_CPU_seconds') or 0 for j in jobs))
        dump(out/'RECEIPT.json',receipt)
    while pending or running:
        stopped=(out/'STOP_DISPATCH').exists()
        occupied={j['cpu'] for j in running.values()}
        for core in args.cpus:
            if not pending or stopped:break
            if core in occupied:continue
            row=pending.pop(0);ident=row['task']['task_id'];command=row['command']+['--cpu',str(core)]
            log=out/'logs'/(ident+'.log');fh=log.open('x')
            info=dict(task_id=ident,case_id=row['task']['case_id'],arm=row['task']['arm'],
                algorithm_label=row['task'].get('material_readout',{}).get('label',row['task']['arm']),
                cpu=core,command=command,log=str(log),output_dir=row['task']['output_dir'],started_utc=utc())
            try:
                p=subprocess.Popen(command,cwd=REPO,env=env,stdout=fh,stderr=subprocess.STDOUT)
                running[p]=dict(row=row,cpu=core,fh=fh,start=time.perf_counter(),info=info)
                emit(dict(event='START',pid=p.pid,**info))
            except Exception as exc:
                fh.write(traceback.format_exc());fh.close()
                info.update(status='LAUNCH_FAILED',error=repr(exc),child_process_CPU_seconds=0.)
                jobs.append(info);emit(dict(event='ERROR',**info))
            save()
        for p,j in list(running.items()):
            code=p.poll()
            if code is None:continue
            info=j['info'];info.update(external_cpu(p));j['fh'].close()
            info.update(returncode=code,finished_utc=utc(),wall_seconds=time.perf_counter()-j['start'])
            path=REPO/info['output_dir']/'RECEIPT.json'
            if path.exists():
                child=json.loads(path.read_text());info.update(child_status=child['status'],
                    child_internal_CPU_seconds=child.get('whole_child_CPU_seconds',child.get('cpu_seconds')),
                    physical_model_cost=child.get('physical_model_cost'),
                    auxiliary_iterations=child.get('auxiliary_iterations'),corrector=child.get('corrector'))
            info['status']='FINISHED' if code==0 and info.get('child_status') in ('SEALED_PUBLIC_PREDICTION','CPU_BUDGET_REACHED') else 'FAILED'
            if info['status']=='FAILED':info['error_log_tail']=log_tail=Path(info['log']).read_text(errors='replace')[-5000:]
            jobs.append(info);del running[p];emit(dict(event='END',**info));save()
        if stopped and not running:break
        if running:time.sleep(.5)
    receipt.update(status=('STOPPED_DISPATCH_CHILDREN_FINISHED' if pending else
        'COMPLETED_WITH_ERRORS' if any(j['status']!='FINISHED' for j in jobs) else 'COMPLETED'),
        not_started_task_ids=[j['task']['task_id'] for j in pending],finished_utc=utc())
    save();events.close()
    return 1 if receipt['status']!='COMPLETED' else 0


if __name__=='__main__':raise SystemExit(main())
