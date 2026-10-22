"""Bounded local process queue; each child pays into the same atomic ledger."""
from pathlib import Path
import sys,os,json,time,subprocess,argparse,shutil
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,write

def schedule(tasks,name,workers=1):
    if workers!=1:raise RuntimeError('memory recovery policy requires one physical worker')
    dest=OUT/'schedules'/name;dest.mkdir(parents=True,exist_ok=True)
    manifest=dest/'TASKS.json'
    if manifest.exists() and json.loads(manifest.read_text())!=tasks:raise RuntimeError('cannot change registered schedule')
    write(manifest,tasks);pending=list(tasks);active=[];done=[];env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONUNBUFFERED='1')
    while pending or active:
        while pending and len(active)<workers:
            args=pending.pop(0);tag='_'.join(args).replace('--','').replace('/','_');log=(dest/(tag+'.log')).open('a')
            operation,parent=args[:2];method=args[args.index('--method')+1];stage=args[args.index('--stage')+1] if '--stage' in args else 'A';suffix=args[args.index('--job-suffix')+1] if '--job-suffix' in args else ''
            job=f'{operation}_{parent}_{method}'+(f'_c{args[args.index("--chain")+1]}' if operation=='reference' else '')+suffix
            closure=OUT/stage/'jobs'/job/'CLOSED.json'
            if closure.exists():
                receipt=json.loads(closure.read_text())
                if receipt['status']=='COMPLETE':
                    log.close();done.append(dict(args=args,exit_code=0,previously_completed=True));continue
                saved=OUT/stage/'anchor_runs'/f'{parent}_{method}{suffix}.json'
                if operation=='anchor' and 'not JSON serializable' in (receipt.get('error') or '') and saved.exists() and saved.with_suffix('.npz').exists():
                    recovery=dict(status='PAID_RESULTS_RECOVERED_WITHOUT_PHYSICAL_REPLAY',original_job_status='ERROR_PRESERVED',saved_report=str(saved),new_physics_RHS=0,reason='Only stdout json.dumps(ndarray) failed after cache/report persistence')
                    write(closure.parent/'POSTPROCESS_RECOVERY.json',recovery);log.close();done.append(dict(args=args,exit_code=0,postprocess_recovery=True));continue
                raise RuntimeError('closed failed job requires explicit diagnosis; no auto physical replay')
            process=subprocess.Popen([sys.executable,str(ROOT/'tools/run_p3l_job.py'),*args],cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
            active.append((args,process,log));print(json.dumps(dict(started=args,PID=process.pid)),flush=True)
        for args,process,log in list(active):
            code=process.poll()
            if code is None:continue
            log.close();active.remove((args,process,log));done.append(dict(args=args,exit_code=code));write(dest/'PROGRESS.json',dict(done=done,running=[dict(args=a,PID=p.pid) for a,p,l in active],pending=pending))
            print(json.dumps(dict(completed=args,exit_code=code)),flush=True)
            if code!=0:
                # Already running independent jobs may finish; no restart/reset.
                for a,p,l in active:p.wait();l.close()
                raise RuntimeError('child failed; preserve ledger and inspect before any recovery')
        if active:time.sleep(.5)
    write(dest/'CLOSED.json',dict(status='COMPLETE',tasks=len(done),all_exit_zero=True))

def main():
    p=argparse.ArgumentParser();p.add_argument('panel',choices=['anchor_pilot','anchors','reference_pilot','reference_revision_pilot','geometry_pilot','mode_reference_pilot','references','crosschecks']);p.add_argument('--workers',type=int,default=1);a=p.parse_args()
    if a.panel=='anchor_pilot':tasks=[['anchor',f'DEV_{i:03d}','--method',method,'--polish'] for i in range(4) for method in ['full_TRF','coarse_TRF']]
    elif a.panel=='anchors':
        selection=json.loads((OUT/'A/BASELINE_SELECTION.json').read_text());method=selection['method'];tasks=[['anchor',f'DEV_{i:03d}','--method',method,'--polish'] for i in range(4,32)]
    else:
        parents=json.loads((OUT/'A/DEEP_REFERENCE_PARENTS.json').read_text())['parents']
        if a.panel=='reference_pilot':parents=parents[:2]
        elif a.panel=='reference_revision_pilot':parents=parents[2:3]
        elif a.panel in ['geometry_pilot','mode_reference_pilot']:parents=parents[2:3]
        elif a.panel=='references':parents=parents[2:]
        if a.panel=='mode_reference_pilot':tasks=[['reference',parent,'--method',method,'--chain',str(c)] for parent in parents for method in ['mala_mode','pcn_mode'] for c in range(4)]
        elif a.panel=='crosschecks':tasks=[['reference',parent,'--method',method,'--chain',str(c),'--warmup','256','--production','512'] for parent in parents[:4] for method in ['pcn','ess'] for c in range(2)]
        elif a.panel=='geometry_pilot':tasks=[['reference',parent,'--method','mala_decoder','--chain',str(c),'--stage','C'] for parent in parents for c in range(4)]
        else:tasks=[['reference',parent,'--method','mala_retracted' if a.panel in ['references','reference_revision_pilot'] else 'mala','--chain',str(c)] for parent in parents for c in range(4)]
    schedule(tasks,a.panel,workers=a.workers)

if __name__=='__main__':main()
