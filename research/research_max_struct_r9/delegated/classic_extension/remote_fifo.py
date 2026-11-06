"""Finite Windows CPU FIFO. Pin one allocated P-core before scientific imports.

Run in a connected SSH session. No services, GPU use, or other process control.
Exact external CPU uses Windows GetProcessTimes on retained Popen handle.
"""
import os
import time
START, CPU_START=time.time(),time.process_time()
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[key]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
import ctypes
from ctypes import wintypes
import json
from pathlib import Path
import subprocess
import sys
import psutil

REPO=Path(__file__).resolve().parents[4]
ROOT=REPO/'research/research_max_struct_r9'
OUT=ROOT/'RESULTS/classic_extension_xinan'
CORES=[0,2,4,6]


def dump(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,allow_nan=False))
    tmp.replace(path)


def process_cpu(handle):
    create,exit,kernel,user=[wintypes.FILETIME() for _ in range(4)]
    api=ctypes.WinDLL('kernel32',use_last_error=True).GetProcessTimes
    api.argtypes=[wintypes.HANDLE]+[ctypes.POINTER(wintypes.FILETIME)]*4
    api.restype=wintypes.BOOL
    if not api(int(handle),ctypes.byref(create),ctypes.byref(exit),ctypes.byref(kernel),ctypes.byref(user)):
        raise ctypes.WinError(ctypes.get_last_error())
    def seconds(value):return ((value.dwHighDateTime<<32)|value.dwLowDateTime)/1e7
    return dict(user_CPU_seconds=seconds(user),kernel_CPU_seconds=seconds(kernel),
                child_process_CPU_seconds=seconds(user)+seconds(kernel),
                process_creation_epoch=seconds(create)-11644473600,
                process_exit_epoch=seconds(exit)-11644473600)


def main():
    psutil.Process().cpu_affinity(CORES)
    manifest=json.loads((ROOT/'configs/classic_extension_main.json').read_text())
    folder=OUT/'queue';folder.mkdir(parents=True,exist_ok=True)
    if (folder/'QUEUE_STARTED.json').exists():raise FileExistsError('Do not overwrite an earlier queue')
    # Same frozen ordering as Mac; four lanes share only immutable input files.
    arms=['cie_single','cie_dual','global_dual','a_material','a_plus','physical_lbfgs','paper_nie','paper_lsie']
    pending=[t for arm in arms for t in manifest['tasks'] if t['arm']==arm]
    dump(folder/'QUEUE_STARTED.json',dict(status='STARTED',UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        epoch_seconds=START,lanes=4,allocated_physical_P_logical_cores=CORES,
        pin_before_numpy=True,threads_per_child=1,GPU=False,tasks=pending,
        external_CPU='Windows GetProcessTimes user+kernel after exit, retained process handle',
        scope='isolated own16-task classic cohort only; no other process modifications',
        platform=dict(system='Windows',python=sys.executable,psutil=psutil.__version__)))
    active={};free=list(CORES);rows=[]
    while pending or active:
        while pending and free:
            task=pending.pop(0);core=free.pop(0)
            if (OUT/task['task_id']).exists():raise FileExistsError('Existing task output '+task['task_id'])
            wrapper=("import os,sys,psutil,runpy;"
                f"psutil.Process().cpu_affinity([{core}]);"
                f"os.environ['CLASSIC_PINNED_CORE']={str(core)!r};"
                "os.environ['CUDA_VISIBLE_DEVICES']='';"
                f"sys.argv=['research/research_max_struct_r9/experiments/classic_extension.py','solve','--task-id',{task['task_id']!r},'--output-root','research/research_max_struct_r9/RESULTS/classic_extension_xinan'];"
                "runpy.run_path(sys.argv[0],run_name='__main__')")
            log=open(folder/(task['task_id']+'.log'),'w')
            started=time.time();p=subprocess.Popen([sys.executable,'-B','-c',wrapper],cwd=REPO,
                stdout=log,stderr=subprocess.STDOUT,env=os.environ.copy())
            active[p.pid]=dict(process=p,core=core,task=task,log=log,epoch_start=started,peak_rss=0)
            print('START',task['task_id'],'PID',p.pid,'P_CORE',core,flush=True)
        for pid,item in list(active.items()):
            p=item['process'];returncode=p.poll()
            if returncode is None:
                try:item['peak_rss']=max(item['peak_rss'],psutil.Process(pid).memory_info().peak_wset)
                except (psutil.NoSuchProcess,psutil.AccessDenied):pass
                continue
            cpu=process_cpu(p._handle)
            item['log'].close()
            dest=OUT/item['task']['task_id']
            if not (dest/'RECEIPT.json').exists() and not (dest/'FAILURE.json').exists():
                dump(dest/'FAILURE.json',dict(status='FAILED',task_id=item['task']['task_id'],
                    error='Worker exited before its own terminal receipt; inspect queue log',returncode=returncode,
                    process_cpu_seconds=cpu['child_process_CPU_seconds'],physical_cost={},
                    physical_count_status='UNKNOWN_IF_IMPORT_OR_EXTERNAL_TERMINATION',GT_loaded=False))
            row=dict(task_id=item['task']['task_id'],returncode=returncode,**cpu,
                epoch_start=item['epoch_start'],epoch_end=time.time(),
                process_wall_seconds=time.time()-item['epoch_start'],
                peak_rss_bytes=item['peak_rss'],peak_rss_semantics='sampled Windows peak_wset; may miss final tail',
                assigned_logical_P_core=item['core'],threads=1,GPU=False,
                external_CPU_semantics='exact GetProcessTimes; not Mac wait4',
                output_terminal='RECEIPT' if (dest/'RECEIPT.json').exists() else 'FAILURE')
            rows.append(row);dump(folder/'TASK_PROCESS_RECEIPTS.json',rows)
            print('FINISH',row['task_id'],returncode,row['child_process_CPU_seconds'],flush=True)
            free.append(item['core']);del active[pid]
        time.sleep(.25)
    dump(folder/'QUEUE_COMPLETE.json',dict(status='COMPLETED',tasks=len(rows),completed=sum(r['returncode']==0 for r in rows),
        failed=sum(r['returncode']!=0 for r in rows),UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        makespan_seconds=time.time()-START,epoch_start=START,epoch_end=time.time(),
        actual_child_CPU_seconds=sum(r['child_process_CPU_seconds'] for r in rows),
        queue_process_CPU_seconds=time.process_time()-CPU_START,physical_P_cores=CORES,
        CPU_cost_semantics='whole16 Windows cohort; no pooling with partial Mac',GPU=False))


if __name__=='__main__':main()
