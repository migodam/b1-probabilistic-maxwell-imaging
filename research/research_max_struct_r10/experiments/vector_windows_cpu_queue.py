"""Three v04 arms, existing Windows CPU-dense/FFT, no CUDA/GT score."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import time
T=time.perf_counter()
from pathlib import Path
import json,sys,subprocess,ctypes,importlib.util
from ctypes import wintypes as w
import psutil
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
OUT=ROOT/'delegated/vector/windows_v04_queue';OUT.mkdir(parents=True,exist_ok=False)
process=psutil.Process();process.cpu_affinity([12,14]);process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
spec=importlib.util.spec_from_file_location('job_helper',REPO/'research/research_max_struct_r9/delegated/3d_extension/windows_resolved.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);job=module.job_object()
os.environ['R10_PARENT_JOB_ACTIVE']='1'
memory=psutil.virtual_memory();load=psutil.cpu_percent(interval=1,percpu=True)
resource=dict(available_RAM_bytes=memory.available,total_RAM_bytes=memory.total,logical_CPU_load={str(i):load[i] for i in (12,14)},
    CPU_affinity=process.cpu_affinity(),priority='BelowNormal',GPU_use=False,job=job)
(OUT/'RESOURCE_RECEIPT.json').write_text(json.dumps(resource,indent=2)+'\n');print(json.dumps(dict(event='CPU_RESOURCE_RECEIPT',**resource)),flush=True)
if memory.available<8_000_000_000:raise RuntimeError('Insufficient current RAM; root review required, no numeric launch')
k=ctypes.WinDLL('kernel32',use_last_error=True);k.GetProcessTimes.argtypes=[w.HANDLE]+[ctypes.POINTER(w.FILETIME)]*4;k.GetProcessTimes.restype=w.BOOL
rows=[]
def run(label,arguments):
    start=time.perf_counter();command=[sys.executable,str(ROOT/'experiments/vector_run.py'),*arguments,'--backend','cpu-dense']
    with (OUT/(label+'.log')).open('w') as log:
        p=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
        (OUT/'ACTIVE.json').write_text(json.dumps(dict(label=label,pid=p.pid,command=command),indent=2)+'\n')
        print(json.dumps(dict(event='START',label=label,pid=p.pid)),flush=True);exit_code=p.wait()
        a,b,c,d=(w.FILETIME() for _ in range(4))
        if not k.GetProcessTimes(w.HANDLE(int(p._handle)),ctypes.byref(a),ctypes.byref(b),ctypes.byref(c),ctypes.byref(d)):raise ctypes.WinError(ctypes.get_last_error())
        cpu=((c.dwHighDateTime<<32)+c.dwLowDateTime+(d.dwHighDateTime<<32)+d.dwLowDateTime)*1e-7
    row=dict(label=label,exit=exit_code,CPU_seconds=cpu,wall_seconds=time.perf_counter()-start,platform='XINAN',
        CPU_source='GetProcessTimes retained child handle after exit',additive=True,scope='Whole child process; nested child/perarm timers not another charge')
    rows.append(row);(OUT/'PROCESS_RECEIPTS.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(dict(event='END',**row)),flush=True)
    if exit_code:raise RuntimeError(label+' failed, no automatic mathematical repair')
try:
    run('generate_v04',['generate','--ids','v04','--generation-sides','12','--output','RESULTS/vector/generation_windows_v04'])
    run('invert_v04',['invert','--queue','configs/vector_remote_v04_queue.json','--max-tasks','3','--output','RESULTS/vector/case_v04_windows_process'])
    status='COMPLETE'
except BaseException as exc:status='FAILED';print(repr(exc),flush=True)
finally:
    cpu=process.cpu_times();(OUT/'CONTROLLER_COST.json').write_text(json.dumps(dict(status=status,platform='XINAN',
        CPU_seconds=cpu.user+cpu.system,wall_seconds=time.perf_counter()-T,rows=rows,additive=True,
        scope='Own finite CPU queue process only; held children separately charged',GPU_use=False,resource=resource),indent=2)+'\n')
if status!='COMPLETE':raise SystemExit(1)
