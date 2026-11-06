"""Entire v08 matched three-arm cohort, one resident CUDA inverse process."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import time
T=time.perf_counter()
from pathlib import Path
import json,sys,subprocess,ctypes,importlib.util
from ctypes import wintypes as w
import psutil
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.src.vector_resources import configure_own_windows_process,CooperativeGate
OUT=ROOT/'delegated/vector/windows_v08_gpu_queue';OUT.mkdir(parents=True,exist_ok=False)
process_config=configure_own_windows_process();process=psutil.Process()
spec=importlib.util.spec_from_file_location('job_helper',REPO/'research/research_max_struct_r9/delegated/3d_extension/windows_resolved.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);job=module.job_object();os.environ['R10_PARENT_JOB_ACTIVE']='1'
gate=CooperativeGate(OUT/'PARENT_GPU_BOUNDARIES.jsonl');gate.wait('before_v08_generation')
(OUT/'RESOURCE_RECEIPT.json').write_text(json.dumps(dict(process=process_config,available_RAM_bytes=psutil.virtual_memory().available,job=job,GPU_memory_policy_bytes=1500000000,cohort='v08_matched12_all3_CUDA'),indent=2)+'\n')
k=ctypes.WinDLL('kernel32',use_last_error=True);k.GetProcessTimes.argtypes=[w.HANDLE]+[ctypes.POINTER(w.FILETIME)]*4;k.GetProcessTimes.restype=w.BOOL
rows=[]
def run(label,args):
    gate.wait('before_'+label);start=time.perf_counter()
    with (OUT/(label+'.log')).open('w') as log:
        command=[sys.executable,str(ROOT/'experiments/vector_run.py'),*args,'--backend','cuda']
        p=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
        (OUT/'ACTIVE.json').write_text(json.dumps(dict(label=label,pid=p.pid,command=command),indent=2)+'\n');print(json.dumps(dict(event='START',label=label,pid=p.pid)),flush=True)
        status=p.wait();a,b,c,d=(w.FILETIME() for _ in range(4))
        if not k.GetProcessTimes(w.HANDLE(int(p._handle)),ctypes.byref(a),ctypes.byref(b),ctypes.byref(c),ctypes.byref(d)):raise ctypes.WinError(ctypes.get_last_error())
        cpu=((c.dwHighDateTime<<32)+c.dwLowDateTime+(d.dwHighDateTime<<32)+d.dwLowDateTime)*1e-7
    row=dict(label=label,exit=status,CPU_seconds=cpu,wall_seconds=time.perf_counter()-start,platform='XINAN',additive=True,CPU_source='GetProcessTimes whole child, not GPU active time',nested_child_timers_not_additive=True)
    rows.append(row);(OUT/'PROCESS_RECEIPTS.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(dict(event='END',**row)),flush=True)
    if status:raise RuntimeError(label+' failed, do not autochange numerical contract')
try:
    run('generate_v08',['generate','--ids','v08','--generation-sides','12','--output','RESULTS/vector/generation_gpu_v08'])
    run('invert_v08',['invert','--queue','configs/vector_gpu_v08_queue.json','--max-tasks','3','--output','RESULTS/vector/case_v08_gpu_process'])
    status='COMPLETE'
except BaseException as exc:status='FAILED';print(repr(exc),flush=True)
finally:
    cpu=process.cpu_times();(OUT/'CONTROLLER_COST.json').write_text(json.dumps(dict(status=status,platform='XINAN',CPU_seconds=cpu.user+cpu.system,
        wall_seconds=time.perf_counter()-T,rows=rows,additive=True,scope='Own finite connected parent only; child scientific CPU and CUDA-event views separate',
        GPU_gate_helper_CPU_seconds=gate.helper_CPU_seconds,GPU_gate_helper_scope='Separate own query child CPU, not included in main psutil CPU',job=job),indent=2)+'\n')
if status!='COMPLETE':raise SystemExit(1)
