"""Whole v08 primary replacement; same public bytes, qualified retained CUDA LU."""
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
from research.research_max_struct_r10.src.vector_resources import configure_own_windows_process
OUT=ROOT/'delegated/vector/windows_v08_dense_primary_queue';OUT.mkdir(parents=True,exist_ok=False)
process_config=configure_own_windows_process();process=psutil.Process()
spec=importlib.util.spec_from_file_location('own_job',REPO/'research/research_max_struct_r9/delegated/3d_extension/windows_resolved.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);job=module.job_object();os.environ['R10_PARENT_JOB_ACTIVE']='1'
for arm in ('sourceCIE','GlobalDual','AplusGlobal'):
    if (ROOT/'RESULTS/vector'/('v08_'+arm)).exists():raise RuntimeError('Old task output not archived; no overwrite')
command=[sys.executable,str(ROOT/'experiments/vector_run.py'),'invert','--queue','configs/vector_gpu_v08_queue.json','--max-tasks','3','--backend','cuda-dense','--output','RESULTS/vector/case_v08_cuda_dense_process']
start=time.perf_counter()
with (OUT/'RUN.log').open('w') as log:
    child=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
    (OUT/'ACTIVE.json').write_text(json.dumps(dict(pid=child.pid,command=command,started_unix=time.time()),indent=2)+'\n');print(json.dumps(dict(event='START_V08_WHOLE_CUDA_DENSE_PRIMARY',pid=child.pid)),flush=True)
    status=child.wait();a,b,c,d=(w.FILETIME() for _ in range(4))
    k=ctypes.WinDLL('kernel32',use_last_error=True);k.GetProcessTimes.argtypes=[w.HANDLE]+[ctypes.POINTER(w.FILETIME)]*4;k.GetProcessTimes.restype=w.BOOL
    if not k.GetProcessTimes(w.HANDLE(int(child._handle)),ctypes.byref(a),ctypes.byref(b),ctypes.byref(c),ctypes.byref(d)):raise ctypes.WinError(ctypes.get_last_error())
    cpu=((c.dwHighDateTime<<32)+c.dwLowDateTime+(d.dwHighDateTime<<32)+d.dwLowDateTime)*1e-7
row=dict(status='COMPLETE' if status==0 else 'FAILED',exit=status,CPU_seconds=cpu,wall_seconds=time.perf_counter()-start,platform='XINAN',additive=True,
    CPU_source='GetProcessTimes full scientific child; replaces nested process/arm CPU',GPU_device_event_not_CPU=True,publicY_regenerated=False,
    scope='One new whole v08 three-arm retained-dense primary; original FFT attempt is separately preserved engineering spend')
(OUT/'PROCESS_COST.json').write_text(json.dumps(row,indent=2)+'\n');own=process.cpu_times();(OUT/'CONTROLLER_COST.json').write_text(json.dumps(dict(CPU_seconds=own.user+own.system,wall_seconds=time.perf_counter()-T,platform='XINAN',additive=True,scope='Own connected finite controller, excludes childCPU',job=job),indent=2)+'\n');print(json.dumps(row),flush=True)
raise SystemExit(status)
