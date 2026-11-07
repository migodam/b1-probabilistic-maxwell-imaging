"""Stop only own pre-CUDA waiting probe; collect retained handle CPU."""
import time,json,ctypes
from ctypes import wintypes as w
from pathlib import Path
import psutil
T=time.perf_counter();C=time.process_time()
ROOT=Path('D:/AI/MAX_STRUCT_R10_VECTOR_20261107/research/research_max_struct_r10');out=ROOT/'RESULTS/vector/queue_probe'
config=json.loads((out/'PROCESS_CONFIGURATION.json').read_text());pid=config['process']['pid'];p=psutil.Process(pid)
assert 'vector_run.py' in ' '.join(p.cmdline()) and 'probe' in p.cmdline()
assert not (out/'SHARED_SETUP.json').exists() and not (ROOT/'RESULTS/vector/generation_initial2').exists()
created=p.create_time();k=ctypes.WinDLL('kernel32',use_last_error=True)
k.OpenProcess.argtypes=[w.DWORD,w.BOOL,w.DWORD];k.OpenProcess.restype=w.HANDLE
h=k.OpenProcess(0x1000|0x100000|1,False,pid)
if not h:raise ctypes.WinError(ctypes.get_last_error())
k.TerminateProcess.argtypes=[w.HANDLE,w.UINT];k.TerminateProcess.restype=w.BOOL
k.WaitForSingleObject.argtypes=[w.HANDLE,w.DWORD];k.WaitForSingleObject.restype=w.DWORD
k.GetProcessTimes.argtypes=[w.HANDLE]+[ctypes.POINTER(w.FILETIME)]*4;k.GetProcessTimes.restype=w.BOOL
if not k.TerminateProcess(h,77):raise ctypes.WinError(ctypes.get_last_error())
k.WaitForSingleObject(h,30000)
a,b,c,d=(w.FILETIME() for _ in range(4))
if not k.GetProcessTimes(h,ctypes.byref(a),ctypes.byref(b),ctypes.byref(c),ctypes.byref(d)):raise ctypes.WinError(ctypes.get_last_error())
cpu=((c.dwHighDateTime<<32)+c.dwLowDateTime+(d.dwHighDateTime<<32)+d.dwLowDateTime)*1e-7
rows=[json.loads(x) for x in (out/'GPU_BOUNDARIES.jsonl').read_text().splitlines()]
receipt=dict(status='OWN_WAITING_PROBE_STOPPED_FOR_ROOT_AUTHORIZED_CPU_RESOURCE_HANDOFF',platform='XINAN',pid=pid,
 CPU_seconds=cpu,wall_seconds=time.time()-created,CPU_source='GetProcessTimes retained process handle after termination',
 GPU_gate_helper_CPU_seconds=sum(x['query_helper_CPU_seconds'] for x in rows if x.get('query_helper_CPU_seconds') is not None),
 gate_queries=len(rows),possible_inflight_query_tail_CPU_seconds=None,private_data_generated=False,physical_RHS=0,
 CUDA_context_created=False,algorithm_failure=False,other_processes_modified=False,
 own_stop_helper_CPU_seconds=time.process_time()-C,own_stop_helper_wall_seconds=time.perf_counter()-T)
(out/'STOPPED_RESOURCE_HANDOFF_COST.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
