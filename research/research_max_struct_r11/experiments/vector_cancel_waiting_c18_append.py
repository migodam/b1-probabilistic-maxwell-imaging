"""Root-conditional engineering cancellation of an owned waiting0FW controller."""
import os,time,json,subprocess,ctypes,datetime
from pathlib import Path
from ctypes import wintypes as w
ENTRY=time.perf_counter();ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'delegated/vector/c18_gpu_append_controller'

def dump(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n')
k=ctypes.WinDLL('kernel32',use_last_error=True)
k.GetProcessTimes.argtypes=[w.HANDLE]+[ctypes.POINTER(w.FILETIME)]*4;k.GetProcessTimes.restype=w.BOOL
k.OpenProcess.argtypes=[w.DWORD,w.BOOL,w.DWORD];k.OpenProcess.restype=w.HANDLE
k.TerminateProcess.argtypes=[w.HANDLE,w.UINT];k.TerminateProcess.restype=w.BOOL
k.WaitForSingleObject.argtypes=[w.HANDLE,w.DWORD];k.CloseHandle.argtypes=[w.HANDLE]
def process_times(handle):
 a,b,c,d=(w.FILETIME() for _ in range(4))
 if not k.GetProcessTimes(handle,ctypes.byref(a),ctypes.byref(b),ctypes.byref(c),ctypes.byref(d)):raise ctypes.WinError(ctypes.get_last_error())
 ticks=lambda t:(t.dwHighDateTime<<32)+t.dwLowDateTime
 return dict(CPU_seconds=(ticks(c)+ticks(d))*1e-7,created_unix=ticks(a)*1e-7-11644473600.)
helperCPU=0.
def cim(filter_text):
 global helperCPU
 command="Get-CimInstance Win32_Process | Where-Object { "+filter_text+" } | Select-Object ProcessId,ParentProcessId | ConvertTo-Json -Compress"
 child=subprocess.Popen(['powershell','-NoProfile','-Command',command],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
 result,error=child.communicate();helperCPU+=process_times(w.HANDLE(int(child._handle)))['CPU_seconds']
 if child.returncode:raise RuntimeError('Read-only owned-process inventory failed:'+error)
 data=json.loads(result) if result.strip() else []
 return data if isinstance(data,list) else [data]
status='FAILED';failure=None
try:
 barrier=json.loads((ROOT/'delegated/vector/confirmation_host_union/ORIGINAL_MAC15_TERMINAL_FOR_C18_RETURN.json').read_text())
 if barrier.get('status')!='ALL_ORIGINAL_MAC15_NATURAL_TERMINAL_PUBLIC_SEALS' or len(set(barrier.get('task_ids',[])))!=15:raise RuntimeError('OriginalMac15 public terminal barrier missing')
 if (ROOT/'delegated/vector/confirmation_process/PROCESS_COST.json').exists():raise RuntimeError('InitialGPU25 already terminal; refuse cancellation race')
 terminal=json.loads((ROOT/'delegated/vector/confirmation_process/TERMINALS.json').read_text())
 if len(terminal['terminals'])>=24:raise RuntimeError('InitialGPU25 too near completion; report instead of cancellation race')
 if any((ROOT/('delegated/vector/'+n)).exists() for n in ('c18_gpu_append_process','c18_gpu_append_inversions')):raise RuntimeError('c18 scientific output already exists')
 if (OUT/'ACTIVE.json').exists() or (OUT/'SCIENTIFIC_PROCESS_COST.json').exists():raise RuntimeError('Append scientific start exists; never interrupt')
 records=cim("$_.Name -eq 'python.exe' -and $_.CommandLine -like '*vector_gpu_c18_append_controller.py*'")
 if len(records)!=1:raise RuntimeError('Require exactly one owned waiting controller')
 pid=int(records[0]['ProcessId'])
 if cim('$_.ParentProcessId -eq '+str(pid)):raise RuntimeError('Owned controller has child; no cancellation')
 handle=k.OpenProcess(0x0400|0x00100000|0x0001,False,pid)
 if not handle:raise ctypes.WinError(ctypes.get_last_error())
 try:
  if (OUT/'ACTIVE.json').exists() or (ROOT/'delegated/vector/confirmation_process/PROCESS_COST.json').exists():raise RuntimeError('Append barrier changed; preserve task')
  if not k.TerminateProcess(handle,1223):raise ctypes.WinError(ctypes.get_last_error())
  k.WaitForSingleObject(handle,10000);measured=process_times(handle)
 finally:k.CloseHandle(handle)
 dump(OUT/'CANCELLATION_ENGINEERING_RECEIPT.json',dict(status='CONTROLLED_CANCEL_WAITING_APPEND_ZERO_SCIENCE',pid=pid,
  UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),authority='Root conditional preauthorization after originalMac15natural completion',
  original_GPU25_interrupted=False,scientific_child_ever_started=False,new_physical_RHS=0,new_inverse_calls=0,GT_access=False,
  originalGPU_scientific_terminals_at_cancel=len(terminal['terminals']),waiting_controller_CPU_seconds=measured['CPU_seconds'],
  waiting_controller_wall_seconds=time.time()-measured['created_unix'],CPU_scope='Whole terminated waiting controller GetProcessTimes; replaces missing nested controllerCPU,not science',
  controlled_exit_code=1223,raw_transport_failure_to_be_retained_and_classified=True))
 status='COMPLETE_ZERO_FW_CANCELLATION'
except BaseException as exc:
 import traceback
 failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
finally:
 dump(OUT/'CANCELLATION_ACTION_COST.json',dict(status=status,failure=failure,own_CPU_seconds=time.process_time(),
  powershell_inventory_child_CPU_seconds=helperCPU,wall_seconds=time.perf_counter()-ENTRY,scope='Only cancellation engineering action,not waiting controllerCPU nor scientificCPU',new_physical_RHS=0,GT_access=False))
print(json.dumps(dict(status=status,failure=failure)),flush=True)
raise SystemExit(0 if status=='COMPLETE_ZERO_FW_CANCELLATION' else 1)
