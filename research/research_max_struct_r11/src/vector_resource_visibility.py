"""R11-only finite retry of own nvidia-smi timeout, before any CUDA work."""
import os,time,json,csv,io,subprocess
from pathlib import Path

def child_cpu(process):
 if os.name!='nt':return None
 import ctypes
 from ctypes import wintypes as w
 a,b,c,d=(w.FILETIME() for _ in range(4));f=ctypes.windll.kernel32.GetProcessTimes
 f.argtypes=[w.HANDLE]+[ctypes.POINTER(w.FILETIME)]*4;f.restype=w.BOOL
 if not f(w.HANDLE(int(process._handle)),ctypes.byref(a),ctypes.byref(b),ctypes.byref(c),ctypes.byref(d)):raise ctypes.WinError()
 return ((c.dwHighDateTime<<32)+c.dwLowDateTime+(d.dwHighDateTime<<32)+d.dwLowDateTime)*1e-7

def install(log,max_attempts=8):
 from research.research_max_struct_r10.src import vector_resources as inherited
 path=Path(log);path.parent.mkdir(parents=True,exist_ok=True)
 def snapshot():
  command=['nvidia-smi','--query-compute-apps=pid,process_name,used_gpu_memory','--format=csv,noheader,nounits']
  start=time.perf_counter();total_cpu=0.;timeouts=0
  for attempt in range(1,int(max_attempts)+1):
   wall=time.perf_counter();process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
   try:stdout,stderr=process.communicate(timeout=10)
   except subprocess.TimeoutExpired:
    process.kill();stdout,stderr=process.communicate();cost=child_cpu(process)
    if cost is None:raise RuntimeError('Timeout helperCPU unavailable; do not substitute zero')
    total_cpu+=cost;timeouts+=1
    with path.open('a') as f:f.write(json.dumps(dict(kind='OWN_NVIDIA_SMI_TIMEOUT_KILLED',unix_time=time.time(),pid=process.pid,attempt=attempt,
      helper_CPU_seconds=cost,helper_wall_seconds=time.perf_counter()-wall,new_physical_RHS=0,CUDA_work_submitted=False))+'\n')
    if attempt==int(max_attempts):raise RuntimeError('GPU visibility unavailable after finite timeout retries; no CUDA work submitted')
    time.sleep(5);continue
   cost=child_cpu(process)
   if cost is None:raise RuntimeError('Valid helperCPU unavailable; do not substitute zero')
   total_cpu+=cost
   if process.returncode:raise RuntimeError('GPU process visibility unavailable:'+stderr.strip())
   rows=[]
   for row in csv.reader(io.StringIO(stdout)):
    if not row or not row[0].strip():continue
    try:pid=int(row[0].strip())
    except ValueError:raise RuntimeError('Unparseable GPU process inventory')
    memory=row[2].strip() if len(row)>2 else None
    rows.append(dict(pid=pid,process_name=row[1].strip(),memory_MiB=int(memory) if memory and memory.isdigit() else None,unknown_memory_is_not_zero=True))
   result=dict(unix_time=time.time(),own_pid=os.getpid(),processes=rows,external=[r for r in rows if r['pid']!=os.getpid()],
    command=command,query_helper_CPU_seconds=total_cpu,query_helper_wall_seconds=time.perf_counter()-start,
    query_attempts=attempt,timeout_attempts=timeouts,visibility='Original compute-app listing; own helper timeout finite retry, no CUDA without valid snapshot')
   if timeouts:
    with path.open('a') as f:f.write(json.dumps(dict(kind='VALID_VISIBILITY_AFTER_TIMEOUT',**result))+'\n')
   return result
 inherited.gpu_process_snapshot=snapshot
 return dict(scope='Only R11 own nvidia-smi TimeoutExpired retry; original valid snapshot parser and gate external-PID policy retained',max_attempts=int(max_attempts),timeout_seconds=10,backoff_seconds=5,no_R10_file_modified=True)
