"""Finite cooperative GPU gate. Inspects only; never changes other jobs.

The gate cannot reserve future GPU availability. A granted idle window is
still required. Active external CUDA PIDs pause the next work unit. Resident
memory while paused must be reported and bounded by the approved contract.
"""
import csv
import io
import json
import os
from pathlib import Path
import subprocess
import time


def gpu_process_snapshot():
    command = ['nvidia-smi', '--query-compute-apps=pid,process_name,used_gpu_memory',
               '--format=csv,noheader,nounits']
    wall=time.perf_counter()
    process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    stdout,stderr=process.communicate(timeout=10)
    helper_cpu=None
    if os.name=='nt':
        import ctypes
        from ctypes import wintypes as w
        creation,exit_time,kernel,user=(w.FILETIME() for _ in range(4))
        function=ctypes.windll.kernel32.GetProcessTimes
        function.argtypes=[w.HANDLE]+[ctypes.POINTER(w.FILETIME)]*4
        function.restype=w.BOOL
        if function(w.HANDLE(int(process._handle)),ctypes.byref(creation),ctypes.byref(exit_time),ctypes.byref(kernel),ctypes.byref(user)):
            helper_cpu=((kernel.dwHighDateTime<<32)+kernel.dwLowDateTime+(user.dwHighDateTime<<32)+user.dwLowDateTime)*1e-7
    if process.returncode:
        raise RuntimeError('GPU process visibility unavailable: '+stderr.strip())
    rows = []
    for row in csv.reader(io.StringIO(stdout)):
        if not row or not row[0].strip():
            continue
        try: pid = int(row[0].strip())
        except ValueError: raise RuntimeError('Unparseable GPU process inventory')
        memory = row[2].strip() if len(row)>2 else None
        rows.append(dict(pid=pid, process_name=row[1].strip(),
            memory_MiB=int(memory) if memory and memory.isdigit() else None,
            unknown_memory_is_not_zero=True))
    return dict(unix_time=time.time(), own_pid=os.getpid(), processes=rows,
        external=[r for r in rows if r['pid'] != os.getpid()],
        command=command,query_helper_CPU_seconds=helper_cpu,query_helper_wall_seconds=time.perf_counter()-wall,
        visibility='Compute-app listing; WDDM perprocess memory may be unavailable')


class CooperativeGate:
    def __init__(self, log, poll_seconds=10, memory_limit_bytes=1500000000, minimum_query_interval_seconds=1.):
        self.log = Path(log); self.poll_seconds = poll_seconds
        self.memory_limit_bytes = memory_limit_bytes
        self.helper_CPU_seconds=0.;self.unmeasured_helper_invocations=0;self.query_count=0
        self.minimum_query_interval_seconds=minimum_query_interval_seconds
        self.last_query_time=-float('inf');self.last_snapshot=None;self.cached_boundary_checks=0
        self.paused_wall_seconds=0.

    def wait(self, boundary, memory_bytes=None):
        """No CUDA work or context is created by this inspection function."""
        if memory_bytes is not None and memory_bytes > self.memory_limit_bytes:
            raise RuntimeError('Own measured allocation exceeds approved memory limit')
        if self.last_snapshot is not None and not self.last_snapshot['external'] and time.perf_counter()-self.last_query_time<self.minimum_query_interval_seconds:
            self.cached_boundary_checks+=1
            return dict(action='ALLOW_RECENT_IDLE_SNAPSHOT',boundary=boundary,
                snapshot_unix_time=self.last_snapshot['unix_time'],maximum_cached_age_seconds=self.minimum_query_interval_seconds)
        announced=False
        while True:
            row = gpu_process_snapshot(); row.update(boundary=boundary, own_memory_bytes=memory_bytes)
            self.last_query_time=time.perf_counter();self.last_snapshot=row
            self.query_count+=1
            if row['query_helper_CPU_seconds'] is None:self.unmeasured_helper_invocations+=1
            else:self.helper_CPU_seconds+=row['query_helper_CPU_seconds']
            row['action'] = 'WAIT_EXTERNAL_GPU' if row['external'] else 'ALLOW_NEXT_WORK_UNIT'
            self.log.parent.mkdir(parents=True, exist_ok=True)
            with self.log.open('a') as handle: handle.write(json.dumps(row)+'\n')
            if not row['external']:
                if announced:print(json.dumps(dict(event='GPU_BOUNDARY_RESUMED',boundary=boundary)),flush=True)
                return row
            if not announced:
                print(json.dumps(dict(event='WAIT_EXTERNAL_GPU',boundary=boundary,external_pids=[r['pid'] for r in row['external']])),flush=True)
                announced=True
            # Finite job waiting, not a daemon/scheduler. Never kill or modify peers.
            start=time.perf_counter();time.sleep(self.poll_seconds);self.paused_wall_seconds+=time.perf_counter()-start


def configure_own_windows_process():
    import psutil
    process = psutil.Process()
    process.cpu_affinity([12, 14])
    process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    return dict(pid=process.pid, affinity=process.cpu_affinity(),
        priority='BelowNormal', scope='Only this authorized process; peer settings unchanged')
