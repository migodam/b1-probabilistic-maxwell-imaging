#!/usr/bin/env python3
"""python tools/run_job.py STAGE LABEL COMMAND...; serial, budgeted jobs."""
import fcntl
import math
import json
import tempfile
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from wslam.accounting import BudgetBook, BudgetExceeded, stage_key, json_default


def _group_cpu(pgid):
    # The final authoritative aggregate is RUSAGE_CHILDREN. This live sample is
    # an additional guard for concurrent descendant CPU before the job exits.
    try:
        output = subprocess.run(['ps','-axo','pgid=,time='],capture_output=True,text=True,check=False).stdout
    except PermissionError:
        # Managed macOS sessions may forbid global process enumeration. These
        # jobs do all physical work in one child process (no detached workers);
        # its RLIMIT_CPU and final RUSAGE_CHILDREN accounting remain active.
        return 0.0
    total = 0.0
    for row in output.splitlines():
        fields = row.split()
        if len(fields) != 2 or fields[0] != str(pgid):
            continue
        value = fields[1]
        days = 0
        if '-' in value:
            day,value = value.split('-',1)
            days = int(day)
        parts = value.split(':')
        seconds = float(parts[-1])
        if len(parts)>1:
            seconds += int(parts[-2])*60
        if len(parts)>2:
            seconds += int(parts[-3])*3600
        total += days*86400+seconds
    return total


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv)<3:
        print(__doc__,file=sys.stderr)
        return 2
    phase,label,command = stage_key(argv[0]),argv[1],argv[2:]
    if not label or Path(label).name != label or label in ('.', '..'):
        raise ValueError('LABEL must be one safe directory name')
    job_dir = ROOT/'results'/'jobs'/label
    job_dir.mkdir(parents=True, exist_ok=True)
    output_path, receipt_path = job_dir/'OUTPUT.log', job_dir/'RECEIPT.json'
    env = os.environ.copy()
    env.update(PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1',BLIS_NUM_THREADS='1',WSLAM_STAGE=phase)
    env['PYTHONPATH'] = os.pathsep.join(map(str,(ROOT/'src',ROOT/'vendor',ROOT.parent/'research-w-continuation'/'.deps')))
    env['MPLCONFIGDIR'] = str(ROOT/'results'/'.matplotlib-cache')
    ledger = Path(env.get('WSLAM_LEDGER',ROOT/'results'/'budget_ledger.json')).resolve()
    # Each project job uses the same lock, avoiding competing CPU allowances.
    ledger.parent.mkdir(parents=True,exist_ok=True)
    with ledger.with_suffix(ledger.suffix+'.job.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        env['WSLAM_LEDGER'] = str(ledger)
        book = BudgetBook(ledger)
        remaining = book.snapshot(phase)[phase]['remaining']['process_cpu_seconds']
        if remaining<=0:
            raise BudgetExceeded(phase+': CPU budget exhausted before launch')
        book.reserve(phase,label,command=command,source='run_job')
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        self_before = time.process_time()
        started = time.monotonic()
        child = None
        peak_cpu = 0.0
        timed_out = False
        exit_code = 1
        failure = None
        output = output_path.open('a', buffering=1)
        output.write('\n'+json.dumps(dict(stage=phase,label=label,command=command,start_unix=time.time()))+'\n')
        output.flush()
        try:
            def set_limits():
                os.setsid()
                _,hard = resource.getrlimit(resource.RLIMIT_CPU)
                cap = max(1,int(math.ceil(remaining)))
                if hard != resource.RLIM_INFINITY:
                    cap = min(cap,hard)
                resource.setrlimit(resource.RLIMIT_CPU,(cap,cap))
            child = subprocess.Popen(command,cwd=ROOT,env=env,preexec_fn=set_limits,stdout=output,stderr=subprocess.STDOUT)
            while child.poll() is None:
                peak_cpu = max(peak_cpu,_group_cpu(child.pid))
                if peak_cpu + (time.process_time()-self_before) >= remaining:
                    timed_out = True
                    os.killpg(child.pid,signal.SIGKILL)
                    break
                time.sleep(0.25)
            exit_code = child.wait()
        except BaseException as exc:
            failure = repr(exc)
            if child is not None:
                try:
                    os.killpg(child.pid,signal.SIGKILL)
                except ProcessLookupError:
                    pass
                child.wait()
            raise
        finally:
            # Stop unjoined descendants before releasing the serial job lock.
            if child is not None:
                try:
                    os.killpg(child.pid,signal.SIGKILL)
                except ProcessLookupError:
                    pass
            after = resource.getrusage(resource.RUSAGE_CHILDREN)
            cpu = after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime + time.process_time()-self_before
            wall = time.monotonic()-started
            output.flush()
            os.fsync(output.fileno())
            output.close()
            receipt = book.event(phase,label,kind='job_completion',command=command,exit_code=exit_code,failure=failure,process_cpu_seconds=max(cpu,peak_cpu),wall_seconds=wall,children_maxrss=after.ru_maxrss,rss_unit='bytes' if sys.platform=='darwin' else 'KiB',rss_semantics='RUSAGE_CHILDREN high-water mark for wrapper lifetime; not additive',live_group_peak_cpu_seconds=peak_cpu,cpu_limit_remaining_at_launch=remaining,cpu_limit_exceeded=timed_out or cpu>=remaining,output_log=str(output_path),receipt_path=str(receipt_path))
            fd, temporary = tempfile.mkstemp(prefix='receipt.',dir=job_dir)
            try:
                with os.fdopen(fd,'w') as stream:
                    json.dump(receipt,stream,indent=2,allow_nan=False,default=json_default)
                    stream.flush()
                    os.fsync(stream.fileno())
                os.replace(temporary,receipt_path)
            finally:
                if os.path.exists(temporary):
                    os.unlink(temporary)
            print(json.dumps(dict(stage=phase,label=label,exit_code=exit_code,receipt=str(receipt_path),output=str(output_path))))
        return 124 if timed_out else exit_code if exit_code>=0 else 128-exit_code


if __name__ == '__main__':
    try:
        sys.exit(main())
    except BudgetExceeded as exc:
        print(str(exc),file=sys.stderr)
        sys.exit(125)
