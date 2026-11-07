"""Finite one-worker first two vector cases; no score and no GPU use.
Wait4 whole-child CPU supersedes nested timers in scientific receipts.
"""
import time,os,subprocess,sys,resource,json
from pathlib import Path
T=time.perf_counter();C=time.process_time();ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'delegated/vector/cpu_first_queue';OUT.mkdir(parents=True,exist_ok=False)
steps=[]
def run(label,arguments):
    start=time.perf_counter();log=(OUT/(label+'.log')).open('w')
    command=[sys.executable,str(ROOT/'experiments/vector_run.py'),*arguments,'--backend','cpu-dense']
    process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
    (OUT/'ACTIVE.json').write_text(json.dumps(dict(label=label,pid=process.pid,command=command,started_unix=time.time()),indent=2)+'\n')
    print(json.dumps(dict(event='START',label=label,pid=process.pid)),flush=True)
    while True:
        pid,status,usage=os.wait4(process.pid,os.WNOHANG)
        if pid:break
        time.sleep(1)
    process.returncode=os.waitstatus_to_exitcode(status);log.close()
    receipt=dict(label=label,command=command,exit=process.returncode,CPU_seconds=usage.ru_utime+usage.ru_stime,
        wall_seconds=time.perf_counter()-start,peak_RSS_bytes=usage.ru_maxrss,platform='Mac',additive=True,
        scope='Whole child wait4, including imports/sharedoperator/allcase arms/export. Internal child/perarm counters are nested.')
    steps.append(receipt);(OUT/'PROCESS_RECEIPTS.json').write_text(json.dumps(steps,indent=2)+'\n')
    print(json.dumps(dict(event='END',**receipt)),flush=True)
    if process.returncode:raise RuntimeError(label+' failed; preserve physics/evidence, do not autochange mathematics')
try:
    run('generate_initial2',['generate','--ids','v02','v09','--output','RESULTS/vector/generation_cpu_initial2'])
    run('case_v02',['invert','--max-tasks','3','--task-offset','0','--output','RESULTS/vector/case_v02_cpu_process'])
    run('case_v09',['invert','--max-tasks','3','--task-offset','3','--output','RESULTS/vector/case_v09_cpu_process'])
    status='COMPLETE'
except BaseException as exc:
    status='FAILED';print(repr(exc),flush=True)
finally:
    (OUT/'CONTROLLER_COST.json').write_text(json.dumps(dict(status=status,platform='Mac',CPU_seconds=time.process_time()-C,
        wall_seconds=time.perf_counter()-T,additive=True,scope='Own finite controller only; wait4 children charged separately, no score',
        steps=steps,workers=1,BLAS_threads_per_child=1,GPU_use=False),indent=2)+'\n')
if status!='COMPLETE':raise SystemExit(1)
