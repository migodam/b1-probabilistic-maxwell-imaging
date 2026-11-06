"""Finite nine-arm Mac extension, starts only after original six children exit."""
import time,os,sys,subprocess,resource,json
from pathlib import Path
T=time.perf_counter();C=time.process_time();ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'delegated/vector/mac_cross_queue';OUT.mkdir(parents=True,exist_ok=False)
steps=[]
def run(label,arguments):
    start=time.perf_counter();command=[sys.executable,str(ROOT/'experiments/vector_run.py'),*arguments,'--backend','cpu-dense']
    log=(OUT/(label+'.log')).open('w');p=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
    (OUT/'ACTIVE.json').write_text(json.dumps(dict(label=label,pid=p.pid,command=command),indent=2)+'\n');print(json.dumps(dict(event='START',label=label,pid=p.pid)),flush=True)
    while True:
        pid,status,usage=os.wait4(p.pid,os.WNOHANG)
        if pid:break
        time.sleep(1)
    p.returncode=os.waitstatus_to_exitcode(status);log.close();row=dict(label=label,exit=p.returncode,CPU_seconds=usage.ru_utime+usage.ru_stime,
        wall_seconds=time.perf_counter()-start,peak_RSS_bytes=usage.ru_maxrss,platform='Mac',additive=True,scope='Whole wait4 scientific child, nested timers not additive')
    steps.append(row);(OUT/'PROCESS_RECEIPTS.json').write_text(json.dumps(steps,indent=2)+'\n');print(json.dumps(dict(event='END',**row)),flush=True)
    if p.returncode:raise RuntimeError(label+' failed; no automatic math/budget change')
try:
    first=ROOT/'delegated/vector/cpu_first_queue/CONTROLLER_COST.json';print(json.dumps(dict(event='WAIT_FIRST_SIX_EXIT')),flush=True)
    while not first.exists():time.sleep(5)
    if json.loads(first.read_text())['status']!='COMPLETE':raise RuntimeError('First six queue not complete; review required')
    run('cross6',['invert','--queue','configs/vector_mac_cross_queue.json','--max-tasks','6','--output','RESULTS/vector/mac_cross6_process'])
    status='COMPLETE'
except BaseException as exc:status='FAILED';print(repr(exc),flush=True)
finally:
    (OUT/'CONTROLLER_COST.json').write_text(json.dumps(dict(status=status,CPU_seconds=time.process_time()-C,wall_seconds=time.perf_counter()-T,
        platform='Mac',additive=True,scope='Own finite waiting/controller only; scientific children separately charged',steps=steps,workers=1,GPU_use=False,GTscore='NOT_RUN'),indent=2)+'\n')
if status!='COMPLETE':raise SystemExit(1)
