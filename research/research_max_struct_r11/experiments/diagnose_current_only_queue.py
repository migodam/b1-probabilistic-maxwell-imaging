"""Finite local two-lane queue with exact wait4 child CPU receipts; no truth access."""
from pathlib import Path
import sys,json,time,subprocess,os,shutil
from datetime import datetime,timezone
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];OUT=ROOT/'delegated/diagnostics/current_only_dev';LOG=OUT/'run_logs'

def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main():
    confirmation=json.loads((ROOT/'delegated/diagnostics/bothblocks_confirmation/QUEUE.json').read_text())
    if confirmation['status'] not in ('COMPLETED','COMPLETED_WITH_ERRORS') or len(confirmation['jobs'])!=80:raise RuntimeError('Do not compete with frozen80 confirmation queue')
    LOG.mkdir(parents=True,exist_ok=True)
    man=ROOT/'configs/twod_dev_current_only_tasks.json';cfg=ROOT/'configs/twod_dev_solver.json';tasks=json.loads(man.read_text())['tasks']
    if (OUT/'CURRENT_ONLY_QUEUE.json').exists():raise FileExistsError('Existing queue audit')
    for src in (ROOT/'src/diagnostic_green_current_only.py',ROOT/'src/material_readouts.py',ROOT/'experiments/diagnose_current_only_run.py',cfg,man):shutil.copy2(src,LOG/src.name)
    state=dict(status='RUNNING',truth_access='NONE',max_lanes=2,host='local Mac CPU BLAS1',jobs=[]);live={};pending=list(tasks);start=time.perf_counter()
    def launch(t):
      log=LOG/(t['case_id']+'.log');handle=log.open('w')
      cmd=[sys.executable,str(ROOT/'experiments/diagnose_current_only_run.py'),'--manifest',str(man),'--config',str(cfg),'--task-id',t['task_id']]
      p=subprocess.Popen(cmd,cwd=REPO,stdout=handle,stderr=subprocess.STDOUT)
      job=dict(task_id=t['task_id'],case_id=t['case_id'],status='RUNNING',pid=p.pid,started_utc=datetime.now(timezone.utc).isoformat(),log=str(log))
      state['jobs'].append(job);live[p.pid]=(p,handle,job,time.perf_counter());dump(OUT/'CURRENT_ONLY_QUEUE.json',state)
    while pending and len(live)<2:launch(pending.pop(0))
    while live:
      pid,status,usage=os.wait4(-1,0)
      p,handle,job,beg=live.pop(pid);p.returncode=os.waitstatus_to_exitcode(status);handle.close()
      job.update(status='FINISHED' if p.returncode==0 else 'FAILED',returncode=p.returncode,child_process_CPU_seconds=usage.ru_utime+usage.ru_stime,user_CPU_seconds=usage.ru_utime,system_CPU_seconds=usage.ru_stime,
        CPU_scope='Exact external wait4 rusage including imports, setup, inversion, serialization and exit',wall_seconds=time.perf_counter()-beg,max_rss=usage.ru_maxrss)
      print(json.dumps(job),flush=True)
      if pending:launch(pending.pop(0))
      dump(OUT/'CURRENT_ONLY_QUEUE.json',state)
    state.update(status='COMPLETED' if all(j['status']=='FINISHED' for j in state['jobs']) else 'COMPLETED_WITH_ERRORS',wall_seconds=time.perf_counter()-start,
      exact_child_CPU_seconds_sum=sum(j['child_process_CPU_seconds'] for j in state['jobs']),controller_CPU_seconds=time.process_time())
    dump(OUT/'CURRENT_ONLY_QUEUE.json',state)
if __name__=='__main__':main()
