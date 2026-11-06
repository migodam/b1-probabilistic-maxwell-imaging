"""Frozen case-cyclic local two-lane 80-task queue; no scorer or truth reads."""
from pathlib import Path
import sys,json,time,subprocess,os,shutil,argparse
from datetime import datetime,timezone
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];OUT=ROOT/'delegated/diagnostics/bothblocks_confirmation';LOG=OUT/'run_logs';SNAP=OUT/'frozen_source'
METHODS=('cie_single','cie_dual','global_dual','fullG19','fullGreenBothBlocks')
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--execute',action='store_true');a=ap.parse_args()
 frozen=json.loads((OUT/'FROZEN_CONTRACT.json').read_text())
 if not a.execute or frozen['status']!='ROOT_FROZEN_BEFORE_GENERATION' or not frozen['execution_authorized']:raise ValueError('Root execution freeze required')
 if (OUT/'QUEUE.json').exists():raise FileExistsError('Queue already exists')
 LOG.mkdir(parents=True,exist_ok=True);SNAP.mkdir(parents=True,exist_ok=False)
 paths=[ROOT/'src/cie_feedback_current.py',ROOT/'src/material_readouts.py',ROOT/'src/twod_budget.py',ROOT/'experiments/twod_run.py',ROOT/'experiments/diagnose_bothblocks_run.py',ROOT/'experiments/diagnose_bothblocks16_generate.py',ROOT/'experiments/diagnose_bothblocks16_queue.py',ROOT/'configs/twod_bothblocks16_generator.json',ROOT/'configs/twod_bothblocks16_solver.json',ROOT/'configs/twod_bothblocks16_scoring.json']
 for p in paths:shutil.copy2(p,SNAP/p.name)
 state=dict(status='GENERATING',truth_access='NONE_CONTROLLER; generator-only private truth construction',max_lanes=2,host='local Mac BLAS1',dispatch='case-index cyclic offset of five frozen arms',jobs=[]);dump(OUT/'QUEUE.json',state);wall=time.perf_counter()
 glog=(OUT/'generation.log').open('w');gp=subprocess.Popen([sys.executable,str(ROOT/'experiments/diagnose_bothblocks16_generate.py'),'--config',str(ROOT/'configs/twod_bothblocks16_generator.json'),'--execute'],cwd=REPO,stdout=glog,stderr=subprocess.STDOUT)
 pid,status,usage=os.wait4(gp.pid,0);gp.returncode=os.waitstatus_to_exitcode(status);glog.close()
 state['generation_job']=dict(returncode=gp.returncode,external_child_CPU_seconds=usage.ru_utime+usage.ru_stime,CPU_scope='Exact external wait4 complete child rusage',max_rss=usage.ru_maxrss)
 if gp.returncode!=0:state['status']='GENERATION_FAILED';dump(OUT/'QUEUE.json',state);raise RuntimeError('Generator failed; preserved outputs; no reroll')
 man=ROOT/'configs/twod_bothblocks16_all_tasks.json';cfg=ROOT/'configs/twod_bothblocks16_solver.json';tasks=json.loads(man.read_text())['tasks'];shutil.copy2(man,SNAP/man.name)
 if len(tasks)!=80 or len({t['case_id'] for t in tasks})!=16:raise ValueError('80-task balance differs')
 for cid in sorted({t['case_id'] for t in tasks}):
  ts=[t for t in tasks if t['case_id']==cid];offset=int(cid[-4:])%5
  if [t['method'] for t in ts]!=list(METHODS[offset:]+METHODS[:offset]):raise ValueError('Cyclic dispatch differs')
 state.update(status='RUNNING',pending=80);pending=list(tasks);live={}
 def launch(t):
  for path in paths:
   if path.read_bytes()!=(SNAP/path.name).read_bytes():raise RuntimeError('Frozen source/config changed: '+str(path))
  both=t['method']=='fullGreenBothBlocks';runner=ROOT/'experiments'/('diagnose_bothblocks_run.py' if both else 'twod_run.py')
  if t['runner']!=runner.name:raise ValueError('Arm runner mismatch')
  log=LOG/(t['task_id']+'.log');handle=log.open('w');cmd=[sys.executable,str(runner),'--manifest',str(man),'--config',str(cfg),'--task-id',t['task_id']]
  p=subprocess.Popen(cmd,cwd=REPO,stdout=handle,stderr=subprocess.STDOUT)
  job=dict(task_id=t['task_id'],case_id=t['case_id'],method=t['method'],runner=runner.name,status='RUNNING',pid=p.pid,started_utc=datetime.now(timezone.utc).isoformat(),log=str(log))
  state['jobs'].append(job);live[p.pid]=(p,handle,job,time.perf_counter());state['pending']=len(pending);dump(OUT/'QUEUE.json',state)
 while pending and len(live)<2:launch(pending.pop(0))
 while live:
  pid,status,usage=os.wait4(-1,0);p,handle,job,beg=live.pop(pid);p.returncode=os.waitstatus_to_exitcode(status);handle.close()
  job.update(status='FINISHED' if p.returncode==0 else 'FAILED',returncode=p.returncode,child_process_CPU_seconds=usage.ru_utime+usage.ru_stime,user_CPU_seconds=usage.ru_utime,system_CPU_seconds=usage.ru_stime,CPU_scope='Exact external wait4 complete child rusage',wall_seconds=time.perf_counter()-beg,max_rss=usage.ru_maxrss)
  done=sum(x['status']!='RUNNING' for x in state['jobs']);print(json.dumps(dict(completed=done,total=80,task_id=job['task_id'],status=job['status'],CPU_seconds=job['child_process_CPU_seconds'])),flush=True)
  if pending:launch(pending.pop(0))
  dump(OUT/'QUEUE.json',state)
 state.update(status='COMPLETED' if all(j['status']=='FINISHED' for j in state['jobs']) else 'COMPLETED_WITH_ERRORS',pending=0,wall_seconds=time.perf_counter()-wall,exact_inversion_child_CPU_seconds_sum=sum(j['child_process_CPU_seconds'] for j in state['jobs']),controller_CPU_seconds=time.process_time(),truth_access='NONE_CONTROLLER; no scoring before all80 terminal')
 dump(OUT/'QUEUE.json',state)
if __name__=='__main__':main()
