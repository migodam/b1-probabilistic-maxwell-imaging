"""Finite two-lane triggered continuation queue with exact external actual CPU."""
from pathlib import Path
from datetime import datetime,timezone
import sys,json,subprocess,os,time,shutil
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];OUT=ROOT/'delegated/diagnostics/bothblocks_confirmation/posthoc_budget_control'
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main():
 manifest=json.loads((OUT/'FROZEN_MANIFEST.json').read_text());pending=[x for x in manifest['pairs'] if x['trigger']]
 if (OUT/'QUEUE.json').exists():raise FileExistsError('Preserve attempts; no duplicate restart')
 logdir=OUT/'run_logs';logdir.mkdir();snap=OUT/'source_snapshot';snap.mkdir()
 for p in (ROOT/'experiments/diagnose_budget_control_run.py',ROOT/'experiments/diagnose_budget_control_queue.py',ROOT/'src/twod_budget.py'):shutil.copy2(p,snap/p.name)
 shutil.copy2(REPO/'research/research_max_struct_r6/src/benchmark_solver.py',snap/'r6_benchmark_solver.py')
 state=dict(status='RUNNING',classification='POSTHOC_BUDGET_CONTROL',requested=len(pending),jobs=[],pending=len(pending),truth_access='NONE');live={};wall=time.perf_counter()
 def launch(r):
  handle=(logdir/(r['pair_id']+'.log')).open('w');p=subprocess.Popen([sys.executable,str(ROOT/'experiments/diagnose_budget_control_run.py'),'--pair-id',r['pair_id']],cwd=REPO,stdout=handle,stderr=subprocess.STDOUT);j=dict(pair_id=r['pair_id'],status='RUNNING',pid=p.pid,started_utc=datetime.now(timezone.utc).isoformat());state['jobs'].append(j);live[p.pid]=(p,handle,j,r,time.perf_counter());state['pending']=len(pending);dump(OUT/'QUEUE.json',state)
 while pending and len(live)<2:launch(pending.pop(0))
 while live:
  pid,status,ru=os.wait4(-1,0);p,h,j,r,beg=live.pop(pid);p.returncode=os.waitstatus_to_exitcode(status);h.close();actual=ru.ru_utime+ru.ru_stime;j.update(status='FINISHED' if p.returncode==0 else 'FAILED',returncode=p.returncode,new_child_actual_CPU_seconds=actual,total_actual_CPU_seconds=r['original_actual_child_CPU_seconds']+actual,target_CPU_seconds=r['target_G19_actual_CPU_seconds'],unused_CPU_seconds=max(0.,r['target_G19_actual_CPU_seconds']-r['original_actual_child_CPU_seconds']-actual),charged_overrun_CPU_seconds=max(0.,r['original_actual_child_CPU_seconds']+actual-r['target_G19_actual_CPU_seconds']),wall_seconds=time.perf_counter()-beg,CPU_scope='Exact external wait4 includes imports/setup/serialization/exit; atomic overrun charged')
  print(json.dumps(j),flush=True)
  if pending:launch(pending.pop(0))
  dump(OUT/'QUEUE.json',state)
 state.update(status='COMPLETED' if all(j['status']=='FINISHED' for j in state['jobs']) else 'COMPLETED_WITH_ERRORS',pending=0,new_child_actual_CPU_seconds_sum=sum(j['new_child_actual_CPU_seconds'] for j in state['jobs']),controller_CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall);dump(OUT/'QUEUE.json',state)
if __name__=='__main__':main()
