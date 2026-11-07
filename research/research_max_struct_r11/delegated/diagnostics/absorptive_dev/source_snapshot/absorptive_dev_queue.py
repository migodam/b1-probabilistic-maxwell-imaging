"""Exactly8 fixed exposed absorptive DEV tasks, one local BLAS1 lane; no GT."""
from pathlib import Path
from datetime import datetime,timezone
import sys,json,subprocess,os,time,shutil
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];OUT=ROOT/'delegated/diagnostics/absorptive_dev'
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main():
 wiring=json.loads((OUT/'WIRING_CHECK.json').read_text())
 if wiring['status']!='PASS':raise RuntimeError('New complex adjoint/material wiring must pass before8run')
 if (OUT/'QUEUE.json').exists():raise FileExistsError('Preserve existing attempts; no restart/rescue')
 man=OUT/'TASKS.json';cfg=OUT/'SOLVER_CONFIG.json';tasks=json.loads(man.read_text())['tasks']
 if len(tasks)!=8 or len({t['case_id'] for t in tasks})!=4 or len({t['material_readout']['label'] for t in tasks})!=2:raise ValueError('Fixed4x2 manifest required')
 snap=OUT/'source_snapshot';paths=[ROOT/'src/absorptive_denominator.py',ROOT/'src/absorptive_adapter.py',ROOT/'src/cie_feedback_current.py',ROOT/'experiments/absorptive_dev_run.py',ROOT/'experiments/twod_run.py',ROOT/'src/twod_budget.py']
 for p in (man,cfg,ROOT/'experiments/absorptive_dev_queue.py',REPO/'research/research_max_struct_r10/src/lowrank_denominator.py'):shutil.copy2(p,snap/p.name)
 logs=OUT/'run_logs';logs.mkdir();state=dict(status='RUNNING',classification='FIXED_EXPOSED_DEV_PASSIVITY_MECHANISM_AUDIT_NOT_CONFIRMATION',requested=8,max_lanes=1,BLAS_threads=1,no_GPU_or_remote_cpu=True,truth_access='NONE',jobs=[]);wall=time.perf_counter()
 for index,t in enumerate(tasks):
  for p in paths:
   if p.read_bytes()!=(snap/p.name).read_bytes():raise RuntimeError('Frozen operator/runner source changed; do not silently continue')
  log=(logs/(t['task_id']+'.log')).open('w');beg=time.perf_counter();p=subprocess.Popen([sys.executable,str(ROOT/'experiments/absorptive_dev_run.py'),'--manifest',str(man),'--config',str(cfg),'--task-id',t['task_id']],cwd=REPO,stdout=log,stderr=subprocess.STDOUT)
  job=dict(task_id=t['task_id'],case_id=t['case_id'],method=t['material_readout']['label'],status='RUNNING',pid=p.pid,started_utc=datetime.now(timezone.utc).isoformat());state['jobs'].append(job);state['pending']=len(tasks)-index-1;dump(OUT/'QUEUE.json',state)
  pid,status,ru=os.wait4(p.pid,0);p.returncode=os.waitstatus_to_exitcode(status);log.close();job.update(status='FINISHED' if p.returncode==0 else 'FAILED',returncode=p.returncode,external_complete_child_CPU_seconds=ru.ru_utime+ru.ru_stime,user_CPU_seconds=ru.ru_utime,system_CPU_seconds=ru.ru_stime,wall_seconds=time.perf_counter()-beg,peak_RSS_bytes=ru.ru_maxrss,CPU_scope='Whole child external wait4 includes imports/setup/physics/output/exit; no scientific child spawned by this worker')
  print(json.dumps(dict(completed=index+1,total=8,task_id=t['task_id'],status=job['status'],CPU_seconds=job['external_complete_child_CPU_seconds'])),flush=True);dump(OUT/'QUEUE.json',state)
 state.update(status='COMPLETED' if all(j['status']=='FINISHED' for j in state['jobs']) else 'COMPLETED_WITH_ERRORS',pending=0,external_child_CPU_seconds_sum=sum(j['external_complete_child_CPU_seconds'] for j in state['jobs']),controller_own_CPU_seconds=time.process_time(),campaign_wall_seconds=time.perf_counter()-wall);dump(OUT/'QUEUE.json',state)
if __name__=='__main__':main()
