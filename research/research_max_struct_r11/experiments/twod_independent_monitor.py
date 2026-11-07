"""Finite campaign monitor with explicit attempt accounting; never GT scoring."""
import json,time,subprocess,sys,os
import twod_dev_transport as t
from twod_recovery import active_task
D=t.ROOT/'delegated/twod/independent48'
def main():
 (D/'ACTIVE_METADATA_MONITOR.json').write_text(json.dumps(dict(pid=os.getpid(),started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),scope='LOCAL_METADATA_ONLY_RECOVERY_AWARE'),indent=2)+'\n')
 previous=-1;previous_failed=0
 master=json.loads((t.ROOT/'configs/twod_independent48_all_tasks.json').read_text())['tasks']
 last_synced=sum((t.REPO/active_task(task)['output_dir']/'RECEIPT.json').exists() for task in master)
 study=t.REMOTE+'/research/research_max_struct_r11'
 while True:
  code="import json,sys,time; sys.path.insert(0,'"+study+"/experiments'); from twod_recovery import queue_paths,active_job,registry; qs=[json.loads(p.read_text()) for p in queue_paths()]; raw=[j for q in qs for j in q['jobs']]; jobs=[j for j in raw if active_job(j)]; reg=registry(); print(json.dumps(dict(status=qs[-1]['status'],completed=len(jobs),failed=sum(j['status']!='FINISHED' for j in jobs),running=qs[-1]['running'],pending=qs[-1].get('pending'),CPU_seconds=sum(j.get('child_process_CPU_seconds') or 0 for j in jobs),task_ids=[j['task_id'] for j in jobs],failed_jobs=[j for j in jobs if j['status']!='FINISHED'],engineering_interrupted_attempts=len(reg.get('interrupted_attempts',[])),zero_physics_startup_failures=len(reg.get('zero_physics_retry02_failures',[])),final_retry03_started=any('final_retry03_queue' in str(p) for p in queue_paths()))))"
  try:x=json.loads(t.ssh('D:/python/python.exe -c "'+code+'"','independent48_monitor_recovery_aware'))
  except Exception as exc:print(json.dumps(dict(event='METADATA_FETCH_FAILED',error=str(exc),physics_interrupted=False)),flush=True);time.sleep(60);continue
  if len(x['task_ids'])!=len(set(x['task_ids'])):raise RuntimeError('Duplicate designated scientific attempt')
  x['observed_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());(t.OUT/'INDEPENDENT48_STATUS.json').write_text(json.dumps(x,indent=2)+'\n')
  if x['completed']!=previous or x['failed']!=previous_failed:
   print(json.dumps({k:v for k,v in x.items() if k!='task_ids'}),flush=True);previous=x['completed']
  if x['failed']>previous_failed:previous_failed=x['failed'];(D/'FAILED_JOB_ALERT.json').write_text(json.dumps(x['failed_jobs'],indent=2)+'\n')
  if x['completed']-last_synced>=12:
   subprocess.run([sys.executable,str(t.ROOT/'experiments/twod_independent_sync.py'),'--finished'],check=True);last_synced=x['completed'];subprocess.run([sys.executable,str(t.ROOT/'experiments/twod_independent_costs.py')],check=True)
  if x['status'] in ('COMPLETED','COMPLETED_WITH_ERRORS') and not x['running']:
   if x['completed']==285 and not x['final_retry03_started']:
    print(json.dumps(dict(event='START_THREE_REGISTERED_ZERO_PHYSICS_RETRY03_ONLY',new_materials=0,GT_access='NONE')),flush=True)
    cmd='D:/python/python.exe '+study+'/experiments/twod_queue.py --execute --receipt-dir '+study+'/delegated/twod/independent48/final_retry03_queue --vendor-path '+t.REMOTE+'/vendor --manifest '+study+'/configs/twod_independent48_final_retry03_tasks.json --config '+study+'/configs/twod_independent48_solver.json --cpus 0 2 8 10 > '+study+'/delegated/twod/independent48/final_retry03_stdout.log 2> '+study+'/delegated/twod/independent48/final_retry03_stderr.log'
    t.ssh(cmd,'run_three_registered_zero_physics_retry03');continue
   if x['completed']!=288:raise RuntimeError('Terminal campaign missing designated scientific tasks')
   subprocess.run([sys.executable,str(t.ROOT/'experiments/twod_independent_sync.py'),'--finished','--include-sealed-truth'],check=True);subprocess.run([sys.executable,str(t.ROOT/'experiments/twod_independent_costs.py')],check=True)
   print(json.dumps(dict(event='ALL288_DESIGNATED_OUTCOMES_SYNCED',engineering_interrupted=4,zero_physics_startup_failures=3,GT_scoring='NOT_RUN')),flush=True);return
  time.sleep(60)
if __name__=='__main__':main()
