"""Finite foreground campaign monitor; no dispatch/retry, no GT scoring."""
import json,time,subprocess,sys
from pathlib import Path
import twod_dev_transport as t
D=t.ROOT/'delegated/twod/independent48'
def main():
 previous=-1;previous_failed=0;last_synced=5
 while True:
  code="import json,pathlib; r=pathlib.Path('"+t.REMOTE+"/research/research_max_struct_r11'); q=[json.loads((r/p).read_text()) for p in ['delegated/twod/independent48/queue/RECEIPT.json','delegated/twod/independent48/interleaved_queue/RECEIPT.json']]; jobs=[j for x in q for j in x['jobs']]; print(json.dumps(dict(status=q[-1]['status'],completed=len(jobs),failed=sum(j['status']!='FINISHED' for j in jobs),running=q[-1]['running'],pending=q[-1].get('pending'),CPU_seconds=sum(j.get('child_process_CPU_seconds') or 0 for j in jobs),task_ids=[j['task_id'] for j in jobs],failed_jobs=[j for j in jobs if j['status']!='FINISHED'])))"
  try:x=json.loads(t.ssh('D:/python/python.exe -c "'+code+'"','independent48_monitor_status'))
  except Exception as exc:
   print(json.dumps(dict(event='MONITOR_METADATA_FETCH_FAILED',error=str(exc),physics_interrupted=False)),flush=True);time.sleep(60);continue
  if len(x['task_ids'])!=len(set(x['task_ids'])):raise RuntimeError('Duplicate task charging; audit required')
  x['observed_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());(t.OUT/'INDEPENDENT48_STATUS.json').write_text(json.dumps(x,indent=2)+'\n')
  if x['completed']!=previous or x['failed']!=previous_failed:
   print(json.dumps({k:v for k,v in x.items() if k!='task_ids'}),flush=True);previous=x['completed']
  if x['failed']>previous_failed:
   previous_failed=x['failed'];(D/'FAILED_JOB_ALERT.json').write_text(json.dumps(x['failed_jobs'],indent=2)+'\n')
   print(json.dumps(dict(event='NEW_FAILED_ATTEMPT_RETAINED',automatic_retry=False,failed_jobs=x['failed_jobs'])),flush=True)
  if x['completed']-last_synced>=12:
   subprocess.run([sys.executable,str(t.ROOT/'experiments/twod_independent_sync.py'),'--finished'],check=True);last_synced=x['completed']
   subprocess.run([sys.executable,str(t.ROOT/'experiments/twod_independent_costs.py')],check=True)
  if x['status'] in ('COMPLETED','COMPLETED_WITH_ERRORS') and not x['running']:
   if x['completed']!=288:raise RuntimeError('Terminal queue missing frozen tasks')
   subprocess.run([sys.executable,str(t.ROOT/'experiments/twod_independent_sync.py'),'--finished','--include-sealed-truth'],check=True)
   subprocess.run([sys.executable,str(t.ROOT/'experiments/twod_independent_costs.py')],check=True)
   print(json.dumps(dict(event='ALL288_TERMINAL_SYNCED',GT_scoring='NOT_RUN',failed=x['failed'])),flush=True);return
  time.sleep(60)
if __name__=='__main__':main()
