"""Public-only independent48 actual-cost compilation, partial or terminal."""
import json,csv,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];D=ROOT/'delegated/twod/independent48';T=ROOT/'delegated/twod/dev_transport'
def read(p):return json.loads(Path(p).read_text())
def main():
 D.mkdir(parents=True,exist_ok=True);jobs={};controllers=[];queues=[]
 for name in ('queue','interleaved_queue'):
  p=D/name/'RECEIPT.json'
  if not p.exists():continue
  q=read(p);queues.append(dict(name=name,status=q['status'],pending=q.get('pending'),running=q.get('running'),completed=len(q['jobs'])))
  controllers.append(dict(queue=name,CPU_seconds=q['controller_CPU_seconds']))
  for j in q['jobs']:
   if j['task_id'] in jobs:raise ValueError('Do not double charge duplicate task attempts without explicit attempt schema')
   jobs[j['task_id']]=dict(queue=name,task_id=j['task_id'],status=j['status'],child_status=j.get('child_status'),CPU_seconds=j.get('child_process_CPU_seconds'),LU=(j.get('physical_model_cost') or {}).get('LU'),complex_RHS=(j.get('physical_model_cost') or {}).get('actual_complex_RHS'))
 transport=[]
 if (T/'TRANSPORT_COSTS.jsonl').exists():
  after_dev=False
  for line in (T/'TRANSPORT_COSTS.jsonl').read_text().splitlines():
   r=json.loads(line)
   if after_dev:transport.append(r)
   if r['label']=='fetch_outputs':after_dev=True
 generation=read(D/'GENERATION_RECEIPT.json') if (D/'GENERATION_RECEIPT.json').exists() else {}
 setup=read(T/'INDEPENDENT48_SOURCE_PACKAGE_RECEIPT.json') if (T/'INDEPENDENT48_SOURCE_PACKAGE_RECEIPT.json').exists() else {}
 sync_receipts=[read(p) for p in D.glob('SYNC_*.json')]
 totals=dict(remote_sync_packaging_CPU_seconds=sum(r['remote_packaging_CPU_seconds'] for r in sync_receipts),local_sync_extraction_CPU_seconds=sum(r['local_extraction_CPU_seconds'] for r in sync_receipts),inversion_child_CPU_seconds=sum(j['CPU_seconds'] or 0 for j in jobs.values()),queue_controller_CPU_seconds=sum(c['CPU_seconds'] for c in controllers),generation_internal_whole_child_CPU_lower_bound_seconds=generation.get('whole_child_CPU_seconds',0),source_packaging_CPU_seconds=setup.get('CPU_seconds',0),local_transport_controller_CPU_seconds=sum(r['controller_CPU_seconds'] for r in transport),local_transport_child_CPU_seconds=sum(r['child_CPU_seconds'] for r in transport))
 status='ALL288_TERMINAL' if len(jobs)==288 and all(not q['running'] for q in queues) and queues[-1]['status'] in ('COMPLETED','COMPLETED_WITH_ERRORS') else 'PARTIAL'
 doc=dict(status=status,queues=queues,tasks_terminal=len(jobs),failed=sum(j['status']!='FINISHED' for j in jobs.values()),totals=totals,known_CPU_lower_bound_seconds=sum(totals.values()),jobs=list(jobs.values()),controllers=controllers,transport=transport,truth_access='NONE',GT_scoring='NOT_RUN',cost_rule='Disjoint actual GetProcessTimes jobs plus each controller lifetime once. Internal objective/aux/corrector timers overlap and are not additive. Explicit STOP_DISPATCH return1 is a controller handoff, not inversion failure.',unknown_costs=['Generation final serialization/process-exit CPU beyond final absolute self timestamp','Remote source extraction, status/probe, final packaging and local extraction if no separate receipt','Metadata preparation/compiler final tail CPU'])
 (D/'COSTS.json').write_text(json.dumps(doc,indent=2)+'\n')
 if jobs:
  with (D/'INVERSION_COSTS.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(next(iter(jobs.values()))));w.writeheader();w.writerows(jobs.values())
 print(json.dumps(dict(status=status,tasks=len(jobs),failed=doc['failed'],known_CPU_lower_bound_seconds=doc['known_CPU_lower_bound_seconds'])))
if __name__=='__main__':main()
