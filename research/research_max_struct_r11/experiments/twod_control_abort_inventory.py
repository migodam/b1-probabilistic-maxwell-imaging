"""Public-only inventory/snapshot after controller interruption; no solver or GT."""
import json,time,shutil,os,psutil
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];D=ROOT/'delegated/twod/independent48';out=D/'control_abort_20261107';out.mkdir(parents=True,exist_ok=False)
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=True)+'\n')
qdir=D/'expanded_queue';q=read(qdir/'RECEIPT.json');tmp=read(qdir/'RECEIPT.json.tmp');manifest=read(ROOT/'configs/twod_independent48_all_tasks.json');tasks={t['task_id']:t for t in manifest['tasks']};events=[json.loads(x) for x in (qdir/'EVENTS.jsonl').read_text().splitlines()];starts={r['task_id']:r for r in events if r['event']=='START'};ends={r['task_id']:r for r in events if r['event']=='END'};inflight=[starts[k] for k in starts if k not in ends]
ours=[]
for p in psutil.process_iter():
 if 'python' not in p.name().lower() or p.pid==os.getpid():continue
 args=p.cmdline()
 if any('MAX_STRUCT_R11_DEV_20261107' in x and 'twod_run.py' in x for x in args):ours.append(dict(pid=p.pid,cmd=args))
if ours:raise RuntimeError('Own physics workers still alive; do not duplicate: '+repr(ours))
for name in ['RECEIPT.json','RECEIPT.json.tmp','EVENTS.jsonl','FROZEN_MANIFEST.json','FROZEN_SOLVER.json','PLAN.json']:
 if (qdir/name).exists():shutil.copy2(qdir/name,out/name)
now=datetime.now(timezone.utc);rows=[]
for e in inflight:
 task=tasks[e['task_id']];path=REPO/task['output_dir'];r=read(path/'RECEIPT.json') if (path/'RECEIPT.json').exists() else {};timed=[];counts=[]
 if path.exists():shutil.copytree(path,out/'attempt01_outputs'/task['task_id'])
 if Path(e['log']).exists():(out/'logs').mkdir(exist_ok=True);shutil.copy2(e['log'],out/'logs'/Path(e['log']).name)
 for p in path.glob('stage_*_joint.json'):
  x=read(p)
  if isinstance(x.get('cpu_seconds'),(int,float)):timed.append(dict(path=str(p),CPU_lower_bound_seconds=x['cpu_seconds']))
  if x.get('physical_model_cost'):counts.append(x['physical_model_cost'])
 lower=max((x['CPU_lower_bound_seconds'] for x in timed),default=None);last_counts=max(counts,key=lambda x:x.get('LU',0)) if counts else None;wall=(now-datetime.fromisoformat(e['started_utc'])).total_seconds()
 rows.append(dict(task_id=task['task_id'],case_id=task['case_id'],attempt_id=task['task_id']+'__engineering_abort01',original_output=task['output_dir'],preserved_output=str((out/'attempt01_outputs'/task['task_id']).relative_to(REPO)),pid=e['pid'],cpu=e['cpu'],started_utc=e['started_utc'],confirmed_dead_by_utc=now.isoformat(),receipt_status=r.get('status','NO_FINAL_RECEIPT'),completed_valid=bool(r.get('status') in ('SEALED_PUBLIC_PREDICTION','CPU_BUDGET_REACHED') and r.get('cap_prediction_status')=='SEALED_PAID_INCUMBENT'),exact_full_child_CPU_seconds=None,CPU_lower_bound_seconds=lower,CPU_upper_observed_wall_bound_seconds=wall,upper_bound_scope='Conservative elapsed wall to same Windows-host dead observation, one-logical-CPU affinity; wall-clock stability assumed, not an exact exitCPU',lower_bound_source=timed,physical_counts_lower_bound=last_counts,cost_status='UNKNOWN_EXACT_NOT_ZERO',engineering_interrupted=True))
# Reproduce the relevant sharing mechanism with a new metadata file, not a scientific output.
target=out/'sharing_probe_target.txt';new=out/'sharing_probe_new.txt';target.write_text('old');new.write_text('new');handle=target.open('r')
try:
 try:new.replace(target);probe=dict(reproduced=False)
 except OSError as exc:probe=dict(reproduced=True,error_type=type(exc).__name__,winerror=getattr(exc,'winerror',None),error=str(exc))
finally:handle.close()
summary=dict(status='ORIGINALS_PRESERVED_PUBLIC_ONLY',observed_utc=now.isoformat(),finished_jobs_retained=5+63+len(q['jobs']),inflight_attempts=rows,own_live_physics=[],receipt_tmp_has_unpromoted_running=tmp['running'],original_failure_exact_stderr='UNAVAILABLE: local strictUTF8 decode failed at byte0xbe offset845 after SSH exit',original_failure_cause='Not conclusively recovered; unpromoted complete receipt.tmp plus no dispatcher continuation is consistent with mutable-file replace failure',sharing_probe=probe,GT_access='NONE',new_physics=0,CPU_seconds=time.process_time())
dump(out/'INVENTORY.json',summary);print(json.dumps(summary,ensure_ascii=True))
