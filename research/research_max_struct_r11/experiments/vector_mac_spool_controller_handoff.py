"""Finite controller-only handoff after c11 whole-case natural completion."""
import os,time,json,sys,subprocess,resource,signal,traceback
from pathlib import Path
ENTRY=time.perf_counter();ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
OUT=ROOT/'delegated/vector/mac_spool_controller_handoff';OUT.mkdir(parents=True,exist_ok=False)
OLD=ROOT/'delegated/vector/mac_confirmation_controller';helper_CPU=0.;status='FAILED';failure=None

def dump(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')
def query(command):
 global helper_CPU
 p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True);text=p.stdout.read();error=p.stderr.read();_,code,r=os.wait4(p.pid,0);p.returncode=os.waitstatus_to_exitcode(code);helper_CPU+=r.ru_utime+r.ru_stime
 return text.strip(),p.returncode
try:
 active=json.loads((OLD/'STATUS.json').read_text())['active'];first=next(a for a in active if a['public_id']=='c11_matched12')
 parent_text,code=query(['/bin/ps','-o','ppid=','-p',str(first['pid'])]);parent=int(parent_text)
 command,_=query(['/bin/ps','-o','command=','-p',str(parent)])
 if 'vector_mac_confirmation_controller.py' not in command:raise RuntimeError('Refuse non-owned controller PID')
 dump(OUT/'PLAN.json',dict(owned_controller_pid=parent,current_c11_worker_pid=first['pid'],current_worker_never_signalled=True,
  boundary='Only after c11 five-arm normal exit0 and parent wait4 record, with no active workers',
  authority='Root Mac25/max2-lane authorization plus future-Mac thin storage spool; resource-policy implementation only',
  second_lane_available_GiB=4,normal_spool_archive_check_required=True,scientific_tasks_precision_budgets_unchanged=True))
 while True:
  path=OLD/'PROCESS_RECEIPTS.json'
  if path.exists():
   jobs=json.loads(path.read_text());job=next((j for j in jobs if j['public_id']=='c11_matched12'),None)
   if job:
    if job['exit']!=0:raise RuntimeError('c11 not normal exit; no controller handoff')
    break
  if (OLD/'CONTROLLER_COST.json').exists():raise RuntimeError('Old controller already ended unexpectedly; no handoff')
  time.sleep(1)
 other=[]
 for group in ('c14_matched12','c15_matched12','c18_matched12','c18_cross18to12'):
  if (ROOT/'delegated/vector/mac_confirmation_processes'/group).exists():other.append(group)
 if other:
  status='HANDOFF_CANCELLED_OLD_CONTROLLER_ALREADY_STARTED_OTHER_OBJECTS'
  dump(OUT/'CANCELLED.json',dict(reason='Do not interrupt any active worker; initial controller continues',objects=other))
 else:
  for task in json.loads((ROOT/'configs/vector_confirmation_queue.json').read_text())['tasks'][:5]:
   d=ROOT/'delegated/vector/inversions'/task['task_id']
   if not all((d/n).exists() for n in ('SEALED.json','COST_PREFIX_SEAL.json','prediction.npz')):raise RuntimeError('c11 full scientific seal barrier missing')
  children,code=query(['/usr/bin/pgrep','-P',str(parent)])
  if children:
   status='HANDOFF_CANCELLED_ACTIVE_CHILD_PRESERVED';dump(OUT/'CANCELLED.json',dict(reason='No signals while any child remains',children=children))
  else:
   os.kill(parent,signal.SIGINT)
   while not (OLD/'CONTROLLER_COST.json').exists():time.sleep(.1)
   old=json.loads((OLD/'CONTROLLER_COST.json').read_text())
   dump(OUT/'ENGINEERING_HANDOFF_CLASSIFICATION.json',dict(status='CONTROLLED_PARENT_DRAIN_HANDOFF_NOT_SCIENCE_FAILURE',
    old_controller_raw_status=old['status'],old_failure=old['failure'],old_cost_path=str(OLD/'CONTROLLER_COST.json'),
    preserved_completed_object='c11_matched12,all5natural terminal',worker_signals_sent=0,root_50_partition_unchanged=True))
   cmd=[sys.executable,str(ROOT/'experiments/vector_mac_confirmation_controller.py'),'--resume','--output-tag','mac_confirmation_spooled_controller','--second-lane-available-gib','4','--require-spool-first-archive']
   with (OUT/'PHASE2.log').open('x') as f:
    child=subprocess.Popen(cmd,cwd=REPO,stdout=f,stderr=subprocess.STDOUT);dump(OUT/'ACTIVE.json',dict(pid=child.pid,command=cmd,started_unix=time.time()))
    _,code,usage=os.wait4(child.pid,0);child.returncode=os.waitstatus_to_exitcode(code)
   dump(OUT/'PHASE2_INCLUSIVE_WAIT4_VIEW.json',dict(exit=child.returncode,CPU_seconds=usage.ru_utime+usage.ru_stime,
    scope='Includes controller and all its children per Macwait4; OVERLAPPING view, do not add to per-object scientific wait4+controller own CPU',additive=False))
   if child.returncode:raise RuntimeError('Resumed CPU resource controller failed; preserve evidence')
   status='COMPLETE'
except BaseException as exc:failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
finally:
 dump(OUT/'PROCESS_COST.json',dict(status=status,failure=failure,CPU_seconds=time.process_time(),query_helper_CPU_seconds=helper_CPU,
  wall_seconds=time.perf_counter()-ENTRY,scope='Own finite handoff helper only; phase2wait4 inclusive is overlapping, scientific objects separately accounted',
  no_GT_access=True,no_new_physical_probe=True,no_current_worker_interruption=True))
print(json.dumps(dict(status=status,failure=failure)),flush=True)
raise SystemExit(0 if status=='COMPLETE' or status.startswith('HANDOFF_CANCELLED') else 1)
