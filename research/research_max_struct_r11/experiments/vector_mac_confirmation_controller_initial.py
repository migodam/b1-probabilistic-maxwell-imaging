"""Root-authorized finite Mac whole-object CPU queue, at most2 lanes, no GT."""
import os,time,json,sys,subprocess,resource,re,traceback
from pathlib import Path
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
 os.environ[name]='1'
ENTRY=time.perf_counter();ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
OUT=ROOT/'delegated/vector/mac_confirmation_controller';OUT.mkdir(parents=True,exist_ok=False)
def dump(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(x,indent=2)+'\n');tmp.replace(p)
def query(command):
 p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True);out=p.stdout.read();err=p.stderr.read();_,status,usage=os.wait4(p.pid,0);p.returncode=os.waitstatus_to_exitcode(status)
 return out,p.returncode,usage.ru_utime+usage.ru_stime
helper_CPU=0.
def available():
 global helper_CPU
 text,code,cpu=query(['/usr/bin/vm_stat']);helper_CPU+=cpu
 if code:raise RuntimeError('Read-only memory inventory unavailable')
 size=int(re.search(r'page size of (\d+) bytes',text).group(1));pages={}
 for line in text.splitlines()[1:]:
  if ':' in line:
   k,v=line.split(':',1)
   try:pages[k]=int(v.strip().rstrip('.'))
   except ValueError:pass
 return size*sum(pages.get(k,0) for k in ('Pages free','Pages inactive','Pages speculative'))
def rss(pid):
 global helper_CPU
 text,code,cpu=query(['/bin/ps','-o','rss=','-p',str(pid)]);helper_CPU+=cpu
 return int(text.strip())*1024 if not code and text.strip() else None
queue=json.loads((ROOT/'configs/vector_confirmation_queue.json').read_text());dispatch=json.loads((ROOT/'configs/vector_confirmation_host_dispatch.json').read_text())
mapping={t['task_id']:t['host'] for t in dispatch['tasks']};groups=[]
for i in range(0,50,5):
 group=queue['tasks'][i:i+5]
 if mapping[group[0]['task_id']]=='MAC_CPU':
  if len(group)!=5 or any(mapping[t['task_id']]!='MAC_CPU' for t in group):raise RuntimeError('Object host partition invalid')
  groups.append(dict(public_id=group[0]['public_id'],offset=i,tasks=[t['task_id'] for t in group]))
if sum(len(g['tasks']) for g in groups)!=25:raise RuntimeError('Mac must cover exactly25 unique tasks')
pending=groups.copy();active={};jobs=[];free_lanes=[0,1];peak_sum_RSS=0;last_monitor=-1e100;known=set();status='FAILED';failure=None
physical_text,_,cost=query(['/usr/sbin/sysctl','-n','hw.memsize']);helper_CPU+=cost;physical_RAM=int(physical_text.strip())
# No new numerical probe: second lane admission conservatively reserves12GiB
# of currently reclaimable memory for the800-request cache and workspace.
second_lane_available_required=12*1024**3;first_lane_available_required=3*1024**3
try:
 while pending or active:
  now=time.perf_counter()
  if now-last_monitor>=20:
   available_bytes=available();values={str(info['pid']):rss(info['pid']) for info in active.values()};summed=sum(v or 0 for v in values.values());peak_sum_RSS=max(peak_sum_RSS,summed);last_monitor=now
   with (OUT/'MEMORY_HISTORY.jsonl').open('a') as f:f.write(json.dumps(dict(elapsed_wall_seconds=now-ENTRY,available_free_inactive_speculative_bytes=available_bytes,active_RSS_bytes=values,active_RSS_sum=summed,active_lanes=len(active)))+'\n')
  while pending and free_lanes:
   needed=first_lane_available_required if not active else second_lane_available_required
   if available_bytes<needed:break
   group=pending.pop(0);lane=free_lanes.pop(0);dest='delegated/vector/mac_confirmation_processes/'+group['public_id']
   command=[sys.executable,str(ROOT/'experiments/vector_run.py'),'--config','configs/vector_confirmation_mac_execution.json','--queue','configs/vector_confirmation_queue.json','--offset',str(group['offset']),'--max-tasks','5','--output',dest]
   path=OUT/(group['public_id']+'.log');handle=path.open('x');child=subprocess.Popen(command,cwd=REPO,stdout=handle,stderr=subprocess.STDOUT)
   info=dict(group,lane=lane,pid=child.pid,command=command,started_wall_seconds=time.perf_counter()-ENTRY,queue_wait_wall_seconds=time.perf_counter()-ENTRY,output=dest,log=str(path.relative_to(REPO)))
   active[child]=(dict(info,pid=child.pid));active[child]['handle']=handle
   print(json.dumps(dict(event='MAC_OBJECT_STARTED',public_id=group['public_id'],lane=lane,pid=child.pid,tasks=5,available_bytes=available_bytes,no_fixed_P_core_affinity=True)),flush=True)
   # Re-read actual availability before admitting another lane.
   available_bytes=available()
  for child,info in list(active.items()):
   terminal=ROOT/info['output']/'TERMINALS.json'
   if terminal.exists():
    for row in json.loads(terminal.read_text())['terminals']:
     if row['task_id'] not in known:
      known.add(row['task_id']);r=json.loads((ROOT/'delegated/vector/inversions'/row['task_id']/'RECEIPT.json').read_text())
      print(json.dumps(dict(event='MAC_TASK_SEALED',completed=len(known),requested=25,task_id=row['task_id'],CPU_seconds=r['CPU_seconds'],wall_seconds=r['wall_seconds'],physical_LU=r['counts']['physical_LU'],Phi=r['physical_Phi'],truth_scoring='NOT_RUN')),flush=True)
   pid,code,usage=os.wait4(child.pid,os.WNOHANG)
   if not pid:continue
   child.returncode=os.waitstatus_to_exitcode(code);info.pop('handle').close();row={k:v for k,v in info.items() if k!='handle'};row.update(exit=child.returncode,whole_child_CPU_seconds=usage.ru_utime+usage.ru_stime,whole_child_wall_seconds=time.perf_counter()-ENTRY-info['started_wall_seconds'],peak_RSS_bytes=usage.ru_maxrss,CPU_scope='Mac wait4 whole scientific child replaces nested process/arm CPU',no_fixed_P_core_affinity=True)
   jobs.append(row);free_lanes.append(info['lane']);del active[child];dump(OUT/'PROCESS_RECEIPTS.json',jobs);print(json.dumps(dict(event='MAC_OBJECT_TERMINAL',**row)),flush=True)
   if child.returncode:raise RuntimeError('Mac object failed; no automatic scientific budget or GT changes')
  dump(OUT/'STATUS.json',dict(status='RUNNING' if pending or active else 'COMPLETE',completed_tasks=len(known),global_frozen_tasks=50,host_assigned_tasks=25,pending=[g['public_id'] for g in pending],active=[{k:v for k,v in a.items() if k!='handle'} for a in active.values()],jobs=jobs,physical_RAM_bytes=physical_RAM,second_lane_available_required_bytes=second_lane_available_required,peak_sum_RSS_observed_bytes=peak_sum_RSS,max_lanes=2,truth_access=False))
  if pending or active:time.sleep(5)
 status='COMPLETE'
except BaseException as exc:failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());print(json.dumps(failure),flush=True)
finally:
 dump(OUT/'CONTROLLER_COST.json',dict(status=status,failure=failure,CPU_seconds=time.process_time(),monitor_helper_CPU_seconds=helper_CPU,wall_seconds=time.perf_counter()-ENTRY,jobs=jobs,max_lanes=2,peak_sum_RSS_observed_bytes=peak_sum_RSS,scope='Own finite Mac controller+read-only memory helpers; scientific wait4 children separately counted; wall overlaps children',no_GPU=True,no_GT_access=True))
raise SystemExit(0 if status=='COMPLETE' else 1)
