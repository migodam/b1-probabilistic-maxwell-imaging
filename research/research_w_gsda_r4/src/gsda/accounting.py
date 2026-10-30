"""Fresh R4 pre-action physics accounting; wrapper CPU is authoritative."""
from pathlib import Path
from contextlib import contextmanager
import json,time,os,fcntl,resource,sys,datetime
ROOT=Path(__file__).resolve().parents[2]
class BudgetExceeded(RuntimeError):pass
def plain(x):
 if hasattr(x,'tolist'):return x.tolist()
 if hasattr(x,'item'):return x.item()
 if isinstance(x,Path):return str(x)
 if isinstance(x,complex):return [x.real,x.imag]
 raise TypeError(type(x).__name__)
def dump(path,obj):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(obj,indent=2,default=plain,allow_nan=False)+'\n')
@contextmanager
def fem_lock(root=ROOT):
 root=Path(root);root.mkdir(parents=True,exist_ok=True)
 with (root/'.fem.lock').open('a') as f:
  fcntl.flock(f,fcntl.LOCK_EX)
  try:yield
  finally:fcntl.flock(f,fcntl.LOCK_UN)
class BudgetBook:
 def __init__(self,root=ROOT,job=None):
  self.root=Path(root);self.path=self.root/'PHYSICS_LEDGER.jsonl';self.clock=json.loads((self.root/'RUN_CLOCK.json').read_text());self.job=job or os.environ.get('GSDA_JOB','unspecified');self.cpu0=time.process_time();self.wall0=time.perf_counter()
 def check(self,new_job=False):
  if (self.root/'RUN_CLOSED.json').exists():raise BudgetExceeded('R4 campaign is closed; new deployment requires a separate explicit budget')
  clock=self.clock;now=time.time();deadline=datetime.datetime.fromisoformat(clock['deadline_utc']).timestamp();heavy=datetime.datetime.fromisoformat(clock['no_new_heavy_after_utc']).timestamp()
  if now>deadline:raise BudgetExceeded('R4 180 minute wall deadline')
  if new_job and now>heavy:raise BudgetExceeded('R4 no-new-heavy cutoff')
  rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
  if rss>clock['single_job_rss_bytes']:raise BudgetExceeded('R4 4GiB job RSS cap')
  return rss
 def append(self,row,path=None):
  with (self.root/'.ledger.lock').open('a') as lock:
   fcntl.flock(lock,fcntl.LOCK_EX)
   with (Path(path) if path else self.path).open('a') as f:f.write(json.dumps(dict(time_unix=time.time(),pid=os.getpid(),job=self.job,**row),default=plain,allow_nan=False)+'\n')
 def reserve(self,phase,label,rhs=0,factorizations=0,**meta):
  self.check();self.append(dict(kind='reservation',phase=phase,label=label,rhs=int(rhs),factorizations=int(factorizations),**meta))
 def event(self,phase='R4',label='',**meta):
  if 'kind' in meta:meta['event_kind']=meta.pop('kind')
  self.append(dict(kind='event',phase=phase,label=label,**meta))
 def failure(self,phase,label,error,**meta):
  r=dict(kind='failure',phase=phase,label=label,error=str(error),**meta);self.append(r);self.append(r,self.root/'FAILURE_LEDGER.jsonl')
 def snapshot(self,pid=None):
  rows=[json.loads(l) for l in self.path.read_text().splitlines()] if self.path.exists() else []
  if pid is not None:rows=[r for r in rows if r.get('pid')==pid]
  else:rows=[r for r in rows if r.get('job')==self.job]
  return dict(rhs=sum(r.get('rhs',0) for r in rows if r['kind']=='reservation'),full_fe=sum(r.get('factorizations',0) for r in rows if r['kind']=='reservation'))
 def receipt(self,status='COMPLETED',**meta):
  d=dict(job=self.job,pid=os.getpid(),status=status,wrapper_process_cpu_seconds=time.process_time()-self.cpu0,wall_seconds=time.perf_counter()-self.wall0,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),**self.snapshot(pid=os.getpid()),**meta)
  dump(self.root/'receipts'/f'{self.job}_{os.getpid()}.json',d);return d
