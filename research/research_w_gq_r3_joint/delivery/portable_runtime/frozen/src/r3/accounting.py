"""Fresh R3 pre-action ledger and shared one-heavy-FEM lock."""
from pathlib import Path
from contextlib import contextmanager
import os,json,time,fcntl,resource,sys,datetime
ROOT=Path(__file__).resolve().parents[2]
class BudgetExceeded(RuntimeError):pass
def plain(x):
    if hasattr(x,'tolist'):return x.tolist()
    if hasattr(x,'item'):return x.item()
    if isinstance(x,Path):return str(x)
    if isinstance(x,complex):return [x.real,x.imag]
    raise TypeError(type(x).__name__)
def utc_epoch(s):return datetime.datetime.fromisoformat(s.replace('Z','+00:00')).timestamp()
@contextmanager
def fem_lock(root=ROOT):
    """Required for model construction and every physical forward/Jacobian action."""
    with (Path(root)/'.fem.lock').open('a') as f:
        fcntl.flock(f,fcntl.LOCK_EX)
        try:yield
        finally:fcntl.flock(f,fcntl.LOCK_UN)
class BudgetBook:
    def __init__(self,root=ROOT):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=True)
        self.path=self.root/'PHYSICS_LEDGER.jsonl';self.cpu0=time.process_time()
        self.clock=json.loads((self.root/'RUN_CLOCK.json').read_text())
        self.deadline=utc_epoch(self.clock['deadline_utc'])
        self.heavy_after=utc_epoch(self.clock['no_new_heavy_after_utc'])
    @contextmanager
    def locked(self):
        with (self.root/'.ledger.lock').open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX)
            try:yield
            finally:fcntl.flock(f,fcntl.LOCK_UN)
    def check(self,new_job=False):
        if time.time()>self.deadline:raise BudgetExceeded('R3 180 minute deadline reached')
        if new_job and time.time()>self.heavy_after:raise BudgetExceeded('R3 new-heavy-job deadline reached')
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
        if rss>self.clock['single_job_rss_bytes']:raise BudgetExceeded('R3 RSS4GiB guard')
        return rss
    def append(self,row):
        with self.path.open('a') as f:f.write(json.dumps(dict(time_unix=time.time(),pid=os.getpid(),**row),default=plain,allow_nan=False)+'\n')
    def reserve(self,phase,label,rhs=0,factorizations=0,**meta):
        self.check()
        with self.locked():self.append(dict(kind='reservation',phase=phase,label=label,rhs=int(rhs),factorizations=int(factorizations),**meta))
    def event(self,phase='R3',label='',**meta):
        with self.locked():self.append(dict(kind='event',phase=phase,label=label or meta.pop('event','event'),**meta))
    def failure(self,phase,label,error,**meta):
        with self.locked():self.append(dict(kind='failure',phase=phase,label=label,error=str(error),**meta))
    def snapshot(self,pid=None):
        rows=[json.loads(l) for l in self.path.read_text().splitlines()] if self.path.exists() else []
        if pid is not None:rows=[r for r in rows if r.get('pid')==pid]
        return dict(rhs=sum(r.get('rhs',0) for r in rows if r['kind']=='reservation'),full_fe=sum(r.get('factorizations',0) for r in rows if r['kind']=='reservation'))
