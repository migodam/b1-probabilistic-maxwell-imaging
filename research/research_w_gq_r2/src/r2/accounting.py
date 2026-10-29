"""Independent phase and campaign hard budgets; reserve before each physics action."""
from pathlib import Path
import json,time,fcntl,contextlib,resource,sys,os
ROOT=Path(__file__).resolve().parents[2]
def plain(x):
    if hasattr(x,'tolist'):return x.tolist()
    if hasattr(x,'item'):return x.item()
    if isinstance(x,complex):return [x.real,x.imag]
    if isinstance(x,Path):return str(x)
    raise TypeError(type(x).__name__)
class BudgetExceeded(RuntimeError):pass
class BudgetBook:
    def __init__(self,root=ROOT):
        self.root=Path(root);self.path=self.root/'PHYSICS_LEDGER.jsonl'
        self.root.joinpath('results').mkdir(exist_ok=True)
        self.cpu0=time.process_time();self.auth=json.loads((self.root/'AUTHORIZATION.json').read_text());self.cap=self.auth['total']
    @contextlib.contextmanager
    def locked(self):
        with (self.root/'results/ledger.lock').open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX)
            try:yield
            finally:fcntl.flock(f,fcntl.LOCK_UN)
    def records(self):return [json.loads(l) for l in self.path.read_text().splitlines() if l] if self.path.exists() else []
    def snapshot(self,phase=None):
        rows=[r for r in self.records() if phase is None or r.get('phase')==phase]
        return dict(rhs=sum(r.get('rhs',0) for r in rows if r['kind']=='reservation'),
            full_fe=sum(r.get('factorizations',0) for r in rows if r['kind']=='reservation'),
            process_cpu_seconds=sum(r.get('process_cpu_seconds',0.) for r in rows if r.get('label')=='job_completion'))
    def append(self,row):
        with self.path.open('a') as f:f.write(json.dumps(dict(time_unix=time.time(),**row),default=plain,allow_nan=False)+'\n')
    def reserve(self,phase,label,rhs=0,factorizations=0,**meta):
        if (rhs or factorizations) and (self.root/'PHYSICS_CLOSED.json').exists():raise BudgetExceeded('physical campaign closed; retained outputs only')
        if phase not in self.auth['caps']:raise ValueError('unknown resource phase '+str(phase))
        if any(isinstance(x,bool) or int(x)!=x or x<0 for x in [rhs,factorizations]):raise ValueError('nonnegative integer charges')
        with self.locked():
            for ph,cap in [(None,self.cap),(phase,self.auth['caps'][phase])]:
                u=self.snapshot(ph)
                if (u['rhs']+rhs>cap['rhs'] or u['full_fe']+factorizations>cap['full_fe'] or
                    u['process_cpu_seconds']+time.process_time()-self.cpu0>cap['process_cpu_seconds']):raise BudgetExceeded(str((ph,u)))
            rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
            if rss>self.auth['rss_bytes']:raise BudgetExceeded('RSS exceeded')
            self.append(dict(kind='reservation',phase=phase,label=label,rhs=int(rhs),factorizations=int(factorizations),**meta))
    def event(self,phase=None,label='',**meta):
        be=meta.pop('event',None);label=label or be or 'event'
        for key in ['rhs','factorizations']:
            if key in meta:meta['observed_'+key]=meta.pop(key)
        with self.locked():self.append(dict(kind='event',phase=phase,label=label,**meta))
    def failure(self,phase,label,error,**meta):
        with (self.root/'FAILURE_LEDGER.jsonl').open('a') as f:f.write(json.dumps(dict(time_unix=time.time(),phase=phase,label=label,error=str(error),**meta),default=plain)+'\n')
