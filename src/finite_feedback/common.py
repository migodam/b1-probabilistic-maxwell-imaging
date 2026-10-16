from pathlib import Path
from collections import Counter
from contextlib import contextmanager
import json,os,time,fcntl
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'research/finite_feedback_fiber_diagnostic_20261016'
LXZ=ROOT.parent/'b1-research-lxz/research/research_lxz_20261016'
def plain(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,dict):return {str(k):plain(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [plain(x) for x in v]
    return v
def dump(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(plain(x),ensure_ascii=False,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def pack(z):return np.concatenate((z.real.T,z.imag.T),axis=1).ravel()
def mass(z):return np.r_[z.real,z.imag]/24
def rms(z):return float(np.linalg.norm(z)/24)

class Book:
    def __init__(self):
        self.path=OUT/'RESOURCE_LEDGER.jsonl';self.counts=Counter();self.context={};self.start=time.perf_counter()
        if (OUT/'RUN_CLOSED.json').exists():raise RuntimeError('diagnostic authorization closed')
        if self.path.exists() and any(json.loads(s).get('event')=='JOB_START' for s in self.path.read_text().splitlines()):raise RuntimeError('existing physical job: explicit recovery required, no reset')
        self.append(dict(event='JOB_START',caps=dict(states=48,RHS=2048,wall_seconds=1200)))
    def append(self,row):
        with self.path.open('a') as f:f.write(json.dumps(plain(dict(utc=time.time(),**row)),allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
    def totals(self):
        rows=[json.loads(s) for s in self.path.read_text().splitlines()]
        return {k:sum(r.get('amount',0) for r in rows if r.get('event')==event) for k,event in [('states','STATE'),('RHS','RHS')]}
    def reserve(self,event,n,**extra):
        with (OUT/'budget.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            if (OUT/'RUN_CLOSED.json').exists() or time.perf_counter()-self.start>1200:raise RuntimeError('diagnostic closed/wall cap')
            key='states' if event=='STATE' else 'RHS';t=self.totals()
            if t[key]+n>{'states':48,'RHS':2048}[key]:raise RuntimeError('diagnostic hard budget cap')
            self.append(dict(event=event,amount=n,cumulative=t[key]+n,context=self.context,**extra))
    def reserve_state(self,material,*,resolution,purpose):self.reserve('STATE',1,resolution=resolution,purpose=purpose)
    def reserve_rhs(self,count,role):self.reserve('RHS',int(count),role=role)
    @contextmanager
    def span(self,label,**counts):
        st=time.perf_counter();self.counts.update(counts);error=None
        try:yield
        except BaseException as e:error=repr(e);raise
        finally:self.append(dict(event='ACTION',label=label,counts=counts,wall_seconds=time.perf_counter()-st,error=error,context=self.context))
