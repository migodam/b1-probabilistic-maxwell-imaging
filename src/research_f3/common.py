from pathlib import Path
from collections import Counter
from contextlib import contextmanager
import json, os, time, fcntl, signal
import numpy as np
from scipy.special import ndtr, ndtri
from finite_feedback.common import plain, dump, pack, mass, rms
from lxz.local2d import MaterialPrior, TAU

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'research/research_f3_quotient_exploration'
OLD=ROOT.parent/'b1-research-f/research/research_f_intrinsic_transport'

class BudgetStop(RuntimeError): pass

class Book:
    """Independent F3 authorization; historical Research F is read-only; reserve before work."""
    def __init__(self,stage):
        self.stage=stage;self.path=OUT/'RESOURCE_LEDGER.jsonl';self.counts=Counter();self.context={}
        path_audit=stage=='F3_PATH_VALIDATION' and (OUT/'configs/PATH_VALIDATION.json').exists()
        if (OUT/'STAGE_CLOSED.json').exists():raise BudgetStop('F3 authorization fully closed')
        if (OUT/'RUN_CLOSED.json').exists() and not ((stage=='F3_RECOVERY_1' and (OUT/'RECOVERY_AUTHORIZATION.json').exists()) or path_audit):raise BudgetStop('Research F3 closed')
        if (OUT/'RUN_FINAL_CLOSED.json').exists() and not path_audit:raise BudgetStop('Research F3 pilot closed')
        rows=self.rows()
        if any(x['event']=='JOB_START' and x['stage']==stage for x in rows):
            raise BudgetStop('stage already attempted: explicit recovery, never reset')
        self.start=time.perf_counter();self.previous=sum(r.get('controlled_wall_seconds',0) for r in rows if r['event']=='JOB_END')
        self.append(dict(event='JOB_START',caps=dict(states=64,RHS=2048,wall_seconds=1200)))
        def stop(a,b):raise BudgetStop('1200 second cumulative physical wall')
        signal.signal(signal.SIGALRM,stop);signal.setitimer(signal.ITIMER_REAL,max(.001,1200-self.previous))
    def rows(self):return [json.loads(x) for x in self.path.read_text().splitlines()] if self.path.exists() else []
    def append(self,r):
        with self.path.open('a') as f:f.write(json.dumps(plain(dict(utc=time.time(),stage=self.stage,**r)),allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
    def totals(self,stage=None):
        rows=[r for r in self.rows() if stage is None or r['stage']==stage]
        return {k:sum(r.get('amount',0) for r in rows if r['event']==e) for k,e in [('states','STATE'),('RHS','RHS')]}
    def reserve(self,event,amount,**meta):
        with (OUT/'budget.lock').open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX)
            t=self.totals();key='states' if event=='STATE' else 'RHS'
            if t[key]+amount>{'states':64,'RHS':2048}[key]:raise BudgetStop('total hard cap')
            if self.context.get('phase')=='F1.5' and event=='STATE' and sum(r.get('amount',0) for r in self.rows() if r['event']=='STATE' and r.get('context',{}).get('phase')=='F1.5')>=16:raise BudgetStop('F1.5 sixteen state allocation')
            if time.perf_counter()-self.start+self.previous>=1200:raise BudgetStop('wall cap')
            self.append(dict(event=event,amount=int(amount),cumulative=t[key]+amount,context=self.context,**meta))
    def reserve_state(self,material,*,resolution,purpose):self.reserve('STATE',1,purpose=purpose,resolution=resolution)
    def reserve_rhs(self,count,role):self.reserve('RHS',count,role=role)
    @contextmanager
    def span(self,label,**counts):
        t=time.perf_counter();error=None;self.counts.update(counts)
        try:yield
        except BaseException as e:error=repr(e);raise
        finally:self.append(dict(event='ACTION',label=label,counts=counts,wall_seconds=time.perf_counter()-t,error=error,context=self.context))
    def end(self,error=None):
        elapsed=time.perf_counter()-self.start;signal.setitimer(signal.ITIMER_REAL,0)
        self.append(dict(event='JOB_END',controlled_wall_seconds=elapsed,error=error,counts=dict(self.counts),totals=self.totals()))

def whiten(z,sigma):return pack(z)*np.sqrt(2)/sigma
