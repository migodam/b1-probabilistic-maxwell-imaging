"""Independent human-authorized192 states/6144 RHS/7200s; append only."""
from pathlib import Path
from contextlib import contextmanager
from collections import Counter
import json,time,os,fcntl
import numpy as np

def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,dict):return {k:plain(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [plain(v) for v in x]
    return x
def dump(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(plain(obj),ensure_ascii=False,indent=2,allow_nan=False)+'\n')
class Limit(RuntimeError):pass

class Book:
    def __init__(self,root):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=True)
        self.path=self.root/'COST_LEDGER.jsonl';self.counts=Counter();self.context={};self.started=None
        self.caps=dict(states=192,rhs=6144,wall=7200)
        self.name=None
    def rows(self):return [json.loads(s) for s in self.path.read_text().splitlines()] if self.path.exists() else []
    def append(self,row):
        with self.path.open('a') as f:
            f.write(json.dumps(plain(dict(utc=time.time(),**row)),allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
    def totals(self):
        rows=self.rows()
        return dict(states=sum(r.get('amount',0) for r in rows if r['event']=='STATE'),rhs=sum(r.get('amount',0) for r in rows if r['event']=='RHS'),wall=sum(r['wall_seconds'] for r in rows if r['event']=='JOB_END'))
    def check(self):
        t=self.totals()
        if t['wall']+(time.perf_counter()-self.started if self.started else 0)>=7200:raise Limit('ResearchT cumulative controlled wall exhausted')
    def reserve(self,event,n,limit,**metadata):
        self.check()
        if (self.root/'RUN_CLOSED.json').exists():raise Limit('ResearchT closed; new solves require new explicit authorization')
        with (self.root/'budget.lock').open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX)
            used=self.totals()['states' if event=='STATE' else 'rhs']
            if used+n>limit:raise Limit(f'ResearchT {event} exhausted ({used}+{n}>{limit})')
            self.append(dict(event=event,amount=n,cumulative=used+n,job=self.name,context=self.context,**metadata))
    def reserve_state(self,material,*,resolution,purpose):self.reserve('STATE',1,192,resolution=resolution,purpose=purpose)
    def reserve_rhs(self,count,role):self.reserve('RHS',int(count),6144,role=role)
    @contextmanager
    def span(self,label,**counts):
        self.check();self.counts.update(counts);t=time.perf_counter();cpu=time.process_time();status='PASS';error=None
        try:yield
        except BaseException as e:status='FAIL';error=repr(e);raise
        finally:self.append(dict(event='ACTION',label=label,counts=counts,wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,status=status,error=error,context=self.context))
    @contextmanager
    def job(self,name):
        if any(r['event']=='JOB_END' and r.get('job')==name for r in self.rows()):raise Limit('completed job cannot be silently repeated')
        self.check();self.name=name;self.started=time.perf_counter();cpu=time.process_time();status='PASS';error=None
        self.append(dict(event='JOB_BEGIN',job=name))
        try:yield
        except BaseException as e:status='FAIL';error=repr(e);raise
        finally:
            elapsed=time.perf_counter()-self.started;self.started=None
            self.append(dict(event='JOB_END',job=name,wall_seconds=elapsed,cpu_seconds=time.process_time()-cpu,status=status,error=error));self.name=None
    def receipt(self):
        r=self.totals();rows=self.rows()
        cumulative=Counter()
        for row in rows:
            if row['event']=='ACTION':cumulative.update(row.get('counts',{}))
        return dict(new_full_state_attempts=r['states'],state_cap=192,complex_full_solver_RHS=r['rhs'],RHS_cap=6144,controlled_execution_wall=r['wall'],wall_cap=7200,GPU=0,old_budget_credit=0,failed_actions=sum(a.get('status')=='FAIL' and a['event']=='ACTION' for a in rows),action_counts=dict(cumulative))
