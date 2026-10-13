"""Only Joint P2 is writable; historical ledgers are never adopted/extended."""
from __future__ import annotations
import json,time,fcntl,os
from pathlib import Path
from contextlib import contextmanager

STAGE_CAPS=dict(mechanism=400,reference=2000,fresh=800,reserve=800)
def dump(path,value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    def numpy_default(x):
        if hasattr(x,'tolist'): return x.tolist()
        if hasattr(x,'item'): return x.item()
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False,default=numpy_default)+'\n')

class Limit(RuntimeError): pass

class Ledger:
    def __init__(self,root):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=True)
        self.path=self.root/'PHYSICS_LEDGER.jsonl';self.contract=self.root/'BUDGET_CONTRACT.json'
        self.stage='mechanism';self.job_start=None;self.job_id=None
        with self.locked():
            n=self.count()
            if self.contract.exists():
                d=json.loads(self.contract.read_text())
                if d['cap']!=4000 or n<d['high_water']:raise Limit('P2 ledger mismatch/truncation')
            else:dump(self.contract,dict(stage='B1-C1-JOINT-P2-20261013',cap=4000,stage_caps=STAGE_CAPS,high_water=n,wall_cap_seconds=7200,authorization='Human approved Joint P2 plan, 2026-10-13',historical_B1=384,historical_C1_P1='read only; no transfer'))
    @contextmanager
    def locked(self):
        with (self.root/'physics.lock').open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX)
            try:yield
            finally:fcntl.flock(f,fcntl.LOCK_UN)
    def rows(self):
        return [json.loads(l) for l in self.path.read_text().splitlines()] if self.path.exists() else []
    def count(self):return sum(r['event']=='BEGIN' for r in self.rows())
    def append(self,row):
        with self.path.open('a') as f:f.write(json.dumps(row,allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
    def charged_wall(self):
        p=self.root/'JOB_LEDGER.jsonl'
        return sum(json.loads(l)['wall_seconds'] for l in p.read_text().splitlines()) if p.exists() else 0.
    def check_time(self):
        running=time.perf_counter()-self.job_start if self.job_start else 0.
        if self.charged_wall()+running>=7200:raise Limit('P2 controlled host wall exhausted')
    def evaluate(self,kind,parameters,fn):
        self.check_time()
        with self.locked():
            rows=self.rows();n=self.count();ns=sum(r['event']=='BEGIN' and r['stage']==self.stage for r in rows)
            if n>=4000 or ns>=STAGE_CAPS[self.stage]:raise Limit(f'P2 budget exhausted: {self.stage} {ns}, total {n}')
            i=n+1;self.append(dict(event='BEGIN',id=i,stage=self.stage,kind=kind,parameters=parameters,job_id=self.job_id,pid=os.getpid()))
            d=json.loads(self.contract.read_text());d['high_water']=i;dump(self.contract,d)
        t=time.perf_counter()
        try:out=fn()
        except BaseException as e:
            with self.locked():self.append(dict(event='END',id=i,status='FAIL',error=type(e).__name__+': '+str(e),wall_seconds=time.perf_counter()-t))
            raise
        info=out[3] if isinstance(out,tuple) and len(out)>3 else {}
        with self.locked():self.append(dict(event='END',id=i,status='PASS',wall_seconds=time.perf_counter()-t,boundary_RHS=info.get('dense_boundary_solves',0)))
        return out
    def cache_hit(self,key):
        with self.locked():self.append(dict(event='CACHE_HIT',new_physics=0,key=key))
    @contextmanager
    def job(self,name):
        p=self.root/'JOB_LEDGER.jsonl'
        if p.exists() and any(json.loads(l)['id']==name for l in p.read_text().splitlines()):raise Limit('completed job id cannot be reused')
        self.check_time();start=time.perf_counter();cpu=time.process_time();n=self.count();self.job_start=start;self.job_id=name;status='PASS';error=None
        try:yield
        except BaseException as e:status='FAIL';error=type(e).__name__+': '+str(e);raise
        finally:
            self.job_start=None;self.job_id=None
            own=sum(a['event']=='BEGIN' and a.get('job_id')==name for a in self.rows())
            row=dict(id=name,status=status,error=error,wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,new_physics=own,global_physics_delta=self.count()-n,accounting='own tagged reservations; global delta may include concurrent independent jobs')
            with p.open('a') as f:f.write(json.dumps(row,ensure_ascii=False)+'\n')
    def summary(self):
        rows=self.rows();beg=[r for r in rows if r['event']=='BEGIN'];end=[r for r in rows if r['event']=='END']
        return dict(stage='B1-C1-JOINT-P2-20261013',cap=4000,used=len(beg),remaining=4000-len(beg),by_stage={s:sum(r['stage']==s for r in beg) for s in STAGE_CAPS},failures=sum(r['status']=='FAIL' for r in end),unfinished=len(beg)-len(end),kernel_wall_seconds=sum(r['wall_seconds'] for r in end),boundary_RHS=sum(r.get('boundary_RHS',0) for r in end),controlled_job_wall_seconds=self.charged_wall(),GPU_training=0,locked_truth_reads=0,historical_B1_unchanged=384)
