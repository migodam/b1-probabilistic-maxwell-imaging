"""Research A only: append-only attempts and cumulative controlled job wall."""
import json, time, os, fcntl
from pathlib import Path
from contextlib import contextmanager
from b1_joint.accounting import dump

class Limit(RuntimeError): pass

class Budget:
    def __init__(self, root, stage='E1'):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=True)
        self.path=self.root/'COST_LEDGER.jsonl';self.stage=stage;self.context={};self.start=None
        self.cap=2304 if stage=='E1' else 384;self.wallcap=1800 if stage=='E1' else 3600
        self.job_id=None
    def rows(self):
        return [json.loads(s) for s in self.path.read_text().splitlines()] if self.path.exists() else []
    def append(self,row):
        with self.path.open('a') as f:
            f.write(json.dumps(dict(stage=self.stage,utc=time.time(),**row),default=lambda x:x.tolist(),allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
    def count(self):return sum(r['event']=='BEGIN' and r['stage']==self.stage for r in self.rows())
    def used_wall(self):return sum(r['wall_seconds'] for r in self.rows() if r['event']=='JOB_END' and r['stage']==self.stage)
    def check(self):
        if self.used_wall()+(time.perf_counter()-self.start if self.start else 0)>=self.wallcap:raise Limit('controlled cumulative wall exhausted')
    def evaluate(self,kind,parameters,fn):
        if (self.root/'RUN_CLOSED.json').exists():raise Limit('Research A scientific run closed; no new forward without fresh explicit authorization')
        self.check()
        with (self.root/'budget.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            n=self.count()
            if n>=self.cap:raise Limit('Research A physical attempt limit')
            i=n+1;self.append(dict(event='BEGIN',attempt=i,kind=kind,parameters=parameters,job=self.job_id,context=self.context))
            fcntl.flock(lock,fcntl.LOCK_UN)
        t=time.perf_counter()
        try:
            result=fn();meta=result[3] if isinstance(result,tuple) and len(result)>3 else {}
        except BaseException as e:
            self.append(dict(event='END',attempt=i,status='FAIL',error=repr(e),wall_seconds=time.perf_counter()-t));raise
        self.append(dict(event='END',attempt=i,status='PASS',wall_seconds=time.perf_counter()-t,boundary_RHS=meta.get('dense_boundary_solves',0)))
        return result
    def cache_hit(self,key):self.append(dict(event='CACHE_HIT',new_physics=0,key=key,context=self.context))
    @contextmanager
    def job(self,name):
        if any(r['event']=='JOB_END' and r.get('job')==name for r in self.rows()):raise Limit('completed job cannot be rerun under same id')
        self.check();self.start=time.perf_counter();self.job_id=name;cpu=time.process_time();status='PASS';error=None
        self.append(dict(event='JOB_BEGIN',job=name))
        try:yield
        except BaseException as e:status='FAIL';error=repr(e);raise
        finally:
            wall=time.perf_counter()-self.start;self.start=None
            self.append(dict(event='JOB_END',job=name,wall_seconds=wall,cpu_seconds=time.process_time()-cpu,status=status,error=error));self.job_id=None
    def summary(self):
        rows=[r for r in self.rows() if r['stage']==self.stage]
        return dict(stage=self.stage,attempts=self.count(),cap=self.cap,remaining=self.cap-self.count(),controlled_execution_wall=self.used_wall(),wall_cap=self.wallcap,failed=sum(r.get('status')=='FAIL' and r['event']=='END' for r in rows),boundary_RHS=sum(r.get('boundary_RHS',0) for r in rows),kernel_wall=sum(r['wall_seconds'] for r in rows if r['event']=='END'),GPU=0,old_P2_budget_adopted=False)
