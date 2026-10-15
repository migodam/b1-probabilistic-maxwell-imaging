"""Research B independent, append-only accounting; never adopts P2 balances."""
from __future__ import annotations
import fcntl, json, os, time
from pathlib import Path
from contextlib import contextmanager

ROOT = Path(__file__).resolve().parents[1]

def dump(path, obj):
    def convert(x):
        if hasattr(x, 'tolist'): return x.tolist()
        if hasattr(x, 'item'): return x.item()
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    p=Path(path); p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=convert,allow_nan=False)+'\n')

class BudgetExceeded(RuntimeError): pass

class Cost:
    def __init__(self, root=ROOT):
        self.root=Path(root); self.results=self.root/'results'; self.results.mkdir(parents=True,exist_ok=True)
        self.physics=self.results/'PHYSICS_LEDGER.jsonl'; self.jobs=self.results/'COMPUTE_LEDGER.jsonl'
        self.stage='E0'; self.job_id=None; self.start=None
        self.overlay=json.loads((self.root/'configs/AUTHORIZATION_OVERLAY.json').read_text())
        authorized=self.overlay.get('E1_E2',{})
        self.execution_authorized=authorized.get('execution_authorized') is True
        attempt_cap=authorized.get('material_frequency_attempt_cap',4096)
        wall_cap=authorized.get('controlled_execution_wall_seconds_cap',7200)
        if (isinstance(attempt_cap,bool) or not isinstance(attempt_cap,int) or attempt_cap<0
                or isinstance(wall_cap,bool) or not isinstance(wall_cap,(int,float))
                or not 0<=wall_cap<float('inf')):
            raise ValueError('Nonnegative finite authorization caps required')
        # A local overlay may tighten the user authorization, never enlarge it.
        self.attempt_cap=min(4096,attempt_cap); self.wall_cap=min(7200,float(wall_cap))
    def rows(self,path):
        return [json.loads(s) for s in path.read_text().splitlines()] if path.exists() else []
    @contextmanager
    def locked(self):
        with (self.results/'accounting.lock').open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX)
            try: yield
            finally: fcntl.flock(f,fcntl.LOCK_UN)
    def append(self,path,row):
        with path.open('a') as f:
            f.write(json.dumps(row,ensure_ascii=False,allow_nan=False)+'\n'); f.flush(); os.fsync(f.fileno())
    def count(self): return sum(r['event']=='BEGIN' for r in self.rows(self.physics))
    def charged_wall(self):
        return sum(r['wall_seconds'] for r in self.rows(self.jobs) if r['stage'] in ('E1','E2'))
    def check_time(self):
        if self.stage in ('E1','E2'):
            if not self.execution_authorized:
                raise BudgetExceeded('Research B E1/E2 execution not authorized by overlay')
            elapsed=time.perf_counter()-self.start if self.start is not None else 0.
            if self.charged_wall()+elapsed>=self.wall_cap:
                raise BudgetExceeded(f'Research B E1/E2 controlled wall cap {self.wall_cap:g}s')
    def evaluate(self,kind,parameters,fn):
        self.check_time()
        with self.locked():
            n=self.count()
            if (self.stage not in ('E1','E2') or not self.execution_authorized
                    or n>=self.attempt_cap or self.job_id is None):
                raise BudgetExceeded('E0/attempt cap/unmetered call rejected before solver')
            self.append(self.physics,dict(event='BEGIN',id=n+1,stage=self.stage,kind=kind,parameters=parameters,job_id=self.job_id))
        start=time.perf_counter()
        try: out=fn()
        except BaseException as e:
            self.append(self.physics,dict(event='END',id=n+1,status='FAIL',error=repr(e),wall_seconds=time.perf_counter()-start)); raise
        info=out[3] if isinstance(out,tuple) and len(out)>3 else {}
        systems=info.get('dense_boundary_solves',0)
        if not systems and parameters.get('implementation')=='recursive':
            systems=2*int(parameters.get('nmax',0))*max(len(parameters.get('radii',[]))-1,0)
        self.append(self.physics,dict(event='END',id=n+1,status='PASS',wall_seconds=time.perf_counter()-start,
            boundary_RHS=systems,factorizations=systems))
        return out
    def cache_hit(self,key):
        self.append(self.physics,dict(event='CACHE_HIT',stage=self.stage,job_id=self.job_id,new_physics=0,key=key))
    @contextmanager
    def job(self,name,stage):
        if any(r['id']==name for r in self.rows(self.jobs)): raise ValueError('Use a new job identity, preserve old receipt')
        self.stage=stage; self.job_id=name; self.start=time.perf_counter(); cpu=time.process_time(); n=self.count(); status='PASS'; error=None
        try:
            try: os.nice(10)
            except OSError: pass
            self.check_time(); yield self
        except BaseException as e: status='FAIL'; error=repr(e); raise
        finally:
            self.append(self.jobs,dict(id=name,stage=stage,status=status,error=error,wall_seconds=time.perf_counter()-self.start,
                cpu_seconds=time.process_time()-cpu,new_maxwell_attempts=self.count()-n,pid=os.getpid(),GPU=0))
            self.start=None; self.job_id=None
            dump(self.results/'RESOURCE_SUMMARY.json',self.summary())
    def summary(self):
        rows=self.rows(self.physics);beg=[r for r in rows if r['event']=='BEGIN']; end=[r for r in rows if r['event']=='END'];jobs=self.rows(self.jobs)
        return dict(authorization_id=self.overlay['authorization_id'],used_attempts=len(beg),cap=self.attempt_cap,remaining_attempts=self.attempt_cap-len(beg),
            by_stage={s:sum(r['stage']==s for r in beg) for s in ('E0','E1','E2')},failed_attempts=sum(r['status']=='FAIL' for r in end),
            unfinished_attempts=len(beg)-len(end),kernel_wall_seconds=sum(r['wall_seconds'] for r in end),
            E1_E2_controlled_wall_seconds=self.charged_wall(),E1_E2_wall_cap=self.wall_cap,
            all_controlled_wall_seconds=sum(r['wall_seconds'] for r in jobs),cpu_seconds=sum(r['cpu_seconds'] for r in jobs),
            boundary_RHS=sum(r.get('boundary_RHS',0) for r in end)+sum(r.get('additional_boundary_RHS',0) for r in rows if r['event']=='ACCOUNTING_CORRECTION'),
            factorizations=sum(r.get('factorizations',0) for r in end)+sum(r.get('additional_factorizations',0) for r in rows if r['event']=='ACCOUNTING_CORRECTION'),
            jobs=len(jobs),failed_jobs=sum(r['status']=='FAIL' for r in jobs),GPU_training=0,historical_ledgers='read only, no balance transfer')
