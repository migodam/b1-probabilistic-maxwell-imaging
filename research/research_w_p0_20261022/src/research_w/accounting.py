"""One append-only physical budget, reserved before any native action."""
from pathlib import Path
from contextlib import contextmanager
import fcntl,json,time,os,resource,uuid

class BudgetStop(RuntimeError): pass

class Book:
    CAPS={"rhs":8000,"factorizations":2000,"cpu_seconds":7200.0}
    PHASE={"search":(.4,3200,800),"endpoint":(.3,2400,600),
           "validation":(.2,1600,400),"repair":(.1,800,200)}
    def __init__(self,root):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=True)
        self.path=self.root/'PHYSICS_LEDGER.jsonl'
    def _read(self,f):
        f.seek(0);return [json.loads(s) for s in f if s.strip()]
    def event(self,**entry):
        with self.path.open('a+') as f:
            fcntl.flock(f,fcntl.LOCK_EX);f.write(json.dumps(dict(time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),pid=os.getpid(),**entry))+'\n');f.flush()
    def reserve(self,phase,label,rhs=0,factorizations=0,**meta):
        if phase not in self.PHASE:raise ValueError(phase)
        with self.path.open('a+') as f:
            fcntl.flock(f,fcntl.LOCK_EX);rows=self._read(f)
            reservations=[r for r in rows if r.get('event')=='reserve']
            totals={k:sum(r.get(k,0) for r in reservations) for k in ['rhs','factorizations']}
            local={k:sum(r.get(k,0) for r in reservations if r['phase']==phase) for k in ['rhs','factorizations']}
            cpu=sum(r.get('charged_cpu_seconds',0) for r in rows)
            limits=dict(rhs=self.PHASE[phase][1],factorizations=self.PHASE[phase][2])
            add=dict(rhs=int(rhs),factorizations=int(factorizations))
            reason=None
            for k in add:
                if totals[k]+add[k]>self.CAPS[k] or local[k]+add[k]>limits[k]:reason=k
            if cpu>=self.CAPS['cpu_seconds']:reason='cpu_seconds'
            deadline=float(os.environ.get('W_JOB_DEADLINE_MONOTONIC','inf'))
            if time.monotonic()>=deadline:reason='in_progress_cpu_wall_deadline'
            entry=dict(event='reserve' if reason is None else 'budget_reject',phase=phase,label=label,**add,**meta,id=uuid.uuid4().hex)
            f.write(json.dumps(entry)+'\n');f.flush()
            if reason:raise BudgetStop(f'{phase}:{reason}')
            return entry['id']
    @contextmanager
    def span(self,label,**meta):
        wall=time.perf_counter();cpu=time.process_time();status='COMPLETE';error=None
        try:yield
        except BaseException as e:
            status='FAILED';error=repr(e);raise
        finally:
            w=time.perf_counter()-wall;c=time.process_time()-cpu
            self.event(event='cost',label=label,status=status,error=error,wall_seconds=w,process_cpu_seconds=c,charged_cpu_seconds=max(w,c),peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,**meta)
    def totals(self):
        if not self.path.exists():return dict(rhs=0,factorizations=0,cpu_seconds=0)
        with self.path.open() as f:rows=self._read(f)
        return {**{k:sum(r.get(k,0) for r in rows if r.get('event')=='reserve') for k in ['rhs','factorizations']},'cpu_seconds':sum(r.get('charged_cpu_seconds',0) for r in rows)}
