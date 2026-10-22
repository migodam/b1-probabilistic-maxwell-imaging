"""Atomic physical reservations; job wrapper bills actual process CPU separately."""
from pathlib import Path
import fcntl, json, os, time, uuid

def encode(value):
    if hasattr(value,'tolist'):return value.tolist()
    raise TypeError(type(value).__name__)

class BudgetStop(RuntimeError): pass

class Book:
    CAPS = dict(rhs=3000, factorizations=1000, cpu_seconds=7200.)
    PHASE=None
    def __init__(self, root):
        self.root=Path(root); self.root.mkdir(parents=True, exist_ok=True)
        self.path=self.root/'PHYSICS_LEDGER.jsonl'
    def event(self, **row):
        with self.path.open('a+') as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            f.write(json.dumps(dict(utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),pid=os.getpid(),**row),default=encode)+'\n'); f.flush()
    def totals(self):
        rows=[json.loads(s) for s in self.path.read_text().splitlines()] if self.path.exists() else []
        return {**{k:sum(r.get(k,0) for r in rows if r.get('event')=='reserve') for k in ('rhs','factorizations')}, 'cpu_seconds':sum(r.get('charged_cpu_seconds',0.) for r in rows)}
    def reserve(self, phase, label, rhs=0, factorizations=0, **meta):
        if min(rhs,factorizations)<0: raise ValueError('negative reservation')
        with self.path.open('a+') as f:
            fcntl.flock(f,fcntl.LOCK_EX); f.seek(0)
            rows=[json.loads(s) for s in f if s.strip()]
            totals={k:sum(r.get(k,0) for r in rows if r.get('event')=='reserve') for k in ('rhs','factorizations')}
            cpu=sum(r.get('charged_cpu_seconds',0.) for r in rows)
            running=time.process_time()-float(os.environ.get('WQ_CPU_START','0'))
            reason=None
            for k,v in dict(rhs=int(rhs),factorizations=int(factorizations)).items():
                if totals[k]+v>self.CAPS[k]:reason=k
                if self.PHASE is not None:
                    local=sum(r.get(k,0) for r in rows if r.get('event')=='reserve' and r.get('phase')==phase)
                    if phase not in self.PHASE:raise ValueError(phase)
                    if local+v>self.PHASE[phase][k]:reason='phase_'+k
            if cpu+running>=self.CAPS['cpu_seconds']: reason='process_cpu_seconds'
            row=dict(event='budget_reject' if reason else 'reserve',utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),pid=os.getpid(),phase=phase,label=label,rhs=int(rhs),factorizations=int(factorizations),id=uuid.uuid4().hex,**meta)
            f.seek(0,2); f.write(json.dumps(row,default=encode)+'\n'); f.flush()
            if reason:raise BudgetStop(reason)
        return row['id']

class SearchBook(Book):
    CAPS=dict(rhs=8000,factorizations=2000,cpu_seconds=7200.)
    PHASE={name:dict(rhs=rhs,factorizations=lu) for name,rhs,lu in [('search',3200,800),('endpoint',2400,600),('validation',1600,400),('repair',800,200)]}
