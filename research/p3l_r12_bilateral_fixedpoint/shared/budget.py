"""Cross-process prepayment guard. No physical operation without reservation.

All reservations remain charged even on exceptions/crashes. Completion records
report success/timing, never refund an attempted solve. flock protects both the
snapshot and append-only ledger; no server, daemon or database is required.
"""
from pathlib import Path
from contextlib import contextmanager
import fcntl, json, os, time, uuid

ROOT=Path(__file__).resolve().parents[1]
COUNTERS=('LU','RHS','RHS576')

class BudgetStop(RuntimeError): pass

def atomic_json(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.'+uuid.uuid4().hex+'.tmp')
    with temp.open('w') as f:
        json.dump(obj,f,ensure_ascii=False,indent=2,allow_nan=False,default=lambda x:x.tolist() if hasattr(x,'tolist') else str(x));f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(temp,path)

class Budget:
    def __init__(self,root=ROOT):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=True)
        self.lock=self.root/'GLOBAL_BUDGET.lock';self.snapshot=self.root/'GLOBAL_BUDGET_STATE.json'
        self.ledger=self.root/'COST_LEDGER.jsonl'
    @contextmanager
    def locked(self):
        with self.lock.open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX)
            yield
            fcntl.flock(f,fcntl.LOCK_UN)
    def initialize(self,deadline,limits=None):
        with self.locked():
            if self.snapshot.exists(): raise RuntimeError('Budget already initialized; no implicit reset')
            s=dict(limits=limits or dict(LU=12000,RHS=300000,RHS576=1024),used=dict.fromkeys(COUNTERS,0),sequence=0,deadline=deadline,initialized=time.time(),stop=False)
            atomic_json(self.snapshot,s)
    def _event(self,s,event,**kw):
        s['sequence']+=1
        row=dict(sequence=s['sequence'],utc=time.time(),pid=os.getpid(),event=event,**kw,used=s['used'].copy())
        # Durable conservative charge FIRST. A crash before ledger append can
        # overcharge an unattempted operation, but can never permit overspending.
        # The operation starts only after BOTH this snapshot and ledger fsync.
        atomic_json(self.snapshot,s)
        with self.ledger.open('a') as f:
            f.write(json.dumps(row,ensure_ascii=False,allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
    def reserve(self,group,LU=0,RHS=0,RHS576=0,**meta):
        request=dict(LU=int(LU),RHS=int(RHS),RHS576=int(RHS576))
        if any(v<0 for v in request.values()):raise ValueError('negative reservation')
        if RHS576>RHS:raise ValueError('RHS576 must also be included in all-order RHS')
        with self.locked():
            s=json.loads(self.snapshot.read_text())
            if s['stop'] or time.time()>=s['deadline']:raise BudgetStop('global STOP or physics deadline')
            if any(s['used'][k]+request[k]>s['limits'][k] for k in COUNTERS):
                self._event(s,'DENIED',group=group,request=request,meta=meta)
                raise BudgetStop('atomic global cap; operation not executed')
            s['used']={k:s['used'][k]+request[k] for k in COUNTERS}
            ticket=uuid.uuid4().hex
            self._event(s,'RESERVE',ticket=ticket,group=group,request=request,meta=meta)
            return ticket
    def complete(self,ticket,group,wall,status='OK',**meta):
        with self.locked():
            s=json.loads(self.snapshot.read_text())
            self._event(s,'COMPLETE',ticket=ticket,group=group,wall_seconds=float(wall),status=status,meta=meta)
    @contextmanager
    def operation(self,group,LU=0,RHS=0,RHS576=0,**meta):
        ticket=self.reserve(group,LU,RHS,RHS576,**meta);ts=time.perf_counter()
        try:yield ticket
        except BaseException as e:
            self.complete(ticket,group,time.perf_counter()-ts,'FAILED',error=repr(e));raise
        else:self.complete(ticket,group,time.perf_counter()-ts)
    def stop(self,reason):
        with self.locked():
            s=json.loads(self.snapshot.read_text());s['stop']=True;self._event(s,'STOP',reason=reason)
    def state(self):
        with self.locked():return json.loads(self.snapshot.read_text())
