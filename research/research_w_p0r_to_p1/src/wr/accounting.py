"""Append-only stage reservations. Historical books are never imported."""
import fcntl, json, os, time, uuid
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
def encode(v):
    if hasattr(v,'tolist'): return v.tolist()
    if isinstance(v,Path): return str(v)
    raise TypeError(type(v).__name__)
def write(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(obj,indent=2,default=encode,allow_nan=False)); temp.replace(path)
class BudgetStop(RuntimeError): pass
class Book:
    def __init__(self,stage):
        self.stage=stage; self.root=ROOT; self.path=ROOT/'PHYSICS_LEDGER.jsonl'
        self.config=json.loads((ROOT/'configs/AUTHORIZATION.json').read_text())
        if (ROOT/'CLOSED.json').exists(): raise BudgetStop('Campaign closed')
        self.start=time.process_time()
    def _rows(self,f):
        f.seek(0); return [json.loads(x) for x in f if x.strip()]
    def totals(self):
        if not self.path.exists(): rows=[]
        else: rows=[json.loads(x) for x in self.path.read_text().splitlines() if x.strip()]
        result={}
        for stage in self.config['stages']:
            rr=[r for r in rows if r.get('stage')==stage]
            result[stage]={k:sum(r.get(k,0) for r in rr if r.get('event')=='reserve') for k in ['rhs','factorizations']}
            result[stage]['cpu_seconds']=sum(r.get('charged_cpu_seconds',0) for r in rr)
        return result
    def event(self,**row):
        row=dict(utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),pid=os.getpid(),stage=self.stage,**row)
        with self.path.open('a+') as f:
            fcntl.flock(f,fcntl.LOCK_EX); f.write(json.dumps(row,default=encode,allow_nan=False)+'\n'); f.flush()
    def reserve(self,phase,label,rhs=0,factorizations=0,**meta):
        if min(rhs,factorizations)<0: raise ValueError('negative budget')
        with self.path.open('a+') as f:
            fcntl.flock(f,fcntl.LOCK_EX); rows=self._rows(f)
            local=[r for r in rows if r.get('stage')==self.stage]
            cap=self.config['stages'][self.stage]; reason=None
            for key,val in [('rhs',rhs),('factorizations',factorizations)]:
                used=sum(r.get(key,0) for r in local if r.get('event')=='reserve')
                if used+val>cap[key]: reason=key
            cpu=sum(r.get('charged_cpu_seconds',0) for r in local)+time.process_time()-self.start
            if cpu>=cap['cpu_seconds']: reason='cpu_seconds'
            row=dict(event='budget_reject' if reason else 'reserve',utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),pid=os.getpid(),stage=self.stage,phase=phase,label=label,rhs=int(rhs),factorizations=int(factorizations),id=uuid.uuid4().hex,**meta)
            if reason: row['reason']=reason
            f.seek(0,2); f.write(json.dumps(row,default=encode,allow_nan=False)+'\n'); f.flush()
            if reason: raise BudgetStop(reason)
        return row['id']
