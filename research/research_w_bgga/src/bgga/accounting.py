"""BGGA scope of the W-SLAM reserve-before-work protocol, with explicit grant.

Append records rather than rewriting a large historical ledger per RHS block.
Only wrapper completion records carry aggregate CPU. Solver events are kernels.
"""
import contextlib, fcntl, json, math, os, time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]

class BudgetExceeded(RuntimeError):pass
class ResourceAuthorizationRequired(BudgetExceeded):pass

def _default(v):
    if hasattr(v,'tolist'):return v.tolist()
    if hasattr(v,'item'):return v.item()
    if isinstance(v,complex):return [v.real,v.imag]
    raise TypeError(type(v).__name__)

class BudgetBook:
    """Same FEM .reserve/.event API, never points at an old W/WSLAM ledger."""
    def __init__(self,root=ROOT):
        self.root=Path(root).resolve();self.path=self.root/'PHYSICS_LEDGER.jsonl'
        self.lock=self.root/'results'/'ledger.lock';self.lock.parent.mkdir(parents=True,exist_ok=True)
        self.declared_completion=False
    def records(self):
        return [json.loads(l) for l in self.path.read_text().splitlines() if l] if self.path.exists() else []
    @contextlib.contextmanager
    def locked(self):
        with self.lock.open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX)
            try:yield
            finally:fcntl.flock(f,fcntl.LOCK_UN)
    def _append(self,rec):
        rec=dict(time_unix=time.time(),**rec)
        with self.path.open('a') as f:f.write(json.dumps(rec,default=_default,allow_nan=False)+'\n')
        return rec
    def snapshot(self,phase):
        keys=('rhs','factorizations','process_cpu_seconds')
        return {k:sum(r.get(k,0) for r in self.records() if r.get('phase')==phase) for k in keys}
    def authorization(self,phase):
        a=json.loads((self.root/'AUTHORIZATION.json').read_text())
        cap=a.get('authorized_caps',{}).get(phase)
        if not a.get('new_physics_authorized') or cap is None:
            raise ResourceAuthorizationRequired(f'{phase}: BGGA new-physics grant absent; old scope is not transferable')
        return cap
    def reserve(self,phase,label,rhs=0,factorizations=0,**meta):
        if any(not isinstance(x,int) or isinstance(x,bool) or x<0 for x in (rhs,factorizations)):
            raise ValueError('Nonnegative integer charges required')
        with self.locked():
            caps=self.authorization(phase);u=self.snapshot(phase)
            auth=json.loads((self.root/'AUTHORIZATION.json').read_text())
            if phase in auth.get('closed_physics_stages',[]) and (rhs or factorizations):
                raise BudgetExceeded(f'{phase}: physical campaign closed; residual quota is not a new authorization')
            if (u['rhs']+rhs>caps['rhs'] or u['factorizations']+factorizations>caps['full_fe']
                    or u['process_cpu_seconds']>=caps['cpu_seconds']):
                raise BudgetExceeded(f'{phase}: budget refuses request rhs={rhs}, FE={factorizations}, used={u}')
            ratio=max([u[k]/limit for k,limit in [('rhs',caps['rhs']),('factorizations',caps['full_fe']),('process_cpu_seconds',caps['cpu_seconds'])] if limit>0]+[0.])
            if ratio>=.7 and not any(r.get('phase')==phase and r.get('label')=='budget_70_percent_checkpoint' for r in self.records()):
                checkpoint=dict(phase=phase,label='budget_70_percent_checkpoint',used=u,caps=caps,
                    policy='90 percent: predeclared necessary verification/closeout only; no new campaign')
                self._append(dict(kind='checkpoint',**checkpoint))
                (self.root/'results'/f'{phase}_70_PERCENT_CHECKPOINT.json').write_text(json.dumps(checkpoint,indent=2)+'\n')
            closing=bool(meta.pop('declared_completion',False) or self.declared_completion)
            if ratio>=.9 and not closing:
                raise BudgetExceeded('90% reached: only explicitly predeclared completion/close-out is allowed')
            if closing:meta['predeclared_closeout']=True
            return self._append(dict(phase=phase,kind='reservation',label=label,rhs=rhs,factorizations=factorizations,**meta))
    def event(self,phase=None,label='',**meta):
        phase=phase or meta.pop('stage',os.environ.get('BGGA_STAGE','PREFLIGHT'))
        for k in ('rhs','factorizations'): 
            if k in meta:meta['observed_'+k]=meta.pop(k)
        cpu=meta.get('process_cpu_seconds',0.)
        if not isinstance(cpu,(int,float)) or not math.isfinite(cpu) or cpu<0:raise ValueError('CPU must be finite and nonnegative')
        with self.locked():return self._append(dict(phase=phase,kind='event',label=label,**meta))
    def dense(self,phase,label,n,**meta):
        return self.event(phase,label,kind_override='dense_core_factorization',dense_dimension=int(n),dense_factorizations=1,**meta)
