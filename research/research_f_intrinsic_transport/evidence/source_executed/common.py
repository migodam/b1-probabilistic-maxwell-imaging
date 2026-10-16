from pathlib import Path
from collections import Counter
from contextlib import contextmanager
import json, os, time, fcntl, signal
import numpy as np
from scipy.special import ndtr, ndtri
from finite_feedback.common import plain, dump, pack, mass, rms
from lxz.local2d import MaterialPrior, TAU

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'research/research_f_intrinsic_transport'

class BudgetStop(RuntimeError): pass

class Book:
    """One authorization, cumulative across F1/F2; reserve before work."""
    def __init__(self,stage):
        self.stage=stage;self.path=OUT/'RESOURCE_LEDGER.jsonl';self.counts=Counter();self.context={}
        if (OUT/'RUN_CLOSED.json').exists():raise BudgetStop('Research F closed')
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
            if self.stage=='F1' and event=='STATE' and self.totals('F1')['states']>=24:raise BudgetStop('F1 24 state allocation')
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

class Prior:
    """Original eight independent uniform coordinates. No changed measure."""
    def __init__(self,metadata):
        self.meta=metadata;self.chart=MaterialPrior(metadata['generator'],metadata['background_id'])
        self.p0=ndtr(self.chart.coordinates())
    def decode(self,p):
        if np.any(np.asarray(p)<=0) or np.any(np.asarray(p)>=1):raise ValueError('open original prior interior')
        return self.chart.decode(ndtri(p))
    def derivative(self,p):
        u=ndtri(p);return self.chart.analytic_derivative(u)/((np.exp(-.5*u*u)/np.sqrt(2*np.pi))[None,:])
    def frame(self,p,C,sigma):
        D=self.derivative(p);DM=np.r_[D.real,D.imag]/24
        Q,sv,Vh=np.linalg.svd(DM,full_matrices=False);rank=int(sum(sv>sv[0]*1e-7))
        E=Vh[:rank].T/sv[:rank];V=C.reshape(256,576)@D
        A=np.column_stack([whiten(V[:,i].reshape(32,8),sigma) for i in range(8)])
        U,s,VT=np.linalg.svd(A@E,full_matrices=False)
        v=E@VT[-1];w=D@v
        if w[np.argmax(abs(w))].real<0:v=-v;w=-w
        return dict(D=D,DM=DM,Q=Q[:,:rank],E=E,A=A,sv=s,material_sv=sv,rank=rank,v=v,w=w)

def save_endpoint(state,p,path,**meta):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('wb') as f:
        np.savez_compressed(f,p=p,chi=state.chi,current=state.current,field=state.field,prediction=state.data,
                            LU=state.factor[0],pivots=state.factor[1],backward_residual=state.residual)
        f.flush();os.fsync(f.fileno())
    dump(path.with_suffix('.json'),meta)

def response_jac(T,state):return np.einsum('qp,ps->qsp',T,state.field).reshape(256,576)
