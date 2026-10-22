from __future__ import annotations
from pathlib import Path
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass
import fcntl,json,os,sys,time,signal
import numpy as np
from scipy.special import ndtr,ndtri
from p1_material.prior import MaterialPrior

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'research/p3l_local_posterior_imaging'
CAPS=dict(A=5800000,B=2500000,C=1500000,D=1500000,E=700000)
TOTAL_RHS=12000000
CPU_CAP=48*3600
GPU_CAP=24*3600

def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,Path):return str(x)
    if isinstance(x,complex):return [x.real,x.imag]
    if isinstance(x,dict):return {str(k):plain(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [plain(v) for v in x]
    return x

def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(plain(obj),ensure_ascii=False,allow_nan=False,indent=2)+'\n');temp.replace(path)

def save(path,**arrays):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.tmp')
    with temp.open('wb') as f:np.savez_compressed(f,**arrays);f.flush()
    temp.replace(path)

def pack(z):return np.concatenate((z.real.T,z.imag.T),axis=1).ravel()

class BudgetStop(RuntimeError):pass

class Book:
    """Atomic campaign reservations. JSONL is authoritative, snapshot recoverable.

    A forward reserves its eight RHS and LU before the unchanged A30 adapter
    starts. The adapter's subsequent full_forward callback consumes that credit,
    so rejected, interrupted and failed attempts cannot disappear or double bill.
    """
    def __init__(self,stage,job,max_wall=7200):
        if stage not in CAPS:raise ValueError(stage)
        if (OUT/'CAMPAIGN_CLOSED.json').exists():raise BudgetStop('P3-L closed')
        auth=json.loads((OUT/'configs/AUTHORIZATION.json').read_text())
        if auth['campaign_id']!='P3L-20261022' or not auth['explicit_user_authorization']:raise BudgetStop('authorization')
        self.stage,self.job=stage,job;self.counts=Counter();self.credit=0;self.context={};self.finished=False
        self.path=OUT/'RESOURCE_LEDGER.jsonl';self.snapshot=OUT/'RESOURCE_LEDGER.json'
        self.jobdir=OUT/stage/'jobs'/job;self.jobdir.mkdir(parents=True,exist_ok=True)
        self.joblock=(self.jobdir/'job.lock').open('a');fcntl.flock(self.joblock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if (self.jobdir/'CLOSED.json').exists():raise BudgetStop('job already closed; new experiment needs registered new job')
        self.control=(OUT/'budget.lock').open('a');self.max_wall=float(max_wall);self.start=time.perf_counter()
        self._update('JOB_START',reserve_wall=self.max_wall)
        signal.signal(signal.SIGALRM,lambda a,b:(_ for _ in ()).throw(BudgetStop('job controlled-wall quota')))
        signal.setitimer(signal.ITIMER_REAL,self.max_wall)

    def _latest(self):
        if not self.path.exists():return dict(sequence=0,stages={s:dict(RHS=0,states=0,LU=0,controlled_wall_seconds=0.,physics_wall_seconds=0.,GPU_seconds=0.) for s in CAPS},active_wall={})
        with self.path.open('rb') as f:
            f.seek(0,2);size=f.tell();f.seek(max(0,size-50000));lines=f.read().splitlines()
        return json.loads(lines[-1])['after']

    def _update(self,event,*,rhs=0,states=0,lu=0,controlled=0.,physics=0.,gpu=0.,reserve_wall=0.,release_wall=False,**meta):
        fcntl.flock(self.control,fcntl.LOCK_EX)
        try:
            d=self._latest();v=d['stages'][self.stage]
            if v['RHS']+rhs>CAPS[self.stage] or sum(x['RHS'] for x in d['stages'].values())+rhs>TOTAL_RHS:raise BudgetStop('RHS cap')
            used=sum(x['controlled_wall_seconds'] for x in d['stages'].values())
            if reserve_wall and used+sum(d['active_wall'].values())+reserve_wall>CPU_CAP:raise BudgetStop('CPU reservation cap')
            if used+controlled>CPU_CAP or sum(x['GPU_seconds'] for x in d['stages'].values())+gpu>GPU_CAP:raise BudgetStop('CPU/GPU cap')
            if reserve_wall:
                if self.job in d['active_wall']:raise BudgetStop('unreconciled previous job; no automatic replay')
                d['active_wall'][self.job]=reserve_wall
            if release_wall:d['active_wall'].pop(self.job,None)
            for k,a in [('RHS',rhs),('states',states),('LU',lu),('controlled_wall_seconds',controlled),('physics_wall_seconds',physics),('GPU_seconds',gpu)]:v[k]+=a
            d['sequence']+=1
            row=dict(sequence=d['sequence'],event=event,stage=self.stage,job=self.job,utc_unix=time.time(),RHS=rhs,states=states,LU=lu,controlled_wall_seconds=controlled,physics_wall_seconds=physics,context=self.context,**plain(meta),after=d)
            with self.path.open('a') as f:f.write(json.dumps(row,allow_nan=False,separators=(',',':'))+'\n');f.flush();os.fsync(f.fileno())
            write(self.snapshot,d)
            return d
        finally:fcntl.flock(self.control,fcntl.LOCK_UN)

    def reserve_state(self,material,*,resolution,purpose):
        if resolution!=24:raise ValueError('only original TM24 authorized')
        self._update('FULL_FORWARD_ATTEMPT',states=1,rhs=8,lu=1,purpose=purpose);self.credit+=8

    def reserve_rhs(self,count,role):
        if role=='full_forward' and self.credit>=count:self.credit-=count;return
        self._update('SOLVER_RHS',rhs=int(count),role=role)

    @contextmanager
    def span(self,label,**counts):
        start=time.perf_counter();self.counts.update(counts)
        try:yield
        finally:self.counts[label+'_seconds']+=time.perf_counter()-start

    @contextmanager
    def physical(self,label):
        start=time.perf_counter()
        try:yield
        finally:self._update('PHYSICAL_OPERATION',physics=time.perf_counter()-start,label=label)

    def close(self,status='COMPLETE',error=None):
        if self.finished:return
        signal.setitimer(signal.ITIMER_REAL,0);self.finished=True
        seconds=time.perf_counter()-self.start
        d=self._update('JOB_END',controlled=seconds,release_wall=True,status=status,error=error,counts=dict(self.counts),unused_forward_credit=self.credit)
        write(self.jobdir/'CLOSED.json',dict(status=status,error=error,controlled_wall_seconds=seconds,counters=dict(self.counts),stage_cumulative=d['stages'][self.stage]))
        fcntl.flock(self.joblock,fcntl.LOCK_UN);self.joblock.close();self.control.close()

@dataclass
class Point:
    q:np.ndarray
    state:object
    key:str
    T:np.ndarray|None=None
    J:np.ndarray|None=None
    gradient:np.ndarray|None=None

class Target:
    def __init__(self,book,*,keep_fields=False):
        self.book=book;self.keep_fields=keep_fields
        self.cfg=json.loads((OUT/'configs/TARGET.json').read_text())
        points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
        self.prior=MaterialPrior(self.cfg,points);self.sigma=self.cfg['noise']['sigma'];self.model=None;self.evaluations=0

    def backend(self):
        if self.model is None:
            p=Path(os.environ.get('B1_P3L_A30_ROOT',str(ROOT.parents[1]/'A30/a30-maxwell-material-manifold')));sys.path.insert(0,str(p/'src'))
            from a30.physics import FrozenMaxwellAdapter
            with self.book.span('paid_background_array_setup'):
                self.model=FrozenMaxwellAdapter.from_background_npz(OUT/'inputs/background/weak.npz',self.book,purpose_prefix='P3L_'+self.book.stage)
        return self.model

    def evaluate(self,q,purpose='target'):
        chi,_=self.prior.decode(q);self.evaluations+=1;key=f'{self.book.job}_{self.evaluations:07d}'
        folder=self.book.jobdir/'endpoints';self.book.context=dict(endpoint=key,purpose=purpose)
        save(self.book.jobdir/'attempts'/f'{key}.npz',q=q,purpose=np.array(purpose))
        try:
            with self.book.physical('full_forward'):
                state=self.backend().full_state(chi,purpose=purpose)
                if state.residual>1e-10:raise RuntimeError('backward residual contract')
            arrays=dict(q=q,chi=chi,prediction=state.data,backward_residual=state.residual)
            if self.keep_fields:arrays.update(current=state.current,field=state.field)
            save(folder/f'{key}.npz',**arrays)
            return Point(np.asarray(q).copy(),state,key)
        except BaseException as e:
            write(self.book.jobdir/'failures'/f'{key}.json',dict(endpoint=key,error=repr(e),purpose=purpose));raise

    def phi(self,p,y):return float(np.vdot(p.state.data-y,p.state.data-y).real/self.sigma**2)
    def white(self,a):return pack(a)*np.sqrt(2)/self.sigma
    def objective(self,p,y):return self.phi(p,y)+.5*float(p.q[:64]@p.q[:64])

    def jacobian(self,p):
        if p.J is None:
            with self.book.physical('receiver_J_32_RHS'):
                p.T=self.backend().receiver_transfer(p.state)
                C=np.einsum('qp,ps->qsp',p.T,p.state.field).reshape(256,576)
                a=C@self.prior.derivative(p.q)
                p.J=np.column_stack([self.white(a[:,k].reshape(32,8)) for k in range(66)])
        return p.J

    def phi_gradient(self,p,y):
        if p.gradient is None:
            model=self.backend();res=p.state.data-y
            with self.book.physical('source_adjoint_gradient_8_RHS'):
                a=model.solve_anchor(model.S.conj().T@res,p.state,adjoint=True,role='likelihood_gradient')
                h=np.sum(np.conj(a)*p.state.field,axis=1)
                p.gradient=2*np.real(self.prior.derivative(p.q).T@h)/self.sigma**2
        return p.gradient

    def persist_cache(self,p,path):
        save(path,q=p.q,chi=p.state.chi,current=p.state.current,field=p.state.field,prediction=p.state.data,LU=p.state.factor[0],pivots=p.state.factor[1],backward_residual=p.state.residual,receiver_T=np.zeros((0,0),complex) if p.T is None else p.T)

    def observation(self,parent):
        with np.load(OUT/f'inputs/online/{parent}.npz',allow_pickle=False) as a:
            if set(a.files)!={'observation','target_id','sigma'} or str(a['target_id'])!=self.cfg['target_id'] or float(a['sigma'])!=self.sigma:raise ValueError('observation metadata')
            return a['observation'].copy()

def to_probit(q):
    x=np.asarray(q).copy();x[64:]=ndtri(np.clip(x[64:],1e-4,1-1e-4));return x

def from_probit(x):
    q=np.asarray(x).copy();q[64:]=ndtr(q[64:]);return q

def probit_derivative(x):return np.r_[np.ones(64),np.exp(-.5*np.asarray(x)[64:]**2)/np.sqrt(2*np.pi)]
