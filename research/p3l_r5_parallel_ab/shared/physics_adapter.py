"""R5 original-P1 Gaussianized decoder and counted TM24 CPU physics.

Use ROOT on sys.path and import this module by its research namespace. No old
Book/Target is imported, and no history is written. All actual solves, including
qualification, use the cumulative R5 reservations. No automatic cache lookup.
"""
from __future__ import annotations
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
from pathlib import Path
from contextlib import contextmanager
from collections import Counter, OrderedDict
from dataclasses import dataclass
import fcntl,json,sys,time,uuid
import numpy as np
from scipy.special import ndtr

ROOT=Path(__file__).resolve().parents[3]
SHARED=Path(__file__).resolve().parent
OUT=SHARED.parent
HISTORY=ROOT.parent/'b1-p3l-r4-conditional-transport'
P1=ROOT/'research/p1_material_target_20261018'
sys.path.insert(0,str(ROOT/'src'))
from p1_material.prior import MaterialPrior
SIGMA=0.00119211798514519

def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.'+uuid.uuid4().hex+'.tmp')
    with tmp.open('w') as f:json.dump(obj,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    tmp.replace(path)

def save(path,**arrays):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.'+uuid.uuid4().hex+'.tmp')
    with tmp.open('wb') as f:np.savez_compressed(f,**arrays);f.flush();os.fsync(f.fileno())
    tmp.replace(path)

def pack(z):
    z=np.asarray(z)
    if z.shape!=(32,8):raise ValueError('source-major complex32x8 required')
    return np.concatenate((z.real.T,z.imag.T),axis=1).ravel()

def white(z):return pack(z)*np.sqrt(2)/SIGMA

def native(x):
    x=np.asarray(x,dtype=float)
    if x.shape!=(66,) or not np.isfinite(x).all():raise ValueError('finite Gaussianized66 required')
    q=x.copy();q[64:]=ndtr(q[64:]);return q

def chain(x):return np.r_[np.ones(64),np.exp(-.5*np.asarray(x)[64:]**2)/np.sqrt(2*np.pi)]

class BudgetStop(RuntimeError):pass

class Ledger:
    def __init__(self,stage,job):
        if stage not in ('A','B'):raise ValueError('stage A or B')
        if not job or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in job):raise ValueError('safe explicit job name')
        self.stage,self.job=stage,job
        self.path=SHARED/'cost_ledger'/'LEDGER.jsonl';self.path.parent.mkdir(parents=True,exist_ok=True)
        self.counts=Counter();self.context={};self.kind=None;self._offset=0
        self._summary=dict(states=0,forward_rhs=0,derivative_rhs=0,cpu_seconds=0.,wall_seconds=0.,actions=0,failures=0);self._active={}
    def rows(self):
        if not self.path.exists():return []
        # A partial final line fails closed; never silently reset or ignore it.
        return [json.loads(line) for line in self.path.read_text().splitlines() if line]
    def append(self,row):
        with (self.path.parent/'journal.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            with self.path.open('a') as f:
                f.write(json.dumps(dict(utc=time.time(),stage=self.stage,job=self.job,**row),allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
    def totals(self,rows=None):
        # Incremental replay of the append-only journal. Restart replays once;
        # subsequent reservations read only new complete records, bounded RAM.
        if rows is not None:raise ValueError('external ledger rows unsupported')
        if self.path.exists():
            if self.path.stat().st_size<self._offset:raise BudgetStop('ledger truncated; reset forbidden')
            with (self.path.parent/'journal.lock').open('a') as lock, self.path.open('r') as f:
                fcntl.flock(lock,fcntl.LOCK_SH)
                f.seek(self._offset)
                for line in f:
                    if not line.endswith('\n'):raise BudgetStop('partial ledger tail; recovery required')
                    r=json.loads(line)
                    if r.get('stage')!=self.stage:continue
                    t=self._summary
                    if r['event']=='RESERVE':
                        for k in ('states','forward_rhs','derivative_rhs'):t[k]+=r.get(k,0)
                    if r['event']=='COMPUTE_BEGIN':self._active[r['id']]=r;t['actions']+=1
                    if r['event']=='COMPUTE_END':
                        self._active.pop(r['id'],None)
                        for k in ('cpu_seconds','wall_seconds'):t[k]+=r[k]
                        t['failures']+=r.get('error') is not None
                self._offset=f.tell()
        t=self._summary.copy()
        elapsed=sum(max(0,time.time()-r['utc']) for r in self._active.values())
        t['cpu_seconds']+=elapsed;t['wall_seconds']+=elapsed
        t['all_rhs']=t['forward_rhs']+t['derivative_rhs'];t['unfinished_actions']=len(self._active)
        return t
    def _check(self,t,states=0,forward_rhs=0,derivative_rhs=0):
        closure=SHARED/'cost_ledger'/'CLOSURE.json'
        if closure.exists() and self.stage in json.loads(closure.read_text()).get('closed_stages',[]):
            raise BudgetStop('R5 stage is closed; no implicit paid replay or budget reset')
        auth=json.loads((SHARED/'cost_ledger/AUTHORIZATION.json').read_text())[self.stage]
        if self.stage=='A':
            if t['states']+states>auth['forward_states'] or t['forward_rhs']+forward_rhs>auth['forward_rhs'] or t['derivative_rhs']+derivative_rhs>auth['derivative_rhs'] or t['all_rhs']+forward_rhs+derivative_rhs>auth['all_rhs']:raise BudgetStop('R5 A cumulative hard cap')
            if max(t['cpu_seconds'],t['wall_seconds'])>=auth['cpu_wall_seconds']:raise BudgetStop('R5 A cumulative CPU or controlled-wall cap')
        else:
            if t['all_rhs']+forward_rhs+derivative_rhs>auth['confirmation_rhs']:raise BudgetStop('R5 B cumulative confirmation RHS cap')
            if max(t['cpu_seconds'],t['wall_seconds'])>=auth['controlled_wall_seconds']:raise BudgetStop('R5 B cumulative CPU or controlled-wall cap')
    def reserve(self,**amount):
        with (self.path.parent/'ledger.lock').open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX);self._check(self.totals(),**amount)
            self.append(dict(event='RESERVE',role=self.kind,**amount,context=self.context))
    def reserve_state(self,material,*,resolution,purpose):
        if self.kind!='forward' or resolution!=24:raise ValueError('state outside public forward action')
        # Check complete forward affordability before LU, then actual RHS is
        # separately reserved immediately before its solve by frozen A30.
        self._check(self.totals(),states=1,forward_rhs=8)
        self.reserve(states=1)
    def reserve_rhs(self,count,role):
        if self.kind is None:raise ValueError('solve outside public action')
        self.reserve(**{('forward_rhs' if self.kind=='forward' else 'derivative_rhs'):int(count)})
    @contextmanager
    def span(self,label,**counts):
        start=time.perf_counter();cpu=time.process_time();err=None;self.counts.update(counts)
        try:yield
        except BaseException as exc:err=repr(exc);raise
        finally:self.append(dict(event='BACKEND_ACTION',label=label,counts=counts,wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,error=err))
    @contextmanager
    def compute(self,kind,purpose):
        # One physics action at a time across A/B, including all setup/arithmetic.
        with (self.path.parent/'physics_worker.lock').open('a') as worker:
            fcntl.flock(worker,fcntl.LOCK_EX)
            self._check(self.totals())
            ident=uuid.uuid4().hex;self.kind=kind;self.context=dict(purpose=purpose,action=ident)
            self.append(dict(event='COMPUTE_BEGIN',id=ident,kind=kind,purpose=purpose))
            cpu=time.process_time();wall=time.perf_counter();err=None
            try:
                # NumPy and SciPy in the frozen local runtime both link Apple's
                # Accelerate; VECLIB_MAXIMUM_THREADS=1 is set before imports.
                yield
            except BaseException as exc:err=repr(exc);raise
            finally:
                self.append(dict(event='COMPUTE_END',id=ident,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,error=err,counters=dict(self.counts)))
                self.kind=None

@dataclass
class Point:
    x:np.ndarray
    q:np.ndarray
    state:object
    key:str
    @property
    def prediction(self):return self.state.data
    @property
    def chi(self):return self.state.chi

class _OneTransfer(OrderedDict):
    def __setitem__(self,k,v):
        super().__setitem__(k,v)
        while len(self)>1:self.popitem(last=False)

class Target:
    def __init__(self,stage='A',job='qualification',*,selection_frozen=None,output_dir=None):
        self.book=Ledger(stage,job);self.stage=stage
        if stage=='B':
            if selection_frozen is None or not Path(selection_frozen).is_file():raise BudgetStop('B physics requires concrete frozen selection/protocol file')
            self.book.append(dict(event='B_FROZEN_PROTOCOL',path=str(Path(selection_frozen).resolve())))
        self.cfg=json.loads((P1/'configs/ONLINE_TARGET.json').read_text())
        points=json.loads((P1/'inputs/geometry.json').read_text())['points']
        self.prior=MaterialPrior(self.cfg,points);self.sigma=SIGMA;self.model=None
        if self.cfg['noise']['sigma']!=SIGMA:raise ValueError('P1 noise drift')
        self.jobdir=Path(output_dir).resolve() if output_dir is not None else SHARED/'jobs'/stage/job
        if not self.jobdir.resolve().is_relative_to(OUT):raise ValueError('R5-only output directory')
        self.jobdir.mkdir(parents=True,exist_ok=True)
    def decoder(self,x):
        with self.book.span('decoder_only_decode',decoder_only_decodes=1):return self.prior.decode(native(x))
    def decoder_jacobian(self,x):
        with self.book.span('decoder_only_derivative',decoder_only_derivative_builds=1):return self.prior.derivative(native(x))*chain(x)[None,:]
    def backend(self):
        if self.model is None:
            p=ROOT.parents[1]/'A30/a30-maxwell-material-manifold'
            sys.path.insert(0,str(p/'src'))
            from a30.physics import FrozenMaxwellAdapter
            background=HISTORY/'research/p3l_r4_conditional_transport/inputs/background/weak.npz'
            self.model=FrozenMaxwellAdapter.from_background_npz(background,self.book,purpose_prefix='R5_'+self.stage)
            self.model._transfers=_OneTransfer()
        return self.model
    def evaluate(self,x,purpose='candidate',*,retain_reference=False):
        x=np.asarray(x,float);q=native(x);chi,_=self.decoder(x)
        key=uuid.uuid4().hex
        save(self.jobdir/'attempts'/f'{key}.npz',x=x,q=q,purpose=np.array(purpose))
        with self.book.compute('forward',purpose):
            state=self.backend().full_state(chi,purpose=purpose)
            if state.residual>1e-10:raise RuntimeError('backward residual >1e-10')
            arrays=dict(x=x,q=q,chi=chi,prediction=state.data,backward_residual=np.array(state.residual))
            if retain_reference:arrays.update(field=state.field,current=state.current,LU=state.factor[0],pivots=state.factor[1])
            save(self.jobdir/'endpoints'/f'{key}.npz',**arrays)
            self.book.append(dict(event='ENDPOINT_PERSISTED',key=key,relative_path=str((self.jobdir/'endpoints'/f'{key}.npz').relative_to(ROOT)),retained_reference=retain_reference))
        return Point(x.copy(),q,state,key)
    def jacobian(self,p,V=None,*,purpose='jacobian'):
        """White derivative wrt Gaussianized x. V66xd direct JVP uses8*d RHS;
        V=None receiver transfer uses32 RHS and yields512x66. No FD physics.
        """
        D=self.decoder_jacobian(p.x)
        if V is not None:
            V=np.asarray(V,float)
            if V.ndim==1:V=V[:,None]
            if V.ndim!=2 or V.shape[0]!=66 or not np.isfinite(V).all():raise ValueError('V finite66xd')
            D=D@V
        with self.book.compute('derivative',purpose):
            model=self.backend()
            if V is None:
                T=model.receiver_transfer(p.state)
                Z=np.einsum('qn,ns,nd->qsd',T,p.state.field,D,optimize=True)
            else:Z=model.anchor_jvp(D,p.state)
            J=np.column_stack([white(Z[:,:,k]) for k in range(D.shape[1])])
            save(self.jobdir/'derivatives'/f'{p.key}_{uuid.uuid4().hex}.npz',x=p.x,V=np.eye(66) if V is None else V,J=J)
        return J
    def jvp(self,p,v,*,purpose='jvp'):return self.jacobian(p,np.asarray(v)[:,None],purpose=purpose)[:,0]
    def phi(self,p,y):
        prediction=p.prediction if isinstance(p,Point) else np.asarray(p)
        return .5*float(white(prediction-np.asarray(y))@white(prediction-np.asarray(y)))
    def white(self,z):return white(z)
    def gradient(self,p,y,*,purpose='adjoint_gradient'):
        y=np.asarray(y)
        if y.shape!=(32,8) or not np.isfinite(y).all():raise ValueError('observation complex32x8')
        with self.book.compute('derivative',purpose):
            m=self.backend();a=m.solve_anchor(m.S.conj().T@(p.prediction-y),p.state,adjoint=True,role='likelihood_gradient')
            h=np.sum(np.conj(a)*p.state.field,axis=1)
            g=2*np.real(self.decoder_jacobian(p.x).T@h)/SIGMA**2
            save(self.jobdir/'derivatives'/f'{p.key}_{uuid.uuid4().hex}.npz',x=p.x,gradient_phi=g)
        return g
    def inverse_action(self,p,rhs,*,adjoint=False,purpose='inverse_action',persist_path=None):
        """Counted generic L^-1 RHS (or L^-H), complex576xn.

        In particular inverse_action(p, eye(576)) costs576 derivative RHS.
        No new state is created; reacquire a fresh paid Point if needed. Optional
        persistence retains the returned operator/reference, not hidden LU caches.
        """
        rhs=np.asarray(rhs,dtype=np.complex128)
        if rhs.ndim not in (1,2) or rhs.shape[0]!=576 or not np.isfinite(rhs).all():raise ValueError('finite complex576xn RHS required')
        if persist_path is not None:
            persist_path=Path(persist_path).resolve()
            if not persist_path.is_relative_to(OUT):raise ValueError('R5-only inverse reference output')
        with self.book.compute('derivative',purpose):
            value=self.backend().solve_anchor(rhs,p.state,adjoint=adjoint,role='R5_inverse_action')
            if persist_path is not None:
                save(persist_path,x=p.x,value=value,adjoint=np.array(adjoint),rhs_columns=np.array(1 if rhs.ndim==1 else rhs.shape[1]),purpose=np.array(purpose))
                self.book.append(dict(event='INVERSE_REFERENCE_PERSISTED',key=p.key,relative_path=str(persist_path.relative_to(ROOT)),rhs_columns=1 if rhs.ndim==1 else rhs.shape[1],adjoint=bool(adjoint)))
        return value

    def load_observation(self,path):
        path=Path(path)
        if any(s.lower() in ('oracle','offline','labels','teacher') or 'locked' in s.lower() for s in path.parts):raise ValueError('online/oracle path isolation')
        with np.load(path,allow_pickle=False) as z:
            if set(z.files)!={'observation','target_id','sigma'}:raise ValueError('online observation schema')
            if str(z['target_id'])!=self.cfg['target_id'] or float(z['sigma'])!=SIGMA:raise ValueError('observation identity/noise')
            y=z['observation'].copy()
        if y.shape!=(32,8) or y.dtype!=np.complex128 or not np.isfinite(y).all():raise ValueError('observation layout')
        return y
