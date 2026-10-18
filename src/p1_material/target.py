"""P1's declared likelihood and paid A30 interface. Importing this is physics-free."""
from pathlib import Path
from contextlib import contextmanager
from collections import Counter
import json,os,time,fcntl,signal,sys
import numpy as np
from .prior import MaterialPrior

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'research/p1_material_target_20261018'


class BudgetDenied(RuntimeError):pass


def pack(z):return np.concatenate((z.real.T,z.imag.T),axis=1).ravel()


def write_json(path,data):
    def plain(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,dict):return {k:plain(x) for k,x in v.items()}
        if isinstance(v,(list,tuple)):return [plain(x) for x in v]
        return v
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(plain(data),ensure_ascii=False,indent=2,allow_nan=False)+'\n');temp.replace(path)


class PhysicsBook:
    """The current authorization is zero; no backend import before this check."""
    def __init__(self,authorization=None):
        self.authorization=authorization or OUT/'configs/EXECUTION_AUTHORIZATION.json'
        self.a=json.loads(Path(self.authorization).read_text());self.counts=Counter();self.start=None;self.context={};self.path=OUT/'PHYSICS_LEDGER.jsonl'
    def rows(self):return [json.loads(x) for x in self.path.read_text().splitlines()] if self.path.exists() else []
    def totals(self):
        rows=self.rows();return dict(states=sum(x.get('amount',0) for x in rows if x['event']=='STATE'),RHS=sum(x.get('amount',0) for x in rows if x['event']=='RHS'))
    def append(self,row):
        with self.path.open('a') as f:f.write(json.dumps(dict(utc=time.time(),context=self.context,**row),allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
    def require(self,states=0,RHS=0):
        cap_s=self.a.get('max_new_Maxwell_states',0);cap_r=self.a.get('max_new_complex_RHS',0)
        if cap_s<=0 or cap_r<=0:raise BudgetDenied('No new P1 physics allocation; zero backend imports/solves permitted')
        if not self.a.get('explicit_user_physics_authorization',False):raise BudgetDenied('Numeric allocation must be accompanied by explicit user authorization')
        target_id=json.loads((OUT/'configs/ONLINE_TARGET.json').read_text())['target_id']
        if self.a.get('target_id')!=target_id or self.a.get('stage')!='P1_physics_qualification_20261018':raise BudgetDenied('allocation belongs to a different target or stage')
        t=self.totals()
        if t['states']+states>cap_s or t['RHS']+RHS>cap_r:raise BudgetDenied('P1 total hard cap; no old-stage balances')
        if self.start is None:
            if (OUT/'PHYSICS_CLOSED.json').exists():raise BudgetDenied('physical stage closed')
            if any(x['event']=='JOB_START' for x in self.rows()):raise BudgetDenied('existing job requires explicit recovery; never reset ledger')
            self.start=time.perf_counter();self.append(dict(event='JOB_START',authorization=self.a))
            def alarm(a,b):raise BudgetDenied('controlled physics wall cap')
            signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,float(self.a['max_controlled_wall_seconds']))
        if time.perf_counter()-self.start>=self.a['max_controlled_wall_seconds']:raise BudgetDenied('P1 physical wall cap')
    def reserve(self,event,amount,**meta):
        with (OUT/'physics.lock').open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX);self.require(states=amount if event=='STATE' else 0,RHS=amount if event=='RHS' else 0)
            self.append(dict(event=event,amount=int(amount),**meta))
    def reserve_state(self,material,*,resolution,purpose):self.reserve('STATE',1,resolution=resolution,purpose=purpose)
    def reserve_rhs(self,count,role):self.reserve('RHS',int(count),role=role)
    @contextmanager
    def span(self,label,**counts):
        start=time.perf_counter();error=None;self.counts.update(counts)
        try:yield
        except BaseException as e:error=repr(e);raise
        finally:self.append(dict(event='ACTION',label=label,counts=counts,wall_seconds=time.perf_counter()-start,error=error))
    def close(self,error=None):
        if self.start is None:return
        signal.setitimer(signal.ITIMER_REAL,0);elapsed=time.perf_counter()-self.start
        self.append(dict(event='JOB_END',wall_seconds=elapsed,totals=self.totals(),counts=dict(self.counts),error=error))
        write_json(OUT/'PHYSICS_CLOSED.json',dict(totals=self.totals(),wall_seconds=elapsed,error=error))


class MaxwellTarget:
    def __init__(self,book=None):
        self.config=json.loads((OUT/'configs/ONLINE_TARGET.json').read_text());points=json.loads((OUT/'inputs/geometry.json').read_text())['points']
        self.prior=MaterialPrior(self.config,points);self.book=book or PhysicsBook();self.model=None;self.serial=0
        self.sigma=self.config['noise']['sigma'];self.target_id=self.config['target_id']
    def _model(self,states=0,RHS=0):
        self.book.require(states,RHS)  # guard BEFORE A30 or its vendored imports
        if self.model is None:
            backend=Path(os.environ.get('B1_P1_A30_ROOT',str(ROOT.parents[1]/'A30/a30-maxwell-material-manifold')))
            sys.path.insert(0,str(backend/'src'));from a30.physics import FrozenMaxwellAdapter
            self.model=FrozenMaxwellAdapter.from_background_npz(OUT/'inputs/background/weak.npz',self.book,purpose_prefix='P1_QUALIFICATION')
        return self.model
    def forward(self,q,purpose='candidate'):
        model=self._model(states=1,RHS=8);chi,_=self.prior.decode(q);self.book.context=dict(purpose=purpose)
        state=model.full_state(chi,purpose=purpose);path=OUT/f'results/physical/endpoints/state_{self.serial:04}.npz';self.serial+=1;path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('wb') as f:
            np.savez_compressed(f,q=q,chi=state.chi,current=state.current,field=state.field,prediction=state.data,LU=state.factor[0],pivots=state.factor[1],backward_residual=state.residual);f.flush();os.fsync(f.fileno())
        self.book.append(dict(event='ENDPOINT_PERSISTED',array=str(path.relative_to(ROOT)),purpose=purpose))
        return state
    def jacobian(self,q,state):
        model=self._model(RHS=32);T=model.receiver_transfer(state)
        C=np.einsum('qp,ps->qsp',T,state.field).reshape(256,576);D=self.prior.derivative(q);Z=C@D
        A=np.column_stack([pack(Z[:,k].reshape(32,8))*np.sqrt(2)/self.sigma for k in range(66)])
        path=OUT/f'results/physical/jacobians/J_{self.serial:04}.npz';self.serial+=1;path.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(path,q=q,T=T,C=C,D=D,A=A)
        return A
    def phi(self,prediction,observation):
        d=np.asarray(prediction)-np.asarray(observation)
        if d.shape!=(32,8):raise ValueError('frozen acquisition shape')
        return float(np.vdot(d,d).real/self.sigma**2)
    def log_likelihood(self,prediction,observation):return -256*np.log(np.pi*self.sigma**2)-self.phi(prediction,observation)
    def log_target(self,q,prediction,observation):return self.prior.logprior(q)+self.log_likelihood(prediction,observation)
    def load_observation(self,parent_id):
        path=OUT/f'inputs/online/{parent_id}.npz'
        if not path.exists():raise FileNotFoundError('P1 actual scattering observation NOT_RUN; no proxy/truth substitution')
        with np.load(path,allow_pickle=False) as z:
            if set(z.files)!={'observation','sigma','target_id'}:raise ValueError('online isolation contract')
            if str(z['target_id'])!=self.target_id or float(z['sigma'])!=self.sigma:raise ValueError('target identity/noise mismatch')
            return z['observation'].copy()
