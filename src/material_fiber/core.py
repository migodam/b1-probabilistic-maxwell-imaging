from __future__ import annotations
from pathlib import Path
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass
import json, os, sys, time, fcntl, signal
import numpy as np
from p1_material.prior import MaterialPrior, structural_change, topology

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'research/p2_material_fiber_20261019'
OLD=ROOT.parent/'b1-p1-material-target/research/p1_material_target_20261018'
CAPS={'A':(600,24000,3600),'B':(1000,34000,7200),'C':(1200,36000,7200),'D':(700,20000,7200),'E':(500,18000,3600)}
TAU=23.785449004057032

def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,Path):return str(x)
    if isinstance(x,complex):return [x.real,x.imag]
    if isinstance(x,dict):return {str(k):plain(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [plain(v) for v in x]
    return x

def write(path,x):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    t=path.with_suffix(path.suffix+'.tmp');t.write_text(json.dumps(plain(x),ensure_ascii=False,indent=2,allow_nan=False)+'\n');t.replace(path)

def save(path,**arrays):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.tmp')
    with temp.open('wb') as f:np.savez_compressed(f,**arrays);f.flush();os.fsync(f.fileno())
    temp.replace(path)

def pack(z):return np.concatenate((z.real.T,z.imag.T),axis=1).ravel()
def mpack(z):return np.r_[z.real,z.imag]/24

class BudgetStop(RuntimeError):pass
class BlockStop(RuntimeError):pass

class Book:
    """Single paid executor; reservation is durable before every native solve."""
    def __init__(self,stage,recovery=None):
        if (OUT/'CAMPAIGN_CLOSED.json').exists():raise BudgetStop('P2 campaign closed; unused allocation cannot be revived')
        self.stage=stage;self.recovery_id=recovery;self.counts=Counter();self.path=OUT/'RESOURCE_LEDGER.jsonl';self.path.parent.mkdir(parents=True,exist_ok=True)
        self.lock=(OUT/'executor.lock').open('a');fcntl.flock(self.lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        self.rows=[json.loads(x) for x in self.path.read_text().splitlines()] if self.path.exists() else []
        self.used={s:Counter() for s in CAPS};self.previous_wall=Counter();self.context={};self.block=None
        for r in self.rows:
            s=r.get('stage')
            if s in self.used and r['event'] in ['STATE','RHS']:self.used[s][r['event']]+=r['amount']
            if s in self.used and r['event']=='JOB_END':self.previous_wall[s]+=r['controlled_wall_seconds']
        self.close_path=OUT/f'{stage}/CLOSED.json'
        if self.close_path.exists():
            if recovery is None:raise BudgetStop('stage already closed; no reset')
            r=json.loads((OUT/f'configs/{stage}_RECOVERY_{recovery}.json').read_text())
            old=json.loads(self.close_path.read_text())
            self.close_path=OUT/f'{stage}/CLOSED_RECOVERY_{recovery}.json'
            permitted=old['status']=='ERROR' or (r.get('baseline_completion',False) and old['status']=='COMPLETE')
            if self.close_path.exists() or not permitted or r['consumed_before']!=dict(self.used[stage]) or not r['same_authorization_and_remaining_caps']:raise BudgetStop('invalid/reused recovery; no reset')
        starts=sum(r['event']=='JOB_START' and r.get('stage')==stage for r in self.rows)
        ends=sum(r['event']=='JOB_END' and r.get('stage')==stage for r in self.rows)
        if starts!=ends:raise BudgetStop('unclosed executor: reconcile recorded attempts before explicit recovery')
        a=json.loads((OUT/'configs/AUTHORIZATION.json').read_text())
        if not a['explicit_user_authorization'] or a['campaign_id']!='P2-material-fiber-20261019':raise BudgetStop('authorization mismatch')
        self.start=time.perf_counter();self.physics_seconds=0.;self.finished=False;self.physical_depth=0
        self.append('JOB_START',authorization='configs/AUTHORIZATION.json',recovery=recovery)
        signal.signal(signal.SIGALRM,lambda a,b:(_ for _ in ()).throw(BudgetStop('controlled stage wall cap')))
        signal.setitimer(signal.ITIMER_REAL,max(.001,CAPS[stage][2]-self.previous_wall[stage]))

    def append(self,event,**kw):
        row=dict(event=event,stage=self.stage,utc_unix=time.time(),context=self.context,**plain(kw))
        with self.path.open('a') as f:f.write(json.dumps(row,ensure_ascii=False,allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())

    def require(self,states=0,RHS=0):
        s=self.stage;u=self.used[s];c=CAPS[s]
        if u['STATE']+states>c[0] or u['RHS']+RHS>c[1]:raise BudgetStop('stage state/RHS hard cap')
        if sum(v['STATE'] for v in self.used.values())+states>4000 or sum(v['RHS'] for v in self.used.values())+RHS>132000:raise BudgetStop('campaign hard cap')
        if self.block:
            bs,br,ns,nr=self.block
            if u['STATE']-bs+states>ns or u['RHS']-br+RHS>nr:raise BlockStop('pre-registered block cap')
        if time.perf_counter()-self.start+self.previous_wall[s]>=c[2]:raise BudgetStop('controlled wall cap')

    def reserve(self,event,n,**meta):
        self.require(states=n if event=='STATE' else 0,RHS=n if event=='RHS' else 0)
        self.append(event,amount=int(n),**meta);self.used[self.stage][event]+=int(n)

    def reserve_state(self,material,*,resolution,purpose):self.reserve('STATE',1,resolution=resolution,purpose=purpose)
    def reserve_rhs(self,count,role):self.reserve('RHS',count,role=role)

    @contextmanager
    def span(self,label,**counts):
        start=time.perf_counter();error=None;self.counts.update(counts)
        try:yield
        except BaseException as e:error=repr(e);raise
        finally:self.append('ACTION',label=label,wall_seconds=time.perf_counter()-start,counts=counts,error=error)

    @contextmanager
    def physical(self,label):
        start=time.perf_counter();outer=self.physical_depth==0;self.physical_depth+=1
        try:yield
        finally:
            dt=time.perf_counter()-start;self.physical_depth-=1
            if outer:self.physics_seconds+=dt
            self.append('PHYSICS_WALL',label=label,wall_seconds=dt,additive=outer)

    @contextmanager
    def scope(self,states,RHS):
        previous=self.block;u=self.used[self.stage];self.block=(u['STATE'],u['RHS'],states,RHS)
        try:yield
        finally:self.block=previous

    def close(self,status='COMPLETE',error=None,final_phase=True):
        if self.finished:return
        self.finished=True;signal.setitimer(signal.ITIMER_REAL,0)
        receipt=dict(status=status,error=error,counts=dict(self.counts),consumed=dict(self.used[self.stage]),controlled_wall_seconds=time.perf_counter()-self.start,physics_operation_wall_seconds=self.physics_seconds)
        self.append('JOB_END',**receipt)
        write(self.close_path if final_phase else OUT/f'{self.stage}/OBSERVATION_JOB_CLOSED.json',receipt)
        fcntl.flock(self.lock,fcntl.LOCK_UN);self.lock.close()


@dataclass
class Point:
    q:np.ndarray
    state:object
    key:str
    T:np.ndarray|None=None
    C:np.ndarray|None=None
    A:np.ndarray|None=None


class Target:
    def __init__(self,book):
        self.book=book;self.cfg=json.loads((OUT/'configs/TARGET.json').read_text())
        points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
        self.prior=MaterialPrior(self.cfg,points);self.sigma=self.cfg['noise']['sigma'];self.target_id=self.cfg['target_id'];self.model=None

    def backend(self):
        if self.model is None:
            p=Path(os.environ.get('B1_P2_A30_ROOT',str(ROOT.parents[1]/'A30/a30-maxwell-material-manifold')));sys.path.insert(0,str(p/'src'))
            from a30.physics import FrozenMaxwellAdapter
            with self.book.physical('saved_background_setup'):
                self.model=FrozenMaxwellAdapter.from_background_npz(OUT/'inputs/background/weak.npz',self.book,purpose_prefix='P2_'+self.book.stage)
        return self.model

    def evaluate(self,q,purpose):
        self.book.require(states=1,RHS=8);chi,_=self.prior.decode(q)
        key=f'{self.book.stage}_{self.book.used[self.book.stage]["STATE"]+1:05d}'
        self.book.context=dict(purpose=purpose,endpoint=key)
        save(OUT/f'{self.book.stage}/attempts/{key}.npz',q=q,chi=chi,purpose=np.array(purpose))
        with self.book.physical('full_forward_and_endpoint_persistence'):
            state=self.backend().full_state(chi,purpose=purpose)
            if state.residual>1e-10:raise RuntimeError('backward residual contract')
            save(OUT/f'{self.book.stage}/endpoints/{key}.npz',q=q,chi=chi,current=state.current,field=state.field,prediction=state.data,backward_residual=state.residual,purpose=np.array(purpose))
        self.book.append('ENDPOINT_PERSISTED',key=key,purpose=purpose)
        return Point(np.asarray(q).copy(),state,key)

    def persist_cache(self,p):
        save(OUT/f'cache/{p.key}.npz',q=p.q,chi=p.state.chi,current=p.state.current,field=p.state.field,prediction=p.state.data,LU=p.state.factor[0],pivots=p.state.factor[1],backward_residual=p.state.residual,receiver_T=np.zeros((0,0),complex) if p.T is None else p.T)

    def load(self,key):
        model=self.backend();path=OUT/f'cache/{key}.npz'
        with np.load(path,allow_pickle=False) as a:d={k:a[k].copy() for k in a.files}
        from a30.physics import _cpu
        chi=d['chi'];L=np.eye(576)-chi[:,None]*model.G
        checks=[np.linalg.norm(L@d['current']-chi[:,None]*model.incident)/max(1e-30,np.linalg.norm(chi[:,None]*model.incident)),np.linalg.norm(d['field']-model.incident-model.G@d['current'])/max(1e-30,np.linalg.norm(d['field'])),np.linalg.norm(model.S@d['current']-d['prediction'])/max(1e-30,np.linalg.norm(d['prediction']))]
        if max(checks)>1e-10 or not np.array_equal(self.prior.decode(d['q'])[0],chi):raise RuntimeError('saved point equation/material mismatch')
        native=_cpu.FullState(chi,d['current'],d['field'],d['prediction'],L,(d['LU'],d['pivots']),float(d['backward_residual']))
        state=model._export_state(native,cached=True,path=path);p=Point(d['q'],state,key)
        if d.get('receiver_T',np.zeros((0,0))).size:
            p.T=d['receiver_T'];model._transfers[state]=p.T;self.compose_J(p)
        self.book.append('CACHE_REUSE',key=key,equation_audit_max=float(max(checks)),incremental_RHS=0,cold_equivalent_forward_RHS=8,cold_equivalent_J_RHS=32 if p.T is not None else 0)
        return p

    def compose_J(self,p):
        p.C=np.einsum('qp,ps->qsp',p.T,p.state.field).reshape(256,576)
        Z=p.C@self.prior.derivative(p.q)
        p.A=np.column_stack([pack(Z[:,k].reshape(32,8))*np.sqrt(2)/self.sigma for k in range(66)])
        return p.A

    def jacobian(self,p):
        if p.A is not None:return p.A
        self.book.require(RHS=32)
        with self.book.physical('receiver_jacobian'):
            p.T=self.backend().receiver_transfer(p.state);self.compose_J(p);self.persist_cache(p)
        return p.A

    def white(self,d):return pack(d)*np.sqrt(2)/self.sigma
    def material_action(self,p,delta):return self.white((p.C@delta).reshape(32,8))
    def phi(self,p,y):return float(np.vdot(p.state.data-y,p.state.data-y).real/self.sigma**2)
    def objective(self,p,y):return self.phi(p,y)+.5*float(p.q[:64]@p.q[:64])
    def qualified(self,p,y):return bool(self.prior.legal(p.q) and 2*self.phi(p,y)<=TAU*TAU and p.state.residual<=1e-10)
    def observation(self,parent):
        with np.load(OUT/f'inputs/online/{parent}.npz',allow_pickle=False) as z:
            if set(z.files)!={'observation','target_id','sigma'} or str(z['target_id'])!=self.target_id or float(z['sigma'])!=self.sigma:raise ValueError('online metadata/keys')
            return z['observation'].copy()

    def record(self,p,parent,method,anchor=None,**kw):
        y=self.observation(parent);aux=self.prior.decode(p.q)[1]
        r=dict(parent=parent,method=method,key=p.key,Phi=self.phi(p,y),objective=self.objective(p,y),white_residual=float(np.sqrt(2*self.phi(p,y))),prior_logdensity=self.prior.logprior(p.q),prior_legal=True,data_feasible=self.qualified(p,y),topology=topology(aux['fraction']),posterior_probability=None,weight_semantics='candidate_only',**kw)
        if anchor is not None:
            structure=structural_change(self.prior,anchor.q,p.q);delta=self.phi(p,y)-self.phi(anchor,y)
            r.update(material_RMS=structure['material_RMS'],structure=structure,DeltaPhi=delta,likelihood_competitive=bool(delta<=2),delta_prior_logdensity=self.prior.logprior(p.q)-self.prior.logprior(anchor.q))
        path=OUT/f'{self.book.stage}/candidates.jsonl'
        with path.open('a') as f:f.write(json.dumps(plain(r),ensure_ascii=False,allow_nan=False)+'\n');f.flush()
        return r
