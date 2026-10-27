"""Paid R6 adapter: original P1/R4 mathematics, independent reservations.

One CPU physics worker, bounded LU ownership. Every attempted state and every
receiver/source adjoint is reserved before the underlying solve. Interrupted
actions fail closed; neither failures nor cache misses refund reservations.
"""
from __future__ import annotations
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_k]='1'
import sys,json,time,uuid,fcntl,resource,signal
from pathlib import Path
from collections import Counter
from contextlib import contextmanager
import numpy as np
sys.dont_write_bytecode=True
OUT=Path(__file__).resolve().parents[1]
REPO=OUT.parents[1]
GAUSSIAN=REPO.parents[1]
R4=REPO.parent/'b1-p3l-r4-conditional-transport'
sys.path.insert(0,str(R4/'src'))
from p3r4.core import Target as OriginalTarget, Point, BudgetStop, plain, pack, to_probit, from_probit, probit_derivative
from p3r4.inversion import reconstruct
from p3r4.anchors import probit_optimize
SIGMA=0.00119211798514519

def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(plain(obj),allow_nan=False,ensure_ascii=False,indent=2)+'\n');temp.replace(path)

def save(path,**arrays):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    with temp.open('wb') as f:np.savez_compressed(f,**arrays)
    temp.replace(path)

def totals():
    p=OUT/'RESOURCE_LEDGER_SNAPSHOT.json'
    if p.exists():return json.loads(p.read_text())
    return dict(sequence=0,RHS=0,forward_RHS=0,derivative_RHS=0,states=0,LU=0,
                physics_wall_seconds=0.,controlled_wall_seconds=0.,CPU_seconds=0.,
                failures=0,stages={},active_actions={},active_jobs={})

class Book:
    def __init__(self,stage,job,parent=None,role='ONLINE_INFERENCE',max_wall=7200):
        if (OUT/'CLOSURE.json').exists():raise BudgetStop('R6 closed')
        if not (OUT/'configs/P0_FREEZE.json').exists():raise BudgetStop('P0 must freeze before paid physics')
        self.stage,self.job,self.context=stage,job,dict(parent=parent,role=role)
        self.jobdir=OUT/'raw/jobs'/job;self.jobdir.mkdir(parents=True,exist_ok=True)
        if (self.jobdir/'START.json').exists():raise RuntimeError('existing job: reconcile, no automatic replay')
        self.lock=(OUT/'shared/physics_worker.lock').open('a')
        fcntl.flock(self.lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        self.counts=Counter();self.credit=0;self.active=None;self.finished=False
        self.started=time.perf_counter();self.cpu0=time.process_time();self.max_wall=float(max_wall)
        if totals()['active_actions'] or totals()['active_jobs']:
            raise RuntimeError('unreconciled interrupted action/job; explicit audit required')
        self._event('JOB_BEGIN',max_wall=self.max_wall)
        write(self.jobdir/'START.json',dict(stage=stage,job=job,pid=os.getpid(),utc=time.time(),role=role))
        signal.signal(signal.SIGALRM,lambda a,b:(_ for _ in ()).throw(BudgetStop('job controlled wall cap')))
        signal.setitimer(signal.ITIMER_REAL,self.max_wall)

    @staticmethod
    def totals():return totals()

    def _event(self,event,*,rhs=0,states=0,lu=0,kind=None,**fields):
        with (OUT/'shared/ledger.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            t=totals();a=json.loads((OUT/'configs/AUTHORIZATION.json').read_text())
            v=t['stages'].setdefault(self.stage,dict(RHS=0,states=0,LU=0,physics_wall_seconds=0.,controlled_wall_seconds=0.))
            if any(x<0 for x in (rhs,states,lu)):raise ValueError('no refunds')
            if t['RHS']+rhs>a['total_complex_RHS'] or v['RHS']+rhs>a['stage_caps'][self.stage]:raise BudgetStop('R6 RHS hard cap')
            if event=='JOB_BEGIN':
                if t['controlled_wall_seconds']+fields['max_wall']>a['controlled_CPU_wall_seconds']:raise BudgetStop('CPU wall reservation')
                t['active_jobs'][self.job]=dict(stage=self.stage,max_wall=fields['max_wall'])
            if event=='PHYSICS_BEGIN':
                if t['active_actions']:raise RuntimeError('physics action already live')
                t['active_actions'][fields['action']]=dict(job=self.job,stage=self.stage)
            if event=='PHYSICS_END':
                t['active_actions'].pop(fields['action'])
                sec=fields['wall_seconds'];t['physics_wall_seconds']+=sec;v['physics_wall_seconds']+=sec
                t['failures']+=int(fields.get('error') is not None)
            if event=='JOB_END':
                t['active_jobs'].pop(self.job)
                sec=fields['wall_seconds'];t['controlled_wall_seconds']+=sec;v['controlled_wall_seconds']+=sec
                t['CPU_seconds']+=fields['CPU_seconds']
            t['RHS']+=rhs;v['RHS']+=rhs;t['states']+=states;v['states']+=states;t['LU']+=lu;v['LU']+=lu
            t['forward_RHS' if kind=='forward' else 'derivative_RHS']+=rhs;t['sequence']+=1
            row=dict(sequence=t['sequence'],event=event,stage=self.stage,job=self.job,utc=time.time(),
                     context=self.context,RHS=rhs,states=states,LU=lu,kind=kind,total_RHS=t['RHS'],**fields)
            with (OUT/'RESOURCE_LEDGER.jsonl').open('a') as stream:
                stream.write(json.dumps(plain(row),allow_nan=False,separators=(',',':'))+'\n');stream.flush()
            write(OUT/'RESOURCE_LEDGER_SNAPSHOT.json',t)

    def reserve_state(self,material,*,resolution,purpose):
        if resolution!=24:raise ValueError('only original TM24')
        self._event('STATE_RESERVED',rhs=8,states=1,lu=1,kind='forward',purpose=purpose)
        self.credit+=8

    def reserve_rhs(self,count,role):
        count=int(count)
        if role=='full_forward':
            if self.credit<count:raise RuntimeError('missing forward reservation')
            self.credit-=count;self.counts['actual_forward_RHS']+=count
        else:
            self._event('DERIVATIVE_RESERVED',rhs=count,kind='derivative',solver_role=role)
            self.counts['actual_derivative_RHS']+=count

    @contextmanager
    def span(self,label,**counts):
        start=time.perf_counter();self.counts.update(counts)
        try:yield
        finally:self.counts[label+'_seconds']+=time.perf_counter()-start

    @contextmanager
    def physical(self,label):
        if self.active is not None:raise RuntimeError('nested physics')
        action=uuid.uuid4().hex;self.active=action;self._event('PHYSICS_BEGIN',action=action,label=label)
        start=time.perf_counter();cpu=time.process_time();err=None
        try:yield
        except BaseException as exc:err=repr(exc);raise
        finally:
            self._event('PHYSICS_END',action=action,label=label,wall_seconds=time.perf_counter()-start,CPU_seconds=time.process_time()-cpu,error=err)
            self.active=None

    def close(self,status='COMPLETE',error=None):
        if self.finished:return
        signal.setitimer(signal.ITIMER_REAL,0)
        receipt=dict(status=status,error=error,wall_seconds=time.perf_counter()-self.started,
                     CPU_seconds=time.process_time()-self.cpu0,peak_RSS_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
                     unused_prepaid_RHS=self.credit,actual_counters=dict(self.counts))
        self._event('JOB_END',**receipt);write(self.jobdir/'CLOSED.json',receipt)
        self.finished=True;fcntl.flock(self.lock,fcntl.LOCK_UN);self.lock.close()

class Target(OriginalTarget):
    def __init__(self,book,*,sigma=SIGMA,mask=None):
        super().__init__(book,keep_fields=False)
        self.sigma=float(sigma)
        self.mask=np.ones((32,8),bool) if mask is None else np.asarray(mask,bool).copy()
        if self.mask.shape!=(32,8) or not self.mask.any():raise ValueError('acquisition mask')
        self._packed_mask=np.concatenate((self.mask.T,self.mask.T),axis=1).ravel()

    def white(self,a):return (pack(a)*np.sqrt(2)/self.sigma)[self._packed_mask]
    def phi(self,p,y):
        r=(p.state.data-y)[self.mask]
        return float(np.vdot(r,r).real/self.sigma**2)

    def evaluate(self,q,purpose='target'):
        # Store lossless coordinates+response per query; the exact deterministic
        # decoder reconstructs material. Final endpoints separately store fields.
        q=np.asarray(q,float);chi,_=self.prior.decode(q);self.evaluations+=1
        key=f'{self.book.job}_{self.evaluations:07d}'
        self.book.context.update(endpoint=key,purpose=purpose)
        try:
            with self.book.physical('full_forward'):
                state=self.backend().full_state(chi,purpose=purpose)
                if state.residual>1e-10:raise RuntimeError('backward residual >1e-10')
            save(self.book.jobdir/'endpoints'/f'{key}.npz',q=q,prediction=state.data,backward_residual=state.residual)
            return Point(q.copy(),state,key)
        except BaseException as exc:
            write(self.book.jobdir/'failures'/f'{key}.json',dict(q=q,error=repr(exc),purpose=purpose));raise

    def persist_cache(self,p,path):
        save(path,q=p.q,chi=p.state.chi,current=p.state.current,field=p.state.field,prediction=p.state.data,
             backward_residual=p.state.residual,J=np.zeros((0,0)) if p.J is None else p.J,
             receiver_T=np.zeros((0,0),complex) if p.T is None else p.T)

    def phi_gradient(self,p,y):
        if self.mask.all():return super().phi_gradient(p,y)
        if p.J is not None:return p.J.T@self.white(p.state.data-y)
        residual=np.where(self.mask,p.state.data-y,0.)
        with self.book.physical('masked_source_adjoint_8_RHS'):
            a=self.backend().solve_anchor(self.backend().S.conj().T@residual,p.state,adjoint=True,role='likelihood_gradient')
            h=np.sum(np.conj(a)*p.state.field,axis=1)
            return 2*np.real(self.prior.derivative(p.q).T@h)/self.sigma**2
