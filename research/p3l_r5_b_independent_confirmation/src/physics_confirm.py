"""Independent paid adapter; unchanged R4 statistical/optimizer mathematics.

Only the ledger, persistence policy and output destinations are new. Historical
Book instances are never constructed. One worker owns physics, with four live
optimizer points and two receiver-transfer references as in R4.
"""
from __future__ import annotations
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_k]='1'
import sys,json,time,uuid,fcntl,resource
from pathlib import Path
from collections import Counter
from contextlib import contextmanager
import numpy as np

OUT=Path(__file__).resolve().parents[1]
REPO=OUT.parents[1]
R4=REPO.parent/'b1-p3l-r4-conditional-transport'
sys.dont_write_bytecode=True
sys.path.insert(0,str(R4/'src'))
from p3r4.core import Target as OriginalTarget, BudgetStop, plain, save, write, pack
from p3r4.inversion import reconstruct
from p3r4.anchors import probit_optimize

class Book:
    """Reserve every attempt before LU, consume exactly eight forward columns.

    A failed reserved action remains paid. Complete states+RHS reservation and
    journal append occur under one file lock. An interrupted action fails closed
    on restart until explicitly reconciled; there is no implicit budget reset.
    """
    def __init__(self,job,parent,role):
        if not (OUT/'configs/FREEZE.json').exists():
            raise RuntimeError('freeze before any paid physics')
        self.job=job;self.stage='B_CONFIRM';self.jobdir=OUT/'raw/physics'/job
        self.jobdir.mkdir(parents=True,exist_ok=True)
        self.counts=Counter();self.context=dict(parent=parent,role=role)
        self.credit=0;self.finished=False;self.active=None
        self.control=(OUT/'raw/physics_worker.lock').open('a')
        if not self._lock():raise RuntimeError('one physical worker only')
        self.started=time.perf_counter();self.cpu0=time.process_time()
        self._event('JOB_BEGIN')
    def _lock(self):
        try:fcntl.flock(self.control,fcntl.LOCK_EX|fcntl.LOCK_NB);return True
        except BlockingIOError:return False
    @staticmethod
    def totals():
        f=OUT/'RESOURCE_LEDGER_SNAPSHOT.json'
        if f.exists():return json.loads(f.read_text())
        return dict(RHS=0,forward_RHS=0,derivative_RHS=0,state_attempts=0,
                    LU_attempts=0,physics_wall_seconds=0.,physics_CPU_seconds=0.,
                    job_wall_seconds=0.,job_CPU_seconds=0.,failures=0,sequence=0,
                    active_actions={},actual_counters={})
    def _event(self,event,*,rhs=0,states=0,lu=0,kind=None,**fields):
        with (OUT/'raw/ledger.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            t=self.totals();auth=json.loads((OUT/'configs/AUTHORIZATION.json').read_text())
            if t['RHS']+rhs>auth['complex_RHS_hard_cap']:
                raise BudgetStop('independent confirmation hard RHS cap')
            if rhs<0 or states<0 or lu<0:raise ValueError('no budget refunds')
            t['RHS']+=rhs;t['state_attempts']+=states;t['LU_attempts']+=lu
            t['forward_RHS' if kind=='forward' else 'derivative_RHS']+=rhs
            t['sequence']+=1
            if event=='PHYSICS_BEGIN':
                if t['active_actions']:raise RuntimeError('unreconciled interrupted physics; stop')
                t['active_actions'][fields['action']]=dict(job=self.job,utc=time.time(),**self.context)
            if event=='PHYSICS_END':
                t['active_actions'].pop(fields['action'])
                t['physics_wall_seconds']+=fields['wall_seconds'];t['physics_CPU_seconds']+=fields['CPU_seconds']
                t['failures']+=int(fields['error'] is not None)
            if event=='JOB_END':
                t['job_wall_seconds']+=fields['wall_seconds'];t['job_CPU_seconds']+=fields['CPU_seconds']
                a=Counter(t['actual_counters']);a.update(dict(self.counts));t['actual_counters']=dict(a)
            row=dict(sequence=t['sequence'],utc=time.time(),event=event,job=self.job,
                     context=self.context,RHS=rhs,state_attempts=states,LU_attempts=lu,kind=kind,
                     cumulative_RHS=t['RHS'],**fields)
            with (OUT/'RESOURCE_LEDGER.jsonl').open('a') as stream:
                stream.write(json.dumps(plain(row),allow_nan=False,separators=(',',':'))+'\n')
                stream.flush();os.fsync(stream.fileno())
            write(OUT/'RESOURCE_LEDGER_SNAPSHOT.json',t)
    def reserve_state(self,material,*,resolution,purpose):
        if resolution!=24:raise ValueError('original TM24 only')
        self._event('FORWARD_ATTEMPT_RESERVED',rhs=8,states=1,lu=1,kind='forward',purpose=purpose)
        self.credit+=8
    def reserve_rhs(self,count,role):
        count=int(count)
        if role=='full_forward':
            if self.credit<count:raise RuntimeError('missing paid forward reservation')
            self.credit-=count;self.counts['actual_full_forward_RHS']+=count
        else:
            self._event('DERIVATIVE_RHS_RESERVED',rhs=count,kind='derivative',solver_role=role)
            self.counts['actual_derivative_RHS']+=count
    @contextmanager
    def span(self,label,**counts):
        start=time.perf_counter();self.counts.update(counts)
        try:yield
        finally:self.counts[label+'_seconds']+=time.perf_counter()-start
    @contextmanager
    def physical(self,label):
        if self.active is not None:raise RuntimeError('nested physical accounting')
        action=uuid.uuid4().hex;self._event('PHYSICS_BEGIN',action=action,label=label)
        self.active=action;start=time.perf_counter();cpu=time.process_time();err=None
        try:yield
        except BaseException as exc:err=repr(exc);raise
        finally:
            self._event('PHYSICS_END',action=action,label=label,wall_seconds=time.perf_counter()-start,
                        CPU_seconds=time.process_time()-cpu,error=err)
            self.active=None
    def close(self,status='COMPLETE',error=None):
        if self.finished:return
        r=resource.getrusage(resource.RUSAGE_SELF)
        self._event('JOB_END',status=status,error=error,wall_seconds=time.perf_counter()-self.started,
                    CPU_seconds=time.process_time()-self.cpu0,peak_RSS_bytes=int(r.ru_maxrss),
                    unused_prepaid_forward_RHS=self.credit,actual_counters=dict(self.counts))
        self.finished=True;fcntl.flock(self.control,fcntl.LOCK_UN);self.control.close()

class Target(OriginalTarget):
    def __init__(self,book):
        super().__init__(book,keep_fields=False)
    def persist_cache(self,p,path):
        # Lossless final current/field/J, without redundant 576x576 LU copies.
        # Optimizers never reload these caches: persistence cannot alter steps.
        save(path,q=p.q,chi=p.state.chi,current=p.state.current,field=p.state.field,
             prediction=p.state.data,backward_residual=p.state.residual,
             J=np.zeros((0,0)) if p.J is None else p.J,
             receiver_T=np.zeros((0,0),complex) if p.T is None else p.T)

def read_online(path):
    with np.load(path,allow_pickle=False) as z:
        if set(z.files)!={'observation','target_id','sigma'}:raise ValueError('online allowlist')
        if (z['observation'].shape!=(32,8) or not np.isfinite(z['observation']).all()
            or not np.iscomplexobj(z['observation'])
            or str(z['target_id'])!='P1-CORRELATED-HOST-MIXTURE-TM24-v1-20261018'
            or float(z['sigma'])!=0.00119211798514519):
            raise ValueError('original data/noise contract')
        return z['observation'].copy()
