"""Exact original TM24 fixed-host ROI physics, with mandatory prepayment.

No GT or digit label is accepted by the online interface. Each paid state is
saved losslessly. A state holds its own L/R/J to avoid duplicate Full-J costs.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[_k]='1'
from pathlib import Path
import importlib.util,json,time
import numpy as np
from scipy.linalg import lu_factor,lu_solve
from scipy.special import ndtr
from .budget import Budget,ROOT,atomic_json

REPO=ROOT.parents[1];B1=REPO.parent
R11=B1/'b1-p3l-r11-stable-digits/research/p3l_r11_stable_structure_digits'
R9=B1/'b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace'
SIGMA=.00119211798514519;HOST=.1+.04j

def load_setup():
    p=ROOT/'inputs/FIXED_HOST_SETUP.npz'
    if not p.exists():p=R11/'inputs/FIXED_HOST_SETUP.npz'
    return dict(np.load(p,allow_pickle=False))

def prior():
    source=ROOT/'inputs/PRIOR_SOURCE.py'
    if not source.exists():source=R11/'inputs/PRIOR_SOURCE.py'
    spec=importlib.util.spec_from_file_location('r12_frozen_prior',source);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    cfg=ROOT/'inputs/PRIOR_CONFIG.json';pts=ROOT/'inputs/PRIOR_POINTS.npz'
    if not cfg.exists():cfg=R11/'inputs/PRIOR_CONFIG.json';pts=R11/'inputs/PRIOR_POINTS.npz'
    return module.MaterialPrior(json.loads(cfg.read_text()),np.load(pts)['points'])

def pack(F):return np.concatenate((F.real.T,F.imag.T),axis=1).ravel()
def white(F):return pack(F)*np.sqrt(2)/SIGMA
def white_jac(J):
    q=J.reshape(8,32,-1)
    return np.concatenate((q.real,q.imag),axis=1).reshape(512,-1)*np.sqrt(2)/SIGMA
def material_rms(a,b):return float(np.linalg.norm(a-b)/24)

class Problem:
    def __init__(self,kind='P1',setup=None):
        self.kind=kind;self.setup=setup or load_setup();self.idx=self.setup['ROI'].astype(int);self.pr=prior()
        self.dimension=66 if kind=='P1' else 258
        self.edge_matrix=None
        if kind!='P1':
            lookup={int(k):i for i,k in enumerate(self.idx)};edges=[]
            for k,i in lookup.items():
                for kk in (k+1 if k%24<23 else -1,k+24):
                    if kk in lookup:edges.append((i,lookup[kk]))
            E=np.zeros((len(edges),256))
            for r,(i,j) in enumerate(edges):E[r,i]=1;E[r,j]=-1
            self.edge_matrix=E
    def regularization(self,x):
        if self.kind=='P1':return .5*float(x@x),x.copy(),np.eye(66)
        a=x[:256];v=self.edge_matrix@a;w=np.sqrt(v*v+.05**2)
        value=.4*np.sum(w-.05)+.01*(a@a)+.05*a.sum()+.5*(x[-2:]@x[-2:])
        grad=np.r_[.4*self.edge_matrix.T@(v/w)+.02*a+.05,x[-2:]]
        H=np.zeros((258,258));H[:256,:256]=.4*(self.edge_matrix.T/w)@self.edge_matrix+.02*np.eye(256);H[-2:,-2:]=np.eye(2)
        return float(value),grad,H
    def project(self,x):
        x=np.asarray(x,float).copy()
        if self.kind!='P1':x[:256]=np.clip(x[:256],0,1)
        return x
    def decode(self,x,derivative=False):
        x=self.project(x);u=ndtr(x[-2:]);du=np.exp(-.5*x[-2:]**2)/np.sqrt(2*np.pi)
        delta=complex(.8+u[0],.03+.15*u[1])
        if self.kind=='P1':
            q=x.copy();q[-2:]=u;mat,aux=self.pr.decode(q);a=aux['fraction'][self.idx]
            D=self.pr.derivative(q)[self.idx]*np.r_[np.ones(64),du][None,:] if derivative else None
        else:
            a=x[:256];mat=np.full(576,HOST,dtype=complex);mat[self.idx]+=delta*a
            D=np.column_stack((delta*np.eye(256),a*du[0],1j*.15*a*du[1])) if derivative else None
        return a,delta,mat,D

class Physics:
    def __init__(self,problem,y,group,out=None,local_rhs_cap=None,local_lu_cap=None):
        self.problem=problem;self.y=np.asarray(y,complex);self.group=group;self.out=Path(out or ROOT/'results'/group/'states')
        self.out.mkdir(parents=True,exist_ok=True);self.budget=Budget();self.calls=0
        s=problem.setup;self.C=s['C'];self.O=s['O'];self.E=s.get('ER',s.get('E'));self.Fh=s['Fh'];self._last=None
        self.rhs=0;self.lu=0;self.local_rhs_cap=local_rhs_cap;self.local_lu_cap=local_lu_cap
        # Explicit restart continues immutable state numbering and paid costs;
        # never overwrite earlier attempts or obtain a fresh local allowance.
        if self.budget.ledger.exists():
            for line in self.budget.ledger.read_text().splitlines():
                row=json.loads(line)
                if row.get('event')=='RESERVE' and row.get('group')==group:
                    self.rhs+=row['request']['RHS'];self.lu+=row['request']['LU']
                    label=row.get('meta',{}).get('state')
                    if label is not None and str(label).isdigit():self.calls=max(self.calls,int(label))
        saved=sorted(self.out.glob('[0-9][0-9][0-9][0-9][0-9].npz'))
        if saved:
            self.calls=max(self.calls,int(saved[-1].stem))
            q=dict(np.load(saved[-1],allow_pickle=False));q.update(prior_legal=True,group=group)
            if problem.kind!='P1':
                rv,rg,rH=problem.regularization(q['x']);q.update(reg_value=rv,reg_gradient=rg,reg_hessian=rH)
            self._last=q
    def evaluate(self,x,need_jac=False):
        x=self.problem.project(x)
        # Exact state cache only: no approximate state reused as a new observation.
        if self._last is not None and np.array_equal(x,self._last['x']) and (not need_jac or 'jac' in self._last):return self._last
        a,delta,material,Dx=self.problem.decode(x,need_jac);d=delta*a
        amount=40 if need_jac else 8
        if self.local_rhs_cap is not None and self.rhs+amount>self.local_rhs_cap:
            from .budget import BudgetStop
            raise BudgetStop('local physical group cap, before operation')
        if self.local_lu_cap is not None and self.lu+1>self.local_lu_cap:
            from .budget import BudgetStop
            raise BudgetStop('local all-order LU group cap, before operation')
        extra_J=need_jac
        self.calls+=1;label=f'{self.calls:05d}'
        with self.budget.operation(self.group,LU=1,RHS=40 if need_jac else 8,order=256,J=need_jac,state=label) as ticket:
            self.rhs+=amount;self.lu+=1
            M=np.eye(256,dtype=complex)-d[:,None]*self.C
            fac=lu_factor(M,check_finite=False);h=lu_solve(fac,d[:,None]*self.E,check_finite=False)
            F=self.Fh+self.O@h;r=white(F-self.y)
            state=dict(x=x.copy(),material=material,occupancy=a,delta=delta,F=F,residual=r,phi=float(.5*r@r),prior_legal=True,ticket=ticket,group=self.group)
            if self.problem.kind!='P1':
                rv,rg,rH=self.problem.regularization(x)
                state.update(reg_value=rv,reg_gradient=rg,reg_hessian=rH)
            if need_jac:
                R=self.E+self.C@h;L=lu_solve(fac,self.O.T,trans=1,check_finite=False).T
                Jm=np.einsum('ri,is->sri',L,R,optimize=True).reshape(256,256)
                Jc=Jm@Dx
                state.update(jac=white_jac(Jc),dmaterial=Dx,L=L,R=R,jac_complex=Jc)
        save={k:v for k,v in state.items() if isinstance(v,np.ndarray) or k in ('phi','delta','reg_value')}
        np.savez_compressed(self.out/f'{label}.npz',**save)
        atomic_json(self.out/f'{label}.json',dict(ticket=ticket,phi=state['phi'],need_jac=need_jac,GT_used=False,utc=time.time()))
        self._last=state
        return state
