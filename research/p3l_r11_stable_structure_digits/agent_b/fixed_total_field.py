"""Fixed total-field Galerkin solver with full residual acceptance and exact fallback.

The state is e=(I-delta*C*diag(a))^-1 E.  a is a shared real occupancy.
J_occ is complex, with observations ordered source first, receiver second.
No data, material decoder, or truth access is performed by this module.
"""
from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[_name]='1'
import json,time
from pathlib import Path
import numpy as np
from scipy.linalg import svd,eigh,lu_factor,lu_solve

SIGMA=.00119211798514519

class CostLedger:
    """Append-only per-attempt ledger; every solved RHS is reserved in advance."""
    def __init__(self,path=None,rhs_cap=30000):
        self.path=Path(path) if path else None
        self.rhs_cap=int(rhs_cap);self.rhs={};self.lu={};self.events=0
        self.green_columns=0;self.start=time.perf_counter();self.cpu=time.process_time()
        if self.path:
            self.path.parent.mkdir(parents=True,exist_ok=True)
            if self.path.exists(): raise FileExistsError(self.path)
        self.event('BEGIN',rhs_cap=self.rhs_cap,BLAS_threads=1,RHS576=0)
    def event(self,event,**kw):
        self.events+=1
        row=dict(sequence=self.events,event=event,utc=time.time(),**kw)
        if self.path:
            with self.path.open('a') as f: f.write(json.dumps(row,allow_nan=False)+'\n')
    def reserve(self,order,rhs,label,lu=1):
        if sum(self.rhs.values())+rhs>self.rhs_cap: raise RuntimeError('CONSERVATIVE_RHS_CAP')
        self.rhs[int(order)]=self.rhs.get(int(order),0)+int(rhs)
        self.lu[int(order)]=self.lu.get(int(order),0)+int(lu)
        self.event('SOLVE_RESERVED',order=int(order),RHS=int(rhs),LU=int(lu),label=label,
                   RHS_conservative_total=sum(self.rhs.values()))
    def summary(self):
        return dict(RHS_by_order=self.rhs,LU_by_order=self.lu,
                    RHS_conservative_total=sum(self.rhs.values()),RHS576=0,
                    Green_columns_applied=self.green_columns,
                    wall_seconds=time.perf_counter()-self.start,CPU_seconds=time.process_time()-self.cpu)
    def close(self,status='COMPLETE',**kw):
        result=dict(status=status,**self.summary(),**kw);self.event('END',**result);return result


class CallbackLedger(CostLedger):
    """Adapter for shared/core.py callable Ledger, without double-charging."""
    def __init__(self,callback,rhs_cap=300000):
        self.callback=callback
        super().__init__(None,rhs_cap)
    def event(self,event,**kw):
        self.events+=1
        row=dict(event=event,**kw)
        if event=='SOLVE_RESERVED':
            row.update(RHS_by_order={str(kw['order']):kw['RHS']})
        self.callback(row)


def block_krylov(operator,E,max_rank,tol=1e-12,ledger=None):
    """Nested block Arnoldi, fixed injection E; finite rank, no silent padding."""
    t=time.perf_counter()
    U,s,_=svd(E,full_matrices=False,check_finite=False)
    V=U[:,s>tol*s[0]][:,:max_rank];block=V.copy();widths=[V.shape[1]]
    opnorm=float(np.linalg.norm(operator,'fro'))
    while V.shape[1]<max_rank and block.shape[1]:
        if ledger: ledger.green_columns+=block.shape[1]
        W=operator@block
        for _ in range(2): W-=V@(V.conj().T@W)
        U,s,_=svd(W,full_matrices=False,check_finite=False)
        if not len(s) or s[0]<=tol*opnorm: break
        block=U[:,s>tol*s[0]][:,:max_rank-V.shape[1]]
        V=np.column_stack((V,block));widths.append(block.shape[1])
    return V,dict(basis_seconds=time.perf_counter()-t,basis_rank=int(V.shape[1]),block_widths=widths,
                  orthogonality_error=float(np.linalg.norm(V.conj().T@V-np.eye(V.shape[1]))))


class LocalDissipativeBound:
    """Current-shape weighted-state coercivity bounds total-field residual.

    r=E-(I-delta C A)e_tilde; P=sqrt(A), B=PCP.
    ||F-F_tilde||F <= ||OP||2 ||P r||F /(Im(delta)/|delta|^2+lambda_min(Im_H B)).
    Numerical eigenvalues/norms are not outward-rounded interval certificates.
    """
    def __init__(self,C,O,sigma=SIGMA,white_gate=.01,potential_gate=.1):
        self.C=np.asarray(C);self.O=np.asarray(O);self.sigma=float(sigma)
        self.white_gate=white_gate;self.potential_gate=potential_gate;self._a=None;self._meta=None
        self.setup_seconds=0.;self.calls=0
    def shape(self,a):
        a=np.asarray(a,float)
        if self._a is not None and np.array_equal(a,self._a): return self._meta,0.
        t=time.perf_counter();p=np.sqrt(a);B=p[:,None]*self.C*p[None,:]
        beta=float(eigh((B-B.conj().T)/(2j),eigvals_only=True,subset_by_index=(0,0),check_finite=False)[0])
        port=float(np.linalg.norm(self.O*p[None,:],2));elapsed=time.perf_counter()-t
        self._a=a.copy();self._meta=dict(beta=beta,port_norm=port)
        self.setup_seconds+=elapsed
        return self._meta,elapsed
    def __call__(self,a,delta,residual,approxF,y=None):
        meta,setup=self.shape(a);self.calls+=1
        b=float(delta.imag/abs(delta)**2+meta['beta']) if delta else float('-inf')
        rn=float(np.linalg.norm(np.sqrt(a)[:,None]*residual))
        bound=float(np.sqrt(2)/self.sigma*meta['port_norm']*rn/b) if b>0 else None
        pot=None if y is None or bound is None else float(np.sqrt(2)/self.sigma*np.linalg.norm(approxF-y)*bound+.5*bound**2)
        accepted=bool(bound is not None and bound<=self.white_gate and (y is None or pot<=self.potential_gate))
        return dict(accepted=accepted,white_bound=bound,potential_bound=pot,b_delta=b,
                    residual_weighted_norm=rn,bound_setup_seconds=setup,**meta)


class TotalFieldSolver:
    def __init__(self,C,O,E,Fh,a_anchor=None,ranks=(32,64,96),basis_kind='anchor',
                 ledger=None,bound_callback=None,sigma=SIGMA,white_gate=.01,potential_gate=.1,
                 exact_jacobian=True):
        self.C=np.asarray(C,complex);self.O=np.asarray(O,complex);self.E=np.asarray(E,complex);self.Fh=np.asarray(Fh,complex)
        self.n=self.C.shape[0];self.nr=self.O.shape[0];self.ns=self.E.shape[1]
        self.ranks=tuple(sorted(set(int(r) for r in ranks)));self.basis_kind=basis_kind
        self.ledger=(ledger if hasattr(ledger,'reserve') else CallbackLedger(ledger)) if ledger is not None else CostLedger();self.exact_jacobian=exact_jacobian
        self.name='FIXED_TOTAL_FIELD_'+basis_kind.upper()
        if max(self.ranks)>self.n: raise ValueError('rank exceeds full state dimension')
        if basis_kind=='anchor':
            if a_anchor is None: raise ValueError('anchor occupancy is required')
            operator=self.C*np.asarray(a_anchor,float)[None,:]
        elif basis_kind=='green': operator=self.C
        else: raise ValueError('basis_kind must be anchor or green')
        self.V,self.setup=block_krylov(operator,self.E,max(self.ranks),ledger=self.ledger)
        self.setup.update(basis_kind=basis_kind,ranks=list(self.ranks),basis_injection='E_fixed_total_field')
        self.bound=bound_callback or LocalDissipativeBound(self.C,self.O,sigma,white_gate,potential_gate)
        self._a=None;self._models={};self.model_setup_seconds=0.;self.evaluations=0
        self.ledger.event('FIXED_BASIS_BUILT',**self.setup)
    def _models_for(self,a):
        if self._a is not None and np.array_equal(a,self._a): return self._models,0.
        t=time.perf_counter();models={}
        # One full Green action up to the maximum usable fixed rank; all nested ranks share it.
        self.ledger.green_columns+=self.V.shape[1]
        CAV=self.C@(a[:,None]*self.V);OAV=self.O@(a[:,None]*self.V)
        for rank in self.ranks:
            if rank>self.V.shape[1]: continue
            V=self.V[:,:rank];CV=CAV[:,:rank]
            models[rank]=(V,CV,OAV[:,:rank],V.conj().T@CV,V.conj().T@self.E)
        self._a=a.copy();self._models=models;dt=time.perf_counter()-t;self.model_setup_seconds+=dt
        self.ledger.event('CURRENT_SHAPE_PROJECTED',available_ranks=list(models),wall_seconds=dt)
        return models,dt
    def _exact(self,a,delta,need_jac,label):
        rhs=self.ns+(self.nr if need_jac else 0)
        self.ledger.reserve(self.n,rhs,label,lu=1);t=time.perf_counter()
        A=np.eye(self.n,dtype=complex)-delta*(self.C*a[None,:]);fac=lu_factor(A,check_finite=False)
        e=lu_solve(fac,self.E,check_finite=False);F=self.Fh+delta*(self.O@(a[:,None]*e))
        Jo=Jd=None
        if need_jac:
            # W=O A_occ (I-delta C A_occ)^-1, via 32 transpose solves sharing the same LU.
            W=lu_solve(fac,(self.O*a[None,:]).T,trans=1,check_finite=False).T
            Q=self.O+delta*(W@self.C) # O(I-delta A_occ C)^-1
            Jo=(delta*np.einsum('ri,is->sri',Q,e)).reshape(self.ns*self.nr,self.n)
            Jd=(Q@(a[:,None]*e)).T.reshape(-1)
        dt=time.perf_counter()-t
        return F,Jo,Jd,dict(exact_seconds=dt,exact_RHS=rhs,exact_LU=1,exact_full_residual_norm=float(np.linalg.norm(self.E-A@e)))
    def evaluate(self,a,delta,need_jac=False,y=None):
        t0=time.perf_counter();before_rhs=self.ledger.rhs.copy();before_lu=self.ledger.lu.copy();a=np.asarray(a,float).reshape(-1);delta=complex(delta)
        if len(a)!=self.n or np.min(a)<0 or np.max(a)>1 or not np.isfinite(a).all():
            raise ValueError('occupancy must be finite in [0,1], one value per ROI cell')
        if y is not None: y=np.asarray(y,complex)
        self.evaluations+=1;attempts=[];setup_time=0.;F=Jo=Jd=None;fallback=False
        if need_jac:
            # The exact state is already necessary for the controlled full derivative.
            # Return its consistent F as well; no artificial reduced solve is charged.
            F,Jo,Jd,ex=self._exact(a,delta,True,'EXACT_FULL_JACOBIAN')
            diag=dict(method='EXACT_JACOBIAN',fallback=False,jacobian_exact=True,chosen_rank=self.n,
                      reduced_derivative_status='NOT_CONTROLLED_NOT_USED',**ex)
        else:
            models,setup_time=self._models_for(a)
            for rank in self.ranks:
                if rank not in models:
                    attempts.append(dict(rank=rank,status='RANK_UNAVAILABLE'));continue
                V,CV,OV,H,R=models[rank];self.ledger.reserve(rank,self.ns,'TOTAL_FIELD_GALERKIN')
                ts=time.perf_counter();fac=lu_factor(np.eye(rank)-delta*H,check_finite=False)
                z=lu_solve(fac,R,check_finite=False);e=V@z;fa=self.Fh+delta*(OV@z)
                residual=self.E-e+delta*(CV@z);self.last_approx_F=fa;self.last_residual=residual;solve_time=time.perf_counter()-ts
                tb=time.perf_counter();gate=self.bound(a,delta,residual,fa,y);bound_time=time.perf_counter()-tb
                attempt=dict(rank=rank,status='ACCEPTED' if gate['accepted'] else 'REJECTED',
                             solve_seconds=solve_time,bound_seconds=bound_time,full_residual_norm=float(np.linalg.norm(residual)),**gate)
                attempts.append(attempt)
                if gate['accepted']: F=fa;break
            if F is None:
                fallback=True;F,Jo,Jd,ex=self._exact(a,delta,False,'EXACT_FORWARD_FALLBACK')
                diag=dict(method='EXACT_FALLBACK',fallback=True,chosen_rank=self.n,jacobian_exact=None,**ex)
            else: diag=dict(method='TOTAL_FIELD_GALERKIN',fallback=False,chosen_rank=attempts[-1]['rank'],jacobian_exact=None)
            diag.update(attempts=attempts,current_shape_projection_seconds=setup_time)
        diag.update(RHS=sum(self.ledger.rhs.values())-sum(before_rhs.values()),
                    LU=sum(self.ledger.lu.values())-sum(before_lu.values()),
                    RHS_by_order={str(k):v-before_rhs.get(k,0) for k,v in self.ledger.rhs.items() if v-before_rhs.get(k,0)},
                    wall_seconds=time.perf_counter()-t0,evaluation=self.evaluations,basis_kind=self.basis_kind,
                    basis_rank=int(self.V.shape[1]),need_jac=bool(need_jac))
        self.ledger.event('EVALUATE_COMPLETE',**{k:v for k,v in diag.items() if k not in ('RHS','LU','RHS_by_order')})
        return dict(F=F,J_occ=Jo,J_delta=Jd,diagnostics=diag)
