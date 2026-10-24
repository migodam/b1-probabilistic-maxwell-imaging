"""Conditional exact-arithmetic decision bounds; no automatic physics solves.

All arrays/callbacks are supplied by the caller. Dense inverse acquisition must
use an explicitly charged callback. Norms are Frobenius upper bounds, never a
random norm estimate. Floating arithmetic here is NOT interval certification.
"""
from dataclasses import dataclass,field
from collections import Counter
import time
import numpy as np
from scipy.sparse.linalg import LinearOperator,gmres


def charged_inverse(callback,n):
    """callback(I)->(inverse, actual_ledger_RHS_increment); no hidden inversion."""
    if not callable(callback) or n<1:raise ValueError('paid callback and dimension required')
    value,charged=callback(np.eye(n,dtype=np.complex128))
    if int(charged)!=n:raise AssertionError('every inverse column must be charged exactly once')
    value=np.asarray(value)
    if value.shape!=(n,n) or not np.isfinite(value).all():raise ValueError('invalid paid inverse')
    return value


def decision_from_phi_interval(phi_lower,phi_upper,base_phi,logu):
    if not np.isfinite(base_phi) or base_phi<0 or not np.isfinite(logu) or logu>0:raise ValueError('known finite basePhi and same logU<=0 required')
    if not 0<=phi_lower<=phi_upper:return 'UNKNOWN'
    # Full66 pCN is reversible for original Gaussianized prior, so only Phi
    # enters this one-step ratio. The same U is reused at all refinements.
    low=min(0.,base_phi-phi_upper);high=min(0.,base_phi-phi_lower)
    if logu<low:return 'ACCEPT'
    if logu>=high:return 'REJECT'
    return 'UNKNOWN'


@dataclass
class DecisionRefresh:
    chi0:np.ndarray
    j0:np.ndarray
    E0:np.ndarray
    G:np.ndarray
    S:np.ndarray
    incident:np.ndarray
    R0:np.ndarray
    sigma:float
    paid_inverse_columns:int
    base_equations_checked:bool
    counts:Counter=field(default_factory=Counter)

    def __post_init__(self):
        n=len(self.chi0)
        if self.paid_inverse_columns!=n or not self.base_equations_checked:raise AssertionError('paid complete inverse and checked exact-anchor equations required')
        for key in ('chi0','j0','E0','G','S','incident','R0'):
            setattr(self,key,np.asarray(getattr(self,key)))
            if not np.isfinite(getattr(self,key)).all():raise ValueError('finite arrays required')
        if self.G.shape!=(n,n) or self.R0.shape!=(n,n) or self.j0.shape!=self.incident.shape or self.E0.shape!=self.j0.shape or self.S.shape[1]!=n:raise ValueError('linear-system shapes')
        if not np.isfinite(self.sigma) or self.sigma<=0:raise ValueError('fixed positive noise sigma required')
        self.R0F=float(np.linalg.norm(self.R0));self.SF=float(np.linalg.norm(self.S));self.sqrt2_sigma=np.sqrt(2)/self.sigma
        L0=np.eye(n)-self.chi0[:,None]*self.G
        self.left_inverse_defect_F=float(np.linalg.norm(np.eye(n)-self.R0@L0))
        self.inverse_residual=float(np.linalg.norm(L0@self.R0-np.eye(n)))
        self.anchor_residual=float(np.linalg.norm(self.chi0[:,None]*self.E0-self.j0))
        if self.anchor_residual>1e-8*max(1.,np.linalg.norm(self.j0)):raise ValueError('supplied inverse/anchor numerical audit failed; not interval certification')

    def bound_constants(self,chi):
        chi=np.asarray(chi)
        if chi.shape!=self.chi0.shape or not np.isfinite(chi).all():raise ValueError('candidate material shape')
        self.counts['bound_constant_evaluations']+=1
        self.counts['bound_norm_matrix_entries']+=self.G.size
        eta=self.left_inverse_defect_F+self.R0F*float(np.linalg.norm((chi-self.chi0)[:,None]*self.G))
        return dict(left_inverse_defect_F=self.left_inverse_defect_F,eta_upper=eta,inverse_upper=self.R0F/(1-eta) if eta<1 else None)

    def inspect(self,chi,current,observation,constants):
        start=time.perf_counter();current=np.asarray(current);observation=np.asarray(observation)
        if current.shape!=self.j0.shape or observation.shape!=(self.S.shape[0],self.j0.shape[1]) or not np.isfinite(current).all() or not np.isfinite(observation).all():raise ValueError('finite current and observation shapes')
        sources=current.shape[1]
        E=self.incident+self.G@current;self.counts['G_matvec_source_columns']+=sources
        rho=chi[:,None]*E-current
        prediction=self.S@current;self.counts['S_matvec_source_columns']+=sources
        rnorm=self.sqrt2_sigma*float(np.linalg.norm(prediction-observation));rhof=float(np.linalg.norm(rho))
        inverse=constants['inverse_upper'];e=None if inverse is None else self.sqrt2_sigma*self.SF*inverse*rhof
        lo=0. if e is None else .5*max(0.,rnorm-e)**2
        hi=float('inf') if e is None else .5*(rnorm+e)**2
        self.counts['inspection_calls']+=1
        return dict(current=current.copy(),field=E,prediction=prediction,residual=rho,
                    rho_F=rhof,residual_white_norm=rnorm,error_white_upper=e,phi_lower=lo,phi_upper=hi,
                    approximate_phi=.5*rnorm*rnorm,eta_upper=constants['eta_upper'],inverse_upper=inverse,
                    bound_scope='conditional exact arithmetic; floating point NOT interval certified',inspection_seconds=time.perf_counter()-start)

    def richardson_decision(self,chi,observation,base_phi,logu,max_refinements=4):
        if not 0<=max_refinements<=4:raise ValueError('frozen maximum4 refinements')
        start=time.perf_counter();before=self.counts.copy();dchi=np.asarray(chi)-self.chi0
        current=self.j0+self.R0@(dchi[:,None]*self.E0)
        self.counts['R0_matvec_source_columns']+=current.shape[1]
        constants=self.bound_constants(chi);history=[];approximate_points=[]
        for iteration in range(max_refinements+1):
            point=self.inspect(chi,current,observation,constants)
            decision=decision_from_phi_interval(point['phi_lower'],point['phi_upper'],base_phi,logu)
            approximate_points.append({k:point[k].copy() for k in ('current','field','prediction','residual')})
            history.append({k:v for k,v in point.items() if k not in ('current','field','prediction','residual')})
            history[-1].update(iteration=iteration,decision=decision)
            if decision!='UNKNOWN' or iteration==max_refinements or constants['inverse_upper'] is None:break
            current=current+self.R0@point['residual'];self.counts['R0_matvec_source_columns']+=current.shape[1]
        return dict(point=point,decision=decision,requires_paid_full_fallback=decision=='UNKNOWN',refinements=iteration,
                    same_logu=logu,history=history,approximate_points=approximate_points,wall_seconds=time.perf_counter()-start,counts=dict(self.counts-before))

    def ordinary_krylov(self,chi,observation,base_phi,logu,rtol=1e-10,restart=8,paid_preconditioner=None):
        """Ordinary GMRES control, at most8 inner steps per source.

        Optional anchored-LU preconditioner must charge each actual solve RHS.

        Iterative actions incur matvec/time/memory costs, not free physical
        full-solver calls. Interval uses separately charged R0 setup; standalone
        Krylov accuracy comparison also reports that setup dependency.
        """
        start=time.perf_counter();before=self.counts.copy();n=len(chi);iterations=[];infos=[]
        def action(v):
            self.counts['krylov_G_matvec_columns']+=1
            return v-chi*(self.G@v)
        op=LinearOperator((n,n),matvec=action,dtype=np.complex128);columns=[]
        preconditioner=None
        if paid_preconditioner is not None:
            if not callable(paid_preconditioner):raise ValueError('explicit charged preconditioner callback')
            def apply_preconditioner(v):
                value,charged=paid_preconditioner(v)
                if int(charged)!=1:raise AssertionError('each GMRES LU application must charge1RHS')
                self.counts['krylov_preconditioner_LU_RHS']+=1
                value=np.asarray(value)
                if value.shape!=(n,) or not np.isfinite(value).all():raise ValueError('preconditioner result')
                return value
            preconditioner=LinearOperator((n,n),matvec=apply_preconditioner,dtype=np.complex128)
        for source in range(self.incident.shape[1]):
            count=[0]
            def callback(value):count[0]+=1
            rhs=chi*self.incident[:,source]
            current,info=gmres(op,rhs,x0=self.j0[:,source],rtol=rtol,atol=0.,restart=restart,maxiter=1,callback=callback,callback_type='pr_norm',M=preconditioner)
            columns.append(current);infos.append(int(info));iterations.append(count[0])
        point=self.inspect(chi,np.column_stack(columns),observation,self.bound_constants(chi))
        decision=decision_from_phi_interval(point['phi_lower'],point['phi_upper'],base_phi,logu)
        return dict(point=point,decision=decision,requires_paid_full_fallback=decision=='UNKNOWN',iterations_per_source=iterations,
                    gmres_info=infos,preconditioner='paid anchored LU' if paid_preconditioner is not None else 'none',rtol=rtol,restart=restart,max_restart_cycles=1,wall_seconds=time.perf_counter()-start,counts=dict(self.counts-before),
                    inverse_bound_setup='same paid completeR0; not free standalone Krylov certification')
