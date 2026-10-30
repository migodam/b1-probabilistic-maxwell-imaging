"""Native fixed-reference contrast sources on the actual FEM quadrature.

The public operators accept physical current density j, not weak forcing or
an observation-space chart. Every inverse action is an accounted sparse solve.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from pathlib import Path
import json
import os
import time
import numpy as np
from scipy.sparse import coo_matrix, diags, save_npz
from scipy.optimize import minimize
from scipy.linalg import lu_factor, lu_solve
from .fem_model import physical_lock, Prediction
from .fwi import FWIOptions, regularization, dump, plain


def quadrature_interpolation(ev, active):
    """P rows follow p.material_points' element-major quadrature ordering."""
    p=ev.p; basis=p.material_basis; ne,nqe=p._material_shape
    selected=np.flatnonzero(active); rowmap=np.full(ne*nqe,-1,np.int64)
    rowmap[selected]=np.arange(len(selected)); rows=[];cols=[];values=[]
    for local,phi in enumerate(basis.basis):
        col=np.repeat(ev._to_free[basis.element_dofs[local]],nqe)
        val=np.asarray(phi[0]).reshape(-1)
        valid=active&(col>=0)&(val!=0)
        rows.append(rowmap[valid]);cols.append(col[valid]);values.append(val[valid])
    return coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(cols))),
                      shape=(len(selected),len(p._free_dofs))).tocsr()


class CurrentGreen:
    """One frequency's L_ref, GD and receiver GS; no all-port dense Green."""
    def __init__(self,model,index,absorbed=None,reference='air'):
        self.model=self.fullmodel=model;self.book=model.book;self.index=index
        self.ev=model.evaluators[index];self.p=self.ev.p;ev=self.ev;p=self.p
        inside,cells=model.cell_data[index]
        self.active=np.asarray(ev._inside_wall|inside);self.indices=np.flatnonzero(self.active)
        self.P=quadrature_interpolation(ev,self.active);self.weights=p.material_weights[self.active].copy()
        if np.any(self.weights<=0):raise ValueError('positive material quadrature weights required')
        self.sqrt_weights=np.sqrt(self.weights);self.nq=len(self.indices)
        self.einc=self.inc=model.incident[index][self.active].copy()
        self.wall=ev._inside_wall[self.active];self.target=inside[self.active]
        self.cells=np.where(self.wall,ev._cells[self.active],model.nwall+cells[self.active]).astype(int)
        self.absorbed=np.zeros(self.nq,complex) if absorbed is None else np.asarray(absorbed,complex).copy()
        if self.absorbed.shape!=(self.nq,):raise ValueError('absorbed quadrature contrast shape')
        self.absorbed_mask=np.abs(self.absorbed)>0
        self.reference=reference
        self._physical_cache=None
        L0=p._background_matrix
        matrix=(L0-p.k**2*(self.P.T@diags(self.weights*self.absorbed)@self.P)).tocsc()
        with physical_lock(self.book):
            self.factor,self.matrix,self.factor_bytes=ev.factor(matrix,label='CSI '+reference+' reference LU')
            # Sixteen receiver adjoints, independent of the four illumination RHS.
            adj=self.factor.solve(ev.R.conj().T.toarray(),trans='H',label='CSI '+reference+' GS receiver adjoint16')
            self.GS=p.k**2*np.asarray(self.P@adj).conj().T*self.weights[None,:]
            del adj
            if np.any(self.absorbed):
                b=p.k**2*(self.P.T@(self.weights[:,None]*self.absorbed[:,None]*self.inc))
                self.uA=self.factor.solve(np.asarray(b),label='CSI absorbed illumination4')
                self.inc=self.einc+np.asarray(self.P@self.uA)
                self.dA=np.asarray(ev.R@self.uA).T
            else:
                self.uA=None;self.dA=np.zeros((4,16),complex)
            self.data_bg=self.dA
            self.book.check()

    def chi(self,z):
        nw=self.model.nwall;nt=self.model.ntarget
        real=np.r_[z[:nw]-1,z[2*nw:2*nw+nt]]
        loss=np.r_[z[nw:2*nw],z[2*nw+nt:]]
        return real[self.cells]+1j*loss[self.cells]/self.p.freq-self.absorbed

    def apply_GD(self,j):
        j=np.asarray(j);vector=j.ndim==1
        if vector:j=j[:,None]
        if j.ndim!=2 or j.shape[0]!=self.nq:raise ValueError('GD current rows must equal nq')
        b=self.p.k**2*(self.P.T@(self.weights[:,None]*j))
        with physical_lock(self.book):
            self.ensure_factor()
            u=self.factor.solve(np.asarray(b),label='CSI '+self.reference+' GD current RHS'+str(j.shape[1]))
        result=np.asarray(self.P@u)
        return result[:,0] if vector else result

    def apply_GDH(self,r):
        r=np.asarray(r);vector=r.ndim==1
        if vector:r=r[:,None]
        if r.ndim!=2 or r.shape[0]!=self.nq:raise ValueError('GDH rows must equal nq')
        with physical_lock(self.book):
            self.ensure_factor()
            u=self.factor.solve(np.asarray(self.P.conj().T@r),trans='H',label='CSI '+self.reference+' GDH current RHS'+str(r.shape[1]))
        result=self.p.k**2*self.weights[:,None]*np.asarray(self.P@u)
        return result[:,0] if vector else result

    def apply_GS(self,j):return np.asarray(self.GS@np.asarray(j)).T
    def apply_GSH(self,d):return np.asarray(self.GS.conj().T@np.asarray(d).T)

    def ensure_factor(self):
        """Caller holds the physical lock; rebuild only after reference switching."""
        if self.factor is None:
            self.factor,_,self.factor_bytes=self.ev.factor(self.matrix,label='CSI '+self.reference+' reopened reference LU')

    def physical_state(self,z):
        """Full material FEM checkpoint, with no iterative-current approximation."""
        state=self.physical_state_with_factor(z)
        return {key:state[key] for key in ('data','current','field','diagnostics')}

    def clear_material_cache(self):
        self._physical_cache=None

    def physical_state_with_factor(self,z):
        """One exact-z full material state/LU, retained only within this case.

        This method acquires the physics lock on a cache miss. Call it outside
        any existing physical_lock, then lock only the requested adjoint solve.
        """
        z=self.model.validate(z)
        if self._physical_cache is not None and np.array_equal(z,self._physical_cache['z']):
            self.book.event(self.model.phase,'CACHE_LU_REUSE',frequency=self.p.freq,
                cache_scope='one exact-z full physical material state',additional_full_fe=0,additional_rhs=0)
            return self._physical_cache
        self.clear_material_cache()
        chi=self.chi(z)+self.absorbed;p=self.p
        A=(p._background_matrix-p.k**2*(self.P.T@diags(self.weights*chi)@self.P)).tocsc()
        b=p.k**2*(self.P.T@(self.weights[:,None]*chi[:,None]*self.einc))
        with physical_lock(self.book):
            factor,A,storage=self.ev.factor(A,label='CSI full FEM checkpoint LU')
            u=factor.solve(np.asarray(b),label='CSI full FEM checkpoint primal4')
            defect=factor.last_defect
            E=self.einc+np.asarray(self.P@u)
            d=self.ev.direct+np.asarray(self.ev.R@u).T
            self.book.check()
        self._physical_cache=dict(z=z.copy(),data=d,current=chi[:,None]*E,field=E,
                    u=u,factor=factor,A=A,matrix=A,
                    diagnostics=dict(frequency=p.freq,factor_bytes=storage,primal_backward=defect))
        return self._physical_cache


class FixedReference:
    def __init__(self,model):
        self.model=model;self.greens=[CurrentGreen(model,i) for i in range(len(model.evaluators))]
        self.original_modes=None;self.mode_setup_records=[]
        self.public_archive_path=None;self.mode_archive_path=None

    def physical_state(self,z):
        states=[g.physical_state(z) for g in self.greens]
        return np.asarray([s['data'] for s in states]),[s['current'] for s in states],[s['field'] for s in states],[s['diagnostics'] for s in states]

    def gradient(self,z,y):
        """Actual full material gradient, reusing cached full FEM LUs at exact z."""
        model=self.model;z=model.validate(z);y=np.asarray(y)
        if y.shape!=(len(self.greens),4,16) or not np.isfinite(y).all():raise ValueError('full weighted reference gradient observations')
        data=[];diagnostics=[];gradient=np.zeros(model.nparam);nw=model.nwall;nt=model.ntarget
        for g,observed in zip(self.greens,y):
            state=g.physical_state_with_factor(z);d=state['data'];factor=state['factor']
            with physical_lock(model.book):
                rhs=np.asarray(g.ev.R.conj().T@((d-observed)/model.sigma**2).T)
                adjoint=factor.solve(rhs,trans='H',label='CSI cached full material weighted adjoint4')
                av=np.asarray(g.P@adjoint)
                products=g.p.k**2*g.weights*np.sum(np.conj(av)*state['field'],axis=1)
                cell_gradient=np.zeros(nw+nt,complex);np.add.at(cell_gradient,g.cells,products)
                real=cell_gradient.real;loss=(1j/g.p.freq*cell_gradient).real
                gradient+=np.r_[real[:nw],loss[:nw],real[nw:],loss[nw:]]
                diagnostic=dict(state['diagnostics'],adjoint_backward=factor.last_defect,
                    gradient_operator='actual P3 products with reused exact-z full material LU')
            data.append(d.copy());diagnostics.append(diagnostic)
        prediction=Prediction(np.asarray(data),diagnostics)
        value=.5*float(np.sum(np.abs(prediction.data-y)**2))/model.sigma**2
        return value,gradient,prediction

    def save_public(self,path):
        """Archive data-independent fixed operators once; LU is rebuilt as metered."""
        path=Path(path);path.mkdir(parents=True,exist_ok=True);self.public_archive_path=str(path.resolve())
        for ii,g in enumerate(self.greens):
            save_npz(path/f'P_f{ii}.npz',g.P);save_npz(path/f'L0_f{ii}.npz',g.matrix)
            save_npz(path/f'R_f{ii}.npz',g.ev.R)
            np.savez_compressed(path/f'operators_f{ii}.npz',GS=g.GS,weights=g.weights,
                incident=g.einc,points=g.p.material_points[g.active],cells=g.cells,wall=g.wall,target=g.target,
                frequency=g.p.freq,material_indices=g.indices)
        dump(path/'PUBLIC_REFERENCE.json',dict(config={k:v for k,v in self.model.config.items() if k!='case_id'},
            data_dependent_arrays=False,truth_access=False,reference='original air/PML L0, actual quadrature P/W, GS0',
            archive=self.public_archive_path,scope='public operator reuse is available equally to ordinary B5 and absorption methods'))

    def prepare_modes(self,rank,path):
        """One paid public GD SVD per frequency, reusable within one batch only."""
        if self.original_modes is not None:
            if any(m['Q'].shape[1]!=rank for m in self.original_modes):raise ValueError('shared GD mode rank mismatch')
            return self.original_modes
        from .spectral import build_modes
        path=Path(path);path.mkdir(parents=True,exist_ok=True);self.mode_archive_path=str(path.resolve())
        self.original_modes=[]
        for ii,g in enumerate(self.greens):
            before=self.model.book.snapshot(pid=os.getpid());modes=build_modes(g,rank=rank)
            kept={key:modes[key] for key in ('Q','CQ','singular_values')};self.original_modes.append(kept)
            after=self.model.book.snapshot(pid=os.getpid())
            row=dict(kind='ORIGINAL_GD_MODE_SETUP',frequency=g.p.freq,rank=rank,
                setup_rhs=after['rhs']-before['rhs'],setup_full_fe=after['full_fe']-before['full_fe'],
                cpu_seconds=modes['cpu_seconds'],wall_seconds=modes['elapsed_seconds'],
                right_modes_fixed=True,shared_across_independent_cold_methods=False)
            self.mode_setup_records.append(row)
            np.savez_compressed(path/f'GD_modes_f{ii}.npz',**kept,
                quadrature_points=g.p.material_points[g.active],quadrature_weights=g.weights,frequency=g.p.freq)
        dump(path/'GD_MODE_SETUP.json',dict(mode_setups=self.mode_setup_records,rank=rank))
        return self.original_modes


@dataclass
class CSIOptions:
    iterations:int=12
    starts:int=2
    checkpoint_every:int=4
    state_weight:float=10000.
    material_maxiter:int=25
    absorption_at:int=4
    reopen_at:int=8
    absorption_fraction:float=.5
    absorption_validation_rtol:float=1e-8
    risk_stability_window:int=2
    risk_max_material_step:float=.05
    risk_max_conditional_sd:float=.25
    risk_max_projected_gradient:float=.1
    risk_max_data_impact:float=1.
    risk_reopen_projected_gradient:float=.1
    risk_reopen_objective_growth:float=1.
    current_initialization:str='observable_residual'
    fixedtik_relative:float=1e-4
    spectral_preconditioner:bool=False
    spectral_rank:int=32
    strong_rank:int=8
    spectral_condition_limit:float=1e10
    target_risk_screen:bool=True
    risk_target_scaled_RMS:float=.1
    physical_reset_at:int=0

    def __post_init__(self):
        if min(self.iterations,self.starts,self.checkpoint_every,self.material_maxiter)<1:raise ValueError('positive CSI counts required')
        if self.state_weight<=0 or not 0<self.absorption_fraction<=1:raise ValueError('CSI fixed weights/fraction')
        if not 0<self.absorption_at<self.reopen_at<=self.iterations:raise ValueError('absorption/reopen order')
        if self.risk_stability_window<1 or min(self.risk_max_material_step,self.risk_max_conditional_sd,self.risk_max_projected_gradient,self.risk_max_data_impact,self.risk_reopen_projected_gradient,self.risk_reopen_objective_growth)<=0:raise ValueError('positive frozen empirical risk thresholds')
        if self.current_initialization not in ('observable_residual','physical'):raise ValueError('unknown current initialization')
        if self.fixedtik_relative<=0:raise ValueError('positive fixed relative Tikhonov weight required')
        if not 0<self.strong_rank<self.spectral_rank or self.spectral_condition_limit<=1:raise ValueError('spectral strong/weak ranks and condition limit')
        if self.risk_target_scaled_RMS<=0:raise ValueError('positive local target-risk scaled RMS threshold required')
        if not 0<=self.physical_reset_at<=self.iterations:raise ValueError('physical reset iteration must be zero or within the frozen run')


class ModalSchurPreconditioner:
    """Classical modal SPD preconditioning; the full current space remains active."""
    def __init__(self,g,Q,CQ,strong_rank=8,condition_limit=1e10):
        self.g=g;self.Q=np.asarray(Q);self.CQ=np.asarray(CQ)
        self.strong_rank=strong_rank;self.condition_limit=condition_limit
        if self.Q.shape!=self.CQ.shape or self.Q.shape[0]!=g.nq:raise ValueError('fixed GD mode dimensions')

    def apply(self,z,r):
        cpu=time.process_time();wall=time.perf_counter();Q=self.Q;n=Q.shape[1];h=self.strong_rank
        chi=self.g.chi(z)
        M=np.eye(n,dtype=complex)-Q.conj().T@(chi[:,None]*self.CQ)
        row=dict(kind='MODAL_SCHUR_PRECONDITIONER',frequency=self.g.p.freq,reference=self.g.reference,
            rank=n,strong_rank=h,weak_rank=n-h,condition_Mq=None,
            condition_limit=self.condition_limit,dense_LU_factorizations=0,dense_condition_SVDs=0,
            dense_LU_count_scope='attempts counted before calls, including failed attempts',
            current_space_dimension=self.g.nq,current_dimension_dropped=0,
            outside_Q_directions_active=True,ordinary_projected_Schur_equivalence=True,
            scope='classical full-current SPD modal preconditioner; no ROM or new-theorem claim')
        result=r
        HH,HL,LH,LL=M[:h,:h],M[:h,h:],M[h:,:h],M[h:,h:]
        def require_condition(matrix,name):
            row['dense_condition_SVDs']+=1
            value=float(np.linalg.cond(matrix));row[name]=value
            if not np.isfinite(value) or value>self.condition_limit:
                raise ValueError(name+' exceeds frozen condition threshold')
        try:
                require_condition(M,'condition_Mq')
                require_condition(LL,'weak_block_condition')
                row['dense_LU_factorizations']+=1;low=lu_factor(LL)
                ilh=lu_solve(low,LH)
                S=HH-HL@ilh
                require_condition(S,'strong_Schur_condition')
                row['dense_LU_factorizations']+=1;strong=lu_factor(S)
                eye=np.eye(n,dtype=complex)
                ilow=lu_solve(low,eye[h:])
                upper=lu_solve(strong,eye[:h]-HL@ilow)
                lower=ilow-ilh@upper
                inverse=np.vstack((upper,lower))
                projected=Q.conj().T@r
                # D=I+Q(Mq^-1 Mq^-H-I)QH is SPD in the full current space.
                result=r-Q@projected+Q@(inverse@(inverse.conj().T@projected))
                positive=float(np.real(np.vdot(r,result)))
                row.update(identity_fallback=False,
                    strong_to_weak_coupling_norm=float(np.linalg.norm(LH)),
                    weak_to_strong_coupling_norm=float(np.linalg.norm(HL)),
                    exact_weak_feedback_norm=float(np.linalg.norm(HL@ilh)),
                    inverse_relative_defect=float(np.linalg.norm(M@inverse-eye)/np.linalg.norm(eye)),
                    weighted_gradient_D_energy=positive)
                if not np.isfinite(result).all() or positive<0:
                    result=r;row.update(identity_fallback=True,fallback_reason='nonfinite or nonpositive numerical SPD action')
        except (ValueError,np.linalg.LinAlgError) as error:
            result=r;row.update(identity_fallback=True,fallback_reason=type(error).__name__+': '+str(error))
        row.update(dense_process_cpu_seconds=time.process_time()-cpu,dense_wall_seconds=time.perf_counter()-wall)
        row={k:(None if isinstance(v,(float,np.floating)) and not np.isfinite(v) else v) for k,v in row.items()}
        self.g.book.event(self.g.model.phase,'ordinary modal Schur SPD current preconditioner',
                          event_kind=row['kind'],**{k:v for k,v in row.items() if k!='kind'})
        return result,row


def observable_residual_initialize(g,j_physical,y,fixedtik_relative=1e-4):
    """Standard GS observable-space seed; no claim of a material tangent or TSOM.

    The 16-dimensional Gram is (GS W^-1/2)(GS W^-1/2)^H. Its cached
    eigendecomposition supplies the weighted GS singular values and a fixed
    relative Tikhonov residual correction; no domain Green matrix is formed.
    """
    g.book.check();cpu=time.process_time();wall=time.perf_counter()
    cache=getattr(g,'_observable_seed_gram',None);new_eigh=cache is None
    if cache is None:
        gram=(g.GS/g.weights[None,:])@g.GS.conj().T
        gram=(gram+gram.conj().T)*.5
        eigenvalues,vectors=np.linalg.eigh(gram)
        eigenvalues=np.maximum(eigenvalues,0.)
        if eigenvalues[-1]<=0:raise RuntimeError('zero GS observable Gram')
        cache=(eigenvalues,vectors);g._observable_seed_gram=cache
    eigenvalues,vectors=cache;regularizer=fixedtik_relative*float(eigenvalues[-1])
    residual=np.asarray(y)-g.ev.direct-g.data_bg-g.apply_GS(j_physical)
    coefficients=(vectors.conj().T@residual.T)/(eigenvalues[:,None]+regularizer)
    receiver_solution=vectors@coefficients
    correction=g.apply_GSH(receiver_solution.T)/g.weights[:,None]
    dense_cpu=time.process_time()-cpu;dense_wall=time.perf_counter()-wall
    current=np.asarray(j_physical)+correction
    field=g.inc+g.apply_GD(current)
    after=np.asarray(y)-g.ev.direct-g.data_bg-g.apply_GS(current)
    row=dict(kind='CURRENT_OBSERVABLE_INITIALIZATION',frequency=g.p.freq,
             reference=g.reference,initialization='observable_residual',
             weighted_GS_singular_values=np.sqrt(eigenvalues[::-1]).tolist(),
             Gram_dimension=int(len(eigenvalues)),dense_eigendecompositions=int(new_eigh),
             dense_eigendecomposition_size=16,dense_solve_rhs=int(residual.shape[0]),
             decomposition_reused=not new_eigh,fixedtik_relative=fixedtik_relative,
             dense_Gram_projection_cpu_seconds=dense_cpu,dense_Gram_projection_wall_seconds=dense_wall,
             tikhonov_lambda=regularizer,data_residual_norm_before=float(np.linalg.norm(residual)),
             data_residual_norm_after=float(np.linalg.norm(after)),
             quadrature_correction_norm=float(np.sqrt(np.sum(g.weights[:,None]*np.abs(correction)**2))),
             cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
             scope='standard GS observable residual initialization; no material-direction legitimacy or full TSOM claim')
    g.book.event(g.model.phase,'standard GS observable current seed',event_kind=row['kind'],
                 **{key:value for key,value in row.items() if key!='kind'})
    return current,field,row


def state_terms(g,z,j,field,alpha):
    residual=j-g.chi(z)[:,None]*field
    value=.5*alpha*float(np.sum(g.weights[:,None]*np.abs(residual)**2))
    prod=-alpha*g.weights*np.sum(np.conj(residual)*field,axis=1)
    nc=g.model.nwall+g.model.ntarget;re=np.zeros(nc);loss=np.zeros(nc)
    np.add.at(re,g.cells,prod.real);np.add.at(loss,g.cells,(1j/g.p.freq*prod).real)
    nw=g.model.nwall
    gradient=np.r_[re[:nw],loss[:nw],re[nw:],loss[nw:]]
    return value,gradient,residual


def current_step(g,z,j,field,y,alpha,sigma):
    chi=g.chi(z);_,_,r=state_terms(g,z,j,field,alpha)
    residual=g.ev.direct+g.dA+g.apply_GS(j)-y
    # Inverse-W preconditioning gives descent in the physical L2 quadrature norm.
    grad=g.apply_GSH(residual)/sigma**2+alpha*(g.weights[:,None]*r-g.apply_GDH(np.conj(chi)[:,None]*g.weights[:,None]*r))
    weighted_gradient=grad/g.sqrt_weights[:,None]
    weighted_gradient[g.absorbed_mask]=0
    preconditioner=getattr(g,'current_preconditioner',None);preconditioner_stats=None
    if preconditioner is not None:
        weighted_direction,preconditioner_stats=preconditioner.apply(z,weighted_gradient)
    else:weighted_direction=weighted_gradient
    weighted_direction[g.absorbed_mask]=0
    direction=-weighted_direction/g.sqrt_weights[:,None]
    # Absorbed cells have no residual current until their material is reopened.
    direction[g.absorbed_mask]=0
    gd=g.apply_GD(direction);gs=g.apply_GS(direction)
    state_direction=direction-chi[:,None]*gd
    numerator=-float(np.real(np.vdot(grad,direction)))
    denominator=float(np.sum(np.abs(gs)**2))/sigma**2+alpha*float(np.sum(g.weights[:,None]*np.abs(state_direction)**2))
    step=0. if denominator<=0 else max(0.,numerator/denominator)
    return j+step*direction,field+step*gd,dict(step=step,directional_derivative=-numerator,
        quadratic_curvature=denominator,modal_Schur=preconditioner_stats)


def invert(model,y,out,options=None,prior_options=None,absorption=False,risk_controlled=False,shared_reference=None):
    o=CSIOptions() if options is None else (CSIOptions(**options) if isinstance(options,dict) else options)
    prior=FWIOptions() if prior_options is None else (FWIOptions(**prior_options) if isinstance(prior_options,dict) else prior_options)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if (out/'ONLINE_SEALED.json').exists():raise RuntimeError('refusing sealed CSI overwrite')
    y=np.asarray(y);nf=len(model.evaluators)
    if y.shape!=(nf,4,16) or not np.isfinite(y).all():raise ValueError('observations [frequency,source,receiver]')
    lo,hi,scale=model.bounds_and_scale();started=time.time();best=None;last=None;fatal=None
    checkpoints=0;checkpoint_requests=0;checkpoint_cache_hits=0;physical_cache=None
    completed=0;absorption_events=[];physical_reset_events=[];ref=None
    shared_boundary_restored=None if shared_reference is None else False
    original_modes=[];mode_setups=[];preconditioner_records=[]
    absorption=absorption or risk_controlled
    method='CSI_EMPIRICAL_RISK_ABSORPTION' if risk_controlled else ('CSI_EXACT_ABSORPTION' if absorption else 'CSI_FIXED_REFERENCE')
    if o.spectral_preconditioner:method+='_MODAL_SCHUR_SPD'
    if o.spectral_preconditioner and not absorption:method+='_NO_ASSIMILATION'
    if risk_controlled:
        baseline_class='B4_EMPIRICAL_RISK_MODAL_SCHUR' if o.spectral_preconditioner else 'B4_EMPIRICAL_RISK'
    elif absorption:
        baseline_class='B3_ABSORPTION_MODAL_SCHUR' if o.spectral_preconditioner else 'B2_HARD_ABSORB_FORCED_REOPEN'
    else:
        baseline_class='B5_ORDINARY_MODAL_SCHUR_NO_ASSIMILATION' if o.spectral_preconditioner else 'B1_NATIVE_CSI'
    reset_control=not absorption and o.physical_reset_at>0
    if reset_control:
        method+='_PHYSICAL_RESET_AT_'+str(o.physical_reset_at)
        baseline_class+='_PHYSICAL_RESET_CONTROL'
    dump(out/'ONLINE_INPUT_CONTRACT.json',dict(method=method,baseline_class=baseline_class,
        truth_access=False,digit_id_access=False,shape_dictionary=False,options=asdict(o),prior=asdict(prior),
        config=model.config,reference_current='j(q,s)=chi(q)*(einc+GD0 j)',
        reference='fixed air/PML original GS0 retained',state_metric='fixed state_weight * quadrature L2',
        data_metric='independent real/imag sigma .01',free_material_parameters=model.nparam,
        checkpoint_policy='evaluate/select material only; currents retained except explicit pre-absorb reset or named B1 physical-reset control; exact-z one-point physical cache available to every CSI method',
        physical_reset_control=reset_control,physical_reset_at=o.physical_reset_at if reset_control else None,
        attribution_limit='absorption gains require paired B1 physical-reset control; pre-absorb reset can itself change optimization',
        optional_target_risk_numerical_failure='KEEP_ACTIVE and continue exact CSI; preserve paid costs and failure; budget/resource/system/program exceptions still stop',
        current_initialization_scope='standard fixed-GS observable residual seed, or physical-current ablation; no material-direction or full TSOM claim',
        spectral_scope='standard full-current SPD modal Schur preconditioning; ordinary projected-Schur B5 is semantically equivalent; no omitted-space truncation' if o.spectral_preconditioner else None,
        public_reference_reuse=shared_reference is not None,
        public_reference_path=getattr(shared_reference,'public_archive_path',None),
        data_dependent_carry_between_cases=False,
        absorption_policy=('B4 preliminary empirical screen then optional paid local full-J target-material tangent screen; conflict-triggered REOPEN; no posterior, nonlinear, or active-box safety certificate' if risk_controlled else 'B2 hard baseline: estimated contrast ranking and scheduled forced REOPEN; not empirical-risk B4') if absorption else None))
    def trace(row):
        with (out/'TRACE.jsonl').open('a') as f:f.write(json.dumps(dict(time_unix=time.time(),**row),default=plain,allow_nan=False)+'\n')
    def install_absorbed_modes(greens,start,iteration):
        if not o.spectral_preconditioner:return
        folder=out/'preconditioner';folder.mkdir(exist_ok=True)
        for ii,g in enumerate(greens):
            cpu=time.process_time();wall=time.perf_counter();before=model.book.snapshot(pid=os.getpid())
            Q=original_modes[ii]['Q']
            # Same fixed original right modes; only the reference GD action changes.
            CQ=g.sqrt_weights[:,None]*g.apply_GD(Q/g.sqrt_weights[:,None])
            g.current_preconditioner=ModalSchurPreconditioner(g,Q,CQ,o.strong_rank,o.spectral_condition_limit)
            after=model.book.snapshot(pid=os.getpid())
            row=dict(kind='MODAL_REFERENCE_REFRESH',start=start,iteration=iteration,frequency=g.p.freq,
                reference=g.reference,rank=o.spectral_rank,original_right_modes_reused=True,
                actual_GD_setup_rhs=after['rhs']-before['rhs'],setup_full_fe=after['full_fe']-before['full_fe'],
                cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
            mode_setups.append(row);trace(row)
            np.savez_compressed(folder/f'absorbed_start{start}_iter{iteration}_f{ii}.npz',CQ=CQ,
                original_mode_file=str(Path(ref.mode_archive_path)/f'GD_modes_f{ii}.npz'),frequency=g.p.freq)
    def checkpoint(z,start,iteration,reason):
        nonlocal best,last,checkpoints,checkpoint_requests,checkpoint_cache_hits,physical_cache
        checkpoint_requests+=1;cache_cpu=time.process_time();cache_wall=time.perf_counter()
        cache_hit=physical_cache is not None and np.array_equal(z,physical_cache['z'])
        if cache_hit:
            data=physical_cache['data'].copy();js=[a.copy() for a in physical_cache['js']]
            fields=[a.copy() for a in physical_cache['fields']];diag=physical_cache['diagnostics']
            checkpoint_cache_hits+=1
        else:
            data,js,fields,diag=ref.physical_state(z);checkpoints+=1
            physical_cache=dict(z=z.copy(),data=data.copy(),js=tuple(a.copy() for a in js),
                fields=tuple(a.copy() for a in fields),diagnostics=diag,
                source_start=start,source_iteration=iteration,source_physical_checkpoint=checkpoints)
        cache_copy_cpu=time.process_time()-cache_cpu;cache_copy_wall=time.perf_counter()-cache_wall
        data_value=.5*float(np.sum(np.abs(data-y)**2))/model.sigma**2
        pv,_,parts=regularization(z,model.nwall,model.target_shape,model.wall_shape,prior)
        state=dict(z=z.copy(),data=data,objective=data_value+pv,data_objective=data_value)
        last=state;io_cpu=time.process_time();io_wall=time.perf_counter()
        np.savez_compressed(out/'endpoint_last.npz',**state,**model.maps(z))
        if best is None or state['objective']<best['objective']:
            best=state;np.savez_compressed(out/'endpoint.npz',**state,**model.maps(z))
            raw=dict(z=z.copy(),data=data,frequencies=np.asarray([g.p.freq for g in ref.greens]),
                     selected_start=start,selected_iteration=iteration,
                     physical_scope='full FEM material checkpoint; original physical j=chi*E',
                     field_scope='total analytic incident plus full FEM scattered field at actual material quadrature')
            for ii,(g,j,field) in enumerate(zip(ref.greens,js,fields)):
                raw[f'j_f{ii}']=j;raw[f'total_field_f{ii}']=field
                raw[f'quadrature_points_f{ii}']=g.p.material_points[g.active]
                raw[f'quadrature_weights_f{ii}']=g.weights
            np.savez_compressed(out/'best_currents.npz',**raw)
        io_cpu=time.process_time()-io_cpu;io_wall=time.perf_counter()-io_wall
        if cache_hit:
            trace(dict(kind='CACHE_REUSE',start=start,iteration=iteration,reason=reason,
                cache='one-point exact material z full FEM state',new_physics_rhs=0,new_physics_full_fe=0,
                source_start=physical_cache['source_start'],source_iteration=physical_cache['source_iteration'],
                source_physical_checkpoint=physical_cache['source_physical_checkpoint'],
                cache_copy_cpu_seconds=cache_copy_cpu,cache_copy_wall_seconds=cache_copy_wall,
                checkpoint_io_cpu_seconds=io_cpu,checkpoint_io_wall_seconds=io_wall))
        trace(dict(kind='FULL_FEM_CHECKPOINT',start=start,iteration=iteration,reason=reason,
                   cache_reuse=cache_hit,physical_checkpoint_solve_count=checkpoints,
                   checkpoint_io_cpu_seconds=io_cpu,checkpoint_io_wall_seconds=io_wall,
                   objective=state['objective'],data_objective=data_value,prior_components=parts,diagnostics=diag))
        return js,fields
    try:
        ref=FixedReference(model) if shared_reference is None else shared_reference
        if ref.model is not model:raise ValueError('shared public reference must belong to the exact same FEM model/book')
        if o.spectral_preconditioner:
            already_public=ref.original_modes is not None
            original_modes=ref.prepare_modes(o.spectral_rank,out/'preconditioner')
            if already_public:
                mode_setups.append(dict(kind='PUBLIC_GD_MODE_REUSE',setup_rhs=0,setup_full_fe=0,
                    original_mode_archive=ref.mode_archive_path,rank=o.spectral_rank,data_dependent_carry=False))
            else:mode_setups.extend(ref.mode_setup_records)
            for row in mode_setups:trace(row)
            for g,modes in zip(ref.greens,original_modes):
                g.current_preconditioner=ModalSchurPreconditioner(g,modes['Q'],modes['CQ'],o.strong_rank,o.spectral_condition_limit)
            if shared_reference is not None:
                dump(out/'PUBLIC_REFERENCE_ROUTE.json',dict(public_reference_path=ref.public_archive_path,
                    mode_archive_path=ref.mode_archive_path,original_modes_reused=True,
                    data_dependent_carry=False,weights='fixed original physical quadrature only'))
        else:
            for g in ref.greens:g.current_preconditioner=None
        from .assimilation import absorb_wall, reopen, risk_probe, risk_selection
        for start in range(o.starts):
            z=model.initial(start);js,fields=checkpoint(z,start,0,'initial full material')
            greens=ref.greens;frozen=np.zeros(model.nparam,bool)
            if o.current_initialization=='observable_residual':
                seed_rows=[]
                for ii,g in enumerate(greens):
                    js[ii],fields[ii],row=observable_residual_initialize(g,js[ii],y[ii],o.fixedtik_relative)
                    row.update(start=start);seed_rows.append(row);trace(row)
                dump(out/'GS_OBSERVABLE_SINGULARS.json',dict(reference='fixed original air/PML GS0',
                    weighting='singular values of GS0 W^-1/2',frequencies=[dict(frequency=row['frequency'],
                    singular_values=row['weighted_GS_singular_values']) for row in seed_rows],
                    standard_initialization=True,material_parameters_remain_free=True))
            else:
                trace(dict(kind='CURRENT_PHYSICAL_INITIALIZATION',start=start,
                           initialization='physical',additional_current_seed_RHS=0))
            history=[z.copy()];anchor_objective=None;is_absorbed=False
            for iteration in range(1,o.iterations+1):
                model.book.check()
                step_records=[]
                for ii,g in enumerate(greens):
                    js[ii],fields[ii],row=current_step(g,z,js[ii],fields[ii],y[ii],o.state_weight,model.sigma)
                    step_records.append(row)
                    if row['modal_Schur'] is not None:preconditioner_records.append(dict(start=start,iteration=iteration,**row['modal_Schur']))
                def material_fun(q):
                    zz=np.clip(q*scale,lo,hi);v,gg,_=regularization(zz,model.nwall,model.target_shape,model.wall_shape,prior)
                    for g,j,field in zip(greens,js,fields):
                        sv,sg,_=state_terms(g,zz,j,field,o.state_weight);v+=sv;gg+=sg
                    return v,gg*scale
                bounds=list(zip(lo/scale,hi/scale))
                for k in np.flatnonzero(frozen):bounds[k]=(z[k]/scale[k],z[k]/scale[k])
                result=minimize(material_fun,z/scale,method='L-BFGS-B',jac=True,bounds=bounds,
                                options=dict(maxiter=o.material_maxiter,maxls=10,ftol=1e-9,gtol=1e-5))
                z=model.validate(result.x*scale)
                history.append(z.copy())
                state_value=sum(state_terms(g,z,j,f,o.state_weight)[0] for g,j,f in zip(greens,js,fields))
                data_value=.5*sum(float(np.sum(np.abs(g.ev.direct+g.dA+g.apply_GS(j)-yy)**2)) for g,j,yy in zip(greens,js,y))/model.sigma**2
                trace(dict(kind='CSI_ITERATION',start=start,iteration=iteration,state_objective=state_value,
                    current_data_objective=data_value,material_optimizer_success=bool(result.success),
                    material_optimizer_stop=str(result.message),material_iterations=int(result.nit),current_steps=step_records,
                    reference=greens[0].reference,free_material_parameters=int((~frozen).sum())))
                completed+=1
                if reset_control and iteration==o.physical_reset_at:
                    before=model.book.snapshot(pid=os.getpid())
                    js,fields=checkpoint(z,start,iteration,'B1 physical reset ablation; no absorption')
                    after=model.book.snapshot(pid=os.getpid())
                    event=dict(kind='PHYSICAL_RESET_CONTROL',start=start,iteration=iteration,
                        absorption=False,all_material_parameters_active=True,
                        free_material_parameters=model.nparam,actual_rhs=after['rhs']-before['rhs'],
                        actual_full_fe=after['full_fe']-before['full_fe'])
                    physical_reset_events.append(event);trace(event)
                if absorption and not risk_controlled and iteration==o.absorption_at:
                    # A full physical state makes the transformation validation independent of CSI feasibility.
                    js,fields=checkpoint(z,start,iteration,'before exact absorption')
                    greens,js,fields,frozen,event=absorb_wall(ref,z,js,fields,o.absorption_fraction,o.absorption_validation_rtol)
                    install_absorbed_modes(greens,start,iteration)
                    event.update(start=start,iteration=iteration);absorption_events.append(event);trace(event)
                if absorption and not risk_controlled and iteration==o.reopen_at:
                    js,fields,event=reopen(ref,greens,js,fields,z,o.absorption_validation_rtol)
                    greens=ref.greens;frozen[:]=False
                    event.update(start=start,iteration=iteration);absorption_events.append(event);trace(event)
                scheduled=iteration%o.checkpoint_every==0 or iteration==o.iterations
                if risk_controlled and scheduled:
                    physical_js,physical_fields=checkpoint(z,start,iteration,'empirical risk checkpoint')
                    probe=risk_probe(model,z,y,prior,reference=ref)
                    if is_absorbed:
                        pg=float(np.max(np.abs(probe['scaled_projected_gradient'][frozen])))
                        growth=float(probe['objective']-anchor_objective)
                        conflict=pg>o.risk_reopen_projected_gradient or growth>o.risk_reopen_objective_growth
                        event=dict(kind='RISK_REOPEN_DECISION',start=start,iteration=iteration,
                            max_absorbed_scaled_projected_gradient=pg,physical_objective_growth=growth,
                            actual_full_FEM_objective=probe['objective'],diagnostics=probe['diagnostics'],
                            decision='REOPEN' if conflict else 'KEEP_ABSORBED')
                        absorption_events.append(event);trace(event)
                        if conflict:
                            js,fields,event=reopen(ref,greens,js,fields,z,o.absorption_validation_rtol)
                            greens=ref.greens;frozen[:]=False;is_absorbed=False
                            event.update(start=start,iteration=iteration,trigger='actual projected-gradient or physical-objective conflict')
                            absorption_events.append(event);trace(event)
                    elif iteration>=o.absorption_at and iteration<o.iterations:
                        selected,event=risk_selection(ref,z,physical_fields,history,probe,prior,o)
                        event.update(start=start,iteration=iteration,actual_full_FEM_objective=probe['objective'],
                                     diagnostics=probe['diagnostics'])
                        absorption_events.append(event);trace(event)
                        if o.target_risk_screen and len(selected):
                            from .risk import target_bias_screen
                            before=model.book.snapshot(pid=os.getpid());cpu=time.process_time();wall=time.perf_counter()
                            try:
                                selected,target_event=target_bias_screen(ref,z,prior,selected,
                                    np.asarray(event['conditional_sd_real']),np.asarray(event['conditional_sd_loss']),
                                    max_scaled_rms=o.risk_target_scaled_RMS)
                            except (ValueError,np.linalg.LinAlgError,FloatingPointError) as error:
                                after=model.book.snapshot(pid=os.getpid());selected=np.empty(0,dtype=int)
                                target_event=dict(kind='TARGET_RISK_SCREEN_FAILED',start=start,iteration=iteration,
                                    decision='KEEP_ACTIVE',optional_risk_failure_fallback=True,
                                    paid_full_J_attempted=True,additional_full_J_computed=False,
                                    error=dict(type=type(error).__name__,message=str(error)),
                                    target_safety_claim=False,
                                    scope='optional local tangent screening incomplete; no target safety claim')
                                model.book.failure(model.phase,'optional local target-risk numerical fallback; KEEP_ACTIVE',error,
                                    paid_rhs=after['rhs']-before['rhs'],paid_full_fe=after['full_fe']-before['full_fe'])
                            except BaseException as error:
                                after=model.book.snapshot(pid=os.getpid())
                                failed=dict(kind='TARGET_RISK_SCREEN_FAILED',start=start,iteration=iteration,
                                    paid_rhs=after['rhs']-before['rhs'],paid_full_fe=after['full_fe']-before['full_fe'],
                                    paid_process_cpu_seconds=time.process_time()-cpu,paid_wall_seconds=time.perf_counter()-wall,
                                    error=dict(type=type(error).__name__,message=str(error)),target_safety_claim=False)
                                absorption_events.append(failed);trace(failed)
                                raise
                            after=model.book.snapshot(pid=os.getpid())
                            target_event.update(start=start,iteration=iteration,
                                paid_rhs=after['rhs']-before['rhs'],paid_full_fe=after['full_fe']-before['full_fe'],
                                paid_process_cpu_seconds=time.process_time()-cpu,paid_wall_seconds=time.perf_counter()-wall,
                                active_box_safe=False,posterior_certificate=False)
                        else:
                            target_event=dict(kind='TARGET_RISK_SCREEN',start=start,iteration=iteration,
                                decision='KEEP_ACTIVE' if not len(selected) else 'TARGET_RISK_DISABLED',
                                reason='no preliminary eligible cells' if not len(selected) else 'target_risk_screen disabled by frozen option',
                                additional_full_J_computed=False,paid_rhs=0,paid_full_fe=0,
                                active_box_safe=False,posterior_certificate=False,
                                scope='TARGET_RISK_NOT_COMPUTED; no target safety claim')
                        absorption_events.append(target_event);trace(target_event)
                        if len(selected):
                            greens,js,fields,frozen,event=absorb_wall(ref,z,physical_js,physical_fields,
                                o.absorption_fraction,o.absorption_validation_rtol,selected=selected)
                            install_absorbed_modes(greens,start,iteration)
                            anchor_objective=probe['objective'];is_absorbed=True
                            event.update(start=start,iteration=iteration,
                                selection='frozen empirical screen and paid local target-material tangent screen' if o.target_risk_screen else 'frozen empirical screen only; target risk not computed',
                                target_risk_computed=bool(target_event.get('additional_full_J_computed')),
                                target_risk_scope='LOCAL_UNCONSTRAINED_GAUSSIAN_TANGENT_WITH_CONDITIONAL_SD_PROXY; not nonlinear, posterior, or active-box safe')
                            absorption_events.append(event);trace(event)
                elif scheduled:
                    # Checkpoints select physical material, never an infeasible current fit.
                    checkpoint(z,start,iteration,'fixed checkpoint schedule')
            dump(out/f'START_{start}.json',dict(start=start,completed_iterations=o.iterations))
    except BaseException as error:
        fatal=dict(type=type(error).__name__,message=str(error));model.book.failure(model.phase,'CSI interrupted',error)
        if isinstance(error,(KeyboardInterrupt,SystemExit)):raise
    finally:
        if shared_reference is not None and ref is not None:
            # Absorbed wall backgrounds and their currents never cross cases.
            for g in locals().get('greens',[]):
                g.clear_material_cache()
                if g not in ref.greens:g.factor=None
            for g in ref.greens:g.clear_material_cache()
            try:
                with physical_lock(model.book):
                    for g in ref.greens:g.ensure_factor()
                shared_boundary_restored=True
                trace(dict(kind='SHARED_CASE_BOUNDARY',original_air_reference_restored=True,
                    previous_material_current_history_retained=False))
            except BaseException as error:
                fatal=dict(type=type(error).__name__,message=str(error))
                model.book.failure(model.phase,'restore shared original air reference',error)
        result=dict(status='COMPLETE_BOUNDED_ATTEMPT' if fatal is None else 'PARTIAL_INTERRUPTED',
            method=method,case_id=model.config.get('case_id'),truth_access=False,template_used=False,
            baseline_class=baseline_class,
            completed_iterations=completed,full_fem_checkpoints=checkpoints,
            checkpoint_requests=checkpoint_requests,checkpoint_cache_hits=checkpoint_cache_hits,
            physical_reset_events=physical_reset_events,physical_reset_control=reset_control,
            error=fatal,absorption_events=absorption_events,
            objective=None if best is None else best['objective'],full_objective_comparable=best is not None,
            selection='minimum full FEM measured objective plus identical frozen FWI TV/spatial prior',
            selected_physical_currents_file='best_currents.npz' if best is not None else None,
            shared_public_reference=shared_reference is not None,
            shared_original_air_reference_restored=shared_boundary_restored,
            wall_seconds=time.time()-started,convergence_claim='bounded attempt; optimizer stops retained',
            original_GS_reference_fixed=True,material_variables_reopened=any(e.get('kind')=='REOPEN' for e in absorption_events))
        result.update(target_risk_screen_enabled=bool(risk_controlled and o.target_risk_screen),
            paid_target_risk_computed=any(e.get('additional_full_J_computed',False) for e in absorption_events),
            optional_risk_failure_fallback=any(e.get('optional_risk_failure_fallback',False) for e in absorption_events),
            optional_risk_numerical_fallback_count=sum(bool(e.get('optional_risk_failure_fallback')) for e in absorption_events),
            target_safety_claim=False,
            target_risk_scope='LOCAL_UNCONSTRAINED_GAUSSIAN_TANGENT_WITH_CONDITIONAL_SD_PROXY; never a posterior, nonlinear, or active-box safety certificate')
        if o.spectral_preconditioner:
            dump(out/'GD_PRECONDITIONER_STATS.json',dict(method='ordinary modal exact weak-Schur SPD preconditioner',
                full_current_space_active=True,dimension_dropped=0,B5_ordinary_projected_Schur_semantically_equivalent=True,
                mode_setups=mode_setups,actions=preconditioner_records,
                dense_LU_factorizations=sum(r['dense_LU_factorizations'] for r in preconditioner_records),
                dense_condition_SVDs=sum(r['dense_condition_SVDs'] for r in preconditioner_records),
                dense_process_cpu_seconds=sum(r['dense_process_cpu_seconds'] for r in preconditioner_records),
                identity_fallback_count=sum(bool(r.get('identity_fallback')) for r in preconditioner_records)))
        dump(out/'ONLINE_SEALED.json',result)
    if last is None:raise RuntimeError('no full FEM checkpoint: '+str(fatal))
    return result
