"""One cross-cell residual exchange after diagonal self absorption.

Separate formula-driven development control. Mature scalar seed, Fourier/DPCS,
global PR, stages, material inner solver and original physical validation are
inherited from cie_block_absorption.py. No original source file is changed.
Gin = diag(Gin) + Goff on fixed disjoint four-cell spatial blocks.
R = (I + diag(a) Goff) u, u=beta*(J-diag(m)Q)/(1+(beta-Gii)m),
a=m/(1+(beta-Gii)m). There is no four-cell inverse or small LU.
"""
from collections import Counter, OrderedDict
from copy import deepcopy
import json
from pathlib import Path
import time
import numpy as np
from scipy.linalg import svd
from scipy.optimize import minimize
from research.research_max_struct_r2.src.cie_global import (
    _source_metric, _inner, _exact_current, fourier_apcs, _regularizer,
    _baseline_readout)


class OneHopAbsorption:
    """Diagonal complex divisions and one fixed off-diagonal block action."""
    gin_mode='onehop_block4'
    def __init__(self,G,shape,cost):
        if len(shape)!=2 or any(int(s)%2 for s in shape):
            raise ValueError('Disjoint 2x2 spatial blocks require even 2D shape')
        self.cost=cost
        grid=np.arange(np.prod(shape)).reshape(shape)
        self.indices=np.array([grid[i:i+2,j:j+2].ravel()
            for i in range(0,shape[0],2) for j in range(0,shape[1],2)])
        self.goff=np.asarray(G)[self.indices[:,:,None],self.indices[:,None,:]].copy()
        self.diag_g=np.diag(np.asarray(G)).copy()
        self.goff[:,np.arange(4),np.arange(4)]=0
        self.n=len(self.diag_g);self.nblocks=len(self.indices)
        cost['block_Green_submatrices_extracted']+=self.nblocks
        cost['block_offdiagonal_nonzero_slots']+=12*self.nblocks
        cost['self_Green_diagonal_entries_extracted']+=self.n

    def parameters(self,m,beta):
        b=float(beta)-self.diag_g
        d=1+b*np.asarray(m)
        a=np.asarray(m)/d
        self.cost['onehop_parameter_builds']+=1
        self.cost['scalar_complex_divisions']+=self.n
        self.cost['onehop_a_scalar_divisions']+=self.n
        return b,d,a

    def off_action(self,rhs,adjoint=False):
        matrix=self.goff.conj().transpose(0,2,1) if adjoint else self.goff
        self.cost['onehop_block_matrix_apply_calls']+=1
        self.cost['onehop_block_matrix_apply_blocks']+=self.nblocks
        self.cost['onehop_block_matrix_complex_RHS']+=self.nblocks*rhs.shape[1]
        self.cost['onehop_adjoint_block_RHS' if adjoint else 'onehop_forward_block_RHS']+=self.nblocks*rhs.shape[1]
        out=np.empty_like(rhs,dtype=complex)
        out[self.indices]=matrix@rhs[self.indices]
        return out

    def scaled_rhs(self,rhs,beta,d):
        self.cost['scalar_complex_divisions']+=rhs.size
        self.cost['onehop_u_scalar_divisions']+=rhs.size
        return float(beta)*rhs/d[:,None]

    def apply(self,m,beta,rhs):
        _,d,a=self.parameters(m,beta)
        u=self.scaled_rhs(rhs,beta,d)
        return u+a[:,None]*self.off_action(u)

    def whole_box_bound(self,box,beta):
        """Exact rowwise maximization of m/|1+b*m| on the declared real box.

        d/dm [m^2/|1+b*m|^2] has sign 1+Re(b)*m for m>0;
        check both endpoints plus -1/Re(b) when inside the interval.
        """
        lo,hi=map(float,box)
        b=float(beta)-self.diag_g
        candidate=np.stack((np.full(self.n,lo),np.full(self.n,hi),
            np.divide(-1.,b.real,out=np.full(self.n,np.nan),where=b.real!=0)))
        legal=np.isfinite(candidate)&(candidate>=lo)&(candidate<=hi)
        ratio=np.where(legal,np.nan_to_num(candidate)/np.abs(1+b*np.nan_to_num(candidate)),-np.inf)
        max_ratio=np.max(ratio,axis=0)
        block_rows=np.sum(np.abs(self.goff),axis=2)
        rows=np.empty(self.n);rows[self.indices]=block_rows
        row_bound=max_ratio*rows
        # Minimum denominator magnitude also certifies the diagonal map is pole-free.
        interior=np.divide(-b.real,np.abs(b)**2,out=np.full(self.n,lo),where=np.abs(b)>0)
        closest=np.clip(interior,lo,hi)
        min_d=np.minimum.reduce([np.abs(1+b*lo),np.abs(1+b*hi),np.abs(1+b*closest)])
        bound=float(np.max(row_bound))
        return dict(beta=float(beta),material_box=[lo,hi],bound_inf=bound,
            sufficient_invertibility_condition=bool(bound<1 and np.min(min_d)>0),
            diagonal_denominator_min_abs=float(np.min(min_d)),
            max_ratio=float(np.max(max_ratio)),max_offdiagonal_row_sum=float(np.max(rows)),
            maximizing_row=int(np.argmax(row_bound)),interior_critical_cells=int(np.sum(legal[2])),
            formula='max_i(max_m_in_box m/abs(1+(beta-Gii)*m) * sum_j!=i_in_block abs(Gij))',
            candidate_rule='box endpoints plus m=-1/Re(beta-Gii) when legal; analytic, no GT',
            scope='Bound<1 suffices for invertibility of I+diag(a)Goff; bound>=1 would be inconclusive, not a failure proof')


def onehop_material_value_gradient(blocks,J,Q,m,beta,Os,rho,reg,cost):
    """Exact real material derivative, with fixed J, Q and source metric."""
    b,d,a=blocks.parameters(m,beta)
    u=blocks.scaled_rhs(J-m[:,None]*Q,beta,d)
    goff_u=blocks.off_action(u)
    R=u+a[:,None]*goff_u
    W=R@Os
    adj_u=W+blocks.off_action(a.conj()[:,None]*W,adjoint=True)
    du=-float(beta)*(Q+b[:,None]*J)/(d[:,None]**2)
    da=1/(d**2)
    cost['scalar_complex_divisions']+=du.size+da.size
    cost['onehop_du_scalar_divisions']+=du.size
    cost['onehop_da_scalar_divisions']+=da.size
    prior,gprior=reg(m)
    cost['material_inner_evaluations']+=1
    cost['material_prior_evaluations']+=1
    cost['source_metric_products']+=1
    value=.5*rho*float(np.vdot(R,W).real)+float(prior)
    gradient=rho*np.real(np.sum(adj_u.conj()*du+W.conj()*goff_u*da[:,None],axis=1))+np.asarray(gprior,float)
    return value,gradient


def onehop_material_readout(blocks,J,Q,m,beta,Os,Ls,rho,bounds,reg,config,cost):
    lo,hi=bounds
    # Preserve mature shared-real scalar CIE readout as seed ONLY. Its scalar
    # state objective is never used for acceptance or inner optimization here.
    seed,_=_baseline_readout(J@Ls,Q@Ls,m,beta,np.ones(J.shape[1]),rho,
                             bounds,reg,False,config,cost)
    cache_values=bool(config.get('cache_material_value_gradient',False))
    value_cache=[]
    def fun(z):
        if not cache_values:
            return onehop_material_value_gradient(blocks,J,Q,z,beta,Os,rho,reg,cost)
        cost['material_value_gradient_requests']+=1
        for old_z,value,gradient in value_cache:
            if np.array_equal(old_z,z):
                cost['material_value_gradient_cache_hits']+=1
                return value,gradient.copy()
        cost['material_value_gradient_cache_misses']+=1
        value,gradient=onehop_material_value_gradient(blocks,J,Q,z,beta,Os,rho,reg,cost)
        value_cache.append((np.asarray(z).copy(),value,gradient.copy()))
        if len(value_cache)>2:value_cache.pop(0)
        return value,gradient
    start_evals=cost['material_inner_evaluations']
    old_value=fun(m)[0];seed_value=fun(seed)[0]
    start=seed if seed_value<old_value else m.copy()
    opt=minimize(fun,start,jac=True,method='L-BFGS-B',bounds=list(zip(lo,hi)),
        options=dict(maxiter=int(config.get('material_inner_maxiter',60)),
                     ftol=float(config.get('material_inner_ftol',1e-12)),
                     gtol=float(config.get('material_inner_gtol',1e-8))))
    chosen=np.asarray(opt.x,float)
    final_value=fun(chosen)[0] if np.all(np.isfinite(chosen)) else np.inf
    fallback=not np.isfinite(final_value) or final_value>min(old_value,seed_value)
    if fallback:chosen=start.copy();final_value=min(old_value,seed_value)
    cost['material_inner_LBFGS_iterations']+=int(opt.nit)
    cost['material_inner_LBFGS_evaluations']+=int(opt.nfev)
    cost['material_inner_fallbacks']+=int(fallback)
    return chosen,dict(material_inner_nit=int(opt.nit),material_inner_nfev=int(opt.nfev),
        material_inner_total_evaluations=cost['material_inner_evaluations']-start_evals,
        material_inner_success=bool(opt.success),material_inner_message=str(opt.message),
        material_old_block_loss=old_value,material_scalar_seed_block_loss=seed_value,
        material_final_block_loss=final_value,material_fallback=fallback)


def run_onehop_cie(model,Y,m0,config):
    """Localized extension of R2 run_cie; only global PR CIE is enabled.

    Same stages/betas/box/prior/DPCS and Fourier masks as mature baseline.
    New config data_metric/state_metric accept full fixed Hermitian SPD source
    metrics; legacy diagonal data_weights/state_weights still accepted. All
    source metrics are held fixed. Scalar seed readout uses mature _readout
    after Cholesky source whitening, followed by the exact onehop-loss solve.
    Global PR uses Frobenius gradients/ratios;
    global stop is ||delta A||F^2/max(||A||F^2,tiny). Unit source rotations
    transport Omega to U* Omega U. Do not recompute per-column normalizers.

    Only `sourcewise_cg=False`, `current_solver='pr'`, and `method='cie'` are
    enabled; dormant inherited branches remain solely for a localized diff.
    Optional receipt_path or cost_callback receives execution costs even on
    failure; original exception is re-raised with cie_execution_receipt attached.
    Parent owns strong-corrector budget and candidate handoff.
    """
    cfg=dict(config);counter=Counter();cfg['_counter']=counter
    if cfg.get('method','cie') != 'cie' or cfg.get('sourcewise_cg',False) or cfg.get('current_solver','pr') != 'pr':
        raise ValueError('This extension is restricted to mature global PR CIE')
    wall=time.perf_counter();cpu=time.process_time()
    before=deepcopy(dict(getattr(model,'cost',{})))
    receipt=dict(status='RUNNING',sourcewise_cg=bool(cfg.get('sourcewise_cg',False)),
                 current_solver=cfg.get('current_solver','pr'),FFT=0,
                 physical_backend='unchanged supplied MaxwellModel',
                 residual='u+diag(a)Goff u; u=beta*(J-diag(m)Q)/d; d=1+(beta-Gii)m; a=m/d',
                 gin_mode='onehop_block4', small_solver='elementwise complex division only; no small LU/inverse',
                 block_geometry='fixed disjoint 2x2 spatial-cell actual Green; zero diagonal in Goff',
                 cache_material_value_gradient=bool(cfg.get('cache_material_value_gradient',False)),
                 prior='unchanged supplied regularizer on physical real material',
                 timer_scope='solver preparation/current/readout/physical validation; excludes parent geometry/data and output serialization')
    result=None
    try:
        result=_run_cie(model,Y,m0,cfg);receipt['status']='COMPLETED'
        return result
    except Exception as exc:
        receipt.update(status='FAILED',error=repr(exc))
        exc.cie_execution_receipt=receipt
        raise
    finally:
        after=dict(getattr(model,'cost',{}))
        delta={k:after.get(k,0)-before.get(k,0) for k in set(before)|set(after)
               if isinstance(after.get(k,0),(int,float)) and isinstance(before.get(k,0),(int,float))}
        receipt.update(wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
                       operation_counts=dict(counter),physical_model_delta=delta,
                       physical_LU=delta.get('LU',0),physical_complex_RHS=delta.get('actual_complex_RHS',0))
        if result is not None:
            result['execution_receipt']=receipt
            result['cost']['cpu_seconds']=receipt['cpu_seconds']
        if cfg.get('receipt_path'):
            path=Path(cfg['receipt_path']);path.parent.mkdir(parents=True,exist_ok=True)
            path.write_text(json.dumps(receipt,indent=2,allow_nan=False))
        if cfg.get('cost_callback'):cfg['cost_callback'](receipt)


def _run_cie(model, Y, m0, config):
    """Return dict(m,candidates,history,cost,last_joint_m,status).

    Required config: stages=[dict(mf=..., beta=..., maxiter=...)],
    dpcs_rank OR dpcs_rtol (noise/calibration-frozen). A stage's optional
    beta_current gives the transported two-beta variant. method='tsom'
    uses ordinary LSIE with Fourier APCS; method='csi' uses unrestricted
    current and ordinary LSIE. All modes use the same real material box.

    Optional config: regularizer(m)->(R,grad_R), physical_whiten(residual),
    material_box=(.02,20), rho=1, data_weights/state_weights (positive source
    weights), max_wall_seconds, budget_check(), max_iterations (global),
    coefficient_tol=1e-3 (global or source-wise SQUARED coefficient change),
    miniter=2, sourcewise_cg=False. Physical whiten defaults to division by
    physical_sigma (default 1). Source normalizations default to the paper.
    Every retained candidate is re-evaluated with the original Maxwell model.
    """
    config = dict(config)
    started = time.perf_counter()
    cost = config['_counter']
    model_before = deepcopy(dict(getattr(model, 'cost', {})))
    G, E, S = [np.asarray(x, complex) for x in (model.G, model.E, model.S)]
    Y, m = np.asarray(Y, complex), np.asarray(m0, float).copy()
    n, ns = E.shape
    shape = tuple(model.shape)
    if G.shape != (n, n) or S.shape[1] != n or Y.shape != (S.shape[0], ns):
        raise ValueError('Inconsistent public operators / observations')
    if m.shape != (n,) or not np.all(np.isfinite(m)) or np.prod(shape) != n:
        raise ValueError('Invalid full real material/grid')
    blocks = OneHopAbsorption(G, shape, cost)
    method = config.get('method', 'cie').lower()
    if method not in ('cie', 'tsom', 'csi'):
        raise ValueError('method must be cie, tsom, or csi')
    stages = [dict(s) for s in config['stages']]
    if not stages:
        raise ValueError('At least one preregistered stage required')
    lo, hi = [np.broadcast_to(np.asarray(x, float), (n,)).copy()
              for x in config.get('material_box', (.02, 20.))]
    if np.any(lo > hi) or np.any(m < lo) or np.any(m > hi):
        raise ValueError('Initial material must belong to the common material box')
    reg = _regularizer(config)
    has_reg = config.get('regularizer') is not None
    rho = float(config.get('rho', 1.))
    if rho <= 0:
        raise ValueError('Positive state penalty required')
    Od, Ld = _source_metric(config, 'data_metric', 'data_weights', ns,
                             1/np.maximum(np.sum(abs(Y)**2, axis=0), 1e-300), cost)
    Os, Ls = _source_metric(config, 'state_metric', 'state_weights', ns,
                             1/np.maximum(np.sum(abs(E)**2, axis=0), 1e-300), cost)
    sourcewise = bool(config.get('sourcewise_cg', False))
    if sourcewise and (not np.allclose(Od,np.diag(np.diag(Od)),rtol=0,atol=0)
                       or not np.allclose(Os,np.diag(np.diag(Os)),rtol=0,atol=0)):
        raise ValueError('Historical sourcewise mode requires diagonal metrics; use sourcewise_cg=False')
    wd, ws = np.diag(Od).real, np.diag(Os).real
    whiten = config.get('physical_whiten', lambda a: a/float(config.get('physical_sigma', 1.)))
    current_solver = config.get('current_solver', 'pr')
    if current_solver not in ('pr', 'exact'):
        raise ValueError('current_solver must be pr or exact')
    metric_info = dict(data_eigenvalues=np.linalg.eigvalsh(Od).tolist(),
                       state_eigenvalues=np.linalg.eigvalsh(Os).tolist(),
                       semantics='fixed throughout run; caller transports full metrics under source mixing')
    history, candidates = [], []
    stage_states = []
    def validate(z, stage, reason):
        state = model.state(z)
        pred = state.get('Y', state.get('y'))
        if pred is None:
            pred = S @ state['j']
        residual = np.asarray(pred)-Y
        nll = .5*float(np.linalg.norm(whiten(residual))**2)
        score = nll+float(reg(z)[0])
        rec = dict(m=z.copy(), stage=stage, reason=reason, physical_nll=nll,
                   objective=score, physical_relative_residual=float(np.linalg.norm(residual)/max(np.linalg.norm(Y), 1e-300)))
        candidates.append(rec)
        cost['physical_validations'] += 1
        return rec
    def expired():
        check = config.get('budget_check')
        if check is not None:
            check()
        return time.perf_counter()-started >= float(config.get('max_wall_seconds', np.inf))
    validate(m, -1, 'initial_incumbent')
    if method == 'csi':
        JD, Vd = np.zeros_like(E), np.empty((n, 0), complex)
        receiver_rank = 0
    else:
        Us, sig, Vh = svd(S, full_matrices=False, check_finite=False)
        cost['receiver_SVD'] += 1
        if 'dpcs_rank' in config:
            receiver_rank = int(config['dpcs_rank'])
        elif 'dpcs_rtol' in config:
            receiver_rank = int(np.count_nonzero(sig > float(config['dpcs_rtol'])*sig[0]))
        else:
            raise ValueError('Freeze dpcs_rank or dpcs_rtol from public noise/calibration')
        if not 0 <= receiver_rank <= len(sig):
            raise ValueError('Invalid DPCS rank')
        if receiver_rank and sig[receiver_rank-1] <= np.finfo(float).eps*max(S.shape)*sig[0]:
            raise ValueError('DPCS attempts to invert a numerically null receiver mode')
        Vd = Vh[:receiver_rank].conj().T
        JD = Vd @ ((Us[:, :receiver_rank].conj().T @ Y)/sig[:receiver_rank, None])
    J = JD.copy()
    old_coefficients = {}
    GJD, SJD = G@JD, S@JD
    cost['G_columns'] += ns
    cost['S_columns'] += ns
    iterations, status = 0, 'SCHEDULE_COMPLETE'
    for stage_index, stage in enumerate(stages):
        if expired():
            status = 'WALL_BUDGET'
            break
        beta_m = float(stage.get('beta_material', stage.get('beta', 1.))) if method == 'cie' else None
        beta_j = float(stage.get('beta_current', beta_m)) if method == 'cie' else None
        if method == 'cie' and (beta_m <= 0 or beta_j <= 0 or np.any(1+beta_m*lo <= 0) or np.any(1+beta_j*lo <= 0)):
            raise ValueError('CIE requires a positive beta on one pole-free material branch')
        if method == 'csi':
            U, basis_info = np.eye(n, dtype=complex), dict(mf=None, fourier_columns=n, apcs_rank=n)
        else:
            U, basis_info = fourier_apcs(shape, stage['mf'], Vd, cost)
        if method == 'csi':
            A = J.copy()
        else:
            frequencies = basis_info.pop('frequencies')
            A = np.array([old_coefficients.get(k, np.zeros(ns, complex)) for k in frequencies])
        # Exact zero-padding carries the original Fourier coefficients between
        # nested masks; it preserves both current and the coefficient metric.
        J = JD + U@A
        GU, SU = G@U, S@U
        cost['G_columns'] += U.shape[1]
        cost['S_columns'] += U.shape[1]
        old_grad = old_direction = None
        stage_stop = 'MAXITER'
        for it in range(int(stage.get('maxiter', 500))):
            if expired() or iterations >= int(config.get('max_iterations', 10**12)):
                status = stage_stop = 'WALL_BUDGET' if expired() else 'ITERATION_BUDGET'
                break
            # Transport through physical m: R_betaJ and R_betaM never reused numerically.
            if beta_j is None:
                H = U-m[:, None]*GU
                T = JD-m[:, None]*(E+GJD)
            else:
                packed = np.concatenate((U-m[:,None]*GU, JD-m[:,None]*(E+GJD)),axis=1)
                transformed = blocks.apply(m,beta_j,packed)
                H, T = transformed[:,:U.shape[1]], transformed[:,U.shape[1]:]
                cost['onehop_current_operator_builds'] += 1
            D, R = SJD-Y+SU@A, T+H@A
            cost['current_objective_evaluations'] += 1
            data_before = .5*_inner(D,D,Od,cost)
            state_before = .5*rho*_inner(R,R,Os,cost)
            grad = SU.conj().T@(D@Od) + rho*H.conj().T@(R@Os)
            cost['source_metric_products'] += 2
            cost['current_gradient_evaluations'] += 1
            if current_solver == 'exact':
                target, exact_info = _exact_current(SU,H,SJD-Y,T,Ld,Ls,rho,config,cost)
                direction = change = target-A
                pr, alpha, restarted = 0., 1., False
            else:
                if old_grad is None:
                    direction = -grad
                    pr = 0.
                else:
                    axis = 0 if sourcewise else None
                    pr = np.maximum(0., np.sum(np.real(grad.conj()*(grad-old_grad)), axis=axis) /
                                    np.maximum(np.sum(abs(old_grad)**2, axis=axis), 1e-300))
                    direction = -grad+old_direction*pr
                axis = 0 if sourcewise else None
                gd = np.sum(np.real(grad.conj()*direction), axis=axis)
                restarted = bool(np.any(gd >= 0))
                if axis == 0:
                    bad = gd >= 0
                    direction[:, bad] = -grad[:, bad]
                    cost['non_descent_source_restarts'] += int(np.count_nonzero(bad))
                elif gd >= 0:
                    direction = -grad
                    cost['non_descent_global_restarts'] += 1
                DD, RR = SU@direction, H@direction
                numerator = -np.sum(np.real(grad.conj()*direction), axis=axis)
                if sourcewise:
                    denom = wd*np.sum(abs(DD)**2, axis=0)+rho*ws*np.sum(abs(RR)**2, axis=0)
                else:
                    denom = _inner(DD,DD,Od,cost)+rho*_inner(RR,RR,Os,cost)
                alpha = np.divide(numerator, denom, out=np.zeros_like(numerator), where=denom > 1e-300)
                change = direction*alpha
                exact_info = {}
            if sourcewise:
                delta = float(np.mean(np.sum(abs(change)**2,axis=0)/np.maximum(np.sum(abs(A)**2,axis=0),1e-300)))
            else:
                delta = float(np.vdot(change,change).real/max(np.vdot(A,A).real,1e-300))
            A = A+change
            J, Q = JD+U@A, E+GJD+GU@A
            data_after_current = .5*_inner(SJD-Y+SU@A,SJD-Y+SU@A,Od,cost)
            state_after_current = .5*rho*_inner(T+H@A,T+H@A,Os,cost)
            m_before_readout = m.copy()
            m, inner = onehop_material_readout(blocks,J,Q,m,beta_m,Os,Ls,rho,
                                               (lo,hi),reg,config,cost)
            cost['material_readout_source_whitening_products'] += 2
            cost['source_metric_products'] += 2
            old_grad, old_direction = grad, direction
            data_f = .5*_inner(SJD-Y+SU@A,SJD-Y+SU@A,Od,cost)
            def state_value(beta):
                e = J-m[:, None]*Q
                if beta is not None:
                    e = blocks.apply(m,beta,e)
                return .5*rho*_inner(e,e,Os,cost)
            history.append(dict(stage=stage_index, iteration=it+1, global_iteration=iterations+1,
                                beta_current=beta_j, beta_material=beta_m, receiver_rank=receiver_rank,
                                **basis_info, coefficient_relative_change_squared=delta,
                                data_term=data_f, state_term_current=state_value(beta_j),
                                state_term_material=state_value(beta_m), regularizer=float(reg(m)[0]),
                                wall=time.perf_counter()-started, sourcewise_cg=sourcewise,
                                current_solver=current_solver, PR_plus=np.asarray(pr).tolist(),
                                line_alpha=np.asarray(alpha).tolist(), direction_restarted=restarted,
                                data_before_current=data_before,state_before_current=state_before,
                                data_after_current=data_after_current,state_after_current=state_after_current,
                                current_step_objective_change=data_after_current+state_after_current-data_before-state_before,
                                readout_material_change=float(np.linalg.norm(m-m_before_readout)),
                                **exact_info, **inner))
            if config.get('iteration_callback'):
                config['iteration_callback'](dict(stage=stage_index,iteration=it+1,
                    m=m.copy(),J=J.copy(),coefficients=A.copy()))
            iterations += 1
            cost['PR_CG_steps' if current_solver=='pr' else 'exact_current_steps'] += 1
            cost['reduced_operator_columns'] += 6*ns
            if it+1 >= int(config.get('miniter', 2)) and delta < float(config.get('coefficient_tol', 1e-3)):
                stage_stop = 'COEFFICIENT_TOL'
                break
        if not np.all(np.isfinite(m)):
            raise FloatingPointError('Nonfinite real material readout')
        validate(m, stage_index, stage_stop)
        stage_state = dict(stage=stage_index,m=m.copy(),J=J.copy(),Q=Q.copy(),
                           coefficients=A.copy(),beta_current=beta_j,beta_material=beta_m,
                           receiver_rank=receiver_rank,basis_info=dict(basis_info))
        stage_states.append(stage_state)
        if config.get('stage_callback'):
            config['stage_callback'](stage_state)
        if method != 'csi':
            old_coefficients = {k: a.copy() for k, a in zip(frequencies, A)}
        history.append(dict(stage=stage_index, event='stage_end', status=stage_stop, **basis_info,
                            beta_current=beta_j, beta_material=beta_m, wall=time.perf_counter()-started))
        if status != 'SCHEDULE_COMPLETE':
            break
    best = min(candidates, key=lambda c: c['objective'])
    cost['wall_seconds'] = time.perf_counter()-started
    model_after = dict(getattr(model, 'cost', {}))
    model_delta = {k: model_after.get(k, 0)-model_before.get(k, 0)
                   for k in set(model_before)|set(model_after)
                   if isinstance(model_after.get(k, 0), (int, float)) and isinstance(model_before.get(k, 0), (int, float))}
    return dict(m=best['m'].copy(), candidates=candidates, history=history,stage_states=stage_states,
                cost=dict(cost, model=model_delta), last_joint_m=m.copy(), status=status,
                method=method, receiver_rank=receiver_rank, source_metrics=metric_info,
                sourcewise_cg=sourcewise, current_solver=current_solver)
