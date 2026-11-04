"""Complete matrix-free CIE/Fourier-TSOM integration for frozen R6 arms.

The alternating algorithm is the R5 factorial/R2 global loop. The changed
computational wiring replaces explicit U/GU/SU/H by native FFT actions and
adjoints; historical solvers remain untouched. S and physical validation stay
dense and unchanged. Native actions omit only the historical roundoff SVD.
"""
from collections import Counter
from copy import deepcopy
from pathlib import Path
import json
import time
import numpy as np
from scipy.linalg import svd
from research.research_max_struct_r2.src.cie_global import _source_metric, _inner, _regularizer, _baseline_readout
from research.research_max_struct_r4.src.cie_onehop import OneHopAbsorption, onehop_material_readout
from research.research_max_struct_r6.src.fft_tsom_ops import NativeFourierAPCS, DiskPulseToeplitzGreen


def run_cie_fft(model,Y,m0,config):
    """Compatible API plus fft_geometry={coords,radius,k} and arm.

    arm='cie_dual' uses the mature shared-real scalar material readout;
    arm='a_material' uses the unchanged onehop material readout. Caller supplies
    exactly the frozen stages/betas/prior/bounds/normalizations and corrector.
    sourcewise_cg retains the legacy default True for cie_dual and False for A.
    No truth is read, and candidates are selected by original physical Phi.
    """
    cfg=dict(config); arm=cfg.get('arm','cie_dual')
    if arm not in ('cie_single','cie_dual','global_dual','a_material','a_plus'): raise ValueError('Unsupported frozen FFT arm')
    cfg['current_transform']='scalar';cfg['material_transform']='onehop' if arm in ('a_material','a_plus') else 'scalar'
    cfg.setdefault('sourcewise_cg',arm.startswith('cie_'))
    cfg['_counter']=Counter()
    wall,cpu=time.perf_counter(),time.process_time();result=None
    before=deepcopy(dict(getattr(model,'cost',{})))
    receipt=dict(status='RUNNING',arm=arm,sourcewise_cg=bool(cfg['sourcewise_cg']),
        current_solver='pr',auxiliary_Green='zero-padded finite disk-pulse Toeplitz FFT',
        Fourier='native coefficients and low-rank DPCS projection; no APCS SVD',
        physical_backend='unchanged supplied MaxwellModel',
        prior='unchanged supplied physical-material prior',
        numerical_equivalence='roundoff-perturbed action; no long-trajectory bitwise guarantee',
        timer_scope='solver setup/current/readout/physical validation; excludes parent assembly and output serialization')
    try:
        result=_run_cie(model,Y,m0,cfg);receipt['status']='COMPLETED'
        return result
    except Exception as exc:
        receipt.update(status='FAILED',error=repr(exc));exc.cie_execution_receipt=receipt
        raise
    finally:
        after=dict(getattr(model,'cost',{}))
        delta={k:after.get(k,0)-before.get(k,0) for k in set(before)|set(after)
            if isinstance(after.get(k,0),(int,float)) and isinstance(before.get(k,0),(int,float))}
        counts=dict(cfg['_counter'])
        fft_calls=(counts.get('Green_kernel_FFT_setup',0)+counts.get('coefficient_FFT_calls',0)+
            counts.get('coefficient_IFFT_calls',0)+2*(counts.get('Green_forward_calls',0)+counts.get('Green_adjoint_calls',0)))
        receipt.update(wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
            FFT=fft_calls,operation_counts=counts,physical_model_delta=delta,
            physical_LU=delta.get('LU',0),physical_complex_RHS=delta.get('actual_complex_RHS',0))
        if result is not None:
            result['execution_receipt']=receipt;result['cost']['cpu_seconds']=receipt['cpu_seconds']
        if cfg.get('receipt_path'):
            path=Path(cfg['receipt_path']);path.parent.mkdir(parents=True,exist_ok=True)
            path.write_text(json.dumps(receipt,indent=2,allow_nan=False))
        if cfg.get('cost_callback'): cfg['cost_callback'](receipt)


run_cie=run_cie_fft


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
    geometry=config['fft_geometry']
    green=DiskPulseToeplitzGreen(geometry['coords'],shape,float(geometry['radius']),float(geometry['k']),cost)
    blocks = OneHopAbsorption(G, shape, cost) if config['material_transform']=='onehop' else None
    method = config.get('method', 'cie').lower()
    if method != 'cie':
        raise ValueError('This FFT integration supports the frozen CIE arms only')
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
    if current_solver != 'pr':
        raise ValueError('FFT path preserves PR; no dense exact-current diagnostic')
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
    GJD, SJD = green.apply(JD), S@JD
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
        U=NativeFourierAPCS(shape,stage['mf'],Vd,cost)
        basis_info=dict(mf=int(stage['mf']),fourier_columns=U.k,apcs_rank=None,
            native_action_without_roundoff_SVD=True)
        frequencies=[tuple(map(int,f)) for f in U.frequencies]
        A=np.array([old_coefficients.get(k,np.zeros(ns,complex)) for k in frequencies])
        cost['peak_retained_FFT_spectrum_bytes']=max(cost['peak_retained_FFT_spectrum_bytes'],green.spectrum.nbytes)
        cost['peak_native_coefficient_bytes']=max(cost['peak_native_coefficient_bytes'],U.k*ns*16)
        cost['explicit_IDFT_columns_allocated']=0
        cost['explicit_GU_columns_allocated']=0
        cost['explicit_H_columns_allocated']=0
        # Exact zero-padding carries the original Fourier coefficients between
        # nested masks; it preserves both current and the coefficient metric.
        J=JD+U.apply(A)
        Q=E+GJD+green.apply(J-JD)
        def SU_action(a):
            cost['S_columns']+=ns
            return S@U.apply(a)
        def SU_adjoint(z):
            cost['S_adjoint_columns']+=ns
            return U.adjoint(S.conj().T@z)
        old_grad = old_direction = None
        stage_stop = 'MAXITER'
        for it in range(int(stage.get('maxiter', 500))):
            if expired() or iterations >= int(config.get('max_iterations', 10**12)):
                status = stage_stop = 'WALL_BUDGET' if expired() else 'ITERATION_BUDGET'
                break
            # Fixed-m CIE residual transform; Fourier/Green are actions.
            p=beta_j/(1+beta_j*m)
            def H_action(a):
                ua=U.apply(a)
                return p[:,None]*(ua-m[:,None]*green.apply(ua))
            def H_adjoint(z):
                v=p.conj()[:,None]*z
                return U.adjoint(v-green.adjoint(m[:,None]*v))
            T=p[:,None]*(JD-m[:,None]*(E+GJD))
            D,R=SJD-Y+SU_action(A),T+H_action(A)
            cost['scalar_current_operator_builds']+=1
            cost['current_objective_evaluations'] += 1
            data_before = .5*_inner(D,D,Od,cost)
            state_before = .5*rho*_inner(R,R,Os,cost)
            grad = SU_adjoint(D@Od) + rho*H_adjoint(R@Os)
            cost['source_metric_products'] += 2
            cost['current_gradient_evaluations'] += 1
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
            DD, RR = SU_action(direction), H_action(direction)
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
            fluctuation=U.apply(A)
            J,Q=JD+fluctuation,E+GJD+green.apply(fluctuation)
            Dafter=SJD-Y+S@fluctuation
            cost['S_columns']+=ns
            Rafter=T+H_action(A)
            data_after_current=.5*_inner(Dafter,Dafter,Od,cost)
            state_after_current=.5*rho*_inner(Rafter,Rafter,Os,cost)
            m_before_readout = m.copy()
            if config['material_transform'] == 'onehop':
                m, inner = config.get('material_readout',onehop_material_readout)(blocks,J,Q,m,beta_m,Os,Ls,rho,
                                                   (lo,hi),reg,config,cost)
            elif sourcewise:
                m, inner = _baseline_readout(J,Q,m,beta_m,ws,rho,
                                             (lo,hi),reg,has_reg,config,cost)
            else:
                m, inner = _baseline_readout(J@Ls,Q@Ls,m,beta_m,np.ones(ns),rho,
                                             (lo,hi),reg,has_reg,config,cost)
            if not sourcewise or config['material_transform']=='onehop':
                cost['material_readout_source_whitening_products'] += 2
                cost['source_metric_products'] += 2
            old_grad, old_direction = grad, direction
            data_f = .5*_inner(Dafter,Dafter,Od,cost)
            def state_value(beta, transform):
                e = J-m[:, None]*Q
                if beta is not None:
                    if transform == 'onehop':
                        e = blocks.apply(m,beta,e)
                    else:
                        e = beta*e/(1+beta*m[:,None])
                return .5*rho*_inner(e,e,Os,cost)
            history.append(dict(stage=stage_index, iteration=it+1, global_iteration=iterations+1,
                                beta_current=beta_j, beta_material=beta_m, receiver_rank=receiver_rank,
                                **basis_info, coefficient_relative_change_squared=delta,
                                data_term=data_f, state_term_current=state_value(beta_j,config['current_transform']),
                                state_term_material=state_value(beta_m,config['material_transform']), regularizer=float(reg(m)[0]),
                                current_transform=config['current_transform'],material_transform=config['material_transform'],
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
                    m=m.copy(),J=J.copy(),coefficients=A.copy(),operation_counts=dict(cost)))
            iterations += 1
            cost['PR_CG_steps' if current_solver=='pr' else 'exact_current_steps'] += 1
            cost['matrix_free_current_steps'] += 1
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
