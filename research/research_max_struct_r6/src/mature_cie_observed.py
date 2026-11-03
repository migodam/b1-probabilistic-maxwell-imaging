"""Observed mature CIE-I / Dual-CIE, with the original sourcewise PR loop.

Derived from research/research_fbg_r3_experiments/src/cie_tsom_baseline.py
(run_cie). The audit diff contains ONLY instrumentation: no update, readout,
Fourier transport, line search, objective, or stopping changes. Helpers are
imported from their original implementation. Historical files stay untouched.
"""
from collections import Counter
from copy import deepcopy
import time
import numpy as np
from scipy.linalg import svd
from research.research_fbg_r3_experiments.src.cie_tsom_baseline import (
    fourier_apcs, _regularizer, _readout)

def run_cie(model, Y, m0, config):
    """Return dict(m,candidates,history,cost,last_joint_m,status).

    Required config: stages=[dict(mf=..., beta=..., maxiter=...)],
    dpcs_rank OR dpcs_rtol (noise/calibration-frozen). A stage's optional
    beta_current gives the transported two-beta variant. method='tsom'
    uses ordinary LSIE with Fourier APCS; method='csi' uses unrestricted
    current and ordinary LSIE. All modes use the same real material box.

    Optional config: regularizer(m)->(R,grad_R), physical_whiten(residual),
    material_box=(.02,20), rho=1, data_weights/state_weights (positive source
    weights), max_wall_seconds, budget_check(), max_iterations (global),
    coefficient_tol=1e-3 (source-wise SQUARED relative coefficient change),
    miniter=2, sourcewise_cg=True. Physical whiten defaults to division by
    physical_sigma (default 1). Source normalizations default to the paper.
    Every retained candidate is re-evaluated with the original Maxwell model.
    """
    config = dict(config)
    started = time.perf_counter()
    cpu_started = time.process_time()
    cost = Counter()
    model_before = deepcopy(dict(getattr(model, 'cost', {})))
    G, E, S = [np.asarray(x, complex) for x in (model.G, model.E, model.S)]
    Y, m = np.asarray(Y, complex), np.asarray(m0, float).copy()
    n, ns = E.shape
    shape = tuple(model.shape)
    if G.shape != (n, n) or S.shape[1] != n or Y.shape != (S.shape[0], ns):
        raise ValueError('Inconsistent public operators / observations')
    if m.shape != (n,) or not np.all(np.isfinite(m)) or np.prod(shape) != n:
        raise ValueError('Invalid full real material/grid')
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
    wd = np.asarray(config.get('data_weights', 1/np.maximum(np.sum(abs(Y)**2, axis=0), 1e-300)), float)
    ws = np.asarray(config.get('state_weights', 1/np.maximum(np.sum(abs(E)**2, axis=0), 1e-300)), float)
    wd, ws = np.broadcast_to(wd, (ns,)), np.broadcast_to(ws, (ns,))
    if np.any(wd <= 0) or np.any(ws <= 0) or not np.all(np.isfinite(np.r_[wd, ws])):
        raise ValueError('Finite positive fixed source weights required')
    whiten = config.get('physical_whiten', lambda a: a/float(config.get('physical_sigma', 1.)))
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
                rj = beta_j*m/(1+beta_j*m)
                H = beta_j*U-rj[:, None]*(GU+beta_j*U)
                T = beta_j*JD-rj[:, None]*(E+GJD+beta_j*JD)
            D, R = SJD-Y+SU@A, T+H@A
            grad = (SU.conj().T@D)*wd[None, :] + rho*(H.conj().T@R)*ws[None, :]
            cost['current_gradient_evaluations'] += 1
            if old_grad is None:
                pr = 0.
                direction = -grad
            else:
                axis = 0 if config.get('sourcewise_cg', True) else None
                pr = np.maximum(0., np.sum(np.real(grad.conj()*(grad-old_grad)), axis=axis) /
                                np.maximum(np.sum(abs(old_grad)**2, axis=axis), 1e-300))
                direction = -grad+old_direction*pr
            axis = 0 if config.get('sourcewise_cg', True) else None
            gd = np.sum(np.real(grad.conj()*direction), axis=axis)
            restarted = bool(np.any(gd >= 0))
            if axis == 0:
                bad = gd >= 0
                direction[:, bad] = -grad[:, bad]
            elif gd >= 0:
                direction = -grad
            DD, RR = SU@direction, H@direction
            numerator = -np.sum(np.real(grad.conj()*direction), axis=axis)
            denom_columns = wd*np.sum(abs(DD)**2, axis=0)+rho*ws*np.sum(abs(RR)**2, axis=0)
            denom = denom_columns if axis == 0 else np.sum(denom_columns)
            alpha = np.divide(numerator, denom, out=np.zeros_like(numerator), where=denom > 1e-300)
            change = direction*alpha
            delta = float(np.mean(np.sum(abs(change)**2, axis=0)/np.maximum(np.sum(abs(A)**2, axis=0), 1e-300)))
            A = A+change
            J, Q = JD+U@A, E+GJD+GU@A
            m, inner = _readout(J, Q, m, beta_m, ws, rho, (lo, hi), reg, has_reg, config, cost)
            old_grad, old_direction = grad, direction
            data_f = .5*float(np.sum(wd[None, :]*abs(SJD-Y+SU@A)**2))
            def state_value(beta):
                e = J-m[:, None]*Q
                if beta is not None:
                    e = beta*e/(1+beta*m[:, None])
                return .5*rho*float(np.sum(ws[None, :]*abs(e)**2))
            history.append(dict(stage=stage_index, iteration=it+1, global_iteration=iterations+1,
                                beta_current=beta_j, beta_material=beta_m, receiver_rank=receiver_rank,
                                **basis_info, coefficient_relative_change_squared=delta,
                                data_term=data_f, state_term_current=state_value(beta_j),
                                state_term_material=state_value(beta_m), regularizer=float(reg(m)[0]),
                                wall=time.perf_counter()-started,
                                cpu_seconds=time.process_time()-cpu_started,
                                sourcewise_cg=bool(config.get('sourcewise_cg', True)),
                                PR_plus=np.asarray(pr).tolist(),
                                line_alpha=np.asarray(alpha).tolist(),
                                direction_restarted=restarted, **inner))
            if config.get('iteration_callback'):
                config['iteration_callback'](dict(stage=stage_index, iteration=it+1,
                    m=m.copy(), J=J.copy(), coefficients=A.copy(),
                    operation_counts=dict(cost)))
            iterations += 1
            cost['PR_CG_steps'] += 1
            cost['reduced_operator_columns'] += 6*ns
            if it+1 >= int(config.get('miniter', 2)) and delta < float(config.get('coefficient_tol', 1e-3)):
                stage_stop = 'COEFFICIENT_TOL'
                break
        if not np.all(np.isfinite(m)):
            raise FloatingPointError('Nonfinite real material readout')
        validate(m, stage_index, stage_stop)
        stage_state = dict(stage=stage_index, m=m.copy(), J=J.copy(), Q=Q.copy(),
                           coefficients=A.copy(), beta_current=beta_j,
                           beta_material=beta_m, receiver_rank=receiver_rank,
                           basis_info=dict(basis_info))
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
    cost['cpu_seconds'] = time.process_time()-cpu_started
    model_after = dict(getattr(model, 'cost', {}))
    model_delta = {k: model_after.get(k, 0)-model_before.get(k, 0)
                   for k in set(model_before)|set(model_after)
                   if isinstance(model_after.get(k, 0), (int, float)) and isinstance(model_before.get(k, 0), (int, float))}
    return dict(m=best['m'].copy(), candidates=candidates, history=history, stage_states=stage_states,
                cost=dict(cost, model=model_delta), last_joint_m=m.copy(), status=status,
                method=method, receiver_rank=receiver_rank)
