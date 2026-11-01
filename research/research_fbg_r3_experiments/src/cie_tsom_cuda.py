"""GPU-resident current block; unchanged CPU real-material CIE readout.

Uses the existing CUDA SharedOperators/ObjectiveLane, never another Maxwell
backend. Importing this module does not initialize CUDA or run an experiment.
"""
from collections import Counter
import importlib.util
import os
from pathlib import Path
import threading
import time

import numpy as np
from scipy.linalg import svd
import torch

from .cie_tsom_baseline import _readout, _regularizer, fourier_apcs

_PREPARE_LOCK = threading.Lock()


def common_material_prior(shape, config=None):
    """Same R4 width-19.5 finite-difference smooth + smoothed-TV prior."""
    cfg = dict(config or {})
    grid = np.arange(np.prod(shape)).reshape(shape)
    left = np.r_[grid[:, :-1].ravel(), grid[:-1].ravel()]
    right = np.r_[grid[:, 1:].ravel(), grid[1:].ravel()]
    smooth, tv, eps = (float(cfg.get(k, v)) for k, v in
                       [('smooth_weight', 1e-4), ('tv_weight', 1e-4), ('tv_epsilon', .01)])
    edges = len(left)
    def prior(m):
        d = (m[right]-m[left])/19.5
        q = np.sqrt(d*d+eps*eps)
        value = smooth*np.mean(d*d)+tv*np.mean(q-eps)
        gd = (2*smooth*d+tv*d/q)/(edges*19.5)
        grad = np.zeros_like(m)
        np.add.at(grad, right, gd)
        np.add.at(grad, left, -gd)
        return float(value), grad
    return prior


def _lane_class(shared):
    if hasattr(shared, 'cie_objective_lane_class'):
        return shared.cie_objective_lane_class
    root = Path(os.environ.get('CUDA_PIPELINE_ROOT',
                    str(Path(__file__).resolve().parents[2]/'research_z_alg_cuda_pipeline')))
    spec = importlib.util.spec_from_file_location('fbg_r3_cie_existing_device', root/'src/device.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.ObjectiveLane


def prepare_cie_shared(shared, config):
    """One immutable public receiver/Fourier preparation per rank/mask schedule.

    Call in the main thread before launching lanes to charge shared setup once.
    The cache contains no observation or material. Native Fourier coefficient
    coordinates are retained even after roundoff-only projected-rank cleanup.
    """
    masks = tuple(sorted(set(int(s['mf']) for s in config['stages'])))
    rank_key = ('rank', int(config['dpcs_rank'])) if 'dpcs_rank' in config else ('rtol', float(config['dpcs_rtol']))
    key = (rank_key, masks)
    with _PREPARE_LOCK:
        cache = getattr(shared, '_cie_public_cache', None)
        if cache is None:
            shared._cie_public_cache = cache = {}
        if key in cache:
            return cache[key], False
        started, counts = time.perf_counter(), Counter()
        host = getattr(shared, 'cie_host_operators', None)
        if host is None:
            torch.cuda.current_stream(shared.device).wait_event(shared.ready)
            S = shared.S.detach().cpu().numpy()
            counts.update(shared_public_D2H_bytes=S.nbytes, shared_public_D2H_returns=1)
        else:
            S = np.asarray(host.get('S', host.get('H')), complex)
        Us, sig, Vh = svd(S, full_matrices=False, check_finite=False)
        counts['receiver_SVD'] += 1
        rank = int(config.get('dpcs_rank', np.count_nonzero(sig > float(config.get('dpcs_rtol', 1e-3))*sig[0])))
        if not 0 <= rank <= len(sig):
            raise ValueError('Invalid frozen DPCS rank')
        Vd = Vh[:rank].conj().T
        def put(a):
            a = np.ascontiguousarray(a)
            counts['shared_basis_H2D_bytes'] += a.nbytes
            return torch.as_tensor(a, device=shared.device)
        stream = torch.cuda.current_stream(shared.device)
        stream.wait_event(shared.ready)
        with torch.no_grad():
            prepared = dict(Vd=put(Vd), Uh=put(Us[:, :rank].conj().T),
                            sigma=put(sig[:rank]), rank=rank, masks={})
            for mf in masks:
                Z, info = fourier_apcs((shared.side, shared.side), mf, Vd, counts)
                frequencies = info.pop('frequencies')
                U = put(Z)
                prepared['masks'][mf] = dict(U=U, GU=shared.G@U, SU=shared.S@U,
                                             frequencies=frequencies, info=info)
                counts['G_columns'] += U.shape[1]
                counts['S_columns'] += U.shape[1]
            prepared['ready'] = torch.cuda.Event()
            prepared['ready'].record(stream)
        # Setup completion is needed once before measuring task wall times.
        prepared['ready'].synchronize()
        prepared['cost'] = dict(counts, wall_seconds=time.perf_counter()-started)
        cache[key] = prepared
        return prepared, True


def run_cie_cuda(shared, Y, m0, config):
    """Same public CIE/TSOM/CSI API and readout as run_cie, CUDA current block.

    Each call owns one existing ObjectiveLane/private stream. All scalar-TM
    operators and coefficient products stay on that stream. One packed J/Q/
    diagnostic download and one material upload occur per alternating step.
    Config may supply a common regularizer; default is exactly the R4 prior.
    Shared preparation is charged separately in shared_setup_cost when new.
    """
    config = dict(config)
    started, cost = time.perf_counter(), Counter()
    torch.cuda.set_device(shared.device)
    method = config.get('method', 'cie')
    if method not in ('cie', 'tsom', 'csi'):
        raise ValueError('method must be cie, tsom or csi')
    m = np.asarray(m0, float).reshape(shared.n).copy()
    lo, hi = [np.broadcast_to(np.asarray(x, float), (shared.n,)).copy()
              for x in config.get('material_box', (.02, 20.))]
    if np.any(m < lo) or np.any(m > hi):
        raise ValueError('Initial material must satisfy the common box')
    config.setdefault('regularizer', common_material_prior((shared.side, shared.side), config))
    reg = _regularizer(config)
    rho = float(config.get('rho', 1.))
    if method == 'csi':
        prepared, created = None, False
    else:
        prepared, created = prepare_cie_shared(shared, config)
    lane = _lane_class(shared)(shared, {'Y': Y}, stage=shared.side, objective='data', config=config)
    context = config.get('_context')
    if context is not None:
        context.setdefault('lanes', []).append(lane)
        context['current_costs'] = cost
    stream = lane.stream
    if prepared is not None:
        stream.wait_event(prepared['ready'])
    shape, ns, n = (shared.side, shared.side), shared.ns, shared.n
    pull_buffer = torch.empty(4*n*ns+2, dtype=torch.float64, pin_memory=True)
    material_host = torch.empty(n, dtype=torch.float64, pin_memory=True)
    material_gpu = torch.empty(n, dtype=torch.float64, device=shared.device)
    history, candidates = [], []
    def validate(z, stage, reason):
        physical = Counter()
        for other in context.get('lanes', []) if context is not None else [lane]:
            physical.update(other.counts)
        matching = lane.last_x is not None and np.array_equal(z, lane.last_x)
        extra_lu, extra_rhs = (0, 0) if matching else (1, ns)
        rhs = sum(physical[k] for k in ('source_RHS', 'adjoint_RHS', 'receiver_RHS', 'tangent_RHS'))
        if physical['LU']+extra_lu > float(config.get('max_lu', np.inf)):
            return 'LU_CAP'
        if rhs+extra_rhs > float(config.get('max_rhs', np.inf)):
            return 'RHS_CAP'
        point = lane.close(z)
        candidate = dict(m=z.copy(), stage=stage, reason=reason,
                         objective=point['ranking_objective'], physical_nll=.5*point['relative_residual']**2,
                         physical_relative_residual=point['relative_residual'],
                         state_relative_residual=point['state_relative_residual'], lu_info=point['lu_info'])
        if point['lu_info'] or not np.isfinite(candidate['objective']):
            raise ArithmeticError('Original Maxwell endpoint validation failed')
        candidates.append(candidate)
        cost['physical_validations'] += 1
        cost['validation_D2H_bytes'] += (n+4)*8
        return None
    def expired():
        check = config.get('budget_check')
        if check is not None:
            check()
        return time.perf_counter()-started >= float(config.get('max_wall_seconds', np.inf))
    first_limit = validate(m, -1, 'initial_incumbent')
    if first_limit is not None:
        raise RuntimeError('Cannot validate initial incumbent within '+first_limit)
    with torch.cuda.stream(stream), torch.no_grad():
        target = lane.target
        wd = torch.as_tensor(config.get('data_weights', 1/np.maximum(np.sum(abs(Y)**2, axis=0), 1e-300)), dtype=torch.float64, device=shared.device)
        ws = torch.as_tensor(config.get('state_weights', np.ones(ns)), dtype=torch.float64, device=shared.device)
        if 'state_weights' not in config:
            ws = 1/shared.E0.abs().square().sum(0).clamp_min(1e-300)
        ws_host = lane._pull(ws)
        cost['state_weights_D2H_bytes'] += ws_host.nbytes
        if method == 'csi':
            JD = torch.zeros_like(shared.E0)
            receiver_rank = 0
        else:
            JD = prepared['Vd']@((prepared['Uh']@target)/prepared['sigma'][:, None])
            receiver_rank = prepared['rank']
        GJD, SJD = shared.G@JD, shared.S@JD
        cost.update(G_columns=ns, S_columns=ns)
        A, old_frequencies = None, []
        iterations, status = 0, 'SCHEDULE_COMPLETE'
        for stage_index, stage in enumerate(config['stages']):
            if expired():
                status = 'WALL_BUDGET'
                break
            bm = float(stage.get('beta_material', stage.get('beta', 1.))) if method == 'cie' else None
            bj = float(stage.get('beta_current', bm)) if method == 'cie' else None
            if method == 'cie' and (bm <= 0 or bj <= 0):
                raise ValueError('Positive CIE beta required')
            if method == 'csi':
                U, GU, SU = shared.eye, shared.G, shared.S
                info = dict(mf=None, fourier_columns=n, apcs_rank=n)
                A = torch.zeros_like(JD) if A is None else A
            else:
                basis = prepared['masks'][int(stage['mf'])]
                U, GU, SU, info = basis['U'], basis['GU'], basis['SU'], basis['info']
                frequencies = basis['frequencies']
                enlarged = torch.zeros((len(frequencies), ns), dtype=torch.complex128, device=shared.device)
                if A is not None:
                    lookup = {key:i for i, key in enumerate(frequencies)}
                    pairs = [(lookup[key], i) for i, key in enumerate(old_frequencies) if key in lookup]
                    if pairs:
                        new_idx, old_idx = zip(*pairs)
                        enlarged[list(new_idx)] = A[list(old_idx)]
                A, old_frequencies = enlarged, frequencies
            old_grad = old_direction = None
            stage_stop = 'MAXITER'
            for it in range(int(stage.get('maxiter', 500))):
                if expired() or iterations >= int(config.get('max_iterations', 10**12)):
                    status = stage_stop = 'WALL_BUDGET' if expired() else 'ITERATION_BUDGET'
                    break
                material_host.copy_(torch.from_numpy(m))
                material_gpu.copy_(material_host, non_blocking=True)
                cost['material_H2D_bytes'] += m.nbytes
                cost['material_H2D_returns'] += 1
                if bj is None:
                    H, T = U-material_gpu[:, None]*GU, JD-material_gpu[:, None]*(shared.E0+GJD)
                else:
                    rj = bj*material_gpu/(1+bj*material_gpu)
                    H = bj*U-rj[:, None]*(GU+bj*U)
                    T = bj*JD-rj[:, None]*(shared.E0+GJD+bj*JD)
                D, R = SJD-target+SU@A, T+H@A
                grad = (SU.mH@D)*wd[None]+rho*(H.mH@R)*ws[None]
                if old_grad is None:
                    direction = -grad
                else:
                    pr = ((grad.conj()*(grad-old_grad)).real.sum(0)/old_grad.abs().square().sum(0).clamp_min(1e-300)).clamp_min(0)
                    direction = -grad+old_direction*pr
                gd = (grad.conj()*direction).real.sum(0)
                direction = torch.where((gd < 0)[None], direction, -grad)
                DD, RR = SU@direction, H@direction
                numerator = -(grad.conj()*direction).real.sum(0)
                denominator = wd*DD.abs().square().sum(0)+rho*ws*RR.abs().square().sum(0)
                alpha = torch.where(denominator > 1e-300, numerator/denominator, torch.zeros_like(numerator))
                change = direction*alpha
                delta = (change.abs().square().sum(0)/A.abs().square().sum(0).clamp_min(1e-300)).mean()
                A = A+change
                J, Q = JD+U@A, shared.E0+GJD+GU@A
                data_value = .5*(wd[None]*(SJD-target+SU@A).abs().square()).sum()
                packed = torch.cat((torch.view_as_real(J).reshape(-1), torch.view_as_real(Q).reshape(-1), delta.reshape(1), data_value.reshape(1)))
                host = lane._pull(packed, pull_buffer)
                cost['current_D2H_bytes'] += host.nbytes
                cost['current_D2H_returns'] += 1
                Jcpu = host[:2*n*ns].view(np.complex128).reshape(n, ns)
                Qcpu = host[2*n*ns:4*n*ns].view(np.complex128).reshape(n, ns)
                m, inner = _readout(Jcpu, Qcpu, m, bm, ws_host, rho, (lo, hi), reg,
                                    config.get('regularizer') is not None, config, cost)
                old_grad, old_direction = grad, direction
                def state_value(beta):
                    e = Jcpu-m[:, None]*Qcpu
                    if beta is not None:
                        e = beta*e/(1+beta*m[:, None])
                    return .5*rho*float(np.sum(ws_host[None]*abs(e)**2))
                history.append(dict(stage=stage_index, iteration=it+1, global_iteration=iterations+1,
                                    beta_current=bj, beta_material=bm, receiver_rank=receiver_rank,
                                    **info, coefficient_relative_change_squared=float(host[-2]),
                                    data_term=float(host[-1]), state_term_current=state_value(bj),
                                    state_term_material=state_value(bm), regularizer=reg(m)[0],
                                    wall=time.perf_counter()-started, **inner))
                iterations += 1
                cost['PR_CG_steps'] += 1
                cost['reduced_operator_columns'] += 6*ns
                if it+1 >= int(config.get('miniter', 2)) and host[-2] < float(config.get('coefficient_tol', 1e-3)):
                    stage_stop = 'COEFFICIENT_TOL'
                    break
            physical_limit = validate(m, stage_index, stage_stop)
            if physical_limit is not None:
                status = stage_stop = physical_limit
            history.append(dict(stage=stage_index, event='stage_end', status=stage_stop,
                                beta_current=bj, beta_material=bm, wall=time.perf_counter()-started, **info))
            if status != 'SCHEDULE_COMPLETE':
                break
    # Only physically evaluated candidates compete. No auxiliary/GT selection.
    best = min(candidates, key=lambda row:row['objective'])
    cost['wall_seconds'] = time.perf_counter()-started
    return dict(m=best['m'].copy(), candidates=candidates, history=history, cost=dict(cost, model=dict(lane.counts)),
                last_joint_m=m.copy(), status=status, method=method, receiver_rank=receiver_rank,
                shared_setup_cost=prepared['cost'] if created else {},
                shared_setup_reused=not created,
                host_result_wait_seconds=lane.wait_seconds,
                validation_scope='original complex128 scalar-TM Maxwell, all cells and sources')
