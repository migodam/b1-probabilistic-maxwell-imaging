"""Ordinary first-order predictor / transverse Gauss--Newton exploration.

No solver, observations, truth, or unmetered Maxwell operation lives here.
The caller supplies ``evaluate(x, need_jac=False)`` and a decoder-only callback.
Weak coordinates are affine functionals held EXACTLY fixed during a corrector.
H1/H2 differ only in the exploration-frame refresh, never the corrector.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from pathlib import Path
import json
import time
import numpy as np
from scipy.linalg import null_space, cho_factor, cho_solve


class LocalBudgetStop(RuntimeError):
    pass


@dataclass
class Config:
    arm: str = "H1"
    weak_rank: int = 16
    layers: int = 3
    proposals_per_layer: int = 4
    predictor_rms: float = .02
    layer_radii: tuple = (.02, .05, .10)
    corrector_states: int = 5
    line_search_trials: int = 4
    local_rhs_cap: int = 5400
    local_lu_cap: int = 240
    reserve_final_rhs: int = 8
    rms_denominator: int = 576
    latent_trust: float = .5
    material_trust: float = .01
    max_latent_radius: float = 4.
    damping: float = 1e-7
    prior_weight: float = 1.
    max_prior_increase: float = 10.
    max_phi_increase: float = 2.
    noise_phi_max: float = 282.8737921622988
    material_radius_relative_tolerance: float = .1
    refresh_layer: int = 1
    adaptive_secant_threshold: float = .25
    adaptive_drift_threshold: float = .25
    adaptive_failure_fraction: float = .5
    adaptive_max_updates: int = 2
    distinct_rms: float = .002
    adaptive_saturation_layers: int = 2
    basis_kind: str = "full"
    legacy_gates_required: bool = False
    full_coordinate_corrector: bool = False
    native_box_retraction: bool = False
    active_set_iterations: int = 4
    seed: int = 20261030


def _plain(x):
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, np.generic): return x.item()
    if isinstance(x, dict): return {k: _plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [_plain(v) for v in x]
    return x


def material_matrix(D, denominator=576):
    D = np.asarray(D)
    return np.concatenate((D.real, D.imag), axis=0) / np.sqrt(denominator)


def material_rms(a, b, denominator=576):
    return float(np.linalg.norm(np.asarray(a)-np.asarray(b))/np.sqrt(denominator))


def canonical_basis(D, denominator=576, rtol=1e-10, atol=0.):
    """Latent lift with material-orthonormal image; removes decoder nullspace."""
    M = material_matrix(D, denominator)
    _, s, Vh = np.linalg.svd(M, full_matrices=False)
    cutoff = max(float(s[0])*max(rtol, np.finfo(float).eps*max(M.shape)), atol) if len(s) else atol
    keep = s > cutoff
    T = Vh[keep].T/s[keep][None, :]
    return T, dict(material_singular=s, material_cutoff=cutoff,
                   effective_material_rank=int(keep.sum()),
                   decoder_gauge_dimension=int(M.shape[1]-keep.sum()))


def weak_frame(state, config, basis_factory=None):
    J, D = np.asarray(state['jac']), np.asarray(state['dmaterial'])
    if basis_factory is not None:
        supplied = basis_factory(state, config)
        W = np.asarray(supplied['W'] if isinstance(supplied, dict) else supplied, float)
        info = supplied.get('info', {}) if isinstance(supplied, dict) else {}
    else:
        T, info = canonical_basis(D, config.rms_denominator)
        if not T.shape[1]: raise ValueError('No effective material tangent')
        K = J@T
        _, singular, Vh = np.linalg.svd(K, full_matrices=True)
        count = min(config.weak_rank, T.shape[1])
        W = T@Vh[-count:][::-1].T
        ev = np.r_[np.zeros(max(0, T.shape[1]-len(singular))), singular[::-1]**2]
        info.update(data_curvature_eigenvalues=ev, weak_extraction='factor_SVD_not_normal_eigenproblem')
    # Reorthogonalize in the material metric, including supplied pair frames.
    M = material_matrix(D, config.rms_denominator)
    gram = (M@W).T@(M@W)
    if not np.allclose(gram, np.eye(W.shape[1]), rtol=1e-7, atol=1e-8):
        _, s, Vh = np.linalg.svd(M@W, full_matrices=False)
        keep = s > max((s[0] if len(s) else 0)*1e-10, 1e-12)
        W = W@Vh[keep].T/s[keep][None, :]
    for i in range(W.shape[1]):
        if W[np.argmax(abs(W[:, i])), i] < 0: W[:, i] *= -1
    C = (M@W).T@M
    return dict(W=W, C=C, info=info, anchor_x=np.asarray(state['x']).copy())


class PaidEvaluator:
    """Conservative equal per-arm cap, in addition to root atomic prepayment.

    A cached anchor supplied by the caller costs zero in this arm. Every new
    requested J state reserves 40 local RHS; every F state reserves 8. Root
    receipts remain authoritative for physical cost. No failed call is refunded.
    """
    def __init__(self, evaluate, cap, reserve_final=8, initial=None, lu_cap=240):
        self.evaluate_raw, self.cap = evaluate, int(cap)
        self.reserve_final = int(reserve_final)
        self.rhs = 0
        self.lu_cap = int(lu_cap)
        self.calls = []
        self.last = initial

    def evaluate(self, x, need_jac=False, final=False):
        if (self.last is not None and np.array_equal(x, self.last['x'])
                and (not need_jac or 'jac' in self.last)):
            return self.last
        amount = 40 if need_jac else 8
        if self.rhs+amount+(0 if final else self.reserve_final) > self.cap:
            raise LocalBudgetStop('Matched local cap, final F slot retained')
        if len(self.calls)+1+(0 if final else 1) > self.lu_cap:
            raise LocalBudgetStop('Matched local LU cap, final F slot retained')
        self.rhs += amount
        row = dict(need_jac=bool(need_jac), local_RHS=amount,
                   cumulative_RHS=self.rhs, status='ATTEMPTED')
        self.calls.append(row)
        try:
            state = self.evaluate_raw(np.asarray(x, float), need_jac=need_jac)
        except BaseException as error:
            row.update(status='FAILED', error=repr(error)); raise
        row.update(status='OK', ticket=state.get('ticket'))
        self.last = state
        return state


def _decode_material(decoder, x):
    """Decoder returns dict(material,legal) or shared Problem.decode tuple."""
    x = np.asarray(x, float)
    owner = getattr(decoder, '__self__', None)
    legal = True
    if owner is not None and hasattr(owner, 'project'):
        legal = bool(np.linalg.norm(owner.project(x)-x) <= 1e-10)
    z = decoder(x)
    if isinstance(z, dict): return np.asarray(z['material']), legal and bool(z.get('prior_legal', True))
    if isinstance(z, tuple): return np.asarray(z[2]), legal
    return np.asarray(z), legal


def predictor(state, direction, config, decoder=None):
    x = np.asarray(state['x'])
    M = material_matrix(state['dmaterial'], config.rms_denominator)
    speed = float(np.linalg.norm(M@direction))
    if speed < 1e-12: return None, dict(status='ZERO_MATERIAL_DIRECTION')
    v = direction/speed
    h = config.predictor_rms
    owner = getattr(decoder, '__self__', None) if decoder is not None else None
    box_retraction = bool(config.native_box_retraction and owner is not None and
                          getattr(owner, 'kind', 'P1') != 'P1' and hasattr(owner, 'project'))
    def retract(z):
        return owner.project(z) if box_retraction else z
    if decoder is not None:
        # Decoder arithmetic alone matches the nonlinear material radius.
        def distance(t):
            material, legal = _decode_material(decoder, retract(x+t*v))
            return material_rms(material, state['material'], config.rms_denominator), legal
        # The local material-unit lift may have huge latent norm. The finite
        # decoder can curve substantially before that linear extrapolation;
        # bracket only inside the declared latent radius and measure it exactly.
        h_limit = config.max_latent_radius/max(float(np.linalg.norm(v)), 1e-300)
        lo, hi = 0., min(h, h_limit)
        dh, legal = distance(hi)
        for _ in range(16):
            if dh >= config.predictor_rms or hi >= h_limit: break
            hi = min(2.*hi, h_limit); dh, legal = distance(hi)
        if dh < config.predictor_rms or not legal:
            return None, dict(status='PREDICTOR_RADIUS_UNAVAILABLE',
                              box_retraction=box_retraction,
                              largest_tested_material_rms=float(dh))
        for _ in range(40):
            mid = .5*(lo+hi); dm, legal = distance(mid)
            if not legal: hi = mid
            elif dm < config.predictor_rms: lo = mid
            else: hi = mid
        h = .5*(lo+hi)
    if np.linalg.norm(h*v) > config.max_latent_radius:
        return None, dict(status='PREDICTOR_LATENT_RADIUS_UNAVAILABLE')
    raw = x+h*v
    xp = retract(raw)
    retraction_norm = float(np.linalg.norm(xp-raw))
    return xp, dict(status='PREDICTOR', scalar=h,
                   latent_step=float(np.linalg.norm(xp-x)),
                   raw_latent_step=float(np.linalg.norm(h*v)),
                   linear_material_rms=float(h), box_retraction=box_retraction,
                   retraction_norm=retraction_norm,
                   leaves_linear_weakspan=bool(retraction_norm > 1e-12),
                   weakspan_scope='box retraction is not constrained to the original linear weak span' if box_retraction else 'unretracted weak ray')


def _objective(state, config, prior_center, residual_shift):
    r = np.asarray(state['residual'])-residual_shift
    if 'reg_value' in state:
        return .5*float(r@r)+float(state['reg_value'])
    z = np.asarray(state['x'])-prior_center
    return .5*float(r@r+config.prior_weight*(z@z))


def transverse_correct(paid, xpred, frame, config, *, decoder=None,
                       prior_center=None, residual_shift=None, reference=None):
    """Dynamic first-order GN with actual-descent backtracking.

    ``corrector_states`` bounds paid Jacobians (default five). Every step is
    in the current effective material tangent AND the one fixed ker(C).
    Trial forwards never update C. Decoder-only feasibility checks happen
    before physical calls; every failed physical attempt remains charged.
    """
    dimension = len(xpred)
    prior_center = np.zeros(dimension) if prior_center is None else np.asarray(prior_center)
    residual_shift = np.zeros_like(reference['residual']) if residual_shift is None else np.asarray(residual_shift)
    C = np.zeros((0, dimension)) if frame is None else frame['C']
    x = np.asarray(xpred, float).copy()
    history, states = [], []
    best, best_obj, stop = None, np.inf, 'J_CAP'
    cholesky_count, cholesky_seconds = 0, 0.
    owner = getattr(decoder, '__self__', None) if decoder is not None else None
    def feasible(z):
        if np.linalg.norm(z-reference['x']) > config.max_latent_radius*(1+1e-10): return False
        if decoder is not None:
            _, legal = _decode_material(decoder, z)
            if not legal: return False
        return bool(np.linalg.norm(C@(z-xpred)) <= 1e-7)
    if not feasible(x):
        return None, dict(status='PREDICTOR_CONSTRAINT_FAIL', history=[], evaluated_states=0,
                          normal_cholesky=0, normal_cholesky_seconds=0.), []
    for iteration in range(config.corrector_states):
        try: state = paid.evaluate(x, need_jac=True)
        except LocalBudgetStop: stop = 'LOCAL_BUDGET_STOP'; break
        states.append(state)
        obj = _objective(state, config, prior_center, residual_shift)
        coordinate_defect = float(np.linalg.norm(C@(state['x']-xpred)))
        same_parameter = np.linalg.norm(state['x']-x) <= 1e-10
        row = dict(iteration=iteration, phi=float(state['phi']), objective=obj,
                   weak_coordinate_defect=coordinate_defect,
                   evaluator_projected=not same_parameter, ticket=state.get('ticket'), trials=[])
        history.append(row)
        if not same_parameter:
            stop = 'PARAMETER_PROJECTION_BREAKS_FIXED_COORDINATE'; break
        if bool(state.get('prior_legal', True)) and feasible(state['x']) and obj < best_obj:
            best, best_obj = state, obj
        T, geometry = canonical_basis(state['dmaterial'], config.rms_denominator)
        correction_lift = np.eye(dimension) if config.full_coordinate_corrector else T
        constraints = C.copy()
        J, M = state['jac'], material_matrix(state['dmaterial'], config.rms_denominator)
        r = state['residual']-residual_shift
        if 'reg_value' in state:
            regularizer_gradient = np.asarray(state['reg_gradient'], float)
            regularizer_hessian = np.asarray(state['reg_hessian'], float)
        else:
            regularizer_gradient = config.prior_weight*(state['x']-prior_center)
            regularizer_hessian = config.prior_weight*np.eye(dimension)
        gradient = J.T@r+regularizer_gradient
        hessian = J.T@J+regularizer_hessian
        step = None
        active_coordinates = set()
        active_limit = max(config.active_set_iterations, dimension+1) if config.native_box_retraction else config.active_set_iterations
        for active_set_iteration in range(active_limit):
            N = null_space(constraints@correction_lift, rcond=1e-10) if constraints.shape[0] else np.eye(correction_lift.shape[1])
            Z = correction_lift@N
            if not Z.shape[1]: break
            gram = Z.T@hessian@Z
            ridge = config.damping*max(float(np.max(np.diag(gram), initial=0.)), 1.)
            ts = time.perf_counter()
            coefficient = -cho_solve(cho_factor(gram+ridge*np.eye(Z.shape[1]), lower=True, check_finite=False), Z.T@gradient, check_finite=False)
            cholesky_count += 1; cholesky_seconds += time.perf_counter()-ts
            step = Z@coefficient
            # A native decoder has a box on its first 256 occupancy coordinates.
            # Add outward active coordinates to the SAME affine constraints;
            # simply clipping a transverse step would destroy fixed weak C.
            if owner is not None and getattr(owner, 'kind', 'P1') != 'P1':
                nbox = dimension-2
                outward_tolerance = 1e-10 if config.native_box_retraction else 0.
                bad = np.flatnonzero(((state['x'][:nbox] <= 1e-10) & (step[:nbox] < -outward_tolerance)) |
                                     ((state['x'][:nbox] >= 1-1e-10) & (step[:nbox] > outward_tolerance)))
                if len(bad):
                    if config.native_box_retraction:
                        bad = np.array([int(j) for j in bad if int(j) not in active_coordinates], dtype=int)
                        if not len(bad): break
                    active_coordinates.update(bad.tolist())
                    constraints = np.vstack((C, np.eye(dimension)[sorted(active_coordinates)]))
                    step = None
                    continue
            break
        if step is None:
            stop = 'NO_FEASIBLE_TRANSVERSE_TANGENT'; break
        scale = min(1., config.latent_trust/max(np.linalg.norm(step), 1e-300),
                    config.material_trust/max(np.linalg.norm(M@step), 1e-300))
        step *= scale
        fraction_max, fraction_scale = 1., 1.
        if config.native_box_retraction and owner is not None and getattr(owner, 'kind', 'P1') != 'P1':
            nbox = dimension-2
            for coordinate in range(nbox):
                if step[coordinate] > 1e-10:
                    fraction_max = min(fraction_max, (1.-state['x'][coordinate])/step[coordinate])
                elif step[coordinate] < -1e-10:
                    fraction_max = min(fraction_max, -state['x'][coordinate]/step[coordinate])
            fraction_max = max(0., float(fraction_max))
            fraction_scale = float(np.nextafter(fraction_max, 0.)) if fraction_max < 1. else 1.
            step *= fraction_scale
        row.update(transverse_rank=int(Z.shape[1]), latent_step=float(np.linalg.norm(step)),
                   linear_material_rms=float(np.linalg.norm(M@step)),
                   linear_weak_defect=float(np.linalg.norm(C@step)),
                   active_set_iterations=active_set_iteration+1,
                   active_set_limit=active_limit, active_coordinates=len(active_coordinates),
                   fraction_to_boundary_max=fraction_max,
                   fraction_to_boundary_scale=fraction_scale, **geometry)
        if np.linalg.norm(step) < 1e-10:
            stop = 'ZERO_TRANSVERSE_STEP'; break
        accepted = False
        for trial in range(config.line_search_trials):
            alpha = .5**trial
            trial_step = alpha*step
            xn = state['x']+trial_step
            record = dict(trial=trial, alpha=alpha)
            row['trials'].append(record)
            if not feasible(xn):
                record['status'] = 'DECODER_OR_LATENT_CONSTRAINT'; continue
            if decoder is not None:
                trial_material, _ = _decode_material(decoder, xn)
                if material_rms(trial_material, state['material'], config.rms_denominator) > config.material_trust*(1+1e-8):
                    record['status'] = 'ACTUAL_MATERIAL_TRUST'; continue
            pred = -float(gradient@trial_step+.5*trial_step@hessian@trial_step)
            record['predicted_reduction'] = pred
            if pred <= 0:
                record['status'] = 'NONPOSITIVE_PREDICTION'; continue
            try: trial_state = paid.evaluate(xn, need_jac=False)
            except LocalBudgetStop:
                stop = 'LOCAL_BUDGET_STOP'; break
            states.append(trial_state)
            trial_obj = _objective(trial_state, config, prior_center, residual_shift)
            actual = obj-trial_obj
            nonlinear = trial_state['residual']-state['residual']-J@trial_step
            response_change = trial_state['residual']-state['residual']
            defect = float(np.linalg.norm(nonlinear)/max(np.linalg.norm(response_change), 1.))
            record.update(actual_reduction=actual, actual_predicted_ratio=actual/pred,
                          secant_relative_defect=defect, phi=float(trial_state['phi']),
                          weak_coordinate_defect=float(np.linalg.norm(C@(trial_state['x']-xpred))),
                          actual_material_rms=material_rms(trial_state['material'], state['material'], config.rms_denominator),
                          ticket=trial_state.get('ticket'))
            row['secant_relative_defect'] = max(row.get('secant_relative_defect', 0.), defect)
            if (actual > 0 and bool(trial_state.get('prior_legal', True)) and
                    np.linalg.norm(trial_state['x']-xn) <= 1e-10 and feasible(trial_state['x'])):
                record['status'] = 'ACCEPTED'
                if trial_obj < best_obj: best, best_obj = trial_state, trial_obj
                x = trial_state['x'].copy(); accepted = True; break
            record['status'] = 'REJECTED_NO_ACTUAL_DESCENT'
        if stop == 'LOCAL_BUDGET_STOP': break
        if not accepted:
            stop = 'LINE_SEARCH_STOP'; break
    return best, dict(status=stop, history=history, evaluated_states=len(states),
                      normal_cholesky=cholesky_count,
                      normal_cholesky_seconds=cholesky_seconds), states


def qualification(state, anchor, config):
    prior_delta = .5*float(state['x']@state['x']-anchor['x']@anchor['x'])
    phi_delta = float(state['phi']-anchor['phi'])
    legal = bool(state.get('prior_legal', True))
    legacy_passed = prior_delta <= config.max_prior_increase and phi_delta <= config.max_phi_increase
    passed = legal and state['phi'] <= config.noise_phi_max and (legacy_passed or not config.legacy_gates_required)
    return bool(passed), dict(prior_increase=prior_delta, phi_increase=phi_delta,
                             legal=legal, legacy_passed=bool(legacy_passed),
                             noise_gate=bool(state['phi'] <= config.noise_phi_max))


def frame_drift(frame, candidate, D, denominator=576):
    M = material_matrix(D, denominator)
    A, B = M@frame['W'], M@candidate['W']
    Qa, _ = np.linalg.qr(A); Qb, _ = np.linalg.qr(B)
    rank = min(Qa.shape[1], Qb.shape[1])
    if not rank: return 1.
    overlap = float(np.linalg.norm(Qa.T@Qb, 'fro')**2)
    return float(np.sqrt(max(0., 1.-overlap/rank)))


def run_arm(evaluate, x0, config=None, *, initial_state=None, decoder=None,
            basis_factory=None, random_draws=None, prior_targets=None,
            residual_targets=None, output=None, callback=None):
    """Run finite heuristic clouds; no posterior/density assertion.

    H0 is ordinary full-coordinate RML when prior_targets/residual_targets are
    supplied. H1 fixes exploration frame. H2 refreshes it exactly once before
    ``refresh_layer``. HA refreshes only on declared feedback thresholds and
    stops after consecutive saturated layers. All use the SAME GN corrector.
    """
    config = config or Config()
    if config.arm not in ('H0', 'H1', 'H2', 'HA'): raise ValueError('Unknown arm')
    paid = PaidEvaluator(evaluate, config.local_rhs_cap, config.reserve_final_rhs, initial_state, config.local_lu_cap)
    anchor = initial_state if initial_state is not None else paid.evaluate(x0, True)
    if 'jac' not in anchor: anchor = paid.evaluate(x0, True)
    frame = weak_frame(anchor, config, basis_factory)
    initial_frame = frame
    rank, dimension = frame['W'].shape[1], len(x0)
    rng = np.random.default_rng(config.seed)
    count = config.layers*config.proposals_per_layer
    draws = rng.normal(size=(count, rank)) if random_draws is None else np.asarray(random_draws)
    if len(draws) < count or draws.shape[1] < rank: raise ValueError('Insufficient common random draws')
    cloud, rows = [anchor], []
    levels, qualified = [0], [qualification(anchor, anchor, config)[0]]
    current = anchor
    stale_layers, refreshes = 0, 0
    previous_failure, previous_secant = 0., 0.
    status = 'DEPTH_CAP'
    for layer in range(config.layers):
        # Required current J: paid only if the previous selected endpoint had F.
        try:
            if 'jac' not in current: current = paid.evaluate(current['x'], True)
        except LocalBudgetStop: status = 'LOCAL_BUDGET_STOP'; break
        candidate = weak_frame(current, config, basis_factory)
        drift = frame_drift(frame, candidate, current['dmaterial'], config.rms_denominator)
        do_refresh = (config.arm == 'H2' and layer == config.refresh_layer) or (
            config.arm == 'HA' and refreshes < config.adaptive_max_updates and layer > 0 and (
                drift >= config.adaptive_drift_threshold or
                previous_secant >= config.adaptive_secant_threshold or
                previous_failure >= config.adaptive_failure_fraction))
        if do_refresh:
            frame = candidate; refreshes += 1
        event = dict(phase='layer', layer=layer+1, action='REFRESH' if do_refresh else 'REUSE',
                     exploration_refreshes=refreshes, projector_drift=drift,
                     previous_secant=previous_secant, previous_failure=previous_failure,
                     cost_RHS=paid.rhs)
        rows.append(event)
        layer_states, layer_success, secants = [], 0, []
        before = [z for z, ok in zip(cloud, qualified) if ok]
        for proposal in range(config.proposals_per_layer):
            index = layer*config.proposals_per_layer+proposal
            layer_config = replace(config, predictor_rms=float(config.layer_radii[min(layer, len(config.layer_radii)-1)])) if config.layer_radii else config
            if config.arm == 'H0':
                center = np.asarray(prior_targets[index]) if prior_targets is not None else np.zeros(dimension)
                shift = np.asarray(residual_targets[index]) if residual_targets is not None else np.zeros_like(anchor['residual'])
                xp = current['x'].copy(); info = dict(status='RML_START')
                correct_frame = None
            else:
                direction = frame['W']@draws[index, :frame['W'].shape[1]]
                xp, info = predictor(current, direction, layer_config, decoder)
                center, shift = np.zeros(dimension), np.zeros_like(anchor['residual'])
                correct_frame = frame
            record = dict(phase='proposal', layer=layer+1, proposal=proposal,
                          candidate_index=len(cloud), predictor=info)
            if xp is None:
                record.update(status='PREDICTOR_UNAVAILABLE', qualified=False); rows.append(record); continue
            endpoint, correction, states = transverse_correct(paid, xp, correct_frame, config,
                decoder=decoder, prior_center=center, residual_shift=shift, reference=anchor)
            record['correction'] = correction
            for hr in correction['history']:
                if 'secant_relative_defect' in hr: secants.append(hr['secant_relative_defect'])
            if endpoint is None:
                record.update(status=correction['status'], qualified=False); rows.append(record)
                if correction['status'] == 'LOCAL_BUDGET_STOP': status='LOCAL_BUDGET_STOP'; break
                continue
            ok, gates = qualification(endpoint, anchor, config)
            layer_success += int(ok)
            cloud.append(endpoint); levels.append(layer+1); qualified.append(ok)
            layer_states.append((endpoint, ok))
            record.update(status='AVAILABLE', qualified=ok, gates=gates,
                          phi=float(endpoint['phi']), cost_RHS=paid.rhs,
                          actual_material_step=material_rms(endpoint['material'], current['material'], config.rms_denominator),
                          anchor_material_distance=material_rms(endpoint['material'], anchor['material'], config.rms_denominator))
            rows.append(record)
            if callback is not None: callback(record, endpoint)
        if status == 'LOCAL_BUDGET_STOP': break
        comparison = list(before)
        distinct = 0
        for z, ok in layer_states:
            if ok and all(material_rms(z['material'], b['material'], config.rms_denominator) > config.distinct_rms for b in comparison):
                distinct += 1
                comparison.append(z)
        stale_layers = stale_layers+1 if distinct == 0 else 0
        previous_failure = 1.-layer_success/max(config.proposals_per_layer, 1)
        previous_secant = max(secants, default=0.)
        rows.append(dict(phase='layer_end', layer=layer+1, action='CONTINUE',
                         distinct_qualified=int(distinct), qualified=layer_success,
                         secant_max=previous_secant, cost_RHS=paid.rhs))
        available = [z for z, ok in layer_states if ok]
        if available: current = min(available, key=lambda z: _objective(z, config, np.zeros(dimension), np.zeros_like(anchor['residual'])))
        if config.arm == 'HA' and stale_layers >= config.adaptive_saturation_layers:
            status = 'ADAPTIVE_SATURATION_STOP'; break
    # Final forward is explicitly available without consuming the reserved slot
    # when its exact state is still cached. This is a receipt, not a new gate.
    selectable = [i for i, ok in enumerate(qualified) if ok]
    selected = min(selectable, key=lambda i: _objective(cloud[i], config, np.zeros(dimension), np.zeros_like(anchor['residual']))) if selectable else 0
    try: paid.evaluate(cloud[selected]['x'], False, final=True)
    except LocalBudgetStop: status += '_FINAL_SLOT_UNAVAILABLE'
    result = dict(status=status, config=asdict(config), truth_read=False,
                  local_RHS=paid.rhs, local_LU=len(paid.calls), local_calls=paid.calls, refreshes=refreshes,
                  trace=rows, states=cloud, layer=levels, qualified=qualified,
                  initial_weak=initial_frame['W'], online_selected_index=selected)
    if output is not None: save_result(output, result)
    return result


def save_result(output, result):
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    cloud = result['states']
    arrays = dict(material=np.array([z['material'] for z in cloud]),
                  materials=np.array([z['material'] for z in cloud]),
                  Phi=np.array([z['phi'] for z in cloud]), phi=np.array([z['phi'] for z in cloud]),
                  legal=np.array([z.get('prior_legal', True) for z in cloud], bool),
                  qualified=np.array(result['qualified'], bool), layer=np.array(result['layer']),
                  params=np.array([z['x'] for z in cloud]), xs=np.array([z['x'] for z in cloud]),
                  prior_cost=np.array([z.get('reg_value', .5*z['x']@z['x']) for z in cloud]),
                  candidate_ids=np.arange(len(cloud)), initial_weak_frame=result['initial_weak'],
                  online_selected_index=np.array(result['online_selected_index']))
    if all('F' in z for z in cloud): arrays['Fs'] = np.array([z['F'] for z in cloud])
    if all('occupancy' in z for z in cloud): arrays['occupancyROI'] = np.array([z['occupancy'] for z in cloud])
    if all('delta' in z for z in cloud): arrays['delta'] = np.array([z['delta'] for z in cloud])
    np.savez_compressed(output/'CLOUD.npz', **arrays)
    with (output/'TRACE.jsonl').open('w') as f:
        for row in result['trace']: f.write(json.dumps(_plain(row), allow_nan=False)+'\n')
    summary = {k: v for k, v in result.items() if k not in ('states', 'trace', 'initial_weak')}
    (output/'SUMMARY.json').write_text(json.dumps(_plain(summary), indent=2, allow_nan=False)+'\n')


def effective_spectrum(J, D, *, denominator=576, fixed_map_derivative=None,
                       weak_functional=None):
    """Cached arithmetic audit; DH must be supplied, never guessed from J."""
    T, info = canonical_basis(D, denominator)
    M = material_matrix(D, denominator)
    K = np.asarray(J)@T
    _, singular, Vh = np.linalg.svd(K, full_matrices=True)
    eigen = np.r_[np.zeros(max(0, T.shape[1]-len(singular))), singular[::-1]**2]
    info.update(data_singular=singular, data_curvature_eigenvalues=eigen,
                effective_parameter_dimension=T.shape[1],
                numerical_first_order_data_nullity=int(T.shape[1]-np.linalg.matrix_rank(K)),
                data_rank_cutoff=float((singular[0] if len(singular) else 0.)*max(K.shape)*np.finfo(float).eps),
                physical_material_metric='||complex delta-material||/sqrt(denominator)',
                denominator=denominator, unit_eigenvalues_not_inferred_from_J=True)
    if fixed_map_derivative is not None:
        H = np.asarray(fixed_map_derivative)
        # Push material tangent forward, then pull back in the orthonormal image.
        effective_H = (M@T).T@M@H@T
        vals = np.linalg.eigvals(effective_H)
        info.update(fixed_map_eigenvalues_real=vals.real,
                    fixed_map_eigenvalues_imag=vals.imag,
                    fixed_map_spectral_radius=float(max(abs(vals), default=0)),
                    material_projection_defect=float(np.linalg.norm(M@H@T-(M@T)@effective_H)),
                    latent_decoder_gauge_component=float(np.linalg.norm(H@T-T@effective_H)))
    if weak_functional is not None:
        C = np.asarray(weak_functional)
        info['deliberately_frozen_weak_coordinates'] = int(np.linalg.matrix_rank(C@T))
    return info
