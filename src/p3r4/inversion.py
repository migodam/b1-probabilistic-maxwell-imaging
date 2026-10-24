"""R4 reuse of unchanged R2 optimizer logic; independent core/ledger."""
"""Independent bounded TRF MAP optimization, with no truth-dependent inputs.

All physical calls go through the supplied Target. This module neither starts
jobs nor owns ledger reservations. The optimized objective is exactly
Phi(q;y) + ||q[:64]||^2/2 on Gaussian64 x Uniform[0,1]^2.
"""
from pathlib import Path
from collections import OrderedDict
import json
import time

import numpy as np
from scipy.optimize import least_squares, minimize
from scipy.stats import chi2

from .core import BudgetStop, Point, Target, plain, write


class _StateCap(RuntimeError):
    pass


def _legal(q):
    q = np.asarray(q, dtype=float)
    if q.shape != (66,) or not np.all(np.isfinite(q)):
        raise ValueError("expected 66 finite native coordinates")
    if np.any(q[64:] < 0) or np.any(q[64:] > 1):
        raise ValueError("uniform properties outside native [0,1] box")
    return q


def projected_gradient(q, gradient, tolerance=1e-8):
    """KKT gradient in the tangent cone of the native property box."""
    q, g = _legal(q), np.asarray(gradient, float).copy()
    if g.shape != (66,) or not np.all(np.isfinite(g)):
        raise ValueError("expected 66 finite objective derivatives")
    for j in (64, 65):
        if (q[j] <= tolerance and g[j] > 0) or (q[j] >= 1-tolerance and g[j] < 0):
            g[j] = 0.
    return g


def _coarse_indices(n):
    # Original MaterialPrior order: offset, then k-major 8x8 DCT, excluding 00.
    modes = [(k,l) for k in range(8) for l in range(8) if k or l]
    return np.array([0]+[i+1 for i,(k,l) in enumerate(modes) if k < n and l < n]
                    +[64,65], dtype=int)


def reconstruct(target: Target, y, parent, method="full_TRF", max_evaluations=96,
                *, initial_point: Point | None = None, polish=False,
                polish_max_evaluations=32, noise_confidence=.95):
    """Return (best objective Point, report), preserving failures and caps.

    The default start is [zeros64,.5,.5]. initial_point is an explicit online
    cache supplied by the caller; its provenance must be appropriate for y.
    `max_evaluations` caps *fresh forward states across every phase*, including
    optional polish. It does not grant a separate RHS/gradient/Jacobian budget.
    With polish=True, up to min(32, max_evaluations//3) fresh-state slots are
    reserved from that same cap. All optimization parameters stay fixed.
    """
    if method not in ("full_TRF", "coarse_TRF"):
        raise ValueError("method must be full_TRF or coarse_TRF")
    if not isinstance(max_evaluations, int) or max_evaluations < 1:
        raise ValueError("positive integer max_evaluations required")
    if not isinstance(polish_max_evaluations, int) or not 0 <= polish_max_evaluations <= 32:
        raise ValueError("polish_max_evaluations must be in [0,32]")
    if not 0 < noise_confidence < 1:
        raise ValueError("noise_confidence must be in (0,1)")
    y = np.asarray(y).copy()
    start = time.perf_counter()
    directory = Path(target.book.jobdir)/"inversion"/str(parent)/method
    directory.mkdir(parents=True, exist_ok=True)
    # A new invocation never overwrites an earlier optimizer's trace.
    trace_path = directory/"trace.jsonl"
    if trace_path.exists():
        raise FileExistsError(f"optimizer trace already exists: {trace_path}")
    trace = trace_path.open("x")
    # Keep the best endpoint separately; historical residual vectors are small
    # but a Point strongly owns a dense LU/matrix. Any genuinely evicted
    # revisit is a new paid state under the same unchanged fresh-state cap.
    cache, residuals = OrderedDict(), {}
    best = None
    best_value = np.inf
    fresh = 0
    jacobian_calls = 0
    source_gradient_calls = 0
    phases = []
    phase = "initialization"
    phase_cap = max_evaluations
    stage_fresh_start = 0
    failure = None
    stopped_by_budget = False

    def log(event, **fields):
        trace.write(json.dumps(plain(dict(event=event, parent=parent, method=method,
            phase=phase, wall_seconds=time.perf_counter()-start,
            fresh_forward_states=fresh, **fields)), allow_nan=False)+"\n")
        trace.flush()

    def residual(p):
        key = p.q.tobytes()
        if key not in residuals:
            r = np.r_[np.asarray(target.white(p.state.data-y), float).ravel(), p.q[:64]]
            if not np.all(np.isfinite(r)):
                raise FloatingPointError("nonfinite residual")
            residuals[key] = r
        return residuals[key]

    def remember(p):
        nonlocal best, best_value
        _legal(p.q)
        cache[p.q.tobytes()] = p
        r = residual(p)
        value = .5*float(r@r)
        if value < best_value:
            best, best_value = p, value
            write(directory/"best_point.json", dict(endpoint=p.key, q=p.q,
                objective=value, phi=target.phi(p,y), fresh_forward_states=fresh))
            log("BEST_OBJECTIVE", endpoint=p.key, objective=value)
        cache.move_to_end(p.q.tobytes())
        while len(cache)>4:
            old=next(iter(cache))
            if best is not None and old==best.q.tobytes():
                cache.move_to_end(old)
                old=next(iter(cache))
            del cache[old]
        return p

    def get(q, purpose):
        nonlocal fresh
        q = _legal(q).copy()
        key = q.tobytes()
        if best is not None and np.array_equal(q,best.q) and key not in cache:
            remember(best)
        hit = key in cache
        if not hit:
            if fresh >= max_evaluations or fresh-stage_fresh_start >= phase_cap:
                raise _StateCap("optimizer fresh-forward-state cap")
            # Increment before invoking Target: failures/attempts cannot hide.
            fresh += 1
            log("EVALUATION_ATTEMPT", q=q, purpose=purpose)
            point = target.evaluate(q, purpose=f"{parent}_{method}_{phase}_{purpose}")
            remember(point)
        point = cache[key]
        log("EVALUATION", endpoint=point.key, cache_hit=hit,
            objective=.5*float(residual(point)@residual(point)), q=q, purpose=purpose)
        return point

    prior_jacobian = np.zeros((64,66))
    prior_jacobian[:,:64] = np.eye(64)

    def full_jacobian(p):
        nonlocal jacobian_calls
        hit = p.J is not None
        if not hit:
            jacobian_calls += 1
        J = np.asarray(target.jacobian(p), float)
        if J.shape != (residual(p).size-64,66) or not np.all(np.isfinite(J)):
            raise ValueError("Target.jacobian shape/finite contract violated")
        log("JACOBIAN", endpoint=p.key, cache_hit=hit,
            receiver_jacobian_calls=jacobian_calls)
        return np.vstack((J, prior_jacobian))

    def kkt(p, allow_new_jacobian=True):
        if p.J is None and not allow_new_jacobian:
            return dict(status="NOT_AVAILABLE_WITHOUT_NEW_JACOBIAN", scaled_KKT=None,
                        raw_projected_gradient_inf=None)
        J = full_jacobian(p)
        g = J.T@residual(p)
        pg = projected_gradient(p.q, g)
        scale = np.sqrt(np.maximum(np.sum(J*J, axis=0),1.))
        return dict(status="COMPUTED", scaled_KKT=float(np.max(np.abs(pg)/scale)),
            raw_projected_gradient_inf=float(np.max(np.abs(pg))),
            objective_gradient=g, projected_gradient=pg,
            active_lower=(p.q[64:] <= 1e-8), active_upper=(p.q[64:] >= 1-1e-8),
            scaling="abs(projected g)/sqrt(max(diag(J_aug.T J_aug),1))",
            box_activity_tolerance=1e-8)

    try:
        q0 = np.r_[np.zeros(64), [.5,.5]]
        if initial_point is not None:
            q0 = _legal(initial_point.q).copy()
            remember(initial_point)
            log("EXPLICIT_ONLINE_POINT", endpoint=initial_point.key)
        reserved = min(polish_max_evaluations, max_evaluations//3) if polish else 0
        trf_cap = max_evaluations-reserved
        schedule = ([("coarse4", _coarse_indices(4), 8),
                     ("coarse6", _coarse_indices(6), 8)] if method == "coarse_TRF" else [])
        schedule.append(("full", np.arange(66), trf_cap))
        for name, indices, nominal_cap in schedule:
            remaining = trf_cap-fresh
            if remaining <= 0:
                break
            phase, phase_cap, stage_fresh_start = name, min(nominal_cap,remaining), fresh
            base = q0.copy() if best is None else best.q.copy()
            lower = np.where(indices >= 64, 0., -np.inf)
            upper = np.where(indices >= 64, 1., np.inf)
            def unpack(x):
                q = base.copy()
                q[indices] = x
                return q
            def fun(x):
                return residual(get(unpack(x), "fun"))
            def jac(x):
                return full_jacobian(get(unpack(x), "jac"))[:,indices]
            iterations = [0]
            def callback(intermediate_result):
                iterations[0] += 1
                log("ITERATION", iteration=iterations[0], q=unpack(intermediate_result.x),
                    cost=float(intermediate_result.cost), nfev=int(intermediate_result.nfev))
            log("PHASE_START", coordinates=indices, fresh_state_limit=phase_cap)
            try:
                result = least_squares(fun, base[indices], jac=jac,
                    bounds=(lower,upper), method="trf", tr_solver="exact",
                    x_scale="jac", ftol=1e-9, xtol=1e-9, gtol=1e-7,
                    max_nfev=phase_cap, callback=callback)
                record = dict(phase=phase, scipy_status=int(result.status),
                    success=bool(result.success), message=str(result.message),
                    nfev=int(result.nfev), njev=int(result.njev),
                    iterations=iterations[0], scipy_optimality=float(result.optimality),
                    terminal_cost=float(result.cost), fresh_states=fresh-stage_fresh_start)
            except _StateCap as exc:
                record = dict(phase=phase, scipy_status=None, success=False,
                    message=str(exc), iterations=iterations[0],
                    fresh_states=fresh-stage_fresh_start, status="STATE_CAP")
            phases.append(record)
            log("PHASE_END", **{k:v for k,v in record.items() if k != "phase"})
        if best is None:
            raise RuntimeError("optimizer produced no finite evaluated point")
        before_polish = kkt(best)
        if polish and reserved > 0 and before_polish["scaled_KKT"] > .001 and fresh < max_evaluations:
            phase, stage_fresh_start = "source_gradient_polish", fresh
            phase_cap = min(reserved, max_evaluations-fresh)
            polish_calls = [0]
            def objective_and_gradient(q):
                nonlocal source_gradient_calls
                if polish_calls[0] >= phase_cap:
                    raise _StateCap("source-gradient polish objective/gradient call cap")
                polish_calls[0] += 1
                p = get(q, "lbfgsb_fun_gradient")
                hit = p.gradient is not None
                if not hit:
                    source_gradient_calls += 1
                g = np.asarray(target.phi_gradient(p,y),float).copy()
                g[:64] += p.q[:64]
                if g.shape != (66,) or not np.all(np.isfinite(g)):
                    raise ValueError("source objective gradient contract violated")
                log("SOURCE_GRADIENT", endpoint=p.key, cache_hit=hit)
                return .5*float(residual(p)@residual(p)), g
            polish_iterations = [0]
            def polish_callback(q):
                polish_iterations[0] += 1
                log("ITERATION", iteration=polish_iterations[0], q=q)
            try:
                result = minimize(objective_and_gradient, best.q.copy(), method="L-BFGS-B",
                    jac=True, bounds=[(None,None)]*64+[(0.,1.)]*2,
                    callback=polish_callback,
                    options=dict(maxfun=phase_cap, maxiter=phase_cap, maxls=20,
                                 ftol=1e-12, gtol=1e-7))
                record = dict(phase=phase, scipy_status=int(result.status),
                    success=bool(result.success), message=str(result.message),
                    nfev=int(result.nfev), njev=int(result.njev),
                    iterations=polish_iterations[0], fresh_states=fresh-stage_fresh_start,
                    objective_gradient_invocations=polish_calls[0])
            except _StateCap as exc:
                record = dict(phase=phase, scipy_status=None, success=False,
                    message=str(exc), status="STATE_CAP", fresh_states=fresh-stage_fresh_start)
                record["objective_gradient_invocations"] = polish_calls[0]
            phases.append(record)
            log("PHASE_END", **{k:v for k,v in record.items() if k != "phase"})
    except Exception as exc:
        failure = dict(type=type(exc).__name__, message=str(exc))
        stopped_by_budget = isinstance(exc, BudgetStop)
        log("FAILURE", **failure)
        if best is None:
            trace.close()
            write(directory/"failure.json",failure)
            raise
    # Never request extra physics after any failed or externally capped call.
    try:
        final_kkt = kkt(best, allow_new_jacobian=failure is None)
    except Exception as exc:
        failure = dict(type=type(exc).__name__, message=str(exc))
        stopped_by_budget = stopped_by_budget or isinstance(exc, BudgetStop)
        log("FINAL_KKT_FAILURE", **failure)
        final_kkt = dict(status="FAILED", scaled_KKT=None, raw_projected_gradient_inf=None)
    r = residual(best)
    data_dimension = int(r.size-64)
    phi = target.phi(best,y)
    threshold = float(chi2.ppf(noise_confidence,data_dimension)/2)
    report = dict(parent=parent, method=method, target="Phi + ||z64||^2/2; native property box",
        endpoint=best.key, q=best.q, objective=best_value, phi=phi,
        fresh_forward_states=fresh, fresh_forward_state_cap=max_evaluations,
        receiver_J32_calls=jacobian_calls, fresh_source_gradient_calls=source_gradient_calls,
        kkt=final_kkt, phases=phases, failure=failure,
        status="BUDGET_STOP" if stopped_by_budget else ("FAILED_BEST_RETAINED" if failure else "OPTIMIZER_RETURNED"),
        noise_ball_qualified=bool(phi <= threshold), noise_ball_phi_threshold=threshold,
        noise_ball_confidence=noise_confidence, noise_ball_degrees_of_freedom=data_dimension,
        noise_ball_is_stopping_rule=False, polish_requested=bool(polish),
        wall_seconds=time.perf_counter()-start,
        cache_provenance="explicit caller-supplied online point" if initial_point is not None else "zero-field midpoint start",
        point_cache_policy="four recent Points with best pinned; evicted revisits count as fresh states within original cap",
        truth_used=False, posterior_reference_claim=False)
    target.persist_cache(best, directory/"best_cache.npz")
    write(directory/"result.json", report)
    log("FINISHED", endpoint=best.key, objective=best_value, status=report["status"])
    trace.close()
    return best, report


# Explicit verb form for callers that do not use the previous reconstruct name.
invert = reconstruct
