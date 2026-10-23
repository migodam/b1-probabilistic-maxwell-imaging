"""Generic real, whitened least-squares optimizers; no physical forward model.

``evaluate(x, jac=True|False) -> (residual, Jacobian_or_None)``.
Callbacks receive one metadata dictionary per event. They must not mutate inputs.
Adapter exceptions are logged; pass budget/authorization exception classes in
``fatal_exceptions`` to stop immediately instead of rejecting the trial.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import OptimizeResult, lsq_linear, minimize


class _CallLimit(RuntimeError):
    pass


class _EvaluationFailure(RuntimeError):
    pass


def _setup(x0, lower, upper, prior_center, prior_lambda):
    x = np.asarray(x0, dtype=float).copy()
    if x.ndim != 1 or not x.size or not np.all(np.isfinite(x)):
        raise ValueError('x0 must be a nonempty finite real vector')
    lo = np.broadcast_to(np.asarray(lower, dtype=float), x.shape).copy()
    hi = np.broadcast_to(np.asarray(upper, dtype=float), x.shape).copy()
    center = np.zeros_like(x) if prior_center is None else np.asarray(prior_center, dtype=float).copy()
    if (np.any(np.isnan(lo)) or np.any(np.isnan(hi)) or np.any(lo > hi)
            or np.any(x < lo) or np.any(x > hi)):
        raise ValueError('invalid bounds or infeasible x0')
    if center.shape != x.shape or not np.all(np.isfinite(center)):
        raise ValueError('invalid prior center')
    if not np.isfinite(prior_lambda) or prior_lambda < 0:
        raise ValueError('prior_lambda must be finite and nonnegative')
    return x, lo, hi, center


def _projected_gradient(x, g, lo, hi):
    # The unit projected-gradient mapping is zero exactly at box KKT points.
    return float(np.linalg.norm(x - np.clip(x - g, lo, hi), ord=np.inf))


class _Evaluator:
    def __init__(self, evaluate, center, lam, max_calls, callback, fatal_exceptions):
        self.evaluate, self.center, self.lam = evaluate, center, lam
        self.max_calls, self.callback = max_calls, callback
        self.fatal_exceptions = tuple(fatal_exceptions)
        self.nfev = self.njev = 0
        self.history = []

    def emit(self, event, **metadata):
        row = dict(event=event, **metadata)
        self.history.append(row)
        if self.callback is not None:
            self.callback(dict(row))

    def call(self, x, *, jac=True, **metadata):
        if self.nfev >= self.max_calls:
            self.emit('call_limit', x=x.tolist(), **metadata)
            raise _CallLimit('max_calls reached')
        self.nfev += 1
        self.njev += int(jac)
        call_id = self.nfev
        self.emit('evaluation_start', call_id=call_id, x=x.tolist(), jac=jac, **metadata)
        try:
            r, j = self.evaluate(x.copy(), jac=jac)
            if np.iscomplexobj(r) or (jac and np.iscomplexobj(j)):
                raise ValueError('evaluate must return real whitened arrays')
            r = np.asarray(r, dtype=float)
            j = np.asarray(j, dtype=float) if jac else None
            if r.ndim != 1 or (jac and j.shape != (r.size, x.size)):
                raise ValueError('residual/Jacobian shape mismatch')
            if not np.all(np.isfinite(r)) or (jac and not np.all(np.isfinite(j))):
                raise ValueError('nonfinite residual/Jacobian')
            delta = x - self.center
            f = float(.5 * np.dot(r, r) + .5 * self.lam * np.dot(delta, delta))
            g = j.T @ r + self.lam * delta if jac else None
            if not np.isfinite(f) or (jac and not np.all(np.isfinite(g))):
                raise ValueError('nonfinite objective/gradient')
        except Exception as exc:
            self.emit('evaluation_failed', call_id=call_id, x=x.tolist(),
                      error_type=type(exc).__name__, error=str(exc), **metadata)
            if isinstance(exc, self.fatal_exceptions):
                raise
            raise _EvaluationFailure(str(exc)) from exc
        self.emit('evaluation_complete', call_id=call_id, x=x.tolist(), fun=f, **metadata)
        return r.copy(), j.copy() if jac else None, f, g.copy() if jac else None, call_id


def _anchor(x, state, lo, hi):
    return dict(x=x.copy(), fun=state[2], call_id=state[4],
                projected_gradient=_projected_gradient(x, state[3], lo, hi))


def _result(x, state, ev, anchors, nit, success, message, lo, hi, **extra):
    return OptimizeResult(x=x.copy(), fun=state[2] if state else np.inf,
                          jac=state[3].copy() if state else None,
                          residual=state[0].copy() if state else None,
                          projected_gradient=(_projected_gradient(x, state[3], lo, hi)
                                              if state else np.inf),
                          nit=nit, nfev=ev.nfev, njev=ev.njev,
                          success=success, message=message, anchors=anchors,
                          history=ev.history, **extra)


def _linear_step(r, j, x, center, lam, lo, hi, radius, l2_radius):
    n = x.size
    free = lo < hi
    step = np.zeros(n)
    if not np.any(free):
        return step
    a, b = j[:, free], -r
    if lam:
        a = np.vstack((a, np.sqrt(lam) * np.eye(n)[:, free]))
        b = np.concatenate((b, -np.sqrt(lam) * (x - center)))
    lower, upper = np.maximum(lo - x, -radius)[free], np.minimum(hi - x, radius)[free]

    def solve(damping):
        aa, bb = a, b
        if damping:
            aa = np.vstack((a, np.sqrt(damping) * np.eye(np.count_nonzero(free))))
            bb = np.concatenate((b, np.zeros(np.count_nonzero(free))))
        fit = lsq_linear(aa, bb, bounds=(lower, upper), tol=1e-12,
                         max_iter=200, lsmr_tol='auto')
        if not fit.success or not np.all(np.isfinite(fit.x)):
            raise RuntimeError('bounded linear subproblem failed: ' + fit.message)
        return fit.x

    p = solve(0.)
    if l2_radius and np.linalg.norm(p) > radius:
        high = 1.
        for _ in range(100):
            candidate = solve(high)
            if np.linalg.norm(candidate) <= radius:
                break
            high *= 4.
        else:
            raise RuntimeError('could not bracket Euclidean radius damping')
        low = 0.
        for _ in range(50):
            mid = (low + high) / 2
            candidate_mid = solve(mid)
            if np.linalg.norm(candidate_mid) > radius:
                low = mid
            else:
                high, candidate = mid, candidate_mid
        p = candidate
    step[free] = p
    return step


def trust_region(evaluate, x0, lower, upper, prior_lambda=0., prior_center=None,
                 max_iter=30, radius=.3, callback=None, *, max_calls=200,
                 gtol=1e-6, gtol_abs=1e-10, min_radius=1e-10, max_radius=10.,
                 acceptance_ratio=.1, max_backtracks=12, l2_radius=True,
                 objective_floor=None, fatal_exceptions=()):
    """Full-J Gauss-Newton with bound-constrained trust-region subproblems.

    Trial points request residuals; provisionally accepted points request full FJ
    through ``evaluate`` before replacing the previous anchor.
    Rejected/failed points never replace the last accepted, fully evaluated anchor.
    Relative convergence: ||x-P_box(x-g)||_inf <= gtol_abs + gtol * pg_initial.
    ``objective_floor`` is an explicitly distinct floor stop, not convergence.
    ``max_calls`` includes failed evaluations; max_iter counts outer attempts.
    Callback exceptions and declared fatal adapter exceptions propagate.
    """
    x, lo, hi, center = _setup(x0, lower, upper, prior_center, prior_lambda)
    if (not 0 < min_radius <= radius <= max_radius or not np.isfinite(max_radius)
            or not 0 < acceptance_ratio < .25 or min(max_iter, max_calls, max_backtracks) < 1
            or min(gtol, gtol_abs) < 0):
        raise ValueError('invalid trust-region controls')
    ev = _Evaluator(evaluate, center, prior_lambda, max_calls, callback, fatal_exceptions)
    anchors, state, nit = [], None, 0
    try:
        state = ev.call(x, iteration=0, trial=0)
    except (_EvaluationFailure, _CallLimit) as exc:
        return _result(x, state, ev, anchors, nit, False, str(exc), lo, hi)
    anchors.append(_anchor(x, state, lo, hi))
    ev.emit('anchor', **{**anchors[-1], 'x': x.tolist()})
    initial_pg = anchors[0]['projected_gradient']
    success, message = False, 'max_iter reached'
    for iteration in range(1, max_iter + 1):
        pg = _projected_gradient(x, state[3], lo, hi)
        if pg <= gtol_abs + gtol * initial_pg:
            success, message = True, 'projected_gradient_converged'
            break
        if objective_floor is not None and state[2] <= objective_floor:
            message = 'objective_floor_reached (not gradient convergence)'
            break
        nit = iteration
        accepted = False
        for trial in range(1, max_backtracks + 1):
            try:
                p = _linear_step(state[0], state[1], x, center, prior_lambda,
                                 lo, hi, radius, l2_radius)
                pred = float(-state[3] @ p - .5 * np.dot(state[1] @ p, state[1] @ p)
                             - .5 * prior_lambda * np.dot(p, p))
            except Exception as exc:
                ev.emit('subproblem_failed', iteration=iteration, trial=trial,
                        radius=radius, error_type=type(exc).__name__, error=str(exc))
                pred = -np.inf
            if not np.isfinite(pred) or pred <= 0:
                ev.emit('trial_rejected', iteration=iteration, trial=trial,
                        radius=radius, reason='nonpositive_or_failed_model', predicted_reduction=None)
                radius *= .25
            else:
                candidate = x + p
                try:
                    next_state = ev.call(candidate, jac=False, iteration=iteration, trial=trial)
                except _CallLimit:
                    message = 'max_calls reached'
                    return _result(x, state, ev, anchors, nit, False, message, lo, hi, radius=radius)
                except _EvaluationFailure as exc:
                    ev.emit('trial_rejected', iteration=iteration, trial=trial,
                            x=candidate.tolist(), radius=radius, reason='evaluation_failed', error=str(exc))
                    radius *= .25
                else:
                    residual_call_id = next_state[4]
                    actual = state[2] - next_state[2]
                    ratio = actual / pred
                    accepted = bool(actual > 0 and ratio >= acceptance_ratio)
                    if accepted:
                        try:
                            full_state = ev.call(candidate, iteration=iteration, trial=trial,
                                                 phase='accepted_refresh')
                        except _CallLimit:
                            ev.emit('trial_rejected', iteration=iteration, trial=trial,
                                    call_id=next_state[4], x=candidate.tolist(), reason='refresh_call_limit')
                            return _result(x, state, ev, anchors, nit, False,
                                           'max_calls reached', lo, hi, radius=radius)
                        except _EvaluationFailure as exc:
                            ev.emit('trial_rejected', iteration=iteration, trial=trial,
                                    call_id=next_state[4], x=candidate.tolist(),
                                    reason='refresh_failed', error=str(exc))
                            radius *= .25
                            if radius < min_radius:
                                return _result(x, state, ev, anchors, nit, False,
                                               'radius_floor_reached (not gradient convergence)',
                                               lo, hi, radius=radius)
                            continue
                        # Acceptance is based on the full objective, even if an adapter
                        # returned inconsistent residual-only and full evaluations.
                        actual = state[2] - full_state[2]
                        ratio = actual / pred
                        accepted = bool(actual > 0 and ratio >= acceptance_ratio)
                        next_state = full_state
                    ev.emit('trial_accepted' if accepted else 'trial_rejected',
                            iteration=iteration, trial=trial, call_id=next_state[4],
                            residual_call_id=residual_call_id,
                            x=candidate.tolist(), radius=radius, actual_reduction=actual,
                            predicted_reduction=pred, ratio=ratio, reason='ratio_test')
                    if ratio < .25:
                        radius *= .25
                    elif ratio > .75 and (np.linalg.norm(p) if l2_radius else np.max(np.abs(p))) >= .8 * radius:
                        radius = min(max_radius, 2 * radius)
                    if accepted:
                        x, state = candidate, next_state
                        anchors.append(_anchor(x, state, lo, hi))
                        ev.emit('anchor', **{**anchors[-1], 'x': x.tolist()})
                        break
            if radius < min_radius:
                message = 'radius_floor_reached (not gradient convergence)'
                return _result(x, state, ev, anchors, nit, False, message, lo, hi, radius=radius)
        if not accepted:
            message = 'max_backtracks reached (not gradient convergence)'
            break
    # Include convergence reached on the final accepted iteration.
    if _projected_gradient(x, state[3], lo, hi) <= gtol_abs + gtol * initial_pg:
        success, message = True, 'projected_gradient_converged'
    return _result(x, state, ev, anchors, nit, success, message, lo, hi, radius=radius)


def lbfgsb(evaluate, x0, lower, upper, prior_lambda=0., prior_center=None,
           max_iter=30, callback=None, *, max_calls=200, gtol=1e-6,
           gtol_abs=1e-10, max_line_search=20, exact_cache=True,
           fatal_exceptions=()):
    """Alternative full-FJ L-BFGS-B; exact consecutive-x reuse is optional.

    Every evaluated line-search trial is ultimately labelled accepted/rejected.
    Adapter failure terminates this alternative at its last accepted anchor.
    scipy's ftol termination alone does not count as gradient convergence.
    """
    x, lo, hi, center = _setup(x0, lower, upper, prior_center, prior_lambda)
    if min(max_iter, max_calls, max_line_search) < 1 or min(gtol, gtol_abs) < 0:
        raise ValueError('invalid L-BFGS-B controls')
    ev = _Evaluator(evaluate, center, prior_lambda, max_calls, callback, fatal_exceptions)
    state, anchors, pending, cache = None, [], [], None
    nit = 0
    try:
        state = ev.call(x, phase='initial')
    except (_EvaluationFailure, _CallLimit) as exc:
        return _result(x, state, ev, anchors, nit, False, str(exc), lo, hi)
    cache = (x.copy(), state)
    anchors.append(_anchor(x, state, lo, hi))
    ev.emit('anchor', **{**anchors[-1], 'x': x.tolist()})
    threshold = gtol_abs + gtol * anchors[0]['projected_gradient']

    def objective(y):
        nonlocal cache
        if exact_cache and cache is not None and np.array_equal(y, cache[0]):
            ev.emit('cache_hit', call_id=cache[1][4], x=y.tolist())
            result = cache[1]
        else:
            result = ev.call(y, phase='line_search')
            cache = (y.copy(), result)
            pending.append((y.copy(), result))
        return result[2], result[3].copy()

    def classify(accepted_x=None):
        for y, data in pending:
            accepted = accepted_x is not None and np.array_equal(y, accepted_x)
            ev.emit('trial_accepted' if accepted else 'trial_rejected',
                    call_id=data[4], x=y.tolist(), fun=data[2], reason='scipy_line_search')
        pending.clear()

    def accepted(y):
        nonlocal x, state, nit
        if cache is None or not np.array_equal(y, cache[0]):
            objective(y)
        classify(y)
        x, state = y.copy(), cache[1]
        nit += 1
        anchors.append(_anchor(x, state, lo, hi))
        ev.emit('anchor', **{**anchors[-1], 'x': x.tolist()})

    raw_message = ''
    try:
        raw = minimize(objective, x, jac=True, method='L-BFGS-B', bounds=list(zip(lo, hi)),
                       callback=accepted, options=dict(maxiter=max_iter, maxfun=max_calls,
                       maxls=max_line_search, gtol=threshold, ftol=0.))
        raw_message = str(raw.message)
    except (_EvaluationFailure, _CallLimit) as exc:
        raw_message = str(exc)
    finally:
        classify()
    success = _projected_gradient(x, state[3], lo, hi) <= threshold
    message = 'projected_gradient_converged' if success else raw_message + ' (not gradient convergence)'
    return _result(x, state, ev, anchors, nit, success, message, lo, hi)


# Spelled-out alias for callers that prefer an explicit method name.
l_bfgs_b = lbfgsb
