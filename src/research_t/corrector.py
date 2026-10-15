"""Bounded material Gauss--Newton corrector for the frozen A30 TM/VIE model.

Importing this file performs no physics and does not load a Maxwell backend.
The supplied model owns all state and RHS reservations.  A newly constructed
receiver transfer costs q adjoint RHS; subsequent material Jacobian products
at that exact state use the paid transfer and field.  Material coordinates are
raw susceptibility Re(chi), Im(chi), not mass-normalized current coordinates.

No truth, topology label, posterior teacher, free nuisance, TV penalty, current
ROM, or A17 reduction is used.  This is an ordinary material-space comparator.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from time import perf_counter
from typing import Any, Mapping

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import LinearConstraint, lsq_linear, minimize
from scipy.sparse.linalg import LinearOperator


@dataclass(frozen=True)
class CorrectorConfig:
    mode: str = "full"
    max_steps: int = 6
    max_full_attempts: int | None = None
    trust_radius: float = .3
    damping: float = 1e-3
    max_backtracks: int = 2
    tau: float = 23.785449004057032
    chart_rank: int = 64
    sigma_convention: str = "complex_rms"
    linear_max_iterations: int = 100
    linear_tolerance: float = 1e-7
    kkt_tolerance: float = 1e-6
    min_step_norm: float = 1e-10
    min_relative_decrease: float = 1e-8
    stop_when_compatible: bool = True
    compatibility_target_factor: float = 1.

    def __post_init__(self):
        if self.mode not in {"full", "multiscale", "multiscale_dct", "multiscale_haar"}:
            raise ValueError("UNKNOWN_MATERIAL_CORRECTOR_MODE")
        if self.trust_radius not in {.3, .6}:
            raise ValueError("PREDECLARED_TRUST_RADII_ARE_0.3_AND_0.6")
        if not 0 <= self.max_backtracks <= 2 or self.max_steps < 0:
            raise ValueError("INVALID_CORRECTOR_ITERATION_POLICY")
        if self.max_full_attempts is not None and self.max_full_attempts < 0:
            raise ValueError("INVALID_FULL_ENDPOINT_ATTEMPT_CAP")
        if not 0 < self.compatibility_target_factor <= 1:
            raise ValueError("INVALID_COMPATIBILITY_HEADROOM")
        if self.damping < 0 or self.tau <= 0 or self.linear_tolerance <= 0:
            raise ValueError("INVALID_CORRECTOR_TOLERANCE")
        if self.sigma_convention not in {"complex_rms", "per_real_component"}:
            raise ValueError("EXPLICIT_SIGMA_CONVENTION_REQUIRED")
        if self.chart_rank < 2 or self.chart_rank % 2:
            raise ValueError("EVEN_REAL_CHART_RANK_REQUIRED")


@dataclass
class JacobianCache:
    state: Any
    model_id: int
    jacobian: np.ndarray
    whitening: np.ndarray
    metadata: dict[str, Any]


@dataclass
class CorrectorResult:
    state: Any
    best_valid_state: Any | None
    attempted_states: list[Any]
    path: list[dict[str, Any]]
    diagnostics: dict[str, Any]
    jacobian_cache: JacobianCache | None = None


def pack(data):
    """A30 data order: source0 [Re(q),Im(q)], then source1, etc."""
    a = np.asarray(data, dtype=np.complex128)
    if a.ndim != 2 or not np.isfinite(a).all():
        raise ValueError("FINITE_Q_BY_S_DATA_REQUIRED")
    return np.concatenate((a.real.T, a.imag.T), axis=1).reshape(-1)


def material_real(chi):
    x = np.asarray(chi, dtype=np.complex128).reshape(-1)
    return np.r_[x.real, x.imag]


def material_complex(x):
    x = np.asarray(x, dtype=float).reshape(-1)
    if len(x) % 2:
        raise ValueError("EVEN_REAL_MATERIAL_COORDINATE_COUNT_REQUIRED")
    n = len(x) // 2
    return x[:n] + 1j * x[n:]


def whitening(sigma, count, convention="complex_rms"):
    sig = np.broadcast_to(np.asarray(sigma, dtype=float), (int(count),)).copy()
    if not np.isfinite(sig).all() or np.any(sig <= 0):
        raise ValueError("FINITE_POSITIVE_NOISE_SCALE_REQUIRED")
    return (np.sqrt(2.) if convention == "complex_rms" else 1.) / sig


def material_bounds(n, protected=None):
    """Absolute real-coordinate bounds, intersected with explicit geometry.

    ``protected`` is a mapping with ``lower`` and ``upper``, each length 2*n
    in Re-all/Im-all order.  Equality bounds are permitted (locked core pixels).
    Core AND guard semantics are determined by the caller and preserved here;
    this function never infers a cavity from background-relative contrast.
    """
    lo = np.r_[np.full(n, -.5), np.zeros(n)]
    hi = np.r_[np.full(n, 4.), np.full(n, .5)]
    if protected is not None:
        if not isinstance(protected, Mapping) or not {"lower", "upper"} <= set(protected):
            raise TypeError("PROTECTED_GEOMETRY_REQUIRES_EXPLICIT_LOWER_UPPER")
        lo = np.maximum(lo, np.broadcast_to(np.asarray(protected["lower"], float), (2*n,)))
        hi = np.minimum(hi, np.broadcast_to(np.asarray(protected["upper"], float), (2*n,)))
    if not np.isfinite(lo).all() or not np.isfinite(hi).all() or np.any(lo > hi):
        raise ValueError("EMPTY_OR_NONFINITE_PROTECTED_MATERIAL_DOMAIN")
    return lo, hi


def linearize(model, state, sigma, *, sigma_convention="complex_rms", cache=None):
    """Exact whitened real material J from one paid receiver-adjoint block.

    J(q,s,p)=T(q,p)*field(p,s); T=S L^-1.  Every source shares the
    same complex material direction.  Cache identity is state-object exact.
    """
    count = np.size(state.data) * 2
    white = whitening(sigma, count, sigma_convention)
    if cache is not None:
        if (cache.state is not state or cache.model_id != id(model)
                or not np.array_equal(cache.whitening, white)):
            raise ValueError("STALE_OR_FOREIGN_MATERIAL_JACOBIAN_CACHE")
        return cache
    tic = perf_counter()
    transfer = model.receiver_transfer(state)
    # (source,receiver,pixel), to preserve each source's Re/Im packing.
    raw = np.einsum("qp,ps->sqp", transfer, state.field)
    re = np.concatenate((raw.real, raw.imag), axis=1).reshape(count, -1)
    im = np.concatenate((-raw.imag, raw.real), axis=1).reshape(count, -1)
    jac = np.concatenate((re, im), axis=1) * white[:, None]
    jac.setflags(write=False)
    meta = dict(real_rows=int(jac.shape[0]), real_material_columns=int(jac.shape[1]),
                receiver_adjoint_RHS_if_transfer_uncached=int(transfer.shape[0]),
                tangent_full_solver_RHS=0,
                construction="T[q,pixel]*field[pixel,source]; shared chi; exact paid state",
                whitened=True, raw_susceptibility_coordinates=True,
                packing="source-major_Re_receivers_then_Im_receivers",
                wall_seconds=perf_counter()-tic)
    return JacobianCache(state, id(model), jac, white, meta)


def multiscale_basis(grid=24, real_rank=64, *, kind="dct"):
    """Fixed representation-only orthonormal chart, no target/truth input."""
    n = int(grid)**2
    rank = int(real_rank)//2
    if not 1 <= rank <= n:
        raise ValueError("MULTISCALE_RANK_OUT_OF_RANGE")
    if kind == "haar":
        # Reuse A30's exact 24->12->6->3 transform if already imported.
        from a30.representations import haar24_transform
        spatial = np.asarray(haar24_transform()[0]).T[:, :rank]
    else:
        yy, xx = np.meshgrid(np.arange(grid)+.5, np.arange(grid)+.5, indexing="ij")
        modes = sorted(((a, b) for a in range(grid) for b in range(grid)),
                       key=lambda ab: (sum(ab), ab))[:rank]
        spatial = np.array([(np.cos(np.pi*a*yy/grid)*np.cos(np.pi*b*xx/grid)).ravel()
                            for a, b in modes]).T
        spatial /= np.linalg.norm(spatial, axis=0)
    z = np.zeros_like(spatial)
    return np.block([[spatial, z], [z, spatial]])


def _projected_kkt(gradient, x, lo, hi, tolerance=1e-9):
    g = np.asarray(gradient, float).copy()
    at_lo, at_hi = x <= lo+tolerance, x >= hi-tolerance
    g[at_lo & (g > 0)] = 0.
    g[at_hi & (g < 0)] = 0.
    fixed = hi-lo <= 1e-13
    g[fixed] = 0.
    return dict(projected_gradient_inf=float(np.max(np.abs(g), initial=0.)),
                active_lower=int(at_lo.sum()), active_upper=int(at_hi.sum()),
                fixed_coordinates=int(fixed.sum()))


def _box_solve(jac, residual, low, high, damping, config):
    """Bounded regularized LS in full material coordinates; LSMR inner solve."""
    fixed = high-low <= 1e-13
    free = ~fixed
    delta = np.zeros(jac.shape[1])
    delta[fixed] = (low[fixed]+high[fixed])/2
    shift = residual + jac[:, fixed] @ delta[fixed]
    scale = max(float(np.linalg.norm(jac, "fro") / np.sqrt(min(jac.shape))), 1e-12)
    a, b = jac[:, free]/scale, -shift/scale
    mu = np.sqrt(damping)
    m, k = a.shape
    if not k:
        return delta, dict(success=True, status="ALL_COORDINATES_FIXED", iterations=0,
                           scale=scale, projected_gradient_inf=0.)
    if damping:
        op = LinearOperator((m+k, k), dtype=float,
                            matvec=lambda v: np.r_[a@v, mu*v],
                            rmatvec=lambda v: a.T@v[:m]+mu*v[m:])
        b = np.r_[b, np.zeros(k)]
    else:
        op = a
    fit = lsq_linear(op, b, bounds=(low[free], high[free]), method="trf",
                     lsq_solver="lsmr", tol=config.linear_tolerance,
                     lsmr_tol=config.linear_tolerance*.1,
                     max_iter=config.linear_max_iterations,
                     lsmr_maxiter=max(100, config.linear_max_iterations*2))
    delta[free] = np.minimum(high[free], np.maximum(low[free], fit.x))
    grad = jac.T @ (residual+jac@delta)/scale**2 + damping*delta
    return delta, dict(success=bool(fit.success), status=str(fit.message),
                       iterations=int(fit.nit), optimality=float(fit.optimality),
                       scale=scale, **_projected_kkt(grad, delta, low, high))


def _chart_solve(jac, residual, low, high, damping, config, basis):
    """Same physical row bounds in a fixed reduced material chart."""
    q = np.asarray(basis, dtype=float)
    if q.shape[0] != jac.shape[1] or not np.isfinite(q).all():
        raise ValueError("INVALID_REAL_MATERIAL_CHART")
    fixed = high-low <= 1e-13
    # All steps originate inside the domain, so fixed-row displacements are0.
    # Eliminate protected equalities rather than give redundant pixel equality
    # rows to a low-rank chart optimizer (which would spuriously fail SLSQP).
    if fixed.any():
        if np.max(np.abs(low[fixed]), initial=0.) > 1e-10:
            raise ValueError("CHART_FIXED_ROWS_REQUIRE_LEGAL_INITIAL_ENDPOINT")
        q = q @ null_space(q[fixed], rcond=1e-12)
    if q.shape[1] == 0:
        return np.zeros(jac.shape[1]), dict(success=True, status="PROTECTED_CHART_ZERO_RANK",
                                          iterations=0, real_chart_rank=0)
    a = jac @ q
    scale = max(float(np.linalg.norm(jac, "fro")/np.sqrt(min(jac.shape))), 1e-12)
    h = (a.T@a)/scale**2 + damping*(q.T@q)
    g = a.T@residual/scale**2
    active = np.linalg.norm(q, axis=1)>1e-14
    # Equalities (core locks) and inequalities are separate for SLSQP.
    eq = active & fixed
    ineq = active & ~fixed
    constraints = []
    # The remaining tiny fixed rows arise only from nullspace roundoff.
    if ineq.any():
        constraints.append(LinearConstraint(q[ineq], low[ineq], high[ineq]))
    fit = minimize(lambda z: float(.5*z@h@z+g@z), np.zeros(q.shape[1]),
                   jac=lambda z: h@z+g, method="SLSQP", constraints=constraints,
                   options=dict(maxiter=config.linear_max_iterations, ftol=config.linear_tolerance))
    delta = q@fit.x
    violation = max(0., float(np.max(low-delta)), float(np.max(delta-high)))
    # A numeric row violation is not repaired by silently changing the chart.
    # Contract a feasible ray only if roundoff is below the declared tolerance.
    if violation > 1e-8:
        return np.zeros_like(delta), dict(success=False, status="CHART_BOUND_VIOLATION",
                                         iterations=int(fit.nit), max_bound_violation=violation)
    delta = _contract_to_bounds(delta, low, high)
    return delta, dict(success=bool(fit.success), status=str(fit.message),
                       iterations=int(fit.nit), real_chart_rank=int(q.shape[1]),
                       max_bound_violation=violation,
                       chart_gradient_inf=float(np.max(np.abs(h@fit.x+g), initial=0.)),
                       KKT_status="SLSQP stationarity; no dual certificate inferred from success")


def _contract_to_bounds(delta, low, high):
    # Low<=0<=high because the current state obeys absolute constraints.
    delta = np.asarray(delta, float).copy()
    tiny = ((delta < low) | (delta > high)) & (np.abs(delta) < 1e-10)
    delta[tiny] = 0.
    pos, neg = delta > 0, delta < 0
    alpha = 1.
    if pos.any(): alpha = min(alpha, float(np.min(high[pos]/delta[pos])))
    if neg.any(): alpha = min(alpha, float(np.min(low[neg]/delta[neg])))
    return max(0., alpha)*(1.-1e-12)*delta


def bounded_compensation_score(model, state, y, sigma, config=None, protected=None, *, cache=None):
    """Primal/dual bounds for this finite edited state's LINEAR box chart.

    With no uniform nonlinear error bound these cannot certify full-Maxwell
    inclusion/exclusion.  The returned feasible delta may initialize correction.
    The dual bound is valid even if the bounded LS optimizer hits its cap.
    """
    config = config or CorrectorConfig()
    linear = linearize(model, state, sigma, sigma_convention=config.sigma_convention, cache=cache)
    x = material_real(state.chi)
    lo, hi = material_bounds(len(state.chi), protected)
    if np.any(x < lo-1e-12) or np.any(x > hi+1e-12):
        raise ValueError("EDITED_STATE_VIOLATES_PROTECTED_DOMAIN")
    low = np.maximum(lo-x, -config.trust_radius)
    high = np.minimum(hi-x, config.trust_radius)
    r = pack(state.data-y)*linear.whitening
    delta, diag = _box_solve(linear.jacobian, r, low, high, 0., config)
    v = r+linear.jacobian@delta
    upper = float(np.linalg.norm(v))
    dual = v/max(upper, 1e-300)
    c = linear.jacobian.T@dual
    lower = max(0., float(dual@r+np.sum(np.where(c >= 0, c*low, c*high))))
    # Dual norm is <=1 by construction; roundoff-only discrepancy is exposed.
    return dict(linear_primal_upper=upper, linear_dual_lower=lower,
                primal_dual_gap=max(0., upper-lower), linear_delta=material_complex(delta),
                raw_residual=float(np.linalg.norm(r)), diagnostics=diag,
                nonlinear_uniform_error_bound=None, nonlinear_feasibility="UNKNOWN",
                scope="FINITE_EDIT_CURRENT_TANGENT_BOX_CHART_ONLY",
                jacobian_cache=linear)


def _append(book, row):
    if hasattr(book, "append"):
        book.append(row)


def correct(model, state, y, sigma, book, config=None, protected=None, *, cache=None):
    """Return attempted full endpoints, best compatible endpoint and receipts.

    ``state`` must already be a paid matching endpoint.  The initial endpoint is
    never re-solved.  max_backtracks=2 means factors 1,1/2,1/4.  Each attempted
    factor invokes model.full_state and therefore reserves one state + s RHS,
    including rejected trials.  Each accepted new state requires a new paid
    receiver transfer before another GN step.  A budget exception terminates
    with already verified results intact and is recorded, not retried for free.
    """
    config = config or CorrectorConfig()
    if isinstance(config, Mapping): config = CorrectorConfig(**config)
    if model.book is not book:
        raise ValueError("CORRECTOR_AND_MODEL_MUST_SHARE_BUDGET")
    y = np.asarray(y, dtype=np.complex128)
    if y.shape != state.data.shape or not np.isfinite(y).all():
        raise ValueError("OBSERVATION_LAYOUT_MISMATCH")
    lo, hi = material_bounds(len(state.chi), protected)
    x = material_real(state.chi)
    if np.any(x < lo-1e-12) or np.any(x > hi+1e-12):
        raise ValueError("INITIAL_STATE_VIOLATES_PROTECTED_MATERIAL_BOUNDS")
    white = whitening(sigma, 2*y.size, config.sigma_convention)
    residual_norm = lambda s: float(np.linalg.norm(pack(s.data-y)*white))
    initial_norm = residual_norm(state)
    current_norm = initial_norm
    best_valid = state if initial_norm <= config.tau else None
    attempted, path = [], []
    stop = "MAX_STEPS"
    damping = config.damping
    last_cache = cache
    if config.mode == "full": basis = None
    else:
        kind = "haar" if config.mode == "multiscale_haar" else "dct"
        basis = multiscale_basis(getattr(model, "grid", 24), config.chart_rank, kind=kind)
    tic = perf_counter()
    for step in range(config.max_steps):
        if config.stop_when_compatible and current_norm <= config.tau*config.compatibility_target_factor:
            stop = "COMPATIBLE_FULL_ENDPOINT"; break
        try:
            if hasattr(book, "check"): book.check()
            last_cache = linearize(model, state, sigma, sigma_convention=config.sigma_convention,
                                   cache=last_cache if last_cache is not None and last_cache.state is state else None)
            x = material_real(state.chi)
            low = np.maximum(lo-x, -config.trust_radius)
            high = np.minimum(hi-x, config.trust_radius)
            r = pack(state.data-y)*white
            if basis is None:
                delta, linear_diag = _box_solve(last_cache.jacobian, r, low, high, damping, config)
            else:
                delta, linear_diag = _chart_solve(last_cache.jacobian, r, low, high, damping, config, basis)
            # Geometry boxes are constraints, not a post hoc image filter.
            delta = _contract_to_bounds(delta, low, high)
            row = dict(step=step, starting_residual=current_norm,
                       trust_radius=config.trust_radius, damping=damping,
                       step_norm=float(np.linalg.norm(delta)), linear=linear_diag,
                       real_material_dofs=len(x), mode=config.mode,
                       predicted_residual=float(np.linalg.norm(r+last_cache.jacobian@delta)),
                       line_search=[], accepted=False)
            if np.linalg.norm(delta) <= config.min_step_norm:
                row["stop"] = "NO_LEGAL_DESCENT_STEP"; path.append(row)
                stop = "NO_LEGAL_DESCENT_STEP"; break
            for backtrack in range(config.max_backtracks+1):
                if config.max_full_attempts is not None and len(attempted) >= config.max_full_attempts:
                    row["stop"] = "FULL_ENDPOINT_ATTEMPT_CAP"
                    break
                alpha = .5**backtrack
                z = x+alpha*delta
                # Clamp only machine-level addition error inside declared boxes.
                z = np.minimum(hi, np.maximum(lo, z))
                candidate = model.full_state(material_complex(z), purpose="ResearchT_material_GN_step_%d_backtrack_%d"%(step, backtrack))
                attempted.append(candidate)
                candidate_norm = residual_norm(candidate)
                receipt = dict(backtrack=backtrack, alpha=alpha, residual=candidate_norm,
                               compatible=candidate_norm<=config.tau,
                               backward_residual=float(candidate.residual),
                               material_box_max_violation=max(0.,float(np.max(lo-z)),float(np.max(z-hi))),
                               status="FULL_ENDPOINT")
                row["line_search"].append(receipt)
                if candidate_norm <= config.tau and (best_valid is None or candidate_norm<residual_norm(best_valid)):
                    best_valid = candidate
                if candidate_norm < current_norm*(1.-config.min_relative_decrease):
                    row["accepted"] = True
                    row["accepted_backtrack"] = backtrack
                    row["actual_reduction"] = current_norm**2-candidate_norm**2
                    linear_residual = r+alpha*last_cache.jacobian@delta
                    predicted = current_norm**2-float(linear_residual@linear_residual)
                    row["actual_to_predicted_reduction"] = row["actual_reduction"]/max(predicted, 1e-300)
                    state, current_norm = candidate, candidate_norm
                    last_cache = None
                    damping = max(config.damping*.1, damping*.5)
                    break
            path.append(row)
            _append(book, dict(event="MATERIAL_GN_STEP", **row))
            if not row["accepted"]:
                # A failed finite-step linearization is a result, not an excuse
                # for hidden extra LM refreshes or asymmetrical fallback.
                stop = row.get("stop", "LINE_SEARCH_FAILED"); break
        except Exception as error:
            failure = dict(step=step, error_type=type(error).__name__, error=str(error),
                           stop="BUDGET_OR_NUMERICAL_FAILURE", accepted=False)
            path.append(failure)
            _append(book, dict(event="MATERIAL_CORRECTOR_FAILURE", **failure))
            stop = "BUDGET_OR_NUMERICAL_FAILURE"; break
    diagnostic = dict(config=asdict(config), initial_residual=initial_norm,
                      final_residual=current_norm, tau=config.tau, stop_reason=stop,
                      attempted_full_endpoints=len(attempted),
                      accepted_steps=sum(bool(r.get("accepted")) for r in path),
                      best_valid_endpoint_found=best_valid is not None,
                      protected_domain_supplied=protected is not None,
                      wall_seconds=perf_counter()-tic,
                      regularization="LM step damping only; no TV or truth penalty",
                      theoretical_status="ordinary bounded nonlinear material GN, no global guarantee")
    _append(book, dict(event="MATERIAL_CORRECTOR_END", **diagnostic))
    return CorrectorResult(state, best_valid, attempted, path, diagnostic, last_cache)
