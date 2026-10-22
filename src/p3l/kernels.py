"""Pure probability kernels on R^64 x [0,1]^2; no solver dependencies.

All targets return the *full* log density with respect to dz du, including
the N(0,I_64) prior. Return -inf outside any additional declared support.
Callables own their costs/caches; counts below count invocations, not RHS.
Freeze all kernel parameters before production. No adaptation is performed.
"""
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.linalg import solve_triangular

LogTarget = Callable[[np.ndarray], float]


def _q(value, property_box=False):
    q = np.array(value, dtype=float, copy=True)
    if q.shape != (66,) or not np.all(np.isfinite(q)):
        raise ValueError("state must have 66 finite real coordinates")
    if property_box and (np.any(q[64:] < 0) or np.any(q[64:] > 1)):
        raise ValueError("property coordinates must lie in [0,1]^2")
    return q


def _log(value):
    result = float(value)
    if np.isnan(result) or np.isposinf(result):
        raise ValueError("target returned NaN or +inf; only finite values/-inf allowed")
    return result


def _gradient(value):
    result = np.asarray(value, dtype=float)
    if result.shape not in ((64,), (66,)) or not np.all(np.isfinite(result)):
        raise ValueError("gradient must contain 64 or 66 finite entries")
    return result.copy()


@dataclass(frozen=True)
class ChainState:
    q: np.ndarray
    log_density: float
    gradient: np.ndarray | None = None

    def __post_init__(self):
        q = _q(self.q)
        q.setflags(write=False)
        lp = _log(self.log_density)
        if not np.isfinite(lp):
            raise ValueError("current state must have finite target density")
        object.__setattr__(self, "q", q)
        object.__setattr__(self, "log_density", lp)
        if self.gradient is not None:
            grad = _gradient(self.gradient)
            grad.setflags(write=False)
            object.__setattr__(self, "gradient", grad)


@dataclass(frozen=True)
class Transition:
    state: ChainState
    metadata: dict

    def __iter__(self):
        yield self.state
        yield self.metadata


def initialize(target: LogTarget, q, gradient=None):
    """Evaluate once, returning (state, invocation-count metadata)."""
    q = _q(q)
    lp = _log(target(q.copy()))
    g = None if gradient is None else _gradient(gradient(q.copy()))
    return ChainState(q, lp, g), dict(target_evaluations=1,
                                     gradient_evaluations=int(gradient is not None))


def reflect_unit_box(u):
    """Repeated reflection, not clipping; accepts arbitrarily many crossings."""
    u = np.asarray(u, dtype=float)
    if not np.all(np.isfinite(u)):
        raise ValueError("cannot reflect nonfinite coordinates")
    return 1.0 - np.abs(1.0 - np.remainder(u, 2.0))


def _scales(value):
    scales = np.broadcast_to(np.asarray(value, dtype=float), (2,)).copy()
    if not np.all(np.isfinite(scales)) or np.any(scales <= 0):
        raise ValueError("fixed diagonal property scales must be positive")
    scales.setflags(write=False)
    return scales


def _mh(current, proposed_q, proposed_lp, ratio, rng, name, gradients=0, grad=None):
    accepted = bool(np.log(max(rng.random(), np.finfo(float).tiny)) < min(0., ratio))
    state = ChainState(proposed_q, proposed_lp, grad) if accepted else current
    return Transition(state, dict(kernel=name, accepted=accepted,
        log_mh_ratio=float(ratio), target_evaluations=1,
        gradient_evaluations=gradients, rejected_proposals=int(not accepted),
        repeated_state=bool(np.array_equal(state.q, current.q)),
        status="EXACT_INVARIANT_KERNEL_NOT_A_MIXING_CERTIFICATE"))


class PCN:
    """Field pCN, optionally combined with symmetric diagonal property reflection.

    The Gaussian-reference proposal ratio is +(||z'||^2-||z||^2)/2.
    property_sd=None leaves properties fixed for a separate PropertyReflection.
    """
    def __init__(self, target: LogTarget, beta=.15, property_sd=None):
        if not np.isfinite(beta) or not 0 < beta <= 1:
            raise ValueError("beta must lie in (0,1]")
        self.target, self.beta = target, float(beta)
        self.property_sd = None if property_sd is None else _scales(property_sd)

    def step(self, current: ChainState, rng):
        _q(current.q, property_box=True)
        q = current.q.copy()
        q[:64] = np.sqrt(1-self.beta**2)*q[:64] + self.beta*rng.normal(size=64)
        if self.property_sd is not None:
            q[64:] = reflect_unit_box(q[64:] + self.property_sd*rng.normal(size=2))
        lp = _log(self.target(q.copy()))
        ratio = lp-current.log_density + .5*(q[:64]@q[:64]-current.q[:64]@current.q[:64])
        return _mh(current, q, lp, ratio, rng, "PCN")


class PropertyReflection:
    """A separate property block; only independent fixed Gaussian scales allowed.

    Coordinatewise reflection of a *correlated* Gaussian is generally not a
    symmetric box proposal, and is deliberately not exposed by this API.
    """
    def __init__(self, target: LogTarget, sd=.04):
        self.target, self.sd = target, _scales(sd)

    def step(self, current: ChainState, rng):
        _q(current.q, property_box=True)
        q = current.q.copy()
        q[64:] = reflect_unit_box(q[64:] + self.sd*rng.normal(size=2))
        lp = _log(self.target(q.copy()))
        return _mh(current, q, lp, lp-current.log_density, rng, "PROPERTY_REFLECTION")


class EllipticalSlice:
    """Gaussian64 conditional elliptical slice; properties remain fixed.

    A finite max_evaluations is a hard computational cap. Hitting it returns
    the repeated current state with SLICE_CAP_UNQUALIFIED. A finite cap can
    change the invariant law; capped runs are not labeled exact references.
    """
    def __init__(self, target: LogTarget, max_evaluations=None):
        if max_evaluations is not None and (not isinstance(max_evaluations, int)
                                            or max_evaluations < 1):
            raise ValueError("max_evaluations must be a positive integer or None")
        self.target, self.max_evaluations = target, max_evaluations

    def step(self, current: ChainState, rng):
        _q(current.q, property_box=True)
        z = current.q[:64]
        nu = rng.normal(size=64)
        log_slice = current.log_density + .5*(z@z) + np.log(max(rng.random(), np.finfo(float).tiny))
        theta = rng.uniform(0., 2*np.pi)
        lower, upper = theta-2*np.pi, theta
        evaluations = 0
        while True:
            q = current.q.copy()
            q[:64] = z*np.cos(theta) + nu*np.sin(theta)
            lp = _log(self.target(q.copy()))
            evaluations += 1
            if lp + .5*(q[:64]@q[:64]) >= log_slice:
                return Transition(ChainState(q, lp), dict(kernel="ELLIPTICAL_SLICE",
                    accepted=True, target_evaluations=evaluations, gradient_evaluations=0,
                    rejected_proposals=evaluations-1, repeated_state=np.array_equal(q, current.q),
                    bracket_width=upper-lower, status="EXACT_INVARIANT_KERNEL_NOT_A_MIXING_CERTIFICATE"))
            if self.max_evaluations is not None and evaluations >= self.max_evaluations:
                return Transition(current, dict(kernel="ELLIPTICAL_SLICE", accepted=False,
                    target_evaluations=evaluations, gradient_evaluations=0,
                    rejected_proposals=evaluations, repeated_state=True,
                    bracket_width=upper-lower, status="SLICE_CAP_UNQUALIFIED"))
            # The original point theta=0 stays in the bracket at every shrink.
            if theta < 0:
                lower = theta
            else:
                upper = theta
            theta = rng.uniform(lower, upper)


class PreconditionedMALA:
    """Fixed SPD MALA in Gaussian coordinates; never reflect a Langevin step.

    gradient(q) is the gradient of the full log target in Gaussian64 (or full
    66 coordinates). Optional basis (dimension,r) must be Euclidean orthonormal.
    q' = q + basis*(h/2*M*basis.T*grad + sqrt(h)*chol(M)*xi).
    Both forward/reverse drifts enter MH. Default dimension=64 keeps properties
    separate. dimension=66 is exclusively for an explicitly declared unbounded
    transformed coordinate target, with all transform Jacobians in the callable.
    ChainState itself only checks finite coordinates, to support this case;
    canonical pCN/reflection/slice additionally enforce the property box.
    """
    def __init__(self, target: LogTarget, gradient, step_size, preconditioner=None,
                 basis=None, dimension=64):
        if not np.isfinite(step_size) or step_size <= 0:
            raise ValueError("step_size must be positive")
        if dimension not in (64, 66):
            raise ValueError("dimension must be 64 (canonical) or 66 (transformed)")
        B = np.eye(dimension) if basis is None else np.array(basis, dtype=float, copy=True)
        if B.ndim != 2 or B.shape[0] != dimension or not 1 <= B.shape[1] <= dimension:
            raise ValueError("basis must have shape (dimension,r), 1 <= r <= dimension")
        if not np.all(np.isfinite(B)) or not np.allclose(B.T@B, np.eye(B.shape[1]), atol=1e-12, rtol=1e-12):
            raise ValueError("basis must be finite and orthonormal")
        r = B.shape[1]
        M = np.eye(r) if preconditioner is None else np.array(preconditioner, dtype=float, copy=True)
        if M.shape != (r,r) or not np.all(np.isfinite(M)) or not np.allclose(M, M.T, atol=1e-12, rtol=1e-12):
            raise ValueError("preconditioner must be a fixed finite symmetric SPD matrix")
        M = .5*(M+M.T)
        L = np.linalg.cholesky(M)
        for a in (B, M, L):
            a.setflags(write=False)
        self.target, self.gradient = target, gradient
        self.step_size, self.basis, self.M, self.L = float(step_size), B, M, L
        self.dimension = dimension

    def step(self, current: ChainState, rng):
        if self.dimension == 64:
            _q(current.q, property_box=True)
        ng = 0
        g = current.gradient
        if g is None:
            g = _gradient(self.gradient(current.q.copy()))
            ng += 1
        h, B, M, L, d = self.step_size, self.basis, self.M, self.L, self.dimension
        if g.size < d:
            raise ValueError("gradient is shorter than MALA dimension")
        forward_drift = .5*h*(M@(B.T@g[:d]))
        delta = forward_drift + np.sqrt(h)*(L@rng.normal(size=B.shape[1]))
        q = current.q.copy()
        q[:d] += B@delta
        lp = _log(self.target(q.copy()))
        if np.isneginf(lp):
            # Gradient need not exist outside hard support; reverse density is
            # irrelevant because the MH numerator is zero.
            out = _mh(current, q, lp, -np.inf, rng, "PRECONDITIONED_MALA", ng)
            return Transition(ChainState(current.q, current.log_density, g), out.metadata)
        gp = _gradient(self.gradient(q.copy()))
        ng += 1
        if gp.size < d:
            raise ValueError("gradient is shorter than MALA dimension")
        reverse_drift = .5*h*(M@(B.T@gp[:d]))
        forward_error = solve_triangular(L, delta-forward_drift, lower=True)
        reverse_error = solve_triangular(L, -delta-reverse_drift, lower=True)
        log_proposal_ratio = (forward_error@forward_error-reverse_error@reverse_error)/(2*h)
        ratio = lp-current.log_density+log_proposal_ratio
        out = _mh(current, q, lp, ratio, rng, "PRECONDITIONED_MALA", ng, gp)
        if not out.metadata["accepted"]:
            out = Transition(ChainState(current.q, current.log_density, g), out.metadata)
        out.metadata["log_proposal_ratio"] = float(log_proposal_ratio)
        return out


def run_chain(target, initial, kernels, draws, rng):
    """Fixed sweep schedule, recording each block transition (including repeats).

    Returns states of shape (draws * n_kernels,66), state log-densities,
    transition metadata, final cache, and initialization invocation counts.
    The initial state is separate from production draws. No thinning occurs.
    """
    if not isinstance(draws, int) or draws < 1 or not kernels:
        raise ValueError("positive draws and a nonempty fixed kernel schedule required")
    if isinstance(initial, ChainState):
        current, init = initial, dict(target_evaluations=0, gradient_evaluations=0)
    else:
        current, init = initialize(target, initial)
    states, logps, metadata = [], [], []
    for sweep in range(draws):
        for block, kernel in enumerate(kernels):
            transition = kernel.step(current, rng)
            current = transition.state
            states.append(current.q.copy())
            logps.append(current.log_density)
            metadata.append(dict(transition.metadata, sweep=sweep, block=block))
    return dict(states=np.asarray(states), log_density=np.asarray(logps),
                transitions=metadata, final_state=current, initialization=init)
