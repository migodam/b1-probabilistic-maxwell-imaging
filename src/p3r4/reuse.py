"""Finite material feedback on caller-owned cached arrays; no backend or LU solve.

The caller alone acquires W = L0^-1 E_I and charges |I| physical RHS.
All solves here are reduced |I|-square algebra, reported separately from physical
RHS. An algebraically exact formula is not a floating-point interval certificate.
"""
from dataclasses import dataclass
import time
import numpy as np
from p3l.decision_bounds import (stability_from_saved_inverse,
    dual_corrected_prediction, same_uniform_decision)


class ReuseUnavailable(ValueError):
    """Caller may choose a charged Full fallback; this module never does so."""


def _array(value, ndim, name):
    value = np.asarray(value)
    if value.ndim != ndim or not np.isfinite(value).all():
        raise ValueError(f'{name}: finite {ndim}-dimensional array required')
    return value


def pack(z):
    z = _array(z, 2, 'complex receiver data')
    return np.concatenate((z.real.T, z.imag.T), axis=1).ravel()


def compose_jacobian(transfer, field, derivative, sigma):
    """J in the coordinates of supplied D; probit chain rule is caller-owned."""
    T, E, D = [_array(a, 2, n) for a, n in
        [(transfer, 'T'), (field, 'E'), (derivative, 'D')]]
    if T.shape[1] != E.shape[0] or D.shape[0] != E.shape[0] or not np.isfinite(sigma) or sigma <= 0:
        raise ValueError('incompatible Jacobian arrays or sigma')
    Z = np.einsum('qn,ns,nk->qsk', T, E, D, optimize=True)
    return np.column_stack([pack(Z[:, :, k])*np.sqrt(2)/sigma for k in range(D.shape[1])])


def jacobian_evolution(T0, T1, E0, E1, D0, D1, sigma):
    """Exact telescoping: propagation + total field + decoder derivative."""
    a = compose_jacobian(np.asarray(T1)-T0, E0, D0, sigma)
    b = compose_jacobian(T1, np.asarray(E1)-E0, D0, sigma)
    c = compose_jacobian(T1, E1, np.asarray(D1)-D0, sigma)
    return dict(propagation=a, total_field=b, decoder_derivative=c, total=a+b+c)


def changed_support(chi0, chi1):
    """Exact stored-value nonzeros: never silently truncate finite increments."""
    a, b = _array(chi0, 1, 'chi0'), _array(chi1, 1, 'chi1')
    if a.shape != b.shape:
        raise ValueError('material shape mismatch')
    return np.flatnonzero(b != a)


@dataclass
class FiniteMaterialCache:
    chi0: np.ndarray
    current0: np.ndarray
    field0: np.ndarray
    prediction0: np.ndarray
    G: np.ndarray
    S: np.ndarray
    incident: np.ndarray
    active: np.ndarray
    W: np.ndarray
    base_id: str
    setup_physical_rhs: int
    setup_owner: str

    def __post_init__(self):
        started = time.perf_counter()
        # Own immutable snapshots, so a caller changing a base cannot stale W.
        names = ('chi0', 'current0', 'field0', 'prediction0', 'G', 'S', 'incident', 'W')
        for name in names:
            a = np.array(getattr(self, name), dtype=complex, copy=True)
            if not np.isfinite(a).all():
                raise ValueError('nonfinite cache array: '+name)
            a.setflags(write=False)
            setattr(self, name, a)
        raw = np.asarray(self.active)
        if raw.ndim != 1 or not np.issubdtype(raw.dtype, np.integer):
            raise ValueError('active must contain integer indices')
        self.active = np.array(raw, dtype=int, copy=True)
        n = len(self.chi0)
        k = len(self.active)
        if self.chi0.shape != (n,) or self.G.shape != (n, n) or self.S.ndim != 2 or self.S.shape[1] != n:
            raise ValueError('cache geometry shape mismatch')
        if self.field0.ndim != 2 or self.current0.shape != self.field0.shape or self.incident.shape != self.field0.shape or self.field0.shape[0] != n or self.prediction0.shape != (self.S.shape[0], self.field0.shape[1]):
            raise ValueError('cache state shape mismatch')
        if len(set(self.active.tolist())) != k or np.any(self.active < 0) or np.any(self.active >= n) or self.W.shape != (n, k):
            raise ValueError('invalid active set or inverse-column layout')
        if not isinstance(self.setup_physical_rhs, (int, np.integer)) or self.setup_physical_rhs < k or not self.base_id or not self.setup_owner:
            raise ValueError('explicit base and paid inverse-column acquisition owner/cost required')
        self.active.setflags(write=False)
        self.K = self.G[self.active]@self.W
        self.O = self.S@self.W
        self.B = self.field0[self.active]
        self.GW = self.G@self.W
        for a in (self.K, self.O, self.B, self.GW):
            a.setflags(write=False)
        self.base_audit = self.audit()
        if any(v > 1e-10 for v in self.base_audit.values() if isinstance(v, float)):
            raise ReuseUnavailable('supplied base or inverse columns fail equation audit')
        self.setup_arithmetic_wall_seconds = time.perf_counter()-started

    def audit(self):
        """Pure matrix residual diagnostics, including supplied inverse columns."""
        L0 = np.eye(len(self.chi0))-self.chi0[:, None]*self.G
        selector = np.eye(len(self.chi0))[:, self.active]
        relative = lambda a, b: float(np.linalg.norm(a-b)/max(np.linalg.norm(b), 1e-30))
        return dict(current=relative(L0@self.current0, self.chi0[:, None]*self.incident),
            total_field=relative(self.incident+self.G@self.current0, self.field0),
            prediction=relative(self.S@self.current0, self.prediction0),
            inverse_columns=relative(L0@self.W, selector),
            qualification='EMPIRICAL_BACKWARD_RESIDUAL_NOT_INTERVAL_CERTIFIED')

    @property
    def bytes(self):
        return sum(getattr(self, n).nbytes for n in ('chi0', 'current0', 'field0', 'prediction0',
            'G', 'S', 'incident', 'active', 'W', 'K', 'O', 'B', 'GW'))

    def evaluate(self, chi, *, max_rank=None, max_condition=None, property_changed=False):
        """Reuse relative to the ORIGINAL base, including all cumulative changes.

        Property changes receive no exemption: any changed cell outside active
        rejects reuse. A returned state has no factorization and cannot be used
        as a new inverse base without acquiring consistent inverse actions.
        """
        started = time.perf_counter()
        chi = _array(chi, 1, 'chi')
        if chi.shape != self.chi0.shape:
            raise ValueError('material shape mismatch')
        changed = changed_support(self.chi0, chi)
        if not np.isin(changed, self.active).all():
            raise ReuseUnavailable('cumulative changed cells outside cached active set; extend paid columns or Full')
        if max_rank is not None and len(changed) > max_rank:
            raise ReuseUnavailable('actual diagonal material increment rank exceeds configured cap')
        d = chi[self.active]-self.chi0[self.active]
        small = np.eye(len(d), dtype=complex)-d[:, None]*self.K
        cond = float(np.linalg.cond(small)) if len(d) else 1.
        if not np.isfinite(cond) or (max_condition is not None and cond > max_condition):
            raise ReuseUnavailable('reduced solve condition monitor requests caller fallback; no stability certificate')
        try:
            c = np.linalg.solve(small, d[:, None]*self.B) if len(d) else np.zeros((0, self.field0.shape[1]), complex)
        except np.linalg.LinAlgError as error:
            raise ReuseUnavailable('reduced system singular') from error
        current = self.current0+self.W@c
        field = self.field0+self.GW@c
        prediction = self.prediction0+self.O@c
        # Residual is full-grid matrix arithmetic, not a physical solve.
        rhs = chi[:, None]*self.incident
        residual = rhs-(current-chi[:, None]*(self.G@current))
        return dict(current=current, field=field, prediction=prediction, coefficients=c,
            backward_residual=float(np.linalg.norm(residual)/max(np.linalg.norm(rhs), 1e-30)),
            material_increment_rank=len(changed), cached_active_size=len(d), property_changed=bool(property_changed),
            small_condition_monitor=cond, rebase_performed=False, cache_base_id=self.base_id,
            qualification='EXACT_FINITE_IDENTITY_FLOAT64_NOT_INTERVAL_CERTIFIED',
            profile=dict(new_Maxwell_states=0, new_physical_RHS=0, new_physical_LU=0,
                reduced_factorizations=int(bool(len(d))), reduced_solve_rhs_columns=self.field0.shape[1] if len(d) else 0,
                cache_bytes=self.bytes, historical_setup_RHS=self.setup_physical_rhs,
                historical_setup_owner=self.setup_owner,
                setup_arithmetic_wall_seconds=self.setup_arithmetic_wall_seconds,
                wall_seconds=time.perf_counter()-started))


def broyden_update(J, step, response_change):
    """Real good-Broyden secant baseline; no claim of derivative integrability."""
    J, h, df = np.asarray(J, float), np.asarray(step, float), np.asarray(response_change, float)
    if J.ndim != 2 or h.shape != (J.shape[1],) or df.shape != (J.shape[0],) or not all(np.isfinite(a).all() for a in (J, h, df)) or h@h <= 0:
        raise ValueError('nonzero finite secant with matching dimensions required')
    return J+np.outer(df-J@h, h)/(h@h)


def secant_decomposition(delta_prediction, J0, step, T0, E0, delta_chi, D0, sigma):
    """Vector remainder sum and cross term; norms must not be added as causes."""
    actual = pack(delta_prediction)*np.sqrt(2)/sigma
    linear = np.asarray(J0)@np.asarray(step)
    material = np.asarray(delta_chi)-np.asarray(D0)@np.asarray(step)
    injection = np.asarray(T0)@(np.asarray(E0)*material[:, None])
    decoder = pack(injection)*np.sqrt(2)/sigma
    maxwell = actual-linear-decoder
    return dict(linear=linear, decoder=decoder, maxwell=maxwell, remainder=decoder+maxwell,
        cross_inner_product=float(decoder@maxwell), actual=actual)


def gn_material_action(J, residual, x, D):
    """Full/frozenJ/Broyden hook in Gaussianized coordinates, H=I+J'J."""
    J, residual, x, D = map(np.asarray, (J, residual, x, D))
    delta = np.linalg.solve(np.eye(J.shape[1])+J.T@J, -(J.T@residual+x))
    return dict(delta=delta, material_action=D@delta)


def line_events(a, b, *, lower=-np.inf, upper=np.inf):
    """All s_i(t)=a_i+b_i*t saturation events, using standardized s."""
    a, b = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float))
    if not np.isfinite(a).all() or not np.isfinite(b).all() or lower >= upper:
        raise ValueError('finite line coefficients and ordered interval required')
    mask = b != 0
    events = np.r_[(-1-a[mask])/b[mask], (1-a[mask])/b[mask]]
    return np.unique(events[(events > lower)&(events < upper)])


def ellipse_events(a, b, c):
    """All roots a*cos(theta)+b*sin(theta)+c=+-1 on [0,2*pi)."""
    a, b, c = np.broadcast_arrays(*[np.asarray(v, float) for v in (a, b, c)])
    if not all(np.isfinite(v).all() for v in (a, b, c)):
        raise ValueError('finite ellipse coefficients required')
    radius = np.hypot(a, b)
    phase = np.arctan2(b, a)
    out = []
    for threshold in (-1., 1.):
        valid = (radius > 0)&(abs(threshold-c) <= radius)
        angle = np.arccos(np.clip((threshold-c)[valid]/radius[valid], -1, 1))
        out.extend(np.mod(phase[valid]+angle, 2*np.pi))
        out.extend(np.mod(phase[valid]-angle, 2*np.pi))
    return np.unique(np.asarray(out))


def compare_to_full(full_prediction, approximate_prediction, observation, sigma):
    """Audit supplied Full/reference arrays; obtaining them is caller-charged.

    This is observed error on this endpoint, never a candidate certificate.
    """
    full, approximate, y = [_array(v, 2, n) for v, n in
        ((full_prediction, 'Full'), (approximate_prediction, 'approximate'), (observation, 'y'))]
    if full.shape != approximate.shape or y.shape != full.shape or not np.isfinite(sigma) or sigma <= 0:
        raise ValueError('matching output shapes and positive sigma required')
    phi = lambda a: float(np.linalg.norm(a-y)**2/sigma**2)
    return dict(prediction_error=float(np.linalg.norm(approximate-full)),
        whitened_prediction_error=float(np.sqrt(2)*np.linalg.norm(approximate-full)/sigma),
        full_Phi=phi(full), approximate_Phi=phi(approximate),
        delta_Phi=phi(approximate)-phi(full),
        qualification='OBSERVED_ENDPOINT_ERROR_NOT_A_BOUND')
