"""Frozen property-spine coordinates; all 66 residual coordinates survive.

T(v, eta) = (v + beta(eta), eta) is a global polynomial shear. Its
Jacobian determinant is identically one. beta is a proposal coordinate,
not a claim about a posterior manifold or an exact conditional mode.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class PolynomialSpine:
    coefficients: np.ndarray
    center: np.ndarray
    scale: float = .5

    def __post_init__(self):
        c = np.array(self.coefficients, dtype=float, copy=True)
        center = np.array(self.center, dtype=float, copy=True)
        if c.shape != (64, 6) or center.shape != (2,):
            raise ValueError('spine requires coefficients 64x6 and center 2')
        if not np.isfinite(c).all() or not np.isfinite(center).all():
            raise ValueError('finite spine coefficients and center required')
        if not np.isfinite(self.scale) or self.scale <= 0:
            raise ValueError('positive finite scale required')
        c.setflags(write=False); center.setflags(write=False)
        object.__setattr__(self, 'coefficients', c)
        object.__setattr__(self, 'center', center)

    def basis(self, eta):
        eta = np.asarray(eta, dtype=float)
        if eta.shape != (2,): raise ValueError('eta must have two coordinates')
        a, b = (eta - self.center) / self.scale
        return np.array([1., a, b, a*a, a*b, b*b])

    def beta(self, eta):
        return self.coefficients @ self.basis(eta)

    def derivative(self, eta):
        a, b = (np.asarray(eta, dtype=float) - self.center) / self.scale
        db = np.array([[0., 0.], [1., 0.], [0., 1.],
                       [2*a, 0.], [b, a], [0., 2*b]]) / self.scale
        return self.coefficients @ db

    @staticmethod
    def _vector(a):
        a = np.array(a, dtype=float, copy=True)
        if a.shape != (66,): raise ValueError('full 66-vector required')
        return a

    def transform(self, r):
        x = self._vector(r); x[:64] += self.beta(x[64:]); return x

    def inverse(self, x):
        r = self._vector(x); r[:64] -= self.beta(r[64:]); return r

    def jacobian(self, r):
        r = self._vector(r)
        D = np.eye(66); D[:64, 64:] = self.derivative(r[64:]); return D

    def pullback_gradient(self, r, gradient_x):
        g = self._vector(gradient_x)
        g[64:] += self.derivative(np.asarray(r)[64:]).T @ g[:64]
        return g

    def log_prior(self, r):
        x = self.transform(r); return -.5 * float(x @ x)

    def prior_gradient(self, r):
        return self.pullback_gradient(r, -self.transform(r))

    @property
    def log_abs_det_jacobian(self):
        return 0.


class SpinePosterior:
    """Full original probit posterior composed with the volume-one shear.

    The underlying Posterior owns physical-x caches and pinning. No physical
    action is hidden here and no N(0,I) prior on r is substituted.
    """
    def __init__(self, posterior, spine):
        if not posterior.probit: raise ValueError('requires full probit posterior')
        self.base, self.spine = posterior, spine

    def point(self, r): return self.base.point(self.spine.transform(r))
    def pin(self, r): self.base.pin(self.spine.transform(r))

    def __call__(self, r):
        if not np.isfinite(r).all(): return -np.inf
        with np.errstate(over='ignore', invalid='ignore'):
            x = self.spine.transform(r)
        if not np.isfinite(x).all(): return -np.inf
        return self.base(x)

    def gradient(self, r):
        with np.errstate(over='ignore', invalid='ignore'):
            x = self.spine.transform(r)
        if not np.isfinite(x).all(): return np.full(66, np.nan)
        return self.spine.pullback_gradient(r, self.base.gradient(x))


def pulled_metric(target, point, physical_x, r, spine):
    """True newly paid full J, then arithmetic GN coordinate pullback."""
    from scipy.linalg import cho_factor, cho_solve
    from .core import probit_derivative
    J = target.jacobian(point) * probit_derivative(physical_x)[None, :]
    Hx = J.T @ J + np.eye(66)
    D = spine.jacobian(r); Hr = D.T @ Hx @ D
    C = cho_solve(cho_factor(Hr, lower=True), np.eye(66))
    return .5 * (C + C.T), Hx, Hr
