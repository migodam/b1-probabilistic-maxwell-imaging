"""Observation-derived full66 GN coordinate maps, without physical execution.

For cached whitened native Jq, Jx=Jq*dq/dx and U.T@U=I+Jx.T@Jx.
Positive-diagonal QR of [Jx;I] computes U without forming a poorly conditioned
normal matrix. This is an invertible coordinate map, not a posterior claim.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import qr, solve_triangular, lu_factor


def covariance_cholesky_from_factor(factor):
    """Fullrank L with L L.T=A A.T, by QR(A.T), no covariance normal matrix."""
    factor = np.asarray(factor)
    if factor.ndim != 3 or factor.shape[-1] != factor.shape[-2] or np.iscomplexobj(factor) or not np.isfinite(factor).all():
        raise ValueError('finite square batch covariance factors required')
    result = np.empty_like(factor, dtype=np.float64)
    for i, a in enumerate(factor):
        upper = qr(a.astype(np.float64).T, mode='r', check_finite=False)[0][:a.shape[0]]
        diagonal = np.diag(upper)
        if not np.isfinite(diagonal).all() or np.any(diagonal == 0):
            raise ValueError('singular Gaussian covariance factor; no rank repair')
        upper = np.where(diagonal >= 0, 1., -1.)[:, None]*upper
        result[i] = upper.T
    return result


def native_to_gaussianized_jacobian(context_x, native_J):
    x, jac = np.asarray(context_x), np.asarray(native_J)
    if x.ndim != 2 or x.shape[1] != 66 or jac.shape != (len(x), 512, 66):
        raise ValueError('context_x[B,66] and native whitened context_J[B,512,66] required')
    if np.iscomplexobj(x) or np.iscomplexobj(jac) or not np.isfinite(x).all() or not np.isfinite(jac).all():
        raise ValueError('finite real observation-derived anchor and native Jacobian required')
    jac = jac.astype(np.float64, copy=True)
    jac[:, :, 64:] *= np.exp(-.5*x[:, None, 64:].astype(np.float64)**2)/np.sqrt(2*np.pi)
    return jac


def fixed_random_rotation(seed=20261025113):
    """One fixed TRAIN-independent real orthogonal control, shared over y."""
    q, upper = qr(np.random.default_rng(seed).normal(size=(66, 66)), mode='economic')
    return q*np.where(np.diag(upper) >= 0, 1., -1.)[None]


class ConditionalCoordinateMap:
    """Linear map of per-y residuals; accepts [B,66] or [B,K,66]."""
    def __init__(self, mode, context_x=None, native_J=None, rotation=None, *, upper=None):
        if mode not in ('none', 'gn', 'random'):
            raise ValueError('coordinate preconditioner must be none/gn/random')
        self.mode = mode
        self.rotation = None
        self.upper = None
        if mode == 'none':
            self.batch = len(context_x) if context_x is not None else (len(upper) if upper is not None else 0)
            self.logabsdet = np.zeros(self.batch)
            return
        if upper is None:
            jac = native_to_gaussianized_jacobian(context_x, native_J)
            self.upper = np.empty((len(jac), 66, 66), np.float64)
            identity = np.eye(66)
            for i, j in enumerate(jac):
                if i and np.array_equal(j, jac[i-1]):
                    self.upper[i] = self.upper[i-1]
                    continue
                u = qr(np.concatenate((j, identity), axis=0), mode='r', check_finite=False)[0][:66]
                signs = np.where(np.diag(u) >= 0, 1., -1.)
                self.upper[i] = signs[:, None]*u
        else:
            self.upper = np.asarray(upper, np.float64)
        if self.upper.ndim != 3 or self.upper.shape[1:] != (66, 66) or not np.isfinite(self.upper).all():
            raise ValueError('finite full66 upper factors required')
        if not np.allclose(self.upper, np.triu(self.upper), rtol=0, atol=1e-12) or np.any(np.diagonal(self.upper, axis1=1, axis2=2) <= 0):
            raise ValueError('positive-diagonal triangular factors required; no H floor beyond I')
        self.batch = len(self.upper)
        self.logabsdet = np.log(np.diagonal(self.upper, axis1=1, axis2=2)).sum(1)
        if mode == 'random':
            rotation = np.asarray(rotation, np.float64)
            if rotation.shape != (66, 66) or not np.isfinite(rotation).all() or not np.allclose(rotation.T@rotation, np.eye(66), atol=1e-10, rtol=0):
                raise ValueError('fixed TRAIN-independent orthogonal rotation66 required')
            self.rotation = rotation
            # |det R|=1 algebraically; retain actual floating-point determinant.
            lu, _ = lu_factor(rotation, check_finite=False)
            logdet = np.log(abs(np.diag(lu))).sum()
            if not np.isfinite(logdet):
                raise ValueError('noninvertible fixed rotation')
            self.logabsdet = self.logabsdet+logdet

    def forward(self, delta):
        delta = np.asarray(delta)
        if delta.ndim not in (2, 3) or delta.shape[0] != self.batch or delta.shape[-1] != 66:
            raise ValueError('aligned residual[B,66] or[B,K,66] required')
        if self.mode == 'none':
            return delta
        if self.rotation is not None:
            delta = delta@self.rotation
        if delta.ndim == 2:
            return np.einsum('bij,bj->bi', self.upper, delta)
        return np.einsum('bij,bkj->bki', self.upper, delta)

    def inverse(self, residual):
        residual = np.asarray(residual)
        if residual.ndim not in (2, 3) or residual.shape[0] != self.batch or residual.shape[-1] != 66:
            raise ValueError('aligned transformed residual[B,66] or[B,K,66] required')
        if self.mode == 'none':
            return residual
        result = np.empty_like(residual, dtype=np.result_type(self.upper.dtype, residual.dtype))
        for i, u in enumerate(self.upper):
            rows = residual[i].reshape(-1, 66)
            value = solve_triangular(u, rows.T, lower=False, check_finite=False).T
            if self.rotation is not None:
                value = value@self.rotation.T
            result[i] = value.reshape(residual[i].shape)
        return result

    def matrix(self):
        if self.mode == 'none':
            return np.broadcast_to(np.eye(66), (self.batch, 66, 66))
        return self.upper if self.rotation is None else self.upper@self.rotation.T

    def transform_base(self, mean_delta, chol, covariance_coordinate='coordinate'):
        """Mean is raw displacement; chol units are explicitly declared.

        coordinate: isotropic/supplied fullrank Gaussian in r space. It induces
        GN/random-aligned raw covariance. raw: transform the supplied raw base.
        """
        mean = self.forward(mean_delta)
        if covariance_coordinate == 'coordinate' or self.mode == 'none':
            return mean, chol
        if covariance_coordinate != 'raw':
            raise ValueError('base covariance coordinate must be coordinate/raw')
        return mean, covariance_cholesky_from_factor(self.matrix()@chol)
