"""Geometry diagnostics and full-rank proposal algebra; no physical calls.

Finite differences/commutators only diagnose supplied fields. No leakage result
claims Lie full rank, noise-ball path certification or posterior mixing.
Covariance/mean rules must be frozen state functions in production MH; history
adaptation needs a separate justification or explicit augmented-state kernel.
"""
import numpy as np


def bracket_from_differentials(X, Y, DX, DY):
    """[X,Y]=DY X-DX Y, with real coordinate derivatives."""
    X, Y, DX, DY = map(lambda a: np.asarray(a, float), (X, Y, DX, DY))
    n = len(X)
    if X.shape != (n,) or Y.shape != (n,) or DX.shape != (n, n) or DY.shape != (n, n) or not all(np.isfinite(v).all() for v in (X, Y, DX, DY)):
        raise ValueError('finite real fields and Jacobians required')
    return DY@X-DX@Y


def material_bracket_diagnostic(bracket, D, weak_directions, *, rtol=1e-10):
    """Absolute material RMS and normal fraction in declared material metric.

    Supply realified D=[Re Dchi;Im Dchi]/sqrt(N), in 64 or full66
    coordinates explicitly; Q is computed on material images, not latent space.
    """
    D, B, W = map(lambda a: np.asarray(a, float), (D, bracket, weak_directions))
    mapped = D@B
    U, singular, _ = np.linalg.svd(D@W, full_matrices=False)
    keep = singular > (singular[0]*rtol if len(singular) else 0)
    Q = U[:, keep]
    normal = mapped-Q@(Q.T@mapped)
    total = float(np.linalg.norm(mapped))
    return dict(material_bracket_amplitude=total, normal_material_amplitude=float(np.linalg.norm(normal)),
        normal_fraction=float(np.linalg.norm(normal)/total) if total else 0.,
        weak_material_rank=int(keep.sum()), coordinate_dimension=D.shape[1],
        full_lie_rank_proven=False, continuous_path_verified=False)


def full_rank_covariance(dimension, *, weak=None, brackets=None, weak_weight=1., bracket_weight=1., floor=.01):
    """Positive floor preserves every latent direction; no pseudodeterminant."""
    if not isinstance(dimension, (int, np.integer)) or dimension <= 0 or not all(np.isfinite(v) for v in (floor, weak_weight, bracket_weight)) or floor <= 0 or min(weak_weight, bracket_weight) < 0:
        raise ValueError('positive dimension/floor and nonnegative weights required')
    C = floor*np.eye(dimension)
    for V, weight in ((weak, weak_weight), (brackets, bracket_weight)):
        if V is None:
            continue
        V = np.asarray(V, float)
        if V.ndim != 2 or V.shape[0] != dimension or not np.isfinite(V).all():
            raise ValueError('finite direction columns required')
        C += weight*(V@V.T)
    return C


def _cholesky(covariance):
    C = np.asarray(covariance, float)
    if C.ndim != 2 or C.shape[0] != C.shape[1] or not np.isfinite(C).all() or not np.allclose(C, C.T, rtol=0, atol=1e-13):
        raise ValueError('finite symmetric full-rank covariance required')
    return np.linalg.cholesky(C)


def gaussian_log_proposal(destination, mean, covariance):
    """Actual Gaussian Lebesgue density, including state-specific determinant."""
    destination, mean = np.asarray(destination, float), np.asarray(mean, float)
    L = _cholesky(covariance)
    if destination.shape != (len(L),) or mean.shape != destination.shape or not np.isfinite(destination).all() or not np.isfinite(mean).all():
        raise ValueError('finite matching Gaussian coordinates required')
    r = np.linalg.solve(L, destination-mean)
    return float(-.5*(len(r)*np.log(2*np.pi)+2*np.log(np.diag(L)).sum()+r@r))


def draw_gaussian(mean, covariance, rng):
    mean = np.asarray(mean, float)
    L = _cholesky(covariance)
    if mean.shape != (len(L),) or not np.isfinite(mean).all():
        raise ValueError('finite mean shape')
    return mean+L@rng.normal(size=len(mean))


def corrected_log_acceptance(x, candidate, *, logtarget_x, logtarget_candidate,
        mean_at_x, covariance_at_x, mean_at_candidate, covariance_at_candidate):
    """Forward/reverse full-rank proposal with caller-supplied Full target.

    No hidden endpoint optimization, clipping or rank-deficient Euler density.
    Means allow either random walk or an explicitly supplied MALA drift.
    """
    forward = gaussian_log_proposal(candidate, mean_at_x, covariance_at_x)
    reverse = gaussian_log_proposal(x, mean_at_candidate, covariance_at_candidate)
    if not np.isfinite(logtarget_x) or np.isnan(logtarget_candidate) or np.isposinf(logtarget_candidate):
        raise ValueError('finite current target and finite/-inf candidate target required')
    return float(min(0., logtarget_candidate-logtarget_x+reverse-forward))


def soft_mobility(J, *, tau=1., floor=.01):
    J = np.asarray(J, float)
    if J.ndim != 2 or not np.isfinite(J).all() or not np.isfinite(tau) or tau <= 0 or not 0 < floor <= 1:
        raise ValueError('finite J, tau>0 and 0<floor<=1 required')
    K = J.T@J
    R = tau*np.linalg.solve(K+tau*np.eye(K.shape[0]), np.eye(K.shape[0]))
    return floor*np.eye(K.shape[0])+(1-floor)*R
