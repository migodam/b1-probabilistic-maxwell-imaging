"""Explicit finite-port formulas, not a new Maxwell solver or a prior decoder.

Physical callers must separately account for every solve. This module is used
only on synthetic matrices during the present zero-physics authorization.
"""
import numpy as np


def q_operator(d, h):
    d = np.asarray(d, complex)
    return np.linalg.solve(np.eye(len(d)) - d[:, None] * h, np.diag(d))


def galerkin(l, v):
    return v @ np.linalg.solve(v.conj().T @ l @ v, v.conj().T)


def exchange(l, v, c, u):
    rv = galerkin(l, v)
    pv = rv @ u
    dressed = c - rv @ l @ c
    schur = c.conj().T @ l @ dressed
    injection = c.conj().T @ (u - l @ pv)
    delta = dressed @ np.linalg.solve(schur, injection)
    return pv, delta, dressed, schur, injection


def finite_state_exchange(l_new, v, c, finite_rhs):
    """Exchange the finite-difference current solve directly, with no port map.

    The same exact Galerkin identity as ordinary current enrichment. Holding
    the proposed material fixed avoids changing its prior/proposal measure.
    """
    return exchange(l_new, v, c, finite_rhs)


def tail_bound(qv, ov, b, eps_h, eps_o):
    theta = float(np.linalg.norm(qv, 2) * eps_h)
    if theta >= 1:
        return dict(theta=theta, bound=None, status='SUFFICIENT_BOUND_UNAVAILABLE')
    value = (eps_o + np.linalg.norm(ov @ qv, 2) * eps_h) * np.linalg.norm(qv @ b) / (1 - theta)
    return dict(theta=theta, bound=float(value), status='FLOATING_NUMERICAL_BOUND_ESTIMATE')


def source_chart(t, b, h):
    e = b + h @ t
    if np.any(e == 0):
        raise ValueError('field zero: no material ratio chart')
    d = t / e
    m = np.eye(len(d)) - d[:, None] * h
    sign, ld = np.linalg.slogdet(m)
    if sign == 0:
        raise ValueError('singular feedback: no invertible source chart')
    derivative = m / e[:, None]
    return d, derivative, float(2 * ld - 2 * np.log(abs(e)).sum())


def row_fit(t, e, block_size=1):
    """One shared complex scalar per voxel, over all components and sources.

    A zero field remains explicitly undetermined, never epsilon-clipped.
    This least-squares fit is not material/prior validity or state validation.
    """
    shape = t.shape
    tt = t.reshape(-1, block_size * shape[1])
    ee = e.reshape(tt.shape)
    energy = np.sum(abs(ee) ** 2, axis=1)
    defined = energy > 0
    d = np.zeros(len(energy), complex)
    d[defined] = np.sum(tt[defined] * ee[defined].conj(), axis=1) / energy[defined]
    residual = tt - d[:, None] * ee
    return d, residual, defined
