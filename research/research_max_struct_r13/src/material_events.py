"""Finite shared-real-material events; no observations, labels, or truth input.

The exact threshold applies to the constitutive normal ray, not to an arbitrary
future DR trajectory. Scalar and cell-major vector fields share this geometry.
"""
import numpy as np


def norm_pair(x, q):
    return float(np.hypot(np.linalg.norm(x), np.linalg.norm(q)))


def finite_threshold(cx, cq, dx, dq, m, scale=6., lower=.02, upper=20., components=1):
    n = len(m)
    shape = (n, components * cx.shape[1])
    cx, cq, dx, dq = [np.asarray(a).reshape(shape) for a in (cx, cq, dx, dq)]
    lo, hi = np.broadcast_to(lower, (n,)), np.broadcast_to(upper, (n,))
    theta = np.arctan(np.asarray(m) / scale)
    sn, co = np.sin(theta)[:, None], np.cos(theta)[:, None]
    v, w = sn * cx + co * cq, co * dx - sn * dq
    a, b = np.sum(abs(v)**2, axis=1), np.sum(abs(w)**2, axis=1)
    gamma = np.real(np.sum(v.conj() * w, axis=1))
    amplitude_normal_defect = np.linalg.norm(sn * dx + co * dq, axis=1)
    lower_active = m <= lo + 1e-9
    upper_active = m >= hi - 1e-9
    interior = ~(lower_active | upper_active)
    t = np.full(n, np.inf)
    nonzero = b > np.finfo(float).tiny
    use = interior & nonzero
    t[use] = np.sqrt(a[use] / b[use])
    for active, other in ((lower_active, np.arctan(hi / scale)),
                          (upper_active, np.arctan(lo / scale))):
        use = active & nonzero
        g = gamma[use] / np.tan(other[use] - theta[use])
        rad = np.hypot(g, np.sqrt(a[use] * b[use]))
        # Avoid cancellation for a large positive boundary correlation.
        pos = g >= 0
        ans = np.empty_like(g)
        ans[pos] = a[use][pos] / (rad[pos] + g[pos])
        ans[~pos] = (rad[~pos] - g[~pos]) / b[use][~pos]
        t[use] = ans
    first = int(np.argmin(t))
    return dict(thresholds=t, first_cell=first, first_t=float(t[first]),
                a=a, b=b, gamma=gamma, interior=interior,
                lower_active=lower_active, upper_active=upper_active,
                amplitude_normal_defect=amplitude_normal_defect,
                interior_stationarity_defect=np.where(interior, abs(gamma), 0.))


def _angle_energy(A, B, C, lo, hi):
    """Maximal projected Gram energy on a closed material-angle interval."""
    theta = .5 * np.arctan2(2 * B, C - A)
    angles = np.stack([lo, hi, theta, theta + np.pi / 2,
                       theta - np.pi / 2, theta + np.pi, theta - np.pi], axis=1)
    valid = (angles >= lo[:, None]) & (angles <= hi[:, None])
    values = A[:, None] * np.sin(angles)**2 + C[:, None] * np.cos(angles)**2
    values += 2 * B[:, None] * np.sin(angles) * np.cos(angles)
    return np.max(np.where(valid, values, -np.inf), axis=1)


def angular_branch_margin(x, q, lower=.02, upper=20., scale=6., components=1,
                          perturbation_radius=0.):
    """Compare low/high halves of the legal material-angle interval.

    If |margin| > 2*rho, its sign is unchanged for every cell-state perturbation
    of Frobenius norm <= radius. rho=2||Z||r+r² bounds the Gram perturbation.
    This certifies branch preference in a declared STATE ball, not recovery or
    a future DR trajectory; the caller must justify that trajectory radius.
    """
    n = x.shape[0] // components
    x, q = x.reshape(n, -1), q.reshape(n, -1)
    A, C = np.sum(abs(x)**2, axis=1), np.sum(abs(q)**2, axis=1)
    B = np.real(np.sum(x.conj() * q, axis=1))
    lo = np.broadcast_to(np.arctan(np.asarray(lower) / scale), (n,))
    hi = np.broadcast_to(np.arctan(np.asarray(upper) / scale), (n,))
    mid = (lo + hi) / 2
    low = _angle_energy(A, B, C, lo, mid)
    high = _angle_energy(A, B, C, mid, hi)
    r = np.broadcast_to(perturbation_radius, (n,))
    rho = 2 * np.sqrt(A + C) * r + r*r
    margin = high - low
    return dict(margin=margin, uncertainty=2*rho,
                robust_high=margin > 2*rho, robust_low=margin < -2*rho)


def event_start(projector, ax, aq, mode, rng=None, overshoot=.001):
    """One globally scaled normal proposal from a point already in A.

    integer: first strictly post-tie integer; fractional: 0.1% beyond event;
    double: t=2 (does not inspect thresholds); random: matched state norm.
    A nonfinite threshold returns the unchanged point with explicit status.
    """
    cx, cq, m = projector.material(ax, aq)
    dx, dq = ax - cx, aq - cq
    tx, tq = projector.affine_ball(cx, cq)
    gap = norm_pair(dx, dq)
    error = norm_pair(tx - ax, tq - aq)
    th = finite_threshold(cx, cq, dx, dq, m, projector.scale,
                          projector.lower, projector.upper, projector.components)
    t = th['first_t']
    meta = dict(first_cell=th['first_cell'], first_t=t if np.isfinite(t) else None,
                gap=gap, pair_defect=error, relative_pair_defect=error/max(gap, 1e-300),
                amplitude_normal_defect=float(np.linalg.norm(th['amplitude_normal_defect'])),
                semantics='Approximate-pair finite state proposal; not exact DR fast-forward')
    if not np.isfinite(t) or gap == 0:
        return ax.copy(), aq.copy(), dict(meta, status='NO_FINITE_EVENT'), th
    if mode == 'double':
        target = 2.
    elif mode in ('integer','random_integer'):
        target = max(1., float(np.floor(t) + 1))
    else:
        target = max(1., (1 + overshoot) * t)
    if mode in ('random','random_integer'):
        rng = np.random.default_rng(13013) if rng is None else rng
        rx = rng.normal(size=ax.shape) + 1j*rng.normal(size=ax.shape)
        rq = rng.normal(size=aq.shape) + 1j*rng.normal(size=aq.shape)
        length = (target - 1) * gap
        factor = length / max(norm_pair(rx, rq), 1e-300)
        x, q = ax + factor*rx, aq + factor*rq
    else:
        x, q = cx + target*dx, cq + target*dq
    meta.update(target_t=target, displacement=norm_pair(x-ax, q-aq), status='PROPOSED')
    return x, q, meta, th
