"""Small real linear/Gaussian audits; no PDE or true-posterior certificate.

Precision inputs are statistical model inputs. Algorithmic GN/LM damping must
never be supplied as a posterior prior. Numerical ranks are tolerance dependent.
"""
from __future__ import annotations

import numpy as np


def _real(value):
    raw = np.asarray(value)
    if np.iscomplexobj(raw):
        raise ValueError("Realify and whiten complex data before this audit")
    out = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError("Inputs must be finite")
    return out


def _sym(value):
    return (value + value.T) / 2


def _symmetric(value):
    value = _real(value)
    if value.ndim != 2 or value.shape[0] != value.shape[1]:
        raise ValueError("Expected square matrix")
    if not np.allclose(value, value.T, rtol=1e-10, atol=1e-12):
        raise ValueError("Expected symmetric statistical precision/metric")
    return _sym(value)


def _spd(value):
    value = _symmetric(value)
    if value.shape[0] and np.linalg.eigvalsh(value)[0] <= 0:
        raise ValueError("Expected positive definite matrix (handle gauge first)")
    return value


def svd_audit(matrix, rtol=None, atol=0.0):
    """Return rank, response/null bases, pseudoinverse and numerical threshold."""
    matrix = _real(matrix)
    if matrix.ndim != 2:
        raise ValueError("Expected matrix")
    if rtol is None:
        rtol = max(matrix.shape, default=0) * np.finfo(float).eps
    if rtol < 0 or atol < 0:
        raise ValueError("Rank tolerances must be nonnegative")
    u, s, vt = np.linalg.svd(matrix, full_matrices=True)
    threshold = max(atol, rtol * (s[0] if s.size else 0.0))
    rank = int(np.count_nonzero(s > threshold))
    v = vt.T
    inverse = (v[:, :rank] / s[:rank]) @ u[:, :rank].T
    return dict(rank=rank, threshold=float(threshold), singular_values=s,
                range_basis=u[:, :rank], parameter_basis=v[:, :rank],
                null_basis=v[:, rank:], pseudoinverse=inverse,
                rank_status="NUMERICAL_TOLERANCE_DEPENDENT")


def schur(a, b, prior_precision=None, *, rtol=None, atol=0.0):
    """Minimize ||Aw+Bu||²+u'prior_precision u without subtractive Schur."""
    a, b = _real(a), _real(b)
    if a.ndim != 2 or b.ndim != 2 or a.shape[0] != b.shape[0]:
        raise ValueError("A and B must have matching observation rows")
    audit = svd_audit(b, rtol, atol)
    c = audit["pseudoinverse"] @ a
    residual = a - b @ c
    data_information = _sym(residual.T @ residual)
    if prior_precision is None:
        prior = np.zeros((b.shape[1], b.shape[1]))
    else:
        prior = _symmetric(prior_precision)
        if prior.shape != (b.shape[1], b.shape[1]):
            raise ValueError("Prior precision has wrong dimension")
    vals, vecs = np.linalg.eigh(prior)
    if vals.size and vals[0] < -1e-12 * max(1.0, np.max(np.abs(vals))):
        raise ValueError("Prior precision must be positive semidefinite")
    root = np.sqrt(np.maximum(vals, 0))[:, None] * vecs.T
    assisted = bool(np.any(vals > 0))
    if assisted:
        aug_b = np.vstack((b, root))
        aug_a = np.vstack((a, np.zeros((b.shape[1], a.shape[1]))))
        aug = svd_audit(aug_b, rtol, atol)
        r = aug_a - aug_b @ (aug["pseudoinverse"] @ aug_a)
        information = _sym(r.T @ r)
    else:
        information = data_information
    return dict(information=information, data_only_information=data_information,
                prior_increment=_sym(information-data_information),
                prior_assisted=assisted, C=c, R=residual, rank_audit=audit)


def dynamic_information(a, b, d, q, *, require_full_rank=False,
                        rtol=None, atol=0.0):
    """D1 or D2; always eliminate ker(B) before motion compensation."""
    base = schur(a, b, rtol=rtol, atol=atol)
    b, d, q = _real(b), _real(d), _spd(q)
    if d.ndim != 2 or d.shape[1] != b.shape[1] or q.shape != (d.shape[0], d.shape[0]):
        raise ValueError("Motion dimensions do not match B/Q")
    audit = base["rank_audit"]
    if require_full_rank and audit["rank"] != b.shape[1]:
        raise ValueError("D1 requires full column rank; use D2")
    c, residual = base["C"], base["R"]
    v, n = audit["parameter_basis"], audit["null_basis"]
    qvals, qvecs = np.linalg.eigh(q)
    qi = (qvecs * qvals**-0.5) @ qvecs.T
    null_motion = svd_audit(qi @ d @ n, rtol, atol)
    un = null_motion["range_basis"]
    projection = np.eye(d.shape[0]) - un @ un.T
    dt = projection @ qi @ d @ v
    ct = projection @ qi @ d @ c
    # In retained SVD coordinates G_r is exactly diag(s_r²).
    gr = np.diag(audit["singular_values"][:audit["rank"]] ** 2)
    denominator = np.eye(d.shape[0]) + dt @ np.linalg.solve(gr, dt.T)
    motion_information = _sym(ct.T @ np.linalg.solve(denominator, ct))
    return dict(information=_sym(residual.T @ residual + motion_information),
                data_only_information=base["data_only_information"],
                motion_information=motion_information, prior_assisted=bool(np.any(motion_information)),
                C=c, R=residual, DC=d @ c, projected_conflict=ct,
                null_motion_projector=projection, rank_audit=audit,
                null_motion_rank_audit=null_motion,
                formula="D1" if audit["rank"] == b.shape[1] else "D2")


def dynamic_information_full_rank(a, b, d, q, **kwargs):
    return dynamic_information(a, b, d, q, require_full_rank=True, **kwargs)


def target_structure(a_blocks, b_blocks, *, structure="shared", d=None, q=None, **kwargs):
    """Shared static nuisance, per-time independent nuisance, or finite motion."""
    if not a_blocks or len(a_blocks) != len(b_blocks):
        raise ValueError("Nonempty matching block lists required")
    aa, bb = [_real(x) for x in a_blocks], [_real(x) for x in b_blocks]
    if any(x.ndim != 2 for x in aa+bb):
        raise ValueError("Per-time arrays must be matrices")
    if any(x.shape[0] != y.shape[0] for x,y in zip(aa,bb)):
        raise ValueError("Per-time observation rows mismatch")
    a = np.vstack(aa)
    if structure == "shared":
        return schur(a, np.vstack(bb), **kwargs)
    rows, cols = sum(x.shape[0] for x in bb), sum(x.shape[1] for x in bb)
    b = np.zeros((rows, cols))
    i = j = 0
    for ab, block in zip(aa, bb):
        if ab.shape[0] != block.shape[0]:
            raise ValueError("Per-time observation rows mismatch")
        b[i:i+block.shape[0], j:j+block.shape[1]] = block
        i += block.shape[0]
        j += block.shape[1]
    if structure == "independent":
        return schur(a, b, **kwargs)
    if structure == "motion" and d is not None and q is not None:
        return dynamic_information(a, b, d, q, **kwargs)
    raise ValueError("Expected shared, independent, or motion with D,Q")


def generalized_retention(a, b, *, material_metric=None, prior_precision=None,
                          rtol=None, atol=0.0):
    """Retention on range(K0); bare principal angles are data-only geometry."""
    a, b = _real(a), _real(b)
    metric = np.eye(a.shape[1]) if material_metric is None else _spd(material_metric)
    if metric.shape != (a.shape[1], a.shape[1]):
        raise ValueError("Material metric dimension mismatch")
    mv, mu = np.linalg.eigh(metric)
    inverse_root = (mu * mv**-0.5) @ mu.T
    am = a @ inverse_root
    audit_a, audit_b = svd_audit(am, rtol, atol), svd_audit(b, rtol, atol)
    ua = audit_a["range_basis"]
    cosines = np.linalg.svd(audit_b["range_basis"].T @ ua, compute_uv=False)
    cosines = np.r_[cosines, np.zeros(ua.shape[1]-cosines.size)]
    cosines = np.clip(cosines, 0, 1)
    angles = np.arccos(cosines)
    va = audit_a["parameter_basis"]
    scale = audit_a["singular_values"][:audit_a["rank"]]
    whitening = va / scale
    result = schur(am, b, prior_precision, rtol=rtol, atol=atol)
    retention = _sym(whitening.T @ result["information"] @ whitening)
    rho = np.linalg.eigvalsh(retention)
    return dict(retention_operator=retention, eigenvalues=rho,
                principal_angles=angles, bare_angle_retention=np.sin(angles)**2,
                angle_sum=float(np.sum(np.sin(angles)**2)),
                retention_trace=float(np.trace(retention)),
                absolute_information=result["information"], K0=am.T @ am,
                material_directions=inverse_root @ whitening,
                wall_rank_audit=audit_a, nuisance_rank_audit=audit_b,
                prior_assisted=result["prior_assisted"],
                zero_response_directions="UNDEFINED_RETENTION_EXCLUDED")


def conditional_innovation(h, wall_dimension, a, b):
    """D3 at a fixed linearization, with no unhandled new free states."""
    h, a, b = _symmetric(h), _real(a), _real(b)
    k = wall_dimension
    if not 0 < k < h.shape[0]:
        raise ValueError("Need nonempty wall and nuisance blocks")
    if np.linalg.eigvalsh(h)[0] < -1e-10 * max(1.0, np.linalg.norm(h, 2)):
        raise ValueError("Old joint H must be positive semidefinite")
    huu, huw = _spd(h[k:, k:]), h[k:, :k]
    cn = a - b @ np.linalg.solve(huu, huw)
    vn = np.eye(a.shape[0]) + b @ np.linalg.solve(huu, b.T)
    old = _sym(h[:k, :k] - h[:k, k:] @ np.linalg.solve(huu, huw))
    increment = _sym(cn.T @ np.linalg.solve(vn, cn))
    return dict(old_information=old, new_information=old+increment,
                increment=increment, conditional_jacobian=cn, conditional_noise=vn,
                increment_eigenvalues=np.linalg.eigvalsh(increment))


def local_gaussian_audit(h, g, c, freeze_indices, *, loss_metric=None,
                         approximate_duals=None):
    """D6/D7/D8 for an exact small Gaussian, never a true-posterior certificate.

    ``h`` is declared statistical precision, excluding algorithmic damping.
    Frozen coordinates are zero increments about the declared frozen point.
    Optional approximate_duals has shape (n_frozen, n_task).
    """
    h, g, c = _spd(h), _real(g), _real(c)
    n = h.shape[0]
    if g.shape != (n,) or c.ndim != 2 or c.shape[1] != n:
        raise ValueError("Gaussian/task dimensions mismatch")
    f = np.asarray(freeze_indices, dtype=int)
    if f.ndim != 1 or len(f) == 0 or len(np.unique(f)) != len(f) or np.any((f < 0) | (f >= n)):
        raise ValueError("Unique valid nonempty frozen indices required")
    u = np.setdiff1d(np.arange(n), f)
    huu, huf = h[np.ix_(u,u)], h[np.ix_(u,f)]
    hfu, hff = h[np.ix_(f,u)], h[np.ix_(f,f)]
    sf = _spd(hff - hfu @ np.linalg.solve(huu, huf))
    rf = g[f] - hfu @ np.linalg.solve(huu, g[u])
    tf = c[:, f] - c[:, u] @ np.linalg.solve(huu, huf)
    mean = np.linalg.solve(h, g)
    pinned_mean = np.zeros(n)
    pinned_mean[u] = np.linalg.solve(huu, g[u])
    covariance = np.linalg.solve(h, np.eye(n))
    pinned_covariance = np.zeros_like(h)
    pinned_covariance[np.ix_(u,u)] = np.linalg.solve(huu, np.eye(len(u)))
    m = np.eye(c.shape[0]) if loss_metric is None else _symmetric(loss_metric)
    if m.shape != (c.shape[0], c.shape[0]):
        raise ValueError("Loss metric dimension mismatch")
    vals, vecs = np.linalg.eigh(m)
    if vals.size and vals[0] < -1e-12:
        raise ValueError("Loss metric must be positive semidefinite")
    mroot = (vecs * np.sqrt(np.maximum(vals, 0))) @ vecs.T
    task_delta = tf @ np.linalg.solve(sf, rf)
    df = mroot @ tf
    exact_duals = np.linalg.solve(sf, df.T)
    duals = exact_duals if approximate_duals is None else _real(approximate_duals)
    if duals.shape != exact_duals.shape:
        raise ValueError("Approximate dual shape mismatch")
    xi = df.T - sf @ duals
    eigenvalues = np.linalg.eigvalsh(sf)
    lower_bound = float(eigenvalues[0])
    beta = np.linalg.norm(xi, axis=0) * np.linalg.norm(rf) / lower_bound
    radius = float(np.linalg.norm(beta))
    error_estimate = duals.T @ rf
    magnitude = np.linalg.norm(error_estimate)
    regret = float(task_delta @ m @ task_delta)
    return dict(scope="LOCAL_GAUSSIAN_EXACT_SMALL_SYSTEM",
                true_posterior_certificate=False, algorithmic_damping_included=False,
                mean=mean, pinned_mean=pinned_mean, covariance=covariance,
                pinned_covariance=pinned_covariance,
                cross_covariance=covariance[np.ix_(f,u)],
                task_mean=c @ mean, pinned_task_mean=c @ pinned_mean,
                task_covariance=c @ covariance @ c.T,
                pinned_task_covariance=c @ pinned_covariance @ c.T,
                omitted_task_covariance=tf @ np.linalg.solve(sf, tf.T),
                S=sf, r=rf, T=tf, objective_gap=float(0.5*rf @ np.linalg.solve(sf, rf)),
                task_regret=regret, task_error=task_delta,
                exact_duals=exact_duals, dual_residuals=xi, dual_beta=beta,
                dual_interval=(float(max(0,magnitude-radius)**2), float((magnitude+radius)**2)),
                lower_bound_source="DENSE_EIGENVALUE_FLOATING_POINT_SMALL_SYSTEM",
                lambda_min=lower_bound, lambda_max=float(eigenvalues[-1]),
                condition_number=float(eigenvalues[-1]/lower_bound),
                gaussian_qualification="Declared unconstrained Gaussian; no physical, mode, tail or PDE qualification")
