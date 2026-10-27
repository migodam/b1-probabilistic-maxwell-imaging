"""Small real/whitened stress diagnostics; no physical or posterior certificate.

The caller supplies A (wall), B (nuisance), and K (future task, in joint
[wall, nuisance] coordinates).  No whitening, material rescaling, finite
differences, reference reads, or truth reads occur here.  ``operator_error``
is an absolute spectral-norm error in the supplied A/B units, not a relative
tolerance.  A caller using one number must bound both blocks with that number.

Data-only geometry and an explicitly declared proper Gaussian prior are
reported separately.  The latter is a local linear Gaussian diagnostic and
must never be interpreted as bare sin-squared geometry or a finite posterior.
"""

from __future__ import annotations

import numpy as np


RESOLVED = "RESOLVED_AT_SUPPLIED_ERROR_RESOLUTION"
UNRESOLVED = "UNRESOLVED"
UNDEFINED_SUPPORT = "UNDEFINED_EMPTY_A_SUPPORT"
PROPER_GAUSSIAN = "PROPER_GAUSSIAN_PRIOR_LOCAL_LINEAR_DIAGNOSTIC"


def _readonly(value):
    result = np.array(value, copy=True)
    result.setflags(write=False)
    return result


def _real_matrix(value, name):
    raw = np.asarray(value)
    if np.iscomplexobj(raw) or raw.dtype.kind not in "fiu":
        raise ValueError(f"{name} requires caller-realified real coordinates")
    result = np.asarray(raw, dtype=float)
    if result.ndim != 2 or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a finite real matrix")
    return result


def _nonnegative_error(value, name):
    raw = np.asarray(value)
    if raw.ndim != 0 or raw.dtype.kind not in "fiu":
        raise ValueError(f"{name} must be one finite nonnegative absolute error")
    result = float(raw)
    if not np.isfinite(result) or result < 0:
        raise ValueError(f"{name} must be one finite nonnegative absolute error")
    return result


def _finite(value, name):
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} exceeds the finite floating-point range")
    return value


def _sym(value):
    return value / 2 + value.T / 2


def _cholesky(value, name):
    try:
        return np.linalg.cholesky(value)
    except np.linalg.LinAlgError as error:
        raise ValueError(f"{name} must be invertible positive definite; no jitter or pseudoinverse is used") from error


def _chol_solve(factor, rhs):
    if factor.shape[0] == 0:
        return np.zeros_like(rhs)
    return _finite(np.linalg.solve(factor.T, np.linalg.solve(factor, rhs)), "linear solve")


def _blocks(a, b):
    a, b = _real_matrix(a, "A"), _real_matrix(b, "B")
    if a.shape[0] == 0 or a.shape[1] == 0 or a.shape[0] != b.shape[0]:
        raise ValueError("A needs nonempty observation/wall axes and B must have matching observation rows")
    return a, b


def svd_rank_audit(matrix, *, operator_error):
    """Threshold singular values by the supplied absolute operator error only.

    ``rank`` counts s > operator_error.  The error envelopes [max(0,s-e),s+e]
    give a lower/upper *effective rank at that same threshold*.  If they differ,
    the rank status is UNRESOLVED.  These are conditional error-resolution
    diagnostics, not a claim about the exact algebraic rank of a true operator.
    No machine-epsilon rank threshold, pseudoinverse, or jitter is introduced.
    """
    matrix = _real_matrix(matrix, "matrix")
    error = _nonnegative_error(operator_error, "operator_error")
    u, singular, vt = np.linalg.svd(matrix, full_matrices=True)
    _finite(singular, "singular values")
    lower = np.maximum(0.0, singular - error)
    with np.errstate(over="ignore"):
        upper = _finite(singular + error, "singular-value error envelope")
    rank = int(np.count_nonzero(singular > error))
    rank_lower = int(np.count_nonzero(lower > error))
    rank_upper = int(np.count_nonzero(upper > error))
    return dict(
        rank=rank, rank_lower=rank_lower, rank_upper=rank_upper,
        status=RESOLVED if rank_lower == rank_upper else UNRESOLVED,
        threshold=error, operator_error=error, threshold_comparison="strict s > operator_error",
        singular_values=_readonly(singular),
        singular_value_lower=_readonly(lower), singular_value_upper=_readonly(upper),
        range_basis=_readonly(u[:, :rank]), parameter_basis=_readonly(vt[:rank].T),
        excluded_parameter_basis=_readonly(vt[rank:].T),
        scope="effective rank at caller-supplied error resolution; not true algebraic rank",
        floating_point_roundoff_certified=False,
    )


def data_only_diagnostics(a, b, *, operator_error):
    """Report absolute information, projected Schur, and A-support geometry.

    Absolute matrices retain all supplied A columns.  The data-only Schur is
    R.T R, R=(I-U_B U_B.T)A, with U_B selected by the error-threshold SVD.
    Generalized retention uses A's retained response support and equals the
    spectrum of (I-U_B U_B.T) restricted to that orthonormal support.  Directions
    outside A's effective support have undefined retention, not zero retention.
    Values remain available when rank is UNRESOLVED, explicitly conditional on
    the recorded threshold ranks.  Projector/covariance error is not certified.
    """
    a, b = _blocks(a, b)
    error = _nonnegative_error(operator_error, "operator_error")
    aa = svd_rank_audit(a, operator_error=error)
    bb = svd_rank_audit(b, operator_error=error)
    ua, ub = aa["range_basis"], bb["range_basis"]
    residual = a - ub @ (ub.T @ a)
    with np.errstate(over="ignore", invalid="ignore"):
        absolute = _finite(_sym(a.T @ a), "wall absolute information")
        schur = _finite(_sym(residual.T @ residual), "data-only Schur information")
    status = UNRESOLVED if UNRESOLVED in (aa["status"], bb["status"]) else RESOLVED
    rank = aa["rank"]
    if rank:
        supported_residual = ua - ub @ (ub.T @ ua)
        retention = _sym(supported_residual.T @ supported_residual)
        cosines_raw = np.linalg.svd(ub.T @ ua, compute_uv=False)
        cosines_raw = np.r_[cosines_raw, np.zeros(rank - cosines_raw.size)]
        # Clipping protects arccos only. Raw cosines and the correction survive.
        cosines = np.clip(cosines_raw, 0.0, 1.0)
        angles = np.arccos(cosines)
        sin_squared = np.sin(angles) ** 2
        retention_eigenvalues = np.linalg.eigvalsh(retention)
        angle_radians_sum = float(angles.sum())
        sin_squared_sum = float(sin_squared.sum())
        retention_trace = float(np.trace(retention))
        geometry_status = status
        clipping = float(np.max(np.abs(cosines - cosines_raw)))
    else:
        retention = np.empty((0, 0))
        retention_eigenvalues = cosines_raw = cosines = angles = sin_squared = np.empty(0)
        angle_radians_sum = sin_squared_sum = retention_trace = None
        geometry_status = UNDEFINED_SUPPORT
        clipping = 0.0
    return dict(
        label="DATA_ONLY_ERROR_THRESHOLD_GEOMETRY", status=status,
        operator_error=error, wall_rank_audit=aa, nuisance_rank_audit=bb,
        wall_absolute_information=_readonly(absolute),
        wall_absolute_eigenvalues=_readonly(np.linalg.eigvalsh(absolute)),
        wall_data_only_schur_information=_readonly(schur),
        wall_data_only_schur_eigenvalues=_readonly(np.linalg.eigvalsh(schur)),
        projected_wall_jacobian=_readonly(residual),
        generalized_retention_operator=_readonly(retention),
        generalized_retention_eigenvalues=_readonly(retention_eigenvalues),
        generalized_retention_status=geometry_status,
        retained_wall_dimension=rank, excluded_wall_dimension=a.shape[1] - rank,
        principal_cosines_raw=_readonly(cosines_raw), principal_cosines=_readonly(cosines),
        principal_angles_radians=_readonly(angles), principal_angle_radians_sum=angle_radians_sum,
        bare_sin_squared=_readonly(sin_squared), bare_sin_squared_sum=sin_squared_sum,
        retention_trace=retention_trace, arccos_clipping_max_abs=clipping,
        principal_angle_convention="one angle per retained A response direction; pad pi/2 when rank(A)>rank(B)",
        excluded_direction_status="UNDEFINED_RETENTION_EXCLUDED",
        numerical_scope="conditional on threshold ranks; operator error retained, no Schur/projector error envelope",
        prior_assisted=False, physical_certificate=False, floating_point_roundoff_certified=False,
    )


def proper_gaussian_diagnostics(a, b, k, *, prior_precision):
    """Use an explicit SPD joint prior in [wall, nuisance] coordinates.

    No singular prior, improper precision, pseudoinverse, or added damping is
    accepted.  K has one row per supplied task coordinate.  K's units/weighting
    are caller-defined; a trace is a sum of variances in those coordinates.
    Error in A/B/K is not propagated into a covariance certificate.
    """
    a, b = _blocks(a, b)
    k = _real_matrix(k, "K")
    nw, nu = a.shape[1], b.shape[1]
    n = nw + nu
    if k.shape[0] == 0 or k.shape[1] != n:
        raise ValueError("K needs nonempty task rows and joint [wall, nuisance] columns")
    prior = _real_matrix(prior_precision, "prior_precision")
    if prior.shape != (n, n) or not np.array_equal(prior, prior.T):
        raise ValueError("prior_precision must be explicitly symmetric with joint [wall, nuisance] shape")
    _cholesky(prior, "proper prior_precision")
    j = np.hstack((a, b))
    with np.errstate(over="ignore", invalid="ignore"):
        precision = _finite(_sym(j.T @ j + prior), "proper Gaussian joint precision")
    factor = _cholesky(precision, "proper Gaussian joint precision")
    covariance = _sym(_chol_solve(factor, np.eye(n)))
    # Nuisance-first Cholesky obtains a Schur factor without subtractive loss.
    order = np.r_[np.arange(nw, n), np.arange(nw)]
    permuted_factor = _cholesky(precision[np.ix_(order, order)], "nuisance-first joint precision")
    schur_factor = permuted_factor[nu:, nu:]
    schur = _finite(_sym(schur_factor @ schur_factor.T), "proper Gaussian wall Schur")
    nuisance_factor = permuted_factor[:nu, :nu]
    compensated_task = k[:, :nw] - k[:, nw:] @ _chol_solve(nuisance_factor, precision[nw:, :nw])
    with np.errstate(over="ignore", invalid="ignore"):
        task_covariance = _finite(_sym(k @ covariance @ k.T), "proper Gaussian task covariance")
        frozen_wall_task_covariance = _finite(_sym(k[:, nw:] @ _chol_solve(nuisance_factor, k[:, nw:].T)),
                                             "wall-conditioned task covariance")
        omitted_task_covariance = _finite(_sym(compensated_task @ _chol_solve(schur_factor, compensated_task.T)),
                                         "wall-omission task covariance")
    return dict(
        label=PROPER_GAUSSIAN, prior_assisted=True,
        prior_precision=_readonly(prior), joint_precision=_readonly(precision),
        joint_covariance=_readonly(covariance),
        wall_marginal_covariance=_readonly(covariance[:nw, :nw]),
        nuisance_marginal_covariance=_readonly(covariance[nw:, nw:]),
        wall_nuisance_cross_covariance=_readonly(covariance[:nw, nw:]),
        wall_schur_information=_readonly(schur), wall_schur_eigenvalues=_readonly(np.linalg.eigvalsh(schur)),
        task_covariance=_readonly(task_covariance),
        task_variance_trace=float(np.trace(task_covariance)),
        wall_conditioned_task_covariance=_readonly(frozen_wall_task_covariance),
        wall_omission_task_covariance=_readonly(omitted_task_covariance),
        wall_compensated_task_derivative=_readonly(compensated_task),
        coordinate_order="[wall, nuisance]", prior_role="explicit statistical prior; never GN/LM damping",
        error_scope="A/B/K operator error is not propagated into this local Gaussian covariance",
        bare_angle_interpretation=False, finite_posterior_certificate=False,
        physical_certificate=False, floating_point_roundoff_certified=False,
    )


def finite_pair_diagnostics(current_difference, task_difference, *, current_error,
                            task_error, task_norm, task_contract_id):
    """Record caller-supplied finite differences and absolute norm errors.

    The two errors are bounds/estimates for the *pair difference*, including
    both endpoints as appropriate. Current coordinates are already whitened
    real measurements. ``task_norm`` is explicitly ``euclidean`` or
    ``max_probe_l2``; the latter requires (probes, real measurements).
    A task/current ratio is withheld unless its denominator remains positive
    after the caller's error. No future point, truth, or reference is read.
    """
    current_raw = np.asarray(current_difference)
    if current_raw.ndim != 1 or current_raw.size == 0:
        raise ValueError("current_difference must be a nonempty real vector")
    current = _real_matrix(current_raw[None, :], "current_difference")[0]
    task = np.asarray(task_difference)
    if np.iscomplexobj(task) or task.dtype.kind not in "fiu" or task.size == 0 or not np.all(np.isfinite(task)):
        raise ValueError("task_difference must contain finite caller-realified real coordinates")
    task = np.asarray(task, dtype=float)
    if not isinstance(task_contract_id, str) or not task_contract_id:
        raise ValueError("nonempty task_contract_id required")
    if task_norm == "euclidean":
        if task.ndim == 0:
            raise ValueError("euclidean task_difference requires an explicit coordinate axis")
        task_magnitude = float(np.linalg.norm(task.ravel()))
    elif task_norm == "max_probe_l2":
        if task.ndim != 2 or task.shape[0] == 0 or task.shape[1] == 0:
            raise ValueError("max_probe_l2 requires (probes, real measurements)")
        task_magnitude = float(np.max(np.linalg.norm(task, axis=1)))
    else:
        raise ValueError("task_norm must be explicitly euclidean or max_probe_l2")
    current_error = _nonnegative_error(current_error, "current_error")
    task_error = _nonnegative_error(task_error, "task_error")
    current_magnitude = float(np.linalg.norm(current))
    _finite([current_magnitude, task_magnitude], "finite-pair norm")
    current_lower, task_lower = max(0.0, current_magnitude-current_error), max(0.0, task_magnitude-task_error)
    with np.errstate(over="ignore"):
        current_upper = float(_finite(current_magnitude+current_error, "current norm upper"))
        task_upper = float(_finite(task_magnitude+task_error, "task norm upper"))
    if current_magnitude == 0:
        ratio_status, ratio, ratio_lower, ratio_upper = "UNDEFINED_ZERO_DENOMINATOR", None, None, None
    elif current_lower == 0:
        ratio_status, ratio, ratio_lower, ratio_upper = "UNRESOLVED_DENOMINATOR", None, None, None
    else:
        with np.errstate(over="ignore"):
            ratios = np.asarray([task_magnitude/current_magnitude, task_lower/current_upper, task_upper/current_lower])
        if np.all(np.isfinite(ratios)):
            ratio_status = "DEFINED_CONDITIONAL_ON_CALLER_ERRORS"
            ratio, ratio_lower, ratio_upper = map(float, ratios)
        else:
            ratio_status, ratio, ratio_lower, ratio_upper = "UNRESOLVED_FLOATING_POINT_RANGE", None, None, None
    return dict(
        label="CALLER_SUPPLIED_FINITE_PAIR_EMPIRICAL_ERROR_DIAGNOSTIC",
        current_difference=_readonly(current), task_difference=_readonly(task),
        current_norm=current_magnitude, task_norm_value=task_magnitude,
        current_error=current_error, task_error=task_error,
        current_norm_lower=current_lower, current_norm_upper=current_upper,
        task_norm_lower=task_lower, task_norm_upper=task_upper,
        task_to_current_ratio=ratio, ratio_lower=ratio_lower, ratio_upper=ratio_upper, ratio_status=ratio_status,
        task_norm=task_norm, task_contract_id=task_contract_id,
        coordinate_assumption="caller-realified, current whitened; task weighting/whitening already applied",
        physical_certificate=False, finite_posterior_certificate=False, floating_point_roundoff_certified=False,
    )


def stress_metrics(a, b, k, *, prior_precision, operator_error):
    """Combined report with mandatory explicit joint prior and absolute error.

    Finite pair information is intentionally a separate caller-supplied input
    through ``finite_pair_diagnostics``; it is never inferred from K or files.
    """
    data = data_only_diagnostics(a, b, operator_error=operator_error)
    gaussian = proper_gaussian_diagnostics(a, b, k, prior_precision=prior_precision)
    return dict(
        status=data["status"], operator_error=data["operator_error"],
        data_only=data, proper_gaussian_prior=gaussian,
        scope="small supplied matrices; no PDE, truth/reference access, class assignment, or scientific gate",
    )
