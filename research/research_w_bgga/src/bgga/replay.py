"""Zero-new-PDE development replay on an explicitly finite response bank.

Online entry points accept only ``(prior, y, oracle)``.  An oracle call returns
one atom's base-TX response; it has no table, truth, reference, or file API.
Every method builds and pays for its own queried responses.  File access and
the later, independent reference calculation belong to the tool layer.

The interval algebra is conditional on the saved discrete responses.  Neither
floating-point roundoff nor an underlying Maxwell error is certified here.
"""

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

from bgga.intervals import (
    LogLikelihoodBounds,
    MATHEMATICAL_INTERVAL,
    finite_posterior_bounds,
)


REPLAY_LABEL = "REPLAY_NOT_NEW_PDE_SPEED_EVIDENCE"
SIGMA_REAL_IMAG = 0.01
M4_UNKNOWN_MASS_LIMIT = 0.02
M4_CHECKPOINTS = (16, 64, 256, 1024, 2048, 2949, 2950)
M5_POTENTIAL_DIAMETER = float(np.arctanh(0.02))
M5_FUTURE_DIAMETER = 0.1
SECANT_SCALES = (1.0, 0.1, 0.3, 0.2)
G0_PARAMETER_NAMES = ("epsinf", "ell", "pose_x", "pose_y", "contrast", "target_loss", "radius")
G0_PARAMETER_SCALES = (1.0, 0.1, 0.3, 0.2, 1.0, 0.1, 0.04)
JOINT_PARTITION_LABEL = "PREDECLARED_JOINT_PARTITION_CELLS_NOT_CONTINUOUS_MODE_DISCOVERY"
ArrayOracle = Callable[[int], np.ndarray]


def _readonly(values, dtype=None):
    result = np.array(values, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


def _logsumexp(values):
    values = np.asarray(values, dtype=float)
    if not values.size:
        return -np.inf
    offset = float(np.max(values))
    if np.isneginf(offset):
        return offset
    return float(offset + np.log(np.exp(values - offset).sum()))


@dataclass(frozen=True)
class SecantQualification:
    """Frozen G0 numerical metadata, never a current-case truth/reference.

    The threshold is max_fixture ||packj(J-J_reference) diag(scales)||_2.
    It is a cross-fixture empirical threshold, not a bank-secants error bound.
    """

    operator_error: Optional[float] = None
    fixture_ids: tuple = ()
    fixture_errors: tuple = ()
    source: str = "G0 derivative qualification unavailable"
    status: str = "NOT_DEFINED"

    def __post_init__(self):
        if self.operator_error is not None:
            if not np.isfinite(self.operator_error) or self.operator_error < 0:
                raise ValueError("finite nonnegative empirical G0 error required")
            if self.status != "EMPIRICAL_CROSS_FIXTURE_THRESHOLD":
                raise ValueError("a supplied threshold must retain its empirical label")


@dataclass(frozen=True)
class FinitePriorGrid:
    """Prior coordinates plus predeclared numerical qualification metadata.

    Atoms need not be ordered by wall: ``atom_by_pose_group`` establishes the
    actual Cartesian lookup, including the historical pose-major ordering.
    """

    prior_weight: np.ndarray
    pose: np.ndarray
    wall_group_id: np.ndarray
    support_id: np.ndarray
    support: np.ndarray
    epsinf: np.ndarray
    ell: np.ndarray
    qualification: SecantQualification = field(default_factory=SecantQualification)
    prior_id: str = "FINITE_PRIOR_ENUM4_V1"

    def __post_init__(self):
        weight = np.asarray(self.prior_weight, dtype=float)
        n = weight.size
        if (weight.ndim != 1 or not n or np.any(~np.isfinite(weight))
                or np.any(weight <= 0)):
            raise ValueError("strictly positive finite atom prior weights required")
        pose = np.asarray(self.pose, dtype=float)
        if pose.shape != (n, 2) or np.any(~np.isfinite(pose)):
            raise ValueError("finite two-coordinate pose per atom required")
        support = np.asarray(self.support)
        if support.shape != (n, 4) or support.dtype.kind != "b":
            raise ValueError("four boolean occupancy cells per atom required")
        for name in ("wall_group_id", "support_id"):
            values = np.asarray(getattr(self, name))
            if values.shape != (n,) or values.dtype.kind not in "iu":
                raise ValueError(f"integer {name} per atom required")
            object.__setattr__(self, name, _readonly(values, np.int64))
        for name in ("epsinf", "ell"):
            values = np.asarray(getattr(self, name), dtype=float)
            if values.shape != (n,) or np.any(~np.isfinite(values)):
                raise ValueError(f"finite {name} per atom required")
            object.__setattr__(self, name, _readonly(values))
        if not isinstance(self.qualification, SecantQualification):
            raise ValueError("explicit SecantQualification metadata required")
        object.__setattr__(self, "prior_weight", _readonly(weight))
        object.__setattr__(self, "pose", _readonly(pose))
        object.__setattr__(self, "support", _readonly(support))
        groups, group_index = np.unique(self.wall_group_id, return_inverse=True)
        poses, pose_index = np.unique(pose, axis=0, return_inverse=True)
        lookup = np.full((len(poses), len(groups)), -1, dtype=np.int64)
        if n != lookup.size:
            raise ValueError("one atom per wall-group/pose pair required")
        for atom, (p, g) in enumerate(zip(pose_index, group_index)):
            if lookup[p, g] != -1:
                raise ValueError("duplicate wall-group/pose pair")
            lookup[p, g] = atom
        if np.any(lookup < 0):
            raise ValueError("incomplete Cartesian wall-group/pose grid")
        for g in range(len(groups)):
            atoms = lookup[:, g]
            for name in ("support_id", "support", "epsinf", "ell"):
                value = getattr(self, name)[atoms]
                if np.any(value != value[0]):
                    raise ValueError(f"wall group changes {name} across poses")
        object.__setattr__(self, "group_ids", _readonly(groups))
        object.__setattr__(self, "group_index", _readonly(group_index))
        object.__setattr__(self, "pose_nodes", _readonly(poses))
        object.__setattr__(self, "pose_index", _readonly(pose_index))
        object.__setattr__(self, "atom_by_pose_group", _readonly(lookup))

    @property
    def size(self):
        return self.prior_weight.size


@dataclass(frozen=True)
class MethodOutput:
    method_id: str
    status: str
    joint_weights: object
    log_joint_weights: object
    query_indices: np.ndarray
    query_log_likelihood: np.ndarray
    details: dict = field(default_factory=dict)
    checkpoints: tuple = ()
    label: str = REPLAY_LABEL
    actual_full_fe: int = 0
    actual_rhs: int = 0
    floating_point_roundoff_certified: bool = False

    @property
    def logical_queries(self):
        return self.query_indices.size


class _OwnQueries:
    def __init__(self, prior, y, oracle):
        if not isinstance(prior, FinitePriorGrid) or not callable(oracle):
            raise ValueError("finite prior and atom-only callable oracle required")
        y = np.asarray(y, dtype=np.complex128)
        if y.shape != (2, 16) or np.any(~np.isfinite(y)):
            raise ValueError("observation must contain base TX0/1 and all 16 RX")
        self.prior, self.y, self.oracle = prior, y, oracle
        self.indices, self.logs = [], []
        self.seen = set()
        self.responses = np.empty((prior.size, 32), dtype=np.complex128)
        self.log_likelihood = np.full(prior.size, np.nan)

    def query(self, index):
        index = int(index)
        if index in self.seen:
            raise ValueError("an exact queried atom must not be queried again")
        if not 0 <= index < self.prior.size:
            raise ValueError("atom index out of range")
        response = np.asarray(self.oracle(index), dtype=np.complex128)
        if response.shape != (2, 16) or np.any(~np.isfinite(response)):
            raise ValueError("oracle must return one finite base-TX 2x16 response")
        self.responses[index] = response.ravel()
        residual = (response - self.y).ravel() / SIGMA_REAL_IMAG
        potential = float(0.5 * (np.dot(residual.real, residual.real)
                                 + np.dot(residual.imag, residual.imag)))
        if not np.isfinite(potential):
            raise ValueError("whitened potential exceeds finite representable range")
        self.log_likelihood[index] = -potential
        self.seen.add(index)
        self.indices.append(index)
        self.logs.append(-potential)

    def full(self):
        for index in range(self.prior.size):
            self.query(index)
        return _exact_posterior(self.prior, self.log_likelihood)

    def output(self, method_id, posterior=None, *, status="COMPLETED", details=None, checkpoints=()):
        return MethodOutput(
            method_id, status,
            None if posterior is None else _readonly(posterior[0]),
            None if posterior is None else _readonly(posterior[1]),
            _readonly(self.indices, np.int64), _readonly(self.logs),
            {} if details is None else details, tuple(checkpoints))


def _exact_posterior(prior, log_likelihood):
    posterior = finite_posterior_bounds(prior.prior_weight, LogLikelihoodBounds.exact(log_likelihood))
    return posterior.joint_weights, posterior.log_joint_weights, posterior.log_evidence_lower


def _group_masses(prior, weights):
    return np.bincount(prior.group_index, weights=weights, minlength=len(prior.group_ids))


def _point_conditional(prior, log_joint, group_index):
    selected = prior.group_index == group_index
    # A medoid may have underflowing displayed mass: use its own finite logs.
    logs = np.full(prior.size, -np.inf)
    log_mass = _logsumexp(log_joint[selected])
    if np.isneginf(log_mass):
        raise ValueError("selected wall has no posterior conditional")
    logs[selected] = log_joint[selected] - log_mass
    return np.exp(logs), logs


def m0_direct(prior, y, oracle):
    own = _OwnQueries(prior, y, oracle)
    q, logs, logz = own.full()
    return own.output("M0", (q, logs), details={
        "definition": "independent full enumeration exact finite Bayes",
        "log_evidence": logz, "posterior_kind": "EXACT_SAVED_FINITE_MODEL"})


def m1_confidence_point(prior, y, oracle):
    own = _OwnQueries(prior, y, oracle)
    q, logs, _ = own.full()
    mass = _group_masses(prior, q)
    group = int(np.argmax(mass))  # group IDs are sorted; exact ties use stable ID.
    return own.output("M1", _point_conditional(prior, logs, group), details={
        "definition": "own posterior MAP wall, retaining target conditional at that wall",
        "selected_wall_group_id": int(prior.group_ids[group]),
        "own_wall_posterior_mass": mass,
        "posterior_kind": "POINT_WALL_COMPRESSION_ABLATION", "safe_by_default": False})


def _secant_indices(axis, value):
    axis = np.asarray(axis)
    here = int(np.flatnonzero(axis == value)[0])
    if axis.size < 2:
        return None
    lo, hi = max(0, here - 1), min(axis.size - 1, here + 1)
    return float(axis[lo]), float(axis[hi]), ("central_secant" if lo < here < hi else "one_sided_secant")


def m2_angle_point(prior, y, oracle):
    own = _OwnQueries(prior, y, oracle)
    qualification = prior.qualification
    if qualification.operator_error is None:
        return own.output("M2", status="NOT_RUN", details={
            "geometry_status": "NOT_DEFINED", "reason": qualification.source,
            "geometry_label": "EMPIRICAL_BANK_SECANT_RETENTION",
            "rank_threshold": None, "safe_by_default": False})
    q, logs, _ = own.full()
    mass = _group_masses(prior, q)
    map_atom = int(np.argmax(logs))
    map_pose = int(prior.pose_index[map_atom])
    lookup = prior.atom_by_pose_group
    group_atoms = lookup[map_pose]
    parameter_lookup = {
        (int(prior.support_id[a]), float(prior.epsinf[a]), float(prior.ell[a])): g
        for g, a in enumerate(group_atoms)}
    pose_lookup = {tuple(p): i for i, p in enumerate(prior.pose_nodes)}
    empirical_error = float(qualification.operator_error)
    rows, legal_groups = [], []
    for group, atom in enumerate(group_atoms):
        row = {"wall_group_id": int(prior.group_ids[group]), "status": "NOT_DEFINED",
               "min_retention": None, "rank_A": None, "rank_B": None, "stencils": []}
        if not np.any(prior.support[atom]):
            row["reason"] = "empty wall has no material secant support"
            rows.append(row)
            continue
        sid = int(prior.support_id[atom])
        same_support = prior.support_id[group_atoms] == sid
        axes = (np.unique(prior.epsinf[group_atoms[same_support]]),
                np.unique(prior.ell[group_atoms[same_support]]),
                np.unique(prior.pose_nodes[:, 0]), np.unique(prior.pose_nodes[:, 1]))
        values = (prior.epsinf[atom], prior.ell[atom], *prior.pose_nodes[map_pose])
        columns = []
        for direction, (axis, value) in enumerate(zip(axes, values)):
            secant = _secant_indices(axis, value)
            if secant is None:
                break
            lo, hi, kind = secant
            if direction < 2:
                material_lo = [float(prior.epsinf[atom]), float(prior.ell[atom])]
                material_hi = material_lo.copy()
                material_lo[direction], material_hi[direction] = lo, hi
                gl = parameter_lookup.get((sid, *material_lo))
                gh = parameter_lookup.get((sid, *material_hi))
                if gl is None or gh is None:
                    break
                il, ih = int(lookup[map_pose, gl]), int(lookup[map_pose, gh])
            else:
                pl, ph = list(prior.pose_nodes[map_pose]), list(prior.pose_nodes[map_pose])
                pl[direction - 2], ph[direction - 2] = lo, hi
                if tuple(pl) not in pose_lookup or tuple(ph) not in pose_lookup:
                    break
                il, ih = int(lookup[pose_lookup[tuple(pl)], group]), int(lookup[pose_lookup[tuple(ph)], group])
            column = (own.responses[ih] - own.responses[il]) / (hi - lo)
            columns.append(np.r_[column.real, column.imag] * (SECANT_SCALES[direction] / SIGMA_REAL_IMAG))
            row["stencils"].append({"parameter": G0_PARAMETER_NAMES[direction],
                "minus_atom": il, "plus_atom": ih, "node_span": hi - lo,
                "kind": kind, "scale": SECANT_SCALES[direction]})
        if len(columns) != 4:
            row["reason"] = "legal finite-grid material/pose neighbors unavailable"
            rows.append(row)
            continue
        matrix = np.column_stack(columns)
        ua, sa, _ = np.linalg.svd(matrix[:, :2], full_matrices=False)
        ub, sb, _ = np.linalg.svd(matrix[:, 2:], full_matrices=False)
        rank_a, rank_b = int(np.sum(sa > empirical_error)), int(np.sum(sb > empirical_error))
        row.update(rank_A=rank_a, rank_B=rank_b, singular_values_A=sa,
                   singular_values_B=sb, rank_threshold=empirical_error)
        if not rank_a:
            row["reason"] = "material secants have zero effective rank at empirical threshold"
            rows.append(row)
            continue
        a_support = ua[:, :rank_a]
        b_support = ub[:, :rank_b]
        projected = a_support - b_support @ (b_support.T @ a_support)
        # Raw squared singular values are retained; no clipping to [0,1].
        retention = np.linalg.svd(projected, compute_uv=False) ** 2
        row.update(status="DEFINED_EMPIRICALLY", retention_spectrum=retention,
                   min_retention=float(np.min(retention)))
        legal_groups.append((group, row["min_retention"]))
        rows.append(row)
    details = {"definition": "point wall selected by finite-bank material/pose secant retention",
        "geometry_label": "EMPIRICAL_BANK_SECANT_RETENTION",
        "legacy_geometry_label": "EMPIRICAL_BANK_DIFFERENCE_GEOMETRY",
        "not_a_local_Jacobian": True, "not_a_uniform_derivative_bound": True,
        "rank_threshold": empirical_error, "threshold_scope": qualification.status,
        "qualification_fixture_ids": qualification.fixture_ids,
        "qualification_fixture_errors": qualification.fixture_errors,
        "secant_parameter_scales": SECANT_SCALES,
        "own_joint_MAP_atom": map_atom, "common_pose_node": prior.pose_nodes[map_pose],
        "group_geometry": rows, "safe_by_default": False,
        "posterior_kind": "POINT_WALL_COMPRESSION_ABLATION"}
    if not legal_groups:
        details["reason"] = "no legal nonempty wall with positive effective material rank"
        return own.output("M2", status="NOT_DEFINED", details=details)
    # No data-dependent tolerance or revised threshold: exact ties use own mass then ID.
    group, _ = min(legal_groups, key=lambda item: (-item[1], -mass[item[0]], int(prior.group_ids[item[0]])))
    details["selected_wall_group_id"] = int(prior.group_ids[group])
    return own.output("M2", _point_conditional(prior, logs, group), details=details)


def m3_task_risk_point(prior, y, oracle):
    own = _OwnQueries(prior, y, oracle)
    q, logs, _ = own.full()
    table = own.responses[prior.atom_by_pose_group]
    joint = q[prior.atom_by_pose_group]
    risk = np.empty(len(prior.group_ids))
    for representative in range(len(prior.group_ids)):
        difference = (table - table[:, representative:representative + 1]) / SIGMA_REAL_IMAG
        squared = np.sum(difference.real ** 2 + difference.imag ** 2, axis=2)
        risk[representative] = float(np.sum(joint * squared))
    selected = int(np.argmin(risk))
    return own.output("M3", _point_conditional(prior, logs, selected), details={
        "definition": "legal wall medoid minimizing own posterior whitened response risk at matched poses",
        "selected_wall_group_id": int(prior.group_ids[selected]), "candidate_risk": risk,
        "candidate_wall_group_ids": prior.group_ids,
        "candidate_count": len(prior.group_ids), "matched_pose_count": len(prior.pose_nodes),
        "risk_computation_uses_own_full_queries": True,
        "posterior_kind": "POINT_WALL_COMPRESSION_ABLATION", "safe_by_default": False})


def _interval_checkpoint(prior, own, count, *, final=False):
    queried = np.zeros(prior.size, dtype=bool)
    queried[own.indices] = True
    lower, upper = np.full(prior.size, -np.inf), np.zeros(prior.size)
    lower[queried] = upper[queried] = own.log_likelihood[queried]
    bounds = LogLikelihoodBounds(lower, upper, MATHEMATICAL_INTERVAL, queried,
        "saved finite-model exact queried values, otherwise Gaussian envelope [0,1]")
    posterior = finite_posterior_bounds(prior.prior_weight, bounds)
    unknown = posterior.event(~queried)
    # The displayed point estimate conditions only on paid queried atoms.  Bounds
    # retain the original prior and every unqueried atom; no true posterior is invented.
    point = _exact_posterior(prior, lower)
    return {"query_count": count, "final": final, "status": "EVALUATED",
        "queried": _readonly(queried), "log_likelihood_lower": bounds.log_lower,
        "log_likelihood_upper": bounds.log_upper,
        "posterior_lower": posterior.lower, "posterior_upper": posterior.upper,
        "log_evidence_lower": posterior.log_evidence_lower,
        "log_evidence_upper": posterior.log_evidence_upper,
        "unknown_mass_lower": unknown.lower, "unknown_mass_upper": unknown.upper,
        "unknown_log_prior_numerator_upper": unknown.log_numerator_upper + unknown.log_mass_offset,
        "point_joint_weights": point[0], "point_log_joint_weights": point[1],
        "point_estimate_kind": ("EXACT_SAVED_FINITE_MODEL" if queried.all()
                                 else "QUERIED_CONDITIONAL_POINT_ESTIMATE"),
        "interval_label": posterior.label, "floating_point_roundoff_certified": False}


def m4_multigreen_interval(prior, y, oracle):
    own = _OwnQueries(prior, y, oracle)
    # All unknown widths and costs equal one.  Sorting once is exactly the
    # preregistered prior*width/cost rule; the order does not depend on y.
    with np.errstate(divide="ignore"):
        log_prior = np.log(prior.prior_weight)
    log_prior -= _logsumexp(log_prior)
    order = np.lexsort((np.arange(prior.size), -log_prior))
    suffix = np.full(prior.size + 1, -np.inf)
    suffix[:-1] = np.logaddexp.accumulate(log_prior[order][::-1])[::-1]
    log_known, checkpoints, trace = -np.inf, [], []
    requested = set(M4_CHECKPOINTS)
    for count, index in enumerate(order, 1):
        own.query(int(index))
        log_known = float(np.logaddexp(log_known, log_prior[index] + own.logs[-1]))
        log_unknown = float(suffix[count])
        log_upper_z = float(np.logaddexp(log_known, log_unknown))
        unknown_upper = 0.0 if np.isneginf(log_unknown) else float(np.exp(log_unknown - log_upper_z))
        trace.append({"query_count": count, "atom_index": int(index),
            "log_priority": float(log_prior[index]), "log_likelihood": own.logs[-1],
            "log_evidence_lower": log_known, "log_evidence_upper": log_upper_z,
            "unknown_mass_upper": unknown_upper})
        stop = unknown_upper <= M4_UNKNOWN_MASS_LIMIT
        if count in requested or stop or count == prior.size:
            checkpoints.append(_interval_checkpoint(prior, own, count, final=stop or count == prior.size))
        if stop:
            break
    final = checkpoints[-1]
    return own.output("M4", (final["point_joint_weights"], final["point_log_joint_weights"]), details={
        "definition": "prior times likelihood interval width per unit logical query cost",
        "query_rule": "PRIOR_LIKELIHOOD_INTERVAL_WIDTH_PER_COST",
        "uniform_cost": 1.0, "unknown_likelihood_interval": [0.0, 1.0],
        "unknown_mass_stop_threshold": M4_UNKNOWN_MASS_LIMIT,
        "unknown_mass_upper_at_stop": final["unknown_mass_upper"],
        "requested_checkpoints": M4_CHECKPOINTS,
        "unreached_checkpoints": [n for n in M4_CHECKPOINTS if n > own.prior.size or n > len(own.indices)],
        "posterior_kind": final["point_estimate_kind"],
        "original_prior_support_retained_in_bounds": True,
        "query_trace": trace, "interval_label": MATHEMATICAL_INTERVAL,
        "scope": "saved finite-model likelihood envelope; roundoff and Maxwell errors uncertified"},
        checkpoints=checkpoints)


def complete_link_wall_clusters(group_ids, potential_distance, future_distance):
    """Stable first-fit complete-link partition; no transitive chain merging."""
    group_ids = np.asarray(group_ids)
    n = group_ids.size
    potential_distance = np.asarray(potential_distance, dtype=float)
    future_distance = np.asarray(future_distance, dtype=float)
    for distance in (potential_distance, future_distance):
        if distance.shape != (n, n) or np.any(~np.isfinite(distance)) or np.any(distance < 0):
            raise ValueError("finite nonnegative pairwise wall-group distances required")
        if not np.array_equal(distance, distance.T):
            raise ValueError("pairwise wall-group distances must be symmetric")
    clusters = []
    for group in np.argsort(group_ids, kind="stable"):
        for members in clusters:
            if (np.all(potential_distance[group, members] <= M5_POTENTIAL_DIAMETER)
                    and np.all(future_distance[group, members] <= M5_FUTURE_DIAMETER)):
                members.append(int(group))
                break
        else:
            clusters.append([int(group)])
    return tuple(tuple(members) for members in clusters)


def m5_quotient(prior, y, oracle):
    own = _OwnQueries(prior, y, oracle)
    own.full()  # All clustering and likelihood computations are deployment work.
    table = own.responses[prior.atom_by_pose_group]
    potentials = -own.log_likelihood[prior.atom_by_pose_group]
    groups = len(prior.group_ids)
    potential_distance, future_distance = np.zeros((groups, groups)), np.zeros((groups, groups))
    for left in range(groups):
        for right in range(left + 1, groups):
            dp = float(np.max(np.abs(potentials[:, left] - potentials[:, right])))
            difference = (table[:, left] - table[:, right]) / SIGMA_REAL_IMAG
            df = float(np.sqrt(np.max(np.sum(difference.real ** 2 + difference.imag ** 2, axis=1))))
            potential_distance[left, right] = potential_distance[right, left] = dp
            future_distance[left, right] = future_distance[right, left] = df
    clusters = complete_link_wall_clusters(prior.group_ids, potential_distance, future_distance)
    cluster_of_group, representative = np.empty(groups, dtype=int), np.empty(groups, dtype=int)
    cluster_rows = []
    for cluster, members in enumerate(clusters):
        cluster_of_group[list(members)] = cluster
        representative[list(members)] = members[0]
        selection = np.ix_(members, members)
        cluster_rows.append({"cluster_id": cluster,
            "wall_group_ids": prior.group_ids[list(members)],
            "representative_wall_group_id": int(prior.group_ids[members[0]]),
            "potential_diameter": float(np.max(potential_distance[selection])),
            "future_max_pose_diameter": float(np.max(future_distance[selection]))})
    representative_atom = prior.atom_by_pose_group[prior.pose_index, representative[prior.group_index]]
    q, logs, logz = _exact_posterior(prior, own.log_likelihood[representative_atom])
    # Quotient cells are (wall cluster, target pose).  The legal representative
    # likelihood is constant on a cell; original prior conditionals lift it.
    atom_cell = prior.pose_index * len(clusters) + cluster_of_group[prior.group_index]
    cell_prior = np.bincount(atom_cell, weights=prior.prior_weight)
    kernel = prior.prior_weight / cell_prior[atom_cell]
    cell_mass = np.bincount(atom_cell, weights=q, minlength=len(cell_prior))
    return own.output("M5", (q, logs), details={
        "definition": "stable complete-link wall quotient using legal representative likelihoods and prior conditional lift",
        "potential_diameter_limit": M5_POTENTIAL_DIAMETER,
        "future_max_pose_diameter_limit": M5_FUTURE_DIAMETER,
        "future_scope": "25 prior pose nodes with fixed known target intrinsic parameters and gamma",
        "G0_other_radius_contrast_task_cluster_qualification": "NOT_RUN_NOT_IN_THIS_TASK_SCOPE",
        "cluster_count": len(clusters), "clusters": cluster_rows,
        "group_ids": prior.group_ids, "group_cluster_id": cluster_of_group,
        "potential_pair_distance": potential_distance, "future_pair_distance": future_distance,
        "atom_quotient_cell": atom_cell, "prior_conditional_kernel": kernel,
        "quotient_cell_mass": cell_mass, "quotient_log_evidence": logz,
        "all_table_queries_attributed_to_this_method": True,
        "joint_error_label": "NOT_CERTIFIED_AWAIT_INDEPENDENT_SCORING",
        "posterior_kind": "REPRESENTATIVE_LIKELIHOOD_PRIOR_CONDITIONAL_LIFT"})


def m6_cache_schur_control(prior, y, oracle):
    own = _OwnQueries(prior, y, oracle)
    q, logs, logz = own.full()
    return own.output("M6", (q, logs), details={
        "definition": "ordinary cache/Schur logical inference control: own exact full table",
        "posterior_kind": "EXACT_SAVED_FINITE_MODEL", "log_evidence": logz,
        "measured_Schur_speed_status": "NOT_RUN", "physics_saving_claim": False})


METHODS = (m0_direct, m1_confidence_point, m2_angle_point, m3_task_risk_point,
           m4_multigreen_interval, m5_quotient, m6_cache_schur_control)
METHOD_DEFINITIONS = {
    "M0": "Independent 2950-query exact finite Bayes.",
    "M1": "Own full-query MAP wall-group; retain its own target conditional and Dirac wall.",
    "M2": "Own full-query legal nonempty wall maximizing empirical finite-secants minimum retention; missing qualification is NOT_RUN.",
    "M3": "Own full-query legal wall medoid minimizing matched-pose posterior whitened response risk.",
    "M4": "Unknown likelihood [0,1]; prior*interval-width/cost queries, stop at unknown-mass upper <=.02.",
    "M5": "Own full-query stable complete-link wall quotient: potential diameter <=atanh(.02), max-pose whitened response diameter <=.1, prior conditional lift.",
    "M6": "Ordinary cache/Schur logical exact full-table control; no measured Schur speed claim.",
}


@dataclass(frozen=True)
class OfflineReference:
    joint_weights: np.ndarray
    log_joint_weights: np.ndarray
    log_likelihood: np.ndarray
    log_evidence: float
    logical_queries: int
    minimum_prior_mass: float
    default_unknown_envelope_floor: float
    default_unknown_envelope_floor_log: float


def offline_reference(prior, y, oracle):
    """Call only after all online outputs are sealed; never feed back to METHODS."""
    own = _OwnQueries(prior, y, oracle)
    q, logs, logz = own.full()
    normalized_prior = _exact_posterior(prior, np.zeros(prior.size))[0]
    pmin = float(np.min(normalized_prior))
    logpmin = float(np.min(np.log(prior.prior_weight))) - _logsumexp(np.log(prior.prior_weight))
    log_floor = float(logpmin - np.logaddexp(logz, logpmin))
    return OfflineReference(_readonly(q), _readonly(logs), _readonly(own.log_likelihood),
        logz, prior.size, pmin, float(np.exp(log_floor)), log_floor)


def _mode_masses(prior, q, *, support_only=False):
    labels = ((int(prior.support_id[i]),) if support_only else
              (int(prior.support_id[i]), "nonnegative" if prior.pose[i, 0] >= 0 else "negative",
               bool(prior.pose[i, 1] >= 2.1)) for i in range(prior.size))
    result = {}
    for label, weight in zip(labels, q):
        result[label] = result.get(label, 0.0) + float(weight)
    return result


def score_joint(prior, output, reference):
    """Offline joint-distribution scores.  No truth/clean-response input exists."""
    if output.joint_weights is None:
        return {"method_id": output.method_id, "status": output.status,
                "joint_TV": None, "coordinate_CDF_sup": None,
                "event_probability_errors": None, "important_joint_modes": None,
                "position_Bayes_regret": None, "ground_truth_loss": None,
                "ground_truth_loss_status": "NOT_SCORED_NO_GT_ACCESS",
                "joint_mode_interpretation": JOINT_PARTITION_LABEL}
    q = np.asarray(output.joint_weights)
    ref = np.asarray(reference.joint_weights)
    if q.shape != (prior.size,) or np.any(~np.isfinite(q)) or np.any(q < 0):
        raise ValueError("valid full joint method weights required")
    if abs(float(q.sum()) - 1.0) > 1e-10:
        raise ValueError("method joint weights must sum to one")
    coordinate = {}
    cdf_details = {}
    for axis, name in enumerate(("x", "y")):
        nodes = np.unique(prior.pose[:, axis])
        q_cdf = np.array([q[prior.pose[:, axis] <= v].sum() for v in nodes])
        r_cdf = np.array([ref[prior.pose[:, axis] <= v].sum() for v in nodes])
        coordinate[name] = float(np.max(np.abs(q_cdf - r_cdf)))
        cdf_details[name] = {"nodes": nodes, "method": q_cdf, "reference": r_cdf}
    events = {}
    for name, mask in (("pose_x>=0", prior.pose[:, 0] >= 0), ("pose_y>=2.1", prior.pose[:, 1] >= 2.1)):
        actual, expected = float(q[mask].sum()), float(ref[mask].sum())
        events[name] = {"method_probability": actual, "reference_probability": expected,
                        "signed_error": actual - expected, "absolute_error": abs(actual - expected)}
    joint_modes, ref_modes = _mode_masses(prior, q), _mode_masses(prior, ref)
    mode_rows = []
    for label in sorted(ref_modes):
        expected, actual = ref_modes[label], joint_modes[label]
        mode_rows.append({"support_id": label[0], "x_sign": label[1], "y>=2.1": label[2],
            "reference_mass": expected, "method_mass": actual,
            "signed_error": actual - expected, "absolute_error": abs(actual - expected),
            "important": expected >= 0.05})
    important = [row for row in mode_rows if row["important"]]
    support_q, support_r = _mode_masses(prior, q, support_only=True), _mode_masses(prior, ref, support_only=True)
    support_rows = [{"support_id": label[0], "reference_mass": support_r[label],
        "method_mass": support_q[label], "absolute_error": abs(support_q[label] - support_r[label])}
        for label in sorted(support_r)]
    mean_q, mean_r = q @ prior.pose, ref @ prior.pose
    ref_risk = float(np.sum(ref * np.sum((prior.pose - mean_r) ** 2, axis=1)))
    regret = float(np.sum((mean_q - mean_r) ** 2))
    limit = 0.05 * ref_risk + 0.005 ** 2
    return {"method_id": output.method_id, "status": output.status,
        "joint_TV": float(0.5 * np.sum(np.abs(q - ref))),
        "coordinate_CDF_sup": coordinate, "coordinate_CDF": cdf_details,
        "event_probability_errors": events, "joint_modes": mode_rows,
        "joint_mode_interpretation": JOINT_PARTITION_LABEL,
        "x_sign_partition": "negative iff x<0; nonnegative iff x>=0; no numerical-zero tolerance",
        "important_joint_modes": important,
        "important_joint_mode_max_absolute_error": max((r["absolute_error"] for r in important), default=None),
        "support_modes_separately": support_rows,
        "posterior_mean_method": mean_q, "posterior_mean_reference": mean_r,
        "position_Bayes_regret": regret, "reference_position_Bayes_risk": ref_risk,
        "position_Bayes_regret_limit": limit, "position_Bayes_regret_minus_limit": regret - limit,
        "ground_truth_loss": None, "ground_truth_loss_status": "NOT_SCORED_NO_GT_ACCESS",
        "floating_point_roundoff_certified": False}


def score_interval_checkpoint(prior, checkpoint, reference):
    """Retain raw signed violations; never clip tiny positive violations away."""
    lower, upper, exact = checkpoint["posterior_lower"], checkpoint["posterior_upper"], reference.joint_weights
    queried = checkpoint["queried"]
    actual_unknown = float(exact[~queried].sum())
    unknown_signed = actual_unknown - checkpoint["unknown_mass_upper"]
    ll_lower_violation = checkpoint["log_likelihood_lower"] - reference.log_likelihood
    ll_upper_violation = reference.log_likelihood - checkpoint["log_likelihood_upper"]
    return {"query_count": checkpoint["query_count"],
        "actual_reference_unknown_mass": actual_unknown,
        "unknown_mass_upper": checkpoint["unknown_mass_upper"],
        "unknown_mass_upper_signed_violation": unknown_signed,
        "atom_lower_signed_violation": lower - exact,
        "atom_upper_signed_violation": exact - upper,
        "maximum_atom_lower_signed_violation": float(np.max(lower - exact)),
        "maximum_atom_upper_signed_violation": float(np.max(exact - upper)),
        "likelihood_lower_signed_violation": ll_lower_violation,
        "likelihood_upper_signed_violation": ll_upper_violation,
        "any_positive_raw_bound_violation": bool(unknown_signed > 0 or np.any(lower - exact > 0)
            or np.any(exact - upper > 0) or np.any(ll_lower_violation > 0) or np.any(ll_upper_violation > 0)),
        "interval_label": checkpoint["interval_label"],
        "floating_point_roundoff_certified": False}
