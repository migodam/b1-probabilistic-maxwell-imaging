"""Finite-prior interval Bayes algebra, conditional on caller-supplied bounds.

Likelihoods omit the common Gaussian normalizing constant and lie in [0, 1].
No function in this module accepts reference weights, truth, or future data.
The mathematical label describes input bounds and finite-model algebra; it is
not a Maxwell stability, continuum, or floating-point roundoff certificate.
An empirical forward error remains EMPIRICAL_INTERVAL through normalization.
"""

from dataclasses import dataclass
from typing import Mapping

import numpy as np


MATHEMATICAL_INTERVAL = "MATHEMATICAL_INTERVAL"
EMPIRICAL_INTERVAL = "EMPIRICAL_INTERVAL"
NOT_CERTIFIED = "NOT_CERTIFIED"
QUERY_RULE = "PRIOR_LIKELIHOOD_INTERVAL_WIDTH_PER_COST"
_LABELS = (MATHEMATICAL_INTERVAL, EMPIRICAL_INTERVAL, NOT_CERTIFIED)


def _readonly(values, dtype=None):
    result = np.array(values, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


def _log_vector(values, name):
    result = np.asarray(values, dtype=float)
    if (result.ndim != 1 or result.size == 0 or np.any(np.isnan(result))
            or np.any(np.isposinf(result))):
        raise ValueError(f"{name} must be a nonempty vector of finite logs or -inf")
    return result


def _logsumexp(values):
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return -np.inf
    maximum = np.max(values)
    if np.isneginf(maximum):
        return -np.inf
    return float(maximum + np.log(np.exp(values - maximum).sum()))


def _normalized_logweights(values, name="weights"):
    weights = np.asarray(values, dtype=float)
    if (weights.ndim != 1 or weights.size == 0 or np.any(~np.isfinite(weights))
            or np.any(weights < 0) or not np.any(weights > 0)):
        raise ValueError(f"{name} must be finite nonnegative weights with positive total")
    with np.errstate(divide="ignore"):
        logs = np.log(weights)
    return logs - _logsumexp(logs)


def _labels(values, size):
    label_dtype = f"<U{max(map(len, _LABELS))}"
    labels = np.asarray(values, dtype=str)
    if labels.ndim == 0:
        labels = np.full(size, labels.item())
    if labels.shape != (size,) or np.any(~np.isin(labels, _LABELS)):
        raise ValueError(f"labels must contain only {_LABELS}")
    return _readonly(labels, dtype=label_dtype)


def _aggregate_label(labels):
    labels = np.asarray(labels)
    if np.any(labels == NOT_CERTIFIED):
        return NOT_CERTIFIED
    if np.any(labels == EMPIRICAL_INTERVAL):
        return EMPIRICAL_INTERVAL
    return MATHEMATICAL_INTERVAL


@dataclass(frozen=True)
class LogLikelihoodBounds:
    """Per-atom likelihood bounds; explicit inputs default to empirical status.

    ``queried`` records oracle access, not correctness. A queried approximate
    likelihood may still have positive interval width and need refinement.
    Arrays are copied and made read-only so registered factors cannot mutate.
    """

    log_lower: np.ndarray
    log_upper: np.ndarray
    labels: object = EMPIRICAL_INTERVAL
    queried: object = None
    provenance: str = "caller-supplied bounds"

    def __post_init__(self):
        lower = _log_vector(self.log_lower, "log_lower")
        upper = _log_vector(self.log_upper, "log_upper")
        if upper.shape != lower.shape or np.any(lower > upper) or np.any(upper > 0):
            raise ValueError("matching log likelihood bounds must satisfy lower <= upper <= 0")
        queried = np.ones(lower.size, dtype=bool) if self.queried is None else np.asarray(self.queried)
        if queried.shape != lower.shape or queried.dtype.kind != "b":
            raise ValueError("queried must be a matching boolean vector")
        object.__setattr__(self, "log_lower", _readonly(lower))
        object.__setattr__(self, "log_upper", _readonly(upper))
        object.__setattr__(self, "labels", _labels(self.labels, lower.size))
        object.__setattr__(self, "queried", _readonly(queried))

    @classmethod
    def unqueried(cls, size):
        """The Gaussian likelihood envelope [0, 1], without deleting atoms."""
        if not isinstance(size, (int, np.integer)) or isinstance(size, (bool, np.bool_)) or size < 1:
            raise ValueError("positive integer atom count required")
        return cls(np.full(size, -np.inf), np.zeros(size), MATHEMATICAL_INTERVAL,
                   np.zeros(size, dtype=bool), "unqueried Gaussian likelihood range [0,1]")

    @classmethod
    def exact(cls, log_likelihood, *, label=MATHEMATICAL_INTERVAL,
              provenance="caller-supplied exact finite-model likelihood"):
        """Treat supplied logs as exact within the declared discrete model."""
        return cls(log_likelihood, log_likelihood, label, provenance=provenance)

    def with_updates(self, indices, bounds):
        """Replace selected atoms, preserving all other atoms and their labels."""
        indices = np.asarray(indices)
        if (indices.ndim != 1 or indices.dtype.kind not in "iu" or indices.size == 0
                or np.any(indices < 0) or np.any(indices >= self.log_lower.size)
                or np.unique(indices).size != indices.size):
            raise ValueError("distinct in-range integer atom indices required")
        if not isinstance(bounds, LogLikelihoodBounds) or bounds.log_lower.size != indices.size:
            raise ValueError("one replacement interval per atom required")
        lower, upper = self.log_lower.copy(), self.log_upper.copy()
        labels, queried = self.labels.copy(), self.queried.copy()
        lower[indices], upper[indices] = bounds.log_lower, bounds.log_upper
        labels[indices], queried[indices] = bounds.labels, bounds.queried
        return LogLikelihoodBounds(lower, upper, labels, queried,
                                   f"{self.provenance}; atom replacement: {bounds.provenance}")


def likelihood_bounds_from_residual(residual_norm, forward_error, *,
                                    label=EMPIRICAL_INTERVAL, provenance="forward-error envelope"):
    """Use 0.5 max(0,r-e)^2 <= Phi <= 0.5(r+e)^2 in whitened space.

    An empirical e is never promoted to a mathematical certificate. Callers
    choosing MATHEMATICAL_INTERVAL own the stability/domain/roundoff assumptions.
    Infinite e is permitted and restores the conservative [0, 1] envelope.
    """
    residual = np.asarray(residual_norm, dtype=float)
    error = np.asarray(forward_error, dtype=float)
    if residual.ndim == 0:
        residual = residual.reshape(1)
    if error.ndim == 0:
        error = np.full(residual.shape, error.item())
    if (residual.ndim != 1 or residual.size == 0 or error.shape != residual.shape
            or np.any(~np.isfinite(residual)) or np.any(residual < 0)
            or np.any(np.isnan(error)) or np.any(error < 0)):
        raise ValueError("finite nonnegative residual norms and nonnegative matching errors required")
    with np.errstate(over="ignore"):
        minimum_residual = np.maximum(0.0, residual - error)
        maximum_residual = residual + error
        phi_lower = (0.5 * minimum_residual) * minimum_residual
        phi_upper = (0.5 * maximum_residual) * maximum_residual
    if np.any(np.isposinf(phi_lower)):
        raise ValueError("residual potential exceeds the representable log range")
    return LogLikelihoodBounds(-phi_upper, -phi_lower, label, provenance=provenance)


def likelihood_bounds_from_potential(potential, potential_error, *,
                                     label=EMPIRICAL_INTERVAL, provenance="potential-error envelope"):
    """Use |Phi-Phi_tilde| <= eta, including Phi >= 0 and likelihood <= 1."""
    potential = np.asarray(potential, dtype=float)
    error = np.asarray(potential_error, dtype=float)
    if potential.ndim == 0:
        potential = potential.reshape(1)
    if error.ndim == 0:
        error = np.full(potential.shape, error.item())
    if (potential.ndim != 1 or potential.size == 0 or error.shape != potential.shape
            or np.any(~np.isfinite(potential)) or np.any(potential < 0)
            or np.any(np.isnan(error)) or np.any(error < 0)):
        raise ValueError("finite nonnegative potentials and nonnegative matching errors required")
    with np.errstate(over="ignore"):
        lower = -(potential + error)
    upper = -np.maximum(0.0, potential - error)
    return LogLikelihoodBounds(lower, upper, label, provenance=provenance)


def combine_likelihood_bounds(factors, size=None):
    """Multiply distinct observation factors once, by adding their log bounds."""
    factors = tuple(factors.values() if isinstance(factors, Mapping) else factors)
    if not factors:
        if size is None:
            raise ValueError("atom count required for an empty factor collection")
        return LogLikelihoodBounds.exact(np.zeros(size), provenance="no observation factors")
    if any(not isinstance(bounds, LogLikelihoodBounds) for bounds in factors):
        raise ValueError("all factors must be LogLikelihoodBounds")
    count = factors[0].log_lower.size
    if (size is not None and count != size) or any(b.log_lower.size != count for b in factors):
        raise ValueError("all factor bounds must match the finite prior")
    lower = np.zeros(count)
    upper = np.zeros(count)
    labels = np.full(count, MATHEMATICAL_INTERVAL)
    queried = np.ones(count, dtype=bool)
    for bounds in factors:
        with np.errstate(over="ignore"):
            next_upper = upper + bounds.log_upper
            next_lower = lower + bounds.log_lower
        if np.any(np.isfinite(upper) & np.isfinite(bounds.log_upper) & np.isneginf(next_upper)):
            raise ValueError("combined upper likelihood log exceeds the representable range")
        lower, upper = next_lower, next_upper
        uncertified = (labels == NOT_CERTIFIED) | (bounds.labels == NOT_CERTIFIED)
        empirical = (labels == EMPIRICAL_INTERVAL) | (bounds.labels == EMPIRICAL_INTERVAL)
        labels[empirical] = EMPIRICAL_INTERVAL
        labels[uncertified] = NOT_CERTIFIED
        queried &= bounds.queried
    return LogLikelihoodBounds(lower, upper, labels, queried, "product of unique observation factors")


def _excluding_logsumexp(values):
    """O(n) complement sums without subtracting nearly equal normalizers."""
    left = np.r_[-np.inf, np.logaddexp.accumulate(values)[:-1]]
    right = np.r_[np.logaddexp.accumulate(values[::-1])[::-1][1:], -np.inf]
    return np.logaddexp(left, right)


def _ratios(log_numerator, log_other):
    numerator, other = np.broadcast_arrays(log_numerator, log_other)
    result = np.zeros(numerator.shape, dtype=float)
    possible = ~np.isneginf(numerator)
    result[possible] = np.exp(numerator[possible] - np.logaddexp(numerator[possible], other[possible]))
    return result


@dataclass(frozen=True)
class EventPosteriorBounds:
    lower: float
    upper: float
    log_numerator_lower: float
    log_numerator_upper: float
    log_complement_lower: float
    log_complement_upper: float
    label: str
    floating_point_roundoff_certified: bool = False
    log_mass_offset: float = 0.0


@dataclass(frozen=True)
class FinitePosteriorBounds:
    """Joint atom bounds and evidence bounds on the original finite prior.

    Marginal interval endpoints need not sum to one. Log mass arrays omit the
    shared ``log_mass_offset``; this preserves prior ratios even with a very
    large common negative likelihood log. Evidence logs include that offset.
    ``joint_weights`` and ``log_joint_weights`` are provided only when all
    supported inputs are exact.
    """

    log_mass_lower: np.ndarray
    log_mass_upper: np.ndarray
    log_evidence_lower: float
    log_evidence_upper: float
    lower: np.ndarray
    upper: np.ndarray
    label: str
    joint_weights: object
    log_joint_weights: object
    floating_point_roundoff_certified: bool = False
    log_mass_offset: float = 0.0

    def event(self, selected):
        selected = np.asarray(selected)
        if selected.shape != self.log_mass_lower.shape or selected.dtype.kind != "b":
            raise ValueError("event must be a boolean mask on the original joint atoms")
        lower = _logsumexp(self.log_mass_lower[selected])
        upper = _logsumexp(self.log_mass_upper[selected])
        other_lower = _logsumexp(self.log_mass_lower[~selected])
        other_upper = _logsumexp(self.log_mass_upper[~selected])
        # Structural zero/full events stay exact even when Z_lower == 0.
        if np.isneginf(upper):
            lo = hi = 0.0
        elif np.isneginf(other_upper):
            lo = hi = 1.0
        else:
            lo = float(_ratios(lower, other_upper))
            hi = float(_ratios(upper, other_lower))
        return EventPosteriorBounds(lo, hi, lower, upper, other_lower, other_upper,
                                    self.label, log_mass_offset=self.log_mass_offset)


def finite_posterior_bounds(prior_weight, likelihood_bounds):
    """Normalize finite-prior mass intervals with event/complement denominators.

    Posterior statements are conditional on positive evidence. An identically
    zero upper evidence is rejected instead of inventing a posterior. A zero
    evidence lower bound is allowed and usually yields broad [0, 1] intervals.
    """
    return _posterior_from_logprior(_normalized_logweights(prior_weight, "prior weights"), likelihood_bounds)


def _posterior_from_logprior(log_prior, likelihood_bounds):
    if not isinstance(likelihood_bounds, LogLikelihoodBounds) or likelihood_bounds.log_lower.shape != log_prior.shape:
        raise ValueError("likelihood bounds must match the finite joint prior")
    supported = ~np.isneginf(log_prior)
    offset = float(np.max(likelihood_bounds.log_upper[supported]))
    if np.isneginf(offset):
        raise ValueError("empty posterior support: evidence upper bound is zero")
    # Center before adding log prior: -1e100 + log(.1) would otherwise erase
    # the prior odds before any logsumexp routine had a chance to normalize.
    lower = log_prior + (likelihood_bounds.log_lower - offset)
    upper = log_prior + (likelihood_bounds.log_upper - offset)
    scaled_evidence_lower, scaled_evidence_upper = _logsumexp(lower), _logsumexp(upper)
    complement_upper = _excluding_logsumexp(upper)
    complement_lower = _excluding_logsumexp(lower)
    atom_lower = _ratios(lower, complement_upper)
    atom_upper = _ratios(upper, complement_lower)
    sole_possible = ~np.isneginf(upper) & np.isneginf(complement_upper)
    atom_lower[sole_possible] = atom_upper[sole_possible] = 1.0
    label = _aggregate_label(likelihood_bounds.labels[supported])
    exact = np.all(lower == upper)
    log_joint = _readonly(upper - scaled_evidence_upper) if exact else None
    joint = _readonly(np.exp(log_joint)) if exact else None
    return FinitePosteriorBounds(_readonly(lower), _readonly(upper), offset + scaled_evidence_lower,
                                offset + scaled_evidence_upper, _readonly(atom_lower), _readonly(atom_upper),
                                label, joint, log_joint, log_mass_offset=offset)


def event_posterior_bounds(prior_weight, likelihood_bounds, selected):
    return finite_posterior_bounds(prior_weight, likelihood_bounds).event(selected)


def query_log_scores(prior_weight, likelihood_bounds, cost=1.0, *, eligible=None):
    """Predeclared rule: log(prior_i * (L_i^+ - L_i^-) / cost_i).

    This is an uncertainty-contribution heuristic, not an optimal-risk theorem.
    It depends only on the prior, current bounds, and declared query costs. The
    log implementation retains ordering even when mass widths would underflow.
    Scores omit a common upper-likelihood offset; their ordering is unchanged.
    """
    log_prior = _normalized_logweights(prior_weight, "prior weights")
    return _query_log_scores_from_logprior(log_prior, likelihood_bounds, cost, eligible=eligible)


def _query_log_scores_from_logprior(log_prior, likelihood_bounds, cost=1.0, *, eligible=None):
    if not isinstance(likelihood_bounds, LogLikelihoodBounds) or likelihood_bounds.log_lower.shape != log_prior.shape:
        raise ValueError("likelihood bounds must match the finite joint prior")
    costs = np.asarray(cost, dtype=float)
    if costs.ndim == 0:
        costs = np.full(log_prior.shape, costs.item())
    if costs.shape != log_prior.shape or np.any(~np.isfinite(costs)) or np.any(costs <= 0):
        raise ValueError("finite positive query costs required")
    if eligible is None:
        eligible = np.ones(log_prior.shape, dtype=bool)
    eligible = np.asarray(eligible)
    if eligible.shape != log_prior.shape or eligible.dtype.kind != "b":
        raise ValueError("eligible must be a matching boolean mask")
    lower, upper = likelihood_bounds.log_lower, likelihood_bounds.log_upper
    width = np.full(log_prior.shape, -np.inf)
    nonzero = (lower < upper) & eligible & ~np.isneginf(log_prior)
    if not np.any(nonzero):
        return _readonly(width)
    offset = np.max(upper[nonzero])
    with np.errstate(divide="ignore"):
        width[nonzero] = (upper[nonzero] - offset) + np.log(-np.expm1(lower[nonzero] - upper[nonzero]))
    scores = log_prior + width - np.log(costs)
    scores[~eligible] = -np.inf
    return _readonly(scores)


def select_query(prior_weight, likelihood_bounds, cost=1.0, *, eligible=None):
    """Select maximum contribution/cost; deterministic first-index tie break.

    Return None if every eligible prior-weighted interval has zero width. This
    does not imply a physical certificate or establish any scientific gate.
    """
    scores = query_log_scores(prior_weight, likelihood_bounds, cost, eligible=eligible)
    index = int(np.argmax(scores))
    return None if np.isneginf(scores[index]) else index


@dataclass(frozen=True)
class IntervalFactor:
    factor_id: str
    observation_id: str
    version: int
    bounds: LogLikelihoodBounds


class FinitePriorIntervals:
    """One current factor per observation, always recomputed from the prior.

    A refresh replaces the old version; it never multiplies a second copy of
    the same observation. Reconstruction also permits old zero-support atoms
    to recover, unlike a log(new)-log(old) correction with -inf entries.
    """

    def __init__(self, prior_weight):
        self.log_prior_weight = _readonly(_normalized_logweights(prior_weight, "prior weights"))
        self.prior_weight = _readonly(np.exp(self.log_prior_weight))
        self._factors = {}
        self._observation_owner = {}

    def set_factor(self, factor_id, bounds=None, *, observation_id=None, version=1):
        if not isinstance(factor_id, str) or not factor_id:
            raise ValueError("nonempty string factor ID required")
        if not isinstance(version, (int, np.integer)) or isinstance(version, (bool, np.bool_)) or version < 1:
            raise ValueError("positive integer version required")
        old = self._factors.get(factor_id)
        if observation_id is None:
            observation_id = old.observation_id if old else factor_id
        if not isinstance(observation_id, str) or not observation_id:
            raise ValueError("nonempty string observation ID required")
        if old and (version <= old.version or observation_id != old.observation_id):
            raise ValueError("replacement must preserve observation ID and increase version")
        owner = self._observation_owner.get(observation_id)
        if owner is not None and owner != factor_id:
            raise ValueError("duplicate likelihood observation ID")
        if bounds is None:
            bounds = LogLikelihoodBounds.unqueried(self.prior_weight.size)
        if not isinstance(bounds, LogLikelihoodBounds) or bounds.log_lower.shape != self.prior_weight.shape:
            raise ValueError("factor intervals must match the finite prior")
        factor = IntervalFactor(factor_id, observation_id, int(version), bounds)
        self._factors[factor_id] = factor
        self._observation_owner[observation_id] = factor_id
        return factor

    def refine_factor(self, factor_id, indices, bounds, *, version=None):
        old = self._factors[factor_id]
        replacement = old.bounds.with_updates(indices, bounds)
        return self.set_factor(factor_id, replacement, observation_id=old.observation_id,
                               version=old.version + 1 if version is None else version)

    def factors(self):
        return tuple(self._factors.values())

    def combined_bounds(self):
        return combine_likelihood_bounds([factor.bounds for factor in self._factors.values()],
                                         self.prior_weight.size)

    def posterior(self):
        return _posterior_from_logprior(self.log_prior_weight, self.combined_bounds())

    def query_log_scores(self, cost=1.0, *, eligible=None):
        return _query_log_scores_from_logprior(self.log_prior_weight, self.combined_bounds(),
                                               cost, eligible=eligible)

    def select_query(self, cost=1.0, *, eligible=None):
        scores = self.query_log_scores(cost, eligible=eligible)
        index = int(np.argmax(scores))
        return None if np.isneginf(scores[index]) else index


def _partition(cluster_ids, size):
    clusters = np.asarray(cluster_ids)
    if clusters.shape != (size,) or clusters.dtype.kind not in "iuUSb":
        raise ValueError("one integer or string cluster ID per joint atom required")
    labels, inverse = np.unique(clusters, return_inverse=True)
    return labels, inverse


@dataclass(frozen=True)
class JointPushforward:
    joint_weights: np.ndarray
    log_joint_weights: np.ndarray
    cluster_labels: np.ndarray
    atom_cluster_index: np.ndarray
    cluster_mass: np.ndarray
    log_cluster_mass: np.ndarray
    conditional_weights: np.ndarray
    log_conditional_weights: np.ndarray
    interpretation: str = "joint pushforward with original within-cluster conditional kernel"

    def lift(self):
        """Restore the original atom joint distribution from retained conditionals."""
        return _readonly(np.exp(self.log_cluster_mass[self.atom_cluster_index] + self.log_conditional_weights))


def joint_pushforward(joint_weights, cluster_ids):
    """Push a joint measure into clusters without factorizing wall and target."""
    logs = _normalized_logweights(joint_weights, "joint weights")
    labels, inverse = _partition(cluster_ids, logs.size)
    mass = np.array([_logsumexp(logs[inverse == index]) for index in range(labels.size)])
    conditional = np.full(logs.shape, -np.inf)
    positive = ~np.isneginf(mass[inverse])
    conditional[positive] = logs[positive] - mass[inverse[positive]]
    return JointPushforward(_readonly(np.exp(logs)), _readonly(logs), _readonly(labels),
                            _readonly(inverse), _readonly(np.exp(mass)), _readonly(mass),
                            _readonly(np.exp(conditional)), _readonly(conditional))


@dataclass(frozen=True)
class CanonicalLift:
    joint_weights: np.ndarray
    log_joint_weights: np.ndarray
    cluster_labels: np.ndarray
    cluster_mass: np.ndarray
    prior_conditional_kernel: np.ndarray
    log_prior_conditional_kernel: np.ndarray
    joint_error_label: str = NOT_CERTIFIED
    interpretation: str = "quotient mass lifted by the original prior conditional kernel"


def canonical_lift(cluster_mass, prior_weight, cluster_ids):
    """Lift quotient masses using the prior's conditional disintegration.

    Cluster masses follow sorted unique cluster IDs, as in joint_pushforward.
    This defines a full atom measure for comparison; it does not certify that
    the original joint posterior was preserved, and never reports zero joint TV.
    """
    prior = joint_pushforward(prior_weight, cluster_ids)
    mass = _normalized_logweights(cluster_mass, "cluster masses")
    if mass.size != prior.cluster_labels.size:
        raise ValueError("one mass per sorted unique cluster ID required")
    if np.any(~np.isneginf(mass) & np.isneginf(prior.log_cluster_mass)):
        raise ValueError("positive quotient mass cannot be lifted into a zero-prior cluster")
    logs = mass[prior.atom_cluster_index] + prior.log_conditional_weights
    return CanonicalLift(_readonly(np.exp(logs)), _readonly(logs), prior.cluster_labels,
                         _readonly(np.exp(mass)), prior.conditional_weights,
                         prior.log_conditional_weights)


@dataclass(frozen=True)
class TaskClusterDiameters:
    cluster_labels: np.ndarray
    diameter: np.ndarray
    member_count: np.ndarray
    task_contract_id: str
    norm: str
    metric: str
    norm_definition: str
    coordinate_assumption: str
    probe_axis: int | None
    measurement_axis: int | None
    input_shape: tuple
    scope: str = "exhaustive finite supplied task points; no future-task or PDE certificate"
    floating_point_roundoff_certified: bool = False


def task_cluster_diameters(task_values, cluster_ids, *, task_contract_id, norm):
    """Compute all-pairs diameter, preventing a neighbor chain from posing as a cover.

    ``norm`` is mandatory; a concatenated L2 norm must not silently stand in
    for a task contract requiring a maximum over probes.

    ``max_probe_l2`` requires already whitened real measurements with shape
    (n_atoms, n_probes, n_real_measurements). For each atom pair it takes the
    L2 norm over measurements separately for each probe, then the maximum
    over probes. This function never performs the whitening itself.

    ``euclidean`` takes L2 over all supplied coordinates concatenated. It is
    weighted L2 only when the caller has already applied the desired weights
    to those coordinates. It also retains the earlier scalar/complex-coordinate
    behavior, but requires the explicit norm choice.

    Coordinates and units define the declared metric. The caller retains this
    contract ID and must reopen/recompute when the task contract changes.
    """
    values = np.asarray(task_values)
    input_shape = values.shape
    if norm not in ("max_probe_l2", "euclidean"):
        raise ValueError("norm must be explicitly max_probe_l2 or euclidean")
    if norm == "max_probe_l2" and values.ndim != 3:
        raise ValueError("max_probe_l2 requires an explicit probe axis: (atoms, probes, real measurements)")
    if norm == "max_probe_l2" and values.dtype.kind == "c":
        raise ValueError("max_probe_l2 requires prewhitened real measurement coordinates")
    if norm == "euclidean" and values.ndim == 1:
        values = values[:, None]
    if values.ndim < 2 or values.size == 0 or values.dtype.kind not in "fciu" or np.any(~np.isfinite(values)):
        raise ValueError("finite task coordinates, with the joint atom axis first, required")
    if not isinstance(task_contract_id, str) or not task_contract_id:
        raise ValueError("nonempty declared task contract ID required")
    labels, inverse = _partition(cluster_ids, values.shape[0])
    values = values.astype(complex if values.dtype.kind == "c" else float)
    if norm == "max_probe_l2":
        metric = "max_probe_l2 on supplied prewhitened real measurements"
        definition = "max_probe ||delta[probe, :]||_2"
        assumption = "caller-prewhitened real measurements; no whitening applied by this function"
        probe_axis, measurement_axis = 1, 2
    else:
        values = values.reshape(values.shape[0], -1)
        metric = "euclidean on supplied coordinates; weightedL2 only if caller-preweighted"
        definition = "||ravel(delta)||_2"
        assumption = "caller-supplied coordinates; any weighting or whitening already applied"
        probe_axis = measurement_axis = None
    diameters = np.zeros(labels.size)
    counts = np.bincount(inverse, minlength=labels.size)
    for index in range(labels.size):
        points = values[inverse == index]
        for point_index in range(points.shape[0] - 1):
            difference = points[point_index + 1:] - points[point_index]
            if norm == "max_probe_l2":
                distances = np.max(np.linalg.norm(difference, axis=2), axis=1)
            else:
                distances = np.linalg.norm(difference, axis=1)
            diameters[index] = max(diameters[index], float(np.max(distances)))
    return TaskClusterDiameters(
        _readonly(labels), _readonly(diameters), _readonly(counts), task_contract_id,
        norm, metric, definition, assumption, probe_axis, measurement_axis, input_shape,
    )
