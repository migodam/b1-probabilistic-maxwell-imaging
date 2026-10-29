"""Finite target likelihood contrasts, conditional on explicit proper priors.

This module has no bank, oracle, truth, reference-file or PDE interface. Strict
labels mean the caller supplied valid envelopes on the declared finite model;
they do not certify Maxwell discretization or floating-point roundoff. Empirical
secants remain empirical after grouping, normalization and TV propagation.
"""

from dataclasses import dataclass
import numpy as np


STRICT = "STRICT_ENVELOPE_CONDITIONAL_ON_SUPPLIED_BOUNDS"
EMPIRICAL = "EMPIRICAL_ENVELOPE_NO_UNIFORM_GUARANTEE"
NO_GUARANTEE = "NO_GUARANTEE"
LABELS = (STRICT, EMPIRICAL, NO_GUARANTEE)


def readonly(value, dtype=None):
    result = np.array(value, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


def logsumexp(values, axis=None):
    values = np.asarray(values, dtype=float)
    if not values.size:
        return -np.inf
    maximum = np.max(values, axis=axis, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        shifted = np.where(np.isfinite(maximum), values - maximum, -np.inf)
        total = maximum + np.log(np.sum(np.exp(shifted), axis=axis, keepdims=True))
    total = np.where(np.isposinf(maximum), np.inf, total)
    total = np.where(np.isneginf(maximum), -np.inf, total)
    if axis is None:
        return float(total.reshape(-1)[0])
    return np.squeeze(total, axis=axis)


def _logs(values):
    with np.errstate(divide="ignore"):
        return np.log(values)


def _log_array(values, name, *, allow_upper_infinity=False):
    values = np.asarray(values, dtype=float)
    if np.any(np.isnan(values)) or (not allow_upper_infinity and np.any(np.isposinf(values))):
        raise ValueError(f"{name} must be finite logs or declared infinities")
    return values


def _label(label):
    if label not in LABELS:
        raise ValueError(f"label must be one of {LABELS}")
    return label


@dataclass(frozen=True)
class FiniteProperPrior:
    """p(X) and explicit p(W|X), with W rows and X columns.

    Finite nonnegative input weights are deliberately normalized into proper
    masses. Every declared target has positive prior mass and a proper wall
    conditional. Zero conditional wall masses are retained, including support
    changes across targets. There is no implicit independent-prior assumption.
    """

    target_mass: np.ndarray
    wall_given_target: np.ndarray
    declaration: str = "explicit finite conditional masses"

    def __post_init__(self):
        target = np.asarray(self.target_mass, dtype=float)
        conditional = np.asarray(self.wall_given_target, dtype=float)
        if (target.ndim != 1 or not target.size or np.any(~np.isfinite(target))
                or np.any(target <= 0)):
            raise ValueError("finite positive target prior weights required")
        if (conditional.ndim != 2 or conditional.shape[1] != target.size or not conditional.shape[0]
                or np.any(~np.isfinite(conditional)) or np.any(conditional < 0)
                or np.any(conditional.sum(axis=0) <= 0)):
            raise ValueError("nonnegative proper wall conditional weights per target required")
        target_scale, conditional_scale = float(target.sum()), conditional.sum(axis=0)
        if not np.isfinite(target_scale) or np.any(~np.isfinite(conditional_scale)):
            raise ValueError("prior normalization exceeds finite range")
        object.__setattr__(self, "target_mass", readonly(target / target_scale))
        object.__setattr__(self, "wall_given_target", readonly(conditional / conditional_scale))
        object.__setattr__(self, "input_target_weight_total", target_scale)
        object.__setattr__(self, "input_wall_conditional_totals", readonly(conditional_scale))

    @classmethod
    def independent(cls, wall_mass, target_mass):
        wall = np.asarray(wall_mass, dtype=float)
        target = np.asarray(target_mass, dtype=float)
        if wall.ndim != 1:
            raise ValueError("wall prior vector required")
        return cls(target, np.repeat(wall[:, None], target.size, axis=1),
                   "explicit independent finite wall and target priors")

    @classmethod
    def from_joint(cls, joint_mass):
        joint = np.asarray(joint_mass, dtype=float)
        if joint.ndim != 2 or np.any(~np.isfinite(joint)) or np.any(joint < 0):
            raise ValueError("finite nonnegative joint mass matrix required")
        target = joint.sum(axis=0)
        if np.any(target <= 0):
            raise ValueError("remove undeclared zero-prior target columns explicitly")
        return cls(target, joint / target, "finite correlated prior disintegrated from explicit joint masses")

    @property
    def joint_mass(self):
        return readonly(self.wall_given_target * self.target_mass)


class AnchorCoverageError(ValueError):
    pass


@dataclass(frozen=True)
class ContrastIdentity:
    log_anchor_wall_mass: np.ndarray
    log_contrast: np.ndarray
    covered_wall: np.ndarray
    log_uncovered_residual: np.ndarray
    log_marginal_likelihood: np.ndarray
    anchor_index: int
    anchor_log_evidence: float
    coverage_complete: bool
    label: str = STRICT
    floating_point_roundoff_certified: bool = False


def target_marginal_contrast(log_likelihood, prior, anchor_index, *, allow_residual=False):
    """Exact finite identity including correlated-prior RN mass ratios.

    Base mass is L(w,x0)p(w|x0), and contrast is the log ratio of
    L(w,x)p(w|x) to that base. For independent priors the conditional masses
    cancel. If another target has contribution outside the anchor support, an
    explicit residual is required; it is never deleted or given infinite odds.
    """
    if not isinstance(prior, FiniteProperPrior):
        raise ValueError("explicit FiniteProperPrior required")
    logs = _log_array(log_likelihood, "likelihood")
    if (logs.shape != prior.wall_given_target.shape
            or not isinstance(anchor_index, (int, np.integer)) or isinstance(anchor_index, (bool, np.bool_))
            or not 0 <= anchor_index < logs.shape[1]):
        raise ValueError("likelihood matrix and legal target anchor required")
    log_conditional = _logs(prior.wall_given_target)
    weighted = logs + log_conditional
    base = weighted[:, int(anchor_index)]
    anchor_z = logsumexp(base)
    if not np.isfinite(anchor_z):
        raise AnchorCoverageError("anchor marginal likelihood must be finite and positive")
    covered = np.isfinite(base)
    contrast = np.full(logs.shape, np.nan)
    contrast[covered] = weighted[covered] - base[covered, None]
    residual = logsumexp(weighted[~covered], axis=0) if np.any(~covered) else np.full(logs.shape[1], -np.inf)
    complete = not np.any(np.isfinite(residual))
    if not complete and not allow_residual:
        raise AnchorCoverageError("anchor misses legal contributing walls; use multiple anchors or retain explicit residual")
    anchored = logsumexp(base[covered, None] + contrast[covered], axis=0)
    marginal = np.logaddexp(anchored, residual)
    return ContrastIdentity(readonly(base), readonly(contrast), readonly(covered), readonly(residual),
        readonly(marginal), int(anchor_index), anchor_z, complete)


@dataclass(frozen=True)
class MarginalEnvelope:
    log_lower: np.ndarray
    log_upper: np.ndarray
    label: str
    provenance: str
    group_labels: object = None
    log_group_anchor_mass: object = None
    floating_point_roundoff_certified: bool = False

    def __post_init__(self):
        lower = _log_array(self.log_lower, "lower")
        upper = _log_array(self.log_upper, "upper", allow_upper_infinity=True)
        if lower.ndim != 1 or not lower.size or upper.shape != lower.shape or np.any(lower > upper):
            raise ValueError("matching nonempty log marginal bounds lower<=upper required")
        object.__setattr__(self, "log_lower", readonly(lower))
        object.__setattr__(self, "log_upper", readonly(upper))
        object.__setattr__(self, "label", _label(self.label))


def group_anchor_envelope(log_anchor_wall_mass, wall_group, representative_indices,
                          representative_contrast, epsilon, *, coverage_complete,
                          label=EMPIRICAL, log_residual_lower=None, log_residual_upper=None,
                          provenance="caller-supplied group contrast envelope"):
    """Sum M_c0 exp(D_rep +/- epsilon); group weights are anchor masses.

    A finite empirical epsilon is not promoted into a strict guarantee. Infinite
    epsilon leaves an uninformative tail. Missing anchor support requires an
    explicitly supplied residual interval, even for an empirical output.
    """
    base = _log_array(log_anchor_wall_mass, "anchor masses")
    groups = np.asarray(wall_group)
    labels, inverse = np.unique(groups, return_inverse=True)
    reps = np.asarray(representative_indices)
    contrasts = _log_array(representative_contrast, "representative contrasts")
    error = np.asarray(epsilon, dtype=float)
    if (base.ndim != 1 or groups.shape != base.shape or reps.shape != (labels.size,)
            or reps.dtype.kind not in "iu" or np.any(reps < 0) or np.any(reps >= base.size)
            or contrasts.ndim != 2 or contrasts.shape[0] != labels.size):
        raise ValueError("one legal representative and contrast row per sorted group required")
    if np.any(inverse[reps] != np.arange(labels.size)):
        raise ValueError("each representative must be an existing wall in its group")
    if error.ndim == 0:
        error = np.full(contrasts.shape, error.item())
    if error.shape != contrasts.shape or np.any(np.isnan(error)) or np.any(error < 0):
        raise ValueError("nonnegative group/target error array required")
    group_base = np.array([logsumexp(base[inverse == c]) for c in range(labels.size)])
    lower, upper = np.full(contrasts.shape, -np.inf), np.full(contrasts.shape, -np.inf)
    active = np.isfinite(group_base)
    with np.errstate(invalid="ignore"):
        dlo, dup = contrasts - error, contrasts + error
    dlo[np.isposinf(error)], dup[np.isposinf(error)] = -np.inf, np.inf
    lower[active] = group_base[active, None] + dlo[active]
    upper[active] = group_base[active, None] + dup[active]
    result_lower, result_upper = logsumexp(lower, axis=0), logsumexp(upper, axis=0)
    if not coverage_complete and (log_residual_lower is None or log_residual_upper is None):
        raise AnchorCoverageError("uncovered anchor support needs an explicit residual interval")
    if log_residual_lower is not None or log_residual_upper is not None:
        if log_residual_lower is None or log_residual_upper is None:
            raise ValueError("both residual endpoints required")
        rlo = _log_array(log_residual_lower, "residual lower")
        rup = _log_array(log_residual_upper, "residual upper", allow_upper_infinity=True)
        if rlo.shape != result_lower.shape or rup.shape != rlo.shape or np.any(rlo > rup):
            raise ValueError("matching residual interval per target required")
        result_lower, result_upper = np.logaddexp(result_lower, rlo), np.logaddexp(result_upper, rup)
    return MarginalEnvelope(result_lower, result_upper, label, provenance, readonly(labels), readonly(group_base))


def _ratio(numerator, other, *, upper):
    if np.isneginf(numerator):
        return 0.0
    if np.isneginf(other):
        return 1.0
    if np.isposinf(numerator):
        return 1.0 if upper else 0.0
    if np.isposinf(other):
        return 0.0
    return float(np.exp(numerator - np.logaddexp(numerator, other)))


@dataclass(frozen=True)
class TargetEventInterval:
    lower: float
    upper: float
    label: str
    floating_point_roundoff_certified: bool = False


@dataclass(frozen=True)
class TargetPosteriorIntervals:
    log_mass_lower: np.ndarray
    log_mass_upper: np.ndarray
    log_evidence_lower: float
    log_evidence_upper: float
    label: str
    log_mass_offset: float = 0.0
    floating_point_roundoff_certified: bool = False

    def event(self, selected):
        selected = np.asarray(selected)
        if selected.shape != self.log_mass_lower.shape or selected.dtype.kind != "b":
            raise ValueError("boolean event on declared target atoms required")
        lo, up = logsumexp(self.log_mass_lower[selected]), logsumexp(self.log_mass_upper[selected])
        clo, cup = logsumexp(self.log_mass_lower[~selected]), logsumexp(self.log_mass_upper[~selected])
        if np.isneginf(up):
            lower = upper = 0.0
        elif np.isneginf(cup):
            lower = upper = 1.0
        else:
            lower, upper = _ratio(lo, cup, upper=False), _ratio(up, clo, upper=True)
        return TargetEventInterval(lower, upper, self.label)


def normalized_target_intervals(target_mass, envelope):
    target = np.asarray(target_mass, dtype=float)
    if (target.shape != envelope.log_lower.shape or np.any(~np.isfinite(target))
            or np.any(target < 0) or not np.any(target > 0)):
        raise ValueError("proper finite target weights matching envelope required")
    logp = _logs(target) - logsumexp(_logs(target))
    maximum_upper = float(np.max(envelope.log_upper))
    if np.isneginf(maximum_upper):
        raise ValueError("evidence upper bound is zero")
    offset = maximum_upper if np.isfinite(maximum_upper) else 0.0
    lower = logp + (envelope.log_lower - offset)
    upper = logp + (envelope.log_upper - offset)
    if np.any(np.isnan(lower)) or np.any(np.isnan(upper)):
        raise ValueError("zero prior combined with infinite upper bound needs explicit removal of zero-prior target atoms")
    zlo, zup = logsumexp(lower), logsumexp(upper)
    if np.isneginf(zup):
        raise ValueError("evidence upper bound is zero")
    return TargetPosteriorIntervals(readonly(lower), readonly(upper), offset + zlo, offset + zup,
                                    envelope.label, log_mass_offset=offset)


def target_posterior(target_mass, log_marginal_likelihood):
    target = np.asarray(target_mass, dtype=float)
    if target.ndim != 1 or np.any(~np.isfinite(target)) or np.any(target < 0) or not np.any(target > 0):
        raise ValueError("finite nonnegative target prior weights with positive total required")
    logp = _logs(target)
    logs = _log_array(log_marginal_likelihood, "marginal likelihood")
    if logp.shape != logs.shape or np.any(np.isnan(logp)):
        raise ValueError("matching nonnegative proper target masses required")
    # Center the likelihood first, retaining prior odds under a large common log.
    offset = float(np.max(logs))
    if not np.isfinite(offset):
        raise ValueError("positive finite target evidence required")
    mass = logp + (logs - offset)
    normalized = mass - logsumexp(mass)
    return readonly(np.exp(normalized)), readonly(normalized)


@dataclass(frozen=True)
class OscillationTVBound:
    log_error_oscillation: float
    conditional_TV_upper: float
    label: str
    guaranteed_for_true_likelihood: bool
    floating_point_roundoff_certified: bool = False


def tv_oscillation_bound(log_error_lower, log_error_upper, *, label=EMPIRICAL):
    """If delta=log(m_true/m_approx) has oscillation <=k, TV<=tanh(k/4).

    A shared constant log error has zero oscillation. For absolute error eta,
    the oscillation is at most 2eta; using tanh(eta/4) would be incorrect.
    """
    lower = _log_array(log_error_lower, "log-error lower")
    upper = _log_array(log_error_upper, "log-error upper", allow_upper_infinity=True)
    if lower.shape != upper.shape or lower.ndim != 1 or not lower.size or np.any(lower > upper):
        raise ValueError("matching ordered log-error endpoints required")
    oscillation = float(np.max(upper) - np.min(lower))
    if np.isnan(oscillation):
        oscillation = np.inf
    tag = _label(label)
    return OscillationTVBound(oscillation, float(np.tanh(oscillation / 4)), tag, tag == STRICT)
