"""Truth-free post-correction archive; no solver, corrector or oracle imports.

Validation flags/receipts must be supplied by the common external validator.
Software acceptance here does not independently establish Maxwell feasibility.
See research/delegated/g_source_audit/ARCHIVE_SPEC.md for frozen semantics.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.ndimage import gaussian_filter, label


METRIC_VERSION = "G-multiscale-support-v1"
MATERIAL_THRESHOLD = .03
STRUCTURE_THRESHOLD = .08
XOR_MIN_CELLS = 6
_FOUR = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
_EIGHT = np.ones((3, 3), dtype=int)


def _frozen_array(value, dtype=np.complex128):
    array = np.asarray(value, dtype=dtype)
    return np.frombuffer(array.tobytes(order="C"), dtype=array.dtype).reshape(array.shape)


def canonical_field(value):
    value = np.asarray(value, dtype=np.complex128)
    if value.shape not in {(576,), (24, 24)}:
        raise ValueError("canonical A30 field must have 576 cells")
    return _frozen_array(value.reshape(24, 24))


def material_distance(first, second):
    delta = canonical_field(first) - canonical_field(second)
    return float(np.sqrt(np.mean((delta.real / 4.5) ** 2 + (delta.imag / .5) ** 2)))


@dataclass(frozen=True)
class StructureFeatures:
    occupancy: np.ndarray
    original_mask: np.ndarray
    version: str = METRIC_VERSION


def structure_features(field, background):
    amplitude = np.abs(canonical_field(field) - canonical_field(background)) / 1.2
    if not np.isfinite(amplitude).all():
        raise ValueError("finite material/background required")
    views = []
    for sigma in (0., 1., 2.):
        smoothed = (amplitude if sigma == 0 else
                    gaussian_filter(amplitude, sigma, mode="constant", cval=0., truncate=4.))
        views.append(np.clip((smoothed - .20) / .10, 0., 1.))
    return StructureFeatures(_frozen_array(np.stack(views), np.float64),
                             _frozen_array(amplitude > .25, np.bool_))


def structural_distance(first: StructureFeatures, second: StructureFeatures):
    if first.version != METRIC_VERSION or second.version != METRIC_VERSION:
        raise ValueError("structure feature version mismatch")
    distance = float(np.sqrt(np.mean((first.occupancy - second.occupancy) ** 2)))
    xor = int(np.count_nonzero(first.original_mask ^ second.original_mask))
    return distance, xor


def structural_distinct(first, second, *, threshold=STRUCTURE_THRESHOLD):
    distance, xor = structural_distance(first, second)
    return bool(distance >= threshold and xor >= XOR_MIN_CELLS)


def topology_summary(features: StructureFeatures):
    """Posthoc grid descriptors only; no hole count enters search/selection."""
    mask = np.asarray(features.original_mask, bool)
    foreground, nforeground = label(mask, _FOUR)
    background, nbackground = label(~mask, _EIGHT)
    fg_sizes = np.bincount(foreground.ravel(), minlength=nforeground + 1)[1:]
    outside = set(np.r_[background[0], background[-1], background[:, 0], background[:, -1]])
    bg_sizes = np.bincount(background.ravel(), minlength=nbackground + 1)
    hole_sizes = [int(bg_sizes[index]) for index in range(1, nbackground + 1)
                  if index not in outside]
    return dict(foreground_connectivity=4, background_connectivity=8,
                components_all=int(nforeground), component_sizes=fg_sizes.tolist(),
                components_ge6=int(np.count_nonzero(fg_sizes >= 6)),
                holes_all=len(hole_sizes), hole_sizes=hole_sizes,
                holes_ge6=sum(size >= 6 for size in hole_sizes),
                scope="POSTHOC_DIGITAL_GRID_NOT_CONTINUUM_TOPOLOGY")


@dataclass(frozen=True, eq=False)
class Candidate:
    candidate_id: str
    field: np.ndarray
    prediction: np.ndarray
    full_residual: float
    full_validated: bool
    solver_validated: bool
    validation_receipt: str | None
    is_base: bool = False

    def __post_init__(self):
        if not isinstance(self.candidate_id, str) or not self.candidate_id.strip():
            raise ValueError("nonempty candidate ID required")
        object.__setattr__(self, "field", canonical_field(self.field))
        object.__setattr__(self, "prediction", _frozen_array(self.prediction))
        object.__setattr__(self, "full_residual", float(self.full_residual))

    @property
    def posterior_weight(self):
        return None

    @property
    def represented_mass(self):
        return None

    @property
    def unmodeled_mass(self):
        return None


@dataclass(frozen=True)
class OfferResult:
    candidate_id: str
    status: str
    active: bool
    new_world_reward: int
    replaced_id: str | None = None
    reason: str | None = None


class OfferBudgetExhausted(RuntimeError):
    """Software endpoint cap only, never a substitute for a paid-solver budget."""


class Archive:
    def __init__(self, background, *, data_threshold, max_offers=None,
                 material_threshold=MATERIAL_THRESHOLD,
                 structure_threshold=STRUCTURE_THRESHOLD, replacement_mode="quality"):
        self.background = canonical_field(background)
        if not np.isfinite(self.background).all():
            raise ValueError("finite declared known background required")
        if not np.isfinite(data_threshold) or data_threshold < 0:
            raise ValueError("declared nonnegative full data threshold required")
        if max_offers is not None and (isinstance(max_offers, bool)
                                       or not isinstance(max_offers, int) or max_offers < 0):
            raise ValueError("max_offers must be a nonnegative integer or None")
        if (not np.isfinite(material_threshold) or material_threshold <= 0
                or not np.isfinite(structure_threshold) or structure_threshold <= 0):
            raise ValueError("strictly positive fixed separation thresholds required")
        if replacement_mode not in {"quality", "first"}:
            raise ValueError("replacement_mode must be 'quality' or 'first'")
        self.data_threshold = float(data_threshold)
        self.material_threshold = float(material_threshold)
        self.structure_threshold = float(structure_threshold)
        self.max_offers = max_offers
        self.replacement_mode = replacement_mode
        self.offers = 0
        self._seen_ids = set()
        self._successes: list[Candidate] = []
        self._active: list[Candidate] = []
        self._frontier: list[Candidate] = []
        self._failures: list[dict] = []
        self._features: dict[str, StructureFeatures] = {}
        self._events: list[OfferResult] = []

    @property
    def success_ledger(self):
        return tuple(self._successes)

    @property
    def active_archive(self):
        return tuple(self._active)

    @property
    def frontier(self):
        return tuple(self._frontier)

    @property
    def failure_ledger(self):
        return tuple(dict(item) for item in self._failures)

    @property
    def event_ledger(self):
        return tuple(self._events)

    def features(self, candidate):
        if candidate.candidate_id in self._features:
            return self._features[candidate.candidate_id]
        return structure_features(candidate.field, self.background)

    def pair(self, first, second):
        material = material_distance(first.field, second.field)
        structure, xor = structural_distance(self.features(first), self.features(second))
        structural = structure >= self.structure_threshold and xor >= XOR_MIN_CELLS
        return dict(dM=material, dS=structure, xor_cells=xor,
                    material_distinct=material >= self.material_threshold,
                    structural_distinct=bool(structural),
                    separated=bool(material >= self.material_threshold or structural))

    def _validation_failure(self, candidate):
        field = candidate.field
        if not np.isfinite(field).all():
            return "NONFINITE_MATERIAL"
        if np.any((field.real < -.5) | (field.real > 4.) | (field.imag < 0.) | (field.imag > .5)):
            return "ILLEGAL_MATERIAL_BOX"
        if candidate.prediction.shape != (32, 8) or not np.isfinite(candidate.prediction).all():
            return "INVALID_A30_FULL_PREDICTION"
        if (not np.isfinite(candidate.full_residual) or candidate.full_residual < 0
                or candidate.full_residual > self.data_threshold):
            return "DATA_INCOMPATIBLE_OR_NONFINITE_RESIDUAL"
        if not candidate.full_validated or not candidate.solver_validated:
            return "MISSING_EXTERNAL_FULL_OR_SOLVER_VALIDATION"
        if not isinstance(candidate.validation_receipt, str) or not candidate.validation_receipt.strip():
            return "MISSING_FULL_VALIDATION_RECEIPT"
        return None

    def offer(self, candidate: Candidate, *, retain_frontier=True):
        if self.max_offers is not None and self.offers >= self.max_offers:
            raise OfferBudgetExhausted("endpoint offer cap exhausted; prior results retained")
        if candidate.candidate_id in self._seen_ids:
            raise ValueError("endpoint receipt ID already submitted; use a new ID for a new attempt")
        self.offers += 1
        self._seen_ids.add(candidate.candidate_id)
        failure = self._validation_failure(candidate)
        if failure:
            self._failures.append(dict(candidate_id=candidate.candidate_id, reason=failure,
                                       full_residual=candidate.full_residual))
            if retain_frontier:
                self._frontier.append(candidate)
            result = OfferResult(candidate.candidate_id, "REJECTED_INVALID_OR_UNVERIFIED", False, 0,
                                 reason=failure)
            self._events.append(result)
            return result
        self._features[candidate.candidate_id] = structure_features(candidate.field, self.background)
        # Historical successes are deliberately never pruned or merged.
        reward = int(all(self.pair(candidate, old)["separated"] for old in self._successes))
        self._successes.append(candidate)
        neighbors = [old for old in self._active if not self.pair(candidate, old)["separated"]]
        if not neighbors:
            self._active.append(candidate)
            result = OfferResult(candidate.candidate_id, "INSERTED", True, reward)
        elif (self.replacement_mode == "quality" and len(neighbors) == 1
              and candidate.full_residual < neighbors[0].full_residual):
            old = neighbors[0]
            others = [item for item in self._active if item.candidate_id != old.candidate_id]
            if all(self.pair(candidate, other)["separated"] for other in others):
                self._active[self._active.index(old)] = candidate
                result = OfferResult(candidate.candidate_id, "REPLACED_LOWER_FULL_RESIDUAL", True,
                                     reward, replaced_id=old.candidate_id)
            else:
                result = OfferResult(candidate.candidate_id, "VALID_DUPLICATE_RETAINED_IN_LEDGER", False, 0)
        else:
            result = OfferResult(candidate.candidate_id, "VALID_DUPLICATE_RETAINED_IN_LEDGER", False, 0)
        self._events.append(result)
        self.assert_pairwise_packing()
        return result

    def assert_pairwise_packing(self):
        for index, first in enumerate(self._active):
            for second in self._active[index + 1:]:
                if not self.pair(first, second)["separated"]:
                    raise AssertionError("active archive violates pairwise separation")
        return True

    def _subpacking(self, criterion):
        result = []
        for candidate in sorted(self._active, key=lambda item: item.candidate_id):
            if all(self.pair(candidate, old)[criterion] for old in result):
                result.append(candidate)
        return tuple(result)

    @property
    def material_representatives(self):
        return self._subpacking("material_distinct")

    @property
    def structural_representatives(self):
        return self._subpacking("structural_distinct")

    def top1(self):
        return min(self._successes, key=lambda item: (item.full_residual, item.candidate_id),
                   default=None)

    def topk(self, k, *, preserve_base=True):
        if isinstance(k, bool) or not isinstance(k, int) or k < 0:
            raise ValueError("K must be a nonnegative integer")
        if not k or not self._successes:
            return tuple()
        candidates = sorted(self._active, key=lambda item: item.candidate_id)
        bases = sorted((item for item in self._successes if item.is_base),
                       key=lambda item: item.candidate_id)
        selected = [bases[0]] if preserve_base and bases else []
        if not selected:
            # Fixed ID seed ensures this coverage selector does not use truth.
            selected = candidates[:1]
        while len(selected) < k:
            remaining = [item for item in candidates
                         if all(self.pair(item, old)["separated"] for old in selected)]
            if not remaining:
                break
            def coverage(candidate):
                pairs = [self.pair(candidate, old) for old in selected]
                return min(pair["dM"] ** 2 + pair["dS"] ** 2 for pair in pairs)
            selected.append(min(remaining, key=lambda item: (-coverage(item), item.candidate_id)))
        return tuple(selected)

    def summary(self):
        return dict(metric_version=METRIC_VERSION, material_threshold=self.material_threshold,
                    structure_threshold=self.structure_threshold, xor_min_cells=XOR_MIN_CELLS,
                    replacement_mode=self.replacement_mode,
                    submitted_endpoints=self.offers, full_validated_endpoints=len(self._successes),
                    active_combined_representatives=len(self._active),
                    active_material_representatives=len(self.material_representatives),
                    active_structural_representatives=len(self.structural_representatives),
                    frontier=len(self._frontier), failures=len(self._failures),
                    historical_new_world_rewards=sum(event.new_world_reward for event in self._events),
                    posterior_weight=None, represented_mass=None, unmodeled_mass=None,
                    probability_status="CANDIDATE_ONLY", physical_validation="EXTERNAL_RECEIPTS_REQUIRED")
