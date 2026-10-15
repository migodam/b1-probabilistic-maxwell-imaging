"""Research T -> G candidate exchange; reading/writing never launches physics.

New physics is possible only through explicitly supplied, budget-enforcing
target operations. This file does not import a Maxwell solver.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Protocol

import numpy as np

from .structure import (absolute_masks, check_protected, material_grid,
                        topology_summary)

SCHEMA_VERSION = "research-t-g/1.0"
_COST_FIELDS = ("new_full_state_attempts", "factorizations", "forward_rhs", "tangent_rhs",
                "adjoint_rhs", "failed_attempts", "setup_seconds", "wall_seconds")


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(_jsonable(value), ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def pack_complex_per_source(channels: np.ndarray) -> np.ndarray:
    """For each source: all receiver Re, then receiver Im; then next source."""
    z = np.asarray(channels, dtype=np.complex128)
    if z.ndim != 2:
        raise ValueError("channels must be complex (receivers, sources)")
    return np.concatenate((z.real.T, z.imag.T), axis=1).reshape(-1)


def whitened_residual(prediction: np.ndarray, observation: np.ndarray, *, sigma: float) -> np.ndarray:
    """Original A30 convention: sigma is complex RMS, W=sqrt(2)/sigma."""
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("positive finite complex RMS sigma required")
    p, y = np.asarray(prediction, complex), np.asarray(observation, complex)
    if p.shape != y.shape:
        raise ValueError("prediction and original observation have different shapes")
    if not np.all(np.isfinite(p)) or not np.all(np.isfinite(y)):
        raise ValueError("non-finite prediction or observation")
    return (np.sqrt(2.0) / sigma) * pack_complex_per_source(p - y)


class BudgetedTarget(Protocol):
    """Implementations must reserve full states/RHS before every attempt."""
    def predict(self, chi: np.ndarray, *, budget: Any, **kwargs: Any) -> Any: ...
    def validate(self, candidate: Any, *, budget: Any, **kwargs: Any) -> Any: ...
    def correct(self, candidate: Any, *, budget: Any, **kwargs: Any) -> Any: ...


class FrozenPhysicsInterface:
    """G-facing facade over the actual A30 model and Research T's sole book.

    Construction does no physics. The caller supplies the already frozen
    model, original observation, and correction configuration. Changing source
    channels/frequency through an acquisition dictionary is explicitly refused.
    """
    def __init__(self, model: Any, book: Any, observation: np.ndarray, sigma: float,
                 *, target_identity: str, observation_id: str, tau: float,
                 fixed_corrector_config: Any,
                 authorization_identity: str,
                 consumer_role: str = "T",
                 producer_ledger_root: str | Path | None = None):
        if model.book is not book:
            raise ValueError("model and interface must share the unique enforcing book")
        if not authorization_identity:
            raise ValueError("explicit current execution authorization identity required")
        if consumer_role not in {"T", "G"}:
            raise ValueError("declare the executing consumer T or G")
        if consumer_role == "G":
            if producer_ledger_root is None or not hasattr(book, "root"):
                raise ValueError("G must declare the read-only producer ledger and its own authorized book")
            if Path(book.root).resolve() == Path(producer_ledger_root).resolve():
                raise ValueError("G cannot borrow the Research T producer budget")
        self.model, self.book = model, book
        self.observation = np.asarray(observation, complex).copy()
        self.observation.setflags(write=False)
        self.sigma, self.tau = float(sigma), float(tau)
        self.target_identity, self.observation_id = target_identity, observation_id
        self.authorization_identity, self.consumer_role = authorization_identity, consumer_role
        self.corrector_config = fixed_corrector_config
        if self.observation.ndim != 2 or self.observation.shape[1] != 8:
            raise ValueError("frozen target requires eight shared-material sources")
        self.last_cost_receipt = None

    def _budget(self, budget: Any, *, may_launch_physics: bool = True) -> None:
        if budget is not self.book:
            raise ValueError("operation must use the same sole explicit enforcing book")
        if may_launch_physics and hasattr(self.book, "root") and (Path(self.book.root)/"RUN_CLOSED.json").exists():
            raise ValueError("producer/consumer run is closed; require a separate real authorization, not a reset")
        self.book.check()

    def _owned(self, state: Any) -> Any:
        if getattr(state, "_owner", None) is not getattr(self.model, "_owner", object()):
            raise ValueError("foreign/synthetic state has no matching frozen-model provenance")
        if np.shape(state.data) != self.observation.shape:
            raise ValueError("state and observation acquisition disagree")
        return state

    @staticmethod
    def _is_endpoint(candidate: Any) -> bool:
        # ndarray also has .data: it is a buffer, never an owned full endpoint.
        return all(hasattr(candidate, key) for key in ("chi", "current", "field", "data"))

    def _cost(self, first_row: int) -> dict[str, Any]:
        rows = self.book.rows()[first_row:]
        actions = [r for r in rows if r.get("event") == "ACTION"]
        rhs = [r for r in rows if r.get("event") == "RHS"]
        role_totals: dict[str, int] = {}
        for row in rhs:
            role = str(row.get("role", "UNCLASSIFIED"))
            role_totals[role] = role_totals.get(role, 0) + int(row["amount"])
        count_totals: dict[str, int] = {}
        for action in actions:
            for k, v in action.get("counts", {}).items():
                count_totals[k] = count_totals.get(k, 0) + int(v)
        # Exact roles/counts are kept; absent roles are zero only for this
        # freshly delimited operation, never for unobserved historical work.
        receipt = {"new_full_state_attempts": sum(r["amount"] for r in rows if r.get("event") == "STATE"),
                   "complex_full_solver_RHS": sum(r["amount"] for r in rhs),
                   "rhs_by_original_role": role_totals,
                   "counts_by_original_action": count_totals,
                   "failed_actions": sum(r.get("status") == "FAIL" for r in actions),
                   "summed_action_wall_seconds": sum(r.get("wall_seconds", 0.) for r in actions),
                   "action_wall_note": "nested action spans can overlap; authoritative cumulative wall is book JOB_END",
                   "authorization_identity": self.authorization_identity,
                   "consumer_role": self.consumer_role,
                   "ledger_row_range": [first_row, first_row+len(rows)]}
        self.last_cost_receipt = receipt
        return receipt

    def predict(self, chi: np.ndarray, *, budget: Any,
                acquisition: dict[str, Any] | None = None) -> Any:
        self._budget(budget)
        allowed = {"target_identity": self.target_identity, "observation_id": self.observation_id}
        if acquisition is not None and (set(acquisition)-set(allowed)
                or any(allowed[k] != v for k, v in acquisition.items())):
            raise ValueError("new acquisition physics is outside the frozen target; no silent solver change")
        first = len(self.book.rows())
        try:
            return self.model.full_state(np.asarray(chi, complex).reshape(-1), purpose="T_G_explicit_full_prediction")
        finally:
            self._cost(first)

    def validate(self, candidate: Any, *, budget: Any) -> dict[str, Any]:
        self._budget(budget, may_launch_physics=not self._is_endpoint(candidate))
        first = len(self.book.rows())
        # A caller-supplied material needs one paid matching full endpoint;
        # an actual owned endpoint can be replayed without repeating physics.
        state = self._owned(candidate) if self._is_endpoint(candidate) else self.predict(candidate, budget=budget)
        x = np.asarray(state.chi, complex)
        residual = float(np.linalg.norm(whitened_residual(state.data, self.observation, sigma=self.sigma)))
        defect = float(np.linalg.norm(state.current-x[:, None]*state.field) / max(np.linalg.norm(state.current), 1e-300))
        legal = bool(np.all((x.real >= -.5-1e-12) & (x.real <= 4+1e-12)
                            & (x.imag >= -1e-12) & (x.imag <= .5+1e-12)))
        backward = float(state.residual)
        numerical_checks = bool(legal and defect < 1e-10 and backward < 1e-10)
        passed = numerical_checks and residual <= self.tau
        return {"status": "FULL_VALIDATED_NUMERICAL" if passed else "REJECTED",
                "state": state, "matching_full_endpoint": True,
                "checks_passed": numerical_checks,
                "observed_residual_norm": residual, "tau": self.tau,
                "material_valid": legal, "constitutive_defect": defect,
                "constitutive_status": "CHECKED_FROM_ACTUAL_OWNED_STATE",
                "backward_residual": backward, "output_error_upper_bound": None,
                "bound_provenance": None,
                "updated_cost_receipt": self._cost(first)}

    def correct(self, candidate: Any, *, budget: Any,
                protected_geometry: dict[str, np.ndarray] | None = None,
                cache: Any = None) -> Any:
        self._budget(budget)
        from .corrector import correct
        first = len(self.book.rows())
        state = self._owned(candidate) if self._is_endpoint(candidate) else self.predict(candidate, budget=budget)
        try:
            return correct(self.model, state, self.observation, self.sigma, self.book,
                           self.corrector_config, protected_geometry, cache=cache)
        finally:
            self._cost(first)

    def compensation_score(self, candidate: Any, *, budget: Any,
                           protected_geometry: dict[str, np.ndarray] | None = None,
                           cache: Any = None) -> Any:
        self._budget(budget)
        from .corrector import bounded_compensation_score
        first = len(self.book.rows())
        state = self._owned(candidate) if self._is_endpoint(candidate) else self.predict(candidate, budget=budget)
        try:
            return bounded_compensation_score(self.model, state, self.observation,
                                               self.sigma, self.corrector_config,
                                               protected_geometry, cache=cache)
        finally:
            self._cost(first)


def new_measurement(target: BudgetedTarget, chi: np.ndarray, *, budget: Any,
                    acquisition: dict[str, Any]) -> Any:
    """Explicit paid operation, never called by bundle I/O or schema checks."""
    if budget is None:
        raise ValueError("new measurement requires an explicit enforcing budget")
    return target.predict(chi, budget=budget, acquisition=acquisition)


def validate_candidate(target: BudgetedTarget, candidate: Any, *, budget: Any,
                       **kwargs: Any) -> Any:
    if budget is None:
        raise ValueError("physics validation requires an explicit enforcing budget")
    return target.validate(candidate, budget=budget, **kwargs)


def correct_candidate(target: BudgetedTarget, candidate: Any, *, budget: Any,
                      **kwargs: Any) -> Any:
    if budget is None:
        raise ValueError("material correction requires an explicit enforcing budget")
    return target.correct(candidate, budget=budget, **kwargs)


def write_bundle(destination: str | Path, *, candidate_id: str,
                 parent_id: str | None, source_revision: str,
                 chi: np.ndarray, observation: np.ndarray,
                 prediction: np.ndarray | None, coordinates: np.ndarray,
                 mass: np.ndarray | float, sigma: float, tau: float,
                 physics: dict[str, Any], observation_id: str, noise_id: str,
                 legal_domain: dict[str, Any], validation: dict[str, Any],
                 cost: dict[str, Any], lineage: dict[str, Any],
                 core: np.ndarray | None = None, guard: np.ndarray | None = None,
                 optional_state: dict[str, np.ndarray] | None = None,
                 nuisance: dict[str, Any] | None = None,
                 minimum_feature_size: float | None = None,
                 minimum_feature_units: str | None = None,
                 artifact_kind: str = "REAL_ENDPOINT_CANDIDATE",
                 shape: tuple[int, int] = (24, 24)) -> Path:
    """Serialize a real caller-supplied endpoint, without recomputation.

    This function can also preserve a rejected endpoint, but cannot promote it
    to validated merely because the saved residual is small. Matching full
    endpoint provenance and numerical checks are supplied by the real target.
    """
    if not candidate_id or not source_revision or not observation_id or not noise_id:
        raise ValueError("nonempty candidate, source revision, observation and noise identities required")
    if artifact_kind not in {"REAL_ENDPOINT_CANDIDATE", "INTERFACE_FIXTURE_NOT_A_PHYSICAL_CANDIDATE"}:
        raise ValueError("undeclared artifact kind")
    missing = [k for k in _COST_FIELDS if k not in cost]
    if missing:
        raise ValueError(f"explicit cost receipt missing {missing}; unknown cost is not zero")
    if any(cost[k] is not None and (not np.isfinite(cost[k]) or cost[k] < 0) for k in _COST_FIELDS):
        raise ValueError("negative/nonfinite cost receipt")
    value = material_grid(chi, shape)
    xy = np.asarray(coordinates, float)
    if xy.shape != (value.size, 2) or not np.all(np.isfinite(xy)):
        raise ValueError("actual grid coordinates required as (grid_cells,2)")
    masses = np.full(value.size, float(mass)) if np.asarray(mass).size == 1 else np.asarray(mass, float).reshape(-1)
    if masses.shape != (value.size,) or not np.all(np.isfinite(masses)) or np.any(masses <= 0):
        raise ValueError("actual positive cell masses required")
    y = np.asarray(observation, complex)
    if y.ndim != 2 or y.shape[1] != 8:
        raise ValueError("frozen A30 acquisition requires (receivers,8 sources)")
    if not np.all(np.isfinite(y)):
        raise ValueError("non-finite original observation")
    if tau <= 0 or not np.isfinite(tau) or sigma <= 0 or not np.isfinite(sigma):
        raise ValueError("finite frozen tau and complex RMS sigma required")
    status = validation.get("status", "NOT_VALIDATED")
    allowed_status = {"NOT_VALIDATED", "REJECTED", "FULL_VALIDATED_NUMERICAL", "FULL_VALIDATED_WITH_BOUND"}
    if status not in allowed_status:
        raise ValueError(f"unknown validation level: {status}")
    p = None if prediction is None else np.asarray(prediction, complex)
    residual = None if p is None else float(np.linalg.norm(whitened_residual(p, y, sigma=sigma)))
    numerical_claim = status in {"FULL_VALIDATED_NUMERICAL", "FULL_VALIDATED_WITH_BOUND"}
    if numerical_claim:
        required_physics = ("backend", "dimension", "field_model", "harmonic_sign",
                            "material_convention", "dispersion_law")
        if any(physics.get(key) is None for key in required_physics):
            raise ValueError("validated bundle lacks an actual declared physics convention")
        if p is None or validation.get("matching_full_endpoint") is not True:
            raise ValueError("validated bundle needs matching real full endpoint provenance")
        if not validation.get("checks_passed", False) or residual > tau + 1e-10 * max(tau, 1):
            raise ValueError("validated bundle contradicts checks or frozen observed-data residual")
    error_bound = validation.get("output_error_upper_bound")
    if status == "FULL_VALIDATED_WITH_BOUND":
        if error_bound is None or not validation.get("bound_provenance") or residual + error_bound > tau:
            raise ValueError("WITH_BOUND requires a justified output bound within tau")
    cc = None if core is None else np.asarray(core, bool).reshape(shape)
    gg = None if guard is None else np.asarray(guard, bool).reshape(shape)
    if (cc is None) != (gg is None):
        raise ValueError("core and guard must be supplied together, or both absent")
    protection = None if cc is None else check_protected(value, cc, gg, shape=shape)
    if numerical_claim and protection is not None and not protection["preserved"]:
        raise ValueError("validated protected edit does not preserve core/guard")
    if optional_state is not None and (not optional_state or any(not np.all(np.isfinite(v)) for v in optional_state.values())):
        raise ValueError("optional state must contain actual finite arrays, never unavailable-field placeholders")
    out = Path(destination)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"refusing to overwrite existing candidate bundle {out}")
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "material.npz", chi=value, coordinates=xy,
                        cell_mass=masses, grid_shape=np.asarray(shape))
    np.savez_compressed(out / "observation.npz", y=y, sigma_complex_rms=sigma,
                        whitening_scale=np.sqrt(2.0) / sigma)
    if p is not None:
        delta = p - y
        block_norms = np.sqrt(2.0) / sigma * np.linalg.norm(delta, axis=0)
        np.savez_compressed(out / "prediction.npz", channels=p,
                            per_source_residual_norm=block_norms,
                            packed_whitened_residual=whitened_residual(p, y, sigma=sigma))
    masks = absolute_masks(value, shape=shape)
    arrays = {"absolute_air": masks["air"], "absolute_material": masks["material"],
              "passive": masks["passive"], "guard_material": masks["guard_material"]}
    if cc is not None:
        arrays.update(protected_core=cc, material_guard=gg)
    np.savez_compressed(out / "structure.npz", **arrays)
    if optional_state is not None:
        np.savez_compressed(out / "optional_state.npz", **optional_state)
    _write_json(out / "cost.json", cost)
    _write_json(out / "lineage.json", lineage)
    _write_json(out / "legal_domain.json", legal_domain)
    _write_json(out / "distance_definition.json", {
        "material_l2": "sqrt(sum(cell_mass*abs(chi_a-chi_b)^2)/sum(cell_mass))/1.2",
        "material_normalization": 1.2,
        "absolute_air": "abs(chi)<=.10; chi=epsilon_r-1",
        "phase_threshold_sensitivity": [.075, .10, .125],
        "significant_phase_changed_cells": 6,
        "structure_components": ["absolute-phase symmetric difference", "enclosed cavity symmetric difference", "2D Betti changes"],
        "material_foreground_connectivity": 4, "air_connectivity": 8,
        "mesh_independent_certificate": False})
    physics_out = dict(physics)
    physics_out.update(geometry_ref="material.npz:coordinates", mass_ref="material.npz:cell_mass")
    vout = dict(validation)
    vout.update(status=status, prediction_ref="prediction.npz:channels" if p is not None else None,
                observed_residual_norm=residual,
                output_error_upper_bound=error_bound,
                bound_provenance=validation.get("bound_provenance"),
                constitutive_defect=validation.get("constitutive_defect"),
                constitutive_status=validation.get("constitutive_status", "NOT_EVALUATED"),
                backward_residual=validation.get("backward_residual"))
    manifest = {"schema_version": SCHEMA_VERSION, "artifact_kind": artifact_kind,
                "candidate_id": candidate_id, "parent_id": parent_id, "source_revision": source_revision,
                "physics": physics_out,
                "data": {"observation_id": observation_id, "noise_id": noise_id,
                         "packing": "each source: receiver Re then receiver Im; concatenate 8 sources",
                         "whitening_ref": "observation.npz:whitening_scale",
                         "sigma_semantics": "complex RMS; real and imaginary SD=sigma/sqrt(2)",
                         "observation_ref": "observation.npz:y", "tau": float(tau),
                         "history_ids": [observation_id]},
                "material": {"array_ref": "material.npz:chi", "shared_across_illuminations": True,
                             "legal_domain_ref": "legal_domain.json"},
                "nuisance": nuisance or {"value": {"real_gain": 1.0}, "legal_domain": {"real_gain": [1.0, 1.0]},
                                         "shared_history_scope": [observation_id], "status": "known/frozen"},
                "structure": {"absolute_air_definition": "abs(chi)<=0.10; chi=epsilon_r-1",
                              "minimum_feature_size": minimum_feature_size, "units": minimum_feature_units,
                              "minimum_feature_status": "DECLARED" if minimum_feature_size is not None else "NOT_DECLARED",
                              "core_ref": "structure.npz:protected_core" if cc is not None else None,
                              "guard_ref": "structure.npz:material_guard" if gg is not None else None,
                              "protection": protection, "betti_numbers": topology_summary(value, shape=shape),
                              "distance_definition_ref": "distance_definition.json"},
                "validation": vout, "cost": cost,
                "lineage": {"array_ref": "lineage.json", "proposal_kind": lineage.get("proposal_kind"),
                            "protected_structure_preserved": protection["preserved"] if protection else None,
                            "verified_path_ref": lineage.get("verified_path_ref")},
                "optional_state_ref": "optional_state.npz" if optional_state is not None else None,
                "optional_state_unavailable_reason": None if optional_state is not None else "not exposed in caller-supplied real state",
                "weight_semantics": "candidate_only", "probability_weight": None,
                "unexplored_structure_status": "UNQUANTIFIED; finite candidate archive is not complete",
                "outside_model_support_mass": None}
    _write_json(out / "manifest.json", manifest)
    return out


def validate_bundle(path: str | Path, *, allow_fixture: bool = False) -> dict[str, Any]:
    """Pure local schema/semantic replay. No physics imports or solver calls."""
    folder = Path(path)
    m = json.loads((folder / "manifest.json").read_text())
    errors: list[str] = []
    if m.get("schema_version") != SCHEMA_VERSION:
        errors.append("unsupported schema version")
    if m.get("artifact_kind") != "REAL_ENDPOINT_CANDIDATE" and not (
            allow_fixture and m.get("artifact_kind") == "INTERFACE_FIXTURE_NOT_A_PHYSICAL_CANDIDATE"):
        errors.append("not a formal real endpoint candidate")
    if m.get("weight_semantics") != "candidate_only" or m.get("probability_weight") is not None:
        errors.append("candidate residuals must not be interpreted as posterior probability")
    try:
        material = np.load(folder / "material.npz", allow_pickle=False)
        chi = material_grid(material["chi"])
        if material["coordinates"].shape != (576, 2) or not np.isfinite(material["coordinates"]).all():
            errors.append("actual 24x24 cell coordinate geometry missing")
        if material["cell_mass"].shape != (576,) or not np.isfinite(material["cell_mass"]).all() or np.any(material["cell_mass"] <= 0):
            errors.append("cell mass missing/nonpositive")
        obs = np.load(folder / "observation.npz", allow_pickle=False)
        if not m["material"]["shared_across_illuminations"]:
            errors.append("each source cannot have a separate material")
        if obs["y"].ndim != 2 or obs["y"].shape[1] != 8:
            errors.append("wrong frozen source count")
        sigma = float(obs["sigma_complex_rms"])
        if not np.isclose(float(obs["whitening_scale"]), np.sqrt(2) / sigma, rtol=1e-14):
            errors.append("whitening convention changed")
        legal = json.loads((folder / "legal_domain.json").read_text())
        re = legal.get("real_chi", [-0.5, 4.0])
        im = legal.get("imag_chi", [0.0, 0.5])
        material_legal = bool(np.all((chi.real >= re[0] - 1e-12) & (chi.real <= re[1] + 1e-12))
                              and np.all((chi.imag >= im[0] - 1e-12) & (chi.imag <= im[1] + 1e-12)))
        status = m["validation"]["status"]
        if status.startswith("FULL_VALIDATED"):
            if not m["validation"].get("matching_full_endpoint") or not m["validation"].get("checks_passed"):
                errors.append("matching full endpoint provenance/checks missing")
            if any(m["physics"].get(k) is None for k in ("backend", "dimension", "field_model", "harmonic_sign", "material_convention", "dispersion_law")):
                errors.append("full endpoint physics convention incomplete")
        if status.startswith("FULL_VALIDATED") and not material_legal:
            errors.append("validated material violates declared legal domain")
        prediction_ref = m["validation"]["prediction_ref"]
        if prediction_ref is not None:
            pred = np.load(folder / "prediction.npz", allow_pickle=False)
            packed = whitened_residual(pred["channels"], obs["y"], sigma=sigma)
            norm = float(np.linalg.norm(packed))
            if not np.allclose(packed, pred["packed_whitened_residual"], rtol=1e-12, atol=1e-12):
                errors.append("packed residual mismatch")
            if not np.isclose(norm, m["validation"]["observed_residual_norm"], rtol=1e-12, atol=1e-12):
                errors.append("saved observed-data residual mismatch")
            if status.startswith("FULL_VALIDATED") and norm > m["data"]["tau"] + 1e-10:
                errors.append("full validated candidate is data-incompatible")
        elif status.startswith("FULL_VALIDATED"):
            errors.append("full validated candidate lacks saved full prediction")
        structure = np.load(folder / "structure.npz", allow_pickle=False)
        if not np.array_equal(structure["absolute_air"], np.abs(chi) <= 0.1):
            errors.append("absolute air mask is not chi-to-vacuum threshold")
        if "protected_core" in structure:
            protected = check_protected(chi, structure["protected_core"], structure["material_guard"])
            if status.startswith("FULL_VALIDATED") and not protected["preserved"]:
                errors.append("protected cavity/guard not retained")
        costs = json.loads((folder / "cost.json").read_text())
        if any(k not in costs for k in _COST_FIELDS):
            errors.append("incomplete real cost receipt")
        if any(costs.get(k) is not None and (not np.isfinite(costs[k]) or costs[k] < 0) for k in _COST_FIELDS):
            errors.append("negative/nonfinite cost receipt")
        if not (folder / m["structure"]["distance_definition_ref"]).is_file():
            errors.append("invented distance definition reference")
        if m.get("optional_state_ref") is not None and not (folder / m["optional_state_ref"]).is_file():
            errors.append("invented optional state reference")
        nuisance = m["nuisance"].get("value", {})
        if nuisance.get("real_gain") != 1.0:
            errors.append("known shared gain=1 frozen target changed")
    except (KeyError, ValueError, OSError, TypeError) as exc:
        errors.append(f"missing or inconsistent bundle content: {exc}")
    return {"schema_version": SCHEMA_VERSION, "candidate_id": m.get("candidate_id"),
            "valid_bundle_semantics": not errors, "errors": errors,
            "new_physics_attempts": 0, "new_full_solver_rhs": 0,
            "scope": "saved-array/schema replay only; does not independently re-solve Maxwell"}
