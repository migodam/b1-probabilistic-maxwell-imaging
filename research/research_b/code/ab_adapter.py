"""Explicit Research A 1.0.0 -> Research B 1.0 import, without solver calls.

The A and B schemas are different.  In particular A does not encode material
units, observation IDs or a gain-sharing graph.  An independently supplied
``import_context`` is therefore required; missing facts are never guessed.
``adapt_a_v1`` returns an in-memory bundle, ``export_bundle`` writes an independent
B package and sidecar.  No A asset, manifest or historical ledger is modified.

The input's validation records are checked for consistency, not recertified by
Maxwell.  A numerical fixture is never evidence that a real A package exists.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any

import numpy as np

ADAPTER_VERSION = "research-a-1.0.0-to-b-1.0-context-1"
CONTRACT_VERSION = "research-b-import-contract-1"
A_EXECUTION_PROFILE = "research-a-native-1.0.0-execution-no-digest-v1"
A_COMMON_EXCHANGE_VERSION = "research-a-common-exchange/1.0.0"
PACKINGS = {"real_only", "all_real_then_all_imag", "interleaved_real_imag"}
MATERIAL_UNITS = {"relative_permittivity": "1", "electric_susceptibility": "1",
                  "absolute_permittivity": "F/m"}
IDENTITY_KEYS = ("hash", "solver_id", "backend_commit", "material_law_id",
                 "time_convention", "dispersion_model", "self_term_convention",
                 "polarizability_map", "boundary_conditions")
TRUTH_KEYS = {"truth", "ground_truth", "truth_material", "truth_mae", "best_of_k_truth"}


class ContractError(ValueError):
    """An import would mix targets, measurements or probability semantics."""


def require(condition: Any, message: str) -> None:
    if not condition:
        raise ContractError(message)


def _walk(value: Any):
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from _walk(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk(child)


def _finite(value: Any, name: str) -> np.ndarray:
    array = np.asarray(value)
    require(array.dtype.kind in "biufc", f"{name}: numeric data required")
    require(np.all(np.isfinite(array)), f"{name}: non-finite data")
    return array


def _load_asset(asset: dict, root: Path) -> np.ndarray | dict | list:
    path = (root / asset["path"]).resolve()
    require(not Path(asset["path"]).is_absolute() and path.is_relative_to(root),
            "Asset must remain inside the read-only source package")
    require(path.is_file(), f"Missing source asset: {asset['path']}")
    if asset["format"] == "json":
        return json.loads(path.read_text(encoding="utf-8"))
    require(asset["format"] in {"npy", "npz"}, "This minimal importer supports npy/npz/json only")
    if asset["format"] == "npz":
        require(asset.get("key"), "npz asset requires an explicit key")
        with np.load(path, allow_pickle=False) as data:
            require(asset["key"] in data.files, "Missing npz array key")
            result = np.array(data[asset["key"]], copy=True)
    else:
        require(asset.get("key") is None, "npy asset cannot have a key")
        result = np.load(path, allow_pickle=False)
    require(list(result.shape) == asset["shape"], "Asset shape declaration differs from data")
    require(np.dtype(result.dtype) == np.dtype(asset["dtype"]), "Asset dtype differs from data")
    return _finite(result, asset["path"])


def _real_vector(array: np.ndarray, packing: str) -> np.ndarray:
    flat = _finite(array, "measurement").reshape(-1)
    if not np.iscomplexobj(flat):
        return flat.astype(np.float64)
    require(packing != "real_only", "Complex measurement cannot use real_only packing")
    if packing == "all_real_then_all_imag":
        return np.concatenate((flat.real, flat.imag)).astype(np.float64)
    return np.column_stack((flat.real, flat.imag)).reshape(-1).astype(np.float64)


def _ref(name: str, data: Any, units: str, provenance: str) -> dict:
    array = np.asarray(data)
    return {"uri": f"bundle://{name}", "format": "json" if array.dtype.kind == "O" else "npy",
            "shape": [] if array.dtype.kind == "O" else list(array.shape),
            "dtype": "json" if array.dtype.kind == "O" else str(array.dtype),
            "units": units, "provenance_id": provenance}


def _uncertainty(source: str) -> dict:
    return {"status": "unquantified", "value": None, "metric": "not_assigned",
            "scope": "Declared source records, not a new Maxwell validation",
            "confidence": None, "source": source}


def _unknown_cost() -> dict:
    return {"measurement_seconds": None, "switching_seconds": None,
            "energy_joules": None, "financial_cost": None, "currency": None,
            "compute_seconds": None, "maxwell_attempts": 0, "factorizations": 0,
            "rhs": 0, "setup_amortization_count": 1, "cost_status": "not_available"}


def _validate_action(action: dict, context: dict, n_groups: set[str]) -> None:
    require(action["id"] and action["noise"]["packing"] == context["noise_packing"],
            "Action noise packing differs from imported measurements")
    require(set(action["nuisance_group_ids"]) == n_groups,
            "Action must share the declared persistent nuisance groups")
    require(action["cost"]["maxwell_attempts"] >= 0, "Negative action cost")


def adapt_a_v1(package: dict, package_root: str | Path, import_context: dict) -> dict:
    """Import actual A data, requiring explicit identity, units/history/sharing.

    Required context keys are demonstrated in tests and documented in the
    delegated artifact.  The target identity is an opaque ID comparison; no
    content digest is computed.  All source asset reads are confined to root.
    """
    p, ctx = deepcopy(package), deepcopy(import_context)
    root = Path(package_root).resolve()
    require(ctx.get("adapter_version") == ADAPTER_VERSION, "Unsupported import-context version")
    require(p.get("version") == "1.0.0", "Unsupported Research A native version")
    require(p.get("status") == "handoff_ready", "A package is not handoff_ready")
    for node in _walk([p, ctx]):
        require(not TRUTH_KEYS.intersection(node), "Truth-derived fields cannot enter online handoff")
        for value in node.values():
            if isinstance(value, float):
                require(np.isfinite(value), "Non-finite manifest value")
    target = p["physics_target"]
    for key in IDENTITY_KEYS:
        require(key in ctx["expected_physics"] and target[key] == ctx["expected_physics"][key],
                f"Physics target mismatch: {key}")
    require(ctx["material_quantity"] in MATERIAL_UNITS and
            ctx["material_units"] == MATERIAL_UNITS[ctx["material_quantity"]],
            "Unknown or inconsistent material quantity/units")
    require(ctx["unit_evidence"] and ctx["observation_units"] == p["noise"]["units"],
            "Observation units mismatch or missing external unit evidence")
    require(ctx["noise_packing"] in PACKINGS and
            ctx["noise_packing"] == p["noise"]["realification_order"],
            "Unsupported or mismatched noise packing; explicit conversion required")
    require(ctx["frequency_units"] == p["acquisition"]["frequency_units"], "Frequency units mismatch")
    require(ctx["frequency_units"] in {"kR", "Hz"}, "Frequency unit conversion is not implemented")
    require(ctx["native_acquisition_id"] == p["acquisition"]["id"], "Native acquisition/history ID mismatch")
    history = ctx["history"]
    obs_ids = [record["observation_id"] for record in history]
    noise_ids = [record["noise_realization_id"] for record in history]
    require(history and len(set(obs_ids)) == len(obs_ids), "Duplicate observation ID")
    require(len(set(noise_ids)) == len(noise_ids), "Duplicate noise realization ID")
    action_ids = [record["action"]["id"] for record in history]
    require(len(set(action_ids)) == len(action_ids), "Duplicate historical action ID")
    frequency_key = "frequency_kR" if ctx["frequency_units"] == "kR" else "frequency_hz"
    native_frequencies = set(p["acquisition"]["frequencies"])
    history_frequencies = {r["action"][frequency_key] for r in history if r["action"][frequency_key] is not None}
    require(native_frequencies == history_frequencies,
            "Context must retain every already conditioned native frequency")
    y = _real_vector(_load_asset(p["observation"], root), ctx["noise_packing"])
    require(len(y) == p["acquisition"]["n_real_measurements"], "Real measurement dimension mismatch")
    indices = [i for record in history for i in record["real_channel_indices"]]
    require(len(indices) == len(set(indices)) and sorted(indices) == list(range(len(y))),
            "Historical channels must form one disjoint full observation, counted once")
    covariance = np.asarray(_load_asset(p["noise"]["covariance"], root), dtype=np.float64)
    require(covariance.shape == (len(y), len(y)), "Covariance dimension mismatch")
    require(np.allclose(covariance, covariance.T, rtol=1e-12, atol=1e-14), "Covariance is not symmetric")
    try:
        chol = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ContractError("Covariance must be strictly positive definite") from exc

    nuisance = ctx["nuisance"]
    require(nuisance["sharing_rule"] == "one_parameter_per_world_across_history",
            "Gain/nuisance must persist within each material world across history")
    require(nuisance["different_worlds_may_use_different_values"] is True,
            "Different material worlds may have different legal gain values")
    require(nuisance["source_evidence"] and nuisance["prior_definition"], "Missing nuisance evidence/prior")
    parameters = nuisance["parameter_names"]
    require(len(parameters) == len(set(parameters)), "Duplicate nuisance parameter")
    require(set(nuisance["group_by_parameter"]) == set(parameters), "Incomplete nuisance sharing map")
    groups = set(nuisance["group_by_parameter"].values())
    for name in parameters:
        bound = nuisance["bounds"][name]
        require(len(bound) == 2 and np.all(np.isfinite(bound)) and bound[0] < bound[1],
                "Invalid legal nuisance bounds")
    for record in history:
        _validate_action(record["action"], ctx, groups)
    available = ctx["available_actions"]
    available_ids = [a["id"] for a in available]
    require(len(available_ids) == len(set(available_ids)), "Duplicate available action ID")
    require(not set(available_ids).intersection(action_ids), "Already observed action offered as new")
    for action in available:
        _validate_action(action, ctx, groups)
        require(not set(action["constituent_ids"]).intersection(action_ids + obs_ids),
                "New action reuses an already conditioned observation")
        if action["kind"] == "repeat":
            require(action["new_independent_noise"], "Repeat requires new independent noise")
        if action["kind"] == "frequency" and action[frequency_key] in native_frequencies:
            raise ContractError("Already conditioned frequency must be an explicit fresh-noise repeat")

    candidates = p["candidates"]
    ids = [c["id"] for c in candidates]
    require(candidates and len(ids) == len(set(ids)) and len(ids) == p["valid_k"] <= p["requested_k"],
            "Invalid candidate IDs/count")
    require(p["base_reconstruction"]["full_validation_passed"], "Base reconstruction lacks full validation")
    base_material = _load_asset(p["base_reconstruction"]["material_field"], root)
    require(isinstance(base_material, np.ndarray) and list(base_material.shape) == ctx["material_shape"],
            "Base material field differs from the declared support")
    event_ids = sorted({e for c in candidates for e in c["structure"]["event_ids"]})
    event_map = {event: i for i, event in enumerate(event_ids)}
    assets: dict[str, Any] = {"observations": y, "covariance": covariance, "base_material": base_material,
        "sharing_graph": {"parameter_groups": nuisance["group_by_parameter"],
                          "observation_ids": obs_ids, "scope": nuisance["sharing_rule"]},
        "nuisance_model": nuisance, "nuisance_bounds": nuisance["bounds"],
        "source_cost": p["total_cost"]}
    b_candidates = []
    native_validation = []
    for c in candidates:
        validation, fit = c["validation"], c["fit"]
        require(validation["passed"] is True and validation["scope"] == "full_discrete_maxwell"
                and validation["all_illuminations_same_material"] is True
                and validation["material_constraints_passed"] is True, "Candidate lacks full Maxwell validation")
        require(validation["solver_id"] == target["solver_id"] and
                validation["physics_target_hash"] == target["hash"] and
                validation["material_law_id"] == target["material_law_id"], "Candidate target/solver mismatch")
        require(fit["noise_model_id"] == p["noise"]["id"] and fit["fit_acceptance"] is True,
                "Candidate fit/noise model mismatch")
        require(fit["nuisance_policy"] != "bounded_profile" or
                fit["likelihood_semantics"] != "proper_marginal", "Profile likelihood is not a proper marginal")
        for name in parameters:
            value = fit["nuisance_parameters"][name]
            lo, hi = nuisance["bounds"][name]
            require(np.isscalar(value) and np.isfinite(value) and lo <= value <= hi,
                    "Nuisance must be one legal scalar per world, not independently fitted per frequency")
        material = _load_asset(c["material_field"], root)
        require(isinstance(material, np.ndarray) and list(material.shape) == ctx["material_shape"],
                "Material field shape does not match the declared unified support")
        if np.iscomplexobj(material):
            sign = 1 if target["time_convention"] == "exp(-i omega t)" else -1
            require(np.all(sign * material.imag >= -1e-12), "Material violates passive sign convention")
        prediction = _real_vector(_load_asset(c["full_predicted_measurement"], root), ctx["noise_packing"])
        require(prediction.shape == y.shape, "Prediction/observation packing or shape mismatch")
        residual = np.linalg.solve(chol, prediction - y)
        l2 = float(np.linalg.norm(residual))
        require(np.isclose(l2, fit["whitened_residual_l2"], rtol=1e-6, atol=1e-9)
                and np.isclose(l2 * l2, fit["squared_residual"], rtol=1e-6, atol=1e-10),
                "Saved predictions and declared joint-history whitened fit disagree")
        require(l2 <= fit["threshold_l2"], "Retained candidate fails its declared data-fit threshold")
        mkey, pkey = f"material_{len(b_candidates)}", f"prediction_{len(b_candidates)}"
        assets[mkey], assets[pkey] = material, prediction
        b_candidates.append({"id": c["id"],
            "material": _ref(mkey, material, ctx["material_units"], p["package_id"]),
            "structure_description": json.dumps(c["structure"], sort_keys=True),
            "event_ids": [event_map[e] for e in c["structure"]["event_ids"]],
            "base_solution": bool(np.array_equal(material, base_material) and
                fit["nuisance_parameters"] == p["base_reconstruction"]["nuisance_parameters"]),
            "full_maxwell_predictions": [_ref(pkey, prediction, ctx["observation_units"], p["package_id"])],
            "prediction_history_id": ctx["data_history_id"],
            "material_distance_to_base": c["material_distance_to_base"],
            "whitened_distance_to_data": l2, "distance_definitions": c["distance_definition_id"],
            "compensation_constraints": json.dumps(c["compensation"], sort_keys=True),
            "physical_feasibility": "Source full-discrete-Maxwell validation record; not recertified by importer",
            "rom": {"used": False, "model_id": None, "valid_material_domain": ctx["material_domain"],
                    "valid_action_ids": [], "response_error": _uncertainty(p["package_id"]),
                    "refresh_rule": "No ROM imported; full predictions only", "setup_and_refresh_cost": _unknown_cost()}})
        native_validation.append({"candidate_id": c["id"], "validation": validation,
                                  "nuisance_parameters": fit["nuisance_parameters"]})

    reliability = p["weight_reliability"]
    kind = {"candidate_only": "none", "learned_proposal": "proposal", "reference": "independent_reference",
            "corrected": "importance_corrected"}.get(reliability)
    require(kind is not None, "Unknown source weight reliability")
    weights = p["weights"]
    allow_bayes = False
    normalized_ref = None
    normalization_status = "No posterior weights; candidates are a finite discovered set"
    if reliability in {"candidate_only", "learned_proposal"}:
        require(weights is None and all(c["posterior_weight"] is None for c in candidates),
                "Candidate/proposal values cannot be promoted to posterior probabilities")
    else:
        require(weights is not None, "Reference/corrected declaration lacks weights")
        w = _finite(weights["values"], "weights").astype(np.float64)
        require(w.shape == (len(candidates),) and np.all(w >= 0), "Invalid probability weights")
        require(np.allclose(w, [c["posterior_weight"] for c in candidates], rtol=1e-9, atol=1e-12),
                "Candidate/global weight mismatch")
        require(np.sum(w) <= 1 + 1e-9, "Probability mass exceeds one")
        assets["source_weights"] = w
        full = (weights["normalization_domain"] == "declared_full_support"
                and np.isclose(w.sum(), 1, rtol=1e-9, atol=1e-12)
                and weights["represented_mass"] == 1 and weights["unrepresented_mass"] == 0)
        if reliability == "reference":
            evidence = ctx.get("reference_evidence", {})
            require(evidence.get("target_id") == ctx["target_id"] and
                    evidence.get("data_history_id") == ctx["data_history_id"] and
                    evidence.get("independent_reference") is True and evidence.get("provenance"),
                    "Independent reference evidence matching target/history is required")
        else:
            require(weights["correction_method"] and
                    all(c["proposal_log_density"] is not None for c in candidates),
                    "Corrected weights require an explicit correction and proposal density")
        allow_bayes = bool(full and p["next_measurement_interface"]["posterior_information_gain_allowed"])
        normalization_status = "Declared full support normalized" if full else "Partial/pool reference; not a full posterior"
        if full:
            normalized_ref = _ref("source_weights", w, "probability", p["package_id"])
    logrefs = {}
    for name, native in (("log_prior", "proper_log_prior"), ("log_proposal", "proposal_log_density")):
        values = [c[native] for c in candidates]
        proper = name != "log_prior" or all(c["prior_score_is_normalized_density"] for c in candidates)
        if proper and all(v is not None for v in values):
            assets[name] = _finite(values, name).astype(np.float64)
            logrefs[name] = _ref(name, assets[name], "log_density", p["package_id"])
        else:
            logrefs[name] = None

    noise = {"packing": ctx["noise_packing"],
             "covariance": _ref("covariance", covariance, ctx["observation_units"] + "^2", p["package_id"]),
             "cross_history_covariance": _ref("covariance", covariance, ctx["observation_units"] + "^2", p["package_id"]),
             "real_channel_variance_convention": ctx["variance_convention"],
             "dwell_scaling": ctx["dwell_scaling"], "reference_id": p["noise"]["id"]}
    b = {"schema_version": "research-b-mhp-1.0", "package_id": p["package_id"] + ":research-b-import-1",
         "created_at": p["created_utc"], "producer": "Research B versioned A import adapter",
         "source_revision": target["backend_commit"], "target_id": ctx["target_id"],
         "data_history_id": ctx["data_history_id"],
         "model": {"physical_equations": "Full discrete vector Maxwell",
                   "solver_id": target["solver_id"], "discretization_id": ctx["discretization_id"],
                   "material_domain": ctx["material_domain"], "material_metric": ctx["material_metric"],
                   "dispersion_model": target["dispersion_model"], "time_convention": target["time_convention"],
                   "reference_frequency": p["acquisition"]["frequencies"][0]},
         "data": {"observations": _ref("observations", y, ctx["observation_units"], p["package_id"]),
                  "acquisition_history": [r["action"] for r in history], "noise": noise},
         "nuisance": {"parameter_names": parameters, "prior_definition": nuisance["prior_definition"],
                      "sharing_graph": _ref("sharing_graph", assets["sharing_graph"], "1", p["package_id"]),
                      "posterior_representation": "None; per-world source fits plus legal persistent compensation",
                      "joint_particles_or_conditional_model": _ref("nuisance_model", nuisance, "1", p["package_id"]),
                      "drift_transition": None,
                      "bounds": _ref("nuisance_bounds", assets["nuisance_bounds"], "1", p["package_id"])},
         "candidates": b_candidates,
         "weights": {"kind": kind, "normalized_weights": normalized_ref, **logrefs,
                     "log_likelihood": None, "base_measure": ctx["weight_base_measure"],
                     "data_history_id": ctx["data_history_id"], "target_id": ctx["target_id"],
                     "normalization_status": normalization_status, "error": _uncertainty(p["package_id"]),
                     "coverage_status": p["unexplored_uncertainty"]["within_support"]},
         "partition": {"version": "a-string-event-to-b-integer-event-1:" + p["package_id"], "events": [
             {"id": value, "definition": ctx["event_definitions"][key], "mass": None,
              "mass_error": _uncertainty(p["package_id"]), "conditional_particles": None,
              "conditional_weights": None, "representative_candidate_ids": [c["id"] for c in candidates
                  if key in c["structure"]["event_ids"]], "is_separated_mode_proven": False}
             for key, value in event_map.items()], "partition_covers_modeled_support": False,
             "displayed_event_ids": list(event_map.values()), "unshown_modeled_mass": None,
             "unshown_modeled_representation": None, "outside_model_mass": None,
             "outside_model_mass_status": p["unexplored_uncertainty"]["outside_support"],
             "unexplored_modes": p["unexplored_uncertainty"]["unexplored_regions"]},
         "available_actions": available,
         "resource_ledger": _ref("source_cost", assets["source_cost"], "counts_and_seconds", p["package_id"]),
         "limitations": ["Native A and B schemas differ; explicit context is retained in contract.json.",
             "No new Maxwell calls; source validation records are not independent recertification.",
             "Candidate absence is not evidence of zero posterior probability.",
             "No posterior is inferred from residuals, candidate counts, or arbitrary equal weights."]}
    contract = {"schema_version": CONTRACT_VERSION, "adapter_version": ADAPTER_VERSION,
                "source_schema_version": p["version"], "source_package_id": p["package_id"],
                "b_schema_compatibility": "research-b-mhp-1.0; original schema unchanged",
                "input_context": ctx, "native_validation_records": native_validation,
                "event_id_map": event_map, "history": history, "observation_ids": obs_ids,
                "noise_realization_ids": noise_ids, "allow_bayesian_risk": allow_bayes,
                "source_weights": weights, "source_probability_reliability": reliability,
                "base_reconstruction": {"material_asset": "base_material.npy",
                    "method": p["base_reconstruction"]["method"],
                    "nuisance_parameters": p["base_reconstruction"]["nuisance_parameters"],
                    "source_full_validation_passed": True},
                "gain_scope": "Each world has its own legal value, persistent across its measurement history",
                "joint_likelihood_rule": "Use each raw observation once in the full joint-history likelihood; do not multiply an old posterior by old data again",
                "status": "REFERENCE_READY" if allow_bayes else "CANDIDATE_ONLY_READY"}
    return {"package": b, "contract": contract, "assets": assets}


def make_recommendation(bundle: dict, action_ids: list[str], *, criterion: str,
                        status: str = "RECOMMEND_MEASURE", risk: dict | None = None,
                        risk_semantics: str | None = None,
                        loss_definition: str = "Not evaluated by adapter") -> dict:
    """Create the original B recommendation format; never manufacture scores."""
    b, contract = bundle["package"], bundle["contract"]
    available = {a["id"]: a for a in b["available_actions"]}
    require(len(action_ids) == len(set(action_ids)) and set(action_ids) <= available.keys(),
            "Unknown/duplicate action, or an old action is being reused as new likelihood")
    seen_tokens = set()
    for action_id in action_ids:
        tokens = {action_id, *available[action_id]["constituent_ids"]}
        require(not tokens.intersection(seen_tokens), "Selected actions duplicate constituent observations")
        seen_tokens.update(tokens)
    if not contract["allow_bayesian_risk"]:
        require("EIG" not in criterion.upper() and "POSTERIOR" not in criterion.upper(),
                "Posterior acquisition requires reliable normalized reference weights")
        require(not risk or risk_semantics in {"minimax_finite_candidates", "declared_design_prior_not_posterior"},
                "Candidate-only risk needs separately declared design/minimax semantics")
        require(not risk or "event_eig" not in risk, "Candidate-only event EIG cannot use posterior semantics")
    require(bool(action_ids) == (status in {"RECOMMEND_MEASURE", "SIMULATION_DESIGN_ONLY"}),
            "Measurement status and action list disagree")
    interval = {"estimate": None, "lower": None, "upper": None, "confidence": None,
                "error_scope": "Not evaluated by the import adapter"}
    cost = _unknown_cost()
    if action_ids:
        selected_cost = [available[a]["cost"] for a in action_ids]
        for key in ("measurement_seconds", "switching_seconds", "energy_joules", "compute_seconds"):
            cost[key] = (sum(c[key] for c in selected_cost)
                         if all(c[key] is not None for c in selected_cost) else None)
        for key in ("maxwell_attempts", "factorizations", "rhs"):
            cost[key] = sum(c[key] for c in selected_cost)
        currencies = {c["currency"] for c in selected_cost}
        if len(currencies) == 1 and None not in currencies and all(c["financial_cost"] is not None for c in selected_cost):
            cost["currency"] = next(iter(currencies))
            cost["financial_cost"] = sum(c["financial_cost"] for c in selected_cost)
        cost["cost_status"] = "estimated" if all(c["cost_status"] != "not_available" for c in selected_cost) else "not_available"
    r = {"schema_version": "research-b-acquisition-1.0", "input_package_id": b["package_id"],
         "input_data_history_id": b["data_history_id"], "status": status, "selected_action_ids": action_ids,
         "criterion": criterion, "criterion_probability_status": "reliable_declared_reference" if contract["allow_bayesian_risk"] else (risk_semantics or "candidate_set_no_posterior"),
         "planning_horizon": len(action_ids), "loss_definition": loss_definition,
         "current_risk": deepcopy(interval), "expected_future_risk": deepcopy(interval),
         "risk_reduction": deepcopy(interval), "event_eig": deepcopy(interval),
         "posterior_event_quantiles": None, "posterior_quantile_definition": "Not evaluated",
         "resolved_ambiguities": [], "remaining_ambiguities": ["Imported finite candidate set; undiscovered worlds remain possible"],
         "calibration_risks": [contract["gain_scope"]], "missing_support_status": b["weights"]["coverage_status"],
         "cost": cost, "utility_error_components": [_uncertainty(b["package_id"])],
         "action_ranking_scope": "Only declared available actions; the adapter performs no ranking",
         "stop_reason": None if action_ids else status,
         "update": {"joint_likelihood_rule": contract["joint_likelihood_rule"],
                    "persistent_nuisance_rule": contract["gain_scope"],
                    "candidate_refresh_rule": "Return new data history to A or an independent reconstruction engine; fixed reweighting cannot improve Best-of-K",
                    "normalization_and_remainder_rule": "Retain declared support, unknown mass and probability reliability",
                    "reuse_or_new_data_ids": contract["observation_ids"] + action_ids,
                    "failure_action": "Preserve old validated candidates and unresolved ambiguity"},
         "limitations": b["limitations"]}
    for key, value in (risk or {}).items():
        require(key in {"current_risk", "expected_future_risk", "risk_reduction", "event_eig"}, "Unknown risk field")
        r[key] = deepcopy(value)
    return r


def updated_data_history(bundle: dict, new_record: dict, *, updated_history_id: str) -> dict:
    """Prepare a handoff after a real/simulated observation, without inference.

    A record must identify target, solver, units, packing, provenance and noise,
    and provide the full old+new covariance.  Distinct IDs alone do not prove
    statistical independence; the supplied covariance and provenance are kept.
    """
    b, contract = bundle["package"], bundle["contract"]
    ctx = contract["input_context"]
    require(updated_history_id != b["data_history_id"] and updated_history_id, "A new history needs a new ID")
    require(new_record["target_id"] == b["target_id"] and new_record["solver_id"] == b["model"]["solver_id"],
            "New measurement target/solver mismatch")
    require(new_record["units"] == ctx["observation_units"] and new_record["packing"] == ctx["noise_packing"],
            "New measurement units/packing mismatch")
    require(new_record["observation_id"] not in contract["observation_ids"], "Duplicate observation likelihood")
    require(new_record["noise_realization_id"] not in contract["noise_realization_ids"], "Duplicate noise realization")
    require(new_record["provenance_id"] not in {r["provenance_id"] for r in contract["history"]},
            "Same observation provenance cannot be renamed and counted twice")
    actions = {a["id"]: a for a in b["available_actions"]}
    require(new_record["action_id"] in actions, "Observation does not correspond to a declared new action")
    require(new_record["nuisance_group_by_parameter"] == ctx["nuisance"]["group_by_parameter"],
            "New observation cannot reset persistent nuisance/gain per frequency")
    require(new_record["source"] in {"simulation", "hardware"} and new_record["noise_evidence"],
            "Missing observation source or noise-dependence evidence")
    new = _real_vector(new_record["values"], ctx["noise_packing"])
    old = bundle["assets"]["observations"]
    joint = _finite(new_record["joint_covariance"], "joint covariance").astype(np.float64)
    require(joint.shape == (len(old) + len(new),) * 2, "Wrong joint-history covariance size")
    require(np.allclose(joint[:len(old), :len(old)], bundle["assets"]["covariance"]),
            "New history changed the already conditioned noise model")
    require(np.allclose(joint, joint.T), "Joint covariance must be symmetric")
    try:
        np.linalg.cholesky(joint)
    except np.linalg.LinAlgError as exc:
        raise ContractError("Joint-history covariance is not positive definite") from exc
    if actions[new_record["action_id"]]["new_independent_noise"]:
        require(new_record.get("independent_new_noise") is True and
                np.allclose(joint[:len(old), len(old):], 0), "Action requires evidenced new independent noise")
    return {"schema_version": "research-b-updated-data-history-1", "source_package_id": b["package_id"],
            "target_id": b["target_id"], "previous_data_history_id": b["data_history_id"],
            "data_history_id": updated_history_id, "old_observation_ids": contract["observation_ids"],
            "new_observation_id": new_record["observation_id"], "action_id": new_record["action_id"],
            "new_record": {k: v for k, v in new_record.items() if k not in {"values", "joint_covariance"}},
            "observations": np.concatenate((old, new)).tolist(), "joint_covariance": joint.tolist(),
            "noise_packing": "history_block_major; each block uses " + ctx["noise_packing"],
            "persistent_nuisance": ctx["nuisance"], "likelihood_rule": contract["joint_likelihood_rule"],
            "candidate_refresh_status": "PENDING_INDEPENDENT_RECONSTRUCTION",
            "fixed_candidate_best_of_k_change_claimed": False}


def export_bundle(bundle: dict, output_dir: str | Path, recommendation: dict | None = None,
                  new_history: dict | None = None) -> dict:
    """Write independent B arrays/metadata; context action URIs stay explicit.

    External geometry/action-noise assets supplied in context are not silently
    copied or declared verified.  The context producer must keep them resolvable.
    """
    root = Path(output_dir).resolve()
    require(not root.exists(), "Export must use a new independent directory")
    root.mkdir(parents=True)
    package = deepcopy(bundle["package"])
    replacements = {}
    for name, data in bundle["assets"].items():
        if isinstance(data, np.ndarray):
            path = root / (name + ".npy")
            np.save(path, data, allow_pickle=False)
        else:
            path = root / (name + ".json")
            path.write_text(json.dumps(data, indent=2, allow_nan=False), encoding="utf-8")
        replacements[f"bundle://{name}"] = path.name
    for node in _walk(package):
        if node.get("uri") in replacements:
            node["uri"] = replacements[node["uri"]]
    for name, value in (("MaterialHypothesisPackage.json", package), ("contract.json", bundle["contract"]),
                        ("AcquisitionRecommendation.json", recommendation), ("updated_data_history.json", new_history)):
        if value is not None:
            (root / name).write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")
    return {"status": bundle["contract"]["status"], "output_dir": str(root),
            "b_schema": "research-b-mhp-1.0", "new_maxwell_attempts": 0}


def waiting_for_a_real_package() -> dict:
    return {"status": "WAITING_FOR_A_REAL_PACKAGE", "real_a_candidates_imported": 0,
            "new_maxwell_attempts": 0, "fixture_is_scientific_evidence": False,
            "required": ["A native 1.0.0 handoff_ready manifest and local assets", "Explicit versioned import context"]}


def adapt_a_execution_v1(package_root: str | Path, *, expected_source: dict,
                         expected_observation: np.ndarray,
                         solver_aliases: dict | None = None,
                         future_action_covariances: dict | None = None) -> dict:
    """Import A's explicit execution profile/CommonMaterialExchange v1.

    ``expected_source`` identifies target_id, solver_id and backend_commit from
    the caller's independently read source contract.  The source target is kept,
    never equated to a differently named B posterior target.  Optional solver
    aliases are naming aliases only.  A raw online observation, supplied by the
    caller from A's observation file, must agree exactly with both manifests.

    Future noise is NOT inferred from old noise: an action is offered only when
    its covariance is explicitly supplied.  This function executes no solver.
    """
    root = Path(package_root).resolve()
    p = json.loads((root / "manifest.json").read_text())
    exchange = json.loads((root / "common_exchange_v1.json").read_text())
    profile = json.loads((root / "EXECUTION_PROFILE.json").read_text())
    require(profile["profile"] == A_EXECUTION_PROFILE and
            exchange["source_profile"] == A_EXECUTION_PROFILE and
            exchange["exchange_version"] == A_COMMON_EXCHANGE_VERSION,
            "Unsupported A execution/exchange profile")
    require(exchange["source_native_version"] == p["version"] == "1.0.0" and
            exchange["package_id"] == p["package_id"], "Native/exchange package mismatch")
    source_target = p["physics_target"]
    for key, native_key in (("target_id", "hash"), ("solver_id", "solver_id"),
                            ("backend_commit", "backend_commit")):
        require(expected_source[key] == source_target[native_key], f"Execution source mismatch: {key}")
    require(exchange["target_identity"] == source_target["hash"] and
            exchange["physics"] == source_target and exchange["noise"] == p["noise"] and
            exchange["acquisition"] == p["acquisition"], "Native/exchange physical contract mismatch")
    require(exchange["probability_semantics"]["status"] == "candidate_only" and
            exchange["probability_semantics"]["weights"] is None and
            p["weight_reliability"] == "candidate_only" and p["weights"] is None,
            "This execution profile import supports candidate-only evidence")
    require(p["noise"]["realification_order"] == "all-real-then-all-imaginary",
            "Unknown execution-profile packing; no inferred alias")
    require(p["noise"]["units"] == "dimensionless scattering amplitude",
            "Unknown execution-profile observation units")
    original = np.asarray(_load_asset(p["observation"], root))
    from_exchange = np.asarray(_load_asset(exchange["observation"], root))
    require(np.array_equal(original, from_exchange) and
            np.array_equal(original, _finite(expected_observation, "independent online observation")),
            "Original/native/exchange observations differ")
    require(not np.iscomplexobj(original) and original.ndim == 1,
            "Execution profile must already contain the declared packed real vector")
    covariance = np.asarray(_load_asset(p["noise"]["covariance"], root))
    grid = _load_asset(source_target["grid"], root)
    support = _load_asset(source_target["material_support"], root)
    nuisance = _load_asset(source_target["nuisance_support"], root)
    require(grid["kind"] == "EXACT_RADIAL_3D_SPHERICAL_MATERIAL" and
            grid["piecewise_constant_exact"] is True and grid["outer_radius"] == 1,
            "Unsupported execution-profile material grid")
    edges = _finite(grid["radial_interval_edges"], "radial grid").astype(np.float64)
    volumes = _finite(grid["interval_normalized_volumes"], "radial volumes").astype(np.float64)
    require(edges[0] == 0 and edges[-1] == 1 and np.all(np.diff(edges) > 0) and
            np.allclose(volumes, np.diff(edges ** 3), rtol=1e-12, atol=1e-14),
            "Radial material volume measure differs from the declared grid")
    require(support["same_material_all_sources"] is True and support["vacuum_epsilon"] == [1., 0.],
            "Unsupported execution-profile material law")
    require(len(exchange["material_worlds"]) == len(p["candidates"]), "Native/exchange candidate count differs")
    physical, latent, effective_geometry = [], [], []
    midpoint = (edges[:-1] + edges[1:]) / 2
    scales = np.array([.45, .18])
    rmax = scales * np.sqrt(1 - 1 / np.sqrt(3.))
    for c, world in zip(p["candidates"], exchange["material_worlds"]):
        require(c["id"] == world["id"] and world["fit"] == c["fit"] and
                world["full_maxwell_validation"] == c["validation"], "Native/exchange candidate mismatch")
        v = c["validation"]
        require(v["constitutive_relative_defect"] is None and
                v["constitutive_relative_defect_status"] == "NOT_EXPOSED_BY_ANALYTIC_MIE_INTERFACE",
                "Unexposed constitutive defect must remain explicitly unknown")
        coord = _finite(_load_asset(c["compensation"]["active_material_coordinates"], root), "physical coordinate")
        require(coord.shape == (4,) and np.all(coord[:2] >= 0) and np.all(coord[:2] <= rmax + 1e-12),
                "Candidate physical coordinates outside original radial support")
        r, h, er, ei = coord
        # Match original P2 from_geometry -> geometry, including its floating
        # point inactive-atom behavior.  A nominal h~1e-11 can encode as z=-1,
        # hence an actually absent ring.  Raw material arrays remain primary.
        z = np.where(coord[:2] > 0, -(1 - (coord[:2] / scales) ** 2) ** -2, 0.)
        z = np.where(coord[:2] == rmax, -3., z)
        active = z < -1.
        effective = scales * np.sqrt(np.maximum(1 - 1 / np.sqrt(np.where(active, -z, 1.)), 0.))
        r_effective, h_effective = effective
        field = _load_asset(c["material_field"], root)
        reconstructed = np.where((midpoint < r_effective) | ((h_effective > 0) & (np.abs(midpoint - .65) < h_effective)),
                                 1 + 0j, complex(er, ei))
        require(np.allclose(field, reconstructed, rtol=1e-12, atol=1e-12),
                "Candidate coordinates do not reconstruct the exact union-layer material")
        require(np.array_equal(field, _load_asset(world["full_material_array"], root)) and
                np.array_equal(_load_asset(c["full_predicted_measurement"], root),
                               _load_asset(world["old_measurement_prediction"], root)),
                "Native/exchange numerical arrays differ")
        physical.append(coord)
        latent.append(np.r_[z, er, ei])
        effective_geometry.append(effective)
    canonical = deepcopy(p)
    canonical["noise"]["realification_order"] = "all_real_then_all_imag"
    history_id = p["package_id"] + ":original-kR1.6-history"
    gain_group = p["package_id"] + ":persistent-real-gain"
    obs_id = p["package_id"] + ":original-observation"
    noise = {"packing": "all_real_then_all_imag",
             "covariance": _ref("covariance", covariance, "dimensionless scattering amplitude^2", p["package_id"]),
             "cross_history_covariance": None, "real_channel_variance_convention": "Covariance of the 32 packed real channels; complex total variance is twice a real-channel variance",
             "dwell_scaling": "Use explicitly provided acquisition covariance; no truth-dependent scaling",
             "reference_id": p["noise"]["id"]}
    def action(action_id, frequency, action_noise, *, old=False):
        cost = _unknown_cost()
        if not old:
            cost.update(maxwell_attempts=p["valid_k"], cost_status="estimated")
        return {"id": action_id, "kind": "frequency", "frequency_kR": frequency, "frequency_hz": None,
                "source_geometry": None, "receiver_geometry": None, "polarization_or_source_weights": None,
                "constituent_ids": [], "feasibility": "simulation_only",
                "feasibility_evidence": "Original full vector-Mie radial simulation contract, not hardware",
                "noise": action_noise, "nuisance_group_ids": [gain_group],
                "new_independent_noise": not old, "cost": cost}
    historical = action("old-kR1.6:" + p["package_id"], 1.6, noise, old=True)
    available, future_assets = [], {}
    native_actions = _load_asset(p["next_measurement_interface"]["allowed_actions"], root)
    for source_action in native_actions:
        aid = source_action["id"]
        if source_action["wavenumber"] in p["acquisition"]["frequencies"]:
            continue
        if aid not in (future_action_covariances or {}):
            continue
        cov = _finite(future_action_covariances[aid], "future action covariance")
        require(cov.shape == covariance.shape and np.allclose(cov, cov.T), "Future covariance shape/symmetry mismatch")
        try:
            np.linalg.cholesky(cov)
        except np.linalg.LinAlgError as exc:
            raise ContractError("Future covariance is not positive definite") from exc
        key = "future_action_covariance_" + str(len(available))
        future_assets[key] = cov
        action_noise = deepcopy(noise)
        action_noise["covariance"] = _ref(key, cov, "dimensionless scattering amplitude^2", p["package_id"])
        available.append(action(aid, source_action["wavenumber"], action_noise))
    event_ids = {e for c in p["candidates"] for e in c["structure"]["event_ids"]}
    ctx = {"adapter_version": ADAPTER_VERSION, "target_id": source_target["hash"],
           "expected_physics": {key: source_target[key] for key in IDENTITY_KEYS},
           "material_quantity": "relative_permittivity", "material_units": "1",
           "material_shape": [len(edges) - 1], "material_domain": "A exact package-specific concentric radial union grid, R=1",
           "material_metric": "volume-weighted componentwise complex-relative-permittivity MAE, fixed scale=1",
           "unit_evidence": grid["array_semantics"], "observation_units": p["noise"]["units"],
           "frequency_units": "kR", "noise_packing": "all_real_then_all_imag",
           "variance_convention": noise["real_channel_variance_convention"], "dwell_scaling": noise["dwell_scaling"],
           "data_history_id": history_id, "native_acquisition_id": p["acquisition"]["id"],
           "discretization_id": p["package_id"] + ":exact-union-radial-grid",
           "history": [{"observation_id": obs_id, "noise_realization_id": obs_id + ":original-noise",
                        "provenance_id": str(root / p["observation"]["path"]),
                        "real_channel_indices": list(range(len(original))), "action": historical}],
           "nuisance": {"parameter_names": ["gain"], "group_by_parameter": {"gain": gain_group},
                        "bounds": {"gain": nuisance["bounds"]},
                        "sharing_rule": "one_parameter_per_world_across_history",
                        "different_worlds_may_use_different_values": True,
                        "source_evidence": exchange["shared_calibration"]["rule"], "prior_definition": nuisance["prior"]},
           "available_actions": available, "weight_base_measure": "No posterior/proposal probability assigned to finite A candidate set",
           "event_definitions": {e: "A post-hoc radial descriptor " + e for e in event_ids}}
    bundle = adapt_a_v1(canonical, root, ctx)
    bundle["assets"].update(radial_interval_edges=edges, radial_volume_weights=volumes,
                            physical_coordinates=np.array(physical), p2_latent_coordinates=np.array(latent),
                            effective_p2_geometry=np.array(effective_geometry),
                            source_native_actions=native_actions,
                            **future_assets)
    fitted_gains = np.asarray([c["fit"]["nuisance_parameters"]["gain"] for c in p["candidates"]], dtype=np.float64)
    applied_predictions = np.stack([bundle["assets"][f"prediction_{i}"] for i in range(len(fitted_gains))])
    unscaled = applied_predictions / fitted_gains[:, None]
    require(np.allclose(unscaled * fitted_gains[:, None], applied_predictions, rtol=1e-14, atol=1e-14),
            "Source gain-applied prediction cannot be consistently unscaled")
    bundle["assets"].update(source_fitted_gains=fitted_gains, unscaled_old_fields=unscaled)
    source_solver = source_target["solver_id"]
    mapped_solver = (solver_aliases or {}).get(source_solver, source_solver)
    bundle["package"]["model"]["solver_id"] = mapped_solver
    bundle["contract"].update(source_profile=A_EXECUTION_PROFILE,
          common_exchange_version=A_COMMON_EXCHANGE_VERSION, original_observation_exactly_preserved=True,
          packing_alias={"source": p["noise"]["realification_order"], "canonical": "all_real_then_all_imag",
                         "permutation_performed": False},
          solver_naming_alias={"source": source_solver, "B": mapped_solver, "physical_equivalence_inferred_from_name": False},
          source_target_id=source_target["hash"], target_identity_renamed=False,
          radial_grid_asset="radial_interval_edges.npy", radial_volume_asset="radial_volume_weights.npy",
          physical_coordinates_asset="physical_coordinates.npy", p2_latent_coordinates_asset="p2_latent_coordinates.npy",
          effective_p2_geometry_asset="effective_p2_geometry.npy",
          coordinate_decoder_roundtrip_max_abs=float(np.max(np.abs(np.array(physical)[:, :2] - np.array(effective_geometry)))),
          coordinate_semantics="Preserve source nominal physical p, source raw field and exact original P2 from_geometry->geometry roundtrip separately; an inactive atom caused by source floating-point decode is not a discovered tiny cavity",
          source_fitted_gains_asset="source_fitted_gains.npy", unscaled_old_fields_asset="unscaled_old_fields.npy",
          prediction_semantics="B prediction_i.npy preserves A's SOURCE-GAIN-APPLIED mean g_source*F at kR1.6; it is not an unscaled field",
          unscaled_old_fields_semantics="F=prediction_i/source_fitted_gain. Use F with one newly fitted legal gain across old+new history; never multiply prediction_i by gain again",
          unscale_roundtrip_max_abs_error=float(np.max(np.abs(unscaled * fitted_gains[:, None] - applied_predictions))),
          source_global_cost_scope="Source GLOBAL ledger snapshot repeated per package; NEVER sum package snapshots",
          source_native_action_ids=[a["id"] for a in native_actions],
          source_native_actions_asset="source_native_actions.json",
          future_actions_without_explicit_noise="NOT_OFFERED; supply separately declared future-action covariance",
          constitutive_relative_defect_status="NOT_EXPOSED_BY_ANALYTIC_MIE_INTERFACE; null retained")
    return bundle


def extend_authorized_b_menu(bundle: dict, frequency_covariances: dict,
                              authorization_id: str) -> dict:
    """Clone a B-only sensing menu under this stage's explicit authorization.

    This reads the independent B overlay and frozen frequency menu.  It does
    not change A's allowed_actions, authorize a solve, reserve budget, or reuse
    historical balance.  The orchestrator still must use B's physical ledger.
    ``frequency_covariances`` maps numeric kR to explicit 32-real covariances.
    """
    stage = Path(__file__).resolve().parents[1]
    overlay = json.loads((stage / "configs/AUTHORIZATION_OVERLAY.json").read_text())
    frozen = json.loads((stage / "configs/E1_FROZEN.json").read_text())
    require(overlay["authorization_id"] == authorization_id and
            overlay["E1_E2"]["execution_authorized"] is True,
            "Missing matching B execution authorization")
    result = deepcopy(bundle)
    b, contract = result["package"], result["contract"]
    require(contract.get("source_profile") == A_EXECUTION_PROFILE and
            contract["original_observation_exactly_preserved"] is True,
            "B extension requires the verified A execution profile")
    expected = contract["input_context"]["expected_physics"]
    require(expected["backend_commit"] == "7ed4faef9bc3dd214d4de689c2203cc383f7495f" and
            expected["solver_id"] == "B1-vector-Mie-tangential-boundary-matching" and
            expected["material_law_id"] == "passive-isotropic-relative-permittivity-mu1" and
            expected["time_convention"] == "exp(-i omega t)" and
            b["model"]["dispersion_model"] == expected["dispersion_model"] and
            expected["dispersion_model"] == "shared frequency-independent complex epsilon over the declared finite band",
            "B extension cannot change solver, material law or shared dispersion target")
    source_actions = result["assets"].get("source_native_actions")
    require(source_actions is not None and any(a["wavenumber"] == 2.4 for a in source_actions),
            "Missing native 2.4 action contract from which to extend")
    prototype = next((deepcopy(a) for a in b["available_actions"] if a["frequency_kR"] == 2.4),
                     deepcopy(b["data"]["acquisition_history"][0]))
    original_ids = [a["id"] for a in b["available_actions"]]
    new_actions = []
    for key, value in frequency_covariances.items():
        frequency = float(key)
        require(np.isfinite(frequency) and frequency in frozen["frequencies"],
                "Frequency is absent from B's frozen authorized menu")
        covariance = _finite(value, "B future covariance").astype(np.float64)
        old_covariance = result["assets"]["covariance"]
        require(covariance.shape == old_covariance.shape and np.allclose(covariance, covariance.T),
                "B future covariance shape/symmetry mismatch")
        try:
            np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError as exc:
            raise ContractError("B future covariance is not positive definite") from exc
        repeat = frequency in {a["frequency_kR"] for a in b["data"]["acquisition_history"]}
        action = deepcopy(prototype)
        action.update(id=("b-repeat-kR-" if repeat else "b-kR-") + f"{frequency:g}",
                      kind="repeat" if repeat else "frequency", frequency_kR=frequency,
                      constituent_ids=[], new_independent_noise=True,
                      feasibility_evidence="B-only simulation extension under " + authorization_id + "; A native menu unchanged")
        name = "B_action_covariance_" + f"{frequency:g}"
        result["assets"][name] = covariance
        action["noise"]["covariance"] = _ref(name, covariance, "dimensionless scattering amplitude^2", b["package_id"])
        action["noise"]["cross_history_covariance"] = None
        action["noise"]["reference_id"] = f"B-{authorization_id}-fixed-anchor-real-noise-kR-{frequency:g}"
        attempts = 0 if repeat else len(b["candidates"])
        action["cost"] = _unknown_cost()
        action["cost"].update(maxwell_attempts=attempts, factorizations=28 * attempts,
                              rhs=28 * attempts, cost_status="estimated")
        _validate_action(action, contract["input_context"],
                         set(contract["input_context"]["nuisance"]["group_by_parameter"].values()))
        new_actions.append(action)
    require(len({a["id"] for a in new_actions}) == len(new_actions), "Duplicate B frequency alias")
    b["available_actions"] = new_actions
    contract.update(source_native_action_ids=[a["id"] for a in source_actions],
        B_menu_authorization_id=authorization_id, original_B_available_action_ids=original_ids,
        B_only_extended_action_ids=[a["id"] for a in new_actions],
        B_menu_extension_statement="These actions belong to B independent frozen authorization, not to A native export. No source manifest changed and no physics was executed.",
        B_menu_cost_scope="Estimated candidate forward predictions only; truth, validation, failures and refreshes must separately enter B ledger. Old-frequency repeat reuses saved fields but requires fresh measurement noise.")
    return result
