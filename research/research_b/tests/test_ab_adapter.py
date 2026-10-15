"""Scientific handoff contract tests.  All input materials here are fixtures.

These tests do not validate a real Research A candidate or run any Maxwell code.
"""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("ab_adapter", ROOT / "code" / "ab_adapter.py")
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)


def schema_check(value, schema):
    """Independent check of every primitive used in the unchanged B schemas."""
    if "anyOf" in schema:
        failures = []
        for branch in schema["anyOf"]:
            try:
                schema_check(value, branch)
                return
            except AssertionError as error:
                failures.append(error)
        raise AssertionError(f"No anyOf branch fits: {failures}")
    if "const" in schema:
        assert value == schema["const"]
    if "enum" in schema:
        assert value in schema["enum"]
    kind = schema.get("type")
    if kind == "object":
        assert isinstance(value, dict)
        assert set(schema.get("required", [])) <= value.keys()
        if schema.get("additionalProperties") is False:
            assert value.keys() <= schema["properties"].keys()
        for key, child in value.items():
            if key in schema.get("properties", {}):
                schema_check(child, schema["properties"][key])
    elif kind == "array":
        assert isinstance(value, list)
        for child in value:
            schema_check(child, schema["items"])
    elif kind == "string":
        assert isinstance(value, str) and len(value) >= schema.get("minLength", 0)
    elif kind == "boolean":
        assert isinstance(value, bool)
    elif kind == "null":
        assert value is None
    elif kind in {"number", "integer"}:
        assert isinstance(value, (int, float)) and not isinstance(value, bool)
        if kind == "integer":
            assert isinstance(value, int)
        assert np.isfinite(value)
        assert value >= schema.get("minimum", -np.inf)
        assert value <= schema.get("maximum", np.inf)
    else:
        assert kind is None, f"Unimplemented schema assertion {schema}"


def fixture(root):
    def asset(name, data):
        np.save(root / name, data)
        # This is inert fixture metadata, never a computed/verified digest.
        return {"path": name, "sha256": "0" * 64, "format": "npy", "key": None,
                "shape": list(data.shape), "dtype": str(data.dtype)}
    y = np.array([1 + 2j, 3 + 4j])
    observation = asset("observation.npy", y)
    cov = asset("covariance.npy", np.eye(4))
    generic = asset("generic.npy", np.ones(1))
    def bref(name, shape, units="1"):
        return {"uri": name, "format": "npy", "shape": shape, "dtype": "float64",
                "units": units, "provenance_id": "fixture_only"}
    noise = {"packing": "all_real_then_all_imag", "covariance": bref("cov.npy", [4, 4]),
             "cross_history_covariance": None, "real_channel_variance_convention": "Real channel variance = 1",
             "dwell_scaling": "fixed", "reference_id": "noise1"}
    def action(name, kind, frequency):
        result = {"id": name, "kind": kind, "frequency_kR": frequency, "frequency_hz": None,
                  "source_geometry": None, "receiver_geometry": None,
                  "polarization_or_source_weights": None, "constituent_ids": [],
                  "feasibility": "simulation_only", "feasibility_evidence": "Interface fixture only",
                  "noise": deepcopy(noise), "nuisance_group_ids": ["episode_gain"],
                  "new_independent_noise": True, "cost": adapter._unknown_cost()}
        result["cost"].update(measurement_seconds=1., compute_seconds=0., cost_status="estimated")
        return result
    error = {"value": None, "norm": "relative_l2", "status": "unknown", "method": "fixture_only"}
    cost = {"new_full_state_frequency_attempts": 0, "vector_rhs": 0, "forward_rhs": 0,
            "adjoint_rhs": 0, "tangent_rhs": 0, "cache_hits": 0, "proposal_count": 2,
            "rejected_count": 0, "failed_count": 0, "wall_seconds": 0., "peak_memory_bytes": None,
            "setup_allocation_seconds": 0., "historical_cost_ledger": None}
    target = {"hash": "fixture_target_id", "material_law_id": "relative_epsilon",
              "time_convention": "exp(-i omega t)", "grid": generic,
              "backend_commit": "fixture_revision", "solver_id": "fixture_solver",
              "self_term_convention": "fixture_self_term", "polarizability_map": "fixture_map",
              "boundary_conditions": "outgoing", "material_support": generic,
              "nuisance_support": generic, "prior_manifest": None, "dispersion_model": "constant_epsilon"}
    candidates = []
    for i, gain in enumerate((.9, 1.1)):
        field = np.array([1. + 0j, 2. + .1j + i, 2. + .1j, 1. + 0j])
        material = asset(f"material{i}.npy", field)
        pred = asset(f"prediction{i}.npy", y + .05 * (i + 1))
        squared = 2 * (.05 * (i + 1)) ** 2
        candidates.append({"id": f"candidate{i}", "material_field": material,
            "material_parameterization_id": "unified_cells_fixture",
            "structure": {"support_volume": 1., "cavity_count": None, "component_count": None,
                "feature_scale": .1, "descriptor_asset": generic, "descriptor_definition_id": "fixture_descriptor",
                "cluster_ids": [], "event_ids": [f"event{i}"], "topology_thresholds": [],
                "topology_persistence_rule": "not_evaluated"},
            "full_predicted_measurement": pred,
            "fit": {"whitened_residual_l2": np.sqrt(squared), "squared_residual": squared,
                "threshold_l2": 1., "noise_model_id": "noise1", "nuisance_policy": "bounded_profile",
                "nuisance_parameters": {"gain": gain}, "log_likelihood": None,
                "likelihood_semantics": "profile_not_marginal", "likelihood_numerical_error": error,
                "fit_acceptance": True},
            "validation": {"passed": True, "scope": "full_discrete_maxwell", "solver_id": "fixture_solver",
                "physics_target_hash": target["hash"], "all_illuminations_same_material": True,
                "constitutive_relative_defect": 0., "solver_relative_residual": 0.,
                "material_constraints_passed": True, "material_law_id": target["material_law_id"],
                "numerical_error": error, "independent_mesh_checked": False,
                "independent_mesh_result": None, "timestamp_utc": "fixture"},
            "compensation": {"parent_candidate_id": None, "edit_family": "fixture",
                "edit_parameters": {}, "material_compensation_parameters": {},
                "nuisance_compensation_parameters": {"gain": gain}, "active_material_coordinates": generic,
                "protected_structure_constraints": [], "path_asset": None},
            "material_distance_to_base": float(i), "measurement_distance_to_base": 0.,
            "distance_definition_id": "fixture_only", "proper_log_prior": None,
            "prior_score_is_normalized_density": False, "proposal_log_density": None,
            "posterior_weight": None, "conditional_samples_asset": None, "cost": cost,
            "rom": {"used": False, "engine": "full_solver", "anchor_id": None, "rank": None,
                "modified_current_rank": None, "refresh_count": 0, "validity": "not_applicable",
                "numerical_error": error, "uniform_certificate_claimed": False}})
    p = {"version": "1.0.0", "package_id": "fixture_not_a_real_handoff", "created_utc": "fixture",
         "status": "handoff_ready", "observation": observation,
         "noise": {"id": "noise1", "realification_order": "all_real_then_all_imag", "covariance": cov,
                   "whitener": generic, "units": "V/m", "model_discrepancy_policy": "none"},
         "acquisition": {"id": "history_native", "manifest": generic, "frequencies": [1.6],
                         "frequency_units": "kR", "n_sources": 1, "n_real_measurements": 4},
         "physics_target": target,
         "base_reconstruction": {"material_field": candidates[0]["material_field"],
             "nuisance_parameters": {"gain": .9}, "method": "fixture_only", "full_validation_passed": True, "cost": cost},
         "requested_k": 2, "valid_k": 2, "candidate_ordering_rule": "fixture_order", "candidates": candidates,
         "pairs": [], "weight_reliability": "candidate_only", "weights": None,
         "unexplored_uncertainty": {"within_support": "unknown", "outside_support": "unquantified",
             "unexplored_regions": ["unknown"], "modeled_residual_status": "unknown",
             "modeled_residual_asset": None, "known_search_restrictions": [], "coverage_guarantee": "none"},
         "next_measurement_interface": {"evaluator_entrypoint": "not_called", "environment_manifest": generic,
             "allowed_actions": generic, "response_cache_manifest": None, "new_action_calls_require_separate_budget": True,
             "posterior_information_gain_allowed": False, "candidate_only_design_semantics": "finite_set_minimax"},
         "total_cost": cost, "upstream_manifest": generic, "failed_attempts_ledger": generic}
    ctx = {"adapter_version": adapter.ADAPTER_VERSION, "target_id": "fixture_target",
           "expected_physics": {k: target[k] for k in adapter.IDENTITY_KEYS},
           "material_quantity": "relative_permittivity", "material_units": "1", "material_shape": [4],
           "material_domain": "unified_fixture_cells", "material_metric": "volume_complex_MAE",
           "unit_evidence": "Explicit synthetic fixture unit declaration", "observation_units": "V/m",
           "frequency_units": "kR", "noise_packing": "all_real_then_all_imag",
           "variance_convention": "Real channel variance = 1", "dwell_scaling": "fixed",
           "data_history_id": "history0", "native_acquisition_id": "history_native", "discretization_id": "fixture_grid",
           "history": [{"observation_id": "obs0", "noise_realization_id": "noise0",
                        "provenance_id": "fixture0", "real_channel_indices": [0, 1, 2, 3],
                        "action": action("initial", "frequency", 1.6)}],
           "nuisance": {"parameter_names": ["gain"], "group_by_parameter": {"gain": "episode_gain"},
               "bounds": {"gain": [.8, 1.2]}, "sharing_rule": "one_parameter_per_world_across_history",
               "different_worlds_may_use_different_values": True,
               "source_evidence": "Synthetic bounded shared-gain contract", "prior_definition": "uniform(.8,1.2)"},
           "available_actions": [action("repeat1", "repeat", 1.6), action("frequency24", "frequency", 2.4)],
           "weight_base_measure": "No probability measure assigned", "event_definitions": {"event0": "fixture_event0", "event1": "fixture_event1"}}
    return p, ctx


class AdapterContractTests(unittest.TestCase):
    def setUp(self):
        scratch = ROOT.parents[1] / "research" / "delegated" / "ab_adapter" / "test_runs"
        scratch.mkdir(parents=True, exist_ok=True)
        self.tmp = tempfile.TemporaryDirectory(dir=scratch)
        self.root = Path(self.tmp.name)
        self.p, self.ctx = fixture(self.root)

    def tearDown(self):
        self.tmp.cleanup()

    def adapt(self):
        return adapter.adapt_a_v1(self.p, self.root, self.ctx)

    def test_candidate_only_not_equal_weight_posterior_and_b_schema_compatible(self):
        bundle = self.adapt()
        self.assertEqual(bundle["package"]["weights"]["kind"], "none")
        self.assertIsNone(bundle["package"]["weights"]["normalized_weights"])
        self.assertFalse(bundle["contract"]["allow_bayesian_risk"])
        np.testing.assert_array_equal(bundle["assets"]["observations"], [1, 3, 2, 4])
        native_material = (self.root / "material0.npy").read_bytes()
        schema_check(bundle["package"], json.loads((ROOT / "source_package/schemas/MaterialHypothesisPackage.schema.json").read_text()))
        rec = adapter.make_recommendation(bundle, ["frequency24"], criterion="robust_composite")
        schema_check(rec, json.loads((ROOT / "source_package/schemas/AcquisitionRecommendation.schema.json").read_text()))
        adapter.export_bundle(bundle, self.root / "independent_export", rec)
        self.assertEqual(native_material, (self.root / "material0.npy").read_bytes())
        self.assertTrue((self.root / "independent_export/contract.json").exists())

    def test_reject_target_and_candidate_solver_mismatch(self):
        self.ctx["expected_physics"]["hash"] = "different_target"
        with self.assertRaisesRegex(adapter.ContractError, "target mismatch"):
            self.adapt()
        self.ctx["expected_physics"]["hash"] = self.p["physics_target"]["hash"]
        self.p["candidates"][0]["validation"]["solver_id"] = "wrong_solver"
        with self.assertRaisesRegex(adapter.ContractError, "target/solver mismatch"):
            self.adapt()

    def test_reject_units_and_wrong_packing(self):
        self.ctx["material_units"] = "S/m"
        with self.assertRaisesRegex(adapter.ContractError, "quantity/units"):
            self.adapt()
        self.ctx["material_units"] = "1"
        self.ctx["noise_packing"] = "interleaved_real_imag"
        with self.assertRaisesRegex(adapter.ContractError, "packing"):
            self.adapt()

    def test_reject_duplicate_observation_and_reused_channels(self):
        self.ctx["history"].append(deepcopy(self.ctx["history"][0]))
        with self.assertRaisesRegex(adapter.ContractError, "Duplicate observation"):
            self.adapt()
        self.ctx["history"] = self.ctx["history"][:1]
        self.ctx["history"][0]["real_channel_indices"] = [0, 1, 1, 3]
        with self.assertRaisesRegex(adapter.ContractError, "counted once"):
            self.adapt()

    def test_reject_old_action_as_extra_independent_evidence(self):
        self.ctx["available_actions"][1]["constituent_ids"] = ["obs0"]
        with self.assertRaisesRegex(adapter.ContractError, "already conditioned"):
            self.adapt()

    def test_native_frequency_cannot_disappear_from_history_then_reappear_as_new(self):
        self.p["acquisition"]["frequencies"] = [1.6, 2.4]
        with self.assertRaisesRegex(adapter.ContractError, "already conditioned native frequency"):
            self.adapt()
        self.p["acquisition"]["frequencies"] = [1.6]
        self.ctx["available_actions"][1]["frequency_kR"] = 1.6
        with self.assertRaisesRegex(adapter.ContractError, "explicit fresh-noise repeat"):
            self.adapt()

    def test_shared_gain_per_world_not_globally_equal_or_per_frequency(self):
        bundle = self.adapt()
        values = [v["nuisance_parameters"]["gain"] for v in bundle["contract"]["native_validation_records"]]
        self.assertEqual(values, [.9, 1.1])
        self.p["candidates"][0]["fit"]["nuisance_parameters"]["gain"] = [0.9, 1.1]
        with self.assertRaisesRegex(adapter.ContractError, "one legal scalar"):
            self.adapt()

    def test_reject_profile_marginal_mislabel_and_inconsistent_saved_fit(self):
        self.p["candidates"][0]["fit"]["likelihood_semantics"] = "proper_marginal"
        with self.assertRaisesRegex(adapter.ContractError, "Profile likelihood"):
            self.adapt()
        self.p["candidates"][0]["fit"]["likelihood_semantics"] = "profile_not_marginal"
        self.p["candidates"][0]["fit"]["whitened_residual_l2"] = .9
        with self.assertRaisesRegex(adapter.ContractError, "whitened fit disagree"):
            self.adapt()

    def test_reject_nonvalidated_or_nonpassive_material(self):
        self.p["candidates"][0]["validation"]["passed"] = False
        with self.assertRaisesRegex(adapter.ContractError, "full Maxwell validation"):
            self.adapt()
        self.p["candidates"][0]["validation"]["passed"] = True
        value = np.load(self.root / "material0.npy")
        value[1] = 2 - .1j
        np.save(self.root / "material0.npy", value)
        with self.assertRaisesRegex(adapter.ContractError, "passive sign"):
            self.adapt()

    def test_proposal_not_posterior_and_no_eig(self):
        self.p["weight_reliability"] = "learned_proposal"
        for i, c in enumerate(self.p["candidates"]):
            c["proposal_log_density"] = -2. - i
        bundle = self.adapt()
        self.assertEqual(bundle["package"]["weights"]["kind"], "proposal")
        self.assertIsNone(bundle["package"]["weights"]["normalized_weights"])
        with self.assertRaisesRegex(adapter.ContractError, "reliable normalized"):
            adapter.make_recommendation(bundle, ["frequency24"], criterion="event_EIG")

    def test_candidate_risk_requires_explicit_design_semantics_and_no_double_bundle(self):
        bundle = self.adapt()
        risk = {"current_risk": {"estimate": .1, "lower": None, "upper": None,
                "confidence": None, "error_scope": "Synthetic finite-candidate minimax example"}}
        with self.assertRaisesRegex(adapter.ContractError, "design/minimax"):
            adapter.make_recommendation(bundle, ["frequency24"], criterion="robust_risk", risk=risk)
        result = adapter.make_recommendation(bundle, ["frequency24"], criterion="robust_risk", risk=risk,
                     risk_semantics="minimax_finite_candidates", loss_definition="Worst finite-candidate material MAE")
        self.assertEqual(result["criterion_probability_status"], "minimax_finite_candidates")
        bundled = deepcopy(bundle["package"]["available_actions"][0])
        bundled.update(id="bundle24", kind="bundle", constituent_ids=["frequency24"])
        bundle["package"]["available_actions"].append(bundled)
        with self.assertRaisesRegex(adapter.ContractError, "duplicate constituent"):
            adapter.make_recommendation(bundle, ["frequency24", "bundle24"], criterion="robust_risk")

    def test_reference_requires_matching_independent_evidence_and_full_support(self):
        self.p["weight_reliability"] = "reference"
        self.p["weights"] = {"values": [.4, .6], "semantics": "weighted_samples",
             "normalization_domain": "declared_full_support", "represented_mass": 1., "unrepresented_mass": 0.,
             "reliability_notes": "fixture_reference_only", "correction_method": None, "effective_sample_size": None}
        for c, w in zip(self.p["candidates"], [.4, .6]):
            c["posterior_weight"] = w
        with self.assertRaisesRegex(adapter.ContractError, "Independent reference evidence"):
            self.adapt()
        self.ctx["reference_evidence"] = {"target_id": "fixture_target", "data_history_id": "history0",
             "independent_reference": True, "provenance": "fixture_only"}
        self.p["next_measurement_interface"]["posterior_information_gain_allowed"] = True
        self.assertTrue(self.adapt()["contract"]["allow_bayesian_risk"])
        self.p["weights"]["normalization_domain"] = "candidate_pool_only"
        self.assertFalse(self.adapt()["contract"]["allow_bayesian_risk"])

    def new_record(self):
        return {"target_id": "fixture_target", "solver_id": "fixture_solver", "units": "V/m",
                "packing": "all_real_then_all_imag", "observation_id": "obs1",
                "noise_realization_id": "noise1", "provenance_id": "fixture1", "action_id": "repeat1",
                "nuisance_group_by_parameter": {"gain": "episode_gain"}, "source": "simulation",
                "noise_evidence": "Independent fixture draw; no experiment inference",
                "values": [1., 2.], "joint_covariance": np.eye(6), "independent_new_noise": True}

    def test_update_counts_once_and_keeps_joint_history_without_best_k_claim(self):
        bundle, record = self.adapt(), self.new_record()
        history = adapter.updated_data_history(bundle, record, updated_history_id="history1")
        self.assertEqual(len(history["observations"]), 6)
        self.assertFalse(history["fixed_candidate_best_of_k_change_claimed"])
        self.assertEqual(history["candidate_refresh_status"], "PENDING_INDEPENDENT_RECONSTRUCTION")
        for field, duplicate in (("observation_id", "obs0"), ("provenance_id", "fixture0"),
                                 ("noise_realization_id", "noise0")):
            bad = deepcopy(record)
            bad[field] = duplicate
            with self.assertRaises(adapter.ContractError):
                adapter.updated_data_history(bundle, bad, updated_history_id="history1")

    def test_update_rejects_gain_reset_and_correlated_repeat(self):
        bundle, record = self.adapt(), self.new_record()
        record["nuisance_group_by_parameter"] = {"gain": "new_gain_at_new_frequency"}
        with self.assertRaisesRegex(adapter.ContractError, "reset persistent"):
            adapter.updated_data_history(bundle, record, updated_history_id="history1")
        record = self.new_record()
        record["joint_covariance"][0, 4] = .1
        record["joint_covariance"][4, 0] = .1
        with self.assertRaisesRegex(adapter.ContractError, "independent noise"):
            adapter.updated_data_history(bundle, record, updated_history_id="history1")

    def test_no_real_package_remains_waiting(self):
        self.assertEqual(adapter.waiting_for_a_real_package()["status"], "WAITING_FOR_A_REAL_PACKAGE")
        self.assertEqual(adapter.waiting_for_a_real_package()["real_a_candidates_imported"], 0)


if __name__ == "__main__":
    unittest.main()
