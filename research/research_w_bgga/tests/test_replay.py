"""Finite supplied-array replay contracts; no Maxwell actions or old data reads.

Written for the parent-owned metered runner. These fixtures exercise query
ownership, the global reference barrier, interval bounds and full joint scores.
They are not scientific qualification of any saved response bank.
"""

from dataclasses import replace
import importlib.util
import inspect
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

import bgga.replay as replay
from bgga.replay import (
    FinitePriorGrid, METHODS, MethodOutput, OfflineReference, REPLAY_LABEL,
    SecantQualification, complete_link_wall_clusters, m0_direct,
    m1_confidence_point, m2_angle_point, m3_task_risk_point,
    m4_multigreen_interval, m5_quotient, m6_cache_schur_control,
    offline_reference, score_interval_checkpoint, score_joint,
)


def grid(poses=((-0.2, 1.9), (0.2, 2.2)), groups=(0, 1, 2), weights=None,
         qualification=None):
    n = len(poses) * len(groups)
    group = np.tile(np.asarray(groups), len(poses))
    occupancy = np.zeros((n, 4), dtype=bool)
    occupancy[group != 0, 0] = True
    return FinitePriorGrid(
        np.ones(n) if weights is None else weights,
        np.repeat(np.asarray(poses, dtype=float), len(groups), axis=0),
        group, (group != 0).astype(int), occupancy,
        np.where(group == 0, 1.0, 1.5 + 2.25 * (group - 1)),
        np.zeros(n), qualification=qualification or SecantQualification())


class IndexOnlyOracle:
    def __init__(self, responses):
        self._responses = responses
        self.calls = []

    def __call__(self, atom_index):
        if not isinstance(atom_index, (int, np.integer)):
            raise AssertionError("oracle received more than one atom index")
        self.calls.append(int(atom_index))
        return self._responses[atom_index].copy()


def responses_for(prior, values=None):
    data = np.zeros((prior.size, 2, 16), dtype=np.complex128)
    data[:, 0, 0] = (0.002 * prior.wall_group_id + 0.003 * prior.pose[:, 0]
                      if values is None else values)
    data[:, 1, 5] = 0.001j * prior.pose[:, 1]
    return data


def manual_reference(weights):
    weights = np.asarray(weights, dtype=float)
    with np.errstate(divide="ignore"):
        logs = np.log(weights)
    return OfflineReference(weights, logs, logs, 0.0, len(weights), 0.1, 0.1, np.log(0.1))


def load_tool(filename="g1_replay.py"):
    path = Path(__file__).resolve().parents[1] / "tools" / filename
    spec = importlib.util.spec_from_file_location("g1_contract_" + path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestReplay(unittest.TestCase):
    def test_online_entry_points_have_only_prior_y_oracle_and_no_file_api(self):
        for method in METHODS:
            self.assertEqual(tuple(inspect.signature(method).parameters), ("prior", "y", "oracle"))
        source = inspect.getsource(replay)
        self.assertNotIn("np.load(", source)
        self.assertNotIn("open(", source)
        self.assertNotIn("pathlib", source)

    def test_atom_lookup_handles_pose_major_and_permuted_order(self):
        prior = grid()
        permutation = np.array([5, 0, 4, 2, 1, 3])
        permuted = FinitePriorGrid(**{name: getattr(prior, name)[permutation]
            for name in ("prior_weight", "pose", "wall_group_id", "support_id", "support", "epsinf", "ell")})
        for p, pose in enumerate(permuted.pose_nodes):
            for g, group_id in enumerate(permuted.group_ids):
                atom = permuted.atom_by_pose_group[p, g]
                np.testing.assert_array_equal(permuted.pose[atom], pose)
                self.assertEqual(permuted.wall_group_id[atom], group_id)
        bad = np.repeat(prior.pose[:1], prior.size, axis=0)
        with self.assertRaisesRegex(ValueError, "wall-group/pose"):
            replace(prior, pose=bad)

    def test_exact_controls_independently_pay_full_queries_and_whiten_complex_noise(self):
        prior = grid(poses=((0.0, 2.1),), groups=(0, 1), weights=np.array([0.25, 0.75]))
        data = np.zeros((2, 2, 16), complex)
        data[1, 0, 0] = 0.01 + 0.02j  # potential=(1^2+2^2)/2=2.5
        y = np.zeros((2, 16), complex)
        first, second = IndexOnlyOracle(data), IndexOnlyOracle(data)
        direct = m0_direct(prior, y, first)
        control = m6_cache_schur_control(prior, y, second)
        self.assertEqual(first.calls, [0, 1])
        self.assertEqual(second.calls, [0, 1])
        np.testing.assert_allclose(direct.query_log_likelihood, [0.0, -2.5])
        expected = np.array([0.25, 0.75 * np.exp(-2.5)])
        np.testing.assert_allclose(direct.joint_weights, expected / expected.sum())
        np.testing.assert_array_equal(direct.joint_weights, control.joint_weights)
        self.assertEqual(control.details["measured_Schur_speed_status"], "NOT_RUN")
        self.assertEqual(direct.label, REPLAY_LABEL)
        self.assertEqual((direct.actual_full_fe, direct.actual_rhs), (0, 0))

    def test_oracle_rejects_outer_channels_or_full_table(self):
        prior = grid()
        y = np.zeros((2, 16), complex)
        with self.assertRaisesRegex(ValueError, "2x16"):
            m0_direct(prior, y, lambda index: np.zeros((4, 16), complex))
        with self.assertRaisesRegex(ValueError, "TX0/1"):
            m0_direct(prior, np.zeros((4, 16), complex), lambda index: y)

    def test_confidence_point_preserves_selected_wall_target_conditional(self):
        prior = grid(weights=np.array([0.1, 0.4, 0.1, 0.1, 0.2, 0.1]))
        data, y = responses_for(prior), np.zeros((2, 16), complex)
        own = m0_direct(prior, y, IndexOnlyOracle(data))
        point = m1_confidence_point(prior, y, IndexOnlyOracle(data))
        selected = prior.wall_group_id == point.details["selected_wall_group_id"]
        self.assertEqual(point.logical_queries, prior.size)
        np.testing.assert_array_equal(point.joint_weights[~selected], 0)
        np.testing.assert_allclose(point.joint_weights[selected],
            own.joint_weights[selected] / own.joint_weights[selected].sum())
        self.assertEqual(point.joint_weights.shape, (prior.size,))
        self.assertFalse(point.details["safe_by_default"])

    def test_task_medoid_can_differ_from_map_wall_and_uses_legal_matched_poses(self):
        # The middle response has only .1 posterior mass, but minimizes squared
        # response distance for .45/.1/.45 at responses 0/1.5/3 noise units.
        target_mass = np.array([0.45, 0.1, 0.45])
        likelihood = np.exp(-np.array([1.125, 0.0, 1.125]))
        prior = grid(poses=((0.0, 2.1),), weights=target_mass / likelihood)
        data = np.zeros((3, 2, 16), complex)
        data[:, 0, 0] = [0.0, 0.015, 0.03]
        y = np.zeros((2, 16), complex)
        y[0, 0] = 0.015
        oracle = IndexOnlyOracle(data)
        point = m3_task_risk_point(prior, y, oracle)
        self.assertEqual(point.details["selected_wall_group_id"], 1)
        self.assertEqual(oracle.calls, [0, 1, 2])
        self.assertEqual(point.details["candidate_count"], 3)
        self.assertTrue(point.details["risk_computation_uses_own_full_queries"])

    def test_interval_stops_with_small_unknown_bound_and_does_not_delete_prior(self):
        prior = grid(poses=((0.0, 2.1),), weights=np.array([0.999, 0.0005, 0.0005]))
        data, y = np.zeros((3, 2, 16), complex), np.zeros((2, 16), complex)
        oracle = IndexOnlyOracle(data)
        output = m4_multigreen_interval(prior, y, oracle)
        self.assertEqual(oracle.calls, [0])
        self.assertEqual(output.details["posterior_kind"], "QUERIED_CONDITIONAL_POINT_ESTIMATE")
        checkpoint = output.checkpoints[-1]
        np.testing.assert_array_equal(checkpoint["queried"], [True, False, False])
        np.testing.assert_array_equal(output.joint_weights, [1.0, 0.0, 0.0])
        self.assertAlmostEqual(checkpoint["unknown_mass_upper"], 0.001)
        self.assertTrue(np.all(checkpoint["posterior_upper"][1:] > 0))
        self.assertTrue(np.all(np.isneginf(checkpoint["log_likelihood_lower"][1:])))
        self.assertTrue(np.all(checkpoint["log_likelihood_upper"][1:] == 0))
        self.assertFalse(checkpoint["floating_point_roundoff_certified"])
        np.testing.assert_array_equal(prior.prior_weight, [0.999, 0.0005, 0.0005])

    def test_interval_uses_log_evidence_and_full_queries_when_default_bound_cannot_shrink(self):
        prior = grid(poses=((0.0, 2.1),), groups=(0, 1))
        data = np.zeros((2, 2, 16), complex)
        data[:, 0, 0] = 0.1  # same phi=50 for all atoms; no float exp underflow shortcut
        y = np.zeros((2, 16), complex)
        output = m4_multigreen_interval(prior, y, IndexOnlyOracle(data))
        reference = offline_reference(prior, y, IndexOnlyOracle(data))
        self.assertEqual(output.logical_queries, 2)
        self.assertEqual(output.details["query_trace"][0]["unknown_mass_upper"], 1.0)
        self.assertAlmostEqual(reference.log_evidence, -50)
        np.testing.assert_allclose(output.joint_weights, [0.5, 0.5])
        self.assertEqual(output.checkpoints[-1]["unknown_mass_upper"], 0)
        # This diagnostic stays near one even though true omitted mass is zero.
        self.assertGreater(reference.default_unknown_envelope_floor, 0.999)

    def test_raw_positive_interval_violation_is_preserved_without_roundoff_clipping(self):
        prior = grid(poses=((0.0, 2.1),), groups=(0, 1))
        data, y = np.zeros((2, 2, 16), complex), np.zeros((2, 16), complex)
        output = m4_multigreen_interval(prior, y, IndexOnlyOracle(data))
        reference = offline_reference(prior, y, IndexOnlyOracle(data))
        checkpoint = dict(output.checkpoints[-1])
        upper = checkpoint["posterior_upper"].copy()
        upper[0] = np.nextafter(reference.joint_weights[0], -np.inf)
        checkpoint["posterior_upper"] = upper
        score = score_interval_checkpoint(prior, checkpoint, reference)
        self.assertTrue(score["any_positive_raw_bound_violation"])
        self.assertGreater(score["maximum_atom_upper_signed_violation"], 0)
        self.assertEqual(score["atom_upper_signed_violation"][0], reference.joint_weights[0] - upper[0])

    def test_complete_link_does_not_chain_merge_and_uses_stable_ids(self):
        # Stable IDs visit coordinate0, then .09, then .18. Pair0--2 fails .1.
        coordinates = np.array([0.09, 0.0, 0.18])
        future = np.abs(coordinates[:, None] - coordinates[None, :])
        clusters = complete_link_wall_clusters([20, 10, 30], np.zeros((3, 3)), future)
        self.assertEqual(clusters, ((1, 0), (2,)))
        self.assertGreater(future[1, 2], replay.M5_FUTURE_DIAMETER)

    def test_quotient_lifts_prior_conditionals_at_each_pose_and_is_scored_for_joint_error(self):
        prior = grid(groups=(0, 1), weights=np.array([0.6, 0.1, 0.1, 0.2]))
        data = np.zeros((4, 2, 16), complex)
        data[:, 0, 0] = 0.0005 * prior.wall_group_id
        y = np.zeros((2, 16), complex)
        output = m5_quotient(prior, y, IndexOnlyOracle(data))
        self.assertEqual(output.logical_queries, 4)
        self.assertEqual(output.details["cluster_count"], 1)
        np.testing.assert_allclose(output.joint_weights, prior.prior_weight)
        cells = output.details["atom_quotient_cell"]
        lifted = (output.details["quotient_cell_mass"][cells]
                  * output.details["prior_conditional_kernel"])
        np.testing.assert_allclose(lifted, output.joint_weights)
        reference = offline_reference(prior, y, IndexOnlyOracle(data))
        self.assertGreater(score_joint(prior, output, reference)["joint_TV"], 0)
        self.assertEqual(output.details["joint_error_label"], "NOT_CERTIFIED_AWAIT_INDEPENDENT_SCORING")

    def test_joint_scores_detect_wall_pose_coupling_when_pose_marginals_match(self):
        prior = grid(groups=(0, 1))
        reference = manual_reference([0.5, 0.0, 0.0, 0.5])
        weights = np.array([0.0, 0.5, 0.5, 0.0])
        with np.errstate(divide="ignore"):
            logs = np.log(weights)
        output = MethodOutput("M1", "COMPLETED", weights, logs, np.arange(4), np.zeros(4))
        score = score_joint(prior, output, reference)
        self.assertEqual(score["coordinate_CDF_sup"], {"x": 0.0, "y": 0.0})
        self.assertEqual(score["joint_TV"], 1.0)
        self.assertEqual(score["position_Bayes_regret"], 0.0)
        self.assertEqual(score["important_joint_mode_max_absolute_error"], 0.5)
        self.assertEqual(score["ground_truth_loss_status"], "NOT_SCORED_NO_GT_ACCESS")

    def test_missing_secant_qualification_and_zero_effective_rank_stay_undefined(self):
        prior = grid()
        oracle = IndexOnlyOracle(responses_for(prior))
        missing = m2_angle_point(prior, np.zeros((2, 16), complex), oracle)
        self.assertEqual((missing.status, missing.logical_queries), ("NOT_RUN", 0))
        self.assertIsNone(missing.joint_weights)
        self.assertEqual(missing.details["geometry_status"], "NOT_DEFINED")
        self.assertEqual(oracle.calls, [])

    def test_secant_retention_uses_legal_neighbors_and_empirical_threshold(self):
        poses = [(x, y) for x in (-0.6, -0.3, 0.0, 0.3, 0.6)
                 for y in (1.7, 1.9, 2.1, 2.3, 2.5)]
        walls = [(e, ell) for e in (1.5, 3.75, 6.0) for ell in (0.0, 0.125, 0.25)]
        n = len(poses) * len(walls)
        group = np.tile(np.arange(1, 10), len(poses))
        parameters = np.tile(np.array(walls), (len(poses), 1))
        occupancy = np.tile([True, False, False, False], (n, 1))
        qualification = SecantQualification(0.05, ("synthetic_fixture",), (0.05,),
            "supplied-array test only", "EMPIRICAL_CROSS_FIXTURE_THRESHOLD")
        prior = FinitePriorGrid(np.ones(n), np.repeat(poses, 9, axis=0), group,
            np.ones(n, int), occupancy, parameters[:, 0], parameters[:, 1], qualification)
        data = np.zeros((n, 2, 16), complex)
        data[:, 0, 0], data[:, 0, 1] = 0.02 * prior.epsinf, 0.2 * prior.ell
        data[:, 0, 2], data[:, 0, 3] = 0.1 * prior.pose[:, 0], 0.15 * prior.pose[:, 1]
        output = m2_angle_point(prior, data[n // 2], IndexOnlyOracle(data))
        self.assertEqual(output.logical_queries, n)
        self.assertEqual(output.status, "COMPLETED")
        self.assertEqual(output.details["geometry_label"], "EMPIRICAL_BANK_SECANT_RETENTION")
        self.assertEqual(output.details["rank_threshold"], 0.05)
        for row in output.details["group_geometry"]:
            self.assertEqual((row["rank_A"], row["rank_B"]), (2, 2))
            self.assertAlmostEqual(row["min_retention"], 1.0)
            self.assertEqual(len(row["stencils"]), 4)
            for stencil in row["stencils"]:
                self.assertTrue(0 <= stencil["minus_atom"] < n)
                self.assertTrue(0 <= stencil["plus_atom"] < n)
        no_rank = replace(prior, qualification=SecantQualification(
            1000.0, ("synthetic_fixture",), (1000.0,), "supplied-array test only",
            "EMPIRICAL_CROSS_FIXTURE_THRESHOLD"))
        undefined = m2_angle_point(no_rank, data[n // 2], IndexOnlyOracle(data))
        self.assertEqual((undefined.status, undefined.logical_queries), ("NOT_DEFINED", n))
        self.assertIsNone(undefined.joint_weights)

    def test_tool_rejects_reference_path_and_global_reference_waits_for_every_case(self):
        tool = load_tool()
        fake_audit = {"roots": {"raw_WSLAM": str(tool.OLD_ONLINE_DIR.parents[3])}, "files": [{"id": identifier,
            "path": str(tool.OLD_ONLINE_DIR / filename), "bytes": 1}
            for identifier, (filename, _) in tool.INPUT_SPECS.items()]}
        self.assertEqual(len(tool.audited_inputs(fake_audit)), 3)
        fake_audit["files"][1]["path"] = "/tmp/inference/POSTERIORS.npz"
        with self.assertRaises(ValueError):
            tool.audited_inputs(fake_audit)
        prior = grid(poses=((0.0, 2.1),))
        data = np.zeros((prior.size, 4, 16), complex)
        observations = np.zeros((2, 1, 4, 16), complex)
        actual_reference = tool.offline_reference
        calls = []
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            (output / "PLAN.json").write_text("{}\n")
            def barrier_reference(prior, y, oracle):
                barrier = json.loads((output / "ALL_ONLINE_OUTPUTS_SEALED.json").read_text())
                self.assertEqual(barrier["case_count"], 2)
                for parent in range(2):
                    directory = output / f"parent_{parent:02d}_repeat_00"
                    for method_id in tool.METHOD_DEFINITIONS:
                        self.assertTrue((directory / (method_id + ".npz")).is_file())
                        self.assertTrue((directory / (method_id + ".json")).is_file())
                calls.append(True)
                return actual_reference(prior, y, oracle)
            with patch.object(tool, "offline_reference", barrier_reference):
                receipt = tool.run_replay(output, prior, data, observations)
            self.assertEqual(len(calls), 2)
            self.assertEqual(receipt["case_method_rows"], 14)
            first = output / "parent_00_repeat_00"
            with np.load(first / "JOINT_METHOD_OUTPUTS.npz", allow_pickle=False) as archive:
                self.assertEqual(archive["joint_weights"].shape, (7, prior.size))
                self.assertTrue(np.isnan(archive["joint_weights"][2]).all())
                self.assertEqual(archive["statuses"][2], "NOT_RUN")

    def test_old_bank_relocation_is_explicit_and_preserves_exact_names_and_identity(self):
        tool = load_tool()
        original = tool.OLD_ONLINE_DIR
        audit = {"roots": {"raw_WSLAM": str(original.parents[3])}, "files": [
            {"id": identifier, "path": str(original / filename), "bytes": 1}
            for identifier, (filename, _) in tool.INPUT_SPECS.items()]}
        default = tool.audited_inputs(audit)
        self.assertEqual(default["FINITE_FULL_RESPONSE_CACHE"]["path"], (original / "FINITE_BANK.npz").resolve())
        relocated_root = Path("/tmp/explicit-old-release/online").resolve()
        relocated = tool.audited_inputs(audit, relocated_root)
        self.assertEqual(relocated["FINITE_FULL_RESPONSE_CACHE"]["path"], relocated_root / "FINITE_BANK.npz")
        self.assertEqual(relocated["FINITE_PRIOR_GRID"]["source_mode"], "EXPLICIT_OLD_BANK_RELOCATION")
        self.assertIn("prior_id", relocated["FINITE_24_DEVELOPMENT_CASES"]["keys"])
        for bad in ("/tmp/locked-bank/online", "/tmp/inference-cache/online", "/tmp/POSTERIORS/online"):
            with self.assertRaises(ValueError):
                tool.audited_inputs(audit, Path(bad))
        # Editing a CACHE_MAP entry is not an authorized relocation, even if an
        # explicit separate root was supplied for the other legitimate sources.
        audit["files"][1]["path"] = "/tmp/unrequested-new-bank/FINITE_BANK.npz"
        for root in (None, relocated_root):
            with self.assertRaises(ValueError):
                tool.audited_inputs(audit, root)
        self.assertEqual(tool.validated_prior_id(np.array("FINITE_PRIOR_ENUM4_V1"), "test"), "FINITE_PRIOR_ENUM4_V1")
        with self.assertRaises(ValueError):
            tool.validated_prior_id(np.array("DIFFERENT_PRIOR"), "test")

    def test_offline_relocation_keeps_truth_sibling_and_truth_cannot_open_before_seal(self):
        replay_tool, tool = load_tool(), load_tool("g1_offline_extra.py")
        original = replay_tool.OLD_ONLINE_DIR
        audit = {"roots": {"raw_WSLAM": str(original.parents[3])}, "files": [
            {"id": identifier, "path": str(original / filename), "bytes": 1}
            for identifier, filename in tool.OLD_INPUT_FILENAMES.items()]}
        audit["files"].append({"id": "FINITE_OFFLINE_TRUTHS_NOISE",
            "path": str(original.parent / "offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz")})
        explicit = Path("/tmp/public-old-release/finite_bank/online").resolve()
        source = tool.resolve_old_sources(audit, explicit)
        truth = source["paths"]["FINITE_OFFLINE_TRUTHS_NOISE"]
        self.assertEqual(truth, explicit.parent / "offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz")
        self.assertTrue(source["explicit_relocation"])
        for bad in ("/tmp/locked/online", "/tmp/inference/online", "/tmp/posteriors-copy/online"):
            with self.assertRaises(ValueError):
                tool.resolve_old_sources(audit, Path(bad))
        with tempfile.TemporaryDirectory() as temporary:
            journal = tool.AccessJournal(Path(temporary))
            journal.old_sources = source
            with patch.object(tool.np, "load", side_effect=AssertionError("NPZ must stay closed")):
                with self.assertRaisesRegex(ValueError, "global online seal"):
                    journal.load_npz(truth, ("truth_index",), truth=True)
            self.assertIsNone(json.loads((Path(temporary) / "SOURCE_ACCESS.json").read_text())["seal_attestation"])
        audit["files"][0]["path"] = "/tmp/new-table/FINITE_GRID.npz"
        with self.assertRaises(ValueError):
            tool.resolve_old_sources(audit)


if __name__ == "__main__":
    unittest.main()
