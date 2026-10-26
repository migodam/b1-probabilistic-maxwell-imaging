"""Synthetic paired tests only; never read physics/truth/model artifacts."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np

spec = importlib.util.spec_from_file_location("statistics_confirm", Path(__file__).with_name("statistics_confirm.py"))
stats = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stats)


class FrozenPairedStatisticsTests(unittest.TestCase):
    def dataset(self, n=512):
        common = np.linspace(1, 4, n)
        base = np.stack([common + x for x in (3, 2, 1, 0.5, 0.8)], axis=1)
        losses = np.repeat(base[:, None, :], 3, axis=1)
        return {
            "parent_ids": np.array([f"p{i}" for i in range(n)]),
            "batch_ids": np.array(["A"] * (n // 2) + ["B"] * (n - n // 2)),
            "seed_ids": np.array([1, 2, 3]), "methods": np.array(stats.METHODS),
            "raw_losses": losses, "clipped_losses": losses * 0.9,
            "split": np.array(["TEST"] * n), "phase": np.array("independent_confirmation"),
        }

    def test_pairing_cancels_shared_parent_variation(self):
        data = self.dataset(12)
        mean_losses = data["raw_losses"].mean(axis=1)
        draws = stats.bootstrap_means(mean_losses, data["batch_ids"])
        effect = stats.effects(draws)
        np.testing.assert_allclose(effect[:, 0], 0.5, atol=1e-14)
        np.testing.assert_allclose(effect[:, 1], 0.2, atol=1e-14)
        np.testing.assert_allclose(effect[:, 3], 0.4, atol=1e-14)
        self.assertGreater(np.std(draws[:, 0]), 0.01)

    def test_stratification_keeps_fixed_batch_weights(self):
        values = np.array([[0.0], [10.0], [10.0], [10.0]])
        draws = stats.bootstrap_means(values, np.array(["A", "B", "B", "B"]))
        np.testing.assert_array_equal(draws, np.full((20000, 1), 7.5))
        with self.assertRaises(ValueError):
            stats.bootstrap_means(values, np.array(["A", "B", "B", "B"]), replicates=19999)

    def test_retention_positive_and_unstable_denominators(self):
        values = self.dataset(12)["raw_losses"].mean(axis=1)
        bootstrap = stats.bootstrap_means(values, np.array(["A"] * 12))
        result = stats.metric_summary(values, bootstrap)
        self.assertAlmostEqual(result["retention_ratio"]["estimate"], 1.2)
        self.assertTrue(result["comparisons"]["retention_80pct"]["confidence_gate_with_stable_denominator"])
        # Zero K32 gain does not permit a retention ratio or retention claim.
        values[:, 2] = values[:, 1]
        bootstrap = stats.bootstrap_means(values, np.array(["A"] * 12))
        result = stats.metric_summary(values, bootstrap)
        self.assertIsNone(result["retention_ratio"]["estimate"])
        self.assertFalse(result["comparisons"]["retention_80pct"]["confidence_gate_with_stable_denominator"])
        self.assertEqual(result["retention_ratio"]["zero_denominator_bootstrap_draws"], 20000)

    def test_loader_split_reorder_seedmean_and_ensemble_are_separate(self):
        data = self.dataset(6)
        data["phase"] = np.array("development")
        data["split"] = np.array(["VAL", "VAL", "VAL", "TRAIN", "TRAIN", "TRAIN"])
        data["raw_losses"][:, :, 1:] += np.array([-0.3, 0.0, 0.3])[None, :, None]
        data["ensemble_raw_losses"] = np.full((6, 5), 0.123)
        order = np.array([4, 3, 2, 1, 0])
        data["methods"] = data["methods"][order]
        for key in ("raw_losses", "clipped_losses"):
            data[key] = data[key][:, :, order]
        data["ensemble_raw_losses"] = data["ensemble_raw_losses"][:, order]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "synthetic.npz"
            np.savez(path, **data)
            with self.assertRaises(ValueError):
                stats.load_losses(path)
            loaded = stats.load_losses(path, split="VAL", phase="development")
        self.assertEqual(len(loaded["parent_ids"]), 3)
        result = stats.analyze(loaded, split="VAL")
        self.assertEqual(result["status"], "DESCRIPTIVE_OUTSIDE_PRIMARY_PHASE_SPLIT")
        self.assertAlmostEqual(result["ensemble_descriptive"]["raw"]["strong_cnn"], 0.123)
        self.assertNotAlmostEqual(result["raw_primary"]["method_means"]["strong_cnn"], 0.123)
        self.assertFalse(result["raw_primary"]["comparisons"]["distilled_vs_strong"]["confirmatory_eligible"])

    def test_confirmatory_shape_and_fail_closed_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "synthetic.npz"
            data = self.dataset(8)
            np.savez(path, **data)
            loaded = stats.load_losses(path)
            with self.assertRaises(ValueError):
                stats.analyze(loaded)
            data["parent_ids"][1] = data["parent_ids"][0]
            np.savez(path, **data)
            with self.assertRaises(ValueError):
                stats.load_losses(path)
            data = self.dataset(8)
            data["raw_losses"][0, 1, 0] += 1
            np.savez(path, **data)
            with self.assertRaises(ValueError):
                stats.load_losses(path)
            data = self.dataset(8)
            data["clipped_losses"][0, 0, 1] = np.nan
            np.savez(path, **data)
            with self.assertRaises(ValueError):
                stats.load_losses(path)

    def test_bonferroni_quantiles_and_centered_pvalues(self):
        draws = np.linspace(-1, 1, 20000)
        bounds = stats.confidence(draws)
        np.testing.assert_allclose(bounds["ci_two_sided_98_75"], [-0.9875, 0.9875])
        self.assertAlmostEqual(bounds["lower_one_sided_98_75"], -0.975)
        p = stats.centered_pvalue(2.0, np.full(20000, 2.0), one_sided=True)
        self.assertAlmostEqual(p["p_bootstrap_centered"], 1 / 20001)
        self.assertAlmostEqual(p["p_bonferroni_4"], 4 / 20001)
        self.assertEqual(stats.centered_pvalue(0.0, draws)["p_bootstrap_centered"], 1.0)

    def test_full_confirmation_has_only_pooled_raw_eligibility_and_json(self):
        data = self.dataset()
        # A and B point in opposite directions for deterministic vs K32.
        data["raw_losses"][:256, :, 3] = data["raw_losses"][:256, :, 2] - 0.6
        data["raw_losses"][256:, :, 3] = data["raw_losses"][256:, :, 2] + 0.2
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "synthetic.npz"
            np.savez(path, **data)
            result = stats.analyze(stats.load_losses(path))
        self.assertEqual(result["status"], "CONFIRMATORY_ANALYSIS")
        self.assertEqual(result["batch_counts"], {"A": 256, "B": 256})
        self.assertTrue(result["batch_consistency"]["deterministic_vs_k32"]["opposite_signs"])
        self.assertTrue(result["raw_primary"]["comparisons"]["deterministic_vs_k32"]["confirmatory_eligible"])
        self.assertFalse(result["clipped_secondary"]["comparisons"]["deterministic_vs_k32"]["confirmatory_eligible"])
        self.assertFalse(result["batches_descriptive"]["A"]["raw"]["comparisons"]["deterministic_vs_k32"]["confirmatory_eligible"])
        self.assertTrue(json.dumps(result, allow_nan=False))

    def test_ratio_never_drops_zero_denominator_draws(self):
        values = self.dataset(12)["raw_losses"].mean(axis=1)
        draws = np.repeat(values.mean(axis=0)[None, :], 20000, axis=0)
        # One zero-denominator tail draw: stable CI still holds, but no ratio.
        draws[0, 2] = draws[0, 1]
        result = stats.metric_summary(values, draws)
        self.assertTrue(result["retention_ratio"]["positive_stable_denominator"])
        self.assertIsNone(result["retention_ratio"]["bounds"])
        self.assertEqual(result["retention_ratio"]["status"], "NOT_ESTIMABLE_EXACT_ZERO_BOOTSTRAP_DENOMINATOR")


if __name__ == "__main__":
    unittest.main()
