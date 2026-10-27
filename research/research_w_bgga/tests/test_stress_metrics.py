"""Small supplied matrices only. No PDE solves or reference/truth file reads."""

import unittest

import numpy as np

from bgga.stress_metrics import (
    PROPER_GAUSSIAN,
    RESOLVED,
    UNDEFINED_SUPPORT,
    UNRESOLVED,
    data_only_diagnostics,
    finite_pair_diagnostics,
    proper_gaussian_diagnostics,
    stress_metrics,
    svd_rank_audit,
)


class TestStressMetrics(unittest.TestCase):
    def test_orthogonal_nuisance_retains_absolute_information(self):
        a = np.array([[2.0, 0], [0, 3], [0, 0]])
        b = np.array([[0.0], [0], [5]])
        result = data_only_diagnostics(a, b, operator_error=0.01)
        self.assertEqual(result["status"], RESOLVED)
        np.testing.assert_allclose(result["wall_absolute_eigenvalues"], [4, 9])
        np.testing.assert_allclose(result["wall_data_only_schur_eigenvalues"], [4, 9])
        np.testing.assert_allclose(result["generalized_retention_eigenvalues"], [1, 1])
        np.testing.assert_allclose(result["principal_angles_radians"], [np.pi/2, np.pi/2])
        self.assertAlmostEqual(result["bare_sin_squared_sum"], 2)
        self.assertAlmostEqual(result["principal_angle_radians_sum"], np.pi)
        self.assertFalse(result["prior_assisted"])

    def test_high_retention_does_not_imply_large_absolute_information(self):
        a = np.array([[1e-4, 0], [0, 2e-4], [0, 0]])
        b = np.array([[0.0], [0], [1]])
        result = data_only_diagnostics(a, b, operator_error=1e-7)
        np.testing.assert_allclose(result["generalized_retention_eigenvalues"], [1, 1])
        np.testing.assert_allclose(result["wall_absolute_eigenvalues"], [1e-8, 4e-8], atol=0)
        self.assertNotIn("class", result)

    def test_confounded_wall_direction_has_zero_data_only_schur(self):
        a = np.array([[4.0, 0], [0, 3], [0, 0]])
        b = np.array([[7.0], [0], [0]])
        result = data_only_diagnostics(a, b, operator_error=0)
        np.testing.assert_allclose(result["wall_absolute_eigenvalues"], [9, 16])
        np.testing.assert_allclose(result["wall_data_only_schur_eigenvalues"], [0, 9])
        np.testing.assert_allclose(result["generalized_retention_eigenvalues"], [0, 1], atol=1e-14)
        np.testing.assert_allclose(result["principal_angles_radians"], [0, np.pi/2], atol=1e-14)
        self.assertAlmostEqual(result["retention_trace"], result["bare_sin_squared_sum"])

    def test_generalized_spectrum_agrees_with_whitened_schur_on_a_support(self):
        a = np.array([[2.0, 1], [1, -2], [0, 1], [1, 0]])
        b = np.array([[1.0], [2], [-1], [1]])
        result = data_only_diagnostics(a, b, operator_error=0)
        audit = result["wall_rank_audit"]
        whitener = audit["parameter_basis"] / audit["singular_values"][:audit["rank"]]
        # Independent direct full-column nuisance elimination in this test.
        schur = a.T @ a - a.T @ b @ np.linalg.solve(b.T @ b, b.T @ a)
        expected = np.linalg.eigvalsh(whitener.T @ schur @ whitener)
        np.testing.assert_allclose(result["generalized_retention_eigenvalues"], expected, atol=1e-14)
        changed = data_only_diagnostics(a @ np.array([[2.0, 1], [0, 0.5]]), b, operator_error=0)
        np.testing.assert_allclose(changed["generalized_retention_eigenvalues"], expected, atol=1e-14)

    def test_rank_threshold_is_absolute_error_without_relative_machine_rank(self):
        audit = svd_rank_audit(np.diag([1e6, 1e-10]), operator_error=1e-12)
        self.assertEqual(audit["rank"], 2)
        self.assertEqual(audit["threshold"], 1e-12)
        self.assertEqual(audit["status"], RESOLVED)
        np.testing.assert_allclose(audit["singular_value_lower"], [1e6-1e-12, 0.99e-10], atol=0)
        boundary = svd_rank_audit(np.diag([3.0, 0.5]), operator_error=0.5)
        self.assertEqual(boundary["rank"], 1)  # Strict threshold comparison.
        self.assertEqual(boundary["status"], UNRESOLVED)

    def test_rank_uncertainty_is_exposed_in_data_only_report(self):
        result = data_only_diagnostics(3*np.eye(2), np.diag([3.0, 0.75]), operator_error=0.5)
        nuisance = result["nuisance_rank_audit"]
        self.assertEqual((nuisance["rank"], nuisance["rank_lower"], nuisance["rank_upper"]), (2, 1, 2))
        self.assertEqual(result["wall_rank_audit"]["status"], RESOLVED)
        self.assertEqual(result["status"], UNRESOLVED)
        self.assertEqual(result["generalized_retention_status"], UNRESOLVED)
        np.testing.assert_allclose(result["wall_data_only_schur_eigenvalues"], [0, 0])

    def test_empty_a_support_has_undefined_retention_sums(self):
        result = data_only_diagnostics(np.zeros((3, 2)), np.zeros((3, 1)), operator_error=0.1)
        self.assertEqual(result["retained_wall_dimension"], 0)
        self.assertEqual(result["excluded_wall_dimension"], 2)
        self.assertEqual(result["generalized_retention_status"], UNDEFINED_SUPPORT)
        self.assertEqual(result["generalized_retention_operator"].shape, (0, 0))
        self.assertIsNone(result["bare_sin_squared_sum"])
        self.assertIsNone(result["principal_angle_radians_sum"])
        self.assertIsNone(result["retention_trace"])
        np.testing.assert_allclose(result["wall_absolute_eigenvalues"], [0, 0])

    def test_proper_prior_gaussian_is_separate_and_preserves_cross_covariance(self):
        a = np.array([[1.0], [2], [0]])
        b = np.array([[2.0], [0], [1]])
        k = np.array([[2.0, -1], [0, 1]])
        prior = np.array([[2.0, 0.25], [0.25, 1.5]])
        result = stress_metrics(a, b, k, prior_precision=prior, operator_error=0.01)
        gaussian = result["proper_gaussian_prior"]
        h = np.hstack((a, b)).T @ np.hstack((a, b)) + prior
        covariance = np.linalg.solve(h, np.eye(2))
        schur = h[:1, :1] - h[:1, 1:] @ np.linalg.solve(h[1:, 1:], h[1:, :1])
        np.testing.assert_allclose(gaussian["joint_covariance"], covariance)
        np.testing.assert_allclose(gaussian["wall_nuisance_cross_covariance"], covariance[:1, 1:])
        np.testing.assert_allclose(gaussian["wall_schur_information"], schur)
        np.testing.assert_allclose(gaussian["task_covariance"], k @ covariance @ k.T)
        np.testing.assert_allclose(gaussian["wall_omission_task_covariance"],
                                   gaussian["task_covariance"]-gaussian["wall_conditioned_task_covariance"], atol=1e-14)
        self.assertEqual(gaussian["label"], PROPER_GAUSSIAN)
        self.assertTrue(gaussian["prior_assisted"])
        self.assertFalse(gaussian["bare_angle_interpretation"])
        self.assertFalse(gaussian["finite_posterior_certificate"])
        self.assertNotIn("bare_sin_squared", gaussian)
        self.assertIn("bare_sin_squared", result["data_only"])

    def test_zero_nuisance_block_needs_no_pseudoinverse(self):
        a, b, k = np.diag([2.0, 3]), np.empty((2, 0)), np.eye(2)
        result = stress_metrics(a, b, k, prior_precision=np.eye(2), operator_error=0.01)
        self.assertEqual(result["data_only"]["nuisance_rank_audit"]["rank"], 0)
        gaussian = result["proper_gaussian_prior"]
        np.testing.assert_allclose(gaussian["wall_schur_eigenvalues"], [5, 10])
        np.testing.assert_allclose(gaussian["wall_conditioned_task_covariance"], np.zeros((2, 2)))
        np.testing.assert_allclose(gaussian["wall_omission_task_covariance"], gaussian["task_covariance"])
        self.assertEqual(gaussian["wall_nuisance_cross_covariance"].shape, (2, 0))

    def test_improper_or_asymmetric_prior_is_rejected_even_if_data_identify_it(self):
        a, b, k = np.array([[1.0], [0]]), np.array([[0.0], [1]]), np.eye(2)
        for prior in (np.diag([1.0, 0]), np.diag([1.0, -1]), np.array([[1.0, 0.1], [0, 1]])):
            with self.subTest(prior=prior):
                with self.assertRaises(ValueError):
                    proper_gaussian_diagnostics(a, b, k, prior_precision=prior)

    def test_finite_pair_explicit_max_probe_task_norm_and_error_envelope(self):
        result = finite_pair_diagnostics([3, 4], [[6, 8], [0, 5]], current_error=0.5,
                                         task_error=2, task_norm="max_probe_l2", task_contract_id="declared-2-probe")
        self.assertEqual(result["current_norm"], 5)
        self.assertEqual(result["task_norm_value"], 10)
        self.assertEqual((result["current_norm_lower"], result["current_norm_upper"]), (4.5, 5.5))
        self.assertEqual((result["task_norm_lower"], result["task_norm_upper"]), (8, 12))
        self.assertEqual(result["task_to_current_ratio"], 2)
        self.assertAlmostEqual(result["ratio_lower"], 8/5.5)
        self.assertAlmostEqual(result["ratio_upper"], 12/4.5)
        self.assertEqual(result["task_contract_id"], "declared-2-probe")
        self.assertFalse(result["physical_certificate"])

    def test_finite_pair_withholds_zero_or_error_unresolved_denominator(self):
        for current, error, status in (([0, 0], 0, "UNDEFINED_ZERO_DENOMINATOR"),
                                       ([0.1, 0], 0.2, "UNRESOLVED_DENOMINATOR")):
            result = finite_pair_diagnostics(current, [2], current_error=error, task_error=0.1,
                                             task_norm="euclidean", task_contract_id="scalar-task")
            self.assertEqual(result["ratio_status"], status)
            self.assertIsNone(result["task_to_current_ratio"])
            self.assertIsNone(result["ratio_lower"])
            self.assertIsNone(result["ratio_upper"])
            self.assertEqual(result["current_error"], error)
            self.assertEqual(result["task_error"], 0.1)

    def test_explicit_real_coordinates_shapes_and_errors_are_required(self):
        with self.assertRaisesRegex(ValueError, "real"):
            data_only_diagnostics(np.eye(2, dtype=complex), np.eye(2), operator_error=0)
        for error in (-1, np.nan, np.inf, True, [0.1]):
            with self.subTest(error=error):
                with self.assertRaises(ValueError):
                    svd_rank_audit(np.eye(2), operator_error=error)
        with self.assertRaisesRegex(ValueError, "current_difference"):
            finite_pair_diagnostics(1, [1], current_error=0, task_error=0,
                                    task_norm="euclidean", task_contract_id="scalar")
        with self.assertRaisesRegex(ValueError, "probes"):
            finite_pair_diagnostics([1], [1, 2], current_error=0, task_error=0,
                                    task_norm="max_probe_l2", task_contract_id="bad-shape")


if __name__ == "__main__":
    unittest.main()
