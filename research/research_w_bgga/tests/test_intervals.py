"""Small finite-model algebra tests. No physics solves or reference-file reads."""

from itertools import product

import numpy as np
import unittest

from bgga.intervals import (
    EMPIRICAL_INTERVAL,
    MATHEMATICAL_INTERVAL,
    NOT_CERTIFIED,
    FinitePriorIntervals,
    LogLikelihoodBounds,
    canonical_lift,
    combine_likelihood_bounds,
    event_posterior_bounds,
    finite_posterior_bounds,
    joint_pushforward,
    likelihood_bounds_from_potential,
    likelihood_bounds_from_residual,
    query_log_scores,
    select_query,
    task_cluster_diameters,
)


def _assert_approx(actual, expected, *, rel=1e-10, abs=1e-14):
    np.testing.assert_allclose(actual, expected, rtol=rel, atol=abs)


class TestIntervals(unittest.TestCase):
    def test_whitened_forward_triangle_envelope_and_empirical_label(self):
        prediction = np.array([[1 + 2j, 2 - 1j], [-1 + 0j, 0 + 1j]])
        observation = np.array([0.5 + 1j, 1 - 1j])
        perturbation = np.array([[0.1 + 0j, -0.2j], [0.05j, -0.1]])
        residual = np.linalg.norm(prediction - observation, axis=1)
        error = np.linalg.norm(perturbation, axis=1)
        bounds = likelihood_bounds_from_residual(residual, error)
        actual = -0.5 * np.linalg.norm(prediction + perturbation - observation, axis=1) ** 2
        assert np.all(bounds.log_lower <= actual)
        assert np.all(actual <= bounds.log_upper)
        assert np.all(bounds.labels == EMPIRICAL_INTERVAL)
        posterior = finite_posterior_bounds([0.4, 0.6], bounds)
        assert posterior.label == EMPIRICAL_INTERVAL
        assert posterior.floating_point_roundoff_certified is False


    def test_potential_envelope_caps_likelihood_and_requires_explicit_label(self):
        bounds = likelihood_bounds_from_potential([0.1, 10000], [1, 2])
        np.testing.assert_allclose(bounds.log_lower, [-1.1, -10002])
        np.testing.assert_allclose(bounds.log_upper, [0, -9998])
        assert bounds.labels[0] == EMPIRICAL_INTERVAL
        mathematical = likelihood_bounds_from_potential([1], [0.1], label=MATHEMATICAL_INTERVAL)
        assert mathematical.labels[0] == MATHEMATICAL_INTERVAL
        unknown = likelihood_bounds_from_residual([5], [np.inf])
        assert np.isneginf(unknown.log_lower[0]) and unknown.log_upper[0] == 0


    def test_unqueried_atoms_keep_support_and_omitted_mass_formula(self):
        bounds = LogLikelihoodBounds.unqueried(3)
        assert not bounds.queried.any()
        initial = finite_posterior_bounds([0.5, 0.3, 0.2], bounds)
        np.testing.assert_allclose(initial.lower, 0)
        np.testing.assert_allclose(initial.upper, 1)
        assert initial.log_evidence_lower == -np.inf
        _assert_approx(initial.log_evidence_upper, 0)
        assert initial.joint_weights is None
        bounds = bounds.with_updates([0], LogLikelihoodBounds.exact(np.log([0.2])))
        omitted = event_posterior_bounds([0.5, 0.3, 0.2], bounds, np.array([False, True, True]))
        assert omitted.lower == 0
        _assert_approx(omitted.upper, 0.5 / (0.1 + 0.5))
        np.testing.assert_array_equal(bounds.queried, [True, False, False])


    def test_event_fraction_bounds_are_sharp_over_small_interval_box(self):
        prior = np.array([0.1, 0.2, 0.3, 0.4])
        lower = np.array([0.01, 0.1, 0.2, 0.05])
        upper = np.array([0.5, 0.9, 0.8, 0.6])
        event = np.array([True, False, True, False])
        bounds = LogLikelihoodBounds(np.log(lower), np.log(upper), MATHEMATICAL_INTERVAL)
        result = event_posterior_bounds(prior, bounds, event)
        probabilities = []
        for corners in product((False, True), repeat=4):
            likelihood = np.where(corners, upper, lower)
            weights = prior * likelihood
            probabilities.append(weights[event].sum() / weights.sum())
        _assert_approx(result.lower, min(probabilities))
        _assert_approx(result.upper, max(probabilities))
        assert result.label == MATHEMATICAL_INTERVAL


    def test_log_normalization_retains_very_small_likelihoods_and_joint_tails(self):
        bounds = LogLikelihoodBounds.exact([-10000, -10001, -11000])
        posterior = finite_posterior_bounds([0.25, 0.75, 1e-200], bounds)
        expected_first = 0.25 / (0.25 + 0.75 * np.exp(-1))
        _assert_approx(posterior.joint_weights[0], expected_first)
        assert np.isfinite(posterior.log_joint_weights[2])
        assert posterior.log_joint_weights[2] < -1000
        _assert_approx(posterior.joint_weights.sum(), 1)
        np.testing.assert_allclose(posterior.lower, posterior.joint_weights)
        np.testing.assert_allclose(posterior.upper, posterior.joint_weights)


    def test_complement_sum_avoids_cancellation_for_tiny_event(self):
        bounds = LogLikelihoodBounds.exact([0, -700])
        posterior = finite_posterior_bounds([0.5, 0.5], bounds)
        event = posterior.event(np.array([False, True]))
        assert event.lower > 0
        _assert_approx(event.lower, np.exp(-700), rel=1e-12, abs=0)
        assert event.upper == event.lower


    def test_huge_common_log_likelihood_preserves_prior_odds(self):
        posterior = finite_posterior_bounds([0.2, 0.8], LogLikelihoodBounds.exact([-1e100, -1e100]))
        np.testing.assert_allclose(posterior.joint_weights, [0.2, 0.8])
        assert posterior.log_mass_offset == -1e100
        _assert_approx(posterior.event(np.array([True, False])).lower, 0.2)
        same_width = LogLikelihoodBounds([-np.inf, -np.inf], [-1e100, -1e100])
        assert select_query([0.2, 0.8], same_width) == 1


    def test_zero_prior_and_structural_events_are_handled_without_zero_division(self):
        posterior = finite_posterior_bounds([0, 1], LogLikelihoodBounds.unqueried(2))
        np.testing.assert_allclose(posterior.lower, [0, 1])
        np.testing.assert_allclose(posterior.upper, [0, 1])
        assert posterior.event(np.array([True, False])).upper == 0
        assert posterior.event(np.array([False, True])).lower == 1
        assert posterior.event(np.array([True, True])).lower == 1
        assert posterior.event(np.array([False, False])).upper == 0
        with self.assertRaisesRegex(ValueError, "evidence upper"):
            finite_posterior_bounds([0.5, 0.5], LogLikelihoodBounds.exact([-np.inf, -np.inf]))


    def test_query_selection_uses_prior_width_and_declared_cost(self):
        bounds = LogLikelihoodBounds([np.log(0.99), -np.inf], [0, 0])
        assert select_query([0.9, 0.1], bounds) == 1
        assert select_query([0.9, 0.1], bounds, cost=[1, 100]) == 0
        assert select_query([0.9, 0.1], bounds, eligible=np.array([True, False])) == 0
        assert select_query([0.5, 0.5], LogLikelihoodBounds.unqueried(2)) == 0


    def test_query_selection_does_not_drop_underflowing_mass_or_requery_exact_atoms(self):
        bounds = LogLikelihoodBounds([-np.inf, -11000], [0, -10000])
        scores = query_log_scores([0, 1e-320], bounds)
        assert np.isneginf(scores[0]) and np.isfinite(scores[1])
        assert select_query([0, 1e-320], bounds) == 1
        assert select_query([0.5, 0.5], LogLikelihoodBounds.exact([-3, -4])) is None


    def test_factor_refresh_replaces_same_observation_and_matches_batch(self):
        prior = np.array([0.4, 0.6])
        state = FinitePriorIntervals(prior)
        old = LogLikelihoodBounds.exact([-1, -4])
        new = LogLikelihoodBounds.exact([-3, -0.5])
        independent = LogLikelihoodBounds.exact([-2, -1])
        state.set_factor("f1", old, observation_id="y1")
        state.set_factor("f2", independent, observation_id="y2")
        state.set_factor("f1", new, observation_id="y1", version=2)
        direct = finite_posterior_bounds(prior, combine_likelihood_bounds([new, independent]))
        np.testing.assert_allclose(state.posterior().log_joint_weights, direct.log_joint_weights)
        assert len(state.factors()) == 2
        assert state.factors()[0].version == 2


    def test_factor_refresh_can_recover_support_and_validates_ownership(self):
        state = FinitePriorIntervals([0.5, 0.5])
        state.set_factor("f", LogLikelihoodBounds.exact([-np.inf, 0]), observation_id="same-data")
        state.set_factor("f", LogLikelihoodBounds.exact([0, -np.inf]), version=2)
        np.testing.assert_allclose(state.posterior().joint_weights, [1, 0])
        with self.assertRaisesRegex(ValueError, "duplicate likelihood"):
            state.set_factor("duplicate", observation_id="same-data")
        with self.assertRaisesRegex(ValueError, "increase version"):
            state.set_factor("f", version=2)
        with self.assertRaisesRegex(ValueError, "preserve observation"):
            state.set_factor("f", observation_id="different-data", version=3)


    def test_state_retains_prior_log_tails_that_normalized_floats_cannot_represent(self):
        state = FinitePriorIntervals([1e300, 1e-300])
        assert state.prior_weight[1] == 0  # Display weights can underflow; source logs must survive.
        state.set_factor("observation")
        assert state.select_query(eligible=np.array([False, True])) == 1
        state.set_factor("observation", LogLikelihoodBounds.exact([-1500, 0]), version=2)
        posterior = state.posterior()
        assert posterior.joint_weights[1] > 0.999
        assert np.all(np.isfinite(posterior.log_joint_weights))


    def test_refinement_keeps_other_atoms_unknown_and_labels_do_not_promote(self):
        state = FinitePriorIntervals([0.1, 0.3, 0.6])
        state.set_factor("observation")
        state.refine_factor("observation", [1], likelihood_bounds_from_potential([2], [0.5]))
        combined = state.combined_bounds()
        np.testing.assert_array_equal(combined.queried, [False, True, False])
        assert np.isneginf(combined.log_lower[0]) and combined.log_upper[0] == 0
        assert state.posterior().label == EMPIRICAL_INTERVAL
        uncertain = LogLikelihoodBounds([-1, -2, -3], [0, 0, 0], NOT_CERTIFIED)
        assert finite_posterior_bounds(state.prior_weight, combine_likelihood_bounds([combined, uncertain])).label == NOT_CERTIFIED
        assert finite_posterior_bounds(state.prior_weight, combine_likelihood_bounds([uncertain, combined])).label == NOT_CERTIFIED
        promoted_input = LogLikelihoodBounds([-2], [0]).with_updates([0], LogLikelihoodBounds.exact([-1]))
        assert promoted_input.labels[0] == MATHEMATICAL_INTERVAL


    def test_joint_pushforward_retains_conditional_wall_target_dependence(self):
        # Atom order is (wall0,target0), (wall0,target1), (wall1,target0), (wall1,target1).
        joint = np.array([0.45, 0.05, 0.05, 0.45])
        pushed = joint_pushforward(joint, [0, 0, 1, 1])
        np.testing.assert_allclose(pushed.cluster_mass, [0.5, 0.5])
        np.testing.assert_allclose(pushed.conditional_weights, [0.9, 0.1, 0.1, 0.9])
        np.testing.assert_allclose(pushed.lift(), joint)
        _assert_approx(np.abs(pushed.lift() - np.full(4, 0.25)).sum() / 2, 0.4)
        with self.assertRaises(ValueError):
            pushed.joint_weights[0] = 0


    def test_canonical_lift_uses_prior_conditional_and_does_not_claim_joint_fidelity(self):
        lifted = canonical_lift([0.75, 0.25], [0.45, 0.05, 0.05, 0.45], [0, 0, 1, 1])
        np.testing.assert_allclose(lifted.joint_weights, [0.675, 0.075, 0.025, 0.225])
        assert lifted.joint_error_label == NOT_CERTIFIED
        with self.assertRaisesRegex(ValueError, "zero-prior cluster"):
            canonical_lift([0.5, 0.5], [1, 0], [0, 1])
        zero = canonical_lift([1, 0], [1, 0], [0, 1])
        np.testing.assert_allclose(zero.joint_weights, [1, 0])


    def test_task_cluster_diameter_checks_full_chain_and_keeps_contract(self):
        geometry = task_cluster_diameters([0, 0.9, 1.8, 4], [0, 0, 0, 1],
                                          task_contract_id="declared-task-space-v1", norm="euclidean")
        np.testing.assert_allclose(geometry.diameter, [1.8, 0])
        np.testing.assert_array_equal(geometry.member_count, [3, 1])
        assert geometry.task_contract_id == "declared-task-space-v1"
        assert geometry.diameter[0] > 1.0  # Adjacent distance <= 1 does not bound the cluster.


    def test_task_diameter_max_probe_l2_differs_from_concatenated_euclidean(self):
        # Both probes have the same unit response difference between the two atoms.
        white_responses = np.array([[[0.0, 0.0], [0.0, 0.0]],
                                    [[1.0, 0.0], [1.0, 0.0]]])
        maximum = task_cluster_diameters(white_responses, [0, 0],
                                         task_contract_id="max-probe-contract", norm="max_probe_l2")
        concatenated = task_cluster_diameters(white_responses, [0, 0],
                                              task_contract_id="concat-coordinate-contract", norm="euclidean")
        np.testing.assert_allclose(maximum.diameter, [1.0])
        np.testing.assert_allclose(concatenated.diameter, [np.sqrt(2.0)])
        self.assertEqual(maximum.norm, "max_probe_l2")
        self.assertEqual(maximum.norm_definition, "max_probe ||delta[probe, :]||_2")
        self.assertEqual((maximum.probe_axis, maximum.measurement_axis), (1, 2))
        self.assertEqual(maximum.input_shape, (2, 2, 2))
        self.assertIn("prewhitened", maximum.coordinate_assumption)
        self.assertEqual(concatenated.norm, "euclidean")
        self.assertIn("weightedL2", concatenated.metric)
        self.assertIsNone(concatenated.probe_axis)


    def test_max_probe_l2_requires_explicit_probe_axis_and_real_measurements(self):
        for without_probe_axis in ([0.0, 1.0], [[0.0, 0.0], [1.0, 1.0]]):
            with self.subTest(shape=np.asarray(without_probe_axis).shape):
                with self.assertRaisesRegex(ValueError, "probe axis"):
                    task_cluster_diameters(without_probe_axis, [0, 0],
                                           task_contract_id="max-probe-contract", norm="max_probe_l2")
        with self.assertRaisesRegex(ValueError, "real measurement"):
            task_cluster_diameters(np.zeros((2, 2, 2), dtype=complex), [0, 0],
                                   task_contract_id="max-probe-contract", norm="max_probe_l2")


    def test_invalid_log_likelihood_bounds_are_rejected(self):
        cases = [([1], [0]), ([float('nan')], [0]), ([-2], [0.1]), ([-3], [-float('inf')])]
        for lower, upper in cases:
            with self.subTest(lower=lower, upper=upper):
                with self.assertRaises(ValueError):
                    LogLikelihoodBounds(lower, upper)


if __name__ == '__main__':
    unittest.main()
