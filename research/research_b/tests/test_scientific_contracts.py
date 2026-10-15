"""Necessary numerical/authorization contracts; no physical solver is invoked.

Synthetic arrays and isolated temporary ledgers are not Maxwell validation.
The independent likelihood oracle is the unchanged Joint P2 implementation.
"""
from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

RESEARCH = Path(__file__).resolve().parents[1]
REPO = RESEARCH.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(RESEARCH / "code"))

from cost import BudgetExceeded, Cost
from sensing import conditional_log_increment, field_l1, normalized
from src.b1_joint.likelihood import conditional_gain, log_likelihood
from src.b1_joint.material import decode, from_geometry

NUMERICAL_METRICS = []


class SharedGainContracts(unittest.TestCase):
    def compare(self, old_channels, new_channels, gain, seed):
        rng = np.random.default_rng(seed)
        candidates, draws = 13, 7
        # Correlated worlds leave nontrivial within-set weights even when the
        # conditional gain is far outside its legal truncation interval.
        f0 = (rng.uniform(.5, 1.2, size=(1, old_channels))
              + rng.normal(scale=.002, size=(candidates, old_channels)))
        fa = (rng.uniform(.4, 1.1, size=(1, new_channels))
              + rng.normal(scale=.01, size=(candidates, new_channels)))
        s0 = np.linspace(.025, .04, old_channels)
        sa = np.linspace(.03, .05, new_channels)
        y0 = gain * f0[4] + rng.normal(scale=.003, size=old_channels)
        ya = gain * fa[4] + rng.normal(scale=.04, size=(draws, new_channels))
        prior = rng.uniform(.1, 1., size=candidates)
        prior /= prior.sum()
        old = log_likelihood(f0, y0, s0)
        joint_fields = np.concatenate((f0, fa), axis=1)
        joint_sd = np.r_[s0, sa]
        direct = np.array([log_likelihood(joint_fields, np.r_[y0, row],
                                          joint_sd) - old for row in ya])
        increment = conditional_log_increment(f0, y0, s0, fa, ya, sa)
        # A draw-only Gaussian normalizer is intentionally omitted by the
        # optimized update. It must cancel after posterior normalization.
        centered_error = float(abs((increment - increment[:, :1])
                                   - (direct - direct[:, :1])).max())
        actual = normalized(np.log(prior)[None, :] + old[None, :] + increment)
        expected = normalized(np.log(prior)[None, :] + old[None, :] + direct)
        posterior_error = float(abs(actual - expected).max())
        m, v = conditional_gain(f0, y0, s0)
        outside_sd = np.maximum((.7-m)/np.sqrt(v), (m-1.3)/np.sqrt(v))
        NUMERICAL_METRICS.append({
            "old_channels": old_channels, "new_channels": new_channels,
            "simulated_unconstrained_gain": gain,
            "max_gain_outside_support_SD": float(outside_sd.max()),
            "centered_log_increment_max_error": centered_error,
            "normalized_posterior_max_error": posterior_error,
        })
        self.assertTrue(np.isfinite(increment).all())
        self.assertLess(centered_error, 2e-8)
        np.testing.assert_allclose(actual, expected, rtol=2e-8, atol=2e-9)
        if gain != 1.:
            self.assertGreater(outside_sd.min(), 5.)

    def test_32_16_1_channel_updates_interior_and_both_truncation_tails(self):
        for old in (32, 16, 1):
            for new in (32, 16, 1):
                for gain in (1., 1.8, .2):
                    with self.subTest(old=old, new=new, gain=gain):
                        self.compare(old, new, gain, 271015 + 10*old + new)

    def test_zero_sensitivity_action_does_not_change_material_weights(self):
        rng = np.random.default_rng(27016)
        f0 = rng.normal(size=(9, 32))
        y0 = rng.normal(size=32)
        w = normalized(rng.normal(size=9))
        increment = conditional_log_increment(f0, y0, .2, np.zeros((9, 1)),
                                              np.array([[.1], [-4.]]), .3)
        np.testing.assert_allclose(normalized(np.log(w) + increment),
                                   np.broadcast_to(w, (2, 9)), atol=2e-14)


def independent_volume_l1(x, y):
    """Integrate the separately decoded piecewise constant materials exactly."""
    a, b = decode(x), decode(y)
    edges = sorted({0., *a["radii"], *b["radii"]})
    total = 0.
    for lo, hi in zip(edges[:-1], edges[1:]):
        midpoint = (lo + hi) / 2
        ea = a["eps_layers"][np.searchsorted(a["radii"], midpoint)]
        eb = b["eps_layers"][np.searchsorted(b["radii"], midpoint)]
        total += (hi**3-lo**3) * (abs(ea.real-eb.real) + abs(ea.imag-eb.imag))
    return total


class MaterialLossContracts(unittest.TestCase):
    def test_vectorized_l1_matches_independent_exact_shell_integral(self):
        rng = np.random.default_rng(27017)
        x = rng.uniform([-3., -3., 2.2, .04], [1., 1., 2.8, .16], (19, 4))
        x = np.vstack((x, from_geometry(0., 0.),
                       from_geometry(.25, 0.), from_geometry(0., .08),
                       from_geometry(.2, .06, 2.3, .15)))
        expected = np.array([[independent_volume_l1(a, b) for b in x] for a in x])
        actual = field_l1(x, x)
        NUMERICAL_METRICS.append({"material_L1_pairs": len(x)**2,
                                  "max_absolute_error": float(abs(actual-expected).max())})
        np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=2e-14)
        np.testing.assert_allclose(actual, actual.T, atol=2e-14)
        np.testing.assert_allclose(np.diag(actual), 0., atol=2e-14)

    def test_solid_cavity_and_material_components_have_declared_scaling(self):
        solid = from_geometry(0., 0., 2.5, .1)
        cavity = from_geometry(.25, 0., 2.5, .1)
        self.assertAlmostEqual(field_l1(solid, cavity)[0, 0], .25**3 * 1.6, 14)
        other = from_geometry(0., 0., 2.7, .15)
        self.assertAlmostEqual(field_l1(solid, other)[0, 0], .25, 14)
        inactive = solid.copy()
        inactive[:2] = [-.5, .8]
        self.assertAlmostEqual(field_l1(solid, inactive)[0, 0], 0., 14)


class CostContracts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="research-b-contract-")
        self.root = Path(self.tmp.name)
        (self.root / "configs").mkdir()
        self.overlay = json.loads((RESEARCH / "configs/AUTHORIZATION_OVERLAY.json").read_text())
        self.nice = patch("cost.os.nice", return_value=0)
        self.nice.start()
        self.write_overlay()

    def tearDown(self):
        self.nice.stop()
        self.tmp.cleanup()

    def write_overlay(self):
        (self.root / "configs/AUTHORIZATION_OVERLAY.json").write_text(
            json.dumps(self.overlay))

    @staticmethod
    def fake_kernel():
        # Accounting fixture only: not an electromagnetic or surrogate solve.
        return (None, None, None, {"dense_boundary_solves": 28})

    def test_e0_rejects_before_callback_or_reservation(self):
        cost = Cost(self.root)
        with cost.job("e0-contract", "E0"):
            with self.assertRaises(BudgetExceeded):
                cost.evaluate("FIXTURE", {}, lambda: self.fail("Callback reached in E0"))
        self.assertEqual(cost.count(), 0)
        self.assertFalse(cost.physics.exists())

    def test_failed_attempt_reserved_before_callback_and_survives_reload(self):
        self.overlay["E1_E2"]["material_frequency_attempt_cap"] = 1
        self.write_overlay()
        cost = Cost(self.root)
        def fail():
            self.assertEqual(cost.count(), 1)
            self.assertEqual(cost.rows(cost.physics)[-1]["event"], "BEGIN")
            raise RuntimeError("Injected accounting fixture failure; no solver")
        with self.assertRaisesRegex(RuntimeError, "Injected accounting"):
            with cost.job("failed-contract", "E1"):
                cost.evaluate("FAILURE_FIXTURE", {}, fail)
        recovered = Cost(self.root)
        with recovered.job("retry-contract", "E2"):
            with self.assertRaises(BudgetExceeded):
                recovered.evaluate("RETRY_FIXTURE", {}, self.fake_kernel)
        self.assertEqual(recovered.count(), 1)
        self.assertEqual(recovered.summary()["failed_attempts"], 1)
        self.assertEqual(recovered.summary()["unfinished_attempts"], 0)
        self.assertEqual(recovered.summary()["remaining_attempts"], 0)

    def test_4096th_reservation_allowed_4097th_rejected(self):
        # Higher overlay cannot enlarge the direct user authorization.
        self.overlay["E1_E2"].update(material_frequency_attempt_cap=10000,
                                     controlled_execution_wall_seconds_cap=10000)
        self.write_overlay()
        cost = Cost(self.root)
        seed = []
        for number in range(1, 4096):
            seed.extend((dict(event="BEGIN", id=number, stage="E1", job_id="synthetic-history"),
                         dict(event="END", id=number, status="PASS", wall_seconds=0.)))
        cost.physics.write_text("".join(json.dumps(row)+"\n" for row in seed))
        with cost.job("last-reservation", "E2"):
            cost.evaluate("FIXTURE", {}, self.fake_kernel)
            with self.assertRaises(BudgetExceeded):
                cost.evaluate("FIXTURE", {}, lambda: self.fail("4097th callback reached"))
        self.assertEqual(cost.count(), 4096)
        self.assertEqual(cost.summary()["remaining_attempts"], 0)
        self.assertEqual(cost.summary()["E1_E2_wall_cap"], 7200.)
        self.assertEqual(cost.summary()["boundary_RHS"], 28)

    def test_lower_overlay_attempt_cap_shared_across_e1_e2(self):
        self.overlay["E1_E2"]["material_frequency_attempt_cap"] = 2
        self.write_overlay()
        cost = Cost(self.root)
        for stage in ("E1", "E2"):
            with cost.job("lower-cap-"+stage, stage):
                cost.evaluate("FIXTURE", {}, self.fake_kernel)
        with cost.job("no-balance-transfer", "E2"):
            with self.assertRaises(BudgetExceeded):
                cost.evaluate("FIXTURE", {}, self.fake_kernel)
        self.assertEqual(cost.summary()["cap"], 2)
        self.assertEqual(cost.count(), 2)

    def test_disabled_authorization_rejected_before_reservation(self):
        self.overlay["E1_E2"]["execution_authorized"] = False
        self.write_overlay()
        cost = Cost(self.root)
        with self.assertRaisesRegex(BudgetExceeded, "not authorized"):
            with cost.job("disabled-overlay", "E1"):
                self.fail("Unauthorized E1 body reached")
        self.assertEqual(cost.count(), 0)

    def test_unregistered_stage_and_unmetered_call_rejected(self):
        cost = Cost(self.root)
        with cost.job("unknown-stage", "E3"):
            with self.assertRaises(BudgetExceeded):
                cost.evaluate("FIXTURE", {}, self.fake_kernel)
        cost.stage = "E1"
        with self.assertRaises(BudgetExceeded):
            cost.evaluate("FIXTURE", {}, self.fake_kernel)
        self.assertEqual(cost.count(), 0)

    def test_lower_overlay_wall_cap_rejects_before_job_body(self):
        self.overlay["E1_E2"]["controlled_execution_wall_seconds_cap"] = .5
        self.write_overlay()
        cost = Cost(self.root)
        # Prior failed jobs and E1/E2 both consume the same wall allocation.
        cost.jobs.write_text(json.dumps(dict(id="prior-failure", stage="E1",
            status="FAIL", wall_seconds=.5, cpu_seconds=.1))+"\n")
        with self.assertRaisesRegex(BudgetExceeded, "wall cap 0.5s"):
            with cost.job("wall-exhausted", "E2"):
                self.fail("Exhausted wall cap body reached")
        self.assertEqual(cost.count(), 0)
        self.assertEqual(cost.summary()["E1_E2_wall_cap"], .5)

    def test_invalid_overlay_caps_do_not_enable_execution(self):
        original = deepcopy(self.overlay)
        for key, value in (("material_frequency_attempt_cap", -1),
                           ("material_frequency_attempt_cap", 2.5),
                           ("controlled_execution_wall_seconds_cap", float("nan"))):
            self.overlay = deepcopy(original)
            self.overlay["E1_E2"][key] = value
            self.write_overlay()
            with self.subTest(key=key, value=value):
                with self.assertRaises(ValueError):
                    Cost(self.root)


if __name__ == "__main__":
    unittest.main()
