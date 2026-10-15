"""Necessary structure/packing/bundle checks. No solver or paid data is read.

Synthetic fixtures exist only inside TemporaryDirectory and are explicitly
labelled NOT_A_PHYSICAL_CANDIDATE. They are never a G discovery handoff.
"""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np

from research_t.structure import (absolute_masks, cavity_birth, check_protected,
                                  material_distances, protected_bounds, topology)
from research_t.handoff import (FrozenPhysicsInterface, new_measurement, pack_complex_per_source,
                               validate_bundle, whitened_residual, write_bundle)


class StructureSemantics(unittest.TestCase):
    def test_absolute_air_is_not_background_increment(self):
        background = np.full((24, 24), 1.2 + .03j)
        self.assertEqual(int((abs(background - background) < .1).sum()), 576)
        self.assertEqual(int(absolute_masks(background)["air"].sum()), 0)
        self.assertEqual(topology(background)["beta_1_material"], 0)

    def test_connected_material_and_enclosed_air(self):
        chi = np.zeros((24, 24), complex)
        chi[3:21, 3:21] = 1.2 + .03j
        chi[8:16, 8:16] = 0
        result = topology(chi)
        self.assertEqual((result["beta_0_material"], result["beta_1_material"]), (1, 1))
        # A one-cell diagonal contact is separate under 4-connected material.
        chi[:] = 0
        chi[3, 3] = chi[4, 4] = 1
        self.assertEqual(topology(chi)["beta_0_material"], 2)

    def test_protected_core_and_guard_cannot_be_filled_or_vented(self):
        base = np.full((24, 24), 1.2 + .03j)
        edit = cavity_birth(base, center_row_col=(12, 12), radius_cells=2,
                            min_feature_radius_cells=2, grid_spacing=1/24)
        self.assertTrue(check_protected(edit.chi, edit.core, edit.guard)["preserved"])
        lower, upper = protected_bounds(edit.core, edit.guard)
        core = np.flatnonzero(edit.core)
        np.testing.assert_array_equal(lower[core], 0)
        np.testing.assert_array_equal(upper[576 + core], 0)
        self.assertTrue(np.all(lower[np.flatnonzero(edit.guard)] >= .3))
        damaged = edit.chi.copy()
        damaged[12, 12] = .5
        self.assertFalse(check_protected(damaged, edit.core, edit.guard)["preserved"])
        vented = edit.chi.copy()
        vented[12, 12:] = 0
        self.assertFalse(check_protected(vented, edit.core, edit.guard)["preserved"])
        with self.assertRaises(ValueError):
            cavity_birth(base, center_row_col=(0, 0), radius_cells=2,
                         min_feature_radius_cells=2, grid_spacing=1/24)

    def test_distance_fixed_units_and_six_cell_novelty(self):
        a = np.full((24, 24), 1.2, complex)
        b = a.copy(); b[9:12, 9:12] = 0
        d = material_distances(a, b)
        self.assertEqual(d["phase_changed_cells"], 9)
        self.assertTrue(d["significant_phase_edit"])
        self.assertAlmostEqual(d["normalized_material_l2"], np.sqrt(9/576))
        b = a.copy(); b[0, 0] = 0
        self.assertFalse(material_distances(a, b)["significant_phase_edit"])


class ExchangeSemantics(unittest.TestCase):
    def test_source_packing_and_complex_rms_noise(self):
        z = np.array([[1+3j, 2+4j], [5+7j, 6+8j]])
        np.testing.assert_array_equal(pack_complex_per_source(z), [1, 5, 3, 7, 2, 6, 4, 8])
        expected = np.sqrt(2) / .2 * np.r_[1, 5, 3, 7, 2, 6, 4, 8]
        np.testing.assert_allclose(whitened_residual(z, 0*z, sigma=.2), expected)

    def test_no_unbudgeted_measurement(self):
        class NoPhysics:
            def predict(self, *args, **kwargs):
                raise AssertionError("a solver must not be touched")
        with self.assertRaises(ValueError):
            new_measurement(NoPhysics(), np.zeros(576), budget=None, acquisition={})

    def test_g_cannot_borrow_producer_budget_or_change_acquisition(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            book = SimpleNamespace(root=root, check=lambda: None)
            model = SimpleNamespace(book=book)
            args = dict(target_identity="frozen", observation_id="original", tau=10,
                        fixed_corrector_config=None, authorization_identity="fixture_only")
            with self.assertRaises(ValueError):
                FrozenPhysicsInterface(model, book, np.ones((32, 8)), .1,
                                       consumer_role="G", producer_ledger_root=root, **args)
            interface = FrozenPhysicsInterface(model, book, np.ones((32, 8)), .1, **args)
            with self.assertRaises(ValueError):
                interface.predict(np.ones(576), budget=book, acquisition={"frequency": 2})
            with self.assertRaises(ValueError):
                interface.predict(np.ones(576), budget=object())
            (root/"RUN_CLOSED.json").write_text("{}")
            with self.assertRaises(ValueError):
                interface.predict(np.ones(576), budget=book)

    def test_saved_array_roundtrip_and_fixture_rejection(self):
        chi = np.full((24, 24), 1.2 + .03j)
        coords = np.stack(np.meshgrid(np.arange(24)/24, np.arange(24)/24), axis=-1).reshape(576, 2)
        obs = np.full((3, 8), .05 + .01j)
        cost = dict(new_full_state_attempts=0, factorizations=0, forward_rhs=0,
                    tangent_rhs=0, adjoint_rhs=0, failed_attempts=0, setup_seconds=0, wall_seconds=0)
        with tempfile.TemporaryDirectory() as tmp:
            path = write_bundle(Path(tmp)/"fixture", candidate_id="SYNTHETIC_INTERFACE_FIXTURE",
                parent_id=None, source_revision="SYNTHETIC_NOT_RUN", chi=chi, observation=obs,
                prediction=obs+.001, coordinates=coords, mass=1/576, sigma=.1, tau=10,
                physics={"backend": "NO_SOLVER_FIXTURE", "dimension": 2}, observation_id="fixture",
                noise_id="fixture", legal_domain={"real_chi": [-.5, 4], "imag_chi": [0, .5]},
                validation={"status": "NOT_VALIDATED"}, cost=cost,
                lineage={"proposal_kind": "fixture_only"}, artifact_kind="INTERFACE_FIXTURE_NOT_A_PHYSICAL_CANDIDATE")
            self.assertFalse(validate_bundle(path)["valid_bundle_semantics"])
            self.assertTrue(validate_bundle(path, allow_fixture=True)["valid_bundle_semantics"])
            manifest = json.loads((path/"manifest.json").read_text())
            self.assertIsNone(manifest["probability_weight"])
            self.assertIsNone(manifest["optional_state_ref"])
            arrays = np.load(path/"prediction.npz")
            np.savez(path/"prediction.npz", channels=arrays["channels"]+1,
                     packed_whitened_residual=arrays["packed_whitened_residual"])
            self.assertFalse(validate_bundle(path, allow_fixture=True)["valid_bundle_semantics"])


if __name__ == "__main__":
    unittest.main()
