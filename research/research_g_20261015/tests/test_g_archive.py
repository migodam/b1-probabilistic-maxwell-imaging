"""Substantive archive invariants on synthetic fields; no Maxwell calls."""
import itertools
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from g_archive import (Archive, Candidate, OfferBudgetExhausted, material_distance,
                       structure_features, structural_distinct, structural_distance,
                       topology_summary)


def field(value=0.):
    return np.full((24, 24), value, np.complex128)


def candidate(name, value=0., residual=2., **kwargs):
    image = field(value) if np.ndim(value) == 0 else value
    defaults = dict(full_validated=True, solver_validated=True,
                    validation_receipt="synthetic-external-receipt:" + name)
    defaults.update(kwargs)
    return Candidate(name, image, np.zeros((32, 8), np.complex128), residual, **defaults)


class ArchiveInvariants(unittest.TestCase):
    def archive(self, **kwargs):
        return Archive(field(), data_threshold=3., **kwargs)

    def test_material_metric_has_fixed_real_and_imag_scales(self):
        self.assertAlmostEqual(material_distance(field(), field(.45 + .05j)), np.sqrt(.02))

    def test_physical_duplicate_survives_in_success_ledger(self):
        archive = self.archive()
        first = candidate("first", .4)
        second = candidate("second", first.field.copy(), residual=2.1)
        self.assertEqual(archive.offer(first).new_world_reward, 1)
        result = archive.offer(second)
        self.assertEqual(result.new_world_reward, 0)
        self.assertEqual(len(archive.success_ledger), 2)
        self.assertEqual(len(archive.active_archive), 1)

    def test_identical_receiver_response_does_not_merge_different_material(self):
        archive = self.archive()
        archive.offer(candidate("a", .1))
        archive.offer(candidate("b", .6))
        self.assertEqual(len(archive.active_archive), 2)
        self.assertEqual(len(archive.material_representatives), 2)

    def test_single_neighbor_lower_residual_replacement_retains_old_success(self):
        archive = self.archive()
        archive.offer(candidate("a", 0., residual=2.))
        archive.offer(candidate("far", .4, residual=2.))
        replacement = archive.offer(candidate("near", .01, residual=1.))
        self.assertEqual(replacement.replaced_id, "a")
        self.assertEqual(replacement.new_world_reward, 0)
        self.assertEqual({item.candidate_id for item in archive.active_archive}, {"near", "far"})
        reinsert = archive.offer(candidate("old-world", 0., residual=.5))
        self.assertEqual(reinsert.new_world_reward, 0)
        self.assertEqual(len(archive.success_ledger), 4)
        archive.assert_pairwise_packing()

    def test_bridge_cannot_single_link_merge_two_separated_worlds(self):
        archive = self.archive()
        archive.offer(candidate("a", 0.))
        archive.offer(candidate("c", .2))
        result = archive.offer(candidate("bridge", .1, residual=.1))
        self.assertFalse(result.active)
        self.assertEqual({item.candidate_id for item in archive.active_archive}, {"a", "c"})
        self.assertEqual(len(archive.success_ledger), 3)

    def test_one_cell_edit_is_not_a_new_structure(self):
        background = field()
        image = field()
        image[10, 10] = 1.2
        first, second = structure_features(background, background), structure_features(image, background)
        self.assertEqual(structural_distance(first, second)[1], 1)
        self.assertFalse(structural_distinct(first, second))
        archive = self.archive()
        archive.offer(candidate("base", background))
        self.assertFalse(archive.offer(candidate("one-cell", image)).active)

    def test_spatially_resolved_hole_is_detected_posthoc_without_hole_input(self):
        solid, edited = field(1.2), field(1.2)
        edited[9:15, 9:15] = 0.
        first, second = structure_features(solid, field()), structure_features(edited, field())
        self.assertTrue(structural_distinct(first, second))
        summary = topology_summary(second)
        self.assertEqual(summary["holes_ge6"], 1)
        self.assertEqual(summary["hole_sizes"], [36])

    def test_negative_complex_contrast_is_not_declared_empty_support(self):
        features = structure_features(field(-.4 + .1j), field())
        self.assertEqual(int(features.original_mask.sum()), 576)

    def test_invalid_and_unverified_never_install(self):
        archive = self.archive()
        archive.offer(candidate("good", .2))
        for item in (candidate("box", 4.1), candidate("residual", .5, residual=4.),
                     candidate("not-full", .6, full_validated=False),
                     candidate("not-solver", .7, solver_validated=False),
                     candidate("no-receipt", .8, validation_receipt=None)):
            self.assertFalse(archive.offer(item).active)
        self.assertEqual(len(archive.success_ledger), 1)
        self.assertEqual(len(archive.failure_ledger), 5)
        self.assertEqual(len(archive.frontier), 5)
        self.assertEqual(archive.top1().candidate_id, "good")

    def test_offer_budget_counts_failures_and_keeps_verified_results(self):
        archive = self.archive(max_offers=2)
        archive.offer(candidate("base", .2, is_base=True))
        archive.offer(candidate("failure", 4.1))
        with self.assertRaises(OfferBudgetExhausted):
            archive.offer(candidate("not-run", .7))
        self.assertEqual(archive.offers, 2)
        self.assertEqual(len(archive.success_ledger), 1)
        self.assertEqual(archive.topk(8)[0].candidate_id, "base")

    def test_topk_preserves_legal_base_and_does_not_pad_duplicates(self):
        archive = self.archive()
        archive.offer(candidate("base", 0., residual=2.5, is_base=True))
        archive.offer(candidate("base-replacement", .01, residual=1.5))
        archive.offer(candidate("far", .5, residual=.5))
        selected = archive.topk(8)
        self.assertEqual([item.candidate_id for item in selected], ["base", "far"])
        self.assertEqual(archive.top1().candidate_id, "far")

    def test_input_and_saved_arrays_are_immutable(self):
        image = field(.3)
        item = candidate("immutable", image)
        image[:] = 3.
        self.assertTrue(np.all(item.field == .3))
        with self.assertRaises(ValueError):
            item.field[0, 0] = 2.
        with self.assertRaises(ValueError):
            item.field.setflags(write=True)

    def test_packing_invariant_holds_across_all_small_input_permutations(self):
        rows = [("a", 0.), ("b", .1), ("c", .2), ("d", .3)]
        for permutation in itertools.permutations(rows):
            archive = self.archive()
            for name, value in permutation:
                archive.offer(candidate(name, value))
            self.assertTrue(archive.assert_pairwise_packing())
            self.assertEqual(len(archive.success_ledger), 4)

    def test_digital_connectivities_are_complementary_and_size_is_posthoc(self):
        image = field()
        image[10, 10] = image[11, 11] = 1.2
        summary = topology_summary(structure_features(image, field()))
        self.assertEqual(summary["components_all"], 2)
        self.assertEqual(summary["components_ge6"], 0)
        self.assertEqual(summary["holes_all"], 0)
        self.assertEqual((summary["foreground_connectivity"], summary["background_connectivity"]), (4, 8))

    def test_probability_and_latent_identity_not_invented(self):
        archive = self.archive()
        item = candidate("plain", .4)
        archive.offer(item)
        for key in ("posterior_weight", "represented_mass", "unmodeled_mass"):
            self.assertIsNone(getattr(item, key))
            self.assertIsNone(archive.summary()[key])

    def test_duplicate_receipt_id_cannot_overwrite_success(self):
        archive = self.archive()
        archive.offer(candidate("same-id", .2))
        with self.assertRaises(ValueError):
            archive.offer(candidate("same-id", .8))
        self.assertEqual(len(archive.success_ledger), 1)
        self.assertTrue(np.all(archive.success_ledger[0].field == .2))

    def test_first_vs_quality_changes_only_representative_replacement(self):
        archives = [self.archive(replacement_mode=mode) for mode in ("first", "quality")]
        for archive in archives:
            archive.offer(candidate("early", .2, residual=2.))
            archive.offer(candidate("better", .21, residual=1.))
        self.assertEqual(archives[0].active_archive[0].candidate_id, "early")
        self.assertEqual(archives[1].active_archive[0].candidate_id, "better")
        self.assertEqual(len(archives[0].success_ledger), len(archives[1].success_ledger))
        self.assertEqual([item.new_world_reward for item in archives[0].event_ledger],
                         [item.new_world_reward for item in archives[1].event_ledger])


if __name__ == "__main__":
    unittest.main()
