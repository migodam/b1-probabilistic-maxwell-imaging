"""Physics/convention regressions; also runnable via unittest without pytest."""
import unittest

import numpy as np

from research_w.analytic import (
    WALL1, WALL2, amplitude_transfer, fixed_plane_transmission, fresnel_slab,
    full_target_reflection, pec_reflection_direct, realify, reflection_recursion,
    scattering, whitened_distance,
)


class AnalyticRegression(unittest.TestCase):
    def test_fresnel_transfer_and_asymmetric_interface_recursion(self):
        f = np.linspace(0.8, 1.2, 1001)
        for wall in (WALL1, WALL2):
            rl, rr, t = scattering(wall, f)
            fr, ft = fresnel_slab(*wall[0], f)
            np.testing.assert_allclose(rl, fr, atol=2e-14, rtol=0)
            np.testing.assert_allclose(t, ft, atol=2e-14, rtol=0)
            np.testing.assert_allclose(rl, rr, atol=2e-14, rtol=0)
        wall = ((1.0, 0.1), (2.3, 0.12), (1.4, 0.08))
        np.testing.assert_allclose(scattering(wall, f)[0], reflection_recursion(wall, f), atol=2e-14, rtol=0)
        np.testing.assert_allclose(scattering(wall[::-1], f)[0], scattering(wall, f)[1], atol=2e-14, rtol=0)

    def test_lossless_unitarity_and_reciprocity(self):
        f = np.linspace(0.8, 1.2, 1001)
        for wall in (WALL1, WALL2, ((1.0, 0.1), (2.3, 0.12), (1.4, 0.08))):
            rl, rr, t = scattering(wall, f)
            s = np.stack((np.stack((rl, t), -1), np.stack((t, rr), -1)), -2)
            np.testing.assert_allclose(s.conj().swapaxes(-1, -2) @ s,
                                       np.broadcast_to(np.eye(2), s.shape), atol=2e-14, rtol=0)
            np.testing.assert_allclose(np.linalg.det(amplitude_transfer(wall, f)), 1, atol=2e-14, rtol=0)

    def test_fixed_reference_plane_and_center_witness(self):
        for wall in (WALL1, WALL2):
            rl, rr, t = scattering(wall, [1.0])
            np.testing.assert_allclose([rl, rr], 0, atol=1e-14, rtol=0)
            np.testing.assert_allclose(t, -1, atol=1e-14, rtol=0)
        t1, t2 = (fixed_plane_transmission(wall, [1.0]) for wall in (WALL1, WALL2))
        np.testing.assert_allclose(t2 / t1, np.exp(1j * np.pi / 6), atol=1e-14, rtol=0)
        np.testing.assert_allclose((t2 / t1) ** 2, np.exp(1j * np.pi / 3), atol=1e-14, rtol=0)
        self.assertAlmostEqual(float(abs(t1 ** 2 - t2 ** 2)[0]), 1.0, places=13)
        # Air-only phase is exp(+ikD) for the stated temporal convention.
        f = np.array([0.91, 1.0, 1.13])
        np.testing.assert_allclose(fixed_plane_transmission((), f), np.exp(2j * np.pi * f), atol=1e-14, rtol=0)

    def test_complete_pec_cascade_and_position_compensation(self):
        f = np.linspace(0.8, 1.2, 1001)
        for wall in (WALL1, WALL2, ((1.0, 0.1), (2.3, 0.12), (1.4, 0.08))):
            y = full_target_reflection(wall, f, 1.5)
            np.testing.assert_allclose(y, pec_reflection_direct(wall, f, 1.5), atol=2e-14, rtol=0)
            np.testing.assert_allclose(abs(y), 1, atol=2e-14, rtol=0)
        y1 = full_target_reflection(WALL1, [1.0], 1.5)
        y2 = full_target_reflection(WALL2, [1.0], 1.5 - 1 / 12)
        np.testing.assert_allclose(y1, y2, atol=1e-14, rtol=0)
        self.assertGreater(abs(full_target_reflection(WALL1, [0.9], 1.5)[0] -
                               scattering(WALL1, [0.9])[2][0] ** 2 * -np.exp(4j * np.pi * 0.9 * 1.25)), 0.1)

    def test_nine_frequency_band_and_real_quadrature_noise(self):
        np.testing.assert_array_equal(realify([1 + 2j, 3 + 4j]), [1, 3, 2, 4])
        self.assertAlmostEqual(whitened_distance([0], [0.01 + 0.01j]), np.sqrt(2))
        f = np.linspace(1 - 0.005 / 2, 1 + 0.005 / 2, 9)
        wall_d = whitened_distance(scattering(WALL1, f)[0], scattering(WALL2, f)[0])
        joint_d = whitened_distance(full_target_reflection(WALL1, f, 1.5),
                                    full_target_reflection(WALL2, f, 1.5 - 1 / 12))
        self.assertAlmostEqual(wall_d, 0.88713775173482, places=11)
        self.assertAlmostEqual(joint_d, 3.0382462624899644, places=10)

    def test_domain_and_reference_plane_checks(self):
        with self.assertRaises(ValueError):
            scattering(((0, 0.2),), [1])
        with self.assertRaises(ValueError):
            scattering(WALL1, [0])
        with self.assertRaises(ValueError):
            fixed_plane_transmission(WALL1, [1], plane=0.1)
        with self.assertRaises(ValueError):
            pec_reflection_direct(WALL1, [1], position=0.1)


if __name__ == "__main__":
    unittest.main()
