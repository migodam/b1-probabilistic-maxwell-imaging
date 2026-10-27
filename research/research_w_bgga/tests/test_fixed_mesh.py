"""Mesh, polygon quadrature, mass/forcing tests only; no PDE or factorization.

These tests were written without being executed by the implementation worker.
The main campaign runner owns any test execution and CPU accounting.
"""
import inspect
import math
import unittest
from unittest.mock import patch

import numpy as np
from skfem import Basis, ElementTriP3, MeshTri

from bgga.fixed_mesh import (
    TARGET_BOX, TARGET_POLYGON_N, TargetAssembler, assemble_target,
    independent_triangle_rule, public_fixed_mesh, public_mesh_key,
)
from wslam.fem_backend import PMLConfig, _validate_mesh


def _polygon_x_moment(vertices, power):
    """Independent boundary formula int x^p dxdy = int_boundary x^(p+1)dy/(p+1)."""
    total = 0.
    for a, b in zip(vertices, np.roll(vertices, -1, axis=0)):
        dx, dy = b-a
        total += dy/(power+1)*sum(
            math.comb(power+1, j)*a[0]**(power+1-j)*dx**j/(j+1)
            for j in range(power+2))
    return total


class TestPublicMesh(unittest.TestCase):
    def test_public_parameters_and_no_overlapping_pslg_segments(self):
        import triangle

        captured = {}
        actual_triangulate = triangle.triangulate

        def capture(data, options):
            captured.update(data)
            return actual_triangulate(data, options)

        with patch.object(triangle, "triangulate", side_effect=capture):
            mesh = public_fixed_mesh(h=.35, grid_refinement=1, target_step=.4)
        self.assertEqual(set(inspect.signature(public_fixed_mesh).parameters),
                         {"h", "grid_refinement", "target_step", "pml"})
        self.assertEqual(mesh.bgga_public_key,
                         public_mesh_key(h=.35, grid_refinement=1, target_step=.4))
        _validate_mesh(mesh, PMLConfig().outer_box)
        vertices = captured["vertices"]
        groups = {}
        for first, second in captured["segments"]:
            a, b = vertices[first], vertices[second]
            # All public PSLG edges are axis aligned.  A real target polygon
            # would introduce oblique edges and fail this check.
            if a[0] == b[0]:
                key, interval = ("v", a[0]), sorted((a[1], b[1]))
            else:
                self.assertEqual(a[1], b[1])
                key, interval = ("h", a[1]), sorted((a[0], b[0]))
            groups.setdefault(key, []).append(interval)
        for intervals in groups.values():
            intervals.sort()
            for a, b in zip(intervals, intervals[1:]):
                self.assertLessEqual(a[1], b[0])
        triangles = mesh.p[:, mesh.t].transpose(2, 1, 0)
        first = triangles[:, 1]-triangles[:, 0]
        second = triangles[:, 2]-triangles[:, 0]
        area = .5*np.abs(first[:, 0]*second[:, 1]-first[:, 1]*second[:, 0]).sum()
        xl, xr, yb, yt = PMLConfig().outer_box
        np.testing.assert_allclose(area, (xr-xl)*(yt-yb), rtol=5e-13)
        self.assertEqual(TARGET_BOX, (-1.2, 1.2, 1.4, 2.8))


class TestTargetQuadrature(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # This is a small integration mesh, not a FEM/PDE backend.  Every DOF
        # is retained so the constant partition of unity is available.
        mesh = MeshTri.init_tensor(np.linspace(-1.4, 1.4, 8),
                                   np.linspace(1.2, 3., 7))
        cls.basis = Basis(mesh, ElementTriP3(), intorder=6)
        cls.free = np.arange(cls.basis.N)
        cls.assembler = TargetAssembler(cls.basis, cls.free, degree=8, chunk_points=211)

    def test_independent_triangle_rule_exact_monomials(self):
        points, weights = independent_triangle_rule(8)
        self.assertTrue(np.all(weights > 0))
        np.testing.assert_allclose(weights.sum(), .5, atol=2e-15)
        for a in range(9):
            for b in range(9-a):
                expected = math.factorial(a)*math.factorial(b)/math.factorial(a+b+2)
                actual = np.dot(weights, points[:, 0]**a*points[:, 1]**b)
                np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-16)

    def test_cut_polygon_area_for_extreme_poses_and_radii(self):
        for pose, radius in (((-1., 1.6), .16), ((1., 2.6), .16),
                             ((-1., 2.6), .06), ((1., 1.6), .06), ((.137, 2.13), .12)):
            with self.subTest(pose=pose, radius=radius):
                q = self.assembler.quadrature(pose, radius)
                expected = .5*TARGET_POLYGON_N*radius**2*np.sin(2*np.pi/TARGET_POLYGON_N)
                np.testing.assert_allclose(q.area, expected, rtol=5e-12, atol=5e-15)
                self.assertTrue(np.all(q.weights > 0))
                self.assertLess(q.area, np.pi*radius**2)
                np.testing.assert_allclose(np.sum(q.points*q.weights[:, None], axis=0),
                                           expected*np.asarray(pose), rtol=5e-12, atol=5e-15)
                self.assertLessEqual(q.candidate_count, len(self.assembler.roi_elements))
                self.assertTrue(set(q.touched_elements).issubset(self.assembler.roi_elements))

    def test_mass_and_direct_constant_linear_forcing(self):
        pose, radius, k, chi = (.137, 2.13), .12, 2.3, 3.+.02j

        def incident(points):
            return np.column_stack((np.ones(len(points)), points[:, 0]+1j*points[:, 1]))

        result = self.assembler.assemble(pose, radius, k, chi, incident, source_count=2)
        mass = result.Vtarget/(k*k*chi)
        self.assertEqual(result.Vtarget.format, "csc")
        self.assertEqual(result.bX.shape, (self.basis.N, 2))
        np.testing.assert_allclose(mass.toarray(), mass.toarray().T, atol=3e-15)
        np.testing.assert_allclose(np.ones(self.basis.N)@mass@np.ones(self.basis.N),
                                   result.target_area, rtol=3e-12, atol=2e-15)
        expected = k*k*chi*result.target_area*np.asarray((1., pose[0]+1j*pose[1]))
        np.testing.assert_allclose(result.bX.sum(axis=0), expected, rtol=5e-12, atol=5e-15)
        dof_field = np.column_stack((np.ones(self.basis.N),
                                    self.basis.doflocs[0]+1j*self.basis.doflocs[1]))
        np.testing.assert_allclose(result.Vtarget@dof_field, result.bX, rtol=2e-11, atol=2e-14)
        active = result.touched_free_indices
        eigenvalues = np.linalg.eigvalsh(mass[active][:, active].toarray().real)
        self.assertGreaterEqual(eigenvalues.min(), -2e-14)
        self.assertEqual(result.touched_rank, len(active))
        self.assertTrue(set(active).issubset(result.public_roi_free_indices))
        outside = np.setdiff1d(np.arange(self.basis.N), active)
        self.assertEqual(result.Vtarget[outside].nnz, 0)
        self.assertTrue(np.all(result.bX[outside] == 0))

    def test_quartic_incident_is_integrated_without_nodal_interpolation(self):
        result = assemble_target(self.basis, self.free, (.137, 2.13), .12, 1., 1.,
                                 lambda points: points[:, 0]**4, assembler=self.assembler)
        # Green's theorem gives an independent polygon moment, while a P3
        # interpolant of x^4 is not the analytic quartic used by the load.
        expected = _polygon_x_moment(result.quadrature.polygon_vertices, 4)
        np.testing.assert_allclose(result.bX.sum(), expected, rtol=2e-11, atol=2e-15)
        interpolated_load = result.Vtarget@(self.basis.doflocs[0]**4)
        self.assertGreater(np.linalg.norm(result.bX[:, 0]-interpolated_load), 1e-9)

    def test_public_ports_survive_a_target_change_and_zero_contrast(self):
        first = self.assembler.assemble((-.9, 1.7), .06, 1., 0., lambda p: np.ones(len(p)))
        second = self.assembler.assemble((.9, 2.5), .16, 1., 1., lambda p: np.ones(len(p)))
        np.testing.assert_array_equal(first.public_roi_free_indices, second.public_roi_free_indices)
        np.testing.assert_array_equal(first.public_roi_global_dofs, second.public_roi_global_dofs)
        self.assertGreater(first.touched_rank, 0)
        self.assertEqual(first.Vtarget_rank, 0)
        self.assertEqual(first.Vtarget.nnz, 0)
        self.assertTrue(np.all(first.bX == 0))
        self.assertNotEqual(set(first.touched_free_indices), set(second.touched_free_indices))

    def test_point_limit_rejects_before_probe_allocation(self):
        limited = TargetAssembler(self.basis, self.free, max_points=10)
        with self.assertRaises(MemoryError):
            limited.quadrature((0., 2.), .12)


if __name__ == "__main__":
    unittest.main()
