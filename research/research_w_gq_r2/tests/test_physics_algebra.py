"""Zero-PDE geometry, accounting and discrete product-kernel identities.

FEM construction and sparse LU are forbidden.  A tiny P3 Basis is used only
for interpolation and cut-integration algebra, with no Helmholtz operator.
"""
from contextlib import ExitStack
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from scipy.sparse import eye
from scipy.spatial import cKDTree
from skfem import Basis, ElementTriP3, MeshTri

from bgga.fixed_mesh import independent_triangle_rule
from r2 import physics


class ProbeBook:
    def __init__(self):
        self.reservations, self.events = [], []

    def reserve(self, phase, label, **metadata):
        self.reservations.append((phase, label, metadata))

    def event(self, phase, label, **metadata):
        self.events.append((phase, label, metadata))


class PhysicsAlgebraTests(unittest.TestCase):
    def setUp(self):
        def forbidden(*args, **kwargs):
            raise AssertionError('zero-PDE tests cannot construct FEM or call sparse LU')
        guards = ExitStack()
        self.addCleanup(guards.close)
        guards.enter_context(patch.object(physics, 'FEM', forbidden))
        guards.enter_context(patch.object(physics, 'splu', forbidden))

    def test_wall_coordinates_and_literal_loss_preserve_public_cells(self):
        x, y = np.linspace(-1.5, 1.5, 7), np.linspace(.5, 1.25, 4)
        eps, loss = np.full((3, 6), 2.), np.full((3, 6), .25)
        eps[0, 0], loss[0, 0] = 1., 0.
        wall = physics.WallGrid(x, y, eps, loss)
        x[1], eps[1, 1], loss[1, 1] = -1.4, 6., 0.
        self.assertEqual(wall.shape, (3, 6))
        self.assertEqual(wall.epsinf[1, 1], 2.)
        for f in (.85, 1., 1.15):
            values = wall.values(np.array(((0., .9), (0., 2.), (-1.49, .51))), f)
            self.assertEqual(values[0].imag, .25/f)
            np.testing.assert_array_equal(values[1:], (1., 1.))
        with self.assertRaises(ValueError):
            wall.epsinf.flags.writeable = True
        for bad_eps, bad_loss in ((1.2, 0.), (6.1, 0.), (1., .01), (2., .251)):
            with self.subTest(eps=bad_eps, ell=bad_loss), self.assertRaises(ValueError):
                physics.WallGrid([-1.5, 1.5], [.5, 1.25], [[bad_eps]], [[bad_loss]])

    def test_atomic_pslg_keeps_every_requested_wall_interface(self):
        x, y = np.array((-1.5, -.4, .2, 1.5)), np.array((.5, .8, 1.25))
        graph, key = physics._public_pslg(x, y, wall_refine=2)
        gx, gy = physics._refined_axis(x, 2), physics._refined_axis(y, 2)
        vertices = graph['vertices']
        segments = {frozenset((tuple(vertices[a]), tuple(vertices[b]))) for a, b in graph['segments']}
        expected = []
        for yi in gy:
            expected.extend(frozenset(((a, yi), (b, yi))) for a, b in zip(gx[:-1], gx[1:]))
        for xi in gx:
            expected.extend(frozenset(((xi, a), (xi, b))) for a, b in zip(gy[:-1], gy[1:]))
        self.assertTrue(set(expected) <= segments)
        self.assertEqual(key[0], physics.MESH_VERSION)
        self.assertEqual(key, physics._public_pslg(x.copy(), y.copy(), wall_refine=2)[1])

    def test_affine_cut_quadrature_partitions_both_orientations(self):
        triangles = np.array((((-1.2, 1.4), (1.2, 1.4), (1.2, 2.8)),
                              ((-1.2, 1.4), (1.2, 2.8), (-1.2, 2.8))))
        for orientation in (0., np.pi/2):
            target = physics.Target((0., 2.1), .12, shape='ellipse', orientation=orientation)
            extent = np.max(np.abs(target.target_vertices()-target.pose), axis=0)
            np.testing.assert_allclose(extent, (.12, .08) if orientation == 0. else (.08, .12), atol=1e-15)
            assembler = object.__new__(physics.AffineTargetAssembler)
            assembler.aspect, assembler.orientation = target.aspect, target.orientation
            assembler._transform = target.transform()
            assembler._triangles, assembler.roi_elements = triangles, np.array((0, 1))
            assembler._lower, assembler._upper = triangles.min(axis=1), triangles.max(axis=1)
            centers = triangles.mean(axis=1)
            assembler._tree = cKDTree(centers)
            assembler._max_radius = float(np.linalg.norm(triangles-centers[:, None], axis=2).max())
            assembler.reference_points, assembler.reference_weights = independent_triangle_rule(8)
            assembler.degree, assembler.max_points, assembler.max_work_bytes = 8, 250000, 128000000
            q = assembler.quadrature(target.pose, target.radius)
            self.assertEqual(set(q.touched_elements), {0, 1})
            self.assertLess(abs(q.area-q.polygon_area), 1e-13)
            self.assertTrue(np.all(q.weights > 0.))

    def test_product_kernel_equals_source_resolved_mass_rhs_duality(self):
        # Independent explicit FE-quadrature construction checks conjugation,
        # source-major row ordering, real arbitrary directions and ell/f.
        rng = np.random.default_rng(2601029)
        probe = rng.normal(size=(9, 5))
        u = rng.normal(size=(5, 2))+1j*rng.normal(size=(5, 2))
        inc = rng.normal(size=(9, 2))+1j*rng.normal(size=(9, 2))
        adjoint = rng.normal(size=(5, 3))+1j*rng.normal(size=(5, 3))
        weights = rng.uniform(.01, .2, size=9)
        cells = np.array((0, 0, 1, 1, 2, 2, 3, 3, 3))
        active = cells != 2
        columns = np.array(((1., .2), (-.4, 1.), (.3, .1), (.8, -.7)))
        total, adj = inc+probe@u, probe@adjoint
        kernel = physics._cell_product_kernel(adj, total, weights, cells, active, 4, 7.)
        product = kernel.transpose(1, 2, 0).reshape(6, 4)@columns
        explicit = []
        for column in columns.T:
            coefficient = np.where(active, column[cells], 0.)
            mass = 7.*probe.T@np.diag(weights*coefficient)@probe
            forcing = 7.*probe.T@(weights[:, None]*coefficient[:, None]*inc)
            explicit.append((adjoint.conj().T@(mass@u+forcing)).T.ravel())
        np.testing.assert_allclose(product, np.column_stack(explicit), rtol=2e-14, atol=2e-14)
        for f in (.85, 1., 1.15):
            projected = 1j/f*product
            scalar = np.vdot(np.arange(6)+1j, projected[:, 0]).real
            direct_scalar = np.vdot((np.arange(6)+1j).reshape(2, 3),
                                    (1j/f*np.column_stack(explicit)[:, 0]).reshape(2, 3)).real
            self.assertAlmostEqual(scalar, direct_scalar, places=12)

    def test_bounded_interpolation_matches_existing_basis_with_boundary_zeros(self):
        mesh = MeshTri.init_tensor(np.array((-1., 0., 1.)), np.array((1.4, 2., 2.8)))
        basis = Basis(mesh, ElementTriP3(), intorder=8)
        free = np.setdiff1d(np.arange(basis.N), basis.get_dofs().all())
        evaluator = object.__new__(physics.Eval)
        evaluator.material_chunk_elements = 3
        evaluator._to_free = np.full(basis.N, -1, dtype=int)
        evaluator._to_free[free] = np.arange(len(free))
        evaluator.p = SimpleNamespace(material_basis=basis, _free_dofs=free,
                                      _material_shape=np.asarray(basis.global_coordinates()).shape[1:])
        rng = np.random.default_rng(2601028)
        u = rng.normal(size=(len(free), 2))+1j*rng.normal(size=(len(free), 2))
        chunked = np.concatenate([fields[0] for _, fields in evaluator._material_chunks(u)], axis=0)
        dofs = np.zeros((basis.N, 2), complex)
        dofs[free] = u
        expected = np.column_stack([np.asarray(basis.interpolate(dofs[:, j])).ravel() for j in range(2)])
        np.testing.assert_allclose(chunked, expected, rtol=0., atol=0.)

    def test_every_factor_action_reserves_before_solve_and_checks_each_rhs(self):
        book = ProbeBook()
        class IdentityAction:
            def solve(self, rhs, trans='N'):
                self_test.assertEqual(book.reservations[-1][2]['rhs'], rhs.shape[1])
                self_test.assertEqual(book.reservations[-1][2]['transpose'], trans)
                return rhs.copy()
        self_test = self
        factor = physics._MeteredFactor(IdentityAction(), eye(2, format='csc', dtype=complex), book, 'E1', 1.)
        rhs = np.array(((1., 1e-20), (0., 0.)), complex)
        for trans in ('N', 'T', 'H'):
            np.testing.assert_array_equal(factor.solve(rhs, trans=trans), rhs)
        self.assertEqual(sum(row[2]['rhs'] for row in book.reservations), 6)
        self.assertEqual(sum(row[2]['factorizations'] for row in book.reservations), 0)
        class WrongAction:
            def solve(self, rhs, trans='N'):
                result = rhs.copy()
                result[0, 1] += 1e-25
                return result
        wrong = physics._MeteredFactor(WrongAction(), eye(2, format='csc', dtype=complex), book, 'E1', 1.)
        with self.assertRaises(RuntimeError):
            wrong.solve(rhs)


if __name__ == '__main__':
    unittest.main(verbosity=2)
