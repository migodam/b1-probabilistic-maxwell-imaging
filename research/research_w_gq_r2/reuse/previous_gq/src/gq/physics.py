"""Ordinary GQ physics interface over the frozen BGGA/W-SLAM P3 backend.

The caller supplies ``book.reserve(phase, label, rhs=..., factorizations=...,
**metadata)`` and ``book.event(phase, label, **metadata)``.  This module knows
neither the book's storage nor its limits.  CPU and process RSS belong to the
outer job wrapper.  Each LU/RHS is reserved before its numerical action.

``vendor`` must precede other BGGA/W-SLAM installations on the import path.
The public mesh, disk assembler, sparse LU, and disk full-Jacobian formulas
come from the pinned implementation.  Ellipses are affine images of the same
inscribed convex 384-gon and are integrated by triangle/polygon intersections.
No target nodal mask enters a matrix or a forcing integral.
"""
from __future__ import annotations

from dataclasses import dataclass, field as dataclass_field, replace
import time

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu
from skfem import asm

from bgga.fixed_mesh import (
    TARGET_POLYGON_N, TargetAssembler, TargetQuadrature,
    _area, _ccw, _clip_convex, _cross, public_fixed_mesh,
)
from bgga.fixed_physics import FixedEvaluator, FixedState
from bgga.operators import relative_backward
from wslam.fem_backend import FEM, _contrast_mass, _forcing, _sparse_bytes, line_green
from wslam.geometry import SHAPE, expand_support, validate_geometry
from wslam.world import polygon_inside, region_indices


BACKWARD_TOLERANCE = 1e-10
ELLIPSE_ASPECT = .65
S1_FREQUENCIES = (.85, 1., 1.15)
PUBLIC_MATERIAL_SHAPES = ((2, 2), (3, 3), SHAPE)
TARGET_GRADIENT_NAMES = ('px', 'py', 'contrast', 'target_loss', 'radius')


def _readonly(value, dtype=None):
    """An owned immutable array, including an immutable underlying buffer."""
    a = np.ascontiguousarray(np.asarray(value, dtype=dtype))
    return np.frombuffer(a.tobytes(), dtype=a.dtype).reshape(a.shape)


S0_TX = _readonly(((-1.5, 0.), (1.5, 0.)), float)
S1_TX = _readonly(((-1.5, 0.), (1.5, 0.), (-2.25, 0.), (2.25, 0.)), float)
RX16 = _readonly(np.column_stack((np.linspace(-2., 2., 16), np.zeros(16))))


def _positive(value, name):
    if isinstance(value, (bool, np.bool_)) or not np.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be positive and finite')
    return float(value)


def _grid_shape(shape):
    shape = tuple(shape)
    if (len(shape) != 2
            or any(isinstance(value, (bool, np.bool_))
                   or not isinstance(value, (int, np.integer)) for value in shape)
            or shape not in PUBLIC_MATERIAL_SHAPES):
        raise ValueError('material grid_shape must be (2,2), (3,3), or (6,24)')
    return tuple(int(value) for value in shape)


def _expanded_material(value, name):
    if np.iscomplexobj(value):
        raise ValueError(f'{name} must contain real coordinates')
    a = np.asarray(value, dtype=float)
    if a.ndim == 0:
        a = np.full(SHAPE, a.item())
    elif a.shape in PUBLIC_MATERIAL_SHAPES:
        a = np.repeat(np.repeat(a, SHAPE[0]//a.shape[0], axis=0),
                      SHAPE[1]//a.shape[1], axis=1)
    else:
        raise ValueError(f'{name} must be scalar or a public 2x2/3x3/6x24 array')
    if not np.isfinite(a).all():
        raise ValueError(f'{name} must be finite')
    return a


@dataclass(frozen=True, eq=False)
class World:
    """Legal, immutable wall/shape/pose/intrinsic coordinates.

    Wall arrays are canonical 6x24 values after public macro expansion; air is
    fixed to epsinf=1 and ell=0.  ``radius`` is the maximum semiaxis, in
    [.08,.16].  The ellipse is axis aligned, with x semiaxis radius and y
    semiaxis .65*radius.  Shape/size, pose, and target material remain separate
    metadata blocks.  No observed data or truth-file loader is present.
    """
    support: object
    epsinf: object
    ell: object
    pose: object
    radius: float = .12
    contrast: float = 3.
    target_loss: float = .02
    shape: str = 'disk'
    aspect: float = ELLIPSE_ASPECT

    def __post_init__(self):
        support = expand_support(self.support)
        eps = _expanded_material(self.epsinf, 'epsinf')
        loss = _expanded_material(self.ell, 'ell')
        pose = np.asarray(self.pose, dtype=float)
        radius = _positive(self.radius, 'radius')
        if not .08 <= radius <= .16:
            raise ValueError('maximum target semiaxis must be in [.08,.16]')
        if self.shape not in ('disk', 'ellipse'):
            raise ValueError("shape must be 'disk' or 'ellipse'")
        if self.shape == 'ellipse' and self.aspect != ELLIPSE_ASPECT:
            raise ValueError('only the frozen axis-aligned ellipse aspect=.65 is supported')
        if (not np.isfinite(self.contrast) or not 0. <= self.contrast <= 8.
                or not np.isfinite(self.target_loss) or not 0. <= self.target_loss <= .25):
            raise ValueError('target requires contrast in [0,8], loss in [0,.25]')
        # The circumscribing disk gives the same legal pose/size bounds and a
        # conservative wall-gap check for the ellipse.  No interface is moved.
        check = validate_geometry(support, eps, loss, pose, radius)
        if not check['legal']:
            raise ValueError('Illegal candidate: '+str(check['reasons']))
        for name, value, dtype in (
            ('support', check['support'], bool), ('epsinf', check['epsinf'], float),
            ('ell', check['ell'], float), ('pose', pose, float),
        ):
            object.__setattr__(self, name, _readonly(value, dtype))
        object.__setattr__(self, 'radius', radius)
        object.__setattr__(self, 'contrast', float(self.contrast))
        object.__setattr__(self, 'target_loss', float(self.target_loss))
        object.__setattr__(self, 'aspect', 1. if self.shape == 'disk' else ELLIPSE_ASPECT)

    @classmethod
    def snapshot(cls, world):
        if isinstance(world, cls):
            return world
        return cls(world.support, world.epsinf, world.ell, world.pose,
                   world.radius, world.contrast, world.target_loss,
                   getattr(world, 'shape', 'disk'), getattr(world, 'aspect', ELLIPSE_ASPECT))

    def cache_key(self):
        return (self.support.tobytes(), self.epsinf.tobytes(), self.ell.tobytes(),
                tuple(self.pose), self.radius, self.contrast, self.target_loss,
                self.shape, self.aspect)

    def occupied(self, points):
        iy, ix, inside = region_indices(points, SHAPE)
        return inside & self.support[iy, ix]

    def wall_values(self, points, f):
        f = _positive(f, 'frequency')
        iy, ix, inside = region_indices(points, SHAPE)
        occupied = inside & self.support[iy, ix]
        value = np.ones(len(points), complex)
        # World validates ell itself.  Preserve ell/f literally, including a
        # legal upper endpoint; do not perturb a coefficient for a subsequent
        # floating-point multiplication used only to re-check the same bound.
        imaginary = self.ell[iy, ix]/f
        value[occupied] = self.epsinf[iy[occupied], ix[occupied]]+1j*imaginary[occupied]
        return value

    def target_vertices(self, n=TARGET_POLYGON_N):
        if isinstance(n, (bool, np.bool_)) or int(n) != n or n < 3:
            raise ValueError('polygon vertex count must be an integer >=3')
        angle = 2*np.pi*np.arange(int(n))/int(n)
        return self.pose+self.radius*np.column_stack((np.cos(angle), self.aspect*np.sin(angle)))

    def target_contains(self, points, n=TARGET_POLYGON_N):
        # This is only a diagnostic coefficient query.  Target assembly uses
        # cut-polygon intersections and never these sampled occupancy values.
        points = np.asarray(points, float)
        local = (points-self.pose)/np.array((1., self.aspect))
        return polygon_inside(local, np.zeros(2), self.radius, n)

    def coefficient(self, f, n=TARGET_POLYGON_N, include_target=True):
        f = _positive(f, 'frequency')
        def value(points):
            result = self.wall_values(points, f)
            if include_target:
                result[self.target_contains(points, n)] = 1+self.contrast+1j*self.target_loss/f
            return result
        return value

    def metadata(self):
        return dict(support=self.support.astype(int).tolist(),
                    epsinf=self.epsinf.tolist(), ell=self.ell.tolist(),
                    target_pose=self.pose.tolist(),
                    target_shape=dict(kind=self.shape, maximum_semiaxis=self.radius,
                                      aspect=self.aspect, polygon_vertices=TARGET_POLYGON_N),
                    target_intrinsic=dict(contrast=self.contrast, loss=self.target_loss),
                    conductivity_convention='imaginary permittivity = ell/f')


class EllipseTargetAssembler(TargetAssembler):
    """Affine 384-gon cut cells; reuse the frozen P3 integral implementation."""
    def quadrature(self, pose, radius):
        center = np.asarray(pose, float)
        radius = _positive(radius, 'radius')
        if (center.shape != (2,) or not np.isfinite(center).all()
                or not -1 <= center[0] <= 1 or not 1.6 <= center[1] <= 2.6
                or not .08 <= radius <= .16):
            raise ValueError('ellipse requires legal pose and maximum semiaxis in [.08,.16]')
        semiaxes = radius*np.array((1., ELLIPSE_ASPECT))
        angle = 2*np.pi*np.arange(TARGET_POLYGON_N)/TARGET_POLYGON_N
        vertices = center+semiaxes*np.column_stack((np.cos(angle), np.sin(angle)))
        near = np.asarray(sorted(self._tree.query_ball_point(center, radius+self._max_radius)),
                          dtype=np.int64)
        if len(near):
            near = near[np.all(self._upper[near] > center-semiaxes, axis=1)
                        & np.all(self._lower[near] < center+semiaxes, axis=1)]
        pieces, touched, count = [], [], 0
        incircle = np.cos(np.pi/TARGET_POLYGON_N)
        for local in near:
            tri = _ccw(self._triangles[local])
            polygon = (tri if np.all(np.linalg.norm((tri-center)/semiaxes, axis=1) <= incircle)
                       else _clip_convex(vertices, tri))
            if _area(polygon) <= 0:
                continue
            element, found = int(self.roi_elements[local]), False
            for j in range(1, len(polygon)-1):
                a, b, c = polygon[0], polygon[j], polygon[j+1]
                determinant = float(_cross(b-a, c-a))
                if determinant <= 0:
                    continue
                count += len(self.reference_weights)
                if count > self.max_points or count*32 > self.max_work_bytes:
                    raise MemoryError('ellipse cut-polygon quadrature point/work-array limit exceeded')
                pieces.append((a.copy(), b.copy(), c.copy(), determinant, element))
                found = True
            if found:
                touched.append(element)
        if not pieces:
            raise RuntimeError('empty ellipse polygon/mesh intersection')
        points = np.empty((count, 2))
        weights = np.empty(count)
        elements = np.empty(count, dtype=np.int64)
        nq = len(self.reference_weights)
        for j, (a, b, c, determinant, element) in enumerate(pieces):
            sl = slice(j*nq, (j+1)*nq)
            points[sl] = a+self.reference_points[:, :1]*(b-a)+self.reference_points[:, 1:]*(c-a)
            weights[sl] = determinant*self.reference_weights
            elements[sl] = element
        area = .5*TARGET_POLYGON_N*np.prod(semiaxes)*np.sin(2*np.pi/TARGET_POLYGON_N)
        return TargetQuadrature(points, weights, elements, np.asarray(touched, dtype=np.int64),
                                float(weights.sum()), float(area), vertices,
                                self.reference_points.copy(), self.reference_weights.copy(),
                                self.degree, len(near))

    def assemble(self, pose, radius, k, chi, incident, *, source_count=None):
        target = super().assemble(pose, radius, k, chi, incident, source_count=source_count)
        metadata = dict(target.metadata)
        metadata.pop('circle_area', None)
        metadata.pop('polygon_relative_circle_area_deficit', None)
        metadata.update(target_shape='axis-aligned ellipse', aspect=ELLIPSE_ASPECT,
                        ellipse_area=np.pi*ELLIPSE_ASPECT*radius**2,
                        polygon_relative_ellipse_area_deficit=(
                            1-np.sin(2*np.pi/TARGET_POLYGON_N)/(2*np.pi/TARGET_POLYGON_N)))
        return replace(target, metadata=metadata)


class _BookAdapter:
    """Normalize the pinned FEM event spelling to the supplied public API."""
    def __init__(self, book):
        if not callable(getattr(book, 'reserve', None)) or not callable(getattr(book, 'event', None)):
            raise TypeError('external book.reserve and book.event are required')
        self.external = book

    def reserve(self, phase, label, rhs=0, factorizations=0, **metadata):
        return self.external.reserve(phase, label, rhs=rhs, factorizations=factorizations, **metadata)

    def event(self, phase=None, label='', **metadata):
        backend_event = metadata.pop('event', None)
        if label and backend_event:
            metadata['backend_event'] = backend_event
        metadata.setdefault('charged_cpu_seconds', 0.)
        return self.external.event(phase, label or backend_event or 'gq_event', **metadata)


def _defect(matrix, solution, rhs):
    """Maximum normwise backward error across RHS, with no denominator jitter."""
    u, b = np.asarray(solution), np.asarray(rhs)
    if u.ndim == 1:
        u, b = u[:, None], b[:, None]
    residual = matrix@u-b
    if not np.isfinite(u).all() or not np.isfinite(residual).all():
        raise RuntimeError('nonfinite GQ linear solution/residual')
    numerator = np.linalg.norm(residual, axis=0)
    denominator = np.linalg.norm(matrix.data)*np.linalg.norm(u, axis=0)+np.linalg.norm(b, axis=0)
    backward = np.full(len(numerator), np.inf)
    np.divide(numerator, denominator, out=backward, where=denominator != 0.)
    backward[(denominator == 0.) & (numerator == 0.)] = 0.
    result = relative_backward(matrix, solution, rhs)
    result.update(backward=float(backward.max()), per_rhs_backward=backward.tolist(),
                  tolerance=BACKWARD_TOLERANCE)
    if not np.isfinite(backward).all() or np.any(backward > BACKWARD_TOLERANCE):
        raise RuntimeError(f'GQ linear backward error exceeds {BACKWARD_TOLERANCE}: {result}')
    return result


def _immutable_sparse(matrix):
    matrix = matrix.tocsc(copy=True)
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    result = csc_matrix((_readonly(matrix.data), _readonly(matrix.indices),
                         _readonly(matrix.indptr)), shape=matrix.shape, copy=False)
    for value in (result.data, result.indices, result.indptr):
        value.flags.writeable = False
    return result


def _immutable_target(target):
    quadrature = replace(target.quadrature, **{
        name: _readonly(getattr(target.quadrature, name))
        for name in ('points', 'weights', 'elements', 'touched_elements', 'polygon_vertices',
                     'reference_points', 'reference_weights')})
    return replace(target, Vtarget=_immutable_sparse(target.Vtarget), bX=_readonly(target.bX),
                   quadrature=quadrature, **{
                       name: _readonly(getattr(target, name))
                       for name in ('touched_global_dofs', 'touched_free_indices',
                                    'public_roi_global_dofs', 'public_roi_free_indices')})


class _CheckedFactor:
    """Use the existing SuperLU factor; validate every direct/adjoint action."""
    def __init__(self, factor, matrix, book, phase):
        self._factor, self._matrix, self._book, self._phase = factor, matrix, book, phase
        self.last_defect = None

    def solve(self, rhs, trans='N'):
        if trans not in ('N', 'T', 'H'):
            raise ValueError('unsupported SuperLU transpose action')
        rhs = np.asarray(rhs, dtype=complex)
        if rhs.ndim not in (1, 2) or rhs.shape[0] != self._matrix.shape[0] or not np.isfinite(rhs).all():
            raise ValueError('linear RHS must have finite matching free-DOF shape')
        result = self._factor.solve(np.asfortranarray(rhs), trans=trans)
        matrix = (self._matrix if trans == 'N' else
                  self._matrix.T if trans == 'T' else self._matrix.conj().T)
        self.last_defect = _defect(matrix, result, rhs)
        self._book.event(self._phase, 'gq_checked_linear_action', transpose=trans,
                         solved_rhs=1 if rhs.ndim == 1 else rhs.shape[1],
                         backward=self.last_defect)
        return result


@dataclass(frozen=True, eq=False)
class PhysicsState:
    data: np.ndarray
    u: np.ndarray
    matrix: object
    rhs: np.ndarray
    factor: object
    target: object
    world: World
    residual: dict
    timings: dict
    epsilon: np.ndarray
    background_epsilon: np.ndarray
    material_points: np.ndarray
    material_weights: np.ndarray
    incident: np.ndarray
    field: np.ndarray
    current: np.ndarray
    scattered_dofs: np.ndarray
    scattered: np.ndarray
    direct: np.ndarray
    target_points: np.ndarray
    target_weights: np.ndarray
    target_incident: np.ndarray
    target_field: np.ndarray
    target_driver_field: np.ndarray
    target_current: np.ndarray
    mode: str
    metadata: dict
    background_u: object = None
    derivative_audit: dict = dataclass_field(default_factory=dict)
    _owner_token: object = dataclass_field(default=None, repr=False)

    @property
    def E(self):
        return self.field


class Evaluator(FixedEvaluator):
    """One public mesh, a generalized illumination, and one immutable state slot.

    ``f=1, tx=S0_TX, rx=RX16`` is S0.  Use one evaluator per supplied frequency
    with ``tx=S1_TX`` for S1; the caller owns ordering/stacking and bookkeeping.
    A cache hit returns the same immutable numerical state without new work.
    Full and Born states share the one slot.  Failed requests clear it first.
    """
    def __init__(self, book, phase='G0', *, f=1., tx=None, rx=None, h=.17,
                 wall_refinement=2, target_step=.125, quadrature_order=None,
                 pml=None, max_dofs=60000, max_factor_bytes=1_500_000_000,
                 target_degree=8, target_chunk_points=4096,
                 target_max_points=250_000, target_max_work_bytes=128_000_000):
        self.external_book = book
        self.book, self.phase = _BookAdapter(book), phase
        self.declared_completion = False
        self._owner_token = object()
        f = _positive(f, 'frequency')
        tx, rx = S0_TX if tx is None else tx, RX16 if rx is None else rx
        self.p = FEM(f, tx, rx, self.book, h=h, order=3, phase=phase,
                     mesh=public_fixed_mesh(h=h, grid_refinement=wall_refinement,
                                            target_step=target_step, pml=pml),
                     pml=pml, quadrature_order=quadrature_order,
                     residual_tolerance=BACKWARD_TOLERANCE,
                     max_dofs=max_dofs, max_factor_bytes=max_factor_bytes)
        self._target_settings = dict(degree=target_degree, chunk_points=target_chunk_points,
                                     max_points=target_max_points,
                                     max_work_bytes=target_max_work_bytes)
        self.ass = TargetAssembler(self.p.basis, self.p._free_dofs, **self._target_settings)
        self._ellipse_ass = None
        self.R = self.p._receiver_probe[:, self.p._free_dofs]
        self.direct = _readonly(line_green(self.p.rx, self.p.tx, self.p.k).T)
        self.incident = _readonly(line_green(self.p.material_points, self.p.tx, self.p.k))
        self._material_points = _readonly(self.p.material_points)
        self._material_weights = _readonly(self.p.material_weights)
        self._key = self._last = None

    def clear_cache(self):
        self._key = self._last = None

    def _state(self, state, *, full=True):
        if not isinstance(state, PhysicsState) or state._owner_token is not self._owner_token:
            raise ValueError('state belongs to a different GQ evaluator')
        if full and state.mode != 'full':
            raise ValueError('full-current derivatives require a full state')

    def _assembler(self, world):
        if world.shape == 'disk':
            return self.ass
        if self._ellipse_ass is None:
            self._ellipse_ass = EllipseTargetAssembler(self.p.basis, self.p._free_dofs,
                                                       **self._target_settings)
        return self._ellipse_ass

    def _target(self, world):
        return self._assembler(world).assemble(
            world.pose, world.radius, self.p.k,
            world.contrast+1j*world.target_loss/self.p.freq,
            lambda points: line_green(points, self.p.tx, self.p.k),
            source_count=len(self.p.tx))

    def background(self, world):
        """The pinned background weak forms with already-validated ell/f.

        GQ validates the real epsinf/ell coordinates before division.  Using
        them directly avoids both an endpoint coefficient adjustment and a
        false rejection from computing f*(ell/f) at ell=.25.
        """
        world = World.snapshot(world)
        p = self.p
        eps = world.wall_values(p.material_points, p.freq)
        mass = asm(_contrast_mass, p.material_basis,
                   contrast=(eps-1).reshape(p._material_shape))
        matrix = (p._background_matrix-p.k**2*mass[p._free_dofs][:, p._free_dofs]).tocsc()
        rhs = np.column_stack([
            asm(_forcing, p.material_basis,
                load=(p.k**2*(eps-1)*self.incident[:, j]).reshape(p._material_shape))[
                    p._free_dofs] for j in range(len(p.tx))])
        return matrix, rhs

    def _factor(self, matrix, label, *, rhs):
        if not np.isfinite(matrix.data).all():
            raise RuntimeError('nonfinite GQ FEM matrix')
        self.book.reserve(self.phase, label, factorizations=1, rhs=rhs,
                          freq=self.p.freq, declared_completion=self.declared_completion)
        factor = splu(matrix, permc_spec='COLAMD')
        storage = (_sparse_bytes(factor.L)+_sparse_bytes(factor.U)
                   + factor.perm_r.nbytes+factor.perm_c.nbytes)
        if storage > self.p.max_factor_bytes:
            raise MemoryError('actual GQ factor storage exceeds backend guard')
        matrix = _immutable_sparse(matrix)
        return _CheckedFactor(factor, matrix, self.book, self.phase), matrix, int(storage)

    def _interpolate(self, free_values):
        free_values = np.asarray(free_values)
        if free_values.ndim == 1:
            free_values = free_values[:, None]
        dofs = np.zeros((self.p.ndofs, free_values.shape[1]), complex)
        dofs[self.p._free_dofs] = free_values
        values = np.column_stack([
            np.asarray(self.p.material_basis.interpolate(dofs[:, j])).reshape(-1)
            for j in range(dofs.shape[1])])
        return values, dofs

    def _target_values(self, target, free_values):
        quadrature, outputs = target.quadrature, []
        for start in range(0, len(quadrature.weights), self.ass.chunk_points):
            end = min(start+self.ass.chunk_points, len(quadrature.weights))
            probe = self.ass._known_element_probes(
                quadrature.points[start:end], quadrature.elements[start:end], self.p._free_dofs)
            outputs.append(np.asarray(probe@free_values))
        return np.concatenate(outputs, axis=0)

    def _finish(self, world, matrix, rhs, factor, target, u, residual, started,
                storage, mode, *, background_u=None, extra_metadata=None):
        scattered_material, dofs = self._interpolate(u)
        field = self.incident+scattered_material
        background_eps = world.wall_values(self.p.material_points, self.p.freq)
        epsilon = world.coefficient(self.p.freq)(self.p.material_points)
        target_incident = line_green(target.quadrature.points, self.p.tx, self.p.k)
        target_field = target_incident+self._target_values(target, u)
        if background_u is None:
            driver_material, target_driver = field, target_field
        else:
            background_material, _ = self._interpolate(background_u)
            driver_material = self.incident+background_material
            target_driver = target_incident+self._target_values(target, background_u)
        chi = world.contrast+1j*world.target_loss/self.p.freq
        current = (background_eps-1)[:, None]*field
        current[world.target_contains(self.p.material_points)] += (
            chi*driver_material[world.target_contains(self.p.material_points)])
        scattered = np.asarray(self.R@u).T
        timings = dict(kernel_cpu_seconds=time.process_time()-started[0],
                       wall_seconds=time.perf_counter()-started[1], factor_bytes=storage)
        metadata = dict(freq=self.p.freq, tx=self.p.tx.tolist(), rx=self.p.rx.tolist(),
                        source_major=True, mode=mode,
                        equation=('(L_W-V_X)u=b_W+b_X' if mode == 'full' else
                                  'L_W u_W=b_W; L_W delta=b_X+V_X u_W; u=u_W+delta'),
                        target_recurrence=(mode == 'full'),
                        target_forcing='analytic incident integral at cut-polygon quadrature',
                        sampled_material_current='diagnostic sample; target cut-current arrays authoritative',
                        target_current_driver=('full total field' if mode == 'full' else
                                               'background total field'),
                        RSS_guard='caller-owned', qualification='NOT_RUN',
                        **(extra_metadata or {}))
        target = _immutable_target(target)
        state = PhysicsState(
            data=_readonly(self.direct+scattered), u=_readonly(u), matrix=matrix,
            rhs=_readonly(rhs), factor=factor, target=target, world=world,
            residual=residual, timings=timings, epsilon=_readonly(epsilon),
            background_epsilon=_readonly(background_eps), material_points=self._material_points,
            material_weights=self._material_weights, incident=self.incident, field=_readonly(field),
            current=_readonly(current), scattered_dofs=_readonly(dofs), scattered=_readonly(scattered),
            direct=self.direct, target_points=target.quadrature.points,
            target_weights=target.quadrature.weights, target_incident=_readonly(target_incident),
            target_field=_readonly(target_field), target_driver_field=_readonly(target_driver),
            target_current=_readonly(chi*target_driver), mode=mode, metadata=metadata,
            background_u=None if background_u is None else _readonly(background_u),
            _owner_token=self._owner_token)
        self.book.event(self.phase, 'gq_physics_state', mode=mode, freq=self.p.freq,
                         backward=residual, ndofs=self.p.ndofs, factor_bytes=storage,
                         kernel_cpu_seconds=timings['kernel_cpu_seconds'])
        return state

    def evaluate(self, world, *, label='gq full world'):
        try:
            world = World.snapshot(world)
        except Exception:
            self.clear_cache()
            raise
        key = ('full', world.cache_key())
        if key == self._key and self._last is not None:
            return self._last
        self.clear_cache()
        started = time.process_time(), time.perf_counter()
        matrix, wall_rhs = self.background(world)
        target = self._target(world)
        matrix, rhs = (matrix-target.Vtarget).tocsc(), wall_rhs+target.bX
        factor, matrix, storage = self._factor(matrix, label, rhs=len(self.p.tx))
        u = factor.solve(rhs)
        state = self._finish(world, matrix, rhs, factor, target, u, factor.last_defect,
                             started, storage, 'full')
        self._key, self._last = key, state
        return state

    def full_current_born(self, world, *, label='gq background Born'):
        """Exact wall propagation and incident forcing, one target insertion.

        The background LU, background source RHS, and Born insertion RHS are
        separately reserved.  The dropped term is V_X*delta, never b_X or
        V_X*u_W.  The residual reported in ``state.residual`` is for the Born
        equation; ``metadata.full_equation_defect`` records its full-equation
        mismatch without claiming that the Born state solves that equation.
        """
        try:
            world = World.snapshot(world)
        except Exception:
            self.clear_cache()
            raise
        key = ('background_born', world.cache_key())
        if key == self._key and self._last is not None:
            return self._last
        self.clear_cache()
        started = time.process_time(), time.perf_counter()
        matrix, wall_rhs = self.background(world)
        target = self._target(world)
        factor, matrix, storage = self._factor(matrix, label+' LU', rhs=0)
        self.book.reserve(self.phase, label+' background RHS', rhs=len(self.p.tx),
                          factorizations=0, freq=self.p.freq)
        background_u = factor.solve(wall_rhs)
        background_defect = dict(factor.last_defect)
        insertion = target.bX+target.Vtarget@background_u
        self.book.reserve(self.phase, label+' target insertion RHS', rhs=len(self.p.tx),
                          factorizations=0, freq=self.p.freq)
        delta = factor.solve(insertion)
        insertion_defect = dict(factor.last_defect)
        u, rhs = background_u+delta, wall_rhs+insertion
        residual = _defect(matrix, u, rhs)
        full_defect = relative_backward(matrix-target.Vtarget, u, wall_rhs+target.bX)
        state = self._finish(
            world, matrix, rhs, factor, target, u, residual, started, storage, 'background_born',
            background_u=background_u,
            extra_metadata=dict(background_backward=background_defect,
                                insertion_backward=insertion_defect,
                                full_equation_defect=full_defect))
        self._key, self._last = key, state
        return state

    born = full_current_born
    background_born = full_current_born

    def derivatives(self, state, *, include_nuisance=True):
        """Delegate the frozen full-Jacobian implementation for disks only.

        A columns are shared wall epsinf/ell.  B columns are
        (px,py,contrast,target_loss,radius,gammaR,gammaI).  Every native solve
        passes through the strict backward-error checker, including wall JVPs.
        """
        self._state(state)
        if state.world.shape != 'disk':
            raise ValueError('native full-Jacobian delegation is disk-only; use objective_gradient for ellipse')
        if include_nuisance and state.world.contrast == 0. and state.world.target_loss == 0.:
            raise ValueError('native zero-chi nuisance Jacobian unavailable; objective_gradient supports zero chi')
        native = FixedState(state.data, state.u, state.matrix, state.rhs, state.factor,
                            state.target, state.world, state.residual, state.timings)
        result = super().derivatives(native, include_nuisance=include_nuisance)
        state.derivative_audit.update(getattr(native, 'derivative_audit', {}))
        return result

    def wall_directional_jvp(self, state, eps_direction, ell_direction=0., *,
                             label='gq source-shared wall directional JVP'):
        """One real material direction, shared by every source; shape (ns,nr)."""
        self._state(state)
        eps = _expanded_material(eps_direction, 'eps direction')
        loss = _expanded_material(ell_direction, 'ell direction')
        iy, ix, inside = region_indices(self.p.material_points, SHAPE)
        occupied = inside & state.world.support[iy, ix]
        coefficient = occupied*(eps[iy, ix]+1j*loss[iy, ix]/self.p.freq)
        rhs = np.column_stack([
            asm(_forcing, self.p.material_basis,
                load=(self.p.k**2*coefficient*state.field[:, j]).reshape(self.p._material_shape))[
                    self.p._free_dofs] for j in range(len(self.p.tx))])
        self.book.reserve(self.phase, label, rhs=len(self.p.tx), factorizations=0,
                          freq=self.p.freq, direction_shared_across_sources=True)
        solved = state.factor.solve(rhs)
        return np.asarray(self.R@solved).T

    wall_jvp = wall_directional_jvp

    def _receiver_adjoint(self, state, label):
        rhs = self.R.conj().T.toarray()
        self.book.reserve(self.phase, label, rhs=len(self.p.rx), factorizations=0,
                          freq=self.p.freq, shared_across_material_cells=True)
        solved = state.factor.solve(rhs, trans='H')
        field, _ = self._interpolate(solved)
        return field

    def _cell_indices(self, world, grid_shape):
        iy, ix, inside = region_indices(self.p.material_points, grid_shape)
        occupied = inside & world.occupied(self.p.material_points)
        return iy*grid_shape[1]+ix, occupied

    def wall_jacobian(self, state, *, grid_shape=SHAPE, include_loss=True):
        """All public-cell columns from nr shared receiver-adjoint RHS.

        Rows flatten (source,receiver).  Columns are row-major epsinf cells,
        followed by row-major ell cells when include_loss=True.  Air-cell
        columns are zero.  The same material direction acts on all sources.
        """
        self._state(state)
        shape = _grid_shape(grid_shape)
        adjoint = self._receiver_adjoint(state, 'gq shared receiver adjoint for wall cells')
        cells, occupied = self._cell_indices(state.world, shape)
        accumulated = np.zeros((int(np.prod(shape)), len(self.p.tx), len(self.p.rx)), complex)
        for start in range(0, len(cells), self.ass.chunk_points):
            end = min(start+self.ass.chunk_points, len(cells))
            active = occupied[start:end]
            products = self.p.k**2*np.einsum(
                'qr,qs,q->qsr', adjoint[start:end].conj(), state.field[start:end],
                self.p.material_weights[start:end], optimize=True)
            np.add.at(accumulated, cells[start:end][active], products[active])
        eps = accumulated.transpose(1, 2, 0).reshape(state.data.size, -1)
        result = np.column_stack((eps, 1j/self.p.freq*eps)) if include_loss else eps
        state.derivative_audit['wall_jacobian'] = dict(
            grid_shape=list(shape), receiver_adjoint_rhs=len(self.p.rx),
            source_shared=True, ordering='source-major rows; row-major epsinf then ell columns',
            backward=state.factor.last_defect)
        return result

    def _weighted_adjoint(self, state, weights, label):
        weights = np.asarray(weights, dtype=complex)
        if weights.shape != state.data.shape or not np.isfinite(weights).all():
            raise ValueError('objective weights must be finite complex (source,receiver) values')
        rhs = np.asarray(self.R.conj().T@weights.T)
        self.book.reserve(self.phase, label, rhs=len(self.p.tx), factorizations=0,
                          freq=self.p.freq, weights_convention='delta_phi=Re(vdot(weights,delta_data))')
        solved = state.factor.solve(rhs, trans='H')
        field, _ = self._interpolate(solved)
        return solved, field

    def _wall_objective(self, state, adjoint_field, grid_shape):
        cells, occupied = self._cell_indices(state.world, grid_shape)
        products = self.p.k**2*self.p.material_weights*np.einsum(
            'qs,qs->q', adjoint_field.conj(), state.field, optimize=True)
        accumulated = np.zeros(int(np.prod(grid_shape)), complex)
        np.add.at(accumulated, cells[occupied], products[occupied])
        return (accumulated.real.reshape(grid_shape),
                (1j/self.p.freq*accumulated).real.reshape(grid_shape))

    def wall_objective_gradient(self, state, complex_residual_weights, *,
                                grid_shape=SHAPE, include_loss=True):
        """Scalar wall gradient using ns weighted-adjoint RHS, eps then ell."""
        self._state(state)
        shape = _grid_shape(grid_shape)
        _, adjoint = self._weighted_adjoint(state, complex_residual_weights,
                                             'gq weighted adjoint for wall objective')
        eps, loss = self._wall_objective(state, adjoint, shape)
        return np.r_[eps.ravel(), loss.ravel()] if include_loss else eps.ravel()

    def objective_gradient(self, state, complex_residual_weights, *, grid_shape=(3, 3)):
        """Shared ns-RHS wall and disk/ellipse target scalar objective gradient.

        The convention is delta_phi=Re(vdot(weights,delta_data)); least-squares
        weights can therefore be residual/sigma**2.  Target order is
        (px,py,contrast,target_loss,radius), with fixed ellipse aspect=.65.
        This method adds no prior, nuisance convention, or optimization step.
        """
        self._state(state)
        shape = _grid_shape(grid_shape)
        solved, adjoint = self._weighted_adjoint(state, complex_residual_weights,
                                                  'gq weighted adjoint for joint objective')
        adjoint_backward = dict(state.factor.last_defect)
        eps_gradient, loss_gradient = self._wall_objective(state, adjoint, shape)
        target_adjoint = self._target_values(state.target, solved)
        volume = self.p.k**2*np.einsum(
            'qs,qs,q->', target_adjoint.conj(), state.target_field,
            state.target_weights, optimize=True)
        w = state.world
        vertices = w.target_vertices()
        edge = np.roll(vertices, -1, axis=0)-vertices
        length = np.linalg.norm(edge, axis=1)
        nodes, weights = np.polynomial.legendre.leggauss(4)
        points = (vertices[:, None, :]+edge[:, None, :]*(nodes[None, :, None]+1)/2).reshape(-1, 2)
        normal = np.repeat(np.column_stack((edge[:, 1], -edge[:, 0]))/length[:, None], 4, axis=0)
        ds = (length[:, None]*weights[None, :]/2).ravel()
        radial_velocity = np.einsum('qi,qi->q', normal, (points-w.pose)/w.radius)
        probe = self.probe(points)
        total = line_green(points, self.p.tx, self.p.k)+probe@state.u
        adjoint_boundary = probe@solved
        boundary = self.p.k**2*(w.contrast+1j*w.target_loss/self.p.freq)*np.einsum(
            'qs,qs->q', adjoint_boundary.conj(), total, optimize=True)
        target_gradient = np.array((
            np.sum(ds*normal[:, 0]*boundary).real,
            np.sum(ds*normal[:, 1]*boundary).real,
            volume.real, (1j/self.p.freq*volume).real,
            np.sum(ds*radial_velocity*boundary).real))
        audit = dict(grid_shape=list(shape), weighted_adjoint_rhs=len(self.p.tx),
                     backward=adjoint_backward, target_shape=w.shape,
                     target_order=list(TARGET_GRADIENT_NAMES),
                     weights_convention='delta_phi=Re(vdot(weights,delta_data))',
                     target_material_integral='same independent cut-polygon quadrature',
                     target_shape_integral='same384gon boundary; degree4 Gaussian edge rule')
        state.derivative_audit['objective_gradient'] = audit
        return dict(epsinf=eps_gradient, ell=loss_gradient, target=target_gradient,
                    target_names=TARGET_GRADIENT_NAMES, audit=audit)


GQEvaluator = Evaluator
GQWorld = World
