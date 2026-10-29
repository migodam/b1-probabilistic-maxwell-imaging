"""R2 interface-fitted P3 adapter over the read-only GQ/BGGA/W-SLAM backend.

The public wall axes, their optional refinement, and the target ROI grid are
the only mesh inputs.  Material values and target parameters never enter mesh
construction.  All numerical LU and RHS actions are reserved before use, and
each solved RHS must meet the inherited 1e-10 normwise backward-error limit.

States retain only the free scattered-field DOFs, receiver data, sparse
operator, forcing, and cut-target assembly.  Material and adjoint fields are
interpolated in bounded element chunks rather than copied into every state.
No qualification result or scientific acceptance is supplied by this module.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
import time

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu
from skfem import MeshTri, asm

from bgga.fixed_mesh import (
    TARGET_BOX, WALL_BOX, TARGET_POLYGON_N, TargetAssembler, TargetQuadrature,
    _area, _ccw, _clip_convex, _cross, _uniform_axis,
)
from bgga.fixed_physics import FixedEvaluator
from gq.physics import (
    BACKWARD_TOLERANCE, RX16, S0_TX, S1_TX, _BookAdapter, _defect,
    _immutable_sparse, _immutable_target, _positive, _readonly,
)
from wslam.fem_backend import (
    FEM, MATERIAL_BOX, PMLConfig, _contrast_mass, _forcing, _sparse_bytes, line_green,
)


MESH_VERSION = 'R2_PUBLIC_WALL_PSLG_V1'
TARGET_NAMES = ('px', 'py', 'contrast', 'loss', 'radius')
ELLIPSE_ASPECT = 2/3


def _axis(value, endpoints, name):
    if np.iscomplexobj(value):
        raise ValueError(f'{name} must be real')
    a = np.asarray(value, dtype=float)
    if (a.ndim != 1 or len(a) < 2 or not np.isfinite(a).all()
            or np.any(np.diff(a) <= 0.) or a[0] != endpoints[0] or a[-1] != endpoints[1]):
        raise ValueError(f'{name} must strictly increase between exact endpoints {endpoints}')
    return _readonly(a)


def _real_array(value, shape, name, *, scalar=False):
    if np.iscomplexobj(value):
        raise ValueError(f'{name} must be real')
    a = np.asarray(value, dtype=float)
    if scalar and a.ndim == 0:
        a = np.full(shape, float(a))
    if a.shape != shape or not np.isfinite(a).all():
        raise ValueError(f'{name} must be finite with shape {shape}')
    return a


@dataclass(frozen=True, eq=False)
class WallGrid:
    """Immutable piecewise-constant wall cells, ordered (y, x).

    ``x`` and ``y`` are cell-boundary axes over WALL_BOX.  Each cell is either
    exactly ambient (epsinf=1, ell=0), or has epsinf in [1.5,6] and ell in
    [0,.25].  The derivative APIs hold this occupancy fixed; ambient columns
    are zero.  No interpolation or macro expansion changes the public cells.
    """
    x: object
    y: object
    epsinf: object
    ell: object

    def __post_init__(self):
        x = _axis(self.x, WALL_BOX[:2], 'wall x')
        y = _axis(self.y, WALL_BOX[2:], 'wall y')
        shape = (len(y)-1, len(x)-1)
        eps = _real_array(self.epsinf, shape, 'epsinf')
        loss = _real_array(self.ell, shape, 'ell')
        air = eps == 1.
        if (np.any((~air) & ((eps < 1.5) | (eps > 6.)))
                or np.any((loss < 0.) | (loss > .25)) or np.any(air & (loss != 0.))):
            raise ValueError('cells require ambient (1,0) or epsinf in [1.5,6], ell in [0,.25]')
        for name, value in (('x', x), ('y', y), ('epsinf', eps), ('ell', loss)):
            object.__setattr__(self, name, _readonly(value))

    @property
    def shape(self):
        return self.epsinf.shape

    @property
    def loss(self):
        return self.ell

    @property
    def support(self):
        return self.epsinf != 1.

    def cache_key(self):
        return (self.x.tobytes(), self.y.tobytes(), self.epsinf.tobytes(), self.ell.tobytes())

    def indices(self, points):
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
            raise ValueError('wall query points must be finite (n,2)')
        ix = np.searchsorted(self.x, points[:, 0], side='right')-1
        iy = np.searchsorted(self.y, points[:, 1], side='right')-1
        # Closed outer box; the final edge belongs to its final cell.  Internal
        # edges use the cell on their positive side, a measure-zero convention.
        inside = ((points[:, 0] >= self.x[0]) & (points[:, 0] <= self.x[-1])
                  & (points[:, 1] >= self.y[0]) & (points[:, 1] <= self.y[-1]))
        return (np.clip(iy, 0, self.shape[0]-1), np.clip(ix, 0, self.shape[1]-1), inside)

    def values(self, points, f):
        f = _positive(f, 'frequency')
        iy, ix, inside = self.indices(points)
        result = np.ones(len(iy), dtype=complex)
        result[inside] = self.epsinf[iy[inside], ix[inside]]+1j*self.ell[iy[inside], ix[inside]]/f
        return result

    wall_values = values

    def metadata(self):
        return dict(x=self.x.tolist(), y=self.y.tolist(), epsinf=self.epsinf.tolist(),
                    ell=self.ell.tolist(), shape=list(self.shape),
                    conductivity_convention='imaginary permittivity = ell/f',
                    ambient='exactly epsinf=1 and ell=0', derivative_support='fixed occupancy')


@dataclass(frozen=True, eq=False)
class Target:
    """Target coordinates only; independent of the wall mesh and wall values.

    The disk or ellipse is the same inscribed 384-gon used for cut integration.
    Ellipses have semiaxes (radius, 2*radius/3), rotated by 0 or pi/2.  Radius is
    the maximum semiaxis, and loss means Im(epsilon_target)=loss/f.
    """
    pose: object
    radius: float = .12
    contrast: float = 3.
    loss: float = .02
    shape: str = 'disk'
    aspect: float = ELLIPSE_ASPECT
    orientation: float = 0.

    def __post_init__(self):
        pose = _real_array(self.pose, (2,), 'target pose')
        if not -1. <= pose[0] <= 1. or not 1.6 <= pose[1] <= 2.6:
            raise ValueError('target pose must lie in [-1,1] x [1.6,2.6]')
        radius = _positive(self.radius, 'target maximum semiaxis')
        if not .08 <= radius <= .16:
            raise ValueError('target maximum semiaxis must be in [.08,.16]')
        if self.shape not in ('disk', 'ellipse'):
            raise ValueError("target shape must be 'disk' or 'ellipse'")
        if self.shape == 'ellipse' and self.aspect != ELLIPSE_ASPECT:
            raise ValueError('R2 ellipse aspect must be exactly 2/3')
        if self.orientation not in (0., np.pi/2):
            raise ValueError('R2 orientation must be 0 or pi/2')
        if (not np.isfinite(self.contrast) or not 0. <= self.contrast <= 8.
                or not np.isfinite(self.loss) or not 0. <= self.loss <= .25):
            raise ValueError('target contrast must be in [0,8] and loss in [0,.25]')
        object.__setattr__(self, 'pose', _readonly(pose))
        object.__setattr__(self, 'radius', radius)
        object.__setattr__(self, 'contrast', float(self.contrast))
        object.__setattr__(self, 'loss', float(self.loss))
        object.__setattr__(self, 'aspect', 1. if self.shape == 'disk' else ELLIPSE_ASPECT)
        object.__setattr__(self, 'orientation', 0. if self.shape == 'disk' else float(self.orientation))

    @property
    def target_loss(self):
        return self.loss

    def cache_key(self):
        return (tuple(self.pose), self.radius, self.contrast, self.loss,
                self.shape, self.aspect, self.orientation)

    def transform(self):
        rotation = (np.eye(2) if self.orientation == 0. else np.array(((0., -1.), (1., 0.))))
        return rotation@np.diag((1., self.aspect))

    def target_vertices(self, n=TARGET_POLYGON_N):
        if isinstance(n, (bool, np.bool_)) or int(n) != n or n < 3:
            raise ValueError('polygon vertex count must be integer >=3')
        angle = 2*np.pi*np.arange(int(n))/int(n)
        return self.pose+self.radius*(np.column_stack((np.cos(angle), np.sin(angle)))@self.transform().T)

    def metadata(self):
        return dict(pose=self.pose.tolist(), radius=self.radius, contrast=self.contrast,
                    loss=self.loss, shape=self.shape, aspect=self.aspect,
                    orientation=self.orientation, polygon_vertices=TARGET_POLYGON_N,
                    conductivity_convention='imaginary permittivity = loss/f')


def _refined_axis(axis, refinement):
    return np.r_[np.concatenate([np.linspace(a, b, refinement+1)[:-1]
                                for a, b in zip(axis[:-1], axis[1:])]), axis[-1]]


def _mesh_settings(wall_x, wall_y, wall_refine, h, targetstep, pml):
    x = _axis(wall_x, WALL_BOX[:2], 'wall x')
    y = _axis(wall_y, WALL_BOX[2:], 'wall y')
    if isinstance(wall_refine, (bool, np.bool_)) or wall_refine not in (1, 2, 3):
        raise ValueError('wall_refine must be 1, 2, or 3')
    pm = PMLConfig() if pml is None else pml
    if not isinstance(pm, PMLConfig):
        raise TypeError('pml must be the pinned W-SLAM PMLConfig')
    return x, y, int(wall_refine), _positive(h, 'h'), _positive(targetstep, 'targetstep'), pm


def _public_pslg(wall_x, wall_y, *, wall_refine=1, h=.17, targetstep=.125, pml=None):
    """Return atomic constrained segments, with no target or material inputs."""
    x, y, refine, h, step, pm = _mesh_settings(wall_x, wall_y, wall_refine, h, targetstep, pml)
    vertices, segments, lookup, edges = [], [], {}, set()
    def vertex(point):
        key = tuple(float(v) for v in point)
        if key not in lookup:
            lookup[key] = len(vertices)
            vertices.append(key)
        return lookup[key]
    def edge(a, b):
        first, second = vertex(a), vertex(b)
        key = tuple(sorted((first, second)))
        if first != second and key not in edges:
            edges.add(key)
            segments.append((first, second))
    def rectangle(box):
        xl, xr, yb, yt = box
        corners = ((xl, yb), (xr, yb), (xr, yt), (xl, yt))
        for a, b in zip(corners, corners[1:]+corners[:1]):
            edge(a, b)
    def grid(gx, gy):
        for yi in gy:
            for a, b in zip(gx[:-1], gx[1:]):
                edge((a, yi), (b, yi))
        for xi in gx:
            for a, b in zip(gy[:-1], gy[1:]):
                edge((xi, a), (xi, b))
    for box in (pm.outer_box, pm.physical_box, MATERIAL_BOX):
        rectangle(box)
    grid(_refined_axis(x, refine), _refined_axis(y, refine))
    grid(_uniform_axis(TARGET_BOX[0], TARGET_BOX[1], step),
         _uniform_axis(TARGET_BOX[2], TARGET_BOX[3], step))
    key = (MESH_VERSION, tuple(x), tuple(y), refine, h, step, TARGET_BOX, WALL_BOX,
           tuple(pm.physical_box), tuple(pm.thickness), tuple(pm.strength), pm.power)
    return dict(vertices=np.asarray(vertices, float), segments=np.asarray(segments, np.int32)), key


def _check_wall_interfaces(mesh, x, y):
    """Reject any triangle crossing the interior of a wall-interface segment."""
    tri = mesh.p[:, mesh.t]
    tolerance = 1000*np.finfo(float).eps*max(1., float(np.abs(mesh.p).max()))
    for axis, cuts, low, high in ((0, x, y[0], y[-1]), (1, y, x[0], x[-1])):
        for cut in cuts:
            crossing = ((tri[axis].min(axis=0) < cut-tolerance)
                        & (tri[axis].max(axis=0) > cut+tolerance))
            selected = tri[:, :, crossing]
            lower = np.full(selected.shape[2], np.inf)
            upper = np.full(selected.shape[2], -np.inf)
            for first, second in ((0, 1), (1, 2), (2, 0)):
                a, b = selected[:, first], selected[:, second]
                denominator = b[axis]-a[axis]
                fraction = np.full(len(denominator), np.nan)
                np.divide(cut-a[axis], denominator, out=fraction, where=denominator != 0.)
                valid = (denominator != 0.) & (fraction >= 0.) & (fraction <= 1.)
                location = a[1-axis]+fraction*(b[1-axis]-a[1-axis])
                lower[valid] = np.minimum(lower[valid], location[valid])
                upper[valid] = np.maximum(upper[valid], location[valid])
            if np.any((upper > low+tolerance) & (lower < high-tolerance)):
                raise ValueError('mesh triangle straddles a public wall interface')


def public_mesh(wall_x, wall_y, *, wall_refine=1, h=.17, targetstep=.125, pml=None):
    """Use the pinned Triangle PSLG recipe with caller-supplied wall axes."""
    import triangle
    x, y, refine, h, step, pm = _mesh_settings(wall_x, wall_y, wall_refine, h, targetstep, pml)
    pslg, key = _public_pslg(x, y, wall_refine=refine, h=h, targetstep=step, pml=pm)
    result = triangle.triangulate(pslg, f'pq28a{h*h/2:.12g}Q')
    if 'triangles' not in result:
        raise RuntimeError('R2 public PSLG triangulation failed')
    mesh = MeshTri(result['vertices'].T, result['triangles'].T)
    _check_wall_interfaces(mesh, x, y)
    mesh.r2_public_key = key
    return mesh


class AffineTargetAssembler(TargetAssembler):
    """Pinned cut-polygon integral with the R2 affine ellipse geometry."""
    def __init__(self, *args, aspect=ELLIPSE_ASPECT, orientation=0., **kwargs):
        # Validate the public geometry without invoking FEM or constructing a target mesh.
        geometry = Target((0., 2.), shape='ellipse', aspect=aspect, orientation=orientation)
        self.aspect, self.orientation = geometry.aspect, geometry.orientation
        self._transform = geometry.transform()
        super().__init__(*args, **kwargs)

    def quadrature(self, pose, radius):
        geometry = Target(pose, radius, shape='ellipse', aspect=self.aspect, orientation=self.orientation)
        center, radius = geometry.pose, geometry.radius
        vertices = geometry.target_vertices()
        extent = radius*np.sum(np.abs(self._transform), axis=1)
        near = np.asarray(sorted(self._tree.query_ball_point(center, radius+self._max_radius)), dtype=np.int64)
        if len(near):
            near = near[np.all(self._upper[near] > center-extent, axis=1)
                        & np.all(self._lower[near] < center+extent, axis=1)]
        pieces, touched, count = [], [], 0
        inverse = np.linalg.inv(radius*self._transform)
        incircle = np.cos(np.pi/TARGET_POLYGON_N)
        for local in near:
            tri = _ccw(self._triangles[local])
            polygon = (tri if np.all(np.linalg.norm((tri-center)@inverse.T, axis=1) <= incircle)
                       else _clip_convex(vertices, tri))
            if _area(polygon) <= 0.:
                continue
            element, found = int(self.roi_elements[local]), False
            for j in range(1, len(polygon)-1):
                a, b, c = polygon[0], polygon[j], polygon[j+1]
                determinant = float(_cross(b-a, c-a))
                if determinant <= 0.:
                    continue
                count += len(self.reference_weights)
                if count > self.max_points or count*32 > self.max_work_bytes:
                    raise MemoryError('R2 ellipse cut-polygon quadrature guard exceeded')
                pieces.append((a.copy(), b.copy(), c.copy(), determinant, element))
                found = True
            if found:
                touched.append(element)
        if not pieces:
            raise RuntimeError('empty R2 ellipse/mesh intersection')
        points, weights = np.empty((count, 2)), np.empty(count)
        elements = np.empty(count, dtype=np.int64)
        nq = len(self.reference_weights)
        for j, (a, b, c, determinant, element) in enumerate(pieces):
            sl = slice(j*nq, (j+1)*nq)
            points[sl] = a+self.reference_points[:, :1]*(b-a)+self.reference_points[:, 1:]*(c-a)
            weights[sl] = determinant*self.reference_weights
            elements[sl] = element
        area = .5*TARGET_POLYGON_N*self.aspect*radius**2*np.sin(2*np.pi/TARGET_POLYGON_N)
        return TargetQuadrature(points, weights, elements, np.asarray(touched, dtype=np.int64),
                                float(weights.sum()), float(area), vertices,
                                self.reference_points.copy(), self.reference_weights.copy(),
                                self.degree, len(near))

    def assemble(self, *args, **kwargs):
        target = super().assemble(*args, **kwargs)
        metadata = dict(target.metadata)
        metadata.pop('circle_area', None)
        metadata.pop('polygon_relative_circle_area_deficit', None)
        radius = float(args[1] if len(args) > 1 else kwargs['radius'])
        metadata.update(target_shape='ellipse', aspect=self.aspect, orientation=self.orientation,
                        ellipse_area=np.pi*self.aspect*radius**2,
                        polygon_relative_ellipse_area_deficit=(
                            1-np.sin(2*np.pi/TARGET_POLYGON_N)/(2*np.pi/TARGET_POLYGON_N)))
        return replace(target, metadata=metadata)


class _MeteredFactor:
    """No direct, transpose, or adjoint action can bypass RHS reservation."""
    def __init__(self, factor, matrix, book, phase, freq):
        self._factor, self._matrix = factor, matrix
        self._book, self._phase, self._freq = book, phase, freq
        self.last_defect = None

    def solve(self, rhs, trans='N', *, label='r2 linear RHS', **metadata):
        rhs = np.asarray(rhs, dtype=complex)
        if (trans not in ('N', 'T', 'H') or rhs.ndim not in (1, 2)
                or rhs.shape[0] != self._matrix.shape[0] or not np.isfinite(rhs).all()
                or (rhs.ndim == 2 and rhs.shape[1] == 0)):
            raise ValueError('linear RHS/transpose must be finite and match the free-DOF operator')
        count = 1 if rhs.ndim == 1 else rhs.shape[1]
        self._book.reserve(self._phase, label, rhs=count, factorizations=0,
                           freq=self._freq, transpose=trans, **metadata)
        result = self._factor.solve(np.asfortranarray(rhs), trans=trans)
        matrix = (self._matrix if trans == 'N' else
                  self._matrix.T if trans == 'T' else self._matrix.conj().T)
        self.last_defect = _defect(matrix, result, rhs)
        self._book.event(self._phase, 'r2_checked_linear_action', transpose=trans,
                         solved_rhs=count, backward=self.last_defect, freq=self._freq)
        return result


@dataclass(frozen=True, eq=False)
class State:
    data: np.ndarray
    u: np.ndarray
    matrix: object
    rhs: np.ndarray
    factor: object
    target: object
    wall: WallGrid
    world: Target
    residual: dict
    timings: dict
    metadata: dict
    derivative_audit: dict = field(default_factory=dict)
    _owner_token: object = field(default=None, repr=False)

    @property
    def target_spec(self):
        return self.world


def _column_basis(shape, selected_cells=None, basis=None):
    """Real cell directions, ncell by ncolumn; no scientific prior is added."""
    size = int(np.prod(shape))
    if selected_cells is not None and basis is not None:
        raise ValueError('choose selected_cells or basis, not both')
    if basis is not None:
        a = np.asarray(basis)
        if np.iscomplexobj(a):
            raise ValueError('wall basis must contain real parameter directions')
        a = np.asarray(a, dtype=float)
        if a.ndim == 3 and a.shape[:2] == shape:
            a = a.reshape(size, a.shape[2])
        if a.ndim != 2 or a.shape[0] != size or not np.isfinite(a).all():
            raise ValueError('wall basis must be finite (ncell,ncolumn) or (ny,nx,ncolumn)')
        return a
    if selected_cells is None:
        return np.eye(size)
    cells = np.asarray(selected_cells)
    if cells.size == 0:
        return np.empty((size, 0))
    if not np.issubdtype(cells.dtype, np.integer):
        raise ValueError('selected_cells must be integer flat indices or (y,x) pairs')
    if cells.ndim == 2 and cells.shape[1] == 2:
        if np.any(cells < 0) or np.any(cells >= np.asarray(shape)):
            raise ValueError('selected wall cell out of range')
        cells = cells[:, 0]*shape[1]+cells[:, 1]
    elif cells.ndim != 1:
        raise ValueError('selected_cells must be flat indices or (y,x) pairs')
    if np.any(cells < 0) or np.any(cells >= size) or len(np.unique(cells)) != len(cells):
        raise ValueError('selected wall cells must be unique and in range')
    return np.eye(size)[:, cells]


def _cell_product_kernel(adjoint, total, weights, cells, active, ncell, k2):
    """Source-major receiver products accumulated into actual public cells."""
    accumulated = np.zeros((ncell, total.shape[1], adjoint.shape[1]), complex)
    products = k2*np.einsum('qr,qs,q->qsr', adjoint.conj(), total, weights, optimize=True)
    np.add.at(accumulated, cells[active], products[active])
    return accumulated


class Eval(FixedEvaluator):
    """A single public wall mesh/frequency and an externally metered state slot.

    ``wall_x``/``wall_y`` are required. ``tx`` is the frozen two-source S0 or
    four-source S1 array; receivers are RX16.  The caller owns phase, stacking,
    CPU/RSS limits, observations, optimizer, and all scientific decisions.
    """
    def __init__(self, book, phase='R2', f=1., tx=None, wall_x=None, wall_y=None,
                 wall_refine=1, h=.17, targetstep=.125, *, pml=None,
                 quadrature_order=None, max_dofs=80000, max_factor_bytes=4_000_000_000,
                 target_degree=8, target_chunk_points=4096, target_max_points=250_000,
                 target_max_work_bytes=128_000_000, material_chunk_elements=128):
        if wall_x is None or wall_y is None:
            raise ValueError('public wall_x and wall_y axes are required')
        x, y, refine, h, step, pm = _mesh_settings(wall_x, wall_y, wall_refine, h, targetstep, pml)
        if (isinstance(material_chunk_elements, (bool, np.bool_))
                or int(material_chunk_elements) != material_chunk_elements or material_chunk_elements < 1):
            raise ValueError('material_chunk_elements must be a positive integer')
        tx = S0_TX if tx is None else np.asarray(tx, dtype=float)
        if not (np.array_equal(tx, S0_TX) or np.array_equal(tx, S1_TX)):
            raise ValueError('R2 tx must be the frozen S0_TX or S1_TX array')
        self.external_book, self.book, self.phase = book, _BookAdapter(book), phase
        self.declared_completion, self._owner_token = False, object()
        self.wall_x, self.wall_y, self.wall_refine = x, y, refine
        self.material_chunk_elements = int(material_chunk_elements)
        self.p = FEM(_positive(f, 'frequency'), tx, RX16, self.book, h=h, order=3,
                     phase=phase, mesh=public_mesh(x, y, wall_refine=refine, h=h,
                                                  targetstep=step, pml=pm),
                     pml=pm, quadrature_order=quadrature_order,
                     residual_tolerance=BACKWARD_TOLERANCE, max_dofs=max_dofs,
                     max_factor_bytes=max_factor_bytes)
        self._target_settings = dict(degree=target_degree, chunk_points=target_chunk_points,
                                     max_points=target_max_points, max_work_bytes=target_max_work_bytes)
        self.ass = TargetAssembler(self.p.basis, self.p._free_dofs, **self._target_settings)
        self._ellipse_assemblers = {}
        self.R = self.p._receiver_probe[:, self.p._free_dofs]
        self.direct = _readonly(line_green(self.p.rx, self.p.tx, self.p.k).T)
        self._to_free = np.full(self.p.ndofs, -1, dtype=np.int64)
        self._to_free[self.p._free_dofs] = np.arange(len(self.p._free_dofs))
        # Cell ownership depends only on the public axes, never on materials.
        template = WallGrid(x, y, np.ones((len(y)-1, len(x)-1)),
                            np.zeros((len(y)-1, len(x)-1)))
        iy, ix, self._inside_wall = template.indices(self.p.material_points)
        self._cells = iy*template.shape[1]+ix
        self._key = self._last = None

    def clear_cache(self):
        self._key = self._last = None

    def _wall(self, wall):
        if not isinstance(wall, WallGrid):
            raise TypeError('wall must be a validated WallGrid')
        if not np.array_equal(wall.x, self.wall_x) or not np.array_equal(wall.y, self.wall_y):
            raise ValueError('wall axes must equal this evaluator public mesh axes')
        return wall

    def _state(self, state):
        if not isinstance(state, State) or state._owner_token is not self._owner_token:
            raise ValueError('state belongs to a different R2 evaluator')

    def _assembler(self, target):
        if target.shape == 'disk':
            return self.ass
        key = (target.aspect, target.orientation)
        if key not in self._ellipse_assemblers:
            self._ellipse_assemblers[key] = AffineTargetAssembler(
                self.p.basis, self.p._free_dofs, aspect=target.aspect,
                orientation=target.orientation, **self._target_settings)
        return self._ellipse_assemblers[key]

    def background(self, wall):
        """Return L_W and its weak forcing b_W; no factorization or RHS solve."""
        wall, p = self._wall(wall), self.p
        contrast = wall.values(p.material_points, p.freq)-1.
        mass = asm(_contrast_mass, p.material_basis, contrast=contrast.reshape(p._material_shape))
        matrix = (p._background_matrix-p.k**2*mass[p._free_dofs][:, p._free_dofs]).tocsc()
        rhs = np.column_stack([
            asm(_forcing, p.material_basis,
                load=(p.k**2*contrast*line_green(p.material_points, p.tx[j:j+1], p.k)[:, 0]
                      ).reshape(p._material_shape))[p._free_dofs]
            for j in range(len(p.tx))])
        return matrix, rhs

    def free_operator(self):
        """Air/PML L0 and exact zero weak source forcing; no numerical action."""
        return self.p._background_matrix, np.zeros((len(self.p._free_dofs), len(self.p.tx)), complex)

    def assemble_target(self, target):
        """Return the sparse V_X and direct analytic b_X cut assembly; zero PDE."""
        if not isinstance(target, Target):
            raise TypeError('target must be a validated Target')
        return self._assembler(target).assemble(
            target.pose, target.radius, self.p.k, target.contrast+1j*target.loss/self.p.freq,
            lambda points: line_green(points, self.p.tx, self.p.k), source_count=len(self.p.tx))

    def _factor(self, matrix, label):
        if not np.isfinite(matrix.data).all():
            raise RuntimeError('nonfinite R2 FEM matrix')
        self.book.reserve(self.phase, label, rhs=0, factorizations=1, freq=self.p.freq,
                          declared_completion=self.declared_completion)
        native = splu(matrix, permc_spec='COLAMD')
        storage = (_sparse_bytes(native.L)+_sparse_bytes(native.U)
                   + native.perm_r.nbytes+native.perm_c.nbytes)
        if storage > self.p.max_factor_bytes:
            raise MemoryError('actual R2 factor storage exceeds backend guard')
        matrix = _immutable_sparse(matrix)
        return _MeteredFactor(native, matrix, self.book, self.phase, self.p.freq), matrix, int(storage)

    def factor(self, matrix, label='r2 external LU'):
        """Meter one caller-assembled LU; returned solve() meters every RHS."""
        return self._factor(csc_matrix(matrix), label)

    def evaluate(self, wall, target, label='r2 full world'):
        """Solve (L_W-V_X)u=b_W+b_X; return source-major total receiver data."""
        try:
            wall = self._wall(wall)
            if not isinstance(target, Target):
                raise TypeError('target must be a validated Target')
            key = (wall.cache_key(), target.cache_key())
        except Exception:
            self.clear_cache()
            raise
        if key == self._key and self._last is not None:
            return self._last
        self.clear_cache()
        cpu, started = time.process_time(), time.perf_counter()
        matrix, rhs = self.background(wall)
        assembly = self.assemble_target(target)
        matrix, rhs = (matrix-assembly.Vtarget).tocsc(), rhs+assembly.bX
        factor, matrix, storage = self._factor(matrix, label+' LU')
        u = factor.solve(rhs, label=label+' primal', declared_completion=self.declared_completion)
        data = self.direct+np.asarray(self.R@u).T
        timings = dict(kernel_cpu_seconds=time.process_time()-cpu,
                       wall_seconds=time.perf_counter()-started, factor_bytes=storage)
        metadata = dict(freq=self.p.freq, tx=self.p.tx.tolist(), rx=self.p.rx.tolist(),
                        source_major=True, equation='(L_W-V_X)u=b_W+b_X',
                        wall_grid_shape=list(wall.shape), wall_refine=self.wall_refine,
                        mesh_version=MESH_VERSION, target=target.metadata(),
                        target_forcing='analytic incident integral at independent cut-polygon quadrature',
                        fields='not retained; bounded material element interpolation chunks',
                        scientific_qualification='NOT_RUN', RSS_guard='caller-owned')
        state = State(_readonly(data), _readonly(u), matrix, _readonly(rhs), factor,
                      _immutable_target(assembly), wall, target, dict(factor.last_defect),
                      timings, metadata, _owner_token=self._owner_token)
        self.book.event(self.phase, 'r2_physics_state', label_of_world=label, freq=self.p.freq,
                         backward=state.residual, ndofs=self.p.ndofs, factor_bytes=storage,
                         kernel_cpu_seconds=timings['kernel_cpu_seconds'])
        self._key, self._last = key, state
        return state

    def _material_chunks(self, *free_values):
        """Interpolate existing P3 basis arrays, keeping only a bounded chunk."""
        arrays = [np.asarray(v)[:, None] if np.asarray(v).ndim == 1 else np.asarray(v)
                  for v in free_values]
        p, basis = self.p, self.p.material_basis
        for a in arrays:
            if a.ndim != 2 or a.shape[0] != len(p._free_dofs):
                raise ValueError('material interpolation requires free-DOF columns')
        nq = p._material_shape[1]
        for begin in range(0, p._material_shape[0], self.material_chunk_elements):
            end = min(begin+self.material_chunk_elements, p._material_shape[0])
            indices = self._to_free[basis.element_dofs[:, begin:end]]
            valid = indices >= 0
            outputs = []
            for a in arrays:
                coefficients = np.zeros((*indices.shape, a.shape[1]), complex)
                coefficients[valid] = a[indices[valid]]
                values = np.zeros((end-begin, nq, a.shape[1]), complex)
                for j, phi in enumerate(basis.basis):
                    values += np.asarray(phi[0])[begin:end, :, None]*coefficients[j, :, None, :]
                outputs.append(values.reshape(-1, a.shape[1]))
            yield slice(begin*nq, end*nq), outputs

    def _direction(self, wall, eps_direction, ell_direction=0.):
        eps = _real_array(eps_direction, wall.shape, 'eps direction', scalar=True)
        loss = _real_array(ell_direction, wall.shape, 'ell direction', scalar=True)
        active = self._inside_wall & wall.support.ravel()[self._cells]
        coefficient = np.zeros(len(self._cells), complex)
        coefficient[active] = (eps.ravel()[self._cells[active]]
                               + 1j*loss.ravel()[self._cells[active]]/self.p.freq)
        return coefficient

    def wall_mass_direction(self, wall, eps_direction, ell_direction=0.):
        """Positive V'=k² int(delta_epsilon phi phi), not A'=-V'; zero PDE."""
        wall, p = self._wall(wall), self.p
        coefficient = self._direction(wall, eps_direction, ell_direction)
        mass = asm(_contrast_mass, p.material_basis, contrast=coefficient.reshape(p._material_shape))
        return (p.k**2*mass[p._free_dofs][:, p._free_dofs]).tocsc()

    def iter_wall_mass_derivatives(self, wall, *, selected_cells=None, basis=None, include_loss=True):
        """Yield each epsinf then ell mass column without storing every matrix."""
        wall = self._wall(wall)
        columns = _column_basis(wall.shape, selected_cells, basis)
        for kind in (('epsinf', 'ell') if include_loss else ('epsinf',)):
            for j in range(columns.shape[1]):
                direction = columns[:, j].reshape(wall.shape)
                yield (kind, j, self.wall_mass_direction(
                    wall, direction if kind == 'epsinf' else 0., direction if kind == 'ell' else 0.))

    def wall_differential_rhs(self, state, eps_direction, ell_direction=0.):
        """V' u+b_W': includes the incident forcing derivative; zero PDE."""
        self._state(state)
        p = self.p
        coefficient = self._direction(state.wall, eps_direction, ell_direction)
        result = self.wall_mass_direction(state.wall, eps_direction, ell_direction)@state.u
        for j in range(len(p.tx)):
            result[:, j] += asm(
                _forcing, p.material_basis,
                load=(p.k**2*coefficient*line_green(p.material_points, p.tx[j:j+1], p.k)[:, 0]
                      ).reshape(p._material_shape))[p._free_dofs]
        return result

    def wall_jvp(self, state, eps_direction, ell_direction=0., *, label='r2 wall JVP'):
        rhs = self.wall_differential_rhs(state, eps_direction, ell_direction)
        solved = state.factor.solve(rhs, label=label, direction_shared_across_sources=True)
        return np.asarray(self.R@solved).T

    wall_directional_jvp = wall_jvp

    def _receiver_adjoint(self, state, label):
        return state.factor.solve(self.R.conj().T.toarray(), trans='H', label=label,
                                  shared_across_material_cells=True)

    def _wall_jacobian(self, state, adjoint, columns, include_loss):
        ncell = int(np.prod(state.wall.shape))
        accumulated = np.zeros((ncell, len(self.p.tx), len(self.p.rx)), complex)
        active = self._inside_wall & state.wall.support.ravel()[self._cells]
        for sl, (scattered, adj) in self._material_chunks(state.u, adjoint):
            total = scattered+line_green(self.p.material_points[sl], self.p.tx, self.p.k)
            accumulated += _cell_product_kernel(adj, total, self.p.material_weights[sl],
                                                self._cells[sl], active[sl], ncell, self.p.k**2)
        eps = (accumulated.transpose(1, 2, 0).reshape(state.data.size, ncell)@columns)
        return np.column_stack((eps, 1j/self.p.freq*eps)) if include_loss else eps

    def wall_jacobian(self, state, *, selected_cells=None, basis=None, include_loss=True):
        """All chosen cell/basis columns using exactly 16 receiver-adjoint RHS.

        Rows flatten (source,receiver); columns are chosen epsinf directions,
        then the same ell directions.  Either 18/144 cells or arbitrary real
        cell-basis columns are supported on the evaluator's public grid.
        """
        self._state(state)
        columns = _column_basis(state.wall.shape, selected_cells, basis)
        if not columns.shape[1]:
            return np.empty((state.data.size, 0), complex)
        adjoint = self._receiver_adjoint(state, 'r2 shared receiver adjoint for wall cells')
        result = self._wall_jacobian(state, adjoint, columns, include_loss)
        state.derivative_audit['wall_jacobian'] = dict(
            columns=int(columns.shape[1]), include_loss=include_loss,
            receiver_adjoint_rhs=len(self.p.rx), source_shared=True,
            ordering='source-major rows; supplied epsinf then ell directions',
            backward=dict(state.factor.last_defect))
        return result

    def _target_names(self, names):
        names = tuple('loss' if name == 'target_loss' else name for name in names)
        if any(name not in TARGET_NAMES for name in names) or len(set(names)) != len(names):
            raise ValueError(f'target differential names must be unique members of {TARGET_NAMES}')
        return names

    def target_differentials(self, state, *, names=TARGET_NAMES):
        """Source-resolved V_X' u+b_X' for five real coordinates; zero PDE.

        Volume material derivatives use the identical cut quadrature.  Shape
        derivatives use the same 384-gon boundary and inherited degree-4 edge
        rule.  Orientation/aspect are fixed hypotheses, not derivative axes.
        """
        self._state(state)
        names, p, target = self._target_names(names), self.p, state.world
        outputs = {name: np.zeros_like(state.u) for name in names}
        quadrature = state.target.quadrature
        material = set(names) & {'contrast', 'loss'}
        if material:
            volume = np.zeros_like(state.u)
            for start in range(0, len(quadrature.weights), self.ass.chunk_points):
                sl = slice(start, min(start+self.ass.chunk_points, len(quadrature.weights)))
                probe = self.ass._known_element_probes(
                    quadrature.points[sl], quadrature.elements[sl], p._free_dofs)
                total = line_green(quadrature.points[sl], p.tx, p.k)+probe@state.u
                volume += p.k**2*np.asarray(probe.T@(quadrature.weights[sl, None]*total))
            if 'contrast' in outputs:
                outputs['contrast'] = volume
            if 'loss' in outputs:
                outputs['loss'] = 1j/p.freq*volume
        shape_names = set(names) & {'px', 'py', 'radius'}
        if shape_names:
            vertices = target.target_vertices()
            edge = np.roll(vertices, -1, axis=0)-vertices
            length = np.linalg.norm(edge, axis=1)
            nodes, weights = np.polynomial.legendre.leggauss(4)
            points = (vertices[:, None, :]+edge[:, None, :]*(nodes[None, :, None]+1)/2).reshape(-1, 2)
            normal = np.repeat(np.column_stack((edge[:, 1], -edge[:, 0]))/length[:, None], 4, axis=0)
            ds = (length[:, None]*weights[None, :]/2).ravel()
            velocities = dict(px=normal[:, 0], py=normal[:, 1],
                              radius=np.einsum('qi,qi->q', normal, (points-target.pose)/target.radius))
            probe = self.probe(points)
            total = line_green(points, p.tx, p.k)+probe@state.u
            scale = p.k**2*(target.contrast+1j*target.loss/p.freq)
            for name in shape_names:
                outputs[name] = scale*np.asarray(probe.T@((ds*velocities[name])[:, None]*total))
        return outputs

    def _target_jacobian(self, state, adjoint, names):
        differentials = self.target_differentials(state, names=names)
        return np.column_stack([(adjoint.conj().T@differentials[name]).T.ravel()
                                for name in names]) if names else np.empty((state.data.size, 0), complex)

    def target_jacobian(self, state, *, names=TARGET_NAMES):
        """Five target columns from a shared 16-RHS receiver adjoint."""
        self._state(state)
        names = self._target_names(names)
        if not names:
            return np.empty((state.data.size, 0), complex)
        adjoint = self._receiver_adjoint(state, 'r2 shared receiver adjoint for target')
        result = self._target_jacobian(state, adjoint, names)
        state.derivative_audit['target_jacobian'] = dict(
            names=list(names), receiver_adjoint_rhs=len(self.p.rx),
            backward=dict(state.factor.last_defect), target_shape=state.world.shape)
        return result

    def jacobians(self, state, *, selected_cells=None, basis=None, include_loss=True,
                  target_names=TARGET_NAMES):
        """Joint wall/target matrices sharing one 16-RHS receiver adjoint."""
        self._state(state)
        columns = _column_basis(state.wall.shape, selected_cells, basis)
        names = self._target_names(target_names)
        if not columns.shape[1] and not names:
            empty = np.empty((state.data.size, 0), complex)
            return empty, empty.copy()
        adjoint = self._receiver_adjoint(state, 'r2 shared receiver adjoint for joint Jacobian')
        A = self._wall_jacobian(state, adjoint, columns, include_loss)
        B = self._target_jacobian(state, adjoint, names)
        state.derivative_audit['joint_jacobian'] = dict(
            wall_columns=int(columns.shape[1]), target_names=list(names),
            receiver_adjoint_rhs=len(self.p.rx), backward=dict(state.factor.last_defect))
        return A, B

    def target_jvp(self, state, direction, *, names=TARGET_NAMES, label='r2 target JVP'):
        self._state(state)
        names = self._target_names(names)
        direction = _real_array(direction, (len(names),), 'target direction')
        differentials = self.target_differentials(state, names=names)
        rhs = np.zeros_like(state.u)
        for name, value in zip(names, direction):
            rhs += value*differentials[name]
        solved = state.factor.solve(rhs, label=label, direction_shared_across_sources=True)
        return np.asarray(self.R@solved).T

    def objective_gradient(self, state, complex_residual_weights, *, selected_cells=None,
                           basis=None, include_loss=True, target_names=TARGET_NAMES):
        """Real scalar gradient from exactly ns weighted-adjoint RHS.

        The caller supplies weights with delta_phi=Re(vdot(weights,delta_data)).
        No whitening, prior, nuisance gain, or optimization step is added here.
        """
        self._state(state)
        weights = np.asarray(complex_residual_weights, dtype=complex)
        if weights.shape != state.data.shape or not np.isfinite(weights).all():
            raise ValueError('objective weights must be finite complex (source,receiver)')
        columns = _column_basis(state.wall.shape, selected_cells, basis)
        names = self._target_names(target_names)
        rhs = np.asarray(self.R.conj().T@weights.T)
        adjoint = state.factor.solve(rhs, trans='H', label='r2 weighted objective adjoint',
                                     weights_convention='delta_phi=Re(vdot(weights,delta_data))')
        accumulated = np.zeros(int(np.prod(state.wall.shape)), complex)
        active = self._inside_wall & state.wall.support.ravel()[self._cells]
        for sl, (scattered, adj) in self._material_chunks(state.u, adjoint):
            total = scattered+line_green(self.p.material_points[sl], self.p.tx, self.p.k)
            products = self.p.k**2*self.p.material_weights[sl]*np.einsum(
                'qs,qs->q', adj.conj(), total, optimize=True)
            mask = active[sl]
            np.add.at(accumulated, self._cells[sl][mask], products[mask])
        eps, loss = accumulated.real, (1j/self.p.freq*accumulated).real
        wall = np.r_[eps@columns, loss@columns] if include_loss else eps@columns
        differentials = self.target_differentials(state, names=names)
        target_gradient = np.asarray([np.vdot(adjoint, differentials[name]).real for name in names])
        audit = dict(weighted_adjoint_rhs=len(self.p.tx), target_names=list(names),
                     backward=dict(state.factor.last_defect), source_shared=True,
                     weights_convention='delta_phi=Re(vdot(weights,delta_data))',
                     shape_derivative='same384gon boundary; degree4 Gaussian edge rule')
        state.derivative_audit['objective_gradient'] = audit
        return dict(wall=wall, epsinf=eps.reshape(state.wall.shape), ell=loss.reshape(state.wall.shape),
                    target=target_gradient, target_names=names, audit=audit)

    def wall_objective_gradient(self, state, complex_residual_weights, *, selected_cells=None,
                                basis=None, include_loss=True):
        return self.objective_gradient(
            state, complex_residual_weights, selected_cells=selected_cells, basis=basis,
            include_loss=include_loss, target_names=())['wall']


Evaluator = Eval
FixedWorld = Target
PhysicsState = State
