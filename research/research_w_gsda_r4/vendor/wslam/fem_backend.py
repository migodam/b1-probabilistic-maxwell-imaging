"""Independent high-order scalar-TM Helmholtz FEM with coordinate-stretch PML.

This module does not import the P0 VIE.  It uses exp(-i omega t), k=2*pi*f,
epsilon=epsilon_infinity+i*conductivity/f, and the unit outgoing line source
g(x,s)=i*H_0^(1)(k*|x-s|)/4.  The material is exactly reference medium outside
D=[-2,2] x [.25,3].  Only the scattered field is discretized:

    (-Delta - k**2 epsilon) u_s = k**2 (epsilon-1) u_inc.

On the Cartesian PML, s_j=1+i*alpha_j*(distance_j/thickness_j)**power,
A=diag(s_y/s_x, s_x/s_y), and B=s_x*s_y.  The weak operator is
integral(A grad(u_s).grad(v) - k**2*B*epsilon*u_s*v).  Positive imaginary
stretch damps outgoing exp(+i*k*r) waves for this time convention.  A zero
Dirichlet condition terminates the outer PML boundary.

P2/P3 continuous Lagrange triangles conform to the rectangular material box
and PML interfaces.  Curved/internal material boundaries are sampled at actual
triangle quadrature points, never represented by center/area-average VIE
cells.  Mesh, quadrature, geometry resolution, and PML convergence must be
qualified by the caller.  No constructor default certifies any of them.

All numerical LU and RHS actions reserve the supplied budget before action.
CPU charging belongs to the caller's job wrapper; this module emits only
uncharged informational timing events.  It has no budget or experiment loop.
"""

from __future__ import annotations

# Adapted from pinned continuation 483d260c6b5b111db9fb814b096685fccb362448.
# New campaign domain, physical box and conductivity range only; see ADAPTER.md.

from dataclasses import dataclass, field as dataclass_field
import time
from typing import Callable

import numpy as np
from scipy.sparse.linalg import splu
from scipy.special import hankel1
from skfem import (Basis, BilinearForm, ElementTriP2, ElementTriP3,
                   LinearForm, MeshTri, asm)


MATERIAL_BOX = (-2.0, 2.0, 0.25, 3.0)


def _four(value, name, *, positive=True):
    """Normalize side parameters in (left, right, bottom, top) order."""
    a = np.asarray(value, dtype=float)
    if a.ndim == 0:
        a = np.repeat(a, 4)
    if a.shape != (4,) or not np.isfinite(a).all():
        raise ValueError(f'{name} must be finite scalar or four side values')
    if positive and np.any(a <= 0):
        raise ValueError(f'{name} must be strictly positive')
    return tuple(float(x) for x in a)


def _coordinates(value, name):
    a = np.asarray(value, dtype=float)
    if a.ndim != 2 or a.shape[1] != 2 or len(a) == 0 or not np.isfinite(a).all():
        raise ValueError(f'{name} requires finite nonempty (n,2) coordinates')
    return a.copy()


def _inside(points, box, *, closed=False):
    x, y = np.asarray(points).T
    lo_x, hi_x, lo_y, hi_y = box
    if closed:
        return ((x >= lo_x) & (x <= hi_x) & (y >= lo_y) & (y <= hi_y))
    return ((x > lo_x) & (x < hi_x) & (y > lo_y) & (y < hi_y))


@dataclass(frozen=True)
class PMLConfig:
    """Explicit PML parameters; these are starting settings, not qualification.

    ``physical_box`` is (xmin,xmax,ymin,ymax) and must strictly enclose the
    material box and all point antennas.  Side parameter tuples are ordered
    (left,right,bottom,top).  ``strength`` is the dimensionless maximum
    imaginary stretch, not a conductivity or an automatic target tolerance.
    """

    physical_box: tuple = (-2.75, 2.75, -0.8, 3.5)
    thickness: object = 1.0
    strength: object = 4.0
    power: int = 3

    def __post_init__(self):
        box = np.asarray(self.physical_box, dtype=float)
        if (box.shape != (4,) or not np.isfinite(box).all()
                or not (box[0] < -2.0 and box[1] > 2.0
                        and box[2] < MATERIAL_BOX[2] and box[3] > MATERIAL_BOX[3])):
            raise ValueError('physical_box must strictly enclose [-2,2]x[.25,3]')
        if isinstance(self.power, bool) or int(self.power) != self.power or self.power < 1:
            raise ValueError('PML power must be a positive integer')
        object.__setattr__(self, 'physical_box', tuple(float(x) for x in box))
        object.__setattr__(self, 'thickness', _four(self.thickness, 'thickness'))
        object.__setattr__(self, 'strength', _four(self.strength, 'strength'))
        object.__setattr__(self, 'power', int(self.power))

    @property
    def outer_box(self):
        xl, xr, yb, yt = self.physical_box
        tl, tr, tb, tt = self.thickness
        return (xl-tl, xr+tr, yb-tb, yt+tt)

    def stretch(self, points):
        """Return s_x,s_y at real points inside the computational rectangle."""
        p = _coordinates(points, 'PML evaluation')
        if not _inside(p, self.outer_box, closed=True).all():
            raise ValueError('PML evaluation lies outside computational box')
        x, y = p.T
        xl, xr, yb, yt = self.physical_box
        tl, tr, tb, tt = self.thickness
        al, ar, ab, at = self.strength
        sx = (1.0 + 1j*(al*(np.maximum(xl-x, 0.0)/tl)**self.power
                         + ar*(np.maximum(x-xr, 0.0)/tr)**self.power))
        sy = (1.0 + 1j*(ab*(np.maximum(yb-y, 0.0)/tb)**self.power
                         + at*(np.maximum(y-yt, 0.0)/tt)**self.power))
        return sx, sy

    def as_dict(self):
        return dict(physical_box=list(self.physical_box),
                    outer_box=list(self.outer_box), thickness=list(self.thickness),
                    strength=list(self.strength), power=self.power,
                    stretch_sign='+i for exp(-i omega t)',
                    termination='homogeneous Dirichlet on outer boundary')


def line_green(points, sources, k):
    """Analytic unit-strength outgoing line Green, shape (point,source)."""
    p = _coordinates(points, 'observation')
    s = _coordinates(sources, 'source')
    distance = np.linalg.norm(p[:, None, :]-s[None, :, :], axis=-1)
    if np.any(distance == 0):
        raise ValueError('Coincident line source and observation is singular')
    return 0.25j*hankel1(0, float(k)*distance)


@BilinearForm(dtype=np.complex128)
def _background(u, v, w):
    return (w.ax*u.grad[0]*v.grad[0] + w.ay*u.grad[1]*v.grad[1]
            - w.k2*w.volume*u*v)


@BilinearForm(dtype=np.complex128)
def _contrast_mass(u, v, w):
    return w.contrast*u*v


@LinearForm(dtype=np.complex128)
def _forcing(v, w):
    return w.load*v


@dataclass
class State:
    """One solved material/frequency and illumination set.

    Receiver arrays are (source,receiver).  ``epsilon`` is sampled at flattened
    material-element quadrature points; ``incident``, ``field`` (total), and
    ``current=(epsilon-1)*field`` are (material_quadrature,source), in the same
    order as ``material_points`` and ``material_weights``.  ``scattered_dofs``
    are (global_FEM_dof,source), including the zero outer boundary DOFs.

    For line sources ``data=direct+scattered``.  For plane-wave qualification
    ``direct=0`` and ``data=scattered``; the analytic plane wave at the receivers
    remains available separately in ``incident_at_receivers``.  This convention
    matches a scattering-data cylinder comparison without subtracting large
    analytic values numerically.
    """

    epsilon: np.ndarray
    material_points: np.ndarray
    material_weights: np.ndarray
    incident: np.ndarray
    field: np.ndarray
    current: np.ndarray
    scattered_dofs: np.ndarray
    data: np.ndarray
    scattered: np.ndarray
    direct: np.ndarray
    incident_at_receivers: np.ndarray
    factor: object
    residual: float
    residuals: np.ndarray
    timings: dict
    illumination: dict
    metadata: dict
    receiver_mode: str = 'probe'
    probe_data: object = None
    probe_scattered: object = None
    source_major: bool = True
    jacobian: object = None
    _backend_token: object = dataclass_field(default=None, repr=False)
    _matrix: object = dataclass_field(default=None, repr=False)

    @property
    def E(self):
        return self.field


def _axis_nodes(breaks, maximum_step):
    b = sorted(set(float(x) for x in breaks))
    pieces = []
    for left, right in zip(b[:-1], b[1:]):
        count = max(1, int(np.ceil((right-left)/maximum_step)))
        pieces.append(np.linspace(left, right, count+1)[:-1])
    return np.r_[np.concatenate(pieces), b[-1]]


def _validate_mesh(mesh, outer_box):
    """Validate rectangular exterior and conformity to the four D boundary cuts.

    Only roundoff-sized coordinate tolerance is used for mesh topology checks;
    it never modifies a point, a coefficient, a source, or the PDE matrix.
    """
    if not isinstance(mesh, MeshTri):
        raise TypeError('mesh must be a scikit-fem MeshTri')
    if mesh.p.shape[0] != 2 or mesh.t.shape[0] != 3 or not np.isfinite(mesh.p).all():
        raise ValueError('mesh must contain finite two-dimensional triangles')
    scale = max(1.0, *np.abs(outer_box))
    tolerance = 1000*np.finfo(float).eps*scale
    actual = (mesh.p[0].min(), mesh.p[0].max(), mesh.p[1].min(), mesh.p[1].max())
    if not np.allclose(actual, outer_box, rtol=0.0, atol=tolerance):
        raise ValueError('mesh extent must equal PMLConfig.outer_box')
    boundary = mesh.p[:, mesh.facets[:, mesh.boundary_facets()]]
    ol, or_, ob, ot = outer_box
    on_side = (np.all(abs(boundary[0]-ol) <= tolerance, axis=0)
               | np.all(abs(boundary[0]-or_) <= tolerance, axis=0)
               | np.all(abs(boundary[1]-ob) <= tolerance, axis=0)
               | np.all(abs(boundary[1]-ot) <= tolerance, axis=0))
    if not np.all(on_side):
        raise ValueError('mesh must cover whole rectangle without internal holes')
    vertices = mesh.p[:, mesh.t]
    for axis, value, other_low, other_high in ((0, MATERIAL_BOX[0], MATERIAL_BOX[2], MATERIAL_BOX[3]),
                                                (0, MATERIAL_BOX[1], MATERIAL_BOX[2], MATERIAL_BOX[3]),
                                                (1, MATERIAL_BOX[2], MATERIAL_BOX[0], MATERIAL_BOX[1]),
                                                (1, MATERIAL_BOX[3], MATERIAL_BOX[0], MATERIAL_BOX[1])):
        crosses = ((vertices[axis].min(axis=0) < value-tolerance)
                   & (vertices[axis].max(axis=0) > value+tolerance))
        selected = vertices[:, :, crosses]
        if selected.shape[2] == 0:
            continue
        lower = np.full(selected.shape[2], np.inf)
        upper = np.full(selected.shape[2], -np.inf)
        for first, second in ((0, 1), (1, 2), (2, 0)):
            a, b = selected[:, first], selected[:, second]
            denominator = b[axis]-a[axis]
            nonzero = denominator != 0.0
            fraction = np.full(len(denominator), np.nan)
            fraction[nonzero] = (value-a[axis, nonzero])/denominator[nonzero]
            valid = nonzero & (fraction >= 0.0) & (fraction <= 1.0)
            location = a[1-axis]+fraction*(b[1-axis]-a[1-axis])
            lower[valid] = np.minimum(lower[valid], location[valid])
            upper[valid] = np.maximum(upper[valid], location[valid])
        if np.any((upper > other_low+tolerance) & (lower < other_high-tolerance)):
            raise ValueError('mesh triangles straddle the rectangular material boundary')
    edges = mesh.p[:, mesh.facets[0]]-mesh.p[:, mesh.facets[1]]
    return float(np.linalg.norm(edges, axis=0).max())


def _factor_memory_estimate(ndofs):
    """Exposed conservative 2D fill heuristic, not a process-RSS guarantee.

    Allow 64*N*ceil(log2(N)) complex factor entries at 24 bytes per entry
    (complex128 value plus 64-bit index), including a conservative fill margin.
    SuperLU's opaque workspaces and existing Python objects are not measured.
    """
    return int(64*ndofs*max(1, int(np.ceil(np.log2(max(2, ndofs)))))*24)


def _sparse_bytes(matrix):
    return int(matrix.data.nbytes+matrix.indices.nbytes+matrix.indptr.nbytes)


class FEM:
    """P2/P3 scattered-field Helmholtz backend with exact single-factor reuse.

    ``epsilon`` in :meth:`solve` must be callable: finite real coordinates
    (n,2) -> complex (n,), or a scalar broadcast to those n points.  It is called
    only within D, directly at triangle quadrature points.  Permittivity is
    exactly one elsewhere.  The frozen physical range is enforced here:
    1<=Re(epsilon)<=9 and 0<=f*Im(epsilon)<=0.25.

    ``h`` bounds Cartesian mesh intervals; triangle diameter can be sqrt(2)*h.
    It is never silently selected from a frequency or an error target.
    ``quadrature_order`` is the triangle integration polynomial order.  The
    default is max(8,2*order+2); a discontinuous curved coefficient still needs
    independent mesh and quadrature checks.

    Repeating ``solve`` with exactly equal quadrature coefficient values reuses
    the same sparse LU, even for a different plane-wave source set.  Otherwise
    one new factorization is reserved and performed.  There is no tolerance,
    coefficient rounding, jitter, least-squares, or pseudoinverse fallback.
    """

    def __init__(self, freq, tx, rx, book, *, h, order=2,
                 quadrature_order=None, pml=None, phase='validation',
                 source_normalization=1.0, residual_tolerance=1e-9,
                 mesh=None, max_dofs=60000, max_factor_bytes=1_500_000_000,
                 receiver_mode='probe'):
        started = time.perf_counter()
        if not np.isfinite(freq) or freq <= 0:
            raise ValueError('freq=f/f0 must be positive and finite')
        if not np.isfinite(h) or h <= 0:
            raise ValueError('h must be positive and finite')
        if isinstance(order, bool) or order not in (2, 3):
            raise ValueError('Only conforming triangular P2 and P3 are supported')
        if book is None or not callable(getattr(book, 'reserve', None)):
            raise TypeError('An external physical-budget book.reserve is required')
        if not np.isfinite(residual_tolerance) or residual_tolerance <= 0:
            raise ValueError('residual_tolerance must be positive and finite')
        if receiver_mode not in ('probe', 'volume'):
            raise ValueError("receiver_mode must be 'probe' or 'volume'")
        self.receiver_mode = receiver_mode
        self.freq = float(freq)
        self.k = 2.0*np.pi*self.freq
        self.h = float(h)
        self.order = int(order)
        self.quadrature_order = (max(8, 2*self.order+2) if quadrature_order is None
                                 else quadrature_order)
        if (isinstance(self.quadrature_order, bool)
                or int(self.quadrature_order) != self.quadrature_order
                or self.quadrature_order < 2*self.order):
            raise ValueError('quadrature_order must be integer >=2*element order')
        self.quadrature_order = int(self.quadrature_order)
        self.pml = PMLConfig() if pml is None else pml
        if not isinstance(self.pml, PMLConfig):
            raise TypeError('pml must be a PMLConfig')
        self.book = book
        self.phase = phase
        if isinstance(max_dofs, bool) or int(max_dofs) != max_dofs or max_dofs <= 0:
            raise ValueError('max_dofs must be a positive integer')
        if not np.isfinite(max_factor_bytes) or max_factor_bytes <= 0:
            raise ValueError('max_factor_bytes must be positive and finite')
        self.max_dofs = int(max_dofs)
        self.max_factor_bytes = int(max_factor_bytes)
        self.residual_tolerance = float(residual_tolerance)
        self.source_normalization = complex(source_normalization)
        if not np.isfinite(self.source_normalization):
            raise ValueError('source_normalization must be finite')
        self.tx = _coordinates(tx, 'tx')
        self.rx = _coordinates(rx, 'rx')
        for antennas, name in ((self.tx, 'tx'), (self.rx, 'rx')):
            if _inside(antennas, MATERIAL_BOX, closed=True).any():
                raise ValueError(f'{name} must be outside closed material domain')
            if not _inside(antennas, self.pml.physical_box).all():
                raise ValueError(f'{name} must be strictly inside physical_box')
        if np.any(np.linalg.norm(self.rx[:, None, :]-self.tx[None, :, :], axis=-1) == 0):
            raise ValueError('Coincident tx/rx has a singular direct line field')

        xl, xr, yb, yt = self.pml.physical_box
        ol, or_, ob, ot = self.pml.outer_box
        custom_mesh = mesh is not None
        if not custom_mesh:
            x = _axis_nodes([ol, xl, MATERIAL_BOX[0], MATERIAL_BOX[1], xr, or_], self.h)
            y = _axis_nodes([ob, yb, MATERIAL_BOX[2], MATERIAL_BOX[3], yt, ot], self.h)
            mesh = MeshTri.init_tensor(x, y)
        self.mesh = mesh
        max_diameter = _validate_mesh(self.mesh, self.pml.outer_box)
        # Exact P2/P3 scalar DOF count from topology, before Basis allocates all
        # global quadrature/basis arrays.  This gate never starts a PDE action.
        expected_dofs = (self.mesh.nvertices+(self.order-1)*self.mesh.nfacets
                         + (self.mesh.nelements if self.order == 3 else 0))
        self.estimated_factor_bytes = _factor_memory_estimate(expected_dofs)
        if (expected_dofs > self.max_dofs
                or self.estimated_factor_bytes > self.max_factor_bytes):
            self._event('backend_memory_reject', ndofs=int(expected_dofs),
                        max_dofs=self.max_dofs,
                        estimated_factor_bytes=self.estimated_factor_bytes,
                        max_factor_bytes=self.max_factor_bytes,
                        basis_constructed=False)
            raise MemoryError(f'FEM pre-Basis gate: {expected_dofs} DOFs '
                              f'(limit {self.max_dofs}), estimated factor '
                              f'{self.estimated_factor_bytes} bytes '
                              f'(limit {self.max_factor_bytes})')
        element = ElementTriP2() if self.order == 2 else ElementTriP3()
        self.basis = Basis(self.mesh, element, intorder=self.quadrature_order)
        self.ndofs = int(self.basis.N)
        self.nelements = int(self.mesh.nelements)
        self._boundary_dofs = self.basis.get_dofs().all()
        self._free_dofs = np.setdiff1d(np.arange(self.ndofs), self._boundary_dofs)
        centroids = self.mesh.p[:, self.mesh.t].mean(axis=1).T
        self._material_elements = np.flatnonzero(_inside(centroids, MATERIAL_BOX))
        self.material_basis = self.basis.with_elements(self._material_elements)
        q = np.asarray(self.material_basis.global_coordinates())
        self._material_shape = q.shape[1:]
        self.material_points = q.reshape(2, -1).T.copy()
        self.material_weights = self.material_basis.dx.reshape(-1).copy()
        self.n_material_quadrature = len(self.material_points)
        # QP arrays follow (element, quadrature); no interpolation of epsilon.
        global_q = np.asarray(self.basis.global_coordinates())
        sx, sy = self.pml.stretch(global_q.reshape(2, -1).T)
        shape = global_q.shape[1:]
        background = asm(_background, self.basis,
                         ax=(sy/sx).reshape(shape),
                         ay=(sx/sy).reshape(shape),
                         volume=(sx*sy).reshape(shape), k2=self.k**2)
        self._background_matrix = background[self._free_dofs][:, self._free_dofs].tocsc()
        self._receiver_probe = self.basis.probes(self.rx.T).tocsr()
        # Exact free-space point-receiver kernel; only the volume integral is
        # quadrature-discretized.  Cache once for either extraction/audit mode.
        self._receiver_green = line_green(self.rx, self.material_points, self.k)
        self._factor = None
        self._matrix = None
        self._epsilon_cache = None
        self._factor_storage_bytes = 0
        self._backend_token = object()
        self.setup_seconds = time.perf_counter()-started
        self.metadata = dict(
            method='independent scattered-field P2/P3 scalar-TM FEM',
            freq=self.freq, k=self.k, order=self.order, h=self.h,
            mesh_kind=('caller-supplied conforming triangle mesh' if custom_mesh
                       else 'Cartesian tensor triangles'),
            h_role=('custom mesh label; use max_triangle_diameter for actual mesh'
                    if custom_mesh else 'upper bound on Cartesian intervals'),
            max_axis_step=(None if custom_mesh
                           else float(max(np.diff(x).max(), np.diff(y).max()))),
            max_triangle_diameter=max_diameter,
            quadrature_order=self.quadrature_order, ndofs=self.ndofs,
            free_dofs=len(self._free_dofs), nelements=self.nelements,
            material_elements=len(self._material_elements),
            material_quadrature_points=self.n_material_quadrature,
            material_area_quadrature=float(self.material_weights.sum()),
            max_dofs=self.max_dofs, max_factor_bytes=self.max_factor_bytes,
            estimated_factor_bytes=self.estimated_factor_bytes,
            background_matrix_bytes=_sparse_bytes(self._background_matrix),
            receiver_mode=self.receiver_mode,
            receiver_green_bytes=int(self._receiver_green.nbytes),
            material_box=list(MATERIAL_BOX), pml=self.pml.as_dict(),
            geometry='rectangular box conforming; internal interfaces quadrature sampled',
            source='analytic i*H0^(1)(k*r)/4; no delta source in FEM mesh',
            time_convention='exp(-i omega t)',
            convergence_status='NOT_QUALIFIED_BY_THIS_MODULE')
        self._event('backend_setup', wall_seconds=self.setup_seconds, **self.metadata)

    def _event(self, event, **metadata):
        emitter = getattr(self.book, 'event', None)
        if callable(emitter):
            emitter(event=event, phase=self.phase, charged_cpu_seconds=0,
                    backend='independent_fem_pml', **metadata)

    def _material(self, epsilon: Callable):
        if not callable(epsilon):
            raise TypeError('epsilon must be a continuous-coefficient callable(points Nx2)')
        value = np.asarray(epsilon(self.material_points.copy()), dtype=np.complex128)
        if value.ndim == 0:
            value = np.full(self.n_material_quadrature, value, dtype=np.complex128)
        if value.shape != (self.n_material_quadrature,):
            raise ValueError('epsilon callable must return scalar or complex (n_points,)')
        if (not np.isfinite(value).all() or np.any(value.real < 1.0)
                or np.any(value.real > 9.0)):
            raise ValueError('Requires finite 1<=Re(epsilon)<=9')
        if np.any(value.imag < 0.0) or np.any(self.freq*value.imag > 0.25):
            raise ValueError('Requires passive 0<=f*Im(epsilon)<=0.25')
        return value.copy()

    def _factor_for(self, epsilon):
        if self._factor is not None and np.array_equal(epsilon, self._epsilon_cache):
            return True, 0.0, 0.0
        # Clear before attempting new work, so a failed assembly/reservation/LU
        # cannot leave a stale factor eligible for a new material.
        self._factor = None
        self._matrix = None
        self._epsilon_cache = None
        self._factor_storage_bytes = 0
        started = time.perf_counter()
        mass = asm(_contrast_mass, self.material_basis,
                   contrast=(epsilon-1.0).reshape(self._material_shape))
        matrix = (self._background_matrix
                  - self.k**2*mass[self._free_dofs][:, self._free_dofs]).tocsc()
        if not np.isfinite(matrix.data).all():
            raise RuntimeError('Nonfinite FEM/PML matrix')
        assembly_seconds = time.perf_counter()-started
        self.book.reserve(self.phase, 'fem_material_frequency_LU',
                          rhs=0, factorizations=1, freq=self.freq,
                          backend='independent_fem_pml', ndofs=self.ndofs,
                          free_dofs=len(self._free_dofs), order=self.order,
                          h=self.h, matrix_nnz=int(matrix.nnz), pml=self.pml.as_dict())
        tic = time.perf_counter()
        factor = splu(matrix, permc_spec='COLAMD')
        factor_seconds = time.perf_counter()-tic
        factor_storage_bytes = (_sparse_bytes(factor.L)+_sparse_bytes(factor.U)
                                + factor.perm_r.nbytes+factor.perm_c.nbytes)
        if factor_storage_bytes > self.max_factor_bytes:
            self._event('backend_memory_reject', ndofs=self.ndofs,
                        actual_factor_storage_bytes=factor_storage_bytes,
                        estimated_factor_bytes=self.estimated_factor_bytes,
                        max_factor_bytes=self.max_factor_bytes,
                        basis_constructed=True, factorization_completed=True)
            raise MemoryError('Actual SuperLU L/U storage exceeds max_factor_bytes')
        self._matrix = matrix
        self._factor = factor
        self._epsilon_cache = epsilon.copy()
        self._factor_storage_bytes = int(factor_storage_bytes)
        self._event('backend_factor', freq=self.freq, wall_seconds=factor_seconds,
                    assembly_seconds=assembly_seconds, matrix_nnz=int(matrix.nnz),
                    ndofs=self.ndofs, free_dofs=len(self._free_dofs),
                    actual_factor_storage_bytes=self._factor_storage_bytes,
                    actual_matrix_bytes=_sparse_bytes(matrix),
                    factor_nnz=int(factor.L.nnz+factor.U.nnz),
                    estimated_factor_bytes=self.estimated_factor_bytes)
        return False, assembly_seconds, factor_seconds

    def _illumination(self, angles, plane_wave_origin, amplitudes):
        if angles is None:
            kind = 'line'
            definition = dict(kind=kind, positions=self.tx.copy())
            count = len(self.tx)
        else:
            a = np.atleast_1d(np.asarray(angles, dtype=float))
            if a.ndim != 1 or not len(a) or not np.isfinite(a).all():
                raise ValueError('angles must be a finite nonempty vector in radians')
            origin = np.asarray(plane_wave_origin, dtype=float)
            if origin.shape != (2,) or not np.isfinite(origin).all():
                raise ValueError('plane_wave_origin must be a finite (2,) coordinate')
            kind = 'plane_wave'
            definition = dict(kind=kind, angles=a.copy(), origin=origin.copy(),
                              directions=np.c_[np.cos(a), np.sin(a)])
            count = len(a)
        amp = np.asarray(1.0 if amplitudes is None else amplitudes, dtype=complex)
        if amp.ndim == 0:
            amp = np.repeat(amp, count)
        if amp.shape != (count,) or not np.isfinite(amp).all():
            raise ValueError('amplitudes must be finite scalar or one value per source')
        definition['amplitudes'] = self.source_normalization*amp
        return definition

    def _incident(self, points, definition):
        if definition['kind'] == 'line':
            values = line_green(points, definition['positions'], self.k)
        else:
            phase = ((np.asarray(points)-definition['origin'])
                     @ definition['directions'].T)
            values = np.exp(1j*self.k*phase)
        return values*definition['amplitudes'][None, :]

    def solve(self, epsilon, *, angles=None, plane_wave_origin=(0.0, 0.0),
              amplitudes=None):
        """Solve one coefficient for all line sources or supplied plane angles.

        ``angles=None`` uses the constructor's tx line sources.  Otherwise each
        angle defines exp(i*k*d(angle).(x-origin)); no cell-average/sinc factor
        is applied.  Plane-wave ``data`` contain scattering only (``direct=0``).
        Both modes reserve their exact number of RHS before forming/solving it.
        """
        started = time.perf_counter()
        definition = self._illumination(angles, plane_wave_origin, amplitudes)
        eps = self._material(epsilon)
        reused, assembly_seconds, factor_seconds = self._factor_for(eps)
        count = len(definition['amplitudes'])
        self.book.reserve(self.phase, 'fem_forward_sources', rhs=count,
                          factorizations=0, freq=self.freq,
                          backend='independent_fem_pml', source_kind=definition['kind'],
                          factor_reused=bool(reused), ndofs=self.ndofs,
                          receiver_mode=self.receiver_mode)
        rhs_started = time.perf_counter()
        incident = self._incident(self.material_points, definition)
        rhs = np.empty((len(self._free_dofs), count), dtype=np.complex128, order='F')
        for j in range(count):
            load = (self.k**2*(eps-1.0)*incident[:, j]).reshape(self._material_shape)
            rhs[:, j] = asm(_forcing, self.material_basis, load=load)[self._free_dofs]
        if not np.isfinite(rhs).all():
            raise RuntimeError('Nonfinite FEM scattered-field RHS')
        rhs_seconds = time.perf_counter()-rhs_started
        tic = time.perf_counter()
        solved = self._factor.solve(rhs)
        solve_seconds = time.perf_counter()-tic
        if not np.isfinite(solved).all():
            raise RuntimeError('Nonfinite FEM/PML solution')
        defects = self._matrix@solved-rhs
        residuals = (np.linalg.norm(defects, axis=0)
                     / np.maximum(np.linalg.norm(rhs, axis=0), np.finfo(float).tiny))
        residual = float(np.max(residuals))
        if not np.isfinite(residual) or residual > self.residual_tolerance:
            self._event('backend_residual_failure', freq=self.freq, residual=residual,
                        tolerance=self.residual_tolerance, rhs=count)
            raise RuntimeError(f'FEM/PML RHS-relative residual {residual} exceeds tolerance')
        dofs = np.zeros((self.ndofs, count), dtype=np.complex128)
        dofs[self._free_dofs] = solved
        probe_scattered = np.asarray(self._receiver_probe@dofs).T
        incident_rx = self._incident(self.rx, definition).T
        direct = (incident_rx.copy() if definition['kind'] == 'line'
                  else np.zeros_like(incident_rx))
        material_scattered = np.column_stack([
            np.asarray(self.material_basis.interpolate(dofs[:, j])).reshape(-1)
            for j in range(count)])
        total = incident+material_scattered
        current = (eps-1.0)[:, None]*total
        scattered = (probe_scattered if self.receiver_mode == 'probe'
                     else (self.k**2*self._receiver_green
                           @ (self.material_weights[:, None]*current)).T)
        timings = dict(setup_seconds=self.setup_seconds,
                       material_assembly_seconds=assembly_seconds,
                       factor_seconds=factor_seconds, rhs_assembly_seconds=rhs_seconds,
                       forward_solve_seconds=solve_seconds,
                       state_seconds=time.perf_counter()-started)
        metadata = dict(self.metadata, factor_reused=bool(reused), rhs=count,
                        actual_factor_storage_bytes=self._factor_storage_bytes,
                        actual_matrix_bytes=_sparse_bytes(self._matrix),
                        source_kind=definition['kind'],
                        receiver_data=('total line-source field' if definition['kind'] == 'line'
                                       else 'scattered plane-wave field; direct explicitly zero'))
        self._event('backend_solve', freq=self.freq, rhs=count,
                    factor_reused=bool(reused), residual=residual,
                    wall_seconds=solve_seconds, rhs_assembly_seconds=rhs_seconds)
        self._event('backend_state', freq=self.freq, rhs=count,
                    residual=residual, timings=timings,
                    source_kind=definition['kind'], jacobian=False,
                    receiver_mode=self.receiver_mode)
        return State(epsilon=eps, material_points=self.material_points.copy(),
                     material_weights=self.material_weights.copy(), incident=incident,
                     field=total, current=current, scattered_dofs=dofs,
                     data=direct+scattered, scattered=scattered, direct=direct,
                     incident_at_receivers=incident_rx, factor=self._factor,
                     residual=residual, residuals=residuals, timings=timings,
                     illumination=definition, metadata=metadata,
                     receiver_mode=self.receiver_mode,
                     probe_data=direct+probe_scattered,
                     probe_scattered=probe_scattered,
                     _backend_token=self._backend_token, _matrix=self._matrix)

    def volume_data(self, state):
        """Return direct+k² integral(g(rx,x)*chi*E_total dx), source-major.

        This postprocesses the same solved state without any additional RHS.
        The analytic free-space receiver Green is exact at quadrature points;
        the material field and volume quadrature retain their numerical errors.
        Plane-wave ``direct`` remains exactly zero, as in :meth:`solve`.
        """
        if state._backend_token is not self._backend_token:
            raise ValueError('State belongs to a different FEM mesh/backend')
        return (state.direct + (self.k**2*self._receiver_green
                @ (state.material_weights[:, None]*state.current)).T)

    def _state_action(self, state, rhs, *, trans, label):
        """Numerical action on an already-reserved RHS and its state's factor."""
        tic = time.perf_counter()
        result = state.factor.solve(np.asfortranarray(rhs), trans=trans)
        elapsed = time.perf_counter()-tic
        if not np.isfinite(result).all():
            raise RuntimeError(f'Nonfinite FEM {label} solution')
        matrix = state._matrix if trans == 'N' else state._matrix.conjugate().T
        residuals = (np.linalg.norm(matrix@result-rhs, axis=0)
                     / np.maximum(np.linalg.norm(rhs, axis=0), np.finfo(float).tiny))
        residual = float(residuals.max())
        self._event('backend_solve', freq=self.freq, label=label,
                    rhs=rhs.shape[1], wall_seconds=elapsed, residual=residual,
                    factor_reused=True, transpose=trans)
        if not np.isfinite(residual) or residual > self.residual_tolerance:
            raise RuntimeError(f'FEM {label} residual {residual} exceeds tolerance')
        return result

    def jvp(self, state, delta_epsilon):
        """Derivative of receiver data along a complex coefficient direction.

        ``delta_epsilon`` is callable(points Nx2)->complex N (or scalar); unlike
        a physical material, its values can have either sign.  It is evaluated
        at the state's actual material quadrature points.  The linearization
        solves L*delta_u=k**2*delta_epsilon*u_total with the existing factor.
        Return shape is (source,receiver).  Exactly n_sources RHS are reserved.
        """
        if state._backend_token is not self._backend_token:
            raise ValueError('State belongs to a different FEM mesh/backend')
        if not callable(delta_epsilon):
            raise TypeError('delta_epsilon must be a coefficient-direction callable')
        direction = np.asarray(delta_epsilon(state.material_points.copy()), dtype=complex)
        if direction.ndim == 0:
            direction = np.full(len(state.material_points), direction, dtype=complex)
        if direction.shape != state.epsilon.shape or not np.isfinite(direction).all():
            raise ValueError('Direction must return finite scalar or complex (n_points,)')
        count = state.field.shape[1]
        self.book.reserve(self.phase, 'fem_material_jvp', rhs=count, factorizations=0,
                          freq=self.freq, backend='independent_fem_pml',
                          factor_reused=True, ndofs=self.ndofs,
                          receiver_mode=state.receiver_mode)
        rhs = np.empty((len(self._free_dofs), count), dtype=complex, order='F')
        for j in range(count):
            load = (self.k**2*direction*state.field[:, j]).reshape(self._material_shape)
            rhs[:, j] = asm(_forcing, self.material_basis, load=load)[self._free_dofs]
        if not np.isfinite(rhs).all():
            raise RuntimeError('Nonfinite FEM coefficient-JVP RHS')
        derivative = self._state_action(state, rhs, trans='N', label='material_jvp')
        dofs = np.zeros((self.ndofs, count), dtype=complex)
        dofs[self._free_dofs] = derivative
        if state.receiver_mode == 'volume':
            sampled = np.column_stack([
                np.asarray(self.material_basis.interpolate(dofs[:, j])).reshape(-1)
                for j in range(count)])
            delta_current = (direction[:, None]*state.field
                             + (state.epsilon-1.0)[:, None]*sampled)
            return (self.k**2*self._receiver_green
                    @ (state.material_weights[:, None]*delta_current)).T
        return np.asarray(self._receiver_probe@dofs).T

    def vjp(self, state, data_cotangent):
        """Real-adjoint coefficient gradient at material quadrature points.

        Input has complex shape (source,receiver).  Output g has complex shape
        (material_quadrature,) and includes triangle integration weights:

          Re(vdot(data_cotangent, jvp(state,delta_epsilon)))
            = Re(vdot(g, delta_epsilon(material_points))).

        Therefore independent real/imaginary sample-coordinate gradients are
        g.real and g.imag.  A low-dimensional material parameter gradient must
        additionally apply that parameterization's chain rule.  This is not a
        continuum density: divide by material_weights only if one is wanted.
        Exactly n_sources conjugate-transpose RHS are reserved.
        """
        if state._backend_token is not self._backend_token:
            raise ValueError('State belongs to a different FEM mesh/backend')
        cotangent = np.asarray(data_cotangent, dtype=complex)
        if cotangent.shape != state.data.shape or not np.isfinite(cotangent).all():
            raise ValueError('data_cotangent must match finite source-major data shape')
        count = state.field.shape[1]
        self.book.reserve(self.phase, 'fem_material_vjp', rhs=count, factorizations=0,
                          freq=self.freq, backend='independent_fem_pml',
                          factor_reused=True, ndofs=self.ndofs,
                          transpose='conjugate', receiver_mode=state.receiver_mode)
        receiver_back = None
        if state.receiver_mode == 'volume':
            receiver_back = self._receiver_green.conjugate().T@cotangent.T
            rhs = np.empty((len(self._free_dofs), count), dtype=complex, order='F')
            for j in range(count):
                # asm(_forcing) supplies the quadrature weights once:
                # B.T [conj(chi) * (k²*weights*g_rx.H*y)].
                load = (self.k**2*np.conjugate(state.epsilon-1.0)
                        * receiver_back[:, j]).reshape(self._material_shape)
                rhs[:, j] = asm(_forcing, self.material_basis, load=load)[self._free_dofs]
        else:
            rhs = np.asarray(self._receiver_probe[:, self._free_dofs].conjugate().T
                             @ cotangent.T)
        if not np.isfinite(rhs).all():
            raise RuntimeError('Nonfinite FEM coefficient-VJP RHS')
        adjoint = self._state_action(state, rhs, trans='H', label='material_vjp')
        dofs = np.zeros((self.ndofs, count), dtype=complex)
        dofs[self._free_dofs] = adjoint
        sampled = np.column_stack([
            np.asarray(self.material_basis.interpolate(dofs[:, j])).reshape(-1)
            for j in range(count)])
        if receiver_back is not None:
            sampled = sampled+receiver_back
        return (self.k**2*state.material_weights
                * np.sum(np.conjugate(state.field)*sampled, axis=1))

    def evaluate(self, state, points, *, total=True):
        """Probe an existing state inside physical_box; no additional PDE solve.

        Return shape (source,point).  ``total=True`` adds the analytic incident
        field for both illumination modes.  A line source point itself cannot
        be probed in total mode.  Values inside PML are deliberately excluded.
        """
        if state._backend_token is not self._backend_token:
            raise ValueError('State belongs to a different FEM mesh/backend')
        p = _coordinates(points, 'field evaluation')
        if not _inside(p, self.pml.physical_box, closed=True).all():
            raise ValueError('Field evaluation must be within physical_box')
        result = np.asarray(self.basis.probes(p.T)@state.scattered_dofs)
        if total:
            result = result+self._incident(p, state.illumination)
        return result.T

    def far_field(self, state, angles):
        """Return F(source,angle) from the existing material total field.

        With direction d=(cos(theta),sin(theta)), this normalization is

          F(theta)=k**2 integral_D (epsilon-1)*u_total(x)*exp(-i*k*d.x) dx,
          u_s(r*d) ~ exp(i*(k*r+pi/4))/sqrt(8*pi*k*r) * F(theta).

        The same FEM material quadrature is used for postprocessing; this is
        not an independent reference or an additional linear solve.
        """
        if state._backend_token is not self._backend_token:
            raise ValueError('State belongs to a different FEM mesh/backend')
        a = np.atleast_1d(np.asarray(angles, dtype=float))
        if a.ndim != 1 or not len(a) or not np.isfinite(a).all():
            raise ValueError('far-field angles must be a finite nonempty vector')
        directions = np.c_[np.cos(a), np.sin(a)]
        kernel = np.exp(-1j*self.k*(directions@state.material_points.T))
        return (self.k**2*kernel@(state.material_weights[:, None]*state.current)).T
