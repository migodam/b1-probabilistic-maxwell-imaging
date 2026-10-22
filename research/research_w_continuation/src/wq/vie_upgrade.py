"""Explicit, optional refinements of the frozen P0 scalar-TM VIE.

The default kernel, source normalization, material law, antenna contract and
source-major data layout are those of ``p0.vie.VIE``.  Nothing in this module
changes qualification tolerances.  ``field_mode='active'`` skips unused air
fields; entries not reconstructed are NaN, while data remain complete.

Two spatial variants are available for controlled qualification:

* cell-average testing along selected axes (pulse trial / averaged test);
* subcell product integration of chi(y)*g(x-y), with a constant total-field
  trial per parent cell.  This is not mean(chi)*integral(g).  The subcell
  material is sampled from the caller's coefficient function, and the Green
  function is integrated over each microcell, including its actual singular
  rectangle when necessary.  It does not claim high-order field accuracy.

All physical LU/RHS actions reserve the supplied Book.  The job wrapper still
owns process-CPU charging, including setup and residual evaluation.
"""
from dataclasses import dataclass, field
import time
import warnings

import numpy as np
from scipy.integrate import quad
from scipy.linalg import LinAlgWarning, lu_factor
from scipy.special import hankel1

from p0.vie import VIE, State, point_green, source_integral, square_self_green
from p0.geometry import material_grid, wall_values, in_query


def _order(value, name, minimum=2):
    if int(value) != value or value < minimum:
        raise ValueError(f'{name} must be integer >= {minimum}')
    return int(value)


def rectangle_self_green(k, h, offset=(0., 0.)):
    """Integral over a square when the observation lies inside or on it.

    ``offset`` is observation minus square center.  Four rectangular quadrants
    about the observation use the same analytic radial antiderivative as P0;
    only a smooth one-dimensional angular quadrature remains.  This includes
    corner/edge singularities of the source microcells.
    """
    offset = np.asarray(offset, float)
    if offset.shape != (2,) or not np.isfinite(offset).all():
        raise ValueError('offset requires two finite coordinates')
    if k <= 0 or h <= 0:
        raise ValueError('positive k and h required')
    tolerance = 32*np.finfo(float).eps*h
    if np.any(np.abs(offset) > h/2+tolerance):
        raise ValueError('observation must be inside the closed square')
    if np.array_equal(offset, np.zeros(2)):
        return square_self_green(k, h)
    offset = np.clip(offset, -h/2, h/2)

    def radial(radius):
        if radius == 0:
            return 0j
        return .25j*(radius*hankel1(1, k*radius)/k + 2j/(np.pi*k*k))

    def integrate(fun, lo, hi):
        if hi <= lo:
            return 0j
        return (quad(lambda t: fun(t).real, lo, hi, epsabs=2e-12,
                     epsrel=2e-12)[0]
                + 1j*quad(lambda t: fun(t).imag, lo, hi, epsabs=2e-12,
                          epsrel=2e-12)[0])

    result = 0j
    for length_x in (h/2-offset[0], h/2+offset[0]):
        for length_y in (h/2-offset[1], h/2+offset[1]):
            if length_x == 0 or length_y == 0:
                continue
            split = np.arctan2(length_y, length_x)
            result += integrate(lambda t: radial(length_x/np.cos(t)), 0., split)
            result += integrate(lambda t: radial(length_y/np.sin(t)), split, np.pi/2)
    return result


def offset_source_integrals(offsets, k, h, quadrature=4):
    """Integrated source squares; all contained singular observations repaired."""
    offsets = np.asarray(offsets, float)
    if offsets.ndim != 2 or offsets.shape[1] != 2:
        raise ValueError('offsets require (n,2) coordinates')
    quadrature = _order(quadrature, 'quadrature')
    contained = np.all(np.abs(offsets) <= h/2+32*np.finfo(float).eps*h, axis=1)
    safe = offsets.copy()
    safe[contained] = (10*h, 10*h)
    values = source_integral(safe, np.zeros((1, 2)), k, h, quadrature)[:, 0]
    for index in np.flatnonzero(contained):
        values[index] = rectangle_self_green(k, h, offsets[index])
    return values


def square_averaged_self_green(k, h):
    """Source-square integral averaged over the identical observation square.

    The difference variable has the triangular overlap weight
    (1-|x|/h)*(1-|y|/h).  Polar integration over one eighth of its square
    support handles the logarithmic singularity without coincident nodes.
    """
    if k <= 0 or h <= 0:
        raise ValueError('positive k and h required')

    def angular(theta):
        cosine, sine = np.cos(theta), np.sin(theta)
        radius = h/cosine

        def radial(r):
            if r == 0:
                return 0j
            overlap = (1-r*cosine/h)*(1-r*sine/h)
            return .25j*hankel1(0, k*r)*r*overlap

        return (quad(lambda r: radial(r).real, 0., radius,
                     epsabs=2e-12, epsrel=2e-12)[0]
                + 1j*quad(lambda r: radial(r).imag, 0., radius,
                          epsabs=2e-12, epsrel=2e-12)[0])

    return 8*(quad(lambda t: angular(t).real, 0., np.pi/4,
                   epsabs=2e-12, epsrel=2e-12)[0]
              + 1j*quad(lambda t: angular(t).imag, 0., np.pi/4,
                        epsabs=2e-12, epsrel=2e-12)[0])


def _observation_rule(axes, order, h):
    axes = tuple(bool(v) for v in axes)
    if len(axes) != 2:
        raise ValueError('field_average_axes must contain (x,y) flags')
    nodes, weights = np.polynomial.legendre.leggauss(_order(order, 'field_quadrature'))
    rules = [(nodes*h/2, weights/2) if enabled
             else (np.zeros(1), np.ones(1)) for enabled in axes]
    return [(np.array([x, y]), wx*wy)
            for x, wx in zip(*rules[0]) for y, wy in zip(*rules[1])]


def integrated_offset_kernel(nx, ny, k, h, quadrature=4, *,
                             source_offset=(0., 0.), source_h=None,
                             field_average_axes=(False, False),
                             field_quadrature=4):
    """Toeplitz source-square integral, optionally averaged over a test cell."""
    dx = np.arange(1-nx, nx)*h
    dy = np.arange(1-ny, ny)*h
    ox, oy = np.meshgrid(dx, dy, indexing='xy')
    offsets = np.c_[ox.ravel(), oy.ravel()]-np.asarray(source_offset)
    source_h = h if source_h is None else float(source_h)
    values = np.zeros(len(offsets), complex)
    for shift, weight in _observation_rule(field_average_axes, field_quadrature, h):
        values += weight*offset_source_integrals(offsets+shift, k, source_h, quadrature)
    if (tuple(bool(v) for v in field_average_axes) == (True, True)
            and source_h == h and np.array_equal(np.asarray(source_offset), np.zeros(2))):
        middle = (ny-1)*(2*nx-1)+nx-1
        values[middle] = square_averaged_self_green(k, h)
    return (k*k*values).reshape(2*ny-1, 2*nx-1)


@dataclass(frozen=True)
class SubcellMaterial:
    """A sampled source-cell coefficient; global active columns are explicit.

    ``contrast[s,j]`` is constant on source microcell s in active parent j.
    ``epsilon`` retains the exact old material_grid values when created by
    ``frozen_subcell_material``.  Small reduction-roundoff differences between
    that array and 1+mean(contrast) are checked rather than repaired.
    """
    epsilon: np.ndarray
    active: np.ndarray
    contrast: np.ndarray
    offsets: np.ndarray
    subsamples: int
    nx: int
    ny: int
    h: float
    frequency: float
    provenance: dict = field(default_factory=dict)


def sample_subcell_material(backend, coefficient, subsamples=4, *,
                           epsilon=None, provenance=None):
    """Sample a shared continuous material law without a physical solve.

    Only active source-cell samples are retained.  ``coefficient(points)`` must
    return complex epsilon at the backend's frequency, including exact air.
    The midpoint locations and loop order match frozen material_grid.
    """
    subsamples = _order(subsamples, 'subsamples', minimum=1)
    positions = (np.arange(subsamples)+.5)/subsamples-.5
    offsets = backend.h*np.array([(a, b) for a in positions for b in positions])

    def evaluate(points):
        value = np.asarray(coefficient(points), complex)
        if value.shape != (len(points),) or not np.isfinite(value).all():
            raise ValueError('coefficient must return finite epsilon per point')
        if (np.any(value.real < 1) or np.any(value.real > 9)
                or np.any(value.imag < 0) or np.any(value.imag*backend.freq > .15)):
            raise ValueError('subcell coefficient violates frozen passive material bounds')
        return value

    if epsilon is None:
        averaged = np.zeros(backend.N, complex)
        for shift in offsets:
            averaged += evaluate(backend.points+shift)
        epsilon = averaged/len(offsets)
    epsilon = np.asarray(epsilon, complex).ravel().copy()
    if epsilon.shape != (backend.N,) or not np.isfinite(epsilon).all():
        raise ValueError('epsilon must be finite and match the grid')
    active = np.flatnonzero(epsilon != 1)
    values = np.empty((len(offsets), len(active)), complex)
    for index, shift in enumerate(offsets):
        values[index] = evaluate(backend.points[active]+shift)-1
    if not np.allclose(1+values.mean(axis=0), epsilon[active], rtol=2e-14, atol=2e-14):
        raise ValueError('provided cell means do not match the subcell material')
    return SubcellMaterial(epsilon, active, values, offsets, subsamples,
                           backend.nx, backend.ny, backend.h, backend.freq,
                           dict(provenance or {}))


def frozen_subcell_material(backend, seed, subsamples=4, probe=None):
    """Use the frozen wall, query clearing, probe override and epsilon(f) law."""
    epsilon = material_grid(backend, seed, subsamples=subsamples, probe=probe)

    def coefficient(points):
        er, sigma = wall_values(points, seed)
        query = in_query(points)
        er = np.where(query, 1., er)
        sigma = np.where(query, 0., sigma)
        if probe is not None:
            center, radius, contrast = probe
            target = np.linalg.norm(points-np.asarray(center), axis=1) <= radius
            er = np.where(target, 1+contrast, er)
            sigma = np.where(target, .02, sigma)
        return er+1j*sigma/backend.freq

    return sample_subcell_material(backend, coefficient, subsamples,
                                   epsilon=epsilon,
                                   provenance=dict(kind='frozen wall_values',
                                                   seed=int(seed), probe=probe))


@dataclass
class UpgradeState(State):
    field_mode: str = 'full'
    source_representation: str = 'cell_mean_current'
    method_metadata: dict = field(default_factory=dict)


class VIEUpgrade(VIE):
    """Drop-in constant-cell VIE plus explicit bounded spatial variants.

    ``solve(epsilon)`` retains the old current formulation and Jacobian.
    ``solve_subcells(material)`` uses total-field unknowns and contrast-weighted
    source integration.  A subcell Jacobian is deliberately unavailable.
    """
    def __init__(self, nx, freq, tx, rx, book, phase='validation', quadrature=4,
                 source_normalization=1., residual_tolerance=1e-10,
                 source_cell_average=False, *, field_mode='full',
                 field_average_axes=(False, False), field_quadrature=4,
                 max_active_cells=4500, max_subcell_kernel_bytes=192*1024**2):
        if field_mode not in ('full', 'active'):
            raise ValueError("field_mode must be 'full' or 'active'")
        self.field_mode = field_mode
        self.field_average_axes = tuple(bool(v) for v in field_average_axes)
        if len(self.field_average_axes) != 2:
            raise ValueError('field_average_axes requires (x,y)')
        self.field_quadrature = _order(field_quadrature, 'field_quadrature')
        self.max_active_cells = _order(max_active_cells, 'max_active_cells', 1)
        if self.max_active_cells > 4500:
            raise ValueError('the frozen 4500 active-cell memory bound cannot be relaxed')
        self.max_subcell_kernel_bytes = int(max_subcell_kernel_bytes)
        if self.max_subcell_kernel_bytes <= 0:
            raise ValueError('positive subcell kernel memory limit required')
        if any(self.field_average_axes) and not source_cell_average:
            raise ValueError('averaged field tests require matching source_cell_average=True')
        super().__init__(nx, freq, tx, rx, book, phase, quadrature,
                         source_normalization, residual_tolerance, source_cell_average)
        extra_started = time.perf_counter()
        if any(self.field_average_axes):
            self.kernel = integrated_offset_kernel(
                self.nx, self.ny, self.k, self.h, self.quadrature,
                field_average_axes=self.field_average_axes,
                field_quadrature=self.field_quadrature)
            self.incident = np.zeros_like(self.incident)
            for shift, weight in _observation_rule(self.field_average_axes,
                                                   self.field_quadrature, self.h):
                self.incident += (weight*self.source_normalization
                                  * point_green(self.points+shift, self.tx, self.k))
        self.setup_seconds += time.perf_counter()-extra_started
        self.book.event(event='upgrade_setup', phase=self.phase,
                        field_mode=self.field_mode,
                        field_average_axes=self.field_average_axes,
                        field_quadrature=self.field_quadrature,
                        source_representation='cell_mean_current',
                        charged_cpu_seconds=0)

    def _factor_matrix(self, matrix, n, label):
        self.book.reserve(self.phase, label, factorizations=1, freq=self.freq,
                          active_cells=n, total_cells=self.N)
        started = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter('error', LinAlgWarning)
            factor = lu_factor(matrix, check_finite=False, overwrite_a=True)
        elapsed = time.perf_counter()-started
        self.book.event(event='backend_factor', phase=self.phase, active_cells=n,
                        wall_seconds=elapsed, charged_cpu_seconds=0)
        return factor, elapsed

    def solve(self, epsilon, jacobian=False):
        # The frozen path is reused verbatim when a complete field/Jacobian is
        # requested.  Only the explicit averaged-test kernel can differ.
        if jacobian and self.field_mode != 'full':
            raise ValueError('full field required for full-domain material Jacobian')
        if self.field_mode == 'full':
            if np.count_nonzero(self._material(epsilon) != 1) > self.max_active_cells:
                raise MemoryError('Active dense-LU bound exceeded; no approximate fallback')
            state = super().solve(epsilon, jacobian)
            return UpgradeState(**state.__dict__, field_mode='full',
                                method_metadata=self._metadata())
        started = time.perf_counter()
        eps = self._material(epsilon)
        chi = eps-1
        active = np.flatnonzero(chi != 0)
        n = len(active)
        if n > self.max_active_cells:
            raise MemoryError('Active dense-LU bound exceeded; no approximate fallback')
        current = np.zeros_like(self.incident)
        factor = None
        factor_seconds = solve_seconds = residual = 0.
        if n:
            matrix = np.empty((n, n), complex, order='F')
            for first in range(0, n, 128):
                block = np.arange(first, min(first+128, n))
                matrix[block] = -chi[active[block], None]*self.green_block(active[block], active)
                matrix[block, block] += 1
            rhs = chi[active, None]*self.incident[active]
            factor, factor_seconds = self._factor_matrix(matrix, n, 'material_frequency_LU')
            current[active], solve_seconds = self._solve(factor, rhs, label='forward_sources')
            residual_sq = 0.
            for first in range(0, n, 128):
                block = np.arange(first, min(first+128, n))
                defect = (current[active[block]]-chi[active[block], None]
                          * (self.green_block(active[block], active)@current[active])-rhs[block])
                residual_sq += float(np.linalg.norm(defect)**2)
            residual = float(np.sqrt(residual_sq)/max(np.linalg.norm(rhs), np.finfo(float).tiny))
            self._check_residual(residual)
        field_values = np.full_like(self.incident, np.nan+1j*np.nan)
        field_values[active] = current[active]/chi[active, None]
        scattered = self.S[:, active]@current[active]
        timings = dict(setup_seconds=self.setup_seconds, factor_seconds=factor_seconds,
                       forward_solve_seconds=solve_seconds, receiver_solve_seconds=0.,
                       state_seconds=time.perf_counter()-started)
        self._state_event(n, residual, timings, 'cell_mean_current')
        return UpgradeState(eps, active, self.incident.copy(), field_values, current,
                            (self.direct+scattered).T, scattered.T, factor, residual,
                            None, timings, field_mode='active',
                            method_metadata=self._metadata())

    def _metadata(self):
        return dict(field_average_axes=list(self.field_average_axes),
                    field_quadrature=self.field_quadrature,
                    source_quadrature=self.quadrature,
                    source_normalization=[self.source_normalization.real,
                                          self.source_normalization.imag],
                    field_complete=self.field_mode == 'full')

    def _check_residual(self, residual):
        if not np.isfinite(residual) or residual > self.residual_tolerance:
            raise RuntimeError(f'VIE backward residual {residual} exceeds tolerance')

    def _state_event(self, n, residual, timings, representation):
        self.book.event(event='backend_state', phase=self.phase, active_cells=n,
                        total_cells=self.N, freq=self.freq, residual=residual,
                        jacobian=False, timings=timings, field_mode=self.field_mode,
                        source_representation=representation, charged_cpu_seconds=0)

    def solve_subcells(self, material, jacobian=False):
        """Solve E = incident + integral(g*chi*E_cell), with actual microcells.

        The retained material memory is O(subsamples²*active); shifted kernels
        are O(subsamples²*N).  Their explicit byte limit is checked before
        allocation or any physical reservation.  This option is intended for
        modest subcell counts; geometry-only 16/32 sampling audits are separate.
        """
        if jacobian:
            raise NotImplementedError('Contrast-weighted subcell Jacobian is not qualified')
        if not isinstance(material, SubcellMaterial):
            raise TypeError('SubcellMaterial required')
        if ((material.nx, material.ny, material.h, material.frequency)
                != (self.nx, self.ny, self.h, self.freq)):
            raise ValueError('subcell material grid/frequency mismatch')
        started = time.perf_counter()
        eps = self._material(material.epsilon)
        active = np.flatnonzero(eps != 1)
        if not np.array_equal(material.active, active):
            raise ValueError('material active columns do not match epsilon')
        n = len(active)
        if n > self.max_active_cells:
            raise MemoryError('Active dense-LU bound exceeded; no approximate fallback')
        samples = material.subsamples**2
        if material.contrast.shape != (samples, n) or material.offsets.shape != (samples, 2):
            raise ValueError('invalid subcell material dimensions')
        if not np.isfinite(material.contrast).all():
            raise ValueError('finite subcell contrasts required')
        if (np.any(material.contrast.real < 0) or np.any(material.contrast.real > 8)
                or np.any(material.contrast.imag < 0)
                or np.any(material.contrast.imag*self.freq > .15)):
            raise ValueError('subcell contrast violates frozen passive material bounds')
        expected_positions = (np.arange(material.subsamples)+.5)/material.subsamples-.5
        expected_offsets = self.h*np.array([(a, b) for a in expected_positions
                                            for b in expected_positions])
        if not np.array_equal(material.offsets, expected_offsets):
            raise ValueError('subcell offsets do not match the declared uniform partition')
        if not np.allclose(1+material.contrast.mean(axis=0), eps[active],
                           rtol=2e-14, atol=2e-14):
            raise ValueError('subcell means do not match epsilon')
        required = samples*(2*self.ny-1)*(2*self.nx-1)*np.dtype(complex).itemsize
        if required > self.max_subcell_kernel_bytes:
            raise MemoryError(f'Subcell kernels need {required} bytes, exceeding explicit '
                              f'limit {self.max_subcell_kernel_bytes}; no silent coarsening')
        kernels = np.empty((samples, 2*self.ny-1, 2*self.nx-1), complex)
        receive = np.zeros((len(self.rx), n), complex)
        micro_h = self.h/material.subsamples
        for index, shift in enumerate(material.offsets):
            kernels[index] = integrated_offset_kernel(
                self.nx, self.ny, self.k, self.h, self.quadrature,
                source_offset=shift, source_h=micro_h,
                field_average_axes=self.field_average_axes,
                field_quadrature=self.field_quadrature)
            receive += (self.k**2*source_integral(self.rx, self.points[active]+shift,
                                                self.k, micro_h, self.quadrature)
                        * material.contrast[index][None, :])

        def weighted_block(rows):
            ry, rx = np.divmod(np.asarray(rows), self.nx)
            cy, cx = np.divmod(active, self.nx)
            dy = ry[:, None]-cy[None, :]+self.ny-1
            dx = rx[:, None]-cx[None, :]+self.nx-1
            result = np.zeros((len(rows), n), complex)
            for index in range(samples):
                result += kernels[index, dy, dx]*material.contrast[index][None, :]
            return result

        factor = None
        factor_seconds = solve_seconds = residual = 0.
        active_field = np.empty((n, len(self.tx)), complex)
        if n:
            matrix = np.empty((n, n), complex, order='F')
            for first in range(0, n, 128):
                block = np.arange(first, min(first+128, n))
                matrix[block] = -weighted_block(active[block])
                matrix[block, block] += 1
            factor, factor_seconds = self._factor_matrix(matrix, n, 'subcell_material_frequency_LU')
            active_field, solve_seconds = self._solve(factor, self.incident[active],
                                                       label='subcell_forward_sources')
            residual_sq = 0.
            for first in range(0, n, 128):
                block = np.arange(first, min(first+128, n))
                defect = (active_field[block]-weighted_block(active[block])@active_field
                          - self.incident[active[block]])
                residual_sq += float(np.linalg.norm(defect)**2)
            residual = float(np.sqrt(residual_sq)/max(np.linalg.norm(self.incident[active]),
                                                     np.finfo(float).tiny))
            self._check_residual(residual)
        if self.field_mode == 'full':
            field_values = self.incident.copy()
            for first in range(0, self.N, 256):
                rows = np.arange(first, min(first+256, self.N))
                field_values[rows] += weighted_block(rows)@active_field
        else:
            field_values = np.full_like(self.incident, np.nan+1j*np.nan)
        field_values[active] = active_field
        current = np.zeros_like(self.incident)
        current[active] = (eps[active, None]-1)*active_field
        scattered = receive@active_field
        timings = dict(setup_seconds=self.setup_seconds, factor_seconds=factor_seconds,
                       forward_solve_seconds=solve_seconds, receiver_solve_seconds=0.,
                       state_seconds=time.perf_counter()-started)
        representation = 'contrast_weighted_subcell_constant_total_field'
        self._state_event(n, residual, timings, representation)
        metadata = dict(self._metadata(), subsamples=material.subsamples,
                        subcell_kernel_bytes=required,
                        material_provenance=material.provenance,
                        current_is_mean_contrast_times_cell_field=True)
        return UpgradeState(eps, active, self.incident.copy(), field_values, current,
                            (self.direct+scattered).T, scattered.T, factor, residual,
                            None, timings, field_mode=self.field_mode,
                            source_representation=representation,
                            method_metadata=metadata)
