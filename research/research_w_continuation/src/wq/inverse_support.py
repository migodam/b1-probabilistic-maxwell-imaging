"""Full-domain pixel coordinates and exact discrete volume-receiver helpers.

No truth support, world generator, optimization algorithm, or experiment loop is
used here.  Parameters are all real contrast pixels followed by all conductivity
pixels, each flattened in (y,x) C order.  The physical material is
epsilon=1+contrast+i*conductivity/f.  Calibration air in the declared query
rectangle is imposed pointwise, independently of the candidate support.
"""

from dataclasses import dataclass
import time

import numpy as np
from skfem import MeshTri, asm

from .fem import PMLConfig, _forcing, _factor_memory_estimate


@dataclass(frozen=True)
class PixelGrid:
    nx: int = 48
    ny: int = 32
    calibration_box: tuple = (-1.2, 1.2, 1.6, 2.8)

    def __post_init__(self):
        for value, name in ((self.nx, 'nx'), (self.ny, 'ny')):
            if isinstance(value, bool) or int(value) != value or value <= 0:
                raise ValueError(f'{name} must be a positive integer')
        box = np.asarray(self.calibration_box, dtype=float)
        if (box.shape != (4,) or not np.isfinite(box).all()
                or not (-2 <= box[0] < box[1] <= 2 and 0 <= box[2] < box[3] <= 3)):
            raise ValueError('calibration_box must be a nonempty rectangle in D')
        object.__setattr__(self, 'nx', int(self.nx))
        object.__setattr__(self, 'ny', int(self.ny))
        object.__setattr__(self, 'calibration_box', tuple(float(x) for x in box))

    @property
    def n_pixels(self):
        return self.nx*self.ny

    @property
    def n_parameters(self):
        return 2*self.n_pixels

    @property
    def x_edges(self):
        return np.linspace(-2.0, 2.0, self.nx+1)

    @property
    def y_edges(self):
        return np.linspace(0.0, 3.0, self.ny+1)

    @property
    def cell_area(self):
        return 12.0/self.n_pixels

    def unpack(self, parameters):
        a = np.asarray(parameters)
        if np.iscomplexobj(a) or a.shape != (self.n_parameters,) or not np.isfinite(a).all():
            raise ValueError('parameters must be a finite real vector of length 2*nx*ny')
        a = np.asarray(a, dtype=float)
        return (a[:self.n_pixels].reshape(self.ny, self.nx).copy(),
                a[self.n_pixels:].reshape(self.ny, self.nx).copy())

    def pack(self, contrast, conductivity):
        a, b = np.asarray(contrast), np.asarray(conductivity)
        shapes = ((self.ny, self.nx), (self.n_pixels,))
        if (a.shape not in shapes or b.shape not in shapes
                or np.iscomplexobj(a) or np.iscomplexobj(b)
                or not np.isfinite(a).all() or not np.isfinite(b).all()):
            raise ValueError('pixel blocks must be finite real (ny,nx) or flat arrays')
        return np.r_[np.asarray(a, dtype=float).ravel(), np.asarray(b, dtype=float).ravel()]

    def indices(self, points):
        """C-order pixel ids; -1 outside closed D, right/top edge in last cell."""
        p = np.asarray(points, dtype=float)
        if p.ndim != 2 or p.shape[1] != 2 or not np.isfinite(p).all():
            raise ValueError('points must have finite shape (n,2)')
        valid = ((p[:, 0] >= -2) & (p[:, 0] <= 2)
                 & (p[:, 1] >= 0) & (p[:, 1] <= 3))
        result = np.full(len(p), -1, dtype=int)
        ix = np.minimum(np.searchsorted(self.x_edges, p[valid, 0], side='right')-1,
                        self.nx-1)
        iy = np.minimum(np.searchsorted(self.y_edges, p[valid, 1], side='right')-1,
                        self.ny-1)
        result[valid] = iy*self.nx+ix
        return result

    def calibration_mask(self, points):
        p = np.asarray(points, dtype=float)
        xl, xr, yb, yt = self.calibration_box
        return ((p[:, 0] >= xl) & (p[:, 0] <= xr)
                & (p[:, 1] >= yb) & (p[:, 1] <= yt))

    def available_areas(self, *, calibration=True):
        """Exact rectangular pixel areas outside forced calibration air."""
        areas = np.full((self.ny, self.nx), self.cell_area)
        if calibration:
            xl, xr, yb, yt = self.calibration_box
            wx = np.maximum(0.0, np.minimum(self.x_edges[1:], xr)
                            - np.maximum(self.x_edges[:-1], xl))
            wy = np.maximum(0.0, np.minimum(self.y_edges[1:], yt)
                            - np.maximum(self.y_edges[:-1], yb))
            areas = np.maximum(0.0, areas-wy[:, None]*wx[None, :])
        return areas


def project_parameters(parameters, grid=None, *, max_area=3.0, calibration=True,
                       conductivity_scale=8.0/0.15):
    """Bound and retain the highest-energy full-domain pixels deterministically.

    Energy is contrast**2+(conductivity_scale*conductivity)**2, with no support
    threshold or truth mask.  At most floor(max_area/cell_area) pixels survive,
    conservatively charging a full pixel even where calibration air cuts it.
    Cells wholly inside calibration air are zeroed before support ranking.
    """
    grid = PixelGrid() if grid is None else grid
    contrast, sigma = grid.unpack(parameters)
    if not np.isfinite(max_area) or not 0 <= max_area <= 12:
        raise ValueError('max_area must be finite and between zero and area(D)=12')
    if not np.isfinite(conductivity_scale) or conductivity_scale <= 0:
        raise ValueError('conductivity_scale must be positive and finite')
    contrast = np.clip(contrast, 0.0, 8.0)
    sigma = np.clip(sigma, 0.0, 0.15)
    if calibration:
        xl, xr, yb, yt = grid.calibration_box
        full_x = (grid.x_edges[:-1] >= xl) & (grid.x_edges[1:] <= xr)
        full_y = (grid.y_edges[:-1] >= yb) & (grid.y_edges[1:] <= yt)
        forced = full_y[:, None] & full_x[None, :]
        contrast[forced] = 0.0
        sigma[forced] = 0.0
    keep = min(grid.n_pixels, int(np.floor(max_area/grid.cell_area)))
    while keep*grid.cell_area > max_area:
        keep -= 1
    energy = (contrast**2+(conductivity_scale*sigma)**2).ravel()
    order = np.argsort(-energy, kind='stable')
    retained = np.zeros(grid.n_pixels, dtype=bool)
    retained[order[:keep]] = True
    contrast.ravel()[~retained] = 0.0
    sigma.ravel()[~retained] = 0.0
    return grid.pack(contrast, sigma)


def occupied_area(parameters, grid=None, *, calibration=True):
    grid = PixelGrid() if grid is None else grid
    contrast, sigma = grid.unpack(parameters)
    active = (contrast != 0.0) | (sigma != 0.0)
    return float(grid.available_areas(calibration=calibration)[active].sum())


def pixel_material(parameters, freq, grid=None, *, calibration=True, max_area=3.0):
    """Return a pointwise coefficient callable after explicit bounds/area checks.

    Inputs are copied and never silently projected.  Call project_parameters
    explicitly when a support/bounds projection is intended.
    """
    grid = PixelGrid() if grid is None else grid
    contrast, sigma = grid.unpack(parameters)
    if not np.isfinite(freq) or freq <= 0:
        raise ValueError('freq must be positive and finite')
    if (np.any(contrast < 0) or np.any(contrast > 8)
            or np.any(sigma < 0) or np.any(sigma > 0.15)):
        raise ValueError('Requires 0<=contrast<=8 and 0<=conductivity<=0.15')
    if not np.isfinite(max_area) or not 0 <= max_area <= 12:
        raise ValueError('max_area must be finite and between zero and twelve')
    active = (contrast != 0.0) | (sigma != 0.0)
    if float(grid.available_areas(calibration=calibration)[active].sum()) > max_area:
        raise ValueError('Pixel material exceeds maximum occupied area')
    values = 1.0+contrast.ravel()+1j*sigma.ravel()/float(freq)

    def epsilon(points):
        points = np.asarray(points, dtype=float)
        ids = grid.indices(points)
        result = np.ones(len(ids), dtype=complex)
        valid = ids >= 0
        result[valid] = values[ids[valid]]
        if calibration:
            result[grid.calibration_mask(points)] = 1.0
        return result

    return epsilon


def _exterior_nodes(left, right, step):
    count = max(1, int(np.ceil((right-left)/step)))
    return np.linspace(left, right, count+1)


def pixel_mesh(grid=None, pml=None, *, exterior_h=0.18, order=3,
               calibration_interfaces=True, max_dofs=60000,
               max_factor_bytes=1_500_000_000):
    """Tensor triangles conforming to every pixel edge and PML interface.

    Optional calibration boundary lines are included by default.  The order is
    used only for the exact pre-mesh DOF/memory gate; pass the same order into
    FEM.  This function neither constructs a FEM basis nor performs a solve.
    """
    grid = PixelGrid() if grid is None else grid
    pml = PMLConfig() if pml is None else pml
    if not np.isfinite(exterior_h) or exterior_h <= 0:
        raise ValueError('exterior_h must be positive and finite')
    if order not in (2, 3) or isinstance(order, bool):
        raise ValueError('order must be two or three')
    ol, or_, ob, ot = pml.outer_box
    xl, xr, yb, yt = pml.physical_box
    x = np.unique(np.r_[grid.x_edges, _exterior_nodes(ol, xl, exterior_h),
                        _exterior_nodes(xl, -2.0, exterior_h),
                        _exterior_nodes(2.0, xr, exterior_h),
                        _exterior_nodes(xr, or_, exterior_h)])
    y = np.unique(np.r_[grid.y_edges, _exterior_nodes(ob, yb, exterior_h),
                        _exterior_nodes(yb, 0.0, exterior_h),
                        _exterior_nodes(3.0, yt, exterior_h),
                        _exterior_nodes(yt, ot, exterior_h)])
    if calibration_interfaces:
        ql, qr, qb, qt = grid.calibration_box
        x = np.unique(np.r_[x, ql, qr])
        y = np.unique(np.r_[y, qb, qt])
    nx, ny = len(x), len(y)
    vertices = nx*ny
    triangles = 2*(nx-1)*(ny-1)
    edges = (nx-1)*ny+(ny-1)*nx+(nx-1)*(ny-1)
    ndofs = vertices+(order-1)*edges+(triangles if order == 3 else 0)
    estimate = _factor_memory_estimate(ndofs)
    if ndofs > max_dofs or estimate > max_factor_bytes:
        raise MemoryError(f'Pixel mesh gate: order={order}, DOFs={ndofs}, '
                          f'estimated factor bytes={estimate}; '
                          f'limits {max_dofs} DOFs/{max_factor_bytes} bytes')
    return MeshTri.init_tensor(x, y)


def volume_pixel_jacobian(fem, state, grid=None, *, calibration=True):
    """Exact pixel Jacobian of the discrete volume-receiver field.

    Return complex (source*receiver,2*nx*ny), source-major, with all real
    contrast parameters followed by all real conductivity parameters.  The
    conductivity block is i/f times the contrast block.  Exactly n_receivers
    transpose RHS are reserved; the state's original LU is reused.
    """
    grid = PixelGrid() if grid is None else grid
    if state._backend_token is not fem._backend_token:
        raise ValueError('State belongs to a different FEM backend')
    ids = grid.indices(state.material_points)
    if np.any(ids < 0):
        raise ValueError('Material quadrature points must all be inside D')
    nr = len(fem.rx)
    ns = state.field.shape[1]
    fem.book.reserve(fem.phase, 'fem_pixel_receiver_transpose', rhs=nr,
                     factorizations=0, freq=fem.freq,
                     backend='independent_fem_pml', receiver_mode='volume',
                     factor_reused=True, pixel_nx=grid.nx, pixel_ny=grid.ny,
                     calibration_air=bool(calibration))
    rhs = np.empty((len(fem._free_dofs), nr), dtype=complex, order='F')
    for r in range(nr):
        load = (fem.k**2*(state.epsilon-1.0)
                * fem._receiver_green[r]).reshape(fem._material_shape)
        rhs[:, r] = asm(_forcing, fem.material_basis, load=load)[fem._free_dofs]
    if not np.isfinite(rhs).all():
        raise RuntimeError('Nonfinite pixel receiver-transpose RHS')
    started = time.perf_counter()
    receive = state.factor.solve(rhs, trans='T')
    seconds = time.perf_counter()-started
    if not np.isfinite(receive).all():
        raise RuntimeError('Nonfinite pixel receiver-transpose solve')
    residuals = (np.linalg.norm(state._matrix.T@receive-rhs, axis=0)
                 / np.maximum(np.linalg.norm(rhs, axis=0), np.finfo(float).tiny))
    residual = float(residuals.max())
    fem._event('backend_solve', label='pixel_receiver_transpose', freq=fem.freq,
               rhs=nr, factor_reused=True, transpose='T', residual=residual,
               wall_seconds=seconds, receiver_mode='volume')
    if not np.isfinite(residual) or residual > fem.residual_tolerance:
        raise RuntimeError(f'Pixel receiver-transpose residual {residual} exceeds tolerance')
    dofs = np.zeros((fem.ndofs, nr), dtype=complex)
    dofs[fem._free_dofs] = receive
    receiving_field = fem._receiver_green.T+np.column_stack([
        np.asarray(fem.material_basis.interpolate(dofs[:, r])).reshape(-1)
        for r in range(nr)])
    weights = fem.k**2*state.material_weights.copy()
    if calibration:
        weights[grid.calibration_mask(state.material_points)] = 0.0
    derivative = np.empty((ns*nr, grid.n_pixels), dtype=complex)
    for s in range(ns):
        source_weight = weights*state.field[:, s]
        for r in range(nr):
            integrand = source_weight*receiving_field[:, r]
            derivative[s*nr+r] = (
                np.bincount(ids, weights=integrand.real, minlength=grid.n_pixels)
                + 1j*np.bincount(ids, weights=integrand.imag, minlength=grid.n_pixels))
    return np.concatenate([derivative, (1j/fem.freq)*derivative], axis=1)


def realify_data(data, *, noise_std=1.0):
    """Flatten source-major, then all Re/all Im; optionally isotropically whiten."""
    if not np.isfinite(noise_std) or noise_std <= 0:
        raise ValueError('noise_std must be positive and finite')
    a = np.asarray(data).ravel(order='C')
    if not np.isfinite(a).all():
        raise ValueError('data must be finite')
    return np.r_[a.real, a.imag]/noise_std


def realify_jacobian(jacobian, *, noise_std=1.0):
    """Stack all Re rows/all Im rows, matching realify_data exactly."""
    if not np.isfinite(noise_std) or noise_std <= 0:
        raise ValueError('noise_std must be positive and finite')
    a = np.asarray(jacobian)
    if a.ndim != 2 or not np.isfinite(a).all():
        raise ValueError('jacobian must be a finite complex matrix')
    return np.concatenate([a.real, a.imag], axis=0)/noise_std
