"""Public target-independent P3 mesh and cut-polygon target assembly.

No FEM operator, factorization, or PDE solve is constructed here.  The caller
owns accounting and qualification.  ``wslam`` must resolve to the frozen copy
under this campaign's ``vendor/`` directory.

The moving disk is the *inscribed regular 384-gon*, not an exact circular
domain.  Each intersection with an existing FEM triangle is triangulated and
integrated independently of the background Basis quadrature.  P3 mass products
are integrated exactly up to floating-point arithmetic when degree >= 6;
analytic incident forcing still needs an independent quadrature check.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.sparse import coo_matrix, csc_matrix
from scipy.spatial import cKDTree
from skfem import ElementTriP3, MeshTri

from wslam.fem_backend import MATERIAL_BOX, PMLConfig


TARGET_BOX = (-1.2, 1.2, 1.4, 2.8)
WALL_BOX = (-1.5, 1.5, 0.5, 1.25)
TARGET_POLYGON_N = 384
MESH_VERSION = "BGGA_PUBLIC_PSLG_V1"


def _positive(value, name):
    if isinstance(value, (bool, np.bool_)) or not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return float(value)


def _settings(h, grid_refinement, target_step, pml):
    h = _positive(h, "h")
    target_step = _positive(target_step, "target_step")
    if isinstance(grid_refinement, (bool, np.bool_)) or grid_refinement not in (1, 2, 3):
        raise ValueError("grid_refinement must be 1, 2, or 3")
    pm = PMLConfig() if pml is None else pml
    if not isinstance(pm, PMLConfig):
        raise TypeError("pml must be the frozen wslam PMLConfig")
    return h, int(grid_refinement), target_step, pm


def public_mesh_key(*, h=.17, grid_refinement=2, target_step=.125, pml=None):
    """Immutable cache key with only public geometry/discretization settings."""
    h, refinement, step, pm = _settings(h, grid_refinement, target_step, pml)
    return (MESH_VERSION, h, refinement, step, TARGET_BOX, WALL_BOX,
            tuple(pm.physical_box), tuple(pm.thickness), tuple(pm.strength), pm.power)


def _uniform_axis(low, high, step):
    return np.linspace(low, high, max(1, int(np.ceil((high-low)/step)))+1)


def public_fixed_mesh(h=.17, grid_refinement=2, target_step=.125, pml=None):
    """One public PSLG mesh for every allowed target pose/radius.

    Preserve every public wall macro interface and add a rectangular grid over
    TARGET_BOX.  Only atomic grid segments are emitted: no long rectangle edge
    overlaps a collection of shorter collinear segments.  No target center,
    radius, material, observation, occupancy, or circle edge enters this mesh.
    ``h`` is Triangle's area-control label, not a certified diameter bound.
    """
    import triangle

    h, refinement, step, pm = _settings(h, grid_refinement, target_step, pml)
    vertices, segments, lookup, edges = [], [], {}, set()

    def vertex(point):
        key = tuple(float(x) for x in point)
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

    def grid(x, y):
        for yi in y:
            for a, b in zip(x[:-1], x[1:]):
                edge((a, yi), (b, yi))
        for xi in x:
            for a, b in zip(y[:-1], y[1:]):
                edge((xi, a), (xi, b))

    rectangle(pm.outer_box)
    rectangle(pm.physical_box)
    rectangle(MATERIAL_BOX)
    grid(np.linspace(WALL_BOX[0], WALL_BOX[1], 24*refinement+1),
         np.linspace(WALL_BOX[2], WALL_BOX[3], 6*refinement+1))
    grid(_uniform_axis(TARGET_BOX[0], TARGET_BOX[1], step),
         _uniform_axis(TARGET_BOX[2], TARGET_BOX[3], step))
    result = triangle.triangulate(
        dict(vertices=np.asarray(vertices, float), segments=np.asarray(segments, np.int32)),
        f"pq28a{h*h/2:.12g}Q")
    if "triangles" not in result:
        raise RuntimeError("public PSLG triangulation failed")
    mesh = MeshTri(result["vertices"].T, result["triangles"].T)
    mesh.bgga_public_key = public_mesh_key(
        h=h, grid_refinement=refinement, target_step=step, pml=pm)
    return mesh


def _cross(a, b):
    return a[..., 0]*b[..., 1]-a[..., 1]*b[..., 0]


def _area(polygon):
    if len(polygon) < 3:
        return 0.
    shifted = polygon-polygon[0]
    return .5*float(np.sum(_cross(shifted, np.roll(shifted, -1, axis=0))))


def _clean_polygon(polygon):
    """Remove roundoff duplicate vertices, never move an interface."""
    polygon = np.asarray(polygon, float).reshape(-1, 2)
    if len(polygon) < 3:
        return np.empty((0, 2))
    tol = 16*np.finfo(float).eps*max(1., float(np.abs(polygon).max()))
    keep = np.linalg.norm(polygon-np.roll(polygon, 1, axis=0), axis=1) > tol
    polygon = polygon[keep]
    return polygon if len(polygon) >= 3 else np.empty((0, 2))


def _clip_convex(subject, clipper):
    """Sutherland-Hodgman intersection; both polygon arguments are CCW."""
    polygon = np.asarray(subject, float)
    for a, b in zip(clipper, np.roll(clipper, -1, axis=0)):
        if len(polygon) < 3:
            return np.empty((0, 2))
        distance = _cross(b-a, polygon-a)
        previous = np.roll(polygon, 1, axis=0)
        dp = np.roll(distance, 1)
        inside = distance >= 0.
        crossing = inside != np.roll(inside, 1)
        counts = inside.astype(np.int64)+crossing.astype(np.int64)
        slots = np.cumsum(counts)-counts
        output = np.empty((int(counts.sum()), 2))
        if np.any(crossing):
            fraction = dp[crossing]/(dp[crossing]-distance[crossing])
            output[slots[crossing]] = previous[crossing]+fraction[:, None]*(
                polygon[crossing]-previous[crossing])
        output[slots[inside]+crossing[inside].astype(np.int64)] = polygon[inside]
        polygon = _clean_polygon(output)
    return polygon


def _ccw(triangle):
    tri = np.asarray(triangle, float)
    return tri if _cross(tri[1]-tri[0], tri[2]-tri[0]) > 0 else tri[[0, 2, 1]]


def independent_triangle_rule(degree=8):
    """Positive Duffy Gauss rule, coordinates (q,2), weights summing to 1/2.

    On the reference triangle (0,0),(1,0),(0,1), x=u, y=(1-u)v,
    and dxdy=(1-u)dudv.  q=ceil((degree+2)/2) Gauss points per axis
    integrate total-degree ``degree`` polynomials.  This uses no Basis.dx or
    background quadrature nodes.
    """
    if isinstance(degree, (bool, np.bool_)) or int(degree) != degree or degree < 6:
        raise ValueError("degree must be an integer >=6 for a P3 mass")
    n = (int(degree)+3)//2
    t, w = np.polynomial.legendre.leggauss(n)
    t, w = .5*(t+1), .5*w
    u, v = np.meshgrid(t, t, indexing="ij")
    weights = (w[:, None]*w[None, :]*(1-u)).ravel()
    return np.column_stack((u.ravel(), ((1-u)*v).ravel())), weights


@dataclass(frozen=True)
class TargetQuadrature:
    points: np.ndarray
    weights: np.ndarray
    elements: np.ndarray
    touched_elements: np.ndarray
    area: float
    polygon_area: float
    polygon_vertices: np.ndarray
    reference_points: np.ndarray
    reference_weights: np.ndarray
    degree: int
    candidate_count: int


@dataclass(frozen=True)
class TargetAssembly:
    Vtarget: csc_matrix
    bX: np.ndarray
    target_area: float
    touched_global_dofs: np.ndarray
    touched_free_indices: np.ndarray
    public_roi_global_dofs: np.ndarray
    public_roi_free_indices: np.ndarray
    touched_rank: int
    Vtarget_rank: int
    quadrature: TargetQuadrature
    metadata: dict


class TargetAssembler:
    """Reusable local triangle index and sparse P3 target assembly.

    ``free_dofs`` is the strictly increasing array of global FEM indices.
    Every returned ``*_free_indices`` indexes this array, hence the free-space
    matrix/RHS.  ROI ports contain the whole public target box, including all
    basis supports intersecting its interior, before a target is supplied.

    The mesh is scanned once to build the ROI index.  A cKDTree query and local
    triangle bounding boxes locate each target; no fine-grid quadrature point
    scan or dense global probe matrix is used.  Point/work-array limits bound
    owned allocations; they do not certify process RSS or Basis storage.
    """
    def __init__(self, basis, free_dofs, *, degree=8, chunk_points=4096,
                 max_points=250_000, max_work_bytes=128_000_000):
        if not isinstance(basis.elem, ElementTriP3):
            raise TypeError("TargetAssembler requires a scalar ElementTriP3 Basis")
        if basis.mesh.p.shape[0] != 2 or basis.mesh.t.shape[0] != 3:
            raise ValueError("only affine two-dimensional triangles are supported")
        free = np.asarray(free_dofs)
        if (free.ndim != 1 or not np.issubdtype(free.dtype, np.integer)
                or np.any(free < 0) or np.any(free >= basis.N)
                or np.any(np.diff(free) <= 0)):
            raise ValueError("free_dofs must be unique sorted global integer indices")
        for value, name in ((chunk_points, "chunk_points"), (max_points, "max_points"),
                            (max_work_bytes, "max_work_bytes")):
            if isinstance(value, bool) or int(value) != value or value <= 0:
                raise ValueError(f"{name} must be a positive integer")
        self.basis, self.free_dofs = basis, free.astype(np.int64, copy=True)
        self.degree = int(degree)
        self.reference_points, self.reference_weights = independent_triangle_rule(degree)
        self.chunk_points, self.max_points = int(chunk_points), int(max_points)
        self.max_work_bytes = int(max_work_bytes)
        self._to_free = np.full(basis.N, -1, dtype=np.int64)
        self._to_free[self.free_dofs] = np.arange(len(free))
        tri = basis.mesh.p[:, basis.mesh.t].transpose(2, 1, 0)
        lo, hi = tri.min(axis=1), tri.max(axis=1)
        xl, xr, yb, yt = TARGET_BOX
        possible = np.flatnonzero((hi[:, 0] > xl) & (lo[:, 0] < xr)
                                  & (hi[:, 1] > yb) & (lo[:, 1] < yt))
        box_polygon = np.array(((xl, yb), (xr, yb), (xr, yt), (xl, yt)))
        roi = [e for e in possible if _area(_clip_convex(_ccw(tri[e]), box_polygon)) > 0]
        self.roi_elements = np.asarray(roi, dtype=np.int64)
        if not len(roi):
            raise ValueError("Basis has no triangles in public target ROI")
        self._triangles = tri[self.roi_elements].copy()
        self._lower, self._upper = lo[self.roi_elements], hi[self.roi_elements]
        centers = self._triangles.mean(axis=1)
        self._max_radius = float(np.linalg.norm(self._triangles-centers[:, None], axis=2).max())
        self._tree = cKDTree(centers)
        global_roi = np.unique(basis.element_dofs[:, self.roi_elements])
        self.public_roi_global_dofs = global_roi[self._to_free[global_roi] >= 0]
        self.public_roi_free_indices = self._to_free[self.public_roi_global_dofs]

    def quadrature(self, pose, radius):
        """Clip the fixed 384-gon against candidate triangles; no PDE action."""
        center = np.asarray(pose, float)
        radius = _positive(radius, "radius")
        if (center.shape != (2,) or not np.isfinite(center).all()
                or not -1 <= center[0] <= 1 or not 1.6 <= center[1] <= 2.6
                or not .06 <= radius <= .16):
            raise ValueError("target requires pose in [-1,1]x[1.6,2.6], radius in [.06,.16]")
        angle = 2*np.pi*np.arange(TARGET_POLYGON_N)/TARGET_POLYGON_N
        disk = center+radius*np.column_stack((np.cos(angle), np.sin(angle)))
        near = np.asarray(sorted(self._tree.query_ball_point(center, radius+self._max_radius)),
                          dtype=np.int64)
        if len(near):
            overlap = np.all(self._upper[near] > center-radius, axis=1) & np.all(
                self._lower[near] < center+radius, axis=1)
            near = near[overlap]
        pieces, touched = [], []
        count = 0
        incircle = radius*np.cos(np.pi/TARGET_POLYGON_N)
        for local in near:
            tri = _ccw(self._triangles[local])
            polygon = (tri if np.all(np.linalg.norm(tri-center, axis=1) <= incircle)
                       else _clip_convex(disk, tri))
            if _area(polygon) <= 0:
                continue
            element = int(self.roi_elements[local])
            found = False
            for j in range(1, len(polygon)-1):
                a, b, c = polygon[0], polygon[j], polygon[j+1]
                determinant = float(_cross(b-a, c-a))
                if determinant <= 0:
                    continue
                count += len(self.reference_weights)
                if count > self.max_points or count*32 > self.max_work_bytes:
                    raise MemoryError("target cut-polygon quadrature point/work-array limit exceeded")
                pieces.append((a.copy(), b.copy(), c.copy(), determinant, element))
                found = True
            if found:
                touched.append(element)
        if not pieces:
            raise RuntimeError("empty target polygon/mesh intersection")
        points = np.empty((count, 2))
        weights = np.empty(count)
        elements = np.empty(count, dtype=np.int64)
        nq = len(self.reference_weights)
        for j, (a, b, c, determinant, element) in enumerate(pieces):
            sl = slice(j*nq, (j+1)*nq)
            points[sl] = a+self.reference_points[:, :1]*(b-a)+self.reference_points[:, 1:]*(c-a)
            weights[sl] = determinant*self.reference_weights
            elements[sl] = element
        exact_polygon_area = .5*TARGET_POLYGON_N*radius**2*np.sin(2*np.pi/TARGET_POLYGON_N)
        return TargetQuadrature(points, weights, elements, np.asarray(touched, dtype=np.int64),
                                float(weights.sum()), float(exact_polygon_area), disk,
                                self.reference_points.copy(), self.reference_weights.copy(),
                                self.degree, len(near))

    def _known_element_probes(self, points, elements, global_support):
        """Basis.probes' sparse evaluation with already-known parent cells.

        The pinned skfem Basis.probes first invokes element_finder, whose
        fallback forms an all-elements by all-query-points temporary array.
        Cut integration already supplies each parent element, so use the
        same mapping/gbasis operations without that locator or its fallback.
        Every temporary here is proportional to the bounded query chunk.
        """
        x = points.T
        local = self.basis.mapping.invF(x[:, :, None], tind=elements)
        values = np.asarray([
            self.basis.elem.gbasis(self.basis.mapping, local, j, tind=elements)[0]
            for j in range(self.basis.Nbfun)]).ravel()
        rows = np.tile(np.arange(len(points), dtype=np.int32), self.basis.Nbfun)
        columns = self.basis.element_dofs[:, elements].ravel()
        return coo_matrix((values, (rows, columns)),
                          shape=(len(points), self.basis.N)).tocsr()[:, global_support]

    def assemble(self, pose, radius, k, chi, incident: Callable, *, source_count=None):
        """Return V=k² chi int(phi phi), bX=k² chi int(phi u_inc).

        ``incident(points Nx2)`` must directly evaluate the analytic incident
        field and return shape (N, sources), or (N,) for one source.  Its values
        are evaluated at the cut quadrature nodes; no nodal interpolation is
        requested or synthesized.  No complex conjugation is used in either
        weak-form integral (the P3 basis itself is real).
        """
        k = _positive(k, "k")
        chi = complex(chi)
        if not np.isfinite(chi) or not callable(incident):
            raise ValueError("chi must be finite and incident must be callable")
        if source_count is not None and (isinstance(source_count, bool)
                or int(source_count) != source_count or source_count < 1):
            raise ValueError("source_count must be a positive integer when supplied")
        quadrature = self.quadrature(pose, radius)
        global_support = np.unique(self.basis.element_dofs[:, quadrature.touched_elements])
        global_support = global_support[self._to_free[global_support] >= 0]
        support = self._to_free[global_support]
        size, nfree = len(support), len(self.free_dofs)
        mass = csc_matrix((size, size), dtype=float)
        force = None
        ns = int(source_count) if source_count is not None else None
        output_estimate = 0
        for start in range(0, len(quadrature.weights), self.chunk_points):
            end = min(start+self.chunk_points, len(quadrature.weights))
            points, weights = quadrature.points[start:end], quadrature.weights[start:end]
            values = np.asarray(incident(points.copy()), dtype=complex)
            if values.ndim == 1:
                values = values[:, None]
            if (values.ndim != 2 or values.shape[0] != len(points) or not values.shape[1]
                    or not np.isfinite(values).all()):
                raise ValueError("analytic incident must return finite (point,source) values")
            if ns is None:
                ns = values.shape[1]
            if values.shape[1] != ns:
                raise ValueError("analytic incident source count changed or disagrees with source_count")
            # Includes retained quadrature, output RHS, two local sparse mass
            # matrices at a conservative dense bound, and chunk work arrays.
            output_estimate = (quadrature.points.nbytes+quadrature.weights.nbytes
                               + quadrature.elements.nbytes+16*nfree*ns
                               + 32*size*size+len(points)*(64*ns+512))
            if output_estimate > self.max_work_bytes:
                raise MemoryError("target assembly owned-array estimate exceeds max_work_bytes")
            if force is None:
                force = np.zeros((size, ns), dtype=complex)
            probe = self._known_element_probes(points, quadrature.elements[start:end], global_support)
            mass = mass+(probe.T@probe.multiply(weights[:, None])).tocsc()
            force += np.asarray(probe.T@(weights[:, None]*values))
        mass.sum_duplicates()
        mass.eliminate_zeros()
        local = mass.tocoo()
        scale = k*k*chi
        Vtarget = coo_matrix((scale*local.data, (support[local.row], support[local.col])),
                             shape=(nfree, nfree), dtype=complex).tocsc()
        Vtarget.eliminate_zeros()
        bX = np.zeros((nfree, ns), dtype=complex)
        bX[support] = scale*force
        # On any positive-area part of a triangle, a P3 polynomial vanishing
        # identically has all ten local coefficients zero.  Thus the exact
        # unscaled Gram matrix has rank equal to this *actual touched* support.
        # This is not a numerical SVD, conditioning, or small-cut qualification.
        touched_rank = len(support)
        metadata = dict(
            target_polygon_n=TARGET_POLYGON_N, target_area=quadrature.area,
            polygon_area=quadrature.polygon_area,
            polygon_partition_area_error=quadrature.area-quadrature.polygon_area,
            circle_area=np.pi*radius**2,
            polygon_relative_circle_area_deficit=1-np.sin(2*np.pi/TARGET_POLYGON_N)/(2*np.pi/TARGET_POLYGON_N),
            quadrature_degree=self.degree, quadrature_points=len(quadrature.weights),
            candidate_triangles=quadrature.candidate_count,
            touched_triangles=len(quadrature.touched_elements),
            touched_ports=touched_rank, public_roi_ports=len(self.public_roi_free_indices),
            rank_meaning="exact-algebraic rank on positive-area touched P3 supports; numerical rank NOT_RUN",
            forcing="direct analytic incident at independent cut-polygon quadrature points",
            locator="one-time triangle ROI index, cKDTree query, local triangle bbox filter",
            sparse_probe="Basis.probes mapping/gbasis formula using known parent cells; no global finder fallback",
            owned_array_estimate_bytes=int(output_estimate), max_points=self.max_points,
            max_work_bytes=self.max_work_bytes, chunk_points=self.chunk_points,
            process_RSS_qualification="NOT_RUN", mass_qualification="NOT_RUN",
            incident_quadrature_qualification="NOT_RUN", PDE_qualification="NOT_RUN")
        return TargetAssembly(Vtarget, bX, quadrature.area, global_support, support,
                              self.public_roi_global_dofs.copy(), self.public_roi_free_indices.copy(),
                              touched_rank, touched_rank if chi != 0 else 0, quadrature, metadata)


def public_roi_ports(basis, free_dofs):
    """Complete fixed target-ROI port indices in the caller's free-DOF order."""
    return TargetAssembler(basis, free_dofs).public_roi_free_indices.copy()


def assemble_target(basis, free_dofs, pose, radius, k, chi, incident, *,
                    assembler=None, source_count=None, **settings):
    """Convenience wrapper; reuse ``TargetAssembler`` across target candidates."""
    if assembler is None:
        assembler = TargetAssembler(basis, free_dofs, **settings)
    elif settings or assembler.basis is not basis or not np.array_equal(assembler.free_dofs, free_dofs):
        raise ValueError("reused assembler must match Basis/free_dofs and owns its settings")
    return assembler.assemble(pose, radius, k, chi, incident, source_count=source_count)
