"""Fit a known mathematical circular probe while preserving every base facet.

This geometry repair uses only the supplied base mesh and circle.  There are
no wall generators, truth labels, coefficient fits, or physical solver calls.
The circle is an inscribed polygon; its approximation error remains explicit.
"""

import time

import numpy as np
from skfem import MeshTri


def _cross(a, b):
    return a[..., 0]*b[..., 1]-a[..., 1]*b[..., 0]


def fit_probe_mesh(base_mesh, center, radius, n=256, *, max_area=0.0162,
                   min_angle=20.0):
    """Return MeshTri conforming to all original facets plus a circle polygon.

    Every base/circle crossing is split before Triangle sees the PSLG.  Shared
    intersection indices and ordered segment splits avoid intersecting input
    constraints.  Triangle may subdivide a preserved facet but cannot remove
    its geometric constraint.  The target interior is meshed, never a hole.
    ``min_angle=None`` omits Triangle's quality refinement; ``max_area=None``
    omits its area refinement.  Neither changes or smooths the input geometry.

    ``mesh.probe_fit_diagnostics`` contains JSON-compatible intersection,
    constraint, mesh-size, and polygon-error notes.  It reports P2/P3 DOF counts;
    the caller must apply its existing DOF/memory gate before constructing FEM.
    No physical budget is consumed here.
    """
    import triangle

    started = time.perf_counter()
    if not isinstance(base_mesh, MeshTri):
        raise TypeError('base_mesh must be a scikit-fem MeshTri')
    if base_mesh.p.shape[0] != 2 or base_mesh.t.shape[0] != 3:
        raise ValueError('base_mesh must have two-dimensional triangular cells')
    base = np.asarray(base_mesh.p.T, dtype=float)
    if not np.isfinite(base).all():
        raise ValueError('base mesh coordinates must be finite')
    c = np.asarray(center, dtype=float)
    if c.shape != (2,) or not np.isfinite(c).all():
        raise ValueError('center must be a finite two-dimensional point')
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError('radius must be positive and finite')
    if isinstance(n, bool) or int(n) != n or n < 8:
        raise ValueError('n must be an integer >=8')
    n = int(n)
    if max_area is not None and (not np.isfinite(max_area) or max_area <= 0):
        raise ValueError('max_area must be positive and finite, or None')
    if min_angle is not None and (not np.isfinite(min_angle) or not 0 < min_angle < 60):
        raise ValueError('min_angle must lie strictly between zero and sixty, or be None')
    if np.any(c-radius <= base.min(axis=0)) or np.any(c+radius >= base.max(axis=0)):
        raise ValueError('circle must lie strictly within the base mesh rectangle')
    facets = np.asarray(base_mesh.facets.T, dtype=np.int64)
    scale = max(1.0, float(np.max(np.abs(base))), float(radius))
    tolerance = 128*np.finfo(float).eps*scale

    # Roundoff-scale spatial bins join independently computed copies of one
    # intersection.  Original coordinates are retained; no grid rounding or
    # physical geometry threshold is applied.
    vertices = []
    bins = {}

    def vertex_id(point):
        point = np.asarray(point, dtype=float)
        key = tuple(int(np.floor(v/tolerance)) for v in point)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for index in bins.get((key[0]+dx, key[1]+dy), ()):
                    if np.linalg.norm(vertices[index]-point) <= tolerance:
                        return index
        index = len(vertices)
        vertices.append(point.copy())
        bins.setdefault(key, []).append(index)
        return index

    base_ids = np.array([vertex_id(p) for p in base], dtype=np.int64)
    merged_base_vertices = len(base)-len(vertices)
    angles = 2*np.pi*np.arange(n)/n
    circle = c+float(radius)*np.c_[np.cos(angles), np.sin(angles)]
    circle_ids = np.array([vertex_id(p) for p in circle], dtype=np.int64)
    base_splits = [[(0.0, int(base_ids[a])), (1.0, int(base_ids[b]))]
                   for a, b in facets]
    circle_splits = [[(0.0, int(circle_ids[j])), (1.0, int(circle_ids[(j+1) % n]))]
                     for j in range(n)]
    pa = base[facets[:, 0]]
    pb = base[facets[:, 1]]
    candidate = np.flatnonzero(np.all(np.maximum(pa, pb) >= c-radius-tolerance, axis=1)
                               & np.all(np.minimum(pa, pb) <= c+radius+tolerance, axis=1))
    p = pa[candidate]
    r = pb[candidate]-p
    q = circle
    s = np.roll(circle, -1, axis=0)-circle
    lr = np.linalg.norm(r, axis=1)
    ls = np.linalg.norm(s, axis=1)
    if np.any(lr == 0) or np.any(ls == 0):
        raise ValueError('Zero-length base or circle segment')
    denominator = _cross(r[:, None, :], s[None, :, :])
    offset = q[None, :, :]-p[:, None, :]
    parallel = (abs(denominator)
                <= 32*np.finfo(float).eps*lr[:, None]*ls[None, :])
    t = np.zeros_like(denominator)
    u = np.zeros_like(denominator)
    np.divide(_cross(offset, s[None, :, :]), denominator, out=t, where=~parallel)
    np.divide(_cross(offset, r[:, None, :]), denominator, out=u, where=~parallel)
    t_tol = tolerance/lr[:, None]
    u_tol = tolerance/ls[None, :]
    hits = (~parallel & (t >= -t_tol) & (t <= 1+t_tol)
            & (u >= -u_tol) & (u <= 1+u_tol))
    intersection_ids = set()
    proper_events = 0
    endpoint_events = 0

    def record(local_base, circle_edge, tb, uc, point):
        index = vertex_id(point)
        intersection_ids.add(index)
        base_splits[int(candidate[local_base])].append((float(np.clip(tb, 0, 1)), index))
        circle_splits[circle_edge].append((float(np.clip(uc, 0, 1)), index))

    for i, j in np.argwhere(hits):
        tb, uc = float(t[i, j]), float(u[i, j])
        # Exact existing endpoint coordinates take priority near an endpoint.
        if abs(tb) <= t_tol[i, 0]:
            point, tb = p[i], 0.0
        elif abs(tb-1) <= t_tol[i, 0]:
            point, tb = p[i]+r[i], 1.0
        elif abs(uc) <= u_tol[0, j]:
            point, uc = q[j], 0.0
        elif abs(uc-1) <= u_tol[0, j]:
            point, uc = q[j]+s[j], 1.0
        else:
            point = 0.5*((p[i]+tb*r[i])+(q[j]+uc*s[j]))
        record(int(i), int(j), tb, uc, point)
        if (t_tol[i, 0] < tb < 1-t_tol[i, 0]
                and u_tol[0, j] < uc < 1-u_tol[0, j]):
            proper_events += 1
        else:
            endpoint_events += 1

    # Collinear overlaps are uncommon, but their endpoints must also become
    # shared splits so the input has no overlapping, unsplit PSLG segments.
    collinear = parallel & (abs(_cross(offset, r[:, None, :]))
                            <= tolerance*lr[:, None])
    collinear_pairs = 0
    for i, j in np.argwhere(collinear):
        rr, ss = float(r[i]@r[i]), float(s[j]@s[j])
        pairs = []
        for uc, point in ((0.0, q[j]), (1.0, q[j]+s[j])):
            tb = float((point-p[i])@r[i]/rr)
            if -t_tol[i, 0] <= tb <= 1+t_tol[i, 0]:
                pairs.append((tb, uc, point))
        for tb, point in ((0.0, p[i]), (1.0, p[i]+r[i])):
            uc = float((point-q[j])@s[j]/ss)
            if -u_tol[0, j] <= uc <= 1+u_tol[0, j]:
                pairs.append((tb, uc, point))
        if pairs:
            collinear_pairs += 1
        for tb, uc, point in pairs:
            record(int(i), int(j), tb, uc, point)

    segments = []
    markers = []
    seen = set()
    base_piece_count = 0
    circle_piece_count = 0
    duplicate_pieces = 0

    def add_pieces(splits, marker):
        nonlocal duplicate_pieces
        ordered = sorted(splits, key=lambda item: item[0])
        indices = []
        for _, index in ordered:
            if not indices or index != indices[-1]:
                indices.append(index)
        count = 0
        for a, b in zip(indices[:-1], indices[1:]):
            if a == b:
                continue
            key = (min(a, b), max(a, b))
            if key in seen:
                duplicate_pieces += 1
                continue
            seen.add(key)
            segments.append((a, b))
            markers.append(marker)
            count += 1
        return count

    for j, splits in enumerate(base_splits):
        base_piece_count += add_pieces(splits, j+1)
    for splits in circle_splits:
        circle_piece_count += add_pieces(splits, len(facets)+1)
    input_vertices = np.asarray(vertices, dtype=float)
    options = ('p'+('' if min_angle is None else f'q{float(min_angle):.12g}')
               + ('' if max_area is None else f'a{float(max_area):.12g}')+'Q')
    result = triangle.triangulate(dict(vertices=input_vertices,
                                       segments=np.asarray(segments, dtype=np.int32),
                                       segment_markers=np.asarray(markers, dtype=np.int32)[:, None]),
                                  options)
    if 'triangles' not in result or not len(result['triangles']):
        raise RuntimeError('Probe-fitting triangulation failed')
    if (len(result['vertices']) < len(input_vertices)
            or not np.array_equal(result['vertices'][:len(input_vertices)], input_vertices)):
        raise RuntimeError('Triangle changed or removed supplied PSLG vertices')
    output_markers = np.asarray(result.get('segment_markers', []), dtype=int).ravel()
    retained = np.unique(output_markers[(output_markers > 0) & (output_markers <= len(facets))])
    if len(retained) != len(facets):
        raise RuntimeError('A preserved base-facet constraint marker is missing')
    mesh = MeshTri(result['vertices'].T, result['triangles'].T)
    cell_points = mesh.p.T[mesh.t.T]
    double_area = abs(_cross(cell_points[:, 1]-cell_points[:, 0],
                             cell_points[:, 2]-cell_points[:, 0]))
    if np.any(double_area == 0):
        raise RuntimeError('Probe-fitting triangulation produced a zero-area cell')
    minimum_angle = np.pi
    maximum_edge_sq = np.zeros(len(cell_points))
    for j in range(3):
        first = cell_points[:, (j+1) % 3]-cell_points[:, j]
        second = cell_points[:, (j+2) % 3]-cell_points[:, j]
        dot = np.einsum('ij,ij->i', first, second)
        minimum_angle = min(minimum_angle, float(np.arctan2(double_area, dot).min()))
        maximum_edge_sq = np.maximum(maximum_edge_sq, np.einsum('ij,ij->i', first, first))
    mesh.probe_fit_diagnostics = dict(
        known_geometry='supplied mathematical circle only; no wall geometry',
        circle_center=c.tolist(), circle_radius=float(radius), circle_polygon_n=n,
        circle_max_sagitta=float(radius*(1-np.cos(np.pi/n))),
        base_vertices=int(base_mesh.nvertices), base_facets=len(facets),
        candidate_base_facets=len(candidate), merged_base_vertices=merged_base_vertices,
        proper_intersection_events=proper_events, endpoint_intersection_events=endpoint_events,
        unique_intersection_vertices=len(intersection_ids), collinear_overlap_pairs=collinear_pairs,
        split_base_segments=base_piece_count, split_circle_segments=circle_piece_count,
        duplicate_segment_pieces_removed=duplicate_pieces,
        input_pslg_vertices=len(input_vertices), input_pslg_segments=len(segments),
        retained_original_facet_markers=len(retained),
        output_vertices=int(mesh.nvertices), output_facets=int(mesh.nfacets),
        output_triangles=int(mesh.nelements),
        p2_dofs=int(mesh.nvertices+mesh.nfacets),
        p3_dofs=int(mesh.nvertices+2*mesh.nfacets+mesh.nelements),
        requested_min_angle=(None if min_angle is None else float(min_angle)),
        requested_max_area=(None if max_area is None else float(max_area)),
        min_triangle_area=float(double_area.min()/2),
        min_triangle_angle_degrees=float(np.degrees(minimum_angle)),
        max_edge_to_min_altitude_ratio=float((maximum_edge_sq/double_area).max()),
        triangle_options=options, coordinate_join_tolerance=float(tolerance),
        circle_interior='meshed medium; no holes supplied',
        guard='caller must apply existing FEM DOF/memory gate before physics',
        wall_seconds=time.perf_counter()-started,
        physical_solves_run=0)
    return mesh
