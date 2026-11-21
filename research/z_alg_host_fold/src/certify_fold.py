"""Arb certificates for the exact-geometry R2 pulse model on this 36-disk host."""

from time import perf_counter, process_time
import json

from flint import acb, acb_mat, arb, arb_mat, ctx
import numpy as np

from host_model import (
    ACTIVE_CENTERS, CENTERS, HOST_CHI, MATERIAL_BOX, RADIUS, RECEIVER_CENTERS, RESULTS,
    SOURCE_CENTERS,
)


PRECISION = 192
HALFWIDTH = "1e-11"
SIGMA_LOWER = "0.26124"


def identity(n):
    return acb_mat([[int(i == j) for j in range(n)] for i in range(n)])


def diagonal(v):
    return acb_mat([[v[i] if i == j else 0 for j in range(len(v))] for i in range(len(v))])


def exact_dyadic_decimal(x):
    """Serialize an exact Arb endpoint without inward decimal rounding."""
    significand, exponent = (int(v) for v in x.man_exp())
    if exponent >= 0:
        return str(significand << exponent)
    return f"{significand * 5**(-exponent)}e{exponent}"


def interval(x):
    return {"lo": exact_dyadic_decimal(x.lower()), "hi": exact_dyadic_decimal(x.upper())}


def complex_interval(x):
    return {"real": interval(x.real), "imag": interval(x.imag)}


def matrix_intervals(A):
    return [[complex_interval(A[i, j]) for j in range(A.ncols())] for i in range(A.nrows())]


def float_interval(x):
    """Readable outer float interval (the proof itself uses Arb)."""
    return {
        "lo": float(np.nextafter(float(x.lower()), -np.inf)),
        "hi": float(np.nextafter(float(x.upper()), np.inf)),
    }


def distance(x, y):
    return norm(a - b for a, b in zip(x, y))


def hankel(order, x):
    return acb(x).bessel_j(order) + acb(0, 1) * acb(x).bessel_y(order)


def pulse_operators():
    """Enclose the decimal geometry, 2π, sqrt distances and the Hankel functions."""
    centers = [[arb(a) for a in point] for point in CENTERS]
    sources = [[arb(a) for a in point] for point in SOURCE_CENTERS]
    receivers = [[arb(a) for a in point] for point in RECEIVER_CENTERS]
    radius, pi = arb(RADIUS), arb.pi()
    k = 2 * pi
    i = acb(0, 1)
    ka = k * radius
    disk_factor = i * pi * ka * acb(ka).bessel_j(1) / 2
    self_pulse = i * pi * ka * hankel(1, ka) / 2 - 1
    pair_distances = [
        (row, col, distance(x, y)) for row, x in enumerate(centers)
        for col, y in enumerate(centers[:row])
    ]
    if not all(d > 2 * radius for _, _, d in pair_distances):
        raise ArithmeticError("Disk disjointness not established.")
    G = acb_mat(len(centers), len(centers))
    for row in range(len(centers)):
        G[row, row] = self_pulse
    for row, col, d in pair_distances:
        G[row, col] = G[col, row] = disk_factor * hankel(0, k * d)
    E0 = acb_mat([[i * hankel(0, k * distance(x, y)) / 4 for y in sources] for x in centers])
    S = acb_mat([[disk_factor * hankel(0, k * distance(x, y)) for y in centers] for x in receivers])
    axes = list(zip(*centers))
    lower = [min(q.lower() for q in axis) - radius for axis in axes]
    upper = [max(q.upper() for q in axis) + radius for axis in axes]
    if not all(any(lo > q or q > hi for q, lo, hi in zip(sensor, lower, upper))
               for sensor in sources + receivers):
        raise ArithmeticError("A transducer is not outside the array bounding rectangle.")
    return G, E0, S


def reduce_host(G, E0, S):
    """Verified A0 solve with three host incident RHS and six coordinate RHS."""
    n, r, ns = len(CENTERS), len(ACTIVE_CENTERS), E0.ncols()
    m0 = [arb(v) for v in HOST_CHI]
    D0 = diagonal(m0)
    A0 = identity(n) - D0 * G
    X = A0.solve(acb_mat([
        [m0[i] * E0[i, s] for s in range(ns)] + [int(i == j) for j in range(r)]
        for i in range(n)
    ]))
    j0 = acb_mat([[X[i, s] for s in range(ns)] for i in range(n)])
    V = acb_mat([[X[i, ns + j] for j in range(r)] for i in range(n)])
    GV, Ehost = G * V, E0 + G * j0
    H = acb_mat([[GV[i, j] for j in range(r)] for i in range(r)])
    B = acb_mat([[Ehost[i, s] for s in range(ns)] for i in range(r)])
    return H, B, S * V, S * j0, A0, j0, V


def pack_rows(A):
    return arb_mat([
        [part(A[s, j]) for j in range(A.ncols())]
        for s in range(A.nrows()) for part in (lambda q: q.real, lambda q: q.imag)
    ])


def response(z, H, B, R, y0):
    D = diagonal(z)
    C = (identity(len(z)) - D * H).inv()
    w = C * D * B
    field, T = B + H * w, R * C
    J = pack_rows(acb_mat([
        [T[0, j] * field[j, s] for j in range(len(z))] for s in range(B.ncols())
    ]))
    return y0 + R * w, J, C, field, T


def real_column(A):
    return [A[i, 0] for i in range(A.nrows())]


def norm(values):
    lower, upper = arb(0), arb(0)
    for x in values:
        lo, hi = x.abs_lower(), x.abs_upper()
        lower += lo * lo
        upper += hi * hi
    return lower.sqrt().union(upper.sqrt())


def complex_norm(A):
    return norm(A[i, j] for i in range(A.nrows()) for j in range(A.ncols()))


def unit(vector):
    length = norm(vector)
    return [x / length for x in vector]


def null_from_minor(J, row, col):
    """From a nonsingular minor, fix the two kernel anchors at one."""
    other_rows = [i for i in range(J.nrows()) if i != row]
    other_cols = [j for j in range(J.ncols()) if j != col]
    minor = arb_mat([[J[i, j] for j in other_cols] for i in other_rows])
    minor_det = minor.det()
    if minor_det.contains(0):
        raise ArithmeticError("Selected rank-five minor not established.")
    v = real_column(minor.solve(arb_mat([[-J[i, col]] for i in other_rows])))
    ell = real_column(minor.transpose().solve(arb_mat([[-J[row, j]] for j in other_cols])))
    v.insert(col, arb(1))
    ell.insert(row, arb(1))
    return ell, v, minor_det, other_rows, other_cols


def positive_ldl_pivots(A):
    """Interval LDL* on one triangle of a Hermitian matrix."""
    n = A.nrows()
    residual = [[A[i, j] for j in range(i + 1)] for i in range(n)]
    pivots = []
    for j in range(n):
        pivot = residual[j][j].real
        if not pivot > 0:
            raise ArithmeticError(f"LDL* positive pivot not established at row {j}: {pivot}")
        pivots.append(pivot)
        for i in range(j + 1, n):
            scale = residual[i][j] / pivot
            for k in range(j + 1, i + 1):
                residual[i][k] -= scale * residual[k][j].conjugate()
    return pivots


def full_matrix(z, G):
    m = z + [arb(v) for v in HOST_CHI[len(z):]]
    return identity(len(m)) - diagonal(m) * G


def sigma_min_enclosure(A, lower, trial_vector):
    """Positive LDL* of A*A−t²I and a verified Rayleigh quotient upper bound."""
    threshold = arb(lower)
    Astar = A.conjugate().transpose()
    pivots = positive_ldl_pivots(Astar * A - identity(A.ncols()) * threshold**2)
    v = acb_mat([[acb(complex(t))] for t in trial_vector])
    upper = (complex_norm(A * v) / complex_norm(v)).upper()
    return threshold, upper, pivots


def full_at_boxes(H, B, R, y0, G, E0, S, z, numeric):
    A = full_matrix(z, G)
    n = A.nrows()
    material = z + [arb(v) for v in HOST_CHI[len(z):]]
    full_y = S * A.solve(diagonal(material) * E0)
    reduced_y = response(z, H, B, R, y0)[0]
    difference = full_y - reduced_y
    if not all(difference[0, j].contains(0) for j in range(difference.ncols())):
        raise ArithmeticError("Independent full/reduced enclosures are inconsistent.")
    A_center = np.eye(n) - np.asarray([float(v.mid()) for v in material])[:, None] * numeric["G"]
    _, values, Vh = np.linalg.svd(A_center)
    floor, upper, pivots = sigma_min_enclosure(A, SIGMA_LOWER, Vh[-1].conjugate())
    return {
        "n_complex_state_degrees_of_freedom": n,
        "sample_sigma_min": float(values[-1]),
        "certified_sigma_min_strictly_exceeds": SIGMA_LOWER,
        "certified_sigma_min_upper": float_interval(upper)["hi"],
        "requested_threshold_0p25": floor > arb("0.25"),
        "positive_ldl_pivots": [float_interval(x) for x in pivots],
        "direct_full_response_enclosure": matrix_intervals(full_y)[0],
        "reduced_response_enclosure": matrix_intervals(reduced_y)[0],
        "independently_evaluated_full_reduced_l2_difference_upper": float_interval(complex_norm(difference))["hi"],
    }


def certify(numeric, H, B, R, y0, G, E0, S):
    axis = int(numeric["axis"])
    centers = [arb(float(x)) for x in numeric["zc"]]
    h = arb(HALFWIDTH)
    a, b = centers.copy(), centers.copy()
    a[axis], b[axis] = a[axis] - h, b[axis] + h
    X = centers.copy()
    X[axis] = a[axis].union(b[axis])
    if not all(arb(MATERIAL_BOX[0]) < z < arb(MATERIAL_BOX[1]) for z in X):
        raise ArithmeticError("Certified interval not entirely inside the material box.")
    ends = [response(z, H, B, R, y0)[1].det() for z in (a, b)]
    if not (ends[0] < 0 < ends[1] or ends[1] < 0 < ends[0]):
        raise ArithmeticError(f"No opposite determinant signs established: {ends}")

    yc, J, C, field, T = response(X, H, B, R, y0)
    n = len(X)
    row = int(np.argmax(abs(numeric["normal"])))
    col = int(np.argmax(abs(numeric["null"])))
    ell, v, minor_det, other_rows, other_cols = null_from_minor(J, row, col)
    vu, lu = unit(v), unit(ell)

    Dv = diagonal(vu)
    second = real_column(pack_rows((2 * T * Dv * H * C * Dv * field).transpose()))
    coefficient = sum((l * k for l, k in zip(lu, second)), arb(0))
    if coefficient.contains(0):
        raise ArithmeticError("Null Hessian pairing not separated from zero.")

    # The derivative of the scalar Schur complement is ell^T J'(t) v
    # even off the singularity, with the anchors of both vectors set to one.
    Da = diagonal([int(j == axis) for j in range(n)])
    dT, dE = T * Da * H * C, H * C * Da * field
    Jprime = pack_rows(acb_mat([
        [dT[0, j] * field[j, s] + T[0, j] * dE[j, s] for j in range(n)]
        for s in range(B.ncols())
    ]))
    slope = (arb_mat([ell]) * Jprime * arb_mat([[q] for q in v]))[0, 0]
    if slope.contains(0):
        raise ArithmeticError("Uniqueness of the critical point on the line not established.")

    numerical_normal = [arb(float(q)) for q in numeric["normal"]]
    normal_error = norm(x - y for x, y in zip(lu, numerical_normal))
    yc_fixed = acb_mat([[acb(complex(q)) for q in numeric["y_critical"]]])
    data_error = complex_norm(yc - yc_fixed)
    sides = [acb_mat([[acb(complex(q)) for q in numeric[name]]]) for name in ("y_minus", "y_plus")]
    side_diffs = [x - yc for x in sides]
    side_projections = [(arb_mat([lu]) * pack_rows(d.transpose()))[0, 0] for d in side_diffs]
    if not side_projections[0] < 0 < side_projections[1]:
        raise ArithmeticError("Finite measurements not separated along the true response normal.")
    return {
        "claim": "One certified corank-one ordinary fold on the specified coordinate segment.",
        "precision_bits": ctx.prec,
        "coordinate_axis_0based": axis,
        "binary64_center": numeric["zc"].tolist(),
        "exact_segment_halfwidth": HALFWIDTH,
        "critical_material_outer_box": [float_interval(x) for x in X],
        "endpoint_real_determinants": [float_interval(x) for x in ends],
        "minor_determinant": float_interval(minor_det),
        "minor_rows_0based": other_rows,
        "minor_cols_0based": other_cols,
        "unit_right_null_intervals": [float_interval(x) for x in vu],
        "unit_left_null_intervals": [float_interval(x) for x in lu],
        "null_hessian_pairing": float_interval(coefficient),
        "schur_complement_derivative_along_segment": float_interval(slope),
        "critical_response_enclosure": matrix_intervals(yc)[0],
        "critical_response_vs_frozen_y_critical_l2_error_upper": float_interval(data_error)["hi"],
        "unit_left_null_vs_frozen_normal_l2_error_upper": float_interval(normal_error)["hi"],
        "true_normal_projection_of_y_minus_and_y_plus_minus_critical_response": [
            float_interval(x) for x in side_projections
        ],
        "true_critical_value_to_frozen_side_observation_l2_intervals": [
            float_interval(complex_norm(x)) for x in side_diffs
        ],
        "full_physical_state": full_at_boxes(H, B, R, y0, G, E0, S, X, numeric),
    }


def main():
    ctx.prec = PRECISION
    start_wall, start_cpu = perf_counter(), process_time()
    numeric = np.load(RESULTS / "host_operators.npz")
    G, E0, S = pulse_operators()
    H, B, R, y0, A0, _, _ = reduce_host(G, E0, S)
    intervals = dict(
        H=matrix_intervals(H), B=matrix_intervals(B), R=matrix_intervals(R),
        y0=matrix_intervals(y0)[0], precision_bits=ctx.prec,
        scope="Ball arithmetic enclosures of the exact decimal geometry, k=2π, outgoing disk-pulse operators and complex host solve.",
        binary64_approximations_file="host_operators.npz",
    )
    (RESULTS / "host_operators_interval.json").write_text(json.dumps(intervals, indent=2) + "\n")

    # G is complex symmetric, so its Hermitian imaginary part equals entrywise Im(G).
    # If Im(G) is positive definite, D^{-1}-G cannot have a null vector for real
    # positive D. This verifies full current regularity throughout the material box.
    ImG = acb_mat([[G[i, j].imag for j in range(G.ncols())] for i in range(G.nrows())])
    pivots = positive_ldl_pivots(ImG)
    host = {
        "host_solve": "Arb complex verified 36x36 solve, 3 incident plus 6 selector RHS.",
        "host_A0_determinant": complex_interval(A0.det()),
        "strictly_positive_imaginary_green_pivots": [float_interval(d) for d in pivots],
        "all_material_box_full_state_regular": True,
        "proof": "Positive definite Im(G) via interval LDL* and positive diagonal D imply D^{-1}−G has no nonzero null vector; A = D(D^{-1}−G).",
        "wall_seconds": perf_counter()-start_wall,
        "cpu_seconds": process_time()-start_cpu,
    }
    (RESULTS / "host_state_certificate.json").write_text(json.dumps(host, indent=2) + "\n")

    start_wall, start_cpu = perf_counter(), process_time()
    proof = certify(numeric, H, B, R, y0, G, E0, S)
    proof["wall_seconds"] = perf_counter()-start_wall
    proof["cpu_seconds"] = process_time()-start_cpu
    (RESULTS / "fold_certificate.json").write_text(json.dumps(proof, indent=2) + "\n")
    print(json.dumps({
        "critical_point": proof["critical_material_outer_box"],
        "null_hessian_pairing": proof["null_hessian_pairing"],
        "critical_response_l2_error": proof["critical_response_vs_frozen_y_critical_l2_error_upper"],
        "determinant_endpoints": proof["endpoint_real_determinants"],
        "certified_sigma_min_bounds": [proof["full_physical_state"][k] for k in (
            "certified_sigma_min_strictly_exceeds", "certified_sigma_min_upper"
        )],
        "smallest_verified_imaginary_G_pivot_lower": min(d["lo"] for d in host["strictly_positive_imaginary_green_pivots"]),
        "fold_wall_seconds": proof["wall_seconds"],
    }, indent=2))


if __name__ == "__main__":
    main()
