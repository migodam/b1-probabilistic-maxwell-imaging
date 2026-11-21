"""One fixed 36-disk patterned host, using the outgoing R2 scalar-TM pulse model."""

from pathlib import Path
from time import perf_counter, process_time
import json

import numpy as np
from scipy.optimize import brentq
from scipy.special import hankel1, jv


RESULTS = Path(__file__).resolve().parents[1] / "results"
RADIUS = "0.04"
ACTIVE_CENTERS = (
    ("-0.29", "-0.15"), ("0.04", "-0.23"), ("0.33", "-0.08"),
    ("-0.22", "0.22"), ("0.09", "0.13"), ("0.38", "0.27"),
)
SOURCE_CENTERS = (("-0.61", "-0.72"), ("0.63", "-0.57"), ("-0.67", "0.42"))
RECEIVER_CENTERS = (("0.17", "0.86"),)
SIDE_X = ("-0.53", "-0.43", "0.49", "0.59")
SIDE_Y = ("-0.30", "-0.15", "0.00", "0.15", "0.30")
RAIL_X = ("-0.32", "-0.17", "-0.02", "0.13", "0.28")
RAIL_Y = ("-0.43", "0.43")
HOST_CENTERS = tuple((x, y) for x in SIDE_X for y in SIDE_Y) + tuple(
    (x, y) for y in RAIL_Y for x in RAIL_X
)
CENTERS = ACTIVE_CENTERS + HOST_CENTERS
WEAK_CHI, STRONG_CHI = "1.4", "2.6"
HOST_CHI = ("0",) * 6 + tuple(
    WEAK_CHI if (row + col) % 2 == 0 else STRONG_CHI
    for col in range(4) for row in range(5)
) + tuple(
    WEAK_CHI if (row + col) % 2 == 0 else STRONG_CHI
    for row in range(2) for col in range(5)
)
T0_START = np.array([
    19.00223026630951, 8.224625770037294, 4.266969667820209,
    2.556933065012014, 19.298009208931692, 25.966536948968322,
])
MATERIAL_BOX = (0.5, 30.0)
TAU = 1e-6


def distances(x, y):
    return np.linalg.norm(x[:, None, :] - y[None, :, :], axis=2)


def geometry():
    """Return the binary64 approximations to the independently specified model."""
    centers = np.asarray(CENTERS, dtype=float)
    sources = np.asarray(SOURCE_CENTERS, dtype=float)
    receivers = np.asarray(RECEIVER_CENTERS, dtype=float)
    return centers, sources, receivers, np.asarray(HOST_CHI, dtype=float)


def pulse_operators():
    centers, sources, receivers, _ = geometry()
    k, radius = 2 * np.pi, float(RADIUS)
    dc = distances(centers, centers)
    offdiag = ~np.eye(len(centers), dtype=bool)
    ds, dr = distances(centers, sources), distances(receivers, centers)
    if min(dc[offdiag]) <= 2 * radius or min(ds.ravel()) <= radius or min(dr.ravel()) <= radius:
        raise ValueError("Overlapping disks or transducer intersects a disk.")
    factor = 0.5j * np.pi * k * radius * jv(1, k * radius)
    G = np.empty_like(dc, dtype=complex)
    G[offdiag] = factor * hankel1(0, k * dc[offdiag])
    np.fill_diagonal(G, 0.5j * np.pi * k * radius * hankel1(1, k * radius) - 1)
    E0 = 0.25j * hankel1(0, k * ds)
    S = factor * hankel1(0, k * dr)
    return G, E0, S


def reduce_host(G, E0, S, m0):
    """Reduce the known-host solve to six additive coordinate susceptibilities."""
    n, r = len(m0), len(ACTIVE_CENTERS)
    A0 = np.eye(n) - m0[:, None] * G
    rhs = np.column_stack((m0[:, None] * E0, np.eye(n, r)))
    X = np.linalg.solve(A0, rhs)
    j0, V = X[:, :E0.shape[1]], X[:, E0.shape[1]:]
    return (G @ V)[:r], (E0 + G @ j0)[:r], S @ V, (S @ j0).ravel(), j0, V


def pack(data):
    """Pack source-major complex values or Jacobian rows in Re, Im order."""
    data = np.asarray(data)
    return np.stack((data.real, data.imag), axis=1).reshape((2 * len(data), *data.shape[1:]))


def reduced_forward(z, H, B, R, y0, derivatives=False):
    A = np.eye(len(z)) - z[:, None] * H
    w = np.linalg.solve(A, z[:, None] * B)
    y = y0 + (R @ w).ravel()
    if not derivatives:
        return y
    C = np.linalg.inv(A)
    field = B + H @ w
    J = pack(field.T * (R @ C))
    return y, J, C, field


def full_forward(z, G, E0, S, m0, derivatives=False):
    """Solve the n-state current equations directly, without host reduction."""
    m = m0.copy()
    m[:len(z)] = z
    A = np.eye(len(m)) - m[:, None] * G
    rhs = m[:, None] * E0
    j = np.linalg.solve(A, rhs)
    y = (S @ j).ravel()
    d = {
        "state_relative_residual": float(np.linalg.norm(A @ j - rhs) / np.linalg.norm(rhs)),
        "sigma_min": float(np.linalg.svd(A, compute_uv=False)[-1]),
    }
    if derivatives:
        T = np.linalg.solve(A.T, S.T).T
        E = E0 + G @ j
        d["J"] = pack(E[:len(z)].T * T[:, :len(z)])
    return y, j, d


def refine_fold(H, B, R, y0, start=T0_START):
    """Try coordinate lines around the public six-disk seed, deterministic order."""
    seed = np.array(start, copy=True)
    search = []

    def determinant(z):
        return np.linalg.det(reduced_forward(z, H, B, R, y0, derivatives=True)[1])

    for axis in (1, 3, 2, 0, 4, 5):
        axis_vector = np.eye(len(seed))[axis]
        lo, hi = MATERIAL_BOX[0] - seed[axis], MATERIAL_BOX[1] - seed[axis]
        offsets = np.unique(np.clip(np.array([-8., -4., -2., -1., -.5, 0., .5, 1., 2., 4., 8.]), lo, hi))
        endpoint_values = [float(determinant(seed + delta * axis_vector)) for delta in offsets]
        changes = [
            [float(a), float(b)] for a, b, da, db in
            zip(offsets[:-1], offsets[1:], endpoint_values[:-1], endpoint_values[1:])
            if da * db < 0
        ]
        search.append({"axis_0based": axis, "offsets": offsets.tolist(),
                       "determinants": endpoint_values, "sign_change_brackets": changes})
        if changes:
            a, b = min(changes, key=lambda bracket: abs(sum(bracket)))
            shift = brentq(lambda delta: determinant(seed + delta * axis_vector), a, b, xtol=6e-16)
            return seed + shift * axis_vector, axis, search
    raise ArithmeticError(f"No sign change on surveyed lines: {search}")


def numerical_artifacts():
    start_wall, start_cpu = perf_counter(), process_time()
    centers, sources, receivers, m0 = geometry()
    G, E0, S = pulse_operators()
    H, B, R, y0, j0, V = reduce_host(G, E0, S, m0)
    zc, axis, search = refine_fold(H, B, R, y0)
    yc, J, C, E = reduced_forward(zc, H, B, R, y0, derivatives=True)
    U, sigma, Vh = np.linalg.svd(J)
    normal, null = U[:, -1], Vh[-1]
    normal *= np.sign(normal[np.argmax(abs(normal))])
    null *= np.sign(null[np.argmax(abs(null))])
    h = pack((2 * R @ C @ np.diag(null) @ H @ C @ np.diag(null) @ E).ravel())
    yfull, jfull, full_metrics = full_forward(zc, G, E0, S, m0, derivatives=True)
    Jfull = full_metrics.pop("J")
    y_nohost, _, nohost = full_forward(zc, G, E0, S, 0 * m0)
    yhost, jhost, host = full_forward(np.empty(0), G, E0, S, m0)
    w = np.linalg.solve(np.eye(len(zc)) - zc[:, None] * H, zc[:, None] * B)
    normal_complex = normal[0::2] + 1j * normal[1::2]
    y_minus, y_plus = yc - TAU * normal_complex, yc + TAU * normal_complex

    RESULTS.mkdir(exist_ok=True)
    np.savez(
        RESULTS / "host_operators.npz", H=H, B=B, R=R, y0=y0,
        zc=zc, axis=axis, normal=normal, null=null, y_critical=yc, y_minus=y_minus,
        y_plus=y_plus, tau=TAU, centers=centers, sources=sources,
        receivers=receivers, radius=float(RADIUS), k=2 * np.pi, G=G, E0=E0, S=S,
        m0=m0, active=np.arange(len(zc)), j0=j0, V=V, material_box=MATERIAL_BOX,
    )
    cell_separation = distances(centers, centers)[np.triu_indices(len(centers), k=1)]
    receipt = {
        "scope": "Finite 2D scalar-TM disk-pulse collocation model, exact outgoing kernel; distances in background wavelengths.",
        "wave_number": "2*pi",
        "time_convention": "e^(-i omega t)",
        "radius": RADIUS,
        "centers": CENTERS,
        "sources": SOURCE_CENTERS,
        "receivers": RECEIVER_CENTERS,
        "active_indices_0based": list(range(len(zc))),
        "known_susceptibilities_and_zero_at_active_sites": HOST_CHI,
        "susceptibility_definition": "epsilon_r - 1",
        "known_material_eps_r": [1 + float(c) for c in HOST_CHI[len(zc):]],
        "active_susceptibility_box": MATERIAL_BOX,
        "minimum_cell_center_separation": float(min(cell_separation)),
        "minimum_disk_edge_separation": float(min(cell_separation) - 2*float(RADIUS)),
        "minimum_source_to_center_separation": float(distances(centers, sources).min()),
        "minimum_receiver_to_center_separation": float(distances(centers, receivers).min()),
        "fold_search": {"source": "R6 T0 historical six-disk start", "seed": T0_START.tolist(),
                        "selected_axis_0based": axis, "survey": search,
                        "fold_center": zc.tolist(), "normal": normal.tolist(),
                        "right_null": null.tolist(), "response_J_singular_values": sigma.tolist(),
                        "unit_vectors_fold_coefficient": float(normal @ h)},
        "full_model_at_fold": full_metrics,
        "separate_direct_full_check": {"response_l2_difference": float(np.linalg.norm(yfull - yc)),
                                       "current_l2_difference": float(np.linalg.norm(jfull - j0 - V @ w)),
                                       "active_jacobian_frobenius_difference": float(np.linalg.norm(Jfull - J))},
        "paired_controls": {"host_only": {**host, "host_only_current_frobenius": float(np.linalg.norm(jhost)),
                                          "host_only_response_l2": float(np.linalg.norm(yhost)),
                                          "host_only_response_vs_y0_l2_difference": float(np.linalg.norm(yhost-y0))},
                            "surrounding_host_removed_at_same_z": {**nohost,
                                "full_vs_host_removed_response_l2_change": float(np.linalg.norm(yfull - y_nohost))}},
        "observation_order": "source-major, real then imaginary",
        "side_offset": TAU,
        "wall_seconds": perf_counter() - start_wall,
        "cpu_seconds": process_time() - start_cpu,
    }
    (RESULTS / "host_model.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"fold_search": receipt["fold_search"], "full_model": full_metrics,
                      "validation": receipt["separate_direct_full_check"], "controls": receipt["paired_controls"]},
                     indent=2))


if __name__ == "__main__":
    numerical_artifacts()
