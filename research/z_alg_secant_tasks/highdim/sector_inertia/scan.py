"""Compare the indefinite sector lower bound with its positive-definite restriction."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from replay import OBSERVATION_RADII, SOURCE, load_operators, read_json

import numpy as np
from scipy.linalg import eigh, eigvalsh
from threadpoolctl import threadpool_limits

CAPS = (.1, .5, 1., 2., 3., 4.)
MULTIPLIERS = (1e-4, 1e-3, 1e-2, 1e-1, 1., 10.)


def adjoint(matrix):
    return matrix.conj().T


def hermitian(matrix):
    return (matrix + adjoint(matrix)) * .5


def clock():
    return perf_counter(), process_time()


def timing(start):
    return dict(zip(("wall_seconds", "cpu_seconds"),
                    (float(t - t0) for t, t0 in zip(clock(), start))))


def load_case(case, radii_by_case):
    with np.load(SOURCE / f"{case}.npz", allow_pickle=False) as archive:
        Y = archive["Y_archive"]
    geometry = load_operators(case)
    return geometry["G"], geometry["S"], geometry["E0"], Y, radii_by_case[case]


def row_for(K, L, gram_y, radii, u, mu, cost):
    start = clock()
    eigenvalues, modes = eigh(K, check_finite=False, driver="evd")
    cost["K_hermitian_eigendecompositions"] += 1
    gap = float(np.min(np.abs(eigenvalues)))
    scale = float(max(abs(eigenvalues[0]), abs(eigenvalues[-1])))
    tolerance = max(np.finfo(float).tiny, 128 * K.shape[0] * np.finfo(float).eps * scale)
    nsource = L.shape[1]
    row = {
        "u": u, "mu": mu, "K_min_eigenvalue": float(eigenvalues[0]),
        "K_max_eigenvalue": float(eigenvalues[-1]),
        "K_min_abs_eigenvalue": gap, "K_singularity_tolerance": float(tolerance),
    }
    if gap <= tolerance:
        return {**row, "status": "numeric_inertia_unresolved", **timing(start)}

    nu = int(np.count_nonzero(eigenvalues < 0))
    below_budget = nu < nsource
    row.update(nu=nu, available_source_directions=max(0, nsource - nu))
    if below_budget:
        # Hermitian eigensolve serves both inertia and all source solves.
        UstarL = adjoint(modes) @ L
        X = modes @ (UstarL / eigenvalues[:, None])
        H = hermitian(gram_y - adjoint(L) @ X)
        h = eigvalsh(H, check_finite=False)[::-1]
        h_squared = float(np.maximum(h[nu:], 0.).sum())
        untrimmed_squared = float(np.maximum(h, 0.).sum())
        residual = float(np.linalg.norm(K @ X - L))
        cost["source_rhs"] += nsource
        cost["H_hermitian_eigendecompositions"] += 1
        row.update(H_eigenvalues_descending=h.tolist(), full_positive_H_trace=untrimmed_squared,
                   KX_minus_L_Frobenius=residual)
    else:
        h_squared = 0.
        cost["K_negative_count_exhausts_sources"] += 1
    d = float(np.sqrt(h_squared))
    row.update(
        status=("positive_definite_numeric" if nu == 0 else
                "indefinite_with_source_directions_numeric" if below_budget else
                "negative_directions_exhaust_sources_numeric"),
        lower_bound_squared_numeric=h_squared, lower_bound_numeric=d,
        exceeds_radius_numeric={ball: bool(d > r) for ball, r in radii.items()},
        margin_numeric={ball: float(d - r) for ball, r in radii.items()},
        **timing(start),
    )
    return row


def brief(row):
    if row is None:
        return None
    return {k: row[k] for k in ("mu", "nu", "status", "lower_bound_numeric",
                               "margin_numeric", "exceeds_radius_numeric")}


def summarize_case(case, radii_by_case, caps=CAPS):
    start = clock()
    G, S, B, Y, radii = load_case(case, radii_by_case)
    n = G.shape[0]
    I = np.eye(n)
    Gh = hermitian(G)
    S2 = hermitian(adjoint(S) @ S)
    SY = adjoint(S) @ Y
    Y2 = hermitian(adjoint(Y) @ Y)
    preparation = timing(start)
    cost = Counter(S_adjoint_S_products=1, S_adjoint_Y_products=1, Y_adjoint_Y_products=1)
    grid = []
    for u in caps:
        C = I - u * Gh
        rows = []
        for mu in MULTIPLIERS:
            K, L = S2 + mu * C, SY + (.5 * mu * u) * B
            rows.append(row_for(K, L, Y2, radii, u, mu, cost))
        positive = [r for r in rows if r.get("nu") == 0]
        indefinite = [r for r in rows if r.get("status") == "indefinite_with_source_directions_numeric"]
        best_pos = max(positive, key=lambda r: r["lower_bound_numeric"], default=None)
        best_ind = max(indefinite, key=lambda r: r["lower_bound_numeric"], default=None)
        best_pos_dist = best_pos["lower_bound_numeric"] if best_pos else 0.
        best_ind_dist = best_ind["lower_bound_numeric"] if best_ind else 0.
        grid.append({
            "u": u, "rows": rows, "best_positive_definite": brief(best_pos),
            "best_indefinite_with_source_directions": brief(best_ind),
            "indefinite_positive_margin_gain": max(0., best_ind_dist - best_pos_dist),
            "new_indefinite_only_exclusion_numeric": {
                name: bool(best_pos_dist <= radius < best_ind_dist) for name, radius in radii.items()
            },
        })
    return {
        "case": case, "state_dimension": n, "receiver_source_shape": list(Y.shape),
        "receiver_operator_norm_Frobenius": float(np.linalg.norm(S)),
        "observation_Frobenius": float(np.linalg.norm(Y)), "radii": radii,
        "preparation": preparation, "grid": grid, "operations": dict(cost), **timing(start),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    parser.add_argument("--case", action="append", help="scan only this historical case")
    parser.add_argument("--caps", nargs="+", type=float, default=CAPS)
    args = parser.parse_args()
    start = clock()
    radii_by_case = read_json(OBSERVATION_RADII)
    case_ids = args.case or radii_by_case
    results = {}
    with threadpool_limits(limits=1, user_api="blas"):
        for case in case_ids:
            result = summarize_case(case, radii_by_case, args.caps)
            results[case] = result
            print(case, "n=", result["state_dimension"], "Y=", result["receiver_source_shape"],
                  "gains", [(g["u"], round(g["indefinite_positive_margin_gain"], 5))
                            for g in result["grid"] if g["indefinite_positive_margin_gain"] > 0], flush=True)
    output = {
        "scope": "Original complex observations and binary64 VIE operators. Numerical scan only.",
        "caps": args.caps, "multipliers": MULTIPLIERS, "cases": results, **timing(start),
    }
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
