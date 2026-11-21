"""Arb proof of a selected indefinite sector bound on the frozen finite model."""
import argparse
import json
from pathlib import Path

from scan import HERE, OBSERVATION_RADII, adjoint, clock, hermitian, load_case, read_json, timing

import numpy as np
from scipy.linalg import LinAlgError, cholesky, eigh, svd
from flint import acb, acb_mat, arb, ctx


def ball(array):
    n, m = array.shape
    return acb_mat(n, m, [acb(float(x.real), float(x.imag)) for x in array.flat])


def star(a):
    return a.conjugate().transpose()


def frobenius_squared(a):
    return sum((v.real * v.real + v.imag * v.imag for v in a.entries()), arb(0))


def gershgorin(a, *, upper):
    sign = 1 if upper else -1
    endpoints = []
    for i in range(a.nrows()):
        off = sum((abs(a[i, j]) for j in range(a.ncols()) if i != j), arb(0))
        edge = a[i, i].real + sign * off
        endpoints.append(edge.upper() if upper else edge.lower())
    return max(endpoints) if upper else min(endpoints)


def safe_float(a, *, upper):
    v = a.upper() if upper else a.lower()
    return float(np.nextafter(float(v), np.inf if upper else -np.inf))


def certify_dark_obstruction(G_np, S_np, G, S, u):
    """Exhibit a dark current with negative sector form, even for ill-conditioned S."""
    started = clock()
    _, sigma, vh = svd(S_np, full_matrices=True, check_finite=False)
    approximate_kernel = adjoint(vh[S_np.shape[0]:, :])
    restricted_green = hermitian(adjoint(approximate_kernel) @ hermitian(G_np) @ approximate_kernel)
    last = restricted_green.shape[0] - 1
    _, vector = eigh(restricted_green, subset_by_index=(last, last), check_finite=False)
    proposal = ball(approximate_kernel @ vector)
    old_precision = ctx.prec
    ctx.prec = max(256, old_precision)
    try:
        # The exact vector enclosed below lies in ker(S). Arb must first
        # establish invertibility of the exact receiver Gram to return a solve.
        St = star(S)
        exact_dark = proposal - St * (S * St).solve(S * proposal)
        a = frobenius_squared(exact_dark)
        b = (star(exact_dark) * (G * exact_dark))[0, 0].real
        gap_upper = (a - arb(str(u)) * b).upper()
        return {
            "status": "Arb_verified_receiver_gram",
            "precision_bits": ctx.prec,
            "receiver_sigma_min_numeric": float(sigma[-1]),
            "receiver_sigma_max_numeric": float(sigma[0]),
            "exact_dark_vector_sector_value_upper": safe_float(gap_upper, upper=True),
            "dark_obstruction_threshold_upper": safe_float(a / b, upper=True) if b > 0 else None,
            "positive_definite_K_impossible_for_every_mu_positive": bool(gap_upper < 0),
            "receiver_gram_dimension": S_np.shape[0],
            **timing(started),
        }
    except (ZeroDivisionError, ValueError) as error:
        return {"status": "Arb_dark_projection_unresolved", "reason": str(error),
                "precision_bits": ctx.prec, **timing(started)}
    finally:
        ctx.prec = old_precision


def prove(case, u, mu, precision, threads):
    ctx.prec, ctx.threads = precision, threads
    started = clock()
    G, S, B, Y, radii = load_case(case, read_json(OBSERVATION_RADII))
    n, s = B.shape
    i_np = np.eye(n)
    start = clock()
    K_np = hermitian(adjoint(S) @ S + mu * (i_np - u * hermitian(G)))
    L_np = adjoint(S) @ Y + (.5 * mu * u) * B
    d, V = eigh(K_np, check_finite=False, driver="evd")
    nu = int(np.count_nonzero(d < 0))
    if not 0 < nu < s:
        raise ValueError(f"Inertia proposal has nu={nu} and {s} sources.")
    target_gap = float(np.min(abs(d)) / 2)
    Vneg = V[:, :nu]
    W = Vneg * np.sqrt(-2 * d[:nu])
    F = cholesky(K_np + W @ adjoint(W) - target_gap * i_np, lower=True, check_finite=False)
    X_np = V @ ((adjoint(V) @ L_np) / d[:, None])
    h_eigenvalues, hv = eigh(hermitian(adjoint(Y) @ Y - adjoint(L_np) @ X_np), check_finite=False)
    h_eigenvalues, hv = h_eigenvalues[::-1], hv[:, ::-1]
    double = timing(start)

    start = clock()
    G_np, S_np = G, S
    G, S, B, Y, I = map(ball, (G, S, B, Y, i_np))
    N, W, F, X = map(ball, (Vneg, W, F, X_np))
    a, b, shift = acb(arb(str(u))), acb(arb(str(mu))), acb(target_gap)
    K = star(S) * S + b * (I - a * (G + star(G)) / 2)
    L = star(S) * Y + (a * b / 2) * B

    # This is the only cubic full-dimensional Arb matrix product.
    residual = K + W * star(W) - shift * I - F * star(F)
    positive_residual = frobenius_squared(residual).sqrt().upper()
    alpha = (shift.real - positive_residual).lower()
    if not alpha > 0:
        raise ArithmeticError("Shifted Cholesky residual cannot certify the positive eigenvalue gap.")
    negative_upper = gershgorin(star(N) * (K * N + shift * N), upper=True)
    if not negative_upper < 0:
        raise ArithmeticError("Projected negative block cannot certify the negative eigenvalue gap.")
    inertia_timing = timing(start)

    start = clock()
    KX = K * X
    defect = frobenius_squared(KX - L).upper()
    correction = (defect / alpha).upper()
    H0 = star(Y) * Y - star(L) * X - star(X) * L + star(X) * KX
    terms = []
    trace_floor = arb(0)
    for j in range(nu + 1, s + 1):
        if h_eigenvalues[j - 1] <= 0:
            break
        Vj = ball(hv[:, :j])
        numerator_floor = gershgorin(star(Vj) * H0 * Vj, upper=False)
        denominator_ceiling = gershgorin(star(Vj) * Vj, upper=True)
        if numerator_floor > 0 and denominator_ceiling > 0:
            term = max(arb(0), (numerator_floor / denominator_ceiling - correction).lower())
        else:
            term = arb(0)
        trace_floor += term
        terms.append({"eigenvalue_index_descending_one_based": j,
                      "eigenvalue_lower_bound": safe_float(term, upper=False) if term > 0 else 0.})
    lower = max(arb(0), trace_floor.lower())
    evaluated = timing(start)
    distance_floor = max(0., safe_float(lower.sqrt(), upper=False))
    dark_control = certify_dark_obstruction(G_np, S_np, G, S, u)
    return {
        "case": case, "u": u, "mu": mu, "status": "Arb_verified",
        "precision_bits": precision, "Arb_threads": threads,
        "state_dimension": n, "source_rhs": s, "verified_negative_eigenvalues": nu,
        "proposed_dyadic_shift": target_gap,
        "positive_shifted_cholesky_defect_Frobenius_upper": safe_float(positive_residual, upper=True),
        "negative_subspace_shifted_gershgorin_upper": safe_float(negative_upper, upper=True),
        "certified_K_absolute_spectral_gap_lower": safe_float(alpha, upper=False),
        "KX_minus_L_Frobenius_upper": safe_float(defect.sqrt(), upper=True),
        "H_operator_correction_upper": safe_float(correction, upper=True),
        "retained_H_eigenvalue_lower_bounds": terms,
        "certified_frobenius_residual_lower": distance_floor,
        "radii": radii,
        "strictly_excludes": {name: bool(lower > arb(radius)**2) for name, radius in radii.items()},
        "positive_definite_obstruction": dark_control,
        "cost": {"double_proposal": double, "Arb_inertia": inertia_timing, "Arb_H": evaluated,
                 "full_dimension_cubic_Arb_products": 1, "double_K_eigendecompositions": 1,
                 "double_shifted_cholesky": 1, "double_H_eigendecompositions": 1, **timing(started)},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan", type=Path, default=HERE / "results.json")
    parser.add_argument("--output", type=Path, default=HERE / "certificate.json")
    parser.add_argument("--case")
    parser.add_argument("--cap", type=float)
    parser.add_argument("--mu", type=float, help="certify the given case and cap directly, with no scan input")
    parser.add_argument("--precision", type=int, default=96)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if args.mu is not None:
        if args.case is None or args.cap is None:
            parser.error("--mu requires both --case and --cap")
        selected = (args.case, args.cap, args.mu)
    else:
        found = []
        for case, data in read_json(args.scan)["cases"].items():
            if args.case and args.case != case:
                continue
            for cap in data["grid"]:
                if args.cap is not None and args.cap != cap["u"]:
                    continue
                if any(cap["new_indefinite_only_exclusion_numeric"].values()):
                    row = cap["best_indefinite_with_source_directions"]
                    strength = max(row["margin_numeric"].values())
                    found.append((data["state_dimension"], -strength, case, cap["u"], row["mu"]))
        selected = min(found)[2:] if found else None
    if selected is None:
        result = {"status": "no_new_indefinite_only_exclusion_on_grid"}
    else:
        case, u, mu = selected
        start = clock()
        try:
            result = prove(case, u, mu, args.precision, args.threads)
        except (ArithmeticError, ValueError, LinAlgError) as error:
            result = {"status": "Arb_unresolved", "case": case, "u": u, "mu": mu,
                      "reason": str(error), **timing(start)}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
