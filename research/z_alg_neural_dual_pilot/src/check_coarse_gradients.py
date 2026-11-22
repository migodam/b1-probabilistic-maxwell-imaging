"""Check analytic derivatives on the complete 32² grid using artificial commands."""

import argparse
from collections import Counter
from pathlib import Path
import json
import platform
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from chart import MU_SCALE, profile, pullback
from coarse import Coarse


HERE = Path(__file__).resolve().parent


def reference(coarse, command, cap, mu, lam):
    """Use an independent Hermitian eigensolve instead of the proposal LDL."""
    blocks = coarse.block_ids[len(mu)]
    z = (cap * mu[blocks] + 1j * lam[blocks]) / 2
    S, G, B = coarse.S, coarse.G, coarse.B
    diagonal_green = z[:, None] * G
    K = S.conj().T @ S + np.diag(mu[blocks]) - diagonal_green - diagonal_green.conj().T
    spectrum, modes = eigh(K, check_finite=False, driver="evd")
    nu = int((spectrum < 0).sum())
    edge = float(spectrum[coarse.s - 1])
    if nu >= coarse.s:
        return 0.0, nu, edge, float(np.min(np.abs(spectrum)))

    L = S.conj().T @ command + z[:, None] * B
    current = modes @ ((modes.conj().T @ L) / spectrum[:, None])
    H = command.conj().T @ command - L.conj().T @ current
    eigenvalues = eigh((H + H.conj().T) / 2, eigvals_only=True)
    bound2 = float(np.clip(eigenvalues[::-1][nu:], 0.0, None).sum())
    return bound2, nu, edge, float(np.min(np.abs(spectrum)))


def verify_case(coarse, rng, Y, name, blocks, cap, lam, smooth, expected_nu):
    # Small spatial variation keeps the derivative test away from symmetry ties.
    a = rng.normal(0, 0.016, size=blocks)
    b = np.full(blocks, (lam / MU_SCALE - 0.5) / 2) + rng.normal(0, 0.009, size=blocks)
    mu, lam = profile(a, b)
    base = coarse.eval(Y, cap, mu, lam, smooth=smooth)
    reference_start = perf_counter(), process_time()
    bound2, nu, edge, state_gap = reference(coarse, Y, cap, mu, lam)
    reference_wall = perf_counter() - reference_start[0]
    reference_cpu = process_time() - reference_start[1]
    assert expected_nu(nu), (name, nu)
    assert base["nu"] == nu
    assert np.isclose(base["d2_exact"], bound2, atol=1e-9, rtol=1e-8), (base, bound2)
    if nu < coarse.s and not smooth:
        assert base["active"] > 0, name
    if nu >= coarse.s:
        assert base["d2"] < 0 and base["d2_exact"] == 0
        assert np.isclose(base["rescue_K_eigenvalue"], edge, atol=1e-11, rtol=1e-8)

    step = 2e-5
    directional = []
    evaluations = [base]
    for variable in ("mu_only", "lambda_only", "joint"):
        d_mu = rng.normal(size=blocks) if variable != "lambda_only" else np.zeros(blocks)
        d_lam = rng.normal(size=blocks) if variable != "mu_only" else np.zeros(blocks)
        scale = np.max(np.abs(np.r_[d_mu, d_lam]))
        d_mu, d_lam = d_mu / scale * MU_SCALE, d_lam / scale * MU_SCALE
        plus = coarse.eval(Y, cap, mu + step * d_mu, lam + step * d_lam, smooth=smooth)
        minus = coarse.eval(Y, cap, mu - step * d_mu, lam - step * d_lam, smooth=smooth)
        evaluations.extend((plus, minus))
        approx = (plus["d2"] - minus["d2"]) / (2 * step)
        analytic = float(base["grad_mu"] @ d_mu + base["grad_lam"] @ d_lam)
        assert base["nu"] == plus["nu"] == minus["nu"]
        assert base["active"] == plus["active"] == minus["active"]
        assert np.isclose(analytic, approx, atol=4e-8, rtol=5e-5), (name, variable, analytic, approx)
        directional.append(dict(
            variable=variable, analytic=analytic, centered_finite_difference=approx,
            abs_error=abs(analytic - approx),
            rel_error=abs(analytic - approx) / max(abs(approx), abs(analytic), 1e-12),
        ))

    da, db = rng.normal(size=(2, blocks)) / np.sqrt(blocks)
    plus_mu, plus_lam = profile(a + step * da, b + step * db)
    minus_mu, minus_lam = profile(a - step * da, b - step * db)
    plus = coarse.eval(Y, cap, plus_mu, plus_lam, smooth=smooth)
    minus = coarse.eval(Y, cap, minus_mu, minus_lam, smooth=smooth)
    evaluations.extend((plus, minus))
    approx = (plus["d2"] - minus["d2"]) / (2 * step)
    ga, gb = pullback(a, base["grad_mu"], base["grad_lam"])
    analytic = float(ga @ da + gb @ db)
    assert base["nu"] == plus["nu"] == minus["nu"]
    assert base["active"] == plus["active"] == minus["active"]
    assert np.isclose(analytic, approx, atol=4e-8, rtol=5e-5), (name, analytic, approx)
    directional.append(dict(
        variable="chart", analytic=analytic, centered_finite_difference=approx,
        abs_error=abs(analytic - approx),
        rel_error=abs(analytic - approx) / max(abs(approx), abs(analytic), 1e-12),
    ))

    operations = Counter()
    for row in evaluations:
        operations.update(row["operations"])
    return dict(
        case=name, cap=cap, blocks=blocks, smooth=smooth, score=base["d2"],
        hard_d2=base["d2_exact"], reference_d2=bound2, nu=nu, active=base["active"],
        minimum_abs_K_eigenvalue=state_gap, source_cutoff_gap=base["source_cutoff_gap"],
        source_zero_gap=base["source_zero_gap"],
        rescue_K_eigenvalue=base.get("rescue_K_eigenvalue"), rescue_K_gap=base.get("rescue_K_gap"),
        directions=directional, dimensionless_centered_step=step, reference_wall_s=reference_wall,
        reference_cpu_s=reference_cpu, operations=dict(operations),
        proposal_wall_s=[row["wall_s"] for row in evaluations],
        proposal_cpu_s=[row["cpu_s"] for row in evaluations],
    )


def verify_optimizer(coarse, Y):
    direct = coarse.optimize(Y, 1.0, blocks=64)
    warm = profile(np.zeros(64), np.zeros(64))
    hybrid = coarse.refine(Y, 1.0, warm, steps=3)
    for limit, prefix in direct["prefixes"].items():
        eligible = [r for r in direct["trace"] if r["elapsed_wall_s"] <= limit]
        assert len(eligible) == prefix["eval_count"]
        assert prefix["best"] in eligible and prefix["best"]["elapsed_wall_s"] <= limit
        assert prefix["best"]["d2"] == max(row["d2"] for row in eligible if row["valid"])
        if not prefix["best"]["d2"]:
            assert prefix["best"]["guide"] == max(row["guide"] for row in eligible if row["valid"])
    assert direct["prefixes"][1.0]["best"]["d2"] <= direct["prefixes"][5.0]["best"]["d2"]
    assert hybrid["eval_count_total"] == 4 and hybrid["trials_total"] == 3
    assert hybrid["best"]["guide"] >= hybrid["trace"][0]["guide"]
    return dict(
        same_dummy_seed=True, cap=1.0, blocks=64,
        direct={key: value for key, value in direct.items() if key != "trace"},
        hybrid={key: value for key, value in hybrid.items() if key != "trace"},
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE.parent / "results" / "coarse_gradient_check.json")
    args = parser.parse_args()
    start = perf_counter(), process_time()
    seed = 20261123
    rng = np.random.default_rng(seed)
    coarse = Coarse()
    assert np.all(np.bincount(coarse.block_ids[16]) == 64)
    assert np.all(np.bincount(coarse.block_ids[64]) == 16)
    Y = .15 * (rng.normal(size=(20, 10)) + 1j * rng.normal(size=(20, 10)))
    with threadpool_limits(limits=1, user_api="blas"):
        cases = [
            verify_case(coarse, rng, Y, "positive-definite hard", 16, .5, .0005, False,
                        lambda nu: nu == 0),
            verify_case(coarse, rng, Y, "indefinite hard", 64, .75, .01, False,
                        lambda nu: 0 < nu < coarse.s),
            verify_case(coarse, rng, Y, "indefinite softened", 16, .75, .01, True,
                        lambda nu: 0 < nu < coarse.s),
            verify_case(coarse, rng, Y, "exhausted-source rescue", 64, 1.0, .0005, True,
                        lambda nu: nu >= coarse.s),
        ]

    doc = dict(
        status="pass", seed=seed, data="deterministic Gaussian complex dummy, no observation files",
        n=coarse.B.shape[0], n_receivers=coarse.r, n_sources=coarse.s,
        python=platform.python_version(), host=platform.platform(),
        setup_cost=coarse.setup_cost, total_wall_s=perf_counter() - start[0],
        total_cpu_s=process_time() - start[1], cases=cases,
    )
    doc["optimizer_prefix_check"] = verify_optimizer(coarse, Y)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(doc, indent=2, allow_nan=False) + "\n")
    for case in cases:
        print(case["case"], "nu", case["nu"], "active", case["active"],
              "hard", case["hard_d2"], "relative-error", max(d["rel_error"] for d in case["directions"]),
              "proposal median ms", np.median(case["proposal_wall_s"]) * 1e3)


if __name__ == "__main__":
    main()
