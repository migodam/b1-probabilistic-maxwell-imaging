"""Sufficient task separation on a full-dimensional manufacturing-tolerance box.

Run after array/array_secant.py. Arb encloses the archived discrete arrays.
The trial real data decoders are chosen by floating-point linear algebra.
All claims for each fixed decoder are evaluated again with real/complex balls.
"""
import os
for var in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[var] = "1"

import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from flint import acb, acb_mat, arb, arb_mat, ctx
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
ctx.prec = 128


def ball_array(x):
    return np.array([[acb(float(c.real), float(c.imag)) for c in row] for row in x],
                    dtype=object)


def from_mat(x):
    return np.array([[x[i, j] for j in range(x.ncols())] for i in range(x.nrows())],
                    dtype=object)


def eye(n):
    return np.array([[arb(int(i == j)) for j in range(n)] for i in range(n)], dtype=object)


def abs_upper(x):
    return np.frompyfunc(lambda t: abs(t).upper(), 1, 1)(x)


def up(x):
    return float(np.nextafter(float(x.upper()), np.inf))


def positive_resolvent(k):
    # A positive w with (I-K)w>0 proves rho(K)<1 for nonnegative K.
    matrix = arb_mat((eye(len(k)) - k).tolist())
    try:
        w = from_mat(matrix.solve(arb_mat([[1]] * len(k))))
        if not all(v > 0 for v in w.ravel()):
            return None
        inverse = from_mat(matrix.inv())
        return abs_upper(inverse)
    except (ValueError, ZeroDivisionError):
        return None


def dressed(raw, size, center, band):
    stem = f"disks{size:02d}__"
    def field(key):
        return ball_array(raw[stem + band + "_" + key])
    g, b, s = field("G"), field("B"), field("S")
    mat = np.concatenate((center, raw[stem + "known_host_chi"]))
    a = acb_mat((eye(len(g)) - mat[:, None] * g).tolist())
    ai = from_mat(a.inv())
    p = (s @ ai)[:, :size]
    h = (g @ ai)[:size, :size]
    e = (b + g @ ai @ (mat[:, None] * b))[:size]
    return p, h, e


def channel_block(dressing, radii, indices, scale):
    p, h, e = dressing
    # dP <= |P|D |H|(I-D|H|)^-1; dE <= (I-|H|D)^-1 |H|D |E|.
    rp = radii[:, None] * abs_upper(h)
    re = abs_upper(h) * radii[None, :]
    qp, qe = positive_resolvent(rp), positive_resolvent(re)
    if qp is None or qe is None:
        return None
    p1 = abs_upper(p) @ rp
    e1 = re @ abs_upper(e)
    dp = p1 @ qp
    de = qe @ e1
    p2 = p1 @ rp @ qp
    e2 = qe @ re @ e1
    jac, hessian, tail = [], [], []
    for a, b in indices:
        j = scale * p[a] * e[:, b]
        first = scale * (p[a][None, :] * h.T * e[:, b][:, None] +
                         p[a][:, None] * h * e[:, b][None, :])
        rest = scale * (p2[a] * abs_upper(e[:, b]) + abs_upper(p[a]) * e2[:, b] +
                        dp[a] * de[:, b])
        jac.append(j)
        hessian.append(first)
        tail.append(rest)
    return np.array(jac), np.array(hessian), np.array(tail)


def realify(a):
    return np.concatenate((np.frompyfunc(lambda v: v.real, 1, 1)(a),
                           np.frompyfunc(lambda v: v.imag, 1, 1)(a)), axis=0)


def real_trial(a):
    return np.vectorize(float, otypes=[float])(a)


def evaluate(raw, case, dress, r, arm):
    q = np.array([arb(int(i in case["roi_indices"])) / len(case["roi_indices"])
                  for i in range(case["n"])], dtype=object)
    prior = 2 * (q @ r)
    n = case["n"]
    base = [(0, 0), (1, 0), (0, 1)]
    pair_scale = 1 / arb(2).sqrt()
    spec = {
        "base": [(0, base, arb(1))],
        "angle": [(0, base + [(0, 2), (1, 2), (0, 3)], pair_scale)],
        "frequency": [(0, base, pair_scale), (1, base, (arb(6) / 5).sqrt() * pair_scale)],
        "rich": [(0, [(a, b) for b in range(4) for a in range(6)], pair_scale)]
    }
    blocks = [channel_block(dress[band], r, mask, scale) for band, mask, scale in spec[arm]]
    out = {"prior_task_diameter": up(prior), "status": "prior_only"}
    if any(block is None for block in blocks):
        out["reason"] = "whole-box absolute Green majorant did not contract"
        return out
    j, hess, rem = [np.concatenate([block[i] for block in blocks], axis=0) for i in range(3)]
    j = realify(j)
    hs = realify(hess).reshape((len(j), n * n))
    j_float = real_trial(j)
    # Noise is a fixed numerical radius. Its source is the archived base measurement.
    case_meta = next(x for x in PAIRS if x["case"] == case["id"] and x["arm"] == arm)
    delta = arb(case_meta["radius"])
    out["delta"] = float(delta)
    u, sv, vh = np.linalg.svd(j_float, full_matrices=False)
    # Trial vectors chosen independently of endpoint responses and of material truth.
    tol = np.finfo(float).eps * max(j.shape) * sv[0]
    modes = sv > tol
    left = (vh[modes].T / sv[modes]) @ u[:, modes].T
    ridge = float(delta) / max(real_trial(r))
    shrink = sv / (sv ** 2 + ridge ** 2)
    trials = [("zero", np.zeros(len(j))), ("least_squares", real_trial(q) @ left),
              ("ridge", (real_trial(q) @ vh.T * shrink) @ u.T)]

    def bounds_for(weights, target):
        # The same real h_j appears in both derivative legs. Combine its signed
        # contribution across measured ports before taking absolute values.
        w = np.vectorize(arb, otypes=[object])(weights)
        jerr = abs_upper(target - w @ j)
        first = abs_upper((w @ hs).reshape((n, n))) @ r
        wabs = np.array([(w[i] ** 2 + w[i + len(rem)] ** 2).sqrt()
                         for i in range(len(rem))], dtype=object)
        return jerr + first + wabs @ rem, sum(x * x for x in w).sqrt()

    table = []
    for name, w in trials:
        er, gain = bounds_for(w, q)
        bias = 2 * (er @ r)
        table.append({"decoder": name, "bias": up(bias), "gain": up(gain),
                      "task_diameter_at_delta": up(bias + 2 * delta * gain)})
    out.update(status="bounded", dual_trials=table, gamma_resolution_task_units=up(prior/2))
    rho = prior / 2
    # The same derivative enclosure yields a multiplicative estimate when a
    # full left inverse is available. Keep the ordinary interval-inverse control.
    if modes.sum() == n:
        target = eye(n)
        rows = [bounds_for(w, t)[0] for w, t in zip(left, target)]
        k = np.array(rows, dtype=object)
        resolvent = positive_resolvent(k)
        if resolvent is not None:
            b = np.vectorize(arb, otypes=[object])(left)
            q_b = q @ b
            tail_weight = (abs_upper(q) @ k) @ resolvent
            gain = sum(v ** 2 for v in q_b).sqrt() + sum(
                tail_weight[i] * sum(x*x for x in b[i]).sqrt() for i in range(n))
            table.append({"decoder": "uniform_left_inverse", "bias": 0.0, "gain": up(gain),
                          "task_diameter_at_delta": up(2 * delta * gain)})
            out["uniform_task_lipschitz_upper"] = up(gain)
    def gamma(item):
        if not item["gain"] or item["decoder"] == "zero":
            return 0.
        v = (1-arb(item["bias"])/rho)/arb(item["gain"])
        return max(0., float(np.nextafter(float(v.lower()), -np.inf)))
    out["task_gamma_lower_at_half_prior"] = max(map(gamma, table))
    out["best_diameter_upper"] = min(up(prior), *(z["task_diameter_at_delta"] for z in table))
    return out


if __name__ == "__main__":
    start, cpu = perf_counter(), process_time()
    data = json.loads((ROOT / "array/results/summary.json").read_text())
    PAIRS = data["pairs"]
    archive = np.load(ROOT / "array/results/physical_arrays.npz")
    rows = []
    with threadpool_limits(limits=1):
        for case in data["cases"]:
            center = np.asarray(case["material_center"])
            dressing = [dressed(archive, case["n"], center, band) for band in ("f0", "f1")]
            # Fixed absolute manufacturing tolerances in chi, intersected
            # symmetrically with the legal box. Each cell varies independently.
            for width in (0.05, 0.25, 0.75):
                radii = np.array([min(arb(width), arb(float(c))-arb(0.5), arb(30)-arb(float(c)))
                                  for c in center], dtype=object)
                for arm in ("base", "angle", "frequency", "rich"):
                    result = evaluate(archive, case, dressing, radii, arm)
                    rows.append(dict(case=case["id"], width=width, arm=arm, **result))
    metadata = {"precision_bits": ctx.prec, "coefficient_scope":
                "exact archived binary64 full G,S,B, known host/centers; sqrt energy scaling in Arb",
                "noise_scope": "archived absolute radius", "wall_s": perf_counter()-start,
                "cpu_s": process_time()-cpu, "rows": rows}
    (ROOT / "bounds/results.json").write_text(json.dumps(metadata, indent=2, allow_nan=False))
    print(len(rows), metadata["wall_s"])
