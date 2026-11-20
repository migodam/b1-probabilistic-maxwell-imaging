"""Measure how the chosen electric-current modes represent a real voxel change."""

import json
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import eigvalsh, lu_solve, qr, svd
from threadpoolctl import threadpool_limits

from cohorts import load_cohort
from current_cycle import (
    OUT,
    gradient_material_direction,
    norm,
    physical,
    prepare,
    tangent_columns,
)


def orthogonal_columns(x):
    return qr(x, mode="economic", check_finite=False)[0]


def capture(q, target):
    return norm(q.conj().T @ target) ** 2 / norm(target) ** 2


def spectrum_and_retention(model, m, data, bank, rng):
    G, V, Q = model["G"], bank["v"], bank["q"]
    current, prediction, factor = physical(model, m)
    field = model["E0"] + G @ current
    power = np.sum(abs(field) ** 2, axis=1)
    light = power > 0.0
    e = field[light] / np.sqrt(power[light, None])

    def propagate(q):
        return orthogonal_columns(q - m[:, None] * (G @ q))

    base = propagate(Q)
    # The projector and the illumination coherence have independent physical origins.
    p = base[light] @ base[light].conj().T
    K = (p * (e.conj() @ e.T)).real
    eigenvalues = eigvalsh(K, check_finite=False)
    d = gradient_material_direction(model, m, current, prediction - data, factor, data)
    t = d[:, None] * field
    c0 = tangent_columns(model, m, field, d, factor)
    q36 = Q[:, :36]
    remainder = c0 - q36 @ (q36.conj().T @ c0)

    def augment(extra):
        u, sigma, _ = svd(extra, full_matrices=False, check_finite=False)
        keep = sigma > max(1.0e-9 * sigma[0], 1.0e-15)
        return orthogonal_columns(np.column_stack((q36, u[:, keep])))

    soft = augment(remainder)
    dark = augment(remainder - V @ (V.conj().T @ remainder))

    forbidden = orthogonal_columns(lu_solve(factor, V, trans=2, check_finite=False))
    best = lu_solve(
        factor, t - forbidden @ (forbidden.conj().T @ t), check_finite=False
    )
    best -= q36 @ (q36.conj().T @ best)
    optimal_dark = augment(best)
    directions = {
        "data_gradient": d,
        "one_voxel": np.eye(1, m.size, m.size // 2 + 2).ravel(),
        "random_real": rng.standard_normal(m.size),
    }
    also = {}
    identity_error = 0.0
    for name, vector in directions.items():
        target = vector[:, None] * field
        if norm(target) == 0.0:
            also[name] = None
            continue
        x = np.sqrt(power[light]) * vector[light]
        direct = capture(base, target)
        identity_error = max(identity_error, abs((x @ K @ x) / (x @ x) - direct))
        also[name] = dict(
            fixed48=direct,
            all_dark_ceiling=1.0 - capture(forbidden, target),
            updated48_dark=capture(propagate(dark), target),
            updated48_soft=capture(propagate(soft), target),
            optimal48_dark=capture(propagate(optimal_dark), target),
        )
    return dict(
        rank=Q.shape[1],
        augmented_rank=soft.shape[1],
        illuminated=int(light.sum()),
        trace=float(np.trace(K)),
        top=float(eigenvalues[-1]),
        min=float(eigenvalues[0]),
        above_half=int(sum(eigenvalues >= 0.5)),
        above_nine_tenths=int(sum(eigenvalues >= 0.9)),
        trace_error=float(abs(np.sum(eigenvalues) - norm(base[light]) ** 2)),
        formula_error=identity_error,
        capture=also,
        eigenvalues=eigenvalues.tolist(),
    )


def main():
    rng = np.random.default_rng(20261120)
    rows = []
    t0, c0 = perf_counter(), process_time()
    for name in ("r10_twod", "r11_crossgrid8", "r12_confirm24"):
        model, cases = load_cohort(name)
        bank = prepare(model)
        with np.load(OUT / ("current_cycle_equal_" + name + ".npz")) as saved:
            estimates = {
                i: x[-1]
                for (_, i, a), x in zip(saved["key"], saved["images"])
                if a == "voxel_only"
            }
        for case in cases:
            # Each anchor is determined only by the archived noisy readings.
            for stage, m in (
                ("background", np.full(len(model["G"]), 0.1)),
                ("voxel_fit", estimates[case["id"]]),
            ):
                row = spectrum_and_retention(model, m, case["Y"], bank, rng)
                rows.append(dict(cohort=name, case=case["id"], stage=stage, **row))
    (OUT / "current_geometry.json").write_text(
        json.dumps(
            dict(wall_s=perf_counter() - t0, cpu_s=process_time() - c0, rows=rows),
            indent=1,
        )
        + "\n"
    )


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
