"""Certify the full-state conditioning margin on the existing physical noise cover."""

import json
from time import perf_counter, process_time

from flint import acb_mat, arb, ctx
import numpy as np

from certify_fold import (
    PRECISION, complex_norm, diagonal, float_interval, full_matrix, identity,
    positive_ldl_pivots, pulse_operators,
)
from host_model import RESULTS
from noisy_cover import path_box


COVER_ARM = "physical_direct_centered"
CENTER_FLOOR = "0.2599"
REQUIRED_FLOOR = "0.25"


def certify_cell(path, order, material_box, G, G_active, shifted_identity):
    loc, den = path_box(path, order)
    center, halfwidth = [], []
    for (left, right), index, subdivisions in zip(material_box, loc, den, strict=True):
        width = (right - left) / int(subdivisions)
        center.append(left + width * (2 * int(index) + 1) / 2)
        halfwidth.append(width / 2)
    A = full_matrix(center, G)
    pivots = positive_ldl_pivots(A.conjugate().transpose() * A - shifted_identity)
    K = A.transpose().solve(G_active.transpose()).transpose()
    relative = complex_norm(diagonal(halfwidth) * K).upper()
    lower = arb(CENTER_FLOOR) * (1 - relative)
    return {
        "path": path,
        "sigma_min_lower": float_interval(lower)["lo"],
        "relative_perturbation_frobenius_upper": float_interval(relative)["hi"],
        "minimum_center_shifted_gram_ldl_pivot_lower": min(float_interval(q)["lo"] for q in pivots),
    }, lower > arb(REQUIRED_FLOOR)


def main():
    ctx.prec = PRECISION
    start_wall, start_cpu = perf_counter(), process_time()
    report = json.loads((RESULTS / "noisy_cover.json").read_text())
    order = report["axis_order_zero_based"]
    material_box = [
        [arb(float(endpoint)) for endpoint in interval]
        for interval in report["scopes"]["physical"]["material_box_binary64"]
    ]
    with np.load(RESULTS / "noisy_cover_cells.npz") as archive:
        inside, undecided = (archive[f"{COVER_ARM}_{kind}"] for kind in ("inside", "undecided"))
        paths = [int(path) for path in np.concatenate((inside, undecided))]
    if len(paths) != len(set(paths)) or not paths:
        raise ValueError("The retained archive must contain distinct nonempty paths.")
    normal_section = json.loads((RESULTS / "normal_section.json").read_text())["global_line"]
    if normal_section["outer_method"] != COVER_ARM or len(paths) != normal_section["outer_cells"]:
        raise ValueError("Cover archive disagrees with the continuous normal-section certificate.")

    G, _, _ = pulse_operators()
    G_active = acb_mat([
        [G[i, j] for j in range(G.ncols())] for i in range(len(material_box))
    ])
    floor = arb(CENTER_FLOOR)
    shifted_identity = identity(G.ncols()) * (floor * floor)
    setup_wall, setup_cpu = perf_counter() - start_wall, process_time() - start_cpu
    certified, unresolved = [], []
    for path in paths:
        try:
            cell, passes = certify_cell(path, order, material_box, G, G_active, shifted_identity)
        except (ArithmeticError, ValueError) as exc:
            unresolved.append({"path": path, "reason": str(exc)})
            continue
        (certified if passes else unresolved).append(cell)

    margin = {
        "scope": "Existing physical_direct_centered outer cover of the fixed 36-state scalar-TM model; each bound applies to the entire material cell.",
        "input_cover": "noisy_cover_cells.npz",
        "cover_arm": COVER_ARM,
        "depth_budget": report["depth_budget"],
        "unique_parent_cells": len(paths),
        "inside_parent_cells": len(inside),
        "undecided_parent_cells": len(undecided),
        "subdivisions_added": 0,
        "precision_bits": ctx.prec,
        "center_sigma_min_strictly_exceeds": CENTER_FLOOR,
        "target_sigma_min_strictly_exceeds": REQUIRED_FLOOR,
        "method": (
            "Interval LDL* proves A_c* A_c − center_floor² I positive definite. "
            "Write A=(I−U D_delta K) A_c with K=U^T G A_c^{-1}. "
            "Interval solve encloses K and the cell halfwidths h give "
            "||U D_delta K||_2 ≤ ||D_h K||_F uniformly. "
            "The whole-cell lower bound is center_floor * (1 − ||D_h K||_F)."
        ),
        "certified_parent_cells": len(certified),
        "unresolved_parent_cells": len(unresolved),
        "whole_outer_cover_certified": not unresolved,
        "minimum_sigma_min_lower_over_certified_cells": min(
            (r["sigma_min_lower"] for r in certified), default=None
        ),
        "maximum_relative_frobenius_bound_over_certified_cells": max(
            (r["relative_perturbation_frobenius_upper"] for r in certified), default=None
        ),
        "tightest_certified_cell": min(certified, key=lambda r: r["sigma_min_lower"], default=None),
        "setup_wall_seconds": setup_wall,
        "setup_cpu_seconds": setup_cpu,
        "wall_seconds": perf_counter() - start_wall,
        "cpu_seconds": process_time() - start_cpu,
        "cells": certified,
        "unresolved": unresolved,
    }
    (RESULTS / "host_noise_margin.json").write_text(json.dumps(margin, indent=2) + "\n")
    print(json.dumps({k: v for k, v in margin.items() if k not in ("cells", "unresolved")}, indent=2))


if __name__ == "__main__":
    main()
