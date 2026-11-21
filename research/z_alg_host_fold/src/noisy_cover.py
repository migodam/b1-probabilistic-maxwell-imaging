"""Outward Bernstein covers of a shared bounded-noise material feasible set."""

import argparse
from fractions import Fraction
from itertools import combinations
import json
import math
from pathlib import Path
import platform
from statistics import median
import subprocess
import tempfile
import time

from flint import acb, acb_mat, arb, ctx
import numpy as np


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
SOURCE = Path(__file__).with_suffix(".cpp")
PHYSICAL_BOX = (0.5, 30.0)
LOCAL_RADIUS = 0.5
NOISE_RADIUS = math.ldexp(5, -22)
PRECISION_BITS = 128
MATERIAL_COUNT = 6
VOLUME_TARGETS = (
    .99, .9, .5, .3, .1, .03, .01, .003, .001, .0003, .0001,
    3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7,
    3e-8, 1e-8, 3e-9, 1e-9, 3e-10, 1e-10, 3e-11, 1e-11,
    3e-12, 1e-12, 3e-13, 1e-13,
)


def exact_complex(value):
    value = complex(value)
    return acb(arb(value.real), arb(value.imag))


def interval(value):
    return arb(value["lo"]).union(arb(value["hi"]))


def complex_interval(value):
    return acb(interval(value["real"]), interval(value["imag"]))


def squared_modulus(value):
    # abs_lower/abs_upper preserve zero as the lower bound across a sign change.
    r, t = value.real, value.imag
    return (
        (r.abs_lower() ** 2 + t.abs_lower() ** 2)
        .union(r.abs_upper() ** 2 + t.abs_upper() ** 2)
    )


def dot(row, column):
    return sum((a * b for a, b in zip(row, column, strict=True)), acb(0))


def load_inputs():
    with np.load(RESULTS / "host_operators.npz") as source:
        saved = {key: source[key] for key in ("H", "B", "R", "y0", "zc", "y_critical", "y_minus", "y_plus")}
    bounded = json.loads((RESULTS / "host_operators_interval.json").read_text())
    unpack = np.vectorize(complex_interval, otypes=[object])
    operators = tuple(unpack(np.asarray(bounded[key], dtype=object)) for key in ("H", "B", "R", "y0"))
    H, B, R, host = operators
    assert H.shape == (6, 6) and B.shape == (6, 3) and R.size == 6 and host.size == 3
    return saved, (H, B, R.reshape(6), host.reshape(3))


def dyadic_triplets(box):
    return [
        (arb(float(lo)), (arb(float(lo)) + arb(float(hi))) / 2, arb(float(hi)))
        for lo, hi in zip(*box, strict=True)
    ]


def quadratic_controls(grid):
    """Convert endpoint/midpoint samples of a degree-two tensor to Bernstein form."""
    for axis in range(grid.ndim):
        controls = np.moveaxis(grid, axis, 0)
        controls[1] = 2 * controls[1] - (controls[0] + controls[2]) / 2
    return grid


def fixed_point(grid, include_upper=False, bits=53):
    """Round both signs downward on one exact binary scale."""
    peak = max(float(abs(value).upper()) for value in grid.flat)
    shift = 0 if peak == 0 else bits - math.frexp(peak)[1]
    scale = arb(2) ** shift

    def coefficients(sign):
        rounded = [(sign * entry * scale).lower().floor().unique_fmpz() for entry in grid.flat]
        if any(value is None or abs(int(value)) >= 2 ** (bits + 2) for value in rounded):
            raise ArithmeticError("Bernstein fixed-point controls exceed certified integer bound")
        return np.asarray([int(value) for value in rounded], dtype="<i8")

    return (coefficients(1), coefficients(-1)) if include_upper else coefficients(1)


def direct_coefficients(box, operators, measured, epsilon, centered):
    """Score sign for the exact six-material additive-noise test."""
    H, B, R, host = operators
    levels = dyadic_triplets(box)
    grid = np.empty((3,) * MATERIAL_COUNT, dtype=object)
    linear = np.empty((8,) + (2,) * MATERIAL_COUNT, dtype=object)
    for idx in np.ndindex(grid.shape):
        z = [levels[i][j] for i, j in enumerate(idx)]
        forward = [[acb(int(i == j)) - H[i, j] * z[j] for j in range(6)] for i in range(6)]
        denominator = acb_mat(forward).det()
        tail = [-R[j] * z[j] for j in range(6)]
        residual = []
        for reading in range(3):
            # This bordered determinant equals det(forward) * (F-D)_reading.
            numerator = acb_mat([
                forward[i] + [B[i, reading]] for i in range(6)
            ] + [tail + [host[reading] - measured[reading]]]).det()
            residual.append(numerator)
        grid[idx] = sum((squared_modulus(x) for x in residual), arb(0)) - epsilon ** 2 * squared_modulus(denominator)
        if centered and all(i != 1 for i in idx):
            corner = tuple(i // 2 for i in idx)
            for j, value in enumerate(residual + [epsilon * denominator]):
                linear[(2 * j,) + corner] = value.real
                linear[(2 * j + 1,) + corner] = value.imag
    score = fixed_point(quadratic_controls(grid), include_upper=True)
    if not centered:
        return score
    # The six-term dot products and squared tests fit within signed 128 bits.
    return (*score, *fixed_point(linear, include_upper=True, bits=38))


def chart_coefficients(box, operators, measured, epsilon, approximate_B):
    """All fifteen R10 minors for three unobserved complex receiver directions."""
    H, B, R, host = operators
    source_rows = max(combinations(range(6), 3), key=lambda ids: abs(np.linalg.det(approximate_B[list(ids)])))
    free_rows = [i for i in range(6) if i not in source_rows]
    inverse = acb_mat([[B[i, a] for a in range(3)] for i in source_rows]).inv()
    right_map = np.asarray([
        [inverse[a, source_rows.index(i)] if i in source_rows else acb(0) for i in range(6)]
        for a in range(3)
    ], dtype=object)
    null_space = [
        [acb(int(j == i)) - dot(B[i], right_map[:, j]) for j in range(6)] for i in free_rows
    ]
    centered = [d - h for d, h in zip(measured, host, strict=True)]
    cell_data = []
    for j, levels in enumerate(dyadic_triplets(box)):
        line = []
        for z in levels:
            u = [acb(int(i == j)) - H[i, j] * z for i in range(6)]
            b = [dot(row, u) for row in null_space]
            source = [dot(row, u) for row in right_map]
            line.append((b, source, dot(centered, source) - z * R[j]))
        cell_data.append(line)

    tensors = []
    for ids in combinations(range(6), 4):
        grid = np.empty((3,) * 4, dtype=object)
        for idx in np.ndindex(grid.shape):
            values = [cell_data[i][j] for i, j in zip(ids, idx, strict=True)]
            minor_col, src, residual = zip(*values, strict=True)
            weights = [
                (-1) ** cut * acb_mat([
                    [minor_col[j][i] for j in range(4) if j != cut] for i in range(3)
                ]).det()
                for cut in range(4)
            ]
            source_coeff = [dot(weights, col) for col in zip(*src, strict=True)]
            grid[idx] = squared_modulus(dot(weights, residual)) - epsilon ** 2 * sum(
                (squared_modulus(x) for x in source_coeff), arb(0)
            )
        stride_axis = [3 - ids.index(i) if i in ids else -1 for i in range(6)]
        tensors.append((np.asarray(stride_axis, dtype=np.int8), fixed_point(quadratic_controls(grid))))
    return tensors, list(source_rows)


def score_material(z, operators, observation):
    H, B, R, host = operators
    forward = acb_mat([[acb(int(i == j)) - H[i, j] * z[j] for j in range(6)] for i in range(6)])
    field = forward.solve(acb_mat(B.tolist()))
    predicted = [host[a] + dot(R, [z[j] * field[j, a] for j in range(6)]) for a in range(3)]
    return sum((squared_modulus(p - y) for p, y in zip(predicted, observation, strict=True)), arb(0))



def score_point(material, operators, observation):
    return score_material([arb(float(q)) for q in material], operators, observation)


def outward_number(ball):
    low = float(ball.lower())
    high = float(ball.upper())
    return [math.nextafter(low, -math.inf), math.nextafter(high, math.inf)]


def seed_real_jacobian(saved):
    z = saved["zc"]
    H, B, R = saved["H"], saved["B"], saved["R"].reshape(6)
    A = np.eye(6) - H * z[None, :]
    fields = np.linalg.solve(A, B)
    dual = R + H.T @ np.linalg.solve(A.T, R * z)
    jacobian = (dual[:, None] * fields).T
    return np.stack((jacobian.real, jacobian.imag), axis=1).reshape(6, 6)


def axis_order(saved):
    return np.argsort(-np.linalg.norm(seed_real_jacobian(saved), axis=0), kind="stable").tolist()


def path_box(path, order):
    widths = [1] * 6
    coords = [0] * 6
    depth = path.bit_length() - 1
    for level in range(depth):
        axis = order[level % 6]
        side = (path >> (depth - level - 1)) & 1
        widths[axis] *= 2
        coords[axis] = 2 * coords[axis] + side
    return coords, widths


def summarize_paths(inside, undecided, order, box):
    retained = np.concatenate([inside, undecided])
    if retained.size == 0:
        return None
    # All endpoints are exact dyadics when decoded at the largest axis split count.
    depth = int(np.max(retained)).bit_length() - 1
    counts = np.array([depth // 6] * 6)
    counts[order[:depth % 6]] += 1
    common = 1 << counts
    coords = np.zeros((len(retained), 6), dtype=np.int64)
    widths = np.ones((len(retained), 6), dtype=np.int64)
    depths = np.frexp(retained.astype(float))[1] - 1  # cell codes use at most 49 bits
    for level in range(depth):
        axis = order[level % 6]
        mask = depths > level
        coords[mask, axis] = 2 * coords[mask, axis] + (
            (retained[mask] >> (depths[mask] - level - 1).astype(np.uint64)) & 1
        ).astype(np.int64)
        widths[mask, axis] *= 2
    scale = common / widths
    lower = np.min(coords * scale, axis=0).astype(int)
    upper = np.max((coords + 1) * scale, axis=0).astype(int)
    hull = []
    for j in range(6):
        lo, hi = [arb(float(box[i, j])) for i in range(2)]
        hull.append([
            outward_number(lo + (hi - lo) * int(lower[j]) / int(common[j]))[0],
            outward_number(lo + (hi - lo) * int(upper[j]) / int(common[j]))[1],
        ])
    return hull


def run_arm(scope, box, arm, operators, measured, saved, order, exe, work, visits, depth, shared_setup):
    path = work / f"{scope}-{arm}.bin"
    frontier_path = work / f"{scope}-{arm}-frontier.bin"
    cpu_start, wall_start = time.process_time(), time.perf_counter()
    source_rows = None
    with path.open("wb") as out:
        out.write(bytes(order))
        if arm.startswith("direct"):
            for tensor in direct_coefficients(box, operators, measured, arb(NOISE_RADIUS), arm == "direct_centered"):
                out.write(tensor.tobytes())
        else:
            tensors, source_rows = chart_coefficients(box, operators, measured, arb(NOISE_RADIUS), saved["B"])
            for axes, tensor in tensors:
                out.write(axes.tobytes())
                out.write(tensor.tobytes())
    build_cpu = time.process_time() - cpu_start
    build_wall = time.perf_counter() - wall_start
    completed = subprocess.run(
        [str(exe), arm, str(path), str(frontier_path), str(visits), str(depth)],
        capture_output=True, check=True, text=True,
    )
    whole_wall = time.perf_counter() - wall_start
    frontier = [json.loads(line) for line in completed.stdout.splitlines()]
    for state in frontier:
        state["cpu_seconds_including_setup"] = build_cpu + shared_setup["cpu_seconds"] + state["cpu_seconds_sieve"]
        state["wall_seconds_including_setup"] = build_wall + shared_setup["wall_seconds"] + state["wall_seconds_sieve"]
    count_inside, count_undecided = np.fromfile(frontier_path, dtype="<u8", count=2).tolist()
    paths = np.fromfile(frontier_path, dtype="<u8", offset=16)
    if count_inside + count_undecided != len(paths):
        raise ArithmeticError("Binary cell archive and counts disagree")
    inside = paths[:count_inside]
    undecided = paths[count_inside:]
    final = frontier[-1]
    if len(inside) != final["inside"] or len(undecided) != final["undecided"]:
        raise ArithmeticError("Frontier changed after final checkpoint")
    fraction = math.fsum(math.ldexp(1., -(int(cell).bit_length() - 1)) for cell in paths)
    if fraction != final["outer_fraction"]:
        raise ArithmeticError("Archived disjoint cells do not recover the outer volume")
    report = {
        "cpu_seconds_polynomial_setup": build_cpu,
        "wall_seconds_polynomial_setup": build_wall,
        "wall_seconds_complete_including_setup_subprocess_and_cell_write": whole_wall + shared_setup["wall_seconds"],
        "frontier": frontier,
        "final_hull_outward": summarize_paths(inside, undecided, order, box),
    }
    if source_rows is not None:
        report["minor_source_rows_zero_based"] = source_rows
    return report, inside, undecided


def median_run(runs, shared_offline, depth):
    first, inside, undecided = runs[0]
    records = [result for result, _, _ in runs]
    for _, a, b in runs[1:]:
        if not np.array_equal(a, inside) or not np.array_equal(b, undecided):
            raise ArithmeticError("Repeated timing runs changed final dyadic cells")
    for key in ("cpu_seconds_polynomial_setup", "wall_seconds_polynomial_setup",
                "wall_seconds_complete_including_setup_subprocess_and_cell_write"):
        values = [record[key] for record in records]
        first[key] = median(values)
        first[key + "_per_trial"] = values
    for index, state in enumerate(first["frontier"]):
        timing_trials = [record["frontier"][index] for record in records]
        algebraic_keys = [key for key in state if "_seconds_" not in key]
        if any(any(trial[key] != state[key] for key in algebraic_keys) for trial in timing_trials):
            raise ArithmeticError("Repeated timing runs changed dyadic refinement")
        for key in ("cpu_seconds_sieve", "wall_seconds_sieve",
                    "cpu_seconds_including_setup", "wall_seconds_including_setup"):
            state[key] = median(trial[key] for trial in timing_trials)
        for clock in ("cpu", "wall"):
            state[clock + "_seconds_with_upstream_certification"] = (
                state[clock + "_seconds_including_setup"] + shared_offline[clock + "_seconds"]
            )
        state["retained_equivalent_cells_at_maximum_depth"] = int(math.ldexp(state["outer_fraction"], depth))
    return first, inside, undecided


def matched_results(arms, baseline):
    a, b = (arms[key]["frontier"] for key in (baseline, "chart"))
    budgets = sorted(set(t["visited"] for t in a) & set(t["visited"] for t in b))

    def small(row):
        return {
            key: row[key] for key in (
                "visited", "retained", "undecided", "inside", "outer_fraction", "inner_fraction",
                "cpu_seconds_including_setup", "wall_seconds_including_setup",
                "cpu_seconds_with_upstream_certification", "wall_seconds_with_upstream_certification",
                "retained_equivalent_cells_at_maximum_depth",
            )
        }

    matched = []
    for volume in VOLUME_TARGETS:
        first = next((t for t in a if t["outer_fraction"] <= volume), None)
        second = next((t for t in b if t["outer_fraction"] <= volume), None)
        ratio = None
        if first and second:
            ratio = first["wall_seconds_including_setup"] / second["wall_seconds_including_setup"]
        matched.append({
            "maximum_outer_fraction": volume,
            "direct": small(first) if first else None,
            "chart": small(second) if second else None,
            "direct_over_chart_wall": ratio,
        })
    budget_records = []
    for budget in budgets:
        if budget not in (100, 1000, 10000, 50000, 100000, 200000, 400000) and budget != budgets[-1]:
            continue
        first = next(t for t in a if t["visited"] == budget)
        second = next(t for t in b if t["visited"] == budget)
        budget_records.append({"visited_cells": budget, "direct": small(first), "chart": small(second)})
    wall_records = []
    for seconds in (.05, .1, .15, .2, .3, .5, 1.):
        first = next((t for t in reversed(a) if t["wall_seconds_including_setup"] <= seconds), None)
        second = next((t for t in reversed(b) if t["wall_seconds_including_setup"] <= seconds), None)
        wall_records.append({
            "wall_seconds_budget": seconds,
            "direct": small(first) if first else None,
            "chart": small(second) if second else None,
        })
    return {"volume_thresholds": matched, "common_box_budgets": budget_records, "common_wall_budgets": wall_records}


def covers_point(material, overlap_paths, order, box):
    # Each exact dyadic point can lie on at most one cell per depth here except at a split.
    # Check closed cell boundaries in both children in that exceptional case.
    coordinates = [
        (Fraction(float(q)) - Fraction(float(lo))) / (Fraction(float(hi)) - Fraction(float(lo)))
        for q, lo, hi in zip(material, *box, strict=True)
    ]
    if any(t < 0 or t > 1 for t in coordinates):
        return None
    codes = {int(p) for p in overlap_paths}
    max_depth = max((p.bit_length() for p in codes), default=1) - 1
    current = [(1, coordinates)]
    for depth in range(max_depth + 1):
        if any(code in codes for code, _ in current):
            return True
        axis = order[depth % 6]
        following = []
        for code, point in current:
            side_value = 2 * point[axis]
            for side in (0, 1):
                if side <= side_value <= side + 1:
                    child = list(point)
                    child[axis] = side_value - side
                    following.append((code * 2 + side, child))
        current = following
    return False


def certified_witnesses(saved, operators, measured, archive, scopes, order):
    points = {"public_fold_seed": saved["zc"]}
    path = RESULTS / "point_fibers.json"
    if path.exists():
        fibers = json.loads(path.read_text())["fibers"]
        for side in ("minus", "plus"):
            for j, root in enumerate(fibers[side]["physical_roots"]):
                points[f"point_fiber_{side}_{j + 1}"] = root["material_midpoint"]
    output = {}
    for name, material in points.items():
        score = score_point(material, operators, measured) / arb(NOISE_RADIUS) ** 2
        inside = bool(score.upper() < 1)
        record = {"material": list(material), "noise_score_enclosure": outward_number(score), "rigorously_feasible": inside}
        record["in_outer_cover"] = {}
        for scope, box in scopes.items():
            for arm in ("direct", "direct_centered", "chart"):
                paths = np.concatenate((archive[f"{scope}_{arm}_inside"], archive[f"{scope}_{arm}_undecided"]))
                in_cover = covers_point(material, paths, order, box)
                if inside and in_cover is False:
                    raise ArithmeticError("Certified physical witness was excluded")
                record["in_outer_cover"][f"{scope}_{arm}"] = in_cover
        output[name] = record
    return output


def main(args):
    ctx.prec = PRECISION_BITS
    load_cpu, load_wall = time.process_time(), time.perf_counter()
    saved, operators = load_inputs()
    measured = [exact_complex(y) for y in saved["y_critical"].reshape(3)]
    epsilon = arb(NOISE_RADIUS)
    order = axis_order(saved)
    scopes = {
        "physical": np.repeat(np.asarray(PHYSICAL_BOX)[:, None], 6, axis=1),
        "local": np.clip(saved["zc"][None, :] + np.array([[-LOCAL_RADIUS], [LOCAL_RADIUS]]), *PHYSICAL_BOX),
    }
    shared_setup = {
        "cpu_seconds": time.process_time() - load_cpu,
        "wall_seconds": time.perf_counter() - load_wall,
    }
    fold = json.loads((RESULTS / "fold_certificate.json").read_text())
    critical = [complex_interval(y) for y in fold["critical_response_enclosure"]]
    critical_distance = sum((squared_modulus(a - b) for a, b in zip(critical, measured, strict=True)), arb(0))
    # A union around a zero lower endpoint can overhang zero by an Arb radius ulp.
    critical_upper = outward_number(critical_distance.upper().sqrt())[1]
    if not critical_upper < NOISE_RADIUS:
        raise ArithmeticError("The certified critical response does not lie within the preset noise set")
    query_distances = {}
    for side in ("minus", "plus"):
        obs = [exact_complex(y) for y in saved["y_" + side].reshape(3)]
        distance = sum((squared_modulus(a - b) for a, b in zip(obs, measured, strict=True)), arb(0)).sqrt()
        if not distance.upper() < epsilon:
            raise ArithmeticError("Side observation is outside the preset noise set")
        query_distances[side] = outward_number(distance)

    output = RESULTS / "noisy_cover.json"
    archive_path = RESULTS / "noisy_cover_cells.npz"
    if args.witnesses_only:
        report = json.loads(output.read_text())
        with np.load(archive_path) as archive:
            report["shared_point_witnesses"] = certified_witnesses(saved, operators, measured, archive, scopes, order)
        output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        return

    upstream_paths = ["host_model.json", "host_state_certificate.json", "fold_certificate.json"]
    upstream_records = [json.loads((RESULTS / path).read_text()) for path in upstream_paths]
    shared_offline = {
        **{key: sum(rec[key] for rec in upstream_records) for key in ("cpu_seconds", "wall_seconds")},
        "recorded_sources": upstream_paths,
    }
    report = {
        "model": "3 complex total-field receiver readings and six real susceptibilities with fixed known host",
        "noise_set": {
            "center_key": "host_operators.npz:y_critical (recorded binary64 center)",
            "center_real_imag": [[y.real, y.imag] for y in saved["y_critical"].reshape(3)],
            "complex_euclidean_radius": NOISE_RADIUS,
            "prescription": "5 * 2^-22, prescribed before both covers",
            "side_observation_distance_enclosures": query_distances,
            "fold_certificate_critical_response_distance_upper": critical_upper,
        },
        "precision_bits": PRECISION_BITS,
        "axis_order_zero_based": order,
        "hierarchy": "Breadth-first shared cyclic bisection, with sensitivity column norms at public seed setting tie order.",
        "visited_budget": args.visits,
        "depth_budget": args.depth,
        "platform": platform.platform(),
        "shared_input_setup": shared_setup,
        "upstream_common_physical_fixture_setup": shared_offline,
        "timing_repetitions": args.repeats,
        "scopes": {},
    }
    archive = {}
    with tempfile.TemporaryDirectory(prefix="host-fold-noisy-") as directory:
        work = Path(directory)
        exe = work / "noisy_cover"
        start = time.perf_counter()
        subprocess.run([args.compiler, "-std=c++20", "-O3", str(SOURCE), "-o", str(exe)], check=True)
        report["shared_optional_compile_wall_seconds"] = time.perf_counter() - start
        for scope, box in scopes.items():
            record = {"material_box_binary64": box.T.tolist(), "arms": {}}
            arm_names = ("direct", "direct_centered", "chart")
            repetitions = {arm: [] for arm in arm_names}
            for trial in range(args.repeats):
                offset = trial % len(arm_names)
                for arm in arm_names[offset:] + arm_names[:offset]:
                    repetitions[arm].append(run_arm(
                        scope, box, arm, operators, measured, saved, order, exe, work, args.visits, args.depth, shared_setup,
                    ))
            for arm in arm_names:
                results, inside, undecided = median_run(repetitions[arm], shared_offline, args.depth)
                print(scope, arm, results["frontier"][-1], flush=True)
                record["arms"][arm] = results
                archive[f"{scope}_{arm}_inside"] = inside
                archive[f"{scope}_{arm}_undecided"] = undecided
            record["matched_frontiers"] = {
                baseline: matched_results(record["arms"], baseline) for baseline in ("direct", "direct_centered")
            }
            report["scopes"][scope] = record
    report["shared_point_witnesses"] = certified_witnesses(saved, operators, measured, archive, scopes, order)
    np.savez_compressed(archive_path, **archive)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--visits", type=int, default=200000)
    parser.add_argument("--depth", type=int, default=48)
    parser.add_argument("--compiler", default="clang++")
    parser.add_argument("--repeats", type=int, choices=range(1, 10), default=3)
    parser.add_argument("--witnesses-only", action="store_true")
    main(parser.parse_args())
