"""Check complete positive task caps by bisecting all six nuisance coordinates."""

import argparse
from collections import Counter, deque
from functools import lru_cache
import json
from pathlib import Path
from time import perf_counter, process_time

from task_caps import ROOT, arb, ctx, dressed, near_face, np, sets, shift, up
from threadpoolctl import threadpool_limits


CASE = "disks08_low"
ARMS = ("base", "angle", "rich")
ETAS = ("0", "0.075", "0.015")
DEFAULT_BUDGET = 4096


def binary(value):
    return arb(float(value))


def round_down(value):
    return float(np.nextafter(float(value.lower()), -np.inf))


def real_norm(values):
    lower, upper = arb(0), arb(0)
    for x in values:
        lo, hi = x.abs_lower(), x.abs_upper()
        lower += lo * lo
        upper += hi * hi
    return lower.sqrt().union(upper.sqrt())


def replay_witness(nominal, center, radii, mask, arm, eta, delta, old,
                   concentrated=None, sign=1):
    # Try one previously frozen nuisance setting. A passing point proves the
    # tested cap is nonempty; a failing point says nothing about all other settings.
    h = sign * radii
    if concentrated is None:
        h[mask] -= sign * eta * arb("0.99")
    else:
        h[concentrated] -= sign * eta * int(mask.sum()) * arb(149) / 150
    h[~mask] = np.array([binary(v) for v in old], dtype=object) - center[~mask]
    if not all(x >= -r and x <= r for x, r in zip(h, radii, strict=True)):
        raise ArithmeticError("The replay point does not belong to the specification box.")
    levels = shift(nominal, h)[0]
    data = [scale * levels[a, b] for _, channels, scale in sets(arm) for a, b in channels]
    distance = real_norm(data)
    gap = sum(radii[mask] - sign * h[mask]) / int(mask.sum())
    return {
        "material_displacement": [float(x) for x in h],
        "task_gap_from_specified_extreme": [round_down(gap), up(gap)],
        "side": sign, "concentrated_roi_index": None if concentrated is None else int(concentrated),
        "residual_interval": [round_down(distance), up(distance)],
        "residual_over_delta_interval": [round_down(distance / delta), up(distance / delta)],
        "proved_inside_cap_and_noise_ball": bool(distance < delta and (gap < eta if eta > 0 else gap == eta)),
    }


def certify_cover(recenter, mask, radii, arm, eta, delta, budget):
    lower = tuple(-float(r) for r in radii[~mask])
    upper = tuple(float(r) for r in radii[~mask])
    for bound, r in zip(upper, radii[~mask], strict=True):
        if binary(bound) != r:
            raise ArithmeticError("Nuisance bounds must be exactly binary-representable.")
    pending = deque([("", lower, upper)])
    leaf_list = []
    decisions = Counter()
    splits = []
    visited = 0
    cache_before = recenter.cache_info()

    while pending:
        path, lo, hi = pending.popleft()
        visited += 1
        mid = tuple((a + b) / 2 for a, b in zip(lo, hi, strict=True))
        side = tuple(b - a for a, b in zip(lo, hi, strict=True))
        local = radii.copy()
        local[~mask] = [binary(s) / 2 for s in side]
        # These dyadic boxes have a small depth here. Check that their endpoints
        # are precisely the centers and radii used by the interval inequality.
        for a, b, c, radius in zip(lo, hi, mid, local[~mask], strict=True):
            if binary(c) - radius != binary(a) or binary(c) + radius != binary(b):
                raise ArithmeticError("Floating subdivision failed to cover the exact box.")
        reason, residual_floor = near_face([recenter(mid)], mask, local, 1, arm, eta, delta)
        decisions[reason] += 1
        if reason == "excluded" or visited + len(pending) + 2 > budget:
            leaf_list.append([path, reason, residual_floor])
            continue
        j = max(range(len(side)), key=side.__getitem__)
        lower_right = lo[:j] + (mid[j],) + lo[j+1:]
        upper_left = hi[:j] + (mid[j],) + hi[j+1:]
        splits.append([path, j])
        pending.append((path + "0", lo, upper_left))
        pending.append((path + "1", lower_right, hi))

    tracked = {path for path, _, _ in leaf_list} | {path for path, _ in splits}
    if "" not in tracked or len(tracked) != visited or any(path + side not in tracked
            for path, _ in splits for side in ("0", "1")):
        raise ArithmeticError("The logged subdivision failed to account for every child.")
    volume_share_not_proved = sum(2.0 ** (-len(path)) for path, reason, _ in leaf_list
                                  if reason != "excluded")
    excluded = all(reason == "excluded" for _, reason, _ in leaf_list)
    cache_after = recenter.cache_info()
    return {
        "status": "all_leaves_excluded" if excluded else "unresolved_complete_cover",
        "visited_boxes": visited, "decisions": dict(decisions),
        "excluded_leaves": sum(reason == "excluded" for _, reason, _ in leaf_list),
        "unresolved_leaves": sum(reason != "excluded" for _, reason, _ in leaf_list),
        "max_depth": max(len(path) for path, _, _ in leaf_list),
        "minimum_lower_bound_all_leaves": min(lower for _, _, lower in leaf_list),
        "unresolved_volume_fraction": volume_share_not_proved,
        "exact_recenters": cache_after.misses - cache_before.misses,
        "reused_recenters": cache_after.hits - cache_before.hits,
        "split_paths_and_zero_based_local_axis": splits,
        "leaf_paths_status_residual_lower_bound": leaf_list,
    }


def run(output, budget):
    start_wall, start_cpu = perf_counter(), process_time()
    source = json.loads((ROOT / "array/results/summary.json").read_text())
    case = next(x for x in source["cases"] if x["id"] == CASE)
    center_float = np.array(case["material_center"])
    center = np.array([binary(x) for x in center_float], dtype=object)
    radii = np.array([min(arb("0.75"), x - arb("0.5"), arb(30) - x) for x in center], dtype=object)
    with np.load(ROOT / "array/results/physical_arrays.npz", allow_pickle=False) as raw:
        mask = raw[f"{CASE}__roi_mask"]
        nominal = dressed(raw, case["n"], center_float, "f0")
    prior = 2 * sum(radii[mask]) / int(mask.sum())
    base_h = np.array([r if inside else arb(0) for r, inside in zip(radii, mask, strict=True)],
                      dtype=object)

    @lru_cache(maxsize=8192)
    def recenter(mid):
        h = base_h.copy()
        h[~mask] = [binary(x) for x in mid]
        return shift(nominal, h)

    rows = []
    brackets = {}
    with np.load(ROOT / "faces/face_parameters.npz", allow_pickle=False) as trial:
        for eta_text in ETAS:
            eta = arb(eta_text)
            for arm in ARMS:
                start, cpu = perf_counter(), process_time()
                delta = binary(next(x["radius"] for x in source["pairs"]
                                    if x["case"] == CASE and x["arm"] == arm))
                old = trial[f"{CASE}__0.75__{arm}__GN_outside_+1"]
                witness = replay_witness(nominal, center, radii, mask, arm, eta, delta, old)
                if witness["proved_inside_cap_and_noise_ball"]:
                    report = {"status": "witnessed_cap_nonempty", "visited_boxes": 0}
                else:
                    report = certify_cover(recenter, mask, radii, arm, eta, delta, budget)
                vertices = []
                if eta > 0 and report["status"] == "unresolved_complete_cover":
                    vertices = [replay_witness(nominal, center, radii, mask, arm, eta, delta, old,
                                               concentrated=i) for i in np.flatnonzero(mask)]
                proved_points = [p for p in [witness, *vertices] if p["proved_inside_cap_and_noise_ball"]]
                conclusion = "witnessed_cap_nonempty" if proved_points else report["status"]
                report.update({
                    "arm": arm, "positive_task_cap_eta": eta_text, "noise_radius": float(delta),
                    "witness_replay": witness,
                    "concentrated_witness_replays": vertices, "positive_cap_conclusion": conclusion,
                    "fixed_reference_task_diameter_upper": up(prior - eta)
                        if report["status"] == "all_leaves_excluded" else up(prior),
                    "wall_seconds": perf_counter() - start, "cpu_seconds": process_time() - cpu,
                })
                interval = brackets.setdefault(arm, [0., up(prior)])
                interval[1] = min(interval[1], report["fixed_reference_task_diameter_upper"])
                if proved_points:
                    candidate = trial[f"{CASE}__0.75__{arm}__GN_outside_-1"]
                    negative = replay_witness(nominal, center, radii, mask, arm, arb(0), delta,
                                              candidate, sign=-1)
                    report["negative_face_replay"] = negative
                    if negative["proved_inside_cap_and_noise_ball"]:
                        cap_cost_upper = min(p["task_gap_from_specified_extreme"][1] for p in proved_points)
                        interval[0] = max(interval[0], round_down(prior - binary(cap_cost_upper)))
                rows.append(report)
                print(arm, eta_text, conclusion, report["visited_boxes"],
                      report.get("unresolved_leaves"), flush=True)
    result = {
        "scope": "fixed observed F(m0), archived binary64 disk-pulse coefficients treated as exact",
        "case": CASE, "precision_bits": ctx.prec, "flint_threads": ctx.threads,
        "halfwidths": [float(x) for x in radii],
        "nuisance_index_order": np.flatnonzero(~mask).tolist(),
        "ROI_index_order": np.flatnonzero(mask).tolist(), "orientation": "positive",
        "prior_task_diameter": [round_down(prior), up(prior)],
        "box_budget_per_test": budget,
        "split_rule": "breadth first, widest nuisance extent, ties choose first local coordinate",
        "certified_fixed_reference_task_diameter_intervals": brackets,
        "rows": rows, "wall_seconds": perf_counter() - start_wall,
        "cpu_seconds": process_time() - start_cpu,
    }
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--max-boxes", type=int, default=DEFAULT_BUDGET)
    cli.add_argument("--output", type=Path, default=ROOT / "bounds/refined_caps.json")
    args = cli.parse_args()
    if args.max_boxes < 1:
        cli.error("--max-boxes must be positive")
    ctx.threads = 1
    with threadpool_limits(limits=1, user_api="blas"):
        run(args.output, args.max_boxes)
