"""Replay whole accepted chords from frozen arrays and summarize physical-area bounds."""
from fractions import Fraction
import csv
import json
from pathlib import Path
from time import perf_counter, process_time

from continuation import (HERE, FROZEN, BOX, Model, Reduction, certify_segment, fro, up, rd,
                          np, threadpool_limits, write_json)

GOAL = 2.0
SCOPES = {
    "": "anchor-response predictor",
    "_ball": "original noisy observation",
    "_broad": "original noisy observation / expanded exterior",
    "_whole_region": "original noisy observation / whole region",
}


def clock_delta(start):
    return dict(wall_s=perf_counter()-start[0], cpu_s=process_time()-start[1])


def safe_mean(x):
    exact = sum(map(Fraction.from_float, x.tolist()), Fraction()) / len(x)
    rounded = float(exact)
    return float(np.nextafter(rounded, -np.inf)) if Fraction.from_float(rounded) > exact else rounded


def upper_cost(change, area):
    magnitude = abs(change)
    n = len(change)
    return dict(L1_m2=up(area * float(magnitude.sum()) * (1 + rd(n))),
                L2_m=up(np.sqrt(area) * fro(magnitude)),
                Linf=up(magnitude.max()))


def cache_name(case, suffix):
    if suffix == "_broad":
        count = 852 if case.startswith("n32") else 512
        return f"{case}_reduction_ext{count}.npz"
    if suffix == "_whole_region":
        count = 512 if case.startswith("n32") else 256
        return f"{case}_reduction_ext{count}_whole_roi.npz"
    return f"{case}_reduction.npz"


def archived_reduction(case, suffix, model):
    with np.load(HERE / cache_name(case, suffix), allow_pickle=False) as cache:
        rec = json.loads(str(cache["receipt"].item()))
        return Reduction(model, cache["x0"], cache["indices"], cache["initial"], cache["transfer"],
                         cache["kernel"], cache["field"], rec["kernel_error_upper"],
                         rec["incident_field_error_upper"], rec["observation_transfer_error_upper"],
                         rec["initial_prediction_error_upper"], rec)


def inspect(paths_file):
    start = perf_counter(), process_time()
    j = json.loads(paths_file.read_text())
    case = j["case"]
    suffix = paths_file.stem.removeprefix(case + "_paths")
    m = Model.read(case)
    r = archived_reduction(case, suffix, m)
    output = []
    support = np.zeros(m.n, bool)
    support[r.indices] = True
    arrays_path = paths_file.with_suffix(".npz")
    with np.load(arrays_path, allow_pickle=False) as arrays:
        roi, patch = arrays["roi_mask"], arrays["patch"]
        task = arrays["task_cells"] if "task_cells" in arrays else patch
        area = j["cell_area_m2"]
        if not np.array_equal(arrays["Y_observed"], m.Y):
            raise ValueError("Stored path does not refer to the original noisy observation")
        for method, p in j["paths"].items():
            materials = arrays[method + "_materials"]
            if np.any(materials < 0.) or np.any(materials > BOX):
                raise ValueError("Stored materials leave the convex material box")
            if not np.array_equal(materials[0], r.x0) or np.any(materials[:, ~support] != r.x0[~support]):
                raise ValueError("Path moves outside the audited support")
            if p["strict_segments"] != len(materials)-1:
                raise ValueError("Segment number differs from stored knots")
            knots = [r.evaluation(node[r.indices], audit=True) for node in materials]
            verified = [certify_segment(r, a, b) for a, b in zip(knots[:-1], knots[1:])]
            if not verified or not all(v["strict"] for v in verified):
                raise ArithmeticError(f"{case}/{method}: archived polyline did not verify")
            peak_field_error = max(v["unsigned_output_error_over_r"] for v in verified)
            means = [safe_mean(node[roi]) for node in materials]
            task_means = [safe_mean(node[task]) for node in materials]
            if np.any(np.diff(task_means) < 0):
                raise ValueError("Archived region means are not monotone")
            costs = [upper_cost(node[~roi] - r.x0[~roi], area) for node in materials]
            idx = int(np.argmax([v["residual_upper_over_r"] for v in verified]))
            mid = (materials[idx] + materials[idx+1]) / 2
            original_check = m.forward(mid, fresh=True)
            reduced_mid = r.evaluation(mid[r.indices])
            reduction_diff = float(np.linalg.norm(reduced_mid-original_check["pred"])/m.radius)
            mid_ratio = original_check["residual_ratio"]
            if mid_ratio >= verified[idx]["residual_upper_over_r"]:
                raise ArithmeticError(f"{case}/{method}: original midpoint exceeds certified enclosure")
            scope = SCOPES[suffix]
            row = dict(case=case, observation_id=j["observation_id"], archive=paths_file.name,
                       method=method, scope=scope, grid=j["grid"], area_target_m2=area * len(task),
                       area_original_roi_m2=area * int(roi.sum()), area_patch_m2=area * len(patch),
                       task_cells=len(task), exterior_support_cells=int(j["ext_support_cells"]),
                       segments=len(verified), attempted_segments=len(p["attempts"]),
                       mean_task_lower=task_means[-1], mean_region_lower=means[-1],
                       mean_patch_lower=safe_mean(materials[-1][patch]),
                       patch_min=float(materials[-1][patch].min()), patch_peak=float(materials[-1][patch].max()),
                       residual_upper_over_r=max(v["residual_upper_over_r"] for v in verified),
                       outward_rounding_solve_and_adjoint_upper_over_r=peak_field_error,
                       maximum_reduced_inverse_defect=max(v["inverse_defect_upper"] for v in verified),
                       original_midpoint_native_ratio=mid_ratio,
                       reduced_original_midpoint_difference_over_r=reduction_diff,
                       terminal_exterior_cost_upper=costs[-1],
                       complete_path_exterior_budget_upper={k:max(v[k] for v in costs) for k in costs[-1]},
                       remaining_to_declared_targets=dict(task_mean=max(0., GOAL-task_means[-1]),
                            all_patch_cells_above_goal=max(0., GOAL-float(materials[-1][patch].min())),
                            material_box_roi_mean_gap=max(0., BOX-means[-1])),
                       stop=p["stop"], method_clocks={k:p[k] for k in ("wall_s", "cpu_s")},
                       reduced_anchor_setup_clocks={k:r.receipt[k] for k in ("wall_s", "cpu_s")},
                       recheck_clocks=clock_delta(start))
            output.append(row)
    return output


def compare_rank_one(case):
    folder, roi = ("stable_profiles", "roi_01") if case.startswith("n32") else ("gn_profiles", "roi_02")
    p = FROZEN / "results" / folder / case / f"{roi}_inner.json"
    item = json.loads(p.read_text())
    peak = item["attainable_peak_to"]
    n = 32 if case.startswith("n32") else 96
    return dict(source=str(p.relative_to(FROZEN)), physical_area_m2=(2 / n)**2,
                peak_to=peak, region_mean_to=peak/item["mask_cells"],
                report_scope="historical one-pixel change from the same vacuum seed")


def inspect_b():
    p = HERE/"B_independent_trajectory.json"
    if not p.exists():
        return None
    from b_control import dark_apertures, layouts
    start = perf_counter(), process_time()
    j = json.loads(p.read_text())
    with np.load(FROZEN/"ports"/"public"/"operators"/"n32_f400_a08.npz",allow_pickle=False) as a:
        op = dict(a)
    _, exc, rec = dark_apertures(op)
    _, _, _, coeff = layouts(op,exc,rec)["B_interleaved"]
    with np.load(p.with_suffix(".npz"),allow_pickle=False) as a:
        model = Model(coeff,a["B_observation"],j["radius_B"],"n32_f400_a08 / B",32)
        x0,ids = a["vacuum"],a["support"]
        info = j["reduction"]
        reduction = Reduction(model,x0,ids,a["F0"],a["D"],a["K"],a["e"],info["kernel_error_upper"],
                    info["incident_field_error_upper"],info["observation_transfer_error_upper"],
                    info["initial_prediction_error_upper"],info)
        knots = [reduction.evaluation(x[ids],audit=True) for x in a["verified_polyline"]]
        certificates = [certify_segment(reduction,u,v) for u,v in zip(knots[:-1],knots[1:])]
        if len(certificates)!=j["continuous"]["strict_segments"] or not all(k["strict"] for k in certificates):
            raise ArithmeticError("Archived B path failed its cold segment replay")
        b_dir = (a["high_patch"]-x0)[~a["roi"]]
        return dict(segments=len(certificates), scope="original archived B acquisition only",
                    strict_segment_max_residual_over_r=max(k["residual_upper_over_r"] for k in certificates),
                    independently_refitted_pair_budget_upper=upper_cost(b_dir,(2/32)**2),
                    exact_chart=j["exact_chart"], **clock_delta(start))


def main():
    start = perf_counter(), process_time()
    files = sorted(HERE.glob("*_paths*.json"))
    rows = []
    for file in files:
        print("replay", file.name, flush=True)
        rows.extend(inspect(file))
    cases = sorted({row["case"] for row in rows})
    exact = {}
    for case in cases:
        source = HERE / f"{case}_exact_chart.json"
        if source.exists():
            exact[case] = json.loads(source.read_text())["exact_constant_field_chart"]
    exact_shape = {case:json.loads((HERE/f"{case}_exact_shape.json").read_text())
                   for case in cases if (HERE/f"{case}_exact_shape.json").exists()}
    write_json(HERE / "summary.json", dict(physical_delta_metrics="equal-area weighted L1/L2/Linf",
        declared_fixed_patch_goal=GOAL, layout=dict(frequency_MHz=400, complex_reads=64,
            pulses=8, available_W_per_pulse=.125, source_column_power_sum_W=1.,
            pulse_dwell_s=None, deterministic_original_Frobenius_ball=True,
            fixed_port_reference_ohm=50., port_count=24), paths=rows, historical_single_pixel={c:compare_rank_one(c) for c in cases},
        exact_constant_field=exact, constant_field_shape_free=exact_shape, original_B=inspect_b(), replay_clock=clock_delta(start)))
    fields = ("case", "archive", "method", "task_cells", "area_target_m2", "exterior_support_cells",
              "segments", "mean_task_lower", "mean_region_lower", "mean_patch_lower", "patch_min",
              "patch_peak", "residual_upper_over_r", "original_midpoint_native_ratio",
              "reduced_original_midpoint_difference_over_r")
    with (HERE/"summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print("rechecked", len(rows), "paths", clock_delta(start), flush=True)


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
