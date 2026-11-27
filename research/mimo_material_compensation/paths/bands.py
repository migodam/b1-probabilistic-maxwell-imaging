"""Export attainable phase bands with a certified full-path exterior budget."""
from fractions import Fraction
import json
from time import perf_counter, process_time

from continuation import HERE, FROZEN, BOX, np, threadpool_limits, write_json
from replay import safe_mean, upper_cost


def as_rational(value):
    exact = sum(map(Fraction.from_float, value.tolist()), Fraction()) / len(value)
    return [exact.numerator, exact.denominator]


def branch(points, task, roi, area, rows, layout, task_label, source):
    ext = ~roi
    result = []
    for i, receipt in enumerate(rows):
        d0 = points[i, ext]-points[0, ext]
        d1 = points[i+1, ext]-points[0, ext]
        one, two = upper_cost(d0, area), upper_cost(d1, area)
        threshold = {k:max(one[k],two[k]) for k in one}
        result.append(dict(layout=layout, task=task_label, source=source,
                           segment=i+1,
                           task_t_exact_start_fraction=as_rational(points[i, task]),
                           task_t_exact_end_fraction=as_rational(points[i+1, task]),
                           task_t_from=float(points[i, task].mean()),
                           task_t_to=safe_mean(points[i+1, task]),
                           sufficient_exterior_budget=threshold,
                           exterior_budget_upper_at_endpoints=[one,two],
                           data_ball_ratio_upper=receipt["residual_upper_over_r"]))
    return result


def main():
    tick, cpu = perf_counter(), process_time()
    items = []
    for path in sorted(HERE.glob("*_paths*.json")):
        if not path.stem.endswith(("_broad","_ball","_whole_region")):
            continue
        j = json.loads(path.read_text())
        with np.load(path.with_suffix(".npz"), allow_pickle=False) as data:
            for name, p in j["paths"].items():
                rows = [r["verification"] for r in p["attempts"] if r["strict"]]
                target = data["task_cells"] if "task_cells" in data else data["patch"]
                label = ("whole ROI shape-free mean" if path.stem.endswith("_whole_region") else
                         "patch shape-free mean" if name=="ball_shape_free" else "uniform connected patch")
                items += branch(data[name+"_materials"],target,data["roi_mask"],j["cell_area_m2"],rows,
                                "A / "+j["observation_id"],label,path.name+" / "+name)
    pairs = []
    for case in ("n32_c01","n96_c01"):
        archived = HERE/f"{case}_paths_ball.npz"
        old = FROZEN/"results"/"extended"/f"{case}.json"
        if not (archived.exists() and old.exists()):
            continue
        rec = json.loads(archived.with_suffix(".json").read_text())
        ext = json.loads(old.read_text())
        high_ratio = ext["upper_residual_over_r"]
        low_ratio = rec["strict_reduction"]["anchor_residual_upper_over_r"]
        with np.load(archived,allow_pickle=False) as ar:
            vacuum, high, roi, patch = [ar[k] for k in
                ("material_anchor","material_old_extended","roi_mask","patch")]
            if np.any(vacuum[roi]) or high[patch].min()<2. or max(low_ratio,high_ratio)>=1.:
                raise ValueError("The archived original A two-endpoint pair violates its physical contract")
            pairs.append(dict(layout="A / "+rec["observation_id"],task="vacuum ROI versus patch floor 2",
                  region_mean_difference_lower=safe_mean(high[roi]),
                  sufficient_exterior_budget=upper_cost((high-vacuum)[~roi],rec["cell_area_m2"]),
                  both_data_ball_ratio_upper=max(low_ratio,high_ratio),
                  source=str(old.relative_to(FROZEN)),status="proven hideable at the given exterior budget"))
    from_b = HERE/"B_independent_trajectory.json"
    if from_b.exists():
        b = json.loads(from_b.read_text())
        with np.load(from_b.with_suffix(".npz"),allow_pickle=False) as a:
            good = [v["verification"] for v in b["continuous"]["attempts"] if v["strict"]]
            items += branch(a["verified_polyline"],a["patch"],a["roi"],(2/32)**2,good,
                            "B / n32_c01 original archived acquisition","uniform connected patch",from_b.name)
            pairs.append(dict(layout="B / n32_c01 original archived acquisition",
                              task="vacuum ROI versus all 43 patch cells >= 2",
                              region_mean_difference_lower=safe_mean(a["high_patch"][a["roi"]]),
                              sufficient_exterior_budget=upper_cost((a["high_patch"]-a["vacuum"])[~a["roi"]],(2/32)**2),
                              both_data_ball_ratio_upper=max(b["independent_endpoint_audits"][k]["B"]["residual_upper_over_r"]
                                for k in ("vacuum","high_patch")),source=from_b.name,
                              status="proven hideable at the given exterior budget"))
    joint_file = HERE/"joint_original_AB.json"
    if joint_file.exists():
        joint = json.loads(joint_file.read_text())
        decision = joint["selected"]
        if set(decision)=={"vacuum","high_patch"}:
            both = all(r["jointly_verified"] for r in decision.values())
            pairs.append(dict(layout="A+B / n32_c01 original archived acquisitions",
                task="vacuum ROI versus all 43 patch cells >= 2",
                both_families_strictly_feasible=both,
                status="proven hideable at recorded budget" if both else "unresolved",
                selected=decision, sufficient_exterior_budget=joint.get("exterior_cost_of_selected_pair") if both else None,
                source=joint_file.name))
    write_json(HERE/"phase_bands.json",dict(
        interpretation="For each closed task interval the given budget bounds the complete segment relative to the registered vacuum anchor. Larger budgets preserve existence. Gaps remain unresolved by this file.",
        calibration="fixed nominal native coefficients; area-weighted differences, not total material mass",
        material_box=[0.,BOX], proved_hideable_segments=items, independently_refitted_endpoint_pairs=pairs,
        wall_s=perf_counter()-tick,cpu_s=process_time()-cpu))


if __name__=="__main__":
    with threadpool_limits(1):
        main()
