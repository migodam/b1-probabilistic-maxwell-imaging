"""Fit from the complete posterior hull midpoint using the same profiled TRF."""
import argparse
from itertools import product
import json
from time import perf_counter, process_time

import numpy as np
from flint import arb, ctx
from scipy.optimize import least_squares
from threadpoolctl import threadpool_limits

from data14 import ACTIONS, HALF_SPLIT_FIRST, ROOT, FrequencyTM, joint_fit, load_cases, measured_reads, ri
from freq_model import SCALE
from material_controls import GRID
from measure_acquire import archived_generator
from measure_recover import FrequencyReuse, solve
from measure_score import compatibility, physical_score
from measure_validation import load_validation_cases, load_validation_reads, private_trials

BUDGET = 780
STAGES = (("first_half", HALF_SPLIT_FIRST), ("both_halves", ACTIONS["split_0.02"]))


def scalar_radius_upper(boxes, point):
    """Directed rounding for the farthest box corner, then an Arb square root."""
    offsets = np.nextafter(np.abs(boxes - point), np.inf).max(axis=1)
    squared = np.nextafter(np.square(offsets), np.inf)
    summation = np.zeros(len(boxes))
    for d in squared.T:
        summation = np.nextafter(summation + d, np.inf)
    with ctx.workprec(128):
        val = (arb(float(summation.max())) / point.size).sqrt()
        return float(np.nextafter(float(val.upper()), np.inf))


def diagnostic(case, model, fit, action, new, boxes, material):
    residual, _, _, _, _ = fit.evaluate(material)
    observed = np.array([model.response(material, d)[0] for d in (0., *action.deltas)])
    old_ok, both_ok = compatibility(case, action, new, observed)
    in_cover = np.any(np.all((boxes[:, 0] <= material) & (material <= boxes[:, 1]), axis=1))
    return dict(material=material.tolist(), old_compatible=old_ok, common_data_compatible=both_ok,
                inside_box_union=bool(in_cover), observed_loss=float(np.linalg.norm(residual / case.sigma) ** 2),
                radius_upper=scalar_radius_upper(boxes, material))


def fit_in_hull(case, action, new, qboxes):
    cpu0, wall0 = process_time(), perf_counter()
    if not len(qboxes):
        return dict(before=None, after=None, state="empty_cover", fit_factors=0,
                    fit_cpu_s=process_time() - cpu0, fit_wall_s=perf_counter() - wall0)
    lo, hi = case.box
    # Conversion is exact: the source coordinates are dyadic integers in [0, GRID].
    boxes = lo + ((hi - lo) / GRID) * qboxes
    low, high = boxes[:, 0].min(axis=0), boxes[:, 1].max(axis=0)
    midpoint = .5 * (low + high)
    model = FrequencyReuse(case.geometry)
    fit = joint_fit(case, model, action, new)
    before = diagnostic(case, model, fit, action, new, boxes, midpoint)
    midpoint_factors = model.factorizations
    free = high > low

    def evaluate(z):
        x = midpoint.copy()
        x[free] = z
        return fit.evaluate(x)

    estimate = midpoint.copy()
    nfev = 0
    if np.any(free):
        n_unique = len({0., *action.deltas})
        result = least_squares(lambda z: ri(evaluate(z)[0]) / SCALE, midpoint[free],
                               jac=lambda z: ri(evaluate(z)[1][:, free]) / SCALE,
                               bounds=(low[free], high[free]), xtol=1e-13, ftol=1e-13,
                               gtol=1e-14, max_nfev=BUDGET // n_unique - 1)
        estimate[free] = result.x
        nfev = int(result.nfev)
    after = diagnostic(case, model, fit, action, new, boxes, estimate)
    return dict(before=before, after=after, state="completed", fit_factors=model.factorizations,
                midpoint_factors=midpoint_factors, hull=[low.tolist(), high.tolist()], fit_nfev=nfev,
                fit_cpu_s=process_time() - cpu0, fit_wall_s=perf_counter() - wall0)


def source_files(validation):
    if validation:
        return load_validation_cases(), load_validation_reads(), ROOT / "results/validation_regions"
    source = ROOT / "results"
    with np.load(source / "measure_new_reads.npz") as all_reads:
        cases = load_cases()
        reads = {case.case_id + "/split_0.02": all_reads[case.case_id + "/split_0.02"] for case in cases}
    return cases, reads, source / "posterior_joint"


def path(name, validation):
    suffix = "_validation" if validation else ""
    return ROOT / "results" / f"measure_region_refine_{name}{suffix}.json"


def fit_all(validation):
    cases, reads, region_dir = source_files(validation)
    records = []
    with threadpool_limits(limits=1):
        for case in cases:
            archive = region_dir / f"{case.case_id}.npz"
            metadata = json.loads(archive.with_suffix(".json").read_text())
            posterior = metadata.get("posterior", metadata)
            upstream = ROOT / "results/material_covers_seeded" / f"{case.case_id}.json"
            old_cpu = (json.loads(upstream.read_text()).get("cpu_s") if upstream.exists()
                       else metadata.get("old", {}).get("cpu_s"))
            with np.load(archive) as covered:
                for stage, action in STAGES:
                    new = measured_reads(reads, case.case_id, action)
                    boxes = covered[stage]
                    selected = fit_in_hull(case, action, new, boxes)
                    t0 = process_time()
                    comparator = solve(case, action, new, [], "ordinary8", BUDGET)
                    comparator["fit_cpu_s"] = process_time() - t0
                    records.append(dict(case_id=case.case_id, cohort=case.cohort, stage=stage,
                                        action=action.name, boxes=len(boxes),
                                        posterior_cpu_s=posterior[stage]["cpu_s"],
                                        old_refinement_cpu_s=old_cpu,
                                        old_seed_cpu_s=metadata.get("r10", {}).get("cpu_s"),
                                        posterior_D_upper=posterior[stage]["D_upper"],
                                        ordinary8=comparator, **selected))
            print(case.case_id, [(r["stage"], r["after"]["radius_upper"] if r["after"] else None)
                                 for r in records[-2:]], flush=True)
    path("predictions", validation).write_text(json.dumps(records, indent=2, allow_nan=False) + "\n")


def score_all(validation):
    records = json.loads(path("predictions", validation).read_text())
    cases, reads, _ = source_files(validation)
    public = {case.case_id: case for case in cases}
    if validation:
        truths = {trial.pair_id: trial.material for trial in private_trials()}
        directions = {}
    else:
        generators = archived_generator()
        truths = {key: generator.truth for key, generator in generators.items()}
        directions = {key: generator.alternatives for key, generator in generators.items()}
    summaries = []
    for row in records:
        name = row["case_id"]
        case, action = public[name], ACTIONS[row["action"]]
        new = measured_reads(reads, name, action)
        truth = truths[name.rsplit("-", 1)[0] if validation else name]
        branches = directions.get(name, np.empty((0, 6)))
        model = FrequencyTM(case.geometry)
        for stage in ("before", "after"):
            point = row[stage]
            if point is not None:
                point.update(physical_score(case, model, action, new,
                                            np.asarray(point["material"]), truth, branches))
        comparison = row["ordinary8"]
        comparison.update(physical_score(case, model, action, new,
                                         np.asarray(comparison["prediction"]), truth, branches))
    for cohort, (stage, _) in product(sorted({r["cohort"] for r in records}), STAGES):
        group = [r for r in records if r["cohort"] == cohort and r["stage"] == stage]
        retained = [r for r in group if r["after"] is not None]
        for name in ("before", "after", "ordinary8"):
            chosen = group if name == "ordinary8" else retained
            results = [r[name] for r in chosen]
            if not results:
                continue
            rmse = np.array([r["rmse"] for r in results])
            fields = dict(cohort=cohort, stage=stage, method=name, total_cases=len(group),
                          nonempty=len(retained), common_compatible=sum(r["common_data_compatible"] for r in results),
                          rmse_median=float(np.median(rmse)), rmse_p90=float(np.quantile(rmse, .9)),
                          rmse_max=float(max(rmse)), rmse_le_point5=int(np.count_nonzero(rmse <= .5)))
            if name != "ordinary8":
                fields["radius_le_point5"] = sum(r["radius_upper"] <= .5 for r in results)
                fields["radius_median"] = float(np.median([r["radius_upper"] for r in results]))
                fields["radius_le_point5_and_compatible"] = sum(
                    r["radius_upper"] <= .5 and r["common_data_compatible"] for r in results)
            if name != "before":
                timing = chosen if name == "after" else results
                fields.update(fit_cpu_median_s=float(np.median([r["fit_cpu_s"] for r in timing])),
                              factor_mean=float(np.mean([r["fit_factors"] for r in timing])))
            summaries.append(fields)
    path("scores", validation).write_text(json.dumps(records, indent=2, allow_nan=False) + "\n")
    path("summary", validation).write_text(json.dumps(summaries, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("fit", "score"))
    parser.add_argument("--validation", action="store_true")
    args = parser.parse_args()
    (fit_all if args.stage == "fit" else score_all)(args.validation)
