"""Compare posterior refinement with and without the intermediate old-domain pass."""
import argparse
from contextlib import nullcontext
import json
from time import perf_counter, process_time

import numpy as np
from threadpoolctl import threadpool_limits

from data14 import ROOT, FrequencyTM, load_cases, measured_reads
from material_controls import GRID
from measure_acquire import archived_generator
from measure_region_refine import STAGES, fit_in_hull
from measure_score import physical_score
from measure_validation import load_validation_cases, load_validation_reads, private_trials
from pivot_experiment import run_arm
from posterior import Joint, group_cover, hull_diam, work

RESULTS = ROOT / "results"
VISITS = 4500
GOAL = .5
FOLDER = RESULTS / "posterior_cost_ablation"


def load(path):
    return json.loads((RESULTS / path).read_text())


def dyadic_seed(boxes, bounds):
    low, high = bounds
    # Same outward grid conversion as the old-domain pass.
    grid = (boxes - low) * (GRID / (high - low))
    lower = np.maximum(0, np.floor(grid[:, 0]).astype(np.int64) - 1)
    upper = np.minimum(GRID, np.ceil(grid[:, 1]).astype(np.int64) + 1)
    return np.stack((lower, upper), axis=1)


def qualifies(point):
    return point is not None and point["common_data_compatible"] and point["radius_upper"] <= GOAL


def apply_stage(case, action, reads, seed):
    t0, c0 = perf_counter(), process_time()
    # The only change from the stored posterior_joint route is the old seed.
    joint = Joint(case, action, reads, numerator=True)
    boxes, stats = work(joint, seed, .085, VISITS)
    diameter = hull_diam(boxes)
    fit = fit_in_hull(case, action, reads, boxes)
    row = dict(**stats, retained=len(boxes), diameter_upper=diameter, energy=action.energy,
               exposure_blocks=action.exposures, fit=fit, certified=qualifies(fit["after"]),
               cpu_s=process_time() - c0, wall_s=perf_counter() - t0)
    return row, boxes


def run(cohorts, validation):
    FOLDER.mkdir(exist_ok=True)
    seeds_meta = None if validation else load("pivot_old.json")["cases"]
    seed_archive = nullcontext(None) if validation else np.load(RESULTS / "pivot_old_regions.npz")
    new_archive = (nullcontext(load_validation_reads()) if validation
                   else np.load(RESULTS / "measure_new_reads.npz"))
    with (seed_archive as seeds, new_archive as observations,
          threadpool_limits(limits=1)):
        for case in load_validation_cases() if validation else load_cases():
            if cohorts and case.cohort not in cohorts:
                continue
            if validation:
                # The independent controls archived the refined old boxes only.
                seed, saved = run_arm(case, "r10_scalar", .70, 4)
            else:
                seed = seeds_meta[case.case_id]["methods"]["r10_scalar"]
            seed_t0, seed_c0 = perf_counter(), process_time()
            if not validation:
                saved = seeds[f"{case.case_id}/r10_scalar"]
            cover = group_cover(dyadic_seed(saved, case.box), min(96, len(saved)))
            seed_setup_wall = perf_counter() - seed_t0
            seed_setup_cpu = process_time() - seed_c0
            stages, retained = {}, {}
            for stage, action in STAGES:
                reads = measured_reads(observations, case.case_id, action)
                stages[stage], retained[stage] = apply_stage(case, action, reads, cover)
                if stages[stage]["certified"]:
                    break
            total_cpu = seed["cpu_s"] + seed_setup_cpu + sum(v["cpu_s"] for v in stages.values())
            total_wall = seed["wall_s"] + seed_setup_wall + sum(v["wall_s"] for v in stages.values())
            record = dict(case_id=case.case_id, cohort=case.cohort, seed_boxes=len(saved),
                          grouped_seed_boxes=len(cover), seed_cpu_s=seed["cpu_s"],
                          seed_wall_s=seed["wall_s"], seed_setup_cpu_s=seed_setup_cpu,
                          seed_setup_wall_s=seed_setup_wall, visit_cap_per_stage=VISITS,
                          target=GOAL, total_cpu_s=total_cpu, total_wall_s=total_wall, stages=stages)
            np.savez_compressed(FOLDER / f"{case.case_id}.npz", **retained)
            (FOLDER / f"{case.case_id}.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
            print(case.case_id, [(name, row["certified"], row["visited"],
                                  round(row["fit"]["after"]["radius_upper"], 4) if row["fit"]["after"] else None)
                                 for name, row in stages.items()], round(total_cpu, 3), flush=True)


def add_metrics(record, case, action, new, truth, alternatives):
    model = FrequencyTM(case.geometry)
    for name in ("before", "after"):
        row = record[name]
        if row is not None:
            row.update(physical_score(case, model, action, new, np.asarray(row["material"]),
                                      truth, alternatives))


def reference_stage(case, reference, seed=None):
    if seed is None:
        validation = load(f"validation_regions/{case}.json")
        seed, old, posterior = (validation[key] for key in ("r10", "old", "posterior"))
    else:
        old = load(f"material_covers_seeded/{case}.json")
        posterior = load(f"posterior_joint/{case}.json")
    selected = {}
    for stage, action in STAGES:
        row, bound = reference[case, stage], posterior[stage]
        selected[stage] = dict(fit={k: row[k] for k in ("before", "after", "fit_factors")},
                               certified=qualifies(row["after"]), visited=bound["visited"],
                               energy=action.energy, exposure_blocks=action.exposures,
                               cpu_s=bound["cpu_s"] + row["fit_cpu_s"],
                               wall_s=bound["wall_s"] + row["fit_wall_s"])
        if selected[stage]["certified"]:
            break
    return dict(total_cpu_s=seed["cpu_s"] + old["cpu_s"] + sum(v["cpu_s"] for v in selected.values()),
                total_wall_s=seed["wall_s"] + old["wall_s"] + sum(v["wall_s"] for v in selected.values()),
                seed_cpu_s=seed["cpu_s"], seed_wall_s=seed["wall_s"],
                intermediate_old_cpu_s=old["cpu_s"], intermediate_old_wall_s=old["wall_s"],
                stages=selected)


def summary(group, key):
    seq = [row[key] for row in group]
    final = [list(v["stages"].values())[-1] for v in seq]
    after = [v["fit"]["after"] for v in final if v["fit"]["after"] is not None]
    rms = [row["rmse"] for row in after]
    return dict(cases=len(group), half_certified=sum(v["stages"]["first_half"]["certified"] for v in seq),
                final_certified=sum(v["certified"] for v in final),
                final_common_compatible=sum(r["common_data_compatible"] for r in after),
                rms_over_radius=sum(r["rmse"] > r["radius_upper"] for r in after),
                rmse_median=float(np.median(rms)) if rms else None,
                rmse_p90=float(np.quantile(rms, .9)) if rms else None,
                rmse_max=max(rms, default=None),
                radius_median=float(np.median([r["radius_upper"] for r in after])) if after else None,
                mean_energy=float(np.mean([r["energy"] for r in final])),
                mean_blocks=float(np.mean([r["exposure_blocks"] for r in final])),
                mean_fit_factors=float(np.mean([
                    sum(s["fit"]["fit_factors"] for s in v["stages"].values()) for v in seq])),
                cpu_median_s=float(np.median([r["total_cpu_s"] for r in seq])),
                wall_median_s=float(np.median([r["total_wall_s"] for r in seq])),
                cpu_total_s=float(sum(r["total_cpu_s"] for r in seq)),
                wall_total_s=float(sum(r["total_wall_s"] for r in seq)))


def score():
    generators = archived_generator()
    reference = {(r["case_id"], r["stage"]): r for r in load("measure_region_refine_scores.json")}
    seed_meta = load("pivot_old.json")["cases"]
    result = []

    def append_case(case, observations, truth, alternatives, ref, seed):
        path = FOLDER / f"{case.case_id}.json"
        if not path.exists():
            return
        record = json.loads(path.read_text())
        for stage, action in STAGES:
            if stage in record["stages"]:
                add_metrics(record["stages"][stage]["fit"], case, action,
                            measured_reads(observations, case.case_id, action), truth, alternatives)
        result.append(dict(case_id=case.case_id, cohort=case.cohort, r10_seed_only=record,
                           two_pass=reference_stage(case.case_id, ref, seed)))

    with np.load(RESULTS / "measure_new_reads.npz") as observations:
        for case in load_cases():
            g = generators[case.case_id]
            append_case(case, observations, g.truth, g.alternatives, reference,
                        seed_meta[case.case_id]["methods"]["r10_scalar"])
    if any(FOLDER.glob("val6-*.json")):
        truth = {trial.pair_id: trial.material for trial in private_trials()}
        reads = load_validation_reads()
        validation_ref = {(r["case_id"], r["stage"]): r
                          for r in load("measure_region_refine_scores_validation.json")}
        for case in load_validation_cases():
            append_case(case, reads, truth[case.case_id.rsplit("-", 1)[0]], np.empty((0, 6)),
                        validation_ref, None)
    summaries = [dict(cohort=cohort, **{name: summary(group, name)
                                       for name in ("r10_seed_only", "two_pass")})
                 for cohort in ("spectral", "natural", "aperture", "low", "high")
                 if (group := [r for r in result if r["cohort"] == cohort])]
    (RESULTS / "posterior_cost_ablation_scores.json").write_text(json.dumps(result, indent=2) + "\n")
    (RESULTS / "posterior_cost_ablation_summary.json").write_text(json.dumps(summaries, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("fit", "score"))
    parser.add_argument("--cohort", action="append", choices=("spectral", "natural", "aperture", "low", "high"))
    parser.add_argument("--validation", action="store_true")
    args = parser.parse_args()
    if args.stage == "fit":
        run(args.cohort, args.validation)
    else:
        score()
