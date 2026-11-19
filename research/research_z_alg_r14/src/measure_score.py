"""Score completed predictions against withheld materials and the shared-error contract."""
import argparse
from collections import defaultdict
from itertools import combinations
import json

import numpy as np

from data14 import ACTIONS, ROOT, FrequencyTM, load_cases, measured_reads
from measure_acquire import archived_generator


def gain_disk(response, data, radius, prior=None):
    """A normalized least-squares sublevel set in u = gain - 1 is exactly a disk."""
    residual = (data - response) / radius
    direction = response / radius
    norm2 = np.vdot(direction, direction).real + (0. if prior is None else prior ** -2)
    if norm2 == 0.:
        return (0., np.inf) if np.linalg.norm(residual) <= 1. else None
    center = np.vdot(direction, residual) / norm2
    minimal = np.linalg.norm(residual - center * direction) ** 2
    if prior is not None:
        minimal += abs(center / prior) ** 2
    if minimal > 1. + 1e-12:
        return None
    return center, np.sqrt(max(0., 1. - minimal) / norm2)


def common_point(disks):
    finite = [(c, r) for c, r in disks if np.isfinite(r)]
    if not finite:
        return True
    candidates = [center for center, _ in finite]
    tol = 32 * np.finfo(float).eps
    for (c1, r1), (c2, r2) in combinations(finite, 2):
        distance = abs(c2 - c1)
        if distance > 0. and abs(r1 - r2) - tol <= distance <= r1 + r2 + tol:
            a = (r1 ** 2 - r2 ** 2 + distance ** 2) / (2 * distance)
            h = np.sqrt(max(0., r1 ** 2 - a ** 2))
            direction = (c2 - c1) / distance
            candidates.extend([c1 + (a + 1j * h) * direction, c1 + (a - 1j * h) * direction])
    return bool(any(all(abs(v - center) <= radius + tol for center, radius in finite)
                    for v in candidates))


def compatibility(case, action, new, predict):
    old = gain_disk(predict[0], case.old, case.additive_radius, case.gain_radius)
    if old is None:
        return False, False
    disks = [old]
    if case.hard_gain_radius is not None:
        disks.append((0., case.hard_gain_radius))
    old_valid = common_point(disks)
    if not action.exposures:
        return old_valid, old_valid
    weight = np.sqrt(action.energies)[:, None]
    extra = gain_disk((weight * predict[1:]).ravel(), (weight * new).ravel(),
                      case.new_noise_radius(action))
    return old_valid, old_valid and extra is not None and common_point([*disks, extra])


def global_ssim(a, b, data_range):
    mu_a, mu_b = np.mean(a), np.mean(b)
    var_a, var_b = np.var(a), np.var(b)
    cov = np.mean((a - mu_a) * (b - mu_b))
    c1, c2 = (.01 * data_range) ** 2, (.03 * data_range) ** 2
    return float((2 * mu_a * mu_b + c1) * (2 * cov + c2)
                 / ((mu_a ** 2 + mu_b ** 2 + c1) * (var_a + var_b + c2)))


def physical_score(case, model, action, new, material, truth, branches):
    predict = np.array([model.response(material, d)[0] for d in (0., *action.deltas)])
    old_valid, joint_valid = compatibility(case, action, new, predict)
    absolute = np.linalg.norm(material - truth)
    truth_high, estimated_high = truth >= 10., material >= 10.
    union = int(np.count_nonzero(truth_high | estimated_high))
    contrast_iou = float(np.count_nonzero(truth_high & estimated_high) / union) if union else 1.
    return dict(
        rmse=float(absolute / np.sqrt(truth.size)), nrmse=float(absolute / np.linalg.norm(truth)),
        iou_contrast_ge_10=contrast_iou, global_cell_ssim=global_ssim(material, truth, np.ptp(case.box)),
        wrong_side=(bool(np.any(np.linalg.norm(branches - material, axis=1) < absolute))
                    if len(branches) else None),
        old_compatible=old_valid, common_data_compatible=joint_valid,
        rms_over_point1=bool(absolute > .1 * np.sqrt(truth.size)),
    )


def summarize(items):
    values = np.asarray([item["rmse"] for item in items])
    relative = np.asarray([item["nrmse"] for item in items])
    wrong = [item["wrong_side"] for item in items if item["wrong_side"] is not None]
    design = [item["design_fallback"] for item in items if item["design_fallback"] is not None]
    return dict(
        cases=len(items), rmse_median=float(np.median(values)),
        rmse_p90=float(np.quantile(values, .9)), rmse_max=float(np.max(values)),
        nrmse_median=float(np.median(relative)), nrmse_p90=float(np.quantile(relative, .9)),
        global_cell_ssim_median=float(np.median([item["global_cell_ssim"] for item in items])),
        contrast_iou_median=float(np.median([item["iou_contrast_ge_10"] for item in items])),
        rms_over_point1=sum(item["rms_over_point1"] for item in items),
        common_compatible=sum(item["common_data_compatible"] for item in items),
        truth_in_contract=sum(item["truth_in_common_contract"] for item in items),
        wrong_side=sum(wrong), wrong_side_denominator=len(wrong),
        fallback=sum(design), fallback_denominator=len(design),
        bandwidth_pct=float(np.mean([item["bandwidth_pct"] for item in items])),
        mean_energy=float(np.mean([item["energy"] for item in items])),
        mean_exposures=float(np.mean([item["exposures"] for item in items])),
        all_factors_mean=float(np.mean([item["total_factors"] for item in items])),
        rhs_mean=float(np.mean([item["total_rhs"] for item in items])),
        design_factors_mean=float(np.mean([item["design_factors"] for item in items])),
        fit_factors_mean=float(np.mean([item["fit_factors"] for item in items])),
        polynomial_paths_mean=float(np.mean([item["polynomial_paths"] for item in items])),
        total_time_median_s=float(np.median([item["total_wall_s"] for item in items])),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default="")
    args = parser.parse_args()
    fixed = json.loads((ROOT / "results" / f"measure_predictions{args.tag}.json").read_text())
    rows = fixed["rows"]
    generators = archived_generator()
    public = {case.case_id: case for case in load_cases()}
    models = {name: FrequencyTM(case.geometry) for name, case in public.items()}
    with np.load(ROOT / "results/measure_new_reads.npz") as reads:
        data = {key: reads[key] for key in reads.files}
    score_cache, truth_cache = {}, {}
    summary = defaultdict(list)
    for row in rows:
        name, action_name = row["case_id"], row["action"]
        case, model, action = public[name], models[name], ACTIONS[action_name]
        generator = generators[name]
        new = measured_reads(data, name, action)
        key = (name, action_name, tuple(row["prediction"]))
        if key not in score_cache:
            score_cache[key] = physical_score(case, model, action, new, np.asarray(row["prediction"]),
                                              generator.truth, generator.alternatives)
        if (name, action_name) not in truth_cache:
            predict = np.array([model.response(generator.truth, d)[0] for d in (0., *action.deltas)])
            truth_cache[name, action_name] = compatibility(case, action, new, predict)[1]
        row.update(score_cache[key], energy=action.energy, exposures=action.exposures,
                   bandwidth_pct=100 * action.bandwidth, max_offset_pct=100 * action.max_offset,
                   truth_in_common_contract=truth_cache[name, action_name],
                   total_factors=row["design_factors"] + row["fit_factors"],
                   total_rhs=row["design_rhs"] + row["fit_rhs"],
                   total_wall_s=row["design_wall_s"] + row["fit_wall_s"])
        title = action_name if row["label"] == "sweep" else row["label"]
        summary[case.cohort, row["label"], title, row["mode"]].append(row)
    all_summary = [{"cohort": key[0], "label": key[1], "action_or_policy": key[2], "mode": key[3],
                    **summarize(items)} for key, items in sorted(summary.items())]
    out = ROOT / "results"
    (out / f"measure_summary{args.tag}.json").write_text(json.dumps(all_summary, indent=2) + "\n")
    (out / f"measure_scores{args.tag}.json").write_text(json.dumps(rows, indent=2) + "\n")
    if not args.tag:
        witnesses = {}
        for name in ("confirm6-00", "confirm6-05"):
            case, model = public[name], models[name]
            generator = generators[name]
            witnesses[name] = {}
            for action_name in ("old_only", "one_minus_0.02", "one_plus_0.02", "split_0.02",
                                "uniform_0.02", "split_0.04", "split_0.08", "uniform_0.08"):
                action = ACTIONS[action_name]
                records = []
                for branch in (generator.truth, *generator.alternatives):
                    response = np.array([model.response(branch, d)[0] for d in (0., *action.deltas)])
                    new = data[f"{name}/{action_name}"]
                    old_ok, joint_ok = compatibility(case, action, data[f"{name}/{action_name}"],
                                                     response)
                    at_one_gain = new - generator.gain * response[1:]
                    old_error = case.old - generator.gain * response[0]
                    old_score = (np.linalg.norm(old_error) / case.additive_radius) ** 2
                    old_score += abs((generator.gain - 1.) / case.gain_radius) ** 2
                    new_score = (float(np.linalg.norm(np.sqrt(action.energies)[:, None] * at_one_gain)
                                       / case.new_noise_radius(action)) if action.exposures else None)
                    records.append(dict(truth_rmse=float(np.linalg.norm(branch - generator.truth) / np.sqrt(6)),
                                        old_compatible=old_ok, common_data_compatible=joint_ok,
                                        old_ellipsoid_at_archived_gain=float(old_score),
                                        new_radius_ratio_at_archived_gain=new_score))
                witnesses[name][action_name] = records
        (out / "measure_witnesses.json").write_text(json.dumps(witnesses, indent=2) + "\n")


if __name__ == "__main__":
    main()
