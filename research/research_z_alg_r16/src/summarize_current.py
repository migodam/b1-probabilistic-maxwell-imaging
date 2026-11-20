"""Case-paired current/voxel outcomes. Reference denominators are shared across methods."""

import argparse
import csv
import json
from collections import defaultdict

import numpy as np

from cohorts import load_cohort
from current_cycle import OUT

SAME = ("r10_twod", "r11_independent48", "r12_confirm24")
COHORTS = (*SAME, "r11_crossgrid8", "r11_crossgrid8_matched", "r12_dev6")


def median(rows, col):
    return float(np.median([r[col] for r in rows]))


def by_arm(rows):
    if not rows:
        return {}
    arms = defaultdict(list)
    for r in rows:
        arms[r["arm"]].append(r)
    control = {r["case"]: r for r in arms["voxel_only"]}
    report = {}
    for arm, items in arms.items():
        ref = [control[r["case"]] for r in items]
        errors = np.array([r["amplitude_NRMSE"] for r in items])
        delta = errors - np.array([r["amplitude_NRMSE"] for r in ref])
        report[arm] = dict(
            n=len(items),
            nrmse=median(items, "amplitude_NRMSE"),
            q90=float(np.quantile(errors, 0.9)),
            iou=median(items, "fixed_support_IoU"),
            rmse=median(items, "RMSE"),
            rel_data=median(items, "full_rel"),
            physical_evals=median(items, "total_full_solves"),
            seconds=median(items, "wall_s"),
            cpu=median(items, "cpu_s"),
            prefix_rel=median(items, "prefix_rel"),
            prefix_calls=median(items, "prefix_full_solves"),
            complex_rhs=median(items, "complex_rhs"),
            success=sum(r["material_success"] for r in items),
            noise2=int(sum(r["noise_rms"] <= 2.0 for r in items)),
            better=int(sum(delta < 0)),
            better_by_01=int(sum(delta < -0.01)),
            worse_by_01=int(sum(delta > 0.01)),
            median_paired_delta=float(np.median(delta)),
            prefix_accepted=median(items, "accepted"),
        )
    return report


def main(extended=False, archive_budget=None):
    cohorts = (*SAME, "r11_crossgrid8") if extended else COHORTS
    tags = ("long",) if extended else ("equal", "tangent", "exactdark", "blend")
    suffix = "_180" if extended else ""
    if archive_budget is not None:
        tags = (f"archive{archive_budget}",)
        cohorts = (*SAME, "r11_crossgrid8")
        suffix = f"_archive{archive_budget}"
    rows = []
    for group in cohorts:
        model, cases = load_cohort(group)
        sigma = {c["id"]: c["sigma"] * np.sqrt(c["Y"].size) for c in cases}
        for tag in tags:
            src = OUT / f"current_scores_{tag}_{group}.csv"
            if not src.exists():
                continue
            with src.open() as f:
                scores = {
                    (r["case"], r["arm"]): r
                    for r in csv.DictReader(f)
                    if r["stage"] == "after_voxel_polish"
                }
            logs = json.loads((OUT / f"current_cycle_{tag}_{group}.json").read_text())[
                "rows"
            ]
            for r in logs:
                work = scores[r["case"], r["arm"]]
                for col in (
                    "RMSE",
                    "amplitude_NRMSE",
                    "fixed_support_IoU",
                    "SSIM_fixed_range",
                ):
                    r[col] = float(work[col])
                r["material_success"] = work["material_success"] == "True"
                r["complex_rhs"] = (
                    r["forward_rhs"]
                    + r["gradient_rhs"]
                    + r["transport_rhs"]
                    + r.get("tangent_rhs", 0)
                )
                r["noise_rms"] = r["full_rel"] * r["initial_norm"] / sigma[r["case"]]
                r.pop("trace")
                rows.append(r)
    groups = {g: by_arm([r for r in rows if r["cohort"] == g]) for g in cohorts}
    groups["samegrid96"] = by_arm([r for r in rows if r["cohort"] in SAME])
    with (OUT / f"current_cases{suffix}.csv").open("w") as f:
        keys = sorted(set().union(*(r.keys() for r in rows)))
        writer = csv.DictWriter(f, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (OUT / f"current_summary{suffix}.json").write_text(
        json.dumps(groups, indent=2) + "\n"
    )
    for g in ("samegrid96", "r11_crossgrid8", "r11_crossgrid8_matched"):
        if g not in groups:
            continue
        print(g)
        for a, r in groups[g].items():
            print(
                f"{a:22s} {r['nrmse']:.4f} {r['iou']:.3f} {r['success']} {r['better']} {r['seconds']:.3f} cpu {r['cpu']:.3f} counts {r['better_by_01']},{r['worse_by_01']} p {r['prefix_accepted']}"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--long", action="store_true")
    parser.add_argument("--archive-budget", type=int, choices=(60, 180))
    args = parser.parse_args()
    main(args.long, args.archive_budget)
