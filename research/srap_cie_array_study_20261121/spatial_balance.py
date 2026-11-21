"""Offline decomposition of paired material improvements. No candidate is reselected."""
from collections import defaultdict
import csv
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
TOL = 1e-10


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    details = []
    for cohort in ("confirmation", "followup"):
        path = RESULTS / f"{cohort}_metrics.csv"
        if not path.exists():
            continue
        grouped = defaultdict(dict)
        for row in read(path):
            if row["budget"] == "30":
                grouped[(row["case"], row["view"])][row["method"]] = row
        for (case, view), arms in grouped.items():
            for label, a, b in (("primary", "dual", "srap"), ("shared_cap", "dual_cap", "srap_cap"),
                                ("cap_Dual", "dual", "dual_cap"), ("cap_SRAP", "srap", "srap_cap")):
                if a not in arms or b not in arms:
                    continue
                base, arm = arms[a], arms[b]
                if not base["rel_L2"] or not arm["rel_L2"]:
                    continue
                def change_sq(key):
                    return float(base[key]) ** 2 - float(arm[key]) ** 2
                details.append(dict(cohort=cohort, case=case, array=base["array"], view=view, comparison=label,
                                    sq_improvement=change_sq("rel_L2"), inside_improvement=change_sq("inside_L2"),
                                    background_improvement=change_sq("background_L2"),
                                    overlap_change=float(arm["half_IoU"])-float(base["half_IoU"]),
                                    base_IoU=float(base["half_IoU"]), arm_IoU=float(arm["half_IoU"]),
                                    base_L2=float(base["rel_L2"]), arm_L2=float(arm["rel_L2"])))
    if not details:
        return
    write(RESULTS / "spatial_balance_objects.csv", details)
    groups = defaultdict(list)
    for row in details:
        groups[tuple(row[k] for k in ("cohort","array","view","comparison"))].append(row)
    out = []
    for key, obs in sorted(groups.items()):
        total, inside, bg, overlap = (np.array([r[k] for r in obs]) for k in
                                    ("sq_improvement","inside_improvement","background_improvement","overlap_change"))
        wins = total > TOL
        out.append(dict(zip(("cohort","array","view","comparison"),key),
                        n=len(obs), total_error_wins=wins.sum(), both_regions_improve=((inside>TOL)&(bg>TOL)).sum(),
                        lower_total_worse_inside=(wins & (inside < -TOL)).sum(),
                        lower_total_no_IoU_gain=(wins & (overlap <= TOL)).sum(),
                        lower_total_both_zero_IoU=sum(r["sq_improvement"]>TOL and r["base_IoU"]==r["arm_IoU"]==0 for r in obs),
                        median_inside_gain=np.median(inside), median_background_gain=np.median(bg),
                        median_IoU_gain=np.median(overlap), max_decomposition_error=np.max(np.abs(total-inside-bg))))
    write(RESULTS / "spatial_balance_summary.csv", out)


if __name__=="__main__":
    main()
