"""Score the full-wave point-template predictions saved prior to opening centre labels."""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def main():
    roster = json.loads((HERE / "public" / "acquisition.json").read_text())["observations"]
    ranked = {}
    for key, spec in roster.items():
        if spec["scene"] != "weak":
            continue
        with np.load(RESULTS / f"{key}__b0.03125__primary.npz") as predictions:
            errors = predictions["reference_error"]
        if not np.all(np.isfinite(errors) & (errors >= 0)):
            raise ArithmeticError(f"Incomplete position sweep for {key}")
        with np.load(HERE / "public" / f"{key}.npz") as obs:
            radius = float(obs["r"])
        side = spec["side"]
        order = np.argsort(errors, kind="stable")
        # All ranking inputs are old measurements and full-wave model predictions.
        ranked[key] = (side, order, errors, radius)
    rows = {}
    for key, (side, order, errors, radius) in ranked.items():
        with np.load(HERE / "offline" / f"{key}.npz") as label:
            cx, cy = map(int, label["centre"])
        truth_rank = int(np.flatnonzero(order == cy * side + cx)[0]) + 1
        rows[key] = dict(side=side, rank_of_generating_centre=truth_rank,
                         best_centre_xy=[int(order[0] % side), int(order[0] // side)],
                         minimum_point_residual_over_radius=float(errors[order[0]] / radius),
                         point_templates_inside_ball=int(np.count_nonzero(errors <= radius)))
    summary = {}
    for side in (32, 96):
        group = [v for v in rows.values() if v["side"] == side]
        ratios = [r["minimum_point_residual_over_radius"] for r in group]
        summary[side] = dict(specimens=len(group), centre_ranked_first=sum(r["rank_of_generating_centre"] == 1 for r in group),
                             cohorts_with_template_in_ball=sum(r["point_templates_inside_ball"] > 0 for r in group),
                             minimum_residual_ratio_range=[min(ratios), max(ratios)])
    report = dict(description="Full-wave midpoint contrast of the specified outline, zero exterior. Uses the same measured data and declared material range. One specific material per centre. Exhaustive original-grid point residuals; rankings do not exclude a class.",
                  scope="Readouts already generated in certified_map.py before opening labels; scoring only, no new acquisition",summary=summary,rows=rows)
    (RESULTS / "nominal_template_control.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
