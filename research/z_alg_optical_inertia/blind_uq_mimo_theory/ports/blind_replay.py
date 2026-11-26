"""Replay sealed image-generated physical ambiguity branches on paired arrays."""
import csv
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np

from model import MAIN_STRENGTH, DiskGreen, Green, operator, network_terms, physical_solve
from sweep import tag
from witness import pair_grid

HERE = Path(__file__).resolve().parent


def main():
    started = perf_counter(), process_time()
    grid = DiskGreen(32, 400)
    op = operator(grid, network_terms(grid, 8, MAIN_STRENGTH), grid.dense())
    green = Green(op)
    rows = []
    for directory in sorted((HERE.parent/"imaging"/"results").glob("n32_c*")):
        receipt_path = directory/"b1.json"
        if not receipt_path.exists():
            continue
        receipt = json.loads(receipt_path.read_text())
        branches = [(roi,v) for roi,v in receipt["witness_fits"].items() if v["within_ball_numerical"]]
        if not branches:
            continue
        oid = directory.name
        with np.load(HERE/"public"/"cases"/f"{oid}.npz") as obs, np.load(directory/"b1_controls.npz") as controls:
            y = obs["Y"]
            y_hold = obs["Y_hold"]
            r = float(obs["radius"])
            hold_r = float(obs["radius_hold"])
            sigma = float(obs["sigma"])
            alternative = controls["physical_chi"].ravel()
        J, cost = physical_solve(alternative, green.apply, op["E0"], green.diagonal, op["G"])
        candidate_r = np.linalg.norm(op["S"]@J-y)/r
        if candidate_r > 1:
            continue
        candidate_hold_r = np.linalg.norm(op["S_hold"]@J-y_hold)/hold_r
        with np.load(directory/"baseline.npz") as base, np.load(directory/"b1_predictions.npz") as candidates:
            masks = base["rois"]
            for roi, _ in branches:
                absent = candidates["x_"+roi].ravel()
                mask = masks[int(roi.rsplit("_",1)[1])-1].ravel()
                if np.max(absent[mask]) != 0 or alternative[mask].mean() == 0:
                    continue
                Ja, ca = physical_solve(absent, green.apply, op["E0"], green.diagonal, op["G"])
                absent_r = np.linalg.norm(op["S"]@Ja-y)/r
                if absent_r > 1:
                    continue
                absent_hold_r = np.linalg.norm(op["S_hold"]@Ja-y_hold)/hold_r
                replay_rows, arrays = pair_grid(alternative, absent, sigma)
                with np.load(HERE/"public"/"cases"/f"{oid}_mimo.npz") as conditions:
                    for row in replay_rows:
                        key = tag(row["frequency_MHz"], row["q"], row["radiation_scale_ohms"], row["calibration_eta"])
                        observed = conditions[key+"_Y"]
                        radius = float(conditions[key+"_radius"])
                        errors = [np.linalg.norm(arrays[key+f"_F_{mode}_region"]-observed)/radius
                                  for mode in ("with", "without")]
                        rows.append(dict(case_id=oid, roi=roi, roi_cells=int(mask.sum()),
                                         physical_roi_mean_difference=float(alternative[mask].mean()),
                                         baseline_with_region_residual_over_r=candidate_r,
                                         baseline_without_region_residual_over_r=absent_r,
                                         baseline_with_region_hold_residual_over_r=candidate_hold_r,
                                         baseline_without_region_hold_residual_over_r=absent_hold_r,
                                         observed_with_region_residual_over_r=errors[0],
                                         observed_without_region_residual_over_r=errors[1],
                                         **row))
                print(oid, roi, float(alternative[mask].mean()), float(candidate_r), float(absent_r), flush=True)
    if not rows:
        raise RuntimeError("No sealed compatible alternatives with distinct region averages")
    with (HERE/"results"/"blind_pair_array_replay.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    cost = dict(seconds_wall=perf_counter()-started[0], seconds_CPU=process_time()-started[1],
                evaluated_rows=len(rows), objects=len({r["case_id"] for r in rows}),
                interval_claim="none here; point distances recomputed with independent bare VIE and exact network")
    (HERE/"results"/"blind_pair_array_replay_cost.json").write_text(json.dumps(cost, indent=2))


if __name__ == "__main__":
    main()
