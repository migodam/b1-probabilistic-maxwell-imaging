"""Restart identical physical inversion arms at four disclosed endpoints and an older bent strip."""

import argparse
import csv
import json

import numpy as np
from threadpoolctl import threadpool_limits

from cohorts import (
    historical_study,
    load_cohort,
    load_r3_curved_bar,
    load_truth,
    load_r3_curved_bar_truth,
)
from current_cycle import OUT, prepare, run_case, use_archived_objective
from score_materials import legacy_r3_material_metrics, material_metrics

ARMS = (
    "voxel_only",
    "lagged48",
    "shared48",
    "shared36_tangent_soft",
    "shared_transport",
)
R10_STARTS = (
    (1, "cie_single"),
    (16, "cie_single"),
    (4, "matrix_green19"),
    (22, "matrix_green19"),
    (4, "cie_single"),
)


def fit(original=False):
    suffix = "_original" if original else ""
    model, cases = load_cohort("r10_twod")
    if original:
        use_archived_objective(model)
    work = []
    for i, method in R10_STARTS:
        case = cases[i]
        path = (
            historical_study(10)
            / "RESULTS"
            / "twod"
            / f"{case['id']}__{method}__noise_ball_99"
            / "prediction.npz"
        )
        with np.load(path) as saved:
            work.append(("r10_twod", model, case, method, saved["m"].ravel()))
    model3, (case3,) = load_r3_curved_bar()
    file3 = historical_study(2) / "predictions" / case3["id"] / "native.npz"
    with np.load(file3) as saved:
        work.append(("r3_curved_bar", model3, case3, "cie_dual", saved["m"].ravel()))
    work.append(
        ("r3_curved_bar", model3, case3, "background", np.full(len(model3["G"]), 0.1))
    )
    keys = []
    images = []
    logs = []
    for group, model, case, start, m in work:
        if original and group != "r10_twod":
            continue
        bank = prepare(model)
        for arm in ARMS:
            x, log = run_case(model, case, bank, arm, 10, 100, 60, initial=m)
            logs.append(
                dict(
                    original_prior=original,
                    cohort=group,
                    start=start,
                    setup_s=bank["setup_s"],
                    **log,
                )
            )
            keys.append([group, case["id"], start, arm])
            images.append([m, *x])
    (OUT / f"current_restarts{suffix}.json").write_text(
        json.dumps(logs, indent=1) + "\n"
    )
    np.savez_compressed(OUT / f"current_restarts{suffix}.npz", key=keys, images=images)


def score(original=False):
    suffix = "_original" if original else ""
    rows = []
    with np.load(OUT / f"current_restarts{suffix}.npz") as saved:
        for (group, cid, start, arm), images in zip(saved["key"], saved["images"]):
            if group == "r3_curved_bar":
                truth, support = load_r3_curved_bar_truth()
            else:
                truth = load_truth(group, cid)
            for stage, m in zip(("start", "prefix", "final"), images):
                if group == "r3_curved_bar":
                    scores = legacy_r3_material_metrics(m, truth, support)
                    nrmse, iou = scores.values()
                else:
                    scores = material_metrics(m, truth)
                    nrmse, iou = scores["amplitude_NRMSE"], scores["fixed_support_IoU"]
                rows.append(
                    dict(
                        cohort=group,
                        case=cid,
                        start=start,
                        arm=arm,
                        stage=stage,
                        nrmse=nrmse,
                        iou=iou,
                    )
                )
    with (OUT / f"current_restart_scores{suffix}.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--score", action="store_true")
    p.add_argument("--original", action="store_true")
    args = p.parse_args()
    with threadpool_limits(1):
        if args.score:
            score(args.original)
        else:
            fit(args.original)
