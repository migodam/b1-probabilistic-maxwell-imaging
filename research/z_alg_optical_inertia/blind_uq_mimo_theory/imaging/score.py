"""Post-seal material scoring, isolated from all candidate-producing commands."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from shared import HERE, OUT, elapsed, now, write_json

import numpy as np
from scipy.ndimage import gaussian_filter, label

SUPPORT_CUTOFF = .25
COMPONENT_LOW_OVERLAP = .1


def materials(path, baseline):
    yield "baseline", baseline
    if (path/"continuation.npz").is_file():
        with np.load(path/"continuation.npz", allow_pickle=False) as z:
            yield "tighter_dual", z["chi"]
    if (path/"b1_predictions.npz").is_file():
        with np.load(path/"b1_predictions.npz", allow_pickle=False) as z:
            yield "b1", z["chi"]
    if (path/"b1_controls.npz").is_file():
        with np.load(path/"b1_controls.npz", allow_pickle=False) as z:
            yield "b1_cie_control", z["cie_chi"]
            yield "b1_wall_control", z["physical_chi"]
    if (path/"b2_predictions.npz").is_file():
        with np.load(path/"b2_predictions.npz", allow_pickle=False) as z:
            for k in z.files:
                if k not in ("dark_modes","random_null_modes"):
                    yield "b2_"+k, z[k]


def metrics(x, truth, threshold, minimum):
    inside = truth > SUPPORT_CUTOFF
    outside = ~inside
    peak = gaussian_filter(x, sigma=x.shape[0]/32.)
    pred = peak >= threshold
    labels, count = label(pred)
    components = [labels == index for index in range(1, count+1)]
    components = [component for component in components if component.sum() >= minimum]
    overlaps = [float(inside[component].mean()) for component in components]
    return dict(
        material_relative_L2=float(np.linalg.norm(x-truth)/np.linalg.norm(truth)),
        inside_relative_L2=float(np.linalg.norm(x[inside]-truth[inside])/np.linalg.norm(truth[inside])),
        background_RMS=float(np.sqrt(np.mean((x[outside]-truth[outside])**2))),
        background_false_positive_fraction=float(pred[outside].mean()),
        segmented_component_count=len(components),
        zero_truth_overlap_components=sum(ov == 0 for ov in overlaps),
        below_tenth_truth_overlap_components=sum(ov < COMPONENT_LOW_OVERLAP for ov in overlaps),
        component_truth_overlap_fractions=overlaps,
        truth_support_iou=float((inside & pred).sum()/max((inside | pred).sum(), 1)),
    )


def score_case(case, directory, path):
    # Caller passes the explicit sealed case ID. No glob or prior load exists.
    with np.load(directory/f"{case}.npz", allow_pickle=False) as secret:
        truth = secret["material"]
    with np.load(path/"baseline.npz", allow_pickle=False) as release:
        x = release["chi"]
        rois = release["rois"]
    b = json.loads((path/"baseline.json").read_text())
    seg = b["segmentation"]
    held = json.loads((path/"held.json").read_text())["predictions"]
    data = {}
    for method, candidate in materials(path, x):
        if candidate.shape != truth.shape:
            raise ValueError("No scoring interpolation between native grids")
        data[method] = metrics(candidate, truth, seg["threshold"], seg["minimum_cells"])
        if method in held:
            data[method].update(held[method])
    roi_table = []
    for i, region in enumerate(rois, 1):
        roi_table.append(dict(id=f"roi_{i:02d}", cells=int(region.sum()),
                              true_mean=float(truth[region].mean()),
                              true_max=float(truth[region].max()),
                              true_support_fraction=float((truth[region] > SUPPORT_CUTOFF).mean())))
    return dict(case=case, side=x.shape[0], method_metrics=data, roi_truth_audit=roi_table)


def score_dataset(manifest, directory, output_stem="material_scores"):
    freeze = json.loads(manifest.read_text())
    if freeze.get("sealed") is not True or not freeze.get("cases"):
        raise ValueError("Require sealed manifest and explicit case IDs before any truth read")
    cases = freeze["cases"]
    stationary = set(freeze.get("stationary", []))
    for case in cases:
        for path in (OUT/case, OUT/case/"stationary") if case in stationary else (OUT/case,):
            required = ("baseline.npz","baseline.json","b1_predictions.npz","b2_predictions.npz",
                        "b1_controls.npz","held.json")
            if not all((path/p).is_file() for p in required):
                raise FileNotFoundError(f"Missing sealed artifact for {path}")
    start = now()
    output = {}
    for case in cases:
        output[case] = {"dual":score_case(case,directory,OUT/case)}
        if case in stationary:
            output[case]["stationary"] = score_case(case,directory,OUT/case/"stationary")
    write_json(OUT/(output_stem+".json"), dict(
        material_support_cutoff=SUPPORT_CUTOFF, low_overlap_fraction=COMPONENT_LOW_OVERLAP,
        write_seal_path=str(manifest.relative_to(HERE.parent)) if manifest.is_relative_to(HERE.parent) else str(manifest),
        objects=output, **elapsed(start)))
    scalar = ("material_relative_L2","inside_relative_L2","background_RMS",
              "background_false_positive_fraction","segmented_component_count",
              "zero_truth_overlap_components","below_tenth_truth_overlap_components",
              "truth_support_iou","held_relative_residual","residual_over_radius")
    with (OUT/(output_stem+".csv")).open("w") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["case","baseline_family","method",*scalar],
                                extrasaction="ignore")
        writer.writeheader()
        for case, families in output.items():
            for family, item in families.items():
                for method, metric in item["method_metrics"].items():
                    writer.writerow(dict(case=case,baseline_family=family,method=method,**metric))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sealed-manifest", type=Path, required=True)
    parser.add_argument("--private-truth-directory", type=Path, required=True)
    parser.add_argument("--output-stem", default="material_scores")
    args = parser.parse_args()
    if Path(args.output_stem).name != args.output_stem:
        raise ValueError("Output stem must be one filename component")
    score_dataset(args.sealed_manifest.resolve(), args.private_truth_directory.resolve(), args.output_stem)
