"""Summarize scored observations with pairing by blind material and acquisition."""

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from research.z_alg_xudong_cie.holdout.data import planned_cases, write_json


METRIC_LABELS = {
    "material_rel_L2": "Material relative L2",
    "inside_L2": "Interior contribution",
    "background_L2": "Air contribution",
    "half_IoU": "Half-peak IoU",
    "native_train": "Native fitting residual",
    "held_ports": "Independent port error",
}
CASE_INDEX = {item["case"]: item for item in planned_cases()}


def row_metrics(row, case_id):
    m = row["material96"]
    train = row["native_fitting"].get(case_id, {})
    ports = row["held_transmit_and_receive"] or {}
    return dict(material_rel_L2=m["rel_L2"], inside_L2=m["inside_L2"],
                background_L2=m["background_L2"], half_IoU=m["half_IoU"],
                native_train=train.get("relative_to_measured"),
                held_ports=ports.get("relative_to_measured"))


def summarize(scored, baseline):
    good, failures = {}, []
    for row in scored["rows"]:
        source = row["source"]
        key = source["case"], source["arm"]
        if row["status"] in ("scored", "material_only"):
            if key in good:
                raise ValueError(f"More than one scored map for {key}.")
            good[key] = row
        if row["status"] != "scored":
            failures.append(row)
    aggregate = defaultdict(list)
    comparisons = defaultdict(list)
    plotting = []
    for (case_id, arm), row in good.items():
        case = CASE_INDEX[case_id]
        group = case["array"], case["noise_relative"], arm
        current = row_metrics(row, case_id)
        aggregate[group].append(current)
        baseline_row = good.get((case_id, baseline))
        if arm == baseline or baseline_row is None:
            continue
        reference = row_metrics(baseline_row, case_id)
        delta = {key: value - reference[key] if value is not None and reference[key] is not None
                 else None for key, value in current.items()}
        comparisons[group].append(delta)
        plotting.append(dict(case=case_id, object_id=case["object_id"], aperture=case["array"],
                             noise=case["noise_relative"], arm=arm, deltas=delta))

    def metrics(items):
        rows = {k: [value[k] for value in items if value[k] is not None] for k in METRIC_LABELS}
        return dict(median={k: float(np.median(v)) if v else None for k, v in rows.items()},
                    counts_by_metric={k: len(v) for k, v in rows.items()})

    def group_summary(groups):
        return [dict(aperture=ap, noise_fraction=noise, arm=arm, n_independent_materials=len(items),
                     **metrics(items)) for (ap, noise, arm), items in sorted(groups.items())]
    counts = defaultdict(int)
    for (_, arm) in good:
        counts[arm] += 1
    by_material = defaultdict(list)
    for pair in plotting:
        key = pair["aperture"], pair["arm"], pair["object_id"]
        by_material[key].append(pair)
    pooled = defaultdict(list)
    for (aperture, arm, _), rows in by_material.items():
        if {row["noise"] for row in rows} == {0.02, 0.10}:
            pooled[(aperture, arm)].append(
                {key: (float(np.mean([row["deltas"][key] for row in rows]))
                       if all(row["deltas"][key] is not None for row in rows) else None)
                 for key in METRIC_LABELS}
            )
    return dict(
        scope="paired comparisons on the same object, aperture and noise level",
        baseline=baseline, scored_observations_by_arm=dict(counts), unresolved=failures,
        absolute_medians=group_summary(aggregate), paired_delta_medians=group_summary(comparisons),
        paired_noise_mean_delta_medians=[
            dict(aperture=ap, arm=arm, n_independent_materials=len(items),
                 **metrics(items))
            for (ap, arm), items in sorted(pooled.items())
        ],
        paired_rows=plotting,
    )


def plot_pairs(summary, path):
    if not summary["paired_rows"]:
        return
    comparisons = summary["paired_rows"]
    arms = sorted({r["arm"] for r in comparisons})
    strata = sorted({(r["aperture"], r["noise"]) for r in comparisons})
    fig, axes = plt.subplots(2, 1, figsize=(9, 6.2), sharex=True, layout="constrained")
    colors = plt.get_cmap("tab10").colors
    for ax, metric in zip(axes, ("material_rel_L2", "held_ports"), strict=True):
        for i, arm in enumerate(arms):
            offset = (i - (len(arms)-1)/2) * min(.18, .75/max(len(arms)-1, 1))
            means = []
            for x, (ap, noise) in enumerate(strata):
                y = [r["deltas"][metric] for r in comparisons
                     if r["arm"] == arm and r["aperture"] == ap and r["noise"] == noise
                     and r["deltas"][metric] is not None]
                if y:
                    jitter = np.linspace(-.04, .04, len(y))
                    ax.scatter(x + offset + jitter, y, color=colors[i % len(colors)], s=12, alpha=.5)
                    means.append((x + offset, float(np.median(y))))
            if means:
                ax.scatter(*np.asarray(means).T, marker="D", s=42,
                           color=colors[i % len(colors)], label=arm)
        ax.axhline(0., color="#525d66", linewidth=1)
        ax.set_ylabel(f"Δ {METRIC_LABELS[metric]}")
        ax.grid(axis="y", color="#dfe5eb")
    axes[0].legend(frameon=False, loc="upper center", bbox_to_anchor=(.5, 1.42),
                  ncols=3, fontsize="small")
    axes[1].set_xticks(range(len(strata)), [f"{ap[2:5]}° · {noise:.0%}" for ap, noise in strata])
    axes[1].set_xlabel("Fitting aperture and noise")
    fig.suptitle(f"Paired differences relative to {summary['baseline']}")
    for extension in (".png", ".pdf"):
        fig.savefig(path.with_suffix(extension), dpi=220, facecolor="white")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scores", type=Path)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args()
    report = summarize(json.loads(args.scores.read_text(encoding="utf-8")), args.baseline)
    write_json(args.output, report)
    if args.plot:
        plot_pairs(report, args.plot)
