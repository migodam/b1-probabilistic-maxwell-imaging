"""Figures for the paired physical model-gap experiment."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/model_gap"
COLOURS = ["#167980", "#ab741d", "#7457a8", "#cc6350", "#495b69"]
LAYOUTS = {"PW360_10x20": "Full ring", "LS_LINE_SAME_10x20": "Same-side line"}
CONDITIONS = {"matched_clean": "32² clean", "matched_noisy": "32² noisy",
              "cross_clean": "96→32 clean", "cross_noisy": "96→32 noisy",
              "historical_Dual_30s_light": "Archived Dual"}


def draw_comparison(gaps, scores):
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.6), layout="constrained")
    for column, (layout, name) in enumerate(LAYOUTS.items()):
        rows = [row for row in gaps if row["array"] == layout]
        ids = [row["observation"].split("__")[0] for row in rows]
        x = np.arange(len(ids))
        axis = axes[0, column]
        keys = ("microgeometry_over_delta", "discretization_over_delta", "total_gap_over_delta",
                "total_gap_plus_noise_over_delta")
        labels = ("Structure averaging", "Grid operator", "Total model gap", "Total gap + archived noise")
        for i, (key, label) in enumerate(zip(keys, labels, strict=True)):
            axis.bar(x + 0.18 * (i - 1.5), [row[key] for row in rows], width=0.16,
                     label=label, color=COLOURS[i])
        axis.axhline(1, color="#394653", ls="--", linewidth=1)
        axis.set(title=f"{name} — coarse reference at the ports", ylabel="Frobenius norm / original noise radius",
                 xticks=x, xticklabels=ids, ylim=(0, 3.15))
        target = axes[1, column]
        local = [row for row in scores if row["array"] == layout]
        for i, (key, label) in enumerate(CONDITIONS.items()):
            subset = ([row for row in local if row["method"] == key] if i == 4 else
                      [row for row in local if row["condition"] == key and row["method"] == "physical_LBFGSB_12"])
            target.plot(x, [row["material_rel_l2"] for row in subset], color=COLOURS[i],
                        alpha=0.9, ls="--" if i == 4 else "-", marker="o", ms=4, label=label)
        target.set(title=f"{name} — material error after 12 steps", ylabel="Relative material ℓ² error",
                   xticks=x, xticklabels=ids)
    for axis in axes.flat:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.2)
        axis.set_axisbelow(True)
    axes[0, 0].legend(frameon=False, fontsize=7.5, ncol=2, loc="upper left")
    axes[1, 0].legend(frameon=False, fontsize=7.5, ncol=2, loc="upper right")
    y_max = max(row["material_rel_l2"] for row in scores)
    axes[1, 0].set_ylim(0, 1.10 * y_max)
    axes[1, 1].set_ylim(0, 1.10 * y_max)
    fig.savefig(ROOT / "figures/model_gap_conditions.png", dpi=175)
    plt.close(fig)


def draw_material(gaps, scores):
    # The pair with the largest fixed-truth model mismatch illustrates the stress case.
    selected = max(gaps, key=lambda x: x["total_gap_over_delta"])["object_id"]
    fig, axes = plt.subplots(2, 6, figsize=(12.0, 4.3), layout="constrained")
    for i, (layout, pretty) in enumerate(LAYOUTS.items()):
        case = f"{selected}__{layout.lower()}"
        with np.load(RESULTS / f"offline/{case}.npz", allow_pickle=False) as model:
            reference = model["material"]
        image = axes[i, 0].imshow(reference.reshape(32, 32), origin="lower", vmin=0, vmax=4.0, cmap="viridis")
        axes[i, 0].set(title=f"{selected} reference", ylabel=pretty)
        for j, (key, title) in enumerate(CONDITIONS.items(), 1):
            selection = [row for row in scores if row["observation"] == case and
                         ((j == 5 and row["method"] == key) or
                          (j < 5 and row["method"] == "physical_LBFGSB_12" and row["condition"] == key))]
            file_suffix = key.lower() if j == 5 else key
            with np.load(RESULTS / f"predictions/{case}__{file_suffix}.npz", allow_pickle=False) as payload:
                material = payload["material"]
            axes[i, j].imshow(material.reshape(32, 32), origin="lower", vmin=0, vmax=4.0, cmap="viridis")
            axes[i, j].set_title(f"{title}\nL2 {selection[0]['material_rel_l2']:.3f}", fontsize=9)
    for axis in axes.flat:
        axis.set(xticks=[], yticks=[])
    fig.colorbar(image, ax=axes, shrink=0.78, pad=0.01, label="Contrast χ (0–4)")
    fig.savefig(ROOT / "figures/model_gap_materials.png", dpi=175)
    plt.close(fig)


if __name__ == "__main__":
    gaps = json.loads((RESULTS / "gaps.json").read_text())["rows"]
    scores = json.loads((RESULTS / "scores.json").read_text())
    draw_comparison(gaps, scores)
    draw_material(gaps, scores)
