"""Plot material arrays and the certified pair separation in the added ports."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from matplotlib.ticker import MultipleLocator
import numpy as np

HERE = Path(__file__).resolve().parent
RESULT = HERE / "results"
ARRAYS = HERE.parent / "array" / "results" / "physical_arrays.npz"
PAIRS = ("disks12_high__1__double", "disks16_high__0.5__double")
TITLES = ("12 disks · high · full amplitude", "16 disks · high · half amplitude")
COLORS = ("#485b70", "#d5933e", "#6c85bb")


def disks_panel(ax, chi, coordinates, task, sign):
    active = len(chi)
    host = PatchCollection([Circle(p, 0.04) for p in coordinates[active:]],
                           facecolor="#d9e1e5", edgecolor="#b5c0c6", linewidth=0.45)
    disks = PatchCollection([Circle(p, 0.04) for p in coordinates[:active]],
                            cmap="viridis", edgecolor="#243d43", linewidth=0.65)
    disks.set_array(chi)
    disks.set_clim(0.5, 30)
    ax.add_collection(host)
    ax.axvline(-0.13, color="#7c5555", linestyle=(0, (3, 3)), linewidth=1, alpha=0.8)
    ax.text(-0.15, 0.45, "ROI", ha="right", va="center", fontsize=8, color="#7c5555")
    ax.add_collection(disks)
    ax.set(xlim=(-0.68, 0.68), ylim=(-0.51, 0.51), aspect="equal")
    ax.set_title(f"{sign} endpoint  ·  ROI mean {task @ chi:.3f}", fontsize=11)
    ax.xaxis.set_major_locator(MultipleLocator(0.4))
    ax.yaxis.set_major_locator(MultipleLocator(0.4))
    ax.tick_params(labelsize=8, length=3, color="#83939b")
    for side in ax.spines:
        ax.spines[side].set_color("#c8d3d6")
    return disks


def channel_panel(ax, pair, audits):
    values = [0]
    for arm in ("angle", "frequency"):
        interval = audits[(pair, arm)]
        values.append((interval["lo"] + interval["hi"]) / 2)
    ax.axvline(1, color="#717a7d", linestyle=(0, (4, 4)), linewidth=1.0)
    ax.barh(range(3), values, color=COLORS, height=0.48)
    ax.scatter([0], [0], color=COLORS[0], s=16, zorder=3)
    for row, value in enumerate(values):
        ax.text(value + 0.045, row, f"{value:.3f}" if value else "exact zero", va="center",
                fontsize=9, color="#213147")
    ax.set_yticks(range(3), ["Base", "+ angle", "+ frequency"])
    ax.tick_params(axis="y", length=0, labelsize=9)
    ax.tick_params(axis="x", length=3, labelsize=8)
    ax.set_xlim(-0.025, 1.70)
    ax.set_ylim(2.56, -0.56)
    ax.set_xlabel("Endpoint distance / 2δ", fontsize=9, labelpad=5)
    ax.set_title("Port separation at 10% noise", fontsize=11)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#c8d3d6")


def main():
    report = json.loads((RESULT / "certified_pair_augmented_ports.json").read_text())
    audits = {(item["baseline_pair"], item["arm"]): item["chord_over_2delta_enclosure"]
              for item in report["results"]}
    plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": "#213147",
                         "axes.labelcolor": "#213147", "xtick.color": "#516371",
                         "ytick.color": "#516371", "axes.titleweight": "normal"})
    figure, axes = plt.subplots(2, 3, figsize=(11.2, 6.6), dpi=160,
                               gridspec_kw={"width_ratios": [1, 1, 1.05]})
    figure.subplots_adjust(left=0.075, bottom=0.13, right=0.97, top=0.83,
                           wspace=0.47, hspace=0.72)
    figure.text(0.075, 0.974, "Different materials, identical base port readings",
                ha="left", va="top", fontsize=19, fontweight="bold")
    with np.load(ARRAYS) as source, np.load(RESULT / "finite_chord_arrays.npz") as materials:
        for row, (pair, title) in enumerate(zip(PAIRS, TITLES, strict=True)):
            case, _, _ = pair.split("__")
            coords = source[f"{case.split('_')[0]}__centers"]
            task = materials[f"{case}__task"]
            for col, sign in enumerate(("minus", "plus")):
                chi = materials[f"{pair}__{sign}"]
                color_key = disks_panel(axes[row, col], chi, coords, task,
                                        "−" if sign == "minus" else "+")
            channel_panel(axes[row, 2], pair, audits)
            label_height = axes[row, 0].get_position().y1 + 0.075
            figure.text(0.075, label_height, title, fontsize=11, fontweight="bold")
    bar_ax = figure.add_axes([0.075, 0.052, 0.48, 0.017])
    bar = figure.colorbar(color_key, cax=bar_ax, orientation="horizontal",
                          ticks=[0.5, 5, 10, 15, 20, 25, 30])
    bar.outline.set_visible(False)
    bar.ax.tick_params(labelsize=8, length=2)
    figure.text(0.59, 0.059, "χ = εᵣ − 1    ·    gray disks: known host",
                va="center", fontsize=9, color="#516371")
    figure.savefig(HERE / "material_pairs.png", facecolor="white")
    plt.close(figure)


if __name__ == "__main__":
    main()
