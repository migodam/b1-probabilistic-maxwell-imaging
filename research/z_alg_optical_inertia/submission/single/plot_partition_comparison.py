"""Plot the matched certificate benchmark from the frozen outward proofs."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parents[1] / "results/theory"
CASES = (("01", "Slotted block"), ("05", "Staggered grooves"))
ROWS = ("4 regions, single coherent SDP", "Joint, uniform weights",
        "16 regions, single coherent", "16 regions, joint matched")
COLORS = ("#273744", "#22548a", "#a34c11", "#086959")
MARKERS = ("o", "D", "^", "s")


def read(name, root=HERE):
    return json.loads((root/name).read_text())


def main():
    old = read("fine_verified.json", PREVIOUS)["cases"]
    local = read("local_moment.json")
    coherent = read("coherent_general_verified.json")["cases"]
    matched = read("matched_joint_replay.json")["cases"]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.45), sharex=True, sharey=True)
    fig.subplots_adjust(left=.337, right=.980, top=.86, bottom=.278, wspace=.21)
    positions = np.array([3, 2, 1, 0])
    for index, (ax, (number, specimen)) in enumerate(zip(axes, CASES, strict=True)):
        case = f"optfresh_{number}__pw360_10x20__n02"
        radius = coherent[case]["radius"]
        assert radius == local["padded_noise_radius"][case]
        uniform = next(row for row in old if row["id"] == case and row["cap"] == .75
                       and row["profile"][3] == .3)
        points = np.array([
            np.sqrt(local["observations"][case]["spatial_SDP_optimal_squared_residual"]),
            uniform["distance_lower"],
            coherent[case]["coherent_lower"],
            matched[case]["joint_trace_distance_lower"],
        ]) / radius
        ax.axvspan(1, 2.42, color="#e8f2ee", alpha=.65, linewidth=0, zorder=0)
        ax.axvline(1, color="#53595d", lw=1.1, ls=(0, (3, 2)), zorder=1)
        for height, value, color, marker in zip(positions, points, COLORS, MARKERS, strict=True):
            ax.hlines(height, 0, value, lw=1.7, color=color, alpha=.81)
            ax.plot(value, height, marker=marker, ms=5.2, color=color, clip_on=False)
            ax.annotate("0" if value == 0 else f"{value:.2f}", (value, height),
                        xytext=(5, 0), textcoords="offset points", va="center",
                        color=color, fontsize=8.3)
        ax.set_xlim(-.04, 2.42)
        ax.set_ylim(-.5, 3.5)
        ax.set_xticks([0, .5, 1, 1.5, 2])
        ax.set_yticks(positions, ROWS)
        ax.set_title(f"({chr(ord('a')+index)}) {number}  {specimen}", pad=5.5, loc="left",
                     fontsize=9.3, fontweight="bold")
        ax.tick_params(axis="y", length=0, pad=7, labelleft=index == 0)
        ax.tick_params(axis="x", length=3, color="#53595d", pad=3)
        for side in ("left", "right", "top"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color("#7b8388")
    fig.text(.66, .137, "Certified residual / error radius", ha="center", fontsize=8.7)
    threshold = Line2D([], [], color="#53595d", lw=1.1, ls=(0, (3, 2)))
    fig.legend([threshold], ["Exclusion threshold for cap 0.75"],
               loc="lower center", bbox_to_anchor=(.66, .0), frameon=False,
               borderaxespad=0.2, handlelength=1.8, fontsize=8.2)
    fig.savefig(HERE/"partition_comparison.pdf")
    fig.savefig(HERE/"partition_comparison.png", dpi=240)


if __name__ == "__main__":
    main()
