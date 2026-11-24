"""Plot the original-noise cardinality tests and the source-count ceilings."""

from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PAIRS = (
    ("r10_twod_0001", "archived_cn99", "R10 ring (24 × 24, threshold 4)", "#006a85"),
    ("optfresh_01__pw360_10x20__n02", "archived_2pct_with_generation_padding",
     "optfresh 01 (96 × 96, threshold 0.5)", "#97393b"),
)


def main():
    report = json.loads((HERE / "area_receipt.json").read_text())["cases"]
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.dpi": 160, "pdf.fonttype": 42})
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.6, 4.3), layout="constrained")
    for case, budget, label, color in PAIRS:
        data = report[case]
        result = data["budgets"][budget]
        part = data["tails"][:6]
        q = [row["maximum_high_pixels_hypothesized"] for row in part]
        ratio = [row["residual_lower"] / result["radius_upper"] for row in part]
        left.plot(q, ratio, color=color, label=label, marker="o", markersize=5)
        last_q = result["proved_high_pixels"] - 1
        left.scatter(q[last_q], ratio[last_q], facecolor=color, edgecolor="white", s=90, zorder=5)
    left.axhline(1, color="#4a5360", linewidth=1, linestyle="--")
    left.text(5, 1.08, "Original noise radius", ha="right", va="bottom", color="#4a5360")
    left.set_xticks(np.arange(6))
    left.set_xlabel("At most q high pixels")
    left.set_ylabel("Certified residual lower bound / noise radius")
    left.set_title("A. Excluding small high-contrast sets", loc="left", fontsize=11, fontweight="bold")
    left.set_yscale("symlog", linthresh=1)
    left.set_yticks([0, .5, 1, 2, 5, 10], labels=["0", "0.5", "1", "2", "5", "10"])
    left.set_ylim(-.12, 20)
    left.legend(loc="upper right", fontsize=8.5, frameon=False)
    left.grid(axis="y", color="#dce2e8", linewidth=.7, alpha=.85)
    categories = ("Certified minimum", "0.5% area target", "Source-count ceiling")
    keys = ("proved_high_pixels", "half_percent_area_target_pixels", "strongest_direct_rank_ceiling")
    for i, (case, budget, _, color) in enumerate(PAIRS):
        row = report[case]["budgets"][budget]
        offset = (-.17 if i == 0 else .17)
        positions = np.arange(3) + offset
        counts = [row[key] for key in keys]
        right.barh(positions, counts, height=.29, color=color, alpha=.9)
        for y, value in zip(positions, counts, strict=True):
            right.annotate(str(value), (value, y), xytext=(4, 0), textcoords="offset points",
                           va="center", fontsize=9, color=color)
    right.set_xlim(0, 54)
    right.set_yticks(np.arange(3), categories)
    right.invert_yaxis()
    right.set_xlabel("High pixels")
    right.set_title("B. Area target on each native grid", loc="left", fontsize=11, fontweight="bold")
    right.grid(axis="x", color="#dce2e8", linewidth=.7, alpha=.85)
    right.set_axisbelow(True)
    fig.savefig(HERE / "area_certificate.png", dpi=220, facecolor="white")
    fig.savefig(HERE / "area_certificate.pdf", facecolor="white")


if __name__ == "__main__":
    main()
