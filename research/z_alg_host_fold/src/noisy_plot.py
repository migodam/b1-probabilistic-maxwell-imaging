"""Render the matched-work frontiers and the final physical coordinate hull."""

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from noisy_cover import RESULTS, ROOT


def main():
    results = json.loads((RESULTS / "noisy_cover.json").read_text())
    colors = {"direct": "#8c95a3", "direct_centered": "#b44438", "chart": "#147f82"}
    styles = {"direct": "--", "direct_centered": "-", "chart": "-"}
    labels = {"direct": "Direct joint score", "direct_centered": "Centered direct", "chart": "Receiver chart"}
    plt.rcParams.update({
        "font.size": 9,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
    })
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.4), gridspec_kw={"width_ratios": [1.1, 1.1, 1.]})
    for ax, (scope, title) in zip(axes[:2], (
        ("physical", "(a) Physical material box"),
        ("local", "(b) Local crop"),
    ), strict=True):
        arms = results["scopes"][scope]["arms"]
        for arm, result in arms.items():
            time = [item["wall_seconds_including_setup"] for item in result["frontier"]]
            volume = [item["outer_fraction"] for item in result["frontier"]]
            ax.loglog(time, volume, color=colors[arm], linestyle=styles[arm], lw=1.9, label=labels[arm])
            ax.scatter([time[-1]], [volume[-1]], color=colors[arm], s=23, zorder=4)
        ax.grid(which="major", alpha=.22, linewidth=.55)
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xlabel("Wall time incl. setup (s)")
        ax.set_ylim(3e-13, 1.45)
        ax.set_xlim(.026, .432)
        ax.set_xticks([.03, .05, .1, .2, .4], ["0.03", "0.05", "0.1", "0.2", "0.4"])
        ax.xaxis.set_minor_formatter(NullFormatter())
    axes[0].set_ylabel("Retained volume / input volume")
    axes[1].legend(loc="upper right", frameon=False, fontsize=8)
    ax = axes[2]
    ax.set_title("(c) Final physical-box hull", loc="left", fontweight="bold")
    for position, arm in ((1, "direct_centered"), (0, "chart")):
        low, high = results["scopes"]["physical"]["arms"][arm]["final_hull_outward"][5]
        ax.hlines(position, low, high, lw=5, color=colors[arm])
        ax.scatter([low, high], [position, position], s=23, color=colors[arm], zorder=3)
    six_values = [v["material"][5] for k, v in results["shared_point_witnesses"].items() if k.startswith("point_fiber_minus")]
    ax.scatter(six_values, [-.7, -.7], marker="D", color="#28282c", s=24, zorder=3)
    ax.set_yticks([1, 0, -.7], ["Centered direct", "Receiver chart", "Minus roots"])
    ax.set_xlim(24.1, 30.25)
    ax.set_ylim(-1.13, 1.55)
    ax.grid(axis="x", alpha=.22, linewidth=.55)
    ax.set_xlabel("Susceptibility in cell 6")
    fig.tight_layout(w_pad=2.1)
    out = ROOT / "figures" / "noisy_cover_frontiers"
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", pad_inches=.08)
    fig.savefig(out.with_suffix(".png"), dpi=180, bbox_inches="tight", pad_inches=.08)


if __name__ == "__main__":
    main()
