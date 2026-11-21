"""Plot the archived mirror certificates and the prospective stopping control."""

from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE/"results"
LABELS = {
    "srap_PW360_10x20": "Full ring, plane waves",
    "srap_PW180_10x20": "Half ring",
    "srap_PW180_FLIP_RX_10x20": "Half ring, flipped Rx",
    "srap_PW180_FLIP_TX_10x20": "Half ring, flipped Tx",
    "srap_PW120_7x14": "120° sparse view",
    "srap_LS_ARC_10x20": "Arc, line sources",
    "srap_LS_LINE_SAME_10x20": "Line, same side",
    "srap_LS_LINE_OPPOSITE_10x20": "Line, opposite sides",
}


def load(name):
    return json.loads((RESULTS/name).read_text())


def main():
    data = load("field_results.json")["rows"]
    summary = load("exposure_prefix.json")["summary"]["prospective_radius"]["all-applicable"]
    plt.rcParams.update({"text.color": "#203246", "axes.labelcolor": "#203246",
                         "font.family": "DejaVu Sans", "font.size": 10})
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12, 5.4), dpi=170,
                                gridspec_kw={"width_ratios": [1.25, 1.]})
    fig.subplots_adjust(left=.181, right=.97, top=.79, bottom=.19, wspace=.52)
    fig.text(.06, .94, "How much of the archive is needed to rule out a mirror pair?",
             fontsize=19, fontweight="bold")
    fig.text(.06, .875, "96² full-wave data → 32² model, 33 independent SRAP materials",
             fontsize=11, color="#4a637a")
    labels = []
    for n, (layout, label) in enumerate(LABELS.items()):
        results = [r["results"][r["largest_verified_mirror"]]["all_data"]["min_EO_over_delta"]
                   for r in data if r["layout"] == layout and r["largest_verified_mirror"]]
        phase = .13*np.sin(.77*np.arange(len(results)))
        ax.scatter(results, n+phase, s=18, alpha=.65, color="#326b90", linewidths=.1)
        labels.append(f"{label} ({len(results)})")
    ax.axvline(1, linestyle=(0,(4,3)), color="#a75337", lw=1.2)
    ax.set(yticks=np.arange(len(labels)), yticklabels=labels, xlim=(.7, 7.05),
           ylim=(len(labels)-.55,-.75), xlabel="min(E, O) / archived joint noise radius")
    ax.set_title("Full-layout mirror statistic", loc="left", fontweight="bold")
    ax.tick_params(axis="y", length=0)
    ax.text(.0, -.225, "162 certified views · 42 views lack compatible port geometry",
            transform=ax.transAxes, color="#4a637a", fontsize=9, ha="left")

    labels = ["Full batch", "Fixed half", "Hindsight fixed", "Adaptive"]
    exposure = [summary["raw_full_exposures"], summary["fixed_half_exposures"],
                summary["hindsight_total_exposures_including_failures"],
                summary["adaptive_exposures_including_failures"]]
    passed = [summary["full_certified"], summary["half_certified"],
              summary["adaptive_certified"], summary["adaptive_certified"]]
    bars = bx.barh(labels, exposure, color=["#91a5b5","#91a5b5","#b18560","#326b90"], height=.56)
    bx.set(xlim=(0,2050), ylim=(3.65,-.65), xlabel="Source exposures across 152 compliant views")
    bx.set_title("Fixed prospective noise cap", loc="left", fontweight="bold")
    bx.tick_params(axis="y", length=0)
    for bar, count, ok in zip(bars, exposure, passed, strict=True):
        bx.text(count+34, bar.get_y()+.5*bar.get_height(), f"{count}  ·  {ok} certified",
                va="center", fontsize=9, color="#203246")
    bx.text(.0, -.225, "162 attempted · 10 exceeded the cap after replay",
            transform=bx.transAxes, color="#4a637a", fontsize=9, ha="left")
    for axis in (ax,bx):
        for side in ("top","right","left"):
            axis.spines[side].set_visible(False)
        axis.spines["bottom"].set_color("#b9c8d4")
    fig.savefig(HERE/"archive_reflection.png", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
