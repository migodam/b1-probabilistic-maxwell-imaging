"""Plot directly certified receiver inertias from the two native controls."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / "figures"
POSITIVE, NEGATIVE = "#17685f", "#a34329"
INK, SECONDARY, RULE = "#172733", "#475c67", "#d9e4e8"


def read(name):
    return json.loads((HERE / name).read_text())


def half_width(receipt):
    return float(receipt["box"].rsplit("_", 1)[1])


def receiver_rows():
    narrow = read("coherent_control.json")
    wide = read("occupied_wide_control.json")
    if narrow["status"] != "CERTIFIED" or wide["status"] != "COMPLETE":
        raise ValueError("Both native controls must have complete proof receipts.")
    receivers = narrow["n_receivers"]
    arms = {arm["arm"]: arm for arm in narrow["arms"]}
    target = narrow["common_target_fluence_interval"]
    arms_in_order = [arms["coherent"], arms["plane_wave"]]
    beams = arms_in_order[0]["number_of_beams"]
    if not all(a["number_of_beams"] == beams and a["target_exact_equal_by_normalization"]
               and a["direct_full_grid_fluence_interval"] == target for a in arms_in_order):
        raise ValueError("Narrow arms need verified matching fluence and exposure counts.")
    narrow_rows = [(label, arm["source_deficit_at_this_acquisition"],
                    arm["receiving_block"]["inertia"])
                   for label, arm in zip(("Coherent eigenbeams", "Bit-reversal plane waves"),
                                         arms_in_order, strict=True)]
    wide_rows = [(f"{r['s']} plane waves", r["h"], r["receive"]["inertia"])
                 for r in sorted(wide["rows"], key=lambda row: row["s"])
                 if r["saturated"] and r["k"] == wide["state_negative"]]
    if len(wide_rows) != len(wide["rows"]):
        raise ValueError("Wide panel requires the verified saturated prefixes.")
    for _, deficit, inertia in narrow_rows + wide_rows:
        if inertia is None or sum(inertia) != receivers or receivers-inertia[0] != deficit:
            raise ValueError("Receiver count needs a direct inertia proof matching the source deficit.")
    return narrow, wide, narrow_rows, wide_rows, receivers, beams, sum(target)/2


def bars(ax, rows, total):
    ax.set_axis_off()
    ax.set(xlim=(-.2, total+.15), ylim=(-.48, 1.48))
    for y, (label, deficit, (pos, neg)) in zip((1, 0), rows, strict=True):
        ax.barh(y, pos, height=.57, color=POSITIVE)
        if neg:
            ax.barh(y, neg, left=pos, height=.57, color=NEGATIVE)
            ax.text(pos+neg/2, y, str(neg), color="white", ha="center", va="center",
                    fontsize=9.4, fontweight="bold")
        ax.text(pos/2, y, str(pos), color="white", ha="center", va="center",
                fontsize=9.4, fontweight="bold")
        ax.text(-1.35, y+.08, label, ha="right", va="center", fontsize=9)
        ax.text(-1.35, y-.20, rf"Source deficit $h={deficit}$",
                ha="right", va="center", fontsize=7.3, color=SECONDARY)


def main():
    narrow, wide, upper, lower, receivers, beams, fluence = receiver_rows()
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "text.color": INK,
        "pdf.fonttype": 42,
    })
    fig = plt.figure(figsize=(7.25, 3.30), facecolor="white")
    for y, bounds, title, subtitle, rows in (
        (.92, (.31, .61, .66, .235),
         f"A   Narrow box, half-width {half_width(narrow):g}",
         f"{beams} illuminations per arm  ·  Total fluence {fluence:.6f}", upper),
        (.50, (.31, .19, .66, .235),
         f"B   Wide box, half-width {half_width(wide):g}",
         rf"Unit-amplitude plane-wave prefixes  ·  Both saturate at $k=\nu={wide['state_negative']}$",
         lower),
    ):
        fig.text(.043, y, title, fontsize=10.3, fontweight="bold", va="center")
        fig.text(.043, y-.058, subtitle, fontsize=8.2, color=SECONDARY, va="center")
        bars(fig.add_axes(bounds), rows, receivers)
    fig.add_artist(plt.Line2D([.043, .968], [.576, .576], transform=fig.transFigure,
                              color=RULE, linewidth=.78))
    fig.legend(
        handles=[Patch(color=POSITIVE, label="Positive receiving directions"),
                 Patch(color=NEGATIVE, label="Negative receiving directions")],
        loc="center left", bbox_to_anchor=(.036, .115), ncol=2, frameon=False,
        handlelength=1.2, handleheight=.8, columnspacing=2.05, fontsize=8.4,
    )
    fig.text(.044, .047, rf"At saturation, $h$ counts missing positive directions out of {receivers}.",
             color=SECONDARY, va="center", fontsize=8.1)
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "native-illumination.pdf", metadata={"Title": "Certified native receiver index"})
    fig.savefig(OUT / "native-illumination.png", dpi=230)
    plt.close(fig)


if __name__ == "__main__":
    main()
