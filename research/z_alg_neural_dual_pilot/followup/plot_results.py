"""Plot the frozen full-ring segment, receiver masks, and same-observation diagnostic."""

from argparse import ArgumentParser
import json
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR
from math import cos, radians, sin
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parent
LEVELS = (3, 5, 9)
COHORTS = {
    "fresh": "Independent follow-up",
    "external8": "External shapes",
    "historical": "Historical stress set",
    "validation": "Development validation",
}
VERDICTS = {
    "PASS": ("Strict pass", "#2265a1"),
    "BELOW": ("Whole segment below", "#227969"),
    "UNRESOLVED": ("Unresolved", "#a75e14"),
}
INK = "#253039"
GRID = "#dce4e9"
ROOF = "#2a8080"
RADIUS = "#a75e14"
MASKS = {
    "full": "Full ring",
    "alternating_even": "Alternating even",
    "alternating_odd": "Alternating odd",
    "arc_0_deg": "Arc from 0°",
    "arc_180_deg": "Arc from 180°",
}


def read_cohorts():
    return {
        name: json.loads((ROOT / "results" / f"{name}_segment_cover.json").read_text())
        for name in COHORTS
    }


def percent_panel(ax, data, title):
    n_views = data["conditions"]
    n_materials = data["independent_objects"]
    view_type = "paired views" if n_views > n_materials else "views"
    ax.set_title(f"{title}\n{n_materials} materials · {n_views} {view_type}", loc="left", pad=10)
    levels = data["summary"]["all_views"]
    for k, (verdict, (_, color)) in enumerate(VERDICTS.items()):
        positions = [j + (k - 1) * 0.23 for j in range(len(LEVELS))]
        counts = [levels[str(level)]["verdicts"][verdict] for level in LEVELS]
        shares = [count / n_views * 100 for count in counts]
        ax.bar(positions, shares, width=0.20, color=color, zorder=3)
        for x, share, count in zip(positions, shares, counts, strict=True):
            ax.annotate(str(count), (x, share), xytext=(0, 3), textcoords="offset points",
                        ha="center", va="bottom", color=color, fontsize=8)
    ax.set_xticks(range(len(LEVELS)), LEVELS)
    ax.set_xlabel("Weight anchors")
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylim(0, 112)
    ax.grid(axis="y", color=GRID, linewidth=0.7)


def profile_panel(ax, row, title):
    nodes = row["nodes"]
    level = row["levels"]["9"]
    pieces = level["pieces"]
    edges = [pieces[0][0], *(right for _, right, _ in pieces)]
    heights = [upper for _, _, upper in pieces]
    lower = level["sampled_table"]["lower"]
    upper = level["continuous_upper"]
    radius = row["radius"]
    precision = Decimal("0.00001")
    lower_label = Decimal.from_float(lower).quantize(precision, rounding=ROUND_FLOOR)
    upper_label = Decimal.from_float(upper).quantize(precision, rounding=ROUND_CEILING)
    ax.set_title(title, loc="left", pad=10)
    ax.stairs(heights, edges, baseline=None, color=ROOF, linewidth=1.35, zorder=3)
    ax.axhline(radius, color=RADIUS, linestyle="--", linewidth=1.15, zorder=2)
    ax.vlines([node["t"] for node in nodes], [node["distance_lower"] for node in nodes],
              [node["distance_upper"] for node in nodes], color=INK, linewidth=1.0, zorder=4)
    ax.scatter([node["t"] for node in nodes],
               [(node["distance_lower"] + node["distance_upper"]) / 2 for node in nodes],
               s=20, color=INK, linewidths=0.7, edgecolors="white", zorder=5)
    ax.set_xlim(-0.025, 1.025)
    ax.set_ylim(-0.003, max(max(heights), radius) * 1.10)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xlabel("Weight coordinate $t$")
    ax.text(0, -0.28, rf"$d_{{\mathrm{{max}}}} \in [{lower_label},\,{upper_label}]$"
            rf"     $r \approx {radius:.5f}$", transform=ax.transAxes, va="top")
    ax.grid(axis="y", color=GRID, linewidth=0.7)


def find_view(data, object_id, noise_pct):
    return next(row for row in data["rows"]
                if row["object_id"] == object_id and row["noise_pct"] == noise_pct)


def save(fig, name):
    out = ROOT / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / f"{name}.png", dpi=220)
    fig.savefig(out / f"{name}.pdf")
    plt.close(fig)


def plot_segment(data):
    fig = plt.figure(figsize=(12, 7.5))
    grid = fig.add_gridspec(2, 4, left=0.075, right=0.98, bottom=0.14, top=0.80,
                            wspace=0.34, hspace=0.69, height_ratios=(1, 1.08))
    fig.suptitle("Strict coverage of one positive-definite 16-block weight segment",
                 x=0.075, y=0.97, ha="left", fontsize=16, weight="bold")
    fig.legend(
        handles=[Patch(facecolor=color, label=label) for label, color in VERDICTS.values()],
        loc="upper left", bbox_to_anchor=(0.069, 0.919), ncols=3, frameon=False, fontsize=9,
    )
    count_axes = [fig.add_subplot(grid[0, j]) for j in range(len(COHORTS))]
    for ax, (name, title) in zip(count_axes, COHORTS.items(), strict=True):
        percent_panel(ax, data[name], title)
    count_axes[0].set_ylabel("Share of views (%)")

    profiles = [fig.add_subplot(grid[1, :2]), fig.add_subplot(grid[1, 2:])]
    profile_panel(profiles[0], find_view(data["fresh"], "ip_port_09", 1),
                  "ip_port_09 · 1% noise · strict pass")
    profile_panel(profiles[1], find_view(data["historical"], "s2d_004", 2),
                  "s2d_004 · 2% noise · unresolved")
    profiles[0].set_ylabel(r"Certificate distance $d(t)$")
    fig.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="none", color=INK, label="Sampled anchor"),
            Line2D([0], [0], color=ROOF, label="Certified subinterval upper bound"),
            Line2D([0], [0], color=RADIUS, linestyle="--", label=r"Noise radius $r$"),
        ],
        loc="lower left", bbox_to_anchor=(0.067, 0.010), ncols=3, frameon=False, fontsize=9,
    )
    save(fig, "conditional_weight_family")


def receiver_panel(ax, coords, mask_geometry, title, retained_modes, universal_zero):
    kept = mask_geometry["0.0"]["kept"]
    ax.add_patch(plt.Circle((0, 0), 1, facecolor="none", edgecolor=GRID, linewidth=1.0))
    ax.scatter(*zip(*coords, strict=True), s=22, facecolors="white", edgecolors="#b5bfc7",
               linewidths=0.9, zorder=2)
    ax.scatter(*zip(*(coords[i] for i in kept), strict=True), s=29, color=VERDICTS["PASS"][1],
               linewidths=0.7, edgecolors="white", zorder=3)
    ax.set_title(f"{title}\n{len(kept)} receivers", pad=9)
    ax.text(0.5, -0.025, str(retained_modes),
            ha="center", va="top", transform=ax.transAxes, fontsize=9)
    if universal_zero:
        ax.text(0.5, -0.145, "Zero for every observation", ha="center", va="top",
                transform=ax.transAxes, fontsize=8, color=ROOF)
    ax.set(xlim=(-1.19, 1.19), ylim=(-1.19, 1.19), aspect="equal")
    ax.set_axis_off()


def plot_aperture(data):
    geometry = json.loads((ROOT / "results" / "masked_geometry.json").read_text())
    continuous = json.loads((ROOT / "results" / "masked_continuous_inertia.json").read_text())
    masked = {
        name: json.loads((ROOT / "results" / f"{name}_masked_ports.json").read_text())
        for name in COHORTS
    }
    angles = geometry["receiver_order"]["receiver_degrees_in_storage_order"]
    coords = [(cos(radians(degrees)), sin(radians(degrees))) for degrees in angles]
    full_example = masked["fresh"]["rows"][0]["masks"]["full"]["points"][0]
    source_count = full_example["retained_source_count"]
    fig = plt.figure(figsize=(11, 6.75))
    grid = fig.add_gridspec(2, 6, left=0.030, right=0.985, bottom=0.15, top=0.79,
                            width_ratios=(1.46, 1, 1, 1, 1, 1),
                            height_ratios=(1.22, 1), hspace=0.32, wspace=0.10)
    fig.suptitle("Passive receiver masks on the same frozen observations",
                 x=0.030, y=0.96, ha="left", fontsize=16, weight="bold")
    fig.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="none", color=VERDICTS["PASS"][1],
                   label="Retained receiver"),
            Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="white",
                   markeredgecolor="#9ba9b5", label="Deleted receiver"),
        ],
        loc="upper left", bbox_to_anchor=(0.024, 0.906), ncols=2, frameon=False, fontsize=9,
    )
    ring_axes = []
    for index, (key, name) in enumerate(MASKS.items(), start=1):
        proof = continuous["profiles"][key] if key != "full" else {
            "negative_count_every_t": 0, "universal_zero_every_Y_every_t": False,
        }
        ax = fig.add_subplot(grid[0, index])
        receiver_panel(ax, coords,
                       geometry["per_mask_geometry"][key], name,
                       source_count - proof["negative_count_every_t"],
                       proof["universal_zero_every_Y_every_t"])
        ring_axes.append(ax)
    ring = ring_axes[0].get_position()
    fig.text(0.030, ring.y0 - 0.025 * ring.height, "Source modes retained\n" r"for all $t \in [0,1]$",
             va="top", fontsize=9)

    label_axis = fig.add_subplot(grid[1, 0])
    label_axis.set_axis_off()
    label_axis.text(0, 1.12, "Exact-zero certificates at the nine sampled weights",
                    va="center", fontsize=10, weight="bold", transform=label_axis.transAxes)
    for index, key in enumerate(MASKS, start=1):
        table_axis = fig.add_subplot(grid[1, index])
        table_axis.set_axis_off()
        for row_index, (cohort, cohort_name) in enumerate(COHORTS.items()):
            y = 0.87 - row_index * 0.245
            summary = masked[cohort]["summary"]["all"]
            denominator = summary["views"] * LEVELS[-1]
            zero_count = summary["masks"][key]["point_exact_zeros"]
            if index == 1:
                n_materials = data[cohort]["independent_objects"]
                label_axis.text(0, y, f"{cohort_name}\n{n_materials} materials · {summary['views']} views",
                                va="center", fontsize=9, transform=label_axis.transAxes)
            table_axis.text(0.5, y, f"{zero_count} / {denominator}", ha="center", va="center",
                            fontsize=10, transform=table_axis.transAxes)
            table_axis.plot([0.04, 0.96], [y - 0.12, y - 0.12], color=GRID, linewidth=0.75,
                            transform=table_axis.transAxes)
    full_passes = sum(masked[name]["summary"]["all"]["masks"]["full"]["view_sampled_passes"]
                      for name in COHORTS)
    restricted_passes = {
        sum(masked[name]["summary"]["all"]["masks"][key]["view_sampled_passes"] for name in COHORTS)
        for key in MASKS if key != "full"
    }
    restricted_range = (f"{min(restricted_passes)}–{max(restricted_passes)}"
                        if len(restricted_passes) > 1 else str(next(iter(restricted_passes))))
    fig.text(0.030, 0.081,
             f"Strictly passing views at any sampled weight   Full ring {full_passes}"
             f"   ·   Per ten-receiver mask {restricted_range}", fontsize=9)
    fig.text(0.030, 0.052,
             f"Each entry counts view–weight pairs. All configurations keep the {source_count} illuminations"
             " and the original absolute noise radius.", fontsize=9)
    save(fig, "passive_receiver_masks")


def plot_diagnostic():
    results = ROOT / "results"
    focus = json.loads((results / "galerkin_historical_focus.json").read_text())
    donor = json.loads((results / "transferred_weight_donor.json").read_text())
    row = find_view(focus, "s2d_004", 2)
    direct = find_view(donor, "s2d_004", 2)
    radius = row["radius"]
    identity = donor["original_individual_identity_check"][0]
    if direct["radius"] != radius or not identity["overlap"]:
        raise ValueError("The direct native and line-focus receipts must refer to the same observation.")
    sampled = row["sampled_native_controls"]["9"]
    intervals = [
        (4, "Native outer chords · 9 nodes", sampled["lower"], sampled["upper"]),
        *((y, f"Galerkin · {pieces} subintervals", row["levels"][str(pieces)]["lower"],
           row["levels"][str(pieces)]["upper"])
          for y, pieces in zip((3, 2, 1), (4, 16, 64), strict=True)),
        (-0.6, "Native 64-block gradient weight", *identity["original"]),
    ]
    transfer = json.loads((results / "transferred_weight_transfer16.json").read_text())
    arms = transfer["summary"]["all_views"]["arms"]
    one = arms["added_weight_alone"]
    controls = arms["full_archived"]
    percentages = sorted({noise for _, noise in controls["passing_views"]})
    noise_label = ", ".join(f"{noise:g}%" for noise in percentages)

    fig, ax = plt.subplots(figsize=(10.4, 4.5))
    fig.subplots_adjust(left=0.254, right=0.738, bottom=0.25, top=0.76)
    fig.text(0.04, 0.94, "A single weight certifies exclusion on the same observation",
             fontsize=15, weight="bold", va="top")
    fig.text(0.04, 0.858, "s2d_004 · 2% noise · 20 receivers × 10 illuminations"
             f" · Material peak cap {donor['cap']:g}", fontsize=9.5)
    ax.set_xlim(0.065, 0.172)
    ax.set_ylim(-1.14, 5.26)
    ax.set_xticks([0.07, 0.09, 0.11, 0.13, 0.15, 0.17])
    ax.set_yticks([y for y, *_ in intervals], [label for _, label, *_ in intervals], fontsize=9)
    ax.tick_params(axis="y", length=0, pad=8)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", color=GRID, linewidth=0.65)
    ax.axvline(radius, color=RADIUS, linestyle="--", linewidth=1.15, zorder=2)
    ax.text(radius + 0.0014, 4.87, f"Noise radius {radius:.8f}", va="center",
            fontsize=8.5, color=RADIUS)
    ax.axhline(0.28, color=GRID, linewidth=0.85)
    ax.set_xlabel("Certificate distance (complex-port Frobenius scale)")

    axis_rows = ax.get_yaxis_transform()
    ax.text(-0.44, 4.87, "Maximum over the\nfixed 16-block weight line",
            transform=axis_rows, fontsize=9, weight="bold", va="center")
    ax.text(-0.44, -0.09, "Observation-selected weight",
            transform=axis_rows, fontsize=9, weight="bold", va="center")
    for x, title in ((1.14, "Strict lower"), (1.41, "Strict upper")):
        ax.text(x, 4.87, title, transform=axis_rows, fontsize=9,
                weight="bold", ha="center", va="center")
    precision = Decimal("0.00000000001")
    for y, _, lower, upper in intervals:
        verdict = "PASS" if lower > radius else "BELOW" if upper < radius else "UNRESOLVED"
        color = VERDICTS[verdict][1]
        mid = (lower + upper) / 2
        ax.errorbar(mid, y, xerr=[[mid - lower], [upper - mid]], fmt="none",
                    linewidth=2, capsize=3, color=color, zorder=4)
        for x, bound, rounding in ((1.14, lower, ROUND_FLOOR), (1.41, upper, ROUND_CEILING)):
            label = Decimal.from_float(bound).quantize(precision, rounding=rounding)
            ax.text(x, y, str(label), transform=axis_rows, fontsize=9,
                    color=color, ha="center", va="center")

    fig.text(0.04, 0.105,
             "Strict exclusion requires the lower bound to exceed the same noise radius.",
             fontsize=9.2)
    fig.text(0.04, 0.054,
             f"Strict exclusions on transfer · {transfer['objects']} materials, {transfer['views']} paired views"
             f" · Frozen gradient weight {len(one['passing_views'])}/{transfer['views']}"
             f" · Five archived weights {len(controls['passing_views'])}/{transfer['views']}"
             f" ({len(controls['passing_objects'])} materials at {noise_label})",
             fontsize=9.2)
    save(fig, "same_observation")


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--aperture", action="store_true", help="also plot passive receiver masks")
    parser.add_argument("--diagnostic", action="store_true", help="plot the same-observation diagnostic")
    args = parser.parse_args()
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "normal",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#a2afb7",
        "text.color": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "pdf.fonttype": 42,
    })
    if args.diagnostic:
        plot_diagnostic()
    if args.aperture or not args.diagnostic:
        data = read_cohorts()
        plot_segment(data)
        if args.aperture:
            plot_aperture(data)


if __name__ == "__main__":
    main()
