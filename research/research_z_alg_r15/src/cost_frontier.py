"""Plot measured reconstruction error against acquisition costs on two cohorts."""

from dataclasses import dataclass
from statistics import mean, median

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from trade_figures import INK, panel_style, read_csv, read_json, save

METHODS = (
    ("fft_tsom", "FFT-TSOM", "#197f88", "o"),
    ("projected_gn", "Projected GN · shared penalty", "#b65537", "^"),
    ("calibrated_gn", "Projected GN · stronger penalty", "#4e5d94", "s"),
)
EXTENSION_NAMES = {
    "pilot_alone": "16 GHz pilot",
    "add_tx": "Add transmitters",
    "switch_rx": "Switch receivers",
    "parallel_rx": "Parallel receivers",
    **{f"add_{f}GHz": f"Add {f} GHz" for f in (6, 12, 20, 24, 30)},
}


@dataclass(frozen=True)
class Point:
    name: str
    reads: int
    shots: int
    chains: int
    rmse: float
    info_gain: float | None = None

    @property
    def objectives(self):
        return self.reads, self.shots, self.chains, self.rmse


def pareto(points):
    """Keep a schedule unless another is no worse on all three costs and RMSE."""
    frontier = set()
    for point in points:
        obj = point.objectives
        if not any(all(a <= b for a, b in zip(other.objectives, obj, strict=True))
                   and any(a < b for a, b in zip(other.objectives, obj, strict=True))
                   for other in points):
            frontier.add(point)
    return frontier


def main_cohort():
    standard = read_csv("inversion_metrics.csv") + read_csv("inversion_metrics_task_design.csv")
    calibrated = read_csv("inversion_metrics_gn_calibrated.csv")
    expected = {r["case"] for r in standard
                if r["split"] == "independent" and r["case"] != "test_empty"}
    by_method = {}
    for key, _, _, _ in METHODS:
        rows = calibrated if key == "calibrated_gn" else [r for r in standard if r["method"] == key]
        grouped = {}
        for row in rows:
            if row["case"] in expected:
                grouped.setdefault(row["arm"], []).append(row)
        points = []
        for arm, data in grouped.items():
            if len(data) != len(expected) or {r["case"] for r in data} != expected:
                raise ValueError(f"Incomplete main cohort for {key}, {arm}")
            first = data[0]
            points.append(Point(arm, int(first["complex_samples"]), int(first["shots"]),
                                int(first["receivers"]), median(float(r["rmse_eps"]) for r in data)))
        by_method[key] = points
    return by_method


def nested_cohort():
    scores = read_csv("nested_scores_holdout72.csv")
    metadata = read_json("nested_holdout72.json")["rows"]
    expected = {r["case"] for r in scores if r["case"] != "blind_08"}
    points = []
    for row in metadata:
        key = row["scheme"]
        data = [r for r in scores if r["scheme"] == key]
        cases = {r["case"] for r in data}
        if len(data) != len(expected) + 1 or cases != (expected | {"blind_08"}):
            raise ValueError(f"Incomplete paired holdout for {key}")
        aggregate = mean(float(r["rmse_eps"]) for r in data if r["case"] in expected)
        chains = 16 if key == "parallel_rx" else 8
        points.append(Point(key, row["samples"], row["shots"], chains, aggregate, row["net_nats"]))
    return points


def plot_main(ax, by_method):
    ax.set_xscale("log", base=2)
    ax.set_xlim(12.4, 660)
    ax.set_ylim(.022, .161)
    ax.set_title("Original cohort · material inversions", color=INK, fontweight="bold",
                 loc="left", fontsize=14, pad=10)
    ax.set_ylabel("Median material RMSE · 8 nonempty objects")
    panel_style(ax)
    read_costs = sorted({(p.reads, p.shots, p.chains) for rows in by_method.values() for p in rows})
    ax.set_xticks([n for n, _, _ in read_costs], [f"{n}\n{shots} / {chains}" for n, shots, chains in read_costs])
    ax.set_xlabel("Complex reads (shots / simultaneous Rx chains)", fontsize=10.5, labelpad=8)

    legend = []
    positions = {
        ("fft_tsom", 16): (0, 8), ("fft_tsom", 64): (0, 8), ("fft_tsom", 512): (0, 8),
        ("projected_gn", 16): (0, -18), ("projected_gn", 64): (-10, 8),
        ("calibrated_gn", 64): (6, -18), ("calibrated_gn", 512): (0, -19),
    }
    for key, label, color, marker in METHODS:
        points = by_method[key]
        keep = pareto(points)
        order = sorted(keep, key=lambda p: p.reads)
        ax.plot([p.reads for p in order], [p.rmse for p in order], linewidth=1.3, color=color,
                alpha=.64, zorder=2)
        for row in points:
            on_frontier = row in keep
            ax.scatter(row.reads, row.rmse, s=65 if on_frontier else 35, marker=marker,
                       facecolor=color if on_frontier else "white", edgecolor=color,
                       alpha=1 if on_frontier else .50, lw=1.3, zorder=4 if on_frontier else 3)
            if on_frontier:
                ax.annotate(f"{row.rmse:.3f}", (row.reads, row.rmse), xytext=positions.get((key, row.reads), (0, 7)),
                            textcoords="offset points", ha="center", va="bottom",
                            fontsize=9.3, color=color, fontweight="bold")
        legend.append(Line2D([0], [0], marker=marker, markersize=5, color=color, label=label, linewidth=1.5))
    ax.legend(handles=legend, loc="upper right", frameon=False, fontsize=9.5, borderaxespad=.5)


def plot_nested(ax, table, points):
    keep = pareto(points)
    y_values = list(range(len(points) - 1, -1, -1))
    ymax, ymin = len(points) + .28, -.65
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(.030, .0525)
    ax.set_title("72-grid holdout · extend the same pilot layout", color=INK,
                 fontweight="bold", loc="left", fontsize=14, pad=10)
    ax.set_yticks(y_values, [EXTENSION_NAMES[p.name] for p in points], fontsize=10.3)
    ax.tick_params(axis="y", length=0, pad=7)
    ax.set_xlabel("Mean material RMSE · 8 nonempty objects", fontsize=10.5, labelpad=8)
    ax.set_xticks((.03, .035, .040, .045, .050))
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis="x", linewidth=.75, color="#dfe7eb")

    table.set_xlim(0, 1)
    table.set_ylim(ymin, ymax)
    table.axis("off")
    columns = (.08, .32, .56, .86)
    for x, label in zip(columns, ("Reads", "Shots", "Rx", "Δ info"), strict=True):
        table.text(x, ymax - .35, label, ha="center", va="center", color=INK, fontsize=9.8,
                   fontweight="bold")

    for y, point in zip(y_values, points, strict=True):
        is_frontier = point in keep
        color = "#126f78" if is_frontier else "#667683"
        if is_frontier:
            ax.axhspan(y - .40, y + .40, color="#eaf4f3", zorder=0)
            table.axhspan(y - .40, y + .40, color="#eaf4f3", zorder=0)
        ax.scatter(point.rmse, y, s=71 if is_frontier else 51, color=color,
                   marker="o" if is_frontier else "D", edgecolor="white", lw=.85, zorder=3)
        ax.annotate(f"{point.rmse:.4f}", (point.rmse, y), xytext=(7, 0),
                    textcoords="offset points", ha="left", va="center", fontsize=9.7,
                    color=color, fontweight="bold" if is_frontier else "normal")
        info = f"{point.info_gain:+.1f}" if point.name != "pilot_alone" else "0.0"
        numbers = (str(point.reads), str(point.shots), str(point.chains), info)
        for x, value in zip(columns, numbers, strict=True):
            table.text(x, y, value, ha="center", va="center", color=color, fontsize=10.1,
                       fontweight="bold" if is_frontier else "normal")


def main():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11, "text.color": INK,
        "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "pdf.fonttype": 42,
    })
    fig = plt.figure(figsize=(14.9, 6.3))
    fig.text(.065, .951, "Acquisition costs and material error",
             fontsize=24, color=INK, fontweight="bold", va="top")
    fig.text(.065, .880,
             "Total radiated energy E = 1 per schedule. Each shot uses E divided by the shot count. "
             "Reference noise is 2%.",
             fontsize=11, color=INK, va="top")
    layout = fig.add_gridspec(1, 2, left=.078, right=.976, top=.755, bottom=.225, width_ratios=(1, 1.25),
                              wspace=.31)
    ax_main = fig.add_subplot(layout[0])
    nested = layout[1].subgridspec(1, 2, width_ratios=(1.64, 1.), wspace=.09)
    ax_nested, table_nested = [fig.add_subplot(nested[j]) for j in range(2)]
    plot_main(ax_main, main_cohort())
    plot_nested(ax_nested, table_nested, nested_cohort())
    fig.text(.078, .082, "Frontiers use reads, shots, simultaneous Rx chains and RMSE. "
             "Filled markers flag the frontier within each method or cohort.",
             fontsize=10.2, color=INK)
    fig.text(.078, .048, "Switched Rx uses eight chains across 16 positions. Δ info is the local Gaussian information gain "
             "over the pilot in nats (known gain, prior SD 0.30).",
             fontsize=10.2, color=INK)
    save(fig, "cost_frontier")


if __name__ == "__main__":
    main()
