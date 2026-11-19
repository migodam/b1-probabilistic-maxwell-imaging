"""Make the paper figure for the bandwidth and computation trade-offs."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from data14 import ROOT


def read(name):
    return json.loads((ROOT / "results" / name).read_text())


COLORS = {
    "blue": "#176187", "green": "#2A786A", "red": "#C25E46",
    "purple": "#7568A1", "black": "#28343B", "gray": "#8A9297",
}
SWEEPS = (
    ("one_minus", "One lower tone", COLORS["blue"], "v"),
    ("one_plus", "One upper tone", COLORS["green"], "^"),
    ("split", "Two symmetric tones", COLORS["red"], "o"),
)
CONTROLS = (
    ("old_only", "root3", "Old only", "gray", "x"),
    ("split_0.02", "single", "Fixed, one start", "blue", "v"),
    ("split_0.02", "ordinary3", "Fixed, 3 ordinary starts", "blue", "^"),
    ("split_0.02", "ordinary8", "Fixed, 8 ordinary starts", "blue", "8"),
    ("split_0.02", "root3", "Fixed, 3 candidate starts", "red", "o"),
    ("uniform_0.02", "root3", "Uniform 4-tone", "green", "D"),
    ("split_0.02", "cumulative3", "Low-to-high cumulative", "purple", "s"),
    ("uniform_0.02", "fbch2", "2-bin cumulative, λ²", "purple", "P"),
    ("candidate_gap", "root3", "Candidate separation", "green", "*"),
    ("fisher", "root3", "Local Fisher", "black", "h"),
    ("random_at_cap", "root3", "Random at cap", "gray", "d"),
    ("region", "root3", "Region rule", "black", "X"),
    ("region", "ordinary8", "Region, 8 ordinary starts", "black", "p"),
)


def plot(rows):
    cfg = {"font.size": 10., "axes.spines.right": False, "axes.spines.top": False,
           "axes.titlepad": 10., "axes.labelpad": 6., "pdf.fonttype": 42}
    with plt.rc_context(cfg):
        fig, axes = plt.subplots(3, 3, figsize=(14, 10.6), layout="constrained")
        titles = {"spectral": "R12 noisy pairs · 36",
                  "natural": "R10 natural targets · 18", "aperture": "R10 finite aperture · 8"}
        for j, (cohort, title) in enumerate(titles.items()):
            band, comp, timer = axes[:, j]
            records = {(r["action_or_policy"], r["mode"]): r for r in rows if r["cohort"] == cohort}
            band.set_title(title, loc="left", fontweight="bold")
            for prefix, label, color, marker in SWEEPS:
                selected = [r for (action, mode), r in records.items()
                            if action.startswith(prefix + "_") and mode == "root3"
                            and (r["label"] == "sweep" or action == "split_0.02")]
                selected.sort(key=lambda r: r["bandwidth_pct"])
                band.plot([r["bandwidth_pct"] for r in selected], [r["rmse_p90"] for r in selected],
                          marker=marker, ms=4.5, lw=1.5, color=color, label=label)
            wide = [records[f"split_{offset}", "root3"] for offset in ("0.02", "0.04", "0.08")]
            band.plot([r["bandwidth_pct"] for r in wide], [r["rmse_p90"] for r in wide],
                      marker="o", ms=4.5, linestyle="--", lw=1.5, color=COLORS["red"],
                      label="Wider-band control")
            ref = records["old_only", "root3"]["rmse_p90"]
            band.axhline(ref, color=COLORS["gray"], lw=1, ls=":", label="Old-only error")
            band.set_xscale("log")
            band.set_yscale("log")
            band.set_xlabel("Relative frequency span (%)")
            if j == 0:
                band.set_ylabel("Material RMSE, 90th percentile")
            for axis, field, name in ((comp, "all_factors_mean", "Mean Maxwell factorizations"),
                                      (timer, "total_time_median_s", "Median design + fit time (s)")):
                frontier = []
                best = float("inf")
                for slug, mode, label, color, marker in CONTROLS:
                    if (slug, mode) not in records or (slug == "region" and axis is comp):
                        continue
                    point = records[slug, mode]
                    x, y = point[field], point["rmse_p90"]
                    axis.plot(x, y, marker=marker, ms=7, ls="", color=COLORS[color], label=label)
                    if slug != "old_only":
                        frontier.append((x, y))
                saved = []
                for x, y in sorted(frontier):
                    if y < best:
                        saved.append((x, y))
                        best = y
                if saved:
                    axis.step(*zip(*saved), where="post", color=COLORS["gray"], alpha=.65,
                              ls="--", lw=1.)
                axis.set_yscale("log")
                axis.set_xlabel(name)
                if j == 0:
                    axis.set_ylabel("Material RMSE, 90th percentile")
            timer.set_xscale("log")
            for axis in (band, comp, timer):
                axis.grid(which="major", color="#DCE0E3", lw=.7, alpha=.75)
                axis.set_axisbelow(True)
        axes[0, 1].legend(loc="lower center", bbox_to_anchor=(.5, 1.17), ncol=5, frameon=False)
        axes[2, 1].legend(loc="upper center", bbox_to_anchor=(.5, -.24), ncol=4, frameon=False)
        for suffix in ("png", "pdf"):
            fig.savefig(ROOT / "figures" / f"measure_pareto.{suffix}", dpi=220, bbox_inches="tight")
        plt.close(fig)


def plot_validation(rows):
    palettes = plt.get_cmap("tab10")
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.4), sharey=True, layout="constrained")
    for ax, kind, title in zip(axes, ("half_first", "full_split"),
                               ("First half exposure", "Both half exposures")):
        selected = {r["case_id"]: r for r in rows if r["label"] == kind}
        for i in range(8):
            lo, hi = [selected[f"val6-{i:02d}-{level}"] for level in ("low", "high")]
            values = [lo["rmse"], hi["rmse"]]
            label = "One geometry and material" if i == 0 else None
            ax.plot([0, 1], values, color=palettes(i), lw=1.4, alpha=.83, label=label)
            for j, row in enumerate((lo, hi)):
                if row["common_data_compatible"]:
                    ax.scatter(j, row["rmse"], color=palettes(i), s=23, zorder=3)
                else:
                    ax.scatter(j, row["rmse"], marker="x", color="#BA332D", s=54, zorder=3,
                               label="Outside common-data contract")
        ax.axhline(.5, color="#525961", ls="--", lw=1., label="0.5 RMSE target")
        ax.set(xlim=(-.26, 1.26), yscale="log", title=title)
        ax.set_xticks([0, 1], ["Low noise\nσ = 3.00 × 10⁻⁷", "High noise\nσ = 6.90 × 10⁻⁵"])
        ax.spines[["right", "top"]].set_visible(False)
        ax.grid(axis="y", which="major", color="#E1E4E7", lw=.7)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Material RMSE")
    axes[1].legend(loc="lower center", fontsize=8., frameon=False)
    for suffix in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"measure_paired_noise.{suffix}", dpi=210, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    results = read("measure_summary.json")
    if (ROOT / "results/measure_summary_region.json").exists():
        results += read("measure_summary_region.json")
    plot(results)
    if (ROOT / "results/validation_scores.json").exists():
        plot_validation(read("validation_scores.json"))
