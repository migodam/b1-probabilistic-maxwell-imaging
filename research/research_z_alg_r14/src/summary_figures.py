"""Draw R14's observed diameter comparison and known-gain continuous-band witness."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
import numpy as np
from scipy.special import j0

from data14 import FrequencyTM, ROOT, load_cases


PAPER = "#fbf8f0"
INK = "#203638"
MUTED = "#536560"
GRID = "#e4ded0"
TRACE = "#87958c"
BLUE = "#345d73"
TEAL = "#1d776b"
BRICK = "#ab5141"
STAGES = ("first_half", "both_halves")
TICKS = (.05, .1, .5, 1., 5., 10., 25.)
RESULTS = ROOT / "results"
STYLE = {
    "figure.facecolor": PAPER, "axes.facecolor": PAPER,
    "savefig.facecolor": PAPER, "axes.edgecolor": GRID,
    "text.color": INK, "axes.labelcolor": INK, "axes.titlecolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.titleweight": "bold", "axes.titlepad": 8,
    "pdf.fonttype": 42,
}


def load(file):
    return json.loads((RESULTS / file).read_text())


def save(fig, name):
    for extension in ("pdf", "png"):
        fig.savefig(ROOT / "figures" / f"{name}.{extension}", dpi=250)
    plt.close(fig)


def diameter_axis(ax, label=False):
    ax.set_yscale("log")
    ax.set_ylim(.029, 36.)
    ax.yaxis.set_major_locator(FixedLocator(TICKS))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.grid(axis="y", color=GRID, linewidth=.7)
    ax.axhline(.5, color=BRICK, ls=(0, (4, 3)), lw=1.05, zorder=1)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", length=0)
    if label:
        ax.set_ylabel("Material diameter upper bound (RMS)")


def histories(ax, forecasts, prefix, title, label=False):
    names = sorted(name for name in forecasts if name.startswith(prefix))
    data = np.array([[forecasts[name]["D_upper"],
                      *[load(f"posterior/{name}.json")[stage]["D_upper"] for stage in STAGES]]
                     for name in names])
    offsets = np.linspace(-.07, .07, len(data))
    for row, shift in zip(data, offsets):
        ax.plot(np.arange(3) + shift, row, color=TRACE, alpha=.25, lw=.8, zorder=2)
    for stage, color, symbol in zip(range(3), (BLUE, TEAL, BRICK), ("o", "o", "s")):
        ax.scatter(stage + offsets, data[:, stage], c=color, marker=symbol, s=17,
                   alpha=.74, edgecolors=PAPER, linewidths=.25, zorder=3)
    passed = (data <= .5).sum(axis=0)
    names = ("Before", "Half · −2%", "Both · ±2%")
    ax.set_xticks(range(3), [f"{name}\n{n}/{len(data)}" for name, n in zip(names, passed)])
    ax.set_xlim(-.25, 2.25)
    ax.set_title(title, loc="left")
    diameter_axis(ax, label=label)


def new_pairs(ax):
    stages = []
    for stage in STAGES:
        paired = []
        for i in range(8):
            paired.append([load(f"validation_regions/val6-{i:02d}-{level}.json")["posterior"][
                stage]["D_upper"] for level in ("low", "high")])
        stages.append(np.asarray(paired))
    for stage, data in enumerate(stages):
        positions = np.array([0., .63]) + stage * 1.70
        jitter = np.linspace(-.053, .053, len(data))
        for row, shift in zip(data, jitter):
            ax.plot(positions + shift, row, color=TRACE, lw=.95, alpha=.58, zorder=2)
        for j, color in enumerate((TEAL, BRICK)):
            ax.scatter(positions[j] + jitter, data[:, j], c=color, s=24,
                       edgecolors=PAPER, linewidths=.35, zorder=3)
    counts = np.asarray([(a <= .5).sum(axis=0) for a in stages]).ravel()
    labels = ("First half · low", "First half · high", "Both halves · low", "Both halves · high")
    ax.set_xticks([0., .63, 1.70, 2.33], [f"{name}\n{c}/8" for name, c in zip(labels, counts)])
    ax.set_xlim(-.27, 2.6)
    ax.set_title("New materials · eight paired noise contrasts", loc="left")
    diameter_axis(ax, label=True)
    ax.annotate("0.5 RMS target", (.965, .5), xycoords=("axes fraction", "data"),
                xytext=(0, 4), textcoords="offset points", ha="right", va="bottom",
                fontsize=9, color=BRICK)


def sequential_energy(ax):
    fixed = next(r for r in load("measure_summary_half.json")
                 if (r["cohort"], r["label"], r["mode"]) == ("spectral", "full_split", "ordinary8"))
    early = next(r for r in load("measure_summary_sequential.json")
                 if (r["cohort"], r["mode"]) == ("spectral", "ordinary8"))
    positions = [1., 0.]
    for y, result, color, name in zip(positions, (fixed, early), (BLUE, TEAL),
                                      ("Fixed two blocks", "Certified early stop")):
        energy = result["mean_energy"]
        ax.text(0, y + .29, name, ha="left", va="bottom", color=INK)
        ax.barh(y, energy, height=.30, color=color, zorder=2)
        ax.text(energy + .03, y, f"{energy:.3g}", ha="left", va="center",
                color=color, weight="bold")
        ax.text(0, y - .26, f"90th-percentile RMSE {result['rmse_p90']:.3g}",
                ha="left", va="top", color=MUTED, fontsize=9)
    ax.set_title("Mean added energy · R12", loc="left")
    ax.set_xlim(0, 1.20)
    ax.set_ylim(-.61, 1.58)
    ax.set_xticks((0, .5, 1.))
    ax.set_xlabel("One unit = fixed two-block exposure")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", color=GRID, linewidth=.7)
    ax.set_axisbelow(True)


def overview():
    forecasts = load("region_actions.json")
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(12.2, 7.5), layout="constrained")
        grid = fig.add_gridspec(2, 6, height_ratios=(1.07, 1.0))
        fig.suptitle("When the next measurement can be skipped", x=.02, ha="left",
                     fontsize=19, weight="bold")
        for j, (prefix, title) in enumerate((("spectral", "R12 noisy reads · 36"),
                                             ("fresh6", "R10 natural cases · 18"),
                                             ("confirm6", "R10 finite aperture · 8"))):
            ax = fig.add_subplot(grid[0, 2*j:2*j+2])
            histories(ax, forecasts, prefix, title, label=j == 0)
            if j:
                ax.tick_params(labelleft=False)
        new_pairs(fig.add_subplot(grid[1, :4]))
        sequential_energy(fig.add_subplot(grid[1, 4:]))
        save(fig, "r14_overview")


def coherent_curve(case, materials, limit):
    model = FrequencyTM(case.geometry)
    positions = np.asarray(case.geometry["sources"])
    separations = np.linalg.norm(positions[:, None] - positions[None, :], axis=-1)
    offsets = np.linspace(-.02, .02, 721)
    ratio = []
    for delta in offsets:
        diff = model.response(materials[0], float(delta))[0] - model.response(materials[1], float(delta))[0]
        gram = j0(case.geometry["k"] * (1. + delta) * separations)
        ratio.append(np.sqrt(3 * np.vdot(diff, np.linalg.solve(gram, diff)).real) / limit)
    return 100. * offsets, np.asarray(ratio)


def continuous_band():
    witness = load("continuous_band/fresh6-16_known_coherent.json")
    case = next(c for c in load_cases() if c.case_id == witness["case_id"])
    materials = np.asarray(witness["materials"])
    limit = witness["threshold_2r"]
    x, ratio = coherent_curve(case, materials, limit)
    intervals = np.asarray(witness["terminal_intervals"])
    mid, width, upper = intervals[:, 0], intervals[:, 1], intervals[:, 3] / limit
    starts = np.column_stack((100 * (mid - width), upper))
    ends = np.column_stack((100 * (mid + width), upper))
    segments = np.stack((starts, ends), axis=1)
    with plt.rc_context(STYLE):
        fig, (left, right) = plt.subplots(1, 2, figsize=(11.7, 4.52),
                                          width_ratios=(.85, 1.3), layout="constrained")
        fig.suptitle("The same readings can fit two materials at known unit gain",
                     x=.02, ha="left", fontsize=16.5, weight="bold")
        offsets = np.arange(1, 7)
        for j, color in enumerate((BLUE, BRICK)):
            left.bar(offsets + (-.18 if j == 0 else .18), materials[j], width=.32,
                     color=color, label=f"Material {'AB'[j]}", zorder=3)
        left.set_title(f"{case.case_id}  ·  separation ≥ {witness['rms_lower']:.4f} RMS", loc="left")
        left.set(xlabel="Material cell", ylabel="Material contrast", xticks=offsets, ylim=(0, 33))
        left.grid(axis="y", color=GRID, linewidth=.7)
        left.set_axisbelow(True)
        left.legend(loc="upper left", frameon=False, fontsize=9)
        right.plot(x, ratio, color=TEAL, label="Full-wave numerical curve", lw=2.0, zorder=2)
        right.add_collection(LineCollection(segments, colors=BLUE, linewidths=1.15,
                                            label=f"Arb upper bounds ({len(intervals)} intervals)", zorder=3))
        right.axhline(1., color=BRICK, lw=1.25, ls=(0, (5, 3)), label="Noise balls touch at 1", zorder=1)
        right.set_title("Any coherent feed on the three original ports", loc="left")
        right.set(xlabel="Frequency offset from the old carrier (%)",
                  ylabel="Worst coherent gap / two noise radii", xlim=(-2., 2.),
                  ylim=(min(.65, float(ratio.min()) - .015), 1.055), xticks=(-2., -1., 0., 1., 2.))
        right.grid(color=GRID, linewidth=.7)
        right.set_axisbelow(True)
        right.legend(loc="upper right", bbox_to_anchor=(1., .865), frameon=False, fontsize=8.8)
        right.annotate(f"Global upper bound {witness['upper_ratio']:.3f}",
                       (.025, upper.max()), xycoords=("axes fraction", "data"),
                       xytext=(0, 7), textcoords="offset points", color=BLUE,
                       fontsize=9.2, va="bottom")
        save(fig, "continuous_band")


if __name__ == "__main__":
    overview()
    continuous_band()
