"""Publication figures for the fixed-energy acquisition and material-witness trade-offs."""

import argparse
import csv
import json
from pathlib import Path
from statistics import mean

import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm, TwoSlopeNorm
from matplotlib.ticker import PercentFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
COUNTS = (1, 2, 4, 8)
RECEIVER_COUNTS = (1, 2, 4, 8, 16)
BUDGETED = ("f1_t8_r8", "f2_t4_r8", "f4_t2_r8", "f8_t1_r8")
SCHEDULE_TEXT = ("1 F × 8 Tx", "2 F × 4 Tx", "4 F × 2 Tx", "8 F × 1 Tx")
WITNESS_SCHEDULES = ("f2_t4_r8", "f4_t2_r8")
SCHEDULE_COLOR = {"f2_t4_r8": "#167e89", "f4_t2_r8": "#b65537"}
INK = "#19252e"
ACCENT = "#a04c22"
NOISE_NOTE = "2% reference noise"
PANEL_HEADING = dict(fontsize=13, color=INK, fontweight="bold")


def read_csv(name):
    with (RESULTS / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def read_json(name):
    with (RESULTS / name).open() as stream:
        return json.load(stream)


def save(fig, name):
    stem = FIGURES / name
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=200, facecolor="white")
    plt.close(fig)


def panel_style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis="y", linewidth=.7, color="#dfe7eb")
    ax.tick_params(length=3)


def mean_rmse(rows, design, variable, value):
    return mean(
        float(r["rmse_eps"]) for r in rows
        if r["arm"] == design and r[variable] == value and
        r["case"].startswith("test_") and r["case"] != "test_empty"
    )


def dof_figure():
    information = {r["design"]: r for r in read_csv("dof_trade.csv")}
    fullwave = read_csv("inversion_metrics.csv")
    sparse = read_csv("sparse_scores.csv")
    fig = plt.figure(figsize=(14.8, 8.2), layout="constrained")
    outer = fig.add_gridspec(3, 1, height_ratios=(.66, 2.55, 2.5))
    heading = fig.add_subplot(outer[0])
    heading.axis("off")
    heading.text(0, .96, "Frequency and aperture change what can be estimated",
                 fontsize=23, color=INK, fontweight="bold", va="top")
    heading.text(0, .38,
                 "Fixed total radiated energy E = 1. "
                 "F × Tx gives the shot count and F × Tx × Rx gives the complex read count.",
                 fontsize=11, color=INK, va="top")
    heading.text(0, .12,
                 "Outlined schedules use 8 shots and 64 complex reads. "
                 "Effective modes use a per-pixel contrast prior SD of 0.30 and 2% reference noise.",
                 fontsize=10.5, color=INK, va="top")

    top = outer[1].subgridspec(1, 6, width_ratios=(1, 1, 1, 1, 1, .15), wspace=.06)
    max_modes = max(float(information[f"f{f}_t{t}_r{r}"]["eff_48"])
                    for f in COUNTS for t in COUNTS for r in RECEIVER_COUNTS)
    normalizer = PowerNorm(.62, 0, np.ceil(max_modes / 20) * 20)
    cmap = plt.get_cmap("Blues")
    for r_idx, receivers in enumerate(RECEIVER_COUNTS):
        ax = fig.add_subplot(top[r_idx])
        arr = np.array([[float(information[f"f{f}_t{t}_r{receivers}"]["eff_48"])
                         for t in COUNTS] for f in COUNTS])
        visual = ax.imshow(arr, cmap=cmap, norm=normalizer, interpolation="nearest")
        for i, (f, values) in enumerate(zip(COUNTS, arr, strict=True)):
            for j, (t, val) in enumerate(zip(COUNTS, values, strict=True)):
                color = "white" if (normalizer(val) > .59) else INK
                ax.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=11, color=color)
                if receivers == 8 and f * t == 8:
                    ax.add_patch(plt.Rectangle((j - .46, i - .46), .92, .92, fill=False,
                                               color=ACCENT, linewidth=2.3))
        ax.set_title(f"{receivers} Rx", **PANEL_HEADING)
        ax.set_xticks(range(4), COUNTS)
        ax.set_xlabel("Tx per frequency")
        ax.set_yticks(range(4), COUNTS if r_idx == 0 else ("",) * 4)
        if r_idx == 0:
            ax.set_ylabel("Frequencies")
        ax.tick_params(length=0)
        ax.spines[:].set_visible(False)
    cax = fig.add_subplot(top[-1])
    cb = fig.colorbar(visual, cax=cax, ticks=(0, 20, 50, 100, 150, 200))
    cb.set_label("Effective material modes", fontsize=11)
    cb.outline.set_visible(False)

    lower = outer[2].subgridspec(1, 3, wspace=.09)
    x = np.arange(4)
    metrics = [information[key] for key in BUDGETED]
    ax_m, ax_l, ax_r = [fig.add_subplot(lower[j]) for j in range(3)]
    ax_m.bar(x, [float(r["eff_48"]) for r in metrics], width=.54, color="#176c86")
    ax_m.bar_label(ax_m.containers[0], fmt="%.1f", padding=3, color=INK, fontsize=11)
    ax_m.set_ylim(0, 112)
    ax_m.set_ylabel("Effective material modes")
    ax_m.set_title("Local image information", **PANEL_HEADING, loc="left")

    loc_specs = (
        ("location_sd_known_model_mm", "Fit amplitude", "#146d7d", "o"),
        ("location_sd_profiled_model_mm", "Fit amplitude and grid discrepancy", "#8f5132", "D"),
    )
    for key, label, color, marker in loc_specs:
        ax_l.plot(x, [1000 * float(r[key]) for r in metrics], label=label,
                  color=color, marker=marker, markersize=6, linewidth=2)
    ax_l.set_ylim(5.75, 20)
    ax_l.set_ylabel("Worst-direction bound (µm)")
    ax_l.set_title("Local flaw position", **PANEL_HEADING, loc="left")
    ax_l.legend(loc="upper right", fontsize=9.3, frameon=False)

    fit_specs = (
        (fullwave, "method", "fft_tsom", "FFT-TSOM", "#556b7e", "s"),
        (fullwave, "method", "projected_gn", "Projected GN", "#ab5936", "D"),
        (sparse, "grid", "48", "Sparse linear", "#136d60", "o"),
    )
    for rows, column, value, label, color, marker in fit_specs:
        ax_r.plot(x, [mean_rmse(rows, design, column, value) for design in BUDGETED],
                  label=label, color=color, marker=marker, markersize=6, linewidth=2)
    ax_r.set_ylim(.020, .156)
    ax_r.set_ylabel("Mean permittivity RMSE")
    ax_r.set_title("Eight nonempty test objects", **PANEL_HEADING, loc="left")
    ax_r.legend(loc="upper right", fontsize=9.3, frameon=False)
    for ax in (ax_m, ax_l, ax_r):
        panel_style(ax)
        ax.set_xticks(x, SCHEDULE_TEXT, fontsize=10)

    save(fig, "dof_trade")


def risk_label(p):
    if p >= .001:
        return f"{100*p:.1f}%"
    if not p:
        return "Below precision"
    mantissa, exponent = f"{100*p:.2e}".split("e")
    return rf"${mantissa}\times10^{{{int(exponent)}}}\%$"


def risk_rows(case, transfer):
    return (
        {case["hard_design"]: case["fine_hard"]["equal_prior_error"],
         case["easy_design"]: case["fine_easy"]["equal_prior_error"]},
        {transfer["hard_design"]: transfer["grid"]["72"]["hard"]["equal_prior_error"],
         transfer["easy_design"]: transfer["grid"]["72"]["easy"]["equal_prior_error"]},
    )


def error_panel(ax, values, grid, show_title, show_labels):
    ax.set_xlim(0, .60)
    ax.set_ylim(-.65, 1.55)
    ax.set_xticks((0, .10, .20, .30, .40, .50))
    ax.xaxis.set_major_formatter(PercentFormatter(1., decimals=0))
    ax.set_yticks((1, 0), ("2 F × 4 Tx", "4 F × 2 Tx") if show_labels else ("", ""))
    ax.tick_params(axis="y", length=0, pad=7)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_axisbelow(True)
    if not show_title:
        ax.set_xlabel("Minimum equal-prior error")
    ax.grid(axis="x", color="#e1e9eb", linewidth=.7)
    ax.axvline(.5, ls="--", linewidth=1.1, color="#889299")
    ax.text(.5, 1.41, "Chance", fontsize=9, ha="center", va="bottom", color="#68747e")
    if show_title:
        ax.set_title(f"{grid} × {grid} full-wave grid", **PANEL_HEADING, loc="left", pad=12)
    for y, design in zip((1, 0), WITNESS_SCHEDULES, strict=True):
        val = values[design]
        color = SCHEDULE_COLOR[design]
        ax.hlines(y, 0, val, lw=4, color=color, zorder=3)
        ax.scatter(val, y, color=color, edgecolor="white", s=65, lw=1.05,
                   zorder=4, clip_on=False)
        ax.annotate(risk_label(val), (val, y), xytext=(8, 0), textcoords="offset points",
                    ha="left", va="center", color=color, fontweight="bold", fontsize=12)


def witness_figure():
    modes = read_json("material_weak_modes.json")["cases"]
    transfers = {(r["hard_design"], r["easy_design"]): r
                 for r in read_json("witness_transfer.json")["cases"]}
    pairs = (WITNESS_SCHEDULES, WITNESS_SCHEDULES[::-1])
    witnesses = [
        (next(case for case in modes if (case["hard_design"], case["easy_design"]) == pair),
         transfers[pair])
        for pair in pairs
    ]

    fig = plt.figure(figsize=(13.8, 7.25))
    outer = fig.add_gridspec(3, 1, left=.055, right=.98, top=.97, bottom=.16,
                             height_ratios=(.78, 2.25, 2.25), hspace=.42)
    heading = fig.add_subplot(outer[0])
    heading.axis("off")
    heading.text(0, .9, "Two material pairs favor different MIMO schedules",
                 fontsize=23, color=INK, fontweight="bold", va="top")
    heading.text(0, .37, "Both schedules span 2–30 GHz using E = 1, 8 shots and 64 complex reads.",
                 fontsize=12, color=INK, va="top")
    heading.text(0, .04,
                 "The two pairs were selected on a 24 × 24 grid and then evaluated unchanged on finer grids.",
                 fontsize=11, color=INK, va="top")

    color_range = max(case["max_cell_difference"] for case, _ in witnesses)
    norm = TwoSlopeNorm(vcenter=0, vmin=-color_range, vmax=color_range)
    for i, (case, transfer) in enumerate(witnesses):
        frame = outer[i + 1].subgridspec(1, 3, width_ratios=(.90, 1.45, 1.45), wspace=.14)
        ax_map, ax_48, ax_72 = [fig.add_subplot(frame[j]) for j in range(3)]
        difference = np.zeros(24**2)
        difference[case["pixel_ids"]] = np.asarray(case["plus"]) - case["minus"]
        displayed = ax_map.imshow(difference.reshape(24, 24), extent=(-12, 12, -12, 12),
                                  origin="lower", cmap="RdBu_r", norm=norm, interpolation="nearest")
        ax_map.set_title(f"Pair {'AB'[i]} material difference", **PANEL_HEADING, loc="left")
        if i == 1:
            ax_map.set_xlabel("x (mm)")
        ax_map.set_ylabel("y (mm)")
        ax_map.set_xticks((-10, 0, 10))
        ax_map.set_yticks((-10, 0, 10))
        ax_map.text(.03, .03, f"Patch RMS  {case['patch_rms']:.3f}", transform=ax_map.transAxes,
                    fontsize=9.5, color=INK, bbox=dict(fc="white", ec="none", pad=1, alpha=.9))
        risk48, risk72 = risk_rows(case, transfer)
        error_panel(ax_48, risk48, 48, show_title=i == 0, show_labels=True)
        error_panel(ax_72, risk72, 72, show_title=i == 0, show_labels=True)
    color_axis = fig.add_axes((.076, .076, .145, .014))
    colorbar = fig.colorbar(displayed, cax=color_axis, orientation="horizontal",
                            ticks=(-.8, -.4, 0., .4, .8), label="Permittivity difference  q₊ − q₋")
    colorbar.outline.set_visible(False)
    fig.text(.35, .073, "2 F uses 2 and 30 GHz. 4 F uses 2, 12, 20 and 30 GHz.",
             fontsize=10.5, color=INK, va="center")
    fig.text(.35, .042, f"Known means, circular Gaussian noise and equal priors · {NOISE_NOTE}",
             fontsize=10.5, color=INK, va="center")
    save(fig, "witness_trade")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure", choices=("dof", "witness", "both"), default="both")
    args = parser.parse_args()
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11.5, "text.color": INK,
        "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "pdf.fonttype": 42,
    })
    if args.figure in ("dof", "both"):
        dof_figure()
    if args.figure in ("witness", "both"):
        witness_figure()


if __name__ == "__main__":
    main()
