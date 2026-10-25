#!/usr/bin/env python3
"""Static scientific figures for W-SLAM.

Figure 1 is a model diagram. Figure 2 contains finite-dimensional analytic
examples, not Maxwell PDE experiments or computational benchmarks.
"""
from pathlib import Path
import os
import tempfile
import json
import numpy as np
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "wslam-matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "mathtext.fontset": "dejavusans",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.facecolor": "white",
})

BLUE = "#245B78"
LIGHT_BLUE = "#E6F0F5"
ORANGE = "#C97927"
LIGHT_ORANGE = "#FCF0E1"
GRAY = "#53616B"
LIGHT_GRAY = "#EFF2F3"
RED = "#B04747"


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"{name}.{ext}", dpi=230, bbox_inches="tight")
    plt.close(fig)


def factor_graph():
    fig, ax = plt.subplots(figsize=(11.2, 5.8))
    ax.set_xlim(-0.3, 11.4)
    ax.set_ylim(-0.95, 6.0)
    ax.axis("off")

    def edge(a, b, color=GRAY, width=1.35):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=width, zorder=1)

    def variable(pos, label, fill=LIGHT_BLUE, radius=0.34, edgecolor=BLUE):
        ax.add_patch(Circle(pos, radius, facecolor=fill, edgecolor=edgecolor,
                            lw=1.7, zorder=3))
        ax.text(*pos, label, ha="center", va="center", fontsize=14, color=edgecolor, zorder=4)

    def factor(pos, label, width=0.52, color=BLUE, face=LIGHT_BLUE):
        ax.add_patch(Rectangle((pos[0] - width / 2, pos[1] - width / 2),
                               width, width, facecolor=face, edgecolor=color, lw=1.6, zorder=3))
        ax.text(*pos, label, ha="center", va="center", fontsize=12, color=color, zorder=4)

    w = (5.55, 5.1)
    wall_prior = (3.55, 5.1)
    edge(wall_prior, w)
    factor(wall_prior, r"$p_W$", width=0.56)
    variable(w, r"$W$", radius=0.41)
    ax.text(6.18, 5.1, "one shared static wall", ha="left", va="center",
            fontsize=11, color=BLUE)

    sx = [1.65, 5.55, 9.45]
    state_y, maxwell_y = 1.65, 3.48
    for i, x in enumerate(sx, start=1):
        f = (x, maxwell_y)
        state = (x, state_y)
        obs = (x - 1.1, maxwell_y)
        edge(w, f, width=1.6)
        edge(f, state)
        edge(obs, f)
        factor(f, rf"$f_{i}$", width=0.58, color=ORANGE, face=LIGHT_ORANGE)
        variable(state, rf"$s_{i}$")
        variable(obs, rf"$y_{i}$", fill=LIGHT_GRAY, radius=0.27, edgecolor=GRAY)
        ax.text(x + 0.44, maxwell_y - 0.15, "Maxwell\nlikelihood",
                ha="left", va="center", fontsize=9.4, color=ORANGE)
        ax.text(x, state_y - 0.64, rf"$s_{i}=(X_{i},q_{i},\gamma_{i})$",
                ha="center", va="center", fontsize=11, color=BLUE)

    initial = (0.55, state_y)
    edge(initial, (sx[0], state_y))
    factor(initial, r"$p_1$", width=0.48)
    for i in range(2):
        pos = ((sx[i] + sx[i + 1]) / 2, state_y)
        edge((sx[i], state_y), pos)
        edge(pos, (sx[i + 1], state_y))
        factor(pos, rf"$\psi_{{{i+1},{i+2}}}$", width=0.56)
        ax.text(pos[0], state_y - 0.50, "state dynamics",
                ha="center", va="center", fontsize=9.2, color=GRAY)

    ax.text(0.0, 5.78, "A shared wall couples every observation", fontsize=14, color="#202A30")
    ax.text(0.0, 0.30,
            r"$f_t=p(y_t\mid W,X_t,q_t,\gamma_t)$"
            "     •     circles: variables     squares: factors     shaded gray: observed data",
            fontsize=10.3, color=GRAY)
    ax.text(0.0, -0.05,
            "Markov state dynamics are shown. A non-Markov calibration prior requires additional temporal factors.\n"
            "Each likelihood appears once; eliminating electromagnetic fields generally yields global map sensitivity.",
            fontsize=9.6, color=GRAY, linespacing=1.55, va="top")
    save(fig, "figure_01_maxwell_factor_graph")


def algebra_examples():
    fig, axes = plt.subplots(1, 3, figsize=(12.1, 4.3), gridspec_kw={"width_ratios": [0.9, 1.0, 1.3]})
    fig.subplots_adjust(left=0.055, right=0.985, top=0.80, bottom=0.23, wspace=0.40)

    ax = axes[0]
    xs = np.arange(3)
    vals = [0.0, 0.0, 2.0]
    ax.bar(xs, vals, width=0.58, color=[LIGHT_BLUE, LIGHT_BLUE, BLUE],
           edgecolor=BLUE, linewidth=1.2)
    ax.scatter(xs[:2], [0.025, 0.025], s=36, facecolors="white", edgecolors=BLUE, zorder=4)
    for x, y in zip(xs, vals):
        ax.text(x, y + 0.10, f"{y:g}", ha="center", color=BLUE, fontsize=12)
    ax.set_ylim(-0.02, 2.5)
    ax.set_xticks(xs, [r"$y_1$", r"$y_2$", r"$y_{1:2}$"])
    ax.set_ylabel("Effective wall information")
    ax.set_title("(a) Shared-state recovery", loc="left", pad=24, fontsize=11.4)
    ax.text(0.5, 1.03, r"$y_1=w+x,\quad y_2=w-x$", transform=ax.transAxes,
            ha="center", fontsize=10.2, color=GRAY)
    ax.yaxis.grid(True, color="#E5EAED", linewidth=0.6)
    ax.set_axisbelow(True)

    ax = axes[1]
    vals = [4/3, 0.0, 0.0]
    bars = ax.bar(np.arange(3), vals, width=0.58, color=[LIGHT_ORANGE, LIGHT_BLUE, LIGHT_BLUE],
                  edgecolor=[RED, BLUE, BLUE], linewidth=1.2)
    bars[0].set_hatch("///")
    ax.scatter([1, 2], [0.018, 0.018], s=36, facecolors="white", edgecolors=BLUE, zorder=4)
    ax.text(0, 4/3 + 0.08, r"$4/3$", ha="center", color=RED, fontsize=12)
    ax.text(1, 0.10, "0", ha="center", color=BLUE, fontsize=12)
    ax.text(2, 0.10, "0", ha="center", color=BLUE, fontsize=12)
    ax.set_xticks([0, 1, 2], ["Naive\n" + r"$G^\dagger$", "Joint\nSchur", "Correct\nquotient"])
    ax.set_ylim(-0.02, 1.68)
    ax.set_ylabel("Effective wall information")
    ax.set_title("(b) False recovery from a pseudoinverse", loc="left", pad=24, fontsize=10.4)
    ax.text(0.5, 1.03, "Unobserved target motion can compensate", transform=ax.transAxes,
            ha="center", fontsize=9.5, color=GRAY)
    ax.yaxis.grid(True, color="#E5EAED", linewidth=0.6)
    ax.set_axisbelow(True)

    ax = axes[2]
    n = np.logspace(0, 4, 500)
    sigma, bias, reopen_at = 0.2, 0.1, 20
    full = 5 * sigma**2 / n
    floor = 9 * bias**2 / 4
    frozen = floor + sigma**2 / (2*n)
    reopened = np.where(n < reopen_at, frozen, full)
    ax.loglog(n, frozen, color=ORANGE, lw=2.1, label="Permanent wrong wall")
    ax.loglog(n, full, color=BLUE, lw=2.1, label="Joint estimation")
    ax.loglog(n, reopened, color="#527B57", lw=1.5, ls=(0, (3, 2)), label="Reopen at n = 20")
    ax.axhline(floor, color=ORANGE, lw=0.9, ls=":")
    ax.text(10500, floor, r"$9b^2/4$", ha="left", va="center", fontsize=9.8, color=ORANGE)
    ax.axvline(reopen_at, color="#84918A", lw=0.8, ls=":", zorder=0)
    ax.set_xlim(1, 1.8e4)
    ax.set_ylim(1.2e-5, 0.32)
    ax.set_xlabel("Independent repeats per design n")
    ax.set_ylabel("Expected squared target error")
    ax.set_title("(c) Irreducible freezing bias", loc="left", pad=24, fontsize=11.4)
    ax.text(0.5, 1.03, r"$y^{(1)}=w+x,\quad y^{(2)}=2w+x$", transform=ax.transAxes,
            ha="center", fontsize=10.2, color=GRAY)
    ax.grid(True, which="major", color="#E5EAED", linewidth=0.6)
    ax.legend(frameon=False, loc="lower left", fontsize=8.3, handlelength=2.8)

    fig.text(0.055, 0.078,
             "Finite-dimensional analytic examples only. No Maxwell PDE solver, imaging accuracy test, or speed benchmark.",
             fontsize=10, color=GRAY)
    fig.text(0.055, 0.028,
             "Panel (c): independent Gaussian noise, σ = 0.2; fixed wall error b = 0.1; all previous data retained when reopening.",
             fontsize=9.5, color=GRAY)
    save(fig, "figure_02_linear_algebra_counterexamples")

    metadata = {
        "scope": "Analytic linear-Gaussian examples, not Maxwell experiments.",
        "shared_state_information": [0, 0, 2],
        "pseudoinverse_false_information": 4/3,
        "correct_information": 0,
        "bias_example": {
            "noise_standard_deviation": sigma,
            "fixed_wall_error": bias,
            "reopen_at_repeats_per_design": reopen_at,
            "full_mse": "5*sigma**2/n",
            "frozen_mse": "9*bias**2/4 + sigma**2/(2*n)",
            "asymptotic_squared_bias": floor
        }
    }
    (FIGURES / "figure_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    factor_graph()
    algebra_examples()
    print(FIGURES)
