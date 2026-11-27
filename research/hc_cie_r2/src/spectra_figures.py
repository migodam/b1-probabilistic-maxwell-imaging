"""Draw the published spectral and finite-perturbation figures from paid arrays."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


VIEWS = {
    "r13_r13_fresh_0001": ("R13 full ring 24²", "#365575", "-"),
    "r13_r13_fresh_0000": ("R13 half ring 24²", "#365575", "--"),
    "srap_s2d_010__pw360_10x20": ("PW 360° 32²", "#137b80", "-"),
    "srap_s2d_010__pw180_10x20": ("PW 180° 32²", "#7565a9", "-"),
    "srap_s2d_010__pw120_10x20": ("PW 120° 32²", "#d4842f", "-"),
    "srap_s2d_010__ls_arc_10x20": ("Line sources to arc 32²", "#429d76", "--"),
    "srap_s2d_010__ls_line_same_10x20": ("Line receive, same side 32²", "#b44264", "-"),
    "srap_s2d_010__ls_line_opposite_10x20": ("Line receive, opposite side 32²", "#b44264", "--"),
}
PARTIAL = ("srap_s2d_010__pw360_10x20", "srap_s2d_010__pw180_10x20",
           "srap_s2d_010__pw120_10x20", "srap_s2d_010__ls_line_same_10x20")


def _styling():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.labelcolor": "#182b3a", "axes.edgecolor": "#afb9bf",
                         "axes.titlecolor": "#182b3a", "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titlesize": 11,
                         "xtick.color": "#38505e", "ytick.color": "#38505e",
                         "grid.alpha": 0.23, "figure.facecolor": "#ffffff"})


def draw_spectra(inputs, figures):
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 7.8), layout="constrained")
    gamma, recv, l_map, rank = axes.flat
    for key, (label, color, style) in VIEWS.items():
        with np.load(inputs / f"{key}.npz", allow_pickle=False) as saved:
            total = 12 if key.startswith("r13_") else 10
            sigma_g = saved["roundoff_cut_gamma_singular_values"][:-1]
            sigma_l = saved[f"roundoff_cut_{total}_sources_reciprocity_l_singular_values"]
            sigma_s = saved["receive_singular_values"]
            gamma.semilogy(np.arange(1, len(sigma_g) + 1) / len(sigma_g), sigma_g / sigma_g[0],
                           color=color, ls=style, lw=1.75, label=label)
            l_map.semilogy(np.arange(1, len(sigma_l) + 1), sigma_l,
                          color=color, ls=style, lw=1.75)
            recv.semilogy(np.arange(1, len(sigma_s) + 1), sigma_s / sigma_s[0],
                          color=color, ls=style, lw=1.75)
    rows = json.loads((inputs / "spectrum.json").read_text())
    by_key = {r["observation"]: [] for r in rows}
    for row in rows:
        if row["mode"] == "roundoff_cut":
            by_key[row["observation"]].append(row)
    representative = ("r13_r13_fresh_0000", "r13_r13_fresh_0001", "srap_s2d_010__pw360_10x20")
    for key in representative:
        _, color, style = VIEWS[key]
        cases = by_key[key]
        rank.plot([r["source_count"] for r in cases], [100 * r["l_over_columns"] for r in cases],
                  color=color, ls=style, marker="o", ms=4)
    for threshold in (1e-2, 1e-4, 1e-6):
        recv.axhline(threshold, color="#8c969b", linestyle=":", linewidth=0.85)
    gamma.set(title="Quadratic curvature Γ on its positive spectrum", xlabel="Fractional mode index",
              ylabel="Singular value / maximum", ylim=(1e-4, 1.3))
    recv.set(title="Receive operator and the three relative cuts", xlabel="Receive mode index",
             ylabel="Singular value / maximum")
    l_map.set(title="Reciprocity L on its reduced row space", xlabel="Reduced mode index",
              ylabel="Singular value in native units")
    rank.set(title="Complex current coordinates constrained by L", xlabel="Number of retained illuminations",
             ylabel="Numerical rank / columns (%)", xticks=[2, 4, 6, 8, 10, 12], ylim=(0, 1.2))
    for axis in axes.flat:
        axis.grid(True, which="major")
    handles, labels = gamma.get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=4, frameon=False, fontsize=9)
    fig.savefig(figures / "spectra_operators.png", dpi=165)
    plt.close(fig)


def draw_partial(inputs, figures):
    rows = json.loads((inputs / "partial.json").read_text())
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 7.8), layout="constrained")
    physical, remainder, truncated, readout = axes.flat
    for key in PARTIAL:
        label, color, style = VIEWS[key]
        group = [row for row in rows if row["observation"] == key]
        scale = [row["scale"] for row in group]
        physical.plot(scale, [r["physical_port_over_delta"] for r in group],
                      label=label, color=color, ls=style, marker="o", ms=3.6)
        remainder.semilogy(scale, [r["full_quadratic_norm"] / r["full_linear_norm"] for r in group],
                           color=color, ls=style, marker="o", ms=3.6)
        truncated.plot(scale, [100 * r["partial"]["0.5"]["full_curvature_fraction"] for r in group],
                       color=color, ls=style, marker="o", ms=3.6)
        readout.semilogy(scale, [r["alternative_physical_readout_port_over_delta"] for r in group],
                        color=color, ls=style, marker="o", ms=3.6)
        readout.semilogy(scale, [r["alternative_full_port_over_delta"] for r in group],
                        color=color, ls=":", alpha=0.75)
    comparison = json.loads((inputs / "feasible_followup.json").read_text())
    group = comparison["contrast"]
    scale = [row["scale"] for row in group]
    refined = dict(color="#152e38", ls="--", marker="D", ms=3, zorder=6)
    physical.plot(scale, [r["physical_port_over_delta"] for r in group],
                  label="PW 360° refined", **refined)
    physical.scatter([1.0], [min(r["physical_port_over_delta"] for r in group)],
                     color="#152e38", marker="D", s=25, zorder=7)
    remainder.semilogy(scale, [r["full_quadratic_norm"] / r["full_linear_norm"] for r in group], **refined)
    truncated.plot(scale, [100 * r["partial"]["0.5"]["full_curvature_fraction"] for r in group], **refined)
    readout.semilogy(scale, [r["alternative_physical_readout_port_over_delta"] for r in group], **refined)
    physical.axhline(1, color="#283c4a", linestyle=":", linewidth=1)
    readout.axhline(1, color="#283c4a", linestyle=":", linewidth=1)
    physical.set(title="Original physical material after amplitude changes",
                 ylabel="Physical port residual / noise radius")
    remainder.set(title="Nonlinear term after exact affine elimination", ylabel="Quadratic / linear residual norm")
    truncated.set(title="Curvature retained by the weakest half of Γ modes",
                  ylabel="Projected / full quadratic norm (%)", ylim=(0, 100))
    readout.set(title="Eliminated current (dotted) and box readout after forward solve",
                ylabel="Port residual / noise radius")
    for axis in axes.flat:
        axis.set_xlabel("Material contrast scale")
        axis.set_xticks([0.5, 0.8, 1.0, 1.2, 1.5])
        axis.grid(True, which="major")
    handles, labels = physical.get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=5, frameon=False, fontsize=9)
    fig.savefig(figures / "spectra_partial_elimination.png", dpi=165)
    plt.close(fig)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "results/spectra")
    parser.add_argument("--figures", type=Path, default=root / "figures")
    args = parser.parse_args()
    args.figures.mkdir(parents=True, exist_ok=True)
    _styling()
    draw_spectra(args.input, args.figures)
    draw_partial(args.input, args.figures)


if __name__ == "__main__":
    main()
