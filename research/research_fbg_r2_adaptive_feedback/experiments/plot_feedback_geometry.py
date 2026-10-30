"""Plot existing development geometry diagnostics; never run a forward solve.

Only the post-box feedback matrices in schur_geometry.json are shown. The
lower panels use one exact public direction at the common 0.1 anchor; they
describe local geometry, not an optimizer basin or fresh-data evidence.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / "figures/feedback_geometry.png"
FEEDBACK = ("zero", "beta_.2", "self", "diagonal", "Schur_diagonal", "block")
LABEL = {"zero": "B = 0", "beta_.2": "fixed β = 0.2", "self": "fitted scalar",
         "diagonal": "fitted diagonal", "Schur_diagonal": "Schur diagonal", "block": "fitted block"}
COLORS = {"zero": "#4c4c4c", "beta_.2": "#e17c20", "self": "#549a7c",
          "diagonal": "#6393c5", "Schur_diagonal": "#a07ac4", "block": "#bd6963"}
TERMS = (("psi_second", "f″ / f"),
         ("positive_weighted_J_term", "positive J term / f"),
         ("physical_residual_second", "physical residual / f"),
         ("signed_Gram_cross", "signed Gram cross / f"),
         ("signed_Gram_second", "signed Gram second / f"))


def main() -> None:
    schur = json.loads((HERE / "results/schur_geometry.json").read_text())
    bias = json.loads((HERE / "results/contrast_bias_geometry.json").read_text())
    high8 = {r["feedback"]: r for r in schur["rows"] if r["case_id"] == "dev_high_digit8"}
    if set(FEEDBACK) - set(high8):
        raise ValueError("Missing measured high-digit-8 post-box feedback rows")
    real = {name: np.asarray(high8[name]["B_diagonal_real"], float) for name in FEEDBACK}
    imag = {name: np.asarray(high8[name]["B_diagonal_imag"], float) for name in FEEDBACK}
    if any(a.size != 576 for a in (*real.values(), *imag.values())):
        raise ValueError("Expected actual 24×24 post-box B diagonals")
    rlim = max(np.max(np.abs(a)) for a in real.values())
    ilim = max(np.max(np.abs(a)) for a in imag.values())
    d_cases = ("dev_high_digit8", "dev_high_ring", "dev_high_curved_bar")
    d_labels = ("Digit 8", "Ring", "Curved bar")
    d_rows = {r["case_id"]: np.asarray(r["d"], float) for r in bias["rows"]}
    if any(d_rows[c].size != 576 for c in d_cases):
        raise ValueError("Expected three measured 24×24 residual-column-energy maps")
    dmax = max(float(np.max(d_rows[c])) for c in d_cases)

    fig = plt.figure(figsize=(17, 20), constrained_layout=False)
    gs = fig.add_gridspec(9, 6, height_ratios=[1]*6 + [1.22, 1.0, 1.25],
                          left=.07, right=.95, bottom=.06, top=.945,
                          hspace=.68, wspace=.35)
    fig.suptitle("Development feedback geometry · measured post-box matrices and one exact direction",
                 fontsize=17, y=.98, fontweight="bold")
    first_real = first_imag = None
    for i, name in enumerate(FEEDBACK):
        row = high8[name]
        a = real[name]
        b = imag[name]
        axr = fig.add_subplot(gs[i, 0:2])
        axi = fig.add_subplot(gs[i, 2:4])
        axs = fig.add_subplot(gs[i, 4:6])
        first_real = axr.imshow(a.reshape(24, 24), cmap="RdBu_r",
                                norm=TwoSlopeNorm(vcenter=0, vmin=-rlim, vmax=rlim), origin="lower")
        # All measured imaginary diagonals may be zero in a degenerate receipt.
        if ilim > 0:
            first_imag = axi.imshow(b.reshape(24, 24), cmap="RdBu_r",
                                    norm=TwoSlopeNorm(vcenter=0, vmin=-ilim, vmax=ilim), origin="lower")
        else:
            axi.imshow(b.reshape(24, 24), cmap="RdBu_r", vmin=-1, vmax=1, origin="lower")
        mag = np.sort(np.hypot(a, b))[::-1]
        nonzero = mag > 0
        if nonzero.any():
            axs.semilogy(np.flatnonzero(nonzero)+1, mag[nonzero], color=COLORS[name], lw=2)
            axs.set_xlim(1, 576)
            axs.set_ylim(1e-5, max(.3, mag.max()*1.25))
        else:
            wording = "rejected to zero" if row["safety"].get("rejected_to_zero") else "all zero: identity baseline"
            axs.text(.5, .5, wording, ha="center", va="center", transform=axs.transAxes,
                     fontsize=10, color="#9c3434" if "rejected" in wording else "#555555")
            axs.set_xlim(1, 576)
            axs.set_ylim(1e-5, .3)
            axs.set_yscale("log")
        axr.set_ylabel(LABEL[name], fontsize=10, fontweight="bold", labelpad=12)
        for ax in (axr, axi):
            ax.set_xticks([])
            ax.set_yticks([])
        axs.grid(alpha=.2, which="both")
        axs.set_xlabel("sorted cell index" if i == len(FEEDBACK)-1 else "", fontsize=8)
        axs.tick_params(labelsize=8)
        if i == 0:
            axr.set_title(f"Re diag(B), common ±{rlim:.3g}", fontsize=11)
            axi.set_title(f"Im diag(B), common ±{ilim:.3g}", fontsize=11)
            axs.set_title("|diag(B)|, post-box spectrum", fontsize=11)
    fig.text(.50, .355, "Only diagonal entries are displayed; fitted block feedback can also have off-diagonal entries.",
             ha="center", fontsize=9, color="#555555")

    d_axes = []
    for j, (case, label) in enumerate(zip(d_cases, d_labels)):
        ax = fig.add_subplot(gs[6, 2*j:2*j+2])
        im = ax.imshow(d_rows[case].reshape(24, 24), cmap="magma", vmin=0, vmax=dmax,
                       origin="lower")
        ax.set_title(f"d: {label} (A = 12)", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        d_axes.append(ax)
    cb = fig.colorbar(im, ax=d_axes, location="right", shrink=.8, pad=.015)
    cb.set_label("residual-column energy, shared raw units", fontsize=9)

    xs = np.arange(len(FEEDBACK))
    axj = fig.add_subplot(gs[7, :])
    jvp = [high8[name]["transformed_projector_jvp_fro"] for name in FEEDBACK]
    axj.bar(xs, jvp, color=[COLORS[name] for name in FEEDBACK], width=.65)
    axj.set_xticks(xs, [LABEL[name] for name in FEEDBACK], rotation=12, fontsize=9)
    axj.set_ylabel("projector JVP Frobenius norm", fontsize=9)
    axj.set_title("Digit 8: exact directional projector rotation", fontsize=11)
    axj.set_ylim(0, max(jvp)*1.22)
    axj.grid(axis="y", alpha=.2)
    axj.set_axisbelow(True)

    axc = fig.add_subplot(gs[8, :])
    width = .15
    term_colors = ("#303030", "#3d8c67", "#b75b51", "#6874b6", "#b28d4f")
    for j, (key, label) in enumerate(TERMS):
        vals = [high8[name][key]/high8[name]["psi"] for name in FEEDBACK]
        axc.bar(xs+(j-2)*width, vals, width=width, label=label, color=term_colors[j])
    axc.axhline(0, color="black", lw=.8)
    axc.set_xticks(xs, [LABEL[name] for name in FEEDBACK], rotation=12, fontsize=9)
    axc.set_ylabel("signed term / raw ψ", fontsize=9)
    axc.set_title("Digit 8: directional curvature terms normalized by raw objective ψ", fontsize=11)
    axc.legend(ncol=5, loc="upper center", bbox_to_anchor=(.5, -.27), fontsize=8, frameon=False)
    axc.grid(axis="y", alpha=.2)
    axc.set_axisbelow(True)
    fig.text(.5, .016, "β = 0.2 improves recovery in Stage 1 despite larger measured projector rotation than B = 0. These local terms are not a basin certificate.",
             ha="center", fontsize=9, color="#333333")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=170)
    plt.close(fig)
    print(OUT)


if __name__ == "__main__":
    main()
