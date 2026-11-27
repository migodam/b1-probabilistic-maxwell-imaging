"""Display change localization on original three-axis reference slices."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/vector_change"
N = 12


def draw():
    receipt = json.loads((RESULTS / "report.json").read_text())
    parent = receipt["config"]["parents"][0]
    with np.load(RESULTS / f"hemi_{parent:02d}.npz", allow_pickle=False) as data:
        ref = data["background"].reshape((N,) * 3)
        actual = data["true_sites"], data["true_change"]
    center = np.unravel_index(actual[0][np.argmax(actual[1])], ref.shape)
    views = [("Reference change", *actual)]
    for array in ("hemi", "plane_grid"):
        with np.load(RESULTS / f"{array}_{parent:02d}.npz", allow_pickle=False) as data:
            for cond in ("diff10", "full10"):
                label = {"hemi": "Hemi", "plane_grid": "Plane grid"}[array]
                label += {"diff10": " · 10% difference", "full10": " · 10% full field"}[cond]
                views.append((label, data[f"{cond}_unknown_sites"], data[f"{cond}_unknown_change"]))
    with np.load(RESULTS / "two_sided.npz", allow_pickle=False) as two_sided:
        prefix = f"plane_grid_{parent:02d}_clean"
        views.insert(3, ("Plane grid · clean, two-sided", two_sided[prefix + "_proposed"],
                         two_sided[prefix + "_proposed_change"]))
    fig, axes = plt.subplots(len(views), 3, figsize=(8.5, 11.6), layout="constrained")
    cmap = plt.get_cmap("coolwarm")
    norm = TwoSlopeNorm(vmin=-1.35, vcenter=0.0, vmax=1.35)
    displayed = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    for row, (label, sites, amplitudes) in enumerate(views):
        for dim, named in enumerate(("YZ", "XZ", "XY")):
            ax = axes[row, dim]
            sl = np.take(ref, center[dim], axis=dim).T
            ax.imshow(sl, cmap="Greys", origin="lower", vmin=0.0, vmax=3.0, interpolation="nearest", alpha=0.63)
            if row == 0:
                ax.set_title(f"{named} at index {center[dim]}", fontsize=10)
            if dim == 0:
                ax.set_ylabel(label, fontsize=8.6)
            for index, value in zip(sites, amplitudes, strict=True):
                if abs(value) <= 1e-6:
                    continue
                xyz = np.unravel_index(index, ref.shape)
                first, second = np.delete(xyz, dim)
                color = cmap(norm(np.clip(value, -1.35, 1.35)))
                off = xyz[dim] - center[dim]
                ax.scatter(first, second, s=66, marker="o" if off else "s",
                           facecolors="none" if off else color, edgecolors=color if off else "#193845",
                           linewidths=1.4, zorder=4)
                adj = f"{off:+d} off plane" if off else f"{value:+.2f}"
                shift = -5 if first >= N - 3 else 5
                vertical = -9 if second >= N - 2 else 5
                ax.annotate(adj, (first, second), xytext=(shift, vertical), textcoords="offset points",
                            ha="right" if shift < 0 else "left", va="top" if vertical < 0 else "bottom",
                            fontsize=6, color="#193845",
                            bbox=dict(fc="white", ec="none", alpha=0.67, pad=1))
            ax.set(xticks=[], yticks=[], xlim=(-0.5, N - 0.5), ylim=(-0.5, N - 0.5))
    fig.colorbar(displayed, ax=axes, pad=0.02, shrink=0.75, label="Contrast change Δχ (clipped at ±1.35 in plot)")
    fig.savefig(ROOT / "figures/vector_change_slices.png", dpi=175)
    plt.close(fig)


if __name__ == "__main__":
    draw()
