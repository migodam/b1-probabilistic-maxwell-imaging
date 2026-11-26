"""Two post-seal examples of removable-by-data but material-positive regions."""
from __future__ import annotations

import json
from shared import HERE, OUT, STUDY

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CASES = ("n32_c02", "n32_c04")


def plot():
    seal = json.loads((OUT/"sealed_32.json").read_text())
    if not seal["sealed"] or not all(c in seal["stationary"] for c in CASES):
        raise RuntimeError("Selected plot labels require both sealed objects")
    fig, axes = plt.subplots(2, 4, figsize=(10.7,6.0), layout="constrained")
    heads = ("Known material", "Stable CIE-I", "Zero-ROI full-wave witness", "Equal-wall full-wave fit")
    for j, title in enumerate(heads):
        axes[0,j].set_title(title, fontsize=10, fontweight="semibold")
    scores = json.loads((OUT/"material_scores32.json").read_text())["objects"]
    for i, case in enumerate(CASES):
        path = OUT/case/"stationary"
        b1 = json.loads((path/"b1.json").read_text())
        selected = int(b1["selected_roi"].rsplit("_",1)[-1])-1
        with np.load(STUDY/"ports/private"/f"{case}.npz", allow_pickle=False) as d:
            truth = d["material"]
        with np.load(path/"baseline.npz", allow_pickle=False) as z:
            initial = z["chi"]
            roi = z["rois"][selected]
        with np.load(path/"b1_predictions.npz", allow_pickle=False) as z:
            cleaned = z["chi"]
        with np.load(path/"b1_controls.npz", allow_pickle=False) as z:
            control = z["physical_chi"]
        values = (truth, initial, cleaned, control)
        keys = (None, "baseline", "b1", "b1_wall_control")
        rows = scores[case]["stationary"]["method_metrics"]
        count = int(np.count_nonzero(truth[roi] > .25))
        for j, (z, key) in enumerate(zip(values, keys)):
            ax = axes[i,j]
            layer = ax.imshow(z, origin="lower", vmin=0., vmax=4.5, cmap="viridis",
                              interpolation="none")
            ax.contour(roi, levels=[.5], colors=["#ff794f"], linewidths=1.25, origin="lower",
                       extent=(-.5,31.5,-.5,31.5))
            ax.set_xticks([])
            ax.set_yticks([])
            if key is None:
                subtitle = f"{case}\n{count}/{int(roi.sum())} ROI cells contain material"
            else:
                m = rows[key]
                subtitle = (f"training {m['residual_over_radius']:.3g} r"
                            f"  ·  held {100*m['held_relative_residual']:.1f}%\n"
                            f"material relative L2 {m['material_relative_L2']:.3f}")
            ax.set_xlabel(subtitle,fontsize=8.3,labelpad=5.8)
    c = fig.colorbar(layer, ax=axes, location="right", shrink=.8, pad=.01)
    c.set_label("Real contrast")
    fig.savefig(HERE/"figures"/"sealed_zero_roi_cases.png",dpi=180)
    fig.savefig(HERE/"figures"/"sealed_zero_roi_cases.pdf")


def plot_native_case():
    case = "n96_c01"
    seal = json.loads((OUT/"sealed_96.json").read_text())
    if not seal["sealed"] or case not in seal["cases"]:
        raise RuntimeError("Selected plot labels require the sealed native object")
    path = OUT/case
    b1 = json.loads((path/"b1.json").read_text())
    selected = int(b1["selected_roi"].rsplit("_",1)[-1])-1
    with np.load(STUDY/"ports/private"/f"{case}.npz", allow_pickle=False) as d:
        truth = d["material"]
    with np.load(path/"baseline.npz", allow_pickle=False) as z:
        initial, roi = z["chi"], z["rois"][selected]
    with np.load(path/"b1_predictions.npz", allow_pickle=False) as z:
        clean = z["chi"]
    with np.load(path/"b1_controls.npz", allow_pickle=False) as z:
        control = z["physical_chi"]
    metrics = json.loads((OUT/"material_scores96.json").read_text())["objects"][case]["dual"]["method_metrics"]
    fig, axes = plt.subplots(1,4,figsize=(10.65,3.5),layout="constrained")
    titles = ("Known material", "Dual-CIE coefficient stop", "Zero-ROI full-wave witness", "Equal-wall full-wave fit")
    for ax, title, z, k in zip(axes,titles,(truth,initial,clean,control),(None,"baseline","b1","b1_wall_control")):
        layer = ax.imshow(z, origin="lower", vmin=0., vmax=4.5, cmap="viridis", interpolation="none")
        ax.contour(roi, levels=[.5], colors=["#ff794f"], linewidths=1.25, origin="lower",
                   extent=(-.5,95.5,-.5,95.5))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title,fontsize=10,fontweight="semibold")
        if k is None:
            caption = f"{case}\nAll {int(roi.sum())} ROI cells are vacuum"
        else:
            m = metrics[k]
            caption = (f"training {m['residual_over_radius']:.3g} r  ·  held "
                       f"{100*m['held_relative_residual']:.1f}%\nmaterial relative L2 "
                       f"{m['material_relative_L2']:.3f}")
        ax.set_xlabel(caption,fontsize=8.4,labelpad=5.5)
    c = fig.colorbar(layer,ax=axes,location="right",shrink=.72,pad=.01)
    c.set_label("Real contrast")
    fig.savefig(HERE/"figures"/"sealed_native_zero_roi.png",dpi=180)
    fig.savefig(HERE/"figures"/"sealed_native_zero_roi.pdf")


if __name__ == "__main__":
    plot()
    if (OUT/"sealed_96.json").is_file():
        plot_native_case()
