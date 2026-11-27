"""Plot the archived strict segments and material maps on common physical scales."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from pathlib import Path
HERE = Path(__file__).resolve().parent

CASES = (("n32_c01", 32), ("n96_c01", 96))
COLORS = {"broad": "#075e88", "fixed": "#678694", "shape": "#b54a2b", "whole": "#986634"}


def read(case, suffix):
    path = HERE / f"{case}_paths{suffix}.json"
    if path.exists():
        return json.loads(path.read_text()), np.load(path.with_suffix(".npz"), allow_pickle=False)
    return None


def main():
    fig = plt.figure(figsize=(12, 5.5), layout="constrained")
    grid = fig.add_gridspec(2, 5, width_ratios=(1., 1., 2.35, 1.5, .12))
    for row, (case, side) in enumerate(CASES):
        files = {key: read(case, suff) for key, suff in [("broad", "_broad"),
                 ("fixed", "_ball"), ("shape", "_ball"), ("whole", "_whole_region")]}
        show, raw = files["broad"] or files["fixed"]
        roi = raw["roi_mask"].reshape(side, side)
        patch = np.zeros(side ** 2, bool)
        patch[raw["patch"]] = True
        patch = patch.reshape(side, side)
        material = raw["ball_min_margin_materials"]
        for col, (x, title) in enumerate(((material[0], "Vacuum region"), (material[-1], "Fixed patch endpoint"))):
            ax = fig.add_subplot(grid[row, col])
            im = ax.imshow(x.reshape(side, side), cmap="cividis", vmin=0, vmax=4.5)
            ax.contour(roi, levels=[.5], colors=["#f7f3df"], linewidths=.75)
            ax.contour(patch, levels=[.5], colors=["#ed8657"], linewidths=.9)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(title if row == 0 else "")
            if col == 0:
                ax.set_ylabel(rf"${side}\times{side}$")
        ax = fig.add_subplot(grid[row, 2])
        ax_cost = fig.add_subplot(grid[row, 3], sharex=ax)
        ax.axhline(1., color="#9f2318", lw=.9, ls="--")
        for key, source in files.items():
            if source is None:
                continue
            j, data = source
            method = "ball_shape_free" if key in ("shape", "whole") else "ball_min_margin"
            materials = data[method + "_materials"]
            ext = ~data["roi_mask"]
            region = materials[:, ~ext].mean(axis=1)
            ext_l2 = np.linalg.norm(materials[:, ext] - materials[[0], :][:, ext], axis=1) * j["cell_area_m2"]**.5
            good = [r["verification"]["residual_upper_over_r"] for r in j["paths"][method]["attempts"] if r["strict"]]
            names = {"broad": "Fixed patch, expanded exterior", "fixed": "Fixed patch, 256 exterior cells",
                     "shape": "Free shape in patch", "whole": "Free shape in entire region"}
            ax.stairs(good, region, color=COLORS[key], label=names[key], baseline=None)
            ax_cost.plot(region, ext_l2, color=COLORS[key], lw=1.3)
            data.close()
        ax.set_ylim(0., 1.09)
        ax.set_ylabel("Full segment residual / radius")
        ax_cost.set_ylabel(r"Exterior weighted $L_2$ [m]")
        if row == 0:
            ax.legend(frameon=False, loc="lower right", fontsize=7.2)
        else:
            ax.set_xlabel("Mean contrast in registered region")
            ax_cost.set_xlabel("Mean contrast in registered region")
        for item in (ax, ax_cost):
            item.xaxis.set_major_locator(MaxNLocator(nbins=4))
            item.grid(True, alpha=.22)
    cax = fig.add_subplot(grid[:, 4])
    fig.colorbar(im, cax=cax, label="Contrast")
    fig.savefig(HERE / "certified_paths.png", dpi=180)
    plt.close(fig)


def b_figure():
    p = HERE / "B_independent_trajectory.json"
    h = HERE / "B_independent_held.json"
    if not (p.exists() and h.exists()):
        return
    report = json.loads(p.read_text())
    hold = json.loads(h.read_text())
    fig, axs = plt.subplots(1, 3, figsize=(10.2, 3.65), layout="constrained")
    with np.load(p.with_suffix(".npz"), allow_pickle=False) as raw:
        roi = raw["roi"].reshape(32,32)
        patch = np.zeros(32**2,bool)
        patch[raw["patch"]] = True
        patch = patch.reshape(32,32)
        for i, (name, label) in enumerate((("vacuum","Vacuum region"), ("high_patch","Connected patch ≥ 2"))):
            ax = axs[i]
            im = ax.imshow(raw[name].reshape(32,32), cmap="cividis",vmin=0.,vmax=4.5)
            ax.contour(roi,levels=[.5],colors=["#f7f3df"],linewidths=.85)
            ax.contour(patch,levels=[.5],colors=["#ed8657"],linewidths=1.0)
            ax.set(xticks=[],yticks=[],title=label)
    fig.colorbar(im, ax=axs[:2], label="Contrast", shrink=.74)
    ax = axs[2]
    names = ("vacuum","high_patch")
    a = report["independent_endpoint_audits"]
    ratios = [[a[k]["B"]["residual_upper_over_r"],
               a[k]["A"]["residual_upper_over_r"],
               hold["scores"][k]["residual_upper_over_r"]] for k in names]
    locations=np.arange(3)
    ax.plot(locations,ratios[0], linestyle="none",marker="o",color=COLORS["broad"],label="Vacuum region")
    ax.plot(locations,ratios[1], linestyle="none",marker="s",color=COLORS["shape"],label="Connected patch ≥ 2")
    ax.axhline(1.,color="#9f2318",lw=.8,ls="--")
    ax.set_yscale("log")
    ax.set(xticks=locations,xticklabels=["Original B","Original A","Held receive arc"],
           ylabel="Certified residual upper / radius",ylim=(.38,46))
    ax.legend(frameon=False,loc="upper left",fontsize=8)
    ax.grid(True,axis="y",alpha=.22)
    fig.savefig(HERE/"B_refit_shapes.png",dpi=180)
    plt.close(fig)


def phase_figure():
    record = HERE/"phase_bands.json"
    if not record.exists():
        return
    log = json.loads(record.read_text())
    known = log["proved_hideable_segments"]
    specs = [
        ("32² A — connected patch","n32", "uniform connected patch", "_broad.json"),
        ("96² A — connected patch","n96", "uniform connected patch", "_broad.json"),
        ("32² B — connected patch","B /", "uniform connected patch", "B_independent_trajectory"),
        ("32² A — whole ROI mean","n32", "whole ROI shape-free mean", "_whole_region.json"),
    ]
    all_pair_costs = [q["sufficient_exterior_budget"]["L2_m"]
                      for q in log["independently_refitted_endpoint_pairs"] if q["sufficient_exterior_budget"]]
    strip_costs = [q["sufficient_exterior_budget"]["L2_m"] for q in known]
    ceiling = max(3.22, 1.09*max([0.]+all_pair_costs+strip_costs))
    fig, axes = plt.subplots(1,4,figsize=(12.3,3.4),sharey=True,layout="constrained")
    for ax, (title, layout, task, source) in zip(axes,specs):
        valid = [q for q in known if layout in q["layout"] and q["task"]==task]
        strips = [q for q in valid if source in q["source"]] or [
            q for q in valid if "_ball.json / ball_min_margin" in q["source"]]
        ax.axvline(2., lw=.75, ls="--", color="#8f3126")
        for row in strips:
            t0,t1 = row["task_t_from"], row["task_t_to"]
            beta = row["sufficient_exterior_budget"]["L2_m"]
            ax.fill_between([t0,t1], [beta,beta], [ceiling,ceiling], color=COLORS["broad"], alpha=.19,lw=0)
            ax.hlines(beta,t0,t1,lw=1.4,color=COLORS["broad"])
        if task=="uniform connected patch":
            pair = next((q for q in log["independently_refitted_endpoint_pairs"]
                if (q["layout"].startswith("B /") if layout=="B /" else
                    q["layout"].startswith("A /") and layout in q["layout"])), None)
            if pair:
                ax.plot(2.,pair["sufficient_exterior_budget"]["L2_m"],marker="D",
                        color=COLORS["shape"],ms=5,label="Pair, patch floor 2")
                ax.legend(frameon=False,loc="lower right",fontsize=6.5)
        ax.set(xlim=(0,2.17),ylim=(0,ceiling),xlabel="Material task amplitude",title=title)
        ax.grid(True,alpha=.2)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    axes[0].set_ylabel(r"Allowed exterior weighted $L_2$ [m]")
    fig.suptitle("Proven feasible bands in the registered noise balls", fontsize=13)
    fig.savefig(HERE/"phase_bands.png",dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
    b_figure()
    phase_figure()
