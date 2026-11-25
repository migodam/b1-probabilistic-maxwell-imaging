"""Render the same-scale material estimates and independent timing measurements."""

import json

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from pool_common import RESULTS
from structured_pool import BRANCHES, CAP, OUT, PRIVATE, PUBLIC, SIDE
from flow_audit import OUTPUT as FLOW_OUT

FIGURES = RESULTS / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)


def payload(p):
    return json.loads(p.read_text())


def write(fig, stem):
    fig.savefig(FIGURES / f"{stem}.png", dpi=180)
    fig.savefig(FIGURES / f"{stem}.pdf")
    plt.close(fig)


def embed(v, active):
    domain = np.zeros(SIDE**2)
    domain[active] = v
    return domain.reshape(SIDE, SIDE)


def structured_plot():
    t = payload(OUT / "timed_pipelines.json")
    scored = payload(OUT / "timed_pipelines_offline.json")
    examined = payload(OUT / "online.json")
    scored_all = payload(OUT / "offline.json")
    with np.load(PUBLIC / "geometry.npz") as data:
        active = data["active"]
    fig, axs = plt.subplots(len(t["cases"]), 3, figsize=(7.8, 7.6), layout="constrained")
    with np.load(OUT / "timed_pipeline_arrays.npz") as saved:
        for j, case in enumerate(t["cases"]):
            with np.load(PRIVATE / f"{case}.npz") as data:
                truth = data["material"]
            vals = [truth, saved[f"{case}__shifted__material"], saved[f"{case}__top_early__material"]]
            colnames = ["True material", "Shifted-box pipeline", "Point ranking to first fit"]
            for i, value in enumerate(vals):
                ax = axs[j, i]
                image = ax.imshow(embed(value, active), cmap="viridis", vmin=0., vmax=CAP,
                                  interpolation="nearest")
                ax.set_xlim(6.5, 18.5)
                ax.set_ylim(18.5, 6.5)
                ax.set(xticks=[], yticks=[])
                if j == 0:
                    ax.set_title(colnames[i], fontsize=11, pad=8)
                if i == 0:
                    ax.set_ylabel(case.split("_")[2][:6], rotation=0, labelpad=28, fontsize=10)
                else:
                    scores = scored[case]["shifted" if i == 1 else "top_early"]
                    ax.text(.5, -.11, f"Material relative L² {100*scores['full_image_relative_l2']:.1f}%",
                            transform=ax.transAxes, ha="center", va="top", fontsize=8)
    fig.colorbar(image, ax=axs, shrink=.87, pad=.02, label="Real contrast")
    fig.suptitle("Material comparisons on the common 49-cell support", fontsize=14)
    write(fig, "structured_materials")

    fig, (ax, timing) = plt.subplots(2, 1, figsize=(10.2, 6.1),
                                     gridspec_kw={"height_ratios": [1., 1.]}, layout="constrained")
    decision = np.zeros((len(t["cases"]), len(BRANCHES)), dtype=int)
    proof_nu, comparisons, open_truth = [], [], []
    for j, (case, item) in enumerate(examined["cases"].items()):
        for i, row in enumerate(item["branch_records"]):
            if row["shifted"]["strict_lower"] > item["radius"]:
                decision[j, i] = 1
                if row["neumann"]["directional_lower"]>item["radius"]:
                    decision[j, i] = 2
                if row["neumann"].get("strict_contraction_lower", 0.) > 1.:
                    comparisons.append((i, j))
                if row["shifted"]["proof"]["negative_inertia"] > 0:
                    proof_nu.append((i, j))
            if scored_all[case]["metrics"][i]["truth_in_box"]:
                open_truth.append((i, j))
    palette = ["#eee8dc", "#66b1ac", "#3b6b98"]
    ax.imshow(decision, cmap=ListedColormap(palette), vmin=0, vmax=2, aspect="auto")
    labels = [x.replace("_", " ") for x in BRANCHES]
    ax.set_xticks(np.arange(len(labels)), labels, rotation=37, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(t["cases"])), [c.split("_")[2][:6] for c in t["cases"]], fontsize=8)
    if comparisons:
        xx, yy = zip(*comparisons)
        ax.scatter(xx, yy, marker="x", c="#173945", linewidth=1.4, s=35)
    if proof_nu:
        xx, yy = zip(*proof_nu)
        ax.scatter(np.array(xx)+.26, np.array(yy)-.20, marker="*", color="#803533", s=48)
    if open_truth:
        xx, yy = zip(*open_truth)
        ax.scatter(xx, yy, marker="o", facecolors="none", edgecolors="#37292c", s=31, linewidth=1.1)
    ax.set_title("Strict whole-box exclusions on the 24² control", loc="left", fontsize=12)
    ax.legend(handles=[Patch(color=c, label=l) for c, l in zip(palette, ("Retained", "Shifted only", "Both"), strict=True)],
              frameon=False, fontsize=8, ncol=3, loc="upper right", bbox_to_anchor=(1., 1.25))
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(0, -.55, "○ truth inside    × Neumann factor > 1    ★ negative inertia > 0",
            transform=ax.transAxes, fontsize=9, va="top")

    logged = ("top_early", "all_cluster", "neumann_directional", "shifted")
    readable = ("Point to first fit", "All clusters", "Directional Neumann", "Shifted box")
    pal = ["#6d6d73", "#b78357", "#3b6b98", "#2e877f"]
    pos = np.arange(len(t["cases"]))
    for method, name, color in zip(logged, readable, pal, strict=True):
        walls = [row["methods"][method]["wall_seconds"] for row in t["cases"].values()]
        offset = (logged.index(method)-1.5)*.16
        timing.bar(pos+offset, walls, width=.15, label=name, color=color)
    timing.set_yscale("log")
    timing.set_ylim(.018, 1.45)
    timing.spines[["top","right"]].set_visible(False)
    timing.set_xticks(pos, [c.split("_")[2][:6] for c in t["cases"]], fontsize=8)
    timing.set_ylabel("Online wall seconds (log scale)")
    timing.set_title("Timed independent runs, including generation, screening and endpoint validation",
                     loc="left", fontsize=11)
    timing.legend(loc="upper center", bbox_to_anchor=(.56, 1.01), frameon=True,
                  facecolor="white", edgecolor="none", ncol=2, fontsize=8)
    write(fig, "structured_exclusions_and_cost")


def flow_plot():
    manifest = payload(FLOW_OUT / "online.json")
    continuation = payload(FLOW_OUT / "continuation.json")["current_pass"]
    fig, axs = plt.subplots(len(manifest["cases"]), 4, figsize=(9.8, 7.3), layout="constrained")
    with np.load(FLOW_OUT / "offline_geometry_scores.npz") as offline, \
            np.load(FLOW_OUT / "continuation_arrays.npz") as later:
        for i, (case, receipt) in enumerate(manifest["cases"].items()):
            with np.load(FLOW_OUT / f"{case}.npz") as saved:
                top = int(np.argmin(receipt["point_residuals"]))
                best = receipt["pipelines"]["all_cluster"]["winner"]
                branches = continuation[case]["refined_box_results"]
                best_after = min(branches, key=lambda branch: branches[branch]["residual"])
                materials = [offline[f"{case}__truth"], saved["samples"][top], saved["refined"][best],
                             later[f"{case}__box{int(best_after):02d}__material"]]
                radius = receipt["radius"]
                residuals = [receipt["point_residuals"][top], receipt["branch_cost"][best]["residual"],
                             branches[best_after]["residual"]]
                captions = [case.split("_")[3].split("__")[0]] + [f"Residual / radius {r/radius:.2f}" for r in residuals]
                names = ["Archived true material", "Best sampled material", "9 evaluations / box", "39 evaluations / box"]
                for j, material in enumerate(materials):
                    curr = axs[i,j]
                    matrix = curr.imshow(material.reshape(32,32), interpolation="nearest",
                                         cmap="RdBu_r", vmin=-2., vmax=2.)
                    if i == 0:
                        curr.set_title(names[j], fontsize=10, pad=6)
                    curr.set(xticks=[], yticks=[])
                    curr.text(.5, -.09, captions[j], ha="center", va="top", fontsize=8, transform=curr.transAxes)
    fig.colorbar(matrix, ax=axs, shrink=.86, label="Signed contrast")
    fig.suptitle("A13 Flow Matching with its archived 32² observation geometry", fontsize=13)
    write(fig, "flow_native_materials")


if __name__ == "__main__":
    structured_plot()
    flow_plot()
