#!/usr/bin/env python3
"""Render figures from the toy replay and independent exact certificates."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np

from toy_models import CurrentProblem, complex_array, frozen_model, wedge
from toy_replay import make_fiber, projection


COLORS = ("#316581", "#b3563e", "#407764", "#674c8b")
ARM_LABELS = [
    ("minors", "Minors LM"), ("real_box", "Real-box elimination"),
    ("trf_gn", "Direct TRF-GN"), ("projected_lm", "Direct projected LM"),
]


def read_json(directory, name):
    return json.loads((directory/(name+".json")).read_text())


def plot_basins(directory, source_dir, fig_dir):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    for column, (name, title) in enumerate([("soft3", "3 disks · soft data"),
                                              ("hard5", "5 disks · exact data fiber")]):
        result = read_json(directory, name)
        model, source = frozen_model(source_dir, name)
        j_true = model.forward(model.truth)
        bad = np.array(result["prescribed_bad"]["position"])
        if name == "soft3":
            problem = CurrentProblem(model, np.eye(len(j_true)), np.zeros_like(j_true),
                                     include_data=True, weight=source["weight"])
        else:
            problem = make_fiber(model, complex_array(source, "N"))
        true = problem.coords(j_true)
        t = np.linspace(-.1, 1.15, 251)
        z = true + t[:,None]*(bad-true)
        ax = axes[0,column]
        for index, (form, label) in enumerate(ARM_LABELS[:2]):
            values = np.array([np.linalg.norm(problem.residual_and_jac(point, form)[0])**2
                               for point in z])
            at_bad = np.linalg.norm(problem.residual_and_jac(bad, form)[0])**2
            ax.plot(t, values/at_bad, label=label.replace(" LM", ""), color=COLORS[index], lw=2)
        if name == "hard5":
            _, (p, *_rest) = projection(problem, j_true)
            def canonical_value(point):
                j = problem.current(point)
                return np.linalg.norm(p@wedge(j, model.field(j)))**2
            ax.plot(t, [canonical_value(point)/canonical_value(bad) for point in z],
                    linestyle="--", label="Curvature projection", color=COLORS[2], lw=2)
        ax.axvline(1., color="#777777", lw=1, ls=":")
        ax.set_yscale("symlog", linthresh=1)
        ax.set(xlabel="Position along current chord (truth 0, original bad point 1)",
               ylabel="Objective ÷ value at original bad point",
               title=title)
        ax.legend(frameon=False)
        ax.grid(alpha=.2)
        ax = axes[1,column]
        for index, (form, label) in enumerate(ARM_LABELS):
            vals = [r[form] for r in result["legally_paired_material_starts"]]
            x = [v["box_material_relative_error"] for v in vals]
            y = [v["physical_data_relative"] for v in vals]
            ax.scatter(x, y, facecolors="none", edgecolors=COLORS[index], alpha=.63,
                       marker=("o", "s", "^", "D")[index], label=label, s=43+10*index)
            ax.scatter(x[0], y[0], marker=("*", "P", "^", "D")[index],
                       edgecolors="white", color=COLORS[index], s=76+10*index)
        ax.set(xlabel="Real material relative error", ylabel="Full-forward port relative residual",
               title="19 paired material starts")
        ax.xaxis.set_major_formatter(PercentFormatter(1))
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.legend(frameon=False, loc="upper left", fontsize=8)
        ax.grid(alpha=.2)
    fig.suptitle("Changing the residual moves the original stationary point", fontsize=15)
    fig.savefig(fig_dir/"toy_basins.png", dpi=160)
    plt.close(fig)


def plot_sweep(directory, fig_dir):
    rows = read_json(directory, "sweep")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.05), layout="constrained")
    for i, (layout, label) in enumerate([("wide", "Wide"), ("one_sided", "One-sided")]):
        selected = [row for row in rows if row["layout"] == layout]
        x = [row["peak"] for row in selected]
        for ax, field, name in [(axes[0], "cond_projected", "Projected inverse"),
                               (axes[1], "material_relative_error", "Material readout"),
                               (axes[1], "physical_data_relative", "Physical ports")]:
            marker = "o" if field in ("cond_projected", "material_relative_error") else "s"
            ls = "-" if marker == "o" else "--"
            ax.loglog(x, [row[field] for row in selected], color=COLORS[i],
                      marker=marker, ls=ls, label=f"{label} · {name}", lw=1.8)
    axes[0].loglog([row["peak"] for row in selected], [row["cond_forward"] for row in selected],
                   color="#555555", marker="s", ls="--", label="Shared forward operator", lw=1.8)
    for ax in axes:
        ax.legend(frameon=False, fontsize=8)
        ax.grid(alpha=.2, which="both")
        ax.set_xlabel("Peak contrast")
    axes[0].set_ylabel("Condition number")
    axes[1].set_ylabel("Relative error")
    fig.suptitle("Five-disk curvature projection · two receiver layouts · noise-free")
    fig.savefig(fig_dir/"toy_contrast.png", dpi=160)
    plt.close(fig)


def plot_feedback(directory, fig_dir):
    values = read_json(directory, "feedback_sweep")
    cert = read_json(directory, "certificates")["tail"]
    x = [v["epsilon"] for v in values]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.05), layout="constrained")
    axes[0].plot(x, [v["curvature_sufficient_bound"] for v in values], color=COLORS[0],
                 label="Sufficient Hessian lower bound", marker="o")
    axes[0].plot(x, [v["sampled_min_hessian"] for v in values], color=COLORS[1],
                 label="Smallest sampled Hessian eigenvalue", marker="s")
    axes[0].set_yscale("symlog", linthresh=.06)
    axes[0].set_yticks([-1000, -100, -10, -1, -.1, 0, .1, 1])
    axes[0].set_ylabel("Minimum curvature")
    axes[1].plot(x, [v["directional_witness"] for v in values], color=COLORS[2],
                 label="Directional curvature at witness point", marker="^")
    axes[1].scatter([.1], [float(cert["negative"]["curvature_display"])], color=COLORS[1],
                    s=110, facecolor="none", linewidth=1.6, label="Exact negative witness")
    axes[1].set_ylabel("Directional second derivative")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("Tail magnitude")
        ax.axhline(0, lw=1, color="#777777", ls=":")
        ax.grid(alpha=.2, which="both")
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Fixed true feedback · a contractive tail can produce negative curvature")
    fig.savefig(fig_dir/"toy_feedback.png", dpi=160)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--results", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    plot_basins(args.results, args.source/"research_work", args.out)
    plot_sweep(args.results, args.out)
    plot_feedback(args.results, args.out)


if __name__ == "__main__":
    main()
