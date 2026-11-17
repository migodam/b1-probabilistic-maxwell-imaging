"""Summarize the physical cell experiments and draw the comparison figure."""
from collections import defaultdict
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
NAMES = dict(ring="Compact ring", arc="Compact arc",
             wide_ring="Expanded ring", wide_arc="Expanded arc")
COLORS = dict(ring="#aa5038", arc="#da8a3b", wide_ring="#174f7d", wide_arc="#338a82")
MAIN_SIGMA = 1e-4


def read():
    raw = json.loads((RESULTS/"physical_raw.json").read_text())
    groups = defaultdict(list)
    with (RESULTS/"physical_individual_rows.csv").open(newline="") as f:
        for row in csv.DictReader(f):
            fields = {"aperture", "arm", "covers_truth", "uninformative"}
            numbers = {k: float(v) if v else None for k, v in row.items() if k not in fields}
            row.update(numbers)
            groups[row["aperture"],row["sigma"],row["arm"]].append(row)
    return raw, groups


def stats(a):
    values = np.array([v for v in a if v is not None])
    if len(values) == 0:
        return None
    return dict(n=len(values), min=float(np.min(values)), q25=float(np.quantile(values,.25)),
                median=float(np.median(values)), q75=float(np.quantile(values,.75)), max=float(np.max(values)))


def summarize(raw, groups):
    out = dict(cpu_seconds_experiment=raw["cpu_seconds"], designed={}, comparisons=[])
    for geom, design in raw["design"].items():
        check = next(v for v in raw["physics_checks"] if v["aperture"] == geom)
        checks = check["physics_rows"]
        out["designed"][geom] = {
            **{k:v for k,v in design.items() if k not in
               ("cells","active","host","sources","receiver","focus","focus_coordinate_matrix","Q_eigenvalues")},
            **{key: stats([r[key] for r in checks]) for key in (
                "target_clean_reciprocal_abs", "target_identity_error", "relative_Born_response_disagreement",
                "max_full_reduced_library_error", "full_VIE_relative_residual", "full_VIE_sigma_min")},
        }
    inversions = defaultdict(list)
    for trial in raw["observations"]:
        inv = trial["interval_inversion"]
        if inv is not None:
            inversions[trial["aperture"],trial["sigma"],trial["arm"]].append(inv)
    for (geom, sigma, arm), rows in groups.items():
        sets = inversions[geom,sigma,arm]
        quantities = (
            "bounded_SNR_dB", "interval_width", "interval_error", "interval_cpu_s",
            "trf_error", "trf_cpu_s", "trf_reduced_lu", "trf_reduced_rhs",
            "witness_spread", "direct_error", "total_energy",
        )
        out["comparisons"].append(
            dict(aperture=geom, sigma=sigma, arms=arm, n=len(rows),
                 covers=sum(v["covers_truth"] == "True" for v in rows),
                 entire_box_hull=sum(v["uninformative"] == "True" for v in rows),
                 entire_box_scalar_set=sum(v["components"] == [[.5,30.]] for v in sets),
                 disconnected_scalar_set=sum(len(v["components"]) > 1 for v in sets),
                 **{key: stats([v[key] for v in rows]) for key in quantities})
        )
    (RESULTS/"physical_summary.json").write_text(json.dumps(out, indent=2)+"\n")
    return out


def geometry_panel(ax, raw):
    d = raw["design"]["wide_ring"]
    arc = np.array(raw["design"]["wide_arc"]["sources"])
    cells = np.array(d["cells"])
    active = np.array(d["active"])
    host = np.array(d["host"])
    sources = np.array(d["sources"])
    receiver = np.array(d["receiver"]).ravel()
    ax.add_patch(Circle((0,0), 1.55, fc="none", ec="#cfd4dc", lw=.9))
    ax.add_patch(Rectangle((-.476,-.476), .952,.952, fc="none", ec="#8e99a6", ls=":", lw=1))
    for position, value in zip(cells, host):
        color = plt.cm.Blues(.22 + .10 * value)
        ax.add_patch(Circle(position, d["radius"], linewidth=.45, ec=color, fc=color if value else "none"))
    ax.scatter(*sources.T, s=18, facecolors="white", edgecolors=COLORS["wide_ring"], lw=1.0, label="Ring feeds")
    ax.scatter(*arc.T, s=13, c=COLORS["wide_arc"], marker="s", label="Arc feeds")
    ax.scatter(*cells[active].T, s=23, marker="+", c="#c26b38", lw=1.4, label="Variable cells")
    ax.scatter(*cells[d["target_cell_index"]], s=110, marker="*", c="#9e3029", label="Target", zorder=5)
    ax.scatter(*receiver, s=38, marker="D", c="#292b2a", label="Receiver", zorder=5)
    ax.set(xlim=(-1.66,1.66), ylim=(-1.66,1.66), xlabel=r"$x/\lambda$", ylabel=r"$y/\lambda$")
    ax.set_aspect("equal")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(.50,-.17), ncol=3, columnspacing=.85, fontsize=8)


def spectrum_panel(ax, raw):
    for geom in NAMES:
        d = raw["design"][geom]
        y = d["incident_singular_values"]
        ax.semilogy(np.arange(1,len(y)+1),y,"o-",color=COLORS[geom],lw=1.7,ms=4,label=NAMES[geom])
    ax.set(xlabel="Singular-value index", ylabel=r"$\sigma_j(BQ^{-1/2})$ in retained modes",
           xticks=np.arange(1,10), ylim=(6e-5,1.15))
    ax.legend(frameon=False, loc="lower left", fontsize=8)


def width_panel(ax, summary):
    for geom in NAMES:
        entries= [d for d in summary["comparisons"] if d["aperture"]==geom and d["sigma"] and d["arms"]=="focused"]
        x = np.array([d["sigma"] for d in entries])
        w = [d["interval_width"] for d in entries]
        ax.loglog(x,[d["median"] for d in w],"o-",color=COLORS[geom],lw=1.7,ms=4,label=NAMES[geom])
        ax.fill_between(x,[d["q25"] for d in w],[d["q75"] for d in w],color=COLORS[geom],alpha=.10)
    ax.axhline(29.5, color="#333333", lw=1, ls=":")
    ax.text(1.12e-6, 32.0, "Entire material box", color="#444444", fontsize=8)
    ax.set(xlabel="Per-shot detector noise bound", ylabel="Two-shot interval width",
           xlim=(7.5e-7,3e-3), ylim=(2e-3,62))


def errors_panel(ax, summary, geom):
    specs = (
        ("focused", "interval_error", "Two-shot interval midpoint", "#174f7d", "o", "-"),
        ("focused", "trf_error", "Same data multistart TRF", "#a44a3f", "d", "--"),
        ("random", "trf_error", "Random two-shot TRF", "#cc8834", "s", "-"),
        ("completion_equal_energy", "trf_error", "9-shot TRF, energy 2", "#25a28e", "^", "-"),
        ("completion_equal_energy", "direct_error", "9-shot direct, energy 2", "#71717a", "x", "--"),
    )
    for arm, key, label, color, marker, style in specs:
        rows = [d for d in summary["comparisons"] if d["aperture"]==geom and d["sigma"] and d["arms"]==arm]
        ax.loglog([d["sigma"] for d in rows], [d[key]["median"] for d in rows], ls=style, marker=marker,
                  color=color, ms=4, lw=1.65,label=label)
    ax.set(xlabel="Per-shot detector noise bound", ylabel="Median absolute target error",
           xlim=(7.5e-7,3e-3))
    ax.legend(frameon=False, loc="lower right", fontsize=7.6)


def material_panel(ax, raw, geom):
    rows = [o for o in raw["observations"] if o["aperture"]==geom and o["sigma"]==MAIN_SIGMA]
    def select(arm):
        return sorted([o for o in rows if o["arm"]==arm], key=lambda row: row["truth"])
    focus = select("focused")
    random = select("random")
    ax.plot([.5,30],[.5,30],color="#71717a",ls=":",lw=1.2)
    for f, r in zip(focus,random):
        truth = f["truth"]
        proposal = f["interval_inversion"]
        if proposal["interval"] is not None:
            ax.vlines(truth-.16, *proposal["interval"], color="#246a9c", lw=1.3, alpha=.78)
            ax.scatter(truth-.16, proposal["estimate"],c="#246a9c",s=12, zorder=4)
        witnesses = r["feasible_profile"]["independently_validated_witnesses"]
        if len(witnesses) == 2:
            ax.vlines(truth+.16, witnesses[0]["target"], witnesses[1]["target"], color="#c98330", lw=1.15, alpha=.7)
            ax.scatter([truth+.16]*2, [v["target"] for v in witnesses], c="#c98330", s=7, zorder=3)
    handles = [
        Line2D([0],[0],color="#246a9c",marker="o",ls="-",lw=1,ms=4,label="Focused outer interval"),
        Line2D([0],[0],color="#c98330",marker="o",ls="-",lw=1,ms=4,label="Random feasible witnesses"),
    ]
    ax.set(xlabel="True target contrast", ylabel="Target contrast", xlim=(0,30.5), ylim=(0,30.5))
    return handles


def draw(raw, summary):
    plt.rcParams.update({
        "font.family":"DejaVu Sans", "axes.spines.top":False, "axes.spines.right":False,
        "axes.titlesize":11, "axes.labelsize":9.5, "xtick.labelsize":8.5, "ytick.labelsize":8.5,
        "pdf.fonttype":42, "axes.labelpad":5,
    })
    fig, axs = plt.subplots(2,3,figsize=(13.2, 8), gridspec_kw=dict(wspace=.34,hspace=.48))
    geometry_panel(axs[0,0],raw)
    spectrum_panel(axs[0,1],raw)
    width_panel(axs[0,2],summary)
    errors_panel(axs[1,0],summary,"wide_ring")
    material_legend = material_panel(axs[1,1],raw,"wide_ring")
    material_panel(axs[1,2],raw,"ring")
    names = ("Known host and apertures","Physical illumination spectrum","Two-shot uncertainty",
             "Expanded ring errors",r"Expanded ring at $10^{-4}$ noise",
             r"Compact ring at $10^{-4}$ noise")
    for number, (ax,name) in enumerate(zip(axs.flat,names)):
        ax.set_title(f"({chr(97+number)}) {name}",loc="left",pad=8,fontweight="bold")
        if number:
            ax.grid(axis="y",alpha=.17)
            ax.set_axisbelow(True)
    fig.subplots_adjust(top=.925,bottom=.11,left=.075,right=.985)
    fig.legend(handles=material_legend, frameon=False, loc="lower right", ncol=2,
               bbox_to_anchor=(.982,.012), fontsize=8.5)
    fig.savefig(FIGURES/"physical_comparison.png",dpi=240,facecolor="white")
    fig.savefig(FIGURES/"physical_comparison.pdf",facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    raw, groups = read()
    summary = summarize(raw,groups)
    draw(raw,summary)
