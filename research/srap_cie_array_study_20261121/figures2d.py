#!/usr/bin/env python3
"""Render planar array designs and scored SRAP-CIE figures without rerunning an inverse."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from html import escape
import json
import os
from pathlib import Path

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle
from matplotlib.ticker import MaxNLocator
import numpy as np

import data2d
import run_study as run

FIGURES = Path(__file__).resolve().parent / "figures"
INK = "#142438"
MUTED = "#536579"
GRID = "#d9e4ed"
DUAL = "#16688c"
SRAP = "#bd511c"
TX = "#c2462d"
RX = "#187fa5"
COLORS = dict(dual=DUAL, srap=SRAP, u11="#626a93", reset="#42916c",
              full_ap="#936b2d", graph="#8f4a85", dual_cap="#4a9ea3", srap_cap="#d17642",
              srap_lazy="#89509e")
METHODS = dict(dual="Dual-CIE", srap="SRAP-CIE", u11="U11 routed AP",
               reset="PR reset", full_ap="Full AP", graph="Graph-compatible AP",
               dual_cap="Dual-CIE capped", srap_cap="SRAP-CIE capped", srap_lazy="SRAP-CIE lazy")
ARRAY_LABELS = {
    "PW360_10x20": ("360° plane wave", "Sparse full view"),
    "PW180_10x20": ("180° plane wave", "Historical half view"),
    "PW120_10x20": ("120° plane wave", "Historical narrow window"),
    "PW180_FLIP_RX_10x20": ("180° plane wave", "Receiver arc flipped"),
    "LS_LINE_SAME_10x20": ("Line source", "Same-side line RX"),
    "LS_LINE_OPPOSITE_10x20": ("Line source", "Opposite-side line RX"),
    "LS_ARC_10x20": ("Line source", "Same-angle arc RX"),
    "PW120_CENTERED_10x20": ("120° plane wave", "Centered on 180° axis"),
    "PW180_FLIP_TX_10x20": ("180° plane wave", "Illumination flipped"),
    "PW180_7x14": ("180° plane wave", "7 TX × 14 RX"),
    "PW120_7x14": ("120° plane wave", "7 TX × 14 RX"),
}
CMAP = plt.get_cmap("viridis")
MATERIAL_NORM = Normalize(0, 6, clip=True)
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "text.color": INK,
    "axes.labelcolor": INK, "axes.edgecolor": GRID, "xtick.color": MUTED,
    "ytick.color": MUTED, "svg.fonttype": "none", "pdf.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white",
})


def save_figure(fig, folder, name, *, vector=True):
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / f"{name}.png", dpi=150, facecolor="white")
    if vector:
        fig.savefig(folder / f"{name}.pdf", facecolor="white")
    plt.close(fig)


def draw_geometry(folder, *, followup=False):
    if followup:
        from data_followup import ARRAYS
        names = ARRAYS
        fig, axes = plt.subplots(1, 5, figsize=(14.4, 4.5))
        fig.subplots_adjust(left=.045, right=.98, bottom=.09, top=.79, wspace=.19)
        title = "Aperture, orientation and illumination controls"
        subtitle = "400 MHz  /  equal channel count  /  fixed receiver geometry for the illumination flip"
    else:
        names = data2d.CONFIRM_ARRAYS
        fig, axes = plt.subplots(2, 4, figsize=(12.8, 7.8))
        fig.subplots_adjust(left=.045, right=.98, bottom=.07, top=.86, wspace=.16, hspace=.35)
        title = "Seven planar arrays at a fixed channel count"
        subtitle = "400 MHz scalar TM  /  10 illuminations × 20 receivers per array"
    fig.text(.045, .925 if followup else .956, title, fontsize=23, weight="bold")
    fig.text(.046, .853 if followup else .913, subtitle,
             fontsize=11, color=MUTED)
    for i, name in enumerate(names):
        ax = axes.flat[i]
        a = data2d.configuration(name)
        ax.set(xlim=(-4.3, 4.3), ylim=(-4.3, 4.3), aspect="equal",
               xticks=(-3.75, 0, 3.75), yticks=(-3.75, 0, 3.75))
        ax.tick_params(labelsize=7, length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.axhline(0, color=GRID, linewidth=.8)
        ax.axvline(0, color=GRID, linewidth=.8)
        ax.add_patch(Rectangle((-1, -1), 2, 2, facecolor="#eaf1f6", edgecolor="#74889b", linewidth=.9))
        if a.tx_kind == "plane_wave" or name == "LS_ARC_10x20":
            ax.add_patch(Circle((0, 0), data2d.RANGE, fill=False,
                                edgecolor=GRID, linestyle=(0, (3, 4)), linewidth=.9))
        ax.scatter(*a.receiver_points.T, color=RX, s=19, edgecolors="white", linewidth=.35, zorder=3)
        if a.tx_kind == "plane_wave":
            directions = data2d.unit_from_degrees(a.tx[:, 0])
            starts = -3.26 * directions
            vectors = .90 * directions
            ax.quiver(*starts.T, *vectors.T, angles="xy", scale_units="xy", scale=1,
                      color=TX, width=.007, headwidth=4, headaxislength=4.5, zorder=4)
        else:
            ax.scatter(*a.tx.T, marker="^", color=TX, s=35,
                       edgecolors="white", linewidth=.45, zorder=4)
        title, description = ARRAY_LABELS[name]
        ax.set_title(title, loc="left", fontsize=11, weight="bold", pad=13)
        ax.text(0, 1.015, description, transform=ax.transAxes, color=MUTED, fontsize=8.6)

    ax = axes.flat[-1]
    ax.axis("off")
    handles = [
        Line2D([], [], linestyle="None", marker=r"$\longrightarrow$", markersize=18,
               color=TX, label="Plane-wave direction d"),
        Line2D([], [], linestyle="None", marker="^", markersize=8,
               color=TX, label="Physical line source"),
        Line2D([], [], linestyle="None", marker="o", markersize=7,
               color=RX, label="Physical receiver"),
        Rectangle((0, 0), 1, 1, facecolor="#eaf1f6", edgecolor="#74889b", label="2 m × 2 m imaging domain"),
    ]
    if followup:
        handles.pop(1)
        note = "The 180° and centered 120°\nwindows share TX and RX axes.\n\nArrows are propagation directions.\nCoordinates are in metres."
    else:
        note = ("Arrows point along propagation.\nPhysical coordinates are in metres.\n\n"
                "The line/arc pair shares the same\nline sources and receiver bearings.")
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, .95),
              frameon=False, labelspacing=1.25, fontsize=9.3, borderpad=0)
    ax.text(0, .46 if followup else .40, note,
            transform=ax.transAxes, color=MUTED, va="top", fontsize=9.4, linespacing=1.5)
    save_figure(fig, folder, "followup_arrays" if followup else "planar_arrays")


def draw_fresnel_masks(folder):
    cases = {c["array"]: c for c in run.cases() if c["cohort"] == "measured_extra"}
    variants = ("FRESNEL_CONTIGUOUS", "FRESNEL_INTERLEAVED")
    if not all(k in cases for k in variants):
        return
    fig, axes = plt.subplots(2, 1, figsize=(10.2, 6.2), sharex=True, sharey=True)
    fig.subplots_adjust(left=.12, right=.96, bottom=.15, top=.79, hspace=.45)
    fig.text(.12, .95, "Fresnel Twin, 8 GHz · matched real-channel masks",
             fontsize=20, weight="bold")
    fig.text(.12, .897, "733 training channels in each mask  /  common holdout: 1,205 TX-side, 390 RX-side",
             color=MUTED, fontsize=10)
    for ax, variant, title in zip(axes, variants, (
            "Continuous half circle", "Interleaved sparse circle"), strict=True):
        c = cases[variant]
        # Read indices from public data. No held-out fields are opened.
        with np.load(run.public_path(c), allow_pickle=False) as z:
            tx, rx = np.asarray(c["train_sources"]) - 1, z["rx_ids"]
        active = np.zeros((18, 360), dtype=bool)
        for source, indices in zip(tx, rx, strict=True):
            active[source, indices[indices >= 0]] = True
        ax.imshow(active, origin="lower", aspect="auto", interpolation="none",
                  extent=(-.5, 359.5, -.5, 17.5), cmap=ListedColormap(["#eaf0f4", RX]))
        ax.set_yticks([0, 4, 8, 12, 16], [f"{s * 20}°" for s in (0, 4, 8, 12, 16)])
        ax.set_title(title, loc="left", fontsize=11, weight="bold", pad=7)
        ax.set_ylabel("TX angle")
        ax.tick_params(length=0)
    axes[-1].set_xticks(np.arange(0, 361, 60), [f"{x}°" for x in range(0, 361, 60)])
    axes[-1].set_xlabel("Receiver angle")
    fig.text(.12, .035, "Blue shows real training channels. Gray includes held-out and unavailable channels.",
             fontsize=9, color=MUTED)
    save_figure(fig, folder, "fresnel_twin_masks")


def draw_receiver_spectra(folder):
    cases = {c["array"]: c for c in run.cases() if c["object_id"] == "s2d_001"}
    if not all(name in cases for name in data2d.CONFIRM_ARRAYS):
        return
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.8), sharey=True)
    fig.subplots_adjust(left=.09, right=.975, bottom=.16, top=.72, wspace=.15)
    fig.text(.09, .927, "Receiver spectra · equal 20-channel receivers", fontsize=20, weight="bold")
    fig.text(.09, .857, "Singular values of the public S maps at 32²  /  normalized within each array",
             color=MUTED, fontsize=10)
    labels = ("360° full", "180° half", "120° narrow", "180° flipped RX",
              "Same-side line RX", "Opposite-side line RX", "Same-angle arc RX")
    palette = (DUAL, SRAP, "#8d4a89", "#e2a06e", DUAL, SRAP, "#50a267")
    styles = ("-", "-", "-", "--", "-", "--", "-")
    for i, (name, label, color, style) in enumerate(zip(
            data2d.CONFIRM_ARRAYS, labels, palette, styles, strict=True)):
        ax = axes[int(i >= 4)]
        with np.load(run.public_path(cases[name]), allow_pickle=False) as z:
            singular = np.linalg.svd(z["S"], compute_uv=False)
        singular /= singular[0]
        rank = np.count_nonzero(singular >= .1)
        ax.semilogy(np.arange(1, singular.size + 1), singular, color=color,
                    linestyle=style, linewidth=1.6, label=f"{label}   ({rank})")
    for ax, title in zip(axes, ("Plane-wave arrays", "External line-source arrays"), strict=True):
        ax.set_title(title, loc="left", fontsize=11, weight="bold", pad=8)
        ax.axhline(.1, color=MUTED, linestyle=":", linewidth=1.2)
        ax.set(xlim=(1, 20), xticks=[1, 5, 10, 15, 20], xlabel="Index")
        ax.grid(color=GRID, linewidth=.7)
        ax.legend(frameon=False, fontsize=8.1, loc="lower left")
    axes[0].set_ylabel(r"$\sigma_i(S)/\sigma_1(S)$")
    fig.text(.09, .035, "Numbers in parentheses count receiver modes at or above the dotted 0.1 threshold.",
             color=MUTED, fontsize=8.8)
    save_figure(fig, folder, "receiver_spectra")


def csv_rows(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def parse_number(v):
    return float(v) if v not in (None, "", "None") else np.nan


def prepare_scored(cohort, required_arms, budget):
    from score_study import expected, prediction_path

    result = run.ROOT / "results"
    path = result / f"{cohort}_metrics.csv"
    if not path.is_file():
        raise SystemExit(f"{cohort}: waiting for sealed scoring file {path.name}")
    cases = [c for c in run.cases() if c["cohort"] == cohort]
    if not cases:
        raise SystemExit(f"{cohort}: case manifest has no matching cases")
    entries = {(r["case"], r["method"], int(r["budget"]), r["view"]): r for r in csv_rows(path)}
    # Every required score and prediction must exist before this function
    # returns to a material plot that can open truth.
    missing = [(c["case"], method, view) for c in cases
               for method in sorted(set(required_arms) | set(expected(c)))
               for view in ("raw", "light")
               if (c["case"], method, budget, view) not in entries
               or not prediction_path(c, method).is_file()]
    if missing:
        raise SystemExit(f"{cohort}: wait for {len(missing)} required scored predictions")
    return cases, entries


def prediction_pair(case, method, budget, directory=None):
    from score_study import prediction_path

    path = (prediction_path(case, method) if directory is None
            else directory / f"{case['case']}__{method}.npz")
    with np.load(path, allow_pickle=False) as z:
        return [z[f"m_{budget}_{view}"].reshape(case["side"], case["side"])
                if f"m_{budget}_{view}" in z else None
                for view in ("raw", "light")]


def materials_truth(case):
    # Caller has verified the scorer completed, including all required cases.
    if case["cohort"] == "legacy":
        from score_study import get_truth
        return get_truth(case, run.load_public(case)).reshape(case["side"], case["side"])
    with np.load(run.DATA / "truth" / f"{case['case']}.npz", allow_pickle=False) as z:
        return z["m"].reshape(case["side"], case["side"])


def material_panel(ax, m, norm):
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    if m is None:
        ax.set_facecolor("#eef2f5")
        ax.text(.5, .5, "No eligible\nsnapshot", color=MUTED,
                ha="center", va="center", transform=ax.transAxes, fontsize=8)
    else:
        ax.imshow(m, origin="lower", cmap=CMAP, norm=norm, interpolation="nearest")


def synthetic_panel_label(row):
    if row is None:
        return ""
    l2, iou = (parse_number(row.get(key)) for key in ("rel_L2", "half_IoU"))
    return f"L2 {l2:.3f}  ·  IoU {iou:.3f}" if np.isfinite(l2) else ""


def measured_panel_label(row):
    if row is None:
        return ""
    rx, tx = (parse_number(row.get(key)) for key in ("held_rx", "held_tx"))
    labels = ([f"TX {tx:.3f}"] if np.isfinite(tx) else []) + ([f"RX {rx:.3f}"] if np.isfinite(rx) else [])
    return "  ·  ".join(labels)


def short_name(case):
    if case["cohort"] == "confirmation":
        line1, line2 = ARRAY_LABELS[case["array"]]
        return f"{line1}\n{line2}"
    if case["cohort"] == "measured_extra":
        return "Twin, 8 GHz\n" + ("Continuous" if "contiguous" in case["case"] else "Interleaved")
    if case["suite"] == "LA_FRESNEL":
        name = case["case"].removeprefix("fresnel_")
        item, freq, mask, side = name.split("_")
        return f"{item.title()}, {freq.replace('GHz', ' GHz')}\n180° {side.upper()}"
    if case["cohort"] in ("development", "followup"):
        line1, line2 = ARRAY_LABELS[case["array"]]
        return f"{line1}\n{line2}"
    pieces = case["case"].split("_")
    aperture = pieces[0].removeprefix("la")
    return f"{aperture}° {' '.join(pieces[1:-1]).replace('rot', 'rot ')}\n$\\epsilon_r$ {pieces[-1].removeprefix('eps').replace('p', '.')}"


def material_page(cases, ledger, arms, budget, title, subtitle, *, shared_truth=None,
                  norm=MATERIAL_NORM, peak_labels=False, prediction_directory=None,
                  views=("raw", "light"), show_object=False):
    synthetic = cases[0]["suite"] == "LA_SYNTH"
    num_rows = len(cases)
    column_titles = (["Reference"] if synthetic else []) + [
        METHODS[arm] for view in views for arm in arms
    ]
    n_panels = len(column_titles)
    fig = plt.figure(figsize=(2.12 + n_panels * 1.65, 1.8 + num_rows * 1.47))
    gs = fig.add_gridspec(
        num_rows, n_panels + 1, width_ratios=[1.45] + [1] * n_panels,
        left=.045, right=.956, top=1 - 1.15 / fig.get_figheight(),
        bottom=.66 / fig.get_figheight(), wspace=.14, hspace=.31,
    )
    fig.text(.045, 1 - .37 / fig.get_figheight(), title, weight="bold", fontsize=18)
    fig.text(.046, 1 - .63 / fig.get_figheight(), subtitle, color=MUTED, fontsize=8.8)
    axes = []
    for i, case in enumerate(cases):
        side = fig.add_subplot(gs[i, 0])
        side.axis("off")
        label_text = f"{case['object_id']}\n{short_name(case)}" if show_object else short_name(case)
        side.text(0, .50, label_text, color=INK, fontsize=8.5, weight="bold",
                  transform=side.transAxes, ha="left", va="center", linespacing=1.4)
        row = [fig.add_subplot(gs[i, j]) for j in range(1, n_panels + 1)]
        axes.append(row)
        j = 0
        if synthetic:
            material_panel(row[0], shared_truth if shared_truth is not None else materials_truth(case), norm)
            j = 1
        models = {arm: prediction_pair(case, arm, budget, prediction_directory) for arm in arms}
        for view in views:
            for arm in arms:
                ax = row[j]
                m = models[arm][int(view == "light")]
                material_panel(ax, m, norm)
                metric = ledger.get((case["case"], arm, budget, view))
                label = synthetic_panel_label(metric) if synthetic else measured_panel_label(metric)
                if peak_labels and m is not None:
                    label = f"{label}\nPeak {np.max(m):.3g}" if label else f"Peak {np.max(m):.3g}"
                ax.set_xlabel(label, labelpad=3, color=MUTED, fontsize=7.3)
                j += 1
    top = axes[0]
    for j, label in enumerate(column_titles):
        j_method = j - int(synthetic)
        view, color = "", INK
        if j_method >= 0:
            view = views[j_method // len(arms)]
            color = COLORS[arms[j_method % len(arms)]]
        top[j].set_title(f"{label}\n{view}" if view else label, loc="left", fontsize=8.8,
                         color=color, pad=5, weight="bold")
    cbax = fig.add_axes([.43, .25 / fig.get_figheight(), .15, .095 / fig.get_figheight()])
    ticks = [0, .25, .5, .75, 1] if norm.vmax == 1 else [0, 2, 4, 6]
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cbax,
                 orientation="horizontal", ticks=ticks)
    cbax.tick_params(labelsize=7, length=2)
    fig.text(.60, .26 / fig.get_figheight(), "Relative permittivity contrast, clipped for display",
             color=MUTED, fontsize=8, va="center")
    return fig


def write_atlas_html(path, items, pdf, budget):
    links = "".join(f'<a href="#{escape(obj)}">{escape(obj)}</a>' for obj in items)
    images = "".join(
        f'<section id="{escape(obj)}"><h2>{escape(obj)}</h2>'
        f'<a href="{escape(obj)}.png"><img loading="lazy" src="{escape(obj)}.png" '
        f'alt="{escape(obj)} 在七种阵列上的材料恢复"></a></section>' for obj in items
    )
    path.write_text(f'''<!doctype html>
<html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>二维确认组材料图集</title><style>
body{{background:#edf3f7;color:#142438;font:16px/1.55 system-ui,sans-serif;margin:0}}
main{{padding:2rem 1rem;max-width:1000px;margin:auto}}
header,section{{padding:1rem;background:white;margin-bottom:1rem;border-radius:8px}}
h1,h2{{line-height:1.15}} h2{{font-size:1.15rem}} nav{{display:flex;flex-wrap:wrap;gap:.4rem}}
nav a{{background:#edf3f7;border-radius:4px;padding:.15rem .45rem}}
a{{color:#17698b}} img{{display:block;width:100%;height:auto}} :target{{outline:2px solid #187fa5}}
</style><main><header><h1>二维确认组材料图集</h1>
<p>24 个材料对象，每个对象 7 种阵列。{budget} 秒预算同时展示直接恢复与轻校正后的结果。全图使用相同的介电对比度色标。</p>
<p><a href="../{escape(pdf)}">浏览整套 PDF</a></p><nav aria-label="对象索引">{links}</nav>
</header>{images}</main></html>
''', encoding="utf-8")


def confirmation_atlas(folder, cases, ledger, budget):
    groups = defaultdict(dict)
    for c in cases:
        groups[c["object_id"]][c["array"]] = c
    if len(groups) != 24 or any(set(g) != set(data2d.CONFIRM_ARRAYS) for g in groups.values()):
        raise RuntimeError("The confirmation atlas requires 24 complete 7-array objects.")
    out = folder / f"confirmation_atlas_{budget}s"
    out.mkdir(parents=True, exist_ok=True)
    pdf_path = folder / f"confirmation_atlas_{budget}s.pdf"
    with PdfPages(pdf_path) as pdf:
        for obj in sorted(groups):
            cs = [groups[obj][a] for a in data2d.CONFIRM_ARRAYS]
            truth = materials_truth(cs[0])
            fig = material_page(cs, ledger, ("dual", "srap"), budget, f"{obj} · material recovery",
                                f"{budget} s  /  400 MHz  /  2 m square, 96² → 32²  /  "
                                "10% complex measurement noise", shared_truth=truth)
            pdf.savefig(fig)
            save_figure(fig, out, obj, vector=False)
    write_atlas_html(out / "index.html", sorted(groups), pdf_path.name, budget)


def line_cap_examples(folder, cases, ledger, budget):
    objects = ("s2d_001", "s2d_012", "s2d_024")
    arrays = ("LS_LINE_SAME_10x20", "LS_LINE_OPPOSITE_10x20", "LS_ARC_10x20")
    index = {(c["object_id"], c["array"]): c for c in cases}
    cs = [index[obj, array] for obj in objects for array in arrays]
    arms = ("dual", "srap", "dual_cap", "srap_cap")
    fig = material_page(cs, ledger, arms, budget, "Line-source arrays · receiver truncation",
                        f"{budget} s light  /  three fixed objects  /  capped modes satisfy "
                        r"$\sigma_i(S)\geq0.1\sigma_1(S)$", views=("light",), show_object=True)
    save_figure(fig, folder, f"line_cap_examples_{budget}s_light")


def legacy_atlas(folder, cases, ledger, budget):
    synth = [c for c in cases if c["suite"] == "LA_SYNTH"]
    fresnel = [c for c in cases if c["suite"] == "LA_FRESNEL"]
    arms = ("dual", "srap", "u11")
    pdf_path = folder / f"legacy_materials_{budget}s.pdf"
    with PdfPages(pdf_path) as pdf:
        slices = (synth[:6], synth[6:], fresnel)
        for j, cs in enumerate(slices, 1):
            measured = cs[0]["suite"] == "LA_FRESNEL"
            title = "Fresnel 2005 · six historical masks" if measured else "Historical synthetic benchmark"
            subtitle = (f"{budget} s  /  independently scored held-out complex fields"
                        if measured else f"{budget} s  /  cross-grid scalar TM  /  cases {6*(j-1)+1}–{6*j} of 12")
            fig = material_page(cs, ledger, arms, budget, title, subtitle)
            pdf.savefig(fig)
            save_figure(fig, folder, f"legacy_materials_{budget}s_p{j}", vector=False)


def legacy_cross_host(folder, budget):
    root = run.ROOT / "results"
    candidates = (
        ("la180_austria_rot90_eps3p5", "Austria, rotated 90°\n$\\epsilon_r$ 3.5"),
        ("la180_core_shell_eps5p0", "Core and shell\n$\\epsilon_r$ 5.0"),
        ("la180_eccentric_shell_eps4p5", "Eccentric shell\n$\\epsilon_r$ 4.5"),
    )
    arms = ("dual", "srap", "u11")
    selected = dict(candidates)
    cases = {c["case"]: c for c in run.cases() if c["cohort"] == "legacy" and c["case"] in selected}
    conditions = (("XINAN", root / "legacy_metrics.csv", run.OUTPUT),
                  ("Mac", root / "legacy_mac_metrics.csv", root / "legacy_mac"))
    groups = []
    for host, score_path, predictions in conditions:
        if not score_path.is_file():
            raise SystemExit(f"{host}: waiting for sealed historical replay scores")
        ledger = {(r["case"], r["method"]): r for r in csv_rows(score_path)
                  if int(r["budget"]) == budget and r["view"] == "light"}
        if any((name, arm) not in ledger or not (predictions / f"{name}__{arm}.npz").is_file()
               for name in selected for arm in arms):
            raise SystemExit(f"{host}: waiting for all selected historical replay predictions")
        groups.append((host, ledger, predictions))

    fig = plt.figure(figsize=(13.7, 6.55))
    fig.text(.048, .944, "Three historical objects replayed on two machines", fontsize=22, weight="bold")
    fig.text(.049, .888, f"{budget} s light  /  same 400 MHz data per object  /  "
             "half-aperture scalar TM", color=MUTED, fontsize=10)
    gs = fig.add_gridspec(3, 8, width_ratios=[1.25] + [1] * 7, left=.048, right=.965,
                          top=.755, bottom=.14, wspace=.195, hspace=.35)
    for i, (case_id, row_name) in enumerate(candidates):
        case = cases[case_id]
        name_ax = fig.add_subplot(gs[i, 0])
        name_ax.axis("off")
        name_ax.text(0, .5, row_name, transform=name_ax.transAxes,
                     fontsize=9, weight="bold", va="center")
        ref_ax = fig.add_subplot(gs[i, 1])
        material_panel(ref_ax, materials_truth(case), MATERIAL_NORM)
        if i == 0:
            ref_ax.set_title("Reference", loc="left", fontsize=9, weight="bold", pad=5)
        for h, (host, ledger, directory) in enumerate(groups):
            for j, arm in enumerate(arms):
                ax = fig.add_subplot(gs[i, 2 + h * 3 + j])
                _, selected_material = prediction_pair(case, arm, budget, directory)
                material_panel(ax, selected_material, MATERIAL_NORM)
                ax.set_xlabel(synthetic_panel_label(ledger[(case_id, arm)]),
                              fontsize=7.5, labelpad=4, color=MUTED)
                if i == 0:
                    ax.set_title(METHODS[arm], loc="left", fontsize=9, color=COLORS[arm],
                                 weight="bold", pad=5)
                    if j == 0:
                        x = ax.get_position().x0
                        fig.text(x, .819, host, fontsize=11, weight="bold")
                        right = gs[i, 2 + h * 3 + 2].get_position(fig).x1
                        fig.lines.append(Line2D([x, right], [.804, .804], transform=fig.transFigure,
                                                color=GRID, linewidth=1.4))
    cbax = fig.add_axes([.42, .055, .15, .016])
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=MATERIAL_NORM, cmap=CMAP), cax=cbax,
                 orientation="horizontal", ticks=[0, 2, 4, 6])
    cbax.tick_params(labelsize=7, length=2)
    fig.text(.59, .062, "Relative permittivity contrast, clipped for display",
             color=MUTED, fontsize=8.5, va="center")
    save_figure(fig, folder, f"legacy_cross_host_{budget}s_light")


def measured_atlas(folder, cases, ledger, budget):
    for suffix, norm in (("", MATERIAL_NORM), ("_low", Normalize(0, 1, clip=True))):
        tab = " · low contrast" if suffix else ""
        fig = material_page(cases, ledger, ("dual", "srap"), budget,
                            "Fresnel Twin, 8 GHz · two training masks" + tab,
                            f"{budget} s / 733 training channels each / common TX and RX holdout / "
                            "no pixel-wise ground truth", norm=norm, peak_labels=True)
        save_figure(fig, folder, f"fresnel_twin_materials_{budget}s{suffix}")


def lazy_fresnel_atlas(folder, budget):
    scored = run.ROOT / "results" / "fresnel_lazy_metrics.csv"
    predictions = run.ROOT / "results" / "lazy_fresnel"
    if not scored.is_file():
        raise SystemExit("Fresnel late-setup figure is waiting for sealed scores.")
    cases = [c for c in run.cases() if c["suite"] == "LA_FRESNEL"]
    arms = ("dual", "srap", "srap_lazy")
    ledger = {(r["case"], r["method"], int(r["budget"]), r["view"]): r for r in csv_rows(scored)}
    if len(cases) != 8 or any(
            (c["case"], arm, budget, view) not in ledger
            or not (predictions / f"{c['case']}__{arm}.npz").is_file()
            for c in cases for arm in arms for view in ("raw", "light")):
        raise SystemExit("Fresnel late-setup figure is waiting for all 24 scored predictions.")
    for scale, norm in (("full", MATERIAL_NORM), ("low", Normalize(0, 1, clip=True))):
        title = "Fresnel 2005 · setup timing comparison" + (" · low contrast" if scale == "low" else "")
        subtitle = (f"{budget} s / eight measured settings / no pixel-wise ground truth / "
                    "TX and RX are held-out errors / peaks use unclipped values")
        fig = material_page(cases, ledger, arms, budget, title, subtitle, norm=norm,
                            peak_labels=True, prediction_directory=predictions)
        save_figure(fig, folder, f"fresnel_lazy_materials_{budget}s_{scale}")


def development_atlas(folder, cases, ledger, budget):
    groups = defaultdict(dict)
    for case in cases:
        groups[case["object_id"]][case["array"]] = case
    names = ("PW180_10x20", "PW180_7x14", "PW120_10x20", "PW120_7x14")
    arms = ("dual", "srap", "dual_cap", "srap_cap")
    with PdfPages(folder / f"development_materials_{budget}s.pdf") as pdf:
        for i, (obj, group) in enumerate(sorted(groups.items()), 1):
            cs = [group[name] for name in names]
            fig = material_page(cs, ledger, arms, budget, obj.replace("_", " "),
                                f"{budget} s  /  two historical angular windows crossed with two densities",
                                shared_truth=materials_truth(cs[0]))
            pdf.savefig(fig)
            save_figure(fig, folder, f"development_materials_{budget}s_p{i}", vector=False)


def scored_examples(cases, ledger, budget, *, eligible=None):
    pairs = []
    for c in cases:
        if c["suite"] != "LA_SYNTH" or (eligible is not None and c["case"] not in eligible):
            continue
        dual, srap = (parse_number(ledger[(c["case"], name, budget, "raw")].get("rel_L2"))
                      for name in ("dual", "srap"))
        if np.isfinite(dual) and np.isfinite(srap) and dual > 0:
            pairs.append(((dual - srap) / dual * 100, c))
    ordered = sorted(pairs, key=lambda pair: (pair[0], pair[1]["case"]))
    return (ordered[-1], ordered[0]) if len(ordered) > 1 else ordered


def route_events(case, arm, budget):
    with np.load(run.output_path(case, arm), allow_pickle=False) as z:
        if "route_json" not in z:
            return []
        events = json.loads(str(z["route_json"]))
    return [e for e in events if e.get("event") == "AP" and e.get("wall", np.inf) <= budget]


def trajectory_plot(folder, cohort, cases, ledger, budget):
    path = run.ROOT / "results" / f"{cohort}_trajectories.csv"
    if not path.is_file():
        raise SystemExit(f"{cohort}: waiting for the scorer's trajectory file {path.name}")
    old = cohort == "legacy"
    examples = scored_examples(cases, ledger, budget, eligible=run.DEV_IDS if old else None)
    if not examples:
        return
    shown = {c["case"] for _, c in examples}
    rows = defaultdict(list)
    for row in csv_rows(path):
        if row["case"] in shown and float(row["seconds"]) <= budget:
            rows[(row["case"], row["method"])].append(row)

    fig, axes = plt.subplots(len(examples), 3, figsize=(13.3, 3.3 * len(examples) + 1.3), squeeze=False)
    fig.subplots_adjust(left=.10, right=.975, bottom=.125, top=.79, wspace=.32, hspace=.57)
    title = "Historical mechanism examples" if old else "Audited error over time · paired extremes"
    fig.text(.10, .955, title, fontsize=20, weight="bold")
    scope = "historical synthetic mechanism-control settings" if old else "the scored synthetic cohort"
    fig.text(.10, .914, f"Highest and lowest raw paired material-L2 change in {scope} at {budget} s.",
             color=MUTED, fontsize=9.1)
    names = ("dual", "srap", "reset", "full_ap", "graph") if old else ("dual", "srap")
    handles = [Line2D([], [], color=COLORS[n], label=METHODS[n],
                      linewidth=1.6 if n in ("dual", "srap") else 1.) for n in names]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.096, .884), ncol=len(names),
               frameon=False, fontsize=8.7)
    for i, (gain, c) in enumerate(examples):
        by_row = axes[i]
        for ax in by_row:
            ax.set_xlim(0, budget + .04)
            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.grid(color=GRID, linewidth=.7)
            if i == len(examples) - 1:
                ax.set_xlabel("Wall time (s)")
        info = f"{c['case'].replace('__', ' · ')}  /  paired L2 change {gain:+.1f}%"
        by_row[0].text(0, 1.15, info, transform=by_row[0].transAxes, fontsize=9.6,
                       weight="bold", clip_on=False)
        for j, (value, label) in enumerate((
                ("rel_L2", "Relative material L2"), ("residual", "Physical data residual"))):
            ax = by_row[j]
            ax.set_title(label, loc="left", fontsize=9.1, color=MUTED, pad=5)
            for name in names:
                samples = sorted(rows.get((c["case"], name), []), key=lambda r: float(r["seconds"]))
                if samples:
                    ax.plot([float(r["seconds"]) for r in samples],
                            [parse_number(r.get(value)) for r in samples],
                            color=COLORS[name], linewidth=1.8 if name in ("dual", "srap") else 1.1,
                            marker=".", markersize=2.4 if name in ("dual", "srap") else 0)
                scored = ledger.get((c["case"], name, budget, "raw"))
                if scored:
                    clock = parse_number(scored.get("available_wall"))
                    level = parse_number(scored.get(value))
                    if np.isfinite(clock) and np.isfinite(level):
                        ax.scatter([clock], [level], s=26, linewidth=.8, edgecolor=COLORS[name],
                                   facecolor="white", zorder=5)
        ax = by_row[2]
        ax.set_title("SRAP intervention ratios", loc="left", fontsize=9.1, color=MUTED, pad=5)
        events = route_events(c["case"], "srap", budget)
        ratios = (
            ("aperture_leak_vs_data", "RX leakage / data", "#20688c"),
            ("commutator_vs_routed_field", "Graph defect / update", "#904c85"),
            ("reinjection_change_error", "Reinjection error / update", "#bd511c"),
        )
        if events:
            for key, label, color in ratios:
                actual = [(float(e["wall"]), float(e[key])) for e in events
                          if e.get(key) is not None and np.isfinite(e[key]) and e[key] > 0]
                if actual:
                    times, values = zip(*actual, strict=True)
                    ax.semilogy(times, values, color=color, linewidth=1.1, marker=".", markersize=3,
                                label=label)
            ax.legend(frameon=False, fontsize=6.6, loc="best", labelspacing=.45)
        else:
            ax.text(.5, .5, "No recorded SRAP intervention", transform=ax.transAxes,
                    color=MUTED, ha="center", va="center", fontsize=8)
    fig.text(.10, .035, "Open markers identify the raw candidate selected by the physical objective. "
             "Intervention ratios measure update geometry.", color=MUTED, fontsize=8.5)
    save_figure(fig, folder, f"{cohort}_trajectory_examples_{budget}s")


def followup_atlas(folder, cases, ledger, budget):
    from data_followup import ARRAYS
    groups = defaultdict(dict)
    for case in cases:
        groups[case["object_id"]][case["array"]] = case
    if len(groups) != 8 or any(set(g) != set(ARRAYS) for g in groups.values()):
        raise RuntimeError("The follow-up atlas requires eight complete four-array objects.")
    with PdfPages(folder / f"followup_materials_{budget}s.pdf") as pdf:
        for obj, group in sorted(groups.items()):
            cs = [group[name] for name in ARRAYS]
            with np.load(run.DATA / "truth" / f"{cs[0]['case']}.npz", allow_pickle=False) as z:
                title = {"continuous_high_low": "Continuous high/low material",
                         "thin_bridge_multibody": "Thin-bridge multibody"}[str(z["family"])]
                truth = z["m"].reshape(cs[0]["side"], cs[0]["side"])
            fig = material_page(cs, ledger, ("dual", "srap"), budget, f"{obj} · {title}",
                                f"{budget} s  /  centered aperture and flipped illumination tests",
                                shared_truth=truth)
            pdf.savefig(fig)
            save_figure(fig, folder / f"followup_atlas_{budget}s", obj, vector=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--section", nargs="+", default=["schematics"], choices=(
        "schematics", "confirmation", "legacy", "measured_extra", "development", "followup",
        "fresnel_lazy", "legacy_cross_host",
    ))
    parser.add_argument("--budget", type=int, choices=(15, 30), default=30)
    parser.add_argument("--output", type=Path, default=FIGURES)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for section in args.section:
        if section == "schematics":
            draw_geometry(args.output)
            draw_geometry(args.output, followup=True)
            draw_fresnel_masks(args.output)
            draw_receiver_spectra(args.output)
            continue
        if section == "fresnel_lazy":
            lazy_fresnel_atlas(args.output, args.budget)
            continue
        if section == "legacy_cross_host":
            legacy_cross_host(args.output, args.budget)
            continue
        arms, painter = {
            "confirmation": (("dual", "srap"), confirmation_atlas),
            "legacy": (("dual", "srap", "u11"), legacy_atlas),
            "measured_extra": (("dual", "srap"), measured_atlas),
            "development": (("dual", "srap", "dual_cap", "srap_cap"), development_atlas),
            "followup": (("dual", "srap"), followup_atlas),
        }[section]
        cs, ledger = prepare_scored(section, arms, args.budget)
        painter(args.output, cs, ledger, args.budget)
        if section in ("confirmation", "legacy"):
            trajectory_plot(args.output, section, cs, ledger, args.budget)
        if section == "confirmation":
            line_cap_examples(args.output, cs, ledger, args.budget)
        print(f"Rendered {section}: {len(cs)} cases at {args.budget}s", flush=True)


if __name__ == "__main__":
    main()
