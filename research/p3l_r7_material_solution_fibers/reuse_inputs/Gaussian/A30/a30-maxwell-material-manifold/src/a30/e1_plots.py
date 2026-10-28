"""One E1 outcome figure from saved corrected-chart local witnesses.

Call ``main(root)`` under the existing A30 accounting wrapper. This module does
not execute at import, recompute F, open labels/locked inputs, or make a gate.
Only compensation_corrected is used; the invalid first-pass chart is excluded.
"""
from __future__ import annotations

from pathlib import Path
import csv
import json
import math
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D


FIXED_CASE = "p250065_weak_a1_snr20"
STEM = "fig08_e1_compensation_outcome"
FAMILIES = ("shell", "asymmetric_piecewise", "smooth", "mixed_scale")
COLORS = {"shell": "#087f8c", "asymmetric_piecewise": "#7c3aed",
          "smooth": "#c2410c", "mixed_scale": "#2563eb"}
LABELS = {"shell": "Shell", "asymmetric_piecewise": "Asymmetric",
          "smooth": "Smooth", "mixed_scale": "Mixed scale"}
MARKERS = {"weak": "o", "ring": "s"}
SOURCES = ("results/a30/compensation_corrected/COMPENSATION_PAIRS.csv",
           "results/a30/compensation_corrected/COMPENSATION_ENDPOINTS.npz",
           "results/a30/hull/HULL_RECONSTRUCTIONS.npz")
LIMITATIONS = (
    "Distances are copied from saved full-F endpoint results; no F evaluation or certificate calculation is performed.",
    "Candidate construction is partial (29/32 in this saved artifact); unavailable rows remain in the raw plotted table.",
    "Points are individual weak/ring conditions and local endpoint witnesses, not global compensated minima or parent-level inference.",
    "Only the corrected chart is displayed; the first-pass invalid chart is not used.",
    "Material panels show magnitude on a fixed 0–1.2 scale; real/imaginary increments and any saturation remain in the raw CSV.",
    "The p250065 weak example is fixed by the task and does not establish a family-wide outcome.",
)


def _path(root, relative):
    path = root / relative
    if not path.resolve().is_relative_to(root) or any(
        part.lower() == "locked" or part.lower().startswith(("locked_id", "locked_ood"))
        for part in path.parts
    ):
        raise ValueError("PLOT_SOURCE_OUTSIDE_DECLARED_ROOT_OR_LOCKED")
    return path


def _finite(value):
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return value if math.isfinite(value) else None


def _usable(row):
    return (str(row.get("valid", "")).lower() in {"true", "1", "yes"}
            and row.get("status") == "COMPLETE"
            and _finite(row.get("nominal_F")) is not None
            and _finite(row.get("compensated_F")) is not None)


def _csv(path, rows):
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main(root):
    """Write PNG, SVG, raw plotted CSV and a concise source/caption JSON."""
    started, cpu = time.perf_counter(), time.process_time()
    root = Path(root).resolve()
    with _path(root, SOURCES[0]).open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    usable = [row for row in rows if _usable(row)]
    if not usable:
        raise ValueError("NO_SAVED_ACTUAL_CORRECTED_E1_DISTANCES")
    fixed = next((row for row in usable if row["case_id"] == FIXED_CASE), None)
    if fixed is None:
        raise ValueError("FIXED_P250065_WEAK_PAIR_UNAVAILABLE")
    with np.load(_path(root, SOURCES[1]), allow_pickle=False) as archive:
        endpoint_names = ("nominal_a", "nominal_b", "compensated_a", "compensated_b")
        endpoints = {name: archive[FIXED_CASE + "__" + name].copy() for name in endpoint_names}
    with np.load(_path(root, SOURCES[2]), allow_pickle=False) as archive:
        case_ids = list(map(str, archive["case_id"]))
        background = archive["background"][case_ids.index(FIXED_CASE)].copy()

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlesize": 11, "axes.titleweight": "semibold",
                         "xtick.labelsize": 9, "ytick.labelsize": 9,
                         "svg.fonttype": "none", "figure.facecolor": "white",
                         "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(13.8, 6.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, width_ratios=(2.2, 1., 1.))
    distance_ax = fig.add_subplot(grid[:, 0])
    image_axes = [fig.add_subplot(grid[row, col]) for row in range(2) for col in (1, 2)]
    plot_rows = []
    all_values = []
    for row in rows:
        nominal, compensated = _finite(row.get("nominal_F")), _finite(row.get("compensated_F"))
        displayed = _usable(row)
        plot_rows.append({"record_type": "distance", "case_id": row["case_id"],
                          "parent_id": row["parent_id"], "family": row["family"],
                          "background_id": row["background_id"], "source_status": row["status"],
                          "source_valid": row["valid"], "displayed": displayed,
                          "nominal_full_F_distance": nominal, "compensated_full_F_distance": compensated,
                          "candidate_id": row.get("candidate_id", ""),
                          "nominal_F_reused_from_paid_run": row.get("nominal_F_reused_from_paid_run", ""),
                          "scope": "CORRECTED_CHART_LOCAL_ENDPOINT_WITNESS"})
        if not displayed:
            continue
        if nominal < 0 or compensated < 0:
            raise ValueError("NEGATIVE_SAVED_WHITENED_DISTANCE")
        distance_ax.scatter(nominal, compensated, s=40, marker=MARKERS[row["background_id"]],
                            color=COLORS[row["family"]], edgecolors="white", linewidths=.5, alpha=.88, zorder=3)
        all_values.extend((nominal, compensated))
    positive = all(value > 0 for value in all_values)
    axis_scale = "log" if positive else "symlog"
    if positive:
        distance_ax.set_xscale("log")
        distance_ax.set_yscale("log")
        low = min(1., min(all_values)) / 1.5
    else:
        distance_ax.set_xscale("symlog", linthresh=1.)
        distance_ax.set_yscale("symlog", linthresh=1.)
        low = 0.
    high = max(3., max(all_values)) * 1.25
    distance_ax.plot([low, high], [low, high], color="#94a3b8", lw=1, ls="--", zorder=1)
    distance_ax.axvline(3., color="#334155", ls=":", lw=1.2)
    distance_ax.axhline(1., color="#334155", ls=":", lw=1.2)
    distance_ax.set_xlim(low, high)
    distance_ax.set_ylim(low, high)
    distance_ax.set_xlabel("Nominal full-F distance (whitened)")
    distance_ax.set_ylabel("Compensated full-F distance (whitened)")
    distance_ax.set_title("Actual saved full-F distances\nIndividual weak/ring local witnesses")
    distance_ax.grid(which="major", color="#e2e8f0", lw=.7)
    distance_ax.text(3., .99, " nominal reference 3", rotation=90, ha="right", va="top", fontsize=8,
                     transform=distance_ax.get_xaxis_transform())
    distance_ax.text(.01, 1., "compensated reference 1", ha="left", va="bottom", fontsize=8,
                     transform=distance_ax.get_yaxis_transform())
    family_handles = [Line2D([], [], marker="o", color=COLORS[f], linestyle="none", markersize=6,
                             label=LABELS[f]) for f in FAMILIES]
    background_handles = [Line2D([], [], marker=MARKERS[bg], color="#475569", linestyle="none", markersize=6,
                                label=bg.capitalize()) for bg in ("weak", "ring")]
    legend = distance_ax.legend(handles=family_handles, loc="upper left", frameon=True,
                                framealpha=.93, edgecolor="none", fontsize=8)
    distance_ax.add_artist(legend)
    distance_ax.legend(handles=background_handles, loc="lower right", frameon=True,
                       framealpha=.93, edgecolor="none", fontsize=8)
    distance_ax.annotate("p250065 weak", xy=(float(fixed["nominal_F"]), float(fixed["compensated_F"])),
                         xytext=(7, -11), textcoords="offset points", fontsize=8, color="#334155")
    plot_rows.extend([
        {"record_type": "protocol_reference", "axis": "nominal_x", "distance": 3.,
         "scope": "ORIGINAL_PROTOCOL_SCREENING_REFERENCE"},
        {"record_type": "protocol_reference", "axis": "compensated_y", "distance": 1.,
         "scope": "ORIGINAL_PROTOCOL_SCREENING_REFERENCE"},
    ])

    image_titles = ("Nominal A", "Nominal B", "Compensated A", "Compensated B")
    norm = Normalize(0., 1.2)
    mappable = None
    for ax, name, title in zip(image_axes, endpoint_names, image_titles):
        value = np.asarray(endpoints[name])
        if value.shape != (576,) or background.shape != (576,) or not np.isfinite(value).all():
            raise ValueError("SAVED_FIXED_ENDPOINT_NOT_RENDERABLE:" + name)
        delta = (value - background).reshape(24, 24)
        mappable = ax.imshow(np.abs(delta), origin="lower", cmap="viridis", norm=norm, interpolation="nearest")
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        pair_distance = float(fixed["nominal_F"] if name.startswith("nominal") else fixed["compensated_F"])
        ax.set_xlabel(f"Actual pair distance {pair_distance:.2f}", fontsize=8)
        if name.endswith("_a"):
            ax.set_ylabel("Fixed p250065 weak cavity pair", fontsize=8)
        for pixel_row in range(24):
            for pixel_col in range(24):
                pixel = delta[pixel_row, pixel_col]
                plot_rows.append({"record_type": "material_pixel", "case_id": FIXED_CASE,
                                  "endpoint": name, "pixel_row": pixel_row, "pixel_column": pixel_col,
                                  "delta_chi_real": float(pixel.real), "delta_chi_imag": float(pixel.imag),
                                  "delta_chi_magnitude": float(abs(pixel)), "colorbar_low": 0., "colorbar_high": 1.2,
                                  "clipped_above_colorbar": bool(abs(pixel) > 1.2),
                                  "actual_pair_full_F_distance": pair_distance,
                                  "protected_cavity_pixels": fixed.get("protected_cavity_pixels", "")})
    fig.colorbar(mappable, ax=image_axes, fraction=.034, pad=.025, label="Actual |delta chi| (fixed 0–1.2)")
    fig.suptitle(f"Corrected E1 local compensation witnesses: {len(usable)}/{len(rows)} constructed\n"
                 "Partial candidate construction; first-pass invalid chart not used", fontsize=13)

    output = _path(root, "results/a30/figures")
    output.mkdir(parents=True, exist_ok=True)
    png, svg, raw = output / (STEM + ".png"), output / (STEM + ".svg"), output / "raw_plot_data" / (STEM + ".csv")
    raw.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=200, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    _csv(raw, plot_rows)
    caption = {
        "figure": STEM, "sources": list(SOURCES), "constructed_conditions": len(usable),
        "registered_conditions": len(rows), "fixed_case_id": FIXED_CASE, "axis_scale": axis_scale,
        "caption_zh": f"修正 chart 下实际保存的 full-F nominal 与 compensated pair 距离。候选构造为 {len(usable)}/{len(rows)}，未构造条件在原表保留；颜色为 family、形状为 weak/ring。3/1 虚线沿用原协议筛选参照，不是普遍可识别性常数。每个点是局部可行 endpoint witness，不声称全局最优补偿。右侧是固定 p250065 weak 空腔 pair 补偿前后的实际 |chi-background|，共用 0–1.2 色标。first-pass invalid chart 未使用。",
        "caption_en": f"Actual saved nominal and compensated full-F pair distances from the corrected chart, with partial construction ({len(usable)}/{len(rows)}). Family is encoded by color and weak/ring by marker. References 3/1 come from the original screening protocol. Points are local endpoint witnesses. The fixed p250065 weak cavity pair shares one actual material-increment scale of 0–1.2. The first-pass invalid chart is not used.",
        "implementation_limitations": list(LIMITATIONS),
        "unavailable_case_ids": [r["case_id"] for r in rows if not _usable(r)],
        "raw_plot_csv": str(raw.relative_to(root)), "Maxwell_calls": 0, "training_updates": 0, "locked_reads": 0,
        "timing": {"wall_seconds": time.perf_counter() - started, "process_cpu_seconds": time.process_time() - cpu,
                   "nested_in_existing_job_wrapper": True},
    }
    caption_path = output / (STEM + ".caption.json")
    caption_path.write_text(json.dumps(caption, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    result = {"png": str(png), "svg": str(svg), "raw_plot_csv": str(raw), "caption": str(caption_path),
              "constructed_conditions": len(usable), "registered_conditions": len(rows)}
    print(json.dumps(result, indent=2))
    return result
