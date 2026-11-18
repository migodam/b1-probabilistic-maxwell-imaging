"""Summarize the R4 weak-direction study and render the two research figures."""
import argparse
import csv
from collections import defaultdict
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from weak_r4 import DEFAULT_ARCHIVE, HERE, R4, archived_material_score, read_npz


def rescore(result, materials, archive):
    sys.path.insert(0, str(archive.resolve()))
    r4 = archive / R4
    rows = json.loads((r4/"data/manifests/evaluation.json").read_text())["records"]
    records = {r["obs_id"]: r for r in rows}
    for c in result["cases"]:
        pid, obs_id = c["parent_id"], c["obs_id"]
        row = records[obs_id]
        truth = read_npz(r4/"data"/row["truth_path"])
        old_gaps = [q["total_over_noise"] for q in c["source_decomposition"] if q["library"]=="repeat"]
        c["old_pair_gap_over_noise"] = float(np.sqrt(np.mean(np.asarray(old_gaps)**2)))
        c["old_pair_relative_data_gap"] = c["old_pair_gap_over_noise"] * row["noise"]
        for selectors in c["selection"].values():
            for label, ids in selectors.items():
                selectors[label] = list(map(int, ids))
        legacy = {
            "allscales_discrepancy": materials[pid+"_anchor"],
            "finite_plus": materials[pid+"_plus"], "finite_minus": materials[pid+"_minus"],
        }
        for name, suffix in (("fulltrf", ""), ("voxel_discrepancy", "_discrepancy"),
                             ("weakreg1e-5", "_weakreg1e-5")):
            legacy[name] = read_npz(r4/"anchors"/f"{obs_id}{suffix}.npz")["m0"]
        for past in c["legacy"]:
            past.update(archived_material_score(legacy[past["name"]], truth))
        c["old_cold_reg"].update(archived_material_score(materials[pid+"_old_cold_reg"], truth))
        for a in c["arms"]:
            a["source_ids"] = list(map(int, a["source_ids"]))
            prefix = f'{pid}_{a["library"]}_{a["selector"]}'
            for k, sol in a["solutions"].items():
                sol.update(archived_material_score(materials[prefix+"_"+k], truth))
    result["config"]["material_scores"] = "Archived GN–NN convention: RMSE and amplitude-normalized RMSE, 256-bin Otsu-IoU against alpha >= .5"


def save_cases(cases, figures):
    fig, axes = plt.subplots(len(cases), 5, figsize=(11, 2.17*len(cases)),
                             layout="constrained", squeeze=False)
    for row, (c, m) in enumerate(cases):
        names = ("Archived material", "Old GN", "Alternative −", "Alternative +", "Difference")
        images = [m[k].reshape((24,24)) for k in ("truth", "anchor", "minus", "plus")]
        lower, upper = min(a.min() for a in images), max(a.max() for a in images)
        for col, (ax, img) in enumerate(zip(axes[row], images)):
            ax.imshow(img, cmap="inferno", interpolation="nearest", vmin=lower, vmax=upper)
            ax.set(xticks=[], yticks=[])
            if row == 0:
                ax.set_title(names[col])
        axes[row, 0].set_ylabel(f'{c["parent_id"]}\nContrast {c["amplitude"]:.2f}')
        axes[row, 0].set_xlabel(f"Material scale {lower:.2f}–{upper:.2f}")
        difference = images[3]-images[2]
        maxdiff = np.abs(difference).max()
        ax = axes[row, 4]
        ax.imshow(difference, cmap="coolwarm", interpolation="nearest", vmin=-maxdiff, vmax=maxdiff)
        ax.set(xticks=[], yticks=[], xlabel=f'RMSE {c["distance_rmse"]:.3f}')
        base, plus, minus = c["residuals"]
        for col, value in zip((1, 2, 3), (base, minus, plus)):
            axes[row, col].set_xlabel(f'Old residual {value*1e3:.3f} × 10⁻³')
        if row == 0:
            ax.set_title(names[4])
    for ext in ("png","pdf"):
        fig.savefig(figures/f"weak_candidates.{ext}", dpi=220)
    plt.close(fig)


def median(xs):
    return float(np.median(xs))


def find_arm(case, library, selector):
    return next(a for a in case["arms"] if a["library"] == library and a["selector"] == selector)


def save_sources(cases, figures):
    fig, (ax, decomp) = plt.subplots(1, 2, figsize=(11, 4.2), layout="constrained",
                                    gridspec_kw={"width_ratios":[1.33,1]})
    x = np.arange(len(cases))
    chosen = "candidate_gap"
    definitions = [
        ("repeat", "equally_spaced", "gap_over_noise", "Old positions repeated", "#766850", "x"),
        ("same_radius", chosen, "gap_over_noise", "New angles, same radius", "#414870", "o"),
        ("near_field", chosen, "field_match_gap_over_noise", "Near field, equal incident norm", "#15837c", "^"),
        ("near_field", chosen, "gap_over_noise", "Near field, equal drive", "#bd5352", "s"),
    ]
    for group, selector, field, name, color, marker in definitions:
        vals = [find_arm(c, group, selector)[field] for c in cases]
        ax.plot(x, vals, color=color, marker=marker, lw=1.2, label=name)
    ax.axhline(1, color=".45", lw=1, ls=":")
    ax.set(xticks=x, xticklabels=[c["parent_id"].split("_p")[-1] for c in cases],
           xlabel="Archived TEST parent", ylabel="Four-source branch data gap / new noise norm",
           title="Source geometry and incident intensity")
    ax.set_ylim(top=max(find_arm(c, "near_field", chosen)["gap_over_noise"] for c in cases)*1.25)
    ax.legend(fontsize=7.4, frameon=False, loc="upper center", ncols=2)
    c = cases[0]
    f = find_arm(c, "same_radius", chosen)
    n = find_arm(c, "near_field", chosen)
    select = [f, n, n]
    legends = ("Same radius", "Near field\nequal drive", "Near field\nequal incident norm")
    keys = ("old_combo_over_joint_new_noise", "perpendicular_over_joint_new_noise", "gap_over_noise")
    values = [[z[k] for k in keys] for z in select]
    near_map = {q["i"]: q for q in c["source_decomposition"]}
    values[-1] = [
        float(np.sqrt(np.mean([near_map[i]["field_match_old_combo_over_noise"]**2 for i in n["source_ids"]]))),
        n["field_match_perpendicular_over_joint_new_noise"], n["field_match_gap_over_noise"],
    ]
    colors = ("#a8a0c4", "#15837c", "#414870")
    for index, name in enumerate(("Old source span", "Outside old source span", "Total gap")):
        decomp.bar(np.arange(3)+.23*(index-1), np.asarray(values)[:,index], .20, color=colors[index],
                   label=name)
    decomp.axhline(1, color=".45", lw=1, ls=":")
    decomp.legend(fontsize=8, frameon=False)
    decomp.set(xticks=np.arange(3), xticklabels=legends,
               ylabel="Finite VIE gap / new noise norm",
               title=f'{c["parent_id"]} projected source fields')
    for ext in ("png","pdf"):
        fig.savefig(figures/f"weak_source_geometry.{ext}", dpi=220)
    plt.close(fig)


def save_summary(cases, path):
    groups = defaultdict(list)
    for c in cases:
        groups["archived", "allscales_discrepancy", "legacy"].append((c["legacy"][0], None))
        groups["old_only", "cold_reg", "cold_reg"].append((c["old_cold_reg"], None))
        for a in c["arms"]:
            for k, sol in a["solutions"].items():
                groups[a["library"], a["selector"], k].append((sol, a))
    names = ("geometry", "selector", "inversion", "n", "median_rmse", "median_nrmse", "median_iou",
             "median_noise_multiple", "discrepancy_n", "median_VIE_lu", "median_complex_rhs",
             "median_seconds", "median_gap_over_noise", "median_field_matched_gap_over_noise")
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=names)
        writer.writeheader()
        for (geometry, selector, inversion), rows in groups.items():
            solutions, arms = zip(*rows)
            row = dict(geometry=geometry, selector=selector, inversion=inversion, n=len(solutions))
            for k in ("rmse", "nrmse", "iou", "noise_multiple"):
                row[f"median_{k}"] = median([s[k] for s in solutions])
            row["discrepancy_n"] = sum(s["noise_multiple"] <= 1.10 for s in solutions)
            if "counts" in solutions[0]:
                row["median_VIE_lu"] = median([s["counts"].get("LU",0) for s in solutions])
                row["median_complex_rhs"] = median([s["counts"].get("actual_complex_RHS",0) for s in solutions])
                row["median_seconds"] = median([s["seconds"] for s in solutions])
            if arms[0] is not None:
                row["median_gap_over_noise"] = median([a["gap_over_noise"] for a in arms])
                row["median_field_matched_gap_over_noise"] = median([a["field_match_gap_over_noise"] for a in arms])
            writer.writerow(row)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", type=Path, default=HERE/"results"/"weak_r4")
    ap.add_argument("--figures", type=Path, default=HERE/"figures")
    ap.add_argument("--rescore-archive", type=Path, nargs="?", const=DEFAULT_ARCHIVE)
    ap.add_argument("--summary-only", action="store_true")
    args = ap.parse_args()
    report = args.results.with_suffix(".json")
    result = json.loads(report.read_text())
    materials = read_npz(args.results.with_suffix(".npz"))
    if args.rescore_archive:
        rescore(result, materials, args.rescore_archive)
        report.write_text(json.dumps(result, indent=2, ensure_ascii=False, default=float)+"\n")
    cases = result["cases"]
    if not args.summary_only:
        args.figures.mkdir(parents=True, exist_ok=True)
        save_cases([(c, {k: materials[f'{c["parent_id"]}_{k}'] for k in ("truth", "anchor", "minus", "plus")})
                    for i in (0, 1, 5) for c in cases[i:i+1]], args.figures)
        save_sources(cases, args.figures)
    budget = result["config"]["supplemental_sources"]
    summary = "weak_summary.csv" if budget == 4 else f"weak_summary_{budget}source.csv"
    save_summary(cases, args.results.with_name(summary))


if __name__ == "__main__":
    main()
