"""Offline scientific figures. Run only after parent-authorized scoring.

Reads offline score products and frozen statistical results. No models, online
decisions, bootstrap recomputation, training, physical solves, or truth sources.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.dont_write_bytecode = True
OUT = Path(__file__).resolve().parents[1]
METHODS = ("physics_probit", "strong_cnn", "k32", "deterministic_matched", "distilled_w128")
LABELS = ("Probit", "StrongCNN", "K32", "Det matched", "Distilled")
COLORS = ("#777777", "#366caa", "#e39b29", "#498c68", "#b95872")
PHASES = ("Outside ROI", "Host plateau", "Transition", "Inclusion plateau")
COMPARISONS = ("deterministic_vs_k32", "distilled_vs_k32", "distilled_vs_strong", "retention_80pct")
COMPARISON_LABELS = ("Det vs K32", "Distilled vs K32", "Distilled vs Strong", "80% retention contrast")


def read_json(path):
    return json.loads(Path(path).read_text())


def csv_rows(path, rows):
    if not rows:
        return
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save(figure, output, name):
    figure.savefig(output / f"{name}.png", dpi=180, bbox_inches="tight")
    figure.savefig(output / f"{name}.pdf", bbox_inches="tight")
    plt.close(figure)


def interval_line(ax, x, y, interval, color, label=None, offset=0, marker="o"):
    lo, hi = interval
    # Percentile intervals need not contain the point; draw segment explicitly.
    ax.plot([lo, hi], [y + offset, y + offset], color=color, lw=2, label=label)
    ax.plot(x, y + offset, marker=marker, color=color, ms=6)


def material_gallery(pack, indices, output, group, scales):
    # One frozen seed is shown; galleries never substitute for seed-mean risk.
    prediction = pack["prediction_raw"][indices, 0, :5]
    truth = pack["truth_chi"][indices]
    images = np.concatenate((truth[:, None], prediction), axis=1)
    fig, axes = plt.subplots(4, 6, figsize=(15.8, 10.3), constrained_layout=True)
    titles = ("Truth",) + LABELS
    for case in range(2):
        for component, component_name in enumerate(("Real", "Imag")):
            for column, title in enumerate(titles):
                ax = axes[2 * case + component, column]
                values = images[case, column].real if component == 0 else images[case, column].imag
                image = ax.imshow(values.reshape(24, 24), origin="lower", cmap="viridis",
                                  vmin=scales[component][0], vmax=scales[component][1])
                if case == 0 and component == 0:
                    ax.set_title(title)
                if column == 0:
                    ax.set_ylabel(f"{group}{case} · {component_name}")
                ax.set_xticks([])
                ax.set_yticks([])
            fig.colorbar(image, ax=axes[2 * case + component, :], shrink=0.78, pad=0.01)
    fig.suptitle(f"Fixed {group} parent indices 0, 1 · raw predictions · seed 2026102601", fontsize=14)
    save(fig, output, f"FIXED_PARENT_MATERIAL_{group}")
    np.savez_compressed(output / f"FIXED_PARENT_MATERIAL_{group}_SOURCE.npz",
        parent_ids=pack["parent_ids"][indices], source_indices=indices, images=images,
        labels=np.array(titles), seed=np.array(2026102601), shared_component_scales=np.array(scales))


def risk_and_effects(statistics, output):
    primary = statistics["raw_primary"]
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2), constrained_layout=True)
    rows = []
    for j, (method, label, color) in enumerate(zip(METHODS, LABELS, COLORS)):
        point = primary["method_means"][method]
        bounds = primary["method_mean_bounds_descriptive"][method]["ci_two_sided_98_75"]
        axes[0].plot([bounds[0], bounds[1]], [j, j], color=color, lw=2)
        axes[0].plot(point, j, "o", color=color)
        rows.append({"kind": "pooled_method_risk", "scope": "pooled_raw_seed_mean", "label": method,
                     "estimate": point, "ci_low": bounds[0], "ci_high": bounds[1]})
    axes[0].set_yticks(range(5), LABELS)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Full-grid complex MSE")
    axes[0].set_title("Pooled risk · method intervals descriptive")
    axes[0].ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    for scope, color, offset in (("pooled", "#222222", -0.2), ("A", "#3976ae", 0), ("B", "#dc9042", 0.2)):
        source = primary if scope == "pooled" else statistics["batches_descriptive"][scope]["raw"]
        for j, name in enumerate(COMPARISONS[:3]):
            value = source["comparisons"][name]
            bounds = value["ci_two_sided_98_75"]
            interval_line(axes[1], value["estimate"], j, bounds, color,
                          label=scope if j == 0 else None, offset=offset)
            rows.append({"kind": "paired_effect", "scope": scope, "label": name,
                         "estimate": value["estimate"], "ci_low": bounds[0], "ci_high": bounds[1]})
    axes[1].axvline(0, color="#aaaaaa", lw=1)
    axes[1].set_yticks(range(3), COMPARISON_LABELS[:3])
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Paired MSE improvement · positive favors candidate")
    axes[1].ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    axes[1].set_title("98.75% CIs · batches are descriptive replication")
    axes[1].legend(frameon=False, loc="best")
    save(fig, output, "RISK_AND_PAIRED_BATCH_EFFECTS")
    csv_rows(output / "RISK_AND_PAIRED_BATCH_EFFECTS_SOURCE.csv", rows)


def phase_figure(pack, compact, output):
    counts = compact["phase_cell_counts"]
    raw_sums = compact["raw_phase_squared_error_sums"][:, :, :5].mean(axis=1)
    roi_count = int(pack["roi_mask"].sum())
    batches = pack["batch_ids"].astype(str)
    fig, axes = plt.subplots(1, 3, figsize=(17.2, 5.0), constrained_layout=True)
    phase_colors = ("#558aaf", "#e1ae48", "#b85d6b")
    rows = []
    bottom = np.zeros(2)
    for p in range(1, 4):
        share = [float(counts[batches == batch, p].mean() / roi_count) for batch in ("A", "B")]
        axes[0].bar([0, 1], share, bottom=bottom, color=phase_colors[p - 1], label=PHASES[p])
        bottom += share
        for batch, value in zip(("A", "B"), share):
            rows.append({"kind": "anchor_roi_share", "scope": batch, "method": "anchor", "phase": PHASES[p], "value": value})
    axes[0].set_xticks([0, 1], ["A", "B"])
    axes[0].set_ylim(0, 1)
    axes[0].set_ylabel("Share of ROI cells")
    axes[0].set_title("Anchor-defined phase proportions")
    axes[0].legend(frameon=False, fontsize=8)
    for j, (method, label, color) in enumerate(zip(METHODS, LABELS, COLORS)):
        values = []
        for p in range(1, 4):
            cells = int(counts[:, p].sum())
            value = float(raw_sums[:, j, p].sum() / cells) if cells else np.nan
            values.append(value)
            rows.append({"kind": "phase_cell_mse", "scope": "pooled_raw_seed_mean", "method": method, "phase": PHASES[p], "value": value})
        axes[1].plot(range(3), values, "o-", color=color, label=label)
    axes[1].set_xticks(range(3), PHASES[1:], rotation=15)
    axes[1].set_ylabel("Complex MSE per phase cell")
    axes[1].set_title("Conditional phase risk · descriptive")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    for k, j in enumerate((2, 3, 4)):
        values = (raw_sums[:, 1] - raw_sums[:, j]).mean(axis=0) / 576
        axes[2].bar(np.arange(4) + (k - 1) * 0.25, values, width=0.24, color=COLORS[j], label=LABELS[j])
        for p, value in enumerate(values):
            rows.append({"kind": "full_grid_gain_contribution_vs_strong", "scope": "pooled_raw_seed_mean", "method": METHODS[j], "phase": PHASES[p], "value": float(value)})
    axes[2].axhline(0, color="#aaaaaa", lw=1)
    axes[2].set_xticks(range(4), PHASES, rotation=20)
    axes[2].set_ylabel("Contribution to full-grid MSE gain vs Strong")
    axes[2].set_title("Phase contributions sum to total gain")
    axes[2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axes[2].legend(frameon=False, fontsize=8)
    save(fig, output, "ANCHOR_PHASE_RISK_AND_GAIN")
    csv_rows(output / "ANCHOR_PHASE_RISK_AND_GAIN_SOURCE.csv", rows)


def retention_figure(statistics, output):
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.4), constrained_layout=True)
    rows = []
    for scope, color, offset in (("raw_primary", "#333333", -0.12), ("clipped_secondary", "#6d8cba", 0.12)):
        source = statistics[scope]
        for j, name in enumerate(COMPARISONS):
            value = source["comparisons"][name]
            bounds = value["ci_two_sided_98_75"]
            interval_line(axes[0], value["estimate"], j, bounds, color,
                          label="Raw primary" if scope == "raw_primary" and j == 0 else "Clipped secondary" if j == 0 else None, offset=offset)
            if name == "retention_80pct":
                axes[0].plot(value["lower_one_sided_98_75"], j + offset, marker=">", color=color, ms=8)
            rows.append({"kind": "paired_contrast", "scope": scope, "label": name, "estimate": value["estimate"],
                         "ci_low": bounds[0], "ci_high": bounds[1], "one_sided_lower": value["lower_one_sided_98_75"],
                         "status": "frozen contrast"})
        ratio = source["retention_ratio"]
        j = 0 if scope == "raw_primary" else 1
        if ratio["status"] == "ESTIMABLE":
            bounds = ratio["bounds"]["ci_two_sided_98_75"]
            interval_line(axes[1], ratio["estimate"], j, bounds, color)
            rows.append({"kind": "retention_ratio_secondary", "scope": scope, "label": "retained_fraction",
                         "estimate": ratio["estimate"], "ci_low": bounds[0], "ci_high": bounds[1],
                         "one_sided_lower": ratio["bounds"]["lower_one_sided_98_75"], "status": ratio["status"]})
        else:
            axes[1].text(0.02, 0.73 - j * 0.4, f"{scope}: ratio unavailable\n{ratio['status']}",
                         transform=axes[1].transAxes, fontsize=8)
    axes[0].set_yticks(range(4), COMPARISON_LABELS)
    axes[0].invert_yaxis()
    axes[0].axvline(0, color="#aaaaaa", lw=1)
    axes[0].set_xlabel("Paired MSE improvement / retention contrast")
    axes[0].ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    axes[0].set_title("98.75% intervals · triangle marks retention lower bound")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].axvline(0.8, color="#aa5555", ls="--", lw=1, label="80% threshold")
    axes[1].set_yticks([0, 1], ["Raw", "Clipped"])
    axes[1].set_ylim(1.6, -0.6)
    axes[1].set_xlabel("Retained incremental gain fraction")
    axes[1].set_title("Ratio interpretation requires a stable positive denominator")
    axes[1].legend(frameon=False)
    save(fig, output, "RETENTION_AND_PAIRED_CONTRASTS")
    csv_rows(output / "RETENTION_AND_PAIRED_CONTRASTS_SOURCE.csv", rows)


def seedmean_ensemble_figure(compact, output):
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.6), constrained_layout=True)
    rows = []
    for ax, kind in zip(axes, ("raw", "clipped")):
        seed_mean = compact[f"{kind}_mse_full"][:, :, :5].mean(axis=(0, 1))
        ensemble = compact[f"ensemble_{kind}_mse_full"][:, 0, :5].mean(axis=0)
        ax.bar(np.arange(5) - 0.18, seed_mean, width=0.35, color="#6087b1", label="Mean of seed losses")
        ax.bar(np.arange(5) + 0.18, ensemble, width=0.35, color="#d49d55", label="Loss of mean prediction")
        ax.set_xticks(range(5), LABELS, rotation=20)
        ax.set_ylabel("Full-grid complex MSE")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.set_title(f"{kind.capitalize()} · separate estimands")
        ax.legend(frameon=False, fontsize=8)
        for method, value, ens in zip(METHODS, seed_mean, ensemble):
            rows.extend([{"loss": kind, "method": method, "estimator": "mean_seed_loss", "mse": float(value)},
                         {"loss": kind, "method": method, "estimator": "ensemble_prediction_loss", "mse": float(ens)}])
    save(fig, output, "RAW_CLIPPED_SEEDMEAN_AND_ENSEMBLE")
    csv_rows(output / "RAW_CLIPPED_SEEDMEAN_AND_ENSEMBLE_SOURCE.csv", rows)


def offline_failure_panels(pack, compact, output):
    loss = compact["raw_mse_full"][:, :, :5].mean(axis=1)
    excess = loss[:, 4] - loss[:, 1]
    selected = np.argsort(-excess, kind="stable")[:2]
    truth = pack["truth_chi"][selected]
    strong = pack["prediction_raw"][selected, :, 1].mean(axis=1)
    distilled = pack["prediction_raw"][selected, :, 4].mean(axis=1)
    images = np.stack((truth, strong, distilled), axis=1)
    scales = [(float(images.real.min()), float(images.real.max())),
              (float(images.imag.min()), float(images.imag.max()))]
    fig, axes = plt.subplots(4, 3, figsize=(8.8, 10.4), constrained_layout=True)
    for case, index in enumerate(selected):
        for component, cname in enumerate(("Real", "Imag")):
            for column, title in enumerate(("Truth", "Strong ensemble", "Distilled ensemble")):
                value = images[case, column].real if component == 0 else images[case, column].imag
                image = axes[2 * case + component, column].imshow(value.reshape(24, 24), origin="lower",
                    cmap="viridis", vmin=scales[component][0], vmax=scales[component][1])
                axes[2 * case + component, column].set_xticks([])
                axes[2 * case + component, column].set_yticks([])
                if case == 0 and component == 0:
                    axes[0, column].set_title(title)
            axes[2 * case + component, 0].set_ylabel(f"Rank {case + 1} · {cname}")
            fig.colorbar(image, ax=axes[2 * case + component, :], shrink=0.75, pad=0.01)
    fig.suptitle("OFFLINE_ERROR_SELECTED\nLargest distilled − Strong raw seed-mean MSE · explanatory only", fontsize=12)
    save(fig, output, "OFFLINE_ERROR_SELECTED_DIAGNOSTICS")
    np.savez_compressed(output / "OFFLINE_ERROR_SELECTED_DIAGNOSTICS_SOURCE.npz",
        selected_indices=selected, parent_ids=pack["parent_ids"][selected], truth=truth,
        strong_ensemble=strong, distilled_ensemble=distilled, seedmean_mse_excess=excess[selected],
        selection=np.array("OFFLINE_ERROR_SELECTED: descending raw seed-mean distilled minus Strong loss"))
    csv_rows(output / "OFFLINE_ERROR_SELECTED_PARENTS.csv", [
        {"rank": j + 1, "parent": str(pack["parent_ids"][i]), "batch": str(pack["batch_ids"][i]),
         "global_index": int(i), "distilled_minus_strong_raw_seedmean_mse": float(excess[i]),
         "selection": "OFFLINE_ERROR_SELECTED"} for j, i in enumerate(selected)])


def optional_cost_figure(root, output):
    path = root / "audit/timing/CONFIRMATION_TIMING.json"
    if not path.exists():
        return {"status": "NOT_RUN_BENCHMARK_UNAVAILABLE"}
    benchmark = read_json(path)
    if benchmark.get("status") != "COMPLETE":
        return {"status": "NOT_RUN_BENCHMARK_INCOMPLETE"}
    rows = []
    categories = (("steady_provider_only", "Provider", "#d29b45"),
                  ("steady_full_predict", "Full prediction", "#4c80ab"),
                  ("load_plus_first_inference", "Load + first prediction", "#9d6c96"))
    batches = sorted(set(row["batch_size"] for row in benchmark["records"]))
    fig, axes = plt.subplots(1, len(batches), figsize=(6.4 * len(batches), 4.9), constrained_layout=True, squeeze=False)
    for ax, batch in zip(axes[0], batches):
        records = [row for row in benchmark["records"] if row["batch_size"] == batch]
        for k, (key, label, color) in enumerate(categories):
            for j, row in enumerate(records):
                timing = row.get(key, {})
                if timing.get("status") == "NOT_APPLICABLE" or "p50_seconds" not in timing:
                    continue
                point = timing["p50_seconds"] * 1000
                low, high = timing["p10_seconds"] * 1000, timing["p90_seconds"] * 1000
                ax.plot([j + (k - 1) * 0.18] * 2, [low, high], color=color, lw=2)
                ax.plot(j + (k - 1) * 0.18, point, "o", color=color, label=label if j == (1 if key == "steady_provider_only" else 0) else None)
                rows.append({"batch_size": batch, "arm": row["arm"], "metric": key,
                             "p10_seconds": timing["p10_seconds"], "p50_seconds": timing["p50_seconds"],
                             "p90_seconds": timing["p90_seconds"], "scope": "inference excludes physical acquisition and anchor construction"})
        ax.set_xticks(range(len(records)), [row["arm"] for row in records], rotation=20)
        ax.set_yscale("log")
        ax.set_ylabel("Milliseconds per whole batch · median, p10–p90")
        ax.set_title(f"CPU1 inference · batch {batch}")
        ax.legend(frameon=False, fontsize=8)
    save(fig, output, "CPU_PROVIDER_AND_PREDICTION_COST")
    csv_rows(output / "CPU_PROVIDER_AND_PREDICTION_COST_SOURCE.csv", rows)
    # No subtraction or extrapolation turns these timings into end-to-end costs.
    return {"status": "PLOTTED", "source": str(path), "end_to_end": "NOT_PLOTTED_NO_ANCHOR_INCLUSIVE_BENCHMARK_IN_THIS_SCHEMA"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=OUT)
    parser.add_argument("--statistics", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    score_root = root / "raw/offline_scores"
    summary = read_json(score_root / "SCORING_SUMMARY.json")
    access = read_json(score_root / "OFFLINE_ACCESS_STATUS.json")
    if summary.get("status") != "OFFLINE_SCORING_COMPLETE_NO_STATISTICAL_DECISION" or access.get("status") != "OFFLINE_SCORING_COMPLETE":
        raise RuntimeError("Plotting requires completed parent-authorized offline scoring")
    statistics_path = args.statistics or root / "RESULTS_AND_STATISTICS.json"
    statistics = read_json(statistics_path)
    statistics = statistics.get("statistics", statistics)
    if statistics.get("status") != "CONFIRMATORY_ANALYSIS" or statistics.get("parents") != 512:
        raise ValueError("Require frozen pooled 512-parent confirmation statistics")
    output = args.output or root / "figures/confirmation"
    if output.exists():
        raise FileExistsError("Figure output already exists; no implicit overwrite")
    with np.load(score_root / "OFFLINE_SCOREPACK.npz", allow_pickle=False) as z:
        pack = {key: np.asarray(z[key]) for key in z.files}
    with np.load(score_root / "PARENT_METRICS.npz", allow_pickle=False) as z:
        compact = {key: np.asarray(z[key]) for key in z.files}
    if list(pack["methods"][:5].astype(str)) != list(METHODS) or list(compact["methods"][:5].astype(str)) != list(METHODS):
        raise ValueError("Canonical method order mismatch")
    if not np.array_equal(pack["parent_ids"], compact["parent_ids"]):
        raise ValueError("Scorepack/metric parent order mismatch")
    batches = pack["batch_ids"].astype(str)
    fixed = {group: np.flatnonzero(batches == group)[:2] for group in ("A", "B")}
    if any(len(np.flatnonzero(batches == group)) != 256 for group in ("A", "B")):
        raise ValueError("Expected exactly256 parents per batch")
    selected = np.concatenate(list(fixed.values()))
    images = np.concatenate((pack["truth_chi"][selected, None], pack["prediction_raw"][selected, 0, :5]), axis=1)
    scales = [(float(images.real.min()), float(images.real.max())), (float(images.imag.min()), float(images.imag.max()))]
    output.mkdir(parents=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "savefig.facecolor": "white"})
    for group in ("A", "B"):
        material_gallery(pack, fixed[group], output, group, scales)
    risk_and_effects(statistics, output)
    phase_figure(pack, compact, output)
    retention_figure(statistics, output)
    seedmean_ensemble_figure(compact, output)
    offline_failure_panels(pack, compact, output)
    cost = optional_cost_figure(root, output)
    captions = """# Confirmation figure captions

**FIXED_PARENT_MATERIAL_A / B.** Fixed within-batch parent indices 0 and 1, chosen before looking at losses. Truth and raw predictions at frozen training seed 2026102601; real and imaginary parts have separate scales shared across both batches, all four selected parents and all methods. These images do not display the seed-mean loss estimator or an ensemble.

**RISK_AND_PAIRED_BATCH_EFFECTS.** Left: pooled raw full-grid complex MSE, averaging each parent's three frozen seed losses first; method intervals are descriptive. Right: paired improvements, positive favoring the named candidate, with Bonferroni four-comparison 98.75% bootstrap intervals. Pooled raw inference is primary; A/B intervals describe replication and cannot be selected as alternative primary tests.

**ANCHOR_PHASE_RISK_AND_GAIN.** Phase is fixed from the online anchor A q0 / halfwidth, using the original prior's ROI. Conditional phase-cell risk divides summed squared error by summed phase-cell count. Phase contributions divide each parent's phase squared-error sum by576 and average; contributions add to total full-grid gain versus StrongCNN. Transition is an anchor-defined prior saturation diagnostic, not a Maxwell-specific mechanism claim. Empty phase groups have no estimated conditional risk.

**RETENTION_AND_PAIRED_CONTRASTS.** The retention contrast is (Strong−Distilled)−0.8(Strong−K32). Intervals are two-sided98.75%; the triangle separately marks the predeclared one-sided98.75% lower contrast bound. The ratio is secondary and is displayed only with the frozen stable-positive-denominator qualification. Raw is primary, clipped secondary. This figure does not change frozen statistical decisions.

**RAW_CLIPPED_SEEDMEAN_AND_ENSEMBLE.** Mean of seed losses and loss of mean prediction are separate estimands. Raw and clipped losses occupy separate panels. Ensemble values are descriptive and never replace primary seed-mean losses.

**OFFLINE_ERROR_SELECTED_DIAGNOSTICS.** The two largest parent raw seed-mean MSE excesses of distilled versus StrongCNN, selected after offline scoring, with stable parent-order tie breaking. Images show mean predictions across seeds for explanation; the selection statistic is mean seed loss. This is explicitly OFFLINE_ERROR_SELECTED, never an online choice rule or a representative sample; excess may be nonpositive if no parent is worse.

**CPU_PROVIDER_AND_PREDICTION_COST, if present.** Measured CPU1 median and p10–p90 repeated-execution timings per whole batch, with batch1 and64 separated. Provider, full prediction and fresh-action load plus first prediction remain distinct. These timings exclude physical acquisition and anchor construction; filesystem/import caches may be warm. They do not establish anchor-inclusive end-to-end runtime or independent training-seed timing variability. Missing timing data produces no substitute cost chart.
"""
    (output / "FIGURE_CAPTIONS.md").write_text(captions)
    receipt = {"status": "FIGURES_CREATED_OFFLINE_REVIEW_REQUIRED", "statistics_source": str(statistics_path),
               "fixed_parent_ids": {g: pack["parent_ids"][indices].astype(str).tolist() for g, indices in fixed.items()},
               "cost_plot": cost, "source_data": "compact CSV/NPZ beside each figure",
               "scientific_decisions": "parent orchestrator; no re-estimation or tuning here"}
    (output / "FIGURE_RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": receipt["status"], "output": str(output)}))


if __name__ == "__main__":
    main()
