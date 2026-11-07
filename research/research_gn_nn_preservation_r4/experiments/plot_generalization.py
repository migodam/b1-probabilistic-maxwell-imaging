"""Post-seal generalization plots from scored scalar tables; never open truth/NPZ.

With no score file, write a NOT_RUN receipt and no scientific figures. Run only
after final TEST scoring. Optional selection figures use VAL scores exclusively.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

STUDY = Path(__file__).resolve().parents[1]
GROUPS = ("new_instance", "unseen_class", "unseen_style")
BOOTSTRAP_SEED = 20261107
BOOTSTRAP_REPETITIONS = 500


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows({k: json.dumps(v) if isinstance(v, (dict, list)) else v
                         for k, v in row.items()} for row in rows)


def number(value):
    if value is None or value == "":
        return None
    value = float(value)
    return value if np.isfinite(value) else None


def flag(value):
    return value is True or str(value).lower() in ("true", "1")


def seed_of(run):
    match = re.search(r"(?:^|_)seed(7|17|27)(?:_|$)", run)
    return int(match.group(1)) if match else None


def family_of(run):
    return re.sub(r"_seed(?:7|17|27)(?=_|$)", "", run)


def short(run):
    return (run.replace("_cascade2", "").replace("_seed", " s")
            .replace("v_physics", "Voxel").replace("j_full", "J-full")
            .replace("j_null", "J-null").replace("g_weak", "G-weak")
            .replace("g_dark", "G-dark").replace("soft_dark", "Soft-dark"))


def read_scores(path, split):
    """Only score JSON/CSV; the JSON sidecar supplies CSV provenance."""
    path = Path(path)
    meta_path = path.with_suffix(".json") if path.suffix.lower() == ".csv" else path
    meta = json.loads(meta_path.read_text())
    if meta.get("schema") != "gn-nn-r4-scored-v1":
        raise ValueError("Require the final R4 score JSON (or its companion CSV)")
    expected_scope = "VAL_MODEL_SELECTION" if split == "val" else "NEW_TEST_AFTER_METHOD_AND_PREDICTION_FREEZE"
    if meta.get("truth_scope") != expected_scope:
        raise ValueError("Score provenance does not match the requested split")
    if path.suffix.lower() == ".csv":
        with path.open(newline="") as handle:
            source = list(csv.DictReader(handle))
    else:
        source = meta["rows"]
    rows = []
    seen = set()
    for raw in source:
        if raw.get("split") != split:
            raise ValueError("Mixed or incorrect score split")
        key = (raw["run"], raw["obs_id"])
        if key in seen:
            raise ValueError("Duplicate run/observation in score table")
        seen.add(key)
        top = raw.get("top1") if isinstance(raw.get("top1"), dict) else None
        anchor = raw.get("anchor") if isinstance(raw.get("anchor"), dict) else {}
        selected_error = number(top.get("NRMSE")) if top else number(raw.get("selected_NRMSE"))
        residual = number(top.get("relative_residual")) if top else number(raw.get("selected_data_relative_residual"))
        row = {k: raw.get(k) for k in ("run", "obs_id", "parent_id", "split", "group", "status")}
        row.update(noise=number(raw.get("noise")), amplitude=number(raw.get("amplitude")),
                   selected_NRMSE=selected_error,
                   anchor_NRMSE=number(anchor.get("NRMSE", raw.get("anchor_NRMSE"))),
                   selected_data_residual=residual,
                   selected=selected_error is not None,
                   not_run=raw.get("status") == "NOT_RUN",
                   quality_joint_success=flag(raw.get("joint_success")),
                   budgeted_joint_success=flag(raw.get("budgeted_joint_success")),
                   total55_joint_success=flag(raw.get("joint_success_total55")))
        gate = .001 if row["noise"] == 0 else .003
        row["data_success"] = bool(row["selected"] and residual is not None and residual <= gate)
        for cost in ("online_total_wall_seconds", "online_total_LU", "online_total_RHS"):
            row[cost] = number(raw.get(cost))
        if row["not_run"] or not row["selected"]:
            for metric in ("data_success", "quality_joint_success", "budgeted_joint_success", "total55_joint_success"):
                row[metric] = False
        rows.append(row)
    return rows, meta


def attach_anchor_diagnostics(rows):
    """Reuse already-scored common-anchor errors; never calculate new GT errors."""
    anchors = {}
    for row in rows:
        value = row["anchor_NRMSE"]
        if value is not None:
            obs = row["obs_id"]
            if obs in anchors and not np.isclose(anchors[obs], value, rtol=1e-8, atol=1e-10):
                raise ValueError("Scores disagree on the common anchor NRMSE")
            anchors[obs] = value
    for row in rows:
        anchor = anchors.get(row["obs_id"])
        row["anchor_NRMSE"] = anchor
        row["diagnostic_anchor_fill"] = not row["selected"] and anchor is not None
        row["diagnostic_NRMSE"] = row["selected_NRMSE"] if row["selected"] else anchor
        row["diagnostic_missing"] = row["diagnostic_NRMSE"] is None


def mean_complete(values):
    return float(np.mean(values)) if values and all(v is not None for v in values) else None


def parent_rows(rows):
    """Collapse correlated clean/noisy observations inside each run and mother."""
    groups = defaultdict(list)
    for row in rows:
        groups[(row["run"], row["parent_id"])].append(row)
    result = []
    for (run, parent), observations in sorted(groups.items()):
        head = observations[0]
        if len({(r["group"], r["amplitude"]) for r in observations}) != 1:
            raise ValueError("A parent's group/amplitude changed across observations")
        out = dict(run=run, parent_id=parent, group=head["group"], amplitude=head["amplitude"],
                   seed=seed_of(run), run_family=family_of(run), observations=len(observations),
                   noise_levels=sorted({r["noise"] for r in observations}),
                   not_run=sum(r["not_run"] for r in observations),
                   selected=sum(r["selected"] for r in observations),
                   anchor_filled=sum(r["diagnostic_anchor_fill"] for r in observations),
                   diagnostic_missing=sum(r["diagnostic_missing"] for r in observations))
        for metric in ("diagnostic_NRMSE", "anchor_NRMSE", "online_total_wall_seconds", "online_total_LU", "online_total_RHS"):
            out[metric] = mean_complete([r[metric] for r in observations])
        for metric in ("data_success", "quality_joint_success", "budgeted_joint_success", "total55_joint_success"):
            out[metric] = float(np.mean([r[metric] for r in observations]))
        out["paired_NRMSE_delta"] = (out["diagnostic_NRMSE"]-out["anchor_NRMSE"]
                                      if out["diagnostic_NRMSE"] is not None and out["anchor_NRMSE"] is not None else None)
        result.append(out)
    return result


def bootstrap(values):
    if not values or any(v is None for v in values):
        return None, None, None
    values = np.asarray(values, float)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    means = values[rng.integers(0, len(values), (BOOTSTRAP_REPETITIONS, len(values)))].mean(1)
    return float(values.mean()), float(np.quantile(means, .025)), float(np.quantile(means, .975))


def summarize(parents, dimension):
    groups = defaultdict(list)
    for row in parents:
        label = row[dimension] if dimension != "all" else "ALL"
        groups[(row["run"], label)].append(row)
    result = []
    for (run, label), mothers in sorted(groups.items(), key=lambda pair: (pair[0][0], str(pair[0][1]))):
        row = dict(run=run, seed=seed_of(run), stratum=label, dimension=dimension,
                   planned_parents=len(mothers), planned_observations=sum(r["observations"] for r in mothers),
                   selected_observations=sum(r["selected"] for r in mothers),
                   not_run_observations=sum(r["not_run"] for r in mothers),
                   anchor_filled_observations=sum(r["anchor_filled"] for r in mothers),
                   unresolved_diagnostic_observations=sum(r["diagnostic_missing"] for r in mothers),
                   parents_with_two_observations=sum(r["observations"] == 2 for r in mothers))
        for metric in ("diagnostic_NRMSE", "anchor_NRMSE", "paired_NRMSE_delta", "data_success", "quality_joint_success", "budgeted_joint_success", "total55_joint_success"):
            estimate, low, high = bootstrap([r[metric] for r in mothers])
            row[metric] = estimate
            row[metric+"_low"] = low
            row[metric+"_high"] = high
        for cost in ("online_total_wall_seconds", "online_total_LU", "online_total_RHS"):
            measured = [r[cost] for r in mothers if r[cost] is not None]
            row[cost] = float(np.mean(measured)) if measured else None
            row[cost+"_measured_parents"] = len(measured)
        result.append(row)
    return result


def pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def error_point(ax, x, row, metric, color, label=None):
    value = row[metric]
    if value is None:
        ax.text(x, .04, "missing", rotation=90, transform=ax.get_xaxis_transform(), fontsize=7)
        return
    low, high = row[metric+"_low"], row[metric+"_high"]
    ax.errorbar(x, value, yerr=[[max(0., value-low)], [max(0., high-value)]],
                fmt="o", color=color, capsize=2, label=label)


def stratum_figure(parents, summary, dimension, labels, runs, output, title):
    plt = pyplot()
    fig, axs = plt.subplots(len(labels), 3, figsize=(max(12, .65*len(runs)+7), 3.8*len(labels)), squeeze=False)
    indexed = {(r["run"], r["stratum"]): r for r in summary}
    colors = plt.get_cmap("tab20")
    for i, label in enumerate(labels):
        for x, run in enumerate(runs):
            row = indexed.get((run, label))
            if row is None:
                for ax in axs[i]:
                    ax.text(x, .04, "no planned rows", rotation=90, transform=ax.get_xaxis_transform(), fontsize=6)
                continue
            mothers = [r for r in parents if r["run"] == run and (dimension == "all" or r[dimension] == label)]
            for mother in mothers:
                value = mother["paired_NRMSE_delta"]
                if value is not None:
                    axs[i, 0].scatter(x, value, s=7, alpha=.18, color=colors(x % 20))
            error_point(axs[i, 0], x, row, "paired_NRMSE_delta", colors(x % 20))
            for offset, metric, color, name in ((-.17, "data_success", "C0", "data"),
                    (0., "quality_joint_success", "C1", "quality joint"),
                    (.17, "total55_joint_success", "C2", "joint / total55 budget")):
                error_point(axs[i, 1], x+offset, row, metric, color, name if x == 0 else None)
            # Both values use all planned mothers; missing top1 is anchor filled.
            error_point(axs[i, 2], x-.10, row, "anchor_NRMSE", "0.5", "common anchor" if x == 0 else None)
            error_point(axs[i, 2], x+.10, row, "diagnostic_NRMSE", "C3", "top1 / anchor-fill diagnostic" if x == 0 else None)
        coverage = [indexed.get((run, label)) for run in runs]
        ticklabels = [short(run)+(f"\np={r['planned_parents']};NR={r['not_run_observations']};fill={r['anchor_filled_observations']}" if r else "\nno plan")
                      for run, r in zip(runs, coverage)]
        for ax in axs[i]:
            ax.set_xticks(range(len(runs)), ticklabels, rotation=65, ha="right", fontsize=7)
            ax.set_xlim(-.6, max(len(runs)-.4, .6))
            ax.grid(axis="y", alpha=.2)
        axs[i, 0].axhline(0, color="k", lw=.7)
        axs[i, 0].set_ylabel(f"{label}\nPaired NRMSE change vs common anchor")
        axs[i, 1].set_ylabel("Fraction of planned observations (parent mean)")
        axs[i, 1].set_ylim(-.04, 1.04)
        axs[i, 2].set_ylabel("NRMSE diagnostic; failed top1 filled by anchor")
        axs[i, 1].legend(fontsize=7)
        axs[i, 2].legend(fontsize=7)
    fig.suptitle(title+"\nClean/noisy averaged within mother; 500 mother-bootstrap draws. NR=NOT_RUN; fill is not a successful output.", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, .94))
    fig.savefig(output, dpi=160)
    plt.close(fig)


def seed_figure(parents, output):
    plt = pyplot()
    families = sorted({r["run_family"] for r in parents if r["seed"] is not None})
    fig, axs = plt.subplots(max(1, len(families)), 2, figsize=(11, max(3.5, 3*len(families))), squeeze=False)
    paired = []
    for i, family in enumerate(families):
        records = [r for r in parents if r["run_family"] == family and r["seed"] is not None]
        by_parent = defaultdict(dict)
        for row in records:
            if row["seed"] in by_parent[row["parent_id"]]:
                raise ValueError("Two run variants collapsed onto one seed family")
            by_parent[row["parent_id"]][row["seed"]] = row
        complete = {parent: rs for parent, rs in by_parent.items() if set(rs) == {7, 17, 27}}
        for parent, records_by_seed in sorted(complete.items()):
            values = [records_by_seed[s]["diagnostic_NRMSE"] for s in (7, 17, 27)]
            if all(v is not None for v in values):
                axs[i, 0].plot([7, 17, 27], values, color="0.5", alpha=.15, lw=.5, marker=".")
            for confirmation in (17, 27):
                a, b = records_by_seed[7], records_by_seed[confirmation]
                paired.append(dict(run_family=family, parent_id=parent, confirmation_seed=confirmation,
                                   paired_NRMSE_change=None if a["diagnostic_NRMSE"] is None or b["diagnostic_NRMSE"] is None else b["diagnostic_NRMSE"]-a["diagnostic_NRMSE"],
                                   paired_joint_change=b["quality_joint_success"]-a["quality_joint_success"]))
        for x, seed in enumerate((7, 17, 27)):
            rs = [r for r in records if r["seed"] == seed]
            values = [r["quality_joint_success"] for r in rs]
            est, low, high = bootstrap(values)
            if est is not None:
                axs[i, 1].errorbar(x, est, yerr=[[max(0., est-low)], [max(0., high-est)]], fmt="o", capsize=3)
            axs[i, 1].text(x, .02, f"p={len(rs)}\nNR={sum(r['not_run'] for r in rs)}", ha="center", fontsize=8)
        axs[i, 0].set(title=f"{short(family)}: matched mothers={len(complete)}", ylabel="NRMSE / anchor-fill diagnostic", xticks=[7, 17, 27], xlabel="Seed; each line is the same mother")
        axs[i, 1].set(title="Seeds reported separately; never independent materials", xticks=[0, 1, 2], xticklabels=["7 primary", "17 confirm", "27 confirm"], ylim=(-.04, 1.04), ylabel="Quality joint success / actual plan")
    if not families:
        axs[0, 0].text(.5, .5, "No seeded runs in score plan", ha="center")
    fig.suptitle("Seed confirmation: clean/noisy averaged inside each mother; incomplete seed families are retained in coverage.csv")
    fig.tight_layout(rect=(0, 0, 1, .96));fig.savefig(output, dpi=150);plt.close(fig)
    return paired


def cost_figure(summary, output):
    plt = pyplot()
    fig, axs = plt.subplots(2, 3, figsize=(15, 8))
    metrics = ("online_total_wall_seconds", "online_total_LU", "online_total_RHS")
    for i, row in enumerate(summary):
        for j, cost in enumerate(metrics):
            x = row[cost]
            if x is None:
                continue
            coverage = f"measured {row[cost+'_measured_parents']}/{row['planned_parents']} mothers"
            for k, quality in enumerate(("diagnostic_NRMSE", "total55_joint_success")):
                if row[quality] is not None:
                    axs[k, j].scatter(x, row[quality], color=f"C{i % 10}", marker="o" if row["seed"] in (None, 7) else "x")
                    axs[k, j].annotate(short(row["run"])+"\n"+coverage, (x, row[quality]), fontsize=6, xytext=(3, 3), textcoords="offset points")
    for j, label in enumerate(("Actual online total wall seconds", "Actual online total LU", "Actual online total complex RHS")):
        for ax in axs[:, j]:
            ax.set_xlabel(label+" (mean of measured parents)");ax.grid(alpha=.2)
        axs[0, j].set_ylabel("NRMSE / anchor-fill diagnostic; complete plan")
        axs[1, j].set_ylabel("Joint success within total55 / all planned mothers");axs[1, j].set_ylim(-.04, 1.04)
    fig.suptitle("Measured online totals from score ledger; cap values are never substituted. Missing cost coverage is explicit.\nTraining and shared setup are listed separately in one_time_costs.csv; x marks seed17/27 confirmation.")
    fig.tight_layout(rect=(0, 0, 1, .93));fig.savefig(output, dpi=160);plt.close(fig)


def one_time_costs(study):
    rows = []
    for path in sorted((study/"models").glob("*/summary.json")):
        data = json.loads(path.read_text())
        rows.append(dict(scope="one_time_training", run=path.parent.name,
                         measured_seconds=number(data.get("training_seconds")), source=str(path),
                         amortization="NONE; separate from per-case online totals"))
    for path in (study/"models/shared_setup.json", study/"predictions/prediction_runs.json"):
        if path.exists():
            data = json.loads(path.read_text())
            value = data.get("seconds", data.get("global_operator_setup_seconds"))
            rows.append(dict(scope="one_time_shared_setup", run="SHARED", measured_seconds=number(value), source=str(path), amortization="NONE"))
    for path in sorted((study/"predictions").glob("*/inference_setup.json")):
        data = json.loads(path.read_text())
        rows.append(dict(scope="one_time_checkpoint_setup", run=path.parent.name,
                         measured_seconds=number(data.get("checkpoint_and_operator_setup_seconds")), source=str(path), amortization="NONE; operator reuse status retained in source"))
    return rows


def training_curves(study, output):
    paths = sorted((study/"models").glob("*/history.json"))
    if not paths:
        return None
    plt = pyplot();fig, ax = plt.subplots(figsize=(11, 6))
    for i, path in enumerate(paths):
        history = json.loads(path.read_text())
        for split, style in (("train", "--"), ("val", "-")):
            values = [(r["epoch"], number(r.get(split, {}).get("material_l1_normalized")), number(r.get(split, {}).get("edge_l1_normalized"))) for r in history]
            values = [(e, l+.1*d) for e, l, d in values if l is not None and d is not None]
            if values:
                ax.plot(*zip(*values), linestyle=style, color=f"C{i % 10}", label=short(path.parent.name)+" "+split.upper())
    ax.set(xlabel="Epoch", ylabel="Material L1 + .1 edge (public .4 scale)", title="TRAIN/VAL only — validation selects checkpoint; no TEST curve")
    ax.legend(fontsize=6, ncol=2);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(output, dpi=150);plt.close(fig)
    return output.name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=STUDY/"results/TEST_SCORES.json")
    parser.add_argument("--val-scores", type=Path, help="Optional VAL-only network/loss/alpha selection figure")
    parser.add_argument("--output", type=Path, default=STUDY/"figures/generalization")
    args = parser.parse_args();args.output.mkdir(parents=True, exist_ok=True)
    if not args.scores.exists():
        write_json(args.output/"STATUS.json", dict(status="NOT_RUN_NO_FINAL_TEST_SCORES", expected_source=str(args.scores), scientific_figures=[], truth_arrays_read=False))
        print("NOT_RUN_NO_FINAL_TEST_SCORES; no scientific figures generated")
        return
    rows, meta = read_scores(args.scores, "test")
    if not rows:
        raise ValueError("Final TEST score table has no planned rows")
    attach_anchor_diagnostics(rows);parents = parent_rows(rows)
    group_summary = summarize(parents, "group");amplitude_summary = summarize(parents, "amplitude");all_summary = summarize(parents, "all")
    write_csv(args.output/"observation_diagnostics.csv", rows)
    write_csv(args.output/"parent_diagnostics.csv", parents)
    write_csv(args.output/"group_summary.csv", group_summary)
    write_csv(args.output/"amplitude_summary.csv", amplitude_summary)
    write_csv(args.output/"cost_summary.csv", all_summary)
    write_csv(args.output/"coverage.csv", [dict(run=r["run"], seed=r["seed"], group=r["stratum"], planned_parents=r["planned_parents"], planned_observations=r["planned_observations"], selected=r["selected_observations"], not_run=r["not_run_observations"], unresolved_diagnostic=r["unresolved_diagnostic_observations"]) for r in group_summary])
    runs = sorted({r["run"] for r in parents});primary = [r for r in runs if seed_of(r) in (None, 7)]
    confirmation = [r for r in runs if seed_of(r) in (17, 27)]
    files = []
    if primary:
        name="generalization_seed7_primary.png"
        stratum_figure(parents, group_summary, "group", GROUPS, primary, args.output/name, "NEW TEST — seed7 primary and actual unseeded baseline plans")
        files.append(name)
        name="amplitude_seed7_primary.png"
        stratum_figure(parents, amplitude_summary, "amplitude", sorted({r["amplitude"] for r in parents}), primary, args.output/name, "NEW TEST — amplitude strata, seed7 primary")
        files.append(name)
    if confirmation:
        name="generalization_seed17_27_confirmation.png"
        stratum_figure(parents, group_summary, "group", GROUPS, confirmation, args.output/name, "NEW TEST — seed17/27 confirmation, reported separately")
        files.append(name)
    name="seed_matched_mothers.png";paired = seed_figure(parents, args.output/name);files.append(name)
    write_csv(args.output/"seed_confirmation_paired.csv", paired)
    name="actual_online_cost_quality_success.png";cost_figure(all_summary, args.output/name);files.append(name)
    write_csv(args.output/"one_time_costs.csv", one_time_costs(STUDY))
    curve = training_curves(STUDY, args.output/"training_validation_curves.png")
    if curve:files.append(curve)
    if args.val_scores:
        val, _ = read_scores(args.val_scores, "val");attach_anchor_diagnostics(val);vp = parent_rows(val)
        vs = summarize(vp, "all");vr = sorted({r["run"] for r in vp})
        if vr:
            name="VAL_ONLY_network_loss_alpha_selection.png"
            stratum_figure(vp, vs, "all", ["ALL"], vr, args.output/name, "VAL MODEL SELECTION ONLY — network / loss / alpha; not generalization evidence")
            files.append(name);write_csv(args.output/"VAL_ONLY_selection_summary.csv", vs)
    write_json(args.output/"STATUS.json", dict(status="GENERATED_FROM_FINAL_SCALAR_SCORES", source=str(args.scores), seal_path=meta.get("seal_path"), scientific_figures=files,
        planned_observations=len(rows), actual_runs=runs, bootstrap_seed=BOOTSTRAP_SEED, bootstrap_repetitions=BOOTSTRAP_REPETITIONS,
        statistical_unit="mother; clean/noisy collapsed within mother; seed runs never pooled as independent materials",
        diagnostic="No top1 => already-scored common anchor fill only; never success. Missing anchor => missing overall diagnostic mean, never selected-subset mean.",
        coverage="Only rows in actual score plan. No unplanned baseline coverage invented. NOT_RUN and FAIL stay in success denominator.",
        cost_scope="actual scored online total wall/LU/RHS; separate unamortized training/setup receipts; missing costs never zero-imputed", truth_arrays_read=False))
    print(json.dumps(dict(status="GENERATED_FROM_FINAL_SCALAR_SCORES", output=str(args.output), figures=len(files))))


if __name__ == "__main__":
    main()
