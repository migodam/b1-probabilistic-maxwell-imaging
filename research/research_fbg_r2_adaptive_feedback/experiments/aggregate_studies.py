"""Offline DEVELOPMENT aggregation for Stage 2 and Stage 3.

Reads only their frozen manifests/receipts and already disclosed R4 truth.
No solver, physical forward, Stage 4 public data, or fresh truth import.
Run only after the corresponding remote receipts have been fetched.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import importlib.util
import itertools
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parents[1]
R4 = HERE.parent / "research_z_alg_r4_diversity"
METHODS2 = ("native_trf", "reciprocal_trf", "data_lbfgs", "gd_all",
            "feedback_beta:0.2", "feedback_qbeta")
METHODS3 = ("native_trf", "gd_all", "feedback_beta:0.2", "feedback_qbeta")
COLORS = {"native_trf": "#4c78a8", "reciprocal_trf": "#72b7b2",
          "data_lbfgs": "#9c755f", "gd_all": "#777777",
          "feedback_beta:0.2": "#f28e2b", "feedback_qbeta": "#b279a2"}


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def write(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def r4_metrics():
    spec = importlib.util.spec_from_file_location("r4_offline_scorer", R4 / "score_results.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.metrics


def old_truth() -> dict:
    cases = read(HERE / "data_contract/stage1_cases.json")["cases"]
    out = {}
    for c in cases:
        ident = Path(c["public_path"]).stem
        with np.load(R4 / "truth_after_freeze" / f"{ident}.npz", allow_pickle=False) as f:
            out[c["case_id"]] = {
                "m": np.asarray(f["m"], float).ravel(),
                "mask": np.asarray(f["mask"], bool).ravel(),
                "amplitude": float(c["amplitude"]),
                "background": float(np.min(f["m"])),
            }
    return out


def cost(receipt: dict | None) -> dict:
    if receipt is None:
        return {"wall_seconds": None, "LU": None, "actual_complex_RHS": None,
                "source_RHS": None, "adjoint_RHS": None, "receiver_RHS": None,
                "tangent_RHS": None}
    result = receipt.get("result") or {}
    counts = result.get("counts", result.get("costcounts")) or {}
    rhs = counts.get("actual_complex_RHS")
    if rhs is None and "source_RHS" in counts:
        rhs = sum(int(counts.get(k, 0)) for k in
                  ("source_RHS", "adjoint_RHS", "receiver_RHS", "tangent_RHS"))
    return {"wall_seconds": receipt.get("wall_s", result.get("task_wall_s")),
            "LU": counts.get("LU"), "actual_complex_RHS": rhs,
            "source_RHS": counts.get("source_RHS"),
            "adjoint_RHS": counts.get("adjoint_RHS"),
            "receiver_RHS": counts.get("receiver_RHS"),
            "tangent_RHS": counts.get("tangent_RHS")}


def cost_sums(rows: list[dict]) -> dict:
    out = {}
    for key in ("wall_seconds", "LU", "actual_complex_RHS", "source_RHS",
                "adjoint_RHS", "receiver_RHS", "tangent_RHS"):
        values = [r["cost"][key] for r in rows if r["cost"].get(key) is not None]
        out[key + "_sum_known"] = float(sum(values)) if values else None
        out[key + "_known_starts"] = len(values)
    out["scope"] = "sum of charged independent task work, not isolated GPU latency or parallel batch wall"
    return out


def scored_task(task: dict, results_dir: Path, truth: dict,
                metric_fn) -> dict:
    path = results_dir / (task["task_id"] + ".json")
    receipt = read(path) if path.exists() else None
    result = (receipt or {}).get("result") or {}
    status = "NOT_RUN" if receipt is None else receipt.get("status", "UNKNOWN")
    material = result.get("m")
    finite = False
    if material is not None:
        try:
            m = np.asarray(material, float).ravel()
            finite = m.size == 576 and bool(np.isfinite(m).all())
        except (TypeError, ValueError):
            finite = False
    row = {k: task.get(k) for k in
           ("task_id", "case_id", "method", "regime", "K_run", "start_index",
            "start_name", "p", "layout", "energy_protocol")}
    row.update({"transport_status": status, "solver_status": result.get("wrapper_status", result.get("status")),
                "stop_reason": result.get("wrapper_stop_reason", result.get("stop_reason")),
                "finite_endpoint": finite, "cost": cost(receipt),
                "ranking_objective": result.get("ranking_objective",
                                    (result.get("top1") or {}).get("ranking_score")),
                "relative_residual": result.get("relative_residual",
                                     (result.get("top1") or {}).get("relative_residual")),
                "all_unknowns_released": result.get("stage") == 24 or
                                         result.get("all_unknowns_released"),
                "feedback_records": [{"stage": f.get("stage"), "kind": f.get("kind"),
                                      "safety": f.get("safety"), "whole_box": f.get("whole_box"),
                                      "tangent_energy_ratio": f.get("tangent_energy_ratio")}
                                     for f in result.get("feedback_records", [])]})
    if finite:
        t = truth[task["case_id"]]
        row.update(metric_fn(m, t["m"], t["mask"], t["amplitude"], t["background"]))
        row["material"] = m  # retained only in memory for pair-RMS calculations
    else:
        row.update({"RMSE": None, "NRMSE": None, "MAE": None, "IoU_Otsu": None,
                    "IoU_known_midpoint_OFFLINE": None})
    limit = .001  # Both development sets are clean only.
    data_ok = row["relative_residual"] is not None and np.isfinite(row["relative_residual"]) and row["relative_residual"] <= limit
    material_ok = finite and row["NRMSE"] <= .15 and row["IoU_Otsu"] >= .7
    row["data_gate"] = "NOT_RUN" if status == "NOT_RUN" else "PASS" if data_ok else "FAIL"
    row["material_gate"] = "NOT_RUN" if status == "NOT_RUN" else "PASS" if material_ok else "FAIL"
    row["joint_gate"] = "NOT_RUN" if status == "NOT_RUN" else "PASS" if data_ok and material_ok else "FAIL"
    return row


def public_row(row: dict) -> dict:
    return {k: v for k, v in row.items() if k != "material"}


def setup_receipt(results_dir: Path) -> dict:
    path = results_dir / "SUMMARY.json"
    if not path.exists():
        return {"status": "NOT_RETURNED", "batch_wall_seconds": None,
                "shared_setup_wall_seconds": None, "shared_setup_counts": None}
    summary = read(path)
    return {"status": "RECORDED", "batch_wall_seconds": summary.get("wall_s"),
            "shared_setup_wall_seconds": summary.get("setup_wall_s"),
            "shared_setup_counts": summary.get("setup_counts"),
            "summary_path": str(path)}


def stage2(manifest_path: Path, results_dir: Path, out_dir: Path, truth: dict,
           metric_fn) -> dict:
    manifest = read(manifest_path)
    tasks = manifest["tasks"]
    if len(tasks) != 360:
        raise ValueError("Actual Stage 2 frozen manifest must have 360 jobs")
    rows = [scored_task(t, results_dir, truth, metric_fn) for t in tasks]
    lookup = {(r["case_id"], r["method"], r["regime"], r["K_run"], r["start_index"]): r
              for r in rows}
    cases = sorted({r["case_id"] for r in rows})
    groups = []
    for case, method in itertools.product(cases, METHODS2):
        for regime, k_run, k in (("per_start", 8, 1), ("per_start", 8, 4),
                                 ("per_start", 8, 8), ("fixed_total", 4, 4),
                                 ("fixed_total", 8, 8)):
            selected = [lookup[(case, method, regime, k_run, i)] for i in range(k)]
            finite = [r for r in selected if r["finite_endpoint"]]
            ranked = [r for r in finite if r["ranking_objective"] is not None and
                      np.isfinite(r["ranking_objective"])]
            online = min(ranked, key=lambda r: (r["ranking_objective"], r["start_index"])) if ranked else None
            oracle = min(finite, key=lambda r: (r["RMSE"], r["start_index"])) if finite else None
            fits = [r for r in selected if r["data_gate"] == "PASS"]
            pair_rms = [float(np.sqrt(np.mean((a["material"]-b["material"])**2)))
                        for a, b in itertools.combinations(fits, 2)]
            groups.append({"case_id": case, "method": method, "regime": regime,
                           "K": k, "source_run_K": k_run,
                           "nested_prefix_of_same_K8_run": regime == "per_start",
                           "planned_starts": k, "finite_starts": len(finite),
                           "failed_starts": sum(r["joint_gate"] == "FAIL" for r in selected),
                           "NOT_RUN_starts": sum(r["joint_gate"] == "NOT_RUN" for r in selected),
                           "start_task_ids": [r["task_id"] for r in selected],
                           "online_top1_start_index": online["start_index"] if online else None,
                           "online_top1_NRMSE": online["NRMSE"] if online else None,
                           "online_top1_IoU_Otsu": online["IoU_Otsu"] if online else None,
                           "online_top1_relative_residual": online["relative_residual"] if online else None,
                           "online_top1_joint_gate": online["joint_gate"] if online else
                                                      "NOT_RUN" if not selected or all(r["joint_gate"] == "NOT_RUN" for r in selected) else "FAIL",
                           "any_joint_at_K_OFFLINE": any(r["joint_gate"] == "PASS" for r in selected),
                           "oracle_best_at_K_OFFLINE": None if oracle is None else
                               {"start_index": oracle["start_index"], "NRMSE": oracle["NRMSE"],
                                "RMSE": oracle["RMSE"], "IoU_Otsu": oracle["IoU_Otsu"],
                                "joint_gate": oracle["joint_gate"]},
                           "data_fit_pair_RMS": {"pair_count": len(pair_rms),
                                "min": min(pair_rms) if pair_rms else None,
                                "max": max(pair_rms) if pair_rms else None,
                                "values": pair_rms,
                                "scope": "only pairs whose separate full nonlinear data gates PASS; distinct images, not proven roots"},
                           "measured_cost": cost_sums(selected)})
    # Fixed-total K=1 uses the already-paid per-start s0, no new task or cost.
    for g in [x for x in groups if x["regime"] == "per_start" and x["K"] == 1]:
        alias = {**g, "regime": "fixed_total", "source_run_K": 8,
                 "K1_alias_of_per_start_s0": True,
                 "nested_prefix_of_same_K8_run": True}
        groups.append(alias)
    conditional = []
    for regime, k_run in (("per_start", 8), ("fixed_total", 4), ("fixed_total", 8)):
        for method in METHODS2:
            for i in range(k_run):
                subset = [lookup[(case, method, regime, k_run, i)] for case in cases]
                conditional.append({"regime": regime, "K_run": k_run,
                                    "method": method, "start_index": i,
                                    "start_name": subset[0]["start_name"],
                                    "planned_cases": len(cases),
                                    "finite_cases": sum(r["finite_endpoint"] for r in subset),
                                    "joint_hit_cases": sum(r["joint_gate"] == "PASS" for r in subset),
                                    "failed_cases": sum(r["joint_gate"] == "FAIL" for r in subset),
                                    "NOT_RUN_cases": sum(r["joint_gate"] == "NOT_RUN" for r in subset)})
    result = {"schema": "fbg-r2-stage2-aggregate-v1", "scope": "DEVELOPMENT_DISCLOSED_R4",
              "manifest": str(manifest_path), "results_dir": str(results_dir),
              "truth_source": "R4/truth_after_freeze; already public DEVELOPMENT truth",
              "planned_jobs": 360, "returned_receipts": sum(r["transport_status"] != "NOT_RUN" for r in rows),
              "setup_and_batch": setup_receipt(results_dir),
              "cost_rule": "per-start K prefixes reuse K8 receipts; do not sum nested K groups as independent jobs. Fixed-total K4/8 are independent runs. Shared setup is separate; sum worker wall is charged work, not GPU batch elapsed.",
              "selection_rule": "online minimum finite original full physical ranking_objective; offline best@K uses disclosed GT only",
              "not_planned_extensions": [
                  {"family": "raw_uniform", "status": "NOT_RUN", "planned_denominator": 0},
                  {"family": "wavelet", "status": "NOT_RUN", "planned_denominator": 0},
                  {"K": 16, "status": "NOT_RUN", "planned_denominator": 0}],
              "tasks": [public_row(r) for r in rows], "groups": groups,
              "conditional_hits_by_initial": conditional}
    write(out_dir / "STAGE2_AGGREGATE.json", result)
    plot_stage2(result, out_dir.parent / "figures")
    return result


def gram_diag(diag_path: Path | None) -> tuple[dict, dict]:
    if diag_path is None or not diag_path.exists():
        return {}, {"status": "NOT_RETURNED", "path": None}
    data = read(diag_path)
    source_rows = data.get("rows", data.get("tasks", [])) if isinstance(data, dict) else data
    by_acquisition = {}
    for row in source_rows if isinstance(source_rows, list) else []:
        if not isinstance(row, dict):
            continue
        selected = {k: v for k, v in row.items()
                    if "gram" in k.lower() or "condition" in k.lower() or
                       k in ("acquisition", "source_row_real_compatibility_residual",
                             "source_row_complex_compatibility_residual", "apparent_real_pass",
                             "interpretation_scope")}
        if row.get("acquisition"):
            by_acquisition[row["acquisition"]] = selected
    return by_acquisition, {"status": "RECORDED", "path": str(diag_path),
                     "rows": len(source_rows) if isinstance(source_rows, list) else None}


def stage3(manifest_path: Path, results_dir: Path, out_dir: Path, truth: dict,
           metric_fn, diag_path: Path | None) -> dict:
    manifest = read(manifest_path)
    tasks = manifest["tasks"]
    if len(tasks) != 64:
        raise ValueError("Actual Stage 3 frozen manifest must have 64 jobs")
    rows = [scored_task(t, results_dir, truth, metric_fn) for t in tasks]
    diag, diag_meta = gram_diag(diag_path)
    for row in rows:
        row["source_Gram_diagnostic"] = diag.get(f"{row['layout']}_p{row['p']}")
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["case_id"], row["layout"], row["p"])].append(row)
    comparison = []
    for (case, layout, p), group in sorted(grouped.items()):
        comparison.append({"case_id": case, "layout": layout, "p": p,
                           "planned_methods": len(group),
                           "method_rows": [public_row(r) for r in group],
                           "joint_pass_methods": [r["method"] for r in group if r["joint_gate"] == "PASS"]})
    result = {"schema": "fbg-r2-stage3-aggregate-v1", "scope": "DEVELOPMENT_DISCLOSED_R4",
              "manifest": str(manifest_path), "results_dir": str(results_dir),
              "truth_source": "R4/truth_after_freeze; already public DEVELOPMENT truth",
              "planned_jobs": 64, "returned_receipts": sum(r["transport_status"] != "NOT_RUN" for r in rows),
              "setup_and_batch": setup_receipt(results_dir),
              "source_Gram_diagnostics": diag_meta,
              "acquisition_diagnostics_public_anchor": diag,
              "acquisition_scope": "two disclosed objects, 4 p values, uniform/random layouts, 4 methods. Changing p changes physical source positions/RHS; layout and p effects remain finite and object-specific.",
              "energy_control": "fixed-total clean is analytic global amplitude scaling; no duplicate inverse task. Fixed-absolute-noise extension NOT_RUN.",
              "p24_failure_policy": "retain all p24 failures and missing receipts, without dropping them from denominators",
              "tasks": [public_row(r) for r in rows], "paired_by_acquisition": comparison}
    write(out_dir / "STAGE3_AGGREGATE.json", result)
    plot_stage3(result, out_dir.parent / "figures")
    return result


def plot_stage2(data: dict, figures: Path) -> None:
    figures.mkdir(parents=True, exist_ok=True)
    # These are the eight actually paid per-start K8 receipts, each scored on
    # three disclosed objects. No invented raw-uniform/wavelet/K16 points.
    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
    for method in METHODS2:
        rows = [r for r in data["conditional_hits_by_initial"]
                if r["regime"] == "per_start" and r["method"] == method]
        rows.sort(key=lambda r: r["start_index"])
        ax.plot([r["start_index"] for r in rows],
                [r["joint_hit_cases"] for r in rows], marker="o", lw=2,
                color=COLORS[method], label=method)
    names = [r["start_name"] for r in data["conditional_hits_by_initial"]
             if r["regime"] == "per_start" and r["method"] == METHODS2[0]]
    ax.set_xticks(range(8), names, rotation=30, ha="right")
    ax.set_yticks(range(4)); ax.set_ylim(-.1, 3.2)
    ax.set_ylabel("Joint hits / 3 planned disclosed cases")
    ax.set_title("Stage 2 DEVELOPMENT · conditional success by public initial state")
    ax.grid(alpha=.2); ax.legend(ncol=3, fontsize=8)
    fig.savefig(figures / "success_by_initials.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True, sharey=True)
    for ax, regime in zip(axes, ("per_start", "fixed_total")):
        for method in METHODS2:
            xs, ys, ks = [], [], []
            for k in (1, 4, 8):
                group = [g for g in data["groups"] if g["regime"] == regime and
                         g["method"] == method and g["K"] == k]
                wall = [g["measured_cost"]["wall_seconds_sum_known"] for g in group]
                xs.append(sum(x for x in wall if x is not None) if any(x is not None for x in wall) else np.nan)
                ys.append(sum(g["online_top1_joint_gate"] == "PASS" for g in group))
                ks.append(k)
            ax.plot(xs, ys, marker="o", lw=2, color=COLORS[method], label=method)
            for x, y, k in zip(xs, ys, ks):
                if np.isfinite(x):
                    ax.annotate(f"K{k}", (x, y), xytext=(3, 2), textcoords="offset points", fontsize=6)
        ax.set_title(regime.replace("_", " "))
        ax.set_xlabel("Cumulative charged solver work: sum worker wall (s) across 3 cases")
        ax.grid(alpha=.2)
    axes[0].set_ylabel("Online top-1 joint successes / 3 cases")
    axes[0].set_yticks(range(4))
    axes[1].legend(ncol=2, fontsize=7)
    fig.suptitle("Stage 2 DEVELOPMENT · success versus paid work (not isolated GPU latency)")
    fig.savefig(figures / "success_vs_budget.png", dpi=150)
    plt.close(fig)


def plot_stage3(data: dict, figures: Path) -> None:
    figures.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True, sharex=True)
    cases = ("dev_high_digit8", "dev_high_ring")
    layouts = ("uniform", "random")
    for i, case in enumerate(cases):
        for j, layout in enumerate(layouts):
            ax = axes[i, j]
            for method in METHODS3:
                group = [r for r in data["tasks"] if r["case_id"] == case and
                         r["layout"] == layout and r["method"] == method]
                group.sort(key=lambda r: r["p"])
                ax.plot([r["p"] for r in group],
                        [np.nan if r["NRMSE"] is None else r["NRMSE"] for r in group],
                        marker="o", lw=1.8, label=method, color=COLORS[method])
            ax.axhline(.15, color="#368a5d", ls="--", lw=1)
            ax.set_title(f"{case} · {layout}")
            ax.set_xticks((3, 6, 12, 24)); ax.set_xlabel("Physical sources p")
            ax.set_ylabel("Material NRMSE (gate also requires IoU and data fit)")
            ax.grid(alpha=.2)
    axes[0, 1].legend(fontsize=7)
    fig.suptitle("Stage 3 DEVELOPMENT · finite physical acquisition outcomes; p24 retained")
    fig.savefig(figures / "stage3_acquisition.png", dpi=150)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=("stage2", "stage3"))
    ap.add_argument("--results-dir", type=Path)
    ap.add_argument("--manifest", type=Path)
    ap.add_argument("--output-dir", type=Path, default=HERE / "results")
    ap.add_argument("--stage3-diag", type=Path,
                    default=HERE / "results/illumination_current_compatibility.json")
    a = ap.parse_args()
    manifest = a.manifest or HERE / "data_contract" / f"{a.stage}_remote.json"
    results_dir = a.results_dir or HERE / "results" / a.stage
    if not results_dir.exists():
        ap.error(f"Receipts not fetched yet: {results_dir}; no partial aggregate written")
    truth, metrics = old_truth(), r4_metrics()
    if a.stage == "stage2":
        result = stage2(manifest, results_dir, a.output_dir, truth, metrics)
    else:
        result = stage3(manifest, results_dir, a.output_dir, truth, metrics, a.stage3_diag)
    print(json.dumps({"stage": a.stage, "planned": result["planned_jobs"],
                      "returned": result["returned_receipts"], "output_dir": str(a.output_dir)}))


if __name__ == "__main__":
    main()
