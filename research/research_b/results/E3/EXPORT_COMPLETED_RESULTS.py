"""Export completed E3 records only. No solver, optimizer or inference imports."""
import csv
import json
import shutil
import sys
from pathlib import Path

BROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BROOT / "code"))
from cost import Cost

E3 = Path(__file__).resolve().parent
METHODS = (
    "A_robust_action_then_profile",
    "A_robust_then_protected_eps_corrector",
    "A_robust_then_independent_prior_multistart_corrector",
)


def save_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+"\n")


def csv_value(value):
    if value is None or isinstance(value, (list, dict, bool)):
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    return value


def write_csv(path, rows, first):
    columns = list(first)
    for row in rows:
        columns.extend(k for k in row if k not in columns)
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: csv_value(v) for k, v in row.items()})


def export():
    source = json.loads((E3 / "RESULTS.json").read_text())
    records = source["records"]
    parents = [p["parent"] for p in source["parent_summaries"]]
    summaries = source["summary"]
    if len(records) != 36 or len(summaries) != 9 or len(parents) != 4:
        raise ValueError("Unexpected completed E3 counts; preserve every method/parent")
    outer = BROOT / "AcquisitionRecommendation/E3"
    if outer.exists():
        raise FileExistsError("Do not overwrite an existing completed E3 handoff export")
    outer.mkdir(parents=True)
    entries = []
    new_counts = {"protected": 0, "prior_multistart": 0}
    total_files, gross_source_candidates = 0, 0
    for parent in parents:
        parent_result = E3 / parent
        for method in METHODS:
            original = parent_result / "AcquisitionRecommendation" / method
            destination = outer / parent / method
            shutil.copytree(original, destination)
            package = json.loads((destination / "MaterialHypothesisPackage.json").read_text())
            history = json.loads((destination / "updated_data_history.json").read_text())
            record = next(r for r in records if r["parent"] == parent and r["method"] == method)
            if package["weights"]["kind"] != "none":
                raise ValueError("E3 handoff cannot be exported as posterior-weighted A candidates")
            if len(history["observations"]) != 64:
                raise ValueError("Completed E3 updated history must retain old32 and new32 once")
            raw = destination / "raw_e3"
            raw.mkdir()
            raw_names = ("A_FULL_ACTION_FIELDS.npz", "SIMULATED_NEW_OBSERVATIONS.npz",
                         "ACTION_RANKINGS.json", "ENDPOINT_TRIALS.jsonl", "PRIOR_MULTISTART_INITIALS.json")
            for filename in raw_names:
                shutil.copy2(parent_result / filename, raw / filename)
            corrector = "corrector" in method
            status = None
            if corrector:
                src_catalogue = parent_result / method
                shutil.copytree(src_catalogue, destination / "completed_corrector")
                status = json.loads((src_catalogue / "CATALOGUE_STATUS.json").read_text())
                kind = "protected" if "protected" in method else "prior_multistart"
                new_counts[kind] += status["installed_endpoints"]
            followups = []
            for completed_method in METHODS[1:]:
                completed_status = json.loads((parent_result / completed_method / "CATALOGUE_STATUS.json").read_text())
                followups.append({
                    "method": completed_method,
                    "relative_handoff": "../"+completed_method,
                    "status": "COMPLETED_FULL_VALIDATED_CORRECTOR",
                    "K_original": completed_status["K_original"],
                    "K_full_validated": completed_status["K_full_validated"],
                    "installed_endpoints": completed_status["installed_endpoints"],
                })
            receipt = {
                "status": "COMPLETED_FULL_VALIDATED_CORRECTOR_HANDOFF" if corrector else "ROBUST_MEASUREMENT_COMPLETED_WITH_TWO_COMPLETED_CORRECTOR_FOLLOWUPS",
                "source_completed_experiment": str(original),
                "source_updated_history_status_preserved": history["candidate_refresh_status"],
                "completed_corrector_followups": followups,
                "completion_explanation": "A source assets and original history fields are preserved verbatim. The generic raw-measurement handoff may say PENDING; the two independently completed B correctors are indexed here and are no longer pending.",
                "parent": parent, "method": method,
                "raw_material_package_K": len(package["candidates"]),
                "candidate_weights_kind": "none", "posterior_assigned": False,
                "measurement_action_id": history["action_id"],
                "old_plus_new_history_real_channels": len(history["observations"]),
                "shared_gain_rule": "Each material world has one bounded gain across complete history; different worlds may use different values",
                "raw_prediction_semantics": "prediction_i.npy is source-gain-applied g_A*F at1.6; unscaled_old_fields.npy is F. raw_e3/A_FULL_ACTION_FIELDS.npz and completed catalogue joint_unscaled_fields contain full unscaled fields.",
                "completed_catalogue_uri": "completed_corrector/FULL_VALIDATED_CATALOGUE.npz" if corrector else None,
                "raw_endpoint_log_uri": "raw_e3/ENDPOINT_TRIALS.jsonl",
                "metrics_from_completed_results": record,
                "scope": "One of four actual A parents, descriptive; no full-space coverage or posterior claim",
                "new_Maxwell_attempts_this_export": 0,
            }
            if status is not None:
                receipt["completed_corrector_counts"] = status
            save_json(destination / "HANDOFF_COMPLETION_RECEIPT.json", receipt)
            count = sum(p.is_file() for p in destination.rglob("*"))
            total_files += count
            gross_source_candidates += len(package["candidates"])
            entries.append(dict(parent=parent, method=method,
                                directory=str(destination.relative_to(BROOT)),
                                selected_MAE=record["selected_MAE"], best_of_K_MAE=record["best_of_K_MAE"],
                                K=record["K"], candidate_weights_kind="none", files=count,
                                completion_receipt="HANDOFF_COMPLETION_RECEIPT.json"))
    parent_rows = []
    for source_record in records:
        row = dict(source_record)
        row["descriptive_parent_count"] = 4
        row["reference_status"] = "UNRESOLVED_AFTER_NEW_HISTORY" if "reference" in row["method"] else "NOT_REQUIRED"
        entry = next((e for e in entries if e["parent"] == row["parent"] and e["method"] == row["method"]), None)
        row["exported_handoff_directory"] = entry["directory"] if entry else None
        parent_rows.append(row)
    write_csv(E3 / "PARENT_METRICS.csv", parent_rows, ("parent", "method"))
    summary_rows = []
    for summary in summaries:
        row = dict(summary)
        group = [r for r in records if r["method"] == row["method"]]
        row["reference_status"] = "UNRESOLVED_AFTER_NEW_HISTORY" if "reference" in row["method"] else "NOT_REQUIRED"
        row["A_posterior_assigned"] = False
        row["total_installed_endpoint_records"] = sum(r.get("new_full_validated_candidates", 0) for r in group)
        row["installed_count_scope"] = "within this method; base-control endpoints reuse protected work and must not be summed as new unique physical work"
        summary_rows.append(row)
    write_csv(E3 / "SUMMARY_METRICS.csv", summary_rows, ("method", "parents", "selected_MAE", "best_of_K_MAE"))
    save_json(outer / "INDEX.json", {
        "status": "REAL_FOUR_A_PARENT_E3_HANDOFFS_EXPORTED", "parents": parents,
        "handoffs": entries, "candidate_weights_kind": "none", "A_posterior_assigned": False,
        "completed_installed_endpoint_records": new_counts,
        "new_Maxwell_attempts_this_export": 0,
    })
    (outer / "README.md").write_text(
        "# E3：四份真实 Research A 材料候选交接\n\n"
        "这里只导出已完成实验，没有重跑 Maxwell、优化或更新。每个 parent 提供 robust 新测量、protected epsilon corrector、三起点连续材料 corrector 三份独立目录。\n\n"
        "每份 `HANDOFF_COMPLETION_RECEIPT.json` 索引两个已经完成的 corrector。原 robust 测量 history 中的通用 `PENDING_INDEPENDENT_RECONSTRUCTION` 字段按原样保留；它不表示后续 corrector 尚未执行。实际完成状态、完整候选目录与负端点日志见外层 receipt 和 `completed_corrector/`。\n\n"
        "A 原材料与 gain-applied 预测、去 gain 的原频响应、旧观测和完整 old+new history 都保留。候选没有 posterior 权重；residual 只作共享有界 gain 的 profile 选择。四个 parent 仅支持描述性结论，不能与 B 的32个 parent 混作一个测试集。\n\n"
        "`raw_e3/` 保留已支付 full fields、真实新增模拟观测和全部端点记录。复制与导出新增物理调用为0。\n",
        encoding="utf-8")
    # The figure is a static research artifact, not a new evaluation.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    labels = {
        "A_base_only_no_measurement": "A base / no new data",
        "A_fixed_candidate_profile_no_measurement": "A candidates / no new data",
        "A_robust_action_then_profile": "Robust action / fixed A catalogue",
        "A_fixed3.2_then_profile": "Fixed kR=3.2 / fixed A catalogue",
        "A_reference_action_then_profile": "Reference action / A set [UNRESOLVED]",
        "A_robust_then_protected_eps_corrector": "Robust / protected epsilon corrector",
        "A_robust_then_independent_prior_multistart_corrector": "Robust / three-start material corrector",
        "declared_continuous_P2_reference_comparator": "Continuous reference [UNRESOLVED]",
        "A_base_only_same_action_eps_corrector": "Base epsilon corrector / same robust action",
    }
    ordered = [s["method"] for s in summaries]
    yy = np.arange(len(ordered))
    colors = ["#c1c9d2" if "no_measurement" in m else "#d9bd8f" if "reference" in m
              else "#67a88d" if "corrector" in m and "base_only" not in m
              else "#a9a1bb" if "base_only_same" in m else "#7e9fbf" for m in ordered]
    fig, axes = plt.subplots(1, 2, figsize=(16.5, 8.8), sharey=True)
    for ax, metric, title in zip(axes, ("selected_MAE", "best_of_K_MAE"),
                                 ("Selected material MAE", "Best-of-K material MAE")):
        means = [next(s[metric] for s in summaries if s["method"] == m) for m in ordered]
        bars = ax.barh(yy, means, color=colors, height=.65, zorder=2)
        for bar, method in zip(bars, ordered):
            if "reference" in method:
                bar.set_hatch("///")
                bar.set_edgecolor("#8f734f")
        maximum = max(r[metric] for r in records)
        ax.set_xlim(0, maximum*1.18)
        for i, method in enumerate(ordered):
            values = [r[metric] for r in records if r["method"] == method]
            ax.scatter(values, i+np.linspace(-.12, .12, len(values)), s=15,
                       color="#293340", alpha=.8, zorder=3)
            ax.text(means[i]+maximum*.018, i+.28, f"{means[i]:.3f}", fontsize=10, color="#293340")
        ax.set_title(title, fontsize=15, loc="left", pad=12)
        ax.set_xlabel("Volume-weighted permittivity error (lower is better)", fontsize=11)
        ax.grid(axis="x", color="#e4e7eb", zorder=0)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
    axes[0].set_yticks(yy, [labels[m] for m in ordered], fontsize=11)
    axes[0].invert_yaxis()
    fig.suptitle("Real A handoff: four parents, descriptive comparison", x=.32,
                 ha="left", fontsize=19, fontweight="bold")
    fig.text(.32, .92, "Bars: four-parent means. Dots: individual parents. No significance claim.", fontsize=11, color="#505967")
    fig.text(.04, .085,
             "Reference-driven rows: new-history audit UNRESOLVED; these are not reliable upper bounds.\n"
             "Fixed-catalogue Best-of-K is unchanged by measurement selection. Only full-validated refresh adds candidates.\n"
             "The base epsilon control uses the same action selected from the full A set; it is not pure-base sensing.",
             fontsize=10.5, linespacing=1.6, color="#444d59")
    fig.subplots_adjust(left=.32, right=.98, bottom=.20, top=.87, wspace=.12)
    fig.savefig(E3 / "E3_MATERIAL_COMPARISON.png", dpi=180, facecolor="white")
    fig.savefig(E3 / "E3_MATERIAL_COMPARISON.pdf", facecolor="white")
    plt.close(fig)
    receipt = {
        "status": "COMPLETED_RESULTS_EXPORTED_NO_RERUN",
        "source_results": str(E3 / "RESULTS.json"),
        "source_parents": parents, "independent_parent_count": 4,
        "method_count": len(summaries), "summary_csv_rows": len(summary_rows),
        "parent_csv_rows": len(parent_rows), "expected_parent_method_rows": 4*9,
        "all_source_parent_method_rows_preserved": len({(r["parent"], r["method"]) for r in parent_rows}) == 36,
        "source_distinct_original_A_candidates": sum(r["K"] for r in records if r["method"] == "A_fixed_candidate_profile_no_measurement"),
        "handoff_bundle_count": len(entries), "handoff_methods_per_parent": list(METHODS),
        "handoff_material_arrays_gross_copies": gross_source_candidates,
        "handoff_files_excluding_outer_index_and_readme": total_files,
        "completed_installed_endpoint_records": new_counts,
        "base_capacity_control_reused_protected_endpoints": sum(r.get("new_full_validated_candidates", 0) for r in records if r["method"] == "A_base_only_same_action_eps_corrector"),
        "new_Maxwell_attempts_this_export": 0, "new_optimizations_this_export": 0,
        "source_E3_callback_reservations": source["local_cost"]["callback_reservations"],
        "A_posterior_assigned": False, "reference_status": "UNRESOLVED_AFTER_NEW_HISTORY",
        "paths": {
            "summary_csv": str(E3 / "SUMMARY_METRICS.csv"),
            "parent_csv": str(E3 / "PARENT_METRICS.csv"),
            "handoff_index": str(outer / "INDEX.json"),
            "figure_png": str(E3 / "E3_MATERIAL_COMPARISON.png"),
            "figure_pdf": str(E3 / "E3_MATERIAL_COMPARISON.pdf"),
        },
        "visual_review": "PENDING_LOCAL_IMAGE_INSPECTION",
        "scope": "Four actual A parents only; no numerical experiment or scientific ranking was rerun",
    }
    save_json(E3 / "SUMMARY_EXPORT.json", receipt)
    return receipt


if __name__ == "__main__":
    cost = Cost()
    with cost.job("e3-results-export-v1", "ANALYSIS"):
        receipt = export()
    job = next(r for r in cost.rows(cost.jobs) if r["id"] == "e3-results-export-v1")
    receipt["analysis_cost_receipt"] = job
    save_json(E3 / "SUMMARY_EXPORT.json", receipt)
    print(json.dumps({k: receipt[k] for k in ("status", "summary_csv_rows", "parent_csv_rows",
                                             "handoff_bundle_count", "new_Maxwell_attempts_this_export")}))
