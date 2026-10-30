"""Merge frozen study receipts into a cost ledger and delivery indexes.

Offline only. This module never opens private Stage 4 truth, imports a solver,
or reruns physics. Run after Stage 4 predictions have been sealed and scored.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
RESULTS = HERE / "results"
AUDIT = HERE / "audit"
STAGES = ("stage1", "stage2", "stage3", "stage4")


def read(path: Path):
    return json.loads(path.read_text())


def write(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def csv_write(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k) for k in columns})


def numeric(value):
    if isinstance(value, bool) or value is None:
        return None
    try:
        x = float(value)
        return x if x == x and abs(x) != float("inf") else None
    except (TypeError, ValueError):
        return None


def canonical_counts(counts: dict | None) -> dict:
    """Count full Maxwell RHS once; expose auxiliary work separately."""
    if not isinstance(counts, dict):
        return {k: None for k in ("physical_LU", "physical_source_RHS", "adjoint_RHS",
                                  "receiver_RHS", "tangent_RHS", "full_complex_RHS",
                                  "joint_source_RHS_auxiliary", "calibration_adjoint_RHS_subset",
                                  "joint_cholesky_auxiliary", "source_Gram_Cholesky_auxiliary",
                                  "Schur_RHS_auxiliary", "Gram_factorizations_auxiliary")}
    source = numeric(counts.get("source_RHS", counts.get("initial_source_RHS")))
    joint = numeric(counts.get("joint_source_RHS"))
    explicit_physical = numeric(counts.get("physical_source_RHS"))
    physical_source = (source - (joint or 0) if source is not None else explicit_physical)
    if (physical_source is None and counts.get("actual_complex_RHS") is not None and
            all(counts.get(k) is None for k in ("adjoint_RHS", "receiver_RHS", "tangent_RHS"))):
        physical_source = numeric(counts.get("actual_complex_RHS"))
    if physical_source is not None and physical_source < 0:
        raise ValueError("joint_source_RHS exceeds source_RHS")
    adjoint = numeric(counts.get("adjoint_RHS", counts.get("data_adjoint_RHS")))
    receiver = numeric(counts.get("receiver_RHS", counts.get("receiver_transpose_RHS")))
    tangent = numeric(counts.get("tangent_RHS"))
    if tangent is None:
        parts = [numeric(counts.get(k)) for k in
                 ("exact_current_tangent_RHS", "directional_tangent_RHS", "tangent_source_RHS")]
        tangent = sum(v for v in parts if v is not None) if any(v is not None for v in parts) else None
    # joint_source_RHS is a subset of source_RHS and uses the full 576-cell
    # system. It is auxiliary to the physical source solve but remains part of
    # the total full-system RHS; include it exactly once.
    components = (physical_source, joint, adjoint, receiver, tangent)
    known = [v for v in components if v is not None]
    # Missing components do not mean zero. For an explicitly returned task
    # count dictionary, absent operation classes are zero by counter contract.
    # Standalone probes may have unknown classes; preserve that distinction.
    total = sum(known) if known else numeric(counts.get("actual_complex_RHS"))
    return {"physical_LU": numeric(counts.get("LU", counts.get("physical_LU", counts.get("full_LU")))),
            "physical_source_RHS": physical_source,
            "adjoint_RHS": adjoint, "receiver_RHS": receiver,
            "tangent_RHS": tangent, "full_complex_RHS": total,
            "joint_source_RHS_auxiliary": joint,
            "calibration_adjoint_RHS_subset": numeric(counts.get("calibration_adjoint_RHS")),
            "joint_cholesky_auxiliary": numeric(counts.get("joint_cholesky")),
            "source_Gram_Cholesky_auxiliary": numeric(counts.get("source_Gram_Cholesky")),
            "Schur_RHS_auxiliary": numeric(counts.get("Schur_RHS")),
            "Gram_factorizations_auxiliary": numeric(counts.get("Gram_factorizations", counts.get("Gram")))}


def task_ledger(stage: str, task: dict, results_dir: Path) -> dict:
    path = results_dir / (task["task_id"] + ".json")
    receipt = read(path) if path.exists() else None
    result = (receipt or {}).get("result") or {}
    counts = result.get("counts", result.get("costcounts"))
    wall = numeric((receipt or {}).get("wall_s", result.get("task_wall_s", result.get("wall_seconds"))))
    cap = numeric((task.get("config") or {}).get("wall_seconds"))
    status = "MISSING" if receipt is None else receipt.get("status", "UNKNOWN")
    return {"category": "INFERENCE_TASK", "stage": stage,
            "task_id": task["task_id"], "case_id": task.get("case_id"),
            "observation_id": task.get("observation_id"),
            "method": task.get("method"), "regime": task.get("regime"),
            "K_run": task.get("K_run"), "start_index": task.get("start_index"),
            "p": task.get("p"), "layout": task.get("layout"),
            "status": status, "solver_status": result.get("wrapper_status", result.get("status")),
            "stop_reason": result.get("wrapper_stop_reason", result.get("stop_reason")),
            "source_path": str(path) if receipt is not None else None,
            "costs_status": "RECORDED" if isinstance(counts, dict) else "UNKNOWN",
            "wall_seconds": wall, "nominal_wall_cap_seconds": cap,
            "wall_cap_overshoot_seconds": max(0., wall-cap) if wall is not None and cap is not None else None,
            "batch_wall_seconds": None,
            **canonical_counts(counts)}


def setup_ledger(stage: str, results_dir: Path) -> dict:
    path = results_dir / "SUMMARY.json"
    if path.exists():
        summary = read(path)
        counts = summary.get("setup_counts")
        status = "RECORDED" if isinstance(counts, dict) else "UNKNOWN"
        wall = numeric(summary.get("setup_wall_s"))
        batch = numeric(summary.get("wall_s"))
    else:
        counts, status, wall, batch = None, "MISSING", None, None
    return {"category": "SHARED_SETUP", "stage": stage, "task_id": None,
            "case_id": None, "observation_id": None, "method": None,
            "regime": None, "K_run": None, "start_index": None, "p": None,
            "layout": None, "status": status, "solver_status": None,
            "stop_reason": None, "source_path": str(path) if path.exists() else None,
            "costs_status": status, "wall_seconds": wall,
            "nominal_wall_cap_seconds": None, "wall_cap_overshoot_seconds": None,
            "batch_wall_seconds": batch, **canonical_counts(counts)}


def auxiliary_row(name: str, source: str, counts: dict | None, wall,
                  status: str = "RECORDED", note: str = "") -> dict:
    return {"category": name, "stage": "development_auxiliary",
            "task_id": None, "case_id": None, "observation_id": None,
            "method": None, "regime": None, "K_run": None, "start_index": None,
            "p": None, "layout": None, "status": status,
            "solver_status": None, "stop_reason": None,
            "source_path": source, "costs_status": "RECORDED" if counts else "UNKNOWN",
            "wall_seconds": numeric(wall), "nominal_wall_cap_seconds": None,
            "wall_cap_overshoot_seconds": None, "batch_wall_seconds": None,
            "note": note, **canonical_counts(counts)}


def auxiliary_ledger() -> list[dict]:
    rows = []
    stage3 = read(HERE / "data_contract/stage3_dev/manifest.json")
    rows.append(auxiliary_row("PHYSICAL_DATA_GENERATION", "data_contract/stage3_dev/manifest.json",
                              stage3.get("generation_costcounts"), stage3.get("generation_wall_seconds"),
                              note="Stage 3 disclosed-object physical generation; 16 LU, 180 source RHS"))
    # This generation receipt was supplied at frozen Stage 4 generation close.
    # Never read private truth or mapping to reconstruct it.
    rows.append(auxiliary_row("PHYSICAL_DATA_GENERATION", "Stage4 generation completion message, 2026-10-30",
                              {"LU": 14, "source_RHS": 168}, .2579298331402242,
                              note="Stage 4 fresh observation generation, 14 LU/168 RHS; private truth unopened"))
    sources = [
        ("SCHUR_GEOMETRY", RESULTS / "schur_geometry.json", "costs", "diagnostic_wall_s"),
        ("CONTRAST_BIAS_GEOMETRY", RESULTS / "contrast_bias_geometry.json", "costs", "wall_s"),
        ("ILLUMINATION_COMPATIBILITY", RESULTS / "illumination_current_compatibility.json", "total_costs", "total_wall_seconds"),
        ("STAGE1_VALIDATOR", AUDIT / "VALIDATED_STAGE1.json", "costs", "wall_s"),
    ]
    for name, path, cost_key, wall_key in sources:
        if path.exists():
            data = read(path)
            rows.append(auxiliary_row(name, str(path.relative_to(HERE)), data.get(cost_key), data.get(wall_key),
                                      note="Separate CPU/validation work; count once"))
        else:
            rows.append(auxiliary_row(name, str(path.relative_to(HERE)), None, None,
                                      status="PENDING", note="Expected auxiliary receipt unavailable"))
    stage0 = AUDIT / "STAGE0_MATH.json"
    rows.append(auxiliary_row("SMALL_MATRIX", str(stage0.relative_to(HERE)), None, None,
                              status=read(stage0).get("status", "UNKNOWN") if stage0.exists() else "PENDING",
                              note="Algebra-only audit; Maxwell LU/RHS and wall NOT_MEASURED"))
    return rows


def delivery_rows(stage1: dict, stage2: dict, stage3: dict, stage4: dict) -> list[dict]:
    rows = []
    for r in stage1["rows"]:
        rows.append({"stage": "stage1", "case_id": r.get("case_id"), "method": r.get("method"),
                     "status": r.get("status"), "NRMSE": r.get("NRMSE"),
                     "IoU_Otsu": r.get("IoU_Otsu"), "relative_residual": r.get("relative_residual"),
                     "joint_gate": "PASS" if r.get("joint_pass") else "FAIL"})
    for g in stage2["groups"]:
        rows.append({"stage": "stage2", "case_id": g["case_id"], "method": g["method"],
                     "regime": g["regime"], "K": g["K"],
                     "status": "MISSING" if g["NOT_RUN_starts"] == g["planned_starts"] else "PARTIAL" if g["NOT_RUN_starts"] else "RECORDED",
                     "NRMSE": g["online_top1_NRMSE"],
                     "IoU_Otsu": g["online_top1_IoU_Otsu"],
                     "relative_residual": g["online_top1_relative_residual"],
                     "joint_gate": g["online_top1_joint_gate"],
                     "any_joint_at_K_OFFLINE": g["any_joint_at_K_OFFLINE"],
                     "oracle_best_at_K_NRMSE_OFFLINE": (g["oracle_best_at_K_OFFLINE"] or {}).get("NRMSE"),
                     "planned_starts": g["planned_starts"], "finite_starts": g["finite_starts"]})
    for r in stage3["tasks"]:
        rows.append({"stage": "stage3", "case_id": r.get("case_id"), "method": r.get("method"),
                     "p": r.get("p"), "layout": r.get("layout"),
                     "status": r.get("transport_status"), "NRMSE": r.get("NRMSE"),
                     "IoU_Otsu": r.get("IoU_Otsu"), "relative_residual": r.get("relative_residual"),
                     "joint_gate": r.get("joint_gate")})
    for r in stage4["rows"]:
        online = r.get("online") or {}
        rows.append({"stage": "stage4", "case_id": r.get("group_id"),
                     "observation_id": r.get("observation_id"), "method": r.get("method"),
                     "noise": r.get("noise"), "layout": r.get("source_layout"),
                     "status": "RECORDED" if online else "MISSING_OR_NONFINITE",
                     "NRMSE": online.get("NRMSE"), "IoU_Otsu": online.get("IoU_Otsu"),
                     "relative_residual": online.get("relative_data_residual"),
                     "joint_gate": online.get("joint_gate", "NOT_RUN"),
                     "planned_starts": r.get("planned_starts"), "finite_starts": r.get("finite_starts"),
                     "oracle_best_at_K_NRMSE_OFFLINE": (r.get("offline_best_at_K_ORACLE") or {}).get("NRMSE")})
    return rows


def actual_task_gates(stage1: dict, stage2: dict, stage3: dict,
                      stage4: dict, probe_scores: list[dict], ledger: list[dict]) -> dict:
    """Join scored endpoints to executed task IDs, never to K-prefix views."""
    score_by_id = {}
    manifests = {stage: read(HERE / "data_contract" / f"{stage}_remote.json")["tasks"]
                 for stage in STAGES}
    stage1_keys = {(r["case_id"], r["method"]): r for r in stage1["rows"]}
    stage4_keys = {(r["observation_id"], r["method"]): r for r in stage4["rows"]}
    if len(stage1_keys) != len(stage1["rows"]) or len(stage4_keys) != len(stage4["rows"]):
        raise ValueError("Stage 1/4 scored task keys are not unique")
    for task in manifests["stage1"]:
        r = stage1_keys[(task["case_id"], task["method"])]
        score_by_id[("stage1", task["task_id"])] = "PASS" if r["joint_pass"] else "FAIL"
    for stage, aggregate in (("stage2", stage2), ("stage3", stage3)):
        for r in aggregate["tasks"]:
            score_by_id[(stage, r["task_id"])] = r["joint_gate"]
    for task in manifests["stage4"]:
        r = stage4_keys[(task["observation_id"], task["method"])]
        score_by_id[("stage4", task["task_id"])] = (r.get("online") or {}).get("joint_gate", "NOT_RUN")
    for r in probe_scores:
        score_by_id[(r["stage"], r["task_id"])] = r["joint_gate"]
    task_rows = [r for r in ledger if r["category"] == "INFERENCE_TASK"]
    ids = [(r["stage"], r["task_id"]) for r in task_rows]
    if len(ids) != len(set(ids)) or set(ids) != set(score_by_id):
        raise ValueError("Executed task IDs and scored task IDs do not match exactly")
    by_stage = {}
    for stage in STAGES + ("schur_probe", "cie_probe"):
        stage_rows = [r for r in task_rows if r["stage"] == stage]
        if not stage_rows:
            continue
        gates = [score_by_id[(r["stage"], r["task_id"])] for r in stage_rows]
        by_stage[stage] = {"planned_tasks": len(gates), "joint_PASS": gates.count("PASS"),
                           "joint_FAIL": gates.count("FAIL"), "joint_NOT_RUN": gates.count("NOT_RUN")}
    gates = [score_by_id[ident] for ident in ids]
    return {"unit": "one executed inference task ID, including each probe task",
            "planned_tasks": len(gates), "joint_PASS": gates.count("PASS"),
            "joint_FAIL": gates.count("FAIL"), "joint_NOT_RUN": gates.count("NOT_RUN"),
            "by_stage": by_stage,
            "failed_task_ids": [{"stage": stage, "task_id": ident} for stage, ident in ids
                                if score_by_id[(stage, ident)] == "FAIL"]}


def _old_probe_truth() -> tuple[dict, object]:
    """Only the already released R4 developmental truth; never Stage 4 private."""
    r4 = HERE.parent / "research_z_alg_r4_diversity"
    cases = read(HERE / "data_contract/stage1_cases.json")["cases"]
    truth = {}
    for case in cases:
        ident = Path(case["public_path"]).stem
        with np.load(r4 / "truth_after_freeze" / f"{ident}.npz", allow_pickle=False) as f:
            m = np.asarray(f["m"], float).ravel()
            truth[case["case_id"]] = (m, np.asarray(f["mask"], bool).ravel(),
                                      float(case["amplitude"]), float(m.min()))
    spec = importlib.util.spec_from_file_location("r4_probe_metrics", r4 / "score_results.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return truth, module.metrics


def _score_probe_task(task: dict, receipt: dict | None, truth: dict, metrics) -> dict:
    result = (receipt or {}).get("result") or {}
    status = "MISSING" if receipt is None else receipt.get("status", "UNKNOWN")
    row = {"stage": "probe", "case_id": task.get("case_id"),
           "method": task.get("method"), "task_id": task["task_id"],
           "status": status, "NRMSE": None, "IoU_Otsu": None,
           "relative_residual": result.get("relative_residual"),
           "data_gate": "NOT_RUN" if status == "MISSING" else "FAIL",
           "material_gate": "NOT_RUN" if status == "MISSING" else "FAIL",
           "joint_gate": "NOT_RUN" if status == "MISSING" else "FAIL",
           "cost": canonical_counts(result.get("counts", result.get("costcounts"))),
           "wall_seconds": (receipt or {}).get("wall_s"),
           "stop_reason": result.get("wrapper_stop_reason", result.get("stop_reason"))}
    material = result.get("m")
    if material is not None:
        try:
            m = np.asarray(material, float).ravel()
            finite = m.size == 576 and bool(np.isfinite(m).all())
        except (TypeError, ValueError):
            finite = False
        if finite:
            t, mask, A, bg = truth[task["case_id"]]
            row.update(metrics(m, t, mask, A, bg))
            residual = numeric(row["relative_residual"])
            data_ok = residual is not None and residual <= .001
            material_ok = row["NRMSE"] <= .15 and row["IoU_Otsu"] >= .7
            row["data_gate"] = "PASS" if data_ok else "FAIL"
            row["material_gate"] = "PASS" if material_ok else "FAIL"
            row["joint_gate"] = "PASS" if data_ok and material_ok else "FAIL"
    return row


def optional_probe_rows(paths: list[Path]) -> tuple[list[dict], list[dict], list[dict]]:
    ledger, summaries, scored = [], [], []
    truth = metrics = None
    for path in paths:
        if not path.exists():
            summaries.append({"path": str(path), "status": "PENDING_OPTIONAL"})
            continue
        if path.is_dir():
            manifest_path = HERE / "data_contract" / f"{path.name}_remote.json"
            if not manifest_path.exists():
                raise FileNotFoundError(f"Probe directory requires frozen manifest: {manifest_path}")
            manifest = read(manifest_path)
            tasks = manifest["tasks"]
            if truth is None:
                truth, metrics = _old_probe_truth()
            ledger.extend(task_ledger(path.name, task, path) for task in tasks)
            ledger.append(setup_ledger(path.name, path))
            rows = []
            for task in tasks:
                receipt_path = path / (task["task_id"] + ".json")
                receipt = read(receipt_path) if receipt_path.exists() else None
                row = _score_probe_task(task, receipt, truth, metrics)
                row["stage"] = path.name
                rows.append(row)
            scored.extend(rows)
            summaries.append({"path": str(path), "status": "RECORDED" if (path / "SUMMARY.json").exists() else "PENDING",
                              "planned_tasks": len(tasks),
                              "received_tasks": sum(r["status"] != "MISSING" for r in rows),
                              "joint_pass": sum(r["joint_gate"] == "PASS" for r in rows),
                              "batch_wall_seconds": read(path / "SUMMARY.json").get("wall_s") if (path / "SUMMARY.json").exists() else None})
            continue
        data = read(path)
        probe_result = data.get("result") if isinstance(data.get("result"), dict) else data
        counts = probe_result.get("counts", probe_result.get("costcounts", probe_result.get("costs")))
        category = "STAGE4_VALIDATOR" if path.name == "VALIDATED_STAGE4.json" else "OPTIONAL_PROBE"
        ledger.append(auxiliary_row(category, str(path), counts,
                                    data.get("wall_s", probe_result.get("wall_s", probe_result.get("wall_seconds"))),
                                    status=data.get("status", "RECORDED"),
                                    note="Separate validator/probe; counted once, no new solve by this merger"))
        summaries.append({"path": str(path), "status": data.get("status", "RECORDED"),
                          "result_fields": [k for k in data if k not in ("m", "truth", "material")]})
    return ledger, summaries, scored


def prerequisites(probes: list[Path]) -> list[str]:
    required = [RESULTS / "stage1_scores.json", RESULTS / "STAGE2_AGGREGATE.json",
                RESULTS / "STAGE3_AGGREGATE.json", RESULTS / "STAGE4_SCORES.json"]
    required += [RESULTS / stage / "SUMMARY.json" for stage in STAGES]
    required += [path / "SUMMARY.json" for path in probes if path.exists() and path.is_dir()]
    return [str(path) for path in required if not path.exists()]


def finalize(probes: list[Path]) -> dict:
    probes = list(dict.fromkeys(path.resolve() for path in probes))
    missing = prerequisites(probes)
    if missing:
        raise RuntimeError("PENDING required inputs; no final delivery written: " + ", ".join(missing))
    stage1 = read(RESULTS / "stage1_scores.json")
    stage2 = read(RESULTS / "STAGE2_AGGREGATE.json")
    stage3 = read(RESULTS / "STAGE3_AGGREGATE.json")
    stage4 = read(RESULTS / "STAGE4_SCORES.json")
    if stage4.get("schema") != "fbg-r2-blind-scores-v2" or not stage4.get("truth_disclosed_only_after_seal"):
        raise ValueError("Stage 4 sealed scorer output required; do not open private truth")
    ledger = []
    stage_counts = {}
    for stage in STAGES:
        manifest = read(HERE / "data_contract" / f"{stage}_remote.json")
        results_dir = RESULTS / stage
        task_rows = [task_ledger(stage, task, results_dir) for task in manifest["tasks"]]
        ledger.extend(task_rows)
        setup = setup_ledger(stage, results_dir)
        ledger.append(setup)
        stage_counts[stage] = {"planned_tasks": len(task_rows),
                               "received_tasks": sum(r["status"] != "MISSING" for r in task_rows),
                               "batch_wall_seconds": setup["batch_wall_seconds"],
                               "shared_setup_wall_seconds": setup["wall_seconds"],
                               "charged_worker_wall_seconds_sum_known": sum(r["wall_seconds"] or 0 for r in task_rows),
                               "physical_LU_sum_known": sum(r["physical_LU"] or 0 for r in task_rows),
                               "full_complex_RHS_sum_known": sum(r["full_complex_RHS"] or 0 for r in task_rows),
                               "unknown_cost_tasks": sum(r["costs_status"] == "UNKNOWN" for r in task_rows)}
    ledger.extend(auxiliary_ledger())
    probe_ledger, probe_summary, probe_scores = optional_probe_rows(probes)
    ledger.extend(probe_ledger)
    write(RESULTS / "DEVELOPMENT_PROBES.json", {
        "schema": "fbg-r2-development-probes-v1",
        "scope": "Previously disclosed R4 truth_after_freeze only; no fresh truth opened",
        "probes": probe_summary, "task_rows": probe_scores,
        "planned_probe_tasks": len(probe_scores),
        "joint_pass": sum(r["joint_gate"] == "PASS" for r in probe_scores),
        "joint_fail": sum(r["joint_gate"] == "FAIL" for r in probe_scores),
        "NOT_RUN": sum(r["joint_gate"] == "NOT_RUN" for r in probe_scores)})
    columns = ["category", "stage", "task_id", "case_id", "observation_id", "method",
               "regime", "K_run", "start_index", "p", "layout", "status", "solver_status",
               "stop_reason", "costs_status", "wall_seconds", "nominal_wall_cap_seconds",
               "wall_cap_overshoot_seconds", "batch_wall_seconds", "physical_LU",
               "physical_source_RHS", "adjoint_RHS", "receiver_RHS", "tangent_RHS",
               "full_complex_RHS", "joint_source_RHS_auxiliary", "calibration_adjoint_RHS_subset",
               "joint_cholesky_auxiliary", "source_Gram_Cholesky_auxiliary",
               "Schur_RHS_auxiliary", "Gram_factorizations_auxiliary", "source_path", "note"]
    csv_write(AUDIT / "COST_LEDGER.csv", ledger, columns)
    by_case = delivery_rows(stage1, stage2, stage3, stage4)
    by_case.extend(probe_scores)
    case_columns = ["stage", "case_id", "observation_id", "method", "regime", "K", "p",
                    "layout", "noise", "status", "NRMSE", "IoU_Otsu", "relative_residual",
                    "joint_gate", "planned_starts", "finite_starts", "any_joint_at_K_OFFLINE",
                    "oracle_best_at_K_NRMSE_OFFLINE"]
    csv_write(RESULTS / "RESULTS_BY_CASE.csv", by_case, case_columns)
    execution_failures = [r for r in ledger if r["category"] == "INFERENCE_TASK" and r["status"] != "OK"]
    scientific_failures = [r for r in by_case if r.get("joint_gate") == "FAIL"]
    pending = [r for r in ledger if r["status"] in ("MISSING", "PENDING", "UNKNOWN")]
    actual_gates = actual_task_gates(stage1, stage2, stage3, stage4, probe_scores, ledger)
    write(HERE / "FAILURES.json", {"schema": "fbg-r2-failures-v1",
                                  "execution_failures": execution_failures,
                                  "scientific_gate_failures_group_views": scientific_failures,
                                  "scientific_gate_failed_actual_task_ids": actual_gates["failed_task_ids"],
                                  "pending_or_unknown_cost_rows": pending,
                                  "interpretation": "Case-table group views include overlapping Stage 2 K prefixes; their 257 FAIL rows are not 257 executed task failures. Actual task gates use unique inference task IDs. Gate FAIL and execution failure are distinct."})
    summary = {"schema": "fbg-r2-results-summary-v1", "status": "RECEIPTS_MERGED_OFFLINE",
               "stages": stage_counts,
               "development_case_counts": {"stage1": 5, "stage2": 3, "stage3": 2},
               "independent_shape_clusters": {"stage1": 3, "stage2": 2, "stage3": 2,
                                              "note": "Stage 1 has five amplitude/case datasets but only ring, curved bar, and digit 8 shape clusters; Stage 2 reuses two of those shapes. These are not iid shape draws."},
               "fresh_independent_objects_stage4": 14,
               "stage4_observations": stage4.get("planned_observations"),
               "stage4_closed_seal": stage4.get("seal"),
               "stage4_by_method": stage4.get("by_method"),
               "stage2_groups": len(stage2.get("groups", [])),
               "stage3_paired_acquisitions": len(stage3.get("paired_by_acquisition", [])),
               "optional_probes": probe_summary,
               "ledger_rows": len(ledger), "by_case_rows": len(by_case),
               "execution_failures": len(execution_failures),
               "actual_task_joint_gates": {k: v for k, v in actual_gates.items() if k != "failed_task_ids"},
               "scientific_gate_failures_group_views": len(scientific_failures),
               "pending_or_unknown_cost_rows": len(pending),
               "task_counter_coverage": {"inference_task_rows": actual_gates["planned_tasks"],
                                         "unknown_task_counter_rows": sum(r["category"] == "INFERENCE_TASK" and r["costs_status"] == "UNKNOWN" for r in ledger)},
               "unmeasured_cost_scope": "Stage 0 small-matrix wall is NOT_MEASURED; import, launch, SSH/transfer, and SM-level allocation costs are outside recorded task and batch counters. Do not describe all costs as known.",
               "cost_rules": [
                   "Full complex RHS = physical source + adjoint + receiver + tangent; calibration_adjoint is a subset of adjoint and excluded from sum.",
                   "For joint probes, joint_source_RHS is a subset of source_RHS: physical_source=source-joint_source, while total full-system RHS includes joint_source once. Auxiliary Cholesky and Schur RHS are separate.",
                   "Only physical Maxwell LU are in physical_LU; small Gram/Cholesky are auxiliary.",
                   "Per-task wall sums are charged solver work and must not replace recorded batch wall; shared setup is separate.",
                   "Nominal wall-cap overshoot is reported; unknown or unrecovered costs remain null.",
                   "Stage 4 scalar labels/material metrics come only from the already sealed scorer output; this script never opens fresh truth."],
               "scientific_interpretation_owner": "parent Codex final review"}
    write(RESULTS / "RESULTS_SUMMARY.json", summary)
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", type=Path, action="append", default=[],
                    help="Optional Schur/CIE probe directory or validator/probe JSON; repeat as needed")
    args = ap.parse_args()
    paths = [path if path.exists() else HERE / path for path in args.probe]
    result = finalize(paths)
    print(json.dumps({"status": result["status"], "ledger_rows": result["ledger_rows"],
                      "by_case_rows": result["by_case_rows"],
                      "execution_failures": result["execution_failures"]}))


if __name__ == "__main__":
    main()
