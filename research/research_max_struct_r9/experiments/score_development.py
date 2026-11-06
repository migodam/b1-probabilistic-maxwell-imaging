"""Public seal, then offline grading of the fixed four-case R9 development set.

No solver imports/calls. `seal` never reads truth; `score` verifies every stamped
public artifact before opening exactly the four selected truth arrays. Failed
terminal attempts remain rows. Nothing executes when this module is imported.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import importlib.util
import json
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
CASES = ("obs_b20d5a55677d4eca", "obs_40ac5b79e8a74ac3",
         "r8_test_0000", "r8_test_0027")
ARMS = ("cie_single", "cie_dual", "global_dual", "a_material", "a_plus",
        "a_plus_sourcewise", "physical_lbfgs")
SUCCESS = {"SEALED_PUBLIC_PREDICTION", "COMPLETED_PUBLIC_PREDICTION"}
TERMINAL = SUCCESS | {"FAILED", "FAILED_R9_ADAPTER"}
DEFAULT_MANIFESTS = ("baseline_development.json", "baseline_historical_fixed.json",
                     "landscape_development.json", "baseline_cost_controls.json")


def load(path):
    return json.loads(Path(path).read_text())


def resolve(path, base=REPO):
    p = Path(path).expanduser()
    return (p if p.is_absolute() else base / p).resolve()


def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as f:
        json.dump(value, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write("\n")


def stamp(path):
    p = Path(path).resolve()
    s = p.stat()
    return dict(path=str(p), bytes=s.st_size, mtime_ns=s.st_mtime_ns)


def public_path(path, base):
    p = resolve(path, base)
    if {"sealed_truth", "data_private", "private", "truth", "gt", "sealed_gt"} & set(p.parts):
        raise ValueError("Private input in public task: " + str(p))
    return p


def public_arrays(files, base):
    """Read only public operators/observations, not labels or initializers."""
    data = {}
    forbidden = {"truth", "alpha", "mask", "yclean", "y_clean", "seed", "family",
                 "gt", "ground_truth", "m_true", "m_truth", "true_material"}
    wanted = {"G", "E", "E0", "S", "Y", "k", "radius", "coords", "centers", "sigma"}
    for file in files:
        with np.load(public_path(file, base), allow_pickle=False) as z:
            if forbidden & {k.lower() for k in z.files}:
                raise ValueError("Private payload in public data: " + str(file))
            for key in wanted & set(z.files):
                canonical = {"E0": "E", "centers": "coords"}.get(key, key)
                if canonical in data and not np.array_equal(data[canonical], z[key]):
                    raise ValueError("Conflicting public field: " + canonical)
                data[canonical] = z[key].copy()
    if not {"G", "E", "S", "Y", "k", "radius", "coords"} <= data.keys():
        raise ValueError("Incomplete public physical contract")
    return data


def same_public(a, b):
    return a.keys() == b.keys() and all(np.array_equal(a[k], b[k]) for k in a)


def metric_module(r8):
    # These pure definitions are exactly those imported by the actual R8 scorer.
    path = r8 / "research/research_max_struct_r7/experiments/score_primary.py"
    spec = importlib.util.spec_from_file_location("r9_offline_r8_metrics", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, path


def frozen_tasks(paths, r8):
    tasks = []
    for path in paths:
        manifest = load(path)
        for source in manifest["tasks"]:
            task = dict(source)
            if task["case_id"] not in CASES:
                raise ValueError("Task outside fixed four-case development set")
            task["_manifest"] = str(path)
            task["_kind"] = "landscape" if "algorithm" in task else "baseline"
            task["_input_base"] = str(resolve(task.get("input_base", manifest.get("input_base", r8))))
            tasks.append(task)
    ids = [t["task_id"] for t in tasks]
    if len(ids) != len(set(ids)) or {t["case_id"] for t in tasks} != set(CASES):
        raise ValueError("Duplicate task or missing frozen development case")
    return tasks


def task_receipt(task):
    output = resolve(task["output_dir"])
    path = output / "RECEIPT.json"
    if not path.exists():
        path = output.with_name(output.name + "__adapter_inputs") / "ADAPTER_PREPARATION.json"
        if not path.exists() or load(path).get("status") != "FAILED_R9_ADAPTER":
            raise RuntimeError("PENDING: " + task["task_id"])
    rec = load(path)
    if rec.get("status") not in TERMINAL:
        raise RuntimeError("PENDING: nonterminal " + task["task_id"])
    if rec.get("truth_access") != "NONE":
        raise ValueError("Prediction receipt does not declare public-only access")
    if path.name == "RECEIPT.json":
        rt = rec.get("task", {})
        for key in ("task_id", "case_id", "algorithm" if task["_kind"] == "landscape" else "arm"):
            if rec.get(key, rt.get(key)) != task[key]:
                raise ValueError("Receipt identity mismatch: " + task["task_id"])
        if task["_kind"] == "landscape" and rec.get("config") != load(task["_manifest"])["landscape_config"]:
            raise ValueError("Changed frozen landscape configuration")
    return output, path, rec


def endpoint_entry(task, output, receipt_path, rec, origin):
    files = [receipt_path]
    if rec["status"] in SUCCESS:
        files.append(output / "prediction.npz")
        with np.load(files[-1], allow_pickle=False) as z:
            if z["m"].size != 576 or not np.isfinite(z["m"]).all() or not np.isfinite(z["Y_prediction"]).all():
                raise ValueError("Invalid sealed prediction: " + task["task_id"])
    # Paid stage/cost outputs are provenance only; never use GT to select stages.
    files += sorted(output.glob("stage_*_joint.json"))
    for name in ("PUBLIC_TASK.json", "FROZEN_CONFIG.json", "NUMERICAL_SCHEDULE.json",
                 "AUXILIARY_COST.json", "R9_DPCS_RULE.json", "CORRECTOR.json",
                 "BUDGET_STATES.npz", "R9_BUDGET_STATES.npz", "BUDGET_MATERIALS.npz", "BUDGET.json"):
        if (output / name).is_file():
            files.append(output / name)
    return dict(task=task, origin=origin, status=rec["status"], output=str(output),
                receipt_path=str(receipt_path), files=[stamp(p) for p in files])


def inherited_entries(r8, public_by_case):
    """Only real R8 uniform endpoints from explicit manifests; no R6/R7 relabel."""
    study = r8 / "research/research_max_struct_r8"
    manifests = [study / "configs/warm_all_manifest.json",
                 study / "experiments/history_warm_manifest.json"]
    entries, available, evidence, seen = [], [], [], set()
    for path in manifests:
        if not path.exists():
            continue
        evidence.append(path)
        for t in load(path)["tasks"]:
            ident = t["case_id"]
            if ident not in CASES or t.get("arm") not in ARMS:
                continue
            if t.get("initialization") not in ("uniform", "uniform0p1") or t.get("initialization_file"):
                continue
            output = resolve(t["output_dir"], r8)
            if output in seen:
                continue
            seen.add(output)
            rp = output / "RECEIPT.json"
            info = dict(case_id=ident, arm=t["arm"], manifest=str(path), output=str(output))
            if not rp.exists():
                available.append(dict(info, status="NOT_AVAILABLE_MISSING_RECEIPT")); continue
            rec = load(rp)
            if rec.get("status") != "SEALED_PUBLIC_PREDICTION":
                available.append(dict(info, status="NOT_AVAILABLE_NOT_SEALED", receipt_status=rec.get("status"))); continue
            rt = rec.get("task", {})
            if (rec.get("task_id"), rec.get("case_id"), rec.get("arm")) != (t["task_id"], ident, t["arm"]):
                raise ValueError("Inherited receipt identity mismatch")
            cfg = rec["configuration"]
            uniform = float(t.get("initial_material", cfg["common"]["host"]))
            if (uniform != .1 or bool(t.get("initialize_physical_current", False))
                    or cfg.get("aplus_options", {}).get("current_anchor")
                    or not same_public(public_arrays(rt["public_files"], r8), public_by_case[ident])):
                available.append(dict(info, status="NOT_AVAILABLE_PUBLIC_OR_INITIALIZATION_MISMATCH")); continue
            task = dict(t, _kind="baseline", _manifest=str(path), _input_base=str(r8))
            entry = endpoint_entry(task, output, rp, rec, "REUSED_R8_UNIFORM_ENDPOINT")
            entry["files"] += [stamp(public_path(p, r8)) for p in t["public_files"]]
            entries.append(entry)
            available.append(dict(info, status="REUSED_EXACT_PUBLIC_INPUTS_UNIFORM0P1", configuration=cfg))
    for ident in CASES:
        for arm in ARMS:
            if not any(e["task"]["case_id"] == ident and e["task"]["arm"] == arm for e in entries):
                available.append(dict(case_id=ident, arm=arm, status="NOT_AVAILABLE_NO_EXACT_R8_UNIFORM_ENDPOINT",
                    note="No R6/R7 result substituted; historical same-code fixed controls are new R9 tasks"))
    return entries, available, evidence


def seal(args):
    r8 = resolve(args.r8_repo)
    paths = [resolve(p) for p in args.manifest]
    tasks = frozen_tasks(paths, r8)
    evidence = set(paths)
    metrics_path = r8 / "research/research_max_struct_r7/experiments/score_primary.py"
    evidence.add(metrics_path)
    entries, publics = [], {}
    # All requested tasks must complete or fail before any truth is opened.
    for task in tasks:
        output, rp, rec = task_receipt(task)
        entries.append(endpoint_entry(task, output, rp, rec, "NEW_R9_" + task["_kind"].upper()))
        base = Path(task["_input_base"])
        pub = public_arrays(task["public_files"], base)
        ident = task["case_id"]
        if ident in publics and not same_public(pub, publics[ident]):
            raise ValueError("Different public acquisition for one development case")
        publics[ident] = pub
        evidence.update(public_path(f, base) for f in task["public_files"])
        if task.get("config"):
            evidence.add(resolve(task["config"]))
        if task.get("initialization_file"):
            evidence.add(public_path(task["initialization_file"], base))
    inherited, availability, inherited_manifests = inherited_entries(r8, publics)
    evidence.update(inherited_manifests)
    for p in (r8 / "research/research_max_struct_r6/data_public/cases.json",
              r8 / "research/research_max_struct_r8/data_public/warm_test/test_manifest.json"):
        evidence.add(p)
    # Preserve archived failed engineering costs without importing queue modules.
    archive = ROOT / "delegated/development_queue/ARCHIVED_FAILED_OUTPUTS.json"
    failed_files = []
    if archive.exists():
        evidence.add(archive)
        for item in load(archive).get("maps", []):
            p = resolve(item["preserved_path"]) / "RECEIPT.json"
            if p.is_file():
                failed_files.append(str(p)); evidence.add(p)
    failed_queue = ROOT / "delegated/development_queue/RECEIPT.json"
    if failed_queue.exists():
        q = load(failed_queue)
        preserve_queue = (q.get("status") == "COMPLETED_WITH_ERRORS"
                          and all(j.get("status") == "FAILED" for j in q.get("jobs", [])))
        if preserve_queue:
            evidence.add(failed_queue)
    else:
        preserve_queue = False
    result = dict(status="ALL_REQUESTED_R9_DEVELOPMENT_PUBLIC_OUTPUTS_SEALED",
        created_utc=datetime.now(timezone.utc).isoformat(), cases=list(CASES),
        manifests=[str(p) for p in paths], r8_repo=str(r8), requested_tasks=len(tasks),
        successful_predictions=sum(e["status"] in SUCCESS for e in entries),
        failed_tasks=sum(e["status"] not in SUCCESS for e in entries),
        entries=entries + inherited, inherited_availability=availability,
        failed_archive=str(archive) if archive.exists() else None, failed_receipts=failed_files,
        failed_queue_receipt=str(failed_queue) if preserve_queue else None,
        evidence=[stamp(p) for p in sorted(evidence)], truth_access="NONE", hash_checks="NONE")
    dump(resolve(args.seal), result)
    print(json.dumps(dict(status=result["status"], tasks=len(tasks), truth_access="NONE")))


def verify_seal(args):
    sealed = load(resolve(args.seal))
    if sealed["status"] != "ALL_REQUESTED_R9_DEVELOPMENT_PUBLIC_OUTPUTS_SEALED":
        raise ValueError("A complete public seal is required")
    if sealed["cases"] != list(CASES):
        raise ValueError("Wrong development case set")
    if [str(resolve(p)) for p in args.manifest] != sealed["manifests"] or str(resolve(args.r8_repo)) != sealed["r8_repo"]:
        raise ValueError("Score inputs differ from the public seal")
    tasks = frozen_tasks([Path(p) for p in sealed["manifests"]], Path(sealed["r8_repo"]))
    actual = [e for e in sealed["entries"] if e["origin"].startswith("NEW_R9_")]
    if len(actual) != len(tasks) or {e["task"]["task_id"] for e in actual} != {t["task_id"] for t in tasks}:
        raise ValueError("Seal does not cover all requested tasks")
    for st in sealed["evidence"] + [p for e in sealed["entries"] for p in e["files"]]:
        if stamp(st["path"]) != st:
            raise RuntimeError("Changed sealed artifact: " + st["path"])
    return sealed


def truth_records(r8):
    """Called only after verify_seal: open exactly four existing label arrays."""
    historical = load(r8 / "research/research_max_struct_r6/sealed_truth/scoring/historical_mapping_PRIVATE.json")
    historical = {r["case_id"]: r for r in historical}
    confirmation = load(r8 / "research/research_max_struct_r6/sealed_truth/confirmation/mapping_PRIVATE.json")
    confirmation = {r["case_id"]: r for r in confirmation["records"]}
    warm = load(r8 / "research/research_max_struct_r8/sealed_truth/warm_test/score_metadata.json")
    warm = {r["case_id"]: r for r in warm["records"]}
    oldpublic = load(r8 / "research/research_max_struct_r6/data_public/cases.json")
    oldpublic_rows = oldpublic.get("records", oldpublic.get("cases"))
    if not isinstance(oldpublic_rows, list):
        raise ValueError("R6 public cases must contain a records/cases list")
    oldpublic = {r["case_id"]: r for r in oldpublic_rows}
    warmpublic = load(r8 / "research/research_max_struct_r8/data_public/warm_test/test_manifest.json")
    warmpublic = {r["case_id"]: r for r in warmpublic["records"]}
    result = {}
    for ident in CASES:
        if ident.startswith("obs_"):
            meta = historical[ident] if ident == CASES[0] else confirmation[ident]
            p = resolve(meta["truth_path"], r8)
            amp = float(oldpublic[ident]["amplitude"])
            same_grid = meta.get("inverse_crime")
        else:
            p = r8 / f"research/research_max_struct_r8/sealed_truth/warm_test/{ident}.npz"
            amp = float(warm[ident]["amplitude"])
            same_grid = warmpublic[ident]["condition"]["generation_grid"] == 24
        with np.load(p, allow_pickle=False) as z:
            result[ident] = dict(m=z["m"].ravel().copy(), mask=z["mask"].ravel().astype(bool),
                amplitude=amp, truth_path=str(p), same_grid=same_grid)
    return result


def score_budgets(entry, receipt, truth, metrics):
    """Grade already-paid incumbents, preserving measured prefix limitations."""
    d, task = Path(entry["output"]), entry["task"]
    records = []
    if task["_kind"] == "landscape":
        path, metadata = d / "BUDGET_MATERIALS.npz", d / "BUDGET.json"
        if not path.exists() or not metadata.exists():
            return []
        raw = load(metadata)
        with np.load(path, allow_pickle=False) as z:
            materials = z["m"].copy()
        if len(raw) != len(materials):
            raise ValueError("Landscape budget arrays/metadata length mismatch")
        shift = receipt.get("pre_solver_CPU_seconds")
        semantics = "EXACT_RECORDED_PRE_SOLVER_PREFIX"
        if shift is None:
            total, inner = receipt.get("process_CPU_seconds"), receipt.get("cost", {}).get("CPU_seconds")
            shift = max(0., total-inner) if total is not None and inner is not None else None
            semantics = "UPPER_SHIFT_including_output_tail" if shift is not None else "UNMEASURED_PREFIX"
        for i, (r, material) in enumerate(zip(raw, materials)):
            c = r.get("physical_cost", {})
            records.append((material, dict(index=i, CPU_seconds=(r["CPU_seconds"]+shift) if shift is not None else None,
                raw_solver_CPU_seconds=r["CPU_seconds"], CPU_shift_seconds=shift, CPU_semantics=semantics,
                original_Phi=r["Phi"], LU=c.get("LU"), complex_RHS=c.get("actual_complex_RHS"),
                stage=r.get("label"), budget_source=str(path), metadata_source=str(metadata))))
    else:
        path = d / "R9_BUDGET_STATES.npz"
        if not path.exists():
            path = d / "BUDGET_STATES.npz"
        if not path.exists():
            return []
        with np.load(path, allow_pickle=False) as z:
            materials, costs = z["m"].copy(), z["costs"].copy()
            columns = [str(k) for k in z["columns"]]
        if len(materials) != len(costs) or costs.ndim != 2 or costs.shape[1] != 4:
            raise ValueError("Inherited budget array shape mismatch")
        for i, (material, c) in enumerate(zip(materials, costs)):
            records.append((material, dict(index=i, CPU_seconds=float(c[0]),
                CPU_semantics="EXISTING_RECORDED_PAID_PREFIX: " + columns[0],
                original_Phi=float(c[3]), LU=float(c[1]), complex_RHS=float(c[2]),
                budget_source=str(path))))
    return [dict(task_id=task["task_id"], case_id=task["case_id"], origin=entry["origin"],
                 method=task.get("algorithm", task.get("arm")), dpcs_rule=task.get("dpcs_rule"),
                 **row, **metrics.material_metrics(material, truth["m"], truth["mask"], truth["amplitude"]))
            for material, row in records]


def score(args):
    begin, cpu = time.perf_counter(), time.process_time()
    sealed = verify_seal(args)  # No private read above this line.
    r8 = Path(sealed["r8_repo"])
    metrics, metric_path = metric_module(r8)
    out = resolve(args.output)
    if out.exists():
        raise FileExistsError(out)
    truths = truth_records(r8)
    rows, stages, budgets = [], [], []
    for e in sealed["entries"]:
        t, d = e["task"], Path(e["output"])
        rec = load(e["receipt_path"]); ident = t["case_id"]; truth = truths[ident]
        landscape = t["_kind"] == "landscape"
        cost = (rec.get("cost") or {}).get("physical", {}) if landscape else (rec.get("physical_model_cost") or {})
        auxcost = (rec.get("cost") or {}).get("auxiliary", {}) if landscape else rec.get("auxiliary_execution_receipt", {}).get("operation_counts", {})
        stage_list = rec.get("stages", []) if landscape else rec.get("auxiliary_candidate_scores", [])
        corrector = next((s for s in reversed(stage_list) if s.get("label") == "original_corrector"), {}) if landscape else (rec.get("corrector") or {})
        prefix_cpu = rec.get("adapter_preparation_cpu_seconds")
        row = dict(task_id=t["task_id"], case_id=ident, origin=e["origin"], cohort="EXPOSED_DEVELOPMENT_ONLY",
            original_source_cohort=t.get("cohort"),
            method=t.get("algorithm", t.get("arm")), initialization=t.get("initialization", "uniform"),
            dpcs_rule=t.get("dpcs_rule"), dpcs_decision=rec.get("dpcs_rule"),
            status=e["status"], material_image_success=False,
            material_nrmse=None, otsu_iou=None, amplitude=truth["amplitude"],
            original_Phi=rec.get("Phi") if landscape else rec.get("physical_objective"),
            relative_residual=rec.get("relative_residual") if landscape else rec.get("physical_relative_residual"),
            cpu_seconds=rec.get("process_CPU_seconds") if landscape else rec.get("r9_total_cpu_seconds", rec.get("cpu_seconds")),
            wall_seconds=rec.get("process_wall_seconds") if landscape else rec.get("r9_total_wall_seconds", rec.get("wall_seconds")),
            nested_solver_cpu_seconds=(rec.get("cost") or {}).get("CPU_seconds") if landscape else rec.get("cpu_seconds"),
            adapter_preparation_cpu_seconds=prefix_cpu, LU=cost.get("LU"), complex_RHS=cost.get("actual_complex_RHS"),
            source_RHS=cost.get("source_RHS"), adjoint_RHS=cost.get("adjoint_RHS"), physical_model_cost=cost,
            auxiliary_cost=auxcost, corrector_stop=corrector.get("stop"), corrector_iterations=corrector.get("iterations"),
            corrector_success=corrector.get("success"), auxiliary_iterations=rec.get("auxiliary_iterations"),
            selected_stage=rec.get("selected_label") if landscape else rec.get("auxiliary_selected_stage"),
            selection_rule=rec.get("online_selection_rule", rec.get("selection")),
            receipt_path=e["receipt_path"], prediction_path=str(d / "prediction.npz") if e["status"] in SUCCESS else None,
            source_manifest=t["_manifest"], configuration=rec.get("config") if landscape else rec.get("configuration"),
            execution_platform=rec.get("execution_platform"), error=rec.get("error"),
            cost_pairing_role=("FRESH_LOCAL_FIXED_RANK_CONTROL" if e["origin"].startswith("NEW_R9_") and t.get("dpcs_rule")=="fixed_001"
                               else "FRESH_LOCAL_TASK" if e["origin"].startswith("NEW_R9_")
                               else "REUSED_HISTORICAL_R8_COST_NOT_FRESH_LOCAL_CONTROL"),
            cost_scope="Inclusive task CPU/wall; nested costs overlap. Reused R8 costs are historical, not new spending.")
        if e["status"] in SUCCESS:
            with np.load(d / "prediction.npz", allow_pickle=False) as z:
                m, yp = z["m"].ravel(), z["Y_prediction"]
            public = public_arrays(t["public_files"], Path(t["_input_base"]))
            row.update(metrics.material_metrics(m, truth["m"], truth["mask"], truth["amplitude"]))
            row.update(metrics.noise_metrics(yp, public, truth["same_grid"] is True))
            if truth["same_grid"] is None:
                row.update(noise_gate_scope="RAW NOISE-ONLY diagnostic; historical grid contract not asserted by mapping",
                           model_mismatch_budget="NOT_ASSERTED_FROM_METADATA")
            row["same_grid_material_noise_success"] = bool(row["material_image_success"] and row["noise_CN99_pass"]) if truth["same_grid"] is True else None
            row["independent_saved_Phi"] = .5 * row["relative_residual"]**2 + metrics.fixed_prior(m)
            row["Phi_receipt_difference"] = row["independent_saved_Phi"] - row["original_Phi"]
        rows.append(row)
        if e["status"] in SUCCESS:
            budgets.extend(score_budgets(e, rec, truth, metrics))
        for i, stage in enumerate(stage_list):
            stages.append(dict(task_id=t["task_id"], case_id=ident, origin=e["origin"], index=i,
                recorded_stage=stage, interpretation="Paid stage metadata only; no GT-based iterate selection"))
    failed = []
    archive = load(sealed["failed_archive"]) if sealed.get("failed_archive") else None
    for path in sealed["failed_receipts"]:
        r = load(path)
        failed.append(dict(receipt_path=path, task=r.get("task"), status=r.get("status"), error=r.get("error"),
            CPU_seconds=r.get("process_CPU_seconds", r.get("cpu_seconds")),
            wall_seconds=r.get("process_wall_seconds", r.get("wall_seconds")),
            physical_cost=r.get("physical_model_cost", (r.get("cost") or {}).get("physical"))))
    new = [r for r in rows if r["origin"].startswith("NEW_R9_")]
    ledger = dict(new_tasks=len(new), CPU_known_seconds=sum(r["cpu_seconds"] for r in new if r["cpu_seconds"] is not None),
        CPU_unmeasured_tasks=[r["task_id"] for r in new if r["cpu_seconds"] is None],
        task_wall_seconds_sum=sum(r["wall_seconds"] for r in new if r["wall_seconds"] is not None),
        wall_scope="Task wall sum is workload, not concurrent queue elapsed time", LU_known=sum(r["LU"] for r in new if r["LU"] is not None),
        complex_RHS_known=sum(r["complex_RHS"] for r in new if r["complex_RHS"] is not None),
        auxiliary_scope="Current Cholesky/RHS and receiver SVD remain separate auxiliary_cost fields, not full-wave LU/RHS",
        archived_failure_ledger=archive, archived_receipt_costs=failed,
        archived_queue_receipt=load(sealed["failed_queue_receipt"]) if sealed.get("failed_queue_receipt") else None,
        archived_measured_CPU_seconds=sum(r["CPU_seconds"] for r in failed if r["CPU_seconds"] is not None),
        failure_cost_limit="Archive index reports28 first import failures and0 physics; baseline child CPU was not measured. Known archived CPU is a lower bound, not0 total.",
        inherited_cost_scope="Original R8 receipt costs retained per row; do not charge again as new R9 solver spending")
    out.mkdir(parents=True)
    dump(out / "SCORES.json", dict(status="SCORED_FROZEN_FOUR_CASE_DEVELOPMENT", rows=rows,
        requested_tasks=sealed["requested_tasks"], new_tasks=len(new), reused_R8_endpoints=len(rows)-len(new),
        metric_source=str(metric_path), normalization="RMSE/declared contrast amplitude; prediction-only256binOtsu",
        material_success="amplitudeNRMSE<=.15 and OtsuIoU>=.7 and legal [.02,20] box; exactly inherited R8 definitions",
        truth_paths_opened=[truths[c]["truth_path"] for c in CASES], new_physical_calls=0))
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with (out / "SCORES.csv").open("x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for r in rows:
            w.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v for k, v in r.items()})
    dump(out / "BUDGET_SCORES.json", dict(rows=budgets,
        selection="Saved original-Phi incumbents only; no interpolation or GT iterate selection",
        CPU_scope="Landscape missing startup prefix uses an explicit conservative upper shift, not an exact prefix",
        new_physical_calls=0))
    dump(out / "STAGES.json", dict(rows=stages, selection="Preserved original-Phi selection; no offline stage reselection"))
    dump(out / "INHERITED_AVAILABILITY.json", sealed["inherited_availability"])
    dump(out / "COST_LEDGER.json", ledger)
    dump(out / "RECEIPT.json", dict(status="SCORED_FROZEN_FOUR_CASE_DEVELOPMENT", cases=list(CASES),
        scoring_CPU_seconds=time.process_time()-cpu, scoring_wall_seconds=time.perf_counter()-begin,
        truth_access="FOUR_SELECTED_EXISTING_TRUTH_ARRAYS_AFTER_COMPLETE_PUBLIC_SEAL", new_physical_calls=0,
        new_LU=0, new_complex_RHS=0, public_seal=stamp(resolve(args.seal))))
    print(json.dumps(dict(status="SCORED", new_tasks=len(new), rows=len(rows), output=str(out))))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=("seal", "score"))
    ap.add_argument("--manifest", action="append", help="Repeat for all frozen manifests; default16+8+12+7 tasks")
    ap.add_argument("--r8-repo", default=str(REPO.parent / "r8_worktree"))
    ap.add_argument("--seal", default=str(ROOT / "audit/DEVELOPMENT_PUBLIC_SEAL.json"))
    ap.add_argument("--output", default=str(ROOT / "RESULTS/development_scoring"))
    args = ap.parse_args()
    args.manifest = args.manifest or [str(ROOT / "configs" / n) for n in DEFAULT_MANIFESTS]
    (seal if args.mode == "seal" else score)(args)


if __name__ == "__main__":
    main()
