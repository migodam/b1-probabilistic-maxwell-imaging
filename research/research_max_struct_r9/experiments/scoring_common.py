"""Offline extension helpers. No solver imports and no truth reads on import."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import json
import time
import numpy as np
from research.research_max_struct_r9.experiments import score_development as dev

ROOT, REPO = dev.ROOT, dev.REPO
load, dump, resolve, stamp = dev.load, dev.dump, dev.resolve, dev.stamp


def tasks_from(manifests, cases):
    result = []
    for path in manifests:
        doc = load(path)
        for task in doc["tasks"]:
            if task["case_id"] not in cases:
                raise ValueError("Task outside declared extension cases")
            result.append(dict(task, _kind="landscape" if "algorithm" in task else "baseline",
                               _manifest=str(path), _input_base=str(REPO)))
    if len({t["task_id"] for t in result}) != len(result):
        raise ValueError("Duplicate frozen task")
    return result


def make_seal(path, manifests, cases, r8, status, extra_evidence=(), portfolio_rule=None):
    if path.exists():
        raise FileExistsError(path)
    tasks = tasks_from(manifests, cases)
    entries, evidence, public_by_case = [], set(manifests) | set(extra_evidence), {}
    evidence.update((Path(__file__), Path(dev.__file__),
                     r8 / "research/research_max_struct_r7/experiments/score_primary.py"))
    for task in tasks:
        output, rp, rec = dev.task_receipt(task)
        entries.append(dev.endpoint_entry(task, output, rp, rec, "NEW_R9_"+task["_kind"].upper()))
        public = dev.public_arrays(task["public_files"], REPO)
        ident = task["case_id"]
        if ident in public_by_case and not dev.same_public(public, public_by_case[ident]):
            raise ValueError("Mismatched public acquisition for one case")
        public_by_case[ident] = public
        evidence.update(dev.public_path(p, REPO) for p in task["public_files"])
        if task.get("config"):
            evidence.add(resolve(task["config"]))
        if task.get("initialization_file"):
            evidence.add(dev.public_path(task["initialization_file"], REPO))
    portfolios = []
    if portfolio_rule is not None:
        evidence.add(portfolio_rule)
        by_id = {e["task"]["task_id"]: e for e in entries}
        for ident in cases:
            for rule in load(portfolio_rule)["portfolios"]:
                members = [ident+"__"+suffix for suffix in rule["members"]]
                candidates = []
                for task_id in members:
                    entry = by_id[task_id]
                    rec = load(entry["receipt_path"])
                    phi = rec.get("Phi") if entry["task"]["_kind"] == "landscape" else rec.get("physical_objective")
                    if entry["status"] in dev.SUCCESS:
                        if phi is None or not np.isfinite(phi):
                            raise ValueError("Nonfinite portfolio objective")
                        candidates.append((phi, task_id))
                portfolios.append(dict(case_id=ident, name=rule["name"], members=members,
                    selected_task_id=min(candidates)[1] if candidates else None,
                    selection="Minimum final original Phi, then task_id; before any truth access",
                    cost="Always pay both full candidates, including rejected/failed work", truth_access="NONE"))
        pp = path.with_name(path.stem+"_PORTFOLIOS.json")
        dump(pp, dict(status="PUBLIC_PHI_SELECTION_FROZEN_BEFORE_GT", records=portfolios))
        evidence.add(pp)
    sealed = dict(status=status, created_utc=datetime.now(timezone.utc).isoformat(),
        cases=list(cases), manifests=[str(p) for p in manifests], requested_tasks=len(tasks),
        r8_repo=str(r8), entries=entries, portfolios=portfolios,
        successful_predictions=sum(e["status"] in dev.SUCCESS for e in entries),
        failed_tasks=sum(e["status"] not in dev.SUCCESS for e in entries),
        evidence=[stamp(p) for p in sorted(evidence)], truth_access="NONE", hash_checks="NONE")
    dump(path, sealed)
    return sealed


def verify(path, manifests, cases, expected_tasks, status):
    sealed = load(path)
    if sealed["status"] != status or sealed["cases"] != list(cases):
        raise ValueError("Wrong extension public seal")
    if sealed["manifests"] != [str(p) for p in manifests]:
        raise ValueError("Changed requested manifest list")
    tasks = tasks_from(manifests, cases)
    if len(tasks) != expected_tasks or len(sealed["entries"]) != expected_tasks:
        raise ValueError("Complete frozen task count required")
    if {t["task_id"] for t in tasks} != {e["task"]["task_id"] for e in sealed["entries"]}:
        raise ValueError("Frozen task identity mismatch")
    for s in sealed["evidence"] + [s for e in sealed["entries"] for s in e["files"]]:
        if stamp(s["path"]) != s:
            raise RuntimeError("Changed sealed public artifact: "+s["path"])
    return sealed


def score_entries(sealed, truths, cohort):
    metrics, metric_path = dev.metric_module(Path(sealed["r8_repo"]))
    rows, stages, budgets, materials = [], [], [], {}
    for entry in sealed["entries"]:
        task, output = entry["task"], Path(entry["output"])
        rec = load(entry["receipt_path"]); ident = task["case_id"]; truth = truths[ident]
        landscape = task["_kind"] == "landscape"
        pcost = (rec.get("cost") or {}).get("physical", {}) if landscape else (rec.get("physical_model_cost") or {})
        acost = (rec.get("cost") or {}).get("auxiliary", {}) if landscape else (rec.get("auxiliary_execution_receipt") or {}).get("operation_counts", {})
        phase = rec.get("stages", []) if landscape else rec.get("auxiliary_candidate_scores", [])
        corr = next((p for p in reversed(phase) if p.get("label")=="original_corrector"), {}) if landscape else (rec.get("corrector") or {})
        row = dict(task_id=task["task_id"], case_id=ident, cohort=cohort,
            method=task.get("algorithm", task.get("arm")), dpcs_rule=task.get("dpcs_rule"),
            dpcs_decision=rec.get("dpcs_rule"), origin=entry["origin"], status=entry["status"],
            initialization=("supplied_public_material" if task.get("initialization_file") else task.get("initialization", "uniform")),
            initialization_file=task.get("initialization_file"), initialization_provenance=task.get("initialization_provenance"),
            amplitude=truth["amplitude"], material_nrmse=None, otsu_iou=None, material_image_success=False,
            original_Phi=rec.get("Phi") if landscape else rec.get("physical_objective"),
            cpu_seconds=rec.get("process_CPU_seconds") if landscape else rec.get("r9_total_cpu_seconds", rec.get("cpu_seconds")),
            wall_seconds=rec.get("process_wall_seconds") if landscape else rec.get("r9_total_wall_seconds", rec.get("wall_seconds")),
            nested_solver_cpu_seconds=(rec.get("cost") or {}).get("CPU_seconds") if landscape else rec.get("cpu_seconds"),
            pre_solver_CPU_seconds=rec.get("pre_solver_CPU_seconds"),
            adapter_preparation_cpu_seconds=rec.get("adapter_preparation_cpu_seconds"),
            LU=pcost.get("LU"), complex_RHS=pcost.get("actual_complex_RHS"),
            physical_model_cost=pcost, auxiliary_cost=acost,
            current_Cholesky=acost.get("current_Cholesky", 0), current_RHS=acost.get("current_RHS", 0),
            corrector_stop=corr.get("stop"), corrector_iterations=corr.get("iterations"), corrector_success=corr.get("success"),
            auxiliary_iterations=rec.get("auxiliary_iterations"),
            selected_stage=rec.get("selected_label") if landscape else rec.get("auxiliary_selected_stage"),
            selection_rule=rec.get("online_selection_rule", rec.get("selection")),
            receipt_path=entry["receipt_path"], prediction_path=str(output/"prediction.npz") if entry["status"] in dev.SUCCESS else None,
            source_manifest=task["_manifest"], configuration=rec.get("config") if landscape else rec.get("configuration"),
            execution_platform=rec.get("execution_platform"), error=rec.get("error"))
        if entry["status"] in dev.SUCCESS:
            with np.load(output/"prediction.npz", allow_pickle=False) as z:
                m, yp = z["m"].ravel(), z["Y_prediction"]
            public = dev.public_arrays(task["public_files"], REPO)
            row.update(metrics.material_metrics(m, truth["m"], truth["mask"], truth["amplitude"]))
            row.update(metrics.noise_metrics(yp, public, truth.get("same_grid") is True))
            if truth.get("same_grid") is None:
                row.update(noise_gate_scope="RAW NOISE diagnostic; historical grid equivalence not asserted")
            row["same_grid_material_noise_success"] = bool(row["material_image_success"] and row["noise_CN99_pass"]) if truth.get("same_grid") is True else None
            row["independent_saved_Phi"] = .5*row["relative_residual"]**2+metrics.fixed_prior(m)
            row["Phi_receipt_difference"] = row["independent_saved_Phi"]-row["original_Phi"]
            materials[task["task_id"]] = m.copy()
            budgets.extend(dev.score_budgets(entry, rec, truth, metrics))
        rows.append(row)
        stages.extend(dict(task_id=task["task_id"], case_id=ident, index=i, recorded_stage=p) for i,p in enumerate(phase))
    return rows, stages, budgets, materials, metric_path


def add_portfolios(sealed, rows, materials):
    by_id = {r["task_id"]: r for r in rows}; result=[]
    for choice in sealed.get("portfolios", []):
        candidates = [by_id[i] for i in choice["members"]]
        chosen = by_id.get(choice["selected_task_id"])
        task_id = choice["case_id"]+"__portfolio__"+choice["name"]
        row = dict(chosen) if chosen is not None else dict(case_id=choice["case_id"], material_image_success=False,
            material_nrmse=None, otsu_iou=None, original_Phi=None)
        row.update(task_id=task_id, method=choice["name"], origin="FROZEN_TWO_BRANCH_PORTFOLIO",
            selected_task_id=choice["selected_task_id"], candidate_task_ids=choice["members"],
            selection_rule=choice["selection"], dpcs_rule="two_frozen_candidates",
            status="SELECTED_PUBLIC_PREDICTION" if chosen is not None else "ALL_CANDIDATES_FAILED")
        for key in ("cpu_seconds", "wall_seconds", "LU", "complex_RHS", "current_Cholesky", "current_RHS"):
            values = [r.get(key) for r in candidates]
            row[key] = sum(values) if all(v is not None for v in values) else None
        row["cost_scope"] = "Sum both full candidates including rejected branch; CPU and task-wall sums, not parallel elapsed time. Nested selected-state costs are not added again."
        for scope in ("physical_model_cost", "auxiliary_cost"):
            row[scope] = {k:sum((r.get(scope) or {}).get(k,0) for r in candidates)
                          for k in set().union(*(set(r.get(scope) or {}) for r in candidates))
                          if all(isinstance((r.get(scope) or {}).get(k,0),(int,float)) for r in candidates)}
        row["candidate_costs"] = [{k:r.get(k) for k in ("task_id","cpu_seconds","wall_seconds","LU","complex_RHS","current_Cholesky","current_RHS")} for r in candidates]
        if chosen is not None:
            materials[task_id] = materials[chosen["task_id"]].copy()
        result.append(row)
    return result


def write_bundle(out, rows, stages, budgets, materials, truths, sealed, metric_path, status, elapsed):
    out.mkdir(parents=True, exist_ok=False)
    dump(out/"SCORES.json", dict(status=status, rows=rows, cases=sealed["cases"],
        metric_source=str(metric_path), truth_paths_opened=[truths[c]["truth_path"] for c in sealed["cases"]],
        material_success="RMSE/amplitude<=.15, prediction-only256binOtsuIoU>=.7, legal[.02,20]",
        new_physical_calls=0))
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (out/"SCORES.csv").open("x", newline="") as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader()
        writer.writerows({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)
    dump(out/"STAGES.json",dict(rows=stages,selection="Recorded paid metadata only; no GT-based stage selection"))
    dump(out/"BUDGET_SCORES.json",dict(rows=budgets,selection="Original-Phi saved incumbent only; exact prefix or explicitly labeled upper shift; no interpolation"))
    regular=[r for r in rows if r["origin"]!="FROZEN_TWO_BRANCH_PORTFOLIO"]
    dump(out/"COST_LEDGER.json",dict(tasks=len(regular),CPU_known_seconds=sum(r["cpu_seconds"] for r in regular if r["cpu_seconds"] is not None),
        CPU_unmeasured_tasks=[r["task_id"] for r in regular if r["cpu_seconds"] is None],
        task_wall_seconds_sum=sum(r["wall_seconds"] for r in regular if r["wall_seconds"] is not None),
        LU_known=sum(r["LU"] for r in regular if r["LU"] is not None),complex_RHS_known=sum(r["complex_RHS"] for r in regular if r["complex_RHS"] is not None),
        current_Cholesky=sum(r["current_Cholesky"] for r in regular), current_RHS=sum(r["current_RHS"] for r in regular),
        scope="Inclusive task totals; nested work overlaps. Portfolio cost reuses already charged tasks, never add portfolio totals to experiment total. Task wall sum is not queue elapsed time."))
    tids=list(materials);ids=sealed["cases"]
    np.savez_compressed(out/"MATERIALS.npz",prediction=np.array([materials[i] for i in tids]),
        truth=np.array([truths[i]["m"] for i in ids]))
    dump(out/"MATERIAL_INDEX.json",dict(task_ids=tids,case_ids=ids,truth_scope="Offline output only after full public seal"))
    dump(out/"RECEIPT.json",dict(status=status,cases=ids,requested_tasks=sealed["requested_tasks"],
        scoring_CPU_seconds=time.process_time()-elapsed[0],scoring_wall_seconds=time.perf_counter()-elapsed[1],
        new_physical_calls=0,new_LU=0,new_complex_RHS=0,
        truth_access="SELECTED_EXTENSION_TRUTH_ONLY_AFTER_COMPLETE_PUBLIC_SEAL"))
