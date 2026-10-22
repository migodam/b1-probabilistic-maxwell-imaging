"""Consolidate saved evidence only; never read truth, checkpoints or JSONL.

This is an evidence index and plotting pass, not a scientific verdict. Existing
reference statuses and FULL/SLICE/PROFILE meanings are preserved verbatim.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import math

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "research/p3l_local_posterior_imaging"


def read(path):
    path = Path(path)
    if path.suffix != ".json" or path.stat().st_size > 24*1024*1024:
        raise ValueError("only bounded JSON summaries may be read")
    return json.loads(path.read_text())


def compact(value):
    """Discard coordinate arrays/traces while retaining scalar evidence rows."""
    omit = {"named_diagnostics", "history", "logs", "gradient", "objective_gradient",
            "projected_gradient", "q", "x", "initialization_x", "initialization_q",
            "source_paths", "chain_sources", "sample_distance_from_mean_rms"}
    if isinstance(value, dict):
        return {k: compact(v) for k, v in value.items() if k not in omit}
    if isinstance(value, list):
        return [compact(v) for v in value] if len(value) <= 256 else dict(omitted_large_list=True, count=len(value))
    return value


def block(relative):
    path = OUT / relative
    if not path.exists():
        return dict(status="NOT_RUN", source=relative)
    try:
        return dict(status="SAVED_EVIDENCE", source=relative, evidence=compact(read(path)))
    except Exception as error:
        return dict(status="READ_FAILED_RETAINED", source=relative, error=repr(error))


def parent_closure():
    path = OUT/"CAMPAIGN_CLOSURE.json"
    if not path.exists():
        return None
    closure = read(path)
    for field in ("scientific_verdict", "closed_at_utc", "probability_reference_status"):
        if not isinstance(closure.get(field), str) or not closure[field].strip():
            raise ValueError(f"parent closure requires nonempty string: {field}")
    for field in ("headline", "evidence_release_url"):
        if field in closure and not isinstance(closure[field], str):
            raise ValueError(f"parent closure {field} must be a string")
    findings = closure.get("findings", [])
    if not isinstance(findings, str) and not (isinstance(findings, list) and all(isinstance(v, str) for v in findings)):
        raise ValueError("parent closure findings must be a string or list of strings")
    return {key: closure[key] for key in ("scientific_verdict", "closed_at_utc", "probability_reference_status",
            "headline", "findings", "evidence_release_url") if key in closure}


def job_cost(stage, job):
    path = OUT / stage / "jobs" / job / "CLOSED.json"
    if not path.exists():
        return dict(status="NOT_RUN", source=str(path.relative_to(OUT)), RHS=None)
    receipt = read(path)
    counters = receipt.get("counters", {})
    # Solver RHS is a per-job saved counter. Never sum cumulative stage values.
    rhs = counters.get("linear_solve_rhs")
    return dict(status=receipt.get("status"), source=str(path.relative_to(OUT)), RHS=rhs,
                forward_RHS=counters.get("full_forward_RHS"),
                receiver_RHS=counters.get("receiver_transfer_adjoint_RHS"),
                controlled_wall_seconds=receipt.get("controlled_wall_seconds"),
                cost_scope="saved per-job linear_solve_rhs; stage cumulative excluded")


def references():
    rows = []
    for stage in ("A", "B", "C"):
        for path in sorted((OUT/stage/"references").glob("*/*/REFERENCE_STATS.json")):
            stats = read(path)
            costs = []
            for source in stats.get("chain_sources", []):
                parts = Path(source["path"]).parts
                if "jobs" in parts:
                    i = parts.index("jobs")
                    costs.append(job_cost(parts[i-1], parts[i+1]))
            known = bool(costs) and all(c["RHS"] is not None for c in costs)
            rows.append(dict(stage=stage, parent=path.parent.parent.name, mechanism=path.parent.name,
                source=str(path.relative_to(OUT)), status=stats.get("status"),
                target_family=stats.get("target_family"),
                source_reference_qualified=stats.get("reference_qualified", False),
                cross_mechanism_status=stats.get("cross_mechanism_status"),
                warmup_removed=stats.get("warmup_removed"), chain_count=stats.get("chain_count"),
                production_draws_per_chain=stats.get("production_draws_per_chain"),
                summaries=compact(stats.get("summaries", {})), gates=stats.get("gates", {}),
                full_grid_mean_rms_mcse=stats.get("full_grid_mean_rms_mcse"),
                cost=dict(RHS=sum(c["RHS"] for c in costs) if known else None, chains=costs),
                interpretation="Source diagnostics only; no automatic qualification from ESS or marginal PIT"))
    return dict(status="SAVED_EVIDENCE" if rows else "NOT_RUN", rows=rows)


def smc():
    rows = []
    for path in sorted((OUT/"A/jobs").glob("smc_*/RESULT.json")):
        source = read(path)
        history = source.get("history", [])
        rows.append(dict(source=str(path.relative_to(OUT)), result=compact(source),
                         final_distinct_prior_ancestors=history[-1].get("distinct_prior_ancestors") if history else None,
                         minimum_distinct_prior_ancestors=min((r["distinct_prior_ancestors"] for r in history if "distinct_prior_ancestors" in r), default=None),
                         cost=job_cost("A", path.parent.name)))
    return dict(status="SAVED_EVIDENCE" if rows else "NOT_RUN", rows=rows,
                semantics="finite-N adaptive SMC; high terminal weight ESS cannot undo ancestry collapse or certify posterior coverage")


def learning():
    rows = []
    for family in ("mlp_diag", "mlp_full", "mlp_mdn", "mlp_flow", "cnn_full"):
        for seed in (2026102270, 2026102271, 2026102272):
            source = f"D/models/{family}/seed{seed}/TRAINING_REPORT.json"
            path = OUT/source
            row = dict(family=family, seed=seed, source=source, status="NOT_RUN")
            if path.exists():
                report = read(path)
                metadata = report.get("metadata", {})
                evaluations = [v for v in report.get("logs", []) if v.get("event") == "COMPLETE_VALIDATION"]
                row.update(status="SAVED_EVIDENCE", metadata=compact(metadata),
                    best_validation_NLL=metadata.get("best_validation_raw_target_nll"),
                    best_update=metadata.get("best_update"),
                    initial_validation_NLL=next((v.get("validation_raw_target_nll") for v in evaluations if v.get("update") == 0), None),
                    final_completed_validation_NLL=evaluations[-1].get("validation_raw_target_nll") if evaluations else None,
                    best_is_initialization=metadata.get("best_update") == 0,
                    complete_validation_evaluations=len(evaluations))
            rows.append(row)
    return dict(neural=dict(expected_runs=15, completed_reports=sum(r["status"] == "SAVED_EVIDENCE" for r in rows),
                            selected_initialization_count=sum(r.get("best_is_initialization", False) for r in rows), rows=rows),
                simple=block("D/simple_density/RESULTS.json"),
                independent_training_review=block("AUDIT/neural_training/ACTUAL_TRAINING_REVIEW.json"),
                interpretation="Validation density fit only; initialization selection and marginal PIT do not qualify a full posterior")


def mean_risks(report):
    rows = report.get("rows", [])
    methods = sorted({name for row in rows for name in row.get("risks", {})})
    result = {}
    for name in methods:
        values = [r.get("risks", {}).get(name, {}).get("full_grid", {}).get("complex_MSE") for r in rows]
        result[name] = dict(parent_count=len(rows), missing_count=sum(v is None for v in values),
                            full_grid_complex_MSE=None if not rows or any(v is None for v in values) else sum(values)/len(values))
    return result


def risk_costs(blocks, anchor_panel, locked):
    points = []
    anchors = anchor_panel.get("evidence", {}).get("parents", [])
    receipts = [job_cost("A", f'anchor_{r["parent"]}_{r["method"]}') for r in anchors]
    base_rhs = sum(r["RHS"] for r in receipts) if len(receipts) == 32 and all(r["RHS"] is not None for r in receipts) else None
    joint_cost = job_cost("B", "joint_probit_gaussian")
    for label in ("fixed_gaussian_risk", "joint_gaussian_risk"):
        means = mean_risks(blocks[label].get("evidence", {}))
        for name, values in means.items():
            if label == "fixed_gaussian_risk":
                cost = base_rhs
            else:
                cost = base_rhs if name == "historical_native_q_anchor" else (None if base_rhs is None or joint_cost["RHS"] is None else base_rhs+joint_cost["RHS"])
                if name == "full_target_IS_unqualified_mean":
                    cost = None  # pilot subset requires its separately declared cost.
            points.append(dict(report=label, estimator=name, **values,
                total_RHS_for_32=None if cost is None else cost,
                mean_RHS_per_parent=None if cost is None else cost/32,
                cost_scope="shared historical-anchor + joint-refinement paid solver costs; zero extra physics for Gaussian integration; not a new equal-cost claim"))
    online_rows = blocks.get("locked_online", {}).get("evidence", {}).get("rows", [])
    for saved in locked:
        for name, values in mean_risks(saved.get("evidence", {})).items():
            phases = {"cheap16": ["cheap16"], "independent_native96": ["independent_native96"],
                      "probit64_coordinate_mode": ["independent_native96", "probit_refinement64"],
                      "model_times_strong_Gaussian_mean": ["cheap16"],
                      "conditional_mode_plugin_mean": ["cheap16"]}.get(name)
            costs = ([r.get("phases", {}).get(phase, {}).get("cost", {}).get("RHS")
                      for r in online_rows for phase in phases] if phases is not None else [])
            total = sum(costs) if len(online_rows) == 128 and costs and all(v is not None for v in costs) else None
            points.append(dict(report="locked_offline_risk", estimator=name, **values,
                 total_RHS_for_128=total, mean_RHS_per_parent=None if total is None else total/128,
                 cost_scope="saved per-parent actual charged phase RHS, including failed attempts; effective fallback risk preserved; not an equal-cost assertion"))
    return dict(points=points, anchor_cost_receipts=receipts, joint_refinement_cost=joint_cost,
                missing_cost_policy="unknown remains null; never infer RHS from ESS/acceptance or stage cumulative")


def charts(report):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    folder = OUT/"figures"
    folder.mkdir(exist_ok=True)
    outputs = []
    def finish(fig, name):
        path = folder/(name+".png")
        fig.savefig(path, dpi=160, bbox_inches="tight"); plt.close(fig)
        outputs.append(str(path.relative_to(OUT)))
    fig, ax = plt.subplots(figsize=(10, 6))
    plotted = 0
    for row in report["risk_vs_RHS"]["points"]:
        x, y = row["mean_RHS_per_parent"], row["full_grid_complex_MSE"]
        if x is not None and y is not None:
            ax.scatter(x, y); ax.annotate(row["report"]+"/"+row["estimator"], (x, y), fontsize=6)
            plotted += 1
    if not plotted:
        ax.text(.5, .5, "Complete risk/cost pair NOT_RUN or unavailable", ha="center", transform=ax.transAxes)
    ax.set(xlabel="Saved solver RHS / parent (DEV and locked labeled; not equal-cost comparison)", ylabel="Full-grid complex MSE",
           title="Saved image risks versus attributable paid RHS; Gaussian references remain approximate")
    finish(fig, "p3l_campaign_risk_vs_RHS")
    fig, ax = plt.subplots(figsize=(10, 6))
    for row in report["references"]["rows"]:
        x = row["cost"]["RHS"]
        y = row["summaries"].get("material_functions", {}).get("minimum_bulk_ess")
        if x is not None and y is not None and x > 0 and y > 0:
            ax.scatter(x, y); ax.annotate(row["stage"]+"/"+row["parent"]+"/"+row["mechanism"], (x, y), fontsize=7)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set(xlabel="Saved per-job solver RHS (all supplied chains)", ylabel="Minimum material-function bulk ESS",
           title="Reference computation versus diagnostics; high ESS alone does not qualify FULL posterior")
    finish(fig, "p3l_campaign_reference_ESS_vs_RHS")
    neural = report["learning"]["neural"]["rows"]
    labels, values, colors = [], [], []
    for row in neural:
        if row.get("best_validation_NLL") is not None:
            labels.append(row["family"]+"/"+str(row["seed"])[-2:]+(" (update0)" if row.get("best_is_initialization") else ""))
            values.append(row["best_validation_NLL"]); colors.append("#b57334")
    simple = report["learning"]["simple"].get("evidence", {})
    for name in ("selected", "unconditional_selected"):
        if name in simple:
            labels.append("simple/"+simple[name]["id"]); values.append(simple[name]["validation_raw_NLL"]); colors.append("#25667d")
    fig, ax = plt.subplots(figsize=(11, max(4, len(labels)*.33)))
    ax.barh(range(len(values)), values, color=colors); ax.set_yticks(range(len(values)), labels); ax.invert_yaxis()
    ax.set(xlabel="Validation raw-coordinate NLL (lower is better)", title="Actual simple and neural density fits; update0 means initialization selected")
    finish(fig, "p3l_campaign_simple_vs_NN_NLL")
    return outputs


def main():
    ledger_block = block("RESOURCE_LEDGER.json")
    ledger = ledger_block.get("evidence", {})
    closure = parent_closure()
    sources = dict(anchors="A/ANCHOR_PANEL.json", fixed_gaussian_risk="B/GAUSSIAN_IMAGE_RISK.json",
        joint_gaussian_risk="B/JOINT_GAUSSIAN_IMAGE_RISK.json", conditional="B/CONDITIONAL_PILOT_RESULTS.json",
        GN_departure="B/GN_DEPARTURE_AUDIT.json", slice3="B/slice3/RESULTS.json",
        decision_cost="C/decision_cost/SUMMARY.json", material_brackets="C/material_brackets/SUMMARY.json",
        geometry="C/jobs/geometry_projector_DEV000_DEV028/RESULTS.json",
        learned_correction="D/jobs/learned_correction_exposed_DEV_s2026102270/SUMMARY.json",
        simple_correction="D/jobs/simple_correction_exposed_DEV/SUMMARY.json",
        learned_correction_audit="AUDIT/learned_correction/RESULTS.json",
        simple_correction_audit="AUDIT/simple_correction/RESULTS.json",
        locked_online="D/LOCKED_ONLINE_IMAGING.json")
    blocks = {name: block(path) for name, path in sources.items()}
    # Read only possible aggregate offline risk reports, never locked arrays.
    locked_paths = ["D/LOCKED_IMAGE_RISK.json", "D/LOCKED_OFFLINE_RISK.json", "D/LOCKED_IMAGING_RISK.json",
                    "D/LOCKED_ESTIMATOR_RISK.json", "D/LOCKED_RISK_RESULTS.json"]
    locked = [block(path) for path in locked_paths if (OUT/path).exists()]
    report = dict(status="CLOSED" if closure else "RUNNING",
        generated_at_utc=datetime.now(timezone.utc).isoformat(),
        final_scientific_judgment=closure["scientific_verdict"] if closure else "PENDING_PARENT_REVIEW",
        parent_authored_closure=closure,
        closure_source="CAMPAIGN_CLOSURE.json" if closure else None,
        active_reservations_present=bool(ledger.get("active_wall")),
        semantics=dict(original_target="Original P1 Gaussian64 × Uniform2, scalar-TM/VIE; unchanged target law",
                       FULL="original full target; source qualification retained, never inferred from PIT or ESS",
                       SLICE="conditioned coordinates/support; no full marginal probability claim",
                       PROFILE="plug-in/optimized coordinates; no volume integration claim",
                       Gaussian="declared Gaussian-reference images, not FULL posterior or qualified teacher",
                       locked="independent estimator-risk confirmation only if separate saved offline report exists"),
        resources=ledger_block, evidence=blocks, references=references(), SMC=smc(), learning=learning(),
        locked_offline_risk=locked if locked else dict(status="NOT_RUN"),
        truth_arrays_read=False, checkpoints_read=False, ledger_journal_read=False,
        new_physics_RHS=0, markdown_modified=False)
    report["risk_vs_RHS"] = risk_costs(blocks, blocks["anchors"], locked)
    report["figures"] = charts(report)
    temporary = OUT/"RESULTS_AND_METRICS.json.tmp"
    temporary.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+"\n")
    temporary.replace(OUT/"RESULTS_AND_METRICS.json")
    print(json.dumps(dict(status=report["status"], references=len(report["references"]["rows"]), figures=report["figures"], new_physics_RHS=0)))


if __name__ == "__main__":
    main()
