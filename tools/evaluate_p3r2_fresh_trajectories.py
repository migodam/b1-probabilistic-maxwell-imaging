"""Post-closure locked trajectory-risk evaluation, scheduled under Book L.

All online estimates and required closures are checked before the first truth
archive is opened. No physics, online changes, or coefficient fitting occurs.
"""
import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"
import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
import numpy as np
from p3r2.trajectory_risk import METHODS
from infer_p3r2_locked import DEFAULT_OUT, resolve_path, save, write
from integrate_p3r2_fresh_trajectories import ALL_ARMS, CORE_ARM, OTHER_ARMS, require_complete, trajectory_contract


def paired_bootstrap(candidate, baseline, indices):
    candidate, baseline = np.asarray(candidate, float), np.asarray(baseline, float)
    difference = candidate-baseline
    replicates = difference[indices].mean(axis=1)
    relative = -replicates/np.maximum(baseline[indices].mean(axis=1), 1e-300)
    result = dict(parents=len(candidate), candidate_mean_MSE=float(candidate.mean()),
        baseline_mean_MSE=float(baseline.mean()), mean_difference=float(difference.mean()),
        relative_risk_reduction=float(1-candidate.mean()/baseline.mean()),
        improved_parent_count=int(np.sum(difference < 0)),
        interpretation="negative difference and positive relative reduction favor trajectory estimator")
    for confidence, name in ((.95, "95"), (.975, "97_5")):
        tail = (1-confidence)/2
        result["difference_CI_"+name] = np.quantile(replicates, [tail, 1-tail]).tolist()
        result["relative_reduction_CI_"+name] = np.quantile(relative, [tail, 1-tail]).tolist()
    result["descriptive_5percent_and_adjusted_CI_criterion"] = bool(
        result["relative_risk_reduction"] >= .05 and result["difference_CI_97_5"][1] < 0)
    return result


def validate_marker(marker, manifest, parents):
    if marker.get("parents") != parents or marker.get("count") != 256:
        raise ValueError("COMPLETE marker must cover all exact256 parents")
    if marker.get("source_revision") != manifest["source_revision"] or float(marker.get("frozen_unix", -1)) != float(manifest["frozen_unix"]):
        raise ValueError("COMPLETE marker belongs to another method freeze")


def counts_from_status(statuses, parents, arm):
    paths = [statuses[p]["arms"][arm] for p in parents]
    rows = [r for path in paths for r in path["selected_endpoint_diagnostics"]]
    paid = [r for row in rows for r in row.get("measured_endpoints", [])]
    paid_with_origins = [(endpoint, row) for row in rows for endpoint in row.get("measured_endpoints", [])]
    for endpoint in paid:
        if not np.isfinite(float(endpoint["Phi"])) or not np.isfinite(float(endpoint["actual_RMS"])):
            raise ValueError("nonfinite paid endpoint diagnostic")
    threshold = float(paths[0]["noise_ball_threshold"])
    return dict(parents=len(parents), expected_selected_endpoints=8*len(parents),
        complete_paths=sum(p["path_complete"] for p in paths),
        paths_with_failures=sum(bool(p["failures"]) for p in paths),
        selected_endpoints_observed=len(rows), selected_noise_ball_count=sum(r["noise_ball"] for r in rows),
        selected_delta_Phi_le_2_count=sum(r["delta_Phi_le_2"] for r in rows),
        selected_prior_support_legal_count=sum(r["native_prior_support_legal"] for r in rows),
        selected_radius_retained_count=sum(r["radius_retained"] for r in rows),
        selected_signed_progress_retained_count=sum(r["progress_retained"] for r in rows),
        selected_noise_ball_and_delta_Phi_le_2_count=sum(r["noise_ball"] and r["delta_Phi_le_2"] for r in rows),
        selected_step_material_RMS=[r["actual_step_material_RMS"] for r in rows],
        selected_anchor_displacement_RMS=[r["anchor_displacement_RMS"] for r in rows],
        selected_Gaussian_prior_penalties=[r["Gaussian_prior_penalty"] for r in rows],
        measured_predictor_corrector_endpoints=len(paid),
        measured_predictor_corrector_noise_ball_count=sum(float(r["Phi"]) <= threshold for r in paid),
        measured_predictor_corrector_delta_Phi_le_2_count=sum(
            float(endpoint["Phi"])-(row["Phi"]-row["delta_Phi_from_anchor"]) <= 2
            for endpoint, row in paid_with_origins),
        measured_predictor_corrector_radius_retained_count=sum(
            row["task_radius"]-1e-10 <= float(endpoint["actual_RMS"]) <= 1.02*row["task_radius"]
            for endpoint, row in paid_with_origins),
        measured_predictor_corrector_signed_progress_retained_count=sum(
            float(endpoint["signed_progress"]) >= .5*row["task_radius"]-1e-9
            for endpoint, row in paid_with_origins),
        measured_predictor_corrector_backward_residual_pass_count=sum(float(r["backward_residual"]) <= 1e-10 for r in paid),
        incremental_step_RHS=sum(row.get("cost", {}).get("linear_solve_rhs", 0) for row in rows),
        noise_ball_threshold=threshold,
        scope="counts of actual paid endpoint metadata; no likelihood assertion for material averages or connecting paths",
        prior_scope="finite native Gaussian64 x Uniform[0,1]^2 support; no event-mass or prior-mass qualification")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--methods", type=Path)
    args = parser.parse_args()
    out = args.root.resolve()
    manifest, parents, mapping, frozen, weights_path, source_closures = trajectory_contract(
        out, args.methods or out/"METHODS_FROZEN_FOR_LOCKED.json")
    point_folder, trajectory_folder = out/"L/estimates", out/"L/trajectory_estimates"
    point_complete = require_complete(point_folder/"COMPLETE.json")
    trajectory_complete = require_complete(trajectory_folder/"COMPLETE.json")
    validate_marker(point_complete, manifest, parents)
    validate_marker(trajectory_complete, manifest, parents)
    if point_complete.get("truth_read") is not False or trajectory_complete.get("truth_read") is not False:
        raise ValueError("both inference stages must explicitly certify truth_read=false")
    if trajectory_complete.get("trajectory_jobs") != manifest["trajectory_jobs"]:
        raise ValueError("trajectory integration differs from frozen registered jobs")
    if manifest.get("include_matched_point176") is not True:
        raise ValueError("matched extra16 ordinary baseline must be frozen for all256")
    generator = resolve_path(manifest.get("generator_closed_json", "L/jobs/data_locked_0_256/CLOSED.json"), out)
    require_complete(generator)
    generated = json.loads((generator.parent/"DATA_MANIFEST.json").read_text())
    if generated.get("group") != "locked" or set(generated.get("parents", [])) != set(parents):
        raise ValueError("completed locked generator must cover exact frozen256")
    context_jobs = manifest.get("context_closed_jsons")
    if not context_jobs:
        raise ValueError("frozen context_closed_jsons required")
    for value in list(context_jobs)+list(manifest.get("additional_online_closures", [])):
        require_complete(resolve_path(value, out))
    destination = out/"L/trajectory_evaluation"
    if (destination/"COMPLETE.json").exists() or (destination/"RESULTS.json").exists():
        raise FileExistsError("locked trajectory evaluation exists; no implicit repeat")
    start_wall, start_cpu = time.monotonic(), time.process_time()
    # Materialize every scored online array and status BEFORE any truth read.
    # Matched control is scored only from the saved frozen L/estimates array;
    # raw MATCHED.npz is never opened by this evaluator.
    online, statuses, point_statuses = {}, {}, {}
    for parent in parents:
        sidecar = json.loads((out/"L/contexts"/f"{parent}.json").read_text())
        if sidecar.get("native_cap") != 96 or sidecar.get("probit_cap") != 64 or sidecar.get("GT_read") is not False:
            raise ValueError("all256 truth-free native96 + probit64 contexts required")
        receipt = json.loads((out/"L/contexts/details"/parent/"matched/MATCHED.json").read_text())
        if receipt.get("extra_fresh_state_cap") != 16 or receipt.get("GT_used") is not False:
            raise ValueError("all256 truth-free matched extra16 receipts required")
        registered = [arm for arm in ALL_ARMS if (parent, arm) in mapping]
        estimates = {}
        for folder, names in ((point_folder, ["probit", "matched_point176"]),
                (trajectory_folder, [f"{arm}__{method}" for arm in registered for method in METHODS])):
            with np.load(folder/f"{parent}.npz", allow_pickle=False) as archive:
                if str(archive["source_revision"]) != manifest["source_revision"] or float(archive["frozen_unix"]) != float(manifest["frozen_unix"]):
                    raise ValueError("parent estimate source freeze mismatch")
                for name in names:
                    value = archive[name].copy()
                    if value.shape != (576,) or not np.isfinite(value).all():
                        raise ValueError("finite online estimate required; no parent exclusions")
                    estimates[name] = value
                if folder == trajectory_folder:
                    if list(archive["registered_arms"]) != registered or list(archive["method_names"]) != list(METHODS):
                        raise ValueError("trajectory parent arm/method schema mismatch")
                    statuses[parent] = json.loads(str(archive["status_json"]))
                else:
                    point_statuses[parent] = json.loads(str(archive["status_json"]))
        if statuses[parent].get("parent") != parent or set(statuses[parent]["arms"]) != set(registered):
            raise ValueError("trajectory online status coverage mismatch")
        online[parent] = estimates
    diagnostics = {arm:counts_from_status(statuses, parents if arm == CORE_ARM else parents[:32], arm) for arm in ALL_ARMS}
    truth_paths = [out/"inputs/offline_only"/f"{parent}.npz" for parent in parents]
    if any(not path.is_file() or path.stat().st_mtime < float(manifest["frozen_unix"]) for path in truth_paths):
        raise ValueError("all256 locked truth files must exist and postdate method freeze")
    destination.mkdir(parents=True, exist_ok=True)
    # This is the first labelled access. All gates and online loads passed.
    risks, truth_opened = {}, []
    for parent, truth_path in zip(parents, truth_paths):
        with np.load(truth_path, allow_pickle=False) as archive:
            truth = archive["chi"].copy()
        truth_opened.append(parent)
        if truth.shape != (576,) or not np.isfinite(truth).all():
            raise ValueError("invalid truth; no parent exclusion allowed")
        risks[parent] = {name:dict(fullgrid_MSE=float(np.mean(abs(estimate-truth)**2)),
            fullgrid_MAE=float(np.mean(abs(estimate-truth)))) for name, estimate in online[parent].items()}
    rng = np.random.default_rng(2026102292)
    core_indices = rng.integers(0, 256, size=(2000, 256))
    secondary_indices = np.random.default_rng(2026102292).integers(0, 32, size=(2000, 32))
    summaries, comparisons, raw = {}, {}, {}
    for arm in ALL_ARMS:
        cohort = parents if arm == CORE_ARM else parents[:32]
        indices = core_indices if arm == CORE_ARM else secondary_indices
        summaries[arm], comparisons[arm] = {}, {}
        for method in METHODS:
            name = f"{arm}__{method}"
            mse = np.array([risks[p][name]["fullgrid_MSE"] for p in cohort])
            mae = np.array([risks[p][name]["fullgrid_MAE"] for p in cohort])
            raw[name+"__MSE"], raw[name+"__MAE"] = mse, mae
            fallback = [statuses[p]["arms"][arm]["methods"][method]["fallback_to_anchor"] for p in cohort]
            summaries[arm][method] = dict(parents=len(cohort), mean_MSE=float(mse.mean()), mean_MAE=float(mae.mean()),
                fallback_count=sum(fallback), fallback_parents=[p for p, used in zip(cohort, fallback) if used])
            comparisons[arm][method] = {}
            for baseline in ("probit", "matched_point176"):
                base = np.array([risks[p][baseline]["fullgrid_MSE"] for p in cohort])
                comparisons[arm][method]["vs_"+baseline] = paired_bootstrap(mse, base, indices)
    # Core curvature also has the same first32 subset as other-arm secondary
    # comparisons, so arm comparisons never mix 256 and32-parent panels.
    first32 = {}
    for arm in ALL_ARMS:
        first32[arm] = {}
        for method in METHODS:
            name = f"{arm}__{method}"
            mse = np.array([risks[p][name]["fullgrid_MSE"] for p in parents[:32]])
            baseline = np.array([risks[p]["probit"]["fullgrid_MSE"] for p in parents[:32]])
            first32[arm][method] = paired_bootstrap(mse, baseline, secondary_indices)
    primary = comparisons[CORE_ARM]["arc_length_trapezoid"]
    save(destination/"RAW_RISKS.npz", parent=np.asarray(parents), secondary_parent=np.asarray(parents[:32]),
        **raw, probit_MSE=np.array([risks[p]["probit"]["fullgrid_MSE"] for p in parents]),
        matched_point176_MSE=np.array([risks[p]["matched_point176"]["fullgrid_MSE"] for p in parents]))
    result = dict(status="COMPLETE", count=256, parents=parents, no_exclusions=True,
        primary_metric="parent mean full-grid complex material MSE",
        primary_arm=CORE_ARM, primary_method="arc_length_trapezoid",
        co_primary_vs_strong_probit=primary["vs_probit"],
        matched_extra16_ordinary_control=primary["vs_matched_point176"],
        co_primary_scope="second prespecified risk estimator alongside separately evaluated fused action; two-sided97.5% CI reported",
        strong_baseline="saved probit point: native96 plus incremental probit64",
        matched_baseline="saved matched_point176: same strong context plus ordinary extra16 fresh-state cap",
        bootstrap=dict(draws=2000, seed=2026102292, unit="paired independent parent", confidence_levels=[.95, .975]),
        summaries=summaries, comparisons=comparisons, same_first32_secondary=first32,
        actual_candidate_diagnostics=diagnostics, per_parent_risks=risks, point_inference_status=point_statuses,
        trajectory_fallbacks_retained=True, source_job_closures=source_closures,
        source_revision=manifest["source_revision"], frozen_unix=manifest["frozen_unix"],
        frozen_weights_path=str(weights_path), truth_opened_only_after_all_gates=True, physics_calls=0,
        coefficient_semantics="convex material estimator coefficients, never posterior probabilities",
        limits=["finite sampled trajectories, not verified connecting paths", "no posterior integration, mixing, or mass claim",
            "convex material average need not equal a single legal decoder state", "estimator averages have not received additional full-forward checks",
            "noise-ball and DeltaPhi counts concern measured endpoints only", "secondary first32 comparisons do not redefine primary"],
        scientific_status="fresh locked material-risk evidence only; final scientific judgment remains with parent")
    write(destination/"RESULTS.json", result)
    lines = ["# Fresh locked trajectory risk", "", "All256 parents retained. Lower full-grid complex MSE is better.", "",
        "| comparison | parents | relative reduction | 95% CI | 97.5% CI |", "|---|---:|---:|---|---|"]
    for name, comparison in (("curvature arc vs strong probit (co-primary)", primary["vs_probit"]),
            ("curvature arc vs matched ordinary extra16", primary["vs_matched_point176"])):
        lines.append(f"| {name} | 256 | {comparison['relative_risk_reduction']:.6g} | {comparison['relative_reduction_CI_95']} | {comparison['relative_reduction_CI_97_5']} |")
    lines.extend(["", "Convex coefficients are estimator weights, not posterior probabilities. Candidate endpoint checks do not certify the averaged material or a continuous path."])
    (destination/"TABLE.md").write_text("\n".join(lines)+"\n")
    write(destination/"COMPLETE.json", dict(status="COMPLETE", count=256, parents=parents,
        truth_archives_opened=truth_opened, source_revision=manifest["source_revision"],
        frozen_unix=manifest["frozen_unix"], physics_calls=0,
        cpu_seconds=time.process_time()-start_cpu, wall_seconds=time.monotonic()-start_wall))


if __name__ == "__main__":
    main()
