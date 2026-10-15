"""Saved-full-field candidate deletion and separate instrument mismatch diagnosis.

No Maxwell solver, training, or physical callback is imported or called here.
Online acquisition receives retained candidates, fields, old data and noise.
The offline simulator alone reads truth, hidden gain and potential outcomes.
"""
from __future__ import annotations
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import chi2

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from src.b1_joint.material import from_geometry, geometry, event_ids
from active_experiment import FREQS, IDS, menu, candidate_scores
from e0_candidate_controls import profile
from sensing import field_l1
from cost import Cost, ROOT, dump

GROUPS = ("full", "central_deleted", "annular_deleted", "combined_deleted")
POLICIES = ("none", "fixed2.4", "fixed3.2", "raw", "robust", "sentinel_first3.2")
SCENARIOS = ("in_model", "instrument_gain_drift")
FAMILIES = ("solid", "central_only", "annular_only", "combined")


def family_ids(x):
    shape = geometry(np.atleast_2d(x))
    return ((shape[:, 0] > 0).astype(int)
            + 2*(shape[:, 1] > 0).astype(int))


def support_masks(x):
    shape = geometry(x)
    return {"full": np.ones(len(x), dtype=bool),
            "central_deleted": shape[:, 0] == 0,
            "annular_deleted": shape[:, 1] == 0,
            "combined_deleted": ~((shape[:, 0] > 0) & (shape[:, 1] > 0))}


def recommend_retained(x, fields, old_y, sigma, calibration_sd):
    """Candidate-only: no truth, future data, prior or posterior arguments."""
    started = time.perf_counter()
    actions = menu(fields, sigma, calibration_sd)
    f0 = fields[:, FREQS.index(1.6)]
    scores = {a.id: candidate_scores(x, f0, old_y,
                                    np.full(32, sigma[FREQS.index(1.6)]),
                                    a.fields, a.sigma) for a in actions}
    return {"picks": {kind: max(scores, key=lambda aid: scores[aid][kind])
                       for kind in ("raw", "robust")},
            "scores": scores,
            "seconds": time.perf_counter()-started,
            "candidate_count": len(x),
            "candidate_weights_kind": "none",
            "recommendation_recomputed_for_retained_set": True}


class OfflineSimulator:
    """Hidden data are used only to return an observation or score a decision."""
    def __init__(self, path):
        with np.load(path, allow_pickle=False) as data:
            self.truth = data["truth"].copy()
            self.fields = data["F_truth"].copy()
            self.gain = data["gain"].copy()
            self.outcomes = data["outcomes"].copy()
            self.calibration = data["calibration"].copy()

    def initial(self, parent, noise):
        return self.outcomes[parent, noise, 0, FREQS.index(1.6)].copy()

    def measure(self, parent, noise, aid, scenario):
        if aid == "reference_calibration":
            return np.atleast_1d(self.calibration[parent, noise, 0]).copy()
        j = IDS.index(aid)
        observed = self.outcomes[parent, noise, 1, j].copy()
        if scenario == "instrument_gain_drift":
            observed += (self.gain[parent, noise] * .03*(FREQS[j]-1.6)
                         * self.fields[parent, j])
        return observed

    def score_catalogue(self, parent, candidates):
        truth = self.truth[parent]
        distances = field_l1(candidates, truth).ravel()
        true_family = int(family_ids(truth)[0])
        events = event_ids(candidates)
        return {"losses": distances,
                "best_of_K_MAE": float(distances.min()),
                "Coverage_at_K_005": bool(distances.min() <= .05),
                "true_family_id": true_family,
                "true_family": FAMILIES[true_family],
                "true_event_id": int(event_ids(truth)),
                "true_family_represented": bool(np.any(family_ids(candidates) == true_family)),
                "true_event_represented": bool(np.any(events == event_ids(truth)))}


def fitted_selection(candidates, fields, y, sd):
    """One-world shared-gain GLRT; neither residual softmax nor a posterior."""
    residual, gain = profile(fields, y, sd)
    chosen = int(np.argmin(residual))
    compatible = residual <= residual[chosen] + 9.
    dof = max(len(y)-1, 1)
    threshold = float(chi2.ppf(.99, dof))
    return {"chosen": chosen, "gain": float(gain[chosen]),
            "chi2": float(residual[chosen]), "dof": dof,
            "threshold": threshold, "sentinel": bool(residual[chosen] > threshold),
            "compatible_candidates": int(compatible.sum()),
            "compatible_families": sorted(np.unique(family_ids(candidates[compatible])).tolist())}


def parent_statistics(rows, seed, bootstrap):
    numeric = ("selected_MAE", "initial_selected_MAE", "best_of_K_MAE",
               "Coverage_at_K_005", "selected_structure_correct",
               "selected_family_correct", "true_family_represented",
               "residual_sentinel", "undetected_missing_family",
               "apparent_wrong_single_family", "decision_accepted", "measurement_blocks")
    by_key = {}
    for row in rows:
        key = (row["scenario"], row["support_group"], row["policy"], row["parent"])
        by_key.setdefault(key, []).append(row)
    parent_rows, lookup = [], {}
    for (scenario, group, policy, parent), children in sorted(by_key.items()):
        if sorted(r["noise"] for r in children) != [0, 1]:
            raise ValueError("Each independent parent requires exactly two noise episodes")
        row = {"scenario": scenario, "support_group": group,
               "policy": policy, "parent": parent,
               **{k: float(np.mean([v[k] for v in children])) for k in numeric}}
        parent_rows.append(row)
        lookup[(scenario, group, policy, parent)] = row
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, 32, size=(bootstrap, 32))
    aggregates, paired_none, paired_full = [], [], []
    for scenario in SCENARIOS:
        for group in GROUPS:
            for policy in POLICIES:
                current = [lookup[(scenario, group, policy, p)] for p in range(32)]
                aggregate = {"scenario": scenario, "support_group": group, "policy": policy,
                             "parents": 32,
                             **{k: float(np.mean([r[k] for r in current])) for k in numeric}}
                child_rows = [r for r in rows if r["scenario"] == scenario
                              and r["support_group"] == group and r["policy"] == policy]
                accepted = [r["selected_MAE"] for r in child_rows if r["decision_accepted"]]
                aggregate["accepted_proposal_MAE"] = float(np.mean(accepted)) if accepted else None
                aggregate["accepted_proposal_MAE_scope"] = (
                    "selective diagnostic only; always compare acceptance rate and unconditional selected_MAE")
                aggregates.append(aggregate)
                def paired(reference_group, reference_policy):
                    delta = np.array([lookup[(scenario, group, policy, p)]["selected_MAE"]
                                      - lookup[(scenario, reference_group, reference_policy, p)]["selected_MAE"]
                                      for p in range(32)])
                    return {"scenario": scenario, "support_group": group, "policy": policy,
                            "reference_group": reference_group, "reference_policy": reference_policy,
                            "paired_MAE_difference": float(delta.mean()),
                            "parent_bootstrap_CI95": np.quantile(delta[indices].mean(1), [.025, .975]).tolist(),
                            "independent_parents": 32, "bootstrap": bootstrap,
                            "negative_means_current_lower_MAE": True}
                if policy != "none":
                    paired_none.append(paired(group, "none"))
                if group != "full":
                    paired_full.append(paired("full", policy))
    return parent_rows, aggregates, paired_none, paired_full


def run(cost):
    started, cpu = time.perf_counter(), time.process_time()
    destination = ROOT / "results/E2_COVERAGE_STRESS.json"
    if destination.exists():
        raise FileExistsError("Preserve completed coverage diagnostic; use a new version instead of overwriting")
    cfg = json.loads((ROOT / "configs/E2_COVERAGE_STRESS_FROZEN.json").read_text())
    with np.load(ROOT / "results/E1_RESPONSE_BANK.npz", allow_pickle=False) as bank:
        physical = np.stack([v.ravel() for v in np.meshgrid(
            *[bank[f"axis_{j}"] for j in range(4)], indexing="ij")], axis=1)
        candidates = np.array([from_geometry(*v) for v in physical])
        fields = np.stack([bank[f"fields_{k:g}"].reshape(-1, 32) for k in FREQS], axis=1)
    with np.load(ROOT / "results/E1_NOISE_ANCHORS.npz", allow_pickle=False) as anchors:
        sigma = anchors["sigma"].copy()
        calibration_sd = float(sigma[1] / np.linalg.norm(anchors["F"][1]))
    masks = support_masks(candidates)
    simulator = OfflineSimulator(ROOT / "results/E1_FRESH_SIMULATOR.npz")
    rows, recommendations = [], []
    for parent in range(32):
        cost.check_time()
        for noise in range(2):
            old_y = simulator.initial(parent, noise)
            old_sd = np.full(32, sigma[1])
            for group in GROUPS:
                retained = masks[group]
                x, ff = candidates[retained], fields[retained]
                # This is rerun after each deletion; no original-full-set picks.
                recommendation = recommend_retained(x, ff, old_y, sigma, calibration_sd)
                recommendations.append(dict(parent=parent, noise=noise, support_group=group,
                                            **recommendation))
                selected_actions = {"none": None, "fixed2.4": "frequency_2.4",
                                    "fixed3.2": "frequency_3.2",
                                    "raw": recommendation["picks"]["raw"],
                                    "robust": recommendation["picks"]["robust"],
                                    "sentinel_first3.2": "frequency_3.2"}
                actions = {a.id: a for a in menu(ff, sigma, calibration_sd)}
                initial = fitted_selection(x, ff[:, 1], old_y, old_sd)
                # Truth metadata is requested only after acquisition is complete.
                offline = simulator.score_catalogue(parent, x)
                for scenario in SCENARIOS:
                    fitted = {None: initial}
                    for aid in dict.fromkeys(selected_actions.values()):
                        if aid is None:
                            continue
                        action = actions[aid]
                        new_y = simulator.measure(parent, noise, aid, scenario)
                        fitted[aid] = fitted_selection(x, np.column_stack((ff[:, 1], action.fields)),
                                                      np.r_[old_y, new_y], np.r_[old_sd, action.sigma])
                    for policy, aid in selected_actions.items():
                        fit = fitted[aid]
                        chosen = fit["chosen"]
                        selected_family = int(family_ids(x[chosen])[0])
                        selected_event = int(event_ids(x[chosen]))
                        missing = not offline["true_family_represented"]
                        accepted = not (policy == "sentinel_first3.2" and fit["sentinel"])
                        rows.append(dict(
                            scenario=scenario, support_group=group, parent=parent, noise=noise,
                            policy=policy, selected_action=aid, K=len(x),
                            selected_MAE=float(offline["losses"][chosen]),
                            initial_selected_MAE=float(offline["losses"][initial["chosen"]]),
                            **{k: v for k, v in offline.items() if k != "losses"},
                            selected_event_id=selected_event, selected_family=FAMILIES[selected_family],
                            selected_family_correct=bool(selected_family == offline["true_family_id"]),
                            selected_structure_correct=bool(selected_event == offline["true_event_id"]),
                            selected_gain=fit["gain"], minimum_joint_profile_chi2=fit["chi2"],
                            residual_sentinel=fit["sentinel"], sentinel_nominal_dof=fit["dof"],
                            sentinel_threshold=fit["threshold"],
                            compatible_candidate_count=fit["compatible_candidates"],
                            compatible_family_ids=fit["compatible_families"],
                            undetected_missing_family=bool(missing and not fit["sentinel"]),
                            apparent_wrong_single_family=bool(missing and len(fit["compatible_families"]) == 1
                                                             and not fit["sentinel"]),
                            decision_accepted=accepted,
                            decision_status="ABSTAIN_PRESERVE_AMBIGUITY" if not accepted else "PROFILE_MATERIAL_PROPOSAL",
                            measurement_blocks=int(aid is not None),
                            history_real_channels=fit["dof"]+1,
                            candidate_weights_kind="none", posterior_probability=None,
                            sentinel_status="heuristic only; finite-set selection and gain boundaries invalidate exact chi-square calibration",
                            new_Maxwell_attempts=0, generator_refresh=False,
                            exposure="existing E1 parents reused for explicitly exposed E2 diagnosis"))
        if parent % 8 == 7:
            print("coverage diagnostic parent", parent+1, "of32; new Maxwell0", flush=True)
    parent_rows, aggregates, paired_none, paired_full = parent_statistics(
        rows, cfg["statistical_seed"], cfg["bootstrap"])
    out = dict(protocol=cfg, status="EXECUTED_SAVED_FULL_FIELD_DIAGNOSIS", new_Maxwell_attempts=0,
               physics_calls_owned_by_this_job=0,
               hardware_measurements=0, GPU_training=0,
               source_results=["E1_RESPONSE_BANK.npz", "E1_FRESH_SIMULATOR.npz", "E1_NOISE_ANCHORS.npz"],
               support_counts={k: int(v.sum()) for k, v in masks.items()},
               rows=rows, recommendations=recommendations, parent_rows=parent_rows,
               aggregates=aggregates, paired_vs_no_measurement=paired_none,
               paired_vs_full_support=paired_full,
               interpretation={
                   "best_of_K": "fixed catalogue: invariant to policy and future instrument gain; deletion changes available support",
                   "sentinel_baseline": "same action/profile as fixed3.2; its distinct output is diagnostic abstention, not a different measurement ranking",
                   "sentinel_selected_MAE": "unconditional proposed-material error still reported even when diagnostic abstention rejects installation",
                   "mismatch": "intentional frequency-dependent instrument gain, not a new physical/material model and never pooled with in-model scores",
                   "no_probability": "raw/robust are finite-set heuristics, profile residual is not a posterior or coverage guarantee",
                   "statistical_unit": "two episodes averaged within each of32 independent parents"},
               body_wall_seconds=time.perf_counter()-started,
               body_cpu_seconds=time.process_time()-cpu,
               body_timer_scope="excludes final JSON export; Cost job includes export and full execution")
    dump(destination, out)
    print(json.dumps({"written": str(destination), "rows": len(rows),
                      "support_counts": out["support_counts"], "new_Maxwell_attempts": 0,
                      "body_wall_seconds": out["body_wall_seconds"]}), flush=True)
    return out


if __name__ == "__main__":
    cost = Cost()
    with cost.job("e2-coverage-stress-v1", "E2"):
        run(cost)
