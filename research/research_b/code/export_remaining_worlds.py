"""Export all saved 625 worlds and final data-compatible masks, no physics.

The mask is a finite-catalogue heuristic, not a posterior credible region.
The observation reader replays existing Simulator potential outcomes only;
the screening function receives no true material, true gain or offline score.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from src.b1_joint.material import from_geometry, geometry, event_ids, EVENT_LABELS
from active_experiment import Simulator, FREQS, IDS
from e0_candidate_controls import profile
from cost import Cost, ROOT, dump

POLICIES = ("fixed_3.2", "robust", "expected_material_risk")
TRAJECTORY_FIELDS = ("parent", "noise", "policy", "step", "selected_action",
                     "history", "candidate_profile_chi2", "measurement_blocks")
FAMILY_LABELS = ("solid", "central_cavity", "annular_cavity", "central_and_annular_cavities")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_trajectory(path):
    """Whitelist online history/fit fields; never pass offline MAE/truth scores."""
    raw = json.loads(path.read_text())
    groups = defaultdict(list)
    for row in raw:
        if row["policy"] in POLICIES:
            online = {key: row[key] for key in TRAJECTORY_FIELDS}
            groups[(online["policy"], online["parent"], online["noise"])].append(online)
    for key, rows in groups.items():
        rows.sort(key=lambda r: r["step"])
        require([r["step"] for r in rows] == list(range(len(rows))), f"Non-contiguous saved steps: {key}")
        require(rows[0]["selected_action"] is None and rows[0]["history"] == ["repeat_1.6"],
                "Original observation history must occur exactly once")
    return groups


def observation_only_simulator(path):
    """Retain only saved observation arrays; truth/gain/score assets unavailable."""
    simulator = Simulator(path)
    observations = {key: simulator.data[key].copy() for key in ("outcomes", "calibration")}
    simulator.data.close()
    simulator.data = observations
    return simulator


def compatible_worlds(fields, observation, sigma):
    """No probabilities: each world gets its own g, shared over the full history."""
    require(fields.ndim == 2 and fields.shape[1] == len(observation) == len(sigma),
            "Joint-history dimensions differ")
    require(np.all(np.isfinite(fields)) and np.all(np.isfinite(observation)) and
            np.all(np.isfinite(sigma)) and np.all(sigma > 0), "Non-finite data or invalid noise")
    require(np.all(np.sum((fields / sigma) ** 2, axis=1) > 0), "Unidentifiable zero-response gain")
    chi2, gain = profile(fields, observation, sigma)
    best = float(chi2.min())
    mask = chi2 <= best + 9.
    require(mask.any() and np.all((gain >= .7) & (gain <= 1.3)), "Invalid bounded-gain compatibility set")
    return chi2, gain, mask, best


def run(cost):
    destination = ROOT / "results"
    npz_path, json_path = destination / "REMAINING_WORLDS.npz", destination / "REMAINING_WORLDS.json"
    require(not npz_path.exists() and not json_path.exists(), "Export exists; preserve it and use an explicitly new revision")
    trajectory = read_trajectory(destination / "E2_POLICY_RESULTS.json")
    parents = sorted({key[1] for key in trajectory})
    noises = sorted({key[2] for key in trajectory})
    require(len(parents) == 32 and len(noises) == 2 and len(trajectory) == 3 * 32 * 2,
            "Expected complete three-policy, 32-parent, two-noise source")
    simulator = observation_only_simulator(destination / "E1_FRESH_SIMULATOR.npz")
    with np.load(destination / "E1_RESPONSE_BANK.npz", allow_pickle=False) as bank:
        axes = [bank[f"axis_{j}"].copy() for j in range(4)]
        physical = np.stack([v.ravel(order="C") for v in np.meshgrid(*axes, indexing="ij")], axis=1)
        fields = np.stack([bank[f"fields_{frequency:g}"].reshape(-1, 32, order="C")
                           for frequency in FREQS], axis=1)
    latent = np.asarray([from_geometry(*p) for p in physical])
    require(latent.shape == (625, 4) and fields.shape == (625, 6, 32), "Not the saved full625 catalogue")
    effective_geometry = geometry(latent)
    families = ((effective_geometry[:, 0] > 0).astype(np.int64) +
                2 * (effective_geometry[:, 1] > 0).astype(np.int64))
    events = event_ids(latent)
    with np.load(destination / "E1_NOISE_ANCHORS.npz", allow_pickle=False) as anchors:
        sigma_per_frequency = anchors["sigma"].copy()
        calibration_sd = float(sigma_per_frequency[1] / np.linalg.norm(anchors["F"][1]))
    shape = (len(POLICIES), len(parents), len(noises), 625)
    chi2_all, gains_all = np.zeros(shape), np.zeros(shape)
    masks = np.zeros(shape, dtype=bool)
    counts = np.zeros(shape[:-1], dtype=np.int64)
    family_counts = np.zeros(shape[:-1] + (4,), dtype=np.int64)
    event_counts = np.zeros(shape[:-1] + (12,), dtype=np.int64)
    history_values = np.zeros(shape[:-1] + (128,))
    history_sigmas = np.zeros_like(history_values)
    history_lengths = np.zeros(shape[:-1], dtype=np.int64)
    history_actions = np.full(shape[:-1] + (4,), "", dtype="U40")
    history_obs_ids = np.full(shape[:-1] + (4,), "", dtype="U96")
    block_lengths = np.zeros(shape[:-1] + (4,), dtype=np.int64)
    records = []
    max_fit_difference = 0.
    steps_verified = 0
    for pi, policy in enumerate(POLICIES):
        for parent_i, parent in enumerate(parents):
            cost.check_time()
            for noise_i, noise in enumerate(noises):
                rows = trajectory[(policy, parent, noise)]
                old = simulator.initial(parent, noise)
                yy = old.copy()
                sd = np.full(32, sigma_per_frequency[FREQS.index(1.6)])
                Fhistory = fields[:, FREQS.index(1.6)].copy()
                actions = ["repeat_1.6"]
                observation_ids = [f"parent-{parent}:noise-{noise}:original-kR1.6:outcome-index-0"]
                lengths = [32]
                for row in rows:
                    step = row["step"]
                    if step:
                        aid = row["selected_action"]
                        require(aid in IDS or aid == "reference_calibration", "Unexpected non-Maxwell/external-context action")
                        observation = simulator.measure(parent, noise, step - 1, aid)
                        if aid == "reference_calibration":
                            added_fields = np.ones((625, 1))
                            added_sd = np.array([calibration_sd])
                            obs_id = f"parent-{parent}:noise-{noise}:calibration-index-{step - 1}"
                        else:
                            frequency_i = IDS.index(aid)
                            added_fields = fields[:, frequency_i]
                            added_sd = np.full(32, sigma_per_frequency[frequency_i])
                            obs_id = f"parent-{parent}:noise-{noise}:{aid}:outcome-index-{step}"
                        Fhistory = np.column_stack((Fhistory, added_fields))
                        yy, sd = np.r_[yy, observation], np.r_[sd, added_sd]
                        actions.append(aid)
                        observation_ids.append(obs_id)
                        lengths.append(len(observation))
                    require(actions == row["history"] and len(set(observation_ids)) == len(observation_ids),
                            "Saved/replayed history differs or repeats the same raw noise observation")
                    chi2, gain, mask, best = compatible_worlds(Fhistory, yy, sd)
                    difference = abs(best - row["candidate_profile_chi2"])
                    max_fit_difference = max(max_fit_difference, difference)
                    require(np.isclose(best, row["candidate_profile_chi2"], rtol=1e-10, atol=1e-9),
                            f"Replayed profile differs from original trajectory: {policy}/{parent}/{noise}/{step}")
                    steps_verified += 1
                index = (pi, parent_i, noise_i)
                chi2_all[index], gains_all[index], masks[index] = chi2, gain, mask
                counts[index] = int(mask.sum())
                family_counts[index] = np.bincount(families[mask], minlength=4)
                event_counts[index] = np.bincount(events[mask], minlength=12)
                history_lengths[index] = len(yy)
                history_values[index][:len(yy)], history_sigmas[index][:len(sd)] = yy, sd
                history_actions[index][:len(actions)] = actions
                history_obs_ids[index][:len(actions)] = observation_ids
                block_lengths[index][:len(lengths)] = lengths
                present_families = np.flatnonzero(family_counts[index])
                present_events = np.flatnonzero(event_counts[index])
                selected = int(np.argmin(chi2))
                records.append({"policy": policy, "parent": parent, "noise": noise,
                    "npz_index": list(index), "catalogue_size_retained": 625,
                    "compatible_count": int(mask.sum()), "selected_catalogue_index": selected,
                    "best_joint_profile_chi2": best, "compatibility_cutoff": best + 9.,
                    "family_ids": present_families.tolist(),
                    "family_labels": [FAMILY_LABELS[j] for j in present_families],
                    "posthoc_event_ids": present_events.tolist(),
                    "multiple_retained_material_worlds": bool(mask.sum() > 1),
                    "unresolved_between_posthoc_families": bool(len(present_families) > 1),
                    "unresolved_between_posthoc_events": bool(len(present_events) > 1),
                    "compatible_physical_coordinate_min": physical[mask].min(axis=0).tolist(),
                    "compatible_physical_coordinate_max": physical[mask].max(axis=0).tolist(),
                    "measurement_history": [{"action": a, "observation_id": oid, "n_real_channels": n}
                        for a, oid, n in zip(actions, observation_ids, lengths)],
                    "new_measurement_blocks_in_original_episode": len(actions) - 1,
                    "acquisition_input": "reference-proxy-assisted action choice" if policy == "expected_material_risk" else "candidate-only or fixed menu",
                    "candidate_probability_kind": "none", "posterior_probability": None,
                    "represented_probability_mass": None, "outside_support_mass": None,
                    "unexplored_material_worlds": "unquantified; mask membership is not a complete-support guarantee"})
    np.savez_compressed(npz_path, policies=np.asarray(POLICIES), parent_ids=np.asarray(parents),
        noise_ids=np.asarray(noises), catalogue_latent=latent, catalogue_physical_coordinates=physical,
        catalogue_effective_geometry=effective_geometry, catalogue_shared_epsilon=latent[:, 2] + 1j * latent[:, 3],
        catalogue_posthoc_family_ids=families, catalogue_posthoc_event_ids=events,
        compatible_mask=masks, joint_profile_chi2=chi2_all, per_world_shared_gain=gains_all,
        compatible_counts=counts, compatible_family_counts=family_counts, compatible_event_counts=event_counts,
        history_values=history_values, history_sigmas=history_sigmas,
        history_valid_real_channels=history_lengths, history_action_ids=history_actions,
        history_observation_ids=history_obs_ids, history_block_real_lengths=block_lengths)
    summary = {"package_version": "research-b-remaining-worlds/1.0.0", "status": "EXPORTED_SAVED_FULL625_CATALOGUE_AND_HEURISTIC_MASKS",
        "array_file": npz_path.name, "catalogue_size": 625, "policies": list(POLICIES), "parents": 32,
        "noise_per_parent": 2, "episode_policy_count": len(records),
        "criterion": "per-world bounded shared-gain joint-history profile chi2 <= best profile chi2 + 9",
        "compatibility_set_status": "HEURISTIC_DATA_COMPATIBILITY_SET_NOT_POSTERIOR_NOT_CONFIDENCE_REGION",
        "gain_bounds": [.7, 1.3], "gain_rule": "Different worlds may have different gains; each world's single gain is shared across all its old/new frequency/calibration blocks",
        "probability_kind": "none", "posterior_probabilities": None, "outside_support_mass": None,
        "candidate_counts_are_probabilities": False, "full_catalogue_retained": True,
        "family_definition": dict(enumerate(FAMILY_LABELS)), "event_definition": dict(enumerate(EVENT_LABELS)),
        "material_definition": "Exact original P2 concentric material decoder from four continuous coefficients; posthoc labels do not condition generation. R=1; relative epsilon dimensionless.",
        "source_forward_bank": {"relative_path": "E1_RESPONSE_BANK.npz", "fields": "full saved 625-by-6-by-32 Maxwell response bank, C flattening; required external evidence asset in same release"},
        "source_trajectory": "E2_POLICY_RESULTS.json", "source_observations": "E1_FRESH_SIMULATOR.npz outcomes/calibration only",
        "trajectory_fields_used": list(TRAJECTORY_FIELDS), "offline_truth_or_true_gain_used_by_filter": False,
        "observation_replay": "Exactly original Simulator.initial and measure(parent,noise,step-1,action) potential outcomes; no new noise drawn; original index0 occurs once",
        "array_axes": {"masks_and_profiles": ["policy", "parent", "noise", "catalogue_world"],
                       "history_values": "Zero padded; slice to history_valid_real_channels. Original and fresh blocks are concatenated once, each block Re16/Im16 or one real calibration channel."},
        "numerical_replay_checks": {"saved_steps_verified": steps_verified, "max_profile_chi2_absolute_difference": max_fit_difference},
        "new_maxwell_attempts": 0, "new_observation_draws": 0, "new_candidate_generation": False,
        "limitations": ["Finite coarse full625 catalogue may miss important worlds.", "Heuristic delta9 mask has no calibrated coverage probability.",
            "Outside-support and undiscovered-family mass is unquantified.", "Different policies observe different action-indexed potential outcomes; they are not extra independent parents.",
            "E2 reuses E1 parents as exposed followup.", "Retaining one event does not prove a globally unique material; original complete catalogue is never discarded."],
        "episodes": records}
    dump(json_path, summary)
    print(json.dumps({"status": summary["status"], "arrays": str(npz_path), "metadata": str(json_path),
        "episodes": len(records), "retained_catalogue": 625, "compatible_count_min_max": [int(counts.min()), int(counts.max())],
        "multiple_families_episodes": sum(r["unresolved_between_posthoc_families"] for r in records),
        "steps_verified": steps_verified, "max_profile_difference": max_fit_difference, "new_maxwell_attempts": 0}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", default="e2-export-remaining-worlds-v1")
    args = parser.parse_args()
    with Cost().job(args.job, "E2") as cost:
        run(cost)
