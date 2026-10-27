"""Read-only audit and additional offline scores for sealed exposed finite24.

No online method is called or modified. The global and every case/method seal
are verified before any old truth archive is opened. Raw bound violations are
reported with their original sign, with no clipping, tolerance or PASS rescue.
"""

from pathlib import Path
import argparse
from datetime import datetime, timezone
import csv
import json
import os
import resource
import sys
import time

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
METHOD_IDS = tuple(f"M{i}" for i in range(7))
CASES = tuple((parent, repeat) for parent in range(12) for repeat in range(2))
LABEL = "REPLAY_NOT_NEW_PDE_SPEED_EVIDENCE"
SIGMA = 0.01
OUTER_REAL_DIMENSION = 64  # 2 TX x16 RX complex =64 real quadratures
RSS_LIMIT_BYTES = 1_000_000_000
GRID_KEYS = ("support", "epsinf", "ell", "pose", "prior_weight", "prior_id")
EXPECTED_PRIOR_ID = "FINITE_PRIOR_ENUM4_V1"
OLD_INPUT_FILENAMES = {"FINITE_PRIOR_GRID": "FINITE_GRID.npz",
    "FINITE_FULL_RESPONSE_CACHE": "FINITE_BANK.npz", "FINITE_24_DEVELOPMENT_CASES": "OBSERVATIONS.npz"}
FEATURE_NAMES = tuple([f"wall_cell_{i}_Re_chi" for i in range(4)]
                      + [f"wall_cell_{i}_Im_chi" for i in range(4)] + ["pose_x", "pose_y"])
CP_ARRAY_KEYS = ("queried", "log_likelihood_lower", "log_likelihood_upper",
                 "posterior_lower", "posterior_upper")


def jsonable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, np.generic):
        return jsonable(value.item())
    if isinstance(value, complex):
        return [value.real, value.imag]
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        return "-Infinity" if value == -np.inf else "Infinity" if value == np.inf else "NaN"
    return value


def write_json(path, value):
    with Path(path).open("x") as stream:
        json.dump(jsonable(value), stream, indent=2, allow_nan=False)
        stream.write("\n")


def write_csv(path, rows):
    if not rows:
        raise ValueError("explicit rows required for CSV output")
    with Path(path).open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: jsonable(value) for key, value in row.items()})


def scalar(value):
    return -np.inf if value == "-Infinity" else np.inf if value == "Infinity" else float(value)


def safe_source(path):
    path = Path(path).resolve()
    if any("locked" in part.lower() or "posteriors" in part.lower()
           or "inference" in part.lower() for part in path.parts):
        raise ValueError("locked/inference/POSTERIORS sources are forbidden")
    return path


def resolve_old_sources(cache_map, old_bank_root=None):
    """Only explicit relocation can replace CACHE_MAP's original old root."""
    original_root = safe_source(Path(cache_map["roots"]["raw_WSLAM"])
                                / "results/e2e3/finite_bank/online")
    selected = original_root if old_bank_root is None else safe_source(old_bank_root)
    entries = {row["id"]: row for row in cache_map["files"]}
    paths, originals = {}, {}
    for identifier, filename in OLD_INPUT_FILENAMES.items():
        original = safe_source(entries[identifier]["path"])
        if original != (original_root / filename).resolve():
            raise ValueError("old source differs from CACHE_MAP's exact original allowlist")
        paths[identifier] = safe_source(selected / filename)
        originals[identifier] = original
    original_truth = safe_source(entries["FINITE_OFFLINE_TRUTHS_NOISE"]["path"])
    if original_truth != (original_root.parent / "offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz").resolve():
        raise ValueError("old truth source is not the original finite-bank sibling offline archive")
    paths["FINITE_OFFLINE_TRUTHS_NOISE"] = safe_source(selected.parent / "offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz")
    originals["FINITE_OFFLINE_TRUTHS_NOISE"] = original_truth
    return {"paths": paths, "original_audited_paths": originals,
        "explicit_relocation": old_bank_root is not None, "requested_online_root": old_bank_root,
        "resolved_online_root": selected,
        "policy": "original CACHE_MAP root by default; explicit online root changes only exact allowlisted filenames; no source search; truth remains a sibling offline archive"}


def validated_prior_id(value, source):
    value = np.asarray(value)
    if value.shape != () or value.dtype.kind not in "US" or str(value.item()) != EXPECTED_PRIOR_ID:
        raise ValueError(f"{source}: expected scalar prior_id={EXPECTED_PRIOR_ID}")
    return EXPECTED_PRIOR_ID


class AccessJournal:
    """Persist access intent before opening a source, including truth intent."""

    def __init__(self, output):
        self.path = output / "SOURCE_ACCESS.json"
        self.events = []
        self.seal_attestation = None
        self.old_sources = None
        self.flush()

    def flush(self):
        self.path.write_text(json.dumps(jsonable({
            "schema": "BGGA_G1_OFFLINE_SOURCE_ACCESS_V1", "label": LABEL,
            "online_source_files_modified": False, "new_truth_labels_created": False,
            "seal_attestation": self.seal_attestation, "events": self.events,
            "old_bank_sources_and_explicit_relocation": self.old_sources,
        }), indent=2, allow_nan=False) + "\n")

    def record(self, action, path=None, keys=None, **metadata):
        self.events.append({"sequence": len(self.events),
            "time_utc": datetime.now(timezone.utc).isoformat(), "action": action,
            "path": None if path is None else str(path), "keys": keys,
            "all_online_outputs_verified_before_access": self.seal_attestation is not None,
            **metadata})
        self.flush()

    def read_json(self, path):
        path = safe_source(path)
        self.record("READ_JSON_INTENT", path)
        return json.loads(path.read_text())

    def load_npz(self, path, keys, *, truth=False):
        path = safe_source(path)
        if self.seal_attestation is None:
            raise ValueError("NPZ payloads require a verified global online seal")
        if truth and (self.old_sources is None
                      or path != self.old_sources["paths"]["FINITE_OFFLINE_TRUTHS_NOISE"]
                      or tuple(keys) != ("truth_index",)):
            raise ValueError("only the fixed old truth_index array may be opened")
        self.record("TRUTH_NPZ_OPEN_INTENT" if truth else "NPZ_OPEN_INTENT", path, list(keys),
                    old_existing_truth_only=truth)
        with np.load(path, allow_pickle=False) as archive:
            values = {key: np.array(archive[key], copy=True) for key in keys}
        self.record("TRUTH_NPZ_READ_COMPLETED" if truth else "NPZ_READ_COMPLETED", path, list(keys),
            array_shapes={key: value.shape for key, value in values.items()},
            numeric_array_bytes=sum(value.nbytes for value in values.values()))
        return values


def verify_sealed_replay(replay, journal):
    """File inventory/metadata barrier only; it does not declare scientific PASS."""
    global_seal = journal.read_json(replay / "ALL_ONLINE_OUTPUTS_SEALED.json")
    if (global_seal.get("all_method_outputs_written") is not True
            or global_seal.get("case_count") != 24
            or global_seal.get("method_count_per_case") != 7
            or global_seal.get("reference_input_to_online_methods") is not False):
        raise ValueError("global seal is incomplete or permits reference action input")
    plan = journal.read_json(replay / "PLAN.json")
    receipt = journal.read_json(replay / "RUN_RECEIPT.json")
    if plan.get("label") != LABEL or receipt.get("label") != LABEL:
        raise ValueError("unexpected replay evidence label")
    metadata = {}
    verified_files = []
    for parent, repeat in CASES:
        directory = replay / f"parent_{parent:02d}_repeat_{repeat:02d}"
        seal_path = directory / "ONLINE_OUTPUTS_SEALED.json"
        seal = journal.read_json(seal_path)
        if (seal.get("parent_index"), seal.get("repeat_index"), seal.get("method_count")) != (parent, repeat, 7):
            raise ValueError("case seal identity/count differs from the fixed cohort")
        required = [directory / "JOINT_METHOD_OUTPUTS.npz", directory / "OFFLINE_REFERENCE.npz",
                    directory / "OFFLINE_REFERENCE.json", directory / "OFFLINE_SCORES.json"]
        for method_id in METHOD_IDS:
            required.extend([directory / (method_id + ".npz"), directory / (method_id + ".json")])
            row = journal.read_json(directory / (method_id + ".json"))
            if row.get("method_id") != method_id or row.get("label") != LABEL:
                raise ValueError("method identity/label differs from sealed replay")
            metadata[(parent, repeat, method_id)] = row
        if any(not path.is_file() for path in required):
            raise ValueError("a sealed method/reference file is missing")
        verified_files.extend([str(seal_path), *map(str, required)])
    attestation = {"status": "SEALED_FILE_INVENTORY_AND_METADATA_VERIFIED",
        "replay": str(replay), "case_count": 24, "method_count_per_case": 7,
        "verified_files": verified_files, "hash_check_performed": False,
        "scope": "presence and recorded identity, not numerical or scientific acceptance",
        "raw_interval_bound_comparison_status_unchanged": receipt.get("raw_interval_bound_comparison_status")}
    journal.seal_attestation = attestation
    journal.record("ALL_ONLINE_METHODS_AND_OFFLINE_REFERENCES_VERIFIED_BEFORE_TRUTH")
    return metadata, receipt


def peak_rss():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def check_memory():
    value = peak_rss()
    if value > RSS_LIMIT_BYTES:
        raise MemoryError(f"offline peak RSS {value} exceeds {RSS_LIMIT_BYTES}")
    return value


def signed_array_stats(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not values.size or np.any(np.isnan(values)):
        raise ValueError("nonempty signed values without NaN required")
    positive = np.flatnonzero(values > 0)
    return {"minimum_signed": float(np.min(values)), "maximum_signed": float(np.max(values)),
            "argmax_atom_index": int(np.argmax(values)), "positive_count": int(positive.size)}, positive


def audit_intervals(replay, output, journal, metadata):
    """Recompute each raw violation from immutable M4/reference payloads."""
    rows, any_positive = [], False
    global_vectors = {name: [] for name in ("atom_lower", "atom_upper", "likelihood_lower", "likelihood_upper")}
    for parent, repeat in CASES:
        directory = replay / f"parent_{parent:02d}_repeat_{repeat:02d}"
        reference = journal.load_npz(directory / "OFFLINE_REFERENCE.npz", ("joint_weights", "log_likelihood"))
        reference_meta = journal.read_json(directory / "OFFLINE_REFERENCE.json")
        q, ll = reference["joint_weights"], reference["log_likelihood"]
        if q.shape != (2950,) or ll.shape != q.shape or np.any(~np.isfinite(q)) or np.any(~np.isfinite(ll)):
            raise ValueError("finite sealed 2950-atom reference required")
        m4_meta = metadata[(parent, repeat, "M4")]
        keys = tuple(row["npz_array_prefix"] + name for row in m4_meta["checkpoints"] for name in CP_ARRAY_KEYS)
        saved = journal.load_npz(directory / "M4.npz", keys)
        arrays, case_rows = {"reference_joint_weights": q, "reference_log_likelihood": ll}, []
        for ordinal, cp in enumerate(m4_meta["checkpoints"]):
            prefix = cp["npz_array_prefix"]
            if prefix != f"checkpoint_{ordinal}_":
                raise ValueError("unexpected checkpoint prefix")
            values = {name: saved[prefix + name] for name in CP_ARRAY_KEYS}
            if any(value.shape != q.shape for value in values.values()) or values["queried"].dtype.kind != "b":
                raise ValueError("checkpoint arrays must match the original atom grid")
            queried = values["queried"]
            if int(queried.sum()) != cp["query_count"]:
                raise ValueError("checkpoint query count and mask differ")
            signed = {"atom_lower": values["posterior_lower"] - q,
                      "atom_upper": q - values["posterior_upper"],
                      "likelihood_lower": values["log_likelihood_lower"] - ll,
                      "likelihood_upper": ll - values["log_likelihood_upper"]}
            actual_unknown = float(q[~queried].sum())
            unknown_lower = scalar(cp["unknown_mass_lower"]) - actual_unknown
            unknown_upper = actual_unknown - scalar(cp["unknown_mass_upper"])
            logz = scalar(reference_meta["logZ_exact"])
            row = {"parent_index": parent, "repeat_index": repeat, "query_count": cp["query_count"],
                "actual_reference_unknown_mass": actual_unknown,
                "unknown_lower_signed_violation": unknown_lower,
                "unknown_upper_signed_violation": unknown_upper,
                "log_evidence_lower_signed_violation": scalar(cp["log_evidence_lower"]) - logz,
                "log_evidence_upper_signed_violation": logz - scalar(cp["log_evidence_upper"]),
                "unknown_likelihood_interval": "[0,1]", "floating_point_roundoff_certified": False}
            for name, vector in signed.items():
                stats, positive = signed_array_stats(vector)
                arrays[prefix + name + "_signed_violation"] = vector
                arrays[prefix + name + "_positive_atom_indices"] = positive
                arrays[prefix + name + "_positive_signed_values"] = vector[positive]
                for stat, value in stats.items():
                    row[name + "_" + stat] = value
                global_vectors[name].append(stats["maximum_signed"])
            # Save source interval values as well as every signed difference.
            arrays.update({prefix + name: value for name, value in values.items()})
            row["any_positive_raw_violation"] = bool(unknown_lower > 0 or unknown_upper > 0
                or row["log_evidence_lower_signed_violation"] > 0
                or row["log_evidence_upper_signed_violation"] > 0
                or any(np.any(vector > 0) for vector in signed.values()))
            any_positive |= row["any_positive_raw_violation"]
            case_rows.append(row)
            rows.append(row)
        np.savez_compressed(output / f"interval_parent_{parent:02d}_repeat_{repeat:02d}.npz", **arrays)
        write_json(output / f"interval_parent_{parent:02d}_repeat_{repeat:02d}.json", case_rows)
        check_memory()
    summary = {"label": LABEL, "case_count": 24, "checkpoint_rows": len(rows),
        "status": "RAW_POSITIVE_VIOLATIONS_PRESENT" if any_positive else "NO_RAW_POSITIVE_OBSERVED",
        "scientific_gate_decision": "PARENT_REVIEW_REQUIRED_NO_PASS_RESCUE",
        "roundoff_certified": False, "clipping": False, "added_tolerance": False,
        "maximum_signed_by_array": {name: max(values) for name, values in global_vectors.items()},
        "maximum_unknown_lower_signed_violation": max(row["unknown_lower_signed_violation"] for row in rows),
        "maximum_unknown_upper_signed_violation": max(row["unknown_upper_signed_violation"] for row in rows),
        "maximum_log_evidence_lower_signed_violation": max(row["log_evidence_lower_signed_violation"] for row in rows),
        "maximum_log_evidence_upper_signed_violation": max(row["log_evidence_upper_signed_violation"] for row in rows),
        "positive_checkpoint_rows": sum(row["any_positive_raw_violation"] for row in rows),
        "interpretation": "signed float comparisons only; magnitude alone does not establish a roundoff proof"}
    write_json(output / "INTERVAL_AUDIT_SUMMARY.json", summary)
    write_csv(output / "INTERVAL_CHECKPOINT_AUDIT.csv", rows)
    return summary


def logsumexp(values):
    offset = float(np.max(values))
    if np.isneginf(offset):
        return offset
    return float(offset + np.log(np.exp(values - offset).sum()))


def feature_table(grid):
    support = np.asarray(grid["support"])
    eps, ell, pose = grid["epsinf"], grid["ell"], grid["pose"]
    if (support.shape != (2950, 4) or support.dtype.kind != "b" or eps.shape != (2950,)
            or ell.shape != (2950,) or pose.shape != (2950, 2)):
        raise ValueError("finite four-cell wall/pose grid required")
    prior_weight = grid["prior_weight"]
    if prior_weight.shape != (2950,) or np.any(~np.isfinite(prior_weight)) or np.any(prior_weight <= 0):
        raise ValueError("positive finite prior weights on the same 2950 atoms required")
    chi = support * ((eps - 1) + 1j * ell)[:, None]
    features = np.column_stack((chi.real, chi.imag, pose))
    if np.any(~np.isfinite(features)):
        raise ValueError("finite physical feature table required")
    return chi, features


def moments(weights, features):
    """Use saved floating weights without clipping/reweighting online results."""
    weights = np.asarray(weights, dtype=float)
    if weights.shape != (2950,) or np.any(~np.isfinite(weights)) or np.any(weights < 0) or not np.any(weights > 0):
        raise ValueError("valid saved 2950-atom joint measure required")
    mean = weights @ features
    centered = features - mean
    covariance = (centered.T * weights) @ centered
    return mean, covariance


def outer_predictive(bank, y_outer, log_weights, counter, owner):
    """Saved TX2/3 predictive density: all 64 real Gaussian constants included.

    This is evaluation of old exposed unassimilated observations. It creates no
    factor, acquisition, new observation, or new PDE response. Every scored
    method/reference is independently charged 2950 offline saved-response rows.
    """
    if bank.shape != (2950, 4, 16) or y_outer.shape != (2, 16) or log_weights.shape != (2950,):
        raise ValueError("frozen outer predictive tensor shapes required")
    if np.any(np.isnan(log_weights)) or np.any(np.isposinf(log_weights)):
        raise ValueError("finite/-inf saved log joint weights required")
    difference = (bank[:, 2:4, :] - y_outer) / SIGMA
    phi = 0.5 * np.sum(difference.real ** 2 + difference.imag ** 2, axis=(1, 2))
    log_constant = -0.5 * OUTER_REAL_DIMENSION * np.log(2 * np.pi * SIGMA ** 2)
    counter[owner] = counter.get(owner, 0) + bank.shape[0]
    return float(logsumexp(log_weights - phi) + log_constant)


def measure_scores(weights, log_weights, chi, features, truth_atom, reference_cov,
                   bank, y_outer, counter, owner):
    mean, covariance = moments(weights, features)
    mean_chi = weights @ chi
    truth_chi, truth_pose = chi[truth_atom], features[truth_atom, 8:]
    delta = covariance - reference_cov
    cross, cross_delta = covariance[:8, 8:], delta[:8, 8:]
    score = {"GT_position_squared_error": float(np.sum((mean[8:] - truth_pose) ** 2)),
        "equal_area_four_cell_wall_chi_MSE": float(np.mean(np.abs(mean_chi - truth_chi) ** 2)),
        "joint_cov_difference_Frobenius": float(np.linalg.norm(delta, "fro")),
        "joint_cov_difference_max_abs_entry": float(np.max(np.abs(delta))),
        "wall_target_cross_cov_difference_Frobenius": float(np.linalg.norm(cross_delta, "fro")),
        "wall_target_cross_cov_difference_max_abs_entry": float(np.max(np.abs(cross_delta))),
        "outer_TX2_3_log_posterior_predictive_density": outer_predictive(bank, y_outer, log_weights, counter, owner),
        "saved_joint_weight_sum": float(np.sum(weights)),
        "saved_log_joint_weight_logsumexp": logsumexp(log_weights),
        "moment_weight_policy": "AS_SAVED_FLOAT_WEIGHTS_NO_RENORMALIZATION",
        "logical_offline_response_atom_queries": 2950, "actual_new_full_FE": 0, "actual_new_RHS": 0}
    arrays = {"joint_mean": mean, "joint_covariance": covariance, "joint_covariance_minus_reference": delta,
        "wall_target_cross_covariance": cross, "wall_target_cross_covariance_minus_reference": cross_delta,
        "posterior_mean_wall_chi": mean_chi}
    return score, arrays


def extra_scores(replay, output, journal, metadata, sources):
    # Fixed old source paths: no source discovery, inference output or new labels.
    grid_path = sources["paths"]["FINITE_PRIOR_GRID"]
    bank_path = sources["paths"]["FINITE_FULL_RESPONSE_CACHE"]
    obs_path = sources["paths"]["FINITE_24_DEVELOPMENT_CASES"]
    truth_path = sources["paths"]["FINITE_OFFLINE_TRUTHS_NOISE"]
    grid = journal.load_npz(grid_path, GRID_KEYS)
    bank_payload = journal.load_npz(bank_path, ("data", "prior_id"))
    obs_payload = journal.load_npz(obs_path, ("data", "prior_id"))
    for name, payload in (("GRID", grid), ("BANK", bank_payload), ("OBSERVATIONS", obs_payload)):
        validated_prior_id(payload["prior_id"], name)
    bank, observations = bank_payload["data"], obs_payload["data"]
    if (bank.shape != (2950, 4, 16) or observations.shape != (12, 2, 4, 16)
            or np.any(~np.isfinite(bank)) or np.any(~np.isfinite(observations))):
        raise ValueError("finite old BANK/OBSERVATIONS tensors required")
    chi, features = feature_table(grid)
    # This event is flushed before opening truth_index. No noise, clean_means,
    # parent IDs or seeds are loaded, and no additional truth labels are made.
    journal.record("OLD_TRUTH_ACCESS_ALLOWED_AFTER_COMPLETE_SEAL_AND_RAW_INTERVAL_AUDIT", truth_path,
                   ["truth_index"], scoring_only=True, never_online_input=True)
    truth = journal.load_npz(truth_path, ("truth_index",), truth=True)["truth_index"]
    if truth.shape != (12,) or truth.dtype.kind not in "iu" or np.any(truth < 0) or np.any(truth >= 2950):
        raise ValueError("exactly twelve existing parent truth atom indices required")
    counter = {owner: 0 for owner in (*METHOD_IDS, "REFERENCE")}
    rows, reference_rows = [], []
    for parent, repeat in CASES:
        directory = replay / f"parent_{parent:02d}_repeat_{repeat:02d}"
        ref = journal.load_npz(directory / "OFFLINE_REFERENCE.npz", ("joint_weights", "log_joint_weights"))
        _, ref_cov = moments(ref["joint_weights"], features)
        ref_score, ref_arrays = measure_scores(ref["joint_weights"], ref["log_joint_weights"], chi,
            features, int(truth[parent]), ref_cov, bank, observations[parent, repeat, 2:4, :], counter, "REFERENCE")
        reference_rows.append({"parent_index": parent, "repeat_index": repeat, "method_id": "REFERENCE",
                               "status": "SEALED_OFFLINE_REFERENCE", **ref_score})
        shape = {"joint_mean": (10,), "joint_covariance": (10, 10),
            "joint_covariance_minus_reference": (10, 10), "wall_target_cross_covariance": (8, 2),
            "wall_target_cross_covariance_minus_reference": (8, 2), "posterior_mean_wall_chi": (4,)}
        saved = {name: np.full((7, *dimensions), np.nan,
            dtype=complex if name == "posterior_mean_wall_chi" else float) for name, dimensions in shape.items()}
        case_rows = []
        for method_index, method_id in enumerate(METHOD_IDS):
            meta = metadata[(parent, repeat, method_id)]
            payload = journal.load_npz(directory / (method_id + ".npz"), ("joint_weights", "log_joint_weights"))
            weights, log_weights = payload["joint_weights"], payload["log_joint_weights"]
            if not weights.size:
                row = {"parent_index": parent, "repeat_index": repeat, "method_id": method_id,
                    "status": meta["status"], **{key: None for key in ref_score},
                    "logical_offline_response_atom_queries": 0, "actual_new_full_FE": 0, "actual_new_RHS": 0,
                    "outer_log_score_minus_reference": None}
            else:
                score, arrays = measure_scores(weights, log_weights, chi, features, int(truth[parent]),
                    ref_cov, bank, observations[parent, repeat, 2:4, :], counter, method_id)
                for name, value in arrays.items():
                    saved[name][method_index] = value
                row = {"parent_index": parent, "repeat_index": repeat, "method_id": method_id,
                    "status": meta["status"], **score,
                    "outer_log_score_minus_reference": (score["outer_TX2_3_log_posterior_predictive_density"]
                        - ref_score["outer_TX2_3_log_posterior_predictive_density"])}
            rows.append(row)
            case_rows.append(row)
        saved.update({"method_ids": np.asarray(METHOD_IDS), "method_statuses": np.asarray([r["status"] for r in case_rows]),
            "feature_names": np.asarray(FEATURE_NAMES), "existing_parent_truth_atom_index": truth[parent],
            "GT_wall_chi_four_cells": chi[truth[parent]], "GT_position": features[truth[parent], 8:]})
        saved.update({"reference_" + name: value for name, value in ref_arrays.items()})
        np.savez_compressed(output / f"extra_parent_{parent:02d}_repeat_{repeat:02d}.npz", **saved)
        write_json(output / f"extra_parent_{parent:02d}_repeat_{repeat:02d}.json", {
            "parent_index": parent, "repeat_index": repeat,
            "truth_source": "existing twelve-parent truth_index; same parent truth for its two old noise repeats",
            "wall_map_scope": "four equal-area wall ROI cells only; known target intrinsic parameters are not scored as a full material image",
            "feature_order": FEATURE_NAMES, "methods": case_rows,
            "reference": reference_rows[-1], "outer_observation_status": "OLD_EXPOSED_UNASSIMILATED_TX2_3",
            "new_likelihood_factor_added": False, "new_acquisition": False, "blind_test": False})
        check_memory()
    write_csv(output / "METHOD_EXTRA_SCORES.csv", rows)
    write_csv(output / "REFERENCE_EXTRA_SCORES.csv", reference_rows)
    numeric_names = ("GT_position_squared_error", "equal_area_four_cell_wall_chi_MSE",
        "joint_cov_difference_Frobenius", "wall_target_cross_cov_difference_Frobenius",
        "outer_TX2_3_log_posterior_predictive_density", "outer_log_score_minus_reference")
    ranges = {}
    for method_id in METHOD_IDS:
        selected = [row for row in rows if row["method_id"] == method_id]
        ranges[method_id] = {"case_count": len(selected), "scored_case_count": sum(
            row["GT_position_squared_error"] is not None for row in selected)}
        for name in numeric_names:
            values = [row[name] for row in selected if row[name] is not None]
            ranges[method_id][name] = None if not values else {
                "minimum": min(values), "maximum": max(values), "mean": float(np.mean(values))}
    summary = {"label": LABEL, "status": "DESCRIPTIVE_OFFLINE_SCORES_NO_NEW_GATE",
        "case_method_rows": len(rows), "reference_case_rows": len(reference_rows), "ranges_by_method": ranges,
        "covariance_coordinate_policy": "raw physical features: eight dimensionless chi coordinates followed by two lambda0 position coordinates; no fitted scaling or dimensionless gate",
        "logical_offline_saved_response_atom_queries_by_owner": counter,
        "logical_offline_saved_response_atom_queries_total": sum(counter.values()),
        "disk_bank_array_load_count": 1, "disk_bank_response_atom_rows_loaded": 2950,
        "old_bank_sources_and_explicit_relocation": sources, "validated_prior_id": EXPECTED_PRIOR_ID,
        "query_accounting": "independent outer-density evaluation charges 2950 cached response rows per scored method/case and reference/case; disk load is separate, do not add it twice",
        "outer_predictive_real_dimension": OUTER_REAL_DIMENSION,
        "outer_predictive_log_Gaussian_constant": float(-0.5 * OUTER_REAL_DIMENSION * np.log(2 * np.pi * SIGMA ** 2)),
        "outer_scope": "old exposed TX2/3 observations, unassimilated in the sealed base-TX replay; prediction only",
        "actual_new_full_FE": 0, "actual_new_RHS": 0, "online_results_modified": False,
        "new_truth_labels_created": False, "new_gate_created": False}
    write_json(output / "EXTRA_SCORE_SUMMARY.json", summary)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", type=Path, default=ROOT / "results/g1_replay_001")
    parser.add_argument("--output", type=Path, default=ROOT / "results/g1_offline_extra")
    parser.add_argument("--old-bank-root", type=Path, default=None,
        help="explicit relocated OLD finite_bank/online directory; truth is derived from its sibling offline directory")
    args = parser.parse_args(argv)
    replay, output = safe_source(args.replay), args.output.resolve()
    results = (ROOT / "results").resolve()
    if not replay.is_relative_to(results) or not output.is_relative_to(results) or output.is_relative_to(replay):
        raise ValueError("use this campaign's sealed replay and a separate fresh result directory")
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "PLAN.json", {"label": LABEL, "scope": "read-only raw-bound audit and extra offline scores",
        "replay": replay, "output": output, "created_utc": datetime.now(timezone.utc).isoformat(),
        "online_inputs_modified": False, "seal_required_before_any_NPZ_and_truth": True,
        "audit_first_before_truth": True, "bound_violation_policy": "raw signed values and all positives, no clip, tolerance or PASS rescue",
        "GRID_keys": GRID_KEYS, "BANK_and_OBSERVATIONS_keys": ["data", "prior_id"], "truth_keys": ["truth_index"],
        "expected_prior_id": EXPECTED_PRIOR_ID,
        "old_bank_root_explicit_relocation": args.old_bank_root,
        "default_old_source_policy": "CACHE_MAP original paths only; no automatic public/local search",
        "truth_path_policy": "selected old online root's parent/offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz; open only after seal and audit",
        "wall_chi_definition": "support_cell*(epsinf-1+i*ell)",
        "wall_material_MSE": "mean over four equal-area cells of abs(posterior mean chi - existing GT chi)^2",
        "joint_cov_feature_order": FEATURE_NAMES, "cross_cov_shape": [8, 2],
        "joint_cov_feature_units": ["dimensionless"] * 8 + ["lambda0"] * 2,
        "moment_weight_policy": "as-saved float joint weights; no renormalization or repair",
        "outer_predictive": {"TX": [2, 3], "RX": list(range(16)), "sigma_each_real_imag": SIGMA,
            "real_dimension": OUTER_REAL_DIMENSION, "all_Gaussian_constants_included": True,
            "status": "OLD_EXPOSED_UNASSIMILATED_OBSERVATIONS", "new_factor_or_acquisition": False},
        "actual_new_full_FE": 0, "actual_new_RHS": 0, "new_gate_created": False,
        "memory_limit_bytes": RSS_LIMIT_BYTES, "floating_point_roundoff_certified": False})
    journal = AccessJournal(output)
    start, wall = time.process_time(), time.monotonic()
    try:
        cache_map = journal.read_json(ROOT / "research/delegated/source_audit/CACHE_MAP.json")
        sources = resolve_old_sources(cache_map, args.old_bank_root)
        journal.old_sources = sources
        journal.record("OLD_BANK_EXPLICIT_RELOCATION_CONFIGURATION", explicit=sources["explicit_relocation"],
                       hash_check_performed=False)
        metadata, original_receipt = verify_sealed_replay(replay, journal)
        audit = audit_intervals(replay, output, journal, metadata)
        scores = extra_scores(replay, output, journal, metadata, sources)
        receipt = {"status": "OFFLINE_EXTRA_COMPLETED_REQUIRES_PARENT_REVIEW", "label": LABEL,
            "process_cpu_seconds": time.process_time() - start, "wall_seconds": time.monotonic() - wall,
            "peak_RSS_bytes": check_memory(), "actual_new_full_FE": 0, "actual_new_RHS": 0,
            "raw_bound_status": audit["status"], "new_scientific_gate_created": False,
            "online_source_files_modified": False, "existing_truth_only_after_seal": True,
            "old_bank_sources_and_explicit_relocation": sources,
            "logical_offline_saved_response_atom_queries_total": scores["logical_offline_saved_response_atom_queries_total"],
            "historical_bank_paid_setup_unchanged": original_receipt.get("historical_bank_paid_setup"),
            "shared_G0_qualification_paid_setup_unchanged": original_receipt.get("shared_G0_qualification_paid_setup")}
        write_json(output / "RUN_RECEIPT.json", receipt)
        print(json.dumps(jsonable({"status": receipt["status"], "output": output,
            "raw_interval_summary": audit, "ranges_by_method": scores["ranges_by_method"]})), flush=True)
        return 0
    except Exception as error:
        write_json(output / "FAILURE.json", {"status": "FAILED", "error_type": type(error).__name__,
            "error": str(error), "actual_new_full_FE": 0, "actual_new_RHS": 0,
            "process_cpu_seconds": time.process_time() - start, "peak_RSS_bytes": peak_rss(),
            "online_source_files_modified": False, "scientific_gate_decision": "PARENT_REVIEW_REQUIRED"})
        raise


if __name__ == "__main__":
    raise SystemExit(main())
