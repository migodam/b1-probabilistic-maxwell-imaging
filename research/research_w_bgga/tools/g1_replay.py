"""Replay exposed finite24 cases. No PDE import/action or hidden-file access.

Run through the parent-owned metered G1_LOGIC wrapper, because this is
zero-new-physics work with its own authorized CPU allowance. The first file in a fresh result directory is PLAN.json.
All seven online outputs for all 24 cases are sealed before any independent
offline reference is calculated. Only this tool can open the audited bank.
"""

from pathlib import Path
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import json
import os
import resource
import sys
import time

for _thread_variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                         "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_thread_variable] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bgga.replay import (
    FinitePriorGrid, G0_PARAMETER_NAMES, G0_PARAMETER_SCALES, METHODS,
    METHOD_DEFINITIONS, JOINT_PARTITION_LABEL, M4_CHECKPOINTS, M5_FUTURE_DIAMETER,
    M5_POTENTIAL_DIAMETER, MethodOutput, REPLAY_LABEL, SECANT_SCALES,
    SIGMA_REAL_IMAG, SecantQualification, offline_reference,
    score_interval_checkpoint, score_joint,
)


PEAK_RSS_LIMIT_BYTES = 1_000_000_000
G0_FIXTURE_IDS = ("strong_center", "right_large", "low_loss_high_contrast")
GRID_KEYS = ("prior_weight", "pose", "wall_group_id", "support_id", "support", "epsinf", "ell")
EXPECTED_PRIOR_ID = "FINITE_PRIOR_ENUM4_V1"
INPUT_SPECS = {
    "FINITE_PRIOR_GRID": ("FINITE_GRID.npz", (*GRID_KEYS, "prior_id")),
    "FINITE_FULL_RESPONSE_CACHE": ("FINITE_BANK.npz", ("data", "prior_id")),
    "FINITE_24_DEVELOPMENT_CASES": ("OBSERVATIONS.npz", ("data", "prior_id")),
}
OLD_ONLINE_DIR = ROOT.parents[2] / "research-w-slam/results/e2e3/finite_bank/online"


def _jsonable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        # JSON has no infinities. NPZ retains the exact IEEE arrays; these strings
        # preserve the distinction from undefined/missing values (JSON null).
        return "-Infinity" if value == -np.inf else "Infinity" if value == np.inf else "NaN"
    return value


def _write_json(path, data):
    with Path(path).open("x") as stream:
        json.dump(_jsonable(data), stream, indent=2, allow_nan=False)
        stream.write("\n")


def _safe_source(path):
    path = Path(path).resolve()
    lowered = [part.lower() for part in path.parts]
    if any("locked" in part or "posteriors" in part or "inference" in part for part in lowered):
        raise ValueError("locked/inference/posterior source paths are forbidden")
    return path


def audited_inputs(cache_map, old_bank_root=None):
    """Original CACHE_MAP sources, or an explicit relocated old online root.

    No local/public fallback search occurs. A relocation changes only the root;
    exact filenames, field allowlists, shape and prior-ID contracts remain.
    """
    entries = {row["id"]: row for row in cache_map["files"]}
    original_root = _safe_source(Path(cache_map["roots"]["raw_WSLAM"])
                                 / "results/e2e3/finite_bank/online")
    relocated = None if old_bank_root is None else _safe_source(old_bank_root)
    result = {}
    for identifier, (filename, keys) in INPUT_SPECS.items():
        row = entries[identifier]
        original = _safe_source(row["path"])
        if original != (original_root / filename).resolve():
            raise ValueError(f"{identifier}: source differs from the fixed audited allowlist")
        path = original if relocated is None else _safe_source(relocated / filename)
        result[identifier] = {"path": path, "keys": keys, "audit_bytes": row["bytes"],
            "original_audited_path": original,
            "source_mode": "CACHE_MAP_ORIGINAL" if relocated is None else "EXPLICIT_OLD_BANK_RELOCATION"}
    return result


def _load_keys(path, keys):
    path = _safe_source(path)
    with np.load(path, allow_pickle=False) as archive:
        return {key: np.array(archive[key], copy=True) for key in keys}


def validated_prior_id(value, source):
    value = np.asarray(value)
    if value.shape != () or value.dtype.kind not in "US" or str(value.item()) != EXPECTED_PRIOR_ID:
        raise ValueError(f"{source}: expected scalar prior_id={EXPECTED_PRIOR_ID}")
    return EXPECTED_PRIOR_ID


def load_qualification():
    """Read only explicitly named new-G0 J/J_reference arrays, never a posterior.

    Require the three predeclared fixtures. Missing/partial G0 qualification
    leaves M2 NOT_RUN rather than substituting machine precision or an invented
    bank-specific error. No current G1 observation enters this threshold.
    """
    directory = ROOT / "results/g0_derivatives"
    metadata = directory / "RESULTS.json"
    if not metadata.is_file():
        return SecantQualification(source="new G0 derivative RESULTS.json absent")
    rows = json.loads(metadata.read_text())
    by_id = {row["case"]: row for row in rows}
    if any(identifier not in by_id or not (directory / (identifier + ".npz")).is_file()
           for identifier in G0_FIXTURE_IDS):
        return SecantQualification(source="new G0 derivative fixtures incomplete; all three required")
    errors = []
    for identifier in G0_FIXTURE_IDS:
        if tuple(by_id[identifier]["names"]) != G0_PARAMETER_NAMES:
            raise ValueError("G0 derivative parameter order differs from frozen seven-column scale contract")
        arrays = _load_keys(directory / (identifier + ".npz"), ("J", "J_reference"))
        jacobian, fitted = arrays["J"], arrays["J_reference"]
        if (jacobian.shape != (32, 7) or fitted.shape != jacobian.shape
                or np.any(~np.isfinite(jacobian)) or np.any(~np.isfinite(fitted))):
            raise ValueError("G0 qualification requires finite base-TX 32x7 complex Jacobians")
        difference = jacobian - fitted
        white = np.concatenate((difference.real, difference.imag), axis=0) / SIGMA_REAL_IMAG
        errors.append(float(np.linalg.norm(white * np.asarray(G0_PARAMETER_SCALES), 2)))
    return SecantQualification(max(errors), G0_FIXTURE_IDS, tuple(errors),
        source="new G0 saved J-J_reference, all seven frozen parameter scales",
        status="EMPIRICAL_CROSS_FIXTURE_THRESHOLD")


def qualification_cost_snapshot():
    """Keep shared paid G0 qualification visible, separate from replay's zero PDE.

    Wrapper job-completion CPU is authoritative and counted once. Physical
    totals below are paid reservations, not an assertion that every reservation
    completed. This does not charge old WSLAM ledgers or transfer their budget.
    """
    path = ROOT / "PHYSICS_LEDGER.jsonl"
    if not path.is_file():
        return {"status": "NOT_DEFINED_LEDGER_ABSENT", "source": str(path),
                "full_FE_paid_reservations": None, "RHS_paid_reservations": None,
                "wrapper_process_cpu_seconds": None}
    rows = [json.loads(line) for line in path.read_text().splitlines() if line]
    rows = [row for row in rows if row.get("phase") == "G0"]
    reservations = [row for row in rows if row.get("kind") == "reservation"]
    completions = [row for row in rows if row.get("label") == "job_completion"]
    return {"status": "RECORDED_G0_SHARED_QUALIFICATION_SETUP", "source": str(path),
        "full_FE_paid_reservations": sum(row.get("factorizations", 0) for row in reservations),
        "RHS_paid_reservations": sum(row.get("rhs", 0) for row in reservations),
        "wrapper_process_cpu_seconds": sum(row.get("process_cpu_seconds", 0.0) for row in completions),
        "scope": "entire new G0 paid setup snapshot; not a per-M2 cost allocation",
        "reservation_count": len(reservations), "wrapper_completion_count": len(completions),
        "attribution": "reused qualification is paid setup, never free speed evidence"}


def _peak_rss_bytes():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(rss if sys.platform == "darwin" else 1024 * rss)


def _check_memory():
    value = _peak_rss_bytes()
    if value > PEAK_RSS_LIMIT_BYTES:
        raise MemoryError(f"replay process peak RSS {value} exceeds {PEAK_RSS_LIMIT_BYTES}")
    return value


def _oracle_factory(bank):
    """Each call creates a plain atom-index-only closure; no table attributes."""
    def make_oracle():
        def oracle(atom_index):
            if not isinstance(atom_index, (int, np.integer)) or isinstance(atom_index, (bool, np.bool_)):
                raise TypeError("oracle accepts one integer atom index")
            if not 0 <= int(atom_index) < bank.shape[0]:
                raise IndexError("atom index outside finite bank")
            return bank[int(atom_index), :2, :].copy()
        return oracle
    return make_oracle


def _save_method(directory, output, process_cpu_seconds):
    n = output.joint_weights.size if output.joint_weights is not None else 0
    arrays = {"joint_weights": np.empty(0) if n == 0 else output.joint_weights,
              "log_joint_weights": np.empty(0) if n == 0 else output.log_joint_weights,
              "query_indices": output.query_indices,
              "query_log_likelihood": output.query_log_likelihood}
    checkpoint_rows = []
    for ordinal, checkpoint in enumerate(output.checkpoints):
        metadata = {}
        for key, value in checkpoint.items():
            if isinstance(value, np.ndarray):
                arrays[f"checkpoint_{ordinal}_{key}"] = value
            else:
                metadata[key] = value
        metadata["npz_array_prefix"] = f"checkpoint_{ordinal}_"
        checkpoint_rows.append(metadata)
    np.savez_compressed(directory / (output.method_id + ".npz"), **arrays)
    _write_json(directory / (output.method_id + ".json"), {
        "method_id": output.method_id, "status": output.status, "label": output.label,
        "logical_query_count": output.logical_queries,
        "actual_new_full_FE": output.actual_full_fe, "actual_new_RHS": output.actual_rhs,
        "deployment_method_process_cpu_seconds_excluding_IO": process_cpu_seconds,
        "floating_point_roundoff_certified": output.floating_point_roundoff_certified,
        "details": output.details, "checkpoints": checkpoint_rows,
        "undefined_joint_representation": "empty per-method vector; NaN full-table row with explicit status" if n == 0 else None,
    })


def _load_sealed_method(directory, method_id):
    metadata = json.loads((directory / (method_id + ".json")).read_text())
    with np.load(directory / (method_id + ".npz"), allow_pickle=False) as archive:
        q, logs = archive["joint_weights"].copy(), archive["log_joint_weights"].copy()
        indices, likelihood = archive["query_indices"].copy(), archive["query_log_likelihood"].copy()
        checkpoints = []
        for row in metadata["checkpoints"]:
            checkpoint = {key: value for key, value in row.items() if key != "npz_array_prefix"}
            prefix = row["npz_array_prefix"]
            checkpoint.update({name[len(prefix):]: archive[name].copy()
                               for name in archive.files if name.startswith(prefix)})
            # Restore the only scalar -inf values needed by scoring.
            for name in ("log_evidence_lower", "log_evidence_upper", "unknown_log_prior_numerator_upper"):
                if checkpoint.get(name) == "-Infinity":
                    checkpoint[name] = -np.inf
            checkpoints.append(checkpoint)
    return MethodOutput(method_id, metadata["status"], None if not q.size else q,
        None if not logs.size else logs, indices, likelihood,
        details=metadata["details"], checkpoints=tuple(checkpoints))


def _summary_row(parent, repeat, output, score):
    cdf = score.get("coordinate_CDF_sup") or {}
    events = score.get("event_probability_errors") or {}
    return {"parent_index": parent, "repeat_index": repeat, "method_id": output.method_id,
        "status": output.status, "logical_queries": output.logical_queries,
        "actual_new_full_FE": 0, "actual_new_RHS": 0,
        "joint_TV": score.get("joint_TV"), "CDF_x_sup": cdf.get("x"), "CDF_y_sup": cdf.get("y"),
        "event_x>=0_absolute_error": events.get("pose_x>=0", {}).get("absolute_error"),
        "event_y>=2.1_absolute_error": events.get("pose_y>=2.1", {}).get("absolute_error"),
        "important_joint_mode_max_absolute_error": score.get("important_joint_mode_max_absolute_error"),
        "position_Bayes_regret": score.get("position_Bayes_regret"),
        "position_Bayes_regret_limit": score.get("position_Bayes_regret_limit"),
        "ground_truth_loss_status": "NOT_SCORED_NO_GT_ACCESS"}


def run_replay(output_dir, prior, bank, observations):
    """Writer barrier is global across all 24 cases, before the first reference."""
    if not (output_dir / "PLAN.json").is_file():
        raise ValueError("write PLAN.json before online replay")
    if (output_dir / "ALL_ONLINE_OUTPUTS_SEALED.json").exists():
        raise ValueError("online outputs are already sealed; use a fresh result directory")
    make_oracle = _oracle_factory(bank)
    timings, method_queries, case_dirs = [], [], []
    online_start = time.process_time()
    for parent in range(observations.shape[0]):
        for repeat in range(observations.shape[1]):
            directory = output_dir / f"parent_{parent:02d}_repeat_{repeat:02d}"
            directory.mkdir()
            case_dirs.append((parent, repeat, directory))
            y = observations[parent, repeat, :2, :]
            # Full seven-row atom table, including explicit undefined rows.
            joint = np.full((len(METHODS), prior.size), np.nan)
            log_joint = np.full_like(joint, np.nan)
            statuses = []
            for method_index, method in enumerate(METHODS):
                started = time.process_time()
                output = method(prior, y, make_oracle())
                cpu = time.process_time() - started
                _check_memory()
                _save_method(directory, output, cpu)
                statuses.append(output.status)
                if output.joint_weights is not None:
                    joint[method_index] = output.joint_weights
                    log_joint[method_index] = output.log_joint_weights
                timings.append({"parent_index": parent, "repeat_index": repeat,
                    "method_id": output.method_id, "process_cpu_seconds_excluding_IO": cpu})
                method_queries.append({"parent_index": parent, "repeat_index": repeat,
                    "method_id": output.method_id, "logical_queries": output.logical_queries})
            np.savez_compressed(directory / "JOINT_METHOD_OUTPUTS.npz", joint_weights=joint,
                log_joint_weights=log_joint, method_ids=np.array(tuple(METHOD_DEFINITIONS)),
                statuses=np.array(statuses), original_atom_indices=np.arange(prior.size))
            _write_json(directory / "ONLINE_OUTPUTS_SEALED.json", {
                "parent_index": parent, "repeat_index": repeat, "method_count": len(METHODS),
                "reference_computed": False, "seal": "write-once files; reference starts after global barrier"})
            print(json.dumps({"stage": "online_outputs_sealed", "parent_index": parent,
                              "repeat_index": repeat, "label": REPLAY_LABEL}), flush=True)
    online_cpu = time.process_time() - online_start
    _write_json(output_dir / "ALL_ONLINE_OUTPUTS_SEALED.json", {
        "case_count": len(case_dirs), "method_count_per_case": len(METHODS),
        "all_method_outputs_written": True, "reference_computed": False,
        "reference_input_to_online_methods": False,
        "logical_queries": method_queries, "method_timings": timings,
        "actual_new_full_FE": 0, "actual_new_RHS": 0})

    # No online method is called below this point. Scoring reads the sealed
    # method artifacts and independently recomputes a reference from raw data.
    offline_start = time.process_time()
    summary = []
    reference_queries = 0
    any_raw_violation = False
    for parent, repeat, directory in case_dirs:
        y = observations[parent, repeat, :2, :]
        reference = offline_reference(prior, y, make_oracle())
        reference_queries += reference.logical_queries
        np.savez_compressed(directory / "OFFLINE_REFERENCE.npz",
            joint_weights=reference.joint_weights, log_joint_weights=reference.log_joint_weights,
            log_likelihood=reference.log_likelihood)
        _write_json(directory / "OFFLINE_REFERENCE.json", {
            "logZ_exact": reference.log_evidence, "p0min": reference.minimum_prior_mass,
            "p0min_over_Zexact_plus_p0min": reference.default_unknown_envelope_floor,
            "log_p0min_over_Zexact_plus_p0min": reference.default_unknown_envelope_floor_log,
            "diagnostic_label": "DEFAULT_[0,1]_UNKNOWN_ENVELOPE_NON_TIGHTENING_DIAGNOSTIC",
            "interpretation": "lower floor on the default unqueried upper envelope when at least one positive-prior atom remains; not a lower bound on true omitted posterior mass",
            "logical_queries": reference.logical_queries, "actual_new_full_FE": 0, "actual_new_RHS": 0,
            "computed_after_all_online_outputs_sealed": True,
            "ground_truth_access": False, "floating_point_roundoff_certified": False})
        scores = []
        for method_id in METHOD_DEFINITIONS:
            output = _load_sealed_method(directory, method_id)
            score = score_joint(prior, output, reference)
            if output.checkpoints:
                checkpoint_scores = []
                for checkpoint in output.checkpoints:
                    comparison = score_interval_checkpoint(prior, checkpoint, reference)
                    any_raw_violation |= comparison["any_positive_raw_bound_violation"]
                    temporary = MethodOutput("M4", "CHECKPOINT_POINT_ESTIMATE",
                        checkpoint["point_joint_weights"], checkpoint["point_log_joint_weights"],
                        output.query_indices[:checkpoint["query_count"]],
                        output.query_log_likelihood[:checkpoint["query_count"]])
                    comparison["queried_conditional_point_estimate_scores"] = score_joint(prior, temporary, reference)
                    checkpoint_scores.append(comparison)
                score["checkpoint_bound_comparisons"] = checkpoint_scores
            scores.append(score)
            summary.append(_summary_row(parent, repeat, output, score))
        _write_json(directory / "OFFLINE_SCORES.json", scores)
        _check_memory()
    _write_json(output_dir / "SCORE_TABLE.json", summary)
    return {"online_phase_process_cpu_seconds_including_output_IO": online_cpu,
        "offline_phase_process_cpu_seconds_including_output_IO": time.process_time() - offline_start,
        "online_method_logical_queries": sum(row["logical_queries"] for row in method_queries),
        "offline_reference_logical_queries": reference_queries,
        "raw_interval_bound_comparison_status": ("VIOLATION_PRESENT" if any_raw_violation
            else "NO_POSITIVE_RAW_VIOLATION_OBSERVED"),
        "scientific_gate_decision": "PARENT_REVIEW_REQUIRED", "case_method_rows": len(summary)}


def main(argv=None):
    qualification=json.loads((ROOT/'results/G0_DECISION.json').read_text())
    if qualification.get('status')!='QUALIFIED_DISCRETE_FULLWAVE_FOR_TESTED_CONTRACT':
        raise RuntimeError('G0 numerical prerequisite has not passed; no G1 replay launch')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-bank-root", type=Path, default=None,
        help="explicit relocated OLD finite_bank/online directory; no automatic source search")
    parser.add_argument("--output", type=Path,
        default=ROOT / "results" / ("g1_replay_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")))
    args = parser.parse_args(argv)
    output_dir = args.output.resolve()
    if not output_dir.is_relative_to((ROOT / "results").resolve()):
        raise ValueError("new replay directory must be within this campaign's results directory")
    audit_path = ROOT / "research/delegated/source_audit/CACHE_MAP.json"
    audit = json.loads(audit_path.read_text())
    inputs = audited_inputs(audit, args.old_bank_root)
    relocation = {"explicit": args.old_bank_root is not None,
        "requested_online_root": args.old_bank_root,
        "resolved_online_root": inputs["FINITE_PRIOR_GRID"]["path"].parent,
        "policy": "CACHE_MAP original paths by default; exact old filenames under supplied root only on explicit relocation; no source autodiscovery"}
    historical_setup = audit["budget_scope"]["physical_bank_paid_setup"]
    output_dir.mkdir(parents=True, exist_ok=False)
    # This is deliberately the first output file and precedes all numerical data
    # loading, qualification-threshold calculation, and online method execution.
    _write_json(output_dir / "PLAN.json", {
        "label": REPLAY_LABEL, "schema": "BGGA_G1_EXPOSED_FINITE24_REPLAY_V1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "data_status": "EXPOSED_DEVELOPMENT_12_PARENTS_X_2_NOISE_REPEATS",
        "methods": METHOD_DEFINITIONS, "method_order": tuple(METHOD_DEFINITIONS),
        "online_core_inputs": ["prior grid and fixed numerical qualification metadata", "y", "atom-index-only oracle"],
        "reference_barrier": "all seven online outputs for all 24 cases written before any reference calculation",
        "channels": {"TX": [0, 1], "RX": list(range(16)), "sigma_each_real_imag": SIGMA_REAL_IMAG},
        "inputs": inputs, "old_bank_relocation": relocation,
        "observation_keys_read": ["data", "prior_id"], "expected_prior_id": EXPECTED_PRIOR_ID,
        "case_identity": "parent_index/repeat_index are positional indices in OBSERVATIONS.data; no ID/truth columns opened",
        "M2": {"geometry_label": "EMPIRICAL_BANK_SECANT_RETENTION",
            "secant_parameter_scales": SECANT_SCALES,
            "G0_parameter_names": G0_PARAMETER_NAMES, "G0_parameter_scales": G0_PARAMETER_SCALES,
            "full_frozen_scale_dictionary": {"epsinf": 1.0, "ell": 0.1,
                "pose_x": 0.3, "pose_y": 0.2, "contrast": 1.0,
                "target_loss": 0.1, "radius": 0.04, "gammaR": 0.005, "gammaI": 0.005},
            "current_secant_nuisance_columns": ["pose_x", "pose_y"],
            "G0_fixture_ids": G0_FIXTURE_IDS,
            "rank_rule": "s > max_fixture ||packj(J-J_reference) diag(frozen_scales)||2",
            "qualification_source": str(ROOT / "results/g0_derivatives"),
            "qualification_NPZ_keys": ["J", "J_reference"],
            "missing_qualification": "NOT_RUN with NOT_DEFINED geometry; no machine threshold",
            "tie_rule": "highest minimum retention, then higher own wall posterior mass, then stable wall-group ID"},
        "M4": {"checkpoints": M4_CHECKPOINTS, "unknown_likelihood_interval": [0.0, 1.0],
            "query_rule": "prior * interval_width / cost; cost=1; exact atoms never requeried",
            "stop_unknown_mass_upper": 0.02, "log_domain_normalization": True},
        "M5": {"potential_complete_link_diameter": M5_POTENTIAL_DIAMETER,
            "future_complete_link_max_pose_whitened_diameter": M5_FUTURE_DIAMETER,
            "future_scope": "25 prior target pose nodes, fixed known intrinsic parameters and gamma",
            "cluster_order_and_representative": "stable wall-group ID, first-fit complete-link, first legal group representative"},
        "scoring": {"joint_TV": True, "coordinate_CDF_sup": ["x", "y"],
            "events": ["pose_x>=0", "pose_y>=2.1"],
            "joint_modes": ["support_id", "x_sign_nonnegative_includes_zero", "y>=2.1"],
            "joint_mode_interpretation": JOINT_PARTITION_LABEL,
            "x_sign_partition": "negative iff x<0; nonnegative iff x>=0; original finite nodes unchanged; no zero tolerance",
            "important_reference_joint_mode_min_mass": 0.05,
            "Bayes_squared_position_regret_limit": ".05*R_ref+(.005)^2",
            "support_modes": "reported separately", "ground_truth_loss": "NOT_SCORED_NO_GT_ACCESS"},
        "new_scientific_gates_created": False,
        "new_actual_full_FE": 0, "new_actual_RHS": 0,
        "historical_bank_paid_setup": historical_setup,
        "historical_setup_policy": "2950 FE/11800 RHS/about9288 CPU seconds remain paid historical setup, not a free online speed result",
        "memory_limit_bytes": PEAK_RSS_LIMIT_BYTES,
        "floating_point_roundoff_certified": False})
    start_cpu, start_wall = time.process_time(), time.monotonic()
    try:
        _write_json(output_dir / "SOURCE_ACCESS.json", {
            "label": REPLAY_LABEL, "access_record_kind": "declared read intents before NPZ loads",
            "old_bank_relocation": relocation, "exact_online_source_arrays": inputs,
            "prior_identity_metadata_keys": ["prior_id"],
            "locked_inference_posteriors_paths_rejected": True,
            "offline_truth_or_posterior_source_access": False, "hash_check_performed": False})
        qualification = load_qualification()
        cost = qualification_cost_snapshot()
        payloads = {identifier: _load_keys(spec["path"], spec["keys"]) for identifier, spec in inputs.items()}
        for identifier, payload in payloads.items():
            validated_prior_id(payload["prior_id"], identifier)
        grid = {key: payloads["FINITE_PRIOR_GRID"][key] for key in GRID_KEYS}
        bank = payloads["FINITE_FULL_RESPONSE_CACHE"]["data"]
        observations = payloads["FINITE_24_DEVELOPMENT_CASES"]["data"]
        prior = FinitePriorGrid(**grid, qualification=qualification, prior_id=EXPECTED_PRIOR_ID)
        if prior.size != 2950 or len(prior.group_ids) != 118 or len(prior.pose_nodes) != 25:
            raise ValueError("audited finite prior must have 2950 atoms =25 poses x118 wall groups")
        if bank.shape != (2950, 4, 16) or observations.shape != (12, 2, 4, 16):
            raise ValueError("audited BANK/OBSERVATIONS shapes differ from the frozen replay")
        if np.any(~np.isfinite(bank)) or np.any(~np.isfinite(observations)):
            raise ValueError("nonfinite audited responses/observations")
        bank.setflags(write=False)
        observations.setflags(write=False)
        _check_memory()
        _write_json(output_dir / "INPUT_RECEIPT.json", {
            "allowlisted_source_arrays": inputs, "bank_shape": bank.shape,
            "old_bank_relocation": relocation, "validated_prior_id": EXPECTED_PRIOR_ID,
            "observations_shape": observations.shape,
            "atom_order": "lookup validated from (pose,wall_group_id); historical pose-major ordering never assumed by reshape",
            "M2_numerical_qualification": asdict(qualification),
            "shared_G0_qualification_paid_setup": cost,
            "historical_bank_paid_setup": historical_setup,
            "reference_or_GT_source_files_opened": False,
            "development_cases_already_exposed": True})
        result = run_replay(output_dir, prior, bank, observations)
        receipt = {"status": "COMPLETED_REQUIRES_PARENT_SCIENTIFIC_REVIEW", "label": REPLAY_LABEL,
            "actual_new_full_FE": 0, "actual_new_RHS": 0, "new_physics_executed": False,
            "process_cpu_seconds_excluding_initial_PLAN_write": time.process_time() - start_cpu,
            "wall_seconds_excluding_initial_PLAN_write": time.monotonic() - start_wall,
            "peak_RSS_bytes": _check_memory(), "memory_limit_bytes": PEAK_RSS_LIMIT_BYTES,
            "historical_bank_paid_setup": historical_setup,
            "shared_G0_qualification_paid_setup": cost,
            "measured_new_PDE_speed_evidence": False, "ground_truth_access": False,
            "new_scientific_gates_created": False, **result}
        _write_json(output_dir / "RUN_RECEIPT.json", receipt)
        print(json.dumps({"status": receipt["status"], "label": REPLAY_LABEL,
                          "output": str(output_dir), "actual_new_full_FE": 0, "actual_new_RHS": 0}), flush=True)
        return 0
    except Exception as error:
        _write_json(output_dir / "FAILURE.json", {
            "status": "FAILED", "label": REPLAY_LABEL, "error_type": type(error).__name__,
            "error": str(error), "actual_new_full_FE": 0, "actual_new_RHS": 0,
            "elapsed_process_cpu_seconds": time.process_time() - start_cpu,
            "peak_RSS_bytes": _peak_rss_bytes(), "scientific_gate_decision": "PARENT_REVIEW_REQUIRED"})
        raise


if __name__ == "__main__":
    raise SystemExit(main())
