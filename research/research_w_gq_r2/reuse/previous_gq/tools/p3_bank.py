"""Bounded P3 old-bank development candidate; online and scorer are separate.

Online reads row_index/data from individual paid_rows only. It never opens
FINITE_BANK, truth, or a saved reference. The mathematical online core receives
only a proper prior/grid, y, and an atom-index callable. The scorer is a later
entry point and requires every online case/method output to be sealed first.
"""

from pathlib import Path
from dataclasses import dataclass
import argparse
import csv
import json
import os
import sys
import time

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gq.contrast import (EMPIRICAL, STRICT, FiniteProperPrior, MarginalEnvelope,
    group_anchor_envelope, logsumexp, normalized_target_intervals, readonly,
    target_posterior, tv_oscillation_bound)

BGGA_COMMIT = "73ba6eec3746548bf21259f627b635678709df00"
WSLAM_COMMIT = "52351aee6ddd41b09851fe5e52d1a8861891ef13"
PRIOR_ID = "FINITE_PRIOR_ENUM4_V1"
LABEL = "EXPOSED_OLD_BANK_DEVELOPMENT_NOT_NEW_PDE_SPEED_EVIDENCE"
METHOD_IDS = ("CONTRAST_P3", "ORDINARY_ADAPTIVE_QUADRATURE_MATCHED_NODES")
GRID_KEYS = ("prior_weight", "pose", "wall_group_id", "support_id", "support", "epsinf", "ell", "prior_id")
DEFAULT_OLD_ROOT = ROOT.parents[2] / "research-w-slam/results/e2e3/finite_bank/online"
FROZEN_CONTROLS = {"BGGA_commit": BGGA_COMMIT, "M4": "all24 queried2950; original F64 bounds failed",
    "M5": "all24 kept118/118 wall groups; no compression", "rerun": False}


def safe_path(path):
    path = Path(path).resolve()
    if any(any(word in part.lower() for word in ("locked", "inference", "posteriors")) for part in path.parts):
        raise ValueError("locked/inference/posterior paths are forbidden")
    return path


def jsonable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, np.generic):
        return jsonable(value.item())
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


def load_keys(path, keys):
    path = safe_path(path)
    with np.load(path, allow_pickle=False) as archive:
        return {key: np.array(archive[key], copy=True) for key in keys}


def prior_identity(value):
    value = np.asarray(value)
    if value.shape != () or str(value.item()) != PRIOR_ID:
        raise ValueError("old finite prior identity differs from the fixed contract")


def emit(book, label, **metadata):
    if book is not None:
        book.event("P3", label, **metadata)


@dataclass(frozen=True)
class BankPrior:
    probability: FiniteProperPrior
    pose_nodes: np.ndarray
    group_ids: np.ndarray
    atom_by_wall_target: np.ndarray
    support_id: np.ndarray
    support: np.ndarray
    epsinf: np.ndarray
    ell: np.ndarray

    @classmethod
    def from_grid(cls, grid):
        prior_identity(grid["prior_id"])
        n = len(grid["prior_weight"])
        poses, px = np.unique(grid["pose"], axis=0, return_inverse=True)
        groups, wx = np.unique(grid["wall_group_id"], return_inverse=True)
        lookup = np.full((len(groups), len(poses)), -1, dtype=int)
        if lookup.size != n:
            raise ValueError("one atom per legal wall/target pair required")
        for atom, (wall, target) in enumerate(zip(wx, px)):
            if lookup[wall, target] >= 0:
                raise ValueError("duplicate wall/target atom")
            lookup[wall, target] = atom
        if np.any(lookup < 0):
            raise ValueError("incomplete Cartesian finite bank")
        representative = lookup[:, 0]
        for key in ("support_id", "support", "epsinf", "ell"):
            for wall in range(len(groups)):
                if np.any(grid[key][lookup[wall]] != grid[key][lookup[wall, 0]]):
                    raise ValueError("wall metadata varies with target")
        joint = np.asarray(grid["prior_weight"], dtype=float)[lookup]
        probability = FiniteProperPrior.from_joint(joint)
        return cls(probability, readonly(poses), readonly(groups), readonly(lookup),
            readonly(grid["support_id"][representative]), readonly(grid["support"][representative]),
            readonly(grid["epsinf"][representative]), readonly(grid["ell"][representative]))

    @property
    def central_target(self):
        x, y = np.unique(self.pose_nodes[:, 0]), np.unique(self.pose_nodes[:, 1])
        central = (x[len(x) // 2], y[len(y) // 2])
        candidates = np.flatnonzero(np.all(self.pose_nodes == central, axis=1))
        if len(candidates) != 1:
            raise ValueError("fixed central target must be a legal prior node")
        return int(candidates[0])

    @property
    def descriptor(self):
        # Frozen material/occupancy coordinates. No y, response or GT enters.
        return np.column_stack((self.support.astype(float), (self.epsinf - 1.5) / 4.5, self.ell / 0.25))


class OwnQueries:
    def __init__(self, prior, y, oracle, book=None):
        self.prior, self.oracle, self.book = prior, oracle, book
        self.y = np.asarray(y, dtype=complex)
        if self.y.shape != (2, 16) or np.any(~np.isfinite(self.y)):
            raise ValueError("base TX0/1 by16RX observation required")
        self.log_likelihood = np.full(prior.atom_by_wall_target.shape, np.nan)
        self.indices, self.responses, self.trace = [], [], []

    def query(self, wall, target, phase):
        if np.isfinite(self.log_likelihood[wall, target]):
            return  # Reuse this method's paid anchor, never query it twice.
        atom = int(self.prior.atom_by_wall_target[wall, target])
        response = np.asarray(self.oracle(atom), dtype=complex)
        if response.shape != (2, 16) or np.any(~np.isfinite(response)):
            raise ValueError("atom oracle returns only one finite base2x16 response")
        difference = (response - self.y).ravel() / 0.01
        phi = float(0.5 * (difference.real @ difference.real + difference.imag @ difference.imag))
        if not np.isfinite(phi):
            raise ValueError("finite saved response potential required")
        self.log_likelihood[wall, target] = -phi
        self.indices.append(atom)
        self.responses.append(response.copy())
        self.trace.append({"query_ordinal": len(self.indices), "atom_index": atom,
            "wall_group_id": int(self.prior.group_ids[wall]), "target_index": target,
            "phase": phase, "log_likelihood": -phi})
        emit(self.book, "old_bank_atom_query", atom_index=atom, query_phase=phase, logical_atom_queries=1,
             actual_full_FE=0, actual_RHS=0)

    def full_wall(self, wall, phase):
        for target in range(len(self.prior.pose_nodes)):
            self.query(wall, target, phase)


def initial_material_nodes(prior, count):
    """Predeclared prior-weighted farthest-point nodes; empty wall excluded."""
    legal = list(np.flatnonzero(np.any(prior.support, axis=1)))
    if not 0 <= count <= len(legal):
        raise ValueError("initial material count exceeds legal nonempty walls")
    wall_mass = prior.probability.wall_given_target[:, prior.central_target]
    chosen = []
    descriptor = prior.descriptor
    while len(chosen) < count:
        def key(wall):
            distance = 1.0 if not chosen else min(np.linalg.norm(descriptor[wall] - descriptor[other]) for other in chosen)
            score = wall_mass[wall] if not chosen else distance * np.sqrt(wall_mass[wall])
            return -float(score), int(prior.group_ids[wall])
        selected = min((wall for wall in legal if wall not in chosen), key=key)
        chosen.append(selected)
    return tuple(chosen)


def representation(prior, own, selected):
    """Only same-support measured secants inform finite empirical errors."""
    selected = tuple(sorted(selected, key=lambda wall: int(prior.group_ids[wall])))
    if not selected:
        raise ValueError("at least one initial material node required")
    anchor = prior.central_target
    descriptor = prior.descriptor
    base = own.log_likelihood[:, anchor] + np.log(prior.probability.wall_given_target[:, anchor])
    contrasts = own.log_likelihood[list(selected)] - own.log_likelihood[list(selected), anchor, None]
    assignment = np.empty(len(prior.group_ids), dtype=int)
    error = np.full(own.log_likelihood.shape, np.inf)
    reasons = []
    for wall in range(len(prior.group_ids)):
        branch = [other for other in selected if prior.support_id[other] == prior.support_id[wall]]
        candidates = branch if branch else selected
        representative = min(candidates, key=lambda other: (
            float(np.linalg.norm(descriptor[wall] - descriptor[other])), int(prior.group_ids[other])))
        assignment[wall] = representative
        if wall in selected:
            error[wall] = 0.0
            reasons.append("QUERIED_WALL_EXACT_SAVED_POTENTIAL")
            continue
        distance = float(np.linalg.norm(descriptor[wall] - descriptor[representative]))
        slopes = []
        for left_index, left in enumerate(branch):
            for right in branch[left_index + 1:]:
                span = float(np.linalg.norm(descriptor[left] - descriptor[right]))
                if span > 0:
                    dl = own.log_likelihood[left] - own.log_likelihood[left, anchor]
                    dr = own.log_likelihood[right] - own.log_likelihood[right, anchor]
                    slopes.append(np.abs(dl - dr) / span)
        if slopes:
            error[wall] = distance * np.max(slopes, axis=0)
            reasons.append("EMPIRICAL_SECANT_EXTRAPOLATION_NO_UNIFORM_TAIL_GUARANTEE")
        else:
            reasons.append("NO_GUARANTEE_NO_SAME_SUPPORT_SECANT")
        error[wall, anchor] = 0.0  # Contrast is exactly zero at the paid anchor.
    representatives = np.unique(assignment)
    group_error = np.array([np.max(error[assignment == wall], axis=0) for wall in representatives])
    representative_d = np.array([own.log_likelihood[wall] - own.log_likelihood[wall, anchor] for wall in representatives])
    empirical = group_anchor_envelope(base, assignment, representatives, representative_d,
        group_error, coverage_complete=True, label=EMPIRICAL,
        provenance="same-support observed finite secants times material descriptor distance; unobserved tails not guaranteed")
    point = group_anchor_envelope(base, assignment, representatives, representative_d,
        np.zeros_like(group_error), coverage_complete=True, label=EMPIRICAL,
        provenance="anchor-weighted representative point approximation, not an exact likelihood")
    # Ordinary piecewise-constant discrete quadrature on exactly the same walls.
    # It also uses every paid exact anchor; no known value is thrown away there.
    ordinary_ll = own.log_likelihood[assignment].copy()
    ordinary_ll[:, anchor] = own.log_likelihood[:, anchor]
    ordinary_marginal = logsumexp(np.log(prior.probability.wall_given_target) + ordinary_ll, axis=0)
    known = np.isfinite(own.log_likelihood)
    strict_lower = logsumexp(np.log(prior.probability.wall_given_target)
        + np.where(known, own.log_likelihood, -np.inf), axis=0)
    strict_upper = logsumexp(np.log(prior.probability.wall_given_target)
        + np.where(known, own.log_likelihood, 0.0), axis=0)
    strict = MarginalEnvelope(strict_lower, strict_upper, STRICT,
        "uninformed Gaussian [0,1] tails from these P3 queries only; not the old M4 query algorithm")
    q, logq = target_posterior(prior.probability.target_mass, point.log_lower)
    ordinary_q, ordinary_logq = target_posterior(prior.probability.target_mass, ordinary_marginal)
    empirical_tv = tv_oscillation_bound(empirical.log_lower - point.log_lower,
        empirical.log_upper - point.log_lower, label=EMPIRICAL)
    strict_tv = tv_oscillation_bound(strict.log_lower - point.log_lower,
        strict.log_upper - point.log_lower, label=STRICT)
    empirical_probability = normalized_target_intervals(prior.probability.target_mass, empirical)
    strict_probability = normalized_target_intervals(prior.probability.target_mass, strict)
    event_intervals = {}
    for name, mask in (("x>=0", prior.pose_nodes[:, 0] >= 0), ("y>=2.1", prior.pose_nodes[:, 1] >= 2.1)):
        e, s = empirical_probability.event(mask), strict_probability.event(mask)
        event_intervals[name] = {"empirical_lower": e.lower, "empirical_upper": e.upper,
            "empirical_label": e.label, "strict_saved_model_lower": s.lower,
            "strict_saved_model_upper": s.upper, "strict_label": s.label,
            "floating_point_roundoff_certified": False}
    masks = np.eye(len(prior.pose_nodes), dtype=bool)
    empirical_atom = [empirical_probability.event(mask) for mask in masks]
    strict_atom = [strict_probability.event(mask) for mask in masks]
    return {"selected_walls": readonly(selected, int), "assignment": readonly(assignment),
        "wall_empirical_contrast_error": readonly(error), "tail_reason": tuple(reasons),
        "representatives": readonly(representatives), "group_empirical_error": readonly(group_error),
        "log_group_anchor_mass": point.log_group_anchor_mass,
        "contrast_log_marginal": point.log_lower, "contrast_target_weights": q,
        "contrast_log_target_weights": logq, "ordinary_log_marginal": readonly(ordinary_marginal),
        "ordinary_target_weights": ordinary_q, "ordinary_log_target_weights": ordinary_logq,
        "empirical_log_marginal_lower": empirical.log_lower, "empirical_log_marginal_upper": empirical.log_upper,
        "strict_log_marginal_lower": strict.log_lower, "strict_log_marginal_upper": strict.log_upper,
        "empirical_target_lower": readonly([interval.lower for interval in empirical_atom]),
        "empirical_target_upper": readonly([interval.upper for interval in empirical_atom]),
        "strict_target_lower": readonly([interval.lower for interval in strict_atom]),
        "strict_target_upper": readonly([interval.upper for interval in strict_atom]),
        "normalized_target_event_intervals": event_intervals,
        "empirical_TV_upper_conditional_only": empirical_tv.conditional_TV_upper,
        "strict_saved_model_TV_upper_roundoff_uncertified": strict_tv.conditional_TV_upper,
        "no_same_support_secant_walls": sum(reason == "NO_GUARANTEE_NO_SAME_SUPPORT_SECANT" for reason in reasons),
        "empirical_interval_label": EMPIRICAL, "strict_interval_label": STRICT,
        "finite_empirical_error_is_not_a_uniform_bound": True,
        "unqueried_tail_guarantee_status": "NO_UNIFORM_GUARANTEE_FOR_UNQUERIED_WALLS",
        "floating_point_roundoff_certified": False}


def adaptive_wall(prior, own, current):
    """Frozen heuristic, not an optimal-risk or safe-query theorem."""
    anchor = prior.central_target
    log_base = own.log_likelihood[:, anchor] + np.log(prior.probability.wall_given_target[:, anchor])
    selected = set(current["selected_walls"].tolist())
    log_target_prior = np.log(prior.probability.target_mass)
    candidates = []
    for wall in range(len(prior.group_ids)):
        if wall in selected:
            continue
        representative = int(current["assignment"][wall])
        contrast = own.log_likelihood[representative] - own.log_likelihood[representative, anchor]
        distance = np.linalg.norm(prior.descriptor[wall] - prior.descriptor[representative])
        error = current["wall_empirical_contrast_error"][wall]
        uncertainty = 1.0 if np.any(np.isinf(error)) else float(np.max(error))
        priority = float(log_base[wall] + np.max(log_target_prior + contrast)
                         + np.log1p(distance) + np.log1p(uncertainty))
        candidates.append((-priority, int(prior.group_ids[wall]), wall))
    return None if not candidates else min(candidates)[2]


def online_core(prior, y, oracle, *, initial_nodes=12, adaptive_nodes=16,
                refine_nodes=0, full_fallback=False, book=None):
    """Bounded access-only core: no files, full bank, GT or reference input."""
    conditional = prior.probability.wall_given_target
    if not np.allclose(conditional, conditional[:, [0]], rtol=1e-12, atol=0):
        raise ValueError("P3 old-bank candidate is frozen to the independent prior; correlated algebra is in gq.contrast")
    if initial_nodes < 1 or min(adaptive_nodes, refine_nodes) < 0:
        raise ValueError("positive initial count and bounded nonnegative refinements required")
    own = OwnQueries(prior, y, oracle, book)
    for wall in range(len(prior.group_ids)):
        own.query(wall, prior.central_target, "ALL_WALL_ANCHORS")
    selected = list(initial_material_nodes(prior, initial_nodes))
    for wall in selected:
        own.full_wall(wall, "PREDECLARED_INITIAL_MATERIAL_NODES")
    current = representation(prior, own, selected)
    snapshots = [{"query_count": len(own.indices), **current}]
    for ordinal in range(adaptive_nodes + refine_nodes):
        wall = adaptive_wall(prior, own, current)
        if wall is None:
            break
        selected.append(wall)
        own.full_wall(wall, "ADAPTIVE16" if ordinal < adaptive_nodes else "EXPLICIT_EXTRA_REFINEMENT")
        current = representation(prior, own, selected)
        snapshots.append({"query_count": len(own.indices), **current})
    if full_fallback:
        for wall in range(len(prior.group_ids)):
            if wall not in selected:
                selected.append(wall)
                own.full_wall(wall, "EXPLICIT_FULL_FALLBACK")
        current = representation(prior, own, selected)
        snapshots.append({"query_count": len(own.indices), **current})
    return {"label": LABEL, "status": ("EXACT_SAVED_FINITE_FULL_REFINEMENT" if len(selected) == len(prior.group_ids)
        else "BOUNDED_EMPIRICAL_CANDIDATE_NO_SAFE_EARLY_STOP_CLAIM"),
        "final": current, "snapshots": tuple(snapshots), "query_indices": readonly(own.indices, int),
        "queried_response_data": readonly(own.responses), "query_trace": own.trace,
        "queried_log_likelihood_matrix": readonly(own.log_likelihood),
        "actual_new_full_FE": 0, "actual_new_RHS": 0,
        "safely_avoided_new_physics": 0,
        "shared_unique_atom_queries": len(own.indices),
        "standalone_query_attribution_per_matched_method": len(own.indices)}


class RowOracle:
    """Tool-owned row provider; the online core sees only its atom callable."""
    def __init__(self, old_online_root):
        self.directory = safe_path(old_online_root) / "paid_rows"
        if not self.directory.is_dir():
            raise FileNotFoundError("MISSING paid_rows: no automatic full-BANK fallback is permitted")
        self.calls = []

    def __call__(self, atom):
        if not isinstance(atom, (int, np.integer)) or not 0 <= atom < 2950:
            raise ValueError("one original atom index required")
        payload = load_keys(self.directory / f"{int(atom):04d}.npz", ("row_index", "data"))
        if int(payload["row_index"]) != atom or payload["data"].shape != (4, 16):
            raise ValueError("paid row identity/shape differs from frozen old schema")
        self.calls.append(int(atom))
        return payload["data"][:2].copy()


def save_case(directory, result, prior):
    arrays = {key: value for key, value in result.items() if isinstance(value, np.ndarray)}
    arrays.update({"final_" + key: value for key, value in result["final"].items() if isinstance(value, np.ndarray)})
    snapshot_metadata = []
    for index, snapshot in enumerate(result["snapshots"]):
        metadata = {key: value for key, value in snapshot.items() if not isinstance(value, np.ndarray)}
        for key, value in snapshot.items():
            if isinstance(value, np.ndarray):
                arrays[f"snapshot_{index}_{key}"] = value
        snapshot_metadata.append(metadata)
    arrays.update(pose_nodes=prior.pose_nodes, wall_group_ids=prior.group_ids,
                  atom_by_wall_target=prior.atom_by_wall_target)
    np.savez_compressed(directory / "ONLINE.npz", **arrays)
    write_json(directory / "ONLINE.json", {key: value for key, value in result.items()
        if key not in ("final", "snapshots") and not isinstance(value, np.ndarray)} | {
        "snapshots": snapshot_metadata,
        "final": {key: value for key, value in result["final"].items() if not isinstance(value, np.ndarray)}})
    for method, prefix in zip(METHOD_IDS, ("contrast", "ordinary")):
        np.savez_compressed(directory / (method + ".npz"),
            target_weights=result["final"][prefix + "_target_weights"],
            log_target_weights=result["final"][prefix + "_log_target_weights"],
            log_marginal_likelihood=result["final"][prefix + "_log_marginal"])
    write_json(directory / "ONLINE_OUTPUTS_SEALED.json", {"method_ids": METHOD_IDS,
        "shared_unique_atom_queries": result["shared_unique_atom_queries"],
        "reference_or_truth_seen": False, "label": LABEL})


def read_online_prior(old_root):
    grid = load_keys(old_root / "FINITE_GRID.npz", GRID_KEYS)
    prior = BankPrior.from_grid(grid)
    if prior.atom_by_wall_target.shape != (118, 25):
        raise ValueError("old finite contract must be118 walls by25 poses")
    observation = load_keys(old_root / "OBSERVATIONS.npz", ("base", "prior_id"))
    prior_identity(observation["prior_id"])
    if observation["base"].shape != (12, 2, 2, 16):
        raise ValueError("old base-only observations must be12x2x2x16")
    return prior, observation["base"]


def online_entry(args, book=None):
    old_root, output = safe_path(args.old_bank_root), Path(args.output).resolve()
    if not output.is_relative_to(ROOT):
        raise ValueError("write only under the local GQ root")
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "PLAN.json", {"label": LABEL, "source_commits": {"BGGA": BGGA_COMMIT, "WSLAM": WSLAM_COMMIT},
        "old_online_root": old_root, "initial_nodes": args.initial_nodes, "adaptive_nodes": args.adaptive_nodes,
        "refine_nodes": args.refine_nodes, "full_fallback": args.full_fallback,
        "fixed_anchor": "central legal prior target, chosen before y", "sigma_real_imag": 0.01,
        "initial_rule": "nonempty walls, prior-weighted farthest point in [4occupancies,(eps-1.5)/4.5,ell/.25], stableID ties",
        "adaptive_rule": "anchor log mass + max_target(log target prior + representative contrast) + log1p(distance) + log1p(empirical error or1 when undefined)",
        "empirical_error_rule": "max observed same-support pair secant times descriptor distance; missing secant gives infinite tail",
        "ordinary_control": "same paid walls, ordinary prior-mass piecewise-constant likelihood quadrature with every exact anchor retained",
        "no_guarantee_tails": True, "safe_stop_claim": False, "frozen_controls_not_rerun": FROZEN_CONTROLS,
        "online_storage": "paid_rows individual row_index/data only; FINITE_BANK unopened",
        "historical_bank_paid_setup": {"full_FE": 2950, "RHS": 11800, "CPU_seconds": 9288.023039999998},
        "new_full_FE": 0, "new_RHS": 0, "local_only": True})
    prior, observations = read_online_prior(old_root)
    initial = initial_material_nodes(prior, args.initial_nodes)
    write_json(output / "PRIOR_PLAN.json", {"written_before_any_oracle_query": True,
        "fixed_anchor_target_index": prior.central_target,
        "fixed_anchor_pose": prior.pose_nodes[prior.central_target],
        "initial_wall_indices": initial, "initial_wall_group_ids": prior.group_ids[list(initial)],
        "planned_unique_atom_queries_without_extra_refinement": (len(prior.group_ids)
            + min(len(prior.group_ids), args.initial_nodes + args.adaptive_nodes) * (len(prior.pose_nodes) - 1)),
        "adaptive_walls_not_chosen_using_reference_or_GT": True})
    oracle = RowOracle(old_root)
    write_json(output / "SOURCE_ACCESS.json", {"old_online_root": old_root,
        "GRID_keys": GRID_KEYS, "OBSERVATIONS_keys": ["base", "prior_id"],
        "oracle_row_keys": ["row_index", "data"], "FINITE_BANK_read_online": False,
        "row_storage_shape": [4, 16], "oracle_delivered_TX": [0, 1],
        "truth_and_reference_payloads_read_online": False, "hashes": False})
    cases, started = [], time.process_time()
    for parent in range(12):
        for repeat in range(2):
            directory = output / f"parent_{parent:02d}_repeat_{repeat:02d}"
            directory.mkdir()
            before = len(oracle.calls)
            result = online_core(prior, observations[parent, repeat], oracle,
                initial_nodes=args.initial_nodes, adaptive_nodes=args.adaptive_nodes,
                refine_nodes=args.refine_nodes, full_fallback=args.full_fallback, book=book)
            save_case(directory, result, prior)
            cases.append({"parent": parent, "repeat": repeat, "queries": len(oracle.calls) - before,
                          "status": result["status"]})
    write_json(output / "ALL_ONLINE_OUTPUTS_SEALED.json", {"case_count": 24, "method_ids": METHOD_IDS,
        "all_method_outputs_written": True, "reference_or_truth_seen": False, "cases": cases,
        "process_CPU_seconds": time.process_time() - started, "actual_new_full_FE": 0, "actual_new_RHS": 0})


def score_entry(args, book=None):
    online, output, old_root = safe_path(args.online), Path(args.output).resolve(), safe_path(args.old_bank_root)
    seal = json.loads((online / "ALL_ONLINE_OUTPUTS_SEALED.json").read_text())
    if (seal.get("case_count") != 24 or seal.get("all_method_outputs_written") is not True
            or seal.get("reference_or_truth_seen") is not False or tuple(seal.get("method_ids", ())) != METHOD_IDS):
        raise ValueError("all online outputs must be sealed before independent scoring")
    for parent in range(12):
        for repeat in range(2):
            directory = online / f"parent_{parent:02d}_repeat_{repeat:02d}"
            case_seal = json.loads((directory / "ONLINE_OUTPUTS_SEALED.json").read_text())
            if tuple(case_seal.get("method_ids", ())) != METHOD_IDS or case_seal.get("reference_or_truth_seen") is not False:
                raise ValueError("case outputs are not sealed for the two declared methods")
            for filename in ("ONLINE.npz", "ONLINE.json", "ONLINE_OUTPUTS_SEALED.json", *(method + ".npz" for method in METHOD_IDS)):
                if not (directory / filename).is_file():
                    raise ValueError("sealed online artifact is missing")
    if not output.is_relative_to(ROOT) or output.is_relative_to(online):
        raise ValueError("independent scorer uses a fresh separate GQ output directory")
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "SOURCE_ACCESS.json", {"all_online_outputs_verified_before_access": True,
        "online_output_root": online, "old_online_root": old_root,
        "truth_path": old_root.parent / "offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz",
        "truth_keys": ["truth_index"], "BANK_keys": ["data", "prior_id"],
        "truth_is_old_exposed_development_only": True, "online_files_modified": False})
    prior, observations = read_online_prior(old_root)
    bank = load_keys(old_root / "FINITE_BANK.npz", ("data", "prior_id"))
    prior_identity(bank["prior_id"])
    if bank["data"].shape != (2950, 4, 16):
        raise ValueError("complete old reference response shape required")
    truth = load_keys(old_root.parent / "offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz", ("truth_index",))["truth_index"]
    if truth.shape != (12,) or truth.dtype.kind not in "iu" or np.any(truth < 0) or np.any(truth >= 2950):
        raise ValueError("existing twelve parent truth atom indices required")
    atom_pose = np.empty((2950, 2))
    for wall in range(118):
        atom_pose[prior.atom_by_wall_target[wall]] = prior.pose_nodes
    rows = []
    for parent in range(12):
        for repeat in range(2):
            directory = online / f"parent_{parent:02d}_repeat_{repeat:02d}"
            response = bank["data"][prior.atom_by_wall_target, :2, :]
            difference = (response - observations[parent, repeat]) / 0.01
            ll = -0.5 * np.sum(difference.real ** 2 + difference.imag ** 2, axis=(2, 3))
            marginal = logsumexp(np.log(prior.probability.wall_given_target) + ll, axis=0)
            reference, logreference = target_posterior(prior.probability.target_mass, marginal)
            mean_ref = reference @ prior.pose_nodes
            risk_ref = float(np.sum(reference * np.sum((prior.pose_nodes - mean_ref) ** 2, axis=1)))
            emit(book, "independent_reference", logical_atom_queries=2950, actual_full_FE=0, actual_RHS=0)
            saved = load_keys(directory / "ONLINE.npz", ("final_assignment", "final_wall_empirical_contrast_error",
                "final_empirical_target_lower", "final_empirical_target_upper", "final_strict_target_lower", "final_strict_target_upper",
                "final_empirical_log_marginal_lower", "final_empirical_log_marginal_upper"))
            actual_d = ll - ll[:, prior.central_target, None]
            signed_error = np.abs(actual_d - actual_d[saved["final_assignment"]]) - saved["final_wall_empirical_contrast_error"]
            interval_errors = {
                "empirical_target_lower_signed_violation": saved["final_empirical_target_lower"] - reference,
                "empirical_target_upper_signed_violation": reference - saved["final_empirical_target_upper"],
                "strict_target_lower_signed_violation": saved["final_strict_target_lower"] - reference,
                "strict_target_upper_signed_violation": reference - saved["final_strict_target_upper"],
                "empirical_marginal_lower_signed_violation": saved["final_empirical_log_marginal_lower"] - marginal,
                "empirical_marginal_upper_signed_violation": marginal - saved["final_empirical_log_marginal_upper"]}
            np.savez_compressed(output / f"reference_parent_{parent:02d}_repeat_{repeat:02d}.npz",
                target_weights=reference, log_target_weights=logreference, log_likelihood=ll,
                actual_contrast=actual_d, empirical_wall_error_signed_violation=signed_error,
                **interval_errors)
            for method in METHOD_IDS:
                q = load_keys(directory / (method + ".npz"), ("target_weights",))["target_weights"]
                if q.shape != (25,) or np.any(~np.isfinite(q)) or np.any(q < 0):
                    raise ValueError("valid sealed target weight vector required")
                mean = q @ prior.pose_nodes
                row = {"parent": parent, "repeat": repeat, "method": method,
                    "target_TV": float(0.5 * np.sum(np.abs(q - reference))),
                    "position_Bayes_regret": float(np.sum((mean - mean_ref) ** 2)),
                    "position_GT_squared_error": float(np.sum((mean - atom_pose[truth[parent]]) ** 2)),
                    "reference_position_Bayes_risk": risk_ref,
                    "reference_GT_squared_error": float(np.sum((mean_ref - atom_pose[truth[parent]]) ** 2)),
                    "max_empirical_contrast_error_signed_violation": float(np.max(signed_error)),
                    "positive_empirical_contrast_error_entries": int(np.sum(signed_error > 0)),
                    "target_weight_sum_as_saved": float(q.sum()),
                    "actual_new_full_FE": 0, "actual_new_RHS": 0}
                for name, vector in interval_errors.items():
                    row["maximum_" + name] = float(np.max(vector))
                    row["positive_count_" + name] = int(np.sum(vector > 0))
                for axis, name in enumerate(("x", "y")):
                    row["CDF_" + name + "_sup"] = max(abs(float(q[prior.pose_nodes[:, axis] <= v].sum()
                        - reference[prior.pose_nodes[:, axis] <= v].sum())) for v in np.unique(prior.pose_nodes[:, axis]))
                for name, mask in (("x>=0", prior.pose_nodes[:, 0] >= 0), ("y>=2.1", prior.pose_nodes[:, 1] >= 2.1)):
                    row["event_" + name + "_signed_error"] = float(q[mask].sum() - reference[mask].sum())
                rows.append(row)
    with (output / "TARGET_RESULTS.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    write_json(output / "SCORE_RECEIPT.json", {"label": LABEL, "case_method_rows": len(rows),
        "reference_logical_atom_queries": 24 * 2950, "actual_new_full_FE": 0, "actual_new_RHS": 0,
        "online_files_modified": False, "scientific_decision": "PARENT_REVIEW_REQUIRED",
        "frozen_M4_M5_not_rerun": FROZEN_CONTROLS})


def main(argv=None, book=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    online = subparsers.add_parser("online")
    online.add_argument("--old-bank-root", type=Path, default=DEFAULT_OLD_ROOT)
    online.add_argument("--output", type=Path, required=True)
    online.add_argument("--initial-nodes", type=int, default=12)
    online.add_argument("--adaptive-nodes", type=int, default=16)
    online.add_argument("--refine-nodes", type=int, default=0)
    online.add_argument("--full-fallback", action="store_true")
    score = subparsers.add_parser("score")
    score.add_argument("--old-bank-root", type=Path, default=DEFAULT_OLD_ROOT)
    score.add_argument("--online", type=Path, required=True)
    score.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if book is None:
        from gq.accounting import BudgetBook
        book = BudgetBook(ROOT)
    return online_entry(args, book) if args.command == "online" else score_entry(args, book)


if __name__ == "__main__":
    main()
