"""Frozen online trajectory risk actions. No truth archives are opened.

The caller accounts CPU under Book L. This module starts no physical work.
Manifest trajectory_jobs entries are {job: C_job_name, parents: [...], arms:
[...]}. Curvature covers the exact 256 cohort; five other arms cover first32.
"""
import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"
import argparse
import json
import re
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
import numpy as np
from scipy.stats import chi2
from p3r2.trajectory_risk import METHODS, apply_frozen, load_frozen
from infer_p3r2_locked import DEFAULT_OUT, frozen_contract, resolve_path, save, write

CORE_ARM = "dynamic_curvature"
OTHER_ARMS = ("frozen_fixed", "dynamic_fixed", "dynamic_bandit", "ordinary_random", "dynamic_mbrl")
ALL_ARMS = (CORE_ARM,)+OTHER_ARMS


def require_complete(path):
    value = json.loads(Path(path).read_text())
    if value.get("status") != "COMPLETE":
        raise ValueError("required COMPLETE closure: "+str(path))
    return value


def trajectory_contract(out, methods):
    manifest, parents, paths, fusion = frozen_contract(out, methods)
    expected = [f"R2_locked_{i:05d}" for i in range(256)]
    if parents != expected:
        raise ValueError("trajectory cohort must be ordered R2_locked_00000..00255")
    if manifest.get("include_matched_point176") is not True:
        raise ValueError("matched extra16 ordinary control must be explicitly frozen")
    weights_path = resolve_path(manifest.get("trajectory_weights_json", "D/trajectory_integration/FROZEN.json"), out)
    record = next((r for r in manifest["frozen_files"] if resolve_path(r["path"], out) == weights_path), None)
    stat = weights_path.stat()
    if record is None or stat.st_size != int(record["size_bytes"]) or stat.st_mtime_ns != int(record["mtime_ns"]) or stat.st_mtime > float(manifest["frozen_unix"]):
        raise ValueError("trajectory weights missing or changed in frozen artifact registry")
    frozen = load_frozen(weights_path)
    if int(frozen["nodes"]) != 9:
        raise ValueError("frozen trajectory action requires anchor plus eight endpoints")
    jobs = manifest.get("trajectory_jobs")
    if not isinstance(jobs, list) or not jobs:
        raise ValueError("frozen trajectory_jobs list required")
    mapping, closures = {}, []
    seen_jobs = set()
    for entry in jobs:
        name = entry.get("job", "")
        if not re.fullmatch(r"[A-Za-z0-9_\-]+", name) or name in seen_jobs:
            raise ValueError("unique safe C job name required")
        seen_jobs.add(name)
        job = (out/"C/jobs"/name).resolve()
        if len(set(entry["parents"])) != len(entry["parents"]) or len(set(entry["arms"])) != len(entry["arms"]):
            raise ValueError("duplicate manifest parent/arm")
        closure = require_complete(job/"CLOSED.json")
        protocol = json.loads((job/"PROTOCOL.json").read_text())
        if protocol.get("parents") != entry["parents"] or protocol.get("arms") != entry["arms"] or protocol.get("steps") != 8:
            raise ValueError("paid C protocol differs from frozen chunk: "+name)
        if protocol.get("split") != "fresh_eval" or protocol.get("GT_used_online") is not False or protocol.get("historical_exposure") is not False:
            raise ValueError("truth-free fresh trajectory protocol required")
        if (job/"PROTOCOL.json").stat().st_mtime < float(manifest["frozen_unix"]):
            raise ValueError("fresh trajectory job predates freeze")
        for parent in entry["parents"]:
            for arm in entry["arms"]:
                key = (parent, arm)
                if key in mapping:
                    raise ValueError("duplicate frozen trajectory parent/arm")
                mapping[key] = job
        closures.append(dict(job=name, closure=closure))
    planned = {(p, CORE_ARM) for p in parents}|{(p, arm) for p in parents[:32] for arm in OTHER_ARMS}
    if set(mapping) != planned:
        raise ValueError("exact256 curvature plus first32 five-arm frozen coverage required")
    return manifest, parents, mapping, frozen, weights_path, closures


def read_material(path):
    with np.load(path, allow_pickle=False) as archive:
        chi = archive["chi"].copy()
    if chi.shape != (576,) or not np.isfinite(chi).all():
        raise ValueError("finite 576-cell online material required")
    return chi


def read_path(out, parent, arm, job, anchor, noise_threshold):
    directory = job/"trajectories"/parent/arm
    trace = directory/"TRACE.jsonl"
    rows = [json.loads(line) for line in trace.read_text().splitlines()] if trace.is_file() else []
    if len(rows) > 8:
        raise ValueError("more than eight planned trajectory decisions")
    materials, diagnostics, failures = [anchor], [], []
    previous = anchor
    previous_endpoint = parent
    for index, row in enumerate(rows):
        if row.get("parent") != parent or row.get("arm") != arm or row.get("step") != index or row.get("GT_used_online") is not False:
            raise ValueError("invalid truth-free ordered online trace")
        # A cache-loaded anchor key is its filename stem, exactly parent.
        if row.get("initial_anchor_endpoint") != parent or row.get("origin_endpoint") != previous_endpoint:
            raise ValueError("trajectory must start from the original online paid anchor and stay consecutive")
        path = resolve_path(row["cache_path"], out)
        if not path.is_relative_to(directory.resolve()) or path.stem != row["endpoint"]:
            raise ValueError("selected endpoint cache outside registered trajectory")
        try:
            with np.load(path, allow_pickle=False) as archive:
                chi, q = archive["chi"].copy(), archive["q"].copy()
            if chi.shape != (576,) or q.shape != (66,) or not np.isfinite(chi).all() or not np.isfinite(q).all():
                raise ValueError("selected endpoint arrays invalid")
            legal = bool(np.all(q[64:] >= 0) and np.all(q[64:] <= 1))
            rms = float(np.linalg.norm(chi-previous)/24)
            radius = float(row["target_RMS"])
            if radius not in (.02, .04):
                raise ValueError("unfrozen material task radius")
            phi, delta_phi = float(row["Phi"]), float(row["delta_Phi_from_anchor"])
            if not np.isfinite([phi, delta_phi, row["signed_progress"]]).all():
                raise ValueError("nonfinite selected endpoint metadata")
            if not legal:
                failures.append(dict(step=index, reason="PRIOR_PROPERTY_SUPPORT_VIOLATION"))
            diagnostics.append(dict(step=index, endpoint=row["endpoint"], Phi=phi,
                delta_Phi_from_anchor=delta_phi, noise_ball=bool(phi <= noise_threshold),
                delta_Phi_le_2=bool(delta_phi <= 2), native_prior_support_legal=legal,
                Gaussian_prior_penalty=.5*float(q[:64]@q[:64]),
                actual_step_material_RMS=rms, task_radius=radius,
                radius_retained=bool(radius-1e-10 <= rms <= 1.02*radius),
                anchor_displacement_RMS=float(np.linalg.norm(chi-anchor)/24),
                signed_progress=float(row["signed_progress"]),
                progress_retained=bool(float(row["signed_progress"]) >= .5*radius-1e-9),
                measured_endpoints=row.get("endpoints", []), cost=row.get("cost", {})))
            materials.append(chi)
            previous = chi
        except (OSError, KeyError, ValueError) as error:
            failures.append(dict(step=index, reason=repr(error)))
            break
        previous_endpoint = row["endpoint"]
    # Never average a truncated or illegal path as if it were complete.
    if failures:
        materials = [anchor]
    return np.stack(materials), dict(parent=parent, arm=arm, source_job=job.name,
        selected_steps_read=len(diagnostics), expected_steps=8, failures=failures,
        path_complete=bool(len(materials) == 9), selected_endpoint_diagnostics=diagnostics,
        noise_ball_threshold=noise_threshold,
        noise_ball_definition="Phi <= chi2.ppf(.95,512)/2; descriptive noise-energy test",
        candidate_scope="selected endpoints and separately retained paid predictor/corrector metadata; no interpolation check")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--methods", type=Path)
    args = parser.parse_args()
    out = args.root.resolve()
    manifest, parents, mapping, frozen, weights_path, closures = trajectory_contract(out, args.methods or out/"METHODS_FROZEN_FOR_LOCKED.json")
    destination = out/"L/trajectory_estimates"
    if (destination/"COMPLETE.json").exists() or any(destination.glob("R2_locked_*.npz")):
        raise FileExistsError("trajectory estimates already exist; no overwrite or implicit replay")
    destination.mkdir(parents=True, exist_ok=True)
    start_wall, start_cpu = time.monotonic(), time.process_time()
    threshold = float(chi2.ppf(.95, 512)/2)
    statuses = []
    for parent in parents:
        anchor = read_material(out/"L/contexts/cache"/f"{parent}.npz")
        arrays, parent_status = {}, dict(parent=parent, arms={})
        registered = [arm for arm in ALL_ARMS if (parent, arm) in mapping]
        for arm in registered:
            stack, status = read_path(out, parent, arm, mapping[(parent, arm)], anchor, threshold)
            status["methods"] = {}
            for method in METHODS:
                estimate, info = apply_frozen(stack, frozen, arm, method=method)
                arrays[f"{arm}__{method}"] = estimate
                status["methods"][method] = info
            parent_status["arms"][arm] = status
        save(destination/f"{parent}.npz", **arrays, parent=np.array(parent),
            registered_arms=np.asarray(registered), method_names=np.asarray(METHODS),
            status_json=np.array(json.dumps(parent_status, allow_nan=False)),
            source_revision=np.array(manifest["source_revision"]), frozen_unix=np.array(manifest["frozen_unix"]))
        write(destination/f"{parent}.json", parent_status)
        statuses.append(parent_status)
    write(destination/"COMPLETE.json", dict(status="COMPLETE", parents=parents, count=256,
        core_arm=CORE_ARM, secondary_parents=parents[:32], secondary_arms=OTHER_ARMS,
        methods=METHODS, trajectory_jobs=manifest["trajectory_jobs"], source_job_closures=closures,
        source_revision=manifest["source_revision"], frozen_unix=manifest["frozen_unix"],
        weights_path=str(weights_path), truth_read=False, physics_calls=0,
        cpu_seconds=time.process_time()-start_cpu, wall_seconds=time.monotonic()-start_wall,
        parent_status=statuses, posterior_integration=False,
        coefficient_semantics="convex material estimator coefficients, never posterior probabilities"))


if __name__ == "__main__":
    main()
