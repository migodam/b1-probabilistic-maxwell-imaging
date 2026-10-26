"""OFFLINE ONLY: score all 512 predictions after validating their global seal.

Do not run until the parent orchestrator authorizes scoring. This module imports
no solver/training/model actions. It never edits prediction or truth inputs.
The first truth open occurs only after BOTH groups and every sealed prediction
have passed content, schema, parent and anchor checks.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import time

import numpy as np

sys.dont_write_bytecode = True

OUT = Path(__file__).resolve().parents[1]
REPO = OUT.parents[1]
R4 = REPO.parent / "b1-p3l-r4-conditional-transport"
METHODS = ("physics_probit", "strong_cnn", "k32", "deterministic_matched", "distilled_w128")
EXPLORATORY = "det_frozen"
SEEDS = (2026102601, 2026102602, 2026102603)
PHASES = ("outside_roi", "host_plateau", "transition", "inclusion_plateau")


def dump(path, value):
    with Path(path).open("w") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def contained_path(raw_path, root, allowed):
    path = Path(raw_path)
    path = (path if path.is_absolute() else root / path).resolve()
    if not path.is_relative_to(allowed.resolve()):
        raise ValueError(f"Seal artifact escapes prediction directory: {path}")
    return path


def artifact_map(entries, root):
    result = {}
    allowed = root / "raw/predictions"
    for entry in entries:
        path = contained_path(entry["path"], root, allowed)
        if path in result:
            raise ValueError(f"Duplicate sealed artifact {path}")
        if not isinstance(entry["bytes"], int) or entry["bytes"] < 0:
            raise ValueError("Invalid sealed byte count")
        digest = str(entry["blake2b"])
        if len(digest) != 128 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError("Require full lowercase BLAKE2b-512 content digest")
        result[path] = {"bytes": entry["bytes"], "blake2b": digest}
    return result


def manifest_parents(root, group):
    manifest = json.loads((root / f"configs/parents_{group}.json").read_text())
    parents = [str(row["parent"]) for row in manifest["parents"]]
    if len(parents) != 256 or len(set(parents)) != 256:
        raise ValueError(f"Frozen parent manifest {group} must have 256 unique IDs")
    if any(Path(parent).name != parent or parent in (".", "..") for parent in parents):
        raise ValueError("Unsafe parent identifier")
    return parents


def require_sealed_parents(value, parents):
    if isinstance(value, int):
        if value != len(parents):
            raise ValueError("Seal parent count mismatch")
    elif isinstance(value, list):
        if len(value) != len(parents) or set(map(str, value)) != set(parents):
            raise ValueError("Seal parent identifier mismatch")
    else:
        raise ValueError("Seal parents must be an integer count or ID list")


def validate_prediction_seals(root):
    """Validate seals and snapshot every prediction before any truth access."""
    prediction_root = root / "raw/predictions"
    global_path = prediction_root / "ALL_PREDICTIONS_SEALED.json"
    global_seal = json.loads(global_path.read_text())
    if global_seal.get("truth_read") is not False or set(global_seal["groups"]) != {"A", "B"}:
        raise ValueError("Global seal must attest truth_read=false and both groups")
    global_entries = artifact_map(global_seal["artifacts"], root)
    union = {}
    parents_by_group = {}
    methods = list(METHODS)
    group_maps = {}
    for group in ("A", "B"):
        parents = manifest_parents(root, group)
        parents_by_group[group] = parents
        expected_seal_path = (prediction_root / group / "ONLINE_COMPLETE.json").resolve()
        seal_path = contained_path(global_seal["groups"][group], root, prediction_root)
        if seal_path != expected_seal_path:
            raise ValueError("Global group seal path mismatch")
        seal = json.loads(seal_path.read_text())
        if seal.get("truth_read") is not False or seal.get("group") != group:
            raise ValueError(f"Group {group} seal missing group/truth-read attestation")
        require_sealed_parents(seal["parents"], parents)
        entries = artifact_map(seal["artifacts"], root)
        required = {(prediction_root / group / f"{method}_s{seed}.npz").resolve()
                    for method in METHODS for seed in (SEEDS[:1] if method == "physics_probit" else SEEDS)}
        optional = {(prediction_root / group / f"{EXPLORATORY}_s{seed}.npz").resolve() for seed in SEEDS}
        if not required.issubset(entries) or set(entries) - required not in (set(), optional):
            raise ValueError(f"Group {group} seal must contain all canonical predictions and zero/all three det_frozen seeds")
        group_maps[group] = entries
        for path, entry in entries.items():
            if path in union:
                raise ValueError("Artifact shared across groups")
            union[path] = entry
    if set(parents_by_group["A"]) & set(parents_by_group["B"]):
        raise ValueError("Parent overlap between A and B")
    if global_entries != union:
        raise ValueError("Global sealed artifacts differ from exact union of A/B seals")
    has_exploratory = [len(group_maps[g]) > 13 for g in ("A", "B")]
    if any(has_exploratory) and not all(has_exploratory):
        raise ValueError("Exploratory deterministic-frozen method must be complete in both groups")
    if all(has_exploratory):
        methods.append(EXPLORATORY)
    # Read/hash/load the same byte snapshot; later changes cannot alter scoring.
    snapshots = {}
    for path, entry in sorted(union.items()):
        content = path.read_bytes()
        if len(content) != entry["bytes"] or hashlib.blake2b(content).hexdigest() != entry["blake2b"]:
            raise ValueError(f"Sealed prediction content changed: {path}")
        with np.load(io.BytesIO(content), allow_pickle=False) as z:
            required = {"parent_ids", "prediction_raw", "prediction_clipped", "fallback_anchor"}
            if not required.issubset(z.files):
                raise ValueError(f"Missing prediction fields: {path}")
            data = {key: np.asarray(z[key]) for key in required}
        group = path.parent.name
        ids = data["parent_ids"].astype(str)
        if ids.shape != (256,) or len(set(ids)) != 256 or set(ids) != set(parents_by_group[group]):
            raise ValueError(f"Prediction parents differ from frozen manifest: {path}")
        index = {parent: j for j, parent in enumerate(ids)}
        order = [index[parent] for parent in parents_by_group[group]]
        for key in ("prediction_raw", "prediction_clipped"):
            a = data[key]
            if a.shape != (256, 576) or not np.isfinite(a).all():
                raise ValueError(f"Prediction must be finite [256,576]: {path}:{key}")
            data[key] = a[order].astype(np.complex128)
        fallback = data["fallback_anchor"]
        if fallback.shape != (256,) or not np.all((fallback == 0) | (fallback == 1)):
            raise ValueError(f"Fallback vector must be boolean [256]: {path}")
        data["fallback_anchor"] = fallback[order].astype(bool)
        snapshots[path] = data
    return parents_by_group, methods, snapshots, global_entries


def load_static_prior():
    # Load only original prior code/config/points, never p3r4.core or a Target.
    spec = importlib.util.spec_from_file_location("confirmation_static_material_prior", R4 / "src/p1_material/prior.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    config = json.loads((R4 / "research/p3l_r4_conditional_transport/configs/TARGET.json").read_text())
    points = json.loads((R4 / "research/p1_material_target_20261018/inputs/geometry.json").read_text())["points"]
    if float(config["field_mean"]) != 0:
        raise ValueError("Frozen anchor phase field assumes original field_mean=0")
    return module.MaterialPrior(config, points)


def load_anchor_phases(root, parents, prior):
    q0 = np.empty((len(parents), 66))
    chi0 = np.empty((len(parents), 576), complex)
    for j, parent in enumerate(parents):
        with np.load(root / f"raw/contexts/{parent}.npz", allow_pickle=False) as z:
            q, chi = np.asarray(z["q0"]), np.asarray(z["chi0"])
        if q.shape != (66,) or chi.shape != (576,) or not np.isfinite(q).all() or not np.isfinite(chi).all():
            raise ValueError(f"Invalid anchor for {parent}")
        q0[j], chi0[j] = q, chi
    field = q0[:, :64] @ prior.A.T / prior.temperature
    roi = np.asarray(prior.mask, bool)
    phase = np.zeros((len(parents), 576), np.int8)
    phase[:, roi] = np.where(field[:, roi] <= -1, 1, np.where(field[:, roi] >= 1, 3, 2))
    counts = np.stack([(phase == p).sum(axis=1) for p in range(4)], axis=1)
    return q0, chi0, field, roi, phase, counts


def metrics(prediction, truth, roi, phase):
    """All arrays [N,S,M,cell]; statistics never pool seed samples as parents."""
    absolute = np.abs(prediction - truth[:, None, None, :])
    squared = absolute ** 2
    sums = np.stack([(squared * (phase == p)[:, None, None, :]).sum(axis=-1) for p in range(4)], axis=-1)
    return {
        "mse_full": squared.mean(axis=-1), "mae_full": absolute.mean(axis=-1),
        "mse_roi": squared[..., roi].mean(axis=-1), "mae_roi": absolute[..., roi].mean(axis=-1),
        "phase_squared_error_sums": sums,
    }


def summarize_metrics(metric, counts, method_names):
    result = {}
    for j, method in enumerate(method_names):
        row = {name: float(value[:, :, j].mean()) for name, value in metric.items() if name != "phase_squared_error_sums"}
        sums = metric["phase_squared_error_sums"][:, :, j].mean(axis=1)
        row["anchor_phase"] = {}
        for p, phase_name in enumerate(PHASES):
            total_cells = int(counts[:, p].sum())
            row["anchor_phase"][phase_name] = {
                "cell_count": total_cells,
                "cell_weighted_mse": float(sums[:, p].sum() / total_cells) if total_cells else None,
                "full_grid_mse_contribution": float(sums[:, p].mean() / 576),
            }
        result[method] = row
    return result


def score(root=OUT, output=None):
    root = Path(root).resolve()
    output = root / "raw/offline_scores" if output is None else Path(output).resolve()
    if output.exists():
        raise FileExistsError("Offline scoring output already exists; no implicit overwrite")
    parents_by_group, method_names, snapshots, sealed_artifacts = validate_prediction_seals(root)
    parents = parents_by_group["A"] + parents_by_group["B"]
    batches = np.array(["A"] * 256 + ["B"] * 256)
    n, m = 512, len(method_names)
    raw = np.empty((n, 3, m, 576), complex)
    clipped = np.empty_like(raw)
    fallback = np.empty((n, 3, m), bool)
    for gi, group in enumerate(("A", "B")):
        section = slice(gi * 256, (gi + 1) * 256)
        for mi, method in enumerate(method_names):
            for si, seed in enumerate(SEEDS):
                source_seed = SEEDS[0] if method == "physics_probit" else seed
                path = (root / f"raw/predictions/{group}/{method}_s{source_seed}.npz").resolve()
                data = snapshots[path]
                raw[section, si, mi] = data["prediction_raw"]
                clipped[section, si, mi] = data["prediction_clipped"]
                fallback[section, si, mi] = data["fallback_anchor"]
    prior = load_static_prior()
    q0, chi0, anchor_field, roi, phase, counts = load_anchor_phases(root, parents, prior)
    if not np.array_equal(raw[:, 0, 0], chi0):
        raise ValueError("Sealed physics baseline raw prediction differs from frozen anchor chi0")
    # All 512 parents/predictions/seals/contexts validated before this point.
    gate = {"status": "ALL_512_PREDICTIONS_SEALED_AND_VALIDATED_BEFORE_GT",
            "truth_read": False, "parents": n, "methods": method_names,
            "artifacts_validated": len(sealed_artifacts), "content_validation": "BLAKE2b-512 and exact bytes",
            "utc_unix": time.time()}
    output.mkdir(parents=True)
    dump(output / "PRE_GT_VALIDATION.json", gate)
    truth = np.empty((n, 576), complex)
    truth_q = np.empty((n, 66))
    access = {"status": "TRUTH_READ_IN_PROGRESS", "truth_open_attempted": True,
              "truth_parents_materialized": 0, "all_prediction_seals_validated": True}
    dump(output / "OFFLINE_ACCESS_STATUS.json", access)
    try:
        for j, parent in enumerate(parents):
            with np.load(root / f"raw/offline_only/{parent}.npz", allow_pickle=False) as z:
                chi, q = np.asarray(z["chi"]), np.asarray(z["q"])
            access["truth_parents_materialized"] = j + 1
            if chi.shape != (576,) or q.shape != (66,) or not np.isfinite(chi).all() or not np.isfinite(q).all():
                raise ValueError(f"Invalid offline truth for {parent}; no parent deletion allowed")
            truth[j], truth_q[j] = chi, q
        access["status"] = "ALL_TRUTH_MATERIALIZED_SCORING_IN_PROGRESS"
        dump(output / "OFFLINE_ACCESS_STATUS.json", access)
    except BaseException as error:
        access.update(status="FAILED_DURING_OFFLINE_TRUTH_ACCESS", error=repr(error))
        dump(output / "OFFLINE_ACCESS_STATUS.json", access)
        raise
    raw_metrics, clipped_metrics = metrics(raw, truth, roi, phase), metrics(clipped, truth, roi, phase)
    # Average frozen predictions before computing ensemble loss, separately.
    ensemble_raw, ensemble_clipped = raw.mean(axis=1, keepdims=True), clipped.mean(axis=1, keepdims=True)
    ensemble_raw_metrics = metrics(ensemble_raw, truth, roi, phase)
    ensemble_clipped_metrics = metrics(ensemble_clipped, truth, roi, phase)
    np.savez_compressed(output / "PARENT_LOSSES.npz", parent_ids=np.array(parents), batch_ids=batches,
        seed_ids=np.array(SEEDS), methods=np.array(METHODS), split=np.array(["TEST"] * n),
        phase=np.array("independent_confirmation"), raw_losses=raw_metrics["mse_full"][:, :, :5],
        clipped_losses=clipped_metrics["mse_full"][:, :, :5],
        ensemble_raw_losses=ensemble_raw_metrics["mse_full"][:, 0, :5],
        ensemble_clipped_losses=ensemble_clipped_metrics["mse_full"][:, 0, :5])
    compact = {"parent_ids": np.array(parents), "batch_ids": batches, "methods": np.array(method_names),
        "seed_ids": np.array(SEEDS), "phase_names": np.array(PHASES), "phase_cell_counts": counts,
        "phase_roi_proportions": counts[:, 1:] / roi.sum(), "fallback_anchor": fallback}
    for label, metric in (("raw", raw_metrics), ("clipped", clipped_metrics),
                          ("ensemble_raw", ensemble_raw_metrics), ("ensemble_clipped", ensemble_clipped_metrics)):
        compact.update({f"{label}_{name}": value for name, value in metric.items()})
        sums = metric["phase_squared_error_sums"]
        denominator = counts[:, None, None, :]
        # Empty phases have zero storage value plus an explicit validity mask;
        # they must never be interpreted as observed zero conditional error.
        compact[f"{label}_phase_cell_mse"] = np.divide(sums, denominator,
            out=np.zeros_like(sums), where=denominator > 0)
        compact[f"{label}_phase_full_grid_mse_contribution"] = sums / 576
    compact["phase_cell_mse_valid"] = counts > 0
    np.savez_compressed(output / "PARENT_METRICS.npz", **compact)
    # Explicit offline-only full truth/prediction pack for later plotting.
    np.savez_compressed(output / "OFFLINE_SCOREPACK.npz", parent_ids=np.array(parents), batch_ids=batches,
        methods=np.array(method_names), seed_ids=np.array(SEEDS), truth_chi=truth, truth_q=truth_q,
        anchor_chi=chi0, anchor_q=q0, anchor_field_over_halfwidth=anchor_field,
        roi_mask=roi, anchor_phase_id=phase, prediction_raw=raw, prediction_clipped=clipped,
        ensemble_prediction_raw=ensemble_raw[:, 0], ensemble_prediction_clipped=ensemble_clipped[:, 0],
        fallback_anchor=fallback)
    with (output / "PARENT_METRICS.csv").open("w", newline="") as stream:
        names = ["parent", "batch", "method", "evaluation", "seed", "fallback_anchor"]
        names += [f"{kind}_{name}" for kind in ("raw", "clipped") for name in ("mse_full", "mae_full", "mse_roi", "mae_roi")]
        names += [f"anchor_{p}_cells" for p in PHASES]
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader()
        for pi, parent in enumerate(parents):
            for mi, method in enumerate(method_names):
                for evaluation, si in [("seed", i) for i in range(3)] + [("seed_mean", None), ("ensemble", None)]:
                    row = {"parent": parent, "batch": str(batches[pi]), "method": method,
                        "evaluation": evaluation, "seed": SEEDS[si] if si is not None else "",
                        "fallback_anchor": int(fallback[pi, si, mi]) if si is not None else int(fallback[pi, :, mi].any())}
                    row.update({f"anchor_{p}_cells": int(counts[pi, j]) for j, p in enumerate(PHASES)})
                    for kind, seedmetric, ensmetric in (("raw", raw_metrics, ensemble_raw_metrics), ("clipped", clipped_metrics, ensemble_clipped_metrics)):
                        for name in ("mse_full", "mae_full", "mse_roi", "mae_roi"):
                            value = ensmetric[name][pi, 0, mi] if evaluation == "ensemble" else seedmetric[name][pi, :, mi].mean() if si is None else seedmetric[name][pi, si, mi]
                            row[f"{kind}_{name}"] = float(value)
                    writer.writerow(row)
    summary = {"status": "OFFLINE_SCORING_COMPLETE_NO_STATISTICAL_DECISION", "truth_read": True,
        "parents": n, "batch_counts": {"A": 256, "B": 256}, "methods": method_names,
        "confirmatory_methods": list(METHODS), "exploratory_methods": [x for x in method_names if x not in METHODS],
        "seeds": list(SEEDS), "physics_seed_policy": "one sealed prediction repeated over three slots",
        "roi_cells": int(roi.sum()), "phase_definition": "anchor-only A@q0[:64]/halfwidth, ROI host<=-1 transition(-1,1) inclusion>=1; outside ROI separate",
        "phase_roi_proportions": dict(zip(PHASES[1:], map(float, counts[:, 1:].mean(axis=0) / roi.sum()))),
        "phase_interpretation": "anchor-defined prior saturation transition diagnostic; not Maxwell-specific mechanism evidence",
        "phase_aggregation": "phase MSE is summed squared error / summed cell count; full-grid contribution is mean parent phase squared-error sum /576",
        "raw_seed_mean": summarize_metrics(raw_metrics, counts, method_names),
        "clipped_seed_mean": summarize_metrics(clipped_metrics, counts, method_names),
        "raw_ensemble": summarize_metrics(ensemble_raw_metrics, counts, method_names),
        "clipped_ensemble": summarize_metrics(ensemble_clipped_metrics, counts, method_names),
        "fallback_counts_parent_seed": {method: int(fallback[:, :, j].sum()) for j, method in enumerate(method_names)},
        "no_physical_solves": True, "no_training": True,
        "statistics_next_step": "parent orchestrator runs frozen statistics_confirm.py on PARENT_LOSSES.npz",
        "utc_unix": time.time()}
    dump(output / "SCORING_SUMMARY.json", summary)
    access["status"] = "OFFLINE_SCORING_COMPLETE"
    dump(output / "OFFLINE_ACCESS_STATUS.json", access)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=OUT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = score(args.root, args.output)
    print(json.dumps({"status": "OFFLINE_SCORING_COMPLETE", "output": str(output)}))


if __name__ == "__main__":
    main()
