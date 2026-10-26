"""Frozen parent-paired R5-B loss analysis. No physical solves or truth access.

NPZ: parent_ids[N], batch_ids[N], seed_ids[3], methods[M],
raw_losses[N,3,M], clipped_losses[N,3,M], split[N], phase scalar.
Optional ensemble_raw_losses[N,M], ensemble_clipped_losses[N,M] are descriptive.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROTOCOL = ROOT / "configs" / "STATISTICAL_PROTOCOL.json"
METHODS = ("physics_probit", "strong_cnn", "k32", "deterministic_matched", "distilled_w128")
COMPARISONS = ("deterministic_vs_k32", "distilled_vs_k32", "distilled_vs_strong", "retention_80pct")


def read_protocol(path=DEFAULT_PROTOCOL):
    with open(path) as f:
        protocol = json.load(f)
    # Fail closed if these frozen core choices change accidentally.
    assert protocol["methods"] == list(METHODS)
    assert protocol["bootstrap"]["minimum_replicates"] >= 20000
    assert protocol["multiplicity"]["family_size"] == 4
    assert protocol["multiplicity"]["family_wise_alpha"] == 0.05
    assert protocol["retention"]["threshold"] == 0.8
    return protocol


def load_losses(path, split="TEST", phase="independent_confirmation"):
    """Only load the allowlisted loss/identifier schema; never pickle objects."""
    required = ("parent_ids", "batch_ids", "seed_ids", "methods", "raw_losses", "clipped_losses", "split", "phase")
    with np.load(path, allow_pickle=False) as z:
        missing = set(required) - set(z.files)
        if missing:
            raise ValueError(f"Missing NPZ keys: {sorted(missing)}")
        data = {key: np.asarray(z[key]) for key in required}
        for key in ("ensemble_raw_losses", "ensemble_clipped_losses"):
            if key in z:
                data[key] = np.asarray(z[key])
    if data["phase"].size != 1 or str(data["phase"].item()) != phase:
        raise ValueError("NPZ phase does not match the requested phase")
    n = len(data["parent_ids"])
    for key in ("parent_ids", "batch_ids", "split"):
        if data[key].shape != (n,):
            raise ValueError(f"{key} must be a one-dimensional parent vector")
    methods = data["methods"].astype(str).tolist()
    if len(methods) != len(METHODS) or set(methods) != set(METHODS):
        raise ValueError(f"Exactly the five frozen methods are required: {METHODS}")
    if data["seed_ids"].shape != (3,) or len(set(data["seed_ids"].astype(str))) != 3:
        raise ValueError("Exactly three distinct frozen seed identifiers are required")
    if len(set(data["parent_ids"].astype(str))) != n:
        raise ValueError("Parents must be unique globally, including across batches")
    order = [methods.index(method) for method in METHODS]
    mask = data["split"].astype(str) == split
    if not np.any(mask):
        raise ValueError(f"No parents for requested split={split}")
    for key in ("raw_losses", "clipped_losses"):
        a = data[key]
        if a.shape != (n, 3, len(METHODS)):
            raise ValueError(f"{key} must have shape [N,3,5]")
        if not np.isfinite(a).all() or (a < 0).any():
            raise ValueError(f"{key} must be finite and nonnegative")
        a = a[mask][:, :, order].astype(np.float64)
        if not np.all(a[:, :, 0] == a[:, :1, 0]):
            raise ValueError("Physics baseline must repeat identical losses across seed slots")
        data[key] = a
    for key in ("ensemble_raw_losses", "ensemble_clipped_losses"):
        if key in data:
            a = data[key]
            if a.shape != (n, len(METHODS)) or not np.isfinite(a).all() or (a < 0).any():
                raise ValueError(f"{key} must be finite nonnegative [N,5]")
            data[key] = a[mask][:, order].astype(np.float64)
    for key in ("parent_ids", "batch_ids", "split"):
        data[key] = data[key][mask].astype(str)
    data["methods"] = np.asarray(METHODS)
    data["phase"] = phase
    return data


def bootstrap_means(values, batch_ids, replicates=20000, rng_seed=2026102705, chunk_size=256):
    """Draw parents within batch; all value columns share each resampling draw."""
    if replicates < 20000:
        raise ValueError("Frozen minimum is 20000 paired parent bootstrap replicates")
    values = np.asarray(values, dtype=np.float64)
    batch_ids = np.asarray(batch_ids)
    if values.ndim != 2 or len(values) != len(batch_ids) or not len(values):
        raise ValueError("Expected nonempty values[N,K] and batch_ids[N]")
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite bootstrap inputs")
    rng = np.random.default_rng(rng_seed)
    groups = [np.flatnonzero(batch_ids == batch) for batch in sorted(set(batch_ids.tolist()))]
    draws = np.empty((replicates, values.shape[1]), dtype=np.float64)
    for start in range(0, replicates, chunk_size):
        count = min(chunk_size, replicates - start)
        accum = np.zeros((count, values.shape[1]))
        for group in groups:
            indices = group[rng.integers(0, len(group), size=(count, len(group)))]
            accum += values[indices].sum(axis=1)
        draws[start:start + count] = accum / len(values)
    return draws


def effects(loss):
    """Positive difference always favors the named candidate; last is retention."""
    _, strong, k32, deterministic, distilled = np.moveaxis(loss, -1, 0)
    numerator, denominator = strong - distilled, strong - k32
    return np.stack((k32 - deterministic, k32 - distilled, numerator,
                     numerator - 0.8 * denominator), axis=-1)


def confidence(draws, alpha=0.0125):
    return {
        "ci_two_sided_98_75": np.quantile(draws, [alpha / 2, 1 - alpha / 2], method="linear").tolist(),
        "lower_one_sided_98_75": float(np.quantile(draws, alpha, method="linear")),
    }


def centered_pvalue(point, draws, one_sided=False):
    centered = draws - point
    exceed = centered >= point if one_sided else np.abs(centered) >= abs(point)
    p = float((1 + np.count_nonzero(exceed)) / (len(draws) + 1))
    return {"p_bootstrap_centered": p, "p_bonferroni_4": min(1.0, 4 * p)}


def metric_summary(losses, bootstrap, include_inference=True):
    """losses are already seed-mean per-parent; bootstrap has same method axis."""
    point = losses.mean(axis=0)
    point_effect, draws_effect = effects(point), effects(bootstrap)
    summaries = {}
    for j, name in enumerate(COMPARISONS):
        bounds = confidence(draws_effect[:, j])
        lo, hi = bounds["ci_two_sided_98_75"]
        direction = "positive" if lo > 0 else "negative" if hi < 0 else "unresolved"
        summaries[name] = {
            "estimate": float(point_effect[j]), **bounds,
            **centered_pvalue(point_effect[j], draws_effect[:, j], one_sided=(j == 3)),
            "adjusted_ci_direction": direction,
            "test": "one_sided_greater" if j == 3 else "two_sided",
            "confidence_gate": bool(bounds["lower_one_sided_98_75"] > 0) if j == 3 else bool(lo > 0 or hi < 0),
            "confirmatory_eligible": include_inference,
        }
    numerator = point[1] - point[4]
    denominator = point[1] - point[2]
    numerator_draws = bootstrap[:, 1] - bootstrap[:, 4]
    denominator_draws = bootstrap[:, 1] - bootstrap[:, 2]
    denom_bounds = confidence(denominator_draws)
    stable = bool(denominator > 0 and denom_bounds["ci_two_sided_98_75"][0] > 0)
    ratio = {
        "numerator_estimate": float(numerator), "denominator_estimate": float(denominator),
        "denominator_bounds": denom_bounds, "positive_stable_denominator": stable,
        "nonpositive_denominator_bootstrap_draws": int(np.count_nonzero(denominator_draws <= 0)),
        "zero_denominator_bootstrap_draws": int(np.count_nonzero(denominator_draws == 0)),
        "estimate": None, "bounds": None,
        "status": "NOT_ESTIMABLE_UNSTABLE_OR_NONPOSITIVE_DENOMINATOR",
        "role": "secondary interpretation; paired linear contrast controls retention inference",
    }
    if stable and not np.any(denominator_draws == 0):
        ratio_draws = numerator_draws / denominator_draws
        ratio.update(estimate=float(numerator / denominator), bounds=confidence(ratio_draws), status="ESTIMABLE")
    elif stable:
        ratio["status"] = "NOT_ESTIMABLE_EXACT_ZERO_BOOTSTRAP_DENOMINATOR"
    retention_gate = bool(stable and summaries["retention_80pct"]["confidence_gate"])
    summaries["retention_80pct"]["confidence_gate_with_stable_denominator"] = retention_gate
    return {
        "method_means": dict(zip(METHODS, map(float, point))),
        "method_mean_bounds_descriptive": {m: confidence(bootstrap[:, j]) for j, m in enumerate(METHODS)},
        "comparisons": summaries, "retention_ratio": ratio,
    }


def analyze(data, protocol=None, replicates=20000, rng_seed=2026102705, split="TEST"):
    protocol = read_protocol() if protocol is None else protocol
    batches = data["batch_ids"]
    counts = {batch: int(np.count_nonzero(batches == batch)) for batch in sorted(set(batches.tolist()))}
    confirmatory = bool(data["phase"] == protocol["phase"] and split == protocol["primary_split"])
    if confirmatory and counts != protocol["independent_batches"]:
        raise ValueError(f"Confirmatory TEST requires exactly A=256,B=256; got {counts}")
    # Concatenate metrics before resampling: identical draws across methods/metrics.
    raw = data["raw_losses"].mean(axis=1)
    clipped = data["clipped_losses"].mean(axis=1)
    draws = bootstrap_means(np.concatenate((raw, clipped), axis=1), batches, replicates, rng_seed)
    result = {
        "protocol_id": protocol["protocol_id"], "phase": data["phase"], "split": split,
        "status": "CONFIRMATORY_ANALYSIS" if confirmatory else "DESCRIPTIVE_OUTSIDE_PRIMARY_PHASE_SPLIT",
        "parents": len(batches), "batch_counts": counts, "seed_ids": data["seed_ids"].tolist(),
        "estimator": protocol["primary_estimator"], "bootstrap_replicates": replicates,
        "bootstrap_rng_seed": rng_seed, "primary_scope": protocol["multiplicity"]["primary_scope"],
        "raw_primary": metric_summary(raw, draws[:, :5], include_inference=confirmatory),
        "clipped_secondary": metric_summary(clipped, draws[:, 5:], include_inference=False),
        "batches_descriptive": {},
    }
    for j, batch in enumerate(sorted(counts)):
        mask = batches == batch
        bd = bootstrap_means(np.concatenate((raw[mask], clipped[mask]), axis=1), batches[mask], replicates, rng_seed + j + 1)
        result["batches_descriptive"][batch] = {
            "raw": metric_summary(raw[mask], bd[:, :5], include_inference=False),
            "clipped": metric_summary(clipped[mask], bd[:, 5:], include_inference=False),
        }
    result["batch_consistency"] = {}
    for name in COMPARISONS:
        points = {batch: result["batches_descriptive"][batch]["raw"]["comparisons"][name]["estimate"] for batch in counts}
        signs = np.sign(list(points.values()))
        result["batch_consistency"][name] = {
            "batch_point_estimates": points,
            "all_positive": bool(np.all(signs > 0)),
            "opposite_signs": bool(np.any(signs > 0) and np.any(signs < 0)),
            "interpretation": "descriptive replication; no batch selection or additional confirmatory claim",
        }
    result["seedwise_losses_descriptive"] = {
        metric: {method: data[key][:, :, j].mean(axis=0).tolist() for j, method in enumerate(METHODS)}
        for metric, key in (("raw", "raw_losses"), ("clipped", "clipped_losses"))
    }
    result["ensemble_descriptive"] = {}
    for metric in ("raw", "clipped"):
        key = f"ensemble_{metric}_losses"
        if key in data:
            result["ensemble_descriptive"][metric] = dict(zip(METHODS, map(float, data[key].mean(axis=0))))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("losses", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--split", default="TEST")
    parser.add_argument("--phase", default="independent_confirmation")
    parser.add_argument("--replicates", type=int, default=20000)
    parser.add_argument("--rng-seed", type=int, default=2026102705)
    args = parser.parse_args()
    data = load_losses(args.losses, split=args.split, phase=args.phase)
    result = analyze(data, replicates=args.replicates, rng_seed=args.rng_seed, split=args.split)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write("\n")
    print(json.dumps({"status": result["status"], "parents": result["parents"], "output": str(args.output)}))


if __name__ == "__main__":
    main()
