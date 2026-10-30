"""Public, truth-free initial states and finite multistart task plans.

This module prepares independent solver starts. It never optimizes or scores GT.
The caller must pass each 576-vector as the adapter's initial material.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

BOUNDS = (0.02, 20.0)
START_NAMES = (
    "constant_0p1", "constant_1", "constant_4", "smooth_random_1",
    "correlated_field_1", "block_field_1", "smooth_random_4",
    "correlated_field_4",
)


def _smooth(x: np.ndarray, passes: int) -> np.ndarray:
    x = np.asarray(x, float)
    for _ in range(passes):
        x = (4*x + np.roll(x, 1, 0) + np.roll(x, -1, 0)
             + np.roll(x, 1, 1) + np.roll(x, -1, 1)) / 8
    return x


def _field(rng: np.random.Generator, mean: float, passes: int,
           spread: float) -> np.ndarray:
    z = _smooth(rng.normal(size=(24, 24)), passes)
    z -= z.mean()
    z /= max(float(z.std()), 1e-12)
    return np.clip(mean + spread*z, *BOUNDS)


def make_starts(seed: int = 20261030) -> np.ndarray:
    """Eight predeclared states, built without labels, observations, or truth."""
    rng = np.random.default_rng(seed)
    block = np.repeat(np.repeat(rng.uniform(0.15, 2.0, (6, 6)), 4, 0), 4, 1)
    starts = [
        np.full((24, 24), .1), np.full((24, 24), 1.),
        np.full((24, 24), 4.), _field(rng, 1., 4, .35),
        _field(rng, 1., 20, .35), _smooth(block, 2),
        _field(rng, 4., 4, 1.2), _field(rng, 4., 20, 1.2),
    ]
    out = np.stack([np.clip(a, *BOUNDS).ravel() for a in starts])
    assert out.shape == (8, 576) and np.isfinite(out).all()
    return out


def split_integer(total: int, count: int) -> list[int]:
    q, r = divmod(total, count)
    return [q + int(i < r) for i in range(count)]


def plan(cases: list[dict], seed: int, wall: float, lu: int,
         rhs: int) -> dict:
    """Plans K=1,4,8 in both separately labeled budget regimes.

    `per_start`: each start gets B, hence total allowance K*B.
    `fixed_total`: the K starts share B, including wall, LU and RHS.
    Plans are descriptive; measured costs must come from execution receipts.
    """
    if wall <= 0 or lu < 8 or rhs < 96:
        raise ValueError("Positive wall and enough counts for eight starts required")
    tasks = []
    for case in cases:
        for k in (1, 4, 8):
            for regime in ("per_start", "fixed_total"):
                walls = [wall]*k if regime == "per_start" else [wall/k]*k
                lus = [lu]*k if regime == "per_start" else split_integer(lu, k)
                rhss = [rhs]*k if regime == "per_start" else split_integer(rhs, k)
                for i in range(k):
                    tasks.append({
                        "task_id": f"{case['case_id']}__K{k}__{regime}__s{i}",
                        "case_id": case["case_id"], "public_path": case["public_path"],
                        "K": k, "regime": regime, "start_index": i,
                        "start_name": START_NAMES[i],
                        "start_npz_key": "starts", "start_row": i,
                        "max_wall_seconds": walls[i], "max_lu": lus[i],
                        "max_complex_rhs": rhss[i],
                        "candidate_selection": "lowest finite public physical ranking_score",
                    })
    return {
        "schema": "fbg-r2-multistart-plan-v1", "status": "DEVELOPMENT",
        "start_seed_public": seed, "start_names": START_NAMES,
        "start_shape": [8, 576], "K_values": [1, 4, 8],
        "per_start_base_budget": {"wall_seconds": wall, "LU": lu, "actual_complex_RHS": rhs},
        "budget_rule": {
            "per_start": "K times base allowance; measured cost not inferred",
            "fixed_total": "one base allowance split across K starts; rounding fixed by task plan",
        },
        "tasks": tasks,
    }


def select_public(receipts: list[dict]) -> dict:
    """Select one top-1 endpoint using solver's public physical score only."""
    eligible = []
    for row in receipts:
        result = row.get("result", row)
        m = result.get("m")
        top = result.get("top1") or {}
        score = top.get("ranking_score")
        if m is None or score is None:
            continue
        material = np.asarray(m, float)
        if material.size == 576 and np.isfinite(material).all() and np.isfinite(score):
            eligible.append((float(score), row.get("task_id", ""), result))
    if not eligible:
        return {"status": "NO_FINITE_ENDPOINT", "m": None, "ranking_score": None}
    score, task_id, result = min(eligible, key=lambda x: (x[0], x[1]))
    return {"status": "FINITE_ENDPOINT", "selected_task_id": task_id,
            "m": result["m"], "ranking_score": score,
            "relative_residual": (result.get("top1") or {}).get("relative_residual"),
            "selection_truth_access": "NONE"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=20261030)
    ap.add_argument("--wall", type=float, default=45.)
    ap.add_argument("--lu", type=int, default=2000)
    ap.add_argument("--rhs", type=int, default=60000)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    starts = make_starts(a.seed)
    np.savez_compressed(a.out / "starts.npz", starts=starts,
                        names=np.asarray(START_NAMES))
    cases = json.loads(a.cases.read_text())["cases"]
    manifest = plan(cases, a.seed, a.wall, a.lu, a.rhs)
    manifest["starts_path"] = "starts.npz"  # relative to multistart_plan.json
    (a.out / "multistart_plan.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"tasks": len(manifest["tasks"]), "starts": 8,
                      "output": str(a.out)}))


if __name__ == "__main__":
    main()
