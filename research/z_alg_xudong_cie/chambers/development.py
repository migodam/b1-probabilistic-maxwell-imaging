"""Audit the sealed Twin development trace at its pre-existing stage endpoints."""

import argparse
import json
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from score import score_material
from training import HERE, elapsed, timers

CASE_ID = "fresnel_twin_8GHz_half180_tx"
TRACE = HERE.parent / "diagnosis" / "trajectories" / CASE_ID


def run(trace, output, terminal_map):
    started = timers()
    receipt = json.loads((trace / "receipt.json").read_text(encoding="utf-8"))
    with np.load(trace / "materials.npz", allow_pickle=False) as src:
        m, step, stage = (src[k] for k in ("materials", "step", "stage"))
    if (receipt.get("truth_access") != "NONE" or not receipt.get("sealed_utc") or
            not np.all(np.diff(step) > 0) or receipt["final_step"] != int(step[-1])):
        raise ValueError("Use the sealed, aligned development trajectory.")
    rows = []
    for tag in np.unique(stage):
        idx = int(np.flatnonzero(stage == tag)[-1])
        report = score_material(CASE_ID, m[idx], map_path=terminal_map if idx == len(step) - 1 else None)
        rows.append(dict(step=int(step[idx]), stage=int(tag), **report))
    report = dict(
        scope="development Twin stage endpoints, fixed by the sealed stage layout; no held-field selection",
        original_trace_status=receipt["stop"], original_trace_sealed_utc=receipt["sealed_utc"],
        original_trace_final_step=receipt["final_step"], rows=rows, **elapsed(started),
    )
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    for row in rows:
        print(row["step"], row["train"]["relative_to_measured"], row["held_transmit"]["relative_to_measured"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace", type=Path, default=TRACE)
    parser.add_argument("--output", type=Path, default=HERE / "dev_twin_stage_scores.json")
    parser.add_argument("--local-terminal-map", type=Path)
    args = parser.parse_args()
    with threadpool_limits(1):
        run(args.trace, args.output, args.local_terminal_map)
