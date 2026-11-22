"""On validation fields only, change which ports condition the frozen network.

The physics objective always sees the original Y. All results are 32-grid
diagnostics. The model, the primary scores, and all native arms stay frozen.
"""

import json
import os
from pathlib import Path
from time import perf_counter, process_time

for variable in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
import torch
from threadpoolctl import threadpool_limits
from coarse import Coarse
from network import load_network, predict

PILOT = Path(__file__).resolve().parents[1]


def main():
    start, cpu = perf_counter(), process_time()
    torch.set_num_threads(1)
    with np.load(PILOT / "data/synthetic_observations.npz", allow_pickle=False) as data:
        val = np.flatnonzero(data["split"] == "validation")
        levels = np.flatnonzero(np.isin(data["noise_fraction"], [0.01, 0.02]))
        ys, rs, ids = (data[key][val] for key in ("Y", "radius", "object_id"))
        fractions = data["noise_fraction"]
    if len(ys) != 12 or len(levels) != 2:
        raise ValueError("Expected 12 separate validation objects at both levels.")
    co = Coarse()
    net = load_network(PILOT / "results/neural_dual_weights.pt")
    rows = []
    for i, oid in enumerate(ids):
        # Same phase scrambling across paired noise levels of a single object.
        rng = np.random.default_rng(np.random.SeedSequence([11943, i]))
        phase = np.exp(2j * np.pi * rng.random((20, 10)))
        for j in levels:
            y = ys[i, j]
            other = ys[(i + 1) % len(ids), j]
            inputs = {
                "matched": y,
                "swapped_same_energy": other
                * (np.linalg.norm(y) / np.linalg.norm(other)),
                "phase_scrambled": y * phase,
            }
            for u in (0.5, 0.75, 1.0):
                for side in (4, 8):
                    for name, source_y in inputs.items():
                        mu, lam = predict(net, source_y, u, side, noise_radius=rs[i, j])
                        result = co.eval(y, u, mu, lam)
                        rows.append(
                            dict(
                                object_id=str(oid),
                                cap=u,
                                noise_pct=round(100 * fractions[j]),
                                side=side,
                                control=name,
                                nu=result["nu"],
                                normalized=float(
                                    result["d2_exact"] / np.vdot(y, y).real
                                ),
                                ratio=float(np.sqrt(result["d2_exact"]) / rs[i, j]),
                            )
                        )
    out = {
        "scope": "validation only, same target energy and radius under all conditions; coarse binary64 diagnostic",
        "wall_s": perf_counter() - start,
        "cpu_s": process_time() - cpu,
        "physics_evaluations": len(rows),
        "rows": rows,
    }
    (PILOT / "results/validation_swaps.json").write_text(
        json.dumps(out, indent=2, allow_nan=False) + "\n"
    )
    for cap in (0.5, 0.75, 1.0):
        print(
            cap,
            {
                arm: round(
                    float(
                        np.mean(
                            [
                                r["normalized"]
                                for r in rows
                                if r["cap"] == cap and r["control"] == arm
                            ]
                        )
                    ),
                    7,
                )
                for arm in inputs
            },
            flush=True,
        )


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        main()
