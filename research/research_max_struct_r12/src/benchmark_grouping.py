"""One existing Pro-toy background projection workload, v1 versus optional v2.

Uses only public operators and stored observations. No inversion or GT scoring.
"""
from pathlib import Path
import importlib.util
import json
import time

import numpy as np

from field_current import DenseFieldCurrentProjectors


def main():
    here = Path(__file__).parent
    source = here.parent / "sources/pro/MAX_STRUCT_RESEARCH/feasibility"
    spec = importlib.util.spec_from_file_location("pro_toy_unmodified", source / "consensus_reflection_toy.py")
    toy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(toy)
    _, G, E, S = toy.setup()
    with np.load(source / "ring_high.npz", allow_pickle=False) as data:
        Y = data["Y"]
    # This is a throughput workload with a declared public tolerance, not a new
    # reconstruction using a changed scientific noise contract.
    delta = .01 * np.linalg.norm(Y)
    backend = toy.Physics(G, E, S, Y)
    _, J, Q = backend.forward(np.full(len(G), .1))
    repeats = 200
    report = {"scope": "one Pro scalar toy fixed-background projection throughput workload",
              "dofs": len(G), "sources": E.shape[1], "receivers": len(S),
              "projection_repeats": repeats, "delta": float(delta),
              "physical_initialization_lu": backend.lu_count,
              "physical_initialization_rhs": backend.rhs_count,
              "timing_limit": "single local run; other research processes may share CPU",
              "variants": {}}
    outputs = []
    for flag in (False, True):
        p = DenseFieldCurrentProjectors(G, E, S, Y, delta, grouped_ports=flag)
        initial_costs = dict(p.costs)
        start_wall, start_cpu = time.perf_counter(), time.process_time()
        for _ in range(repeats):
            result = p.affine_ball(J / p.scale, Q)
        wall, cpu = time.perf_counter() - start_wall, time.process_time() - start_cpu
        outputs.append(result)
        report["variants"]["grouped_v2" if flag else "frozen_v1"] = {
            "setup_cpu_seconds": p.setup_cpu_seconds,
            "setup_wall_seconds": p.setup_wall_seconds,
            "projection_cpu_seconds": cpu, "projection_wall_seconds": wall,
            "costs": {key: value - initial_costs[key] for key, value in p.costs.items()}}
    differences = [np.linalg.norm(a - b) / max(1., np.linalg.norm(a))
                   for a, b in zip(*outputs)]
    assert max(differences) < 1e-11
    report["relative_v1_v2_output_difference"] = [float(v) for v in differences]
    report["cpu_speedup_this_workload"] = (
        report["variants"]["frozen_v1"]["projection_cpu_seconds"]
        / report["variants"]["grouped_v2"]["projection_cpu_seconds"])
    (here / "projection_grouping_benchmark.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
