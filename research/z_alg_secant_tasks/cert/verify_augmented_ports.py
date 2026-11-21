"""Audit the certified baseline material pairs on the two matched acquisitions."""

import argparse
import json
from pathlib import Path
from time import perf_counter, process_time

from finite_chord import HERE
from array_secant import make_models, save_json
from interval_physics import IntervalModel, arb, ball_norm, bounds, point


def run(source, output):
    start_wall, start_cpu = perf_counter(), process_time()
    input_data = json.loads(source.read_text())
    models = {}
    results = []
    for case in input_data["cases"]:
        n = case["n"]
        if n not in models:
            models[n] = IntervalModel(make_models(n)[0].centers, n, with_frequency=True)
        physical = models[n]
        signal = None
        for row in case["arms"]:
            certificate = row.get("certificate", {})
            if certificate.get("status") != "certified_exact_root_in_chart":
                continue
            if signal is None:
                base = certificate["noise_10pct_radius"]
                signal = point(base["lo"]).union(point(base["hi"]))
            minus, plus = [
                [point(x["lo"]).union(point(x["hi"])) for x in certificate[f"material_{sign}_enclosure"]]
                for sign in ("minus", "plus")
            ]
            for arm in ("angle", "frequency"):
                ym, _ = physical.response(minus, arm)
                yp, _ = physical.response(plus, arm)
                distance = ball_norm(a - b for a, b in zip(yp, ym, strict=True))
                ratio = distance / (2 * arb(2).sqrt() * signal)
                decision = "unresolved"
                if ratio > arb(1):
                    decision = "pair_excluded"
                if ratio < arb(1):
                    decision = "common_midpoint_exists"
                results.append({"baseline_pair": row["id"], "arm": arm,
                                "roi_gap_enclosure": certificate["task_gap_plus_minus"],
                                "chord_over_2delta_enclosure": bounds(ratio), "decision": decision})
    save_json(output, {"model_inputs": source.name, "results": results,
                       "interval_state_solves": sum(model.solve_count for model in models.values()),
                       "wall_seconds": perf_counter() - start_wall,
                       "cpu_seconds": process_time() - start_cpu})


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--input", type=Path, default=HERE / "results" / "finite_chord.json")
    cli.add_argument("--output", type=Path, default=HERE / "results" / "certified_pair_augmented_ports.json")
    args = cli.parse_args()
    run(args.input, args.output)
