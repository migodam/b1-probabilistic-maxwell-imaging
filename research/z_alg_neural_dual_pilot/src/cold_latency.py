"""Compare sealed 64-block proposals using an independent native verifier per row."""

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import platform
from statistics import median

from holdout import observations
from native import Native
from projected_upper import NativeRestriction
from score_batch import PROPOSALS, RESULTS, atomic_json, clock, elapsed, evaluate, find_open, inputs, read_new_lines
from threadpoolctl import threadpool_limits


CASES = (
    ("optfresh_01", 2), ("optfresh_03", 2), ("optfresh_05", 2),
    ("s2d_008", 10), ("s2d_012", 10), ("s2d_016", 10),
)
ARMS = ("neural", "hybrid", "gradient1", "gradient5")
CAP = 0.5
BLOCK_SIDE = 8


def sealed_trials():
    selection = set(CASES)
    fields = {
        (oid, noise): (radius, y, caps)
        for oid, noise, radius, y, caps in observations()
        if (oid, noise) in selection
    }
    if set(fields) != selection:
        raise ValueError("The selected original primary views were not found.")

    for index, (oid, noise) in enumerate(CASES):
        path = PROPOSALS / f"{oid}.jsonl"
        rows, end = read_new_lines(path, 0)
        if end != path.stat().st_size or not rows or rows[-1].get("_complete") is not True:
            raise ValueError(f"The selected proposal file {path} is not sealed.")
        records = [
            row for row in rows if row.get("arm") in ARMS and row.get("noise") == noise
            and row.get("side") == BLOCK_SIDE and row.get("u") == CAP
        ]
        by_arm = {record["arm"]: record for record in records}
        if set(by_arm) != set(ARMS) or len(records) != len(ARMS):
            raise ValueError(f"Expected one sealed proposal for each of four arms on {oid}.")
        offset = index % len(ARMS)
        for arm in ARMS[offset:] + ARMS[:offset]:
            record = by_arm[arm]
            values = inputs(record, oid, {noise: fields[(oid, noise)]})
            yield path, record, values


def evaluate_cold(record, values):
    native_start = clock()
    verifier = Native()
    native_setup = elapsed(native_start)

    restriction_start = clock()
    restriction = NativeRestriction(16)
    restricted_setup = elapsed(restriction_start)

    outcome = evaluate(record, values, restriction, verifier)
    complete = elapsed(native_start)
    cache = verifier.cache_summary()
    operations = cache["operation_totals"]
    if any(operations[key] for key in ("profile_hits", "weyl_hits", "witness_hits")):
        raise ArithmeticError("A cold verifier unexpectedly found a stored multiplier proof.")
    return dict(
        **outcome, native_setup=native_setup, restriction_setup=restricted_setup,
        native_setup_and_cache=cache, restriction_setup_detail=restriction.setup,
        complete_native_wall=complete["wall"], complete_native_cpu=complete["cpu"],
        end_to_end_wall=record["proposal_wall"] + complete["wall"],
        end_to_end_cpu=record["proposal_cpu"] + complete["cpu"],
    )


def summary(rows):
    successful_pairs = {"neural": [], "hybrid": []}
    controls = []
    for oid, noise in CASES:
        group = {row["proposal"]["arm"]: row for row in rows if row["proposal"]["oid"] == oid}
        baseline = min(
            (row for arm, row in group.items() if arm.startswith("gradient") and row["passed"]),
            key=lambda row: row["end_to_end_wall"], default=None,
        )
        if baseline is None:
            continue
        controls.append(dict(
            oid=oid, noise=noise, arm=baseline["proposal"]["arm"], wall=baseline["end_to_end_wall"],
            cpu=baseline["end_to_end_cpu"], lower=baseline["lower"], radius=baseline["proposal"]["radius"],
        ))
        for arm, pairs in successful_pairs.items():
            proposal = group.get(arm)
            if proposal is None or not proposal["passed"]:
                continue
            pairs.append(dict(
                oid=oid, noise=noise, direct_arm=baseline["proposal"]["arm"],
                direct_wall=baseline["end_to_end_wall"], direct_cpu=baseline["end_to_end_cpu"],
                direct_lower=baseline["lower"], comparison_wall=proposal["end_to_end_wall"],
                comparison_cpu=proposal["end_to_end_cpu"], comparison_lower=proposal["lower"],
                radius=proposal["proposal"]["radius"],
                speedup=baseline["end_to_end_wall"] / proposal["end_to_end_wall"],
            ))
    by_arm = {}
    for arm in ARMS:
        selected = [row for row in rows if row["proposal"]["arm"] == arm]
        if selected:
            by_arm[arm] = dict(
                n=len(selected), statuses=dict(Counter(row["status"] for row in selected)),
                passed=sum(row["passed"] for row in selected),
                median_end_to_end_wall=median(row["end_to_end_wall"] for row in selected),
                median_end_to_end_cpu=median(row["end_to_end_cpu"] for row in selected),
                median_cold_native_wall=median(row["complete_native_wall"] for row in selected),
            )
    comparisons = {
        arm: dict(n=len(pairs), pairs=pairs,
                  median_speedup=median(r["speedup"] for r in pairs) if pairs else None,
                  n_at_least_5x=sum(r["speedup"] >= 5.0 for r in pairs))
        for arm, pairs in successful_pairs.items()
    }
    return dict(by_arm=by_arm, paired_fastest_passing_direct=controls, comparisons=comparisons)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=RESULTS / "cold_latency.json")
    args = parser.parse_args()
    start = clock()
    freeze = find_open(False, 0)
    schedule = list(sealed_trials())
    prep = elapsed(start)
    if len(schedule) != len(CASES) * len(ARMS):
        raise ValueError("The preselected independent comparison must contain exactly 24 sealed rows.")
    setup_records = {
        path.name: json.loads(path.read_text()) for path in sorted(PROPOSALS.glob("setup_*.json"))
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    doc = dict(
        status="in_progress", expected_rows=len(schedule), rows=[],
        case_order=CASES, base_arm_order=ARMS, cap=CAP, blocks=BLOCK_SIDE**2,
        schedule="one worker; contiguous object groups; rotate initial arm by object index",
        cache_scope="fresh Native and NativeRestriction and all fine geometry per score; OS file cache allowed",
        timing_contract="Original frozen proposal cost plus fresh verifier construction and native evaluation.",
        proposal_cpu_note="The gradient1 CPU entry is the recorded CPU through the last 1 s completed trial.",
        proposal_setup_scope=(
            "Original combined model/checkpoint, coarse geometry and readout setup is listed "
            "separately because its numerator includes both networks and shared inputs."
        ),
        proposal_worker_setup=setup_records,
        sealed_at=freeze.get("time_utc"), machine=platform.platform(),
        workers=1, requested_blas_threads=1, read_and_prepare=prep,
    )
    previous_elapsed = {"wall": 0., "cpu": 0.}
    completed_ids = set()
    if args.output.is_file():
        old = json.loads(args.output.read_text())
        if old["sealed_at"] != doc["sealed_at"]:
            raise ValueError("The stored proposals belong to another training freeze.")
        doc["rows"] = old["rows"]
        completed_ids = {row["proposal"]["id"] for row in doc["rows"]}
        previous_elapsed = old["elapsed"]
    with threadpool_limits(limits=1):
        for position, (source, record, values) in enumerate(schedule, 1):
            if record["id"] in completed_ids:
                continue
            load = os.getloadavg()
            result = evaluate_cold(record, values)
            row = dict(
                turn=position, proposal_file=str(source.relative_to(RESULTS.parent)),
                proposal=record, system_load_at_start=load, **result,
            )
            doc["rows"].append(row)
            doc.update(elapsed=elapsed(start), summary=summary(doc["rows"]))
            atomic_json(args.output, doc)
            print(
                position, record["oid"], record["arm"], row["status"],
                "pass" if row["passed"] else "", "cold wall", round(row["complete_native_wall"], 3),
                "end-to-end", round(row["end_to_end_wall"], 3), flush=True,
            )
    if len(doc["rows"]) != 24 or len({row["proposal"]["id"] for row in doc["rows"]}) != 24:
        raise ValueError("The cold latency replay did not produce 24 unique rows.")
    doc.update(status="complete", elapsed={unit: previous_elapsed[unit]+v for unit,v in elapsed(start).items()})
    atomic_json(args.output, doc)


if __name__ == "__main__":
    main()
