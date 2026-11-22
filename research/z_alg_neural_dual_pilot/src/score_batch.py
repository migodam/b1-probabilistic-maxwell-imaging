"""Replay sealed pilot proposals on the original 96² observations, shard by object."""

import argparse
from collections import Counter
import json
import math
import os
from pathlib import Path
import struct
import sys
from time import perf_counter, process_time, sleep

sys.dont_write_bytecode = True

import numpy as np


PILOT = Path(__file__).resolve().parents[1]
RESULTS = PILOT / "results"
PROPOSALS = RESULTS / "proposals"
RECEIPTS = RESULTS / "native"
OPEN = RESULTS / "holdout_open.json"
DONE = RESULTS / "proposals_done.json"
TIMING_LABELS = {
    "profile_operator_setup",
    "verify_cached_geometry",
    "native_binary64_inertia_proposal",
    "verify_native_negative_subspace",
    "build_regional_witness",
    "verify_Fourier_lower_envelope",
}


def clock():
    return perf_counter(), process_time()


def elapsed(start):
    return {"wall": perf_counter() - start[0], "cpu": process_time() - start[1]}


def atomic_json(path, payload):
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def read_new_lines(path, offset):
    """Consume complete flushed JSON lines, leaving an unfinished tail unread."""
    if not path.exists():
        return [], offset
    with path.open("rb") as stream:
        stream.seek(offset)
        raw = stream.read()
    end = raw.rfind(b"\n") + 1
    return [json.loads(line) for line in raw[:end].splitlines() if line], offset + end


def load_receipts(path):
    rows, offset = read_new_lines(path, 0)
    if path.exists() and path.stat().st_size != offset:
        with path.open("r+b") as handle:
            handle.truncate(offset)
    values = {row["id"]: row for row in rows}
    if len(values) != len(rows):
        raise ValueError(f"Repeated receipt IDs in {path.name}.")
    return values


def append_receipt(path, receipt):
    text = json.dumps(receipt, allow_nan=False, separators=(",", ":"))
    with path.open("a", encoding="utf-8") as stream:
        stream.write(text + "\n")
        stream.flush()
        os.fsync(stream.fileno())


def shard_pair(text):
    try:
        index, count = map(int, text.split("/"))
        if count <= 0 or not 0 <= index < count:
            raise ValueError(text)
        return index, count
    except ValueError as error:
        raise argparse.ArgumentTypeError("Specify a zero-based shard, for example 0/2.") from error


def find_open(watch, poll):
    while True:
        if OPEN.is_file():
            try:
                return json.loads(OPEN.read_text())
            except json.JSONDecodeError:
                if not watch:
                    raise
        if not watch:
            raise SystemExit("The official holdout-open receipt is required before reading observations.")
        sleep(poll)


def inputs(record, oid, observations):
    if record["oid"] != oid:
        raise ValueError("Proposal object ID disagrees with its source filename.")
    side = record["side"]
    if side not in (4, 8):
        raise ValueError("Expected a 4x4 or 8x8 profile.")
    noise = int(record["noise"])
    if noise != record["noise"]:
        raise ValueError("Nominal noise label must be an integer percentage.")
    radius, y, caps = observations[noise]
    physical_radius = float(record["radius"])
    if physical_radius != radius:
        raise ValueError("Proposal radius differs from the frozen observation radius.")
    u = float(record["u"])
    if u not in caps:
        raise ValueError("Proposal cap is outside this observation's frozen comparison.")
    mu = np.asarray(record["mu"], dtype=np.float64) / 9.
    lam = np.asarray(record["lam"], dtype=np.float64) / 9.
    if mu.shape != (side * side,) or lam.shape != mu.shape:
        raise ValueError("Expected flat row-major sector and optical profiles.")
    if not np.isfinite(mu).all() or not np.isfinite(lam).all() or np.min(mu) <= 0.:
        raise ValueError("Physical sector weights must be finite positive and optical weights finite.")
    key = (y.tobytes(order="C"), struct.pack("!ddi", u, radius, side),
           mu.tobytes(order="C"), lam.tobytes(order="C"))
    return y, u, mu, lam, radius, key


def base_row(record):
    return {label: record.get(label) for label in (
        "id", "oid", "noise", "u", "side", "arm", "radius", "proposal_wall",
        "proposal_cpu", "coarse_d2",
    )} | {"coarse_nu": record.get("nu")}


def profile_charge(result, costs):
    """Track reusable state-only costs and charge them to an independent profile."""
    native = result.get("native")
    if native is None:
        return {"wall": 0., "cpu": 0.}
    weights = native["profile"]
    key = (weights["cap"], tuple(weights["weights"]), tuple(weights["optical"]))
    timing = native["timing"]
    current = {unit: sum(timing[label][f"{unit}_seconds"] for label in TIMING_LABELS
                         if label in timing) for unit in ("wall", "cpu")}
    proof = native.get("proof")
    if proof is not None and proof.get("cache_hit"):
        previous = costs.get(key, current)
        return {unit: max(0., previous[unit] - current[unit]) for unit in current}
    costs[key] = current
    return {"wall": 0., "cpu": 0.}


def evaluate(record, values, restriction, verifier):
    y, u, mu, lam, radius, _ = values
    started = clock()
    row = {"upper_check": None, "native": None, "lower": 0., "upper": None, "nu": None,
           "nu_lower": 0, "status": "UNKNOWN", "passed": False, "score_reuse_id": None}
    try:
        upper_started = clock()
        try:
            row["upper_check"] = restriction.check(y, u, mu, lam, radius)
        except (ArithmeticError, RuntimeError, ValueError) as error:
            upper_elapsed = elapsed(upper_started)
            row["upper_check"] = {"status": "INCONCLUSIVE", "upper": None, "nu_lower": None,
                                  "explanation": str(error), "wall_s": upper_elapsed["wall"],
                                  "cpu_s": upper_elapsed["cpu"]}
        upper = row["upper_check"]
        if upper["status"] in ("PROVED_ZERO", "PROVED_BELOW_RADIUS"):
            row.update(upper=upper["upper"], status=upper["status"],
                       nu_lower=upper["nu_lower"] or 0)
        else:
            native = verifier.score(y, u, mu, lam, radius=radius, nu_hint=None)
            row.update(native=native, lower=native["lower"], upper=native["upper"],
                       nu=native["nu"], nu_lower=native["nu_lower"], status=native["status"],
                       passed=native["strictly_excludes"])
            if native["lower"] == 0. and native["upper"] is not None and native["upper"] <= radius:
                if row["status"] != "PROVED_ZERO":
                    row["status"] = "PROVED_BELOW_RADIUS"
    except (ArithmeticError, RuntimeError, ValueError) as error:
        row["explanation"] = str(error)
    spent = elapsed(started)
    row.update(proof_wall=spent["wall"], proof_cpu=spent["cpu"])
    return row


def failed_record(record, error, started):
    spent = elapsed(started)
    row = base_row(record)
    row.update(upper_check=None, native=None, lower=0., upper=None, nu=None, nu_lower=0,
               status="UNKNOWN", passed=False, score_reuse_id=None, explanation=str(error),
               proof_wall=spent["wall"], proof_cpu=spent["cpu"],
               charged_proof_wall=spent["wall"], charged_proof_cpu=spent["cpu"],
               score_reuse_saved_wall=0., score_reuse_saved_cpu=0., profile_reuse_wall=0.,
               profile_reuse_cpu=0.)
    return row


def process_proposal(record, source_id, observed, previous, cache, profile_costs, rejection_counts,
                     restriction, verifier):
    started = clock()
    try:
        values = inputs(record, source_id, observed)
    except (KeyError, TypeError, ValueError, OverflowError) as error:
        rejection_counts[type(error).__name__] += 1
        return previous or failed_record(record, error, started)
    key = values[-1]
    if previous is not None:
        if previous.get("score_reuse_id") is None:
            profile_charge(previous, profile_costs)
        cache.setdefault(key, previous)
        return previous
    if key in cache:
        reused = cache[key]
        # The inputs, radius and mapped fine profiles are exactly equal. Preserve
        # the full proof record, plus a distinct independent-arm and reuse charge.
        proof = {name: reused.get(name) for name in (
            "upper_check", "native", "lower", "upper", "nu", "nu_lower", "status", "passed",
            "explanation",
        )}
        actual = elapsed(started)
        proof.update(
            score_reuse_id=reused["id"], proof_wall=actual["wall"], proof_cpu=actual["cpu"],
            charged_proof_wall=reused["charged_proof_wall"],
            charged_proof_cpu=reused["charged_proof_cpu"],
            score_reuse_saved_wall=max(0., reused["proof_wall"] - actual["wall"]),
            score_reuse_saved_cpu=max(0., reused["proof_cpu"] - actual["cpu"]),
            profile_reuse_wall=0., profile_reuse_cpu=0.,
        )
    else:
        proof = evaluate(record, values, restriction, verifier)
        extra = profile_charge(proof, profile_costs)
        proof.update(
            charged_proof_wall=proof["proof_wall"] + extra["wall"],
            charged_proof_cpu=proof["proof_cpu"] + extra["cpu"],
            score_reuse_saved_wall=0., score_reuse_saved_cpu=0.,
            profile_reuse_wall=extra["wall"], profile_reuse_cpu=extra["cpu"],
        )
    row = base_row(record) | proof
    if key not in cache:
        cache[key] = row
    return row


def run(args):
    started = clock()
    find_open(args.watch, args.poll)
    after_open = clock()
    from holdout import observations
    from native import Native
    from projected_upper import NativeRestriction

    frozen = {}
    for oid, noise, radius, y, caps in observations(supplement=args.supplement):
        per_object = frozen.setdefault(oid, {})
        if noise in per_object:
            raise ValueError(f"Duplicate frozen observation {oid}, noise {noise}.")
        per_object[noise] = (radius, np.ascontiguousarray(y, dtype=np.complex128), caps)
    index, count = args.shard
    names = [name for position, name in enumerate(sorted(frozen)) if position % count == index]
    if args.object is not None:
        if not set(args.object) <= set(names):
            raise ValueError("The requested object is absent from the frozen shard.")
        names = [name for name in names if name in args.object]
    read_setup = elapsed(after_open)
    RECEIPTS.mkdir(parents=True, exist_ok=True)
    output = {name: RECEIPTS / f"{name}.jsonl" for name in names}
    persisted = {name: load_receipts(path) for name, path in output.items()}
    offsets = dict.fromkeys(names, 0)
    seen = {name: set() for name in names}
    seen_ids = set()
    complete = set()
    cache, profile_costs, rejection_counts = {}, {}, Counter()
    verifier = Native()
    restriction = NativeRestriction(16)
    print(f"Native shard {index}/{count} has {len(names)} objects.", flush=True)
    scope = args.object[0] if args.object and len(args.object) == 1 else f"{index}-of-{count}"
    summary = RECEIPTS / f"worker-{scope}.json"

    def write_summary(finished):
        rows = [row for name in names for row in persisted[name].values()]
        totals = {field: sum(r.get(field, 0.) for r in rows) for field in (
            "proof_wall", "proof_cpu", "charged_proof_wall", "charged_proof_cpu",
            "score_reuse_saved_wall", "score_reuse_saved_cpu", "profile_reuse_wall", "profile_reuse_cpu",
        )}
        payload = {
            "shard": [index, count],
            "objects": names,
            "complete_objects": sorted(complete),
            "finished": finished,
            "read_holdout_and_import_setup": read_setup,
            "restriction_setup": restriction.setup,
            "native_setup_and_cache": verifier.cache_summary(),
            "proof_status_counts": dict(Counter(row["status"] for row in rows)),
            "proof_rows": len(rows),
            "proof_totals": totals,
            "receipt_rejections": dict(rejection_counts),
            "elapsed": elapsed(started),
        }
        atomic_json(summary, payload)

    while True:
        aggregate_before_scan = DONE.is_file()
        any_new = False
        for name in names:
            if name in complete:
                continue
            rows, offsets[name] = read_new_lines(PROPOSALS / f"{name}.jsonl", offsets[name])
            object_finished = False
            for proposal in rows:
                if proposal.get("_complete") is True:
                    if object_finished:
                        raise ValueError(f"Repeated completion marker for {name}.")
                    object_finished = True
                    continue
                if object_finished:
                    raise ValueError(f"Proposal appended after completion for {name}.")
                identifier = proposal["id"]
                if not isinstance(identifier, str) or not identifier or identifier in seen_ids:
                    raise ValueError(f"Malformed or repeated proposal identifier in {name}.")
                seen[name].add(identifier)
                seen_ids.add(identifier)
                previous = persisted[name].get(identifier)
                result = process_proposal(
                    proposal, name, frozen[name], previous, cache, profile_costs,
                    rejection_counts, restriction, verifier
                )
                if previous is None:
                    any_new = True
                    append_receipt(output[name], result)
                    persisted[name][identifier] = result
            if object_finished:
                complete.add(name)
                any_new = True
                print(f"{name}: {len(persisted[name])} native receipts complete.", flush=True)
        if aggregate_before_scan:
            # By this scan's start all final appends were already complete.
            complete.update(name for name in names if (PROPOSALS / f"{name}.jsonl").is_file())
        for name in complete:
            if seen[name] != set(persisted[name]):
                raise ValueError(f"Native receipts and sealed proposal IDs disagree for {name}.")
        finished = len(complete) == len(names)
        if any_new or finished or not args.watch or aggregate_before_scan:
            write_summary(finished)
        if finished:
            print(f"Native shard {index}/{count} complete.", flush=True)
            return 0
        if not args.watch or aggregate_before_scan:
            print(f"Native shard {index}/{count} is waiting for {len(names)-len(complete)} proposal files.",
                  flush=True)
            return 2 if aggregate_before_scan else 0
        if not any_new:
            sleep(args.poll)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shard", type=shard_pair, default=(0, 1))
    parser.add_argument("--object", action="append")
    parser.add_argument("--watch", action="store_true")
    parser.add_argument("--supplement", action="store_true", help="Replay the two independent low-noise SRAP fields per object.")
    parser.add_argument("--poll", type=float, default=2.)
    config = parser.parse_args()
    if config.supplement:
        PROPOSALS = RESULTS / 'supplement' / 'proposals'
        RECEIPTS = RESULTS / 'supplement' / 'native'
        DONE = RESULTS / 'supplement' / 'proposals_done.json'
    if not math.isfinite(config.poll) or config.poll <= 0.:
        parser.error("--poll must be positive.")
    raise SystemExit(run(config))
