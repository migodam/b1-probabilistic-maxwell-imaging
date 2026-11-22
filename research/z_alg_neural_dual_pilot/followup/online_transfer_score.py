"""Strict native replay and paired comparison of the frozen adaptive Transfer16 arm."""

import argparse
from collections import Counter
import gzip
import json
from statistics import median

from online_transfer_propose import CAP, FROZEN, NOISES, PILOT, RESULTS, identifier, observations
from score_batch import append_receipt, atomic_json, clock, elapsed, load_receipts, process_proposal, shard_pair


def proposals(workers):
    rows, setup = [], []
    frozen = json.loads(FROZEN.read_text())
    for rank in range(workers):
        path = RESULTS / f"online_transfer_proposals-{rank}-of-{workers}.json.gz"
        with gzip.open(path, "rt") as stream:
            data = json.load(stream)
        if data["freeze_utc"] != frozen["freeze_utc"] or data["shard"] != [rank, workers]:
            raise ValueError("Proposal receipt is unrelated to this frozen run.")
        rows.extend(data.pop("rows"))
        setup.append(data)
    indexed = {row["id"]: row for row in rows}
    if len(indexed) != len(rows) or len(rows) != 48:
        raise ValueError("All 48 unique proposal receipts are required.")
    return indexed, setup


def native_verdict(lower, upper, radius):
    return "PASS" if lower > radius else "BELOW_OR_EQUAL" if upper is not None and upper <= radius else "UNKNOWN"


def verify(shard, proposal_workers):
    from flint import ctx
    from native import Native
    from projected_upper import NativeRestriction
    from recover_unknown import RecoveryNative, eligible, examine
    from threadpoolctl import threadpool_limits

    started = clock()
    cases = observations()
    sealed, _ = proposals(proposal_workers)
    rank, workers = shard
    names = sorted(cases)[rank::workers]
    path = RESULTS / f"online_transfer_native-{rank}-of-{workers}.jsonl"
    old = load_receipts(path)
    ctx.prec, ctx.threads = 112, 1
    primary, recovery = Native(threads=1), RecoveryNative()
    recovery.threads = 1
    reuse, profile_costs, rejected = {}, {}, Counter()
    recovery_reuse, recovery_costs = {}, {}
    with threadpool_limits(limits=1):
        restriction = NativeRestriction()
        setup = elapsed(started)
        for oid in names:
            for noise in NOISES:
                record = sealed[identifier(oid, noise)]
                if record["id"] in old:
                    continue
                proof = process_proposal(
                    record, oid, cases[oid], None, reuse, profile_costs, rejected, restriction, primary,
                )
                append_receipt(path, {"id": record["id"], "primary": proof, "recovery": None, "stage": "primary"})
                print(record["id"], "primary", proof["status"], proof["lower"], proof["upper"], flush=True)
        initial = load_receipts(path)
        recovery_path = RESULTS / f"online_transfer_recovery-{rank}-of-{workers}.jsonl"
        recovered = load_receipts(recovery_path)
        for oid in names:
            for noise in NOISES:
                record = sealed[identifier(oid, noise)]
                key = record["id"]
                if key in recovered:
                    continue
                proof = initial[key]["primary"]
                if eligible(proof) and record["mu"] is not None:
                    extra = examine(record, proof, oid, cases[oid], recovery, recovery_reuse, recovery_costs)
                    append_receipt(recovery_path, {"id": key, "recovery": extra})
                    print(key, "recovery", extra["status"], extra["lower"], extra["upper"], flush=True)
    atomic_json(RESULTS / f"online_transfer_worker-{rank}-of-{workers}.json", {
        "shard": shard, "objects": names, "setup": setup, "restriction_setup": restriction.setup,
        "native_cache": primary.cache_summary(), "recovery_cache": recovery.cache_summary(),
        "primary_rejection_counts": dict(rejected), **elapsed(started),
    })


def read_score_shards(workers):
    primary, recovery, metadata = {}, {}, []
    for rank in range(workers):
        suffix = f"{rank}-of-{workers}"
        report = RESULTS / f"online_transfer_worker-{suffix}.json"
        metadata.append(json.loads(report.read_text()))
        primary.update(load_receipts(RESULTS / f"online_transfer_native-{suffix}.jsonl"))
        recovery.update(load_receipts(RESULTS / f"online_transfer_recovery-{suffix}.jsonl"))
    if len(primary) != 48 or not set(recovery) <= set(primary):
        raise ValueError("Primary and recovery replay has not completed for all 48 views.")
    return primary, recovery, metadata


def paired(row, archived, proposal, extra):
    radius = proposal["radius"]
    primary = row["primary"]
    if radius != archived["radius"] or radius != primary["radius"] or (extra and radius != extra["radius"]):
        raise ValueError("Adaptive and fixed comparators do not use the same matched radius.")
    lo = max(primary["lower"], extra["lower"]) if extra else primary["lower"]
    bounds = [r["upper"] for r in (primary, extra) if r is not None and r["upper"] is not None]
    hi = min(bounds) if bounds else None
    if hi is not None and lo > hi:
        raise ArithmeticError("The first and second native proofs are inconsistent.")
    control = archived["arms"]["full_archived"]
    online = {
        unit: proposal[f"paid_{unit}"] + primary[f"proof_{unit}"] + (extra[f"proof_{unit}"] if extra else 0.)
        for unit in ("wall", "cpu")
    }
    fixed = {name: control[name] for name in ("lower", "upper", "best_profile", "wall_seconds", "cpu_seconds")}
    return {
        "id": proposal["id"], "object_id": proposal["oid"], "noise_pct": proposal["noise"], "radius": radius,
        "initial": {name: proposal["initial"][name] for name in ("mu", "lam", "d2")},
        "selected_coarse_weights": {name: proposal[name] for name in ("mu", "lam")},
        "selected_native_weights": {name: [v / 9 for v in proposal[name]] if proposal[name] is not None else None
                                    for name in ("mu", "lam")},
        "coarse_d2_diagnostic": proposal["coarse_d2"], "coarse_nu_diagnostic": proposal["nu"],
        "selection": {name: proposal[name] for name in (
            "paid_wall", "paid_cpu", "evals", "overrun_evals", "optimizer_stop", "last_eligible_wall",
            "best_completed_wall",
        )},
        "native": {"lower": lo, "upper": hi, "verdict": native_verdict(lo, hi, radius),
                   "wall_seconds": online["wall"], "cpu_seconds": online["cpu"]},
        "primary": primary, "recovery": extra,
        "fixed_archived_five": fixed | {"verdict": native_verdict(fixed["lower"], fixed["upper"], radius)},
    }


def status_summary(rows):
    def keys(predicate):
        return [[r["object_id"], r["noise_pct"]] for r in rows if predicate(r)]

    method, comparator = "native", "fixed_archived_five"
    timing = ("wall_seconds", "cpu_seconds")
    passes = {name: keys(lambda r: r[name]["verdict"] == "PASS") for name in (method, comparator)}
    return {
        "views": len(rows), "independent_objects": len({r["object_id"] for r in rows}),
        "verdicts": {name: dict(Counter(r[name]["verdict"] for r in rows)) for name in (method, comparator)},
        "pass_views": passes,
        "pass_objects": {name: sorted({key[0] for key in subset}) for name, subset in passes.items()},
        "gain_over_fixed_bound": keys(lambda r: r[method]["verdict"] == "PASS"
                                      and r[comparator]["verdict"] != "PASS"),
        "loss_vs_fixed_pass": keys(lambda r: r[comparator]["verdict"] == "PASS"
                                   and r[method]["verdict"] == "BELOW_OR_EQUAL"),
        "fixed_pass_adaptive_unknown": keys(lambda r: r[comparator]["verdict"] == "PASS"
                                           and r[method]["verdict"] == "UNKNOWN"),
        "all_adaptive_unknown": keys(lambda r: r[method]["verdict"] == "UNKNOWN"),
        "median_online": {name: {unit: median(r[name][unit] for r in rows) for unit in timing}
                          for name in (method, comparator)},
        "total_online": {name: {unit: sum(r[name][unit] for r in rows) for unit in timing}
                         for name in (method, comparator)},
    }


def assemble(workers, proposal_workers):
    begin = clock()
    frozen = json.loads(FROZEN.read_text())
    first, second, worker_reports = read_score_shards(workers)
    sealed, proposal_reports = proposals(proposal_workers)
    fixed_path = RESULTS / "transferred_weight_transfer16.json"
    fixed = json.loads(fixed_path.read_text())
    by_view = {(r["object_id"], r["noise_pct"]): r for r in fixed["rows"]}
    if fixed["cap"] != CAP or len(by_view) != len(sealed):
        raise ValueError("Mismatch against the saved fixed comparison.")
    rows = [
        paired(first[key], by_view[proposal["oid"], proposal["noise"]], proposal,
               second[key]["recovery"] if key in second else None)
        for key, proposal in sorted(sealed.items())
    ]
    groups = {"all_views": rows} | {
        f"{noise}_pct": [r for r in rows if r["noise_pct"] == noise] for noise in NOISES
    }
    ledger = json.loads((PILOT / "followup/transfer16_offline.json").read_text())
    low = {r["object_id"] for r in ledger["rows"] if r["under_cap"]}
    low_rows = [r for r in rows if r["object_id"] in low]
    high_rows = [r for r in rows if r["object_id"] not in low]
    if len(low_rows) != 12 or any(r["native"]["verdict"] == "PASS" for r in low_rows):
        raise ArithmeticError("Low-cap controls disagree with the independent sealed material audit.")
    groups.update(low_cap_controls=low_rows, above_cap_objects=high_rows)
    output = {
        "scope": frozen["scope"], "protocol": FROZEN.name, "cap": CAP,
        "composition": "48 matched views; candidate generation and native scoring complete before scalar label lookup",
        "fixed_control": fixed_path.name, "fixed_prep": fixed["preparation_cost_by_arm"]["full_archived"],
        "proposal_workers": proposal_reports, "verification_workers": worker_reports,
        "after_scoring_offline_low_cap_ids": sorted(low),
        "summary": {key: status_summary(group) for key, group in groups.items()},
        "rows": rows, **elapsed(begin),
    }
    atomic_json(RESULTS / "online_transfer_paired.json", output)
    print(json.dumps(output["summary"]["all_views"]), flush=True)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--shard", type=shard_pair, default=(0, 1))
    cli.add_argument("--proposal-workers", type=int, default=4)
    cli.add_argument("--assemble", type=int, metavar="VERIFICATION_WORKERS")
    args = cli.parse_args()
    if args.assemble:
        assemble(args.assemble, args.proposal_workers)
    else:
        verify(args.shard, args.proposal_workers)
