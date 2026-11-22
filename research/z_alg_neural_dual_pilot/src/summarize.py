"""Combine native certificates, fixed controls, and all paid proposal and replay costs."""

import argparse
from collections import Counter, defaultdict
import csv
import json
from math import sqrt
from pathlib import Path
from statistics import median

BASE = Path(__file__).resolve().parents[1]
RESULT = BASE / "results"
ARMS = ("gradient1", "gradient5", "neural", "hybrid", "no_feature")


def read_json_lines(path):
    if not path.exists():
        return []
    raw = path.read_bytes().split(b"\n")
    return [json.loads(x) for x in raw[:-1] if x.strip()]


def median_of(values):
    values = list(values)
    return median(values) if values else None


def scope(oid):
    return "Fresh6" if oid.startswith("optfresh_") else "SRAP22"


def classify(lo, hi, r):
    if lo > r:
        return "PASS"
    if hi == 0.0:
        return "ZERO"
    if hi is not None and hi <= r:
        return "BELOW"
    return "UNKNOWN"


def learned_rows(study, directory):
    out = []
    for path in sorted((directory / "native").glob("*.jsonl")):
        improvements = {
            r["id"]: r for r in read_json_lines(directory / "recovery" / path.name)
        }
        for r in read_json_lines(path):
            rid = r["id"]
            hi, lo, radius = r["upper"], r["lower"], r["radius"]
            nu, nlow = r["nu"], r["nu_lower"]
            first_decision = classify(lo, hi, radius)
            repair = improvements.get(rid)
            fee = {"wall": 0.0, "cpu": 0.0}
            extra_independent = dict(fee)
            if repair:
                for unit in fee:
                    fee[unit] = repair["proof_" + unit]
                    extra_independent[unit] = repair.get(
                        "charged_proof_" + unit, fee[unit]
                    )
                if repair.get("upper") is not None:
                    hi = min(hi, repair["upper"]) if hi is not None else repair["upper"]
                lo = max(lo, repair.get("lower") or 0.0)
                nu = repair["nu"] if repair.get("nu") is not None else nu
                nlow = max(nlow or 0, repair.get("nu_lower") or 0)
            decision = classify(lo, hi, radius)
            if hi is not None and lo > hi:
                raise ValueError(f"Incompatible independent intervals on {rid}")
            proposal_wall, proposal_cpu = r["proposal_wall"], r["proposal_cpu"]
            proof = r.get("native") or {}
            state = proof.get("proof") or {}
            base = dict(
                study=study,
                cohort=scope(r["oid"]),
                oid=r["oid"],
                noise=r["noise"],
                u=r["u"],
                side=r["side"],
                arm=r["arm"],
                radius=radius,
                lower=lo,
                upper=hi,
                decision=decision,
                passed=int(decision == "PASS"),
                coarse_d=sqrt(r["coarse_d2"]),
                coarse_nu=r.get("coarse_nu"),
                nu=nu,
                nu_lower=nlow,
                proposal_wall=proposal_wall,
                proposal_cpu=proposal_cpu,
                proof_wall=r["proof_wall"] + fee["wall"],
                proof_cpu=r["proof_cpu"] + fee["cpu"],
                charged_wall=proposal_wall
                + r["charged_proof_wall"]
                + extra_independent["wall"],
                charged_cpu=proposal_cpu
                + r["charged_proof_cpu"]
                + extra_independent["cpu"],
                base_decision=first_decision,
                recovery_attempt=int(bool(repair)),
                recovery_wall=fee["wall"],
                recovery_cpu=fee["cpu"],
                state_method=state.get("method"),
                score_reuse_id=r.get("score_reuse_id"),
                prefilter=(r.get("upper_check") or {}).get("status"),
                id=rid,
            )
            out.append(base)
    return out


def fixed_rows(study, directory):
    file = directory / "fixed_baselines.json"
    if not file.exists():
        return [], {}
    source = json.loads(file.read_text())
    out = []
    shared = {
        (p["cap"], p["profile_id"]): p["shared_timing"] for p in source["profiles"]
    }
    for k, family in (("rows", False), ("portfolio", True)):
        for p in source[k]:
            arm = "fixed_portfolio" if family else p["profile_id"]
            time_key = "online_" if family else ""
            decision = classify(p["lower"], p["upper"], p["padded_radius"])
            out.append(
                dict(
                    study=study,
                    cohort=scope(p["object_id"]),
                    oid=p["object_id"],
                    noise=p["noise_pct"],
                    u=p["cap"],
                    side=4 if arm.startswith("fixed_16") else 0,
                    arm=arm,
                    radius=p["padded_radius"],
                    lower=p["lower"],
                    upper=p["upper"],
                    decision=decision,
                    passed=int(decision == "PASS"),
                    coarse_d=None,
                    coarse_nu=None,
                    nu=p.get("nu"),
                    nu_lower=p.get("nu_lower"),
                    proposal_wall=0.0,
                    proposal_cpu=0.0,
                    proof_wall=p[time_key + "wall_seconds"],
                    proof_cpu=p[time_key + "cpu_seconds"],
                    charged_wall=p[time_key + "wall_seconds"],
                    charged_cpu=p[time_key + "cpu_seconds"],
                    base_decision=decision,
                    recovery_attempt=0,
                    recovery_wall=0.0,
                    recovery_cpu=0.0,
                    state_method="offline_fixed_state",
                    score_reuse_id=None,
                    prefilter=None,
                    id=p["condition"] + "__" + arm,
                )
            )
    return out, {
        "whole_replay_wall": source["wall_seconds"],
        "whole_replay_cpu": source["cpu_seconds"],
        "fixed_profile_setup": [
            dict(cap=u, arm=arm, **cost) for (u, arm), cost in shared.items()
        ],
    }


def strata(rows, group_keys):
    bins = defaultdict(list)
    for row in rows:
        bins[tuple(row[key] for key in group_keys)].append(row)
    output = []
    for key, group in sorted(bins.items()):
        counts = Counter(x["decision"] for x in group)
        has_coarse = [x for x in group if x["coarse_d"] is not None]
        coarse_pass = [x for x in has_coarse if x["coarse_d"] > x["radius"]]
        agreed = [x for x in group if x["passed"]]
        output.append(
            dict(
                zip(group_keys, key),
                records=len(group),
                materials=len({x["oid"] for x in group}),
                **{k: counts[k] for k in ("PASS", "BELOW", "ZERO", "UNKNOWN")},
                coarse_pass=len(coarse_pass),
                coarse_pass_fine_yes=sum(x["passed"] for x in coarse_pass),
                coarse_pass_proved_lost=sum(
                    x["decision"] in ("BELOW", "ZERO") for x in coarse_pass
                ),
                coarse_pass_fine_unknown=sum(
                    x["decision"] == "UNKNOWN" for x in coarse_pass
                ),
                coarse_weak_fine_pass=sum(
                    x["passed"] and x["coarse_d"] <= x["radius"] for x in has_coarse
                ),
                proven_at_least_ten_neg=sum(
                    (x["nu_lower"] or 0) >= 10 or (x["nu"] or 0) >= 10 for x in group
                ),
                exact_fewer_than_ten=sum(
                    x["nu"] is not None and x["nu"] < 10 for x in group
                ),
                exact_nu=sum(x["nu"] is not None for x in group),
                hot_wall_median=median_of(
                    x["proposal_wall"] + x["proof_wall"] for x in group
                ),
                charged_wall_median_passing=median_of(
                    x["charged_wall"] for x in agreed
                ),
                proposal_median=median_of(x["proposal_wall"] for x in group),
                normalized_margin_median=median_of(
                    x["lower"] / x["radius"] for x in group
                ),
                object_ids_passing=sorted({x["oid"] for x in agreed}),
            )
        )
    return output


def paired(rows):
    by_key = {(r["oid"], r["noise"], r["u"], r["side"], r["arm"]): r for r in rows}
    pairs = []
    for r in rows:
        if r["arm"] not in ("neural", "hybrid") or not r["passed"]:
            continue
        key = r["oid"], r["noise"], r["u"]
        candidates = [
            by_key.get((*key, r["side"], a)) for a in ("gradient1", "gradient5")
        ]
        both = [a for a in candidates if a is not None and a["passed"]]
        fixed = by_key.get((*key, 0, "fixed_portfolio"))
        pairs.append(
            dict(
                study=r["study"],
                cohort=r["cohort"],
                oid=r["oid"],
                noise=r["noise"],
                u=r["u"],
                side=r["side"],
                arm=r["arm"],
                seconds=r["charged_wall"],
                vs_gradient1=candidates[0]["charged_wall"] / r["charged_wall"]
                if candidates[0] and candidates[0]["passed"]
                else None,
                vs_gradient5=candidates[1]["charged_wall"] / r["charged_wall"]
                if candidates[1] and candidates[1]["passed"]
                else None,
                vs_fastest_passing_gradient=min(
                    (g["charged_wall"] for g in both), default=0.0
                )
                / r["charged_wall"]
                if both
                else None,
                vs_cached_fixed=fixed["charged_wall"] / r["charged_wall"]
                if fixed and fixed["passed"]
                else None,
            )
        )
    return pairs


def write_csv(file, rows):
    if not rows:
        return
    with file.open("w", newline="", encoding="utf-8") as stream:
        wr = csv.DictWriter(stream, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)


def main(preview):
    all_rows, fixed_cost, ledgers = [], {}, {}
    for dataset, directory, expected in (
        ("primary", RESULT, 980),
        ("supplement", RESULT / "supplement", 440),
    ):
        dynamic = learned_rows(dataset, directory)
        controls, fixed_cost[dataset] = fixed_rows(dataset, directory)
        recovery_workers = [
            json.loads(p.read_text())
            for p in sorted((directory / "recovery").glob("worker-*.json"))
        ]
        if len({r["id"] for r in dynamic}) != len(dynamic) or (
            not preview and len(dynamic) != expected
        ):
            raise ValueError(
                f"{dataset} expected {expected} unique physical receipts, got {len(dynamic)}"
            )
        if not preview and len(controls) != 4 * (expected // 10):
            raise ValueError(f"{dataset} has incomplete fixed baselines")
        if not preview:
            completed = {oid for w in recovery_workers if w["finished"] for oid in w["complete_objects"]}
            eligible = {
                r["id"] for path in (directory / "native").glob("*.jsonl")
                for r in read_json_lines(path)
                if r["status"] == "UNKNOWN" and r["lower"] <= 0
                and (r["upper"] is None or r["upper"] > r["radius"])
            }
            repaired = {r["id"] for r in dynamic if r["recovery_attempt"]}
            if completed != {r["oid"] for r in dynamic} or eligible != repaired:
                raise ValueError(f"{dataset} is awaiting complete matched recovery receipts")
        all_rows.extend(dynamic + controls)
        ledger = {
            "recovery_workers": recovery_workers,
            "native_workers": [
                json.loads(p.read_text())
                for p in sorted((directory / "native").glob("worker-*.json"))
            ],
            "proposal_workers": [
                json.loads(p.read_text())
                for p in sorted((directory / "proposals").glob("worker_*.json"))
            ],
            "proposal_setup": [
                json.loads(p.read_text())
                for p in sorted((directory / "proposals").glob("setup_*.json"))
            ],
            "diagnostic_wall": 0.0,
            "diagnostic_cpu": 0.0,
            "diagnostic_coarse_evals": 0,
        }
        for p in (directory / "proposals").glob("*.jsonl"):
            for row in read_json_lines(p):
                if "diagnostic_wall" in row:
                    for name in ("diagnostic_wall", "diagnostic_cpu"):
                        ledger[name] += row[name]
                    ledger["diagnostic_coarse_evals"] += 1
        ledgers[dataset] = ledger
    all_rows.sort(key=lambda r: (r["oid"], r["noise"], r["u"], r["side"], r["arm"]))
    suffix = "_preview" if preview else ""
    write_csv(RESULT / f"cases{suffix}.csv", all_rows)
    pairs = paired(all_rows)
    write_csv(RESULT / f"paired_timing{suffix}.csv", pairs)
    summary = dict(
        account="Averages are descriptive paired repeated fields. Object counts are distinct within each stratum.",
        material_count=len({r["oid"] for r in all_rows}),
        measured_field_count=len({(r["oid"], r["noise"]) for r in all_rows}),
        condition_count=len({(r["oid"], r["noise"], r["u"]) for r in all_rows}),
        dynamic_count=sum(r["arm"] in ARMS for r in all_rows),
        fixed_individual_count=sum(
            r["arm"] not in (*ARMS, "fixed_portfolio") for r in all_rows
        ),
        strata=strata(all_rows, ("cohort", "noise", "u", "side", "arm")),
        arm_totals=strata(all_rows, ("side", "arm")),
        costs=ledgers,
        fixed_setup=fixed_cost,
    )
    (RESULT / f"summary{suffix}.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    )
    print({k: v for k, v in summary.items() if k.endswith("count")})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview", action="store_true")
    main(parser.parse_args().preview)
