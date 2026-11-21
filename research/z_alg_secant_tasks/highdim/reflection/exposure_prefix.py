"""Replay a geometry-ordered stopping rule on archived reflection-closed shots."""

import argparse
from collections import defaultdict
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from flint import arb

from replay_symmetry import HERE, lift_complex, read_json, rounded_ball, squared_length, write


def orbits(symmetry):
    paired = dict(zip(symmetry["source_order"], symmetry["source_image_order"], strict=True))
    seen = set()
    groups = []
    for source in sorted(paired):
        if source in seen:
            continue
        partner = paired[source]
        if paired.get(partner) != source:
            raise ValueError("The source schedule is not closed under reflection.")
        orbit = tuple(sorted({source, partner}))
        groups.append(orbit)
        seen.update(orbit)
    return paired, groups


def scenario(y, plan, delta):
    sources, groups = orbits(plan)
    rec = plan["receiver_order"]
    re_image = plan["receiver_image_order"]
    total, record, d2, e2, o2 = 0, [], delta*delta, arb(0), arb(0)
    first = None
    half = None
    n_shots = y.shape[1]
    half_budget = n_shots//2
    for group in groups:
        new = lift_complex(y[np.ix_(rec, group)])
        reflected = lift_complex(y[np.ix_(re_image, [sources[j] for j in group])])
        e2 += squared_length([(p+q)/2 for p, q in zip(new, reflected, strict=True)])
        o2 += squared_length([(p-q)/2 for p, q in zip(new, reflected, strict=True)])
        total += len(group)
        passes = bool(e2 > d2 and o2 > d2)
        status = {"exposures": total, "passes": passes, "source_group": list(group),
                  "E_squared": rounded_ball(e2), "O_squared": rounded_ball(o2)}
        record.append(status)
        if passes and first is None:
            first = total
        if total <= half_budget:
            half = status
    if not record:
        raise ValueError("At least one complete source group is required.")
    return {
        "available_original_exposures": n_shots,
        "max_closed_exposures": total,
        "receivers_per_exposure": len(rec),
        "raw_full_exposures": n_shots,
        "whole_closed_set_certifies": record[-1]["passes"],
        "fixed_half_budget": half_budget,
        "fixed_half_actual_exposures": 0 if half is None else half["exposures"],
        "fixed_half_certifies": bool(half and half["passes"]),
        "adaptive_exposures": total if first is None else first,
        "adaptive_certifies": first is not None,
        "first_certifying_exposures": first,
        "prefix_trace": record,
    }


def aggregate(rows, case_key):
    buckets = defaultdict(list)
    for row in rows:
        if case_key in row and row.get("axis") is not None:
            buckets[row["cohort"]].append(row)
            buckets["all-applicable"].append(row)
    summary = {}
    for group, cluster in buckets.items():
        usable = [r for r in cluster if r[case_key]["within_actual_contract"] is not False]
        runs = [r[case_key] for r in usable]
        early = [r for r in runs if r["adaptive_certifies"]]
        pulse = sum(x["available_original_exposures"] for x in runs)
        by_layout = defaultdict(list)
        for row in usable:
            by_layout[row["layout"]].append(row[case_key])
        hindsight = {}
        for layout, vals in by_layout.items():
            prefix = max((a["first_certifying_exposures"] for a in vals if a["adaptive_certifies"]),
                         default=0)
            hindsight[layout] = {"views": len(vals), "fixed_prefix": prefix,
                                  "total_exposures": prefix*len(vals),
                                  "certifies": sum(a["adaptive_certifies"] for a in vals),
                                  "failures": sum(not a["adaptive_certifies"] for a in vals)}
        attempted = [r[case_key] for r in cluster]
        all_pulses = sum(x["available_original_exposures"] for x in attempted)
        all_adaptive = sum(x["adaptive_exposures"] for x in attempted)
        summary[group] = {
            "geometry_applicable_views": len(cluster),
            "within_contract_views": len(runs),
            "outside_actual_contract_views": len(cluster)-len(runs),
            "unqualified_apparent_adaptive_certifications": sum(
                r[case_key]["adaptive_certifies"] for r in cluster
                if r[case_key]["within_actual_contract"] is False),
            "full_certified": sum(r["whole_closed_set_certifies"] for r in runs),
            "half_certified": sum(r["fixed_half_certifies"] for r in runs),
            "adaptive_certified": len(early),
            "adaptive_no_certificate": len(runs)-len(early),
            "raw_full_exposures": pulse,
            "closed_full_exposures": sum(x["max_closed_exposures"] for x in runs),
            "fixed_half_exposures": sum(x["fixed_half_actual_exposures"] for x in runs),
            "adaptive_exposures_including_failures": sum(x["adaptive_exposures"] for x in runs),
            "adaptive_fraction_of_raw_full": sum(x["adaptive_exposures"] for x in runs)/pulse if pulse else None,
            "adaptive_successful_exposure_minmax": (
                [min(x["adaptive_exposures"] for x in early), max(x["adaptive_exposures"] for x in early)]
                if early else None),
            "hindsight_min_fixed_prefix_by_layout_for_same_successes": hindsight,
            "hindsight_total_exposures_including_failures": sum(x["total_exposures"] for x in hindsight.values()),
            "all_geometry_attempts_cost_including_out_of_contract": {
                "total_attempted_views": len(attempted), "full_raw_exposures": all_pulses,
                "fixed_half_exposures": sum(x["fixed_half_actual_exposures"] for x in attempted),
                "adaptive_exposures": all_adaptive,
                "adaptive_fraction_of_raw_full": all_adaptive/all_pulses if all_pulses else None,
                "hindsight_prefix_exposures_applied_to_every_attempt": sum(
                    hindsight.get(r["layout"], {}).get("fixed_prefix", 0) for r in cluster),
                "adaptive_exposures_spent_on_out_of_contract": sum(
                    x["adaptive_exposures"] for x in attempted if x["within_actual_contract"] is False),
            }
        }
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=HERE/"results")
    args = p.parse_args()
    start, cpu = perf_counter(), process_time()
    plan = read_json(args.output/"SELECTION.json")
    rows = []
    with np.load(args.output/"complex_readouts.npz", allow_pickle=False) as raw:
        for case in plan["view_manifest"]:
            ident = case["id"]
            geo = plan["geometries"][case["layout"]]
            syms = geo["symmetries"]
            axis = next((x for x in ("mirror_x", "mirror_y")
                         if syms[x]["status"] == "verified_closed_permutation"), None)
            row = {k: case[k] for k in ("id", "cohort", "layout", "object")}
            row["axis"] = axis
            y = raw[f"{ident}__Y_noisy"]
            if case["cohort"].startswith("SRAP"):
                prechosen = arb("0.05")*arb(y.size).sqrt()
                row["prospectively_fixed_radius_interval"] = rounded_ball(prechosen)
            if axis is not None:
                delta = arb(float(raw[f"{ident}__proof_delta_upper"]))
                row["archived_radius"] = scenario(y, syms[axis], delta)
                row["archived_radius"]["within_actual_contract"] = None
                if case["cohort"].startswith("SRAP"):
                    # Compliance is checked offline, after the stopping policy
                    # has been replayed with only its fixed prospective radius.
                    prefix = scenario(y, syms[axis], prechosen)
                    contract_check = bool(arb(float(raw[f"{ident}__joint_delta"])) <= prechosen)
                    prefix["within_actual_contract"] = contract_check
                    row["prospective_radius"] = prefix
            else:
                row["geometry_reason"] = {
                    x: syms[x]["status"] for x in ("mirror_x", "mirror_y")
                }
            if case["cohort"].startswith("SRAP"):
                row["archived_actual_noise_inside_prospective_radius"] = bool(
                    arb(float(raw[f"{ident}__joint_delta"])) <= prechosen)
            rows.append(row)
    result = {
        "contract": "fixed x reflection (otherwise y if x has no nontrivial closure); source orbits ordered by minimum original source ID",
        "prefix": "Only each newly closed source orbit enters every decision; preserve physical source and original source power, use all preassigned closed receivers.",
        "comparators": "full raw acquisition; fixed budget floor(raw source count/2), without splitting any source orbit; optional hindsight fixed prefix",
        "archived_radius": "original joint radius known conditionally from historical archive; Austria uses a conservative cap on the generator radius",
        "prospective_SRap_radius": "same deterministic absolute cap .05 sqrt(M_full), chosen before any Y prefix; actual archived noise checked after replay",
        "interpretation": "archive replay, no real hardware energy measurement or identification of general non-mirror ambiguities",
        "classical_identical_prefix_baseline": "sigma_min([Y, QY])/sqrt(2) = min(E, O) on every closed mirror prefix; ordering, stopping and exposure counts coincide",
        "counts": {
            "views": len(rows), "without_usable_geometry": sum(r["axis"] is None for r in rows),
            "SRAP_whole_archive_within_prospective_cap": sum(r.get("archived_actual_noise_inside_prospective_radius", False)
                                                          for r in rows),
            "SRAP_whole_archive_outside_prospective_cap": sum(r.get("archived_actual_noise_inside_prospective_radius") is False
                                                           for r in rows),
        },
        "summary": {"archived_radius": aggregate(rows, "archived_radius"),
                    "prospective_radius": aggregate(rows, "prospective_radius")},
        "rows": rows, "new_full_wave_rhs": 0, "new_hardware_exposures": 0,
        "wall_seconds": perf_counter()-start, "cpu_seconds": process_time()-cpu,
    }
    write(args.output/"exposure_prefix.json", result)
    for radius, stats in result["summary"].items():
        print(radius, stats.get("all-applicable"), flush=True)


if __name__ == "__main__":
    main()
