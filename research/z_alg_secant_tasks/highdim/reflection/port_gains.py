"""Bound mirror-pair ambiguity under independent receive or transmit gains."""

import argparse
from collections import defaultdict
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from flint import arb

from replay_symmetry import HERE, channels, lift_complex, read_json, rounded_ball, squared_length, write


def block_lambda(first, second):
    """Smallest Gram eigenvalue of two complex column vectors."""
    if len(first) <= 1:
        return arb(0)
    a, b = lift_complex(first), lift_complex(second)
    p, q = squared_length(a), squared_length(b)
    c = sum(x.conjugate()*y for x, y in zip(a, b, strict=True))
    d = p-q
    discriminant = (d*d+4*(c.real*c.real+c.imag*c.imag)).nonnegative_part().sqrt()
    return ((p+q-discriminant)/2).nonnegative_part()


def block_svd(first, second):
    if len(first) <= 1:
        return 0.
    s = np.linalg.svd(np.column_stack((first, second)), compute_uv=False)
    return float(s[-1]**2)


def losses(y, reflected):
    """All block sums; the independent SVD is an ordinary numerical control."""
    cases = {"common_scalar": (y.reshape(1, -1), reflected.reshape(1, -1)),
             "receive_independent": (y, reflected),
             "transmit_independent": (y.T, reflected.T)}
    result = {}
    for name, (f, z) in cases.items():
        distance = (sum(block_lambda(a, b) for a, b in zip(f, z, strict=True))/2).sqrt()
        reference = np.sqrt(sum(block_svd(a, b) for a, b in zip(f, z, strict=True))/2)
        if abs(float(distance)-reference) > 1e-10 * max(1., np.linalg.norm(y)):
            raise ArithmeticError("Arb and ordinary rank-one replay disagree.")
        result[name] = distance, float(reference)
    for name in ("receive_independent", "transmit_independent"):
        if result[name][0] > result["common_scalar"][0]:
            raise ArithmeticError("Free port gains cannot strengthen the common scalar certificate.")
    return result


def outcomes(values, radius, inside):
    return {
        name: {
            "D_interval": [0., 0.] if value == 0 else rounded_ball(value),
            "classical_block_svd_D": numeric,
            "D_over_delta_interval": [0., 0.] if value == 0 else rounded_ball(value/radius),
            "D_minus_delta_interval": rounded_ball(value-radius),
            "threshold_exceeded": bool(value > radius),
            "certified_with_noise_contract": bool(value > radius and inside is True),
        }
        for name, (value, numeric) in values.items()
    }


def inside_interval(noise, radius):
    if noise <= radius:
        return True
    if noise > radius:
        return False
    return None


def noise_bound(case, y, raw, nominal, proof):
    if case["cohort"] == "R10":
        clean = raw[f"{case['id']}__saved_clean"]
        exact_res = [x-z for x, z in zip(lift_complex(y), lift_complex(clean), strict=True)]
        return squared_length(exact_res).sqrt(), "saved original noisy and clean observations"
    if case["cohort"] == "Austria":
        signal = squared_length(lift_complex(y)).sqrt()
        low = signal/(arb(11)+arb("1e-11"))
        return low.lower().union(proof), "original 10% full-matrix noise law; triangle bounds"
    return nominal, "actual archived SRAP joint noise norm"


def statistic_for(y, plan, budgets):
    a, b = channels(y, plan), channels(y, plan, partner=True)
    values = losses(a, b)
    table = {name: outcomes(values, radius, inside) for name, (radius, inside) in budgets.items()}
    return {"measured_ports": int(a.size), "all_original_ports": int(y.size),
            "closed_sources": int(a.shape[1]), "closed_receivers": int(a.shape[0]), "budgets": table}


def summarize(rows, choices):
    gathered = defaultdict(list)
    for row in rows:
        for scope, field in choices:
            if field not in row:
                continue
            gathered[(scope, row["cohort"])].append(row[field])
            if row["cohort"].startswith("SRAP/"):
                gathered[(scope, "SRAP_all_eligible")].append(row[field])
    result = []
    for (scope, cohort), group in sorted(gathered.items()):
        for radius in ("archived", "prospective"):
            trials = [trial for trial in group if radius in trial["budgets"]]
            if not trials:
                continue
            verified = [t for t in trials if t["noise_within_budget"][radius] is True]
            counts = {}
            for control in ("common_scalar", "receive_independent", "transmit_independent"):
                conditional = [t["budgets"][radius][control] for t in trials]
                pass_set = [t for t in conditional if t["certified_with_noise_contract"]]
                counts[control] = {
                    "ratio_range": [min(x["D_over_delta_interval"][0] for x in conditional),
                                    max(x["D_over_delta_interval"][1] for x in conditional)],
                    "threshold_exceeds": sum(x["threshold_exceeded"] for x in conditional),
                    "certified": len(pass_set),
                    "within_contract_without_certificate": len(verified)-len(pass_set),
                    "apparent_exceeds_outside_or_uncertain_contract": sum(
                        x["threshold_exceeded"] and not x["certified_with_noise_contract"] for x in conditional),
                    "minimum_certified_margin": min((x["D_minus_delta_interval"][0] for x in pass_set),
                                                    default=None),
                }
            result.append({"use": scope, "cohort": cohort, "radius": radius, "tested_views_or_turns": len(trials),
                           "within_budget": len(verified), "outside_budget": sum(
                               t["noise_within_budget"][radius] is False for t in trials),
                           "uncertain_budget": sum(t["noise_within_budget"][radius] is None for t in trials),
                           "controls": counts})
    return result


def run(directory):
    start, cpu = perf_counter(), process_time()
    plan = read_json(directory/"SELECTION.json")
    inputfile = directory/"complex_readouts.npz"
    report = []
    rotation_rows = []
    with np.load(inputfile, allow_pickle=False) as raw:
        for case in plan["view_manifest"]:
            case_id = case["id"]
            geometry = plan["geometries"][case["layout"]]
            y = raw[f"{case_id}__Y_noisy"]
            delta = arb(float(raw[f"{case_id}__joint_delta"]))
            conservative = arb(float(raw[f"{case_id}__proof_delta_upper"]))
            prospective = arb(".05")*arb(y.size).sqrt()
            noise, description = noise_bound(case, y, raw, delta, conservative)
            actual = {name: inside_interval(noise, cap) for name, cap in (
                ("archived", conservative), ("prospective", prospective))}
            if case["cohort"] == "Austria":
                actual["archived"] = True  # The conservative radius is the validated envelope.
            budgets = {"archived": (conservative, actual["archived"]),
                       "prospective": (prospective, actual["prospective"])}
            record = {key: case[key] for key in ("id", "object", "cohort", "layout")}
            record.update({
                "archive_nominal_delta": float(delta), "archive_proof_delta": float(conservative),
                "prospective_radius_interval": rounded_ball(prospective),
                "actual_joint_noise_interval": rounded_ball(noise),
                "noise_definition": description,
                "primary_axis": geometry["largest_verified_mirror"],
                "single_source_axis": geometry["primary_single_source_mirror"],
                "geometry_status": {}, "all_qualified_actions": {},
            })
            for action, band in geometry["symmetries"].items():
                record["geometry_status"][action] = band["status"]
                if band["status"] != "verified_closed_permutation":
                    continue
                measured = statistic_for(y, band, budgets)
                measured["noise_within_budget"] = actual
                record["all_qualified_actions"][action] = measured
                if action == record["primary_axis"]:
                    record["primary_mirror"] = measured
                if action.startswith("turn_"):
                    rotation_rows.append({"cohort": case["cohort"], action: measured})
            if record["single_source_axis"] is not None:
                d = geometry["symmetries"][record["single_source_axis"]]
                source = d["single_source"]["source"]
                indexed = np.ix_(d["receiver_order"], [source])
                reflected = np.ix_(d["receiver_image_order"], [source])
                values = losses(y[indexed], y[reflected])
                assert values["receive_independent"][0] == 0
                single = {name: outcomes(values, radius, valid) for name, (radius, valid) in budgets.items()}
                record["single_source"] = {
                    "source_index": source, "retained_receivers": len(d["receiver_order"]),
                    "noise_within_budget": actual, "budgets": single,
                }
            report.append(record)
    counts = summarize(report, (("primary_mirror", "primary_mirror"), ("fixed_single_source", "single_source")))
    counts.extend(summarize(rotation_rows, tuple((name, name) for name in ("turn_90", "turn_180", "turn_270"))))
    write(directory/"port_gains.json", {
        "input_scope": "Frozen archived observations, selected transforms and noise radii only. No material labels or new field solves.",
        "threshold": "D > joint delta on the reflection-closed product mask, using original full noise contract after crop",
        "per_port_controls": "Independent complex gains on the named axis only, separately per candidate material. The other port axis stays coherent.",
        "single_source_limit": "The receive-independent loss is identically zero for one source. The transmit-independent bound becomes the common-scalar row bound.",
        "physical_budget": "Prospective delta=.05*sqrt(M_full); numerically re-evaluated outside SRAP and annotated separately.",
        "proof": "128-bit Arb Gram eigenvalues and separate ordinary SVD replay on every retained block",
        "count": len(report), "total_new_physical_rhs": 0, "material_arrays_read": 0,
        "wall_seconds": perf_counter()-start, "cpu_seconds": process_time()-cpu,
        "summaries": counts, "rows": report,
    })
    for row in counts:
        print(row["use"], row["cohort"], row["radius"], row["within_budget"],
              {name: control["certified"] for name, control in row["controls"].items()}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE/"results")
    run(parser.parse_args().output)
