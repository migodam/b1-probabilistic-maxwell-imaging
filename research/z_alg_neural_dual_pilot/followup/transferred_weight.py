"""Replay one frozen indefinite native weight and the archived fixed portfolios."""

import argparse
from collections import defaultdict
import json
from statistics import median
from types import SimpleNamespace
import sys

sys.dont_write_bytecode = True

import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

import jensen_port as jp


FROZEN = jp.HERE / "TRANSFERRED_WEIGHT_FROZEN.json"
RESULTS = jp.HERE / "results"
CACHE = RESULTS / "transferred_weight_common.json"
ADDED = "new_indefinite"
ANCHORS = ("fixed_16_0", "fixed_16_2")
PILOT = ("global_real", *ANCHORS)
FULL = ("global_real", "global_optical", "quadrant_optical", *ANCHORS)
PORTFOLIOS = {
    "added_weight_alone": (ADDED,), "positive_anchors": ANCHORS, "pilot_fixed": PILOT,
    "full_archived": FULL, "anchors_plus_added": (*ANCHORS, ADDED),
    "pilot_plus_added": (*PILOT, ADDED), "full_plus_added": (*FULL, ADDED),
}
TIMES = ("wall_seconds", "cpu_seconds")
ALIAS = {
    "global_real_075": "global_real", "global_optical_075": "global_optical",
    "quadrant_optical_075": "quadrant_optical", "fixed_4x4_0_075": "fixed_16_0",
    "fixed_4x4_2_075": "fixed_16_2",
}


def archived_row(path, identifier):
    with path.open() as source:
        found = [row for row in map(json.loads, source) if row.get("id") == identifier]
    if len(found) != 1:
        raise ValueError(f"Expected one row for {identifier} in {path.name}.")
    return found[0]


def read_freeze():
    frozen = json.loads(FROZEN.read_text())
    identifier, = frozen["dictionary_new_entries"]
    proposal = archived_row(jp.PILOT / frozen["source_ids"]["proposals"], identifier)
    individual = archived_row(jp.PILOT / frozen["source_ids"]["native"], identifier)
    native, physical, proof = individual["native"], frozen["physical_coefficients"], frozen["proof"]
    if (
        native["profile"] != physical or native["proof"] != proof or not individual["passed"]
        or individual["score_reuse_id"] is not None
        or (proposal["side"], proposal["u"], proposal["arm"]) != (8, jp.CAP, "gradient5")
        or [individual["lower"], individual["upper"]] != frozen["donor_original_interval"]
        or proposal["radius"] != individual["radius"] or individual["radius"] != frozen["donor_radius"]
        or physical["weights"] != [x / 9 for x in proposal["mu"]]
        or physical["optical"] != [x / 9 for x in proposal["lam"]]
        or physical["cap"] != jp.CAP or proof["nu"] != proof["nu_lower"] or proof["nu"] != 5
        or min(proof["positive_gap"], proof["negative_gap"]) <= 0.
        or proof["method"] != "native_Fourier_and_negative_subspace"
    ):
        raise ValueError("Frozen parameters or certified state gaps disagree with the admitted donor row.")
    return frozen


def profiles(frozen):
    archived = json.loads(jp.archived.ORIGINAL_ONLINE.read_text())
    if not archived["sealed"] or archived["native_dimension"] != 9216 or archived["joint_channel_shape"] != [20, 10]:
        raise ValueError("Original fixed comparator has a different geometry or is unsealed.")
    by_label = {entry["label"]: entry for entry in archived["profiles"] if entry["cap"] == jp.CAP}
    choices = {entry["label"]: entry for entry in jp.original.choices() if entry["cap"] == jp.CAP}
    if set(by_label) != set(choices) or set(choices) != set(ALIAS):
        raise ValueError("The complete archived fixed same-cap portfolio differs from the frozen set.")

    chosen = {
        "label": ADDED, "cap": jp.CAP, "block_side": 8, "nu": frozen["proof"]["nu"],
        "parameters": frozen["physical_coefficients"],
        "positive_gap_lower": frozen["proof"]["positive_gap"],
        "negative_gap_lower": frozen["proof"]["negative_gap"],
        "inherited_proof": frozen["source_ids"]["native"],
        "state_proof_work": frozen["original_individual_verification"],
        "coefficient_convention": "admitted native binary64 physical coefficients",
    }
    specs = {ADDED: chosen}
    for label, key in ALIAS.items():
        row, choice = by_label[label], choices[label]
        keys = ("label", "case", "block_side", "cap", "nu", "parameters", "positive_gap_lower")
        if any(row[name] != choice[name] for name in keys):
            raise ValueError(f"Archived parameters and current provenance for {label} disagree.")
        negative = row["negative_proof"].get("negative_gap_lower")
        if not row["positive_gap_lower"] > 0. or (row["nu"] and not (negative and negative > 0.)):
            raise ValueError(f"Uncertified spectral gap in {label}.")
        specs[key] = choice | {
            "archive_label": label, "label": key, "negative_gap_lower": negative,
            "archive_negative_proof": row["negative_proof"],
            "archive_previous_common_response": row["shared_stationary_work"],
            "coefficient_convention": "exact represented scale times archive normalized weights"
            if "profile_index" in choice else "archived binary64 physical coefficients",
        }
    return specs


def prepare(frozen, specs):
    start = jp.started()
    if CACHE.is_file():
        stored = json.loads(CACHE.read_text())
        if stored["profiles"] != specs or stored["dictionary_new_entries"] != frozen["dictionary_new_entries"]:
            raise ValueError("The cached state responses do not belong to the frozen profiles.")
    else:
        stored = {
            "format": "compact m and residual Gram; outward dyadic ACB intervals",
            "dictionary_new_entries": frozen["dictionary_new_entries"], "profiles": specs,
            "precision_bits": jp.PRECISION, "interval_threads": ctx.threads, "blas_threads": 1,
            "responses": {}, "setup": {"geometry_read": None, "partitions": {}},
        }
    responses = stored["responses"]
    if set(responses) != set(specs):
        begin = jp.started()
        arrays = jp.archived.native_geometry()
        stored["setup"]["geometry_read"] = jp.elapsed(begin)
        if arrays[0].shape != (192, 192) or arrays[1].shape != (20, 9216) or arrays[2].shape != (9216, 10):
            raise ValueError("Unexpected archive geometry dimensions.")
        geometries = {}
        for name, spec in specs.items():
            if name in responses:
                continue
            side = spec["block_side"]
            if side not in geometries:
                begin = jp.started()
                geometries[side] = jp.archived.Geometry(*arrays, block_side=side)
                stored["setup"]["partitions"][str(side)] = jp.elapsed(begin)
            geometry = geometries[side]
            begin = jp.started()
            costs = geometry.cost.copy()
            action = jp.original.saved_operator(geometry, spec)
            operator_setup = jp.elapsed(begin)
            m, defect, work = jp.original.exact_common_response(action, spec["nu"])
            responses[name] = {
                "m": jp.encode_complex_matrix(m), "defect": jp.encode_complex_matrix(defect),
                "stationary_work": work, "operator_setup": operator_setup,
                "operations": dict(geometry.cost - costs), **jp.elapsed(begin),
            }
            jp.atomic_json(CACHE, stored)
            print(name, "shared native work", round(work["wall_seconds"], 3),
                  "s, common RHS", work["shared_rhs"], flush=True)
    models = {}
    for name, spec in specs.items():
        response = responses[name]
        scale = spec["parameters"]["mu" if "profile_index" in spec else "scale"]
        models[name] = (
            SimpleNamespace(w=SimpleNamespace(scale=scale)), spec,
            jp.decode_complex_matrix(response["m"]), jp.decode_complex_matrix(response["defect"]),
        )
    return models, stored, jp.elapsed(start)


def read_transfer16():
    with np.load(jp.PILOT / "data/transfer16.npz", allow_pickle=False) as content:
        if set(content.files) != {"object_id", "Y", "radius", "noise_pct"}:
            raise ValueError("Transfer16 online inputs must be exactly IDs, noisy Y, radius and noise labels.")
        oids, y, radius, noise = (content[name] for name in ("object_id", "Y", "radius", "noise_pct"))
    if oids.shape != (16,) or y.shape != (16, 3, 20, 10) or radius.shape != (16, 3):
        raise ValueError("Transfer16 must contain 16 objects, 3 noise views, 20 receivers, and 10 sources.")
    if len(set(oids.tolist())) != 16 or not np.all(np.broadcast_to(noise, radius.shape) == [1, 2, 10]):
        raise ValueError("Transfer16 must have distinct objects with ordered paired 1%, 2%, and 10% views.")
    return [(str(oids[i]), (1, 2, 10)[j], y[i, j], radius[i, j]) for i, j in np.ndindex(radius.shape)]


def views(cohort):
    if cohort in ("donor", "historical"):
        archived = jp.historical_observations()
        if cohort == "donor":
            return [entry for entry in archived if entry[:2] == ("s2d_004", 2)]
        return archived
    if cohort == "transfer16":
        return read_transfer16()
    files = {"prior32": "followup_observations.npz", "external8": "followup_external8.npz"}
    return jp.read_observations(jp.PILOT / "data" / files[cohort])


def certify_view(view, models):
    oid, noise, y, radius = view
    radius = float(radius)
    if y.shape != (jp.PORTS, jp.SOURCES) or not np.isfinite(y).all() or not np.isfinite(radius) or radius < 0:
        raise ValueError(f"Invalid noisy port/radius for {oid}.")
    scores, arms = {}, {}
    for name, (operator, spec, m, defect) in models.items():
        start = jp.started()
        result = jp.original.bound_source(
            y, operator, spec["nu"], spec["positive_gap_lower"], spec["negative_gap_lower"],
            m, defect, dict(object_id=oid, noise_pct=noise, padded_radius=radius),
        )
        lo, hi = result["joint_distance_lower"], result["joint_distance_upper"]
        if not 0. <= lo <= hi or not np.isfinite(hi):
            raise ArithmeticError(f"Invalid native interval for {oid} {noise}%, {name}.")
        scores[name] = {"lower": lo, "upper": hi, "strictly_excludes": lo > radius,
                        "original_native_receipt": result, **jp.elapsed(start)}
    for name, members in PORTFOLIOS.items():
        winner = max(members, key=lambda item: scores[item]["lower"])
        lo, hi = scores[winner]["lower"], max(scores[item]["upper"] for item in members)
        arms[name] = {
            "lower": lo, "upper": hi, "best_profile": winner, "strictly_excludes": lo > radius,
            "verdict": "EXCLUDES" if lo > radius else "BELOW_OR_EQUAL" if hi <= radius else "UNRESOLVED",
            **{unit: sum(scores[item][unit] for item in members) for unit in TIMES},
        }
    return {"object_id": oid, "noise_pct": noise, "radius": radius, "profiles": scores, "arms": arms}


def identity(row):
    return [row["object_id"], row["noise_pct"]]


def comparison(rows, test, base):
    lost, gained, strict_better, strict_worse = [], [], [], []
    changes, relative = [], []
    for row in rows:
        trial, control = (row["arms"][arm] for arm in (test, base))
        case = identity(row)
        if trial["strictly_excludes"] and not control["strictly_excludes"]:
            gained.append(case)
        if control["strictly_excludes"] and not trial["strictly_excludes"]:
            lost.append(case)
        if trial["lower"] > control["upper"]:
            strict_better.append(case)
        if control["lower"] > trial["upper"]:
            strict_worse.append(case)
        changes.append(trial["lower"] - control["lower"])
        if control["lower"] > 0.:
            relative.append(trial["lower"] / control["lower"] - 1.)
    return {
        "gained_views": gained, "lost_views": lost,
        "objects_with_gained_view": sorted({key[0] for key in gained}),
        "objects_with_lost_view": sorted({key[0] for key in lost}),
        "objects_gained_any_view_exclusion": sorted(
            {r["object_id"] for r in rows if r["arms"][test]["strictly_excludes"]}
            - {r["object_id"] for r in rows if r["arms"][base]["strictly_excludes"]}),
        "objects_lost_any_view_exclusion": sorted(
            {r["object_id"] for r in rows if r["arms"][base]["strictly_excludes"]}
            - {r["object_id"] for r in rows if r["arms"][test]["strictly_excludes"]}),
        "strictly_larger_intervals": strict_better, "strictly_smaller_intervals": strict_worse,
        "median_lower_change": median(changes),
        "median_relative_lower_change_given_positive_baseline": median(relative) if relative else None,
        "views_with_positive_baseline": len(relative),
    }


def summarize(rows):
    arms = {}
    for name in PORTFOLIOS:
        values = [row["arms"][name] for row in rows]
        passing = [identity(row) for row, r in zip(rows, values, strict=True) if r["strictly_excludes"]]
        arms[name] = {
            "passing_views": passing, "passing_objects": sorted({case[0] for case in passing}),
            "views_by_verdict": {verdict: sum(r["verdict"] == verdict for r in values)
                                 for verdict in ("EXCLUDES", "BELOW_OR_EQUAL", "UNRESOLVED")},
            "median_online": {unit: median(r[unit] for r in values) for unit in TIMES},
            "total_online": {unit: sum(r[unit] for r in values) for unit in TIMES},
        }
    pairs = [(test, base) for test, base in (
        ("added_weight_alone", "positive_anchors"), ("added_weight_alone", "full_archived"),
        ("anchors_plus_added", "positive_anchors"), ("pilot_plus_added", "pilot_fixed"),
        ("full_plus_added", "full_archived"),
    )]
    return {"arms": arms, "comparisons": {f"{test}_vs_{base}": comparison(rows, test, base)
                                        for test, base in pairs}}


def costs(stored, frozen):
    results = {}
    selection = frozen["selection_and_freeze_metering"]
    proposal = frozen["original_proposal_cost"]
    proof = frozen["original_individual_verification"]
    fixed_geometry = stored["setup"]
    for name, members in PORTFOLIOS.items():
        responses = [stored["responses"][item] for item in members]
        geometries = {
            str(stored["profiles"][item]["block_side"]) for item in members
        }
        replay = {
            unit: fixed_geometry["geometry_read"][unit]
            + sum(fixed_geometry["partitions"][side][unit] for side in geometries)
            + sum(response[unit] for response in responses) for unit in TIMES
        }
        results[name] = {
            "native_common_rhs": sum(response["stationary_work"]["shared_rhs"] for response in responses),
            "matched_response_preparation": replay,
            "with_original_donor_proposal_verification_and_metered_freeze": {
                unit: replay[unit] + (proposal[unit] + proof[unit] + selection[unit] if ADDED in members else 0.)
                for unit in TIMES
            },
        }
    return results


def score(cohort, models, stored, frozen, cache_use):
    start = jp.started()
    own_views = views(cohort)
    n_objects, n_views = {"donor": (1, 1), "historical": (28, 84), "prior32": (32, 96),
                         "external8": (8, 24), "transfer16": (16, 48)}[cohort]
    if (len({row[0] for row in own_views}), len(own_views)) != (n_objects, n_views):
        raise ValueError(f"Unexpected frozen cohort size for {cohort}.")
    if len({row[:2] for row in own_views}) != n_views:
        raise ValueError(f"Repeated view IDs in {cohort}.")
    rows = [certify_view(view, models) for view in own_views]
    original_checks = []
    for row in rows:
        if identity(row) == ["s2d_004", 2]:
            if row["radius"] != frozen["donor_radius"]:
                raise ValueError("The transferred donor reload has a different radius.")
            original_lo, original_hi = frozen["donor_original_interval"]
            replay = row["profiles"][ADDED]
            if max(original_lo, replay["lower"]) > min(original_hi, replay["upper"]):
                raise ArithmeticError("Replay and original donor certificates disagree.")
            original_checks.append({"original": [original_lo, original_hi],
                                    "common_replay": [replay["lower"], replay["upper"]], "overlap": True})
    subsets = {"all_views": rows}
    subsets.update({f"{noise:g}_pct": [row for row in rows if row["noise_pct"] == noise]
                    for noise in sorted({row["noise_pct"] for row in rows})})
    output = {
        "cohort": cohort, "scope": frozen["cohorts"][cohort], "objects": n_objects, "views": n_views,
        "donor_is_included_in_historical_only_for_paired_accounting": cohort == "historical",
        "contract": frozen["physical_geometry"], "cap": frozen["cap"], "strict_threshold": "lower > radius",
        "frozen": jp.relative_receipt(FROZEN), "cache": jp.relative_receipt(CACHE),
        "profile_sets": PORTFOLIOS, "cache_load_or_prepare": cache_use,
        "preparation_cost_by_arm": costs(stored, frozen),
        "common_prep_cost_ledger": stored["setup"],
        "original_individual_identity_check": original_checks,
        "summary": {name: summarize(subset) for name, subset in subsets.items()}, "rows": rows,
        **jp.elapsed(start),
    }
    jp.atomic_json(RESULTS / f"transferred_weight_{cohort}.json", output)
    print(cohort, {name: len(record["passing_views"]) for name, record
                   in output["summary"]["all_views"]["arms"].items()}, flush=True)


def offline_audit():
    """Attach scalar generator flags only after all online result files exist."""
    profiles_and_cases = {}
    label_contracts = {
        "prior32": ("offline_scores.json", "generation.json"),
        "external8": ("external8_offline_scores.json", "external8_generation.json"),
        "transfer16": ("transfer16_offline.json", "transfer16_generation.json"),
    }
    frozen_online = {
        name: json.loads((RESULTS / f"transferred_weight_{name}.json").read_text())
        for name in ("donor", "historical", *label_contracts)
    }
    for name, (label_file, generator_file) in label_contracts.items():
        label_data = json.loads((jp.HERE / label_file).read_text())
        labels = label_data["rows"] if isinstance(label_data, dict) else label_data
        label_by_id = {row.get("object_id", row.get("oid")): row for row in labels}
        generator = json.loads((jp.HERE / generator_file).read_text())
        generator_by_id = {row.get("object_id", row.get("oid")): row for row in generator["rows"]}
        rows = frozen_online[name]["rows"]
        if set(label_by_id) != set(generator_by_id) or set(label_by_id) != {r["object_id"] for r in rows}:
            raise ValueError(f"Offline IDs disagree with scored inputs for {name}.")
        levels = generator.get("noise_pct", [1, 2, 10])
        checks, high, low = [], [], []
        for row in rows:
            oid, noise = identity(row)
            label, audit = label_by_id[oid], generator_by_id[oid]
            under = label["under_cap"]
            if under != (label["peak"] <= jp.CAP):
                raise ValueError(f"Offline cap label disagrees with reported scalar peak for {oid}.")
            index = levels.index(noise)
            upper = audit["truth_distance_upper"][index] if "truth_distance_upper" in audit else None
            slack = audit["radius_slack_lower"][index]
            if slack <= 0. or (upper is not None and upper > row["radius"]):
                raise ValueError(f"The recorded true object's fit is unqualified for {oid}, {noise}%.")
            passed = [key for key, v in row["profiles"].items() if v["strictly_excludes"]]
            if under and passed:
                raise ArithmeticError(f"A below-cap in-ball control was incorrectly excluded for {oid}, {noise}%.")
            checks.append({
                "object_id": oid, "noise_pct": noise, "under_cap": under, "peak": label["peak"],
                "radius": row["radius"], "truth_distance_upper": upper, "truth_radius_slack_lower": slack,
                "passing_physical_profiles": passed,
            })
            (low if under else high).append(row)
        profiles_and_cases[name] = {
            "scalar_labels": label_file, "generator_receipt": generator_file,
            "low_control_objects": sorted({r["object_id"] for r in low}),
            "low_control_views": len(low), "high_objects": sorted({r["object_id"] for r in high}),
            "high_views": len(high), "minimum_truth_radius_slack_lower": min(
                rec["truth_radius_slack_lower"] for rec in checks),
            "low_control_results": summarize(low) if low else None,
            "above_cap_results": summarize(high) if high else None, "view_checks": checks,
        }
    old = frozen_online["historical"]["rows"]
    historic = {(r["object_id"], r["noise_pct"]): r for r in old}
    original = json.loads(jp.archived.ORIGINAL_ONLINE.read_text())
    original_checks = []
    for profile in original["profiles"]:
        label = profile["label"]
        if label not in ALIAS:
            continue
        for expected in profile["rows"]:
            case_id = expected["object_id"], 10
            row = historic[case_id]
            interval = [expected["joint_distance_lower"], expected["joint_distance_upper"]]
            candidate = row["profiles"][ALIAS[label]]
            if row["radius"] != expected["padded_radius"] or max(interval[0], candidate["lower"]) > min(
                interval[1], candidate["upper"]
            ):
                raise ArithmeticError(f"Archival proof replay disagrees for {label}, {case_id}.")
            original_checks.append({"case": list(case_id), "profile": ALIAS[label],
                                    "archived_interval": interval,
                                    "replay_interval": [candidate["lower"], candidate["upper"]]})
    no_donor = [r for r in old if r["object_id"] != "s2d_004"]
    jp.atomic_json(RESULTS / "transferred_weight_offline_audit.json", {
        "scope": "separate scalar label and same-radius archive audit after scoring",
        "control_cohorts": profiles_and_cases,
        "archived_srap10_interval_checks": original_checks,
        "historical_without_any_donor_material_view": {
            "objects": len({r["object_id"] for r in no_donor}), "views": len(no_donor),
            "summary": summarize(no_donor),
        },
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    retrospective = ("donor", "historical", "prior32", "external8")
    parser.add_argument("cohort", choices=("prepare", "all", "retrospective", "transfer16", "offline_audit", *retrospective))
    choice = parser.parse_args().cohort
    if choice == "offline_audit":
        offline_audit()
        return
    ctx.prec, ctx.threads = jp.PRECISION, 1
    frozen = read_freeze()
    specs = profiles(frozen)
    with threadpool_limits(limits=1):
        models, stored, load_cost = prepare(frozen, specs)
        if choice == "prepare":
            return
        selected = {"all": ("transfer16", *retrospective), "retrospective": retrospective}
        for cohort in selected.get(choice, (choice,)):
            score(cohort, models, stored, frozen, load_cost)


if __name__ == "__main__":
    main()
