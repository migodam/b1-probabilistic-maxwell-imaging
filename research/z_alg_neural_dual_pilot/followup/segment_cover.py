"""Strict nested 3, 5, and 9-node covers of the fixed positive multiplier segment."""

import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
from statistics import median
import sys

sys.dont_write_bytecode = True

import numpy as np
from flint import acb, acb_mat, arb, ctx
from threadpoolctl import threadpool_limits

import jensen_port as jp


CACHE = jp.HERE / "results/segment_enclosures.json"
FROZEN = jp.HERE / "SEGMENT_COVER_FROZEN.json"
LEVELS = (3, 5, 9)
SUBDIVISIONS = 4
PROTOCOL = {
    "cap": jp.CAP,
    "anchor_indices": [0, 2],
    "levels": list(LEVELS),
    "new_node_order": [.25, .75, .125, .375, .625, .875],
    "subdivisions_per_interval": SUBDIVISIONS,
    "nodes": {str(n): [i / (n - 1) for i in range(n)] for n in LEVELS},
    "roof": "nearest outer chords; min of upper endpoint maxima on each piece, then max across pieces",
    "lower": "largest outward lower of the exact sampled native source heights",
    "verdicts": {"PASS": "sampled lower > radius",
                 "BELOW": "certified continuous upper < radius",
                 "UNRESOLVED": "neither strict inequality certified"},
    "comparators": "best of two endpoints, fixed midpoint, and best table at identical preparation",
    "precision_bits": jp.PRECISION,
    "threads": jp.THREADS,
    "development": "12 pre-existing independent validation objects at the fixed 2 percent view",
    "fresh": "32 independent follow-up objects with paired 1, 2, and 10 percent views",
    "historical": "28 exposed historical stress objects, 84 views including two anchor donors",
    "external8": "8 independent pre-existing thin-bridge and continuous shapes with newly simulated full-ring observations",
    "inputs": "Frozen noisy complex port readings and radius only. No scoring labels or material arrays.",
}


def freeze_protocol():
    if FROZEN.is_file():
        old = json.loads(FROZEN.read_text())
        if old["protocol"] != PROTOCOL:
            raise ValueError("Existing frozen segment-cover protocol differs from this code.")
    else:
        jp.atomic_json(FROZEN, {"frozen_utc": datetime.now(timezone.utc).isoformat(), "protocol": PROTOCOL})


class DyadicProfile(jp.SavedProfile):
    """Affine interpolation of represented normalized sector and optical entries."""

    def __init__(self, geometry, left, right, t):
        if left["mu"] != right["mu"] or t not in jp.GRID:
            raise ValueError("The two fixed endpoints must share the physical scale and dyadic grid.")
        weight = arb(t)

        def mix(name):
            return [(1 - weight) * arb(float(a)) + weight * arb(float(b))
                    for a, b in zip(np.asarray(left[name]).flat, np.asarray(right[name]).flat, strict=True)]

        self.exact_real = mix("sector")
        self.exact_optical = mix("optical")
        trial = {
            "mu": left["mu"],
            "sector": [float(value) for value in self.exact_real],
            "optical": [float(value) for value in self.exact_optical],
        }
        super().__init__(geometry, trial, jp.CAP)

    def exact_scaled_coefficients(self):
        z = [acb(arb(jp.CAP) * a, b) / 2 for a, b in zip(self.exact_real, self.exact_optical, strict=True)]
        return [self.exact_real[i] for i in self.g.blocks], [z[i] for i in self.g.blocks]


def setup_counts(base, extra, new_receipts):
    timing_keys = ("wall_seconds", "cpu_seconds")
    rows = {}
    for n in LEVELS:
        nodes = set(PROTOCOL["nodes"][str(n)])
        paid = [receipt for t, receipt in new_receipts.items() if t in nodes]
        rows[str(n)] = {"native_rhs": 30 * n, "new_rhs_relative_to_3": 30 * (n - 3)}
        for key in timing_keys:
            previous = sum(base["profile_work"][name][key] for name in ("anchor_0", "anchor_2", "native_midpoint"))
            propagation = sum(receipt[key] for receipt in paid)
            extra_geometry = extra[key] if len(paid) else 0.
            rows[str(n)][key] = base["geometry_setup"][key] + previous + propagation + extra_geometry
    return rows


def load_cache(path):
    if not path.is_file():
        return {}, None
    receipt = json.loads(path.read_text())
    if receipt["format"] != "native-exact-dyadic-segment-v1" or receipt["protocol"] != PROTOCOL:
        raise ValueError("The stored native response cache belongs to a different segment protocol.")
    models = {record["t"]: jp.PortModel(record["label"], jp.decode_complex_matrix(record["response"]),
                                        record["receipt"])
              for record in receipt["records"]}
    return models, receipt


def prepare(path=CACHE):
    """Prepare the six additional outward responses and reuse the frozen three."""
    freeze_protocol()
    ctx.prec, ctx.threads = jp.PRECISION, jp.THREADS
    path = Path(path)
    models, persisted = load_cache(path)
    if set(models) == set(jp.GRID):
        return models, persisted["setup"]
    if not jp.DEFAULT_CACHE.is_file():
        raise FileNotFoundError("The three prerequisite exact response enclosures are unavailable.")
    original = jp.prepare(jp.DEFAULT_CACHE)
    for weight, base in zip((0., .5, 1.), (original.anchors[0], original.midpoint, original.anchors[1]), strict=True):
        models.setdefault(weight, base)
    started = jp.started()
    profiles = jp.archived_anchors()
    with threadpool_limits(limits=1):
        gbegin = jp.started()
        geometry = jp.archived.Geometry(*jp.archived.native_geometry(), block_side=4)
        new_geometry = jp.elapsed(gbegin)
        geometry_ledger = new_geometry if persisted is None else {
            key: persisted["setup"]["additional_geometry_setup"][key] + new_geometry[key]
            for key in new_geometry
        }
        gap = min(p.archived_inertia["positive_gap"] for p in profiles)
        for t in PROTOCOL["new_node_order"]:
            if t in models:
                continue
            begin = jp.started()
            action = DyadicProfile(geometry, profiles[0].specification, profiles[1].specification, t)
            action_setup = jp.elapsed(begin)
            spec = {"anchor_indices": [0, 2], "t": t, "gap_method": "exact affine convex combination",
                    "coefficient_convention": "exact represented normalized products and dyadic interpolation"}
            model = jp.enclose_response(action, gap, f"native_t_{t:g}", spec)
            model.receipt["action_setup"] = action_setup
            model.receipt.update(jp.elapsed(begin))
            models[t] = model
            new = {q: v.receipt for q, v in models.items() if q not in (0., .5, 1.)}
            setup = {
                "predecessor_cache": jp.relative_receipt(jp.DEFAULT_CACHE),
                "additional_geometry_setup": geometry_ledger,
                "preparation_by_level": setup_counts(original.receipt, geometry_ledger, new),
                "node_work": {str(q): m.receipt for q, m in models.items()},
                "reused_initial_rhs": 90,
                "additional_rhs_completed": 30 * len(new),
                "new_physical_forward_rhs": 0,
                "latest_invocation_since_geometry_start": jp.elapsed(started),
            }
            payload = {"format": "native-exact-dyadic-segment-v1", "protocol": PROTOCOL, "setup": setup,
                       "records": [{"t": q, "label": model.label, "receipt": model.receipt,
                                    "response": jp.encode_complex_matrix(model.response)}
                                   for q, model in sorted(models.items())]}
            jp.atomic_json(path, payload)
            print(f"native t={t:g}, {model.receipt['wall_seconds']:.3f}s, {len(models)}/9", flush=True)
    return models, setup


def cover(nodes, spectra, heights, cache, subs=SUBDIVISIONS):
    """Return certified roof caps on every computational subinterval."""
    hit = 0

    def extrapolation(a, b, x):
        nonlocal hit
        if x == a or x == b:
            return spectra[x]["distance_upper"]
        if a < x < b:
            raise ArithmeticError("An interior chord gives a lower, not an upper.")
        key = a, b, x
        if key not in cache:
            z = (arb(x) - arb(a)) / (arb(b) - arb(a))
            roof = (1 - z) * heights[a] + z * heights[b]
            cache[key] = jp.source_bounds(roof, upper=True)["distance_upper"]
        else:
            hit += 1
        return cache[key]

    pieces = []
    for i, (a, b) in enumerate(zip(nodes[:-1], nodes[1:], strict=True)):
        lines = []
        if i > 0:
            lines.append((nodes[i - 1], a))
        if i + 2 < len(nodes):
            lines.append((b, nodes[i + 2]))
        for j in range(subs):
            low, high = a + j / subs * (b - a), a + (j + 1) / subs * (b - a)
            upper = min(max(extrapolation(*line, low), extrapolation(*line, high)) for line in lines)
            pieces.append([low, high, upper])
    return max(piece[2] for piece in pieces), pieces, hit


def evaluate(models, y, radius, *, object_id, noise):
    begin = jp.started()
    y, radius = np.asarray(y, dtype=np.complex128), float(radius)
    if y.shape != (jp.PORTS, jp.SOURCES) or not np.isfinite(y).all() or not np.isfinite(radius) or radius < 0:
        raise ValueError("Invalid measured 20×10 complex observation or radius.")
    states, spectra, cache, levels = {}, {}, {}, {}
    t0 = jp.started()
    yb = jp.ball(y)
    eye = acb_mat(jp.SOURCES, jp.SOURCES)
    for i in range(jp.SOURCES):
        eye[i, i] = 1
    xi = jp.original.vertical(yb, eye)
    xstar = jp.star(xi)
    ygram = jp.star(yb) * yb
    preprocessing = jp.elapsed(t0)
    for n in LEVELS:
        nodes = PROTOCOL["nodes"][str(n)]
        t0 = jp.started()
        for t in nodes:
            if t not in states:
                states[t] = ygram - xstar * models[t].response * xi
                spectra[t] = jp.source_bounds(states[t], upper=True)
        table_time = jp.elapsed(t0)
        t0 = jp.started()
        count_before = len(cache)
        raw_upper, pieces, repeats = cover(nodes, spectra, states, cache)
        roof_time = jp.elapsed(t0)
        best_t = max(nodes, key=lambda t: (spectra[t]["distance_lower"], -t))
        best_inner_t = max(nodes[1:-1], key=lambda t: (spectra[t]["distance_lower"], -t))
        endpoint_upper = max(spectra[0.]["distance_upper"], spectra[1.]["distance_upper"])
        best = spectra[best_t]["distance_lower"]
        upper = min([raw_upper] + [row["continuous_upper"] for row in levels.values()])
        if best > upper:
            raise ArithmeticError("A certified exact sample exceeds the continuous upper bound.")
        table_upper = max(spectra[t]["distance_upper"] for t in nodes)
        midpoint = spectra[.5]
        interior_lower = spectra[best_inner_t]["distance_lower"]
        verdict = "PASS" if best > radius else "BELOW" if upper < radius else "UNRESOLVED"
        levels[str(n)] = {
            "verdict": verdict, "best_t": best_t, "best_interior_t": best_inner_t,
            "best_interior_lower": interior_lower,
            "sampled_table": {"lower": best, "upper": table_upper, "passed": best > radius},
            "continuous_upper": upper, "continuous_upper_without_previous_levels": raw_upper,
            "optimality_gap_upper": jp.high(arb(upper) - arb(best)) if upper > best else 0.,
            "best_interior_strictly_above_endpoints": interior_lower > endpoint_upper,
            "best_interior_strictly_above_fixed_midpoint": interior_lower > midpoint["distance_upper"],
            "new_crossing_over_endpoints": best > radius and endpoint_upper <= radius,
            "new_crossing_over_midpoint": best > radius and midpoint["distance_upper"] <= radius,
            "pieces": pieces, "new_extrapolated_endpoint_proofs": len(cache) - count_before,
            "reused_extrapolated_endpoint_proofs": repeats,
            "incremental_table_time": table_time, "incremental_roof_time": roof_time,
        }
    for k in ("wall_seconds", "cpu_seconds"):
        table, roof = preprocessing[k], 0.
        for row in levels.values():
            table += row["incremental_table_time"][k]
            roof += row["incremental_roof_time"][k]
            row.setdefault("cumulative_table_time", {})[k] = table
            row.setdefault("cumulative_table_and_cover_time", {})[k] = table + roof
    return {"object_id": object_id, "noise_pct": noise, "radius": radius,
            "nodes": [{"t": t, **record} for t, record in sorted(spectra.items())],
            "best_endpoint_lower": max(spectra[t]["distance_lower"] for t in (0., 1.)),
            "best_endpoint_upper": max(spectra[t]["distance_upper"] for t in (0., 1.)),
            "fixed_midpoint_lower": spectra[.5]["distance_lower"],
            "fixed_midpoint_upper": spectra[.5]["distance_upper"],
            "levels": levels, "preprocessing": preprocessing, "evaluation": jp.elapsed(begin)}


def summarize(rows):
    groups = {"all_views": rows}
    for noise in sorted({row["noise_pct"] for row in rows}):
        groups[f"{noise:g}_percent"] = [row for row in rows if row["noise_pct"] == noise]
    output = {}
    for name, group in groups.items():
        output[name] = {}
        for n in LEVELS:
            key = str(n)
            levels = [row["levels"][key] for row in group]
            counts = Counter(level["verdict"] for level in levels)
            output[name][key] = {
                "conditions": len(group),
                "objects": len({r["object_id"] for r in group}),
                "verdicts": {state: counts[state] for state in ("PASS", "BELOW", "UNRESOLVED")},
                "fixed_endpoint_passes": sum(row["best_endpoint_lower"] > row["radius"] for row in group),
                "fixed_midpoint_passes": sum(row["fixed_midpoint_lower"] > row["radius"] for row in group),
                "median_optimum_gap_upper": median(row["optimality_gap_upper"] for row in levels),
                "median_table_wall_seconds": median(row["cumulative_table_time"]["wall_seconds"]
                                                    for row in levels),
                "median_table_and_cover_wall_seconds": median(row["cumulative_table_and_cover_time"]["wall_seconds"]
                                                              for row in levels),
            }
    named = {}
    for n in LEVELS:
        matches = []
        for row in rows:
            record = row["levels"][str(n)]
            flags = ("best_interior_strictly_above_endpoints", "best_interior_strictly_above_fixed_midpoint",
                     "new_crossing_over_endpoints", "new_crossing_over_midpoint")
            if not any(record[flag] for flag in flags):
                continue
            matches.append({
                "object_id": row["object_id"], "noise_pct": row["noise_pct"], "radius": row["radius"],
                "best_interior_t": record["best_interior_t"], "best_interior_lower": record["best_interior_lower"],
                "best_endpoint_upper": row["best_endpoint_upper"], "fixed_midpoint_upper": row["fixed_midpoint_upper"],
                **{flag: record[flag] for flag in flags},
            })
        named[str(n)] = matches
    return output, named


def run(cohort, path=None):
    initial = jp.started()
    models, setup = prepare()
    read_setup = jp.elapsed(initial)
    if cohort == "historical":
        views = jp.historical_observations()
        source = "src/holdout.py historical stress inputs"
    else:
        views = jp.read_observations(path, validation=cohort == "validation")
        source = jp.relative_receipt(path)
    expected = {"validation": (12, 12), "fresh": (32, 96), "historical": (28, 84), "external8": (8, 24)}
    if (len({row[0] for row in views}), len(views)) != expected[cohort]:
        raise ValueError(f"Frozen {cohort} pairing does not match the registered scope.")
    rows = []
    with threadpool_limits(limits=1):
        for oid, noise, y, radius in views:
            rows.append(evaluate(models, y, radius, object_id=oid, noise=noise))
    summary, named = summarize(rows)
    output = jp.HERE / "results" / f"{cohort}_segment_cover.json"
    jp.atomic_json(output, {
        "scope": PROTOCOL["development" if cohort == "validation" else cohort],
        "input": source, "cohort": cohort, "conditions": len(rows),
        "independent_objects": len({row["object_id"] for row in rows}),
        "protocol": jp.relative_receipt(FROZEN), "setup_cache": jp.relative_receipt(CACHE),
        "preparation_by_level": setup["preparation_by_level"], "setup_read": read_setup,
        "summary": summary, "named_improvements": named, "rows": rows, **jp.elapsed(initial),
    })
    print(f"{cohort}: {len(rows)} views -> {jp.relative_receipt(output)}", flush=True)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("cohort", choices=("prepare", "validation", "fresh", "historical", "external8"))
    cli.add_argument("--observations", type=Path)
    args = cli.parse_args()
    if args.cohort == "prepare":
        prepare()
        return
    defaults = {"validation": "synthetic_observations.npz", "fresh": "followup_observations.npz",
                "external8": "followup_external8.npz"}
    path = args.observations or (jp.PILOT / "data" / defaults[args.cohort] if args.cohort in defaults else None)
    run(args.cohort, path)


if __name__ == "__main__":
    main()
