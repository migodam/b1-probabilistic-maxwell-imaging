"""Cover the fixed positive multiplier segment using a native Galerkin trial space."""

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import json
from statistics import median
import sys

sys.dont_write_bytecode = True

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import solve, svd
from scipy.sparse.linalg import LinearOperator, cg
from threadpoolctl import threadpool_limits

import jensen_port as jp
import segment_cover as existing


RESULTS = jp.HERE / "results"
FROZEN = jp.HERE / "GALERKIN_FROZEN.json"
CACHE = RESULTS / "galerkin_enclosures.json.gz"
LEVELS = (4, 16, 64)
KEY_VIEW = ("dualpilot_validation_000", 1.)
HISTORICAL_VIEW = ("s2d_004", 2.)
PROTOCOL = {
    "cap": jp.CAP, "anchor_indices": [0, 2], "partitions": list(LEVELS),
    "native_cg": {"rtol": 3e-13, "maxiter": 700, "rhs_per_endpoint": jp.SIZE},
    "basis": "SVD of separately normalized endpoint trial columns, cutoff 1e-11 relative to largest",
    "upper": "one center Galerkin trial per equal dyadic subinterval, reused at both edges; max of convex positive trace at edges",
    "nesting": "min of the current raw piece upper and its already certified parent upper",
    "lower": "maximum of the two previously paid strict endpoint lowers, without any interior node",
    "verdicts": {"PASS": "endpoint lower > radius", "BELOW": "continuous upper < radius",
                 "UNRESOLVED": "otherwise"},
    "precision_bits": jp.PRECISION, "arb_threads": 1, "blas_threads": 1,
    "validation_focus": list(KEY_VIEW), "historical_focus": list(HISTORICAL_VIEW),
    "cohorts": {"validation_comparator": "twelve validation objects at 2%",
                "fresh": "all 32 objects, 96 views", "external8": "all 8 objects, 24 views"},
    "comparison": "previous same-view 3/5/9 results; focus validation 1% replay uses saved matrices only",
}


def freeze():
    if FROZEN.is_file():
        if json.loads(FROZEN.read_text())["protocol"] != PROTOCOL:
            raise ValueError("The Galerkin protocol differs from its frozen record.")
    else:
        jp.atomic_json(FROZEN, {"frozen_utc": datetime.now(timezone.utc).isoformat(), "protocol": PROTOCOL})


def endpoint_trials(action):
    """Recreate the archived numerical common solves and retain the trial columns."""
    start = jp.started()
    geom, scale = action.g, action.w.scale
    rhs = np.column_stack((geom.sensor_star / scale, action.z[:, None] * geom.incident))
    preconditioner = LinearOperator((geom.n, geom.n), matvec=action.precondition, dtype=complex)
    columns, iterations = [], []
    for index, vector in enumerate(rhs.T):
        steps = []
        trial, status = cg(action.operator, vector, M=preconditioner, atol=0.,
                           rtol=PROTOCOL["native_cg"]["rtol"],
                           maxiter=PROTOCOL["native_cg"]["maxiter"],
                           callback=lambda _: steps.append(1))
        if status:
            raise ArithmeticError(f"Native common CG column {index} failed with status {status}.")
        columns.append(trial)
        iterations.append(len(steps))
    return np.column_stack(columns), dict(
        rhs=len(columns), iterations=iterations, **jp.elapsed(start),
    )


def compressed_endpoints(geometry, basis, actions):
    """Enclose V*AV and V*(T/alpha) from one set of exact Green actions."""
    trial = jp.ball(basis)
    propagated = geometry.exact_green(trial)
    observed = geometry.sensor_ball * trial
    observed_star = jp.star(observed)
    scale = arb(actions[0].w.scale)
    common = observed_star * observed / scale
    multiply_rows = jp.original.mult_rows
    matrices = []
    for action in actions:
        mu, z = action.exact_scaled_coefficients()
        z_star_v = multiply_rows(trial, [entry.conjugate() for entry in z])
        optical = jp.star(z_star_v) * propagated
        compressed = common + jp.star(trial) * multiply_rows(trial, mu) - optical - jp.star(optical)
        forcing = observed_star / scale
        incident = jp.star(z_star_v) * geometry.incident_ball
        compact_rhs = acb_mat(trial.ncols(), jp.SIZE, [
            value for i in range(trial.ncols()) for value in (
                [forcing[i, j] for j in range(jp.PORTS)]
                + [incident[i, j] for j in range(jp.SOURCES)]
            )
        ])
        matrices.append((compressed, compact_rhs))
    return matrices


def mix(left, right, t):
    return (1 - t) * left + t * right


def roofs(matrices, scale):
    """Construct two compact affine trial forms for every dyadic subinterval."""
    (a0, b0), (a1, b1) = matrices
    amid = [jp.hermitian(jp.mid(x)) for x in (a0, a1)]
    bmid = [jp.mid(x) for x in (b0, b1)]
    levels, timings = {}, {}
    for size in LEVELS:
        start = jp.started()
        pieces = []
        for index in range(size):
            center = (index + .5) / size
            coefficient = solve(mix(*amid, center), mix(*bmid, center), assume_a="pos",
                                check_finite=False)
            c = jp.ball(coefficient)
            cstar = jp.star(c)
            q = []
            for a, b in matrices:
                cross = cstar * b
                q.append(scale * (cstar * a * c - cross - jp.star(cross)))
            a, b = arb(index) / size, arb(index + 1) / size
            pieces.append((mix(*q, a), mix(*q, b)))
        levels[size] = pieces
        timings[str(size)] = dict(projected_solves=size, compact_trial_forms=2 * size, **jp.elapsed(start))
    return levels, timings


def encode(matrix):
    return {"shape": [matrix.nrows(), matrix.ncols()], "entries": jp.encode_complex_matrix(matrix)}


def decode(record):
    return acb_mat(*record["shape"], [
        acb(jp.decode_real(real), jp.decode_real(imag)) for real, imag in record["entries"]
    ])


def prepare():
    profiles = jp.archived_anchors()
    if CACHE.is_file():
        start = jp.started()
        with gzip.open(CACHE, "rt") as stream:
            raw = json.load(stream)
        if raw["protocol"] != PROTOCOL or raw["profiles"] != [p.specification for p in profiles]:
            raise ValueError("The Galerkin cache differs from the frozen protocol or archived anchors.")
        levels = {int(key): [(decode(a), decode(b)) for a, b in values]
                  for key, values in raw["trial_forms"].items()}
        return levels, {**raw["setup"], "read_compact_cache": jp.elapsed(start)}

    begin = jp.started()
    if profiles[0].specification["mu"] != profiles[1].specification["mu"]:
        raise ValueError("The anchor operators must have a common exact scale.")
    start = jp.started()
    geometry = jp.archived.Geometry(*jp.archived.native_geometry(), block_side=4)
    actions = [profile.action(geometry) for profile in profiles]
    load_geometry = jp.elapsed(start)
    columns, native_solves = [], []
    for profile, action in zip(profiles, actions, strict=True):
        vectors, work = endpoint_trials(action)
        columns.append(vectors)
        native_solves.append(dict(index=profile.index, **work))
        print(f"Anchor {profile.index} trial solves: {work['wall_seconds']:.3f} s", flush=True)
    start = jp.started()
    stacked = np.column_stack(columns)
    stacked /= np.linalg.norm(stacked, axis=0)
    modes, values, _ = svd(stacked, full_matrices=False, check_finite=False)
    keep = values > 1.e-11 * values[0]
    basis = np.ascontiguousarray(modes[:, keep])
    build_basis = dict(rank=basis.shape[1], singular_values=values.tolist(), **jp.elapsed(start))
    start = jp.started()
    matrices = compressed_endpoints(geometry, basis, actions)
    compression = jp.elapsed(start)
    levels, compact = roofs(matrices, arb(actions[0].w.scale))
    setup = {
        "load_geometry": load_geometry, "native_solves": native_solves,
        "build_basis": build_basis, "compression": compression, "compact_levels": compact,
        "native_rhs": sum(entry["rhs"] for entry in native_solves),
        "numerical_green_columns": geometry.cost["native_binary64_green_columns"],
        "interval_green_columns": geometry.cost["native_arb_green_columns"],
        **jp.elapsed(begin),
    }
    start = jp.started()
    payload = {
        "protocol": PROTOCOL, "profiles": [p.specification for p in profiles], "setup": setup,
        "native_gap_proofs": [p.archived_inertia for p in profiles],
        "projected_endpoint_K_scaled_and_rhs_scaled": [
            [encode(a), encode(b)] for a, b in matrices
        ],
        "trial_forms": {str(n): [[encode(a), encode(b)] for a, b in values]
                        for n, values in levels.items()},
    }
    temporary = CACHE.with_suffix(".tmp")
    with gzip.open(temporary, "wt", compresslevel=3) as stream:
        json.dump(payload, stream, separators=(",", ":"), allow_nan=False)
    temporary.replace(CACHE)
    setup["write_compact_cache"] = jp.elapsed(start)
    print(f"Geometry and reduced forms built: rank {basis.shape[1]}, {setup['wall_seconds']:.3f} s", flush=True)
    return levels, setup


def paid_endpoint_models():
    raw = json.loads(jp.DEFAULT_CACHE.read_text())
    profiles = jp.archived_anchors()
    if (raw["format"], raw["cap"], raw["precision_bits"], raw["portfolio_order"]) != (
        "dyadic-arb-balls-v1", jp.CAP, jp.PRECISION, [p.index for p in profiles]
    ):
        raise ValueError("The previously paid endpoint source enclosures differ from this trial.")
    endpoints = []
    for profile, entry in zip(profiles, raw["records"][:2], strict=True):
        if entry["receipt"]["profile"]["specification"] != profile.specification:
            raise ValueError("Stored endpoint source enclosure has a different physical profile.")
        endpoints.append(jp.decode_complex_matrix(entry["response"]))
    setup = raw["setup"]
    paid = setup.get("reusable_anchor_setup")
    if paid is None:
        paid = {
            unit: setup["geometry_setup"][unit] + sum(record["receipt"][unit] for record in raw["records"][:2])
            for unit in ("wall_seconds", "cpu_seconds")
        }
    return endpoints, paid


def online_upper(xistar, xi, gram, form):
    return jp.source_bounds(gram + xistar * form * xi, upper=True)["distance_upper"]


def evaluate(levels, paid_endpoints, observation):
    oid, noise, y, radius = observation
    radius = float(radius)
    begin = jp.started()
    y = np.asarray(y, dtype=np.complex128)
    if y.shape != (jp.PORTS, jp.SOURCES) or not np.isfinite(y).all() or not np.isfinite(radius) or radius < 0.:
        raise ValueError(f"Bad fixed observation or radius for {oid}.")
    yb = jp.ball(y)
    eye = acb_mat(jp.SOURCES, jp.SOURCES)
    for i in range(jp.SOURCES):
        eye[i, i] = 1
    xi = jp.original.vertical(yb, eye)
    xistar = jp.star(xi)
    gram = jp.star(yb) * yb
    preprocess = jp.elapsed(begin)
    start = jp.started()
    endpoint_bounds = [jp.source_bounds(gram - xistar * q * xi, upper=True) for q in paid_endpoints]
    lower = max(proof["distance_lower"] for proof in endpoint_bounds)
    endpoint_timing = jp.elapsed(start)
    rows, previous = {}, None
    for size, pieces in levels.items():
        start = jp.started()
        raw_bounds = [[online_upper(xistar, xi, gram, form) for form in pair] for pair in pieces]
        piece_bounds = [max(pair) for pair in raw_bounds]
        if previous is not None:
            ratio = size // len(previous)
            piece_bounds = [min(value, previous[i // ratio]) for i, value in enumerate(piece_bounds)]
        upper = max(piece_bounds)
        if lower > upper:
            raise ArithmeticError(f"Strict endpoint lower exceeds the trial upper for {oid} at {noise} percent.")
        rows[str(size)] = {
            "verdict": "PASS" if lower > radius else "BELOW" if upper < radius else "UNRESOLVED",
            "lower": lower, "upper": upper, "raw_upper": max(max(pair) for pair in raw_bounds),
            "gap_upper": jp.high(arb(upper) - arb(lower)) if upper > lower else 0.,
            "piece_bound": piece_bounds, "uncombined_edge_uppers": raw_bounds,
            "new_port_contractions": 2 * size, **jp.elapsed(start),
        }
        previous = piece_bounds
    return {
        "object_id": oid, "noise_pct": noise, "radius": radius, "endpoint_bounds": endpoint_bounds,
        "levels": rows, "preprocessing": preprocess, "contract_paid_endpoints": endpoint_timing,
        **jp.elapsed(begin),
    }


def focus_validation():
    with np.load(jp.PILOT / "data/synthetic_observations.npz", allow_pickle=False) as data:
        index = data["object_id"].tolist().index(KEY_VIEW[0])
        noise = data["noise_fraction"].tolist().index(KEY_VIEW[1] / 100)
        if data["split"][index] != "validation":
            raise ValueError("The focus observation is outside its declared validation split.")
        return (*KEY_VIEW, data["Y"][index, noise], data["radius"][index, noise])


def read_views(cohort):
    if cohort == "validation_focus":
        return [focus_validation()]
    if cohort == "historical_focus":
        return [row for row in jp.historical_observations() if row[:2] == HISTORICAL_VIEW]
    if cohort == "historical":
        return jp.historical_observations()
    files = {"validation": "synthetic_observations.npz", "fresh": "followup_observations.npz",
             "external8": "followup_external8.npz"}
    return jp.read_observations(jp.PILOT / "data" / files[cohort], validation=cohort == "validation")


def existing_controls(cohort, views):
    if cohort == "validation_focus":
        models, receipt = existing.load_cache(existing.CACHE)
        if set(models) != set(jp.GRID):
            raise ValueError("The reference 9-node enclosures have not been prepared.")
        oid, noise, y, radius = views[0]
        row = existing.evaluate(models, y, radius, object_id=oid, noise=noise)
        return [row], receipt["setup"]["preparation_by_level"], "archived node enclosures replayed"
    archived_cohort = "historical" if cohort == "historical_focus" else cohort
    path = RESULTS / f"{archived_cohort}_segment_cover.json"
    old = json.loads(path.read_text())
    return old["rows"], old["preparation_by_level"], jp.relative_receipt(path)


def attach_controls(rows, archived_rows):
    by_key = {(r["object_id"], r["noise_pct"]): r for r in archived_rows}
    for row in rows:
        prior = by_key[row["object_id"], row["noise_pct"]]
        if row["radius"] != prior["radius"]:
            raise ValueError(f"The control uses a different radius for {row['object_id']}.")
        row["sampled_native_controls"] = {}
        for n, control in prior["levels"].items():
            old_lower, old_upper = control["sampled_table"]["lower"], control["continuous_upper"]
            if any(max(r["lower"], old_lower) > min(r["upper"], old_upper)
                   for r in row["levels"].values()):
                raise ArithmeticError(f"Galerkin and node-bound intervals do not overlap for {row['object_id']}.")
            row["sampled_native_controls"][n] = {
                "lower": old_lower, "upper": old_upper, "verdict": control["verdict"],
                "cumulative_table_time": control["cumulative_table_time"],
                "cumulative_table_and_cover_time": control["cumulative_table_and_cover_time"],
                "intersects_all_galerkin_levels": True,
            }


def summarize(rows, field, levels, proposed):
    out = {}
    for level in map(str, levels):
        level_rows = [row[field][level] for row in rows]
        levels_used = [n for n in levels if n <= int(level)]
        if proposed:
            keys = ("preprocessing", "contract_paid_endpoints")
            costs = {unit: median(
                sum(row[key][unit] for key in keys)
                + sum(row[field][str(n)][unit] for n in levels_used)
                for row in rows
            ) for unit in ("wall_seconds", "cpu_seconds")}
        else:
            costs = {unit: median(record["cumulative_table_and_cover_time"][unit] for record in level_rows)
                     for unit in ("wall_seconds", "cpu_seconds")}
        counts = Counter(record["verdict"] for record in level_rows)
        out[level] = {
            "verdict_counts": {key: counts[key] for key in ("PASS", "BELOW", "UNRESOLVED")},
            "median_continuous_upper": median(record["upper"] for record in level_rows),
            "median_cumulative_online_time": costs,
        }
    return out


def run_cohort(cohort, levels, setup, endpoint_models, paid_cost):
    begin = jp.started()
    views = read_views(cohort)
    expected = {"validation_focus": 1, "historical_focus": 1, "historical": 84, "validation": 12,
                "fresh": 96, "external8": 24}
    if len(views) != expected[cohort] or len({view[:2] for view in views}) != len(views):
        raise ValueError(f"Unexpected frozen coverage for {cohort}.")
    rows = [evaluate(levels, endpoint_models, view) for view in views]
    archived_rows, archived_cost, control_source = existing_controls(cohort, views)
    attach_controls(rows, archived_rows)
    result = {
        "cohort": cohort, "conditions": len(rows), "independent_objects": len({r["object_id"] for r in rows}),
        "scope": "retrospective replay of all previously exposed historical views"
                 if cohort == "historical" else "frozen diagnostic cohort",
        "protocol": jp.relative_receipt(FROZEN), "setup_cache": jp.relative_receipt(CACHE),
        "galerkin_new_setup": setup, "previously_paid_anchor_lower_setup": paid_cost,
        "previously_paid_anchor_lower_rhs": 60,
        "sampled_native_preparation_by_level": archived_cost, "sampled_native_control_source": control_source,
        "galerkin_summary": summarize(rows, "levels", LEVELS, proposed=True),
        "sampled_native_control_summary": summarize(rows, "sampled_native_controls", existing.LEVELS, proposed=False),
        "rows": rows, **jp.elapsed(begin),
    }
    jp.atomic_json(RESULTS / f"galerkin_{cohort}.json", result)
    if len(rows) == 1:
        row = rows[0]
        print(cohort, row["object_id"], "radius", row["radius"],
              [(n, r["lower"], r["upper"], r["verdict"]) for n, r in row["levels"].items()], flush=True)
    else:
        print(cohort, {n: r["verdict_counts"] for n, r in result["galerkin_summary"].items()}, flush=True)


def paired_summary():
    """Summarize the retrospective matched comparison without rerunning proofs."""
    cohorts = ("validation", "historical", "fresh", "external8")
    measures = ("wall_seconds", "cpu_seconds")
    all_rows, summaries = [], {}
    for cohort in cohorts:
        result = json.loads((RESULTS / f"galerkin_{cohort}.json").read_text())
        controls = json.loads((RESULTS / f"{cohort}_segment_cover.json").read_text())
        rows = result["rows"]
        all_rows.extend(rows)
        midpoints = {(r["object_id"], r["noise_pct"]): r["fixed_midpoint_lower"] for r in controls["rows"]}
        summaries[cohort] = {
            "views": len(rows), "objects": result["independent_objects"],
            "exact_zero_by_galerkin_level": {
                str(n): sum(row["levels"][str(n)]["upper"] == 0. for row in rows) for n in LEVELS
            },
            "exact_zero_by_sampled_native_level": {
                str(n): sum(row["sampled_native_controls"][str(n)]["upper"] == 0. for row in rows)
                for n in existing.LEVELS
            },
            "individual_fixed_anchor_passes": {
                str(anchor): sum(row["endpoint_bounds"][index]["distance_lower"] > row["radius"] for row in rows)
                for index, anchor in enumerate(PROTOCOL["anchor_indices"])
            },
            "fixed_midpoint_passes": sum(midpoints[row["object_id"], row["noise_pct"]] > row["radius"]
                                         for row in rows),
            "nine_node_passes_missed_by_best_fixed_anchor": [
                [row["object_id"], row["noise_pct"]] for row in rows
                if row["sampled_native_controls"]["9"]["verdict"] == "PASS"
                and row["levels"]["4"]["lower"] <= row["radius"]
            ],
            "m4_median_incremental_roof": {unit: median(r["levels"]["4"][unit] for r in rows)
                                           for unit in measures},
            "galerkin": result["galerkin_summary"],
            "sampled_native": result["sampled_native_control_summary"],
        }
    if len(all_rows) != 216 or len({(r["object_id"], r["noise_pct"]) for r in all_rows}) != 216:
        raise ValueError("The archived comparison is not paired across all 216 views.")
    setup = result["galerkin_new_setup"]
    paid = result["previously_paid_anchor_lower_setup"]
    new_m4 = {unit: setup[unit] - sum(setup["compact_levels"][str(n)][unit] for n in LEVELS[1:])
              for unit in measures}
    jp.atomic_json(RESULTS / "galerkin_paired_summary.json", {
        "scope": "216 paired archived cases; the separate validation_000 at 1 percent is an extra focus",
        "cohorts": summaries, "views": len(all_rows),
        "m4_median_incremental_roof": {unit: median(r["levels"]["4"][unit] for r in all_rows)
                                       for unit in measures},
        "galerkin": summarize(all_rows, "levels", LEVELS, proposed=True),
        "sampled_native": summarize(all_rows, "sampled_native_controls", existing.LEVELS, proposed=False),
        "m4_new_preparation": {**new_m4, "native_rhs": setup["native_rhs"]},
        "m4_preparation_including_archived_lower": {
            **{unit: new_m4[unit] + paid[unit] for unit in measures},
            "native_rhs_including_archived_lower": setup["native_rhs"] + result["previously_paid_anchor_lower_rhs"],
        },
        "sampled_native_preparation_by_level": result["sampled_native_preparation_by_level"],
        "setup_cache": jp.relative_receipt(CACHE),
    })


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cohorts = ("validation_focus", "historical_focus", "validation", "fresh", "external8")
    cli.add_argument("cohort", choices=("freeze", "all", "historical", "paired_summary", *cohorts))
    requested = cli.parse_args().cohort
    freeze()
    if requested == "freeze":
        return
    if requested == "paired_summary":
        paired_summary()
        return
    ctx.prec, ctx.threads = jp.PRECISION, PROTOCOL["arb_threads"]
    with threadpool_limits(limits=1):
        levels, setup = prepare()
        paid, paid_cost = paid_endpoint_models()
        for cohort in cohorts if requested == "all" else (requested,):
            run_cohort(cohort, levels, setup, paid, paid_cost)


if __name__ == "__main__":
    main()
