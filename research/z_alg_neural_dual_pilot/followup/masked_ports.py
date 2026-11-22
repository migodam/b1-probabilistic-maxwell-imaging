"""Outward receiver-mask state/source replay from the fixed full-ring responses."""

import argparse
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
from statistics import median
import sys

sys.dont_write_bytecode = True

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator, cg
from scipy.special import hankel1, jv
from threadpoolctl import threadpool_limits

import jensen_port as jp
import segment_cover as segment


RESULTS = jp.HERE / "results"
PROTOCOL = jp.HERE / "MASKED_PORTS_FROZEN.json"
MASKS = {
    "full": tuple(range(20)),
    "alternating_even": tuple(range(0, 20, 2)),
    "alternating_odd": tuple(range(1, 20, 2)),
    "arc_0_deg": tuple(range(10)),
    "arc_180_deg": tuple(range(10, 20)),
}
COHORTS = {
    "validation": ("synthetic_observations.npz", 12, 12),
    "fresh": ("followup_observations.npz", 32, 96),
    "external8": ("followup_external8.npz", 8, 24),
    "historical": (None, 28, 84),
}
SANITY = (("full", 0.), ("alternating_even", .5), ("arc_0_deg", .5))
RULE = {
    "cap": jp.CAP, "t_nodes": list(jp.GRID),
    "physical_arrays": {key: list(values) for key, values in MASKS.items()},
    "input_contract": "native Y port rows only, all ten original illuminations; keep the original common absolute radius",
    "blocks": "same exact dyadic weights on the published fixed16_0 to fixed16_2 path",
    "D_sign": "verify both congruence signs with Arb; retain source eigenvalues after the top nu terms",
    "unresolved": "UNKNOWN if signs or the compact solve cannot be enclosed",
    "direct_96_squared_checks": [{"mask": mask, "t": t} for mask, t in SANITY],
    "direct_96_squared_observation": "lexically lowest fresh object ID at one percent, selected before scoring",
    "allocation": "every mask and weight is evaluated on the same frozen noisy view, no mask selection",
}


@dataclass
class MaskModel:
    response: acb_mat | None
    receipt: dict


def freeze():
    if PROTOCOL.is_file():
        if json.loads(PROTOCOL.read_text())["protocol"] != RULE:
            raise ValueError("The saved mask protocol differs from this run.")
    else:
        jp.atomic_json(PROTOCOL, {"frozen_utc": datetime.now(timezone.utc).isoformat(), "protocol": RULE})


def select(a, rows, columns):
    return acb_mat(len(rows), len(columns), [a[i, j] for i in rows for j in columns])


def identity(size):
    m = acb_mat(size, size)
    for i in range(size):
        m[i, i] = 1
    return m


def frobenius_upper(matrix):
    # The ball sum of squares can cross zero due to interval dependency.
    # Only its outward nonnegative endpoint is needed before taking a root.
    squared_upper = jp.frobenius_squared(matrix).upper()
    if not squared_upper.is_finite():
        return squared_upper
    return max(arb(0), squared_upper).sqrt().upper()


def verified_receiver_order():
    begin = jp.started()
    coordinates_path = jp.archived.PAPER / "synthetic/public/pw360_10x20.npz"
    with np.load(coordinates_path, allow_pickle=False) as f:
        points = f["rx_m"]
    _, stored_s, _ = jp.archived.native_geometry()
    if points.shape != (jp.PORTS, 2):
        raise ValueError("Unexpected frozen receiver coordinates.")
    degrees = np.mod(np.degrees(np.arctan2(points[:, 1], points[:, 0])), 360.)
    degrees[np.isclose(degrees, 360., rtol=0, atol=1e-8)] = 0.
    order = np.argsort(degrees, kind="stable")
    if not np.allclose(degrees[order], np.arange(jp.PORTS) * 18., rtol=0, atol=1e-10):
        raise ValueError("The frozen geometry is not the expected 20-position angular ring.")
    # Check that the coordinate row order is also S's order on the 96² grid.
    step = 2. / 96
    coord = (np.arange(96) + .5) * step - 1.
    yy, xx = np.meshgrid(coord, coord, indexing="ij")
    ctr = np.column_stack((xx.ravel(), yy.ravel()))
    distance = np.linalg.norm(points[:, None] - ctr, axis=2)
    k = 2 * np.pi * 400e6 / 3e8
    cell_radius = step / np.sqrt(np.pi)
    recreated = .5j * np.pi * k * cell_radius * jv(1, k * cell_radius) * hankel1(0, k * distance)
    discrepancy = float(np.max(abs(recreated - stored_s)))
    if discrepancy > 1e-13:
        raise ValueError("Receiver positions cannot be matched to the original native S row order.")
    return tuple(int(k) for k in order), {
        "coordinate_source": jp.relative_receipt(coordinates_path),
        "stored_native_s_source": jp.relative_receipt(jp.archived.FINE),
        "receiver_degrees_in_storage_order": [round(float(x), 10) for x in degrees],
        "indices_in_angular_order": order.tolist(),
        "numeric_rebuilt_s_max_entry_error": discrepancy,
        **jp.elapsed(begin),
    }


def signed_inertia(d):
    if d.nrows() == 0:
        return {"nu": 0, "positive": 0, "absolute_gap_lower": None, "status": "RESOLVED"}
    values, basis = eigh(jp.hermitian(jp.mid(d)), check_finite=False, driver="evd")
    count = int(np.count_nonzero(values < 0.))
    receipt = {"nu": None, "positive": None, "absolute_gap_lower": None,
               "numeric_eigenvalues": values.tolist(), "status": "UNKNOWN"}
    floors = []
    for sign, vectors in ((-1, basis[:, :count]), (1, basis[:, count:])):
        if vectors.shape[1] == 0:
            continue
        v = jp.ball(vectors)
        gram_up = jp.gershgorin(jp.star(v) * v, upper=True)
        congruence = jp.star(v) * d * v
        floor = jp.gershgorin(congruence, upper=False) if sign > 0 else -jp.gershgorin(congruence, upper=True)
        if not (floor > 0. and gram_up > 0.):
            receipt["reason"] = "A proposed sign subspace could not be enclosed."
            return receipt
        edge = (floor / gram_up).lower()
        floors.append(edge)
        receipt["positive_gap_lower" if sign > 0 else "negative_gap_lower"] = jp.low(edge)
    return receipt | {"status": "RESOLVED", "nu": count, "positive": d.nrows() - count,
                      "absolute_gap_lower": jp.low(min(floors))}


def mask_q(response, kept):
    start = jp.started()
    missing = tuple(i for i in range(jp.PORTS) if i not in kept)
    coordinates = kept + tuple(range(jp.PORTS, jp.SIZE))
    if not missing:
        q = response if kept == tuple(range(jp.PORTS)) else select(response, coordinates, coordinates)
        return MaskModel(q, {"kept": kept, "deleted": missing, "nu": 0, "status": "RESOLVED",
                             "inverse_status": "full_reference_no_update", **jp.elapsed(start)})
    d = identity(len(missing)) - select(response, missing, missing)
    proof = signed_inertia(d)
    if proof["status"] == "UNKNOWN":
        return MaskModel(None, proof | {"kept": kept, "deleted": missing, **jp.elapsed(start)})
    nu = proof["nu"]
    if nu >= jp.SOURCES:
        q = None
        inverse_status = "unneeded_all_source_modes_discarded"
    else:
        coupling = select(response, coordinates, missing)
        base = select(response, coordinates, coordinates)
        try:
            q = base + coupling * d.solve(jp.star(coupling))
            q = (q + jp.star(q)) / 2
        except (ValueError, ZeroDivisionError, ArithmeticError) as error:
            return MaskModel(None, proof | {"status": "UNKNOWN", "kept": kept, "deleted": missing,
                                           "inverse_status": "interval_solve_failed", "reason": str(error),
                                           **jp.elapsed(start)})
        inverse_status = "enclosed"
    return MaskModel(q, proof | {"kept": kept, "deleted": missing,
                                 "inverse_status": inverse_status, **jp.elapsed(start)})


def retain_source(h, nu):
    midpoint = jp.hermitian(jp.mid(h))
    if not np.isfinite(midpoint).all():
        return {"distance_lower": 0., "distance_upper": None, "status": "UNKNOWN",
                "reason": "Indeterminate compact source form."}
    vals, vec = eigh(midpoint, check_finite=False, driver="evd")
    vals, vec = vals[::-1], vec[:, ::-1]
    v = jp.ball(vec)
    gram = jp.star(v) * v
    g_lo, g_hi = jp.gershgorin(gram, upper=False), jp.gershgorin(gram, upper=True)
    if not g_lo > 0:
        return {"distance_lower": 0., "distance_upper": None, "status": "UNKNOWN",
                "reason": "Source trial basis is not certified invertible."}
    c = jp.star(v) * h * v
    for i, value in enumerate(vals):
        c[i, i] -= arb(float(value))
    defect = frobenius_upper(c)
    if not defect.is_finite():
        return {"distance_lower": 0., "distance_upper": None, "status": "UNKNOWN",
                "reason": "No finite source congruence defect enclosure."}
    intervals = []
    for value in vals[nu:]:
        a = max(arb(0), (arb(float(value)) - defect).lower())
        b = max(arb(0), (arb(float(value)) + defect).upper())
        intervals.append((a / g_hi, b / g_lo))
    lower = sum((lo for lo, _ in intervals), arb(0))
    upper = sum((hi for _, hi in intervals), arb(0))
    # Sylvester inertia and generalized min-max make these inequalities
    # valid for each positive sorted eigenvalue after the same nu entries.
    endpoint_l = jp.low(lower.sqrt()) if lower > 0. else 0.
    endpoint_u = jp.high(upper.sqrt()) if upper > 0. else 0.
    return {"distance_lower": endpoint_l, "distance_upper": endpoint_u,
            "retained_source_count": jp.SOURCES - nu,
            "retained_positive_eigenvalue_bounds": [
                [jp.low(a) if a > 0. else 0., jp.high(b) if b > 0. else 0.] for a, b in intervals
            ],
            "numeric_source_spectrum_descending": vals.tolist(),
            "source_congruence_defect_upper": jp.high(defect),
            "status": "PROVED_ZERO" if endpoint_u == 0. else "INTERVAL"}


def pointwise(model, y):
    begin = jp.started()
    kept = model.receipt["kept"]
    assert y.shape == (len(kept), jp.SOURCES)
    nu = model.receipt["nu"]
    if model.receipt["status"] == "UNKNOWN":
        return {"status": "UNKNOWN", "nu": nu, "distance_lower": 0., "distance_upper": None,
                **jp.elapsed(begin)}
    if nu >= jp.SOURCES:
        return {"status": "PROVED_ZERO", "nu": nu, "distance_lower": 0., "distance_upper": 0.,
                "retained_source_count": 0, **jp.elapsed(begin)}
    yb = jp.ball(y)
    xi = jp.original.vertical(yb, identity(jp.SOURCES))
    h = jp.star(yb) * yb - jp.star(xi) * model.response * xi
    return retain_source(h, nu) | {"nu": nu, **jp.elapsed(begin)}


def prepare():
    freeze()
    ctx.prec, ctx.threads = jp.PRECISION, jp.THREADS
    begin = jp.started()
    models, original = segment.load_cache(segment.CACHE)
    if set(models) != set(jp.GRID):
        raise ValueError("The nine fixed positive full-ring nodes must already be enclosed.")
    read = jp.elapsed(begin)
    order, order_proof = verified_receiver_order()
    start = jp.started()
    masks = {
        name: {t: mask_q(model.response, tuple(order[a] for a in ids)) for t, model in models.items()}
        for name, ids in MASKS.items()
    }
    record = {
        "full_node_cache": jp.relative_receipt(segment.CACHE),
        "source_setup_predecessor_9_nodes": original["setup"]["preparation_by_level"]["9"],
        "geometry_read": read,
        "receiver_order": order_proof,
        "pointwise_mask_setup": jp.elapsed(start),
        "extra_native_rhs_for_all_mask_production": 0,
        "per_mask_geometry": {name: {str(t): model.receipt for t, model in table.items()}
                              for name, table in masks.items()},
        **jp.elapsed(begin),
    }
    jp.atomic_json(RESULTS / "masked_geometry.json", record)
    return masks, record


def views_for(cohort):
    name, objects, conditions = COHORTS[cohort]
    views = (jp.historical_observations() if name is None
             else jp.read_observations(jp.PILOT / "data" / name, validation=cohort == "validation"))
    if len(views) != conditions or len({row[0] for row in views}) != objects:
        raise ValueError(f"Frozen scope of {cohort} differs from the registered mask protocol.")
    return views


def propagated_rows(cohort):
    cover = RESULTS / f"{cohort}_segment_cover.json"
    raw = json.loads(cover.read_text())
    return {(r["object_id"], r["noise_pct"]): r for r in raw["rows"]}


def evaluate(masks, views, propagated):
    rows = []
    for oid, noise, y, radius in views:
        radius = float(radius)
        previous = propagated[(oid, noise)]
        if previous["radius"] != radius:
            raise ValueError("Full-ring propagation must use the same frozen absolute radius.")
        full_continuous = previous["levels"]["9"]["continuous_upper"]
        old_points = {record["t"]: record for record in previous["nodes"]}
        row = {
            "object_id": oid, "noise_pct": noise, "radius": radius,
            "passive_all_masks_continuous_propagation": {
                "upper": full_continuous,
                "exact_zero": full_continuous == 0.,
                "strictly_below_common_radius": full_continuous < radius,
                "method": "bordered_positive_count_monotonicity_under_receiver_contraction",
            },
            "masks": {},
        }
        for name, table in masks.items():
            y_kept = y[list(table[0.].receipt["kept"])]
            values = []
            for t, model in table.items():
                bound = pointwise(model, y_kept)
                bound.update(t=t, passed=bound["distance_lower"] > radius)
                if name == "full":
                    old = old_points[t]
                    if (bound["distance_upper"] is None or bound["distance_lower"] > old["distance_upper"]
                            or old["distance_lower"] > bound["distance_upper"]):
                        raise ArithmeticError("Full algebra check disagrees with the prior native node interval.")
                    bound["overlap_with_original_full_node_interval"] = True
                else:
                    full_upper = next(v["distance_upper"] for v in row["masks"]["full"]["points"] if v["t"] == t)
                    if bound["distance_lower"] > full_upper:
                        raise ArithmeticError("A masked lower endpoint exceeds the full-ring upper endpoint.")
                if bound["distance_lower"] > full_continuous:
                    raise ArithmeticError("A pointwise lower bound exceeds the certified full-segment upper bound.")
                values.append(bound)
            best = max(values, key=lambda record: (record["distance_lower"], -record["t"]))
            row["masks"][name] = {
                "points": values, "best_point_t": best["t"], "best_point_lower": best["distance_lower"],
                "best_point_passed": best["passed"],
                "max_point_upper": (None if any(v["distance_upper"] is None for v in values)
                                    else max(v["distance_upper"] for v in values)),
            }
        rows.append(row)
    return rows


def summarize(rows):
    summary = {}
    bands = {"all": rows} | {
        f"{noise:g}_percent": [r for r in rows if r["noise_pct"] == noise]
        for noise in sorted({r["noise_pct"] for r in rows})
    }
    for name, band in bands.items():
        summary[name] = {
            "views": len(band),
            "universal_strict_below_via_full_segment": sum(r["passive_all_masks_continuous_propagation"]
                                                          ["strictly_below_common_radius"] for r in band),
            "universal_exact_zero_via_full_segment": sum(r["passive_all_masks_continuous_propagation"]
                                                        ["exact_zero"] for r in band),
            "masks": {},
        }
        for mask in MASKS:
            all_points = [p for r in band for p in r["masks"][mask]["points"]]
            summary[name]["masks"][mask] = {
                "view_sampled_passes": sum(r["masks"][mask]["best_point_passed"] for r in band),
                "point_passes": sum(p["passed"] for p in all_points),
                "point_exact_zeros": sum(p["status"] == "PROVED_ZERO" for p in all_points),
                "point_positive_lowers": sum(p["distance_lower"] > 0. for p in all_points),
                "numerical_unknown": sum(p["status"] == "UNKNOWN" for p in all_points),
                "nu_counts": dict(Counter(str(p["nu"]) for p in all_points)),
                "median_point_wall_seconds": median(p["wall_seconds"] for p in all_points),
                "median_point_cpu_seconds": median(p["cpu_seconds"] for p in all_points),
            }
    return summary


def independent_d(full, sensor_block, normalized_gap):
    """Independent 96² 10-RHS solve used once for the direct signed reference."""
    started = jp.started()
    g = full.g
    scale = arb(float(full.w.scale))
    s = jp.ball(sensor_block)
    ell = jp.star(s) / scale
    pre = LinearOperator((g.n, g.n), matvec=full.precondition, dtype=np.complex128)
    columns, iterations = [], []
    for l in (sensor_block.conj().T / full.w.scale).T:
        track = []
        x, info = cg(full.operator, l, M=pre, rtol=3e-13, atol=0., maxiter=700,
                     callback=lambda _: track.append(1))
        if info:
            raise ArithmeticError(f"Independent native deleted-row solve failed, info={info}.")
        columns.append(x)
        iterations.append(len(track))
    x = jp.ball(np.column_stack(columns))
    ax = full.exact_action(x)
    defect = ax - ell
    cap = scale * (jp.star(defect) * defect) / arb(normalized_gap)
    center = scale * (jp.star(ell) * x + jp.star(x) * ell - jp.star(x) * ax)
    size = len(sensor_block)
    q = acb_mat(size, size)
    bounds = [max(arb(0), cap[i, i].real.upper()) for i in range(size)]
    for i, di in enumerate(bounds):
        q[i, i] = acb(center[i, i].real + arb(0).union(di))
        for j in range(i + 1, size):
            eps = (di * bounds[j]).sqrt().upper()
            enclosure = (-eps).union(eps)
            q[i, j] = (center[i, j] + center[j, i].conjugate()) / 2 + acb(enclosure, enclosure)
            q[j, i] = q[i, j].conjugate()
    proof = signed_inertia(identity(size) - q)
    if proof["status"] != "RESOLVED":
        raise ArithmeticError("The independent D sign proof did not resolve.")
    gap_d = arb(proof["absolute_gap_lower"])
    g_normal = arb(normalized_gap)
    # ||(K_kept/scale)^{-1}|| <= 1/g + ||S_deleted||^2/(scale g^2 gap_D).
    gap = 1 / (1 / g_normal + jp.frobenius_squared(s) / (scale * g_normal**2 * gap_d))
    return proof | {"detector_rhs": size, "CG_iterations": iterations,
                    "direct_mask_normalized_absolute_gap_lower": jp.low(gap),
                    "stationary_defect_upper": jp.high(frobenius_upper(defect)),
                    **jp.elapsed(started)}


def direct_checks(masks):
    start = jp.started()
    oid, noise, y, radius = min((r for r in views_for("fresh") if r[1] == 1.),
                                key=lambda row: row[0])
    profiles = jp.archived_anchors()
    spectrum, sensor, illumination = jp.archived.native_geometry()
    g = min(p.archived_inertia["positive_gap"] for p in profiles)
    full = jp.archived.Geometry(spectrum, sensor, illumination, block_side=4)
    direct_rows = []
    for name, t in SANITY:
        stage = jp.started()
        model = masks[name][t]
        kept = list(model.receipt["kept"])
        missing = list(model.receipt["deleted"])
        subset = y[kept]
        geometry = full if not missing else jp.archived.Geometry(spectrum, sensor[kept], illumination, block_side=4)
        action = segment.DyadicProfile(geometry, profiles[0].specification, profiles[1].specification, t)
        receipt = {}
        if name == "full":
            n, absolute_gap = 0, g
            receipt["state_method"] = "archived_independent_full_state_PSD_proof"
        elif name == "alternating_even":
            receipt["native_inertia_numeric_proposal"] = action.modes(jp.SOURCES)
            receipt["direct_negative_subspace_upper"], _ = action.ensure_negative(jp.SOURCES, return_gap=True)
            # The direct strict 10-dimensional negative subspace equals the
            # rank of the downdate, so the known positive full operator gives
            # the complementary count. No stationary source solve is needed.
            n, absolute_gap = len(missing), None
            receipt["state_method"] = "direct_negative_subspace_plus_downdate_rank_upper"
        else:
            full_action = segment.DyadicProfile(full, profiles[0].specification, profiles[1].specification, t)
            sign = independent_d(full_action, sensor[missing], g)
            n, absolute_gap = sign["nu"], sign["direct_mask_normalized_absolute_gap_lower"]
            receipt["independent_deleted_row_gap"] = sign
            receipt["native_inertia_numeric_proposal"] = action.modes(n)
            receipt["direct_negative_subspace_upper"], _ = action.ensure_negative(n, return_gap=True)
            receipt["state_method"] = "independent_deleted_row_stationarity_plus_direct_negative_subspace"
        if n >= jp.SOURCES:
            lower = upper = 0.
            receipt["direct_masked_stationary_rhs"] = 0
        else:
            trial, vec, vals, numeric = action.numerical_distance(subset, n)
            result = action.ensure_source(subset, n, absolute_gap, trial, vec, vals,
                                          negative_gap=absolute_gap)
            lower, upper = result["distance_lower"], result["distance_upper"]
            receipt.update({"direct_masked_stationary_rhs": jp.SOURCES,
                            "original_native_stationary_numeric": numeric, "original_native_Arb": result})
        compact = pointwise(model, subset)
        agree = (compact["distance_upper"] is not None and n == compact["nu"]
                 and lower <= compact["distance_upper"] and compact["distance_lower"] <= upper)
        if not agree:
            raise ArithmeticError(f"Independent native/compact signed source checks disagree in {name}.")
        direct_rows.append({"mask": name, "t": t, "object_id": oid, "noise_pct": noise,
                            "radius": float(radius), "nu": n, "native_direct_interval": [lower, upper],
                            "compact_result": compact, "verified_overlap_and_inertia_agreement": agree,
                            "native_computation": receipt, **jp.elapsed(stage)})
        print("Direct", name, "nu", n, "native interval", [lower, upper], flush=True)
        jp.atomic_json(RESULTS / "masked_direct_checks.json", {"rows": direct_rows,
                       "physics_rhs_in_production_masks": 0, **jp.elapsed(start)})
    return direct_rows


def affine_between(points, a, b, x):
    w = (arb(x) - arb(a)) / (arb(b) - arb(a))
    return (1 - w) * points[a] + w * points[b]


def common_sign_margin(left, right, dimension, *, positive):
    if dimension == 0:
        return None
    numerical = jp.hermitian(jp.mid(left) + jp.mid(right)) / 2
    _, candidates = eigh(numerical, check_finite=False, driver="evd")
    guess = candidates[:, -dimension:] if positive else candidates[:, :dimension]
    basis = jp.ball(guess)
    gram_upper = jp.gershgorin(jp.star(basis) * basis, upper=True)
    if not gram_upper > 0.:
        return None
    signed = [jp.star(basis) * matrix * basis for matrix in (left, right)]
    if positive:
        gap = min(jp.gershgorin(m, upper=False) for m in signed)
    else:
        gap = -max(jp.gershgorin(m, upper=True) for m in signed)
    if not gap > 0.:
        return None
    return jp.low(gap / gram_upper)


def prove_continuous_array(d, expected_neg):
    """Use matrix concavity of D(t) and certify every angular-geometry interval."""
    start = jp.started()
    nodes = sorted(d)
    piece_cache = {}

    def affine(a, b, x):
        key = (a, b, x)
        if key not in piece_cache:
            piece_cache[key] = affine_between(d, a, b, x)
        return piece_cache[key]

    expected_pos = d[nodes[0]].nrows() - expected_neg
    history = []
    for nsub in (4, 8, 16, 32, 64):
        pieces = []
        for index, (a, b) in enumerate(zip(nodes[:-1], nodes[1:], strict=True)):
            upper_chords = []
            if index > 0:
                upper_chords.append((nodes[index - 1], a))
            if index + 2 < len(nodes):
                upper_chords.append((b, nodes[index + 2]))
            for j in range(nsub):
                low = a + (b - a) * j / nsub
                high = a + (b - a) * (j + 1) / nsub
                negatives = []
                for p, q in upper_chords:
                    left, right = affine(p, q, low), affine(p, q, high)
                    margin = common_sign_margin(left, right, expected_neg, positive=False)
                    if margin is None and expected_pos == 0:
                        # A full negative cone is convex, so these two full
                        # rank sign checks need not share an eigenbasis.
                        endpoints = [signed_inertia(m) for m in (left, right)]
                        if all(r["nu"] == expected_neg for r in endpoints):
                            margin = min(r["negative_gap_lower"] for r in endpoints)
                    if margin is not None and margin > 0.:
                        negatives.append((margin, (p, q)))
                if expected_pos:
                    positive = common_sign_margin(affine(a, b, low), affine(a, b, high),
                                                  expected_pos, positive=True)
                else:
                    positive = None
                record = {"interval": [low, high], "lower_chord": [a, b],
                          "negative_margin_lower": None,
                          "positive_margin_lower": positive,
                          "negative_directions_lower": 0,
                          "positive_directions_lower": expected_pos if positive is not None and positive > 0. else 0}
                if negatives:
                    negative, chord = max(negatives)
                    record.update(negative_margin_lower=negative, upper_chord=chord,
                                  negative_directions_lower=expected_neg)
                record["resolved"] = record["negative_directions_lower"] + record["positive_directions_lower"] == d[a].nrows()
                pieces.append(record)
        unresolved = [p["interval"] for p in pieces if not p["resolved"]]
        history.append({"subdivisions_per_node_interval": nsub, "unresolved": unresolved})
        if not unresolved:
            break
    neg = [x["negative_margin_lower"] for x in pieces if x["negative_margin_lower"] is not None]
    pos = [x["positive_margin_lower"] for x in pieces if x["positive_margin_lower"] is not None]
    negative_check = all(p["negative_directions_lower"] == expected_neg for p in pieces)
    positive_check = all(p["positive_directions_lower"] == expected_pos for p in pieces)
    return {
        "status": "RESOLVED" if not unresolved else "UNKNOWN",
        "negative_count_every_t": expected_neg if not unresolved else None,
        "at_least_negative_every_t": expected_neg if negative_check else 0,
        "at_least_positive_every_t": expected_pos if positive_check else 0,
        "universal_zero_every_Y_every_t": negative_check and expected_neg >= jp.SOURCES,
        "worst_certified_negative_margin_lower": min(neg) if negative_check and neg else None,
        "worst_certified_positive_margin_lower": min(pos) if positive_check and pos else None,
        "adaptive_search_geometry_only": history, "pieces": pieces,
        "n_unique_dyadic_chord_evals": len(piece_cache),
        **jp.elapsed(start),
    }


def continuous_geometry(masks, setup):
    started = jp.started()
    models, _ = segment.load_cache(segment.CACHE)
    output = {}
    for name in MASKS:
        if name == "full":
            continue
        missing = tuple(masks[name][0.].receipt["deleted"])
        d = {t: identity(len(missing)) - select(full.response, missing, missing) for t, full in models.items()}
        counts = {model.receipt["nu"] for model in masks[name].values()}
        if None in counts or len(counts) != 1:
            output[name] = {"status": "UNKNOWN", "reason": "No shared strict sampled inertia."}
        else:
            output[name] = prove_continuous_array(d, counts.pop())
        print("Continuous", name, output[name]["status"], output[name].get("negative_count_every_t"),
              output[name].get("universal_zero_every_Y_every_t"), flush=True)
    record = {
        "contract": jp.relative_receipt(PROTOCOL),
        "geometry_setup": setup["pointwise_mask_setup"],
        "original_nine_exact_nodes": jp.relative_receipt(segment.CACHE),
        "extra_native_rhs": 0, "observations_opened": 0,
        "method": ("D concave from affine positive K and constant deleted receiver block. "
                   "Outer-chord upper proves a common negative subspace on each interval; "
                   "the inner-chord lower proves a common positive subspace."),
        "profiles": output,
        **jp.elapsed(started),
    }
    jp.atomic_json(RESULTS / "masked_continuous_inertia.json", record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-direct", action="store_true",
                        help="Use only compact per-mask calculations, with zero extra native RHS.")
    parser.add_argument("--continuous-only", action="store_true", help="Prove geometry-only signed counts before any Y.")
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        masks, setup = prepare()
        if args.continuous_only:
            continuous_geometry(masks, setup)
            return
        for cohort in COHORTS:
            start = jp.started()
            rows = evaluate(masks, views_for(cohort), propagated_rows(cohort))
            summary = summarize(rows)
            record = {"cohort": cohort, "scope": COHORTS[cohort],
                      "frozen_contract": jp.relative_receipt(PROTOCOL),
                      "mask_geometry_receipt": jp.relative_receipt(RESULTS / "masked_geometry.json"),
                      "passive_full_segment_receipt": jp.relative_receipt(RESULTS / f"{cohort}_segment_cover.json"),
                      "summary": summary, "rows": rows, **jp.elapsed(start)}
            jp.atomic_json(RESULTS / f"{cohort}_masked_ports.json", record)
            print(cohort, len(rows), [(name, info["view_sampled_passes"], info["point_exact_zeros"],
                                      info["point_positive_lowers"], info["numerical_unknown"])
                                     for name, info in summary["all"]["masks"].items()], flush=True)
        if not args.skip_direct:
            direct_checks(masks)


if __name__ == "__main__":
    main()
