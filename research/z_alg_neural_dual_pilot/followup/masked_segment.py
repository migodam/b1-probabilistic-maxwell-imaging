"""Enclose zero certificates over the fixed segment for the two contiguous arcs."""

import json

from flint import acb_mat, ctx
import numpy as np
from threadpoolctl import threadpool_limits

import masked_ports as mp


jp = mp.jp
segment = mp.segment
ARCS = ("arc_0_deg", "arc_180_deg")
LEVELS = (4, 8, 16, 32)
OUTPUT = mp.RESULTS / "masked_segment.json"


def bordered_forms(models, y, kept, deleted):
    """Build F(t) using the full positive-state response and a fixed embedding."""
    n_deleted = len(deleted)
    y_kept = jp.ball(y[list(kept)])
    y_gram = jp.star(y_kept) * y_kept
    embedding = acb_mat(jp.SIZE, n_deleted + jp.SOURCES)
    constant = acb_mat(n_deleted + jp.SOURCES, n_deleted + jp.SOURCES)
    for col, row in enumerate(deleted):
        embedding[row, col] = 1
        constant[col, col] = 1
    for j in range(jp.SOURCES):
        embedding[jp.PORTS + j, n_deleted + j] = -1
        for i, row in enumerate(kept):
            embedding[row, n_deleted + j] = -y_kept[i, j]
        for i in range(jp.SOURCES):
            constant[n_deleted + i, n_deleted + j] = y_gram[i, j]
    adjoint = jp.star(embedding)
    return {t: constant - adjoint * model.response * embedding for t, model in models.items()}


def continuous_cover(forms):
    """Each successful proof uses the same rank-ten subspace at both roof endpoints."""
    begin = jp.started()
    nodes = sorted(forms)
    cache = {}
    resolved, history = [], []
    pending = []
    for i, (low, high) in enumerate(zip(nodes[:-1], nodes[1:], strict=True)):
        chords = []
        if i:
            chords.append((nodes[i - 1], low))
        if i + 2 < len(nodes):
            chords.append((high, nodes[i + 2]))
        pending.append((low, high, chords))

    for nsub in LEVELS:
        parts = nsub if nsub == LEVELS[0] else 2
        attempts, congruences, next_pending = 0, 0, []
        layer_begin = jp.started()
        for a, b, chords in pending:
            for j in range(parts):
                low = a + (b - a) * j / parts
                high = a + (b - a) * (j + 1) / parts
                attempts += 1
                margins = []
                for p, q in chords:
                    endpoints = []
                    for x in (low, high):
                        if p < x < q:
                            raise ArithmeticError("An interior chord cannot bound F from above.")
                        key = (p, q, x)
                        if key not in cache:
                            cache[key] = mp.affine_between(forms, p, q, x)
                        endpoints.append(cache[key])
                    congruences += 1
                    margin = mp.common_sign_margin(*endpoints, jp.SOURCES, positive=False)
                    if margin is not None and margin > 0.:
                        margins.append((margin, (p, q)))
                if margins:
                    margin, chord = max(margins)
                    resolved.append({"interval": [low, high], "upper_chord": chord,
                                     "negative_directions_lower": jp.SOURCES,
                                     "negative_margin_lower": margin})
                else:
                    next_pending.append((low, high, chords))
        pending = next_pending
        history.append({
            "subdivisions_per_node_interval": nsub,
            "attempted_intervals": attempts,
            "endpoint_pair_congruences": congruences,
            "unresolved": [[a, b] for a, b, _ in pending],
            **jp.elapsed(layer_begin),
        })
        if not pending:
            break

    unresolved = [{"interval": [a, b], "attempted_upper_chords": chords} for a, b, chords in pending]
    pieces = sorted(resolved + unresolved, key=lambda row: row["interval"][0])
    edges = [row["interval"] for row in pieces]
    if (edges[0][0] != nodes[0] or edges[-1][1] != nodes[-1]
            or any(left[1] != right[0] for left, right in zip(edges[:-1], edges[1:], strict=True))):
        raise ArithmeticError("The certified and unresolved pieces must partition the whole segment.")
    return {
        "exact_zero_every_t": not unresolved,
        "status": "PROVED_ZERO" if not unresolved else "UNKNOWN",
        "subdivision_history": history,
        "proved_pieces": sorted(resolved, key=lambda row: row["interval"][0]),
        "unresolved_pieces": unresolved,
        "worst_proved_margin_lower": min((r["negative_margin_lower"] for r in resolved), default=None),
        "unique_affine_roof_values": len(cache),
        **jp.elapsed(begin),
    }


def index_rows(rows):
    indexed = {(row["object_id"], row["noise_pct"]): row for row in rows}
    if len(indexed) != len(rows):
        raise ValueError("Duplicate rows in a frozen observation receipt.")
    return indexed


def run_cohort(name, models, geometry):
    begin = jp.started()
    masked_path = mp.RESULTS / f"{name}_masked_ports.json"
    masked = index_rows(json.loads(masked_path.read_text())["rows"])
    full = mp.propagated_rows(name)
    observations = mp.views_for(name)
    if set(masked) != set(full) or set(full) != {(oid, noise) for oid, noise, _, _ in observations}:
        raise ValueError(f"The frozen observation IDs differ for {name}.")

    rows = []
    for oid, noise, y, radius in observations:
        full_row, masked_row = full[(oid, noise)], masked[(oid, noise)]
        radius = float(radius)
        if (y.shape != (jp.PORTS, jp.SOURCES) or not np.isfinite(y).all()
                or not np.isfinite(radius) or radius < 0.
                or radius != full_row["radius"] or radius != masked_row["radius"]):
            raise ValueError(f"Invalid or mismatched original observation for {name}/{oid}/{noise}.")
        full_upper = full_row["levels"]["9"]["continuous_upper"]
        if full_upper != masked_row["passive_all_masks_continuous_propagation"]["upper"]:
            raise ArithmeticError("The two receipts disagree about the inherited continuous upper bound.")
        arc_results = {}
        for arc in ARCS:
            arc_started = jp.started()
            sampled = masked_row["masks"][arc]
            point_t = {point["t"] for point in sampled["points"]}
            if point_t != set(models) or sampled["max_point_upper"] != 0.:
                raise ArithmeticError("Expected the already verified zero sampled-node upper bounds.")
            if full_upper == 0.:
                result = {"exact_zero_every_t": True, "status": "PROVED_ZERO",
                          "method": "inherited_full_ring_zero_by_receiver_contraction"}
            else:
                coords = geometry[arc]["0.0"]
                forms = bordered_forms(models, y, coords["kept"], coords["deleted"])
                result = {"method": "rank_ten_negative_of_concave_border", **continuous_cover(forms)}
            upper = 0. if result["exact_zero_every_t"] else full_upper
            arc_results[arc] = result | {
                "max_sampled_distance_upper": sampled["max_point_upper"],
                "continuous_distance_upper": upper,
                "strictly_below_original_radius": upper < radius,
                "total_including_forms": jp.elapsed(arc_started),
            }
        rows.append({"object_id": oid, "noise_pct": noise, "radius": radius,
                     "full_ring_nine_node_continuous_upper": full_upper,
                     "full_ring_nine_node_verdict": full_row["levels"]["9"]["verdict"],
                     "arcs": arc_results})
    inherited = sum(row["full_ring_nine_node_continuous_upper"] == 0. for row in rows)
    return {
        "views": len(rows),
        "materials": len({row["object_id"] for row in rows}),
        "full_ring_receipt": jp.relative_receipt(mp.RESULTS / f"{name}_segment_cover.json"),
        "arc_sampled_receipt": jp.relative_receipt(masked_path),
        "zero_views_inherited_by_each_arc": inherited,
        "views_checked_directly_on_each_arc": len(rows) - inherited,
        "arcs": {
            arc: {
                "new_exact_zero_views": sum(r["arcs"][arc]["exact_zero_every_t"]
                                             and r["full_ring_nine_node_continuous_upper"] > 0. for r in rows),
                "total_exact_zero_views": sum(r["arcs"][arc]["exact_zero_every_t"] for r in rows),
                "strictly_below_original_radius": sum(r["arcs"][arc]["strictly_below_original_radius"]
                                                       for r in rows),
                "unresolved_views": [
                    [r["object_id"], r["noise_pct"]] for r in rows if not r["arcs"][arc]["exact_zero_every_t"]
                ],
            }
            for arc in ARCS
        },
        "rows": rows,
        **jp.elapsed(begin),
    }


def main():
    begin = jp.started()
    ctx.prec, ctx.threads = jp.PRECISION, jp.THREADS
    with threadpool_limits(limits=1):
        models, cache = segment.load_cache(segment.CACHE)
        if set(models) != set(jp.GRID):
            raise ValueError("The nine frozen full-ring interval responses are required.")
        loaded_cache = jp.elapsed(begin)
        geometry_path = mp.RESULTS / "masked_geometry.json"
        geometry = json.loads(geometry_path.read_text())["per_mask_geometry"]
        continuous_path = mp.RESULTS / "masked_continuous_inertia.json"
        previous = json.loads(continuous_path.read_text())["profiles"]
        for arc in ARCS:
            if (previous[arc]["negative_count_every_t"] != jp.SOURCES - 1
                    or set(geometry[arc]["0.0"]["kept"]) != set(mp.MASKS[arc])):
                raise ValueError(f"Unexpected receiver geometry or continuous inertia for {arc}.")
        cohorts = {}
        for name in mp.COHORTS:
            cohorts[name] = run_cohort(name, models, geometry)
            print(name, cohorts[name]["arcs"], flush=True)
    inherited = sum(r["zero_views_inherited_by_each_arc"] for r in cohorts.values())
    checked = sum(r["views_checked_directly_on_each_arc"] for r in cohorts.values())
    jp.atomic_json(OUTPUT, {
        "policy": {
            "cap": jp.CAP, "receiver_masks": ARCS, "nodes": jp.GRID, "subdivisions": LEVELS,
            "refinement": "Bisect only unresolved subintervals at successive frozen levels.",
            "criterion": "One strictly negative rank-ten compression of both endpoints of an exterior chord.",
            "precision_bits": jp.PRECISION, "flint_threads": jp.THREADS,
            "input": "Only original full-ring noisy Y, same receiver rows, all sources, and original absolute radius.",
        },
        "response_cache": jp.relative_receipt(segment.CACHE),
        "receiver_geometry_receipt": jp.relative_receipt(geometry_path),
        "previous_continuous_inertia_receipt": jp.relative_receipt(continuous_path),
        "predecessor_native_rhs_and_cost": cache["setup"]["preparation_by_level"]["9"],
        "new_native_rhs": 0,
        "gt_reads": 0,
        "response_cache_load": loaded_cache,
        "view_arc_combinations_inherited": inherited * len(ARCS),
        "view_arc_combinations_checked": checked * len(ARCS),
        "cohorts": cohorts,
        **jp.elapsed(begin),
    })


if __name__ == "__main__":
    main()
