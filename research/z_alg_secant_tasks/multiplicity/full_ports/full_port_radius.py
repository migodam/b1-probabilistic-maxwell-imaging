"""Bound common full-matrix noise radii from already certified root enclosures."""

from itertools import combinations
import json
from math import inf, nextafter
from pathlib import Path

from flint import acb, arb, ctx


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
ctx.prec = 192
ctx.threads = 1


def real_ball(endpoints):
    lo, hi = (arb(float(x)) for x in endpoints)
    return lo.union(hi)


def complex_ball(box):
    return acb(real_ball(box["real"]), real_ball(box["imag"]))


def complex_center(box):
    endpoints = [[arb(float(v)) for v in box[key]] for key in ("real", "imag")]
    return acb(*(sum(row)/2 for row in endpoints))


def outward(box):
    return [nextafter(float(box.lower()), -inf), nextafter(float(box.upper()), inf)]


def complex_outward(z):
    return {"real": outward(z.real), "imag": outward(z.imag)}


def ball_norm(values):
    lower, upper = arb(0), arb(0)
    for value in values:
        lo, hi = value.abs_lower(), value.abs_upper()
        lower += lo * lo
        upper += hi * hi
    return lower.sqrt().union(upper.sqrt())


def envelope(values, order):
    """Enclose min/max of the same list of unknown quantities."""
    return order(x.lower() for x in values).union(order(x.upper() for x in values))


def bound_cases(cases, local_target, root_key, cross_key, source, group):
    results = []
    for case in cases:
        roots = case[root_key]
        cross = [tuple(x) for x in case[cross_key]]
        if len(roots) not in (4, 8) or any(x["status"] != "certified_simple_real_root" for x in roots):
            raise ValueError("Expected a complete set of previously certified nominal roots.")
        local_size = len(local_target)
        count = len(roots[0]["all_complex_port_intervals"])
        columns = len(roots[0]["all_complex_port_intervals"][0])
        expected_cross = [(i, j) for i in range(count) for j in range(columns) if j//local_size != i]
        if cross != expected_cross or columns != count*local_size or len(roots) != 2**count:
            raise ValueError("Cross-port ordering or number of certified modes changed.")
        matrices = [root["all_complex_port_intervals"] for root in roots]
        cross_ball = [[complex_ball(row[i][j]) for i, j in cross] for row in matrices]
        cross_mean = [
            sum(complex_center(row[i][j]) for row in matrices)/len(roots)
            for i, j in cross
        ]

        # In all root certificates the local data are the same exact y0.
        # Thus they contribute zero to differences but enter the full signal
        # normalization as separately measured repeated complex samples.
        common_signal_norm = ball_norm(local_target*count + cross_mean)
        center_distances = [
            ball_norm(a-b for a, b in zip(x, cross_mean, strict=True)) for x in cross_ball
        ]
        pair_list = []
        for (a, x), (b, z) in combinations(enumerate(cross_ball), 2):
            norm = ball_norm(u-v for u, v in zip(x, z, strict=True))
            pair_list.append((roots[a]["branch"], roots[b]["branch"], norm))
        farthest = envelope([row[2] for row in pair_list], max)
        nearest = envelope([row[2] for row in pair_list], min)
        radius_lower = farthest.lower()/2
        radius_upper = envelope(center_distances, max).upper()
        certain_bound_interval = [
            nextafter(float(radius_lower.lower()), -inf),
            nextafter(float(radius_upper.upper()), inf),
        ]
        # Do not infer the optimal center. This is a diameter lower bound
        # and a constructive upper bound furnished by one admissible center.
        normalized = [(radius_lower / common_signal_norm).lower(),
                      (radius_upper / common_signal_norm).upper()]
        normalized_bounds = [nextafter(float(normalized[0]), -inf),
                             nextafter(float(normalized[1]), inf)]
        cross_index = dict(zip(cross, cross_mean, strict=True))
        common_array = [[cross_index[(i, j)] if (i, j) in cross_index else local_target[j % local_size]
                         for j in range(columns)] for i in range(count)]
        delta_10pct = common_signal_norm / 10

        results.append({
            "group": group, "id": case["id"], "certified_root_source": source,
            "root_count": len(roots), "matrix_shape_receiver_by_shot": [count, columns],
            "local_measured_complex_count": count*local_size, "cross_measured_complex_count": len(cross),
            "common_full_scattered_signal_frobenius": outward(common_signal_norm),
            "fixed_common_reading_full_matrix": [[complex_outward(z) for z in row] for row in common_array],
            "centroid_radius_interval": outward(envelope(center_distances, max)),
            "optimal_all_roots_covering_radius_lower_upper": certain_bound_interval,
            "optimal_radius_relative_to_fixed_common_signal_lower_upper": normalized_bounds,
            "optimal_radius_percent_signal_lower_upper": [
                nextafter(float((100*normalized[0]).lower()), -inf),
                nextafter(float((100*normalized[1]).upper()), inf),
            ],
            "strict_10_percent_common_observation_contains_all_roots": bool(delta_10pct > radius_upper),
            "strict_10_percent_no_other_ball_contains_all_roots": bool(radius_lower > delta_10pct),
            "strict_10_percent_all_pairs_have_disjoint_noise_balls": bool(nearest > 2*delta_10pct),
            "strict_10_percent_some_pair_has_a_common_noise_ball": bool(nearest < 2*delta_10pct),
            "minimum_pairwise_port_frobenius_distance": outward(nearest),
            "minimum_pair_midpoint_radius_percent_common_signal": outward(50*nearest/common_signal_norm),
            "maximum_pairwise_port_frobenius_distance": outward(farthest),
            "root_to_fixed_common_reading": [
                {"branch": root["branch"], "distance": outward(distance)}
                for root, distance in zip(roots, center_distances, strict=True)
            ],
            "pairwise_full_matrix_distances": [
                {"branch_a": a, "branch_b": b, "distance": outward(norm)} for a, b, norm in pair_list
            ],
        })
        print(group, case["id"], certain_bound_interval,
              [100*v for v in normalized_bounds], "minimum pair", outward(nearest), flush=True)
    return results


def main():
    main_name = "results/krawczyk.json"
    conditioned_name = "conditioning/results.json"
    earlier = json.loads((PARENT/main_name).read_text())
    changed = json.loads((PARENT/conditioned_name).read_text())
    historical_local = [acb(arb(float(re)), arb(float(im))) for re, im in earlier["exact_target_complex"]]
    conditioned_local = [complex_ball(changed["displaced_pair"]["exact_target"])]
    results = bound_cases(earlier["results"], historical_local, "root_certificates",
                          "cross_ports_receiver_and_source_zero_based", main_name, "historical_six_disk_cell")
    results.extend(bound_cases(changed["coupled_certificates_at_first_complete_spacing"],
                               conditioned_local, "certificates", "cross_ports_receiver_and_source_indices",
                               conditioned_name, "symmetry_broken_two_disk_cell"))
    result = {
        "scope": "Only these certified roots; fixed deterministic common center, no direct incident field",
        "data": "full receiver-by-shot complex scattered matrix, one sample per acquired port",
        "center": "exact common local target repeated by component; arithmetic mean of saved cross-field box midpoints",
        "radius_upper": "maximum distance to the specified central matrix, enclosed by Arb",
        "radius_lower": "half the largest pairwise Frobenius distance lower bound, valid for any center",
        "reference": "Frobenius norm of the constructed complete scattered matrix, counting each repeated local observation once",
        "precision_bits": ctx.prec, "results": results,
    }
    (HERE/"results.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
