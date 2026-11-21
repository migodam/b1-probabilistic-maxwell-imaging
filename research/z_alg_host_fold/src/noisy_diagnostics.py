"""Independent inner witnesses, feasible bridge and single-minor relaxation gap."""

from fractions import Fraction
from itertools import combinations
import json
import math
import time

from flint import acb, acb_mat, arb, arb_mat, ctx
import numpy as np
from scipy.optimize import root

from noisy_cover import (
    NOISE_RADIUS, PRECISION_BITS, RESULTS, dot, exact_complex, load_inputs,
    outward_number, score_material, seed_real_jacobian, squared_modulus,
)

CIRCUITS = list(combinations(range(6), 4))
RADIUS = arb(NOISE_RADIUS)


def approximate_forward(z, data):
    H, B, R, host = data["H"], data["B"], data["R"].reshape(6), data["y0"]
    field = np.linalg.solve(np.eye(6) - H * z[..., None, :], np.broadcast_to(B, (*z.shape[:-1], 6, 3)))
    return host + np.einsum("i,...i,...ia->...a", R, z, field)


def preconditioned_bases(B, approximate_B):
    pivots = list(max(combinations(range(6), 3), key=lambda ids: abs(np.linalg.det(approximate_B[list(ids)]))))
    free = [i for i in range(6) if i not in pivots]
    numerical = B.dtype != object
    inverse = np.linalg.inv(B[pivots, :]) if numerical else np.asarray(acb_mat(B[pivots, :].tolist()).inv().tolist())
    right = np.zeros((3, 6), dtype=B.dtype)
    right[:, pivots] = inverse
    null = np.eye(6, dtype=B.dtype)[free] - B[free] @ right
    return right, null


def approximate_minor_scores(z, data):
    right, null = preconditioned_bases(data["B"], data["B"])
    A = np.eye(6) - data["H"] * z[:, None, :]
    X, N = right @ A, null @ A
    centered = data["y_critical"] - data["y0"]
    tail = z * data["R"].reshape(6)
    values = []
    for ids in CIRCUITS:
        lam = np.array([
            (-1) ** j * np.linalg.det(N[:, :, [i for i in ids if i != omitted]])
            for j, omitted in enumerate(ids)
        ]).T
        alpha = np.einsum("nai,ni->na", X[:, :, ids], lam)
        residual = alpha @ centered - np.sum(tail[:, ids] * lam, axis=1)
        values.append(abs(residual) ** 2 / (NOISE_RADIUS ** 2 * np.sum(abs(alpha) ** 2, axis=1)))
    return np.max(values, axis=0)


def certified_minor_scores(z, operators, observed, approximate_B):
    H, B, R, host = operators
    right, null = preconditioned_bases(B, approximate_B)
    A = np.array([[acb(int(i == j)) - H[i, j] * z[j] for j in range(6)] for i in range(6)], dtype=object)
    X, N = right @ A, null @ A
    centered = [y - y0 for y, y0 in zip(observed, host, strict=True)]
    tail = [z[i] * R[i] for i in range(6)]
    ratios = []
    for ids in CIRCUITS:
        weights = [
            (-1) ** j * acb_mat(N[:, [i for i in ids if i != omitted]].tolist()).det()
            for j, omitted in enumerate(ids)
        ]
        alpha = [dot(X[i, ids], weights) for i in range(3)]
        residual = dot(centered, alpha) - dot([tail[i] for i in ids], weights)
        denominator = RADIUS ** 2 * sum((squared_modulus(a) for a in alpha), arb(0))
        ratios.append(squared_modulus(residual) / denominator)
    return ratios


def find_gap(data, operators, observed):
    _, sigma, vh = np.linalg.svd(seed_real_jacobian(data))
    rng = np.random.default_rng(7206)
    direction = rng.normal(size=(512, 5))
    direction /= np.linalg.norm(direction, axis=1, keepdims=True)
    candidates = data["zc"] + ((direction / sigma[:5]) @ vh[:5]) * (1.2 * NOISE_RADIUS)
    direct_score = np.sum(abs(approximate_forward(candidates, data) - data["y_critical"]) ** 2, axis=1) / NOISE_RADIUS ** 2
    minor_scores = approximate_minor_scores(candidates, data)
    candidate_indices = np.where((direct_score > 1.01) & (minor_scores < .99))[0]
    if not len(candidate_indices):
        return {"status": "No gap found in the deterministic tangential diagnostic"}
    idx = candidate_indices[np.argmin(minor_scores[candidate_indices])]
    z = candidates[idx]
    for power in range(18, 34):
        rad = math.ldexp(1., -power)
        box = [arb(float(x), rad) for x in z]
        physical = score_material(box, operators, observed) / RADIUS ** 2
        minors = certified_minor_scores(box, operators, observed, data["B"])
        if physical.lower() > 1 and all(value.upper() < 1 for value in minors):
            return {
                "status": "Nonempty material cube is excluded by full noise test and passes all 15 separate chart tests",
                "material_center_binary64": z.tolist(),
                "material_half_width": rad,
                "complete_noise_score_enclosure": outward_number(physical),
                "largest_minor_noise_score_upper": max(outward_number(a)[1] for a in minors),
                "all_minor_noise_score_enclosures": [outward_number(a) for a in minors],
                "candidate_count": len(candidates),
                "approximate_eligible_count": len(candidate_indices),
                "seed": 7206,
            }
    return {"status": "Numerical candidate found; attempted material cubes not certified"}


def inner_cubes(witnesses, operators, observed):
    records = {}
    for name, details in witnesses.items():
        seed = details["material"]
        for power in range(15, 40):
            half_width = math.ldexp(1., -power)
            score = score_material([arb(float(x), half_width) for x in seed], operators, observed) / RADIUS ** 2
            if score.upper() < 1:
                enclosure = outward_number(score)
                enclosure[0] = max(0., enclosure[0])
                records[name] = {
                    "material_center_binary64": seed,
                    "material_half_width": half_width,
                    "noise_score_enclosure": enclosure,
                    "material_volume_exact_dyadic": (2 * half_width) ** 6,
                }
                break
    for a, b in combinations(records.values(), 2):
        if np.linalg.norm(
            np.array(a["material_center_binary64"]) - b["material_center_binary64"], ord=np.inf
        ) <= a["material_half_width"] + b["material_half_width"]:
            raise ArithmeticError("Certified inner cubes must be disjoint for the volume sum")
    return {
        "cubes": records,
        "physical_volume_lower": math.fsum(record["material_volume_exact_dyadic"] for record in records.values()),
    }


def score_polynomials(z, operators, observed):
    H, B, R, host = operators
    A = [[acb(int(i == j)) - H[i, j] * z[j] for j in range(6)] for i in range(6)]
    D = acb_mat(A).det()
    tail = [-R[j] * z[j] for j in range(6)]
    P = [acb_mat([
        A[i] + [B[i, a]] for i in range(6)
    ] + [tail + [host[a] - observed[a]]]).det() for a in range(3)]
    return sum((squared_modulus(x) for x in P), arb(0)), RADIUS ** 2 * squared_modulus(D)


def prove_segment(a, b, operators, observed, nodes, bernstein_inverse):
    start = [arb(float(q)) for q in a]
    delta = [arb(float(q)) - p for q, p in zip(b, start, strict=True)]
    grid = [score_polynomials([p + t * d for p, d in zip(start, delta, strict=True)], operators, observed) for t in nodes]
    controls = bernstein_inverse * arb_mat(grid)
    low = min(controls[i, 1].lower() for i in range(len(grid)))
    return max(controls[i, 0].upper() for i in range(len(grid))) / low if low > 0 else arb("inf")


def prove_bridge(data, operators, observed, witnesses):
    center = data["zc"]
    branches = [(key, raw["material"]) for key, raw in witnesses.items() if key.startswith("point_fiber_minus")]
    if len(branches) != 2:
        return {"status": "Two certified minus readings are needed for this diagnostic"}
    material_axis = int(np.argmax(abs(np.linalg.svd(seed_real_jacobian(data))[2][-1])))
    free_axes = [i for i in range(6) if i != material_axis]
    transverse = np.linalg.svd(seed_real_jacobian(data))[0][:, :5].T
    degree = 12
    nodes = [arb(i) / degree for i in range(degree + 1)]
    inverse = arb_mat([
        [math.comb(degree, j) * t ** j * (1 - t) ** (degree - j) for j in range(degree + 1)]
        for t in nodes
    ]).inv()
    output = {}
    for branch, endpoint in branches:
        end = np.array(endpoint)
        for segments in (8, 16, 32, 64):
            vertex = [center]
            for fraction in np.linspace(0., 1., segments + 1)[1:-1]:
                trial = center + fraction * (end - center)

                def residual(free_values):
                    q = trial.copy()
                    q[free_axes] = free_values
                    err = approximate_forward(q, data) - data["y_critical"]
                    return transverse @ np.stack((err.real, err.imag), axis=1).reshape(6)

                solved = root(residual, trial[free_axes], method="hybr", options={"xtol": 1e-11})
                trial[free_axes] = solved.x
                vertex.append(trial)
            vertex.append(end)
            bounds = [prove_segment(a, b, operators, observed, nodes, inverse) for a, b in zip(vertex[:-1], vertex[1:], strict=True)]
            if all(b.upper() < 1 for b in bounds):
                output[branch] = {
                    "status": "certified",
                    "vertex_binary64": [z.tolist() for z in vertex],
                    "piece_count": segments,
                    "maximum_noise_score_upper_on_segments": max(outward_number(b)[1] for b in bounds),
                }
                break
        if branch not in output:
            output[branch] = {"status": "segments unresolved"}
    return {
        "status": "both branches certified" if all(path["status"] == "certified" for path in output.values()) else "incomplete",
        "material_axis_fixed_in_transverse_solve": material_axis,
        "polynomial_degree_for_each_rational_score_enclosure": degree,
        "paths_from_public_fold_seed": output,
    }


def root_box_in_cube(cube, bounds):
    radius = Fraction(cube["material_half_width"])
    for q, (lo, hi) in zip(cube["material_center_binary64"], bounds, strict=True):
        lower, upper = Fraction(float(q)) - radius, Fraction(float(q)) + radius
        if not lower < Fraction(float(lo)) <= Fraction(float(hi)) < upper:
            return False
    return True


def link_certified_roots(cubes, bridge, scope_bounds):
    fold = json.loads((RESULTS / "fold_certificate.json").read_text())
    fibers = json.loads((RESULTS / "point_fibers.json").read_text())["fibers"]
    critical = [[x["lo"], x["hi"]] for x in fold["critical_material_outer_box"]]
    included = {"critical_root": root_box_in_cube(cubes["public_fold_seed"], critical)}
    for i, root in enumerate(fibers["minus"]["physical_roots"], 1):
        key = f"point_fiber_minus_{i}"
        included[key] = root_box_in_cube(cubes[key], root["material_intervals"])
    box = np.asarray(scope_bounds)
    in_local = all(
        np.all(q > box[:, 0]) and np.all(q < box[:, 1])
        for path in bridge["paths_from_public_fold_seed"].values() for q in path.get("vertex_binary64", [])
    )
    cubes_in_local = all(
        Fraction(float(lo)) < Fraction(float(q)) - Fraction(cube["material_half_width"])
        and Fraction(float(q)) + Fraction(cube["material_half_width"]) < Fraction(float(hi))
        for cube in cubes.values()
        for q, (lo, hi) in zip(cube["material_center_binary64"], box, strict=True)
    )
    bridge["certificate_root_boxes_within_inner_cubes"] = included
    bridge["all_cubes_and_bridge_pieces_within_local_crop"] = bool(in_local and cubes_in_local)
    bridge["all_three_certified_roots_in_same_noisy_path_component"] = bool(
        bridge["status"] == "both branches certified" and all(included.values()) and in_local and cubes_in_local
    )


def main():
    ctx.prec = PRECISION_BITS
    start_cpu, start_wall = time.process_time(), time.perf_counter()
    data, operators = load_inputs()
    observed = [exact_complex(y) for y in data["y_critical"]]
    cover = json.loads((RESULTS / "noisy_cover.json").read_text())
    witnesses = cover["shared_point_witnesses"]
    result = {
        "precision_bits": PRECISION_BITS,
        "inner_cubes": inner_cubes(witnesses, operators, observed),
        "relaxation_gap": find_gap(data, operators, observed),
        "bridge": prove_bridge(data, operators, observed, witnesses),
    }
    link_certified_roots(result["inner_cubes"]["cubes"], result["bridge"],
                         cover["scopes"]["local"]["material_box_binary64"])
    result["cpu_seconds"] = time.process_time() - start_cpu
    result["wall_seconds"] = time.perf_counter() - start_wall
    (RESULTS / "noisy_diagnostics.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print("Inner cubes", len(result["inner_cubes"]["cubes"]), "score gap", result["relaxation_gap"]["status"], "bridge", result["bridge"]["status"], flush=True)


if __name__ == "__main__":
    main()
