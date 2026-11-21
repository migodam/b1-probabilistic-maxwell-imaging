"""Replicate the historical regular pair in a fully coupled outgoing disk model."""

import os
for var in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[var] = "1"

import argparse
from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.special import hankel1, jv
from threadpoolctl import threadpool_limits


HERE = Path(__file__).resolve().parent
LEGACY = HERE.parents[1] / "delegated/z_alg_r5_delivery/host_publication/research/z_alg_host_fold/results"
FIXED = HERE / "historical_inputs.json"
DISTANCES = (20, 200, 2000, 20000, 200000)
SIZES = (2, 3)
TOL = 3e-14
ITERATION_LIMIT = 40
HALVING_LIMIT = 18


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def freeze_inputs(historical, output):
    """Copy small, immutable inputs once; historical roots are preserved as seeds."""
    old = json.loads((historical / "point_fibers.json").read_text())
    fiber = old["fibers"]["minus"]
    if fiber["legal_count"] != 2:
        raise ValueError("The historical observation must carry two certified legal roots.")
    roots = sorted(enumerate(fiber["physical_roots"]), key=lambda item: item[1]["material_midpoint"][0])
    with np.load(historical / "host_operators.npz", allow_pickle=False) as arrays:
        source = {
            "reference_url": "https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/z_alg_host_fold/",
            "archived_npz": "host_operators.npz",
            "archived_root_fiber": "point_fibers.json / fibers.minus.physical_roots",
            "geometry_decimal": {
                key: [[str(float(x)) for x in coord] for coord in arrays[key]]
                for key in ("centers", "sources", "receivers")
            },
            "disk_radius_decimal": str(float(arrays["radius"])),
            "known_material_contrast_decimal": [str(float(v)) for v in arrays["m0"]],
            "active_indices": arrays["active"].tolist(),
            "legal_material_range": arrays["material_box"].tolist(),
            "local_observation_binary64": [[float(v.real), float(v.imag)] for v in arrays["y_minus"]],
            "component_seed_A": roots[0][1]["material_midpoint"],
            "component_seed_B": roots[1][1]["material_midpoint"],
            "component_seed_A_index_in_archive": roots[0][0],
            "component_seed_B_index_in_archive": roots[1][0],
            "historical_root_records": [item for _, item in roots],
        }
    write_json(output, source)
    return source


def pair_distance(local_x, local_y, component_x, component_y, spacing):
    delta = local_x[:, None, :] - local_y[None, :, :]
    delta[:, :, 0] += spacing * (component_x[:, None] - component_y[None, :])
    return np.hypot(delta[:, :, 0], delta[:, :, 1])


def l2(x):
    return float(np.linalg.norm(x))


def pack(x):
    return np.concatenate((x.real, x.imag), axis=0)


def encoded(x):
    return [[[float(a.real), float(a.imag)] for a in row] for row in x]


class FullModel:
    """Every source sees every disk; every receiver evaluates every source."""

    def __init__(self, parameters, components, R):
        self.components = components
        self.spacing = R
        self.cost = Counter()
        self.xy = {}
        local = {}
        unit_labels = {}
        for key, value in parameters["geometry_decimal"].items():
            sub = np.asarray(value, dtype=float)
            local[key] = np.tile(sub, (components, 1))
            unit_labels[key] = np.repeat(np.arange(components), len(sub))
            self.xy[key] = local[key].copy()
            self.xy[key][:, 0] += R * unit_labels[key]
        labels = unit_labels["centers"]
        point_to_disks = pair_distance(local["centers"], local["centers"], labels, labels, R)
        disks_to_sources = pair_distance(local["centers"], local["sources"], labels, unit_labels["sources"], R)
        receiver_to_disks = pair_distance(local["receivers"], local["centers"], unit_labels["receivers"], labels, R)
        self.fixed = np.tile(np.asarray(parameters["known_material_contrast_decimal"], dtype=float), components)
        self.active = np.array([j + i*len(self.fixed)//components for i in range(components)
                                for j in parameters["active_indices"]])
        self.lo, self.hi = parameters["legal_material_range"]
        radius = float(parameters["disk_radius_decimal"])
        k = 2*np.pi
        diameter = 2 * radius
        self.id = np.eye(len(self.fixed))
        nonsingular = ~self.id.astype(bool)
        if min(point_to_disks[nonsingular]) <= diameter or min(disks_to_sources.ravel()) <= radius \
                or min(receiver_to_disks.ravel()) <= radius:
            raise ValueError("Physical geometry is invalid.")
        u = 0.5j * np.pi * k * radius * jv(1, k*radius)
        self.G = np.empty_like(point_to_disks, dtype=complex)
        self.G[nonsingular] = u * hankel1(0, k*point_to_disks[nonsingular])
        np.fill_diagonal(self.G, 0.5j * np.pi * k * radius * hankel1(1, k*radius) - 1)
        self.B = 0.25j * hankel1(0, k*disks_to_sources)
        self.S = u * hankel1(0, k*receiver_to_disks)
        self.local = np.arange(components)[:, None] == unit_labels["sources"]
        self.target = np.tile(np.array([complex(*p) for p in parameters["local_observation_binary64"]]), components)
        block_offdiag = labels[:, None] != labels[None, :]
        cross_disks_source = labels[:, None] != unit_labels["sources"][None, :]
        cross_receiver_disks = unit_labels["receivers"][:, None] != labels[None, :]

        def minimum_gap(distances, selected, width):
            entries = distances[selected]
            return float(np.min(entries) - width) if len(entries) else None

        self.geometry = {
            "components": components, "spacing_reference_wavelengths": R,
            "minimum_cross_disk_edge_gap": minimum_gap(point_to_disks, block_offdiag, diameter),
            "minimum_disk_to_cross_source_gap": minimum_gap(disks_to_sources, cross_disks_source, radius),
            "minimum_disk_to_cross_receiver_gap": minimum_gap(receiver_to_disks, cross_receiver_disks, radius),
            "maximum_receiver_x_extent": float(np.ptp(self.xy["receivers"][:, 0])),
            "off_block_green_frobenius": l2(self.G[block_offdiag]),
        }

    def solve(self, z, differentiate=False):
        m = self.fixed.copy()
        m[self.active] = z
        A = self.id - m[:, None] * self.G
        lu = lu_factor(A)
        self.cost["full_state_lu"] += 1
        RHS = m[:, None] * self.B
        currents = lu_solve(lu, RHS)
        self.cost["full_source_rhs"] += RHS.shape[1]
        y = self.S @ currents
        f = pack(y[self.local] - self.target)
        if not differentiate:
            return y, f
        E = (self.B + self.G @ currents)[self.active]
        P = lu_solve(lu, self.S.T, trans=1).T[:, self.active]
        self.cost["full_receiver_transpose_rhs"] += self.components
        J = pack((P[:, None, :] * E.T[None, :, :])[self.local])
        return y, f, J, A, RHS, currents


def newton(model, initial):
    start, cpu = perf_counter(), process_time()
    before = model.cost.copy()
    z = initial.copy()
    y, f, J, A, RHS, currents = model.solve(z, differentiate=True)
    sigma_at_seed = float(np.linalg.svd(J, compute_uv=False)[-1])
    initial_norm, initial_ports = l2(f), y.copy()
    trace = []
    outcome = "newton_step_limit"
    for i in range(ITERATION_LIMIT + 1):
        residual = l2(f)
        if residual <= TOL:
            outcome = "numerical_root"
            break
        if i == ITERATION_LIMIT:
            break
        try:
            update = np.linalg.solve(J, -f)
            model.cost["newton_small_real_lu"] += 1
        except np.linalg.LinAlgError:
            outcome = "singular_newton_matrix"
            break
        row = {"residual_before": residual, "unscaled_step_l2": l2(update),
               "outside_box_trials": 0, "nondecreasing_trials": 0}
        for halvings in range(HALVING_LIMIT + 1):
            scale = 2.0**(-halvings)
            trial = z + scale * update
            if np.any((trial <= model.lo) | (trial >= model.hi)):
                row["outside_box_trials"] += 1
                continue
            _, trial_res = model.solve(trial)
            if l2(trial_res) >= residual:
                row["nondecreasing_trials"] += 1
                continue
            z = trial
            y, f, J, A, RHS, currents = model.solve(z, differentiate=True)
            row["accepted_step_scale"] = scale
            break
        else:
            trace.append(row)
            outcome = "box_or_descent_backtracking_limit"
            break
        trace.append(row)
    sv = np.linalg.svd(J, compute_uv=False)
    overall_diagnostics = {
        "status": outcome, "material_at_seed": initial.tolist(), "material_at_stop": z.tolist(),
        "initial_local_residual": initial_norm, "initial_local_residual_times_R": model.spacing*initial_norm,
        "final_local_residual": l2(f), "J_initial_sigma_min": sigma_at_seed,
        "J_final_sigma_min": float(sv[-1]), "J_final_sigma_max": float(sv[0]),
        "min_material_clearance": float(min(np.min(z - model.lo), np.min(model.hi - z))),
        "state_relative_residual": l2(A @ currents - RHS)/l2(RHS),
        "state_sigma_min": float(np.linalg.svd(A, compute_uv=False)[-1]),
        "material_displacement_from_seed_l2": l2(z-initial),
        "source_receiver_ports_order": "receiver-major complex",
        "all_port_response_at_seed": encoded(initial_ports), "all_port_response_at_stop": encoded(y),
        "newton_trace": trace, "cost": dict(model.cost-before),
        "wall_seconds": perf_counter()-start, "cpu_seconds": process_time()-cpu,
    }
    return z, y, J, overall_diagnostics


def root_gaps(roots, local):
    gaps = []
    for (first, z, y), (second, cz, cy) in combinations(roots, 2):
        gaps.append({"first": first, "second": second,
                     "material_l2": l2(z-cz), "local_port_l2": l2((y-cy)[local]),
                     "cross_port_l2": l2((y-cy)[~local])})
    return gaps


def run(parameters, output, sizes=SIZES, distances=DISTANCES):
    started, cpu = perf_counter(), process_time()
    output.mkdir(exist_ok=True, parents=True)
    seeds = {name: np.asarray(parameters[f"component_seed_{name}"]) for name in "AB"}
    arrays, rows = {}, []
    for components in sizes:
        for R in distances:
            current = FullModel(parameters, components, R)
            label = f"L{components}_R{R}"
            arrays[f"{label}__local_selection"] = current.local
            arrays[f"{label}__target_local"] = current.target
            for kind, xy in current.xy.items():
                arrays[f"{label}__{kind}"] = xy
            roots, attempt = [], []
            for signs in product("AB", repeat=components):
                key = "".join(signs)
                z, y, J, row = newton(current, np.concatenate([seeds[s] for s in signs]))
                row["branch_seed"] = key
                row["nearest_historical_branches_at_stop"] = "".join(
                    min(seeds, key=lambda k: l2(zi-seeds[k])) for zi in np.split(z, components))
                attempt.append(row)
                if row["status"] == "numerical_root":
                    roots.append((key, z, y))
                prefix = f"{label}__{key}"
                arrays[f"{prefix}__material"] = z
                arrays[f"{prefix}__all_ports"] = y
                arrays[f"{prefix}__local_jacobian"] = J
            gaps = root_gaps(roots, current.local)
            unique = []
            for name, material, _ in roots:
                if not any(l2(material-old) <= 1e-6 for _, old in unique):
                    unique.append((name, material))
            rows.append({"id": label, **current.geometry, "root_attempts": attempt,
                         "converged_nominees": len(roots), "distinct_nominees": len(unique),
                         "unique_branch_representatives": [name for name, _ in unique],
                         "pairwise_gaps_between_converged_nominees": gaps,
                         "cost": dict(current.cost)})
            print(label, len(roots), len(unique), f"{min(a['initial_local_residual'] for a in attempt):.4e}",
                  f"{max(a['initial_local_residual'] for a in attempt):.4e}", flush=True)
    slopes = []
    for L in sizes:
        series = [row for row in rows if row["components"] == L][-3:]
        if len(series) != 3:
            continue
        spacing = np.log([row["spacing_reference_wavelengths"] for row in series])
        for together in zip(*(row["root_attempts"] for row in series), strict=True):
            slopes.append({"components": L, "seed": together[0]["branch_seed"],
                           "fitted_slope_initial_local_residual_vs_R": float(np.polyfit(
                               spacing, np.log([item["initial_local_residual"] for item in together]), 1)[0])})
        slopes.append({"components": L, "fitted_slope_off_block_Green_vs_R": float(np.polyfit(
            spacing, np.log([row["off_block_green_frobenius"] for row in series]), 1)[0])})
    np.savez_compressed(output/"materials_and_ports.npz", **arrays)
    write_json(output/"numeric_ladder.json", {
        "physical_model": "historical 2D outgoing scalar-TM circular pulse, full all-to-all matrix",
        "distance_unit": "background wavelength at fixed k0=2pi", "cadence": distances,
        "components_translated_along_x_axis": sizes,
        "independent_shots": "each component's 3 sources are switched on separately",
        "observation": "same historical y_minus, exactly the saved binary64 data, tiled onto each local receiver",
        "cross_ports": "all are forward-calculated and left out of the solve",
        "stationary_test": "binary64 residual <= 3e-14, all materials strictly legal",
        "duplicate_report_threshold_material_l2": 1e-6, "rows": rows, "last_three_distance_log_fits": slopes,
        "wall_seconds": perf_counter()-started, "cpu_seconds": process_time()-cpu,
    })


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--inputs", type=Path, default=FIXED)
    cli.add_argument("--historical", type=Path, default=LEGACY)
    cli.add_argument("--output", type=Path, default=HERE/"results")
    cli.add_argument("--sizes", nargs="+", type=int, default=SIZES, choices=SIZES)
    cli.add_argument("--spacings", nargs="+", type=int, default=DISTANCES)
    args = cli.parse_args()
    inputs = json.loads(args.inputs.read_text()) if args.inputs.exists() else freeze_inputs(args.historical, args.inputs)
    with threadpool_limits(limits=1, user_api="blas"):
        run(inputs, args.output, args.sizes, args.spacings)
