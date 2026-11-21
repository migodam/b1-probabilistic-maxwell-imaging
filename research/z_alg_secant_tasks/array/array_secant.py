"""Deterministic finite-chord experiment in a known outgoing disk-pulse host."""

import os

for thread_var in ("VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                   "OMP_NUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[thread_var] = "1"

import argparse
import csv
from dataclasses import dataclass
import json
from pathlib import Path
import platform
from time import perf_counter, process_time

import numpy as np
from scipy import __version__ as scipy_version
from scipy.linalg import lu_factor, lu_solve
from scipy.special import hankel1, jv
from threadpoolctl import threadpool_info, threadpool_limits


HERE = Path(__file__).resolve().parent
LOWER, UPPER = 0.5, 30.0
RADIUS = 0.04
K0 = 2 * np.pi
FREQUENCY_RATIO = 1.2
SOURCE_ANGLES = (158, 116, 72, 32)
RECEIVER_ANGLES = (102, 54, 14, 77, 129, 146)
ACTIVE_CENTERS = (
    (-0.290, -0.301), (0.094, 0.099), (-0.112, 0.306), (0.323, -0.106),
    (-0.281, 0.111), (0.103, -0.296), (-0.101, -0.108), (0.307, 0.294),
    (-0.301, -0.096), (0.312, 0.116), (-0.102, -0.311), (0.083, 0.298),
    (-0.292, 0.302), (0.099, -0.104), (-0.116, 0.104), (0.320, -0.293),
)
SIDE_X = (-0.53, -0.43, 0.49, 0.59)
SIDE_Y = (-0.30, -0.15, 0.00, 0.15, 0.30)
RAIL_X = (-0.32, -0.17, -0.02, 0.13, 0.28)
RAIL_Y = (-0.43, 0.43)

ARMS = ("base", "angle", "frequency", "rich")
CHANNEL_COUNT = dict(zip(ARMS, (3, 6, 6, 24), strict=True))
LEVELS = (("quarter", 0.25), ("half", 0.5), ("match", 1.0))
SEEDS = {
    (8, "low"): 2026112181, (8, "high"): 2026112182,
    (12, "low"): 20261121121, (12, "high"): 20261121122,
    (16, "low"): 20261121161, (16, "high"): 20261121162,
}
EXPOSURE = 1 / np.sqrt(2)


@dataclass
class PhysicalCost:
    factorizations: int = 0
    forward_rhs: int = 0
    dressed_rhs: int = 0
    transpose_rhs: int = 0
    odd_resolvents: int = 0


@dataclass
class Dress:
    response: np.ndarray
    P: np.ndarray
    Gd: np.ndarray
    E: np.ndarray
    diagnostics: dict

    @property
    def jacobian(self):
        return self.P[:, None, :] * self.E.T[None, :, :]


def distances(x, y):
    return np.linalg.norm(x[:, None, :] - y[None, :, :], axis=2)


def polar_locations(radius, angles):
    theta = np.deg2rad(angles)
    return radius * np.column_stack((np.cos(theta), np.sin(theta)))


def known_host():
    side_centers = [(x, y) for x in SIDE_X for y in SIDE_Y]
    rail_centers = [(x, y) for y in RAIL_Y for x in RAIL_X]
    side_chi = [1.4 if (col + row) % 2 == 0 else 2.6
                for col in range(len(SIDE_X)) for row in range(len(SIDE_Y))]
    rail_chi = [1.4 if (col + row) % 2 == 0 else 2.6
                for row in range(len(RAIL_Y)) for col in range(len(RAIL_X))]
    return np.array(side_centers + rail_centers), np.array(side_chi + rail_chi)


class Model:
    """Full current model at one physical frequency, at unit source amplitude."""

    def __init__(self, active, k, sources, receivers, cost=None):
        host, self.host_chi = known_host()
        self.n = len(active)
        self.k = k
        self.centers = np.vstack((active, host))
        self.sources = sources
        self.receivers = receivers
        self.cost = PhysicalCost() if cost is None else cost
        self.identity = np.eye(len(self.centers))

        separation = distances(self.centers, self.centers)
        offdiag = ~self.identity.astype(bool)
        ds = distances(self.centers, sources)
        dr = distances(receivers, self.centers)
        self.min_geometry_gaps = {
            "disk_edge_to_disk_edge": float(np.min(separation[offdiag]) - 2 * RADIUS),
            "disk_edge_to_source": float(np.min(ds) - RADIUS),
            "disk_edge_to_receiver": float(np.min(dr) - RADIUS),
        }
        if min(self.min_geometry_gaps.values()) <= 0:
            raise ValueError("Disks overlap or meet a transducer.")
        # The diagonal is the analytic outgoing integral over the source disk.
        factor = 0.5j * np.pi * k * RADIUS * jv(1, k * RADIUS)
        self.G = np.empty_like(separation, dtype=complex)
        self.G[offdiag] = factor * hankel1(0, k * separation[offdiag])
        np.fill_diagonal(self.G, 0.5j * np.pi * k * RADIUS * hankel1(1, k * RADIUS) - 1)
        self.B = 0.25j * hankel1(0, k * ds)
        self.S = factor * hankel1(0, k * dr)

    def _state(self, z):
        material = np.concatenate((z, self.host_chi))
        current_system = self.identity - material[:, None] * self.G
        rhs = material[:, None] * self.B
        lu = lu_factor(current_system)
        self.cost.factorizations += 1
        self.cost.forward_rhs += rhs.shape[1]
        current = lu_solve(lu, rhs)
        singular = np.linalg.svd(current_system, compute_uv=False)
        diagnostics = {
            "current_relative_residual": float(np.linalg.norm(current_system @ current - rhs)
                                                / np.linalg.norm(rhs)),
            "sigma_min_current_system": float(singular[-1]),
            "condition_current_system": float(singular[0] / singular[-1]),
        }
        return self.S @ current, current, lu, diagnostics

    def state(self, z):
        response, _, _, diagnostics = self._state(z)
        return response, diagnostics

    def dress(self, z):
        response, current, lu, diagnostics = self._state(z)
        self.cost.dressed_rhs += self.n
        self.cost.transpose_rhs += len(self.receivers)
        current_from_active = lu_solve(lu, self.identity[:, :self.n])
        Gd = self.G[:self.n] @ current_from_active
        P = lu_solve(lu, self.S.T, trans=1).T[:, :self.n]
        E = self.B[:self.n] + self.G[:self.n] @ current
        return Dress(response, P, Gd, E, diagnostics)


def make_models(n, cost=None):
    if n not in (8, 12, 16):
        raise ValueError("The active disk count must be 8, 12, or 16.")
    active = np.asarray(ACTIVE_CENTERS[:n])
    sources = polar_locations(0.95, SOURCE_ANGLES)
    receivers = polar_locations(1.07, RECEIVER_ANGLES)
    return (Model(active, K0, sources, receivers, cost),
            Model(active, FREQUENCY_RATIO * K0, sources[:2], receivers[:2], cost))


def take_channels(response, arm):
    """Select a single unweighted channel block from receiver-by-source data."""
    if arm == "base":
        return response[[0, 1, 0], [0, 0, 1]]
    if arm == "angle":
        return response[[0, 1, 0], [2, 2, 3]]
    if arm == "rich":
        return response.swapaxes(0, 1).reshape((CHANNEL_COUNT["rich"], *response.shape[2:]))
    raise ValueError(f"No single-frequency channel block for {arm}.")


def split_observation(f0_response, f1_response):
    """Form the measured complex vectors; preserves trailing derivative dimensions."""
    base = take_channels(f0_response, "base")
    angle = take_channels(f0_response, "angle")
    frequency = take_channels(f1_response, "base")
    return {
        "base": base,
        "angle": EXPOSURE * np.concatenate((base, angle)),
        "frequency": EXPOSURE * np.concatenate((base, np.sqrt(FREQUENCY_RATIO) * frequency)),
        "rich": EXPOSURE * take_channels(f0_response, "rich"),
    }


def pack(z):
    return np.concatenate((z.real, z.imag), axis=0)


def null_projection(J, task):
    U, singular, Vh = np.linalg.svd(J, full_matrices=True)
    threshold = max(J.shape) * np.finfo(float).eps * singular[0]
    rank = int(np.count_nonzero(singular > threshold))
    row_space = Vh[:rank]
    coeff = row_space @ task
    invisible = task - row_space.T @ coeff
    return invisible, {
        "rank": rank,
        "singular_values": singular.tolist(),
        "sensitivity_to_null_relative": float(np.linalg.norm(invisible) / np.linalg.norm(task)),
        "linear_visible_task_noise_gain": float(np.linalg.norm(U[:, :rank] @ (coeff / singular[:rank]))),
    }


def paired_box_limit(center, direction):
    clearance = np.minimum(center - LOWER, UPPER - center)
    return float(np.min(clearance[direction != 0] / np.abs(direction[direction != 0])))


def odd_prediction(t, v, dressed, cost):
    K = dressed.Gd * v[None, :]
    PD = dressed.P * v[None, :]
    linear = PD @ dressed.E
    quadratic = PD @ K @ dressed.E
    cubic = PD @ K @ K @ dressed.E
    Z = np.eye(len(v)) - t**2 * (K @ K)
    cost.odd_resolvents += 1
    exact = 2 * t * PD @ np.linalg.solve(Z, dressed.E)
    return np.stack((2 * t * linear, t**2 * quadratic, 2 * t**3 * cubic, exact), axis=2)


def noise_vectors(seed):
    generator = np.random.default_rng(seed + 1_000_000)
    standard = generator.standard_normal((2, CHANNEL_COUNT["rich"]))
    complex_noise = standard[0] + 1j * standard[1]
    directions = {}
    for arm, count in CHANNEL_COUNT.items():
        prefix = complex_noise[:count]
        directions[arm] = prefix / np.linalg.norm(prefix)
    return directions


def norm(value):
    return float(np.linalg.norm(value))


def status(tested, reference):
    return bool(tested <= reference * (1 + 8 * np.finfo(float).eps))


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def run_case(n, strength, models, arrays, cost):
    start_wall, start_cpu = perf_counter(), process_time()
    cost_before = vars(cost).copy()
    seed = SEEDS[(n, strength)]
    case = f"disks{n:02d}_{strength}"
    center_range = (1.2, 3.6) if strength == "low" else (9., 23.)
    center = np.random.default_rng(seed).uniform(*center_range, n)
    mask = models[0].centers[:n, 0] < -0.13
    roi_task = mask / mask.sum()
    background_task = (~mask) / (~mask).sum()
    full_task = np.ones(n) / n
    reference = [model.dress(center) for model in models]
    host = [model.state(np.zeros(n))[0] for model in models]
    center_data = split_observation(*(dressed.response for dressed in reference))
    known_host_data = split_observation(*host)
    jacobian_complex = split_observation(*(dress.jacobian for dress in reference))
    per_channel_rms = 0.1 * norm(center_data["base"] - known_host_data["base"]) / np.sqrt(CHANNEL_COUNT["base"])
    base_noise_to_total = 0.1 * norm(center_data["base"]) / np.sqrt(CHANNEL_COUNT["base"])

    jacobians = {}
    for arm in ARMS:
        sigma = per_channel_rms * np.sqrt(CHANNEL_COUNT[arm])
        J = pack(jacobian_complex[arm])
        _, jacobian = null_projection(J, roi_task)
        jacobian["linear_visible_task_halfwidth"] = sigma * jacobian["linear_visible_task_noise_gain"]
        jacobian["linear_full_task_halfwidth"] = (
            jacobian["linear_visible_task_halfwidth"]
            if jacobian["sensitivity_to_null_relative"] < 1e-10 else None
        )
        jacobians[arm] = jacobian
        arrays[f"{case}__{arm}__J"] = J
        arrays[f"{case}__{arm}__center"] = center_data[arm]
        arrays[f"{case}__{arm}__host"] = known_host_data[arm]

    projected, _ = null_projection(pack(jacobian_complex["base"]), roi_task)
    if norm(projected) <= 20 * np.finfo(float).eps * norm(roi_task):
        raise ArithmeticError(f"No task gradient in the baseline null space for {case}.")
    directions = {"null": projected / norm(projected), "roi": roi_task / norm(roi_task)}
    null_limit = paired_box_limit(center, directions["null"])
    ordinary_limit = paired_box_limit(center, directions["roi"])
    amplitude = 0.9 * min(null_limit, ordinary_limit)
    radii = [
        (name, level, fraction * amplitude)
        for name in directions for level, fraction in LEVELS
    ]
    radii.append(("null", "self", 0.9 * null_limit))
    noise = noise_vectors(seed)

    arrays[f"{case}__material_center"] = center
    arrays[f"{case}__roi_mask"] = mask
    for name, vector in directions.items():
        arrays[f"{case}__direction_{name}"] = vector

    cases = []
    for name, level, t in radii:
        v = directions[name]
        pair = f"{case}__{name}_{level}"
        minus, plus = center - t * v, center + t * v
        if np.min([minus, plus]) < LOWER or np.max([minus, plus]) > UPPER:
            raise ArithmeticError(f"Illegal endpoint material for {pair}.")
        data_minus, diag_minus = zip(*(model.state(minus) for model in models), strict=True)
        data_plus, diag_plus = zip(*(model.state(plus) for model in models), strict=True)
        minus_observation = split_observation(*data_minus)
        plus_observation = split_observation(*data_plus)
        prediction = split_observation(*(odd_prediction(t, v, d, cost) for d in reference))
        arrays[f"{pair}__material_minus"] = minus
        arrays[f"{pair}__material_plus"] = plus
        for frequency, (yminus, yplus) in enumerate(zip(data_minus, data_plus, strict=True)):
            arrays[f"{pair}__f{frequency}_minus"] = yminus
            arrays[f"{pair}__f{frequency}_plus"] = yplus

        material_delta = plus - minus
        for arm in ARMS:
            prefix = f"{pair}__{arm}"
            yminus, yplus = minus_observation[arm], plus_observation[arm]
            c = center_data[arm]
            clean_delta = yplus - yminus
            midpoint = 0.5 * (yplus + yminus)
            linear, quadratic, cubic, exact = prediction[arm].T
            delta = per_channel_rms * np.sqrt(CHANNEL_COUNT[arm])
            self_delta = 0.1 * norm(c - known_host_data[arm])
            total_delta = base_noise_to_total * np.sqrt(CHANNEL_COUNT[arm])
            observed = yplus + delta * noise[arm]
            residual_mid = max(norm(yplus - midpoint), norm(yminus - midpoint))
            full_diff = norm(clean_delta)
            rms_delta = norm(material_delta) / np.sqrt(n)
            task_gap = float(roi_task @ material_delta)
            row = {
                "case": case, "n": n, "strength": strength, "direction": name,
                "level": level, "arm": arm, "amplitude_material_l2": t,
                "roi_minus": float(roi_task @ minus), "roi_plus": float(roi_task @ plus),
                "roi_task_gap": task_gap,
                "nonroi_mean_gap": float(background_task @ material_delta),
                "all_mean_gap": float(full_task @ material_delta),
                "material_rms_gap": rms_delta,
                "per_pair_minimax_roi_lb": abs(task_gap) / 2 if status(residual_mid, delta) else None,
                "per_pair_minimax_rms_lb": rms_delta / 2 if status(residual_mid, delta) else None,
                "roi_rms_material_gap": norm(material_delta[mask]) / np.sqrt(mask.sum()),
                "nonroi_rms_material_gap": norm(material_delta[~mask]) / np.sqrt((~mask).sum()),
                "full_wave_chord_norm": full_diff,
                "linear_chord_norm": norm(linear),
                "cubic_chord_norm": norm(cubic),
                "single_end_quadratic_norm": norm(quadratic),
                "single_end_even_actual_norm": norm(midpoint - c),
                "single_end_minus_center_norm": norm(yminus - c),
                "single_end_plus_center_norm": norm(yplus - c),
                "linear_prediction_rel_error": norm(clean_delta - linear) / full_diff,
                "cubic_prediction_rel_error": norm(clean_delta - linear - cubic) / full_diff,
                "identity_absolute_error": norm(clean_delta - exact),
                "identity_rel_error": norm(clean_delta - exact) / full_diff,
                "reference_total_scattered_norm": norm(c),
                "reference_host_subtracted_norm": norm(c - known_host_data[arm]),
                "radius": delta,
                "radius_local_acquisition": self_delta,
                "radius_baseline_total_scattered": total_delta,
                "chord_over_twice_radius": residual_mid / delta,
                "chord_over_twice_radius_local_acquisition": residual_mid / self_delta,
                "chord_over_twice_radius_baseline_total_scattered": residual_mid / total_delta,
                "even_over_radius": norm(midpoint - c) / delta,
                "common_midpoint": status(residual_mid, delta),
                "noisy_positive_residual_over_radius": norm(observed - yplus) / delta,
                "noisy_negative_residual_over_radius": norm(observed - yminus) / delta,
                "noisy_negative_compatible": status(norm(observed - yminus), delta),
                "minimum_endpoint_current_system_sigma": min(
                    diagnostic["sigma_min_current_system"] for diagnostic in (*diag_minus, *diag_plus)),
                "maximum_endpoint_state_relative_residual": max(
                    diagnostic["current_relative_residual"] for diagnostic in (*diag_minus, *diag_plus)),
            }
            cases.append(row)
            arrays[f"{prefix}__minus"] = yminus
            arrays[f"{prefix}__plus"] = yplus
            arrays[f"{prefix}__midpoint"] = midpoint
            arrays[f"{prefix}__noisy_positive"] = observed
            arrays[f"{prefix}__predictions"] = prediction[arm]

    return {
        "id": case, "n": n, "strength": strength, "seed": seed,
        "material_center": center.tolist(),
        "center_allowed_drawing_interval": center_range,
        "roi_indices": np.flatnonzero(mask).tolist(),
        "jacobians": jacobians, "null_box_limit": null_limit,
        "ordinary_box_limit": ordinary_limit, "matched_amplitude": amplitude,
        "center_states": [d.diagnostics for d in reference],
        "physical_cost": {name: val - cost_before[name] for name, val in vars(cost).items()},
        "wall_seconds": perf_counter() - start_wall, "cpu_seconds": process_time() - start_cpu,
    }, cases


def run(output):
    output.mkdir(parents=True, exist_ok=True)
    arrays, cases, pair_records, geometries = {}, [], [], {}
    cost = PhysicalCost()
    start_wall, start_cpu = perf_counter(), process_time()
    for n in (8, 12, 16):
        models = make_models(n, cost)
        geom = f"disks{n:02d}"
        geometries[geom] = [model.min_geometry_gaps for model in models]
        arrays[f"{geom}__centers"] = models[0].centers
        arrays[f"{geom}__known_host_chi"] = models[0].host_chi
        for i, model in enumerate(models):
            arrays[f"{geom}__f{i}_sources"] = model.sources
            arrays[f"{geom}__f{i}_receivers"] = model.receivers
            arrays[f"{geom}__f{i}_G"] = model.G
            arrays[f"{geom}__f{i}_B"] = model.B
            arrays[f"{geom}__f{i}_S"] = model.S
        for strength in ("low", "high"):
            case, records = run_case(n, strength, models, arrays, cost)
            cases.append(case)
            pair_records.extend(records)

    summary = {
        "schema": "ZALG_HOST_ARRAY_SECANT_V1",
        "model": "2D outgoing scalar-TM disk pulse at full host state",
        "length_unit": "reference vacuum wavelength", "disk_radius": RADIUS,
        "chi_box": [LOWER, UPPER], "wavenumbers": [K0, FREQUENCY_RATIO * K0],
        "source_angles_deg": SOURCE_ANGLES, "receiver_angles_deg": RECEIVER_ANGLES,
        "channel_count_by_arm": CHANNEL_COUNT,
        "noise_reference": "10% of baseline host-subtracted center RMS per complex channel",
        "raw_prediction_columns": ["first_order_chord", "single_end_quadratic",
                                   "third_order_chord", "odd_resolvent_chord"],
        "platform": platform.platform(), "python": platform.python_version(),
        "numpy": np.__version__, "scipy": scipy_version, "threadpools": threadpool_info(),
        "blas_thread_environment": {var: os.environ[var] for var in (
            "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "OMP_NUM_THREADS", "BLIS_NUM_THREADS")},
        "wall_seconds_through_physical_diagnostics": perf_counter() - start_wall,
        "process_seconds_through_physical_diagnostics": process_time() - start_cpu,
        "physical_cost": vars(cost), "geometry_gaps": geometries,
        "cases": cases, "pairs": pair_records,
    }
    np.savez_compressed(output / "physical_arrays.npz", **arrays)
    save_json(output / "summary.json", summary)
    with (output / "pairs.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=pair_records[0].keys())
        writer.writeheader()
        writer.writerows(pair_records)
    print(f"{len(cases)} centers, {len(pair_records)} paired diagnostic rows, "
          f"{summary['wall_seconds_through_physical_diagnostics']:.3f} s wall, "
          f"{summary['process_seconds_through_physical_diagnostics']:.3f} s CPU")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "results")
    args = parser.parse_args()
    with threadpool_limits(limits=1, user_api="blas"):
        run(args.output)
