"""Find finite task-changing material pairs and certify them at the physical ports."""

import argparse
from pathlib import Path
import sys
from time import perf_counter, process_time

HERE = Path(__file__).resolve().parent
ARRAY = HERE.parent / "array"
sys.path.insert(0, str(ARRAY))

from array_secant import (
    EXPOSURE, FREQUENCY_RATIO, LOWER, UPPER, make_models, norm, null_projection, pack,
    save_json, split_observation, take_channels,
)

import numpy as np
from scipy.linalg import lu_factor, lu_solve, qr
from threadpoolctl import threadpool_limits

from interval_physics import IntervalModel, certify_chord, interval_environment


FRACTIONS = (0.25, 0.5, 1.0)
TOLERANCE = 3e-14
STEP_LIMIT = 24
HALVING_LIMIT = 16


class NumericalModel:
    def __init__(self, model, angular=False, cost=None):
        n = model.n
        m = np.r_[np.zeros(n), model.host_chi]
        self.angular = angular
        self.source_count = 4 if angular else 2
        B, S = model.B[:, :self.source_count], model.S[:2]
        rhs = np.column_stack((m[:, None] * B, np.eye(len(m), n)))
        X = np.linalg.solve(np.eye(len(m)) - m[:, None] * model.G, rhs)
        host, V = X[:, :self.source_count], X[:, self.source_count:]
        self.H = model.G[:n] @ V
        self.B = B[:n] + model.G[:n] @ host
        self.R = S @ V
        self.host = S @ host
        self.identity = np.eye(n)
        self.cost = ({"lu": 0, "forward_rhs": 0, "transpose_rhs": 0, "correction_solves": 0}
                     if cost is None else cost)

    def complex_response(self, material):
        factorization = lu_factor(self.identity - material[:, None] * self.H)
        self.cost["lu"] += 1
        self.cost["forward_rhs"] += self.source_count
        w = lu_solve(factorization, material[:, None] * self.B)
        field = self.B + self.H @ w
        self.cost["transpose_rhs"] += 2
        observed = lu_solve(factorization, self.R.T, trans=1).T
        response = self.host + self.R @ w
        jacobian = observed[:, None, :] * field.T[None, :, :]
        options = ("base", "angle") if self.angular else ("base",)
        return tuple(np.concatenate([take_channels(x, kind) for kind in options], axis=0)
                     for x in (response, jacobian))

    def response(self, material):
        response, jacobian = self.complex_response(material)
        return response, pack(jacobian)


class NumericalAcquisition:
    def __init__(self, models, arm):
        self.arm = arm
        self.baseline = NumericalModel(models[0], angular=arm == "angle")
        self.cost = self.baseline.cost
        self.second = NumericalModel(models[1], cost=self.cost) if arm == "frequency" else None

    def response(self, material):
        response, jacobian = self.baseline.complex_response(material)
        if self.second is not None:
            more, derivative = self.second.complex_response(material)
            power = np.sqrt(FREQUENCY_RATIO)
            response = np.concatenate((response, power * more))
            jacobian = np.concatenate((jacobian, power * derivative))
        gain = 1 if self.arm == "base" else EXPOSURE
        return gain * response, gain * pack(jacobian)


def nuisance_basis(roi_mask, J):
    n = len(roi_mask)
    roi = np.flatnonzero(roi_mask)
    nonroi = np.flatnonzero(~roi_mask)
    coordinates = np.eye(n, dtype=int)
    candidate = np.column_stack((coordinates[:, nonroi],
                                 coordinates[:, roi[1:]] - coordinates[:, roi[0], None]))
    _, _, pivots = qr(J @ candidate, pivoting=True, mode="economic")
    data_dimension = J.shape[0]
    U = candidate[:, pivots[:data_dimension]]
    if np.any(U[roi_mask].sum(axis=0) != 0):
        raise ArithmeticError("The nuisance chart changes the ROI task.")
    singular = np.linalg.svd(J @ U, compute_uv=False)
    if U.shape[1] < data_dimension or singular[-1] <= data_dimension * np.finfo(float).eps * singular[0]:
        raise ArithmeticError("No nuisance chart with full data rank.")
    return U, {
        "roi_anchor": int(roi[0]), "selected_column_indices": pivots[:data_dimension].tolist(),
        "center_chart_jacobian_singular_values": singular.tolist(),
    }


def legal(z):
    return bool(np.all((LOWER < z) & (z < UPPER)))


def correction_solve(model, z, base, U, mode, center_y):
    started, c0 = perf_counter(), process_time()
    initial_cost = model.cost.copy()
    u = np.zeros(U.shape[1])

    def evaluate(x):
        h = base + U @ x
        plus = z + h
        if not legal(plus):
            return None
        yplus, Jplus = model.response(plus)
        if mode == "single":
            return pack(yplus - center_y), Jplus @ U
        yminus, Jminus = model.response(z - h)
        return pack((yplus - yminus) / 2), ((Jplus + Jminus) / 2) @ U

    trace = []
    f, J = evaluate(u)
    initial_residual = norm(f)
    initial_correction = None
    outcome = "step_limit"
    for iteration in range(STEP_LIMIT + 1):
        residual = norm(f)
        if residual <= TOLERANCE:
            outcome = "numerical_root"
            break
        if iteration == STEP_LIMIT:
            break
        try:
            step = np.linalg.solve(J, -f)
            model.cost["correction_solves"] += 1
        except np.linalg.LinAlgError:
            outcome = "singular_chart_at_iterate"
            break
        if iteration == 0:
            initial_correction = norm(U @ step)
        row = {"residual_before": residual, "material_step_l2_unscaled": norm(U @ step),
               "box_rejections": 0, "residual_rejections": 0}
        for halvings in range(HALVING_LIMIT + 1):
            next_u = u + 2.0**(-halvings) * step
            next_h = base + U @ next_u
            if not legal(z + next_h) or (mode == "double" and not legal(z - next_h)):
                row["box_rejections"] += 1
                continue
            next_f, next_J = evaluate(next_u)
            if norm(next_f) >= residual:
                row["residual_rejections"] += 1
                continue
            u, f, J = next_u, next_f, next_J
            row["accepted_step_fraction"] = 2.0**(-halvings)
            break
        else:
            trace.append(row)
            outcome = "box_backtracking_limit" if row["box_rejections"] else "residual_backtracking_limit"
            break
        trace.append(row)

    return u, J, {
        "status": outcome, "initial_residual": initial_residual,
        "final_residual": norm(f), "initial_linearized_correction_l2": initial_correction,
        "correction_material_l2": norm(U @ u),
        "chart_jacobian_sigma_min": float(np.linalg.svd(J, compute_uv=False)[-1]),
        "newton_trace": trace,
        "cost": {k: v - initial_cost[k] for k, v in model.cost.items()},
        "wall_seconds": perf_counter() - started, "cpu_seconds": process_time() - c0,
    }


def check_full(model, material, probability_radius, center_data):
    y, diagnostic = model.state(material)
    selected = take_channels(y, "base")
    return selected, {
        "material_min": float(np.min(material)), "material_max": float(np.max(material)),
        "material_legal": legal(material),
        "distance_to_F0": norm(selected - center_data),
        "distance_to_F0_over_10pct_radius": norm(selected - center_data) / probability_radius,
        "full_physics": diagnostic,
    }


def augmented_counterexample(source, output):
    """Fix F_arm(center) and the two prescribed signed task shifts."""
    start_wall, start_cpu = perf_counter(), process_time()
    arrays, cases = {}, []
    with np.load(source) as inputs:
        models = make_models(16)
        interval = IntervalModel(models[0].centers, 16, with_frequency=True)
        for mode in ("low", "high"):
            case = f"disks16_{mode}"
            z = inputs[f"{case}__material_center"]
            roi = inputs[f"{case}__roi_mask"]
            task = roi / roi.sum()
            full_ref = split_observation(*(item.state(z)[0] for item in models))
            for arm in ("base", "angle", "frequency"):
                key = f"{case}__{arm}"
                numeric = NumericalAcquisition(models, arm)
                y0, J0 = numeric.response(z)
                projected, local = null_projection(J0, task)
                v = projected / norm(projected)
                U, chart = nuisance_basis(roi, J0)
                arrays[f"{key}__direction"] = v
                arrays[f"{key}__U"] = U
                arrays[f"{key}__centre"] = z
                arrays[f"{key}__J0"] = J0
                arrays[f"{key}__full_F0"] = full_ref[arm]
                endpoints = []
                full_predictions = []
                materials = []
                for sign in (-1, 1):
                    target = "-0.10" if sign == -1 else "0.10"
                    t = float(target) / float(task @ v)
                    name = f"{key}__{'minus' if sign == -1 else 'plus'}"
                    if not legal(z + t * v):
                        endpoints.append({"id": name, "status": "prescribed_perturbation_outside_box"})
                        continue
                    u, J, result = correction_solve(numeric, z, t * v, U, "single", y0)
                    endpoint = z + t * v + U @ u
                    full_y = split_observation(*(item.state(endpoint)[0] for item in models))[arm]
                    material_audit = {"minimum": float(endpoint.min()), "maximum": float(endpoint.max()),
                                      "within_box": legal(endpoint)}
                    result.update({"id": name, "t": t, "exact_prescribed_task_offset": target,
                                   "full_response_error": norm(full_y - full_ref[arm]),
                                   "full_material_audit": material_audit})
                    if result["status"] == "numerical_root":
                        result["certificate"] = certify_chord(
                            interval, z, t, v, U, u, J, roi, arm=arm, single=True, task_target=target,
                        )
                    arrays[f"{name}__u"] = u
                    arrays[f"{name}__material"] = endpoint
                    arrays[f"{name}__full_response"] = full_y
                    full_predictions.append(full_y)
                    materials.append(endpoint)
                    endpoints.append(result)
                both = all(e.get("certificate", {}).get("status") == "certified_exact_root_in_chart"
                           for e in endpoints)
                cases.append({
                    "id": key, "arm": arm, "endpoints": endpoints, "both_endpoints_certified": both,
                    "certified_task_gap": "0.20" if both else None,
                    "full_state_response_gap": norm(full_predictions[0] - full_predictions[1])
                        if len(full_predictions) == 2 else None,
                    "full_material_rms_gap": norm(materials[0] - materials[1]) / np.sqrt(16)
                        if len(materials) == 2 else None,
                    "local": local, "chart": chart,
                })
                print(key, "both sides certified" if both else "two-sided certificate incomplete", flush=True)
    result = {"schema": "ZALG_HOST_AUGMENTED_FIXED_REFERENCE_V1",
              "interval_environment": interval_environment(), "cases": cases,
              "wall_seconds": perf_counter() - start_wall,
              "cpu_seconds": process_time() - start_cpu}
    save_json(output / "enhanced_counterexample.json", result)
    np.savez_compressed(output / "enhanced_counterexample_arrays.npz", **arrays)


def run(inputs, output):
    start_wall, start_cpu = perf_counter(), process_time()
    output.mkdir(exist_ok=True, parents=True)
    arrays = {}
    cases = []
    with np.load(inputs) as source:
        for n in (8, 12, 16):
            full, _ = make_models(n)
            numeric = NumericalModel(full)
            interval = IntervalModel(full.centers, n)
            for strength in ("low", "high"):
                key = f"disks{n:02d}_{strength}"
                z = source[f"{key}__material_center"]
                v = source[f"{key}__direction_null"]
                roi = source[f"{key}__roi_mask"]
                y0, J0 = numeric.response(z)
                full_y0 = take_channels(full.state(z)[0], "base")
                full_host = take_channels(full.state(np.zeros(n))[0], "base")
                delta = 0.1 * norm(full_y0 - full_host)
                u_basis, chart_diagnostic = nuisance_basis(roi, J0)
                q = roi / roi.sum()
                tmax = 0.9 * np.min(np.minimum(z - LOWER, UPPER - z) / np.abs(v))
                arrays[f"{key}__centre"] = z
                arrays[f"{key}__v"] = v
                arrays[f"{key}__U"] = u_basis
                arrays[f"{key}__J0"] = J0
                arrays[f"{key}__task"] = q
                arrays[f"{key}__full_F0"] = full_y0
                arrays[f"{key}__full_host"] = full_host
                arms = []
                for fraction in FRACTIONS:
                    t = float(fraction * tmax)
                    base = t * v
                    for mode in ("double", "single"):
                        name = f"{key}__{fraction:g}__{mode}"
                        u, J, result = correction_solve(numeric, z, base, u_basis, mode, y0)
                        h = base + u_basis @ u
                        plus, minus = z + h, z - h
                        yp, plus_check = check_full(full, plus, delta, full_y0)
                        ym, minus_check = check_full(full, minus, delta, full_y0)
                        result.update({
                            "id": name, "mode": mode, "fraction": fraction, "t": t,
                            "prescribed_roi_shift_single_end": float(t * (q @ v)),
                            "material_pair_roi_gap": float(q @ (plus - minus)),
                            "material_pair_rms_gap": norm(plus - minus) / np.sqrt(n),
                            "full_plus_minus_field_gap": norm(yp - ym),
                            "full_plus_minus_field_gap_over_2delta": norm(yp - ym) / (2 * delta),
                            "full_plus": plus_check, "full_minus": minus_check,
                        })
                        if mode == "double" and result["status"] == "numerical_root":
                            result["certificate"] = certify_chord(interval, z, t, v, u_basis, u, J, roi)
                        arrays[f"{name}__u"] = u
                        arrays[f"{name}__plus"] = plus
                        arrays[f"{name}__minus"] = minus
                        arrays[f"{name}__full_yplus"] = yp
                        arrays[f"{name}__full_yminus"] = ym
                        arms.append(result)
                cases.append({
                    "id": key, "n": n, "strength": strength, "asymmetric_chart": chart_diagnostic,
                    "full_vs_reduced_F0_gap": norm(full_y0 - y0), "noise_radius": delta,
                    "center_Jv_norm": norm(J0 @ v), "t_max": float(tmax), "arms": arms,
                })
                print(key, [(a["fraction"], a["mode"], a["status"],
                             a.get("certificate", {}).get("status")) for a in arms], flush=True)
    result = {
        "schema": "ZALG_HOST_FINITE_CHORD_CERT_V1", "source_initial_values": str(inputs.name),
        "interval_environment": interval_environment(), "normal_equation_tolerance": TOLERANCE,
        "newton_step_limit": STEP_LIMIT, "backtracking_halving_limit": HALVING_LIMIT,
        "cases": cases,
        "wall_seconds": perf_counter() - start_wall, "cpu_seconds": process_time() - start_cpu,
    }
    save_json(output / "finite_chord.json", result)
    np.savez_compressed(output / "finite_chord_arrays.npz", **arrays)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ARRAY / "results" / "physical_arrays.npz")
    parser.add_argument("--output", type=Path, default=HERE / "results")
    parser.add_argument("--enhanced", action="store_true",
                        help="Only certify the matched fixed-reference counterexamples after augmentation.")
    args = parser.parse_args()
    with threadpool_limits(limits=1, user_api="blas"):
        if args.enhanced:
            args.output.mkdir(exist_ok=True, parents=True)
            augmented_counterexample(args.input, args.output)
        else:
            run(args.input, args.output)
