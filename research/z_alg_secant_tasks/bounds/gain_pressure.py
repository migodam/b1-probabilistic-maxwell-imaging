"""Calibrated and unknown shared-complex-gain checks on saved task witnesses."""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

from refine_caps import CASE, ROOT, arb, binary, ctx, np, real_norm, round_down, up
from box_modulus import acb
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(ROOT / "faces"))
from face_modulus import FrozenVIE, acb_channels, exact_face


ARMS = ("base", "angle", "frequency", "rich")
GAIN_LIMITS = ("0.01", "0.05")


def interval(value):
    return [round_down(value), up(value)]


def inner(x, y):
    return sum(a.conjugate() * b for a, b in zip(x, y, strict=True))


def data_mid(values):
    return np.array([complex(x) for x in values])


def choose_feasible_gain(change, limit):
    """Emit a dyadic point slightly inside the allowed complex gain disc."""
    radius = float(limit) * (1 - 64 * np.finfo(float).eps)
    if abs(change) > radius:
        change *= radius / abs(change)
    gain = 1 + change
    out = acb(float(gain.real), float(gain.imag))
    if not abs(out - 1) < limit:
        raise ArithmeticError("The proposed gain failed the directional-rounding disc check.")
    return out


def verify_gain(gain, xs, y, delta, limit):
    distance = [real_norm(gain * xj - yj for xj, yj in zip(x, y, strict=True)) for x in xs]
    return {
        "gain": [float(gain.real), float(gain.imag)], "gain_deviation": interval(abs(gain - 1)),
        "within_allowed_gain_disc": bool(abs(gain - 1) <= limit),
        "residual_over_delta": [interval(d / delta) for d in distance],
        "all_noise_compatible": bool(abs(gain - 1) <= limit and all(d < delta for d in distance)),
    }


def single_fit(x, y, delta, limit):
    L = inner(x, x).real
    difference = y - x
    correlation = inner(x, difference)
    change = complex(correlation / L)
    gain = choose_feasible_gain(change, limit)
    result = verify_gain(gain, [x], y, delta, limit)

    # The true residual minimum over the complex gain disc has this closed form.
    excess = abs(correlation) / L - limit
    nonnegative_excess = excess.nonnegative_part()
    sq = inner(difference, difference).real - abs(correlation)**2 / L + L * nonnegative_excess**2
    best_norm = sq.nonnegative_part().sqrt()
    result["exact_best_residual_over_delta"] = interval(best_norm / delta)
    result["candidate_excluded_for_all_permitted_gains"] = bool(best_norm > delta)
    return result


def shared_fit(xs, y, delta, limit):
    # Minimize t with both squared residuals <= t delta^2. Each residual
    # constraint and the gain disc are convex in real and imaginary gain.
    radius = float(limit)
    denominator = float(delta)
    fields = np.array([data_mid(x) for x in xs]) / denominator
    shift = fields - data_mid(y)[None, :] / denominator

    def residuals(v):
        return shift + radius * complex(v[0], v[1]) * fields

    def constraints(v):
        squared = np.sum(np.abs(residuals(v))**2, axis=1)
        return np.r_[1 - v[0]**2 - v[1]**2, v[2] - squared]

    def jacobian(v):
        corr = 2 * radius * np.sum(fields.conj() * residuals(v), axis=1)
        return np.vstack(([-2 * v[0], -2 * v[1], 0], np.column_stack(
            (-corr.real, -corr.imag, np.ones(len(xs))))))

    initial_t = float(np.max(np.sum(np.abs(shift)**2, axis=1)))
    optimum = minimize(lambda v: v[2], [0., 0., initial_t],
                       jac=lambda v: np.array([0., 0., 1.]),
                       constraints={"type": "ineq", "fun": constraints, "jac": jacobian},
                       method="SLSQP", options={"ftol": 2e-13, "maxiter": 80})
    gain = choose_feasible_gain(radius * complex(*optimum.x[:2]), limit)
    checked = verify_gain(gain, xs, y, delta, limit)
    checked.update({"floating_optimizer_status": int(optimum.status),
                    "floating_optimizer_iterations": int(optimum.nit)})
    return checked


def cap_point(row, config, center, radii, mask, outer_plus):
    h = radii.copy()
    eta = arb(row["positive_task_cap_eta"])
    concentrated = config["concentrated_roi_index"]
    if concentrated is None:
        h[mask] -= eta * arb("0.99")
    else:
        h[concentrated] -= eta * int(mask.sum()) * arb(149) / 150
    h[~mask] = np.array([binary(x) for x in outer_plus], dtype=object) - center[~mask]
    if any(a < -r or a > r for a, r in zip(h, radii, strict=True)):
        raise ArithmeticError("An archived cap point fails the specification box.")
    return center + h


def material_candidates(arm, center, radii, mask, face_input, caps):
    prefix = f"{CASE}__0.75__{arm}__"
    result = {}
    for label, sign in (("minus", -1), ("plus", 1)):
        outer_options = {"roi_only": np.array([float(c) for c in center[~mask]])}
        for name in ("LS", "GN"):
            k = f"{prefix}{name}_outside_{sign:+}"
            if k in face_input:
                outer_options[name] = face_input[k]
        for method, outer in outer_options.items():
            values = np.array([float(x) for x in center])
            values[~mask] = outer
            result[f"{label}_{method}_extreme"] = (sign, exact_face(values, center, radii, mask, sign))
    for row in caps:
        if row["arm"] != arm or row["positive_task_cap_eta"] == "0":
            continue
        eta = row["positive_task_cap_eta"]
        options = [row["witness_replay"], *row["concentrated_witness_replays"]]
        outer = face_input[f"{prefix}GN_outside_+1"]
        for item in options:
            label = item["concentrated_roi_index"]
            suffix = "balanced" if label is None else f"slot{label}"
            result[f"plus_cap{eta}_{suffix}"] = (1, cap_point(row, item, center, radii, mask, outer))
    return result


def cap_margin(rows, arm, y, delta):
    proven = next((row for row in rows if row["arm"] == arm
                   and row["positive_task_cap_eta"] == "0.015"
                   and row["status"] == "all_leaves_excluded"), None)
    if proven is None:
        return None
    certified_floor = binary(proven["minimum_lower_bound_all_leaves"])
    mid_norm = real_norm(y)
    excess = certified_floor - delta
    safe_limit = excess / (certified_floor + mid_norm)
    return {
        "known_gain_cap_residual_floor": float(certified_floor),
        "known_gain_cap_residual_floor_minus_radius": interval(excess),
        "reference_full_scattered_norm": interval(mid_norm),
        "guaranteed_gain_radius_threshold_open": interval(safe_limit),
        "guaranteed_gain_radius_threshold_percent_open": interval(100 * safe_limit),
        "analytic_lower_bound_with_unknown_gain": {
            test: interval((1 - arb(test)) * certified_floor - arb(test) * mid_norm)
            for test in GAIN_LIMITS
        },
    }


def run(output):
    start_wall, start_cpu = perf_counter(), process_time()
    reference = json.loads((ROOT / "array/results/summary.json").read_text())
    archive_caps = json.loads((ROOT / "bounds/refined_caps.json").read_text())
    caps = archive_caps["rows"]
    c = np.array(next(row["material_center"] for row in reference["cases"] if row["id"] == CASE))
    z = np.array([binary(v) for v in c], dtype=object)
    costs = Counter()
    results = []
    with np.load(ROOT / "array/results/physical_arrays.npz") as source, np.load(
            ROOT / "faces/face_parameters.npz") as face_points:
        mask = source[f"{CASE}__roi_mask"]
        radii = np.array([binary(r) for r in archive_caps["halfwidths"]], dtype=object)
        n_roi = int(mask.sum())
        prior = 2 * sum(radii[mask]) / n_roi
        model = FrozenVIE(source, len(c), costs)
        center_fields = model.certify(z, "reference")
        for arm in ARMS:
            y = acb_channels(center_fields, arm)
            delta = binary(next(row["radius"] for row in reference["pairs"]
                                if row["case"] == CASE and row["arm"] == arm))
            choices = material_candidates(arm, z, radii, mask, face_points, caps)
            xs, materials = {}, {}
            for label, (sign, state) in choices.items():
                fields = model.certify(state, "candidate")
                xs[label] = acb_channels(fields, arm)
                task_change = (sum(state[mask]) - sum(z[mask])) / n_roi
                endpoint_gap = sum(radii[mask]) / n_roi - sign * task_change
                materials[label] = {
                    "side": sign, "material": [interval(v) for v in state],
                    "task_change_relative_to_center": interval(task_change),
                    "task_gap_from_extreme": interval(endpoint_gap),
                    "nominal_gain_residual_over_delta": interval(real_norm(xs[label] - y) / delta),
                    "gains": {},
                }
            separate = {}
            shared = {}
            for gain in GAIN_LIMITS:
                limit = arb(gain)
                for label, d in materials.items():
                    d["gains"][gain] = single_fit(xs[label], y, delta, limit)

                checked_pairs = []
                for minus, (sminus, zm) in choices.items():
                    if sminus != -1:
                        continue
                    for plus, (splus, zp) in choices.items():
                        if splus != 1:
                            continue
                        task_gap = sum(zp[mask] - zm[mask]) / n_roi
                        # If either point failed all possible gains, a shared gain
                        # for this exact material pair is impossible.
                        if any(materials[label]["gains"][gain]["candidate_excluded_for_all_permitted_gains"]
                               for label in (minus, plus)):
                            check = {"all_noise_compatible": False, "reason": "single_point_excluded_for_all_gains"}
                        else:
                            check = shared_fit([xs[minus], xs[plus]], y, delta, limit)
                        check.update({"minus": minus, "plus": plus, "task_gap": interval(task_gap)})
                        checked_pairs.append(check)
                certified = [d for d in checked_pairs if d["all_noise_compatible"]]
                individual = [d for d in checked_pairs if all(
                    materials[label]["gains"][gain]["all_noise_compatible"]
                    for label in (d["minus"], d["plus"])
                )]
                best_individual = max(individual, key=lambda d: d["task_gap"][0], default=None)
                best_shared = max(certified, key=lambda d: (d["task_gap"][0],
                    -max(bound[1] for bound in d["residual_over_delta"])), default=None)
                separate[gain] = None if best_individual is None else {
                    "minus": best_individual["minus"], "plus": best_individual["plus"],
                    "task_gap": best_individual["task_gap"],
                    "distinct_material_hypotheses_both_noise_compatible": True,
                    "minus_gain": materials[best_individual["minus"]]["gains"][gain]["gain"],
                    "plus_gain": materials[best_individual["plus"]]["gains"][gain]["gain"],
                }
                shared[gain] = {"best_verified_pair": best_shared, "all_tested_pairs": checked_pairs}
            results.append({
                "arm": arm, "delta": float(delta), "reference_full_scattered_norm": interval(real_norm(y)),
                "prior_task_diameter": interval(prior), "positive_cap_guarantee": cap_margin(caps, arm, y, delta),
                "candidates": materials, "best_verified_pair_independent_gains": separate,
                "shared_gain_results": shared,
            })
            print(arm, [[gain,
                         None if separate[gain] is None else separate[gain]["task_gap"],
                         None if shared[gain]["best_verified_pair"] is None
                            else shared[gain]["best_verified_pair"]["task_gap"]] for gain in GAIN_LIMITS], flush=True)
    result = {
        "case": CASE, "scope": "archived coefficients, fixed nominal observation F_arm(m0), global scalar complex gain per acquisition",
        "gain_multiplies": "full target-plus-host scattered field before comparison; direct incident path absent",
        "gain_radii": GAIN_LIMITS, "wall_seconds": perf_counter() - start_wall,
        "cpu_seconds": process_time() - start_cpu, "precision_bits": ctx.prec,
        "flint_threads": ctx.threads, "full_model_cost": dict(costs), "arms": results,
    }
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, default=ROOT / "bounds/gain_pressure.json")
    args = cli.parse_args()
    ctx.threads = 1
    with threadpool_limits(limits=1, user_api="blas"):
        run(args.output)
