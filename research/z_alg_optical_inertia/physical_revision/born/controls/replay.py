#!/usr/bin/env python3
"""Paired native-grid weak-medium and exactly Born-compatible controls."""
from __future__ import annotations

from contextlib import contextmanager
import json
from pathlib import Path
import platform
import sys
from time import perf_counter, process_time

HERE = Path(__file__).resolve().parent
BORN = HERE.parent
sys.path[:0] = [str(BORN), str(BORN.parent / "modes")]
import run as bound

from flint import acb, acb_mat, arb, ctx
import numpy as np
from scipy.linalg import solve
from threadpoolctl import threadpool_limits

from analyze import Geometry, green_audit, physical_currents
from core.certificate import ball, frobenius_squared, safe_float

SEED = 20261124183
ARMS = ("ultraweak_ellipse_vie", "modest_paired_disks_vie", "modest_paired_disks_born")
NOISE = (0.01, 0.02)


@contextmanager
def clock(timing, label):
    wall, cpu = perf_counter(), process_time()
    yield
    timing[label] = dict(wall_seconds=perf_counter()-wall, cpu_seconds=process_time()-cpu)


def profiles(xy):
    x, y = xy.T
    q = ((x + .21) / .37)**2 + ((y - .16) / .24)**2
    ellipse = 1e-6 * np.maximum(1 - q, 0.)
    disk_a = (x + .29)**2 + (y + .12)**2 <= .32**2
    disk_b = (x - .31)**2 + (y - .19)**2 <= .22**2
    disks = .04 * np.maximum(disk_a, .6*disk_b)
    return ellipse, disks


def noisy_fields(clean, arm):
    rng = np.random.default_rng(np.random.SeedSequence([SEED, arm]))
    noise = rng.standard_normal(clean.shape) + 1j*rng.standard_normal(clean.shape)
    noise *= np.linalg.norm(clean) / np.linalg.norm(noise)
    for fraction in NOISE:
        y = clean + fraction * noise
        radius = np.linalg.norm(y-clean) * (1 + 8 * y.size * np.finfo(float).eps)
        radius = np.nextafter(radius, np.inf) + 1e-6*np.linalg.norm(y)
        yield fraction, y, float(radius)


def prepare_bound(s, b):
    complex_a = np.einsum("pi,ij->pji", s, b).reshape(s.shape[0]*b.shape[1], s.shape[1])
    a = np.ascontiguousarray(np.r_[complex_a.real, complex_a.imag])
    gram = a @ a.T
    interior = solve(gram + 1e-9*np.linalg.norm(gram, 2)*np.eye(len(gram)),
                     -a.sum(axis=1), assume_a="pos")
    margins = a.T @ interior
    bs = [bound.mat(s[:, i:i+384].conj().T) for i in range(0, len(b), 384)]
    bb = [bound.mat(b[i:i+384].conj()) for i in range(0, len(b), 384)]
    return a, interior, margins, bs, bb


def inverse_action_error(g, contrast, current_fft, s, b, symbol_norm_upper):
    chi = [arb(float(z)) for z in contrast]
    q = arb(float(contrast.max())) * symbol_norm_upper
    denom = 1-q
    if not denom > 0:
        raise ArithmeticError(f"No strict contraction bound: {q!s}")
    rest = b + g.exact_green(current_fft)
    for i, value in enumerate(chi):
        for j in range(rest.ncols()):
            rest[i, j] *= value
    defect = frobenius_squared(current_fft-rest).sqrt()
    # In the supported pixel space the cell operator is a compression of a
    # unitary FFT symbol. The Frobenius bound on S is sufficient for this radius.
    field_error = frobenius_squared(s).sqrt() * defect / denom
    return q.upper(), field_error.upper()


def diagonal_times(field, diagonal):
    return acb_mat(
        field.nrows(), field.ncols(),
        [value * field[i,j] for i, value in enumerate(diagonal) for j in range(field.ncols())],
    )


def residual_upper(y, prediction):
    return safe_float(frobenius_squared(y-prediction).sqrt(), upper=True)


def explicit_fit(y, born, peak, radius, cap, rigorous_y, rigorous_born):
    ray = float(np.vdot(born, y).real / np.vdot(born, born).real)
    alpha = max(0., ray) if cap is None else float(np.clip(ray, 0., cap/peak))
    if cap is not None and not (arb(alpha)*arb(peak)).upper() <= cap:
        raise ArithmeticError("Fallback least-squares coefficient exceeds the proposed cap")
    upper = residual_upper(rigorous_y, rigorous_born * arb(alpha))
    return dict(real_nonnegative_scalar=alpha, fitted_peak=alpha*peak,
                residual_upper=upper, residual_upper_over_radius=upper/radius,
                within_radius=upper <= radius)


def main():
    timing, numeric, raw = {}, {}, {}
    start = perf_counter(), process_time()
    ctx.prec, ctx.threads = 80, 1
    with clock(timing, "prepare_inputs"):
        with np.load(bound.STUDY / "synthetic/public/fine_geometry.npz", allow_pickle=False) as f:
            spec, s, b = [f[k] for k in ("spectrum", "PW360_10x20__S", "PW360_10x20__B")]
        g = Geometry(spec, s, b)
        xy, _, diagonal, _ = green_audit(g)
        weak, modest = profiles(xy)
        comparisons = []
        s_exact, b_exact = ball(s), ball(b)
        symbol_norm = max(abs(acb(complex(z))).upper() for z in spec.ravel())
    with clock(timing, "native_full_wave_and_interval_check"):
        for arm, chi in zip(ARMS[:2], (weak, modest), strict=True):
            j, iterations, defects = physical_currents(g, chi, diagonal)
            j_exact = ball(j)
            trial_full = s_exact * j_exact
            q_upper, tail = inverse_action_error(g, chi, j_exact, s_exact, b_exact, symbol_norm)
            born_exact = s_exact * diagonal_times(b_exact, [arb(float(z)) for z in chi])
            born = s @ (chi[:, None] * b)
            full = s @ j
            generator_clean_error = safe_float(
                frobenius_squared(ball(full)-trial_full).sqrt() + tail, upper=True)
            true_full_to_born = frobenius_squared(trial_full-born_exact).sqrt()
            numeric[arm] = dict(
                chi=chi, full=full, born=born, born_exact=born_exact, full_exact_trial=trial_full,
                full_exact_tail=tail,
                fit_evidence=dict(
                    peak_contrast=float(chi.max()), support_area_m2=float(np.mean(chi > 0)*4),
                    current_max_relative_defect=float(defects.max()), current_gmres_iterations=iterations,
                    uniform_symbol_norm_upper=safe_float(symbol_norm, upper=True),
                    material_contraction_q_upper=safe_float(q_upper, upper=True),
                    generated_full_field_absolute_error_upper=generator_clean_error,
                    born_vs_full_clean_relative_error=float(np.linalg.norm(born-full)/np.linalg.norm(full)),
                    born_vs_full_clean_absolute_error_interval=[
                        max(0., safe_float(true_full_to_born-tail, upper=False)),
                        safe_float(true_full_to_born+tail, upper=True),
                    ],
                ),
            )
            raw.update({f"{arm}__chi": chi.reshape(96,96), f"{arm}__clean": full,
                        f"{arm}__born": born, f"{arm}__physical_current": j})
        born_arm = numeric[ARMS[1]].copy()
        born_arm["fit_evidence"] = dict(born_arm["fit_evidence"])
        born_arm["fit_evidence"]["generated_born_field_absolute_error_upper"] = residual_upper(
            ball(born_arm["born"]), born_arm["born_exact"])
        numeric[ARMS[2]] = born_arm
    with clock(timing, "prepare_original_born_test"):
        setup = prepare_bound(s, b)
        matrix, interior, margins, bs, bb = setup
    # No generator-specific profile, shape, or fit is passed to the exclusion test.
    with clock(timing, "data_only_born_test"):
        for index, arm in enumerate(ARMS, 1):
            data = numeric[arm]
            clean = data["born"] if index == 3 else data["full"]
            for fraction, y, radius in noisy_fields(clean, index):
                noise_tag = f"pct{int(100 * fraction)}"
                name = f"{arm}__{noise_tag}"
                row = dict(arm=arm, observation=name, nominal_noise_fraction=fraction,
                           radius=radius, noisy_to_generated_clean=float(np.linalg.norm(y-clean)),
                           clean_field_frobenius=float(np.linalg.norm(clean)), tests={})
                raw[name+"__Y"] = y
                for label, cap in (("cone",None),("box4",4.)):
                    method_timing = {}
                    with clock(method_timing, "direction"):
                        estimate, vector = bound.polar_fit(
                            matrix, bound.realvector(y), interior, margins, cap, 64)
                        direction = bound.complexmatrix(vector, y.shape)
                    with clock(method_timing, "interval"):
                        lower, maxg = bound.interval_bounds(direction, y, bs, bb, cap)
                    raw[name+"__"+label+"_dual"] = direction
                    row["tests"][label] = dict(
                        search_estimate=float(estimate), lower=float(lower), lower_over_radius=float(lower/radius),
                        max_cell_pair_upper=float(maxg), direct_test_excludes=bool(lower>radius),
                        timing=method_timing,
                    )
                comparisons.append(row)
    with clock(timing, "explicit_offline_upper_fits"):
        for row in comparisons:
            arm, name, radius = (row[k] for k in ("arm", "observation", "radius"))
            data = numeric[arm]
            y = raw[name+"__Y"]
            y_exact = ball(y)
            peak = data["fit_evidence"]["peak_contrast"]
            row["known_born_model_residual_upper_over_radius"] = (
                residual_upper(y_exact, data["born_exact"]) / radius)
            row["known_full_model_residual_upper_over_radius"] = safe_float(
                frobenius_squared(y_exact-data["full_exact_trial"]).sqrt()+data["full_exact_tail"], upper=True
            ) / radius
            for label, cap in (("cone", None), ("box4", 4.)):
                fit = explicit_fit(y, data["born"], peak, radius, cap, y_exact, data["born_exact"])
                test = row["tests"][label]
                direct = test["direct_test_excludes"]
                by_subset = label == "box4" and row["tests"]["cone"]["direct_test_excludes"]
                if (direct or by_subset) and fit["within_radius"]:
                    raise ArithmeticError(f"Inconsistent interval result on {name}")
                state = ("excluded" if direct or by_subset else
                         "explicit_Born_fit_in_ball" if fit["within_radius"] else "undecided")
                test.update(fit=fit, decision=state, exclusion_inherited_from_cone=by_subset,
                            strongest_proved_lower_over_radius=max(
                                test["lower_over_radius"],
                                row["tests"]["cone"]["lower_over_radius"] if by_subset else 0.
                            ))
            print(name, [(label, d["decision"], d["lower_over_radius"],
                          d["fit"]["residual_upper_over_radius"])
                         for label,d in row["tests"].items()], flush=True)
    with clock(timing, "write_raw_arrays"):
        np.savez_compressed(HERE / "arrays.npz", **raw)
    result = dict(seed=SEED, setup=dict(
        grid=96, width_m=2., frequency_hz=400e6, illuminations=10, receivers=20,
        independent_shapes=2, arms=3, paired_noise_fractions=NOISE, bound_iterations=64,
        precision_bits=ctx.prec, controls_have_same_fixed_G_S_B_as_original=True,
        tuning_on_control_outcomes=False, host=platform.platform(),
    ), design={arm: data["fit_evidence"] for arm,data in numeric.items()},
        observations=comparisons, timing=timing, wall_seconds=perf_counter()-start[0],
        cpu_seconds=process_time()-start[1], native_actions=dict(g.cost))
    (HERE / "results.json").write_text(json.dumps(result, allow_nan=False, indent=2) + "\n")


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
