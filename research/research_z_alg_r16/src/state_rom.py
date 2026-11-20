"""FFT-TSOM current-state ROM and direct full-wave material inversion."""
import os

# Accelerate and OpenBLAS both read their settings before the first BLAS call.
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import lu_factor, lu_solve, qr, svd
from scipy.optimize import minimize

from cohorts import load_cohort

ROOT = Path(__file__).resolve().parents[1]
RANKS = (36, 80, 144)
RECEIVER_CUTOFF = 1e-2
BACKGROUND = 0.1
SMOOTH_WEIGHT = 1e-3
RIDGE_WEIGHT = 1e-5
REFRESH_AT = 36
EXTRA_RANK = 12
DEFAULT_ITERATIONS = 120
COHORTS = ("r10_twod", "r11_independent48", "r11_crossgrid8", "r11_bothblocks16")


def timestamp():
    return np.array([perf_counter(), process_time()])


def time_since(start):
    return timestamp() - start


def timing(time_pair):
    return dict(wall_s=float(time_pair[0]), cpu_s=float(time_pair[1]))


def fft_dictionary(shape, size):
    """Inverse-DFT atoms ordered by radial spatial frequency."""
    n_y, n_x = shape
    signed_y = np.rint(np.fft.fftfreq(n_y) * n_y).astype(int)
    signed_x = np.rint(np.fft.fftfreq(n_x) * n_x).astype(int)
    modes = sorted(
        ((int(y), int(x)) for y in signed_y for x in signed_x),
        key=lambda yx: (yx[0] ** 2 + yx[1] ** 2, yx[0], yx[1]),
    )[:size]
    spectra = np.zeros((n_y, n_x, size), dtype=np.complex128)
    for column, (y, x) in enumerate(modes):
        spectra[y % n_y, x % n_x, column] = 1
    return np.fft.ifft2(spectra, axes=(0, 1), norm="ortho").reshape(n_y * n_x, size)


def regularizer(material, shape):
    grid = material.reshape(shape)
    dy = np.diff(grid, axis=0)
    dx = np.diff(grid, axis=1)
    centered = grid - BACKGROUND
    w_s = SMOOTH_WEIGHT / material.size
    w_r = RIDGE_WEIGHT / material.size
    value = 0.5 * (w_s * (np.sum(dy**2) + np.sum(dx**2)) + w_r * np.sum(centered**2))
    gradient = w_r * centered
    gradient[:-1] -= w_s * dy
    gradient[1:] += w_s * dy
    gradient[:, :-1] -= w_s * dx
    gradient[:, 1:] += w_s * dx
    return value, gradient.ravel()


class CurrentModel:
    """Identity trial/test bases for Full; orthonormal Galerkin bases for ROM."""

    def __init__(self, G, S, E0, shape, U=None, GU=None, SU=None):
        self.G, self.S, self.E0 = G, S, E0
        self.shape = shape
        self.U = U
        self.n, self.ns = E0.shape
        if U is None:
            self.rank = self.n
            self.SU = S
        else:
            self.rank = U.shape[1]
            self.Uh = np.ascontiguousarray(U.conj().T)
            self.GU = G @ U if GU is None else GU
            self.SU = S @ U if SU is None else SU
        self.SUh = np.ascontiguousarray(self.SU.conj().T)

    def forward(self, material):
        m = material[:, None]
        rhs = m * self.E0
        if self.U is None:
            lhs = np.multiply(self.G, -m, order="F")
        else:
            lhs = np.asfortranarray(-(self.Uh @ (m * self.GU)))
            rhs = self.Uh @ rhs
        lhs.flat[:: self.rank + 1] += 1.0
        factor = lu_factor(lhs, overwrite_a=True, check_finite=False)
        z = lu_solve(factor, rhs, check_finite=False)
        return z, self.SU @ z, factor

    def field(self, z):
        return self.E0 + (self.G @ z if self.U is None else self.GU @ z)

    def lift(self, z):
        return z if self.U is None else self.U @ z

    def objective(self, material, Y, scale):
        z, prediction, factor = self.forward(material)
        residual = prediction - Y
        w = lu_solve(factor, self.SUh @ residual / scale, trans=2, check_finite=False)
        gradient = np.sum((self.lift(w).conj() * self.field(z)).real, axis=1)
        reg, reg_gradient = regularizer(material, self.shape)
        data_objective = 0.5 * np.vdot(residual, residual).real / scale
        return data_objective + reg, gradient + reg_gradient

    def enrich(self, material):
        z, _, _ = self.forward(material)
        residual = material[:, None] * self.field(z) - self.lift(z)
        complement = residual - self.U @ (self.Uh @ residual)
        candidate, singular, _ = svd(complement, full_matrices=False, check_finite=False)
        k = min(EXTRA_RANK, int(np.count_nonzero(singular > 1e-10 * singular[0])))
        new = candidate[:, :k]
        new -= self.U @ (self.Uh @ new)
        new, _ = qr(new, mode="economic", check_finite=False)
        updated = CurrentModel(
            self.G, self.S, self.E0, self.shape, np.column_stack([self.U, new]),
            np.column_stack([self.GU, self.G @ new]),
            np.column_stack([self.SU, self.S @ new]),
        )
        return updated, singular[:k].tolist()


def build_models(model):
    start = timestamp()
    G = np.asarray(model["G"], dtype=np.complex128, order="F")
    S = np.asarray(model["S"], dtype=np.complex128)
    E0 = np.asarray(model["E0"], dtype=np.complex128)
    shape = tuple(model["shape"])
    full = CurrentModel(G, S, E0, shape)
    common = time_since(start)
    models = {"full": full}
    preparation = {"full": common}
    start = timestamp()
    _, singular, Vh = svd(S, full_matrices=False, check_finite=False)
    kept = int(np.count_nonzero(singular > RECEIVER_CUTOFF * singular[0]))
    V = np.ascontiguousarray(Vh[:kept].conj().T)
    Vh = V.conj().T
    receiver_time = time_since(start)
    for rank in RANKS:
        start = timestamp()
        dictionary = fft_dictionary(shape, rank - kept)
        perpendicular = dictionary - V @ (Vh @ dictionary)
        Q, _ = qr(perpendicular, mode="economic", check_finite=False)
        key = f"tsom{rank}"
        models[key] = CurrentModel(G, S, E0, shape, np.column_stack([V, Q]))
        preparation[key] = common + receiver_time + time_since(start)
    start = timestamp()
    models["fft80"] = CurrentModel(G, S, E0, shape, fft_dictionary(shape, 80))
    preparation["fft80"] = common + time_since(start)
    return models, preparation, singular.tolist(), kept


def score_optimization(arm, observation, initial, bounds, max_iterations):
    scale = np.vdot(observation, observation).real
    trace = []
    state = None

    def objective(material):
        nonlocal state
        value, gradient = arm.objective(material, observation, scale)
        state = (value, gradient)
        return value, gradient

    def record(material):
        value, gradient = state
        lo, hi = bounds
        projected_gradient = material - np.clip(material - gradient, lo, hi)
        trace.append([len(trace) + 1, float(value), float(np.max(abs(projected_gradient)))])

    start = timestamp()
    result = minimize(
        objective, initial, jac=True, method="L-BFGS-B", bounds=[bounds] * initial.size,
        callback=record,
        options={"maxiter": max_iterations, "gtol": 1e-7, "ftol": 1e-11, "maxls": 40},
    )
    elapsed = time_since(start)
    return result.x, {
        "termination": str(result.message),
        "iterations": int(result.nit),
        "evaluations": int(result.nfev),
        "objective": float(result.fun),
        "rank": arm.rank,
        "timing": timing(elapsed),
        "trace": trace,
    }, elapsed


def terminal_readout(arm, full, material, observation):
    start = timestamp()
    z, prediction, _ = arm.forward(material)
    J = arm.lift(z)
    rhs = material[:, None] * arm.E0
    defect = J - material[:, None] * arm.field(z)
    defect_rel = np.linalg.norm(defect) / np.linalg.norm(rhs)
    rom_time = time_since(start)
    exact_time = np.zeros(2)
    exact_prediction = prediction
    if arm.U is not None:
        start = timestamp()
        _, exact_prediction, _ = full.forward(material)
        exact_time = time_since(start)
    scale = np.linalg.norm(observation)
    return {
        "model_data_residual": float(np.linalg.norm(prediction - observation) / scale),
        "true_data_residual": float(np.linalg.norm(exact_prediction - observation) / scale),
        "rom_to_full_data_difference": float(np.linalg.norm(prediction - exact_prediction) / scale),
        "state_defect_rel_rhs": float(defect_rel),
        "reduced_terminal_readout": timing(rom_time),
        "exact_terminal_resolve": timing(exact_time),
        "full_terminal_exact_solve": int(arm.U is not None),
    }, rom_time + exact_time


def recovered(arm, full, observation, bounds, max_iterations, setup, base=None):
    initial = np.full(arm.n, np.clip(BACKGROUND, *bounds))
    stages = []
    enrich_time = np.zeros(2)
    new_singular = []
    if base is None:
        material, stage, optimize_time = score_optimization(
            arm, observation, initial, bounds, max_iterations,
        )
        stages.append(stage)
    else:
        first_material, first_stage, first_time = base
        stages.append(first_stage)
        start = timestamp()
        arm, new_singular = arm.enrich(first_material)
        enrich_time = time_since(start)
        material, stage, second_time = score_optimization(
            arm, observation, first_material, bounds, max_iterations - REFRESH_AT,
        )
        stages.append(stage)
        optimize_time = first_time + second_time

    readout, terminal_time = terminal_readout(arm, full, material, observation)
    evals = int(sum(s["evaluations"] for s in stages))
    total_warm = optimize_time + enrich_time
    data = {
        "rank_per_illumination": arm.rank,
        "complex_current_unknowns_before_elimination": arm.rank * arm.ns,
        "real_current_dof_before_elimination": 2 * arm.rank * arm.ns,
        "real_material_unknowns": arm.n,
        "nonlinear_iterations": int(sum(s["iterations"] for s in stages)),
        "objective_evaluations": evals,
        "warm_factorizations": evals + int(base is not None),
        "warm_forward_rhs": arm.ns * (evals + int(base is not None)),
        "warm_adjoint_rhs": arm.ns * evals,
        "warm_optimize": timing(optimize_time),
        "warm_online_enrichment": timing(enrich_time),
        "warm_reconstruction": timing(total_warm),
        "shared_operator_and_method_setup": timing(setup),
        "cold_reconstruction": timing(total_warm + setup),
        "warm_including_terminal_validation": timing(total_warm + terminal_time),
        "cold_including_terminal_validation": timing(total_warm + setup + terminal_time),
        "online_enrichment_singular_values": new_singular,
        "stages": stages,
        **readout,
    }
    return material, data


def check_gradient(arm, observation, material=None):
    rng = np.random.default_rng(1616)
    if material is None:
        material = 1.1 + 0.1 * rng.standard_normal(arm.n)
    direction = rng.standard_normal(arm.n)
    direction /= np.linalg.norm(direction)
    scale = np.vdot(observation, observation).real
    _, gradient = arm.objective(material, observation, scale)
    step = 1e-4
    plus, _ = arm.objective(material + step * direction, observation, scale)
    minus, _ = arm.objective(material - step * direction, observation, scale)
    fd = (plus - minus) / (2 * step)
    via_adjoint = gradient @ direction
    return {
        "central_difference": float(fd),
        "adjoint_directional_derivative": float(via_adjoint),
        "abs_error": float(abs(fd - via_adjoint)),
        "rel_error": float(abs(fd - via_adjoint) / max(abs(fd), abs(via_adjoint), 1e-14)),
    }


def run_cohort(name, max_iterations, max_cases, output_tag):
    start = timestamp()
    model, cases = load_cohort(name)
    data_load = time_since(start)
    if max_cases is not None:
        cases = cases[:max_cases]
    bounds = tuple(map(float, model["bounds"]))
    models, setups, singular, kept = build_models(model)
    full = models["full"]
    case_results = {}
    materials = {}
    gradient_checks = {}
    for i, case in enumerate(cases):
        case_id = case["id"]
        Y = np.asarray(case["Y"], dtype=np.complex128)
        if i == 0:
            gradient_checks.update({key: check_gradient(arm, Y) for key, arm in models.items()})
        case_results[case_id] = {"sigma": float(case["sigma"]), "roles": {}}
        results = case_results[case_id]["roles"]
        for key, arm in models.items():
            materials[f"{case_id}::{key}"], results[key] = recovered(
                arm, full, Y, bounds, max_iterations, setups[key],
            )
        initial = np.full(full.n, np.clip(BACKGROUND, *bounds))
        base = score_optimization(models["tsom80"], Y, initial, bounds, REFRESH_AT)
        materials[f"{case_id}::tsom80_enriched"], results["tsom80_enriched"] = recovered(
            models["tsom80"], full, Y, bounds, max_iterations, setups["tsom80"], base=base,
        )
        if i == 0:
            enriched, _ = models["tsom80"].enrich(base[0])
            gradient_checks["tsom80_enriched"] = check_gradient(enriched, Y, base[0])
        print(name, case_id, "exact data", " ".join(
            f"{key}={value['true_data_residual']:.4f}" for key, value in results.items()
        ), flush=True)

    saved = {
        "cohort": name,
        "grid_shape": list(model["shape"]),
        "material_bounds": list(bounds),
        "illuminations": full.ns,
        "receivers": int(full.S.shape[0]),
        "receiver_cutoff_relative": RECEIVER_CUTOFF,
        "receiver_rank": kept,
        "receiver_singular_values": singular,
        "uniform_initial_material": float(np.clip(BACKGROUND, *bounds)),
        "smooth_weight_per_pixel": SMOOTH_WEIGHT,
        "ridge_weight_per_pixel": RIDGE_WEIGHT,
        "iterations_per_inversion": max_iterations,
        "enrichment_split_iteration": REFRESH_AT,
        "enrichment_max_added_rank": EXTRA_RANK,
        "shared_cohort_archive_load": timing(data_load),
        "gradient_checks": gradient_checks,
        "cases": case_results,
    }
    # Write material predictions before GT is accessed by the separate score action.
    out = ROOT / "results" / f"state_{name}{output_tag}"
    np.savez_compressed(out.with_suffix(".npz"), **materials)
    out.with_suffix(".json").write_text(json.dumps(saved, indent=2, allow_nan=False) + "\n")


def score_cohort(name, output_tag):
    from cohorts import load_truth
    out = ROOT / "results" / f"state_{name}{output_tag}"
    saved = json.loads(out.with_suffix(".json").read_text())
    with np.load(out.with_suffix(".npz")) as material:
        for case_id, case in saved["cases"].items():
            truth = np.asarray(load_truth(name, case_id)).ravel()
            target = truth > 1.0
            for key, result in case["roles"].items():
                prediction = material[f"{case_id}::{key}"]
                rmse = float(np.sqrt(np.mean((prediction - truth) ** 2)))
                mask = prediction > 1.0
                union = np.count_nonzero(target | mask)
                result["material_rmse"] = rmse
                result["material_amplitude_nrmse"] = float(rmse / (truth.max() - BACKGROUND))
                result["material_support_iou"] = float(np.count_nonzero(target & mask) / union if union else 1.0)
    out.with_suffix(".json").write_text(json.dumps(saved, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("run", "score"))
    parser.add_argument("cohorts", nargs="*")
    parser.add_argument("--iterations", type=int, default=DEFAULT_ITERATIONS)
    parser.add_argument("--max-cases", type=int)
    parser.add_argument("--tag", default="")
    args = parser.parse_args()
    if args.action == "run" and args.iterations <= REFRESH_AT:
        parser.error(f"--iterations must exceed {REFRESH_AT}")
    for cohort in args.cohorts or COHORTS:
        if args.action == "run":
            run_cohort(cohort, args.iterations, args.max_cases, args.tag)
        else:
            score_cohort(cohort, args.tag)
