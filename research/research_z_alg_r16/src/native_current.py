"""Short-budget rerun of the R10 native global-dual FFT-CIE loop."""
import os

os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
import json
import sys
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import qr, svd, svdvals

from cohorts import historical_study, load_cohort

ROOT = Path(__file__).resolve().parents[1]
CASE_COHORT = "r10_twod"
ARMS = ("native_fft", "small_fft", "small_swap")
ADDITION_ARMS = ("small_add_fourier", "small_add_residual")
NATIVE_MF = (3, 6, 10, 12)
SMALL_MF = (3, 4, 5, 6)
ITERATIONS = (24, 40, 56, 80)
SWAP_COUNT = 6


def load_native():
    """Use the same R8 implementation and R10 DPCS ranking called by R10."""
    sys.path.insert(0, str(historical_study(10).parents[1]))
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from research.research_max_struct_r8.src import cie_fft
    from research.research_max_struct_r8.src.benchmark_solver import arm_config
    from research.research_max_struct_r10.src.twod_adapter import noise_ball_rank
    return MaxwellModel, cie_fft, arm_config, noise_ball_rank


def times():
    return (perf_counter(), process_time())


def elapsed(start):
    return dict(zip(("wall_s", "cpu_s"), (float(x - t) for x, t in zip(times(), start))))


def radial_frequency(mode, shape):
    signed = [min(k, n - k) for k, n in zip(mode, shape)]
    return sum(k**2 for k in signed), mode


class HybridAPCS:
    """Retained native Fourier columns and added projected residual columns."""

    def __init__(self, fft, keep, extra, cost):
        self.fft, self.keep, self.extra, self.cost = fft, keep, extra, cost
        self.k = int(np.count_nonzero(keep)) + extra.shape[1]
        self.frequencies = np.vstack([fft.frequencies[keep], [(-1, i) for i in range(extra.shape[1])]])

    def apply(self, coefficients):
        packed = np.zeros((self.fft.k, coefficients.shape[1]), dtype=coefficients.dtype)
        size = np.count_nonzero(self.keep)
        packed[self.keep] = coefficients[:size]
        self.cost["residual_basis_apply_rhs"] += coefficients.shape[1]
        return self.fft.apply(packed) + self.extra @ coefficients[size:]

    def adjoint(self, field):
        self.cost["residual_basis_adjoint_rhs"] += field.shape[1]
        return np.vstack([self.fft.adjoint(field)[self.keep], self.extra.conj().T @ field])


class ResidualReplacement:
    """Build a fixed residual block after the common first native CIE stage."""

    def __init__(self, native_factory, model, Y, cfg, replace=True):
        self.native = native_factory
        self.model, self.Y, self.cfg = model, Y, cfg
        self.replace = replace
        self.anchor = None
        self.basis = None
        self.last_basis = None
        self.records = {}

    def record_first_stage(self, state):
        if state["stage"] == 0:
            self.anchor = (state["m"], state["J"])

    def __call__(self, shape, mf, Vd, cost):
        fft = self.native(shape, mf, Vd, cost)
        if mf == SMALL_MF[0]:
            return fft
        if self.basis is None:
            self.prepare(shape, fft, Vd, cost)
        keep = np.array([tuple(pair) not in self.drop for pair in fft.frequencies])
        self.last_basis = HybridAPCS(fft, keep, self.basis, cost)
        return self.last_basis

    def prepare(self, shape, second_fft, Vd, cost):
        start = times()
        material, J = self.anchor
        self.drop = set()
        if self.replace:
            first_fft = self.native(shape, SMALL_MF[0], Vd, cost)
            first_modes = {tuple(x) for x in first_fft.frequencies}
            new_modes = {tuple(x) for x in second_fft.frequencies} - first_modes
            ranked = sorted(new_modes, key=lambda mode: radial_frequency(mode, shape), reverse=True)
            self.drop = set(ranked[:SWAP_COUNT])
        final_fft = self.native(shape, SMALL_MF[-1], Vd, cost)
        keep = np.array([tuple(pair) not in self.drop for pair in final_fft.frequencies])
        dark = final_fft.apply(np.eye(final_fft.k, dtype=complex)[:, keep])
        orthogonal, _ = qr(np.column_stack([Vd, dark]), mode="economic", check_finite=False)

        G, S, E = self.model.G, self.model.S, self.model.E
        M = material[:, None]
        beta = float(self.cfg["stages"][1]["beta_current"])
        p = beta / (1 + beta * M)
        d = (S @ J - self.Y) * self.cfg["data_weights"]
        e = (p**2) * (J - M * (E + G @ J)) * self.cfg["state_weights"]
        gradient = S.conj().T @ d + self.cfg.get("rho", 1.0) * (e - G.conj().T @ (M * e))
        projected = gradient - orthogonal @ (orthogonal.conj().T @ gradient)
        left, singular, _ = svd(projected, full_matrices=False, check_finite=False)
        self.basis = left[:, :SWAP_COUNT]
        self.records = {
            "setup_timing": elapsed(start),
            "dropped_fourier_modes": [list(map(int, p)) for p in sorted(self.drop)],
            "candidate_singular_values": singular.tolist(),
            "retained_singular_energy": float(np.sum(singular[:SWAP_COUNT] ** 2) / np.sum(singular ** 2)),
            "receiver_overlap": float(np.linalg.norm(Vd.conj().T @ self.basis)),
            "retained_dark_overlap": float(np.linalg.norm(dark.conj().T @ self.basis)),
        }
        cost["basis_prep_dense_G_rhs"] += self.model.ns
        cost["basis_prep_dense_G_adjoint_rhs"] += self.model.ns
        cost["basis_prep_dense_S_rhs"] += self.model.ns
        cost["basis_prep_dense_S_adjoint_rhs"] += self.model.ns
        cost["basis_prep_retained_dark_QR"] += 1
        cost["basis_prep_gradient_SVD"] += 1
        cost["basis_prep_column_count"] += dark.shape[1]


class FourierAddition:
    """Keep the native mask and append the six nearest modes of the next mask."""

    def __init__(self, native_factory):
        self.native = native_factory
        self.last_basis = None
        self.previous_modes = set()
        self.records = {}

    def __call__(self, shape, mf, Vd, cost):
        native = self.native(shape, mf, Vd, cost)
        modes = {tuple(pair) for pair in native.frequencies}
        if mf != SMALL_MF[0]:
            expanded = self.native(shape, mf + 1, Vd, cost)
            candidates = {tuple(pair) for pair in expanded.frequencies} - modes
            added = sorted(candidates, key=lambda mode: radial_frequency(mode, shape))[:SWAP_COUNT]
            modes.update(added)
            self.records[str(mf)] = [list(map(int, p)) for p in added]
            expanded.frequencies = np.array(sorted(modes))
            expanded.k = len(modes)
            native = expanded
        if not self.previous_modes <= modes:
            raise RuntimeError("New Fourier basis must retain previous-stage modes")
        self.previous_modes = modes
        self.last_basis = native
        return native


def rank_readout(basis):
    start = times()
    singular = svdvals(basis.apply(np.eye(basis.k, dtype=complex)), check_finite=False)
    tol = np.finfo(float).eps * max(24 * 24, basis.k) * singular[0]
    return {
        "effective_dark_rank": int(np.count_nonzero(singular > tol)),
        "coordinates": basis.k,
        "smallest_singular_value": float(singular[-1]),
        "numerical_rank_tolerance": float(tol),
        "offline_rank_check_timing": elapsed(start),
    }


def cfg_for(arm_config, Y, model, model_data, dpcs_rank, arm):
    native = arm_config(
        "global_dual", {"host": .1, "material_box": model_data["bounds"], "dpcs_rank": dpcs_rank},
        Y, model.E, [{"mf": 12, "maxiter": 1500}],
    )
    mask = NATIVE_MF if arm == "native_fft" else SMALL_MF
    native["stages"] = [
        dict(stage, mf=mf, maxiter=maxiter)
        for stage, mf, maxiter in zip(native["stages"], mask, ITERATIONS)
    ]
    native.update(
        fft_geometry={key: model_data[key] for key in ("coords", "radius", "k")},
        regularizer=lambda m: model.reg.evaluate(m)[:2],
        physical_sigma=float(np.linalg.norm(Y)),
    )
    return native


def run_arm(solver, model_type, model_data, Y, dpcs_rank, arm, arm_config):
    start = times()
    model = model_type({**model_data, "grid_shape": model_data["shape"]})
    config = cfg_for(arm_config, Y, model, model_data, dpcs_rank, arm)
    initial = np.full(model.n, 0.1)
    factory = solver.NativeFourierAPCS
    replacement = None
    if arm in ("small_swap", "small_add_residual"):
        replacement = ResidualReplacement(factory, model, Y, config, replace=arm == "small_swap")
        config["stage_callback"] = replacement.record_first_stage
    elif arm == "small_add_fourier":
        replacement = FourierAddition(factory)
    if replacement is not None:
        solver.NativeFourierAPCS = replacement
    try:
        answer = solver.run_cie_fft(model, Y, initial, config)
    finally:
        solver.NativeFourierAPCS = factory
    total_time = elapsed(start)
    trace = answer["history"]
    selected = min(answer["candidates"], key=lambda row: row["objective"])
    stage_stats = [
        {
            "iterations": sum(row.get("stage") == i and "iteration" in row for row in trace),
            "mf": config["stages"][i]["mf"],
            "outer_budget": config["stages"][i]["maxiter"],
            "coordinate_rank": stage["basis_info"]["fourier_columns"],
            "true_data_residual": float(candidate["physical_relative_residual"]),
            "physical_objective": float(candidate["objective"]),
            "last_data_term": float(last["data_term"]),
            "last_cie_state_term": float(last["state_term_current"]),
        }
        for i, (stage, candidate, last) in enumerate(zip(
            answer["stage_states"], answer["candidates"][1:],
            [trace[index - 1] for index, row in enumerate(trace) if row.get("event") == "stage_end"],
        ))
    ]
    receipt = answer["execution_receipt"]
    record = {
        "timing": total_time,
        "solver_timing": {
            "wall_s": receipt["wall_seconds"],
            "cpu_s": receipt["cpu_seconds"],
        },
        "full_wave_candidate_lu": receipt["physical_LU"],
        "full_wave_candidate_rhs": receipt["physical_complex_RHS"],
        "selected_stage": int(selected["stage"]),
        "true_data_residual": float(selected["physical_relative_residual"]),
        "physical_objective": float(selected["objective"]),
        "stages": stage_stats,
        "operations": receipt["operation_counts"],
    }
    if replacement is not None:
        if arm == "small_add_fourier":
            record["additional_fourier_modes"] = replacement.records
        else:
            record["basis_replacement"] = replacement.records
        if arm in ADDITION_ARMS:
            record["final_basis_rank"] = rank_readout(replacement.last_basis)
    return answer["m"], record


def run(max_cases, tag, arms=ARMS):
    model_type, solver, arm_config, rank_rule = load_native()
    start = times()
    model_data, cases = load_cohort(CASE_COHORT)
    load_time = elapsed(start)
    records = {}
    predictions = {}
    for case in cases[:max_cases]:
        case_id, Y, sigma = (case[k] for k in ("id", "Y", "sigma"))
        start = times()
        rank = rank_rule(model_data["S"], Y, sigma)["rank"]
        setup = elapsed(start)
        case_results = {}
        records[case_id] = {
            "dpcs_rank": rank, "DPCS_rank_timing": setup, "sigma": sigma,
            "relative_rms_noise": float(sigma * np.sqrt(Y.size) / np.linalg.norm(Y)),
            "arms": case_results,
        }
        for arm in arms:
            prediction, record = run_arm(solver, model_type, model_data, Y, rank, arm, arm_config)
            predictions[f"{case_id}::{arm}"] = prediction
            case_results[arm] = record

        stages0 = [case_results[arm]["stages"][0] for arm in arms]
        agreement = max(s["true_data_residual"] for s in stages0) - min(s["true_data_residual"] for s in stages0)
        records[case_id]["identical_initial_stage_data_residual_max_difference"] = agreement
        print(case_id, "full-wave residual",
              " ".join(f"{arm}={case_results[arm]['true_data_residual']:.4f}" for arm in arms),
              "CPU " + "/".join(f"{case_results[arm]['timing']['cpu_s']:.2f}" for arm in arms), flush=True)

    out = ROOT / "results" / f"native_current{tag}"
    np.savez_compressed(out.with_suffix(".npz"), **predictions)
    result = {
        "cohort": CASE_COHORT,
        "archive_load_timing": load_time,
        "schedule_short_budget": list(ITERATIONS),
        "native_mf": list(NATIVE_MF), "small_mf": list(SMALL_MF),
        "swap_count": SWAP_COUNT,
        "cases": records,
    }
    out.with_suffix(".json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


def score(tag):
    from cohorts import load_truth
    from score_materials import material_metrics
    out = ROOT / "results" / f"native_current{tag}"
    result = json.loads(out.with_suffix(".json").read_text())
    with np.load(out.with_suffix(".npz")) as predictions:
        for case_id, row in result["cases"].items():
            truth = load_truth(CASE_COHORT, case_id)
            row["truth_peak"] = float(truth.max())
            for arm, record in row["arms"].items():
                prediction = predictions[f"{case_id}::{arm}"]
                record["materials"] = material_metrics(prediction, truth, (24, 24), (0.02, 20.0))
    out.with_suffix(".json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("run", "score", "run_additions", "score_additions"))
    parser.add_argument("--max-cases", type=int)
    parser.add_argument("--tag", default="")
    args = parser.parse_args()
    additions = args.action.endswith("_additions")
    tag = args.tag or ("_additions" if additions else "")
    if args.action.startswith("run"):
        run(args.max_cases, tag, ADDITION_ARMS if additions else ARMS)
    else:
        score(tag)
