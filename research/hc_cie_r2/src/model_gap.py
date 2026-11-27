"""Reconstruct native SRAP model error and paired clean/noisy ablations."""
import argparse
import csv
import json
import os
from pathlib import Path
import sys
import time

import spectra as sp  # initialize single-thread native BLAS before NumPy is imported
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import minimize


HERE = Path(__file__).resolve().parents[1]
CONFIG = json.loads((HERE / "configs/model_gap.json").read_text())


def write_json(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf8")


def write_csv(path, rows):
    keys = list(dict.fromkeys(key for row in rows for key in row if not isinstance(row[key], dict)))
    with path.open("w", encoding="utf8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def original_sources(workspace):
    srap = workspace / "research/srap_cie_array_study_20261121"
    os.environ.setdefault("LA_WORKSPACE", str(workspace))
    sys.path.insert(0, str(srap))
    import data2d
    return srap, data2d


def plan(data):
    permutation = np.random.default_rng(data.FAMILY_SEED).permutation(np.repeat(np.arange(6), 4))
    objects = sorted((int(np.flatnonzero(permutation == family)[0] + 1), family) for family in range(6))
    return [(f"s2d_{index:03d}", index, family) for index, family in objects]


def select_ports(matrix, operators, name):
    return matrix[np.ix_(operators.rx_ids[name], operators.tx_ids[name])]


def prepare(workspace, output):
    srap, data = original_sources(workspace)
    observations, offline = output / "public", output / "offline"
    observations.mkdir(exist_ok=True)
    offline.mkdir(exist_ok=True)
    start, cpu = time.perf_counter(), time.process_time()
    operators = data.SyntheticOperators(96, data.CONFIRM_ARRAYS)
    results, forward_costs = [], []
    fine_step = operators.fine // data.SIDE
    assert operators.fine == fine_step * data.SIDE
    coarse_green = None
    for object_id, index, family in plan(data):
        clock, cpu0 = time.perf_counter(), time.process_time()
        material_fine = data.random_material(operators.fine_centers, index, family)
        material = data.area_average(material_fine, operators.fine, data.SIDE)
        block_lift = np.repeat(np.repeat(material.reshape(data.SIDE, data.SIDE), fine_step, axis=0),
                               fine_step, axis=1).ravel()
        fine_clean, actual_cost = data.history.fine_solve(material_fine, operators.fine_e, operators.fine_s,
                                                         operators.green.apply, operators.green.diagonal)
        fine_lift, lift_cost = data.history.fine_solve(block_lift, operators.fine_e, operators.fine_s,
                                                       operators.green.apply, operators.green.diagonal)
        forward_costs.append(dict(object_id=object_id, full_fine=actual_cost, lifted_fine=lift_cost,
                                  wall_s=time.perf_counter() - clock, cpu_s=time.process_time() - cpu0))
        for name in CONFIG["layouts"]:
            array_index = data.CONFIRM_ARRAYS.index(name)
            key = data.synthetic_name("confirmation", object_id, name)
            source = srap / f"data/public/{key}.npz"
            with np.load(source, allow_pickle=False) as saved:
                delta, observed, b, s = float(saved["delta"]), saved["Y"], saved["E"], saved["S"]
                canonical = {k: saved[k] for k in ("centers", "radius", "k", "E", "S")}
                meta = json.loads(str(saved["case_json"]))
            # The historical truth check is isolated in data preparation. Its archive is never used by solve().
            with np.load(srap / f"data/truth/{key}.npz", allow_pickle=False) as saved:
                assert np.array_equal(material, saved["m"])
                assert str(saved["family"]) == data.FAMILY_NAMES[family]
            if coarse_green is None:
                coarse_green = data.history.PaperModel(canonical, meta, data.history.Clock(1800)).G
                np.savez_compressed(observations / "coarse_green.npz", G=coarse_green)
            assert np.array_equal(b, operators.coarse_e[:, operators.tx_ids[name]])
            assert np.array_equal(s, operators.coarse_s[operators.rx_ids[name]])
            rng_seed = np.random.SeedSequence([data.NOISE_SEED, index, array_index, 0])
            rng = np.random.default_rng(rng_seed)
            direction = rng.standard_normal(observed.shape) + 1j * rng.standard_normal(observed.shape)
            noise = direction * (delta / np.linalg.norm(direction))
            cross_clean = observed - noise
            cross_check = select_ports(fine_clean, operators, name)
            lifted = select_ports(fine_lift, operators, name)
            independent_noise = data.history.complex_noise(cross_check, data.NOISE_FRACTION, rng_seed)
            replay_error = np.linalg.norm(cross_clean - cross_check) / np.linalg.norm(cross_clean)
            noise_error = np.linalg.norm(noise - independent_noise) / delta
            # Use the exact saved noisy field as primary and validate it against an independent fine solve.
            if replay_error > 5e-7 or noise_error > 5e-7:
                raise RuntimeError(f"{key}: historical noise/port reconstruction disagreement ({replay_error}, {noise_error})")
            matched = s @ sp._current(coarse_green, b, material)
            geometry_gap = cross_clean - lifted
            discretization_gap = lifted - matched
            total_gap = geometry_gap + discretization_gap
            norms = {"noise": np.linalg.norm(noise), "microgeometry": np.linalg.norm(geometry_gap),
                     "discretization": np.linalg.norm(discretization_gap), "total_gap": np.linalg.norm(total_gap),
                     "total_gap_plus_noise": np.linalg.norm(observed - matched)}
            full = np.vdot(total_gap, noise).real
            results.append(dict(
                observation=key, object_id=object_id, family=data.FAMILY_NAMES[family], array=name,
                source=str(source.relative_to(workspace)), native_fine_side=meta["fine"], coarse_side=32,
                delta=float(delta), native_noise_replay_relative=float(noise_error),
                fine_clean_replay_relative=float(replay_error),
                historical_noisy_replay_relative=float(np.linalg.norm(observed - cross_check - independent_noise) /
                                                       np.linalg.norm(observed)),
                reference_contrast_peak=float(material.max()),
                **{f"{key}_over_delta": float(value / delta) for key, value in norms.items()},
                gap_noise_cosine=float(full / (norms["total_gap"] * delta)),
                gap_noise_interference_over_delta2=float(2 * full / delta**2),
                structural_discrete_cross_term_over_delta2=float(2 * np.vdot(geometry_gap, discretization_gap).real / delta**2),
            ))
            np.savez_compressed(observations / f"{key}.npz", B=b, S=s, delta=delta, matched_clean=matched,
                                matched_noisy=matched + noise, cross_clean=cross_clean, cross_noisy=observed)
            np.savez_compressed(offline / f"{key}.npz", material=material, noise=noise, clean_resolved=cross_check,
                                lifted=lifted, microgeometry=geometry_gap, discretization=discretization_gap,
                                total_gap=total_gap)
            print("prepared", key, f"gap/noise={norms['total_gap'] / delta:.3f}", flush=True)
    write_csv(output / "gaps.csv", results)
    write_json(output / "gaps.json", dict(rows=results, forward_costs=forward_costs, setup=operators.cost,
                                        wall_s=time.perf_counter() - start, cpu_s=time.process_time() - cpu,
                                        coarse_truth_forward_rhs=10 * len(results),
                                        original_material_and_noise_seeds=dict(material=data.MATERIAL_SEED,
                                                                                family=data.FAMILY_SEED,
                                                                                noise=data.NOISE_SEED)))


def sample_score(estimate, truth):
    error = estimate - truth
    threshold = 0.1 * truth.max()
    interior = truth > threshold
    return dict(material_rel_l2=float(np.linalg.norm(error) / np.linalg.norm(truth)),
                interior_rel_l2=float(np.linalg.norm(error[interior]) / np.linalg.norm(truth[interior])),
                background_error_over_truth_l2=float(np.linalg.norm(error[~interior]) / np.linalg.norm(truth)),
                truth_interior_cells=int(interior.sum()), estimate_peak=float(estimate.max()))


def native_solve(g, s, b, y, common_yn):
    objective_calls = 0
    t0, c0 = time.perf_counter(), time.process_time()
    width = len(b)
    previous = None
    history = []

    def objective(m):
        nonlocal objective_calls, previous
        objective_calls += 1
        system = -m[:, None] * g
        system.flat[::width + 1] += 1
        lu = lu_factor(system, overwrite_a=True, check_finite=False)
        j = lu_solve(lu, m[:, None] * b, check_finite=False)
        predicted = s @ j
        error = predicted - y
        adjoint = lu_solve(lu, s.conj().T @ error, trans=2, check_finite=False)
        grad = np.sum((adjoint.conj() * (b + g @ j)).real, axis=1) / common_yn**2
        last = float(np.linalg.norm(error))
        previous = dict(residual=last, ev=objective_calls)
        return 0.5 * (last / common_yn)**2, grad

    def accept(_m):
        history.append(previous.copy())

    fit = minimize(objective, np.full(width, CONFIG["initial_material"]), jac=True, method="L-BFGS-B",
                   bounds=[tuple(CONFIG["box"])] * width, callback=accept,
                   options=dict(maxiter=CONFIG["iterations"], maxfun=CONFIG["function_cap"],
                                ftol=1e-13, gtol=1e-8))
    prediction = s @ sp._current(g, b, fit.x)
    return fit.x, prediction, dict(iterations=int(fit.nit), evaluations=objective_calls, history=history,
                                   stop=str(fit.message), normalized_objective=float(fit.fun),
                                   initial_residual=float(np.linalg.norm(y)),
                                   residual=float(np.linalg.norm(prediction - y)),
                                   forward_rhs=(objective_calls + 1) * b.shape[1],
                                   adjoint_rhs=objective_calls * b.shape[1],
                                   wall_s=time.perf_counter() - t0, cpu_s=time.process_time() - c0)


def solve(output):
    source = output / "public"
    prediction = output / "predictions"
    prediction.mkdir(exist_ok=True)
    with np.load(source / "coarse_green.npz", allow_pickle=False) as saved:
        g = saved["G"]
    collected = []
    for path in sorted(source.glob("s2d_*.npz")):
        with np.load(path, allow_pickle=False) as inputs:
            g_data = {name: inputs[name] for name in ("B", "S", *CONFIG["conditions"])}
            common_y_norm = np.linalg.norm(inputs["cross_noisy"])
            delta = float(inputs["delta"])
        for condition in CONFIG["conditions"]:
            material, pred, stats = native_solve(g, g_data["S"], g_data["B"], g_data[condition], common_y_norm)
            np.savez_compressed(prediction / f"{path.stem}__{condition}.npz", material=material, Y_pred=pred)
            stats.update(observation=path.stem, condition=condition,
                         residual_over_delta=float(stats["residual"] / delta))
            collected.append(stats)
            write_json(output / "solve_trace.json", collected)
            print("solve", path.stem, condition, f"res/δ={stats['residual_over_delta']:.3f}", flush=True)


def score(workspace, output):
    srap = workspace / "research/srap_cie_array_study_20261121"
    gaps = json.loads((output / "gaps.json").read_text())["rows"]
    trace = json.loads((output / "solve_trace.json").read_text())
    assert len(trace) == len(gaps) * len(CONFIG["conditions"])
    assert all((output / "predictions" / f"{r['observation']}__{r['condition']}.npz").is_file() for r in trace)
    results = []
    with np.load(output / "public/coarse_green.npz", allow_pickle=False) as saved:
        g = saved["G"]
    for gap in gaps:
        case = gap["observation"]
        with np.load(output / f"offline/{case}.npz", allow_pickle=False) as saved:
            truth = saved["material"]
        with np.load(output / f"public/{case}.npz", allow_pickle=False) as saved:
            b, s, y, delta = saved["B"], saved["S"], saved["cross_noisy"], float(saved["delta"])
            native_clean, matched_clean = saved["cross_clean"], saved["matched_clean"]
        for condition in CONFIG["conditions"]:
            with np.load(output / f"predictions/{case}__{condition}.npz", allow_pickle=False) as saved:
                material, pred = saved["material"], saved["Y_pred"]
            stats = next(r for r in trace if r["observation"] == case and r["condition"] == condition)
            results.append(dict(observation=case, family=gap["family"], array=gap["array"],
                                method="physical_LBFGSB_12", condition=condition,
                                residual_over_delta=stats["residual_over_delta"],
                                archival_port_residual_over_delta=float(np.linalg.norm(pred - y) / delta),
                                fine_clean_residual_over_delta=float(np.linalg.norm(pred - native_clean) / delta),
                                matched_clean_residual_over_delta=float(np.linalg.norm(pred - matched_clean) / delta),
                                **sample_score(material, truth)))
        source = srap / f"results/predictions/{case}__dual.npz"
        started, cpu = time.perf_counter(), time.process_time()
        with np.load(source, allow_pickle=False) as saved:
            material = saved["m_30_light"]
        pred = s @ sp._current(g, b, material)
        results.append(dict(observation=case, family=gap["family"], array=gap["array"],
                            method="historical_Dual_30s_light", condition="cross_noisy",
                            relative_source=str(source.relative_to(workspace)),
                            residual_over_delta=float(np.linalg.norm(pred - y) / delta),
                            archival_port_residual_over_delta=float(np.linalg.norm(pred - y) / delta),
                            fine_clean_residual_over_delta=float(np.linalg.norm(pred - native_clean) / delta),
                            matched_clean_residual_over_delta=float(np.linalg.norm(pred - matched_clean) / delta),
                            audit_wall_s=time.perf_counter() - started, audit_cpu_s=time.process_time() - cpu,
                            audit_rhs=b.shape[1], **sample_score(material, truth)))
        np.savez_compressed(output / f"predictions/{case}__historical_dual_30s_light.npz", material=material, Y_pred=pred)
    write_csv(output / "scores.csv", results)
    write_json(output / "scores.json", results)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "solve", "score", "all"))
    parser.add_argument("--workspace", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, default=HERE / "results/model_gap")
    args = parser.parse_args()
    root, output = args.workspace.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.stage in ("prepare", "all"):
        prepare(root, output)
    if args.stage in ("solve", "all"):
        solve(output)
    if args.stage in ("score", "all"):
        score(root, output)


if __name__ == "__main__":
    main()
