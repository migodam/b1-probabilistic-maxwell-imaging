"""Reconstruct every independent NDT geometry with both material solvers."""

import argparse
import csv
from dataclasses import replace
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.ndimage import label, maximum_filter
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from src.designs import RECOVERY_ARMS, Design, all_designs, read_arm
from src.inversion import (
    MaterialPenalty, acquisition_models, energy, fft_tsom, fullwave, projected_gn,
)
from src.ndt import BOUNDS, grid, host, read_data


ROOT = Path(__file__).resolve().parent
SUPPORT_THRESHOLD = 0.20
PENALTY = MaterialPenalty(tv=0.0055, amplitude=0.0035)
TSOM_STAGES = ((3, 10), (6, 14), (9, 60))
TSOM_STATE_WEIGHT = 0.32
TSOM_SINGULAR_CUT = 0.13
GN_ITERATIONS = 24
GN_CALIBRATION_SCALES = (1, 2, 4, 8, 16, 32, 64)


def _centroid(x, xy, selected):
    mass = abs(x[selected])
    if mass.size:
        return mass @ xy[selected] / mass.sum()
    return None


def _paired_peaks(estimate, centers, xy):
    """Score whether both true holes yield nearby distinct material minima."""
    amp = -estimate.reshape(24, 24)
    maxima = (amp >= SUPPORT_THRESHOLD) & (amp == maximum_filter(amp, size=3))
    components, count = label(maxima, structure=np.ones((3, 3)))
    candidate = []
    for index in range(1, count + 1):
        members = components.ravel() == index
        candidate.append((amp.ravel()[members].max(), xy[members].mean(axis=0)))
    candidate.sort(key=lambda peak: peak[0], reverse=True)
    apart = float(np.linalg.norm(centers[0] - centers[1]))
    stats = {"true_separation_mm": apart, "two_peaks_resolved": False,
             "paired_peak_mean_error_mm": None}
    if len(candidate) < 2:
        return stats
    distances = cdist(centers, np.stack([peak[1] for peak in candidate[:2]]))
    row, col = linear_sum_assignment(distances)
    best = distances[row, col]
    stats["paired_peak_mean_error_mm"] = float(best.mean())
    stats["two_peaks_resolved"] = bool(best.max() <= min(1.5, apart / 2))
    return stats


def pair_centers(path=None):
    if path is None:
        path = ROOT / "results" / "acquisition.json"
    with path.open() as source:
        metadata = json.load(source)
    return {
        scene["id"]: np.array([[item["x_mm"], item["y_mm"]] for item in scene["primitives"]])
        for scene in metadata["case_definitions"]
        if (scene["id"].startswith("test_pair") or scene.get("spacing_mm") is not None)
        and len(scene["primitives"]) == 2
    }


def score(estimate, truth, xy, true_pair_centers=None):
    predicted = abs(estimate) >= SUPPORT_THRESHOLD
    actual = abs(truth) >= SUPPORT_THRESHOLD
    union = (actual | predicted).sum()
    correct_sign = (estimate * truth) > 0
    signal = actual & predicted & correct_sign
    center_true = _centroid(truth, xy, actual)
    center_pred = _centroid(estimate, xy, predicted)
    denom = energy(truth[actual])
    metrics = {
        "rmse_eps": float(np.sqrt(np.mean((estimate - truth) ** 2))),
        "support_iou": float(signal.sum() / union) if union and np.any(actual) else None,
        "centroid_error_mm": (float(np.linalg.norm(center_true - center_pred))
                              if center_true is not None and center_pred is not None else None),
        "contrast_retention": float(estimate[actual] @ truth[actual] / denom) if denom else None,
        "mean_absolute_contrast_error": (float(np.mean(abs(estimate[actual] - truth[actual])))
                                          if denom else None),
        "peak_abs_eps": float(np.max(abs(estimate))),
        "min_eps_delta": float(np.min(estimate)),
        "max_eps_delta": float(np.max(estimate)),
        "false_support_cells": int(np.sum(predicted & ~actual)),
        "false_support_fraction": float(np.mean(predicted & ~actual)),
        "detected": bool(predicted.any()),
    }
    if true_pair_centers is not None:
        metrics.update(_paired_peaks(estimate, true_pair_centers, xy))
    return metrics


def evaluate(models, observed, clean, estimate, noise_scale, material_factors=None):
    predicted = fullwave(models, estimate, material_factors=material_factors)
    err_rms = float(np.sqrt(np.mean(abs(predicted - observed) ** 2)))
    clean_norm = np.linalg.norm(clean)
    return predicted, {
        "measured_relative_fullwave_residual": float(np.linalg.norm(predicted - observed) /
                                                      np.linalg.norm(observed))
                                              if np.any(observed) else None,
        "clean_relative_fullwave_residual": (float(np.linalg.norm(predicted - clean) / clean_norm)
                                             if clean_norm else None),
        "fullwave_residual_rms": err_rms,
        "fullwave_residual_in_noise_sd": float(err_rms / noise_scale) if noise_scale else None,
    }


def _write_figures(data_pack, arms, methods, results, rows, output_suffix, *, case_titles=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_arms = tuple(a for a in ("f1_t8_r8", "f2_t4_r8", "f4_t2_r8", "f8_t1_r8",
                                   "f4_t8_r16", "oneside_f2_t4_r8", "high_t8_r8",
                                   "loc_nominal") if a in arms)
    all_labels = list(data_pack["case_ids"])
    if case_titles is None:
        cases_by_id = ("test_hole", "test_pair_close", "test_pair_wide", "test_slot",
                       "test_mixed", "test_empty")
        labels = ("Hole", "3.0 mm pair", "6.2 mm pair", "Thin slot", "Mixed contrast", "Intact")
        case_titles = zip(cases_by_id, labels)
    arm_labels = {"f1_t8_r8": "16 GHz", "f2_t4_r8": "2 / 30 GHz", "f4_t2_r8": "4 frequencies",
                  "f8_t1_r8": "8 frequencies", "f4_t8_r16": "4 f × 8 Tx × 16 Rx",
                  "oneside_f2_t4_r8": "2 / 30 GHz, one side", "high_t8_r8": "30 GHz",
                  "loc_nominal": "28 / 30 GHz, designed"}
    method_labels = {"fft_tsom": "FFT-TSOM", "projected_gn": "Projected GN"}
    scenario = {
        "_holdout_constant": "constant dielectric plate",
        "_holdout_debye": "known Debye plate",
        "_holdout_debye_wrong_host": "Debye data with a constant host model",
    }.get(output_suffix, "finite dielectric plate")
    if not figure_arms:
        return
    case_labels = [(case, title) for case, title in case_titles
                   if np.isfinite(results[arms.index(figure_arms[0]), all_labels.index(case), 0]).all()]
    if not case_labels:
        return
    row_dict = {(x["arm"], x["case"], x["method"]): x for x in rows}
    for method in methods:
        fig, ax = plt.subplots(len(case_labels), len(figure_arms) + 1,
                               figsize=(2.05 * (len(figure_arms) + 1), 2.05 * len(case_labels)),
                               squeeze=False, layout="constrained")
        method_idx = methods.index(method)
        for row, (case, label_name) in enumerate(case_labels):
            case_idx = all_labels.index(case)
            truth = data_pack["truth"][case_idx]
            ax[row, 0].set_ylabel(label_name, fontsize=11)
            images = [truth] + [results[arms.index(a), case_idx, method_idx]
                                for a in figure_arms]
            for col, image in enumerate(images):
                panel = ax[row, col]
                im = panel.imshow(image, extent=(-12, 12, -12, 12), origin="lower",
                                  cmap="RdBu_r", vmin=-1.5, vmax=1.5,
                                  interpolation="nearest")
                if col:
                    metrics = row_dict[(figure_arms[col - 1], case, method)]
                    iou = metrics["support_iou"]
                    detail = f"IoU {iou:.2f}" if iou is not None else (
                        f"{metrics['false_support_cells']} flagged cells")
                    if "two_peaks_resolved" in metrics:
                        detail += "  ·  " + ("two peaks" if metrics["two_peaks_resolved"] else "unresolved")
                    panel.text(0.02, 0.02, detail, ha="left", va="bottom", fontsize=8,
                               transform=panel.transAxes,
                               bbox={"fc": "white", "alpha": 0.88, "ec": "none", "pad": 2})
                panel.set_xticks([])
                panel.set_yticks([])
        for i, a in enumerate(("Truth",) + tuple(arm_labels[x] for x in figure_arms)):
            ax[0, i].set_title(a, fontsize=10)
        fig.suptitle(f"{method_labels[method]} · {scenario}", fontsize=14)
        fig.colorbar(im, ax=ax, shrink=.5, label="Relative permittivity deviation")
        fig.savefig(ROOT / "figures" / f"reconstruction_{method}{output_suffix}.png", dpi=175)
        plt.close(fig)


def _plot_dense_regularization(data_pack, arms, methods, results, rows):
    dense_arm = "f4_t8_r16"
    baseline_path = ROOT / "results" / "inversion_materials.npz"
    if dense_arm not in arms or not baseline_path.exists() or "projected_gn" not in methods:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with np.load(baseline_path, allow_pickle=False) as baseline:
        baseline_materials = baseline["estimates"][list(baseline["arms"]).index(dense_arm)]
        baseline_methods = list(baseline["methods"])
    calibrated_materials = results[arms.index(dense_arm), :, methods.index("projected_gn")]
    baseline_report = json.loads((ROOT / "results" / "inversion_metrics.json").read_text())
    method_metrics = {(row["case"], row["method"]): row for row in baseline_report["rows"]
                      if row["arm"] == dense_arm}
    calibration_metrics = {row["case"]: row for row in rows
                           if row["arm"] == dense_arm and row["method"] == "projected_gn"}
    desired = {"test_hole": "Hole", "test_pair_close": "3.0 mm pair",
               "test_slot": "Thin slot", "test_mixed": "Mixed contrast", "test_empty": "Intact"}
    scenes = [(i, c) for i, c in enumerate(data_pack["case_ids"])
              if c in desired and c in calibration_metrics]
    if not scenes:
        return
    fig, axes = plt.subplots(len(scenes), 4, figsize=(8.5, 10.0), layout="constrained", squeeze=False)
    first_method = baseline_methods.index("fft_tsom")
    second_method = baseline_methods.index("projected_gn")
    for line, (i, c) in enumerate(scenes):
        metrics = [None, method_metrics[(c, "fft_tsom")], method_metrics[(c, "projected_gn")],
                   calibration_metrics[c]]
        material = [data_pack["truth"][i], baseline_materials[i, first_method],
                    baseline_materials[i, second_method], calibrated_materials[i]]
        axes[line, 0].set_ylabel(desired[c], fontsize=10)
        for j, (panel, image, score_row) in enumerate(zip(axes[line], material, metrics)):
            im = panel.imshow(image, extent=(-12, 12, -12, 12), origin="lower", cmap="RdBu_r",
                              vmin=-1.5, vmax=1.5, interpolation="nearest")
            panel.set_xticks([])
            panel.set_yticks([])
            if score_row is not None:
                if c == "test_empty":
                    caption = f"{score_row['false_support_cells']} flagged cells"
                else:
                    e = score_row["rmse_eps"]
                    r = score_row["measured_relative_fullwave_residual"]
                    caption = f"RMSE {e:.3f}\nData error {100 * r:.1f}%"
                panel.text(.02, .02, caption, transform=panel.transAxes, fontsize=7.5,
                           bbox={"fc": "white", "alpha": .88, "ec": "none", "pad": 2})
    for panel, title in zip(axes[0], ("Truth", "FFT-TSOM", "Projected GN",
                                      "Projected GN, calibrated")):
        panel.set_title(title, fontsize=10)
    fig.suptitle("Four-frequency full-array reconstruction with cross-grid data", fontsize=12)
    fig.colorbar(im, ax=axes, shrink=.56, label="Relative permittivity deviation")
    fig.savefig(ROOT / "figures" / "reconstruction_dense_regularization.png", dpi=180)
    fig.savefig(ROOT / "figures" / "reconstruction_dense_regularization.pdf")
    plt.close(fig)


def _penalty_at_scale(scale):
    return replace(PENALTY, tv=PENALTY.tv * scale, amplitude=PENALTY.amplitude * scale)


def _calibrate_gn(models, pack, design, same_grid, noise_fraction, reference, material_factors):
    started = perf_counter()
    development_cases = np.flatnonzero(pack["splits"] == "development")
    results = {}
    for scale in GN_CALIBRATION_SCALES:
        rmse = {}
        for index in development_cases:
            y, _ = read_arm(pack, index, design, same_grid=same_grid,
                            sigma_fraction=noise_fraction)
            estimate, _, _ = projected_gn(
                models, y, reference, _penalty_at_scale(scale), iterations=GN_ITERATIONS,
                material_factors=material_factors,
            )
            rmse[str(pack["case_ids"][index])] = float(
                np.sqrt(np.mean((estimate.reshape(24, 24) - pack["truth"][index]) ** 2))
            )
        results[scale] = rmse
    chosen = min(results, key=lambda s: np.mean(list(results[s].values())))
    return {"scale": chosen, "development_rmse_eps": results,
            "calibration_seconds": perf_counter() - started}


def _write_report(suffix, report):
    path = ROOT / "results" / f"inversion_metrics{suffix}.json"
    with path.open("w") as output:
        json.dump(report, output, indent=2)
    rows = report["rows"]
    columns = list(dict.fromkeys(key for item in rows for key in item))
    with path.with_suffix(".csv").open("w") as output:
        writer = csv.DictWriter(output, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def _score_saved_predictions(pack, truth_metadata, suffix):
    xy_mm = grid(24) * 1000
    pair_centers_by_id = pair_centers(truth_metadata)
    result_path = ROOT / "results" / f"inversion_materials{suffix}.npz"
    with np.load(result_path, allow_pickle=False) as raw:
        arms, methods, cases = map(list, (raw["arms"], raw["methods"], raw["case_ids"]))
        estimates = raw["estimates"]
    if cases != list(pack["case_ids"]):
        raise ValueError("Saved predictions and truth have different case ordering")
    path = ROOT / "results" / f"inversion_metrics{suffix}.json"
    with path.open() as source:
        report = json.load(source)
    for row in report["rows"]:
        scene = row["case"]
        idx = cases.index(scene)
        individual = estimates[arms.index(row["arm"]), idx, methods.index(row["method"])].ravel()
        row.update(score(individual, pack["truth"][idx].ravel(), xy_mm, pair_centers_by_id.get(scene)))
    report["scored_after_predictions"] = True
    _write_report(suffix, report)
    with truth_metadata.open() as source:
        scenes = json.load(source)["case_definitions"]
    shape_titles = {
        "single_air_hole": "Round hole",
        "elliptical_air_hole": "Elliptical hole",
        "thin_slot": "Thin slot",
        "mixed_air_hole_and_dielectric_load": "Mixed contrast",
        "three_air_defects": "Three holes",
        "bent_slot_unseen_geometry": "Bent slot",
        "empty": "Intact",
    }
    titles = []
    for scene in scenes:
        spacing = scene.get("spacing_mm")
        name = shape_titles.get(scene["shape"], scene["shape"].replace("_", " "))
        titles.append((scene["id"], f"{spacing:.2f} mm pair" if spacing is not None else name))
    _write_figures(pack, tuple(arms), tuple(methods), estimates, report["rows"], suffix,
                   case_titles=titles)


def main():
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("--arms", nargs="+")
    args.add_argument("--cases", nargs="+")
    args.add_argument("--design-json", type=Path)
    args.add_argument("--methods", nargs="+", choices=("fft_tsom", "projected_gn"))
    args.add_argument("--calibrate-gn", action="store_true")
    args.add_argument("--gn-penalty-json", nargs="+", type=Path)
    args.add_argument("--data", type=Path)
    args.add_argument("--geometry", type=Path)
    args.add_argument("--host-model", choices=("constant", "debye"), default="constant")
    args.add_argument("--score-blind", action="store_true")
    args.add_argument("--same-grid", action="store_true")
    args.add_argument("--sigma-fraction", type=float, default=.02)
    args.add_argument("--suffix", default="")
    args.add_argument("--development-only", action="store_true")
    opts = args.parse_args()
    if opts.development_only and opts.cases:
        args.error("--development-only and --cases cannot be combined")
    if opts.score_blind and (not opts.data or not opts.suffix):
        args.error("--score-blind requires --data and --suffix")
    suffix = opts.suffix or ("_same_grid" if opts.same_grid else "")
    pack = read_data(opts.data)
    geometry = opts.geometry or (opts.data.with_suffix(".json") if opts.data else None)
    if opts.score_blind:
        _score_saved_predictions(pack, geometry, suffix)
        return
    observed_only = "truth" not in pack
    pair_centers_by_id = {} if observed_only else pair_centers(geometry)
    if observed_only and opts.calibrate_gn:
        args.error("Calibrate on development data before reading holdout observations")
    fixed_scales = {}
    for path in opts.gn_penalty_json or ():
        with path.open() as source:
            for arm, item in json.load(source)["gn_development_calibration"].items():
                fixed_scales[arm] = item["scale"]
    ids = list(pack["case_ids"])
    design_options = all_designs()
    extra = {}
    if opts.design_json:
        with opts.design_json.open() as design_file:
            for specification in json.load(design_file)["designs"]:
                design = Design(specification["name"], tuple(specification["freqs"]),
                                tuple(specification["tx"]), tuple(specification["rx"]))
                extra[design.name] = design
        design_options.update(extra)
    arms = tuple(opts.arms or extra or RECOVERY_ARMS)
    ref = float(pack["reference_amplitude"])
    methods = tuple(opts.methods or ("fft_tsom", "projected_gn"))
    fields = np.full((len(arms), len(ids), len(methods), 24, 24), np.nan)
    rows, all_traces, calibrations = [], [], {}
    debye = opts.host_model == "debye"
    if debye:
        from src.ndt import dielectric_factor, host_debye
    for j, arm in enumerate(arms):
        d = design_options[arm]
        host.cache_clear()
        if debye:
            host_debye.cache_clear()
        material_factors = [dielectric_factor(f) for f in d.freqs] if debye else None
        acquisition_models.cache_clear()
        started = perf_counter()
        models = acquisition_models(24, d, debye=debye)
        setup_time = perf_counter() - started
        xy = models[0].xy * 1000
        gn_scale = fixed_scales[arm] if fixed_scales and "projected_gn" in methods else 1
        if opts.calibrate_gn and "projected_gn" in methods:
            comparison = _calibrate_gn(models, pack, d, opts.same_grid, opts.sigma_fraction, ref,
                                       material_factors)
            calibrations[arm] = comparison
            gn_scale = comparison["scale"]
            print(f"{arm} projected_gn calibration chose scale {gn_scale}", flush=True)
        gn_penalty = _penalty_at_scale(gn_scale)
        for i, case in enumerate(ids):
            if opts.cases and case not in opts.cases:
                continue
            if opts.development_only and pack["splits"][i] != "development":
                continue
            y, sigma = read_arm(pack, i, d, same_grid=opts.same_grid,
                                sigma_fraction=opts.sigma_fraction)
            clean, _ = read_arm(pack, i, d, same_grid=opts.same_grid, sigma_fraction=0.0)
            true = None if observed_only else pack["truth"][i].ravel()
            for k, method in enumerate(methods):
                t0 = perf_counter()
                if method == "fft_tsom":
                    result = fft_tsom(models, y, ref, sigma, PENALTY, stages=TSOM_STAGES,
                                      relative_cut=TSOM_SINGULAR_CUT, state_tradeoff=TSOM_STATE_WEIGHT,
                                      material_factors=material_factors)
                else:
                    result = projected_gn(models, y, ref, gn_penalty, iterations=GN_ITERATIONS,
                                          material_factors=material_factors)
                solver_time = perf_counter() - t0
                x, trace, snapshots = result[:3]
                fields[j, i, k] = x.reshape(24, 24)
                _, residuals = evaluate(models, y, clean, x, sigma, material_factors)
                row = {
                    "arm": arm, "case": case, "split": str(pack["splits"][i]), "method": method,
                    "frequency_count": len(d.freqs), "transmitters": len(d.tx),
                    "receivers": len(d.rx), "complex_samples": d.samples, "shots": d.shots,
                    "shot_energy": d.shot_energy, "noise_sd": sigma,
                    "shared_host_seconds": setup_time, "solver_seconds": solver_time,
                    "end_to_end_seconds": setup_time + solver_time,
                    "accepted_iterations": len(trace),
                    **residuals,
                }
                if true is not None:
                    row.update(score(x, true, xy, pair_centers_by_id.get(case)))
                if method == "fft_tsom":
                    row["deterministic_current_ranks"] = result[3]
                    if trace:
                        row["fft_tsom_source_data_residual_rms"] = ref * np.sqrt(2 * trace[-1][3])
                        row["fft_tsom_weighted_state_loss"] = trace[-1][4]
                else:
                    row["gn_penalty_scale"] = gn_scale
                if len(snapshots):
                    for checkpoint in sorted({0, len(trace)//3, 2*len(trace)//3, len(trace)-1}):
                        estimate = snapshots[checkpoint]
                        _, checkpoint_fit = evaluate(models, y, clean, estimate, sigma, material_factors)
                        all_traces.append({
                            "arm": arm, "case": case, "method": method,
                            "accepted_iteration": checkpoint + 1,
                            "solver_seconds": trace[checkpoint][1],
                            "end_to_end_seconds": setup_time + trace[checkpoint][1],
                            "objective": trace[checkpoint][2],
                            "state_loss": trace[checkpoint][4],
                            "rmse_eps": (float(np.sqrt(np.mean((estimate - true) ** 2)))
                                         if true is not None else None),
                            **checkpoint_fit,
                        })
                rows.append(row)
                residual = row["measured_relative_fullwave_residual"]
                fit_label = f"{residual:.3f}" if residual is not None else "n/a"
                rmse_label = f" RMSE {row['rmse_eps']:.3f}" if true is not None else ""
                print(f"{arm} {case} {method} {solver_time:.2f}s{rmse_label} res {fit_label}",
                      flush=True)

    results_dir = ROOT / "results"
    truth_fields = {} if observed_only else {"truth": pack["truth"]}
    np.savez_compressed(results_dir / f"inversion_materials{suffix}.npz",
                        estimates=fields, case_ids=pack["case_ids"], splits=pack["splits"],
                        arms=arms, methods=methods, **truth_fields)
    _write_report(suffix, {
        "same_grid": opts.same_grid, "sigma_fraction": opts.sigma_fraction,
        "host_model": opts.host_model, "source_data": opts.data.name if opts.data else None,
        "support_threshold_eps": SUPPORT_THRESHOLD, "material_bounds": BOUNDS,
        "penalty": vars(PENALTY), "gn_iteration_cap": GN_ITERATIONS,
        "tsom_stages": TSOM_STAGES, "tsom_state_tradeoff": TSOM_STATE_WEIGHT,
        "tsom_relative_singular_cut": TSOM_SINGULAR_CUT,
        "tsom_projected_signal_noise_cut": 2.75, "gn_development_calibration": calibrations,
        "rows": rows, "trace": all_traces,
    })
    if not observed_only:
        _write_figures(pack, arms, methods, fields, rows, suffix)
        if opts.calibrate_gn:
            _plot_dense_regularization(pack, arms, methods, fields, rows)


if __name__ == "__main__":
    main()
