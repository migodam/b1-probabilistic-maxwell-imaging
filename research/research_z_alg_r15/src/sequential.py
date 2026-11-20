"""Two-stage frequency selection using only a noisy 16 GHz pilot."""

import argparse
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from threadpoolctl import threadpool_limits

from designs import Design, spaced
from ndt import Host, dielectric_factor, grid, host, host_debye
from port_online import PortFields
from sparse_features import estimate, prepare, real_stack, remove_shared_gain

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from run_inversion import pair_centers, score
from src.inversion import fullwave
from src.roi_refine import region


OUT = ROOT / "results"
PACKS = ("holdout", "holdout_debye")
PILOT = 16.0
MENU = (6., 12., 20., 24., 30.)
FALLBACK = 20.0
TX, RX = spaced(8, 4), spaced(16, 8)
SHOT_ENERGY = 1 / 8
PRIOR_SD = .30
RANDOM_SEED = 15112030
FIXED_NAMES = tuple(f"fixed_{f:g}" for f in MENU)
ARMS = ("pilot_only", "active", "uniform_random", *FIXED_NAMES)


def aperture():
    return Design("pilot_and_followup_ports", (PILOT,), TX, RX)


def selected(data, i, f, sd):
    """Reveal one frequency on the programmed physical ports."""
    selection = np.ix_(RX, TX)
    return (data["clean"][i, f][selection] + sd * data["noise"][i, f][selection]).ravel()


def fit(prepared, y, blank, sd):
    target = real_stack(remove_shared_gain(y[:, None], blank)) * (np.sqrt(2) / sd)
    prediction, evidence, threshold = estimate(*prepared, target)
    return prediction[0].reshape(24, 24), float(evidence[0]), float(threshold)


def local_tangent(model, f, ids, q, debye):
    restricted = Host(model.Eh[np.ix_(ids, TX)], model.Sh[np.ix_(RX, ids)],
                      model.Gh[np.ix_(ids, ids)], model.xy[ids])
    multiplier = dielectric_factor(f) if debye else 1.
    _, j = fullwave([restricted], q, jacobian=True, material_factors=[multiplier])
    return j


def posterior_trace(j, blank, sd):
    j = real_stack(remove_shared_gain(j, blank)) * (np.sqrt(2) / sd)
    precision = j.T @ j
    precision.flat[:: len(precision) + 1] += 1 / PRIOR_SD**2
    c = cho_factor(precision, check_finite=False)
    covariance = cho_solve(c, np.eye(len(precision)), check_finite=False)
    return float(np.trace(covariance))


def choose_second(pilot_image, hosts, blanks, sd, debye):
    """The target readout for the second frequency is outside this interface."""
    ids, _ = region(pilot_image)
    if not len(ids):
        return FALLBACK, len(ids), {}, None
    q = pilot_image.ravel()[ids]
    j0 = local_tangent(hosts[PILOT], PILOT, ids, q, debye)
    trace0 = posterior_trace(j0, blanks[PILOT], sd)
    gains = {}
    for f in MENU:
        jf = local_tangent(hosts[f], f, ids, q, debye)
        stacked = np.vstack((j0, jf))
        blank = np.concatenate((blanks[PILOT], blanks[f]))
        gains[f] = (trace0 - posterior_trace(stacked, blank, sd)) / trace0
    return max(gains, key=gains.get), len(ids), gains, trace0


def run_pack(name):
    start, cpu = perf_counter(), process_time()
    with np.load(OUT / f"{name}_observations.npz", allow_pickle=False) as pack:
        fields = ("case_ids", "frequencies_ghz", "clean", "noise", "sigma_2pct", "background")
        data = {k: pack[k] for k in fields}
    debye = name == "holdout_debye"
    sd = float(data["sigma_2pct"]) / np.sqrt(SHOT_ENERGY)
    fi = {float(f): i for i, f in enumerate(data["frequencies_ghz"])}
    designs = (PILOT, *MENU)
    d = aperture()
    blanks = {f: data["background"][fi[f]][np.ix_(RX, TX)].ravel() for f in designs}
    host_function = host_debye if debye else host
    fine_tangent, hosts, prep, preparation = {}, {}, {}, {}
    for f in designs:
        tic = perf_counter()
        fine_tangent[f] = PortFields(48, f, debye).tangent(d)
        after_fine = perf_counter()
        hosts[f] = host_function(24, f)
        preparation[f] = dict(fine_tangent_seconds=after_fine - tic,
                              coarse_host_seconds=perf_counter() - after_fine)
    tic = perf_counter()
    pilot_fit = prepare(fine_tangent[PILOT], blanks[PILOT], sd, True)
    pilot_design_prepare = perf_counter() - tic
    for f in MENU:
        tic = perf_counter()
        blank = np.concatenate((blanks[PILOT], blanks[f]))
        prep[f] = (prepare(np.vstack((fine_tangent[PILOT], fine_tangent[f])), blank, sd, True), blank)
        preparation[f]["fixed_design_prepare_seconds"] = perf_counter() - tic
    random = np.random.default_rng(RANDOM_SEED).choice(MENU, size=len(data["case_ids"]))
    images = {name: [] for name in ARMS}
    rows = []
    for case_index, label in enumerate(data["case_ids"]):
        pilot = selected(data, case_index, fi[PILOT], sd)
        pilot_start = perf_counter()
        pilot_image, evidence, critical = fit(pilot_fit, pilot, blanks[PILOT], sd)
        after_pilot = perf_counter()
        f_active, nroi, gains, trace0 = choose_second(pilot_image, hosts, blanks, sd, debye)
        after_choice = perf_counter()
        images["pilot_only"].append(pilot_image)

        def after_second(f):
            y = np.concatenate((pilot, selected(data, case_index, fi[f], sd)))
            prepared, blank = prep[f]
            return fit(prepared, y, blank, sd)[0]

        # The adaptive image is formed before opening any other candidate readout.
        adaptive_image = after_second(f_active)
        after_active = perf_counter()
        images["active"].append(adaptive_image)
        fixed_fits, fixed_seconds = {}, {}
        for f, arm in zip(MENU, FIXED_NAMES):
            tic = perf_counter()
            fixed_fits[f] = after_second(f)
            fixed_seconds[f] = perf_counter() - tic
            images[arm].append(fixed_fits[f])
        if not np.array_equal(adaptive_image, fixed_fits[f_active]):
            raise RuntimeError("Adaptive and fixed solvers differ on identical readouts")
        images["uniform_random"].append(fixed_fits[random[case_index]])
        rows.append(dict(case=str(label), chosen_frequency_ghz=f_active, random_frequency_ghz=float(random[case_index]),
                         ROI_pixels=nroi, pilot_evidence=evidence, pilot_threshold=critical,
                         pilot_covariance_trace=trace0, relative_covariance_reduction=gains,
                         pilot_reconstruction_seconds=after_pilot - pilot_start,
                         selection_seconds=after_choice - after_pilot, active_final_reconstruction_seconds=after_active - after_choice,
                         fixed_final_reconstruction_seconds=fixed_seconds))
        print(f"{name} {label}: ROI {nroi:3d}, second {f_active:g} GHz", flush=True)
    timing = dict(host_and_operator_seconds=preparation, pilot_operator_prepare_seconds=pilot_design_prepare,
                  elapsed_seconds=perf_counter() - start, cpu_seconds=process_time() - cpu)
    return data["case_ids"], np.array([images[arm] for arm in ARMS]), rows, timing


def predict():
    fitted = [run_pack(name) for name in PACKS]
    config = dict(pilot_frequency_ghz=PILOT, candidates_ghz=MENU, blank_fallback_frequency_ghz=FALLBACK,
                  receivers=list(RX), transmitters=list(TX), shots=8, complex_readouts=64,
                  unit_total_energy=1., prior_pixel_sd=PRIOR_SD, uniform_random_seed=RANDOM_SEED,
                  selector="Local A-optimal conditional trace reduction for exact 24-cell full-wave Jacobian at the 16-GHz pilot, profiled shared complex blank-field gain.",
                  material_solver="Shared full-grid signed sparse solver, 48-grid intact-plate tangent area-aggregated to 24-grid material variables.")
    result = dict(protocol=config, data_sets={name: dict(cases=a[2], timing=a[3]) for name, a in zip(PACKS, fitted)})
    np.savez_compressed(OUT / "sequential_predictions.npz", packs=PACKS, arms=ARMS,
                        case_ids=np.stack([item[0] for item in fitted]), images=np.stack([item[1] for item in fitted]))
    (OUT / "sequential_choices.json").write_text(json.dumps(result, indent=2) + "\n")


def score_saved():
    rows, summaries = [], {}
    xy = 1000 * grid(24)
    with np.load(OUT / "sequential_predictions.npz", allow_pickle=False) as z:
        packs, arms, predictions, labels = (z[key] for key in ("packs", "arms", "images", "case_ids"))
    choices = json.loads((OUT / "sequential_choices.json").read_text())
    for name, all_images, labels0 in zip(packs, predictions, labels):
        with np.load(OUT / f"{name}.npz", allow_pickle=False) as scoring:
            truth = scoring["truth"]
            if not np.array_equal(labels0, scoring["case_ids"]):
                raise ValueError("Prediction/truth case ordering differs")
        centers = pair_centers(OUT / f"{name}.json")
        arm_summary = {}
        for arm, stack in zip(arms, all_images):
            records = []
            for label, image, material in zip(labels0, stack, truth):
                metrics = score(image.ravel(), material.ravel(), xy, centers.get(label))
                records.append(dict(pack=str(name), arm=str(arm), case=str(label), empty=not np.any(material), **metrics))
            nonempty = [row for row in records if not row["empty"]]
            pairs = sorted((row for row in records if "true_separation_mm" in row),
                           key=lambda r: r["true_separation_mm"])
            arm_summary[str(arm)] = dict(nonempty_mean_rmse=float(np.mean([r["rmse_eps"] for r in nonempty])),
                                         nonempty_median_rmse=float(np.median([r["rmse_eps"] for r in nonempty])),
                                         near_pair=pairs[0], wide_pair=pairs[-1],
                                         empty_detected=next(r["detected"] for r in records if r["empty"]))
            rows.extend(records)
        design_records = choices["data_sets"][str(name)]["cases"]
        frequencies = [r["chosen_frequency_ghz"] for r in design_records]
        arms_fixed = [arm_summary[key] for key in FIXED_NAMES]
        summaries[str(name)] = dict(arms=arm_summary,
                                   expected_uniform_random_mean_rmse=float(np.mean([r["nonempty_mean_rmse"] for r in arms_fixed])),
                                   expected_uniform_random_near_resolution=float(np.mean([r["near_pair"]["two_peaks_resolved"] for r in arms_fixed])),
                                   choice_counts={str(f): frequencies.count(f) for f in MENU},
                                   median_extra_pilot_selection_seconds=float(np.median([
                                       r["pilot_reconstruction_seconds"] + r["selection_seconds"] for r in design_records
                                   ])))
    (OUT / "sequential_scores.json").write_text(json.dumps(dict(summary=summaries, records=rows), indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score", action="store_true", help="Use separate material truth to score saved predictions")
    options = parser.parse_args()
    with threadpool_limits(1):
        if options.score:
            score_saved()
        else:
            predict()
