"""Score follow-up experiments through the frozen inverse algorithm on new geometries."""

import argparse
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from threadpoolctl import threadpool_limits

from acquire import COARSE, fft_field, material, relative
from designs import Design
from holdout import cases, forward_batch
from ndt import Air, BACKGROUND, Host, PLATE_PERMITTIVITY, grid, host
from port_online import PortFields
from sequential import PILOT, RX, SHOT_ENERGY, TX, fit, posterior_trace
from sparse_features import estimate, prepare, real_stack, remove_shared_gain

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from run_inversion import score
from src.inversion import fullwave
from src.roi_refine import region


OUT = ROOT / "results"
GEOMETRY_SEEDS = (20261121, 20261122)
DATA_NOISE_SEEDS = (20261123, 20261124)
FORECAST_SEED = 20261125
RANDOM_SEED = 20261126
FREQUENCIES = (12., PILOT, 24., 30.)
FINE = 72
QUADRATURE = 8
PILOT_PORTS = Design("pilot", (PILOT,), TX, RX)
ACTIONS = {
    "fixed_12": Design("fixed_12", (12.,), TX, RX),
    "fixed_24": Design("fixed_24", (24.,), TX, RX),
    "fixed_30": Design("fixed_30", (30.,), TX, RX),
    "fixed_angle": Design("fixed_angle", (PILOT,), (1, 3, 5, 7), RX),
}
ARMS = ("pilot_only", "forecast_risk", "fisher", "uniform_random", *ACTIONS)
FALLBACK = "fixed_angle"


def write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def output(kind, seed, extension="npz"):
    return OUT / f"pilot_risk_{kind}_{seed}.{extension}"


def generate():
    with np.load(OUT / "acquisition.npz", allow_pickle=False) as original:
        sigma = .02 * float(original["reference_amplitude"])
    for seed, noise_seed in zip(GEOMETRY_SEEDS, DATA_NOISE_SEEDS):
        wall, cpu = perf_counter(), process_time()
        definitions = cases(seed)
        fine = np.array([material(item["primitives"], FINE, QUADRATURE) for item in definitions])
        cells_per_side = FINE // COARSE
        truth = fine.reshape(len(fine), COARSE, cells_per_side, COARSE, cells_per_side).mean(axis=(2, 4))
        clean, blank, diagnostics = [], [], []
        for f in FREQUENCIES:
            total, _, background, _, diag = forward_batch(f, truth, fine, PLATE_PERMITTIVITY, 1., host)
            clean.append(total[:, RX, :])
            blank.append(background[list(RX)])
            diagnostics.append(diag)
        clean = np.stack(clean, axis=1)
        rng = np.random.default_rng(noise_seed)
        noise = (rng.standard_normal(clean.shape) + 1j * rng.standard_normal(clean.shape)) / np.sqrt(2)
        ids = [item["id"] for item in definitions]
        np.savez_compressed(output("observations", seed), case_ids=ids, freq_ghz=FREQUENCIES,
                            clean=clean, noise=noise, blank=blank, sigma_2pct=sigma, receiver_ids=RX)
        np.savez_compressed(output("truth", seed), case_ids=ids, truth=truth, truth_fine=fine)
        write_json(output("truth", seed, "json"), dict(
            case_definitions=definitions, shape_seed=seed, noise_seed=noise_seed,
            fine_grid_side=FINE, geometry_subpixels_per_cell=QUADRATURE**2,
            wall_seconds=perf_counter() - wall, cpu_seconds=process_time() - cpu, diagnostics=diagnostics,
        ))
        print(f"Generated {seed}, {len(ids)} objects, {perf_counter() - wall:.2f} s", flush=True)


def local_tangent(calibrated, design, selected_ids, x):
    selected = Host(calibrated.Eh[np.ix_(selected_ids, design.tx)],
                    calibrated.Sh[np.ix_(design.rx, selected_ids)],
                    calibrated.Gh[np.ix_(selected_ids, selected_ids)], calibrated.xy[selected_ids])
    return fullwave([selected], x, jacobian=True)[1]


def fisher_choice(seed_image, host24, blanks, sd):
    ids, _ = region(seed_image)
    if not len(ids):
        return FALLBACK, {}
    q = seed_image.ravel()[ids]
    base = local_tangent(host24[PILOT], PILOT_PORTS, ids, q)
    before = posterior_trace(base, blanks["pilot"], sd)
    difference = {}
    for key, design in ACTIONS.items():
        candidate = local_tangent(host24[design.freqs[0]], design, ids, q)
        joint = np.vstack((base, candidate))
        after = posterior_trace(joint, np.r_[blanks["pilot"], blanks[key]], sd)
        difference[key] = (before - after) / before
    return max(difference, key=difference.get), difference


def prediction_from_proxy(ports, design, ids, image):
    block = ports[design.freqs[0]].extract(design, ids)
    fine = np.repeat(image.ravel()[ids], 4)
    return fullwave([block], fine)[0].ravel()


def risk_choice(seed_image, y0, ports, prepared, blanks, sd, noise):
    ids = np.flatnonzero(seed_image.ravel())
    if not len(ids):
        return FALLBACK, {}, 0
    risks = {}
    for key, design in ACTIONS.items():
        forecast = prediction_from_proxy(ports, design, ids, seed_image)
        y = np.vstack((np.broadcast_to(y0[:, None], noise.shape), forecast[:, None] + sd * noise))
        projected = remove_shared_gain(y, np.r_[blanks["pilot"], blanks[key]])
        predicted = estimate(*prepared[key], real_stack(projected) * (np.sqrt(2) / sd))[0]
        error = predicted - seed_image.ravel()[None, :]
        risks[key] = float(np.mean(np.sqrt(np.mean(error**2, axis=1))))
    return min(risks, key=risks.get), risks, len(ids)


def compare_forecast_projection(ports, image):
    """Check the 24 variables tie to the correct four fine cells inside the VIE."""
    selected = ACTIONS["fixed_30"]
    ids = np.flatnonzero(image.ravel())
    sub = prediction_from_proxy(ports, selected, ids, image)
    air = Air(48, selected.freqs[0])
    op = air.convolution()
    bg, _ = fft_field(air, op, np.full(48**2, BACKGROUND))
    m = BACKGROUND + image.repeat(2, axis=0).repeat(2, axis=1).ravel()[:, None]
    field, _ = fft_field(air, op, m, bg)
    reference = air.S @ (m * field - BACKGROUND * bg)
    return relative(sub, reference[np.ix_(selected.rx, selected.tx)].ravel())


def prepare_models():
    prep_start = perf_counter()
    ports, host24, costs = {}, {}, {}
    for f in FREQUENCIES:
        started = perf_counter()
        ports[f] = PortFields(48, f, False)
        mid = perf_counter()
        host24[f] = host(24, f)
        costs[f] = dict(port48_seconds=mid - started, host24_seconds=perf_counter() - mid)
    return ports, host24, dict(frequencies=costs, total_seconds=perf_counter() - prep_start)


def predict_pack(seed, ports, host24, check_projection):
    wall, cpu = perf_counter(), process_time()
    with np.load(output("observations", seed), allow_pickle=False) as source:
        pack = {name: source[name] for name in source.files}
    sd = float(pack["sigma_2pct"]) / np.sqrt(SHOT_ENERGY)
    fi = {f: i for i, f in enumerate(pack["freq_ghz"])}
    if tuple(pack["receiver_ids"]) != RX:
        raise ValueError("Data and design receiver coordinates differ")
    labels = pack["case_ids"]
    names = {"pilot": PILOT_PORTS, **ACTIONS}
    blanks = {key: pack["blank"][fi[d.freqs[0]]][:, d.tx].ravel() for key, d in names.items()}
    linear = {key: ports[d.freqs[0]].tangent(d) for key, d in names.items()}
    pilot = prepare(linear["pilot"], blanks["pilot"], sd, True)
    prep = {key: prepare(np.vstack((linear["pilot"], linear[key])), np.r_[blanks["pilot"], blanks[key]], sd, True)
            for key in ACTIONS}

    rng = np.random.default_rng(FORECAST_SEED)
    samples = (rng.standard_normal((32, 2)) + 1j * rng.standard_normal((32, 2))) / np.sqrt(2)
    forecast_draws = np.concatenate((samples, -samples), axis=1)
    random = np.random.default_rng(np.random.SeedSequence([RANDOM_SEED, seed])).choice(list(ACTIONS), len(labels))

    def measure(key, index):
        design = names[key]
        f = fi[design.freqs[0]]
        return (pack["clean"][index, f][:, design.tx] + sd * pack["noise"][index, f][:, design.tx]).ravel()

    records = []
    projection_error = None
    images = {arm: [] for arm in ARMS}
    for i, case_id in enumerate(labels):
        tick = perf_counter()
        y0 = measure("pilot", i)
        initial, evidence, threshold = fit(pilot, y0, blanks["pilot"], sd)
        after_initial = perf_counter()
        risk, risks, support = risk_choice(initial, y0, ports, prep, blanks, sd, forecast_draws)
        after_risk = perf_counter()
        fisher, scores = fisher_choice(initial, host24, blanks, sd)
        after_fisher = perf_counter()
        images["pilot_only"].append(initial)

        def fit_revealed(action):
            y = np.r_[y0, measure(action, i)]
            return fit(prep[action], y, np.r_[blanks["pilot"], blanks[action]], sd)[0]

        # Both choices are made before reading any post-pilot target signal.
        images["forecast_risk"].append(fit_revealed(risk))
        after_risk_reconstruction = perf_counter()
        images["fisher"].append(fit_revealed(fisher))
        if check_projection and projection_error is None and support:
            projection_error = compare_forecast_projection(ports, initial)

        fixed = {action: fit_revealed(action) for action in ACTIONS}
        for key, image in fixed.items():
            images[key].append(image)
        images["uniform_random"].append(fixed[random[i]])
        for arm, key in (("forecast_risk", risk), ("fisher", fisher)):
            if not np.array_equal(images[arm][-1], fixed[key]):
                raise RuntimeError("Candidate and fixed solve disagree on the same revealed measurement")
        records.append(dict(
            case=str(case_id), forecast_choice=risk, fisher_choice=fisher, random_choice=str(random[i]),
            pilot_evidence=evidence, pilot_threshold=threshold, proxy_nonzero_pixels=support,
            fine_proxy_forward_subregion_states=4 * support,
            forecast_rmse=risks, fisher_fractional_trace_reduction=scores,
            pilot_seconds=after_initial - tick, forecast_selection_seconds=after_risk - after_initial,
            fisher_selection_seconds=after_fisher - after_risk,
            forecast_final_reconstruction_seconds=after_risk_reconstruction - after_fisher,
        ))
        print(f"{seed} {case_id}: proxy {support} cells, risk {risk}, Fisher {fisher}", flush=True)
    np.savez_compressed(output("predictions", seed), arms=ARMS, case_ids=labels,
                        images=np.asarray([images[arm] for arm in ARMS]))
    return dict(cases=records, wall_seconds=perf_counter() - wall, cpu_seconds=process_time() - cpu,
                exact_proxy_forward_vs_fft_relative_error=projection_error)


def predict():
    wall, cpu = perf_counter(), process_time()
    ports, host24, setup = prepare_models()
    poses = {key: dict(frequency_ghz=d.freqs[0], transmitters=d.tx, receivers=d.rx)
             for key, d in {"pilot": PILOT_PORTS, **ACTIONS}.items()}
    batches = {seed: predict_pack(seed, ports, host24, seed == GEOMETRY_SEEDS[0]) for seed in GEOMETRY_SEEDS}
    write_json(OUT / "pilot_risk_choices.json", dict(
        measurement_plan=poses, shots=8, complex_readouts=64, energy=1., fallback=FALLBACK,
        forecast_noise_seed=FORECAST_SEED, forecast_draws=4, antithetic_pairs=2, random_seed=RANDOM_SEED,
        setup=setup, batches=batches, wall_seconds=perf_counter() - wall, cpu_seconds=process_time() - cpu,
    ))


def score_saved():
    wall, cpu = perf_counter(), process_time()
    all_summaries, all_scores = {}, []
    with (OUT / "pilot_risk_choices.json").open() as f:
        decisions = json.load(f)["batches"]
    for seed in GEOMETRY_SEEDS:
        with np.load(output("truth", seed), allow_pickle=False) as t:
            truth, labels = (t[key] for key in ("truth", "case_ids"))
        with np.load(output("predictions", seed), allow_pickle=False) as r:
            if not np.array_equal(labels, r["case_ids"]):
                raise ValueError("Prediction/truth order differs")
            images = dict(zip(r["arms"], r["images"]))
        meta = json.loads(output("truth", seed, "json").read_text())["case_definitions"]
        centers = {d["id"]: np.array([[p["x_mm"], p["y_mm"]] for p in d["primitives"]])
                   for d in meta if "spacing_mm" in d}
        results = {}
        for arm, photos in images.items():
            rows = []
            for cid, image, reference in zip(labels, photos, truth):
                s = score(image.ravel(), reference.ravel(), grid(24) * 1000, centers.get(cid))
                rows.append(dict(seed=seed, arm=str(arm), case=str(cid), empty=not np.any(reference), **s))
            details = [r for r in rows if not r["empty"]]
            pairs = sorted((r for r in details if "true_separation_mm" in r), key=lambda r: r["true_separation_mm"])
            results[str(arm)] = dict(
                mean_rmse=float(np.mean([r["rmse_eps"] for r in details])),
                median_rmse=float(np.median([r["rmse_eps"] for r in details])),
                near_pair=pairs[0], wide_pair=pairs[-1],
                empty_false_detection=next(r["detected"] for r in rows if r["empty"]),
            )
            all_scores.extend(rows)
        choice_counts = {}
        for mode in ("forecast", "fisher", "random"):
            # Nonempty subset is for reporting only and loaded after predictions.
            picked = [record[mode + "_choice"] for record, ref in zip(decisions[str(seed)]["cases"], truth) if np.any(ref)]
            choice_counts[mode] = {key: picked.count(key) for key in ACTIONS}
        all_summaries[seed] = dict(arms=results, nonempty_choice_counts=choice_counts,
                                   random_expected_mean_rmse=float(np.mean([results[a]["mean_rmse"] for a in ACTIONS])))
    write_json(OUT / "pilot_risk_scores.json", dict(batches=all_summaries, cases=all_scores,
                                                    wall_seconds=perf_counter() - wall,
                                                    cpu_seconds=process_time() - cpu))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=("generate", "predict", "score"))
    args = p.parse_args()
    with threadpool_limits(1):
        {"generate": generate, "predict": predict, "score": score_saved}[args.stage]()
