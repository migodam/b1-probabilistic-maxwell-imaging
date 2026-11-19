"""Generate, fit, and score eight paired held-out six-material noise controls."""
import argparse
from dataclasses import dataclass
import json

import numpy as np
from threadpoolctl import threadpool_limits

from data14 import (ACTIONS, BOX, GAIN_SD, HALF_SPLIT_FIRST, R10_OLD_EPS, R12, ROOT, SIGMA,
                    FrequencyTM, OldCase, cn_radius, complex_reads, encode_complex, measured_reads)

SEED = 1491625
N_PAIRS = 8
LOW_SIGMA = SIGMA
HIGH_SIGMA = R10_OLD_EPS / cn_radius(1., 3)
DEFAULT = ROOT / "results/validation_public.json"
FULL_SPLIT = ACTIONS["split_0.02"]


@dataclass
class PrivateTrial:
    pair_id: str
    geometry: dict
    material: np.ndarray
    gain: complex
    standard_noise: np.ndarray
    overlap_rejections: int


def physical_layout(geometry, rng):
    rejections = 0
    radius = geometry["radius"]
    while True:
        coordinates = {key: (np.asarray(geometry[key]) + rng.normal(0., sd, np.shape(geometry[key])))
                       for key, sd in (("centers", .003), ("sources", .006), ("receivers", .006))}
        centers = coordinates["centers"]
        separations = np.linalg.norm(centers[:, None] - centers[None, :], axis=-1)
        np.fill_diagonal(separations, np.inf)
        nodes = np.concatenate((coordinates["sources"], coordinates["receivers"]))
        node_gaps = np.linalg.norm(nodes[:, None] - centers[None, :], axis=-1)
        if np.min(separations) > 2 * radius and np.min(node_gaps) > radius:
            return dict(geometry, **{key: value.tolist() for key, value in coordinates.items()}), rejections
        rejections += 1


def private_trials():
    layouts = json.loads((R12 / "results/spectral_new_layouts.json").read_text())
    for i, key in enumerate(np.random.SeedSequence(SEED).spawn(N_PAIRS)):
        geo_rng, material_rng, noise_rng, gain_rng = map(np.random.default_rng, key.spawn(4))
        geometry, rejected = physical_layout(layouts[i % len(layouts)]["geometry"], geo_rng)
        material = material_rng.uniform(*BOX, size=6)
        noise = (noise_rng.normal(size=(3, 3)) + 1j * noise_rng.normal(size=(3, 3))) / np.sqrt(2.)
        amount = GAIN_SD * np.sqrt(gain_rng.uniform())
        gain = 1. + amount * np.exp(2j * np.pi * gain_rng.uniform())
        yield PrivateTrial(f"val6-{i:02d}", geometry, material, gain, noise, rejected)


def generate():
    cases, rejections = [], {}
    for trial in private_trials():
        model = FrequencyTM(trial.geometry)
        clean = np.array([model.response(trial.material, delta)[0]
                          for delta in (0., *FULL_SPLIT.deltas)])
        for level, sigma in (("low", LOW_SIGMA), ("high", HIGH_SIGMA)):
            noisy = trial.gain * clean + sigma * trial.standard_noise / np.sqrt([1., .5, .5])[:, None]
            cases.append(dict(case_id=f"{trial.pair_id}-{level}", geometry=trial.geometry, sigma=sigma,
                              old_ri=encode_complex(noisy[0])))
            cases[-1]["split_0.02_ri"] = encode_complex(noisy[1:])
        rejections[trial.pair_id] = trial.overlap_rejections
    DEFAULT.write_text(json.dumps(dict(cases=cases), indent=2, allow_nan=False) + "\n")
    settings = dict(seed=SEED, independent_geometry_material_pairs=N_PAIRS,
                    levels=dict(low=LOW_SIGMA, high=HIGH_SIGMA),
                    material_box=BOX, material_sampling="independent coordinate uniform",
                    jitter_sd=dict(centers=.003, sources=.006, receivers=.006),
                    layout_order="cycle the six R12 layouts once, then repeat layouts 0 and 1",
                    overlap_rejections=rejections, gain_radius=GAIN_SD)
    (ROOT / "results/validation_design.json").write_text(json.dumps(settings, indent=2) + "\n")


def load_validation_cases(path=DEFAULT):
    cases = []
    for row in json.loads(path.read_text())["cases"]:
        old = complex_reads(row["old_ri"])
        sigma = row["sigma"]
        level = row["case_id"].rsplit("-", 1)[1]
        old.flags.writeable = False
        cases.append(OldCase(row["case_id"], level, row["geometry"], old, sigma, GAIN_SD,
                             np.sqrt(2.) * cn_radius(sigma, old.size), np.sqrt(2.) * GAIN_SD,
                             hard_gain_radius=GAIN_SD))
    return tuple(cases)


def load_validation_reads(path=DEFAULT):
    return {row["case_id"] + "/split_0.02": complex_reads(row["split_0.02_ri"])
            for row in json.loads(path.read_text())["cases"]}


def fit():
    from measure_recover import solve
    cases, reads = load_validation_cases(), load_validation_reads()
    records = []
    with threadpool_limits(limits=1):
        for case in cases:
            for action, label in ((HALF_SPLIT_FIRST, "half_first"), (FULL_SPLIT, "full_split")):
                result = solve(case, action, measured_reads(reads, case.case_id, action),
                               np.empty((0, 6)), "ordinary8", 780)
                records.append(dict(case_id=case.case_id, cohort=case.cohort, label=label,
                                    action=action.name, mode="ordinary8", **result))
    path = ROOT / "results/validation_predictions.json"
    path.write_text(json.dumps(records, indent=2, allow_nan=False) + "\n")


def score():
    from measure_score import compatibility, physical_score, summarize
    pred = json.loads((ROOT / "results/validation_predictions.json").read_text())
    private = {t.pair_id: t for t in private_trials()}
    cases = {c.case_id: c for c in load_validation_cases()}
    reads = load_validation_reads()
    no_alternatives = np.empty((0, 6))
    for row in pred:
        case = cases[row["case_id"]]
        trial = private[case.case_id.rsplit("-", 1)[0]]
        action = ACTIONS[row["action"]]
        new = measured_reads(reads, case.case_id, action)
        model = FrequencyTM(case.geometry)
        response = np.array([model.response(trial.material, delta)[0] for delta in (0., *action.deltas)])
        _, truth_valid = compatibility(case, action, new, response)
        row.update(physical_score(case, model, action, new, np.asarray(row["prediction"]),
                                  trial.material, no_alternatives))
        row.update(energy=action.energy, exposures=action.exposures,
                   bandwidth_pct=100 * action.bandwidth, truth_in_common_contract=truth_valid,
                   total_factors=row["fit_factors"], total_rhs=row["fit_rhs"], design_factors=0,
                   polynomial_paths=0, design_fallback=None, total_wall_s=row["fit_wall_s"])
    stats = []
    for level in ("low", "high"):
        for label in ("half_first", "full_split"):
            selected = [row for row in pred if (row["cohort"], row["label"]) == (level, label)]
            stats.append(dict(level=level, action=label, **summarize(selected),
                              rms_le_point1=sum(r["rmse"] <= .1 for r in selected),
                              rms_le_point5=sum(r["rmse"] <= .5 for r in selected),
                              rms_le_one=sum(r["rmse"] <= 1. for r in selected)))
    out = ROOT / "results"
    (out / "validation_scores.json").write_text(json.dumps(pred, indent=2) + "\n")
    (out / "validation_summary.json").write_text(json.dumps(stats, indent=2) + "\n")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("stage", choices=("generate", "fit", "score"))
    {"generate": generate, "fit": fit, "score": score}[cli.parse_args().stage]()
