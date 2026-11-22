"""Replay the native verifier checks without entering the holdout data path."""

import argparse
import json
from pathlib import Path

import numpy as np

from native import Native


PILOT = Path(__file__).resolve().parents[1]


def dummy_profiles():
    rng = np.random.default_rng(2783)
    y = (rng.standard_normal((20, 10)) + 1j*rng.standard_normal((20, 10))) * .15
    yield "global_4x4_psd", y, .2, .5, np.full(16, .0025), np.zeros(16)
    yield "global_4x4_indefinite", y, .2, .75, np.full(16, .0001), np.full(16, 1.e-5)
    yield "global_4x4_cache_reuse", y, .2, .75, np.full(16, .0001), np.full(16, 1.e-5)
    yield ("global_4x4_weyl_reuse", y, .2, .75, np.full(16, .0001),
           1.e-5*(1 + .001*rng.normal(size=16)))
    yield ("nonuniform_4x4_indefinite", y, .2, .75, np.full(16, .0001),
           1.e-5*(1 + .25*rng.normal(size=16)))
    yield ("nonuniform_8x8_signed_psd", y, .2, .5, .0025*(1 + .02*rng.normal(size=64)),
           1.e-7*rng.normal(size=64))
    yield "full_cap_8x8_zero", y, .2, 1., np.full(64, .0001), np.zeros(64)
    yield ("full_cap_8x8_shared_zero", y, .2, 1., .0001*(1 + .06*rng.normal(size=64)),
           .00004*rng.normal(size=64))


def training_profiles():
    source = PILOT / "data/synthetic_observations.npz"
    with np.load(source, allow_pickle=False) as file:
        index = file["object_id"].tolist().index("dualpilot_train_000")
        noise = file["noise_fraction"].tolist().index(.02)
        assert file["split"][index] == "train"
        y, radius = file["Y"][index, noise], float(file["radius"][index, noise])
    rng = np.random.default_rng(9021)
    yield "fixed_uniform_4x4", y, radius, .5, np.full(16, .0025), np.zeros(16)
    yield ("varying_signed_8x8", y, radius, .5, .0025*(1 + .004*rng.normal(size=64)),
           .000002*rng.normal(size=64))
    yield "fixed_uniform_4x4", y, radius, .75, np.full(16, .0001), np.full(16, 1.e-5)
    yield ("varying_signed_8x8", y, radius, .75, .0001*(1 + .005*rng.normal(size=64)),
           .00001*(1 + .02*rng.normal(size=64)))
    yield "fixed_uniform_4x4", y, radius, 1., np.full(16, .0001), np.zeros(16)
    yield ("dense_variation_witness_reuse", y, radius, 1., .0001*(1 + .07*rng.normal(size=16)),
           .000025*rng.normal(size=16))


def main(training):
    native = Native()
    label = "dualpilot_train_000, 2% paired noise" if training else "Nonphysical dummy Y, seed 2783"
    name = "native_training_check.json" if training else "native_verifier_check.json"
    records = []
    candidates = training_profiles() if training else dummy_profiles()
    for case, y, radius, u, mu, lam in candidates:
        result = native.score(y, u, mu, lam, radius=radius)
        records.append(dict(case=case, **result))
        elapsed = result["timing"]["total"]["wall_seconds"]
        print(case, u, result["status"], result["nu_lower"], result["lower"],
              result["upper"], result["radius"], f"{elapsed:.3f}s", flush=True)
        payload = dict(input=label, rows=records, cache=native.cache_summary())
        (PILOT / "results" / name).write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training", action="store_true")
    main(parser.parse_args().training)
