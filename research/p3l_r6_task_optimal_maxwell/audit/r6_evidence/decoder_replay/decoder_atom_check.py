#!/usr/bin/env python3
"""Decoder-only prior diagnostic; zero Maxwell calls and no posterior inference.

Reproduce with: python decoder_atom_check.py
Uses exact committed prior.py/config/geometry copied to decoder_sources/.
Scope: independent original untruncated Gaussian64 x Uniform2 prior.
Replicate SE and ranges are diagnostics, not rigorous probability intervals.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import time

import numpy as np
import scipy
from scipy.special import ndtr, ndtri
from scipy.stats import qmc
import os
# Single thread requested by launch environment; threadpoolctl unavailable.


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "decoder_sources"
COMMIT = "ad71a86ef3c28bf48049199f733080a0f2b95e33"
SOURCE_PATHS = {
    "prior.py": "src/p1_material/prior.py",
    "ONLINE_TARGET.json": "research/p1_material_target_20261018/configs/ONLINE_TARGET.json",
    "geometry.json": "research/p1_material_target_20261018/inputs/geometry.json",
}


def summarize(values):
    values = np.asarray(values, dtype=float)
    return {
        "replicate_values": values.tolist(),
        "mean": float(values.mean()),
        "replicate_min": float(values.min()),
        "replicate_max": float(values.max()),
        "replicate_sample_SD": float(values.std(ddof=1)),
        "replicate_SE_of_mean": float(values.std(ddof=1) / np.sqrt(len(values))),
        "qualification": "NUMERICAL_DIAGNOSTIC_NOT_RIGOROUS_INTERVAL",
    }


def main():
    start = time.perf_counter()
    cfg = json.loads((SRC / "ONLINE_TARGET.json").read_text())
    geometry = json.loads((SRC / "geometry.json").read_text())
    spec = importlib.util.spec_from_file_location("committed_decoder_prior", SRC / "prior.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    prior = module.MaterialPrior(cfg, geometry["points"])
    roi = np.flatnonzero(prior.mask)
    A = prior.A[roi]
    offset = A[:, 0]
    assert np.all(offset == cfg["offset_SD"])
    a0 = float(cfg["offset_SD"])
    assert a0 > 0
    tau = float(cfg["phase_transition_halfwidth"])
    mean = float(cfg["field_mean"])
    assert mean == 0.0  # Needed only for exact host/pure symmetry claim.
    B = A[:, 1:]
    assert B.shape == (256, 63)

    # A direct source-level check at two known complete plateau witnesses.
    q_host = np.r_[[-2.0], np.zeros(63), [0.37, 0.61]]
    q_pure = np.r_[[2.0], np.zeros(63), [0.37, 0.61]]
    chi_host, info_host = prior.decode(q_host)
    chi_pure, info_pure = prior.decode(q_pure)
    checks = {
        "host_fraction_all_zero": bool(np.all(info_host["fraction"] == 0)),
        "host_chi_all_host": bool(np.all(chi_host == prior.host)),
        "pure_ROI_fraction_all_one": bool(np.all(info_pure["fraction"][roi] == 1)),
        "outside_ROI_fraction_zero": bool(np.all(info_pure["fraction"][~prior.mask] == 0)),
        "host_decoder_J_max_abs": float(np.abs(prior.derivative(q_host)).max()),
        "pure_decoder_J_numerical_rank": int(np.linalg.matrix_rank(
            np.r_[prior.derivative(q_pure).real, prior.derivative(q_pure).imag]
        )),
    }
    assert all(checks[k] for k in checks if isinstance(checks[k], bool))
    assert checks["host_decoder_J_max_abs"] == 0
    assert checks["pure_decoder_J_numerical_rank"] == 2

    levels = []
    # First two fixed groups followed by a pre-fixed 16x65536 extension.
    # The extension was authorized after observing 3.37% difference between the
    # first two group means, to address that concrete numerical risk.
    # Every scramble is independent; host/pure pairing is NOT independent replication.
    plan = [(12, 8), (13, 8), (16, 16)]
    for power, repetitions in plan:
        host, pure, seeds = [], [], []
        for replicate in range(repetitions):
            seed = 2026102700 + 100 * power + replicate
            seeds.append(seed)
            engine = qmc.Sobol(d=63, scramble=True, bits=52, seed=seed)
            u = engine.random_base2(power)
            # No endpoint clipping: abort rather than silently change the prior.
            assert np.all((u > 0) & (u < 1))
            z_rest = ndtri(u)
            # Bounded decoder-only audit adaptation; largest b is 4096x256 float64.
            host_sum = pure_sum = 0.0
            for start_batch in range(0, len(z_rest), 4096):
                b = mean + z_rest[start_batch:start_batch+4096] @ B.T
                host_sum += ndtr((-tau - b.max(axis=1)) / a0).sum()
                pure_sum += ndtr((b.min(axis=1) - tau) / a0).sum()
            host.append(float(host_sum / len(z_rest)))
            pure.append(float(pure_sum / len(z_rest)))
        symmetric = (np.array(host) + np.array(pure)) / 2
        levels.append({
            "sobol_points_per_replicate": 2**power,
            "replicates": repetitions,
            "independent_scramble_seeds": seeds,
            "host": summarize(host),
            "pure": summarize(pure),
            "symmetry_pooled_host_equals_pure": summarize(symmetric),
        })

    sources = []
    for name, path in SOURCE_PATHS.items():
        sources.append({
            "repository_path": path,
            "commit": COMMIT,
            "url": f"https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/{COMMIT}/{path}",
            "local_blake2b": hashlib.blake2b((SRC / name).read_bytes()).hexdigest(),
        })
    result = {
        "status": "DECODER_ONLY_PRIOR_NUMERICAL_DIAGNOSTIC",
        "scope": "Original untruncated Gaussian64 x Uniform2 prior; not frozen 2D submodel",
        "target_id": cfg["target_id"],
        "prior_id": cfg["prior_id"],
        "sources": sources,
        "dimensions": {"latent": 66, "Gaussian_field": 64, "integrated_offset": 1,
                       "remaining_Gaussian_QMC": 63, "grid_cells": 576, "ROI_cells": int(len(roi))},
        "parameters": {"offset_SD": a0, "field_SD": cfg["field_SD"],
                       "field_mean": mean, "tau": tau, "correlation_rolloff_length": cfg["correlation_rolloff_length"]},
        "formula_host": "E_b Phi((-tau-max_i b_i)/a0), b_i=field_mean+(A_ROI[:,1:] z_rest)_i",
        "formula_pure": "E_b Phi((min_i b_i-tau)/a0)",
        "analytic_positive_mass_proof": "For every finite z_rest the finite ROI maximum/minimum is finite, so each conditional Gaussian CDF is strictly positive. Its expectation is positive. The independent U2 property coordinates integrate to one. With zero field mean, z_rest and -z_rest have the same law, hence host and pure masses are exactly equal.",
        "material_measure_interpretation": {
            "all_host": "Rank-0 material atom at the all-host chi vector; atom prior mass estimated here.",
            "all_pure_ROI": "Rank-2 continuous material layer because both uniform material property coordinates remain free; not a point atom.",
            "posterior": "NOT_COMPUTED; host posterior mass is p_host*likelihood(host)/full66_evidence, which is unknown.",
        },
        "checks": checks,
        "levels": levels,
        "fine_minus_coarse_symmetry_mean": float(levels[1]["symmetry_pooled_host_equals_pure"]["mean"] - levels[0]["symmetry_pooled_host_equals_pure"]["mean"]),
        "extension_minus_original_fine_symmetry_mean": float(levels[2]["symmetry_pooled_host_equals_pure"]["mean"] - levels[1]["symmetry_pooled_host_equals_pure"]["mean"]),
        "extension_rationale": "A pre-fixed independent 16x65536 scramble group was added after the two initial groups differed by 3.37%; no optional convergence stopping. All levels are retained.",
        "cost": {"new_Maxwell_states": 0, "new_Maxwell_RHS": 0,
                 "Gaussian_Rao_Blackwell_points": sum(repetitions*2**power for power,repetitions in plan),
                 "wall_seconds": time.perf_counter()-start},
        "versions": {"numpy": np.__version__, "scipy": scipy.__version__},
        "limitations": [
            "Finite scrambled-QMC replicate variability is a numerical diagnostic, not a rigorous enclosure or posterior certificate.",
            "No directed rounding or full deterministic quadrature error bound was used.",
            "No Maxwell forward, Jacobian, likelihood, posterior normalization, posterior mass, or high-SNR experiment was run.",
            "The frozen two-dimensional prior mass 0.22086048 is not used anywhere in this computation.",
        ],
    }
    (ROOT / "decoder_atom_results.json").write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n")
    print(json.dumps({"checks": checks, "levels": levels,
                      "fine_minus_coarse_symmetry_mean": result["fine_minus_coarse_symmetry_mean"],
                      "cost": result["cost"]}, indent=2))


if __name__ == "__main__":
    main()
