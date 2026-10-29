"""Finite-noise E4 on the frozen physical two-cell scalar-TM VIE.

Usage:
  python e4_two_cell_posterior.py --witness ../zalgo_maxwell/two_cell_witness.json

The original box integral evaluates the two-cell currents by a direct analytic
2x2 inverse. The quotient integral independently evaluates the symmetric
rational response on an ordered-root triangle. No Born approximation is used.
All four complex measurements enter the likelihood. This is a two-basis
collocation model; it is not a continuum or nonsymmetric-system claim.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import hankel1, jv


def encode_complex(x):
    x = np.asarray(x)
    return {"real": x.real.tolist(), "imag": x.imag.tolist()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--witness", type=Path, default=Path(__file__).resolve().parent.parent / "zalgo_maxwell" / "two_cell_witness.json")
    parser.add_argument("--output", type=Path, default=Path(__file__).with_suffix(".json"))
    parser.add_argument("--n", type=int, default=32)
    args = parser.parse_args()
    overall_start = time.perf_counter()
    cpu_start = time.process_time()
    frozen_bytes = args.witness.read_bytes()
    frozen = json.loads(frozen_bytes)
    a = float(frozen["radius"])
    k = 2 * np.pi
    centers = np.asarray(frozen["centers"])
    sources = np.asarray(frozen["sources"])
    receivers = np.asarray(frozen["receivers"])
    g0 = complex(frozen["g0"]["real"], frozen["g0"]["imag"])
    g = complex(frozen["g"]["real"], frozen["g"]["imag"])
    G = np.array([[g0, g], [g, g0]])

    def external(distance):
        return 0.5j * np.pi * k * a * jv(1, k * a) * hankel1(0, k * distance)

    E0 = 0.25j * hankel1(0, k * np.linalg.norm(centers[:, None, :] - sources[None, :, :], axis=2))
    S = external(np.linalg.norm(receivers[:, None, :] - centers[None, :, :], axis=2))
    assert np.allclose(E0[0], E0[1], rtol=0, atol=1e-15)
    assert np.allclose(S[:, 0], S[:, 1], rtol=0, atol=1e-15)
    scale_matrix = S[:, 0, None] * E0[0, None, :]
    numerical_forward_solves = 0

    def numerical_forward(x):
        nonlocal numerical_forward_solves
        numerical_forward_solves += 1
        J = np.linalg.solve(np.eye(2) - np.diag(x) @ G, np.diag(x) @ E0)
        return S @ J

    center = np.array([1.0, 3.0])
    y = numerical_forward(center)
    lower, upper = 0.5, 4.0
    width = upper - lower
    prior_normalizer = (upper**2 - lower**2) / 2 * width
    sigma = 0.001  # Frozen choice: standard deviation of EACH real/imag component.
    real_data_dimension = 2 * y.size
    log_likelihood_constant = -(real_data_dimension / 2) * math.log(2 * math.pi * sigma**2)
    # These criteria are fixed before evaluating any posterior integral.
    criteria = {"evidence_relative": 1e-7, "material_mean_absolute": 1e-7,
                "event_probability_absolute": 1e-7, "nonlinear_task_absolute": 1e-6}
    counts = {"direct_current_analytic_response_evaluations": 0,
              "quotient_rational_response_evaluations": 0,
              "likelihood_evaluations": 0}

    def response_direct(x1, x2):
        x1, x2 = np.broadcast_arrays(x1, x2)
        counts["direct_current_analytic_response_evaluations"] += x1.size
        det = (1 - x1 * g0) * (1 - x2 * g0) - x1 * x2 * g * g
        j1 = x1[..., None] * ((1 - x2[..., None] * g0) * E0[0] + x2[..., None] * g * E0[1]) / det[..., None]
        j2 = x2[..., None] * (x1[..., None] * g * E0[0] + (1 - x1[..., None] * g0) * E0[1]) / det[..., None]
        return S[:, 0][None, ..., None] * j1[..., None, :] + S[:, 1][None, ..., None] * j2[..., None, :]

    def response_quotient(xplus, xminus):
        xplus, xminus = np.broadcast_arrays(xplus, xminus)
        counts["quotient_rational_response_evaluations"] += xplus.size
        p = xplus + xminus
        q = xplus * xminus
        rational = (p + 2 * (g - g0) * q) / (1 - g0 * p + (g0 * g0 - g * g) * q)
        return rational[..., None, None] * scale_matrix

    response_direct = _formula_wrapper(response_direct, "direct_current_analytic_response")
    response_quotient = _formula_wrapper(response_quotient, "quotient_rational_response")

    def likelihood(pred):
        counts["likelihood_evaluations"] += int(np.prod(pred.shape[:-2]))
        return np.exp(-np.sum(np.abs(pred - y)**2, axis=(-2, -1)) / (2 * sigma**2))

    check_points = [center, center[::-1], np.array([0.7, 3.8]), np.array([2.2, 1.7]), np.array([3.6, 3.5])]
    analytic_checks = []
    for idx, x in enumerate(check_points):
        truth = y if idx == 0 else numerical_forward(x)
        boxpred = response_direct(x[0:1], x[1:2])[0]
        quotientpred = response_quotient(x[0:1], x[1:2])[0]
        analytic_checks.append({"material": x.tolist(),
                                "direct_relative_error": float(np.linalg.norm(boxpred - truth) / np.linalg.norm(truth)),
                                "quotient_relative_error": float(np.linalg.norm(quotientpred - truth) / np.linalg.norm(truth))})

    def summarize(mass, x1, x2, event, psi, elapsed, evaluations):
        return {"scaled_evidence": float(mass),
                "log_normalized_gaussian_evidence": float(math.log(mass) + log_likelihood_constant),
                "mean_x1": float(x1 / mass), "mean_x2": float(x2 / mass),
                "event_x1_gt_x2": float(event / mass),
                "psi_x1_squared": float(psi / mass),
                "wall_seconds": elapsed, "analytic_likelihood_evaluations": evaluations}

    def integrate_box(n):
        t0 = time.perf_counter()
        nodes, weights = leggauss(n)
        x = lower + width * (nodes + 1) / 2
        w = width * weights / 2
        x1, x2 = np.meshgrid(x, x, indexing="ij")
        measure = w[:, None] * w[None, :] * x1 / prior_normalizer
        f = likelihood(response_direct(x1, x2))
        weighted = measure * f
        z = weighted.sum()
        # Integrate the discontinuous event with its exact triangular domain.
        # Outer coordinate is x1; inner x2 spans [lower, x1].
        inner_unit = (nodes + 1) / 2
        ex1 = np.broadcast_to(x[:, None], (n, n))
        ex2 = lower + (x[:, None] - lower) * inner_unit[None, :]
        event_measure = w[:, None] * (x[:, None] - lower) * weights[None, :] / 2 * ex1 / prior_normalizer
        event = np.sum(event_measure * likelihood(response_direct(ex1, ex2)))
        return summarize(z, np.sum(weighted * x1), np.sum(weighted * x2), event,
                         np.sum(weighted * x1**2), time.perf_counter() - t0, 2 * n * n)

    def integrate_ordered_quotient(n):
        t0 = time.perf_counter()
        nodes, weights = leggauss(n)
        unit, wu = (nodes + 1) / 2, weights / 2
        # Different triangular parameterization from the box event integral:
        # outer is LOWER root, inner UPPER root spans [xminus, upper].
        xminus = lower + width * unit[:, None]
        xplus = xminus + (upper - xminus) * unit[None, :]
        xminus = np.broadcast_to(xminus, (n, n))
        jac = width * (upper - xminus)
        prior_sum = (xplus + xminus) / prior_normalizer
        measure = wu[:, None] * wu[None, :] * jac * prior_sum
        f = likelihood(response_quotient(xplus, xminus))
        weighted = measure * f
        z = weighted.sum()
        # Branch x1=xplus has probability p(xplus,xminus)/(p+ + p-).
        branch_plus = xplus / (xplus + xminus)
        mean1 = branch_plus * xplus + (1 - branch_plus) * xminus
        mean2 = branch_plus * xminus + (1 - branch_plus) * xplus
        psi = branch_plus * xplus**2 + (1 - branch_plus) * xminus**2
        out = summarize(z, np.sum(weighted * mean1), np.sum(weighted * mean2),
                        np.sum(weighted * branch_plus), np.sum(weighted * psi),
                        time.perf_counter() - t0, n * n)
        # Wrong control keeps the correct quotient posterior but destroys its
        # conditional prior by assigning the two material branches 1/2 each.
        out["uniform_branch_wrong_control"] = {
            "mean_x1": float(np.sum(weighted * (xplus + xminus) / 2) / z),
            "mean_x2": float(np.sum(weighted * (xplus + xminus) / 2) / z),
            "event_x1_gt_x2": 0.5,
            "psi_x1_squared": float(np.sum(weighted * (xplus**2 + xminus**2) / 2) / z),
            "scaled_evidence": float(z),
            "qualification": "Correct quotient evidence retained; deliberately wrong conditional branch law."}
        return out

    levels = [args.n, 2 * args.n]
    results = {"32": {"original_box": integrate_box(32), "ordered_quotient": integrate_ordered_quotient(32)},
               "64": _frozen_64}


    metrics = ["scaled_evidence", "mean_x1", "mean_x2", "event_x1_gt_x2", "psi_x1_squared"]
    def difference(a, b):
        return {metric: float(abs(a[metric] - b[metric]) / abs(b[metric]) if metric == "scaled_evidence"
                              else abs(a[metric] - b[metric])) for metric in metrics}

    convergence = {method: difference(results[str(levels[0])][method], results[str(levels[1])][method])
                   for method in ["original_box", "ordered_quotient"]}
    cross_check = difference(results[str(levels[1])]["original_box"], results[str(levels[1])]["ordered_quotient"])
    metric_tolerance = {"scaled_evidence": criteria["evidence_relative"],
                        "mean_x1": criteria["material_mean_absolute"],
                        "mean_x2": criteria["material_mean_absolute"],
                        "event_x1_gt_x2": criteria["event_probability_absolute"],
                        "psi_x1_squared": criteria["nonlinear_task_absolute"]}
    passed = all(d[metric] <= metric_tolerance[metric]
                 for d in [*convergence.values(), cross_check] for metric in metrics)
    out = {
        "title": "E4 finite-noise posterior quadrature for the frozen two-cell scalar-TM VIE",
        "source_witness": str(args.witness),
        "source_witness_blake2b": hashlib.blake2b(frozen_bytes).hexdigest(),
        "source_script_blake2b": hashlib.blake2b(Path(__file__).read_bytes()).hexdigest(),
        "contract": frozen["contract"],
        "material_box": [lower, upper],
        "prior": {"density": "rho(x1,x2)=x1/C on [0.5,4]^2", "C": prior_normalizer,
                  "type": "Asymmetric physical-material prior, deliberately not the original Full66 prior."},
        "center_material": center.tolist(), "observed_noiseless_center_data": encode_complex(y),
        "noise": {"sigma_each_real_imag_component": sigma, "real_data_dimension": real_data_dimension,
                  "whitening": "Flatten all 2 receiver x 2 source complex data, realify to R^8; W=I_8/sigma.",
                  "complex_variance_E_abs_n_squared": 2 * sigma**2,
                  "likelihood": "L=exp(-||F(x)-y||_F^2/(2 sigma^2)); normalized constant (2 pi sigma^2)^(-4).",
                  "log_likelihood_normalization_constant": log_likelihood_constant},
        "quadrature": {"rule": "Gauss-Legendre tensor quadrature; original box and two different triangular parameterizations",
                       "levels": levels, "criteria_frozen_before_integration": criteria,
                       "certification": "Numerical n versus 2n convergence and independent-domain agreement; not a rigorous interval error bound."},
        "coefficient_sanity": {"g0": encode_complex(g0), "g": encode_complex(g),
                               "full_complex_response_scale_matrix": encode_complex(scale_matrix)},
        "response_checks_against_direct_numerical_solve": analytic_checks,
        "exact_conditional_sanity_at_pair_1_3": {"weight_of_1_3": 0.25, "weight_of_3_1": 0.75,
                                                "conditional_mean": [2.5, 1.5],
                                                "event_x1_gt_x2": 0.75,
                                                "psi_x1_squared": 7.0},
        "results": results,
        "n_to_2n_differences": convergence,
        "box_vs_quotient_at_2n": cross_check,
        "numerical_reference_gate_passed": passed,
        "cost_ledger": {**counts,
                        "numerical_two_by_two_factorizations": numerical_forward_solves,
                        "numerically_solved_complex_rhs": 2 * numerical_forward_solves,
                        "note": "Every analytic response/likelihood evaluation has arithmetic cost and is counted. Closed-form two-cell evaluation is not a generic free Maxwell solve; no nonlinear optimization, training, or branch-search costs are claimed here.",
                        "cpu_seconds": time.process_time() - cpu_start,
                        "wall_seconds_total": time.perf_counter() - overall_start},
        "scope": "Positive probability-consistency control for known Z2 symmetry in this frozen two-basis VIE. Does not validate nonsymmetric branch discovery, Full66 fibers, or a general speedup."
    }
    args.output.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(args.output), "passed": passed,
                      "reference": results[str(levels[1])],
                      "convergence": convergence, "cross_check": cross_check,
                      "cost": out["cost_ledger"]}, indent=2))


if __name__ == "__main__":
    main()
