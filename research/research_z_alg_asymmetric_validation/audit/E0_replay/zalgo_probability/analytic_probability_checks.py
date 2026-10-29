"""Reproduce analytic_checks.json from the Z-ALG probability toy checks.

These are low-dimensional analytic probability examples, not additional
Maxwell experiments. NumPy and SciPy are the only external dependencies.

Usage:
    python analytic_probability_checks.py
    python analytic_probability_checks.py --output /absolute/path/checks.json
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad


def compute_checks():
    # Q,D independent Uniform[1,2], conditioned on QD=2.
    # Conditional Q density is 1/(q log(2)).
    qnorm = quad(lambda q: 1 / q, 1, 2, epsabs=1e-13)[0]
    eq = quad(lambda q: q / (q * math.log(2)), 1, 2, epsabs=1e-13)[0]

    # H(x)=x^3-3x, uniform material prior on [-2,2], exact observation u=0.
    roots = np.array([-math.sqrt(3), 0, math.sqrt(3)])
    weights = 1 / np.abs(3 * roots**2 - 3)
    weights = weights / weights.sum()

    # Attainment of the sharp posterior TV bound for unnormalized likelihood
    # ratio w in [exp(-eta), exp(eta)]. The base law has its high-w atom with
    # probability t; its normalizer is sqrt(a*b).
    eta = 0.7
    a, b = math.exp(-eta), math.exp(eta)
    t = math.sqrt(a) / (math.sqrt(a) + math.sqrt(b))
    normalizer = (1 - t) * a + t * b
    tv = 0.5 * ((1 - t) * abs(a / normalizer - 1)
                + t * abs(b / normalizer - 1))

    # A separate unit-box probability toy: rho(x1,x2)=2*x1 on [0,1]^2.
    # Invariant coordinates are s=x1+x2, q=x1*x2. Ordered gap g=|x1-x2|
    # removes the 1/sqrt(s^2-4q) density singularity; its density becomes s.
    # This unit-box toy must not be confused with the physical E4 box [0.5,4]^2.
    mass = quad(
        lambda s: quad(lambda g: s, 0, min(s, 2 - s), epsabs=1e-13)[0],
        0, 2, epsabs=1e-13,
    )[0]
    branch_probability = quad(
        lambda s: quad(lambda g: (s + g) / 2, 0, min(s, 2 - s), epsabs=1e-13)[0],
        0, 2, epsabs=1e-13,
    )[0]

    return {
        "QD_conditional_normalization": qnorm,
        "QD_E_Q": eq,
        "QD_independent_product_mean": eq**2,
        "QD_true_product": 2,
        "cubic_roots": roots.tolist(),
        "cubic_coarea_weights": weights.tolist(),
        "cubic_event_abs_x_lt_1": float(weights[1]),
        "cubic_second_moment": float(np.dot(weights, roots**2)),
        "sharp_TV_eta": eta,
        "TV_extremal": tv,
        "bound_tanh_half_eta": math.tanh(eta / 2),
        "mirror_sum_product_prior_mass": mass,
        "mirror_prior_event_x1_gt_x2": branch_probability,
        "mirror_direct_event_reference": 2 / 3,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("analytic_checks.json"))
    args = parser.parse_args()
    result = compute_checks()
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
