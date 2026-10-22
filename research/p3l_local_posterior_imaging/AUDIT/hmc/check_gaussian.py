"""Analytic, full 66-D HMC controls; no physics, data, or training imports."""
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "src"))
from p3l.hmc import PreconditionedHMC
from p3l.kernels import ChainState


class FixedRandom:
    def __init__(self, xi, uniform):
        self.xi, self.uniform = xi, uniform

    def normal(self, size):
        assert size == 66
        return self.xi.copy()

    def random(self):
        return self.uniform


def main():
    rng = np.random.default_rng(20261022)
    d = 66
    # A non-diagonal matrix distinguishes L^-T from the incorrect L^-1 draw.
    S = np.diag(np.linspace(.4, 2., d))
    S[0, 1] = S[1, 0] = .24
    precision = np.linalg.inv(S)
    mu = np.linspace(-.5, .5, d)
    target = lambda q: -.5 * (q - mu) @ precision @ (q - mu)
    gradient = lambda q: -precision @ (q - mu)
    q0 = mu + rng.normal(size=d)
    xi = rng.normal(size=d)
    h, steps = .31, 5
    kernel = PreconditionedHMC(target, gradient, h, S, steps)
    initial = ChainState(q0, target(q0), gradient(q0))
    result = kernel.step(initial, FixedRandom(xi, 1e-300))
    # Independent matrix construction for a quadratic Hamiltonian.
    A = np.block([[np.eye(d), np.zeros((d, d))],
                  [-.5 * h * precision, np.eye(d)]])
    B = np.block([[np.eye(d), h * S],
                  [np.zeros((d, d)), np.eye(d)]])
    leapfrog = np.linalg.matrix_power(A @ B @ A, steps)
    flip = np.diag(np.r_[np.ones(d), -np.ones(d)])
    involution = flip @ leapfrog
    p0 = np.linalg.solve(np.linalg.cholesky(S).T, xi)
    endpoint = involution @ np.r_[q0 - mu, p0]
    q1, p1 = endpoint[:d] + mu, endpoint[d:]
    ratio = target(q1) - target(q0) + .5 * (p0 @ S @ p0 - p1 @ S @ p1)
    assert result.metadata["accepted"]
    np.testing.assert_allclose(result.state.q, q1, atol=2e-14, rtol=2e-14)
    np.testing.assert_allclose(result.metadata["log_mh_ratio"], ratio, atol=2e-13)
    np.testing.assert_allclose(involution @ endpoint, np.r_[q0 - mu, p0], atol=2e-14)
    # Each triangular shear has determinant one; the flip has absolute
    # determinant one. This checks volume without an unnecessary dense det.
    assert np.prod(np.diag(A)) == np.prod(np.diag(B)) == 1.
    assert abs(np.prod(np.diag(flip))) == 1.
    reverse_ratio = target(q0) - target(q1) + .5 * (p1 @ S @ p1 - p0 @ S @ p0)
    np.testing.assert_allclose(reverse_ratio, -ratio, atol=2e-14)
    assert result.metadata["target_evaluations"] == 1
    assert result.metadata["gradient_evaluations"] == steps
    assert result.metadata["leapfrog_steps"] == steps
    center = ChainState(mu, target(mu), gradient(mu))
    threshold_kernel = PreconditionedHMC(target, gradient, .7, S, 2)
    probe = threshold_kernel.step(center, FixedRandom(xi, 1e-300))
    probability = np.exp(probe.metadata["log_mh_ratio"])
    assert 0. < probability < 1.
    assert threshold_kernel.step(center, FixedRandom(xi, probability / 2)).metadata["accepted"]
    assert not threshold_kernel.step(center, FixedRandom(xi, (1 + probability) / 2)).metadata["accepted"]
    # A deliberately poor trajectory must finish before rejecting, keeping g0.
    bad = PreconditionedHMC(target, gradient, 3., S, 3)
    rejected = bad.step(ChainState(q0, target(q0)), FixedRandom(xi, .999999))
    assert not rejected.metadata["accepted"]
    assert rejected.metadata["leapfrog_steps"] == 3
    assert rejected.metadata["large_energy_error"]
    assert rejected.metadata["gradient_evaluations"] == 4
    np.testing.assert_array_equal(rejected.state.q, q0)
    np.testing.assert_array_equal(rejected.state.gradient, gradient(q0))
    # Numerical nonfiniteness must never receive an exact-kernel certificate.
    fail = PreconditionedHMC(target, lambda q: np.full(d, np.nan), .1, S, 2)
    failed = fail.step(initial, FixedRandom(xi, .5))
    assert failed.metadata["divergence"]
    assert failed.metadata["status"].startswith("NUMERICAL_FAILURE")
    assert failed.metadata["gradient_evaluations"] == 1
    assert failed.metadata["attempted_leapfrog_steps"] == 1
    assert failed.metadata["leapfrog_steps"] == 0
    # Stationary moments test all coordinates and the nonzero cross covariance.
    stationary = PreconditionedHMC(target, gradient, .22, S, 7)
    state = ChainState(mu, target(mu), gradient(mu))
    draws, acceptances = [], 0
    for i in range(14500):
        out = stationary.step(state, rng)
        state = out.state
        acceptances += out.metadata["accepted"]
        if i >= 500:
            draws.append(state.q)
    draws = np.asarray(draws)
    mean = draws.mean(axis=0)
    covariance = np.cov(draws.T)
    max_mean_se = float(np.max(np.abs(mean - mu) / np.sqrt(np.diag(S) / len(draws))))
    max_variance_relative_error = float(np.max(np.abs(np.diag(covariance) / np.diag(S) - 1)))
    cross_covariance_error = float(abs(covariance[0, 1] - S[0, 1]))
    assert max_mean_se < 5.5
    assert max_variance_relative_error < .07
    assert cross_covariance_error < .025
    assert acceptances / 14500 > .8
    report = dict(status="PASS_ANALYTIC_GAUSSIAN_ONLY", dimension=d,
                  samples=len(draws), acceptance_rate=acceptances / 14500,
                  max_mean_standard_errors=max_mean_se,
                  max_variance_relative_error=max_variance_relative_error,
                  cross_covariance_error=cross_covariance_error,
                  matrix_endpoint_max_error=float(np.max(np.abs(result.state.q - q1))),
                  log_mh_ratio_error=float(abs(result.metadata["log_mh_ratio"] - ratio)),
                  deterministic_acceptance_threshold=float(probability),
                  reversibility_max_error=float(np.max(np.abs(involution @ endpoint - np.r_[q0 - mu, p0]))),
                  rejection_metadata=rejected.metadata,
                  numerical_failure_metadata=failed.metadata,
                  scope="Analytic Gaussian algebra and stationary moments only; no Maxwell/posterior/mixing claim")
    # Strict JSON preserves the explicit undefined numerical-error fields as null.
    report["numerical_failure_metadata"] = {
        k: (None if isinstance(v, float) and not np.isfinite(v) else v)
        for k, v in failed.metadata.items()}
    output = Path(__file__).with_name("gaussian_results.json")
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
