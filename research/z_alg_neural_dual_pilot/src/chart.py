"""Shared data-independent multiplier parameterization."""

import numpy as np
from scipy.special import expit


MU_SCALE = 1.0e-3
MU_FLOOR = 0.05
MU_OFFSET = float(np.log(np.expm1(1.0 - MU_FLOOR)))
LAMBDA_BIAS = 0.5
LAMBDA_SLOPE = 2.0


def profile(a, b):
    """Convert signed coordinates to physical block multipliers."""
    mu = MU_SCALE * (MU_FLOOR + np.logaddexp(0.0, MU_OFFSET + np.asarray(a)))
    lam = MU_SCALE * (LAMBDA_BIAS + LAMBDA_SLOPE * np.asarray(b))
    return mu, lam


def coordinates(mu, lam):
    """Invert the chart for a proposed warm start above the mu floor."""
    target = np.asarray(mu, dtype=float) / MU_SCALE - MU_FLOOR
    if np.any(target <= 0):
        raise ValueError("Warm start mu must be above the chart floor.")
    a = target + np.log(-np.expm1(-target)) - MU_OFFSET
    b = (np.asarray(lam, dtype=float) / MU_SCALE - LAMBDA_BIAS) / LAMBDA_SLOPE
    return a, b


def pullback(a, grad_mu, grad_lam):
    """Transform gradients with respect to physical weights to chart coordinates."""
    return MU_SCALE * expit(MU_OFFSET + a) * grad_mu, MU_SCALE * LAMBDA_SLOPE * grad_lam
