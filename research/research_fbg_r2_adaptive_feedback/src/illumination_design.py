"""Physical disk-pulse illumination geometry and source-energy bookkeeping.

The incident field is exactly the R4 scalar-TM formula, .25j H0^(1)(k r).
Changing source positions changes physical data; right-invertible source mixing
of a fixed collection does not create independent information.
"""
from __future__ import annotations

import numpy as np

SOURCE_RADIUS = 1.6


def source_positions(p: int, design: str, old_sources: np.ndarray,
                     seed: int = 20261030) -> np.ndarray:
    if p not in (3, 6, 12, 24):
        raise ValueError("p must be 3, 6, 12, or 24")
    old = np.asarray(old_sources, float)
    if old.shape != (12, 2):
        raise ValueError("Expected twelve original source positions")
    if design == "uniform":
        if p <= 12:
            return old[np.arange(0, 12, 12 // p)].copy()
        angles = 2*np.pi*np.arange(24)/24
        positions = SOURCE_RADIUS*np.column_stack((np.cos(angles), np.sin(angles)))
        # Exact old positions are retained at even indices; odd indices are new.
        positions[::2] = old
        return positions
    if design == "random":
        rng = np.random.default_rng(seed)
        # One jittered 12-angle geometry; p=3,6 are prefixes of the same 12.
        angles12 = 2*np.pi*(np.arange(12) + rng.uniform(-.20, .20, 12))/12
        order = rng.permutation(12)
        angles12 = angles12[order]
        pos12 = SOURCE_RADIUS*np.column_stack((np.cos(angles12), np.sin(angles12)))
        if p <= 12:
            return pos12[:p].copy()
        # Twelve additional, genuinely new physical source locations.
        angles_extra = 2*np.pi*(np.arange(12) + .5 + rng.uniform(-.20, .20, 12))/12
        extra = SOURCE_RADIUS*np.column_stack((np.cos(angles_extra), np.sin(angles_extra)))
        return np.vstack((pos12, extra))
    raise ValueError("design must be uniform or random")


def incident_field(centers: np.ndarray, sources: np.ndarray, k: float) -> np.ndarray:
    from scipy.special import hankel1
    r = np.linalg.norm(np.asarray(centers)[:, None, :] - np.asarray(sources)[None, :, :], axis=2)
    if np.any(r <= 0):
        raise ValueError("Source coincides with a cell center")
    return .25j*hankel1(0, float(k)*r)


def energy_scale(p: int, protocol: str, reference_p: int = 12) -> float:
    if protocol == "per_source":
        return 1.
    if protocol == "fixed_total":
        return float(np.sqrt(reference_p/p))
    raise ValueError("protocol must be per_source or fixed_total")


def source_mixing(E0: np.ndarray, Y: np.ndarray, M: np.ndarray):
    """Invertible column mixing; require whitening in comparisons.

    If independent source noise has covariance sigma^2 I before mixing,
    the mixed source covariance is sigma^2 M^H M. Residual @ inv(M)
    restores the original residual and is one valid whitening transform.
    """
    E0, Y, M = map(np.asarray, (E0, Y, M))
    p = E0.shape[1]
    if Y.shape[1] != p or M.shape != (p, p):
        raise ValueError("Incompatible source dimensions")
    if not np.isfinite(np.linalg.cond(M)) or np.linalg.cond(M) > 1e8:
        raise ValueError("Mixing matrix must be stably invertible")
    return E0@M, Y@M


def whitened_relative_residual(pred_mixed: np.ndarray, obs_mixed: np.ndarray,
                               M: np.ndarray) -> float:
    inverse = np.linalg.inv(np.asarray(M))
    err = (np.asarray(pred_mixed) - np.asarray(obs_mixed))@inverse
    ref = np.asarray(obs_mixed)@inverse
    return float(np.linalg.norm(err)/max(np.linalg.norm(ref), 1e-300))
