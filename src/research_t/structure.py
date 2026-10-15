"""Truth-free absolute-material topology and protected raster cavity edits.

This module contains no forward solver, observation, truth, or budget access.
The raster topology convention is 4-connected foreground / 8-connected air.
Susceptibility is chi=epsilon_r-1, so air is near chi=0, not near a background.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.ndimage import binary_dilation, binary_erosion, distance_transform_edt, label
from scipy.optimize import linear_sum_assignment

AIR_THRESHOLDS = (0.075, 0.1, 0.125)
FOREGROUND_CONNECTIVITY = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
AIR_CONNECTIVITY = np.ones((3, 3), dtype=int)


def material_grid(chi: np.ndarray, shape: tuple[int, int] = (24, 24)) -> np.ndarray:
    value = np.asarray(chi, dtype=np.complex128)
    if value.size != np.prod(shape):
        raise ValueError(f"chi has {value.size} values, expected grid {shape}")
    value = value.reshape(shape)
    if not np.all(np.isfinite(value)):
        raise ValueError("non-finite material")
    return value


def absolute_masks(chi: np.ndarray, *, threshold: float = 0.1,
                   shape: tuple[int, int] = (24, 24)) -> dict[str, np.ndarray]:
    if threshold <= 0:
        raise ValueError("absolute air threshold must be positive")
    value = material_grid(chi, shape)
    air = np.abs(value) <= threshold
    return {"air": air, "material": ~air,
            "passive": value.imag >= -1e-12,
            "guard_material": (value.real >= 0.3) & (value.imag >= -1e-12)}


def topology(chi: np.ndarray, *, threshold: float = 0.1,
             shape: tuple[int, int] = (24, 24)) -> dict[str, Any]:
    masks = absolute_masks(chi, threshold=threshold, shape=shape)
    foreground, beta0 = label(masks["material"], FOREGROUND_CONNECTIVITY)
    air_labels, air_count = label(masks["air"], AIR_CONNECTIVITY)
    border = np.concatenate((air_labels[0], air_labels[-1],
                             air_labels[:, 0], air_labels[:, -1]))
    exterior_ids = np.unique(border[border != 0])
    all_ids = np.arange(1, air_count + 1)
    enclosed_ids = np.setdiff1d(all_ids, exterior_ids)
    cavity_mask = np.isin(air_labels, enclosed_ids)
    cavities = []
    for component_id in enclosed_ids:
        rc = np.argwhere(air_labels == component_id)
        cavities.append({"label": int(component_id), "cells": len(rc),
                         "centroid_row_col": rc.mean(axis=0).tolist(),
                         "bounding_box_row_col": [rc.min(axis=0).tolist(),
                                                   rc.max(axis=0).tolist()]})
    return {"threshold": float(threshold), "beta_0_material": int(beta0),
            "beta_1_material": int(len(enclosed_ids)),
            "air_components": int(air_count),
            "exterior_air_components": int(len(exterior_ids)),
            "enclosed_air_cavities": cavities, "air_cells": int(masks["air"].sum()),
            "material_cells": int(masks["material"].sum()),
            "foreground_labels": foreground, "air_labels": air_labels,
            "enclosed_air_mask": cavity_mask,
            "convention": "2D material foreground 4-connected; air/background 8-connected"}


def topology_summary(chi: np.ndarray, *, shape: tuple[int, int] = (24, 24),
                     thresholds: tuple[float, ...] = AIR_THRESHOLDS) -> dict[str, Any]:
    """Threshold persistence is a sensitivity diagnostic, not persistent homology."""
    rows = []
    for threshold in thresholds:
        value = topology(chi, threshold=threshold, shape=shape)
        rows.append({k: v for k, v in value.items() if not isinstance(v, np.ndarray)})
    return {"absolute_material_convention": "chi=epsilon_r-1; air iff abs(chi)<=threshold",
            "thresholds": list(thresholds), "threshold_sensitivity": rows,
            "cavity_count_stable_across_thresholds": len({r["beta_1_material"] for r in rows}) == 1,
            "continuous_topology_certificate": False}


def material_distances(a: np.ndarray, b: np.ndarray, *, mass: np.ndarray | float = 1.0,
                       normalization: float = 1.2,
                       shape: tuple[int, int] = (24, 24),
                       air_threshold: float = 0.1,
                       min_phase_changed_cells: int = 6) -> dict[str, Any]:
    """Fixed-unit distances. Normalization never depends on either endpoint."""
    if normalization != 1.2:
        raise ValueError("Research T frozen material L2 normalization is 1.2")
    aa, bb = material_grid(a, shape), material_grid(b, shape)
    weights = np.broadcast_to(np.asarray(mass, float).reshape(shape)
                              if np.asarray(mass).size == np.prod(shape)
                              else np.asarray(mass, float), shape)
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0):
        raise ValueError("positive finite grid mass required")
    delta, denominator = aa - bb, float(weights.sum())
    ma = absolute_masks(aa, threshold=air_threshold, shape=shape)
    mb = absolute_masks(bb, threshold=air_threshold, shape=shape)
    phase_difference = ma["air"] ^ mb["air"]
    ta, tb = topology(aa, threshold=air_threshold, shape=shape), topology(bb, threshold=air_threshold, shape=shape)
    ca, cb = ta["enclosed_air_mask"], tb["enclosed_air_mask"]
    boundary_a = ma["material"] & ~binary_erosion(ma["material"], FOREGROUND_CONNECTIVITY, border_value=0)
    boundary_b = mb["material"] & ~binary_erosion(mb["material"], FOREGROUND_CONNECTIVITY, border_value=0)
    hausdorff = None
    if boundary_a.any() and boundary_b.any():
        hausdorff = float(max(distance_transform_edt(~boundary_b)[boundary_a].max(),
                              distance_transform_edt(~boundary_a)[boundary_b].max()))
    centres_a = np.asarray([c["centroid_row_col"] for c in ta["enclosed_air_cavities"]])
    centres_b = np.asarray([c["centroid_row_col"] for c in tb["enclosed_air_cavities"]])
    matches = []
    if len(centres_a) and len(centres_b):
        centre_distances = np.linalg.norm(centres_a[:, None, :] - centres_b[None, :, :], axis=-1)
        ia, ib = linear_sum_assignment(centre_distances)
        matches = [{"cavity_a": int(i), "cavity_b": int(j),
                    "centroid_distance_cells": float(centre_distances[i, j])} for i, j in zip(ia, ib)]
    return {"mass_weighted_complex_l2": float(np.sqrt(np.sum(weights * abs(delta) ** 2) / denominator)),
            "normalized_material_l2": float(np.sqrt(np.sum(weights * abs(delta) ** 2) / denominator) / normalization),
            "mass_weighted_complex_l1": float(np.sum(weights * abs(delta)) / denominator),
            "real_mae": float(np.sum(weights * abs(delta.real)) / denominator),
            "imag_mae": float(np.sum(weights * abs(delta.imag)) / denominator),
            "phase_changed_cells": int(phase_difference.sum()),
            "phase_symmetric_difference_mass_fraction": float(np.sum(weights * phase_difference) / denominator),
            "enclosed_cavity_symmetric_difference_cells": int((ca ^ cb).sum()),
            "beta_0_change": tb["beta_0_material"] - ta["beta_0_material"],
            "beta_1_change": tb["beta_1_material"] - ta["beta_1_material"],
            "material_boundary_hausdorff_cells": hausdorff,
            "cavity_centroid_minimum_cost_matches": matches,
            "unmatched_cavities_a": len(centres_a) - len(matches),
            "unmatched_cavities_b": len(centres_b) - len(matches),
            "significant_phase_edit": bool(phase_difference.sum() >= min_phase_changed_cells),
            "min_phase_changed_cells": int(min_phase_changed_cells),
            "normalization": normalization, "air_threshold": air_threshold}


@dataclass(frozen=True)
class ProtectedEdit:
    chi: np.ndarray
    core: np.ndarray
    guard: np.ndarray
    metadata: dict[str, Any]


def cavity_birth(chi: np.ndarray, *, center_row_col: tuple[float, float], radius_cells: float,
                 min_feature_radius_cells: float, guard_width_cells: int = 1,
                 grid_spacing: float, shape: tuple[int, int] = (24, 24),
                 guard_min_real: float = 0.3,
                 material_real_bounds: tuple[float, float] = (-0.5, 4.0),
                 material_imag_bounds: tuple[float, float] = (0.0, 0.5)) -> ProtectedEdit:
    """Make an explicit raster air core plus a closed material guard.

    The guard adjustment is part of the paid proposed material, not a free
    correction. A core/guard collision with the grid exterior is rejected.
    """
    if radius_cells < min_feature_radius_cells or min_feature_radius_cells <= 0:
        raise ValueError("cavity radius violates frozen physical minimum feature")
    if guard_width_cells < 1 or grid_spacing <= 0:
        raise ValueError("positive guard width and grid spacing required")
    if not (material_real_bounds[0] <= 0 <= material_real_bounds[1]
            and material_imag_bounds[0] <= 0 <= material_imag_bounds[1]
            and material_real_bounds[0] <= guard_min_real <= material_real_bounds[1]):
        raise ValueError("global material bounds exclude the intended air core or guard")
    value = material_grid(chi, shape).copy()
    rr, cc = np.indices(shape)
    core = (rr - center_row_col[0]) ** 2 + (cc - center_row_col[1]) ** 2 <= radius_cells ** 2
    if core.sum() < 6:
        raise ValueError("sub-significant raster core: fewer than six cells")
    guard = binary_dilation(core, AIR_CONNECTIVITY, iterations=guard_width_cells) & ~core
    protected = core | guard
    if np.any(protected[0]) or np.any(protected[-1]) or np.any(protected[:, 0]) or np.any(protected[:, -1]):
        raise ValueError("protected cavity/guard touches exterior boundary")
    value[core] = 0.0
    value[guard] = np.maximum(value[guard].real, guard_min_real) + 1j * np.maximum(value[guard].imag, 0.0)
    if (np.any(value.real < material_real_bounds[0] - 1e-12)
            or np.any(value.real > material_real_bounds[1] + 1e-12)
            or np.any(value.imag < material_imag_bounds[0] - 1e-12)
            or np.any(value.imag > material_imag_bounds[1] + 1e-12)):
        raise ValueError("proposal violates material bounds")
    return ProtectedEdit(value, core, guard,
                         {"kind": "protected_cavity_birth", "center_row_col": list(center_row_col),
                          "radius_cells": float(radius_cells),
                          "radius_physical": float(radius_cells * grid_spacing),
                          "min_feature_radius_cells": float(min_feature_radius_cells),
                          "min_feature_radius_physical": float(min_feature_radius_cells * grid_spacing),
                          "guard_width_cells": int(guard_width_cells), "guard_min_real_chi": guard_min_real,
                          "core_cells": int(core.sum()), "guard_cells": int(guard.sum()),
                          "truth_used": False, "raster_pilot_only": True})


def protected_bounds(core: np.ndarray, guard: np.ndarray, *, shape: tuple[int, int] = (24, 24),
                     real_bounds: tuple[float, float] = (-0.5, 4.0),
                     imag_bounds: tuple[float, float] = (0.0, 0.5),
                     guard_min_real: float = 0.3) -> tuple[np.ndarray, np.ndarray]:
    """Return [Re(flat chi), Im(flat chi)] bounds for a shared real corrector."""
    core, guard = np.asarray(core, bool).reshape(shape), np.asarray(guard, bool).reshape(shape)
    if np.any(core & guard):
        raise ValueError("core and guard overlap")
    lower = np.r_[np.full(np.prod(shape), real_bounds[0]), np.full(np.prod(shape), imag_bounds[0])]
    upper = np.r_[np.full(np.prod(shape), real_bounds[1]), np.full(np.prod(shape), imag_bounds[1])]
    n = int(np.prod(shape))
    idx = np.flatnonzero(core)
    lower[idx] = upper[idx] = lower[n + idx] = upper[n + idx] = 0.0
    lower[np.flatnonzero(guard)] = max(real_bounds[0], guard_min_real)
    lower[n + np.flatnonzero(guard)] = max(imag_bounds[0], 0.0)
    if np.any(lower > upper):
        raise ValueError("empty protected legal domain")
    return lower, upper


def check_protected(chi: np.ndarray, core: np.ndarray, guard: np.ndarray, *,
                    shape: tuple[int, int] = (24, 24), guard_min_real: float = 0.3,
                    tolerance: float = 1e-10) -> dict[str, Any]:
    value = material_grid(chi, shape)
    core, guard = np.asarray(core, bool).reshape(shape), np.asarray(guard, bool).reshape(shape)
    if not core.any() or not guard.any() or np.any(core & guard):
        return {"preserved": False, "reason": "empty or overlapping core/guard"}
    core_error = float(np.max(np.abs(value[core])))
    guard_re = float(np.min(value[guard].real))
    guard_im = float(np.min(value[guard].imag))
    air_topology = topology(value, shape=shape)
    enclosed = air_topology["enclosed_air_mask"]
    enclosed_core = bool(np.all(enclosed[core]))
    ok = core_error <= tolerance and guard_re >= guard_min_real - tolerance and guard_im >= -tolerance and enclosed_core
    return {"preserved": bool(ok), "core_max_abs_chi": core_error,
            "guard_min_real_chi": guard_re, "guard_min_imag_chi": guard_im,
            "all_core_enclosed": enclosed_core,
            "reason": None if ok else "core, guard contrast, passivity, or enclosure failed"}
