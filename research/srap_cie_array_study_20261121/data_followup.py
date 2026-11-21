#!/usr/bin/env python3
"""Generate four continuous materials and four connected thin-bridge materials."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform

import data2d as base
import numpy as np
from scipy.ndimage import label

ARRAYS = ("PW180_10x20", "PW120_10x20", "PW120_CENTERED_10x20", "PW180_FLIP_TX_10x20")
MATERIAL_SEED = 20261121031
ORDER_SEED = 20261121032


def bump(x, y, center, radii, transition):
    rho = ((x - center[0]) / radii[0]) ** 2 + ((y - center[1]) / radii[1]) ** 2
    t = np.clip((1 - rho) / transition, 0, 1)
    return t**2 * (3 - 2 * t)


def continuous_material(x, y, rng):
    support = bump(x, y, (0, 0), rng.uniform([.66, .57], [.76, .68]), .54)
    low = rng.uniform(.63, 1.08) + .23 * (x + .3) + .11 * y
    a = np.array([-.22, .19]) + rng.uniform(-.08, .08, 2)
    b = np.array([.27, -.16]) + rng.uniform(-.065, .065, 2)
    c = np.array([-.07, -.27]) + rng.uniform(-.08, .08, 2)
    high_a = rng.uniform(2.45, 3.35) * bump(x, y, a, rng.uniform(.21, .32, 2), 1)
    high_b = rng.uniform(1.95, 3.0) * bump(x, y, b, rng.uniform(.18, .29, 2), .75)
    middle = rng.uniform(.55, 1.25) * bump(x, y, c, rng.uniform(.21, .29, 2), 1)
    m = support * (low + high_a + high_b + middle)
    if m.max() < 2.8 or not np.any((m > .4) & (m < 1.2)):
        raise RuntimeError("Frozen smooth material lacks the intended contrast range")
    return m


def bridge_material(x, y, rng, variant, fine_side):
    centers = np.array(((-.44, -.18), (.36, .23), (.32, -.46)))
    centers = centers[:2 + variant % 2] + rng.uniform(-.06, .06, (2 + variant % 2, 2))
    radii = rng.uniform(.15, .21, (len(centers), 2))
    peak = rng.uniform(3.1, 4.1)
    levels = peak * np.linspace(1, .77, len(centers))
    lobes = np.maximum.reduce([
        level * base.ellipse(x, y, *center, *radius)
        for level, center, radius in zip(levels, centers, radii, strict=True)
    ])
    half_width = rng.uniform(.037, .050)
    bends = rng.uniform(-.10, .10, (len(centers) - 1, 2))
    links = [
        base.history.stroke(x, y, [start, .5 * (start + end) + bend, end], half_width)
        for start, end, bend in zip(centers[:-1], centers[1:], bends, strict=True)
    ]
    m = np.maximum(lobes, .78 * peak * np.logical_or.reduce(links))
    if label((m > 0).reshape(fine_side, fine_side))[1] != 1:
        raise RuntimeError("Frozen thin-bridge material has disconnected fine-grid support")
    if label((lobes > 0).reshape(fine_side, fine_side))[1] != len(centers):
        raise RuntimeError("Frozen thin-bridge material has lobes that touch without a bridge")
    return m


def material(centers, index, family, variant, fine_side):
    rng = np.random.default_rng(np.random.SeedSequence([MATERIAL_SEED, index]))
    angle = rng.uniform(-np.pi, np.pi)
    translation = rng.uniform(-.055, .055, 2)
    x, y = base.rotated_xy(centers, angle, translation).T
    if family == 0:
        return continuous_material(x, y, rng)
    return bridge_material(x, y, rng, variant, fine_side)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=base.HERE / "data_followup")
    root = parser.parse_args().output
    for folder in ("public", "truth", "heldout"):
        (root / folder).mkdir(parents=True, exist_ok=True)
    manifest = root / "cases.json"
    costs = root / "generation_cost.json"
    cases = {entry["case"]: entry for entry in json.loads(manifest.read_text())} if manifest.exists() else {}
    receipt = json.loads(costs.read_text()) if costs.exists() else {
        "platform": platform.platform(), "numpy": np.__version__, "threads_per_process": 1,
        "setups": {}, "forwards": {}, "measured": {},
    }
    assignment = np.random.default_rng(ORDER_SEED).permutation(8)
    jobs = [(i, j, f"Cx{i:02d}") for i, j in enumerate(assignment, 1)]
    pending = [(i, j, obj) for i, j, obj in jobs if not all(
        base.present(root, cases, base.synthetic_name("followup", obj, a), "truth") for a in ARRAYS
    )]
    if not pending:
        return
    ops = base.SyntheticOperators(96, ARRAYS)
    receipt["setups"].setdefault("followup/96", []).append(ops.cost)
    families = ("continuous_high_low", "thin_bridge_multibody")
    for i, j, obj in pending:
        family, variant = divmod(j, 4)
        chi = material(ops.fine_centers, i, family, variant, ops.fine)
        base.synthetic(root, cases, receipt, ops, "followup", obj, 100 + i, chi, families[family])


if __name__ == "__main__":
    main()
