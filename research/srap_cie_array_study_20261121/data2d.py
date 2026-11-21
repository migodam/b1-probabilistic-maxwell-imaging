#!/usr/bin/env python3
"""Generate the SRAP planar scalar-TM data without running inversion."""
from __future__ import annotations

import os

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import argparse
import json
from pathlib import Path
import platform
import sys
import tempfile
import time
from dataclasses import dataclass

import numpy as np
from scipy.special import hankel1

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
WORKSPACE = Path(os.environ.get("LA_WORKSPACE", HERE.parents[1]))
LEGACY = WORKSPACE / "research" / "limited_aperture_routed_dr_20261120"
if (LEGACY / "run_la_study.py").is_file():
    sys.path.insert(0, str(LEGACY))
    import run_la_study as history
else:
    import legacy_la as history

from research.research_max_struct_r8.src.fresnel_adapter import incident_fit, load_file, operators
from research.research_max_struct_r8.src.paper_physics import EqualAreaGreen

FRESNEL_2005_DIR = Path(os.environ.get("FRESNEL_2005_DIR", str(history.FRESNEL_DIR))).expanduser().resolve()
FREQUENCY = 400e6
LIGHT_SPEED = 3e8  # Historical synthetic convention; Fresnel uses exact c.
WIDTH = 2.0
SIDE = 32
RANGE = 3.75
MATERIAL_SEED = 20261121021
FAMILY_SEED = 20261121022
NOISE_SEED = 20261121023
NOISE_FRACTION = 0.10
FAMILY_NAMES = (
    "notched_eccentric_shell",
    "offset_elliptical_rings",
    "elliptical_triplet",
    "offset_angular_arms",
    "eccentric_layered_core",
    "hollow_polygon_and_bars",
)
CONFIRM_ARRAYS = (
    "PW360_10x20",
    "PW180_10x20",
    "PW120_10x20",
    "PW180_FLIP_RX_10x20",
    "LS_LINE_SAME_10x20",
    "LS_LINE_OPPOSITE_10x20",
    "LS_ARC_10x20",
)
DEV_ARRAYS = ("PW180_10x20", "PW120_10x20", "PW180_7x14", "PW120_7x14")
DEV_CASES = ("la120_austria_eps3p5", "la120_eccentric_shell_eps4p0")


@dataclass
class Array:
    name: str
    receiver_points: np.ndarray
    tx_kind: str
    tx: np.ndarray
    aperture: int | None


def unit_from_degrees(angle):
    rad = np.deg2rad(angle)
    return np.stack((np.cos(rad), np.sin(rad)), axis=-1)


def line_coordinates(n, y):
    x = RANGE * (2 * (np.arange(n) + 0.5) / n - 1)
    return np.c_[x, np.full(n, y)]


def configuration(name):
    if name.startswith("PW"):
        ns, nr = (int(s) for s in name.rsplit("_", 1)[1].split("x"))
        aperture = int(name[2:5])
        if aperture == 360:
            tx_degrees = np.arange(ns) * 360.0 / ns
            rx_degrees = np.arange(nr) * 360.0 / nr
        else:
            tx_limit, rx_limit = {180: (162, 171), 120: (108, 117)}[aperture]
            tx_degrees = np.linspace(0, tx_limit, ns)
            rx_degrees = np.linspace(0, rx_limit, nr)
        if "CENTERED" in name:
            tx_degrees += 27
            rx_degrees += 27
        if "FLIP_TX" in name:
            tx_degrees += 180
        if "FLIP_RX" in name:
            rx_degrees += 180
        return Array(name, RANGE * unit_from_degrees(rx_degrees), "plane_wave",
                     tx_degrees[:, None], aperture)
    tx = line_coordinates(10, RANGE)
    rx = line_coordinates(20, -RANGE if name == "LS_LINE_OPPOSITE_10x20" else RANGE)
    if name == "LS_ARC_10x20":
        rx = RANGE * rx / np.linalg.norm(rx, axis=1)[:, None]
    return Array(name, rx, "line_source", tx, None)


def rotated_xy(points, angle, translation):
    c, s = np.cos(angle), np.sin(angle)
    return (points - translation) @ np.array([[c, -s], [s, c]])


def ellipse(x, y, cx, cy, a, b):
    return ((x - cx) / a) ** 2 + ((y - cy) / b) ** 2 <= 1


def random_material(points, index, family):
    rng = np.random.default_rng(np.random.SeedSequence([MATERIAL_SEED, index]))
    angle = rng.uniform(-np.pi, np.pi)
    translation = rng.uniform(-0.075, 0.075, 2)
    scale = rng.uniform(0.85, 1.07, 2)
    x, y = (rotated_xy(points, angle, translation) / scale).T
    contrast = rng.uniform(2.4, 4.1)
    aspect = rng.uniform(0.8, 1.20)
    shift = rng.uniform(-0.065, 0.065, 5)

    if family == 0:
        hole = ellipse(x, y, 0.09 + shift[0], -0.10 + shift[1], 0.32, 0.29)
        outer = ellipse(x, y, 0, 0, 0.58, 0.53 * aspect)
        angle = np.arctan2(y, x)
        notch = np.abs(angle - (0.9 + 4 * shift[2])) < rng.uniform(0.12, 0.30)
        island = ellipse(x, y, 0.10 + shift[0], -0.15 + shift[1],
                         0.10 + shift[3] / 2, 0.11 + shift[4] / 2)
        return contrast * (outer & ~hole & ~notch) + 0.72 * contrast * island

    if family == 1:
        first = ellipse(x, y, -0.32 + shift[0], 0.28 + shift[1], 0.29, 0.25 * aspect)
        first &= ~ellipse(x, y, -0.34 + shift[0], 0.29 + shift[1], 0.15 + shift[4], 0.13)
        second = ellipse(x, y, 0.34 + shift[2], -0.29 + shift[3], 0.28, 0.30 / aspect)
        second &= ~ellipse(x, y, 0.32 + shift[2], -0.31 + shift[3], 0.13, 0.15 + shift[4])
        return contrast * np.maximum(first.astype(float), rng.uniform(0.72, 0.98) * second)

    if family == 2:
        first = ellipse(x, y, -0.36 + shift[0], 0.22 + shift[1], 0.20, 0.23 * aspect)
        second = ellipse(x, y, 0.34 + shift[2], 0.23 + shift[3], 0.21 * aspect, 0.20)
        third = ellipse(x, y, 0.02 + shift[4], -0.35 + shift[1], rng.uniform(0.16, 0.27), 0.25)
        return contrast * np.maximum.reduce([first.astype(float),
                                               rng.uniform(0.68, 0.95) * second,
                                               rng.uniform(0.72, 0.97) * third])

    if family == 3:
        first = history.stroke(x, y, [[-0.62, -0.20 + shift[0]], [-0.20, -0.20],
                                      [-0.01 + shift[1], 0.09], [0.54, 0.09 + shift[2]]],
                               rng.uniform(0.09, 0.15))
        second = history.stroke(x, y, [[-0.45, 0.32 + shift[3]], [0.02, 0.32],
                                       [0.29, 0.47 + shift[4]]], rng.uniform(0.07, 0.12))
        third = ellipse(x, y, 0.39 + shift[1], -0.32 + shift[2], 0.19, 0.13 * aspect)
        return contrast * np.maximum.reduce([first.astype(float),
                                               rng.uniform(0.55, 0.94) * second,
                                               rng.uniform(0.75, 0.96) * third])

    if family == 4:
        outer = ellipse(x, y, 0, 0, 0.62, 0.51 * aspect)
        cavity = ellipse(x, y, shift[0], shift[1], 0.40, 0.34 * aspect)
        core = ellipse(x, y, shift[2], -0.07 + shift[3], 0.22 + shift[4], 0.18 / aspect)
        return rng.uniform(1.3, 2.4) * (outer & ~cavity & ~core) + contrast * core

    box = (np.abs(x + 0.23 + shift[0]) < 0.29) & (np.abs(y + 0.06 + shift[1]) < 0.32 * aspect)
    box &= ~((np.abs(x + 0.23 + shift[2]) < 0.14) & (np.abs(y + 0.06 + shift[1]) < 0.16))
    bar1 = (np.abs(x - 0.38 + shift[3]) < 0.12) & (np.abs(y - 0.10 + shift[4]) < 0.32)
    bar2 = (np.abs(x - 0.22 + shift[1]) < 0.18) & (np.abs(y + 0.48 + shift[2]) < 0.08)
    return contrast * np.maximum.reduce([box.astype(float), 0.82 * bar1, 0.66 * bar2])


def area_average(m, fine, side):
    """Average fine piecewise-constant cells over coarse cell areas."""
    edges = np.linspace(0, fine, side + 1)
    coarse_left = edges[:-1, None]
    coarse_right = edges[1:, None]
    fine_left = np.arange(fine)[None, :]
    overlap = np.maximum(0, np.minimum(coarse_right, fine_left + 1) -
                         np.maximum(coarse_left, fine_left))
    weights = overlap / (fine / side)
    return (weights @ m.reshape(fine, fine) @ weights.T).ravel()


def line_wave(centers, tx, k):
    distance = np.linalg.norm(centers[:, None] - tx[None, :], axis=2)
    return 0.25j * hankel1(0, k * distance)


class SyntheticOperators:
    """One Green setup and one solve per distinct illumination of an object."""

    def __init__(self, fine, names):
        wall, cpu = time.perf_counter(), time.process_time()
        arrays = [configuration(name) for name in names]
        self.arrays = arrays
        self.fine = fine
        self.k = 2 * np.pi * FREQUENCY / LIGHT_SPEED
        self.coarse_centers = history.grid(SIDE, WIDTH)
        self.fine_centers = history.grid(fine, WIDTH)
        self.coarse_radius = WIDTH / SIDE / np.sqrt(np.pi)
        self.fine_radius = WIDTH / fine / np.sqrt(np.pi)
        self.green = EqualAreaGreen(self.fine_centers, (fine, fine), self.fine_radius, self.k)

        self.tx_ids = {}
        coarse_e, fine_e = [], []
        source_gain = {}
        offset = 0
        for kind in dict.fromkeys(a.tx_kind for a in arrays):
            same_kind = [a for a in arrays if a.tx_kind == kind]
            tx, ids = np.unique(np.concatenate([a.tx for a in same_kind]), axis=0,
                                return_inverse=True)
            for a, index in zip(same_kind, np.split(ids, np.cumsum(
                    [len(a.tx) for a in same_kind])[:-1]), strict=True):
                self.tx_ids[a.name] = index + offset
            if kind == "plane_wave":
                directions = unit_from_degrees(tx[:, 0])
                coarse_e.append(np.exp(1j * self.k * (self.coarse_centers @ directions.T)))
                fine_e.append(np.exp(1j * self.k * (self.fine_centers @ directions.T)))
                source_gain[kind] = np.ones(len(tx))
            else:
                reference_e = line_wave(history.grid(128, WIDTH), tx, self.k)
                gain = 1 / np.sqrt(np.mean(np.abs(reference_e) ** 2, axis=0))
                coarse_e.append(line_wave(self.coarse_centers, tx, self.k) * gain)
                fine_e.append(line_wave(self.fine_centers, tx, self.k) * gain)
                source_gain[kind] = gain
            offset += len(tx)
        self.coarse_e = np.concatenate(coarse_e, axis=1)
        self.fine_e = np.concatenate(fine_e, axis=1)
        self.source_gain = np.concatenate(list(source_gain.values()))
        rx, ids = np.unique(np.concatenate([a.receiver_points for a in arrays]), axis=0,
                            return_inverse=True)
        self.rx_ids = dict(zip(names, np.split(ids, np.cumsum(
            [len(a.receiver_points) for a in arrays])[:-1]), strict=True))
        self.coarse_s = operators(
            self.coarse_centers, None, rx, WIDTH / SIDE, self.k, dense=False,
        )["S"]
        self.fine_s = operators(
            self.fine_centers, None, rx, WIDTH / fine, self.k, dense=False,
        )["S"]
        self.cost = {
            "wall_s": time.perf_counter() - wall,
            "cpu_s": time.process_time() - cpu,
            "green_kernel_fft_setups": 1,
            "fine_source_columns": self.fine_e.shape[1],
            "fine_receiver_rows": self.fine_s.shape[0],
            "fine_side": fine,
        }


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                     delete=False) as file:
        json.dump(obj, file, indent=2, ensure_ascii=False)
        file.write("\n")
    Path(file.name).replace(path)


def write_npz(path, **arrays):
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as file:
        np.savez(file, **arrays)
    Path(file.name).replace(path)


def public_case(root, meta, **arrays):
    write_npz(root / "public" / f"{meta['case']}.npz",
              case_json=np.array(json.dumps(meta, ensure_ascii=False)), **arrays)


def synthetic_name(cohort, object_id, array):
    prefix = "dev__" if cohort == "development" else ""
    return f"{prefix}{object_id}__{array.lower()}"


def present(root, cases, name, scoring_folder):
    return (name in cases and (root / "public" / f"{name}.npz").is_file() and
            (root / scoring_folder / f"{name}.npz").is_file())


def synthetic(root, cases, receipt, ops, cohort, object_id, index, chi, family=None):
    wall, cpu = time.perf_counter(), time.process_time()
    clean, cost = history.fine_solve(chi, ops.fine_e, ops.fine_s,
                                    ops.green.apply, ops.green.diagonal)
    coarse = area_average(chi, ops.fine, SIDE)
    receipt["forwards"][f"{cohort}/{object_id}"] = dict(
        cost, family=family, fine_side=ops.fine,
        wall_s=time.perf_counter() - wall, cpu_s=time.process_time() - cpu,
    )
    for array_index, array in enumerate(ops.arrays):
        case = synthetic_name(cohort, object_id, array.name)
        tx, rx = ops.tx_ids[array.name], ops.rx_ids[array.name]
        selected = clean[np.ix_(rx, tx)]
        if cohort == "development" and array.name == "PW120_7x14":
            seed = 20261121000 + index  # Exact previous run seed.
        else:
            seed = np.random.SeedSequence([NOISE_SEED, index, array_index,
                                           int(cohort == "development")])
        noise = history.complex_noise(selected, NOISE_FRACTION, seed)
        meta = dict(case=case, suite="LA_SYNTH", side=SIDE, width=WIDTH,
                    object_id=object_id, array=array.name, cohort=cohort,
                    frequency=FREQUENCY, fine=ops.fine, noise_relative=NOISE_FRACTION,
                    train_channels=selected.size, tx_count=len(tx), rx_count=len(rx),
                    tx_kind=array.tx_kind, nominal_aperture_deg=array.aperture,
                    tx_plane_wave_propagation_deg=None, tx_line_source_m=None,
                    rx_m=array.receiver_points.tolist(),
                    source_normalization="unit_mean_incident_power_on_fixed_domain")
        if array.tx_kind == "plane_wave":
            meta["tx_plane_wave_propagation_deg"] = array.tx[:, 0].tolist()
        else:
            meta["tx_line_source_m"] = array.tx.tolist()
            meta["tx_line_source_gain"] = ops.source_gain[tx].tolist()
        score = dict(m=coarse, object_id=np.array(object_id))
        if family is not None:
            score["family"] = np.array(family)
        write_npz(root / "truth" / f"{case}.npz", **score)
        public_case(root, meta, E=ops.coarse_e[:, tx], S=ops.coarse_s[rx],
                    Y=selected + noise, centers=ops.coarse_centers, k=ops.k,
                    radius=ops.coarse_radius, delta=np.linalg.norm(noise))
        cases[case] = meta

    log_update(root, cases, receipt)
    print(f"Prepared {object_id}: {len(ops.arrays)} arrays, "
          f"{cost['forward_rhs']} RHS, {receipt['forwards'][f'{cohort}/{object_id}']['wall_s']:.2f} s",
          flush=True)


def log_update(root, cases, receipt):
    write_json(root / "generation_cost.json", receipt)
    write_json(root / "cases.json", sorted(cases.values(), key=lambda x: x["case"]))


def run_synthetic(root, cases, receipt, cohort):
    if cohort == "development":
        by_case = {case["case"]: case for case in history.all_specs()}
        for object_id in DEV_CASES:
            names = DEV_ARRAYS
            if all(present(root, cases, synthetic_name(cohort, object_id, a), "truth") for a in names):
                continue
            spec = by_case[object_id]
            ops = SyntheticOperators(spec["fine"], names)
            receipt["setups"][f"development/{object_id}"] = ops.cost
            chi = history.material(spec, ops.fine_centers)
            synthetic(root, cases, receipt, ops, cohort, object_id, spec["index"], chi)
        return

    names = CONFIRM_ARRAYS
    permutation = np.random.default_rng(FAMILY_SEED).permutation(np.repeat(np.arange(6), 4))
    missing = [
        (index, family, f"s2d_{index:03d}")
        for index, family in enumerate(permutation, 1)
        if not all(present(root, cases, synthetic_name(cohort, f"s2d_{index:03d}", a), "truth")
                   for a in names)
    ]
    if missing:
        ops = SyntheticOperators(96, names)
        receipt["setups"].setdefault("confirmation/96", []).append(ops.cost)
        for index, family, object_id in missing:
            chi = random_material(ops.fine_centers, index, family)
            synthetic(root, cases, receipt, ops, cohort, object_id, index, chi,
                      FAMILY_NAMES[family])


def fresnel_selection(y):
    contiguous = np.arange(9)
    interleaved = np.arange(0, 18, 2)
    receiver_ids = np.arange(y.shape[0])
    finite = np.isfinite(y)
    upper = receiver_ids <= 180
    even = receiver_ids % 2 == 0
    continuous_rx = [np.flatnonzero(finite[:, source] & upper) for source in contiguous]
    interleaved_rx = []
    for source, equal in zip(interleaved, continuous_rx, strict=True):
        candidates = np.flatnonzero(finite[:, source] & even)
        # Use the same physical angular ordering for every original source.
        selected = np.round(np.linspace(0, len(candidates) - 1, len(equal))).astype(int)
        if len(np.unique(selected)) != len(equal):
            raise RuntimeError("Interleaved mask has too few real channels")
        interleaved_rx.append(candidates[selected])

    joint_hold_tx = np.setdiff1d(np.arange(18), np.union1d(contiguous, interleaved))
    joint_seen_tx = np.intersect1d(contiguous, interleaved)
    held_rx = y[:, joint_seen_tx].copy()
    held_rx[upper | even] = np.nan + 1j * np.nan
    return ((contiguous, continuous_rx), (interleaved, interleaved_rx),
            joint_hold_tx, joint_seen_tx, held_rx)


def run_fresnel(root, cases, receipt):
    variants = ("FRESNEL_CONTIGUOUS", "FRESNEL_INTERLEAVED")
    names = tuple(f"fresnel_twin_8GHz__{array.lower().removeprefix('fresnel_')}"
                  for array in variants)
    if all(present(root, cases, name, "heldout") for name in names):
        return
    wall, cpu = time.perf_counter(), time.process_time()
    frequency, side, width = 8e9, 60, 0.2
    k = 2 * np.pi * frequency / history.C0
    path = FRESNEL_2005_DIR / "FoamTwinDielTM.exp"
    _, y, inc = load_file(path, frequency / 1e9)
    sources = 1.67 * unit_from_degrees(np.arange(y.shape[1]) * 360 / y.shape[1])
    receivers = 1.67 * unit_from_degrees(np.arange(y.shape[0]))
    centers = history.grid(side, width)
    ops = operators(centers, sources, receivers, width / side, k, dense=False)
    e, calibration, _ = incident_fit(inc, sources, receivers, centers, k)
    continuous, interleaved, held_tx, common_tx, y_held_rx = fresnel_selection(y)
    held = dict(E=e[:, held_tx], S=ops["S"], Y=y[:, held_tx],
                E_train_hemi=e[:, common_tx], Y_train_opposite_rx=y_held_rx,
                source_ids=held_tx, source_ids_train_hemi=common_tx)
    for name, variant, (tx, receivers_by_source) in zip(
            names, variants, (continuous, interleaved), strict=True):
        n_columns = len(tx)
        n_rows = max(map(len, receivers_by_source))
        valid = np.zeros((n_rows, n_columns), dtype=bool)
        rx_ids = np.full((n_columns, n_rows), -1, dtype=int)
        y_train = np.zeros((n_rows, n_columns), dtype=complex)
        for j, (source_id, rx) in enumerate(zip(tx, receivers_by_source, strict=True)):
            rx_ids[j, :len(rx)] = rx
            valid[:len(rx), j] = True
            y_train[:len(rx), j] = y[rx, source_id]
            if source_id in common_tx and np.isfinite(y_held_rx[
                    rx, np.searchsorted(common_tx, source_id)]).any():
                raise RuntimeError("Shared held receive side overlaps training")
        if np.intersect1d(held_tx, tx).size or not np.all(np.isfinite(y_train[valid])):
            raise RuntimeError("Fresnel held Tx or finite scatter channels are inconsistent")
        meta = dict(
            case=name, suite="LA_FRESNEL", side=side, width=width,
            object_id="FoamTwinDielTM", array=variant, cohort="measured_extra",
            frequency=frequency, fine=0, source_specific=True, ragged_rx=True,
            train_sources=(tx + 1).tolist(), rx_radius_m=1.67, tx_radius_m=1.67,
            train_channels=int(valid.sum()), held_channels=int(np.isfinite(held["Y"]).sum()),
            held_rx_channels=int(np.isfinite(y_held_rx).sum()), data_ball_ratio=0.15,
            calibration_median_validation=float(np.median([
                r["validation_relative"] for r in calibration])),
        )
        write_npz(root / "heldout" / f"{name}.npz", **held)
        public_case(root, meta, E=e[:, tx], S=ops["S"], Y=y_train, centers=centers,
                    k=k, radius=ops["radius"], delta=0.15 * np.linalg.norm(y_train[valid]),
                    rx_ids=rx_ids, valid_mask=valid)
        cases[name] = meta
    receipt["measured"]["Fresnel_Twin_8GHz"] = {
        "wall_s": time.perf_counter() - wall, "cpu_s": time.process_time() - cpu,
        "path": str(path.relative_to(WORKSPACE)) if path.is_relative_to(WORKSPACE) else str(path),
        "empty_target_source_fits": len(calibration), "fine_physics_RHS": 0,
        "train_channels_each": cases[names[0]]["train_channels"],
        "held_source_channels_common": cases[names[0]]["held_channels"],
        "held_receive_channels_common": cases[names[0]]["held_rx_channels"],
    }
    log_update(root, cases, receipt)
    print(f"Prepared Fresnel pair: {cases[names[0]]['train_channels']} train, "
          f"{cases[names[0]]['held_channels']} shared held Tx, "
          f"{cases[names[0]]['held_rx_channels']} shared held Rx channels", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "data")
    parser.add_argument("--cohort", choices=("development", "confirmation", "fresnel", "all"),
                        default="all")
    args = parser.parse_args()
    root = args.output
    for folder in ("public", "truth", "heldout"):
        (root / folder).mkdir(parents=True, exist_ok=True)
    manifest = root / "cases.json"
    costs = root / "generation_cost.json"
    cases = {entry["case"]: entry for entry in json.loads(manifest.read_text())} if manifest.exists() else {}
    receipt = json.loads(costs.read_text()) if costs.exists() else {
        "platform": platform.platform(), "numpy": np.__version__, "threads_per_process": 1,
        "setups": {}, "forwards": {}, "measured": {},
    }
    if args.cohort in ("development", "all"):
        run_synthetic(root, cases, receipt, "development")
    if args.cohort in ("fresnel", "all"):
        run_fresnel(root, cases, receipt)
    if args.cohort in ("confirmation", "all"):
        run_synthetic(root, cases, receipt, "confirmation")


if __name__ == "__main__":
    main()
