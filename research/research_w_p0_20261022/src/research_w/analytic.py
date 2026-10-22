"""Independent 1D planar Maxwell regression for Research W P0-A.

Lengths are in lambda0, frequencies in f0, fields use exp(-i omega t).
The implementation uses traveling electric-field amplitudes and interface
continuity, independently of the source package's characteristic matrices.
Both exterior media have n=1. This module makes no 2D imaging claim.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

WALL1 = ((2.0, 0.25),)
WALL2 = ((3.0, 1.0 / 6.0),)
BANDWIDTHS = (0.0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)


def _inputs(layers, frequencies):
    layers = tuple((float(n), float(d)) for n, d in layers)
    f = np.atleast_1d(np.asarray(frequencies, dtype=float))
    if f.ndim != 1 or not np.all(np.isfinite(f)) or np.any(f <= 0):
        raise ValueError("frequencies must be a finite positive 1D array")
    if any(not math.isfinite(n + d) or n <= 0 or d < 0 for n, d in layers):
        raise ValueError("lossless layers require finite n>0 and d>=0")
    return layers, f


def _interface(n_left, n_right):
    """Map [forward, backward] amplitudes from left to right at an interface."""
    ratio = n_left / n_right
    return 0.5 * np.array([[1 + ratio, 1 - ratio],
                           [1 - ratio, 1 + ratio]], dtype=complex)


def amplitude_transfer(layers, frequencies):
    """Map amplitudes at the left physical face to the right physical face."""
    layers, f = _inputs(layers, frequencies)
    total = np.broadcast_to(np.eye(2, dtype=complex), (len(f), 2, 2)).copy()
    n_previous = 1.0
    for n, thickness in layers:
        total = _interface(n_previous, n) @ total
        phase = 2 * np.pi * f * n * thickness
        total[:, 0, :] *= np.exp(1j * phase)[:, None]
        total[:, 1, :] *= np.exp(-1j * phase)[:, None]
        n_previous = n
    return _interface(n_previous, 1.0) @ total


def scattering(layers, frequencies):
    """Return (r_left, r_right, t), referenced at the two physical faces."""
    transfer = amplitude_transfer(layers, frequencies)
    denominator = transfer[:, 1, 1]
    return (-transfer[:, 1, 0] / denominator,
            transfer[:, 0, 1] / denominator, 1 / denominator)


def fresnel_slab(n, thickness, frequencies):
    """Closed-form interface multiple-reflection sum for one slab."""
    _, f = _inputs(((n, thickness),), frequencies)
    rho = (1 - n) / (1 + n)
    propagation = np.exp(2j * np.pi * f * n * thickness)
    denominator = 1 - rho ** 2 * propagation ** 2
    reflection = rho * (1 - propagation ** 2) / denominator
    transmission = (1 - rho ** 2) * propagation / denominator
    return reflection, transmission


def reflection_recursion(layers, frequencies):
    """Independent reflection recursion, including every interface return."""
    layers, f = _inputs(layers, frequencies)
    reflection = np.zeros(len(f), dtype=complex)
    n_right = 1.0
    for index in range(len(layers) - 1, -1, -1):
        n, thickness = layers[index]
        r_exit = (n - n_right) / (n + n_right)
        reflection = (r_exit + reflection) / (1 + r_exit * reflection)
        reflection *= np.exp(4j * np.pi * f * n * thickness)
        n_right = n
    r_entrance = (1 - n_right) / (1 + n_right)
    return (r_entrance + reflection) / (1 + r_entrance * reflection)


def fixed_plane_transmission(layers, frequencies, plane=1.0):
    layers, f = _inputs(layers, frequencies)
    thickness = sum(d for _, d in layers)
    if plane < thickness:
        raise ValueError("fixed reference plane must lie behind the wall")
    return scattering(layers, f)[2] * np.exp(2j * np.pi * f * (plane - thickness))


def full_target_reflection(layers, frequencies, position, rho=-1.0):
    """Complete wall/load cascade: direct wall return plus all round trips."""
    layers, f = _inputs(layers, frequencies)
    thickness = sum(d for _, d in layers)
    if position < thickness:
        raise ValueError("target must lie behind the wall")
    r_left, r_right, transmission = scattering(layers, f)
    load = rho * np.exp(4j * np.pi * f * (position - thickness))
    return r_left + transmission ** 2 * load / (1 - r_right * load)


def pec_reflection_direct(layers, frequencies, position):
    """Solve E(position)=0 directly, independently of the cascade formula."""
    layers, f = _inputs(layers, frequencies)
    gap = position - sum(d for _, d in layers)
    if gap < 0:
        raise ValueError("PEC must lie behind the wall")
    transfer = amplitude_transfer((*layers, (1.0, gap)), f)
    # a_PEC+b_PEC=0, with a_left=1. This includes the full multilayer field.
    return -(transfer[:, 0, 0] + transfer[:, 1, 0]) / (
        transfer[:, 0, 1] + transfer[:, 1, 1])


def realify(data):
    data = np.asarray(data, dtype=complex).ravel()
    return np.concatenate((data.real, data.imag))


def whitened_distance(a, b, sigma_each_real_quadrature=0.01):
    sigma = float(sigma_each_real_quadrature)
    if sigma <= 0 or not math.isfinite(sigma):
        raise ValueError("quadrature noise standard deviation must be positive")
    return float(np.linalg.norm(realify(np.asarray(a) - np.asarray(b))) / sigma)


def _defect(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b))))


def run_regression(output_dir, source_results=None):
    """Write deterministic numerical evidence and CPU timing, without sampling."""
    start_wall, start_cpu = time.perf_counter(), time.process_time()
    started = datetime.now(timezone.utc).isoformat()
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    sigma, nfreq, plane = 0.01, 9, 1.0
    target1, target2 = 1.5, 1.5 - 1.0 / 12.0
    center_r1, _, center_t1 = scattering(WALL1, [1.0])
    center_r2, _, center_t2 = scattering(WALL2, [1.0])
    h1 = fixed_plane_transmission(WALL1, [1.0], plane) ** 2
    h2 = fixed_plane_transmission(WALL2, [1.0], plane) ** 2
    y1 = full_target_reflection(WALL1, [1.0], target1)
    y2 = full_target_reflection(WALL2, [1.0], target2)
    tests = {
        "max_center_reflection_abs": max(_defect(center_r1, 0), _defect(center_r2, 0)),
        "center_task_response_difference_abs": _defect(h1, h2),
        "center_task_phase_difference_rad": float(np.angle(h2 / h1)[0]),
        "joint_center_data_difference_abs": _defect(y1, y2),
        "joint_target_position_difference_wavelength": abs(target1 - target2),
        "center_slab_transmission_minus_one_defect": max(_defect(center_t1, -1), _defect(center_t2, -1)),
    }
    dense = np.linspace(0.8, 1.2, 1001)
    configs = (WALL1, WALL2, ((1.0, 0.1), (2.3, 0.12), (1.4, 0.08)))
    unitarity, fresnel, recursion, direct_pec, reciprocity = [], [], [], [], []
    arrays = {"dense_frequencies": dense, "center_h1": h1, "center_h2": h2,
              "center_y1": y1, "center_y2": y2}
    for index, layers in enumerate(configs):
        rl, rr, t = scattering(layers, dense)
        matrix = np.stack((np.stack((rl, t), axis=-1),
                           np.stack((t, rr), axis=-1)), axis=-2)
        unitarity.append(_defect(matrix.conj().swapaxes(-1, -2) @ matrix, np.eye(2)))
        reciprocity.append(_defect(np.linalg.det(amplitude_transfer(layers, dense)), 1))
        recursion.append(_defect(rl, reflection_recursion(layers, dense)))
        direct = pec_reflection_direct(layers, dense, target1)
        cascade = full_target_reflection(layers, dense, target1)
        direct_pec.append(_defect(direct, cascade))
        arrays.update({f"config{index}_r_left": rl, f"config{index}_r_right": rr,
                       f"config{index}_transmission": t, f"config{index}_pec_direct": direct,
                       f"config{index}_pec_cascade": cascade})
        if len(layers) == 1:
            fr, ft = fresnel_slab(*layers[0], dense)
            fresnel.extend((_defect(fr, rl), _defect(ft, t)))
    tests.update({"maximum_lossless_unitarity_defect": max(unitarity),
                  "maximum_independent_fresnel_difference": max(fresnel),
                  "maximum_interface_recursion_difference": max(recursion),
                  "maximum_direct_pec_cascade_difference": max(direct_pec),
                  "maximum_transfer_determinant_defect": max(reciprocity)})
    rows = []
    for index, bandwidth in enumerate(BANDWIDTHS):
        f = np.linspace(1 - bandwidth / 2, 1 + bandwidth / 2, nfreq)
        r1, r2 = scattering(WALL1, f)[0], scattering(WALL2, f)[0]
        joint1 = full_target_reflection(WALL1, f, target1)
        joint2 = full_target_reflection(WALL2, f, target2)
        task1, task2 = (fixed_plane_transmission(wall, f, plane) ** 2 for wall in (WALL1, WALL2))
        distance = whitened_distance(r1, r2, sigma)
        rows.append({"fractional_bandwidth": bandwidth, "wall_distance": distance,
                     "wall_optimal_equal_prior_test_error": 0.5 * math.erfc(distance / (2 * math.sqrt(2))),
                     "joint_shifted_target_distance": whitened_distance(joint1, joint2, sigma),
                     "task_response_rms_difference": float(np.linalg.norm(task1 - task2) / math.sqrt(nfreq))})
        arrays.update({f"band{index}_frequencies": f, f"band{index}_wall1": r1,
                       f"band{index}_wall2": r2, f"band{index}_joint1": joint1,
                       f"band{index}_joint2": joint2, f"band{index}_task1": task1,
                       f"band{index}_task2": task2})
    payload = {"status": "executed_independent_analytic_1d_regression_only",
               "method": "traveling-wave interface matrices; Fresnel sum; direct PEC boundary",
               "reference_planes": {"common_front": 0.0, "fixed_transmission_plane": plane,
                                    "wall_scattering": "physical faces"},
               "noise_sigma_each_real_quadrature": sigma, "number_of_frequencies": nfreq,
               "wall1_n_d": WALL1, "wall2_n_d": WALL2,
               "target1_position": target1, "target2_position": target2,
               "tests": tests, "bandwidth_sweep": rows,
               "limitations": ["Infinite planar lossless nonmagnetic layers at normal incidence.",
                               "No 2D/3D fullwave, inverse imaging, posterior mass, or speedup measured.",
                               "Bandwidth grids are not nested; zero bandwidth means repeated same-frequency samples."]}
    if source_results is not None:
        source = json.loads(Path(source_results).read_text())
        differences = {key: abs(tests[key] - value) for key, value in source["tests"].items()}
        sweep = [{key: abs(row[key] - reference[key]) for key in row}
                 for row, reference in zip(rows, source["bandwidth_sweep"], strict=True)]
        payload["archived_package_comparison"] = {"source": str(Path(source_results).resolve()),
             "test_scalar_absolute_differences": differences, "sweep_absolute_differences": sweep,
             "maximum_sweep_absolute_difference": max(value for row in sweep for value in row.values())}
    np.savez_compressed(out / "arrays.npz", **arrays)
    (out / "results.json").write_text(json.dumps(payload, indent=2) + "\n")
    receipt = {"started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
               "wall_seconds": time.perf_counter() - start_wall,
               "process_cpu_seconds": time.process_time() - start_cpu,
               "python_executable": sys.executable, "numpy_version": np.__version__,
               "platform": platform.platform(), "device": "CPU",
               "timing_scope": "run_regression numerical checks and artifact writes; excludes Python imports",
               "threads": {key: os.environ.get(key) for key in
                           ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
               "new_2d_or_3d_fullwave_calls": 0, "material_search_calls": 0}
    (out / "cpu_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-results", type=Path)
    args = parser.parse_args()
    result = run_regression(args.output, args.source_results)
    print(json.dumps({"output": str(args.output.resolve()), "tests": result["tests"],
                      "narrow_band": result["bandwidth_sweep"][2]}, indent=2))


if __name__ == "__main__":
    main()
