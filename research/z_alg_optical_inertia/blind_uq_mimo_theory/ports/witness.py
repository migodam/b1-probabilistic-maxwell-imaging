"""Nonlinear two-material resolution mechanism, separate from the blind cohort."""
import argparse
import csv
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import lu_factor, lu_solve, solve

from model import BANDS, COUNTS, STRENGTHS, DiskGreen, network_terms
from model import operator, physical_solve, port_point_matrices
from sweep import ERROR_RELATIVE, change_network, clean_transfer, normalized_reactance_error, tag, union_positions

HERE = Path(__file__).resolve().parent
WEAK_COUNT, BASE_BAND, BASE_STRENGTH = 6, 400, 48.0


def planted_material(disk):
    x, y = disk.centers.T
    body = np.hypot((x+0.34)/0.69, (y-0.08)/0.46)
    dent = np.hypot((x+0.17)/0.26, (y-0.30)/0.17)
    component = ((x-0.58)/0.23)**2 + ((y+0.54)/0.24)**2 < 1
    main = (body < 1) & (dent > 0.9)
    region = ((x-0.49)/0.14)**2 + ((y-0.47)/0.13)**2 < 1
    m = np.where(main, 3.35 + 0.38*np.sin(3*x+6*y), 0)
    m = np.maximum(m, np.where(component, 2.45 + 0.3*np.cos(6*x), 0))
    return m, region


class RestrictedInverse:
    def __init__(self, data, active):
        self.active = active
        self.G = data["G"][np.ix_(active, active)]
        self.E = data["E0"][active]
        self.S = data["S"][:, active]
        self.lu = 0
        self.rhs = 0

    def __call__(self, m, grad=False):
        M = np.eye(len(m)) - m[:, None] * self.G
        LU = lu_factor(M, overwrite_a=True, check_finite=False)
        J = lu_solve(LU, m[:, None]*self.E, check_finite=False)
        self.lu += 1
        self.rhs += self.E.shape[1]
        F = self.S @ J
        if not grad:
            return F
        receiver = lu_solve(LU, self.S.T, trans=1, check_finite=False).T
        self.rhs += self.S.shape[0]
        field = self.E + self.G @ J
        derivative = (receiver[:, None, :] * field.T[None, :, :]).reshape(F.size, len(m))
        return F, np.concatenate((derivative.real, derivative.imag), axis=0)


def fit_no_patch(disk, m0, region, strength):
    active = np.flatnonzero(m0)
    net = network_terms(disk, WEAK_COUNT, BASE_STRENGTH)
    data = operator(disk, net, disk.dense())
    no_patch = RestrictedInverse(data, active)
    with_patch_model = RestrictedInverse(data, np.flatnonzero((m0 > 0) | region))
    with_patch = m0 + strength * region
    target = with_patch_model(with_patch[with_patch_model.active])
    x = m0[active].copy()
    trace = []
    level = np.linalg.norm(target)
    for _ in range(18):
        value, jac = no_patch(x, grad=True)
        diff = (target-value).ravel()
        loss = np.linalg.norm(diff)
        trace.append(float(loss/level))
        if loss < 1e-5*level:
            break
        rhs = np.r_[diff.real, diff.imag]
        gram = jac @ jac.T
        ridge = (1e-7+min(1e-4, loss*0.001/level)) * np.max(np.diag(gram))
        step = jac.T @ solve(gram + ridge*np.eye(len(gram)), rhs, assume_a="pos")
        improved = False
        for scale in (1.0, 0.5, 0.25, 0.125, 0.0625):
            trial = np.clip(x + scale*step, 0, 4.5)
            if np.linalg.norm(no_patch(trial) - target) < loss:
                x = trial
                improved = True
                break
        if not improved:
            break
    absent = m0.copy()
    absent[active] = x
    return with_patch, absent, trace, dict(
        LU=no_patch.lu+with_patch_model.lu, complex_rhs=no_patch.rhs+with_patch_model.rhs,
        compensated_cells=len(active), roi_cells=int(region.sum()),
        roi_material_difference=strength,
        pair_nrmse=float(np.linalg.norm(with_patch-absent)/np.linalg.norm(with_patch)),
        compensated_material_difference_RMS=float(np.linalg.norm(absent-m0)/np.sqrt(len(active))))


def pair_grid(material0, material1, sigma):
    positions, geometry_indices = union_positions()
    rows, arrays = [], {}
    for band in BANDS:
        disk = DiskGreen(32, band)
        dense = disk.dense()
        src, rcv, _ = port_point_matrices(disk, positions)
        fft = lambda a: disk.kernel.apply(a)
        j0, cost0 = physical_solve(material0, fft, src, disk.kernel.diagonal, dense)
        j1, cost1 = physical_solve(material1, fft, src, disk.kernel.diagonal, dense)
        common = [rcv @ j0, rcv @ j1]
        for q in COUNTS:
            subports = geometry_indices[q]
            submat = [resp[np.ix_(subports, subports)] for resp in common]
            for strength in STRENGTHS:
                net = network_terms(disk, q, strength)
                direction = normalized_reactance_error(net["positions"])
                nominal_a, nominal_b = [clean_transfer(v, net) / np.sqrt(q) for v in submat]
                for level in ERROR_RELATIVE:
                    actual = change_network(net, level*direction)
                    true_a, true_b = [clean_transfer(v, actual) / np.sqrt(q) for v in submat]
                    select = np.ix_(net["rx"], np.arange(q))
                    y_a, y_b = true_a[select], true_b[select]
                    distance = np.linalg.norm(y_a-y_b)
                    nominal_distance = np.linalg.norm(nominal_a[select]-nominal_b[select])
                    noise_radius = 1.5 * sigma * q
                    # Benchmark using all-terminal reference transfers supplied
                    # as known; do not infer this penalty from masked data.
                    beta = (2*level+level**2)*max(np.linalg.norm(true_a), np.linalg.norm(true_b))
                    noisy_midpoint = 0.5*(y_a+y_b)
                    key = tag(band, q, strength, level)
                    arrays[key+"_midpoint_Y"] = noisy_midpoint
                    arrays[key+"_F_with_region"] = y_a
                    arrays[key+"_F_without_region"] = y_b
                    held = np.ix_(net["hold"], np.arange(q))
                    arrays[key+"_F_with_region_hold"] = true_a[held]
                    arrays[key+"_F_without_region_hold"] = true_b[held]
                    rows.append(dict(frequency_MHz=band, q=q, radiation_scale_ohms=strength,
                                     calibration_eta=level, readout_frobenius_separation=distance,
                                     nominal_readout_frobenius_separation=nominal_distance,
                                     noise_radius=noise_radius, calibration_radius_full=beta,
                                     actual_pair_sep_over_two_r=float(distance/(2*noise_radius)),
                                     nominal_pair_sep_over_two_inflated_r=float(nominal_distance/(2*(noise_radius+beta))),
                                     both_present_in_shared_ball_numeric=bool(distance < 2*noise_radius),
                                     radius_kind="fixed absolute stochastic-free hard radius; instrument rms from 8x8 mechanism object",
                                     solve_cost_left=cost0["seconds_wall"], solve_cost_right=cost1["seconds_wall"]))
    return rows, arrays


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--strength", type=float, default=0.75)
    args = p.parse_args()
    started = perf_counter(), process_time()
    disk = DiskGreen(32, BASE_BAND)
    base, region = planted_material(disk)
    with_patch, absent, fit_trace, cost = fit_no_patch(disk, base, region, args.strength)
    net = network_terms(disk, 8, BASE_STRENGTH)
    nominal = operator(disk, net, disk.dense())
    forward = RestrictedInverse(nominal, np.flatnonzero(with_patch))
    sigma = 0.02 * np.linalg.norm(forward(with_patch[forward.active])) / 8
    rows, arrays = pair_grid(with_patch, absent, sigma)
    np.savez_compressed(HERE/"results"/"mechanism_witness.npz", with_region=with_patch.reshape(32, 32),
                        without_region=absent.reshape(32, 32), roi=region.reshape(32, 32), reference_sigma=sigma,
                        **arrays)
    (HERE/"results"/"mechanism_witness.json").write_text(json.dumps(dict(
        cohort="separate designed resolution mechanism; no blind object used",
        construction=f"optimized at {WEAK_COUNT}x{WEAK_COUNT}, {BASE_BAND} MHz",
        patch_amplitude=args.strength, trace_relative=fit_trace, cost=cost,
        calibration="same additive reciprocal perturbation for two materials",
        equal_available_power_W=1, single_band_rows=len(rows),
        seconds_wall=perf_counter()-started[0], seconds_CPU=process_time()-started[1]), indent=2))
    with (HERE/"results"/"mechanism_witness.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        if row["calibration_eta"] in (0, 0.005, 0.025) and row["radiation_scale_ohms"] == BASE_STRENGTH:
            print(row["frequency_MHz"], row["q"], row["calibration_eta"],
                  row["actual_pair_sep_over_two_r"], row["nominal_pair_sep_over_two_inflated_r"])


if __name__ == "__main__":
    main()
