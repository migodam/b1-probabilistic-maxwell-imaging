"""Paired fixed-aperture array and additive calibration sensitivity.

Generation reads sealed labels. The optimization agents get only the released
complex measurements and nominal models, not this script's private transfers.
"""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import solve

from model import BANDS, COUNTS, R, STRENGTHS, DiskGreen, antenna_forward, aperture, network_terms
from model import operator, physical_solve, port_point_matrices, power_report

HERE = Path(__file__).resolve().parent
ERROR_RELATIVE = (0.0, 0.005, 0.025, 0.06)


def normalized_reactance_error(positions):
    """Fixed reciprocal parasitic reactance. The spectral norm is exactly one."""
    distance = np.linalg.norm(positions[:, None] - positions[None, :], axis=-1)
    angles = np.arctan2(positions[:, 1], positions[:, 0])
    real_symmetric = np.cos(3*(angles[:, None]+angles[None, :])) * np.exp(-distance/0.75)
    np.fill_diagonal(real_symmetric, 0)
    return 1j * real_symmetric / np.linalg.norm(real_symmetric, ord=2)


def change_network(network, normalized_error):
    shift = R * normalized_error
    out = dict(network)
    out["matrix"] = network["matrix"] + shift
    out["antenna"] = network["antenna"] + shift
    out["inv"] = np.linalg.solve(out["matrix"], np.eye(len(shift)))
    return out


def union_positions():
    geometry = np.concatenate([aperture(count) for count in COUNTS])
    canonical = geometry.round(11)
    _, index, table = np.unique(canonical, axis=0, return_index=True, return_inverse=True)
    positions = geometry[index]
    blocks = np.split(table, np.cumsum([2*q+8 for q in COUNTS])[:-1])
    return positions, dict(zip(COUNTS, blocks, strict=True))


def clean_transfer(bare_response, network):
    g = network["gain"]
    loaded = network["matrix"] - 1j * g[:, None] * bare_response * g[None, :]
    # Dimensionless, all-terminal, 2R times difference of the impedance inverses.
    return 2 * R * (network["inv"] - np.linalg.solve(loaded, np.eye(len(g))))


def noise_direction(case_id, band, count, hold):
    token = 867992 + 1001*sum((i+1)*ord(v) for i, v in enumerate(case_id))
    rng = np.random.default_rng(token + 500*band + count + 914503*hold)
    shape = (8 if hold else count, count)
    draw = rng.normal(size=shape) + 1j*rng.normal(size=shape)
    return draw / np.linalg.norm(draw)


def off_diagonal_impedance(antenna):
    return antenna - np.diag(np.diag(antenna))


def response_row(disk, case_id, count, coupling, error, bare_response, network, sigma_ref, sigma_hold):
    start = perf_counter()
    err = error * normalized_reactance_error(network["positions"])
    actual_net = change_network(network, err)
    nominal_transfer = clean_transfer(bare_response, network)
    actual_transfer = clean_transfer(bare_response, actual_net)
    f_nominal, f_actual = nominal_transfer/np.sqrt(count), actual_transfer/np.sqrt(count)
    rx, hold, tx = network["rx"], network["hold"], np.arange(count)
    clean = f_actual[np.ix_(rx, tx)]
    y_nominal = f_nominal[np.ix_(rx, tx)]
    clean_hold = f_actual[np.ix_(hold, tx)]
    n = sigma_ref * count * noise_direction(case_id, disk.band, count, False)
    n_hold = sigma_hold * np.sqrt(8*count) * noise_direction(case_id, disk.band, count, True)
    raw_y, currents, z_loaded = antenna_forward(bare_response, actual_net)
    power = power_report(currents, z_loaded, actual_net["drive"])
    norm = np.linalg.norm
    a_empty = network["matrix"] / R
    b_loaded = z_loaded / R + np.eye(len(err)) - err
    factored_nominal = ((np.eye(len(err)) + solve(b_loaded, err)) @ actual_transfer @
                        (np.eye(len(err)) + solve(a_empty.T, err.T).T))
    eps_rel = np.linalg.norm(err, ord=2)
    # Full-transfer benchmark for the theorem. A visible block has no access
    # to the norms unless separately calibrated/measured.
    full_bound = (2*eps_rel+eps_rel**2) * norm(actual_transfer) / np.sqrt(count)
    absolute_block_bound = 4 * eps_rel
    # The adjoint-compatible discretization retains a true all-source 8-16 s cap.
    s_nom = np.linalg.svd(network["inv"][rx] @ network["port_readout"], compute_uv=False)
    b_nom = np.linalg.svd(network["port_excitation"] @ (network["inv"] @ network["drive"]), compute_uv=False)
    absolute_difference = float(norm(clean-y_nominal))
    if absolute_difference > min(full_bound, absolute_block_bound) + 1e-12:
        raise ArithmeticError(f"Calibration containment failed for {case_id}, {count}, {coupling}, {error}")
    if norm(raw_y[rx] - clean) > 1e-12:
        raise ArithmeticError(f"Physical amplitude normalization mismatch on {case_id}")
    if power["energy_balance_max_W"] > 1e-12 or power["impedance_min_hermitian_eigen_ohms"] < 0:
        raise ArithmeticError(f"Non-passive loaded network on {case_id}: {power}")
    return dict(
        case_id=case_id, side=disk.shape[0], frequency_MHz=disk.band, q=count, radiation_scale_ohms=coupling,
        normalized_calibration_error=eps_rel, calibration_error_ohms=R*eps_rel,
        off_diagonal_Z_operator_ohms=float(norm(off_diagonal_impedance(network["antenna"]), ord=2)),
        loaded_total_condition=float(np.linalg.cond(z_loaded + R*np.eye(len(z_loaded)))),
        loaded_resistance_min_ohm=power["impedance_min_hermitian_eigen_ohms"],
        loaded_reciprocity_relative=power["reciprocity_relative"],
        accepted_total_W=sum(power["driven_accepted_W"]),
        reflected_total_W=sum(power["reflected_W"]),
        antenna_total_W=sum(power["antenna_W"]), other_termination_total_W=sum(power["other_termination_W"]),
        energy_balance_max_W=power["energy_balance_max_W"],
        transfer_identity_max_abs=float(np.max(abs(factored_nominal - nominal_transfer))),
        radius_noise=float(1.5*norm(n)), radius_hold_noise=float(1.5*norm(n_hold)),
        Y_train_calibration_difference=absolute_difference,
        Y_train_clean_norm=float(norm(clean)), Y_hold_clean_norm=float(norm(clean_hold)),
        Y_train_nominal_minus_observed=float(norm(y_nominal-clean-n)),
        Y_train_noise_norm=float(norm(n)), Y_train_calibration_bound_full=float(full_bound),
        Y_train_calibration_bound_absolute=float(absolute_block_bound),
        full_transfer_Fro=float(norm(actual_transfer)), full_transfer_spectral=float(norm(actual_transfer, ord=2)),
        receiver_numeric_rank=int(np.count_nonzero(s_nom / s_nom[0] > 1e-10)),
        illumination_numeric_rank=int(np.count_nonzero(b_nom / b_nom[0] > 1e-10)),
        receiver_10pct_modes=int(np.count_nonzero(s_nom / s_nom[0] > 0.1)),
        illumination_10pct_modes=int(np.count_nonzero(b_nom / b_nom[0] > 0.1)),
        seconds_wall=perf_counter()-start,
    ), clean+n, clean_hold+n_hold, nominal_transfer, actual_transfer, power


def tag(band, q, coupling, err):
    return f"f{band}_a{q:02}_c{int(coupling):03}_e{int(round(1000*err)):03}"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--side", type=int, required=True, choices=[32, 96])
    p.add_argument("--items", type=int, required=True, nargs="+")
    args = p.parse_args()
    side = args.side
    start_all = perf_counter(), process_time()
    positions, geometry_indices = union_positions()
    rows, costs = [], {}
    all_data = defaultdict(dict)
    private_transfer = defaultdict(dict)
    port_powers = defaultdict(dict)
    for band in BANDS:
        assembly = perf_counter()
        disk = DiskGreen(side, band)
        bare = disk.dense() if side == 32 else None
        src, rcv, _ = port_point_matrices(disk, positions)
        networks = {}
        for q in COUNTS:
            for coupling in STRENGTHS:
                network = network_terms(disk, q, coupling)
                networks[q, coupling] = network
                op_id = f"n{side}_f{band}_a{q:02}_c{int(coupling):03}"
                if q == 8 and coupling == 48 and band == 400:
                    continue
                data = operator(disk, network)
                np.savez_compressed(HERE / "public/operators" / f"{op_id}.npz", **data)
        costs[f"f{band}_setup"] = dict(wall_seconds=perf_counter()-assembly, all_point_ports=len(positions))
        bare_action = lambda a: disk.kernel.apply(a[:, None])[:, 0] if a.ndim == 1 else disk.kernel.apply(a)
        for item in args.items:
            case_id = f"n{side}_c{item:02}"
            with np.load(HERE / "private" / f"{case_id}.npz") as label:
                m = label["material"].ravel()
            with np.load(HERE / "public/cases" / f"{case_id}.npz") as obs:
                sigma, sigma_hold = float(obs["sigma"]), float(obs["sigma_hold"])
            current, cost = physical_solve(m, bare_action, src, disk.kernel.diagonal, bare)
            costs[f"{case_id}_f{band}"] = cost
            bare_response = rcv @ current
            for q in COUNTS:
                indices = geometry_indices[q]
                response = bare_response[np.ix_(indices, indices)]
                for strength in STRENGTHS:
                    for eta in ERROR_RELATIVE:
                        key = tag(band, q, strength, eta)
                        row, y, y_hold, nominal, actual, power = response_row(
                            disk, case_id, q, strength, eta, response, networks[q, strength], sigma, sigma_hold)
                        rows.append(row)
                        all_data[case_id][key+"_Y"] = y
                        all_data[case_id][key+"_Y_hold"] = y_hold
                        all_data[case_id][key+"_radius"] = row["radius_noise"]
                        all_data[case_id][key+"_radius_hold"] = row["radius_hold_noise"]
                        private_transfer[case_id][key+"_full_true_transfer"] = actual
                        if eta == 0:
                            private_transfer[case_id][key+"_full_nominal_transfer"] = nominal
                        port_powers[case_id][key] = power
            print(case_id, band, "MHz", len(positions), "unique source positions",
                  f"state residual max {max(cost['state_residual']):.2e}, wall {cost['seconds_wall']:.2f}s", flush=True)
    for case_id, arr in all_data.items():
        np.savez_compressed(HERE / "public/cases" / f"{case_id}_mimo.npz", **arr)
        np.savez_compressed(HERE / "private" / f"{case_id}_mimo_transfer.npz", **private_transfer[case_id])
        (HERE/"private"/f"{case_id}_mimo_power.json").write_text(json.dumps(port_powers[case_id]))
    with (HERE/"results"/f"array_grid_n{side}.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    costs["end_to_end"] = dict(seconds_wall=perf_counter()-start_all[0], seconds_CPU=process_time()-start_all[1])
    (HERE/"results"/f"array_grid_n{side}_cost.json").write_text(json.dumps(costs, indent=2))


if __name__ == "__main__":
    main()
