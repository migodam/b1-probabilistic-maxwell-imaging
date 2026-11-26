"""Hold self-impedance, radiative resistance and all geometry fixed."""
import csv
from pathlib import Path

import numpy as np

from model import COUNTS, DiskGreen, antenna_forward, network_terms, physical_solve, power_report
from sweep import change_network, clean_transfer, off_diagonal_impedance

HERE = Path(__file__).resolve().parent


def main():
    rows = []
    with np.load(HERE/"results"/"mechanism_witness.npz") as w:
        a, b, sigma = w["with_region"].ravel(), w["without_region"].ravel(), float(w["reference_sigma"])
    disk = DiskGreen(32, 400)
    dense = disk.dense()
    for q in COUNTS:
        net = network_terms(disk, q, 48.0)
        data = [physical_solve(m, disk.kernel.apply, net["excitation"], disk.kernel.diagonal, dense)
                for m in (a,b)]
        response = [net["readout"] @ rays for rays, _ in data]
        reactive_off_diagonal = off_diagonal_impedance(net["antenna"]).imag
        for scale in (0.0, 1.0, 2.0, 4.0):
            tuned = change_network(net, 1j * (scale-1) * reactive_off_diagonal / 50)
            current_a, impedance_a = antenna_forward(response[0], tuned)[1:]
            current_b, impedance_b = antenna_forward(response[1], tuned)[1:]
            power_a = power_report(current_a, impedance_a, tuned["drive"])
            power_b = power_report(current_b, impedance_b, tuned["drive"])
            full_a, full_b = [clean_transfer(p, tuned) for p in response]
            block = (full_a-full_b)[np.ix_(tuned["rx"], np.arange(q))]/np.sqrt(q)
            rows.append(dict(q=q, reactive_offdiag_multiplier=scale,
                             reactive_offdiag_2ohms=scale*np.linalg.norm(reactive_off_diagonal, 2),
                             pair_distance_over_two_r=np.linalg.norm(block)/(3*sigma*q),
                             loaded_condition_1=np.linalg.cond(impedance_a+50*np.eye(len(impedance_a))),
                             loaded_condition_2=np.linalg.cond(impedance_b+50*np.eye(len(impedance_b))),
                             accepted_W1=sum(power_a["driven_accepted_W"]),
                             accepted_W2=sum(power_b["driven_accepted_W"]),
                             min_resistance_ohms=min(power_a["impedance_min_hermitian_eigen_ohms"],
                                                     power_b["impedance_min_hermitian_eigen_ohms"]),
                             shared_solve_time=data[0][1]["seconds_wall"]+data[1][1]["seconds_wall"]))
    with (HERE/"results"/"coupling_control.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for d in rows:
        print(d["q"], d["reactive_offdiag_multiplier"], d["pair_distance_over_two_r"], d["loaded_condition_1"])


if __name__ == "__main__":
    main()
