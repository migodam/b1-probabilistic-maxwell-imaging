"""Release a full-port synthetic calibration diagnostic for the public mechanism."""
import argparse
from pathlib import Path
import json

import numpy as np

from model import DiskGreen, antenna_forward, network_terms, physical_solve
from sweep import change_network, clean_transfer, normalized_reactance_error

HERE = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--eta", type=float, default=0.025)
    p.add_argument("--tag", default="")
    args = p.parse_args()
    w = np.load(HERE/"results"/"mechanism_witness.npz")
    disk = DiskGreen(32, 400)
    nominal = network_terms(disk, 8, 48)
    eta = args.eta
    direction = normalized_reactance_error(nominal["positions"])
    true = change_network(nominal, eta*direction)
    mat = w["with_region"].ravel()
    J, cost = physical_solve(mat, disk.kernel.apply, nominal["excitation"], disk.kernel.diagonal, disk.dense())
    bare = nominal["readout"]@J
    f_nominal, f_true = [clean_transfer(bare, a) for a in [nominal, true]]
    _, _, loaded_true = antenna_forward(bare, true)
    full_incident_normalization = 2 * 50.0
    m = true["matrix"]
    empty_full = np.eye(len(m)) - full_incident_normalization * true["inv"]
    loaded_full = np.eye(len(m)) - full_incident_normalization * np.linalg.solve(
        loaded_true+50*np.eye(len(m)), np.eye(len(m)))
    np.savez_compressed(HERE/"public"/f"calibration_example{args.tag}.npz", F_true=f_true, F_nominal=f_nominal,
                        empty_scattering_full=empty_full, loaded_scattering_full=loaded_full,
                        receive_index=true["rx"], source_index=np.arange(true["drive"].shape[1]),
                        reference_ohms=50.0, normalized_additive_error=eta*direction,
                        antenna_impedance_empty_nominal=nominal["antenna"],
                        antenna_impedance_empty_true=true["antenna"], antenna_impedance_loaded_true=loaded_true)
    data = dict(mode="separate numerical hypothetical all-terminal probe; not part of blind masked input",
                frequency_MHz=400, terminal_count=len(m), error_ohms=eta*50,
                full_Fro=float(np.linalg.norm(f_true)),
                full_spectral=float(np.linalg.norm(f_true, ord=2)),
                matrix_convention="Dimensionless scattering matrix; loaded - empty = 2R (M_empty^-1 - M_loaded^-1)",
                consistency_error=float(np.linalg.norm(loaded_full-empty_full-f_true)), physical_cost=cost)
    (HERE/"results"/f"calibration_example{args.tag}.json").write_text(json.dumps(data, indent=2))


if __name__ == "__main__":
    main()
