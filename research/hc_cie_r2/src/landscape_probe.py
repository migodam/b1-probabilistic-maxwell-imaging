"""Compare original physical centered second differences inside the unclipped box."""
import argparse
from pathlib import Path
import sys
import time

from imaging_run import R13_TREE, json_write, public_input
from imaging_score import read, write_csv

import numpy as np


def run(sealed, output):
    start_batch = time.perf_counter()
    sys.path.insert(0, str(R13_TREE))
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel

    seal = read(sealed / "PREDICTION_SEAL.json")
    if not seal["without_truth"] or seal["expected"] != len(seal["jobs"]):
        raise ValueError("Incomplete public seal for fixed directions")
    entries = []
    for job in seal["jobs"]:
        if job["arm"] != "curvature":
            continue
        start_cpu = time.process_time()
        directory = sealed / job["directory"]
        trace = read(directory / "TRACE.json")
        with np.load(directory / "PREDICTIONS.npz", allow_pickle=False) as z:
            u = z["direction"]
            archived = z["material"][0]
        public, material, _, paths = public_input(job["case"], job["endpoint"])
        if not np.array_equal(material, archived):
            raise ValueError("Historical input changed after seal")
        Y = public["Y"]
        y_scale = np.vdot(Y, Y).real
        lower, upper = float(public.get("lower", .02)), float(public.get("upper", 20.))
        model = MaxwellModel(public)
        center = model.state(material)["Y"]
        r0 = center - Y
        b = model.reg.D @ material
        db = model.reg.D @ u
        r_value = model.reg.evaluate(material)[0]
        base = .5*np.vdot(r0, r0).real / y_scale + r_value
        both = np.abs(u) > 0
        half_width = float(np.min(np.minimum(upper - material[both], material[both] - lower) / np.abs(u[both])))
        curvature = trace["diagnostics"]["free_min_eigenvalue"]
        source_rows = []
        for fraction in [.25, .5]:
            a = fraction * half_width
            delta_ports = []
            finite = {}
            for sign in [-1, 1]:
                m = material + sign*a*u
                if m.min() < lower or m.max() > upper:
                    raise ValueError("Unclipped derivative probe must stay in the material box")
                observed = model.state(m)["Y"]
                res = observed - Y
                phi = .5*np.vdot(res, res).real/y_scale + model.reg.evaluate(m)[0]
                delta_ports.append(observed-center)
                finite[sign] = float(phi)
            dm, dp = delta_ports
            q = model.reg.epsilon**2
            edge_step = a * db
            base_root = np.sqrt(b**2+q)
            tv_plus = (2*b*edge_step+edge_step**2)/(np.sqrt((b+edge_step)**2+q)+base_root)
            tv_minus = (-2*b*edge_step+edge_step**2)/(np.sqrt((b-edge_step)**2+q)+base_root)
            stable_prior = 2*model.reg.smooth*np.mean(edge_step**2) + model.reg.tv*np.mean(tv_plus+tv_minus)
            stable_data = (np.vdot(r0, dm+dp).real+.5*(np.vdot(dm, dm).real+np.vdot(dp, dp).real))/y_scale
            actual = finite[-1]+finite[1]-2*base
            stable = stable_data + stable_prior
            predicted = curvature*a*a
            row = dict(case=job["case"], endpoint=job["endpoint"], common_straight_radius=half_width,
                       radius_fraction=fraction, radius=a, lambda_min=curvature, Phi_base=base,
                       Phi_minus=finite[-1], Phi_plus=finite[1], ordinary_central_difference=actual,
                       stable_central_difference=stable, predicted=float(predicted),
                       relative_stable_error=float(abs(stable-predicted)/max(abs(predicted), 1e-100)),
                       state_full_LU=model.cost["LU"], full_RHS=model.cost["actual_complex_RHS"],
                       cpu_through_probe=time.process_time()-start_cpu, observation=paths[-2])
            source_rows.append(row)
        entries.extend(source_rows)
    output.mkdir(parents=True, exist_ok=False)
    write_csv(output/"FREE_WINDOW.csv", entries)
    json_write(output/"RECEIPT.json", dict(n_cases=len(entries)//2, radius_fractions=[.25, .5],
               no_truth_read=True, full_lu=sum(row["state_full_LU"] for row in entries[1::2]),
               total_cpu=sum(row["cpu_through_probe"] for row in entries[1::2]),
               total_wall=time.perf_counter()-start_batch))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sealed", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.sealed, args.out)


if __name__ == "__main__":
    main()
