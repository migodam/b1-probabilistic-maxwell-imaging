"""Anchor-free bihomogeneous POLSYS_PLP control on held-out balanced cases."""

import argparse
from collections import Counter
from importlib.metadata import version
from math import comb
from time import perf_counter, process_time

from branch_discovery import Experiment, HI, LO, MATERIAL_TOL, ROOT, rms
import json
import numpy as np
import pypolsys


def coefficients(p):
    u, sigma, _ = np.linalg.svd(p.b, full_matrices=True)
    rank = np.sum(sigma > 1e-12 * sigma[0])
    w = p.obs @ np.linalg.pinv(p.b, rcond=1e-12)
    v = p.r + w @ p.h
    q = u[:, rank:].conj().T * max(np.linalg.norm(w), 1e-12)
    qh = q @ p.h
    d = len(q)
    if 2 * d != p.n:
        raise ValueError("This experiment needs balanced data.")
    # z and zet are independent complex variables. Real materials require zet = conj(z).
    # 2i F = (w + z Q)(conj(v) + zet conj(QH)) - (conj(w) + zet conj(Q))(v + z QH).
    constant = (w * v.conj()).imag[None, :]
    linear_z = (q * v.conj() - w.conj() * qh) / (2j)
    cross = (q[:, None, :] * qh[None, :, :].conj() - qh[:, None, :] * q[None, :, :].conj()) / (2j)
    c = np.vstack((constant, linear_z, linear_z.conj(), cross.reshape(d*d, -1))).T
    c = c / np.linalg.norm(c, axis=1)[:, None]
    eye = np.eye(2*d, dtype=np.int32)
    degrees = np.vstack((np.zeros(2*d, dtype=np.int32), eye,
                         (eye[:d, None, :] + eye[None, d:, :]).reshape(d*d, 2*d)))
    return w, q, c, degrees


def run(geometry, observation):
    wall, cpu = perf_counter(), process_time()
    p = Experiment(geometry, observation)
    w, q, c, degrees = coefficients(p)
    d = len(q)
    partitions = pypolsys.utils.make_mh_part(p.n, [list(range(1, d+1)), list(range(d+1, 2*d+1))])
    pypolsys.polsys.init_poly(p.n, np.full(p.n, c.shape[1], dtype=np.int32), c.reshape(-1),
                             np.tile(degrees, (p.n, 1)))
    pypolsys.polsys.init_partition(*partitions)
    setup = process_time() - cpu
    n_paths = pypolsys.polsys.solve(1e-9, 1e-14, np.array(0.), 0)
    continuation_cpu = process_time() - cpu - setup
    solve_status = int(pypolsys.polsys.solve_status)
    paths = []
    roots = []
    for i, (status, evaluations, endpoint) in enumerate(zip(pypolsys.polsys.path_status,
                                                            pypolsys.polsys.num_jac_eval,
                                                            pypolsys.polsys.myroots.T)):
        classification = "tracking_failed"
        row = {"path": i + 1, "status": int(status), "jacobian_evaluations": int(evaluations)}
        if status >= 11 and status % 10 == 1:
            h = endpoint[-1]
            z, zet = endpoint[:d], endpoint[d:2*d]
            row["homogenizing_coordinate_abs"] = float(abs(h))
            if not np.all(np.isfinite(endpoint)) or abs(h) <= 1e-8:
                classification = "at_infinity"
            else:
                monomials = np.r_[1., z, zet, np.outer(z, zet).reshape(-1)]
                row["polynomial_residual_norm"] = float(np.linalg.norm(c @ monomials))
                err = np.max(np.abs(zet - z.conj())) / max(1., np.max(np.abs(endpoint[:-1])))
                row["conjugacy_defect_relative"] = float(err)
                classification = "finite_nonreal"
                if err <= 1e-7:
                    wl = w + ((z + zet.conj()) / 2.) @ q
                    vl = p.r + wl @ p.h
                    classification = "finite_zero_node"
                    if np.min(np.abs(vl)) > 1e-15:
                        material = wl / vl
                        row["material_imaginary_max"] = float(np.max(np.abs(material.imag)))
                        row["material"] = material.real.tolist()
                        lo, hi = np.min(material.real), np.max(material.real)
                        if row["material_imaginary_max"] > 1e-7:
                            classification = "finite_material_not_real"
                        elif lo < LO - 1e-10 or hi > HI + 1e-10:
                            classification = "finite_real_outside_box"
                        else:
                            validated = p.validate(material.real)
                            classification = "finite_failed_full_VIE"
                            if validated is not None:
                                row["full_VIE"] = validated
                                classification = "accepted_duplicate"
                                if all(rms(material.real - root) > MATERIAL_TOL for root in roots):
                                    classification = "accepted_distinct"
                                    roots.append(material.real)
        row["classification"] = classification
        paths.append(row)
    pypolsys.polsys.cleanup_pol()
    pypolsys.polsys.cleanup_par()
    return {
        "expected_bihomogeneous_bound": comb(2*d, d), "paths_attempted": int(n_paths),
        "solve_status": solve_status, "path_classifications": dict(Counter(v["classification"] for v in paths)),
        "path_jacobian_evaluations": sum(p["jacobian_evaluations"] for p in paths),
        "construction_cpu_seconds": setup, "tracking_cpu_seconds": continuation_cpu,
        "total_cpu_seconds": process_time() - cpu, "total_wall_seconds": perf_counter() - wall,
        "validated_roots": [r.tolist() for r in roots], "paths": paths,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", help="Target ID(s), default all held-out six-cell targets")
    parser.add_argument("--tag", default="")
    args = parser.parse_args()
    case_file = json.loads((ROOT / "results/branch_holdout.json").read_text())
    observed = np.load(ROOT / "results/branch_holdout.npz")
    names = args.case or [k for k in case_file["cases"] if k.startswith("fresh6-")]
    solved = {}
    for name in names:
        geometry = case_file["cases"][name]["geometry"]
        solved[name] = run(geometry, observed[f"{name}/y"])
        print(name, solved[name]["path_classifications"], f"{solved[name]['total_cpu_seconds']:.4f} CPU seconds",
              flush=True)
    # Compare predictions only after the anchor-free solver has returned for every case.
    for name, result in solved.items():
        roots = np.asarray(result["validated_roots"])
        truth = observed[f"{name}/truth"]
        result["best_truth_RMS"] = min((float(rms(m - truth)) for m in roots), default=None)
        if f"{name}/anchor" in observed:
            anchor = observed[f"{name}/anchor"]
            result["matches_anchor"] = any(rms(m - anchor) <= MATERIAL_TOL for m in roots)
            result["method_coverage"] = {}
            for method in case_file["cases"][name]["methods"]:
                conditional_roots = np.vstack((anchor, observed[f"{name}/{method}"]))
                result["method_coverage"][method] = {
                    "conditional_roots": len(conditional_roots),
                    "homotopy_roots_matched": sum(any(rms(m - q) <= MATERIAL_TOL for q in conditional_roots) for m in roots),
                    "conditional_roots_matched": sum(any(rms(m - q) <= MATERIAL_TOL for m in roots) for q in conditional_roots),
                }
    (ROOT / "results" / f"homotopy_control{args.tag}.json").write_text(json.dumps({
        "package": "pypolsys / ACM TOMS 801 POLSYS_PLP", "version": version("pypolsys"),
        "library": "https://github.com/nennigb/pypolsys", "paper": "https://doi.org/10.1145/347837.347885",
        "cohort": "Held-out balanced clean Maxwell VIE", "cases": solved,
    }, indent=2))
