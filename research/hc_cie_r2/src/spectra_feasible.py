"""Follow up with a bounded physical anchor refined against the original observation."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import time

import spectra as sp
import numpy as np
from scipy.linalg import lu_factor, lu_solve, svd
from scipy.optimize import minimize


def audit_candidates(observation, design, follow):
    g, s, b, y = observation.g, observation.s, observation.b, observation.y
    bounds = [(observation.lower, observation.upper)] * len(b)
    norm_y_squared = np.linalg.norm(y)**2

    def objective(material):
        system = -material[:, None] * g
        system.flat[::len(b) + 1] += 1
        factors = lu_factor(system, overwrite_a=True, check_finite=False)
        j = lu_solve(factors, material[:, None] * b, check_finite=False)
        error = s @ j - y
        adjoint = lu_solve(factors, s.conj().T @ error, trans=2, check_finite=False)
        gradient = np.sum((adjoint.conj() * (b + g @ j)).real, axis=1) / norm_y_squared
        return float(0.5 * np.linalg.norm(error)**2 / norm_y_squared), gradient

    fit = minimize(objective, observation.base, jac=True, method="L-BFGS-B", bounds=bounds,
                   options=dict(maxiter=follow["iteration_cap"], maxfun=follow["function_cap"],
                                ftol=1e-12, gtol=1e-9))
    u, receive, vh = svd(s, full_matrices=True, check_finite=False)
    rank = np.count_nonzero(receive > follow["relative_receive_cut"] * receive[0])
    dark = vh[rank:].conj().T
    j0 = (vh[:rank].conj().T / receive[:rank]) @ (u[:, :rank].conj().T @ y)
    g_off = (g + g.T) * 0.5
    np.fill_diagonal(g_off, 0)
    gamma = sp.gamma_spectrum(g_off, dark, design["action_seed"])
    model = sp._reciprocal_model(dark, b)
    design = dict(design, contrast_scales=follow["contrast_scales"],
                  save_contrast_scales=follow["save_contrast_scales"])
    stats, arrays = sp._contrast_partial(replace(observation, base=fit.x), dark, model,
                                        j0, b, y, g_off, gamma, design)
    # The same-data minimum-norm affine candidate is a separate control.
    affine = arrays["affine_current"]
    field = b + g @ affine
    material = np.clip(np.sum((field.conj() * affine).real, axis=1) / np.sum(abs(field)**2, axis=1),
                       *bounds[0])
    forward_affine = sp._current(g, b, material)
    info = dict(
        physical_material_optimization=dict(iterations=fit.nit, evaluations=fit.nfev,
                                            status=fit.status, stop=str(fit.message), loss=fit.fun),
        unrefined_seed=observation.provenance, reduced_rank=int(rank),
        affine_only_port_over_delta=float(np.linalg.norm(s @ affine - y) / observation.delta),
        affine_only_physical_readout_port_over_delta=float(np.linalg.norm(s @ forward_affine - y) / observation.delta),
        affine_only_closure=float(np.linalg.norm(affine - material[:, None] * field) / np.linalg.norm(affine)),
        physical_audit_forward_solves=2 * len(stats) + 1,
        physical_audit_forward_rhs=(2 * len(stats) + 1) * b.shape[1],
        optimizer_forward_rhs=fit.nfev * b.shape[1], optimizer_adjoint_rhs=fit.nfev * b.shape[1],
    )
    arrays.update(initial_material=observation.base, affine_only_box_material=material)
    return stats, arrays, info


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, default=Path.cwd())
    parser.add_argument("--spectra-output", type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve().parents[1]
    folder = args.spectra_output or here / "results/spectra"
    design = json.loads((here / "configs/spectra.json").read_text())
    follow = json.loads((here / "configs/spectra_feasible.json").read_text())
    baseline = json.loads((folder / "partial.json").read_text())
    best = min((row for row in baseline if row["scale"] == 1), key=lambda r: r["physical_port_over_delta"])
    clock, cpu = time.perf_counter(), time.process_time()
    obs = next(o for o in sp.load_observations(args.workspace.resolve(), design) if o.key == best["observation"])
    stats, arrays, info = audit_candidates(obs, design, follow)
    np.savez_compressed(folder / "feasible_followup.npz", **arrays)
    sp.write_csv(folder / "feasible_followup.csv", stats)
    info.update(observation=obs.key, unrefined_port_over_delta=best["physical_port_over_delta"],
                wall_s=time.perf_counter() - clock, cpu_s=time.process_time() - cpu, contrast=stats)
    (folder / "feasible_followup.json").write_text(json.dumps(info, ensure_ascii=False, indent=2) + "\n")
    print(obs.key, info["physical_material_optimization"], flush=True)
    for r in stats:
        print(r["scale"], r["physical_port_over_delta"], r["alternative_physical_readout_port_over_delta"],
              r["full_quadratic_norm"] / r["full_linear_norm"], flush=True)


if __name__ == "__main__":
    main()
