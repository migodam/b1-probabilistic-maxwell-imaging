"""Free-face curvature and equal-norm finite-recovery interventions, without truth."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys
import time
import zlib

from imaging_run import R13_TREE, json_write, public_input

import numpy as np
from scipy.linalg import eigh, lu_solve, solve
from scipy.optimize import minimize
from scipy.stats import chi2
from threadpoolctl import threadpool_info


def prior_terms(model, material):
    reg = model.reg
    delta = reg.D @ material
    eps2 = reg.epsilon**2
    exact_edge_curvature = (2 * reg.smooth +
                            reg.tv * eps2 / (delta**2 + eps2)**1.5) / reg.edges
    value, grad, _ = reg.evaluate(material)
    hess = (reg.D.T @ reg.D.multiply(exact_edge_curvature[:, None])).toarray()
    return value, grad, hess


def material_derivatives(model, state, Y, material, full):
    """All contractions use the original unwhitened relative complex data term."""
    norm2 = np.vdot(Y, Y).real
    residual = state["Y"] - Y
    u = state["Q"]
    with model._physics.operation("MATERIAL_RECEIVER", rhs=model.nr, receiver_transpose_RHS=model.nr):
        t = lu_solve(state["lu"], model.S.T, trans=1, check_finite=False).T
    jac = np.einsum("ri,is->rsi", t, u).reshape(residual.size, model.n)
    jr = np.ascontiguousarray(jac.real)
    ji = np.ascontiguousarray(jac.imag)
    jac_gram = (jr.T @ jr + ji.T @ ji) / norm2
    _, prior_gradient, prior_hessian = prior_terms(model, material)
    gradient = (jr.T @ residual.real.ravel() + ji.T @ residual.imag.ravel()) / norm2 + prior_gradient
    gauss_newton = jac_gram + prior_hessian
    if not full:
        return gradient, gauss_newton, None, None, None

    with model._physics.operation("MATERIAL_INTERNAL", rhs=model.n, material_propagation_RHS=model.n):
        internal = lu_solve(state["lu"], model.G.T, trans=1, check_finite=False).T
    weighted_ports = t.T @ residual.conj() @ u.T
    cross = (weighted_ports * internal).real / norm2
    data_curvature = cross + cross.T
    exact_hessian = gauss_newton + data_curvature
    return gradient, gauss_newton, exact_hessian, data_curvature, prior_hessian


def gradient_check(model, material, Y, hessian, lower, upper, rng):
    direction = rng.standard_normal(len(material))
    direction[(material <= lower + .004) | (material >= upper - .004)] = 0
    norm = np.linalg.norm(direction)
    if norm == 0:
        return None
    direction /= norm
    steps = [.004, .002]
    errors = []
    exact = hessian @ direction
    for step in steps:
        forward = model.data_value_gradient(material + step * direction, Y)[1]
        backward = model.data_value_gradient(material - step * direction, Y)[1]
        fd = (forward - backward) / (2 * step)
        errors.append(float(np.linalg.norm(fd - exact) / max(np.linalg.norm(exact), 1e-100)))
    return {"steps": steps, "relative_Hv_errors": errors}


def direction_diagnostics(direction, free):
    squares = direction**2
    participation = float(1 / (squares @ squares))
    iy, ix = np.indices((24, 24), dtype=float)
    coords = np.column_stack((iy.ravel(), ix.ravel()))
    center = squares @ coords
    radius = np.sqrt(np.sum(squares * np.sum((coords - center)**2, axis=1)))
    ordered = np.sort(squares)[::-1]
    return {"participation_pixels": participation,
            "participation_fraction_free": participation / np.count_nonzero(free),
            "mass90_pixels": int(np.searchsorted(ordered.cumsum(), .90) + 1),
            "energy_rms_radius_cells": float(radius),
            "weighted_center_rc": center.tolist()}


def ray(material, v, norm, lower, upper):
    """Find the requested realized radius along a clipped, monotone-norm ray."""
    moves = v != 0
    span = np.where(v > 0, upper - material, material - lower)
    max_norm = float(np.linalg.norm(span[moves]))
    if norm > max_norm:
        return None
    lo, hi = 0., norm
    work = np.clip(material + hi * v, lower, upper)
    while np.linalg.norm(work - material) < norm and hi < 1e12:
        hi *= 2
        work = np.clip(material + hi * v, lower, upper)
    if np.linalg.norm(work - material) < norm:
        return None
    for _ in range(64):
        mid = .5 * (lo + hi)
        proposal = np.clip(material + mid * v, lower, upper)
        if np.linalg.norm(proposal - material) < norm:
            lo = mid
        else:
            hi = mid
    return np.clip(material + hi * v, lower, upper), hi


def run_one(task):
    cfg, case, endpoint, arm = task
    cpu0, wall0 = time.process_time(), time.perf_counter()
    out = Path(cfg["output"]) / f"{case}__{arm}"
    out.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(R13_TREE))
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel

    public, start, y_saved, paths = public_input(case, endpoint)
    model = MaxwellModel(public)
    Y = public["Y"]
    scale = float(np.linalg.norm(Y))
    delta = float(public.get("delta", float(public.get("sigma", 0.)) *
                    np.sqrt(chi2.ppf(.99, 2 * Y.size) / 2)))
    lower, upper = float(public.get("lower", .02)), float(public.get("upper", 20.))
    cap = cfg["budget_cpu_seconds"]
    rows, materials, ports = [], [], []
    winner = None
    new_seed = None
    direction = np.zeros(len(start))
    tol = cfg["active_tol"]
    face_lo, face_hi = start <= lower + tol, start >= upper - tol
    free = ~(face_lo | face_hi)
    diagnostics = dict(free_cells=int(free.sum()), lower_cells=int(face_lo.sum()), upper_cells=int(face_hi.sum()))
    cost = Counter()
    started = time.process_time()
    state = model.state(start)
    disagreement = float(np.linalg.norm(y_saved - state["Y"]) / max(np.linalg.norm(y_saved), 1e-100))

    class BudgetEnd(Exception):
        """The next action would start outside the shared paid prefix."""

    def check():
        if time.process_time() - cpu0 >= cap:
            raise BudgetEnd

    def audit(m, physical, label, details=None):
        nonlocal winner
        clock = time.process_time() - cpu0
        ym = physical["Y"]
        resid = float(np.linalg.norm(ym - Y))
        d = .5 * (resid / scale)**2
        r = model.reg.evaluate(m)[0]
        if not np.all(np.isfinite(m)) or m.min() < lower - 1e-12 or m.max() > upper + 1e-12:
            raise ValueError("An original-Maxwell candidate is outside the material box")
        index = len(rows)
        rows.append(dict(i=index, stage=label, cpu=clock, wall=time.perf_counter() - wall0,
                         Phi=d+r, D=d, R=r, relative_port_error=resid/scale,
                         residual_over_delta=resid/delta if delta else None,
                         paid=clock<=cap, setup_counts=dict(cost), original_counts=dict(model.cost),
                         **(details or {})))
        materials.append(m.copy())
        ports.append(ym.copy())
        if clock <= cap and (winner is None or rows[winner]["Phi"] > d+r):
            winner = index
        return index

    if disagreement > 1e-7:
        raise ValueError("Old port prediction disagrees with the original operator")
    audit(start, state, "endpoint")
    rng = np.random.default_rng(cfg["rng_seed"] + zlib.crc32(case.encode()))
    hessian = gradient = data_curvature = prior_hessian = gn = None
    direction_cpu = 0.
    termination = None
    unique_rays = []
    try:
        if arm in {"curvature", "gauss_newton"}:
            gradient, gn, hessian, data_curvature, prior_hessian = material_derivatives(
                model, state, Y, start, full=arm == "curvature")
            cost["jacobian_material_columns"] += model.n
            cost["jac_gram"] += 1
        if np.any(free):
            if arm == "curvature":
                evals, evecs = eigh(hessian[np.ix_(free, free)], check_finite=False)
                direction[free] = evecs[:, 0]
                cost["eigensystem_free_dimension"] = int(free.sum())
                diagnostics.update(
                    free_min_eigenvalue=float(evals[0]), free_max_eigenvalue=float(evals[-1]),
                    free_negative_eigenvalues=int(np.count_nonzero(evals < -cfg["negative_eigenvalue_tol"])),
                    free_min_is_negative=bool(evals[0] < -cfg["negative_eigenvalue_tol"]),
                    curvature_from_GN=float(direction @ gn @ direction),
                    curvature_from_residual=float(direction @ data_curvature @ direction),
                    curvature_from_prior=float(direction @ prior_hessian @ direction),
                    gradient_projection_inf=float(np.linalg.norm(start-np.clip(start-gradient, lower, upper), np.inf)),
                    free_gradient_inf=float(np.linalg.norm(gradient[free], np.inf)),
                    lower_release_count=int(np.count_nonzero(face_lo & (gradient < -1e-8))),
                    upper_release_count=int(np.count_nonzero(face_hi & (gradient > 1e-8))),
                    Hv_check=gradient_check(model, start, Y, hessian, lower, upper, rng),
                )
            elif arm == "gauss_newton":
                reduced = gn[np.ix_(free, free)]
                ridge = max(1e-4*np.trace(reduced)/len(reduced), 1e-14)
                direction[free] = solve(reduced + ridge*np.eye(len(reduced)), -gradient[free],
                                        assume_a="pos", check_finite=False)
                cost["GN_dimension"] = int(free.sum())
                diagnostics["damping"] = float(ridge)
            elif arm == "random":
                direction[free] = rng.standard_normal(free.sum())
        if np.linalg.norm(direction) > 0:
            direction /= np.linalg.norm(direction)
            diagnostics.update(direction_diagnostics(direction, free))
            direction_cpu = time.process_time() - started
            for sign in [-1, 1]:
                v = sign * direction
                bound_dist = np.where(v > 0, upper - start, start - lower)
                nonzero = np.abs(v) > 1e-14
                diagnostics[f"straight_step_max_sign_{sign}"] = float(np.min(bound_dist[nonzero]/np.abs(v[nonzero])))
                for target in cfg["radius_L2"]:
                    check()
                    candidate = ray(start, v, target, lower, upper)
                    if candidate is None:
                        cost["unreachable_radii"] += 1
                        continue
                    trial, a = candidate
                    if any(np.array_equal(trial, previous) for previous in unique_rays):
                        cost["duplicates_removed"] += 1
                        continue
                    unique_rays.append(trial)
                    change = trial - start
                    norm = np.linalg.norm(change)
                    details = {"radius": target, "signed_direction": sign,
                               "actual_radius": float(norm), "ray_parameter": a,
                               "direction_cosine": float(change @ v / norm),
                               "n_changed": int(np.count_nonzero(np.abs(change) > 1e-10)),
                               "n_newly_clipped": int(np.count_nonzero(
                                   free & ((trial <= lower + tol) | (trial >= upper - tol)))),
                               "linear_term": None if gradient is None else float(gradient @ change),
                               "quad_term": None if hessian is None else .5*float(change @ hessian @ change)}
                    idx = audit(trial, model.state(trial), "jump", details)
                    if rows[idx]["paid"] and (new_seed is None or rows[new_seed]["Phi"] > rows[idx]["Phi"]):
                        new_seed = idx
        seed = start if new_seed is None else materials[new_seed]

        def objective(m):
            check()
            value, grad = model.data_value_gradient(m, Y)
            audit(m, model.state(m), "polish")
            check()
            return 1e6*value, 1e6*grad

        result = minimize(objective, seed, method="L-BFGS-B", jac=True,
                          bounds=[(lower, upper)]*model.n,
                          options=dict(maxiter=6000, maxfun=24000, ftol=1e-12, gtol=.01, maxls=30))
        termination = str(result.message)
    except BudgetEnd:
        termination = "PROCESS_CPU_LIMIT"

    np.savez_compressed(out / "PREDICTIONS.npz", material=np.asarray(materials), direction=direction,
                        free=free, selected_port=ports[winner], endpoint_port=ports[0],
                        selected_index=np.asarray(winner))
    json_write(out / "TRACE.json", dict(case=case, endpoint=endpoint, arm=arm, without_truth=True,
               diagnostics=diagnostics, input_paths=paths, selected=winner, selected_jump=new_seed,
               total_cpu=time.process_time()-cpu0, total_wall=time.perf_counter()-wall0,
               direction_cpu=direction_cpu, termination=termination, endpoint_port_disagreement=disagreement,
               configuration={k:v for k,v in cfg.items() if k!="output"},
               threadpool=threadpool_info(), cost=dict(model.cost), auxiliary_cost=dict(cost), chain=rows))
    return dict(case=case, endpoint=endpoint, arm=arm, directory=out.name,
                selected=winner, completed_audits=len(rows), final_phi=rows[winner]["Phi"],
                free_min_eigenvalue=diagnostics.get("free_min_eigenvalue"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2, choices=(1, 2))
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    out = args.out.absolute()
    out.mkdir(parents=True, exist_ok=False)
    cfg["output"] = str(out)
    tasks = []
    for i, (case, endpoint) in enumerate(cfg["cases"]):
        arms = cfg["arms"][i % len(cfg["arms"]):] + cfg["arms"][:i % len(cfg["arms"])]
        tasks.extend((cfg, case, endpoint, arm) for arm in arms)
    started = time.perf_counter()
    jobs = []
    with ProcessPoolExecutor(args.workers) as pool:
        for row in pool.map(run_one, tasks):
            jobs.append(row)
            print(row, flush=True)
    json_write(out / "PREDICTION_SEAL.json",
               dict(id=cfg["id"], batch_wall=time.perf_counter()-started,
                    without_truth=True, workers=args.workers, expected=len(tasks), jobs=jobs))


if __name__ == "__main__":
    main()
