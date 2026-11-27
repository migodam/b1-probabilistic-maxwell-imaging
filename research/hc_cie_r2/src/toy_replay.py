#!/usr/bin/env python3
"""Replay the frozen disk-cell witnesses, and run paired material-box controls."""
import argparse
import json
import os
from pathlib import Path
import platform
from time import perf_counter, process_time
from scipy.optimize import root

import numpy as np
import scipy

from toy_models import (
    BOX_UPPER, MAX_EVALS, RANK_RTOL, SUCCESS_RTOL, CurrentProblem, MaterialProblem,
    Model, complex_array, complex_json, disk_operators, frozen_model, null_space, pack,
    wedge,
)


def make_fiber(model, n_basis=None):
    if n_basis is None:
        n_basis = null_space(model.receiver, rcond=RANK_RTOL)
    return CurrentProblem(model, n_basis, np.linalg.pinv(model.receiver, rcond=RANK_RTOL)
                          @model.data)


def audit_case(problem, position, truth_current, form, solve_info=None):
    f, a = problem.residual_and_jac(position, form)
    h = (problem.analytic_hessian(position) if form == "minors"
         else problem.hessian_box(position))
    m = problem.model
    return dict(position=position.tolist(), objective=float(.5*f@f),
                gradient_norm=float(np.linalg.norm(a.T@f)),
                hessian_eigenvalues=np.linalg.eigvalsh(h).tolist(),
                **m.audit(problem.current(position), truth_current),
                solver=solve_info)


def projection(fiber, truth_current, projector=None):
    t0, c0 = perf_counter(), process_time()
    n, g = fiber.basis, fiber.model.green
    gn = g@n
    gamma = (n[:, :0] if n.shape[1] == 1
             else (n[:, 0]*gn[:, 1]-n[:, 1]*gn[:, 0])[:, None])
    zero = np.zeros(2*np.prod(fiber.shape))
    b, a = fiber.residual_and_jac(zero)
    nrows = len(b)//2
    b = b[:nrows]+1j*b[nrows:]
    a = a[:nrows, :len(zero)//2]+1j*a[nrows:, :len(zero)//2]
    p = np.eye(len(g), dtype=complex)-gamma@np.linalg.pinv(gamma, rcond=RANK_RTOL)
    if projector is not None:
        p = projector
    lin = p@a
    sv = np.linalg.svd(lin, compute_uv=False)
    rank = np.count_nonzero(sv > RANK_RTOL*sv[0])
    solution = np.linalg.lstsq(lin, -p@b, rcond=RANK_RTOL)[0]
    elapsed, cpu = perf_counter()-t0, process_time()-c0
    j = fiber.current(pack(solution))
    pj = n@n.conj().T
    gram = pj*(g@pj@g.conj().T) - (pj@g.conj().T)*(g@pj)
    measures = dict(
        operator_rank=int(rank), expected_rank=int(a.shape[1]), numerical_cutoff=RANK_RTOL,
        gamma_rank=int(gamma.shape[1]), singular_values=sv.tolist(),
        retained_rank=int(np.linalg.matrix_rank(p, tol=RANK_RTOL)), alpha_squared=float(sv[-1]**2),
        nonlinear_term_remaining=float(np.linalg.norm(p@gamma)),
        gram_identity_error=float(np.linalg.norm(gram-gamma@gamma.conj().T)),
        current_absolute_error=float(np.linalg.norm(j-truth_current)),
        raw_material_absolute_error=float(np.linalg.norm(fiber.model.box_readout(j)[0]-fiber.model.truth)),
        **fiber.model.audit(j, truth_current), elapsed_seconds=elapsed, cpu_seconds=cpu,
    )
    return measures, (p, a, b, gamma)


def random_projection_controls(fiber, truth_current, starts, reference, seed=23041):
    rng = np.random.default_rng(seed)
    p_reference, a, b, gamma = reference
    controls = []
    for index in range(32):
        n = len(fiber.model.green)
        u = rng.normal(size=n) + 1j*rng.normal(size=n)
        u /= np.linalg.norm(u)
        p = np.eye(n) - np.outer(u, u.conj())
        shortcut, _ = projection(fiber, truth_current, projector=p)
        t0, c0 = perf_counter(), process_time()
        a0 = np.linalg.lstsq(p@a, -p@b, rcond=RANK_RTOL)[0].reshape(2, 2)
        a1 = np.linalg.lstsq(p@a, -p@gamma, rcond=RANK_RTOL)[0].reshape(2, 2)
        mixed = (a0[0, 0]*a1[1, 1]+a1[0, 0]*a0[1, 1]
                 -a0[0, 1]*a1[1, 0]-a1[0, 1]*a0[1, 0])
        scalar_roots = np.roots([np.linalg.det(a1), mixed-1, np.linalg.det(a0)])
        root_records = []
        for d in scalar_roots:
            z = pack(a0+d*a1)
            j = fiber.current(z)
            root_records.append(dict(
                det_closure_error=float(abs(np.linalg.det(a0+d*a1)-d)),
                projected_residual=float(np.linalg.norm(p@wedge(j, fiber.model.field(j)))),
                **fiber.model.audit(j, truth_current)))
        exact_reduced_seconds = perf_counter()-t0
        exact_reduced_cpu = process_time()-c0
        nonlinear = {}
        for label, start in starts.items():
            z, cost = fiber.solve(start, "minors", projector=p)
            nonlinear[label] = dict(**fiber.model.audit(fiber.current(z), truth_current), **cost)
        controls.append(dict(seed=seed, index=index, removed_direction=complex_json(u),
                             relative_projection_distance=float(np.linalg.norm(p-p_reference)),
                             invalid_linear_shortcut=shortcut, nonlinear_projected_lm=nonlinear,
                             quadratic_scalar_roots=root_records,
                             reduced_seconds=exact_reduced_seconds, reduced_cpu_seconds=exact_reduced_cpu))
    return controls


def derivative_checks(problem, midpoint, direction, material, x, h):
    f, j = problem.residual_and_jac(midpoint)
    _, jbox = problem.residual_and_jac(midpoint, "real_box")
    hess = problem.analytic_hessian(midpoint)
    box_hess = problem.hessian_box(midpoint)
    s = 0.17
    fg = lambda m: problem.residual_and_jac(m)[0]
    directed_difference = fg(midpoint+s*direction)-fg(midpoint-s*direction)
    jac_errors = []
    for step in (1e-4, 1e-5, 1e-6):
        xp, xm = midpoint+step*direction, midpoint-step*direction
        jac_fd = (fg(xp)-fg(xm))/(2*step)
        hess_fd = (problem.gradient(xp)-problem.gradient(xm))/(2*step)
        f_plus = problem.residual_and_jac(xp, "real_box")[0]
        f_minus = problem.residual_and_jac(xm, "real_box")[0]
        g_plus = problem.gradient(xp, "real_box")
        g_minus = problem.gradient(xm, "real_box")
        delta = step*h
        mat_jac = material.jac(x)
        mat_plus = material.residual(x+delta)
        mat_minus = material.residual(x-delta)
        hessian_fd = (.5*mat_plus@mat_plus + .5*mat_minus@mat_minus
                      -material.residual(x)@material.residual(x))/step**2
        jac_errors.append(dict(
            step=step,
            current_jacobian_error=float(np.linalg.norm(jac_fd-j@direction)),
            current_hessian_error=float(np.linalg.norm(hess_fd-hess@direction)),
            real_box_jacobian_error=float(np.linalg.norm((f_plus-f_minus)/(2*step)-jbox@direction)),
            real_box_hessian_error=float(np.linalg.norm((g_plus-g_minus)/(2*step)-box_hess@direction)),
            physical_jacobian_error=float(np.linalg.norm((mat_plus-mat_minus)/(2*step)-mat_jac@h)),
            physical_directional_hessian_error=float(abs(hessian_fd-material.objective_hessian(x, h))),
        ))
    column = np.arange(1, problem.basis.shape[1]+1)+0.2j
    rank_one = np.outer(column, [0.8-0.2j, -1.0j])
    r1, t = pack(rank_one), 0.25
    recover = fg(midpoint+t*r1)+fg(midpoint-t*r1)-2*f
    return dict(
        midpoint_secant_error=float(np.linalg.norm(directed_difference-2*s*(j@direction))),
        rank_one_centered_second_difference=float(np.linalg.norm(recover)),
        physical_finite_response_error=material.finite_response(x, x+1.4*h),
        real_box_raw_readout_at_diagnostic_point=complex_json(problem.model.box_readout(
            problem.current(midpoint))[0]), finite_differences=jac_errors,
    )


def run_case(directory, case):
    model, source = frozen_model(directory, case)
    physical = MaterialProblem(model)
    truth_j = model.forward(model.truth)
    case_seed = 200311 if case == "soft3" else 200312
    rng = np.random.default_rng(case_seed)
    if case == "soft3":
        fiber = make_fiber(model)
        problem = CurrentProblem(model, np.eye(len(truth_j)), np.zeros_like(truth_j),
                                 include_data=True, weight=source["weight"])
        seed_current = complex_array(source, "J_bad")
    else:
        n = complex_array(source, "N")
        fiber = make_fiber(model, n)
        problem = fiber
        seed_current = complex_array(source, "J0")+n@complex_array(source, "Z_bad")

    refined, refinements = problem.refine_bad(problem.coords(seed_current))
    bad_current = problem.current(refined)
    bad = audit_case(problem, refined, truth_j, "minors", refinements)
    real_box_at_bad = audit_case(problem, refined, truth_j, "real_box")
    projected, operators = projection(fiber, truth_j)
    diagnostic = None
    if case == "hard5":
        truth_position = fiber.coords(truth_j)
        _, derivative = fiber.residual_and_jac(truth_position)
        gamma_norm = np.linalg.norm(operators[3])
        sigma_min = np.linalg.svd(derivative, compute_uv=False)[-1]
        diagnostic = dict(
            lambda_at_truth=float(sigma_min),
            quadratic_bound=float(gamma_norm/2),
            sufficient_onepoint_radius=float(sigma_min/gamma_norm),
            distance_original_bad_to_truth=float(np.linalg.norm(refined-truth_position)),
        )

    x_bad = model.box_readout(bad_current)[1]
    start_materials = [("bad_box_readout", x_bad), ("zero", np.zeros(len(truth_j))),
                       ("box_midpoint", np.full(len(truth_j), BOX_UPPER/2))]
    start_materials += [(f"uniform_{i:02d}", rng.uniform(0, BOX_UPPER, len(truth_j)))
                        for i in range(16)]
    projector_starts = {}
    starts = []
    for label, chi0 in start_materials:
        t0, c0 = perf_counter(), process_time()
        j0 = model.forward(chi0)
        position = problem.coords(j0)
        reconstructed = problem.current(position)
        conversion_seconds = perf_counter()-t0
        conversion_cpu = process_time()-c0
        if label in ("bad_box_readout", "box_midpoint"):
            projector_starts[label] = fiber.coords(j0)
        result = dict(
            start=label, chi_initial=chi0.tolist(),
            current_conversion_distance=float(np.linalg.norm(reconstructed-j0)),
            current_conversion_port_error=float(np.linalg.norm(model.receiver@reconstructed-model.data)),
            current_conversion_seconds=conversion_seconds,
            current_conversion_cpu_seconds=conversion_cpu, current_conversion_lu=1,
            current_conversion_rhs=int(model.source.shape[1]),
        )
        for form in ("minors", "real_box"):
            u, info = problem.solve(position, form)
            result[form] = audit_case(problem, u, truth_j, form, info)
        for method in ("trf_gn", "projected_lm"):
            x, info = physical.solve(chi0, method)
            audit_start, audit_cpu_start = perf_counter(), process_time()
            j = model.forward(x)
            result[method] = dict(**model.audit(j, truth_j), solver=info, material=x.tolist(),
                                  endpoint_forward_lu=1, endpoint_forward_rhs=int(model.source.shape[1]),
                                  endpoint_with_readout_audit_seconds=perf_counter()-audit_start,
                                  endpoint_with_readout_audit_cpu_seconds=process_time()-audit_cpu_start)
        starts.append(result)
    same_bad = {}
    for form in ("minors", "real_box"):
        u, info = problem.solve(refined, form)
        same_bad[form] = audit_case(problem, u, truth_j, form, info)
    box_stationary = root(lambda z: problem.gradient(z, "real_box"),
                          np.asarray(same_bad["real_box"]["position"]),
                          jac=problem.hessian_box, options={"xtol": 1e-11})
    refined_box = audit_case(problem, box_stationary.x, truth_j, "real_box",
                             dict(root_flag=bool(box_stationary.success), message=box_stationary.message))
    refined_box["distance_from_lm_endpoint"] = float(np.linalg.norm(
        box_stationary.x-np.asarray(same_bad["real_box"]["position"])))
    physical_stationary = physical.refine_box_stationary(
        np.asarray(starts[0]["trf_gn"]["material"]), truth_j)
    projector_starts["curated_bad"] = fiber.coords(bad_current)
    direction = rng.normal(size=refined.size)
    direction /= np.linalg.norm(direction)
    mat_dir = rng.normal(size=len(truth_j))
    mat_dir /= np.linalg.norm(mat_dir)
    geom = np.linalg.norm(model.points[:,None]-model.points[None,:], axis=-1)
    np.fill_diagonal(geom, np.inf)
    min_sep = float(np.min(geom))
    reconstructed_g, reconstructed_s = disk_operators(
        model.points, source["k"], model.radius, source["receiver_angles"], 2.)
    incident = np.asarray(source["incident_angles"])
    incident_directions = np.column_stack([np.cos(incident), np.sin(incident)])
    reconstructed_b = np.exp(1j*source["k"]*(model.points@incident_directions.T))
    return dict(
        case=case, input_kind="previously exposed witness",
        operator_checks=dict(
            green_geometry_error=float(np.linalg.norm(reconstructed_g-model.green)),
            receiver_geometry_error=float(np.linalg.norm(reconstructed_s-model.receiver)),
            incident_geometry_error=float(np.linalg.norm(reconstructed_b-model.source)),
            fiber_kernel_residual=float(np.linalg.norm(model.receiver@fiber.basis)),
            fiber_basis_isometry_error=float(np.linalg.norm(
                fiber.basis.conj().T@fiber.basis-np.eye(fiber.basis.shape[1]))),
            true_forward_data_absolute_error=float(np.linalg.norm(model.receiver@truth_j-model.data)),
        ),
        cell_radius=model.radius, min_cell_center_distance=min_sep, disks_disjoint=bool(min_sep>2*model.radius),
        truth_material=model.truth.tolist(), case_seed=case_seed, source_count=int(model.source.shape[1]),
        receiver_count=int(model.receiver.shape[0]), forward_condition=float(np.linalg.cond(
            np.eye(len(truth_j))-model.truth[:, None]*model.green)),
        box=[0., BOX_UPPER], max_solver_nfev=MAX_EVALS,
        source_loss_weight=problem.weight,
        prescribed_bad=bad, box_objective_at_prescribed_bad=real_box_at_bad,
        same_bad_current_starts=same_bad, box_stationary_from_bad=refined_box,
        direct_stationary_from_legal_bad_readout=physical_stationary, legally_paired_material_starts=starts,
        projection=projected, local_truth_diagnostics=diagnostic, random_equal_retained_rank=(
            random_projection_controls(fiber, truth_j, projector_starts, operators)
            if case == "hard5" else []),
        derivative_checks=derivative_checks(problem, refined+.13*direction, direction,
                                            physical, np.full(len(truth_j), 12.0), mat_dir),
        recovery_counts={key: int(sum(row[key]["box_material_relative_error"] < SUCCESS_RTOL
                                  and row[key]["physical_data_relative"] < SUCCESS_RTOL
                                  and row[key]["current_relative_error"] < SUCCESS_RTOL
                                  for row in starts))
                         for key in ("minors", "real_box", "trf_gn", "projected_lm")},
        data_audit_cohort="same two incident waves and receivers; independent holdout absent",
        multistart_selected_by_physical_residual={
            key: {"start": row["start"], "physical_data_relative": row[key]["physical_data_relative"],
                  "material_relative_error": row[key]["box_material_relative_error"],
                  "elapsed_seconds_all_starts": sum(arm[key]["solver"]["elapsed_seconds"] for arm in starts),
                  "cpu_seconds_all_starts": sum(arm[key]["solver"]["cpu_seconds"] for arm in starts),
                  "nfev_all_starts": sum(arm[key]["solver"]["nfev"] for arm in starts),
                  "njev_all_starts": sum(arm[key]["solver"]["njev"] for arm in starts)}
            for key in ("minors", "real_box", "trf_gn", "projected_lm")
            for row in [min(starts, key=lambda r: r[key]["physical_data_relative"])]
        },
    )


def run_sweep():
    records = []
    points = np.array([[-.18, -.10], [.02, -.17], [.17, -.01], [.08, .16], [-.14, .14]])
    angles = np.array([.22, 1.67])
    k, radius = 2*np.pi, .025
    b = np.exp(1j*k*(points@np.column_stack([np.cos(angles), np.sin(angles)]).T))
    for name, directions in [("wide", [.13, 2.15, 4.52]), ("one_sided", [.12, .48, 1.04])]:
        g, s = disk_operators(points, k, radius, directions, 1.7)
        for peak in (1, 5, 10, 25, 50, 100):
            chi = peak*np.array([.7, 1., .85, .6, .9])
            a = np.eye(len(chi))-chi[:, None]*g
            j = np.linalg.solve(a, chi[:, None]*b)
            model = Model(g, b, s, s@j, chi, points, radius)
            metrics, _ = projection(make_fiber(model), j)
            # The [0,30] diagnostic box is specific to the basin experiments.
            # Audit these exact-noise-free equations without clipping the stated contrasts.
            raw = np.array(metrics["raw_material"]["real"])+1j*np.array(metrics["raw_material"]["imag"])
            records.append(dict(layout=name, peak=peak, cond_forward=float(np.linalg.cond(a)),
                                cond_projected=metrics["singular_values"][0]/metrics["singular_values"][-1],
                                rank=metrics["operator_rank"], rank_gamma=metrics["gamma_rank"],
                                physical_data_relative=float(np.linalg.norm(s@model.forward(raw)-s@j)
                                                             /np.linalg.norm(s@j)),
                                material_relative_error=metrics["raw_material_relative_error"],
                                current_relative_error=metrics["current_relative_error"],
                                max_imaginary_material=float(np.max(abs(raw.imag)))))
    return records


def save(path, values):
    path.write_text(json.dumps(values, indent=2, allow_nan=False)+"\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True, help="unmodified attachment root")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    t0, c0 = perf_counter(), process_time()
    for name in ("soft3", "hard5"):
        result = run_case(args.source/"research_work", name)
        save(args.out/(name+".json"), result)
        print(name, result["recovery_counts"], flush=True)
    save(args.out/"sweep.json", run_sweep())
    save(args.out/"replay_environment.json", dict(
        platform=platform.platform(), processor=platform.machine(), python=platform.python_version(),
        numpy=np.__version__, scipy=scipy.__version__, wall_seconds=perf_counter()-t0,
        cpu_seconds=process_time()-c0,
        threads={key: os.environ.get(key) for key in (
            "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
            "MKL_NUM_THREADS", "BLIS_NUM_THREADS")},
    ))


if __name__ == "__main__":
    main()
