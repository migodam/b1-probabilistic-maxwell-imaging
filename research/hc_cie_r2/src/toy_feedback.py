#!/usr/bin/env python3
"""Floating-point retained-tail replay, alongside the exact toy certificates."""
import argparse
from itertools import product
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.optimize import least_squares, minimize

from toy_models import pack, real_columns, unpack
from toy_replay import save


TAU_MAX = 15/7
FEEDBACK_VECTOR = np.array([[1., .8, .6]])
INCIDENT = np.array([[1., .4+.7j], [.6+.1j, 1.], [.3+.8j, -.5+.6j]])
RECEIVER = np.array([[1., 0., .35], [0., 1., -.2]])
AUGMENTED_RECEIVER = np.vstack([RECEIVER, FEEDBACK_VECTOR])
TRUE_CONTRAST = np.array([3., 8., 13.])
FEEDBACK_OPERATOR = (.15+.1j)*FEEDBACK_VECTOR.T


def toy(epsilon):
    tail = 1j*epsilon
    g = FEEDBACK_OPERATOR@FEEDBACK_VECTOR+(-.4+tail)*np.eye(3)

    def physical(chi):
        return np.linalg.solve(np.eye(3)-chi[:, None]*g, chi[:, None]*INCIDENT)

    j_true = physical(TRUE_CONTRAST)
    data = RECEIVER@j_true
    q_true = FEEDBACK_VECTOR@j_true
    return tail, g, physical, data, q_true


def inner(tau, q, t, data):
    # In this particular toy T=tI. These elementwise formulas keep its full value.
    drive = INCIDENT+FEEDBACK_OPERATOR@q
    inverse = 1/(1-t*tau)
    current = (tau*inverse)[:, None]*drive
    residual = AUGMENTED_RECEIVER@current-np.vstack([data, q])
    dj = AUGMENTED_RECEIVER.T[:, :, None]*(inverse**2)[:, None, None]*drive[:, None, :]
    d2j = AUGMENTED_RECEIVER.T[:, :, None]*(2*t*inverse**3)[:, None, None]*drive[:, None, :]
    jac = real_columns(dj.reshape(3, -1).T)
    hessian = jac.T@jac+np.diag((d2j.reshape(3,-1)@residual.ravel().conj()).real)
    return residual, jac, hessian, current


def run_tail(epsilon):
    t0, c0 = perf_counter(), process_time()
    t, g, physical, data, q_true = toy(epsilon)
    calls = inner_evals = inner_failures = 0

    def solve_inner(q):
        nonlocal calls, inner_evals, inner_failures
        calls += 1

        def fg(tau):
            f, jac, *_ = inner(tau, q, t, data)
            r = pack(f)
            return .5*r@r, jac.T@r

        sol = minimize(fg, np.full(3, TAU_MAX/2), jac=True, method="L-BFGS-B",
                       bounds=[(0, TAU_MAX)]*3, options=dict(ftol=1e-15, gtol=1e-10,
                                                            maxiter=160, maxls=30))
        inner_evals += sol.nfev
        inner_failures += int(not sol.success)
        return sol.x

    def reduced(q_real):
        q = unpack(q_real, (1,2))
        tau = solve_inner(q)
        return pack(inner(tau, q, t, data)[0])

    scale = max(.2, np.linalg.norm(data)/3)
    starts = [np.zeros(4)] + [scale*np.array(v) for v in product((-1, 1), repeat=4)]
    multistart = []
    for index, start in enumerate(starts):
        sol = least_squares(reduced, start, xtol=2e-10, ftol=2e-10, gtol=2e-10, max_nfev=120)
        q = unpack(sol.x, (1,2))
        tau = solve_inner(q)
        chi = tau/(1-.4*tau)
        multistart.append(dict(
            index=index, initial_q_realified=start.tolist(), material=chi.tolist(),
            optimizer_status=int(sol.status), outer_reported_nfev=int(sol.nfev),
            augmented_residual=float(np.linalg.norm(inner(tau, q, t, data)[0])),
            physical_data_relative=float(np.linalg.norm(RECEIVER@physical(chi)-data)
                                         / np.linalg.norm(data)),
            material_relative=float(np.linalg.norm(chi-TRUE_CONTRAST)/np.linalg.norm(TRUE_CONTRAST)),
        ))

    alpha = TAU_MAX*epsilon
    drive = INCIDENT+FEEDBACK_OPERATOR@q_true
    a_f, a_inf = np.linalg.norm(drive), np.linalg.norm(drive, axis=1).max()
    obs_norm = np.linalg.norm(AUGMENTED_RECEIVER, 2)
    linear_columns = (AUGMENTED_RECEIVER.T[:, :, None]*drive[:, None, :]).reshape(3,-1).T
    gamma = np.linalg.svd(real_columns(linear_columns), compute_uv=False)[-1]
    perturb = obs_norm*(alpha*a_inf/(1-alpha)+alpha*a_f/(1-alpha)**2)
    c2 = 2*obs_norm*epsilon*a_f/(1-alpha)**3
    rho = obs_norm*TAU_MAX*a_f/(1-alpha)+np.linalg.norm(np.vstack([data, q_true]))
    rng = np.random.default_rng(47)
    sampled = [np.linalg.eigvalsh(inner(TAU_MAX*rng.random(3), q_true, t, data)[2])[0]
               for _ in range(400)]
    tau = np.array([2.1, 2.1, .05])
    h = np.array([-.5, .1, .86])
    f, jac, hs, j = inner(tau, q_true, t, data)
    dh = 1e-5*h
    fp, jp, *_ = inner(tau+dh, q_true, t, data)
    fm, jm, *_ = inner(tau-dh, q_true, t, data)
    hess_err = np.linalg.norm((jp.T@pack(fp)-jm.T@pack(fm))/(2e-5)-hs@h)
    # Evaluate a genuinely unclosed auxiliary state before comparing to full scattering.
    q_probe = q_true+.11*np.array([[.3+.7j, -.2+.5j]])
    j_aux = inner(tau, q_probe, t, data)[-1]
    closure = FEEDBACK_VECTOR@j_aux-q_probe
    correction = np.linalg.solve(
        np.eye(3)-tau[:, None]*(t*np.eye(3)+FEEDBACK_OPERATOR@FEEDBACK_VECTOR),
        tau[:, None]*(FEEDBACK_OPERATOR@closure))
    chi_probe = tau/(1-.4*tau)
    physical_difference = physical(chi_probe)-j_aux
    return dict(
        epsilon=epsilon, model_rank=int(np.linalg.matrix_rank(g)), alpha=float(alpha),
        gamma0=float(gamma), error_bound=float(perturb), c2_bound=float(c2), rho_bound=float(rho),
        curvature_sufficient_bound=float(max(gamma-perturb,0)**2-rho*c2),
        positive_linear_margin=bool(gamma>perturb),
        sampled_min_hessian=float(min(sampled)), directional_witness=float(h@hs@h),
        q_scope="Frozen q_true for certificate and derivatives; the 17 optimization starts use data only.",
        derivative_checks=dict(
            jacobian_error=float(np.linalg.norm(pack(fp-fm)/(2e-5)-jac@h)),
            hessian_error=float(hess_err),
            finite_closure_response_error=float(np.linalg.norm(physical_difference-correction)),
            probe_closure_norm=float(np.linalg.norm(closure)),
            probe_physical_difference_norm=float(np.linalg.norm(physical_difference)),
        ),
        success_count=sum(r["material_relative"]<1e-5 for r in multistart),
        inner_calls=calls, inner_evaluations=inner_evals, inner_optimizer_failures=inner_failures,
        multistart=multistart, elapsed_seconds=perf_counter()-t0, cpu_seconds=process_time()-c0,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    results = []
    for epsilon in (.0001, .001, .01, .05, .1, .3):
        result = run_tail(epsilon)
        print(epsilon, result["sampled_min_hessian"], result["success_count"], flush=True)
        results.append(result)
    save(args.out/"feedback_sweep.json", results)


if __name__ == "__main__":
    main()
