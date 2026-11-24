#!/usr/bin/env python3
"""Find and certify feasible low-peak materials on the frozen R14 acquisitions."""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from warnings import catch_warnings, simplefilter

import numpy as np
from scipy.linalg import LinAlgError, LinAlgWarning, lu_factor, lu_solve
from scipy.optimize import Bounds, minimize
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
OPTICS = RESEARCH / "z_alg_optical_inertia"
TRANSFER = OPTICS / "validation/max_struct"
R14 = RESEARCH / "delegated/max_struct_r1/r14_worktree/research/research_max_struct_r14"
ARCHIVE = R14 / "RESULTS/independent24_v1"
SOURCE = R14 / "data_public/independent24"
sys.path[:0] = [str(OPTICS), str(OPTICS / "submission/diagnostic")]

from audit import dense_audit
from core.inertia import clock, timing
from flint import ctx


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path, record):
    pending = path.with_suffix(".partial")
    pending.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    pending.replace(path)


class Forward:
    def __init__(self, G, S, B, Y):
        self.G, self.S, self.E0, self.Y = G, S, B, Y
        self.n = G.shape[0]
        self.operations = Counter()

    def view(self, m):
        operator = -m[:, None] * self.G
        operator.flat[::self.n + 1] += 1.
        self.operations["physical_factorizations"] += 1
        with catch_warnings():
            simplefilter("error", LinAlgWarning)
            lu = lu_factor(operator, check_finite=False, overwrite_a=True)
        self.operations["physical_source_rhs"] += self.E0.shape[1]
        J = lu_solve(lu, m[:, None] * self.E0, check_finite=False)
        field = self.S @ J
        return field, lu, J

    def linearize(self, lu, current, residual):
        self.operations["physical_receiver_rhs"] += self.S.shape[0]
        self.operations["dense_exact_jacobians"] += 1
        receive_green = lu_solve(lu, self.S.T, trans=1, check_finite=False).T
        total = self.E0 + self.G @ current
        derivative = (receive_green[:, None, :] * total.T[None, :, :]).reshape(-1, self.n)
        jacobian = np.concatenate((derivative.real, derivative.imag), axis=0)
        target = np.concatenate((residual.real.ravel(), residual.imag.ravel()))
        return jacobian, target


def load(model_file):
    provenance = read(TRANSFER / "public_provenance.json")
    entries = provenance["observations"]
    with np.load(model_file, allow_pickle=False) as raw:
        cases = raw["cases"].tolist()
        arrays = [raw[key] for key in ("G", "S", "E0", "Y", "delta", "centers", "k")]
    G, S, B, Y, radius, centers, wave_number = arrays
    if cases != [entry["case"] for entry in entries] or Y.shape != (7, 48, 12):
        raise ValueError("The archive must have seven matched cases with all 48×12 ports.")
    for index, entry in enumerate(entries):
        with np.load(SOURCE / f"{cases[index]}.npz", allow_pickle=False) as raw:
            for key, value in (("G", G), ("S", S), ("E0", B), ("Y", Y[index]), ("delta", radius[index])):
                if not np.array_equal(raw[key], value):
                    raise ValueError(f"Original array mismatch: {cases[index]} {key}")
            if (int(raw["inverse_side"]), int(raw["forward_side"])) != (24, 24):
                raise ValueError("The data must be matched 24² disks.")
    return entries, G, S, B, Y, radius, centers, float(wave_number)


def history(model, case, floor, ceiling):
    prior = []
    rows = []
    for source in sorted(ARCHIVE.glob(f"{case}__*/primary.npz")):
        provenance = read(source.parent / "WRAPPER_PROVENANCE.json")
        if provenance.get("truth_access") != "NONE" or provenance.get("case_id") != case:
            rows.append(dict(method=source.parent.name, status="unqualified_provenance"))
            continue
        with np.load(source, allow_pickle=False) as archive:
            m = archive["m"]
        tol = 32 * np.finfo(float).eps * ceiling
        if m.shape != (model.n,) or not np.isrealobj(m) or (
                float(m.min()) < floor - tol or float(m.max()) > ceiling + tol):
            rows.append(dict(method=source.parent.name, status="outside_original_box"))
            continue
        m = np.clip(m, floor, ceiling)
        norm = float(np.linalg.norm(model.view(m)[0] - model.Y))
        receipt = read(source.parent / "RECEIPT.json")
        archived_budget = receipt.get("protocol", {}).get("CPU_cap",
                                                         receipt.get("task", {}).get("CPU_cap"))
        row = dict(method=source.parent.name, source=str(source.relative_to(RESEARCH)),
                   status="numeric_replay", peak=float(m.max()), residual=norm,
                   archived_actual_cpu_s=receipt.get("whole_process_CPU_seconds"),
                   archived_budget_cpu_s=archived_budget)
        rows.append(row)
        prior.append((row, m))
    return prior, rows


def solve_cap(model, initial, cap, floor, target, iterations):
    """Projected Gauss–Newton with boxed quadratic steps and physical backtracking."""
    m = np.clip(initial, floor, cap)
    trail = []
    damping = 3.e-5
    view = model.view(m)
    stop = "iteration_limit"
    for _ in range(iterations):
        residual = view[0] - model.Y
        resid = float(np.linalg.norm(residual))
        trail.append(resid)
        if resid <= target:
            stop = "numerically_in_ball"
            break
        jacobian, error = model.linearize(*view[1:], residual)
        gradient = jacobian.T @ error
        hessian = jacobian.T @ jacobian
        scale = float(max(np.max(hessian.diagonal()), np.finfo(float).tiny))
        gradient /= scale
        hessian /= scale
        hessian.flat[::model.n + 1] += damping

        def value_and_grad(step):
            curvature = hessian @ step
            return float(step @ (gradient + .5 * curvature)), gradient + curvature

        step = minimize(value_and_grad, np.zeros(model.n), method="L-BFGS-B", jac=True,
                        bounds=Bounds(floor - m, cap - m),
                        options=dict(maxiter=120, ftol=1.e-12, gtol=1.e-8))
        model.operations["boxed_quadratic_objective_calls"] += step.nfev
        model.operations["boxed_quadratic_iterations"] += step.nit
        accepted = False
        for fraction in (1., .5, .25, .125, .0625):
            next_m = np.clip(m + fraction * step.x, floor, cap)
            next_view = model.view(next_m)
            next_res = float(np.linalg.norm(next_view[0] - model.Y))
            if next_res < resid - max(1.e-12, resid * 1.e-8):
                m = next_m
                view = next_view
                damping = max(2.e-8, damping * .5)
                accepted = True
                break
        if not accepted:
            damping *= 10
            if damping >= .1:
                stop = "no_descent_within_box"
                break
    resid = float(np.linalg.norm(view[0] - model.Y))
    return m, view[0], dict(status=stop, cap=float(cap), residual=resid,
                            satisfies_numeric_target=resid <= target, residual_path=trail)


def run_case(entry, G, S, B, Y, radius, attempts, iterations):
    t0 = clock()
    radius = float(radius)
    case = entry["case"]
    floor = float(entry["lower"])
    ceiling = float(entry["upper"])
    model = Forward(G, S, B, Y)
    target = float(radius * (1 - 1.e-6))
    prior, history_rows = history(model, case, floor, ceiling)
    eligible = [(row, m) for row, m in prior if row["residual"] <= target]
    if eligible:
        seed, initial = min(eligible, key=lambda item: (item[0]["peak"], item[0]["residual"]))
        provenance = seed
    elif prior:
        seed, initial = min(prior, key=lambda item: item[0]["residual"])
        provenance = seed
    else:
        initial = np.full(model.n, floor)
        provenance = dict(method="constant_original_floor")
    seeded = timing(t0)
    started = clock()
    stages = []
    candidates = []
    failed_cap = None
    current = initial
    numerical = np.linalg.norm(model.view(current)[0] - Y)
    if numerical <= target:
        candidates.append((current, model.view(current)[0]))
    else:
        current, field, first = solve_cap(model, current, ceiling, floor, target, 2 * iterations)
        stages.append(first)
        if first["satisfies_numeric_target"]:
            candidates.append((current, field))
    print(case, "initial eligible", len(eligible), "peak", float(current.max()),
          "residual", float(numerical), flush=True)
    if candidates:
        for _ in range(attempts):
            high = float(current.max())
            if failed_cap is not None and high - failed_cap < .0015:
                break
            cap = floor + .86 * (high - floor) if failed_cap is None else (high + failed_cap) / 2
            build_start = clock()
            old_counts = model.operations.copy()
            proposal, field, stage = solve_cap(model, current, cap, floor, target, iterations)
            stage["operations"] = dict(model.operations - old_counts)
            stage.update(timing(build_start))
            stages.append(stage)
            if stage["satisfies_numeric_target"]:
                current = proposal
                candidates.append((proposal, field))
            else:
                failed_cap = cap
            print(case, "cap", cap, "res/r", stage["residual"] / radius,
                  "solver_time", stage["wall_seconds"], flush=True)
    optimization = timing(started)
    tried = []
    winner = None
    materials = {}
    for m, field in reversed(candidates):
        started = clock()
        if not np.all((m >= floor) & (m <= ceiling)):
            raise ValueError("Final candidate outside the original material box.")
        audit = dense_audit(SimpleNamespace(G=G, S=S, E0=B), m, Y, field)
        audit.update(dict(peak=float(m.max()), minimum=float(m.min()), radius=float(radius),
                          accepted=audit["full_wave_distance_hi"] < radius, **timing(started)))
        tried.append(audit)
        if audit["accepted"]:
            winner = audit
            materials = dict(material=m, online_endpoint_field=field, initial=initial)
            break
    if winner is None:
        print(case, "no strictly accepted candidate", flush=True)
    else:
        print(case, "accepted", winner["peak"], winner["full_wave_distance_hi"], radius, flush=True)
    return dict(case=case, stated_contract=dict(side=24, receivers=48, sources=12, lower=floor,
                          ceiling=ceiling, radius=float(radius), original_input=entry["source"]),
                truth_access="NONE", status="strictly_feasible" if winner else "unresolved",
                history_endpoints=history_rows, warm_start=provenance, history_time=seeded,
                peak_trials=stages, optimization_time=optimization,
                local_failed_cap=failed_cap, local_failure_is_global_exclusion=False,
                accepted=winner, final_audit_trials=tried, operations=dict(model.operations),
                **timing(t0)), materials


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--attempts", default=12, type=int)
    parser.add_argument("--iterations", default=16, type=int)
    parser.add_argument("--case", action="append")
    parser.add_argument("--output", default=HERE / "online", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    ctx.prec, ctx.threads = 128, 1
    t0 = clock()
    entries, G, S, B, Y, radii, _, k = load(TRANSFER / "public_inputs.npz")
    summary = dict(shared_geometry=dict(wave_number=k, model="matched saved disk-pulse scalar TM",
                   k_source="public_inputs.npz", grid_side=24,
                   noise_contract="original complex 1% RMS, CN99 per saved observation"),
                   online_setup=timing(t0), observations={})
    for index, entry in enumerate(entries):
        case = entry["case"]
        if args.case and case not in args.case:
            continue
        try:
            result, arrays = run_case(entry, G, S, B, Y[index], radii[index], args.attempts, args.iterations)
        except (ArithmeticError, LinAlgError, LinAlgWarning, ValueError) as error:
            result, arrays = dict(case=case, status="unresolved", error=str(error)), {}
        if arrays:
            np.savez_compressed(args.output / f"{case}_feasible.npz", **arrays)
        write(args.output / f"{case}.json", result)
        summary["observations"][case] = result
        write(args.output / "online.json", {**summary, **timing(t0)})


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
