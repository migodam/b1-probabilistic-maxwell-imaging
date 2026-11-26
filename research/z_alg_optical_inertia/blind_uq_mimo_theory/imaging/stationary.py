"""Box-constrained stationary single-beta CIE-I on native Fourier coordinates."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from shutil import copy2

from shared import OUT, elapsed, now, write_arrays, write_json, col_energy
from cie import BASE, CIE, CurrentState, UPPER, segment
from run import read_baseline

import numpy as np
from scipy.optimize import Bounds, minimize
from threadpoolctl import threadpool_limits

BETA = 0.5
WINDOW = 64
REQUIRED_WINDOWS = 3
MATERIAL_WINDOW_RTOL = 0.002
AUXILIARY_WINDOW_RTOL = 0.0002
PG_REDUCTION = 0.005
PG_ABSOLUTE = 2e-6
MAX_UPDATES = 6000


class Stable(Exception):
    pass


class Objective:
    def __init__(self, model, initial, beta=BETA, scales=None):
        self.model = model
        self.beta = beta
        self.cie = CIE(model)
        self.stage = len(BASE["fourier_corner_width"])-1
        self.U = self.cie.stage(self.stage)
        self.a_shape = (self.U.k, model.ns)
        self.a_len = int(np.prod(self.a_shape))
        a = self.U.transport(initial.coefficients, initial.frequencies)
        j, q = self.cie.auxiliary_current(CurrentState(initial.x, a, self.U.frequencies, self.stage))
        z_max = beta*UPPER/(1+beta*UPPER)
        self.a_scale = 1/(beta * np.sqrt(self.cie.ws))
        self.z_scale = 1/np.sqrt(np.maximum(np.sum(abs(q+beta*j)**2 * self.cie.ws, axis=1),
                                            .01/model.n))
        if scales is not None:
            self.a_scale, self.z_scale = scales
        self.lower = np.r_[np.full(2*self.a_len, -np.inf), np.zeros(model.n)]
        self.upper = np.r_[np.full(2*self.a_len, np.inf), z_max/self.z_scale]
        self.initial = self.pack(a / self.a_scale, beta*initial.x/(1+beta*initial.x)/self.z_scale)
        self.latest = None

    @staticmethod
    def pack(a, z):
        return np.r_[a.real.ravel(), a.imag.ravel(), z]

    def unpack(self, variables):
        a_re, a_im, z = np.split(variables, (self.a_len, 2*self.a_len))
        a = (a_re + 1j*a_im).reshape(self.a_shape) * self.a_scale
        z = z*self.z_scale
        x = z/(self.beta*(1-z))
        return a, z, np.clip(x, 0., UPPER)

    def evaluate(self, variables):
        model, cie = self.model, self.cie
        beta = self.beta
        a, z, x = self.unpack(variables)
        v = self.U.apply(a)
        j = cie.JD + v
        q = cie.Q0 + model.green.apply(v)
        d = cie.D0 + model.S @ v
        r = beta*j - z[:, None]*(q+beta*j)
        r_w = r*cie.ws
        reg, regg = model.reg.evaluate(x)[:2]
        f_data = .5*float(cie.wd @ col_energy(d))
        f_state = .5*float(cie.ws @ col_energy(r))
        value = f_data + f_state + reg
        g_a = self.U.adjoint(model.SH @ (d*cie.wd) + beta*(1-z[:, None])*r_w -
                             model.green.adjoint(z[:, None]*r_w))
        g_z = -np.sum(((q+beta*j).conj()*r_w).real, axis=1) + regg/(beta*(1-z)**2)
        grad = self.pack(g_a*self.a_scale, g_z*self.z_scale)
        projected = variables - np.clip(variables-grad, self.lower, self.upper)
        self.latest = (variables.copy(), value, grad, x,
                       dict(data_term=f_data, state_term=f_state, material_prior=float(reg),
                            projected_inf=float(np.max(np.abs(projected))),
                            projected_l2=float(np.linalg.norm(projected))))
        model.cost["stationary_evaluations"] += 1
        return value, grad

    def directional_check(self):
        # Only interior material coordinates move in this deterministic check.
        # Its direction is generated from public state, with no held data.
        rng = np.random.default_rng(1421)
        z = self.initial
        d = rng.standard_normal(z.size)
        d[(z < self.lower+1e-4) | (z > self.upper-1e-4)] = 0.
        d /= np.linalg.norm(d)
        _, g = self.evaluate(z)
        eps = 4e-6
        p, _ = self.evaluate(z+eps*d)
        m, _ = self.evaluate(z-eps*d)
        numerical = (p-m)/(2*eps)
        analytical = float(g @ d)
        absolute_defect = abs(numerical-analytical)
        if absolute_defect > 1e-7*(1+abs(analytical)):
            raise ArithmeticError(f"Single-beta CIE objective gradient check failed: {absolute_defect}")
        return dict(approximate=numerical, analytical=analytical, defect=absolute_defect)


def run_stationary(case, extend=False):
    started = now()
    case_path, model, _, _ = read_baseline(case)
    stage = len(BASE["fourier_corner_width"])-1
    path = case_path / "stationary"
    path.mkdir(exist_ok=True)
    previous = None
    scales = None
    prefix_path = path if extend else case_path
    if extend:
        previous = json.loads((path/"baseline.json").read_text())
        if previous["stable"]:
            return
        suffix = f"cap_{previous['updates']:05d}"
        copy2(path/"baseline.json",path/(suffix+".json"))
        copy2(path/"baseline.npz",path/(suffix+".npz"))
    with np.load(prefix_path/"baseline.npz", allow_pickle=False) as data:
        x = data["chi"] if extend else data[f"chi_{stage}"]
        initial = CurrentState(x.ravel(), data[f"coefficients_{stage}"],
                               data[f"frequencies_{stage}"], stage)
        if extend:
            scales = (data["coefficient_scale"], data["material_scale"])
    objective = Objective(model, initial, scales=scales)
    gradient_check = objective.directional_check()
    initial_value = objective.evaluate(objective.initial)[0]
    initial_gradient = (objective.latest[-1]["projected_inf"] if previous is None
                        else previous["initial_projected_inf"])
    projected_ceiling = max(PG_ABSOLUTE, PG_REDUCTION*initial_gradient)
    stable_windows = 0
    trace, windows = [], []
    last_x, last_value = initial.x.copy(), initial_value
    final_v = objective.initial.copy()

    def callback(intermediate_result):
        nonlocal stable_windows, last_x, last_value, final_v
        final_v = intermediate_result.x.copy()
        if objective.latest is None or not np.array_equal(objective.latest[0], final_v):
            objective.evaluate(final_v)
        _, f, _, x, info = objective.latest
        update = len(trace) + 1
        trace.append(dict(update=update, auxiliary=float(f), **info, **elapsed(started)))
        if update % WINDOW:
            return
        difference = float(np.linalg.norm(x-last_x)/max(np.linalg.norm(last_x),1e-300))
        decrease = (last_value-f)/max(abs(last_value),1e-300)
        is_stable = (difference < MATERIAL_WINDOW_RTOL and abs(decrease) < AUXILIARY_WINDOW_RTOL
                     and info["projected_inf"] < projected_ceiling)
        stable_windows = stable_windows + 1 if is_stable else 0
        last_x, last_value = x.copy(), f
        windows.append(dict(update=update, material_relative_change=difference, auxiliary_relative_decrease=decrease,
                            stable_windows=stable_windows, projected_inf=info["projected_inf"]))
        if update % (5*WINDOW) == 0 or is_stable:
            print(case, "single-beta", update, "aux", round(f, 8), "pgrad", info["projected_inf"],
                  "window_dx", difference, "window_df", decrease, "stable", stable_windows, flush=True)
        if stable_windows >= REQUIRED_WINDOWS:
            raise Stable("CONSECUTIVE_STATIONARY_WINDOWS")

    try:
        ans = minimize(objective.evaluate, objective.initial, method="L-BFGS-B", jac=True,
                       bounds=Bounds(objective.lower, objective.upper), callback=callback,
                       options=dict(maxiter=MAX_UPDATES, maxfun=2*MAX_UPDATES, maxcor=15,
                                    ftol=1e-15, gtol=1e-10, maxls=30))
        final_v, status = ans.x, str(ans.message)
    except Stable as stop:
        status = str(stop)
    final_a, _, final_x = objective.unpack(final_v)
    final_aux = objective.evaluate(final_v)[0]
    stable = stable_windows >= REQUIRED_WINDOWS
    physics_start = model.audit(initial.x)
    physics_end = model.audit(final_x)
    masks, bg, smooth, roi_info = segment(final_x.reshape(model.shape))
    state = model.forward(final_x)
    previous_trace = [] if previous is None else previous["trace"]
    previous_windows = [] if previous is None else previous["all_windows"]
    start_ordinal = len(previous_trace)
    for row in [*trace,*windows]:
        row["update"] += start_ordinal
    total = elapsed(started)
    counters = Counter(model.cost)
    if previous is not None:
        summary_initial_value = previous["initial_auxiliary"]
        physics_start = previous["physics_start"]
        for key in total:
            total[key] += previous[key]
        for key, old in previous["counts"].items():
            counters[key] = max(old, counters[key]) if key.startswith("max_") else old+counters[key]
    else:
        summary_initial_value = initial_value
    summary = dict(case=case, side=model.side, beta=BETA,
                   objective="Single beta CIE-I reduced Fourier current and shared real material with common R8 prior",
                   stable=stable, status=status, updates=start_ordinal+len(trace),
                   initial_auxiliary=summary_initial_value,
                   final_auxiliary=final_aux, final=objective.latest[-1], initial_projected_inf=initial_gradient,
                   gradient_check=gradient_check, physics_start=physics_start, original_physics=physics_end,
                   segmentation=roi_info, all_windows=previous_windows+windows,
                   trace=previous_trace+trace, counts=dict(counters),
                   extension_counts=dict(model.cost) if previous else None,
                   initial_budget_stop=previous["status"] if previous else None,
                   extension_start_auxiliary=initial_value if previous else None,
                   criteria=dict(window=WINDOW, required_windows=REQUIRED_WINDOWS,
                                 material_window_relative_change=MATERIAL_WINDOW_RTOL,
                                 auxiliary_window_relative_decrease=AUXILIARY_WINDOW_RTOL,
                                 projected_gradient_ceiling=projected_ceiling,
                                 max_updates=MAX_UPDATES), **total)
    write_json(path/"baseline.json", summary)
    write_arrays(path/"baseline.npz", chi=final_x.reshape(model.shape), roi=masks[0], rois=masks,
                 background_tiles=bg, smooth=smooth, current=state["j"], train_prediction=state["pred"],
                 stage=stage, **{f"coefficients_{stage}":final_a,f"frequencies_{stage}":objective.U.frequencies},
                 material_scale=objective.z_scale, coefficient_scale=objective.a_scale, beta=BETA)
    print(case, "single-beta", status, "updates", start_ordinal+len(trace), "train/r",
          physics_end["residual_over_radius"], "stable", stable, flush=True)


if __name__ == "__main__":
    parse = argparse.ArgumentParser()
    parse.add_argument("cases", nargs="+")
    parse.add_argument("--extend", action="store_true")
    args = parse.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        for case in args.cases:
            run_stationary(case, extend=args.extend)
