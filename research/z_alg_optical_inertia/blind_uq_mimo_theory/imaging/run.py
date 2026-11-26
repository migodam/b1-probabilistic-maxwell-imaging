"""Blind baseline, matched algorithm arms, and separate held-port scoring."""
from __future__ import annotations

import argparse
from pathlib import Path
import json
import sys
from time import perf_counter

from shared import (CASE_DIR, FREEZE, OPERATOR_DIR, OUT, Model, elapsed, energy,
                    now, write_arrays, write_json)
from cie import CIE, CurrentState, Penalty, UPPER, compute_k_modes, run_baseline, segment

import numpy as np
from scipy.optimize import Bounds, minimize
from threadpoolctl import threadpool_limits


def get_case_path(case):
    result = OUT / case
    result.mkdir(parents=True, exist_ok=True)
    return result


def baseline(case):
    path = get_case_path(case)
    started = now()
    model = Model.read(case)
    setup = elapsed(started)
    x, states, report = run_baseline(model)
    masks, tiles, smooth, roi_info = segment(x.reshape(model.shape))
    baseline_audit = model.audit(x)
    st = model.forward(x)
    raw = dict(chi=x.reshape(model.shape), roi=masks[0], rois=masks,
               background_tiles=tiles, smooth=smooth, current=st["j"],
               train_prediction=st["pred"], stage=int(report["selected"]["stage"]))
    for state in states:
        raw[f"chi_{state.stage}"] = state.x.reshape(model.shape)
        raw[f"coefficients_{state.stage}"] = state.coefficients
        raw[f"frequencies_{state.stage}"] = state.frequencies
    report.update(case=case, side=model.side, operator_id=model.operator_id,
                  radius=model.radius, segmentation=roi_info, original_physics=baseline_audit,
                  setup=setup, end_to_end=elapsed(started), end_to_end_counts=dict(model.cost))
    write_json(path / "baseline.json", report)
    write_arrays(path / "baseline.npz", **raw)
    print(case, "baseline", [r["updates"] for r in report["stages"]],
          [r["stop"] for r in report["stages"]], "roi", [int(m.sum()) for m in masks],
          "train/r", baseline_audit["residual_over_radius"], flush=True)


def read_baseline(case, stationary=False):
    path = get_case_path(case)
    beta = None
    if stationary:
        path /= "stationary"
        report = json.loads((path / "baseline.json").read_text())
        if not report["stable"]:
            raise RuntimeError(f"{case}: stationary baseline did not qualify")
        beta = float(report["beta"])
    model = Model.read(case)
    with np.load(path / "baseline.npz", allow_pickle=False) as data:
        stage = int(data["stage"])
        if stage < 0:
            a = np.empty((0, model.ns), complex)
            freqs = np.empty((0, 2), int)
        else:
            a = data[f"coefficients_{stage}"]
            freqs = data[f"frequencies_{stage}"]
        start = CurrentState(data["chi"].ravel(), a, freqs, max(stage, 0), beta)
        rois = data["rois"].reshape(-1, model.n)
    return path, model, start, rois


def continuation(case):
    path, model, _, _ = read_baseline(case)
    started = now()
    stage = len(FREEZE["baseline"]["fourier_corner_width"]) - 1
    with np.load(path / "baseline.npz", allow_pickle=False) as data:
        initial = CurrentState(data[f"chi_{stage}"].ravel(), data[f"coefficients_{stage}"],
                               data[f"frequencies_{stage}"], stage)
    original = model.audit(initial.x)
    cie = CIE(model)
    tighter = .5 * FREEZE["baseline"]["relative_coefficient_rms_stop_strictly_below"]
    final, report, audited = cie.advance(initial, stage,
                                          FREEZE["baseline"]["maximum_updates_per_stage"],
                                          stop=True, coefficient_rms_tol=tighter)
    original_aux = json.loads((path / "baseline.json").read_text())["stages"][-1]["trace"][-1]
    ending_aux = report["trace"][-1]
    report.update(case=case, original=original, tightened=audited[-1][1],
                  original_aux_objective=original_aux["data_term"] + original_aux["state_term"],
                  tightened_aux_objective=ending_aux["data_term"] + ending_aux["state_term"],
                  relative_material_change=float(np.linalg.norm(final.x-initial.x)/
                                                 max(np.linalg.norm(initial.x),1e-300)),
                  end_to_end=elapsed(started), total_counts=dict(model.cost))
    write_arrays(path / "continuation.npz", chi=final.x.reshape(model.shape),
                 coefficients=final.coefficients, frequencies=final.frequencies)
    write_json(path / "continuation.json", report)
    print(case, "tighter", report["stop"], report["updates"],
          [original["residual_over_radius"], report["tightened"]["residual_over_radius"]],
          "relative dx", report["relative_material_change"], flush=True)


class FitStop(Exception):
    pass


def physical_fit(model, seed, upper, *, max_updates, target=None, wall_budget=None):
    """Shared real-material nonlinear objective and exact-solve adjoint.

    The discrepancy target is an early exit available only to existence
    search. The control spends its paid budget reducing the same data loss.
    """
    model.reset()
    started = now()
    trial = np.minimum(seed, upper)
    best = trial.copy()
    best_loss = np.inf
    records, accepted = [], []

    def calculate(x):
        nonlocal best, best_loss
        if records and wall_budget is not None and perf_counter() - started[0] >= wall_budget:
            raise FitStop("WALL_BUDGET")
        value, grad = model.data_grad(x)
        state = model.forward(x)
        records.append(dict(loss=value, residual_over_radius=state["residual_ratio"],
                            state_relative_defect=state["state_relative_defect"],
                            gradient_inf=float(np.max(np.abs(grad))), **elapsed(started)))
        if value < best_loss:
            best, best_loss = x.copy(), value
        if target is not None and state["residual_ratio"] <= target:
            raise FitStop("DISCREPANCY_TARGET")
        return value, grad

    def callback(x):
        accepted.append(dict(evaluation=len(records), **elapsed(started)))

    try:
        opt = minimize(calculate, trial, method="L-BFGS-B", jac=True,
                       bounds=Bounds(np.zeros_like(upper), upper), callback=callback,
                       options=dict(maxiter=max_updates, maxls=25, ftol=1e-13, gtol=1e-9))
        status = str(opt.message)
    except FitStop as stop:
        status = str(stop)
    state = model.forward(best)
    info = model.audit(best)
    info.update(stop=status, accepted_updates=len(accepted), attempted_evaluations=len(records),
                trace=records, accepted_trace=accepted, counts=dict(model.cost), **elapsed(started))
    return best, state["j"], info


def paired_cie(model, seed, updates, *, upper=None, penalty=None):
    model.reset()
    started = now()
    c = CIE(model)
    _, info, candidates = c.advance(seed, len(FREEZE["baseline"]["fourier_corner_width"])-1,
                                    updates, box_upper=upper, penalty=penalty,
                                    audit_stride=FREEZE["B2"]["audit_every_updates"])
    info.update(global_cost=dict(model.cost), cie_setup=c.setup, end_to_end=elapsed(started))
    return info, candidates


def b1(case, stationary=False):
    path, model, initial, rois = read_baseline(case, stationary)
    started = now()
    cfg = FREEZE["B1"]
    results = {}
    arrays = {}
    feasible = []
    for j, region in enumerate(rois):
        tag = f"roi_{j+1:02d}"
        upper = np.where(region, 0., UPPER)
        x, current, out = physical_fit(model, initial.x, upper,
                                       max_updates=cfg["witness_updates_per_roi"],
                                       target=cfg["witness_target_radius_ratio"])
        out.update(fixed_cells=int(np.count_nonzero(region)), max_fixed_contrast=float(np.max(x[region])),
                   roi_contrast_mass=float(initial.x[region].sum()),
                   proposed_feasible=bool(out["residual_over_radius"] < cfg["checked_q_preferred_margin"]))
        results[tag] = out
        arrays[f"x_{tag}"] = x.reshape(model.shape)
        arrays[f"j_{tag}"] = current
        if out["proposed_feasible"]:
            feasible.append((out["roi_contrast_mass"], out["residual_over_radius"], j))
        print(case, "B1 witness", "stationary" if stationary else "dual", tag,
              "train/r", out["residual_over_radius"], "cost", out["wall_s"], flush=True)
        write_json(path / "b1_fit_progress.json", dict(sealed=False, fits=results, **elapsed(started)))
    fits = elapsed(started)
    output = initial.x.copy()
    final_info = model.audit(output, source="baseline_fallback")
    chosen = None
    constrained = None
    if feasible:
        _, _, j = min(feasible)
        chosen = f"roi_{j+1:02d}"
        upper = np.where(rois[j], 0., UPPER)
        witness = arrays[f"x_{chosen}"].ravel()
        witness_audit = model.audit(witness, source="full_wave_witness")
        witness_seed = CurrentState(witness, initial.coefficients, initial.frequencies,
                                    initial.stage, initial.common_beta)
        constrained, candidates = paired_cie(model, witness_seed, cfg["cie_review_updates"], upper=upper)
        eligible = [(x, dict(a, source="constrained_cie")) for x, a in candidates
                    if a["residual_over_radius"] < cfg["checked_q_preferred_margin"]]
        eligible.append((witness, witness_audit))
        output, final_info = min(eligible,key=lambda v: v[1]["physical_objective"])
    paid = elapsed(started)
    final_fields = model.forward(output)["j"]
    write_arrays(path / "b1_predictions.npz", chi=output.reshape(model.shape), j=final_fields, **arrays)
    # The saved witness material/current precede all comparison and held-port reading.
    u_info, u_candidates = paired_cie(model, initial, cfg["cie_review_updates"])
    u_trial, u_selected = min(u_candidates + [(initial.x, model.audit(initial.x, step=0))],
                              key=lambda a: a[1]["physical_objective"])
    full_x, full_j, full = physical_fit(model, initial.x, np.full(model.n, UPPER),
                                       max_updates=2000, wall_budget=paid["wall_s"])
    write_arrays(path / "b1_controls.npz", cie_chi=u_trial.reshape(model.shape),
                 physical_chi=full_x.reshape(model.shape), physical_j=full_j)
    write_json(path / "b1.json", dict(case=case, side=model.side, selected_roi=chosen,
                                     selected=final_info, existence="numerical original operator",
                                     witness_fits=results, witness_fits_timing=fits,
                                     constrained_cie=constrained, compared_paid_budget=paid,
                                     matched_cie=u_info, matched_cie_selected=u_selected,
                                     matched_wall_full_wave=full, **elapsed(started)))
    write_json(path / "b1_fit_progress.json", dict(sealed=True, fits=results, **elapsed(started)))
    print(case, "B1", "stationary" if stationary else "dual", chosen, "train/r", final_info["residual_over_radius"],
          "wall control/r", full["residual_over_radius"], flush=True)


def b2(case, stationary=False):
    path, model, initial, _ = read_baseline(case, stationary)
    started = now()
    cfg = FREEZE["B2"]
    dark, random, modes = compute_k_modes(model)
    c = CIE(model)
    aux_j, aux_q = c.auxiliary_current(initial)
    beta = initial.common_beta if stationary else FREEZE["baseline"]["beta_current"][-1]
    p = beta / (1 + beta * initial.x)
    aux_r = p[:, None] * (aux_j - initial.x[:, None]*aux_q)
    target = .1 * float(c.ws @ np.sum(abs(aux_r)**2,axis=0))
    modes["initial_target_quadratic_energy_before_half"] = target
    incumbent = model.audit(initial.x, step=0)
    masks = (("plain", "noop"), ("ridge", None), ("random", random), ("dark", dark))
    reports, arrays, audited = {}, {}, {}
    for name, basis in masks:
        if isinstance(basis, str):
            penalty = None
        else:
            q_unscaled = Penalty(1., basis)
            scale = min(1., target/max(float(c.ws @ q_unscaled.quadratic(aux_j)), 1e-300))
            penalty = Penalty(scale, basis)
        report, candidates = paired_cie(model, initial, cfg["cie_updates"], penalty=penalty)
        x, candidate = min(candidates + [(initial.x, incumbent)], key=lambda a: a[1]["physical_objective"])
        arrays[name] = x.reshape(model.shape)
        arrays[name + "_last"] = candidates[-1][0].reshape(model.shape)
        st = model.forward(x)
        report["selected"] = candidate
        report["actual_physical_k_negative_energy_fraction"] = (
            energy(dark.conj().T @ st["j"])/max(energy(st["j"]), 1e-300))
        report["penalty_initial_energy_before_half"] = (
            None if penalty is None else float(c.ws @ penalty.quadratic(aux_j)))
        report["penalty_weight"] = None if penalty is None else penalty.weight
        report["with_selected_replay_cost"] = dict(model.cost)
        reports[name] = report
        audited[name] = candidates
    mode_prepaid = {name: modes["wall_s"] if name in ("dark","random") else 0. for name in reports}
    common_wall = min(report["end_to_end"]["wall_s"] + mode_prepaid[name]
                      for name, report in reports.items())
    for name, report in reports.items():
        prefix = [(x, audit) for x, audit in audited[name]
                  if audit["wall_s"] + report["cie_setup"]["wall_s"] + mode_prepaid[name] <= common_wall]
        x, selected = min(prefix + [(initial.x, incumbent)], key=lambda a: a[1]["physical_objective"])
        arrays[name+"_time"] = x.reshape(model.shape)
        report["matched_wall_selected"] = selected
        report["spectral_setup_charged_wall_s"] = mode_prepaid[name]
    physical_initial = model.forward(initial.x)["j"]
    modes["physical_initial_energy_fraction"] = energy(dark.conj().T @ physical_initial)/max(energy(physical_initial), 1e-300)
    write_arrays(path / "b2_predictions.npz", dark_modes=dark, random_null_modes=random, **arrays)
    write_json(path / "b2.json", dict(case=case, side=model.side, modes=modes, common_extra_wall_s=common_wall,
                                     arms=reports, **elapsed(started)))
    print(case, "B2", "stationary" if stationary else "dual", {k: round(v["selected"]["residual_over_radius"], 4)
          for k,v in reports.items()}, flush=True)


def held(case, stationary=False):
    """Evaluate immutable prediction arrays on observed held rows, after sealing."""
    path, model, initial, rois = read_baseline(case, stationary)
    with np.load(CASE_DIR / f"{case}.npz", allow_pickle=False) as data:
        Y, radius = data["Y_hold"], float(data["radius_hold"])
    with np.load(OPERATOR_DIR / f"{model.operator_id}.npz", allow_pickle=False) as data:
        S = data["S_hold"]
    started = now()
    matrices = {"baseline": initial.x}
    if (path / "continuation.npz").exists():
        with np.load(path / "continuation.npz", allow_pickle=False) as data:
            matrices["tighter_dual"] = data["chi"].ravel()
    if (path / "b1_predictions.npz").exists():
        with np.load(path / "b1_predictions.npz", allow_pickle=False) as data:
            matrices["b1"] = data["chi"].ravel()
            matrices.update({k: data[k].ravel() for k in data.files if k.startswith("x_roi_")})
        with np.load(path / "b1_controls.npz", allow_pickle=False) as data:
            matrices["b1_cie_control"] = data["cie_chi"].ravel()
            matrices["b1_wall_control"] = data["physical_chi"].ravel()
    if (path / "b2_predictions.npz").exists():
        with np.load(path / "b2_predictions.npz", allow_pickle=False) as data:
            matrices.update({"b2_"+k: data[k].ravel() for k in data.files
                             if k not in ("dark_modes","random_null_modes")})
    predictions, report, checked = {}, {}, []
    for name, x in matrices.items():
        duplicate = next((old_name for old_name, old_x in checked if np.array_equal(x, old_x)), None)
        if duplicate is not None:
            predictions[name] = predictions[duplicate].copy()
            report[name] = dict(report[duplicate], reused_exact_material=duplicate)
            continue
        state = model.forward(x, fresh=True)
        v = S @ state["j"]
        predictions[name] = v
        residual = np.linalg.norm(v-Y)
        report[name] = dict(held_relative_residual=float(residual/np.linalg.norm(Y)),
                            held_residual_over_radius=float(residual/radius),
                            **model.audit(x))
        checked.append((name, x))
    write_arrays(path / "held_predictions.npz", **predictions)
    write_json(path / "held.json", dict(case=case, predictions=report, counts=dict(model.cost),
                                       **elapsed(started)))
    print(case, "held", {k: round(v["held_relative_residual"], 4)
          for k,v in report.items()}, flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=("baseline","continuation","b1","b2","held"))
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--stationary", action="store_true")
    cli = ap.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        for oid in cli.cases:
            if cli.phase in ("b1", "b2", "held"):
                globals()[cli.phase](oid, stationary=cli.stationary)
            else:
                globals()[cli.phase](oid)
