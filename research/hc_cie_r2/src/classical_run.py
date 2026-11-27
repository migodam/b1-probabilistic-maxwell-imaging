"""Matched original-physics continuation versus an exact at-most-two-pixel jump."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

from imaging_run import REPO, json_write

import numpy as np
from scipy.linalg import lu_solve
from scipy.optimize import minimize


TOP = Path(os.environ.get("HC_CIE_R2_TOP_ROOT", REPO / "research/top_venue_killtest_20261120"))
RUNTIME = TOP / ".runs"


def inherited():
    if "hc_top_venue" in sys.modules:
        return sys.modules["hc_top_venue"]
    spec = importlib.util.spec_from_file_location("hc_top_venue", TOP / "run_killtest.py")
    source = importlib.util.module_from_spec(spec)
    os.environ.setdefault("MAXWELL_KILLTEST_WORKSPACE", str(REPO))
    sys.modules[spec.name] = source
    spec.loader.exec_module(source)
    return source


def load_start(case, method, cfg):
    context = RUNTIME / "public" / f"{case}.npz"
    past = RUNTIME / "pred" / f"{case}.npz"
    key = method + cfg["parent_suffix"]
    with np.load(context, allow_pickle=False) as z:
        obs = {name: z[name] for name in z.files}
    with np.load(past, allow_pickle=False) as z:
        m = z[key]
        saved = z[f"{key}_Y"]
        rows = json.loads(str(z["rows_json"]))
    receipt = next(r for r in rows if r["key"] == key)
    return obs, json.loads(str(obs["case_json"])), m, saved, receipt


def finite_two(model, m, physical, Y, config):
    _, grad = model.data_value_gradient(m, Y)
    lo, hi = config["material_box"]
    pg = m - np.clip(m-grad, lo, hi)
    indices = np.argsort(-np.abs(pg), kind="stable")[:2]
    directions = np.sign(-pg[indices])
    rhs = np.zeros((len(m), len(indices)))
    rhs[indices, np.arange(len(indices))] = 1.
    if physical["lu"] is None:
        dressed_current = rhs
    else:
        dressed_current = lu_solve(physical["lu"], rhs, check_finite=False)
        model.clock.counts["complex_RHS"] += rhs.shape[1]
        model.clock.counts["intervention_columns_RHS"] += rhs.shape[1]
    local_green = model.G[indices, :] @ dressed_current
    local_ports = model.ports.S @ dressed_current
    field = physical["Q"][indices]
    norm2 = np.vdot(Y, Y).real
    results = []
    attempted = set()
    for amplitude in config["amplitudes"]:
        for chosen in [[1, 0], [0, 1], [1, 1]]:
            change = np.clip(m[indices]+directions*amplitude*chosen, lo, hi)-m[indices]
            signature = tuple(change)
            if signature in attempted or not np.any(change):
                continue
            attempted.add(signature)
            w = np.linalg.solve(np.eye(len(indices))-change[:, None]*local_green, change[:, None]*field)
            model.clock.counts["local_resolvent_LU"] += 1
            model.clock.counts["local_resolvent_RHS"] += Y.shape[1]
            predicted = local_ports @ w
            if model.ports.mask is not None:
                predicted = np.take_along_axis(predicted, model.ports.mask.T, axis=0)
            y = physical["Y"] + predicted
            r = y-Y
            trial = m.copy()
            trial[indices] += change
            phi = .5*np.vdot(r, r).real/norm2 + model.reg.evaluate(trial)[0]
            results.append((float(phi), trial, y, dict(indices=indices.tolist(), changes=change.tolist())))
    return results


def run_one(args):
    case, parent, arm, cfg = args
    start_cpu, start_wall = time.process_time(), time.perf_counter()
    out = Path(cfg["output"]) / f"{case}__{parent}__{arm}"
    out.mkdir(parents=True, exist_ok=False)
    top = inherited()
    public, meta, initial, saved_y, receipt = load_start(case, parent, cfg)
    model = top.PaperModel(public, meta, top.Clock(np.inf))
    measured = "rx_ids" in public
    cap = cfg["budget_cpu_s"]["measured" if measured else "synthetic"]
    lo, hi = cfg["material_box"]
    mats, ports, chain = [], [], []
    Y = public["Y"]
    norm = np.linalg.norm(Y)
    best = None
    candidates = []
    phase = None

    class BudgetEnd(Exception):
        """An operation cannot begin outside the paid CPU prefix."""

    def check():
        if time.process_time()-start_cpu >= cap:
            raise BudgetEnd

    def record(m, state, stage):
        nonlocal best
        cpu = time.process_time()-start_cpu
        residual = float(np.linalg.norm(state["Y"]-Y)/norm)
        d = .5*residual**2
        r = float(model.reg.evaluate(m)[0])
        if np.min(m) < lo-1e-12 or np.max(m) > hi+1e-12 or not np.isfinite(m).all():
            raise ValueError("Invalid audited two-pixel material or continuation point")
        index = len(chain)
        chain.append(dict(index=index, phase=stage, cpu=cpu, wall=time.perf_counter()-start_wall,
                          phi=d+r, D=d, R=r, relative_port_error=residual, paid=cpu<=cap,
                          counts=dict(model.clock.counts)))
        mats.append(m.copy())
        ports.append(state["Y"].copy())
        if chain[-1]["paid"] and (best is None or chain[best]["phi"] > d+r):
            best = index
        return chain[-1]

    base = model.state(initial)
    archived_error = float(np.linalg.norm(base["Y"]-saved_y)/max(np.linalg.norm(saved_y), 1e-100))
    if archived_error > 1e-7:
        raise ValueError(f"Saved observation and current original operator disagree: {archived_error:.3g}")
    record(initial, base, "endpoint")
    formula_check = None
    try:
        if arm == "finite_two_pixels":
            check()
            screened = finite_two(model, initial, base, Y, cfg)
            for phi, _, _, action in screened:
                candidates.append(dict(predicted_phi=phi, **action))
            if screened:
                chosen = min(screened, key=lambda a: a[0])
                check()
                full = model.state(chosen[1])
                row = record(chosen[1], full, "jump")
                formula_check = dict(forecast_phi=chosen[0], full_phi=row["phi"],
                                     port_relative_error=float(np.linalg.norm(full["Y"]-chosen[2])/
                                                               max(np.linalg.norm(full["Y"]), 1e-100)),
                                     **chosen[3])
            phase = "FINITE_JUMP_DONE"
        else:
            objective_scale = cfg.get("continuation_scale", 1e6)
            def objective(m):
                check()
                val, grad = model.data_value_gradient(m, Y)
                record(m, model.state(m), "continuation")
                check()
                return objective_scale*val, objective_scale*grad

            result = minimize(objective, initial, method="L-BFGS-B", jac=True,
                              bounds=[(lo, hi)]*len(initial), options=dict(
                                  maxiter=6000, maxfun=24000, ftol=1e-12,
                                  gtol=1e-8*objective_scale, maxls=cfg.get("continuation_maxls", 30)))
            phase = str(result.message)
    except BudgetEnd:
        phase = "PROCESS_CPU_LIMIT"
    if best is None:
        raise RuntimeError("Original endpoint audit did not finish within the declared prefix.")
    calc_cpu = time.process_time()-start_cpu
    np.savez_compressed(out/"AUDITS.npz", material=np.asarray(mats), selected_index=np.asarray(best),
                        selected_port=ports[best], endpoint_port=ports[0])
    report = dict(case=case, parent=parent, arm=arm, best=best, without_truth_or_holdout=True,
                  budget_cpu=cap, online_calculation_cpu=calc_cpu, total_cpu=time.process_time()-start_cpu,
                  total_wall=time.perf_counter()-start_wall,
                  stop=phase, full_counts=dict(model.clock.counts), chain=chain,
                  historical_receipt=receipt, parent_port_relative_discrepancy=archived_error,
                  formula_check=formula_check, local_candidates=candidates)
    json_write(out/"TRACE.json", report)
    return dict(case=case, parent=parent, arm=arm, directory=out.name, selected=best,
                original=chain[0]["relative_port_error"], final=chain[best]["relative_port_error"],
                candidate_count=len(chain))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, choices=(1, 2), default=2)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    output = args.out.absolute()
    output.mkdir(parents=True, exist_ok=False)
    cfg["output"] = str(output)
    started = time.perf_counter()
    tasks = [(case, method, arm, cfg) for case in cfg["observation_ids"] for method in cfg["parent_methods"]
             for arm in cfg["arms"]]
    rows = []
    with ProcessPoolExecutor(args.workers) as ex:
        for row in ex.map(run_one, tasks):
            rows.append(row)
            print(row, flush=True)
    json_write(output/"PREDICTION_SEAL.json", dict(id=cfg["id"], expected=len(tasks), jobs=rows,
               workers=args.workers, batch_wall=time.perf_counter()-started, no_truth_or_holdout_read=True))


if __name__ == "__main__":
    main()
