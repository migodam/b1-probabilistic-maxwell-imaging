"""Single CUDA context, immutable shared operators, FIFO tasks and stream lanes.

Run as a module or directly: ``python src/executor.py run manifest.json --lanes 2``.
The main thread writes receipts while lane threads immediately take their next task.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import importlib
import json
import os
from pathlib import Path
import queue
import time
import traceback

for _variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_variable] = "1"


class BudgetStop(RuntimeError):
    pass


def _json_safe(value):
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if hasattr(value, "tolist"):
        return _json_safe(value.tolist())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float):
        import math
        return value if math.isfinite(value) else None
    return value


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(_json_safe(value), ensure_ascii=False, indent=2), encoding="utf-8")
    temporary.replace(path)


def _load_npz(path):
    import numpy as np
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def _resolve(root, path):
    path = Path(path)
    return path if path.is_absolute() else root / path


def solve_task(lane, task):
    """Optimize one lane. No process-wide synchronization or I/O in callbacks."""
    import numpy as np
    from scipy.optimize import least_squares, minimize

    cfg = task.get("config", {})
    method = {"data_lbfgs": "data", "globaldual": "dual"}.get(task["method"], task["method"])
    stage = int(task.get("stage", cfg.get("stage", 24)))
    n = int(getattr(lane, "n_coeffs", stage * stage))
    lower, upper = float(cfg.get("lower", .02)), float(cfg.get("upper", 20.))
    start = time.perf_counter()
    wall_limit = float(cfg.get("wall_seconds", float("inf")))
    x0 = np.broadcast_to(np.asarray(cfg.get("initial_coefficients", cfg.get("initial_material", .1)), float), (n,)).copy()
    x0 = np.clip(x0, lower, upper)
    best_coeffs, best_value = x0.copy(), float("inf")
    trace = []

    def guard():
        if time.perf_counter() - start >= wall_limit:
            raise BudgetStop("WALL_SECONDS")

    def evaluate(coeffs):
        nonlocal best_coeffs, best_value
        guard()
        value, gradient = lane.evaluate(coeffs)
        value = float(value)
        physical = float(lane.last_physical_objective)
        if physical < best_value:
            best_coeffs, best_value = np.asarray(coeffs).copy(), physical
        trace.append({"evaluation": len(trace) + 1, "seconds": time.perf_counter() - start,
                      "objective": value, "physical_objective": physical})
        return value, np.asarray(gradient, float)

    status, reason, fit = "COMPLETED_NUMERICAL", None, None
    try:
        if method in ("data", "dual"):
            width = upper - lower

            def scaled(z):
                value, gradient = evaluate(lower + width * z)
                return value, width * gradient

            fit = minimize(scaled, (x0 - lower) / width, jac=True, method="L-BFGS-B",
                           bounds=[(0., 1.)] * n,
                           options={"maxiter": int(cfg.get("maxiter", 2000)),
                                    "maxfun": int(cfg.get("maxfun", 4000)),
                                    "ftol": float(cfg.get("ftol", 1e-12)),
                                    "gtol": float(cfg.get("gtol", 1e-8)),
                                    "maxls": int(cfg.get("maxls", 25))})
        elif method in ("native_trf", "reciprocal_trf"):
            if not (hasattr(lane, "residual") and hasattr(lane, "jacobian")):
                raise NotImplementedError("TRF requires ObjectiveLane.residual(coeffs) and .jacobian(coeffs); Full-J interface unavailable")
            reciprocal = method == "reciprocal_trf"
            cached_x, cached_r = None, None

            def residual(x):
                nonlocal best_coeffs, best_value, cached_x, cached_r
                guard()
                coeffs = 1. / x if reciprocal else x
                cached_r = np.asarray(lane.residual(coeffs), float)
                cached_x = np.asarray(x).copy()
                value = .5 * float(cached_r @ cached_r)
                if value < best_value:
                    best_coeffs, best_value = np.asarray(coeffs).copy(), value
                trace.append({"evaluation": len(trace) + 1, "seconds": time.perf_counter() - start, "objective": value})
                return cached_r

            def jacobian(x):
                guard()
                if cached_x is None or not np.array_equal(x, cached_x):
                    residual(x)
                coeffs = 1. / x if reciprocal else x
                jac = np.asarray(lane.jacobian(coeffs), float)
                return jac * (-coeffs ** 2)[None, :] if reciprocal else jac

            fit = least_squares(residual, 1. / x0 if reciprocal else x0, jac=jacobian,
                                method="trf", bounds=(1. / upper, 1. / lower) if reciprocal else (lower, upper),
                                x_scale="jac", tr_solver="exact", max_nfev=int(cfg.get("maxfun", 4000)),
                                ftol=float(cfg.get("trf_ftol", 1e-10)), xtol=float(cfg.get("trf_xtol", 1e-10)),
                                gtol=float(cfg.get("trf_gtol", 1e-8)))
        else:
            raise ValueError("Unknown method " + str(method))
    except BudgetStop as exc:
        status, reason = "BUDGET_STOP", str(exc)

    result = lane.close(best_coeffs)
    result.update({"wrapper_method": method, "optimizer_status": None if fit is None else int(fit.status),
                   "optimizer_success": None if fit is None else bool(fit.success),
                   "optimizer_message": None if fit is None else str(fit.message),
                   "optimizer_nfev": None if fit is None else int(fit.nfev),
                   "wrapper_status": status, "wrapper_stop_reason": reason,
                   "optimizer_wall_s": time.perf_counter() - start,
                   "best_physical_objective": best_value, "evaluation_trace": trace,
                   "_coefficients": best_coeffs,
                   "selection": "minimum evaluated original data objective plus common regularizer; no GT",
                   "stage": stage, "truth_access": "NONE_PUBLIC_INPUT_ONLY"})
    return result


def _resample(coeffs, old_stage, new_stage):
    """Endpoint-aligned bilinear interpolation, matching the frozen basis rule."""
    import numpy as np
    source = np.asarray(coeffs).reshape(old_stage, old_stage)
    old_axis, new_axis = np.linspace(0., 1., old_stage), np.linspace(0., 1., new_stage)
    horizontal = np.asarray([np.interp(new_axis, old_axis, row) for row in source])
    return np.asarray([np.interp(new_axis, old_axis, row) for row in horizontal.T]).T.ravel()


def solve_multiscale(shared, public, task, device_module, shared_ready):
    import torch
    cfg = dict(task.get("config", {}))
    stages = list(map(int, cfg.get("stages", [task.get("stage", cfg.get("stage", 24))])))
    method = task["method"]
    started = time.perf_counter()
    total_wall = float(cfg.get("wall_seconds", float("inf")))
    coefficients, previous_stage, stage_results = None, None, []
    for stage in stages:
        remaining = total_wall - (time.perf_counter() - started)
        if remaining <= 0 and stage_results:
            break
        stage_cfg = {**cfg, "stage": stage, "wall_seconds": max(0., remaining)}
        if coefficients is not None:
            stage_cfg["initial_coefficients"] = _resample(coefficients, previous_stage, stage)
        objective = "dual" if method == "dual" or (method == "globaldual" and stage == stages[-1]) else "data"
        phase_method = "data" if method == "globaldual" and stage != stages[-1] else method
        if method == "globaldual" and stage == stages[-1]:
            # Fixed-work allocation permits a throughput comparison without a
            # concurrency-dependent change to the numerical experiment budget.
            fraction = float(cfg.get("dual_fraction", .45))
            stage_cfg["maxiter"] = max(1, int(cfg.get("maxiter", 100) * fraction))
            stage_cfg["maxfun"] = max(1, int(cfg.get("maxfun", 200) * fraction))
            stage_cfg["wall_seconds"] = max(0., remaining * fraction)
        lane = device_module.ObjectiveLane(shared, public, stage=stage, objective=objective, config=stage_cfg)
        lane.stream.wait_event(shared_ready)
        with torch.cuda.stream(lane.stream):
            result = solve_task(lane, {**task, "method": phase_method, "stage": stage, "config": stage_cfg})
        coefficients = result.pop("_coefficients")
        previous_stage = stage
        stage_results.append(result)
        del lane
        if method == "globaldual" and stage == stages[-1]:
            remaining = total_wall - (time.perf_counter() - started)
            tail_cfg = {**cfg, "stage": stage, "wall_seconds": max(0., remaining),
                        "initial_coefficients": coefficients,
                        "maxiter": max(1, int(cfg.get("maxiter", 100) * (1. - fraction))),
                        "maxfun": max(1, int(cfg.get("maxfun", 200) * (1. - fraction)))}
            tail = device_module.ObjectiveLane(shared, public, stage=stage, objective="data", config=tail_cfg)
            with torch.cuda.stream(tail.stream):
                corrected = solve_task(tail, {**task, "method": "data", "stage": stage, "config": tail_cfg})
            coefficients = corrected.pop("_coefficients")
            corrected["globaldual_proposal_phase"] = result
            stage_results[-1] = corrected
            result = corrected
            del tail
        if result["wrapper_status"] == "BUDGET_STOP":
            break
    from collections import Counter
    phases = []
    for row in stage_results:
        if "globaldual_proposal_phase" in row:
            phases.append(row["globaldual_proposal_phase"])
        phases.append(row)
    selected = min(phases, key=lambda row: row["ranking_objective"])
    final = dict(selected)
    all_counts = Counter()
    for row in phases:
        all_counts.update(row["counts"])
    final.update({"configured_stages": stages, "completed_stages": [r["stage"] for r in stage_results],
                  "stage_results": stage_results, "task_optimize_wall_s": time.perf_counter() - started,
                  "counts": dict(all_counts),
                  "cumulative_optimizer_evaluations": sum(len(r["evaluation_trace"]) for r in phases),
                  "all_material_candidates": [{k: row[k] for k in ("m", "ranking_objective", "relative_residual", "stage")} for row in phases],
                  "maxiter_maxfun_scope": "per stage; shared task wall limit",
                  "all_stages_reached": len(stage_results) == len(stages)})
    if len(stage_results) != len(stages):
        final.update({"wrapper_status": "BUDGET_STOP", "wrapper_stop_reason": "WALL_SECONDS"})
    return final


def run(manifest_path, lanes=None):
    """Execute the finite FIFO manifest in one process, with no batch barrier."""
    import torch
    import scipy.optimize  # Load BLAS libraries before applying one global limit.
    manifest_path = Path(manifest_path).resolve()
    root = manifest_path.parent
    manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    if "vendor_path" in manifest:
        import sys
        sys.path.insert(0, str(_resolve(root, manifest["vendor_path"])))
    from threadpoolctl import threadpool_limits
    lane_count = int(lanes if lanes is not None else manifest.get("lanes", 1))
    if lane_count not in (1, 2, 4):
        raise ValueError("lanes must be 1, 2, or 4")
    output = _resolve(root, manifest.get("output_dir", "results"))
    output.mkdir(parents=True, exist_ok=True)
    module_name = manifest.get("device_module", f"{__package__}.device" if __package__ else "device")
    device_module = importlib.import_module(module_name)
    torch.cuda.set_device(int(manifest.get("device", 0)))
    torch.set_num_threads(1)
    started = time.perf_counter()
    shared_public = _load_npz(_resolve(root, manifest["operators_path"]))
    shared = device_module.SharedOperators(shared_public)
    shared_ready = torch.cuda.Event()
    shared_ready.record(torch.cuda.current_stream())
    # The finite observation set is small; preload it once while GPU setup runs.
    # This keeps disk reads out of the task handoff and avoids four reloads per Y.
    public_cache = {task["public_path"]: None for task in manifest["tasks"]}
    for path in public_cache:
        public_cache[path] = _load_npz(_resolve(root, path))
    init_seconds = time.perf_counter() - started
    hybrid = None
    if manifest.get("svd_backend") == "hybrid":
        module_prefix = device_module.__name__.rsplit(".", 1)[0]
        HybridSVD = importlib.import_module(module_prefix + ".hybrid_svd").HybridSVD
        hybrid = HybridSVD()
        hybrid.install()
    incoming, outgoing = queue.Queue(), queue.Queue()
    for index, task in enumerate(manifest["tasks"]):
        incoming.put((index, task, time.perf_counter()))

    def worker(slot):
        torch.cuda.set_device(int(manifest.get("device", 0)))
        try:
            while True:
                try:
                    index, task, enqueued = incoming.get_nowait()
                except queue.Empty:
                    return
                began = time.perf_counter()
                ident = task["task_id"]
                outgoing.put({"event": "started", "slot": slot, "task_id": ident, "sequence": index,
                              "seconds": began - started, "queue_wait_s": began - enqueued})
                try:
                    public = public_cache[task["public_path"]]
                    result = solve_multiscale(shared, public, task, device_module, shared_ready)
                    # close() returns host-ready JSON; releasing a task never waits on other lanes.
                    receipt = {"event": "result", "task_id": ident, "sequence": index, "slot": slot,
                               "status": "OK", "result": result,
                               "optimizer_evaluations": result["cumulative_optimizer_evaluations"]}
                except Exception as exc:
                    receipt = {"event": "result", "task_id": ident, "sequence": index, "slot": slot,
                               "status": "UNSUPPORTED_METHOD" if isinstance(exc, NotImplementedError) else "FAIL",
                               "error": str(exc), "traceback": traceback.format_exc(),
                               "costs_status": "INCOMPLETE_EXCEPTION", "costs": None,
                               "optimizer_evaluations": 0}
                receipt.update({"task": task, "worker_wall_s": time.perf_counter() - began,
                                "queue_wait_s": began - enqueued, "finished_seconds": time.perf_counter() - started,
                                "pid": os.getpid()})
                outgoing.put(receipt)
                incoming.task_done()
        finally:
            outgoing.put({"event": "lane_stopped", "slot": slot})

    receipts, events, stopped = [], [], 0
    # BLAS limits apply once around all workers; nested per-task contexts would race.
    with threadpool_limits(limits=1), ThreadPoolExecutor(max_workers=lane_count, thread_name_prefix="cuda-lane") as pool:
        futures = [pool.submit(worker, slot) for slot in range(lane_count)]
        while stopped < lane_count:
            event = outgoing.get()
            if event["event"] == "lane_stopped":
                stopped += 1
            elif event["event"] == "started":
                events.append(event)
            else:
                task = event["task"]
                result_path = _resolve(root, task["result_path"]) if "result_path" in task else output / (task["task_id"] + ".json")
                write_json(result_path, event)
                receipts.append({key: event[key] for key in ("task_id", "sequence", "slot", "status", "worker_wall_s", "queue_wait_s", "finished_seconds", "optimizer_evaluations")})
        for future in futures:
            future.result()
    if hybrid is not None:
        hybrid.restore()
    summary = {"schema": "cuda-pipeline-executor-v1", "utc": datetime.now(timezone.utc).isoformat(),
               "pid": os.getpid(), "lanes": lane_count, "expected": len(manifest["tasks"]),
               "completed": len(receipts), "shared_operator_imports": 1,
               "shared_init_wall_s": init_seconds, "wall_s": time.perf_counter() - started,
               "observations_preloaded": len(public_cache), "shared_setup_counts": shared.setup_counts,
               "svd_backend": "cpu" if hybrid is None else hybrid.summary(),
               "status_counts": {status: sum(r["status"] == status for r in receipts) for status in sorted({r["status"] for r in receipts})},
               "receipts": receipts, "starts": events,
               "cumulative_optimizer_evaluations": sum(r.get("optimizer_evaluations", 0) for r in receipts),
               "queue_policy": "FIFO; each lane immediately takes next task; main thread writes receipts",
               "synchronization": "shared-ready event; lane stream context; host returns synchronize only their lane",
               "budget_note": "shared task wall guards run between callbacks; an active CUDA solve cannot be interrupted; physical costs are recorded by device lanes"}
    write_json(output / "runner_summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("run")
    command.add_argument("manifest")
    command.add_argument("--lanes", type=int, choices=(1, 2, 4))
    args = parser.parse_args()
    summary = run(args.manifest, args.lanes)
    print(json.dumps({key: summary[key] for key in ("completed", "expected", "lanes", "wall_s", "status_counts")}, indent=2))


if __name__ == "__main__":
    main()
