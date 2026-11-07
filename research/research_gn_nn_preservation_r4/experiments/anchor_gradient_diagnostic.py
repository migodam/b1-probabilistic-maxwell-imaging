"""Offline endpoint gradients from public anchors; no optimization/Full-J/GT.

Prepared for a parent-coordinated single-CPU run. Default scope is TRAIN/VAL.
The original anchor/solver receipts and stopping statuses are never rewritten.
"""
from __future__ import annotations

import os
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from threadpoolctl import threadpool_limits, threadpool_info

STUDY = Path(__file__).resolve().parents[1]
LABEL = "OFFLINE_RECOMPUTED_NOT_SOLVER_REPORTED"
FIELDS = (
    "obs_id", "parent_id", "split", "noise", "diagnostic_origin", "status", "error",
    "original_solver_status", "original_stop_reason", "original_final_stage",
    "normalized_data_objective", "relative_data_residual", "gradient_L2", "gradient_Linf",
    "unit_step_box_gradient_mapping_L2", "unit_step_box_gradient_mapping_Linf",
    "box_Frank_Wolfe_gap", "lower_active_cells", "upper_active_cells", "active_box_fraction",
    "active_tolerance", "LU", "actual_complex_RHS", "adjoint_complex_RHS", "forward_complex_RHS",
    "io_wall_seconds", "LU_wall_seconds", "adjoint_wall_seconds", "wall_seconds", "cpu_seconds",
    "anchor_npz", "anchor_receipt",
)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def endpoint_gradient(operators, record, anchor_root, row):
    """The parent-verified data-objective adjoint, with one LU and 12 RHS."""
    n, ns = 576, 12
    path = anchor_root/f"{record['obs_id']}.npz"
    sidecar = path.with_suffix(".json")
    row.update(anchor_npz=str(path), anchor_receipt=str(sidecar))
    tick = time.perf_counter()
    with np.load(path, allow_pickle=False) as archive:
        m = np.asarray(archive["m0"], np.float64).reshape(n)
        j = np.asarray(archive["j0"], np.complex128).reshape(n, ns)
        q = np.asarray(archive["Q0"], np.complex128).reshape(n, ns)
        y = np.asarray(archive["Y"], np.complex128).reshape(48, ns)
    original = json.loads(sidecar.read_text())
    row.update(original_solver_status=original.get("status"),
               original_stop_reason=original.get("stop_reason"),
               original_final_stage=original.get("final_stage"))
    lower = float(original.get("config", {}).get("lower", .02))
    upper = float(original.get("config", {}).get("upper", 20.))
    row["io_wall_seconds"] = time.perf_counter()-tick
    scale2 = float(np.vdot(y, y).real)
    if scale2 <= 0:
        raise ValueError("Nonzero public Y is required for the declared normalized objective")
    residual = operators["S"]@j-y
    objective = .5*float(np.vdot(residual, residual).real)/scale2
    rhs = operators["S"].conj().T@residual/scale2
    matrix = np.eye(n, dtype=np.complex128)-m[:, None]*operators["G"]
    tick = time.perf_counter()
    row["LU"] += 1
    try:
        factor = lu_factor(matrix, check_finite=False)
    finally:
        row["LU_wall_seconds"] = time.perf_counter()-tick
    tick = time.perf_counter()
    row["actual_complex_RHS"] += ns
    row["adjoint_complex_RHS"] += ns
    try:
        adjoint = lu_solve(factor, rhs, trans=2, check_finite=False)
    finally:
        row["adjoint_wall_seconds"] = time.perf_counter()-tick
    gradient = np.real(np.sum(adjoint.conj()*q, axis=1))
    if not np.isfinite(gradient).all():
        raise FloatingPointError("Nonfinite recomputed endpoint gradient")
    mapping = m-np.clip(m-gradient, lower, upper)
    vertex = np.where(gradient >= 0, lower, upper)
    gap = float(np.dot(gradient, m-vertex))
    tolerance = 1e-8
    lower_count = int(np.count_nonzero(m <= lower+tolerance))
    upper_count = int(np.count_nonzero(m >= upper-tolerance))
    row.update(normalized_data_objective=objective,
               relative_data_residual=float(np.sqrt(2*objective)),
               gradient_L2=float(np.linalg.norm(gradient)),
               gradient_Linf=float(np.linalg.norm(gradient, np.inf)),
               unit_step_box_gradient_mapping_L2=float(np.linalg.norm(mapping)),
               unit_step_box_gradient_mapping_Linf=float(np.linalg.norm(mapping, np.inf)),
               box_Frank_Wolfe_gap=gap, lower_active_cells=lower_count,
               upper_active_cells=upper_count, active_box_fraction=(lower_count+upper_count)/n,
               active_tolerance=tolerance, status="COMPLETED")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--splits", nargs="+", choices=("train", "val", "test"), default=["train", "val"])
    parser.add_argument("--anchor-root", type=Path, default=STUDY/"anchors")
    parser.add_argument("--output", type=Path, help="New output directory; default is a timestamped audit subdirectory")
    parser.add_argument("--limit", type=int, help="Optional finite prefix of the public online manifest; labelled in the receipt")
    args = parser.parse_args()
    if "test" in args.splits and not (STUDY/"audit/METHOD_FREEZE.json").exists():
        raise RuntimeError("TEST endpoint diagnostics require METHOD_FREEZE; no TEST truth is ever used")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    records = json.loads((STUDY/"data/manifests/online.json").read_text())["records"]
    records = [r for r in records if r["split"] in args.splits]
    if args.limit is not None:
        records = records[:args.limit]
    if not records:
        raise ValueError("No public observations match this finite plan")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    output = args.output or STUDY/"audit/anchor_gradient_diagnostic"/stamp
    output.mkdir(parents=True, exist_ok=False)
    started, cpu_started = time.perf_counter(), time.process_time()
    setup_start = time.perf_counter()
    with np.load(STUDY/"data/operators.npz", allow_pickle=False) as archive:
        operators = {name: np.asarray(archive[name], np.complex128) for name in ("G", "S")}
    setup_seconds = time.perf_counter()-setup_start
    rows, interrupted = [], False
    receipt_path = output/"RECEIPT.json"

    def receipt(status):
        counts = Counter()
        for row in rows:
            counts.update({key: row[key] for key in ("LU", "actual_complex_RHS", "adjoint_complex_RHS", "forward_complex_RHS")})
        return dict(schema="gn-nn-r4-offline-anchor-gradient-v1", status=status,
                    diagnostic_origin=LABEL, start_utc=stamp, splits=args.splits,
                    planned_observations=len(records), attempted_observations=len(rows),
                    planned_parents=len({r["parent_id"] for r in records}), limit=args.limit,
                    statuses=dict(Counter(r["status"] for r in rows)), costcounts=dict(counts),
                    operation_counter_semantics="Actual dispatched LU/RHS calls, including failed/interrupted calls; no assumed successful completion of a failed solve",
                    wall_seconds=time.perf_counter()-started, cpu_seconds=time.process_time()-cpu_started,
                    operator_setup_seconds=setup_seconds,
                    sum_observation_wall_seconds=sum(r["wall_seconds"] for r in rows),
                    BLAS_threads=1, CPU_workers=1, threadpools=threadpool_info(),
                    objective="0.5 ||S j0 - Y||_F^2 / ||Y||_F^2; original data objective only",
                    derivative="L=I-diag(m0)G; lambda=L^{-*}S*(Sj0-Y)/||Y||^2; g=Re sum_s conj(lambda) Q0",
                    gradient_mapping="m0 - clip(m0 - g, .02, 20), unit step in original real material coordinates",
                    Frank_Wolfe_gap="g dot (m0-v), v_i=lower if g_i>=0 else upper",
                    activity_tolerance=1e-8, physics_dtype="complex128", forward_RHS=0,
                    truth_access="NONE; online manifest, operators, public anchor m0/j0/Q0/Y and scalar original receipt only",
                    optimization_performed=False, Full_J_constructed=False,
                    formula_or_FD_validation_performed=False,
                    original_anchor_and_solver_receipts_modified=False,
                    interpretation="Offline diagnostic only; does not replace solver-reported optimality or alter original stop status; no GT recovery or causal overfit conclusion",
                    csv=str(output/"PER_OBSERVATION.csv"))

    write_json(receipt_path, receipt("RUNNING"))
    try:
        with (output/"PER_OBSERVATION.csv").open("w", newline="") as handle, threadpool_limits(limits=1):
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            for index, record in enumerate(records):
                tick, cpu_tick = time.perf_counter(), time.process_time()
                row = {key: record[key] for key in ("obs_id", "parent_id", "split", "noise")}
                row.update(diagnostic_origin=LABEL, status="RUNNING", error=None, LU=0,
                           actual_complex_RHS=0, adjoint_complex_RHS=0, forward_complex_RHS=0,
                           io_wall_seconds=0., LU_wall_seconds=0., adjoint_wall_seconds=0.)
                try:
                    endpoint_gradient(operators, record, args.anchor_root, row)
                except KeyboardInterrupt:
                    row.update(status="INTERRUPTED", error="KeyboardInterrupt; dispatched cost retained")
                    interrupted = True
                except Exception as exc:
                    row.update(status="FAILED", error=f"{type(exc).__name__}: {exc}")
                row.update(wall_seconds=time.perf_counter()-tick, cpu_seconds=time.process_time()-cpu_tick)
                rows.append(row);writer.writerow(row);handle.flush()
                if interrupted:
                    break
                if (index+1) % 20 == 0 or index+1 == len(records):
                    write_json(receipt_path, receipt("RUNNING"))
                    print(json.dumps(dict(completed_plan_entries=index+1, planned=len(records), statuses=dict(Counter(r["status"] for r in rows)))), flush=True)
    finally:
        final = "INTERRUPTED" if interrupted or len(rows)<len(records) else ("COMPLETED_WITH_FAILURES" if any(r["status"]!="COMPLETED" for r in rows) else "COMPLETED")
        write_json(receipt_path, receipt(final))
    print(json.dumps(dict(status=final, output=str(output), diagnostic_origin=LABEL)), flush=True)


if __name__ == "__main__":
    main()
