"""Finite public-only CPU anchor queue, one fresh subprocess per observation.

BLAS=1; at most two workers on macOS and four on Windows. No persistent service.
Public online manifest required, so evaluator-only metadata never enters workers.
"""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import multiprocessing as mp
import os
from pathlib import Path
import subprocess
import sys
import time

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"
STUDY = Path(__file__).resolve().parents[1]
REPOSITORY = STUDY.parents[1]
sys.path.insert(0, str(REPOSITORY))


def launch_job(task):
    command, output = task
    tick = time.perf_counter()
    finished = subprocess.run(command, capture_output=True, text=True, env=os.environ.copy())
    return dict(output_path=str(output), exit_code=finished.returncode,
                child_process_wall_seconds=time.perf_counter()-tick,
                stdout=finished.stdout[-1000:], stderr=finished.stderr[-3000:])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=STUDY/"data")
    parser.add_argument("--output", type=Path, default=STUDY/"anchors")
    parser.add_argument("--reference", type=Path, default=STUDY/"cache/reference_jacobian.npz")
    parser.add_argument("--split", action="append", choices=["train", "val", "test", "dev"])
    parser.add_argument("--obs-id", action="append")
    parser.add_argument("--noise", type=float, choices=[0., 1e-3], help="Optional paired-observation filter")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--weakreg", type=float, default=0.)
    parser.add_argument("--worker-obs", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--worker-output", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--core-mask", type=int, help="Optional Windows affinity mask; affects this worker only")
    args = parser.parse_args()
    from research.research_gn_nn_completion_r3.src.anchors import generate_anchor, prepare_reference, dump
    if args.worker_obs is not None:
        if args.core_mask and sys.platform == "win32":
            import ctypes
            ctypes.windll.kernel32.SetProcessAffinityMask(ctypes.windll.kernel32.GetCurrentProcess(), args.core_mask)
        receipt = generate_anchor(args.data/"operators.npz", args.worker_obs, args.worker_output, args.reference, args.weakreg)
        print(json.dumps({key: receipt[key] for key in ("status", "solve_wall_seconds", "feature_wall_seconds", "final_stage", "final_relative_residual", "counts")}))
        return
    if args.workers < 1 or args.workers > (4 if sys.platform == "win32" else 2):
        parser.error("At most 4 Windows or 2 local CPU workers")
    if not args.split and not args.obs_id:
        parser.error("Select --split or --obs-id explicitly; default never launches the full dataset")
    tick = time.perf_counter()
    shared_setup = prepare_reference(args.data/"operators.npz", args.reference)
    manifest = json.loads((args.data/"manifests/online.json").read_text(encoding="utf-8"))
    allowed = {"obs_id", "parent_id", "split", "noise", "public_path"}
    records = []
    for record in manifest["records"]:
        if set(record) != allowed:
            raise ValueError("Online manifest contains unexpected metadata")
        if args.split and record["split"] not in args.split: continue
        if args.obs_id and record["obs_id"] not in args.obs_id: continue
        if args.noise is not None and record["noise"] != args.noise: continue
        records.append(record)
    if args.limit: records = records[:args.limit]
    suffix = "_weakreg1e-5" if args.weakreg else ""
    tasks, skipped = [], []
    for record in records:
        output = args.output/(record["obs_id"]+suffix+".npz")
        if output.exists() and output.with_suffix(".json").exists():
            skipped.append(record["obs_id"])
            continue
        command = [sys.executable, str(Path(__file__).resolve()), "--data", str(args.data.resolve()),
                   "--reference", str(args.reference.resolve()), "--worker-obs", str((args.data/record["public_path"]).resolve()),
                   "--worker-output", str(output.resolve()), "--weakreg", str(args.weakreg)]
        if args.core_mask:
            command.extend(["--core-mask", str(args.core_mask)])
        tasks.append((command, output))
    completions = []
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("spawn")) as pool:
        futures = [pool.submit(launch_job, task) for task in tasks]
        for future in as_completed(futures):
            result = future.result()
            completions.append(result)
            print(json.dumps(result), flush=True)
    costs, stops = Counter(), Counter()
    solve_seconds = feature_seconds = 0.
    for result in completions:
        receipt_path = Path(result["output_path"]).with_suffix(".json")
        if result["exit_code"] == 0 and receipt_path.exists():
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            costs.update(receipt["counts"])
            stops[receipt["stop_reason"] or receipt["status"]] += 1
            solve_seconds += receipt["solve_wall_seconds"]
            feature_seconds += receipt["feature_wall_seconds"]
        else:
            stops["FAILED_COSTS_UNKNOWN"] += 1
    batch = dict(status="COMPLETED" if all(r["exit_code"] == 0 for r in completions) else "COMPLETED_WITH_FAILURES",
                 workers=args.workers, selected=len(records), launched=len(tasks), skipped_existing=skipped,
                 wall_seconds=time.perf_counter()-tick, shared_reference_setup=shared_setup,
                 aggregate_model_counts=dict(costs), stop_histogram=dict(stops),
                 sum_solver_wall_seconds=solve_seconds, sum_feature_wall_seconds=feature_seconds,
                 jobs=completions, online_truth_access="NONE", weakreg=args.weakreg)
    batch_id = "_".join(args.split or ["selected"])+suffix+"_"+str(time.time_ns())
    dump(args.output/("BATCH_"+batch_id+".json"), batch)
    print(json.dumps({key: batch[key] for key in ("status", "selected", "launched", "wall_seconds")}))


if __name__ == "__main__":
    main()
