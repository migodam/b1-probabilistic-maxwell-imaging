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
    try: worker_result=json.loads(finished.stdout.strip().splitlines()[-1])
    except (ValueError,IndexError): worker_result={}
    return dict(worker_action=worker_result.get("worker_action","UNKNOWN_LEGACY"), output_path=str(output), exit_code=finished.returncode,
                child_process_wall_seconds=time.perf_counter()-tick,
                stdout=finished.stdout[-1000:], stderr=finished.stderr[-3000:])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=STUDY/"data")
    parser.add_argument("--output", type=Path, default=STUDY/"anchors")
    parser.add_argument("--reference", type=Path, default=STUDY/"data/reference_jacobian.npz")
    parser.add_argument("--split", action="append", choices=["train", "val", "test", "dev"])
    parser.add_argument("--obs-id", action="append")
    parser.add_argument("--id-file", type=Path, help="JSON list of public observation IDs for a finite disjoint batch")
    parser.add_argument("--noise", type=float, choices=[0., 1e-3], help="Optional paired-observation filter")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--discrepancy", action="store_true", help="Separate discrepancy baseline, never primary NN anchor")
    parser.add_argument("--discrepancy-all-scales", action="store_true", help="Separate standard discrepancy control at every8/12/24 scale")
    parser.add_argument("--full-trf", action="store_true", help="Disable discrepancy stopping for explicit mature control")
    parser.add_argument("--weakreg", type=float, default=0.)
    parser.add_argument("--worker-obs", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--worker-output", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--core-mask", type=int, help="Optional Windows affinity mask; affects this worker only")
    args = parser.parse_args()
    if args.discrepancy_all_scales and (args.full_trf or args.weakreg):
        parser.error("--discrepancy-all-scales is an independent unregularized control")
    solver_config=dict(discrepancy_stop=bool((args.discrepancy or args.discrepancy_all_scales) and not args.full_trf),
                       discrepancy_each_scale=bool(args.discrepancy_all_scales))
    from research.research_gn_nn_preservation_r4.src.anchors import generate_anchor, prepare_reference, dump
    if args.worker_obs is not None:
        # Fresh subprocesses recheck canonical completion after queue submission.
        # JSON is committed last when a local result is imported remotely.
        output=args.worker_output
        if output.exists() and output.with_suffix(".json").exists():
            receipt=json.loads(output.with_suffix(".json").read_text(encoding="utf-8"))
            expected=solver_config["discrepancy_stop"]
            if receipt.get("weakreg",0.)!=args.weakreg or bool(receipt.get("config",{}).get("discrepancy_stop",False))!=expected or bool(receipt.get("config",{}).get("discrepancy_each_scale",False))!=solver_config["discrepancy_each_scale"]:
                raise ValueError("Existing canonical receipt has incompatible public solver configuration")
            print(json.dumps(dict(worker_action="REUSED_COMPLETE_PAIR",status="REUSED_COMPLETE_PAIR",new_physics_counts={},canonical_origin=receipt.get("execution_origin",{}))))
            return
        reservations=args.output/"LOCAL_RESERVED_IDS.json"
        # Actual worker --output defaults to STUDY/anchors, matching canonical
        # paths for both the live remote queue and local finite batch.
        if sys.platform=="win32" and reservations.exists() and not args.weakreg and not args.discrepancy and not args.discrepancy_all_scales and not args.full_trf:
            reserved=json.loads(reservations.read_text(encoding="utf-8"))
            if output.stem in reserved.get("obs_ids",[]):
                print(json.dumps(dict(worker_action="DEFERRED_TO_LOCAL",status="DEFERRED_TO_LOCAL",new_physics_counts={})))
                return
        if args.core_mask and sys.platform == "win32":
            import ctypes
            ctypes.windll.kernel32.SetProcessAffinityMask(ctypes.windll.kernel32.GetCurrentProcess(), args.core_mask)
        receipt = generate_anchor(args.data/"operators.npz", args.worker_obs, args.worker_output, args.reference, args.weakreg, solver_config)
        print(json.dumps(dict(worker_action="COMPUTED_NEW", **{key: receipt[key] for key in ("status", "solve_wall_seconds", "feature_wall_seconds", "final_stage", "final_relative_residual", "counts")})))
        return
    if args.workers < 1 or args.workers > (4 if sys.platform == "win32" else 2):
        parser.error("At most 4 Windows or 2 local CPU workers")
    if args.id_file:
        args.obs_id = json.loads(args.id_file.read_text())
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
    suffix = "_discrepancy_allscales" if args.discrepancy_all_scales else ("_weakreg1e-5" if args.weakreg else ("_discrepancy" if args.discrepancy else ("_fulltrf" if args.full_trf else "")))
    tasks, skipped = [], []
    for record in records:
        output = args.output/(record["obs_id"]+suffix+".npz")
        if output.exists() and output.with_suffix(".json").exists():
            skipped.append(record["obs_id"])
            continue
        command = [sys.executable, str(Path(__file__).resolve()), "--data", str(args.data.resolve()),
                   "--reference", str(args.reference.resolve()), "--worker-obs", str((args.data/record["public_path"]).resolve()),
                   "--worker-output", str(output.resolve()), "--weakreg", str(args.weakreg)]
        if args.full_trf: command.append("--full-trf")
        if args.discrepancy: command.append("--discrepancy")
        if args.discrepancy_all_scales: command.append("--discrepancy-all-scales")
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
    costs, stops, canonical_costs, actions = Counter(), Counter(), Counter(), Counter()
    solve_seconds = feature_seconds = 0.
    for result in completions:
        receipt_path = Path(result["output_path"]).with_suffix(".json")
        action=result.get("worker_action","UNKNOWN_LEGACY");actions[action]+=1
        if action=="DEFERRED_TO_LOCAL":
            stops["DEFERRED_TO_LOCAL_NO_PHYSICS"]+=1
            continue
        if result["exit_code"] == 0 and receipt_path.exists():
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            canonical_costs.update(receipt["counts"])
            if action!="REUSED_COMPLETE_PAIR":costs.update(receipt["counts"])
            stops[receipt["stop_reason"] or receipt["status"]] += 1
            if action!="REUSED_COMPLETE_PAIR":
                solve_seconds += receipt["solve_wall_seconds"]
                feature_seconds += receipt["feature_wall_seconds"]
        else:
            stops["FAILED_COSTS_UNKNOWN"] += 1
    batch = dict(status="COMPLETED" if all(r["exit_code"] == 0 for r in completions) else "COMPLETED_WITH_FAILURES",
                 workers=args.workers, selected=len(records), launched=len(tasks), skipped_existing=skipped,
                 wall_seconds=time.perf_counter()-tick, shared_reference_setup=shared_setup,
                 aggregate_model_counts=dict(costs), canonical_case_counts=dict(canonical_costs), worker_actions=dict(actions), stop_histogram=dict(stops),
                 sum_solver_wall_seconds=solve_seconds, sum_feature_wall_seconds=feature_seconds,
                 jobs=completions, online_truth_access="NONE", weakreg=args.weakreg, **solver_config)
    batch_id = "_".join(args.split or ["selected"])+suffix+"_"+str(time.time_ns())
    dump(args.output/("BATCH_"+batch_id+".json"), batch)
    print(json.dumps({key: batch[key] for key in ("status", "selected", "launched", "wall_seconds")}))


if __name__ == "__main__":
    main()
