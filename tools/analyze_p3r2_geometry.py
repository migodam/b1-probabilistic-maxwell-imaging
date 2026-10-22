"""Post-closure offline trajectory risk/cost summary with explicit GT routing.

No Target or solver is constructed. A caller must select one split and its
truth directory; all online jobs must already have COMPLETE paid closures.
This evaluator does not certify a posterior or unobserved connecting paths.
"""
from pathlib import Path
import argparse
import json
import os
import sys

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from p3r2.core import OUT, Book, write


def risk(chi, truth):
    delta = np.asarray(chi)-np.asarray(truth)
    return dict(complex_MSE=float(np.mean(abs(delta)**2)), complex_MAE=float(np.mean(abs(delta))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", required=True)
    parser.add_argument("--trajectory-jobs", nargs="+", type=Path, required=True)
    parser.add_argument("--split", choices=("train", "validation", "exposed_eval", "fresh_eval"), required=True)
    parser.add_argument("--truth-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output is None:
        args.output = OUT/"C/analysis"/args.job/"RESULTS.json"
    if args.output.exists():
        raise FileExistsError("offline analysis output already exists")
    # Entire online manifests and traces are loaded before truth is opened.
    trajectories, sources = [], []
    seen = set()
    for job in args.trajectory_jobs:
        protocol = json.loads((job/"PROTOCOL.json").read_text())
        closed = json.loads((job/"CLOSED.json").read_text())
        if protocol.get("split") != args.split:
            raise ValueError("split-specific analysis cannot mix source splits")
        if closed.get("status") != "COMPLETE":
            raise ValueError("all source online jobs must close COMPLETE before any GT read")
        if args.split == "train" and any(p not in {f"DEV_{k:03d}" for k in range(16)} for p in protocol["parents"]):
            raise ValueError("train split provenance")
        if args.split == "fresh_eval" and any(p.startswith(("DEV_", "LOCKED_")) for p in protocol["parents"]):
            raise ValueError("historical exposed parents cannot become fresh evaluation")
        sources.append(dict(job=str(job), protocol=protocol, paid_closure=closed))
        for parent in protocol["parents"]:
            for arm in protocol["arms"]:
                key = (parent, arm)
                if key in seen:
                    raise ValueError("duplicate parent/arm comparison; use separate replicate analysis")
                seen.add(key)
                directory = job/"trajectories"/parent/arm
                trace = directory/"TRACE.jsonl"
                rows = [json.loads(line) for line in trace.read_text().splitlines()] if trace.is_file() else []
                if any(row.get("GT_used_online") is not False or row["parent"] != parent or row["arm"] != arm for row in rows):
                    raise ValueError("online trace GT/split provenance")
                trajectories.append(dict(parent=parent, arm=arm, directory=directory, rows=rows,
                    requested_steps=protocol["steps"]))
    book = Book("C", args.job, max_wall=600, analysis_only=True)
    status, error = "COMPLETE", None
    try:
        truth = {}
        for parent in sorted({row["parent"] for row in trajectories}):
            with np.load(args.truth_directory/f"{parent}.npz", allow_pickle=False) as archive:
                truth[parent] = archive["chi"].copy()
            if truth[parent].shape != (576,) or not np.isfinite(truth[parent]).all():
                raise ValueError("truth material array schema")
        summaries = []
        for trajectory in trajectories:
            parent, arm, rows = trajectory["parent"], trajectory["arm"], trajectory["rows"]
            if not rows:
                summaries.append(dict(parent=parent, arm=arm, status="NO_MEASURED_STEPS", final_risk=None,
                    requested_steps=trajectory["requested_steps"], completed_steps=0))
                continue
            cache = trajectory["directory"]/"cache"
            anchor_path = cache/f"{rows[0]['initial_anchor_endpoint']}.npz"
            with np.load(anchor_path, allow_pickle=False) as archive:
                anchor_risk = risk(archive["chi"], truth[parent])
            offline_rows = []
            for row in rows:
                with np.load(row["cache_path"], allow_pickle=False) as archive:
                    score = risk(archive["chi"], truth[parent])
                offline_rows.append(dict(step=row["step"], endpoint=row["endpoint"], risk=score,
                    Phi=row["Phi"], delta_Phi_from_anchor=row["delta_Phi_from_anchor"],
                    actual_step_RMS=row["actual_step_RMS"], anchor_displacement_RMS=row["anchor_displacement_RMS"],
                    signed_progress=row["signed_progress"], target_RMS=row["target_RMS"],
                    saturation_fraction=row["saturation_fraction"],
                    saturation_strata_changed_from_anchor=row["saturation_strata_changed_from_anchor"],
                    projector_rotation_step=row["projector_rotation_step"],
                    projector_rotation_from_anchor=row["projector_rotation_from_anchor"],
                    Maxwell_secant_relative_error=row["Maxwell_secant_relative_error"], cost=row["cost"]))
            final_risk = offline_rows[-1]["risk"]
            summary = dict(parent=parent, arm=arm,
                status="COMPLETE" if len(rows) == trajectory["requested_steps"] else "INCOMPLETE_RETAINED",
                requested_steps=trajectory["requested_steps"], completed_steps=len(rows),
                anchor_risk=anchor_risk, final_risk=final_risk,
                final_minus_anchor_MSE=final_risk["complex_MSE"]-anchor_risk["complex_MSE"],
                cost={k:sum(r["cost"].get(k, 0) for r in rows) for k in ("full_forward_calls", "linear_solve_rhs")},
                radius_or_progress_violations=sum(r["actual_step_RMS"] < r["target_RMS"]-1e-10
                    or r["actual_step_RMS"] > 1.02*r["target_RMS"]
                    or r["signed_progress"] < .5*r["target_RMS"]-1e-9 for r in rows),
                steps=offline_rows)
            summaries.append(summary)
        by_arm = {}
        for arm in sorted({row["arm"] for row in summaries}):
            panel = [row for row in summaries if row["arm"] == arm]
            complete = all(row["status"] == "COMPLETE" for row in panel)
            by_arm[arm] = dict(parents=len(panel), complete_trajectories=sum(row["status"] == "COMPLETE" for row in panel),
                mean_final_MSE=float(np.mean([row["final_risk"]["complex_MSE"] for row in panel])) if complete else None,
                mean_final_MAE=float(np.mean([row["final_risk"]["complex_MAE"] for row in panel])) if complete else None,
                mean_final_minus_anchor_MSE=float(np.mean([row["final_minus_anchor_MSE"] for row in panel])) if complete else None,
                aggregate_status="COMPLETE_PANEL" if complete else "NULL_INCOMPLETE_PANEL",
                measured_RHS=sum(row.get("cost", {}).get("linear_solve_rhs", 0) for row in panel),
                measured_states=sum(row.get("cost", {}).get("full_forward_calls", 0) for row in panel))
        result = dict(status="COMPLETE", split=args.split, truth_directory=str(args.truth_directory),
            GT_used="OFFLINE_POST_ONLINE_CLOSURE_ONLY", trajectories=summaries, by_arm=by_arm,
            source_jobs=sources, historical_exposure=args.split != "fresh_eval",
            candidate_probabilities=None, posterior_claim=False, continuous_path_verified=False,
            interpretation="material-risk scores compare measured candidate trajectories; no mass, mixing, or mode-coverage claim",
            cost_note="step totals exclude shared anchor/cache setup; authoritative source paid closures included")
        write(args.output, result)
        write(book.jobdir/"RESULTS.json", dict(output=args.output, split=args.split,
            parents=len(truth), trajectories=len(summaries), by_arm=by_arm, physics_RHS=0))
        print(json.dumps(dict(status="COMPLETE", output=str(args.output), parents=len(truth), arms=list(by_arm))), flush=True)
    except BaseException as exc:
        status, error = "ERROR", repr(exc)
        raise
    finally:
        book.close(status, error)


if __name__ == "__main__":
    main()
