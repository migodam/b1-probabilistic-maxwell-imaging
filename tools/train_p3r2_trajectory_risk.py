"""D3 train-only trajectory-risk action and separate exposed-only offline audit.

No physics execution. Training reads only DEV_000..007 truth after COMPLETE
C_train_r0 closure. Optional evaluation reads only DEV_024..031 truth after
its own COMPLETE closure and after freezing the exported estimator.
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
from p3r2.core import OUT, Book, save, write
from p3r2.trajectory_risk import METHODS, REGULARIZATION, apply_frozen, fit_convex_weights, load_frozen

OLD = ROOT.parent/"b1-p3l-local-posterior/research/p3l_local_posterior_imaging"
TRAIN = [f"DEV_{k:03d}" for k in range(8)]
EXPOSED = [f"DEV_{k:03d}" for k in range(24, 32)]
ARMS = ["frozen_fixed", "dynamic_fixed", "dynamic_curvature", "dynamic_bandit", "ordinary_random"]


def load_online(job, split, parents):
    job = Path(job).resolve()
    protocol = json.loads((job/"PROTOCOL.json").read_text())
    closure = json.loads((job/"CLOSED.json").read_text())
    if closure.get("status") != "COMPLETE":
        raise ValueError("COMPLETE paid online closure required before reading any truth")
    if protocol.get("split") != split or protocol.get("parents") != parents or protocol.get("arms") != ARMS or protocol.get("steps") != 8:
        raise ValueError("predeclared eight-parent/five-arm/eight-step source protocol required")
    stacks, metadata = {}, []
    for parent in parents:
        # The original paid anchor, rather than an endpoint chosen using GT.
        anchor_path = OLD/"A/anchors"/f"{parent}.npz"
        with np.load(anchor_path, allow_pickle=False) as archive:
            anchor = archive["chi"].copy()
        for arm in ARMS:
            trace = job/"trajectories"/parent/arm/"TRACE.jsonl"
            rows = [json.loads(line) for line in trace.read_text().splitlines()] if trace.is_file() else []
            if any(row.get("GT_used_online") is not False or row["parent"] != parent or row["arm"] != arm for row in rows):
                raise ValueError("online GT provenance failure")
            materials = [anchor]
            for index, row in enumerate(rows):
                if row["step"] != index:
                    raise ValueError("ordered consecutive path required")
                if index == 0 and row["initial_anchor_endpoint"] != parent:
                    raise ValueError("original paid anchor endpoint provenance")
                path = Path(row["cache_path"]).resolve()
                if not path.is_relative_to(job):
                    raise ValueError("trajectory endpoint must be in its paid source job")
                with np.load(path, allow_pickle=False) as archive:
                    materials.append(archive["chi"].copy())
            stacks[(parent, arm)] = np.stack(materials)
            metadata.append(dict(parent=parent, arm=arm, nodes=len(materials), original_anchor=str(anchor_path),
                source_trace=str(trace), paid_RHS=sum(row["cost"]["linear_solve_rhs"] for row in rows)))
    return stacks, dict(job=str(job), protocol=protocol, paid_closure=closure, paths=metadata)


def load_truth(parents):
    if not set(parents) <= set(TRAIN+EXPOSED):
        raise ValueError("this tool cannot read fresh or locked GT")
    result = {}
    for parent in parents:
        with np.load(OLD/"inputs/offline_only"/f"{parent}.npz", allow_pickle=False) as archive:
            result[parent] = archive["chi"].copy()
        if result[parent].shape != (576,) or not np.isfinite(result[parent]).all():
            raise ValueError("finite full-grid truth required")
    return result


def evaluate(stacks, truth, parents, frozen, *, destination=None):
    rows = []
    for parent in parents:
        for arm in ARMS:
            for method in METHODS:
                estimate, info = apply_frozen(stacks[(parent, arm)], frozen, arm, method=method)
                delta = estimate-truth[parent]
                row = dict(parent=parent, arm=arm,
                    full_grid_complex_MSE=float(np.mean(abs(delta)**2)),
                    full_grid_complex_MAE=float(np.mean(abs(delta))), **info)
                rows.append(row)
                if destination is not None:
                    save(Path(destination)/parent/arm/f"{method}.npz", chi=estimate)
    summary = []
    for arm in ARMS:
        anchor = [row for row in rows if row["arm"] == arm and row["method"] == "anchor"]
        base = float(np.mean([row["full_grid_complex_MSE"] for row in anchor]))
        for method in METHODS:
            selected = [row for row in rows if row["arm"] == arm and row["method"] == method]
            mse = float(np.mean([row["full_grid_complex_MSE"] for row in selected]))
            summary.append(dict(arm=arm, method=method, parents=len(selected), mean_MSE=mse,
                mean_MAE=float(np.mean([row["full_grid_complex_MAE"] for row in selected])),
                MSE_change_percent_vs_anchor=100*(mse/base-1),
                improved_parent_count=sum(row["full_grid_complex_MSE"] < ref["full_grid_complex_MSE"] for row, ref in zip(selected, anchor)),
                fallback_count=sum(row["fallback_to_anchor"] for row in selected)))
    return dict(rows=rows, summary=summary, parents=parents, posterior_claim=False,
        coefficients_are_probabilities=False, physics_RHS=0)


def write_table(path, results, title):
    lines = [f"# {title}", "", "Complex full-grid risk; lower is better. Convex coefficients are estimator weights, not posterior probabilities.",
        "", "| arm | estimator | MSE | MAE | MSE change vs anchor | improved parents | fallback |",
        "|---|---|---:|---:|---:|---:|---:|"]
    for row in results["summary"]:
        lines.append(f"| {row['arm']} | {row['method']} | {row['mean_MSE']:.9g} | {row['mean_MAE']:.9g} | {row['MSE_change_percent_vs_anchor']:+.3f}% | {row['improved_parent_count']}/8 | {row['fallback_count']} |")
    lines.extend(["", "Eight train and eight historically exposed evaluation parents are small panels. Exposed evaluation is not blind.",
        "No posterior density, posterior mean, mode coverage, or path compatibility claim follows from this estimator.",
        "A convex material average is a risk action and may not equal a single legal decoder output.",
        "Incomplete paths fall back to the original paid anchor and remain in the risk table."])
    Path(path).write_text("\n".join(lines)+"\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", default="D3_trajectory_integration")
    parser.add_argument("--train-job", type=Path, default=OUT/"C/jobs/C_train_r0")
    parser.add_argument("--exposed-job", type=Path, default=OUT/"C/jobs/C_exposed_eval_r0")
    parser.add_argument("--output-directory", type=Path, default=OUT/"D/trajectory_integration")
    parser.add_argument("--evaluate-exposed", action="store_true")
    parser.add_argument("--evaluate-only", action="store_true")
    parser.add_argument("--reuse-frozen", action="store_true",
        help="rebuild missing reports from an existing frozen export without refitting or overwriting it")
    args = parser.parse_args()
    destination = args.output_directory
    if args.evaluate_only and not args.evaluate_exposed:
        parser.error("evaluate-only requires --evaluate-exposed")
    if not args.evaluate_only and not args.reuse_frozen and (destination/"FROZEN.json").exists():
        raise FileExistsError("never overwrite frozen trajectory-risk weights")
    if args.evaluate_exposed and (destination/"EXPOSED_RESULTS.json").exists():
        raise FileExistsError("exposed evaluation already exists")
    book = Book("D", args.job, max_wall=600, analysis_only=True)
    status, error = "COMPLETE", None
    try:
        destination.mkdir(parents=True, exist_ok=True)
        for name in ("src/p3r2/trajectory_risk.py", "tools/train_p3r2_trajectory_risk.py"):
            path = book.jobdir/"source"/name; path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes((ROOT/name).read_bytes())
        if not args.evaluate_only:
            stacks, source = load_online(args.train_job, "train", TRAIN)
            # Truth phase starts only after the complete online inputs are loaded.
            truth = load_truth(TRAIN)
            weights, fits = {}, {}
            existing = load_frozen(destination/"FROZEN.json") if args.reuse_frozen else None
            if existing is not None and (existing["train_parents"] != TRAIN or list(existing["arms"]) != ARMS):
                raise ValueError("existing frozen train panel mismatch")
            for arm in ARMS:
                panel = np.stack([stacks[(parent, arm)] for parent in TRAIN])
                if panel.shape != (8, 9, 576):
                    raise ValueError("training requires all eight complete paths; no GT-dependent filtering")
                if existing is None:
                    w, report = fit_convex_weights(panel, np.stack([truth[parent] for parent in TRAIN]))
                else:
                    w = np.array(existing["arms"][arm], float)
                    estimates = np.einsum("n,pnc->pc", w, panel)
                    mse = float(np.mean(abs(estimates-np.stack([truth[parent] for parent in TRAIN]))**2))
                    e0 = np.zeros(9); e0[0] = 1.
                    report = dict(reused_frozen_without_refit=True, train_full_grid_complex_MSE=mse,
                        regularization=REGULARIZATION, regularization_penalty=float(REGULARIZATION*np.dot(w-e0, w-e0)))
                weights[arm], fits[arm] = w, report
            frozen = dict(schema_version=1, kind="constant_convex_trajectory_risk_action", training_split="train",
                train_parents=TRAIN, nodes=9, arms=weights, regularization=REGULARIZATION,
                objective="mean_parent mean_fullgrid_cell abs(sum_node(w*chi)-GT)^2 + 1e-4*||w-e_anchor||^2; simplex coefficients",
                frozen_before_fresh_evaluation=True, historical_train_exposure=True,
                coefficient_semantics="constant convex material-estimator coefficients, not probabilities",
                baselines={"anchor":"original paid anchor", "last_point":"eighth measured endpoint",
                    "equal_steps":"equal coefficients over anchor plus eight ordered endpoints (nine nodes)",
                    "arc_length_trapezoid":"material-RMS segment-length trapezoid coefficients, normalized by total length"},
                limitations=["8 training parents only", "exposed evaluation is historically exposed, not blind",
                    "risk action, not posterior integration", "material average may not be a single legal decoder state",
                    "missing arm or incomplete future path falls back to anchor"], source_online_job=str(args.train_job.resolve()))
            if existing is None:
                write(destination/"FROZEN.json", frozen)
            frozen = load_frozen(destination/"FROZEN.json")
            train_results = evaluate(stacks, truth, TRAIN, frozen)
            write(destination/"TRAIN_RESULTS.json", dict(**train_results, fit_reports=fits, source=source))
            write_table(destination/"TRAIN_TABLE.md", train_results, "D3 training risk (8 parents)")
        else:
            frozen = load_frozen(destination/"FROZEN.json")
        if args.evaluate_exposed:
            # Separate post-freeze closure/read phase; no fit or coefficient
            # update is permitted using these evaluation labels.
            stacks, source = load_online(args.exposed_job, "exposed_eval", EXPOSED)
            truth = load_truth(EXPOSED)
            exposed = evaluate(stacks, truth, EXPOSED, frozen, destination=destination/"exposed_estimates")
            write(destination/"EXPOSED_RESULTS.json", dict(**exposed, split="historically_exposed_evaluation_not_blind",
                frozen_path=str(destination/"FROZEN.json"), source=source))
            write_table(destination/"EXPOSED_TABLE.md", exposed, "D3 historically exposed evaluation (8 parents; not blind)")
        write(book.jobdir/"RESULTS.json", dict(status="COMPLETE", frozen_path=destination/"FROZEN.json",
            exposed_evaluated=args.evaluate_exposed, physics_RHS=0, GT_used="TRAIN_AND_OPTIONAL_EXPOSED_OFFLINE_ONLY",
            fresh_locked_GT_read=False, independent_of_D1_and_D2=True))
        print(json.dumps(dict(status="COMPLETE", output=str(destination), physics_RHS=0, exposed_evaluated=args.evaluate_exposed)), flush=True)
    except BaseException as exc:
        status, error = "ERROR", repr(exc)
        write(book.jobdir/"ERROR.json", dict(error=error, partial_outputs_preserved=True))
        raise
    finally:
        book.close(status, error)


if __name__ == "__main__":
    main()
