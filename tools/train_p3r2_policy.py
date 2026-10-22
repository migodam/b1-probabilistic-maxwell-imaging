"""OFFLINE train-only contextual bandit reward fitting; zero Maxwell calls.

Recorded choices are supervised with material-MSE improvements per measured
step and a frozen physical cost penalty. This is offline reward regression;
it is not a posterior density, an off-policy value certificate, or full RL.
"""
from pathlib import Path
import argparse
import json
import os
import sys
import time

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from p3r2.core import OUT, Book, write
from p3r2.geometry import FEATURES, load_policy, policy_score

TRAIN = {f"DEV_{k:03d}" for k in range(16)}


def fit_model(features, rewards, ridge):
    X, y = np.asarray(features, float), np.asarray(rewards, float)
    mean, scale = X.mean(axis=0), X.std(axis=0)
    mean[0], scale[0] = 0., 1.
    scale[scale < 1e-8] = 1.
    Z = (X-mean)/scale
    penalty = ridge*np.eye(Z.shape[1]); penalty[0, 0] = 0.
    coefficients = np.linalg.solve(Z.T@Z+penalty, Z.T@y)
    return mean, scale, coefficients


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", required=True)
    parser.add_argument("--trajectory-jobs", nargs="+", required=True, type=Path)
    parser.add_argument("--kind", choices=("bandit", "model_based"), default="bandit")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--horizon", type=int, choices=(2, 3), default=3)
    parser.add_argument("--discount", type=float, default=.9)
    parser.add_argument("--truth-directory", type=Path, default=OUT/"inputs/offline_only")
    parser.add_argument("--ridge", type=float, default=10.)
    parser.add_argument("--cost-penalty", type=float, default=.001,
        help="reward penalty per 80 physical RHS, frozen before fitting")
    parser.add_argument("--min-rows", type=int, default=128)
    args = parser.parse_args()
    if args.output is None:
        args.output = OUT/"C/policy"/("FROZEN_POLICY.json" if args.kind == "bandit" else "FROZEN_MODEL_POLICY.json")
    if args.ridge <= 0 or args.cost_penalty < 0 or not 0 < args.discount <= 1:
        parser.error("positive ridge and nonnegative cost penalty required")
    if args.output.exists():
        raise FileExistsError("never overwrite a frozen policy")
    # Load online-generated records before opening any truth file. Check whole
    # source jobs: silently filtering eval rows would conceal a leaked source.
    records = []
    for job in args.trajectory_jobs:
        protocol = json.loads((job/"PROTOCOL.json").read_text())
        if protocol.get("split") != "train" or not set(protocol["parents"]) <= TRAIN:
            raise ValueError("training accepts only predeclared train-only jobs")
        closed = json.loads((job/"CLOSED.json").read_text())
        if closed.get("status") != "COMPLETE":
            raise ValueError("training source paid job must close COMPLETE")
        for trace in sorted((job/"trajectories").glob("*/*/TRACE.jsonl")):
            for line in trace.read_text().splitlines():
                row = json.loads(line)
                if row["parent"] not in TRAIN or row.get("GT_used_online") is not False:
                    raise ValueError("training trace split/GT provenance")
                if row["feature_names"] != list(FEATURES):
                    raise ValueError("training feature schema")
                records.append(row)
    unique = {(r["parent"], r["arm"], r["step"], r["endpoint"]) for r in records}
    if len(unique) != len(records):
        raise ValueError("duplicate source transition")
    parents = sorted({r["parent"] for r in records})
    if len(records) < args.min_rows or len(parents) < 4:
        raise ValueError("need >=128 train transitions and >=4 distinct development parents by default")
    book = Book("C", args.job, max_wall=600, analysis_only=True)
    status, error = "COMPLETE", None
    try:
        truth = {}
        for parent in parents:
            with np.load(args.truth_directory/f"{parent}.npz", allow_pickle=False) as archive:
                truth[parent] = archive["chi"].copy()
            if truth[parent].shape != (576,):
                raise ValueError("offline truth material layout")
        labels, X = [], []
        # All endpoint reads are paid caches; no Target construction exists.
        endpoint_paths = {}
        for job in args.trajectory_jobs:
            endpoint_paths.update({p.stem:p for p in (job/"endpoints").glob("*.npz")})
            endpoint_paths.update({p.stem:p for p in (job/"trajectories").glob("*/*/cache/*.npz")})
        for row in records:
            parent = row["parent"]
            endpoint = Path(row["cache_path"])
            origin = endpoint_paths.get(row["origin_endpoint"])
            if origin is None:
                origin = OUT/"A/anchors"/f"{parent}.npz"
            with np.load(endpoint, allow_pickle=False) as a, np.load(origin, allow_pickle=False) as b:
                endpoint_risk = float(np.mean(abs(a["chi"]-truth[parent])**2))
                origin_risk = float(np.mean(abs(b["chi"]-truth[parent])**2))
            # Relative local risk improvement has a stable material scale; no
            # likelihood softmax or candidate probabilities are manufactured.
            penalty = args.cost_penalty*row["cost"].get("linear_solve_rhs", 0)/80
            reward = (origin_risk-endpoint_risk)/max(origin_risk, 1e-8)-penalty
            X.append(row["features"])
            labels.append(dict(parent=parent, arm=row["arm"], step=row["step"],
                endpoint=row["endpoint"], reward=reward, origin_MSE=origin_risk, endpoint_MSE=endpoint_risk,
                physical_RHS=row["cost"].get("linear_solve_rhs", 0), GT_used="OFFLINE_TRAIN_REWARD_ONLY"))
        X = np.asarray(X, float); y = np.array([v["reward"] for v in labels])
        mean, scale, coefficients = fit_model(X, y, args.ridge)
        model = dict(kind="ridge_contextual_bandit" if args.kind == "bandit" else "linear_model_based_policy", schema_version=1, feature_names=list(FEATURES),
            mean=mean, scale=scale, coefficients=coefficients, ridge=args.ridge,
            training_split="train", frozen_before_evaluation=True, fit_rows=len(X), fit_parents=parents,
            reward_definition=f"relative per-step complex material MSE reduction minus {args.cost_penalty} * physical RHS/80",
            no_posterior_claim=True)
        transition_report = None
        if args.kind == "model_based":
            # Only real consecutive train transitions with the same measured
            # endpoint/origin join are used. No fabricated transition targets.
            by_origin = {row["origin_endpoint"]:row for row in records}
            transitions = [(row, by_origin[row["endpoint"]]) for row in records
                if row["endpoint"] in by_origin and by_origin[row["endpoint"]]["parent"] == row["parent"]
                and by_origin[row["endpoint"]]["arm"] == row["arm"]
                and by_origin[row["endpoint"]]["step"] == row["step"]+1]
            if len(transitions) < 32:
                raise ValueError("model-based fitting needs >=32 measured consecutive train transition pairs")
            context_indices = [9, 10, 11, 12]
            TX = np.array([pair[0]["features"] for pair in transitions], float)
            TY = np.array([pair[1]["features"] for pair in transitions], float)[:, context_indices]
            Z = (TX-mean)/scale
            penalty = args.ridge*np.eye(Z.shape[1]); penalty[0, 0] = 0.
            transition_coefficients = np.linalg.solve(Z.T@Z+penalty, Z.T@TY)
            action_counts = {str(action):sum(int(row["action"]) == action for row in records)
                for action in sorted({int(row["action"]) for row in records})}
            limits = ["offline behavior coverage is sparse and nonrandom",
                "future geometry/action bank is locally frozen; only online context evolves in rollout",
                "future radius response uses a local quadratic scaling approximation",
                "future decoder feasibility is checked only at the next real step",
                "bootstrapped planning inherits learned transition and reward error; no off-policy value guarantee"]
            model.update(transition_coefficients=transition_coefficients, context_indices=context_indices,
                context_min=TY.min(axis=0), context_max=TY.max(axis=0), horizon=args.horizon,
                discount=args.discount, step_increment=float(np.median([
                    pair[1]["features"][2]-pair[0]["features"][2] for pair in transitions])),
                transition_pairs=len(transitions), behavior_action_counts=action_counts,
                planning_limits=limits)
            transition_report = dict(pairs=len(transitions), context_feature_names=[FEATURES[i] for i in context_indices],
                train_context_MSE=np.mean((Z@transition_coefficients-TY)**2, axis=0),
                behavior_action_counts=action_counts, horizon=args.horizon, discount=args.discount,
                limits=limits)
        # Honest leave-parent-out prediction errors, with no hyperparameter
        # selection on exposed eval/validation parents and no OPE claim.
        loo = []
        for parent in parents:
            test = np.array([row["parent"] == parent for row in labels])
            m, s, c = fit_model(X[~test], y[~test], args.ridge)
            pred = ((X[test]-m)/s)@c
            loo.append(dict(parent=parent, rows=int(test.sum()), reward_MSE=float(np.mean((pred-y[test])**2))))
        write(args.output, model)
        policy = load_policy(args.output)
        pred = policy_score(policy, X)
        write(book.jobdir/"OFFLINE_LABELS.json", dict(labels=labels, GT_online=False))
        write(book.jobdir/"RESULTS.json", dict(policy_path=args.output, rows=len(X), parents=parents,
            train_reward_MSE=float(np.mean((pred-y)**2)), leave_parent_out=loo,
            GT_used="OFFLINE_TRAIN_ONLY", online_policy_has_no_GT=True,
            off_policy_value_claim=False, posterior_claim=False,
            finite_horizon_model_based_policy=args.kind == "model_based", transition_model=transition_report,
            PPO_executed=False,
            limitation="limited behavior coverage; online learned-arm evaluation still required"))
        print(json.dumps(dict(status="COMPLETE", policy=str(args.output), rows=len(X), parents=len(parents))), flush=True)
    except BaseException as exc:
        status, error = "ERROR", repr(exc)
        raise
    finally:
        book.close(status, error)


if __name__ == "__main__":
    main()
