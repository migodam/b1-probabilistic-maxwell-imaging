"""Explicit Track C paid entry point. Import and --help perform no physics."""
from pathlib import Path
import argparse
import json
import os
import sys

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from p3r2.core import OUT, Book, Target, write
from p3r2.geometry import ARMS, FEATURES, load_policy, run_trajectory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", required=True)
    parser.add_argument("--parents", nargs="+", default=[f"DEV_{k:03d}" for k in range(4)])
    parser.add_argument("--arms", nargs="+", choices=ARMS, default=list(ARMS[:4]))
    parser.add_argument("--steps", type=int, default=8)
    parser.add_argument("--max-wall", type=float, default=600.)
    parser.add_argument("--weak-dimension", type=int, default=4)
    parser.add_argument("--corrector-iterations", type=int, default=35)
    parser.add_argument("--policy", type=Path)
    parser.add_argument("--model-policy", type=Path)
    parser.add_argument("--seed", type=int, default=2026102203)
    parser.add_argument("--anchor-directory", type=Path, default=OUT/"A/anchors")
    parser.add_argument("--split", choices=("train", "validation", "exposed_eval", "fresh_eval"), default="train")
    args = parser.parse_args()
    if not 1 <= args.steps <= 32 or not 0 < args.max_wall <= 1200:
        parser.error("steps must be 1..32 and controlled wall must be <=1200s")
    if args.split == "train" and any(p not in [f"DEV_{k:03d}" for k in range(16)] for p in args.parents):
        parser.error("train parents predeclared DEV_000..DEV_015")
    if args.split == "validation" and any(p not in [f"DEV_{k:03d}" for k in range(16, 24)] for p in args.parents):
        parser.error("validation parents predeclared DEV_016..DEV_023")
    if args.split == "exposed_eval" and any(p not in [f"DEV_{k:03d}" for k in range(24, 32)] for p in args.parents):
        parser.error("exposed eval parents predeclared DEV_024..DEV_031")
    if args.split == "fresh_eval" and any(p.startswith(("DEV_", "LOCKED_")) for p in args.parents):
        parser.error("historically exposed DEV/LOCKED parents cannot be relabeled fresh")
    if args.split != "train" and "dynamic_bandit" in args.arms and args.policy is None:
        parser.error("evaluation bandit requires a frozen train-only policy")
    if args.split != "train" and "dynamic_mbrl" in args.arms and args.model_policy is None:
        parser.error("evaluation model-based arm requires a separate frozen train-only transition model")
    if len(set(args.parents)) != len(args.parents) or len(set(args.arms)) != len(args.arms):
        parser.error("duplicate parents/arms prohibited")
    policy = load_policy(args.policy) if args.policy is not None else None
    model_policy = load_policy(args.model_policy) if args.model_policy is not None else None
    if policy is not None and policy["kind"] != "ridge_contextual_bandit":
        parser.error("--policy expects the bandit export")
    if model_policy is not None and model_policy["kind"] != "linear_model_based_policy":
        parser.error("--model-policy expects the model-based export")
    # Validate online inputs and paid caches before opening a paid job.
    for parent in args.parents:
        if not (args.anchor_directory/f"{parent}.npz").is_file():
            raise FileNotFoundError(f"missing paid anchor for {parent}")
        if not (OUT/f"inputs/online/{parent}.npz").is_file():
            raise FileNotFoundError(f"missing online observation for {parent}")
    book = Book("C", args.job, max_wall=args.max_wall)
    status, error = "COMPLETE", None
    results = []
    try:
        protocol = dict(job=args.job, stage="C", split=args.split, parents=args.parents,
            arms=args.arms, steps=args.steps, radii=[.02, .04], feature_names=FEATURES,
            weak_dimension=args.weak_dimension, Gaussian_properties="native fixed u",
            task_constraints="actual per-step RMS radius and signed material progress >= half radius",
            direction_metric="unit material RMS; Jdata^T Jdata + original Gaussian I",
            predictor_corrector=True, corrector="strong complement only; no return to origin",
            full_states_max=2*args.steps*len(args.arms)*len(args.parents),
            fresh_endpoint_RHS_max=80*args.steps*len(args.arms)*len(args.parents),
            anchor_J_RHS_max=32*len(args.parents), max_wall=args.max_wall,
            policy_path=str(args.policy) if args.policy else None,
            model_policy_path=str(args.model_policy) if args.model_policy else None,
            GT_used_online=False, posterior_samples=False, interpolated_path_verified=False,
            historical_exposure=args.split != "fresh_eval")
        write(book.jobdir/"PROTOCOL.json", protocol)
        for name in ("tools/run_p3r2_geometry.py", "tools/train_p3r2_policy.py", "src/p3r2/geometry.py", "src/p3r2/core.py"):
            destination = book.jobdir/"source"/name
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes((ROOT/name).read_bytes())
        target = Target(book, keep_fields=True)
        for parent_index, parent in enumerate(args.parents):
            anchor = target.load_cache(args.anchor_directory/f"{parent}.npz",
                owner="P3R2 C reused y-derived paid P1-compatible anchor; historical exposure recorded")
            for arm in args.arms:
                result = run_trajectory(target, parent, anchor, arm=arm,
                    directory=book.jobdir/"trajectories"/parent/arm, steps=args.steps,
                    weak_dimension=args.weak_dimension, policy=model_policy if arm == "dynamic_mbrl" else policy,
                    seed=args.seed+parent_index, max_corrector_iterations=args.corrector_iterations)
                results.append(result)
                write(book.jobdir/"RESULTS.json", dict(protocol=protocol, trajectories=results,
                    status="IN_PROGRESS", GT_used_online=False, candidate_probabilities=None))
                print(json.dumps(dict(parent=parent, arm=arm, status=result["status"],
                    completed_steps=result["completed_steps"], cost=result["cost"])), flush=True)
        write(book.jobdir/"RESULTS.json", dict(protocol=protocol, trajectories=results,
            status="COMPLETE", GT_used_online=False, candidate_probabilities=None))
    except BaseException as exc:
        status, error = "ERROR", repr(exc)
        write(book.jobdir/"ERROR.json", dict(status=status, error=error, partial_outputs_preserved=True))
        raise
    finally:
        book.close(status, error)


if __name__ == "__main__":
    main()
