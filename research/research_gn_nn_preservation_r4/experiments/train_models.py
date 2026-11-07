"""Sequential R4 trainer; fixed operators/stats reused across from-scratch arms."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import time

STUDY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STUDY/"src"))
from networks import ARMS, OperatorCache, choose_device, fit_feature_stats
from train import train_model, load_anchor, released_training_cohort
from predict import predict_anchors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-ready", action="store_true")
    parser.add_argument("--data-root", type=Path, default=STUDY/"data")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--operators", type=Path)
    parser.add_argument("--anchors", type=Path, default=STUDY/"anchors")
    parser.add_argument("--output", type=Path, default=STUDY/"models")
    parser.add_argument("--arms", nargs="+", choices=ARMS, default=list(ARMS[:-1]))
    parser.add_argument("--seeds", type=int, nargs="+", default=[7])
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--patience", type=int, default=15)
    parser.add_argument("--device", default="mps")
    parser.add_argument("--cascade-passes", type=int, choices=(2,), default=2)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--alpha", type=float, default=.3)
    parser.add_argument("--loss", choices=("material", "const", "internal"), default="material")
    parser.add_argument("--cache-ram-mb", type=int, default=1536)
    parser.add_argument("--predict-val", action="store_true", help="Immediately write public VAL predictions per trained arm")
    args = parser.parse_args()
    if not args.data_ready:
        parser.error("--data-ready is required after complete TRAIN/VAL anchor release")
    manifest = args.manifest or args.data_root/"manifests"/"evaluation.json"
    operators_path = args.operators or args.data_root/"operators.npz"
    cohort = released_training_cohort(args.data_root, manifest, args.anchors)
    train, val = cohort["train"], cohort["val"]
    args.output.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    device = choose_device(args.device)
    operators = OperatorCache(operators_path).to(device)
    stats = fit_feature_stats(load_anchor(args.anchors/f"{r['obs_id']}.npz") for r in train)
    setup = {"seconds": time.perf_counter()-start, "operator_uploads": 1,
             "stats_fitted_on": "TRAIN_public_anchors_only", "train_count": len(train),
             "projectors": operators.projector_metadata, "device": str(device)}
    setup["invocation_unix_ns"] = time.time_ns()
    text = json.dumps(setup, indent=2)+"\n"
    (args.output/f"shared_setup_{setup['invocation_unix_ns']}.json").write_text(text)
    (args.output/"shared_setup.json").write_text(text)
    results = []
    for seed in args.seeds:
        if seed not in (7,17,27):
            raise ValueError("Frozen R4 seeds are 7/17/27")
        for arm in args.arms:
            suffix = f"_a{args.alpha:g}" if arm == "soft_dark" else ""
            suffix += f"_{args.loss}" if args.loss != "material" else ""
            output = args.output/f"{arm}_seed{seed}_cascade2{suffix}"
            if (output/"best.pt").exists():
                raise FileExistsError(f"Refusing to overwrite prior trained model: {output}")
            result = train_model(args.data_root, manifest, operators_path, output, arm,
                                 seed, args.epochs, args.batch_size, args.patience,
                                 args.device, args.cascade_passes, args.anchors, args.num_workers,
                                 cache_ram_mb=args.cache_ram_mb, alpha=args.alpha, loss_name=args.loss,
                                 operators=operators, feature_stats=stats, cohort=cohort)
            if args.predict_val:
                pred_dir = STUDY/"predictions"/output.name
                receipts = predict_anchors(output/"best.pt", [args.anchors/f"{r['obs_id']}.npz" for r in val],
                                           operators_path, pred_dir, args.device, operators)
                result["public_val_predictions"] = {"output": str(pred_dir), "count": len(receipts), "truth_access": "NONE"}
                print(json.dumps({"status": "ARM_TRAINED_VAL_PREDICTED", "arm": arm, "seed": seed,
                                  "checkpoint": str(output/"best.pt"), "predictions": str(pred_dir)}), flush=True)
            results.append(result)
            (args.output/f"training_runs_{args.loss}_seed{'_'.join(map(str,args.seeds))}.json").write_text(json.dumps(results, indent=2)+"\n")


if __name__ == "__main__":
    main()
