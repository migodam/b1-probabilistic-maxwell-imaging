"""Sequential model launcher; requires an explicit data-ready acknowledgement."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

STUDY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STUDY/"src"))
from networks import ARMS
from train import train_model


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-ready", action="store_true", help="Data/anchors are released for training")
    parser.add_argument("--data-root", type=Path, default=STUDY/"data")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--operators", type=Path)
    parser.add_argument("--anchors", type=Path, default=STUDY/"anchors")
    parser.add_argument("--output", type=Path, default=STUDY/"models")
    parser.add_argument("--arms", nargs="+", choices=ARMS, default=list(ARMS))
    parser.add_argument("--seeds", type=int, nargs="+", default=[7])
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--patience", type=int, default=15)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--cascade-passes", type=int, choices=(1, 2), default=1)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--dev-only", action="store_true", help="Use explicitly released debugtrain/debugval only")
    args = parser.parse_args()
    if not args.data_ready:
        parser.error("--data-ready is required; training must wait for released data/anchors")
    manifest = args.manifest or args.data_root/"manifests"/"evaluation.json"
    operators = args.operators or args.data_root/"operators.npz"
    results = []
    for seed in args.seeds:
        for arm in args.arms:
            if args.cascade_passes == 2 and arm not in ("v_physics", "j_null"):
                continue
            suffix = "_cascade2" if args.cascade_passes == 2 else ""
            output = args.output/f"{arm}_seed{seed}{suffix}"
            results.append(train_model(args.data_root, manifest, operators, output, arm,
                                        seed, args.epochs, args.batch_size, args.patience,
                                        args.device, args.cascade_passes, args.anchors, args.num_workers, args.dev_only))
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/"training_runs.json").write_text(json.dumps(results, indent=2)+"\n")


if __name__ == "__main__":
    main()
