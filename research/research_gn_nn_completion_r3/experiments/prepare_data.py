"""Local full-VIE grouped data generation, never an inference entry point."""
import argparse
import json
from pathlib import Path
import sys

STUDY = Path(__file__).resolve().parents[1]
REPOSITORY = STUDY.parents[1]
sys.path.insert(0, str(REPOSITORY))
from research.research_gn_nn_completion_r3.src.data import prepare_dataset, expand_training_dataset, summarize_anchor_receipts, dump_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=STUDY/"data")
    parser.add_argument("--operators", type=Path, default=STUDY.parent/"research_z_alg_r4_diversity/public_data/operators.npz")
    parser.add_argument("--seed", type=int, default=20261107)
    parser.add_argument("--train", type=int, default=256)
    parser.add_argument("--val", type=int, default=64)
    parser.add_argument("--test", type=int, default=64)
    parser.add_argument("--dev", type=int, default=8)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--expand-train", type=int, choices=[1024, 2048],
                      help="Append training parents to an existing dataset; preserve all old files and entries")
    mode.add_argument("--summarize-anchors", type=Path,
                      help="Read-only receipt summary for an existing anchor directory; no generation")
    parser.add_argument("--summary-output", type=Path, help="Optional JSON destination for the read-only summary")
    args = parser.parse_args()
    if args.summarize_anchors is not None:
        receipt = summarize_anchor_receipts(args.output, args.summarize_anchors)
        if args.summary_output:
            dump_json(args.summary_output, receipt)
        print(json.dumps(receipt))
        return
    if args.expand_train is not None:
        receipt = expand_training_dataset(args.output, args.expand_train)
    else:
        receipt = prepare_dataset(args.output, args.operators, args.seed, args.train, args.val, args.test, args.dev)
    print(json.dumps({key: value for key, value in receipt.items() if key != "forward_receipts"}))


if __name__ == "__main__":
    main()
