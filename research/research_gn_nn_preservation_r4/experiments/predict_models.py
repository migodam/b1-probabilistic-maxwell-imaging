"""Public-only R4 bulk inference from explicit anchor lists; never reads labels."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import time

STUDY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STUDY/"src"))
from networks import OperatorCache, choose_device
from predict import predict_anchors


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--checkpoints", nargs="+", type=Path, required=True)
    p.add_argument("--anchor-list", type=Path, required=True, help="JSON list of public anchor paths")
    p.add_argument("--operators", type=Path, default=STUDY/"data"/"operators.npz")
    p.add_argument("--output", type=Path, default=STUDY/"predictions")
    p.add_argument("--device", default="mps")
    p.add_argument("--recurrent", action="store_true", help="Use the recurrent checkpoint family; keeps the same public-only anchor interface")
    args = p.parse_args()
    paths = [Path(v) for v in json.loads(args.anchor_list.read_text())]
    if not paths or not all(v.exists() and v.suffix == ".npz" for v in paths):
        raise ValueError("Explicit nonempty existing public anchor list required")
    device = choose_device(args.device)
    t = time.perf_counter()
    operators = OperatorCache(args.operators).to(device)
    setup_seconds = time.perf_counter()-t
    summaries=[]
    if args.recurrent:
        from predict_recurrent import predict_recurrent_anchors
    for checkpoint in args.checkpoints:
        out = args.output/checkpoint.parent.name
        receipts = (predict_recurrent_anchors(checkpoint, paths, args.operators, out, args.device, operators=operators)
                    if args.recurrent else predict_anchors(checkpoint, paths, args.operators, out, args.device, operators))
        summaries.append({"model": checkpoint.parent.name, "output": str(out), "count":len(receipts), "truth_access":"NONE"})
    args.output.mkdir(parents=True,exist_ok=True)
    receipt=json.dumps({"global_operator_setup_seconds":setup_seconds,
        "operator_uploads":1,"recurrent_family":args.recurrent,"runs":summaries},indent=2)+"\n"
    (args.output/f"prediction_runs_{time.time_ns()}.json").write_text(receipt)
    (args.output/"prediction_runs.json").write_text(receipt)


if __name__ == "__main__":
    main()
