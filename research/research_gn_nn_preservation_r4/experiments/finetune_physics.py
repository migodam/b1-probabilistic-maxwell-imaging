"""Optional matched R4 fine-tuning; explicit parent VAL/time authorization required."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

STUDY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STUDY/"src"))
from train import finetune_model


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--checkpoint", type=Path, required=True, help="One R4 best checkpoint shared by both arms")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--data-root", type=Path, default=STUDY/"data")
    p.add_argument("--manifest", type=Path, default=STUDY/"data/manifests/evaluation.json")
    p.add_argument("--operators", type=Path, default=STUDY/"data/operators.npz")
    p.add_argument("--anchors", type=Path, default=STUDY/"anchors")
    p.add_argument("--contract", type=Path, default=STUDY/"audit/OPTIONAL_FT_CONTRACT.json")
    p.add_argument("--device", default="mps")
    p.add_argument("--authorized", action="store_true", help="Parent has passed the frozen corrected-VAL and remaining-time gate")
    p.add_argument("--authorization", type=Path, help="Optional parent gate receipt recorded with both arms")
    args = p.parse_args()
    if not args.authorized:
        p.error("Actual optional fine-tuning awaits the parent's corrected-VAL/time gate authorization")
    authorization = json.loads(args.authorization.read_text()) if args.authorization else None
    result = finetune_model(args.checkpoint, args.data_root, args.manifest, args.operators, args.output,
                            args.anchors, args.contract, args.device, authorized=True,
                            authorization=authorization)
    print(json.dumps({"status": result["status"], "output": str(args.output),
                      "modes": [r["mode"] for r in result["runs"]],
                      "pair_wall_seconds": result["pair_wall_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
