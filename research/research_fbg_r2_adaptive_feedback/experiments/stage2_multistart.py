"""Plan finite K=1/4/8 multistart jobs, then select by public physics score."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from src.multistart import make_starts, plan, select_public, START_NAMES
import numpy as np


def write_new(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as f:
        json.dump(data, f, indent=2, allow_nan=False)
        f.write("\n")


def create(a) -> None:
    cases = json.loads(a.cases.read_text())["cases"]
    chosen = set(a.case_ids)
    cases = [c for c in cases if c["case_id"] in chosen]
    if {c["case_id"] for c in cases} != chosen:
        raise ValueError("Unknown Stage 2 case ID")
    a.output.mkdir(parents=True, exist_ok=True)
    starts = make_starts(a.seed)
    path = a.output / "starts.npz"
    if path.exists():
        raise FileExistsError(path)
    np.savez_compressed(path, starts=starts, names=np.asarray(START_NAMES))
    manifest = plan(cases, a.seed, a.wall, a.lu, a.rhs)
    manifest["starts_path"] = "starts.npz"  # relative to plan.json for transport
    write_new(a.output / "plan.json", manifest)
    print(json.dumps({"planned_tasks": len(manifest["tasks"]), "output": str(a.output)}))


def select(a) -> None:
    manifest = json.loads(a.plan.read_text())
    grouped = {}
    for t in manifest["tasks"]:
        key = (t["case_id"], t["K"], t["regime"])
        grouped.setdefault(key, []).append(t)
    rows = []
    for (case_id, k, regime), tasks in grouped.items():
        receipts = []
        for t in tasks:
            path = a.receipts / (t["task_id"] + ".json")
            if path.exists():
                r = json.loads(path.read_text())
                r["task_id"] = t["task_id"]
                receipts.append(r)
        selected = select_public(receipts)
        rows.append({"case_id": case_id, "K": k, "regime": regime,
                     "planned_starts": k, "received": len(receipts),
                     "missing": k-len(receipts), "selected": selected,
                     "selection_uses_truth": False,
                     "cost_accounting": "Sum measured independent start receipts; never infer observed cost from allowances."})
    write_new(a.output, {"schema": "fbg-r2-multistart-selection-v1",
                         "status": "DEVELOPMENT_PUBLIC_SCORE_SELECTION",
                         "rows": rows})
    print(json.dumps({"groups": len(rows), "output": str(a.output)}))


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="command", required=True)
    p = sub.add_parser("plan")
    p.add_argument("--cases", type=Path, default=HERE / "data_contract/stage1_cases.json")
    p.add_argument("--case-ids", nargs="+", default=["dev_moderate_ring", "dev_high_digit8", "dev_high_ring"])
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--seed", type=int, default=20261030)
    p.add_argument("--wall", type=float, default=45.)
    p.add_argument("--lu", type=int, default=2000)
    p.add_argument("--rhs", type=int, default=60000)
    s = sub.add_parser("select")
    s.add_argument("--plan", type=Path, required=True)
    s.add_argument("--receipts", type=Path, required=True)
    s.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    if a.command == "plan":
        create(a)
    else:
        select(a)


if __name__ == "__main__":
    main()
