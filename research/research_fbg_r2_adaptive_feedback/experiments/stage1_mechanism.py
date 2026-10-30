"""Freeze a small disclosed-data mechanism set and paired run list."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", type=Path,
                    default=HERE / "data_contract/stage1_cases.json")
    ap.add_argument("--output", type=Path,
                    default=HERE / "data_contract/stage1_plan.json")
    ap.add_argument("--wall", type=float, default=45.)
    ap.add_argument("--lu", type=int, default=2000)
    ap.add_argument("--rhs", type=int, default=60000)
    a = ap.parse_args()
    cases = json.loads(a.cases.read_text())["cases"]
    if len(cases) != 5 or len({c["case_id"] for c in cases}) != 5:
        raise ValueError("Stage 1 requires exactly five preselected unique cases")
    tasks = []
    for c in cases:
        if not (HERE / c["public_path"]).exists():
            raise FileNotFoundError(HERE / c["public_path"])
        for arm in ("frozen_baseline", "adaptive_feedback"):
            tasks.append({"task_id": c["case_id"] + "__" + arm,
                          "case_id": c["case_id"], "arm": arm,
                          "public_path": c["public_path"],
                          "budget": {"wall_seconds": a.wall, "LU": a.lu,
                                     "actual_complex_RHS": a.rhs}})
    plan = {
        "schema": "fbg-r2-stage1-plan-v1", "status": "DEVELOPMENT_DISCLOSED_R4",
        "selection_rule": "fixed five clean cases recorded before new stage 1 runs",
        "task_count": len(tasks), "tasks": tasks,
        "mechanism_fields_required": ["initial_residual", "feedback_trigger",
             "feedback_update_count", "accepted_update_count", "full_physical_residual",
             "ranking_score", "all_unknowns_released", "LU", "actual_complex_RHS",
             "wall_seconds", "stop_reason"],
        "interpretation": "Paired finite cases diagnose feedback behavior; no blind or global claim.",
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(plan, indent=2) + "\n")
    print(json.dumps({"tasks": len(tasks), "output": str(a.output)}))


if __name__ == "__main__":
    main()
