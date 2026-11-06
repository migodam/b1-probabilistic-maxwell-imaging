"""Run one parent-frozen R9 public task through the inherited R8 wrapper.
No launch occurs on import. A fresh output directory is mandatory.
"""
import os
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[key] = "1"
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from baseline_adapter import execute, RULES


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True,
                        help="R9 manifest with an explicit public task and solver config")
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--config", help="Frozen config override; relative paths use the R9 repository root")
    parser.add_argument("--output-dir", help="Fresh output, confined to the R9 study")
    parser.add_argument("--r8-repo", help="Actual inherited checkout; default is sibling r8_worktree")
    parser.add_argument("--dpcs-rule", choices=RULES,
                        help="Default task.dpcs_rule, otherwise fixed_001")
    parser.add_argument("--prepare-only", action="store_true",
                        help="Validate public inputs and freeze adapter inputs; no inverse solve")
    args = parser.parse_args()
    receipt, out = execute(args)
    print(json.dumps(dict(task_id=args.task_id, status=receipt["status"], output_dir=str(out),
        truth_access="NONE", scoring="NOT_RUN",
        r9_total_cpu_seconds=receipt.get("r9_total_cpu_seconds", receipt.get("cpu_seconds")),
        r9_total_wall_seconds=receipt.get("r9_total_wall_seconds", receipt.get("wall_seconds")))), flush=True)
    return 1 if receipt["status"].startswith("FAILED") else 0


if __name__ == "__main__":
    raise SystemExit(main())
