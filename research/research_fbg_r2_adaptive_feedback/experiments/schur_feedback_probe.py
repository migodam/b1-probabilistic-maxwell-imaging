"""Q9 exploratory Schur feedback runner; the frozen adaptive core is untouched.

Usage: python -m research.research_fbg_r2_adaptive_feedback.experiments.schur_feedback_probe MANIFEST
The manifest is passed unchanged to adaptive_solver.main. Only the named
``feedback_schur_diag`` arm is intercepted; every other arm uses its original
constructor. The script performs no GPU work at import time.
"""
from pathlib import Path
import sys

PROJECT = Path(__file__).resolve().parents[1]
if str(PROJECT) not in sys.path:
    sys.path.insert(0, str(PROJECT))

from src import adaptive_solver
from src.schur_feedback_lane import make_schur_feedback


def main():
    original = adaptive_solver.make_feedback

    def patched(lane, x, kind, cfg):
        if kind == "schur_diag":
            return make_schur_feedback(lane, x, kind, cfg)
        return original(lane, x, kind, cfg)

    adaptive_solver.make_feedback = patched
    try:
        adaptive_solver.main()
    finally:
        adaptive_solver.make_feedback = original


if __name__ == "__main__":
    main()
