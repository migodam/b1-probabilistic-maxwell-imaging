"""Run independent checks and an unchanged source copy; retain CPU receipts."""
import importlib.util
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
from research_w import analytic
import numpy as np


def run_logged(name, command):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    wall = time.perf_counter()
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                               env={**os.environ, "PYTHONPATH": str(ROOT / "src")})
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    elapsed = time.perf_counter() - wall
    (HERE / f"{name}.stdout.txt").write_text(completed.stdout)
    (HERE / f"{name}.stderr.txt").write_text(completed.stderr)
    return {"argv": command, "returncode": completed.returncode,
            "wall_seconds": elapsed,
            "child_user_cpu_seconds": after.ru_utime - before.ru_utime,
            "child_system_cpu_seconds": after.ru_stime - before.ru_stime,
            "child_process_cpu_seconds": after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime}


def max_scalar_difference(first, second):
    differences = {key: abs(first["tests"][key] - value) for key, value in second["tests"].items()}
    sweep = [{key: abs(row[key] - ref[key]) for key in ref}
             for row, ref in zip(first["bandwidth_sweep"], second["bandwidth_sweep"], strict=True)]
    return {"tests_absolute_differences": differences, "sweep_absolute_differences": sweep,
            "max_test_difference": max(differences.values()),
            "max_sweep_difference": max(value for row in sweep for value in row.values())}


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    rerun_dir = HERE / "package_rerun"
    rerun_dir.mkdir(exist_ok=True)
    copied = rerun_dir / "analytic_witness.py"
    shutil.copyfile(ROOT / "source_package" / "analytic_witness.py", copied)
    receipts = {}
    receipts["independent_regression"] = run_logged("independent", [sys.executable, "-m", "research_w.analytic",
        "--output", str(ROOT / "results" / "p0a"), "--source-results",
        str(ROOT / "source_package" / "analytic_results" / "results.json")])
    receipts["unit_tests"] = run_logged("tests", [sys.executable, "-m", "unittest", "discover", "-s", "tests",
                                              "-p", "test_analytic.py", "-v"])
    receipts["unchanged_package_copy"] = run_logged("package_rerun", [sys.executable, str(copied)])
    if any(receipt["returncode"] != 0 for receipt in receipts.values()):
        (HERE / "execution_receipt.json").write_text(json.dumps(receipts, indent=2) + "\n")
        raise RuntimeError("Execution failed; inspect retained logs")
    original = json.loads((ROOT / "source_package" / "analytic_results" / "results.json").read_text())
    fresh = json.loads((rerun_dir / "analytic_results" / "results.json").read_text())
    independent = json.loads((ROOT / "results" / "p0a" / "results.json").read_text())
    spec = importlib.util.spec_from_file_location("copied_witness", copied)
    witness = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(witness)
    array_differences = []
    dense = np.linspace(0.8, 1.2, 1001)
    for label, layers in (("wall1", analytic.WALL1), ("wall2", analytic.WALL2),
                          ("asymmetric", ((1.0, 0.1), (2.3, 0.12), (1.4, 0.08)))):
        fields = (("r_left", analytic.scattering(layers, dense)[0], witness.scattering(layers, dense)[0]),
                  ("r_right", analytic.scattering(layers, dense)[1], witness.scattering(layers, dense)[1]),
                  ("t", analytic.scattering(layers, dense)[2], witness.scattering(layers, dense)[2]),
                  ("fixed_plane_t", analytic.fixed_plane_transmission(layers, dense),
                   witness.fixed_plane_transmission(layers, dense)),
                  ("pec_y", analytic.full_target_reflection(layers, dense, 1.5),
                   witness.full_target_reflection(layers, dense, 1.5)))
        for quantity, ours, theirs in fields:
            array_differences.append({"configuration": label, "quantity": quantity,
                                      "frequency_count": len(dense),
                                      "maximum_absolute_difference": float(np.max(abs(ours - theirs)))})
    comparison = {"archived_vs_fresh_package": max_scalar_difference(fresh, original),
                  "independent_vs_fresh_package": max_scalar_difference(independent, fresh),
                  "dense_complex_array_comparisons": array_differences,
                  "maximum_dense_array_difference": max(row["maximum_absolute_difference"] for row in array_differences),
                  "source_package_modified": False,
                  "scope": "1D regression reproducibility only; no scientific novelty or 2D verdict"}
    (ROOT / "results" / "p0a" / "package_comparison.json").write_text(json.dumps(comparison, indent=2) + "\n")
    receipts["runtime"] = {"python": sys.executable, "numpy": np.__version__, "device": "CPU",
        "threads": {key: os.environ.get(key) for key in
                   ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
        "new_2d_or_3d_fullwave_calls": 0, "source_copy_provenance": "shutil.copyfile, unchanged bytes, no hash computed"}
    (HERE / "execution_receipt.json").write_text(json.dumps(receipts, indent=2) + "\n")
    print(json.dumps({"receipts": receipts, "maximum_dense_array_difference": comparison["maximum_dense_array_difference"],
                      "maximum_sweep_difference": comparison["independent_vs_fresh_package"]["max_sweep_difference"]}, indent=2))


if __name__ == "__main__":
    main()
