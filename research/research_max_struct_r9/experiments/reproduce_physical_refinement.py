"""Isolated wrapper around the executed three-start physical experiment.

No physics, budget, source, or candidate-selection changes. Historical runner
and scoring artifacts stay untouched. Run on the documented Windows CUDA host;
score may run locally after the entire output directory has been transferred.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("run", "score"))
    parser.add_argument("--output-dir", default="reproduction/physical_refinement")
    args = parser.parse_args()
    out = (ROOT / args.output_dir).resolve()
    out.relative_to((ROOT / "reproduction").resolve())
    if args.action == "run":
        spec = importlib.util.spec_from_file_location(
            "executed_refinement", ROOT / "experiments/run_original_physical_refinement.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.OUT = out
        module.main()
        return

    cpu, wall = time.process_time(), time.perf_counter()
    import numpy as np
    seal = json.loads((out / "SEALED_ALL_THREE.json").read_text())
    if len(seal["entries"]) != 3:
        raise ValueError("All three completed predictions are required before scoring")
    for entry in seal["entries"]:
        for record in entry["files"]:
            path = out / entry["start"] / Path(record["path"].replace("\\", "/")).name
            if path.stat().st_size != record["bytes"]:
                raise ValueError("Incomplete sealed prediction transfer")
    target = out / "REPRODUCTION_SCORES.json"
    if target.exists():
        raise FileExistsError(target)
    # First private material read, after all three sealed predictions exist.
    with np.load(ROOT / "RESULTS/original_vector3d/sealed_truth/truth.npz") as data:
        truth = data["inverse_material"]
    with np.load(ROOT / "RESULTS/original_vector3d/cohorts/matched_grid/public/observations.npz") as data:
        observed = data["Y"]
    rows = []
    for entry in seal["entries"]:
        name = entry["start"]
        with np.load(out / name / "prediction.npz") as data:
            material, prediction = data["material"], data["prediction"]
        receipt = json.loads((out / name / "RECEIPT.json").read_text())
        rows.append(dict(start=name,
            contrast_relative_L2=float(np.linalg.norm(material-truth)/np.linalg.norm(truth)),
            paper_Err=float(np.sqrt(np.mean((material-truth)**2/(1+truth)**2))),
            Phi=float(np.mean(np.sum(abs(prediction-observed)**2,axis=0)/np.sum(abs(observed)**2,axis=0))),
            stop=receipt["stop"], selected_state_id=receipt["selected_state_id"],
            CPU_seconds=receipt["CPU_seconds"], wall_seconds=receipt["wall_seconds"]))
    target.write_text(json.dumps(dict(
        classification="REPRODUCTION_OF_PUBLIC_DEVELOPMENT_NOT_NEW_CONFIRMATION",
        selection="original runner public Phi; no truth selection",
        rows=rows, scoring_CPU_seconds=time.process_time()-cpu,
        scoring_wall_seconds=time.perf_counter()-wall,
        scoring_physical_RHS=0), indent=2)+"\n")
    print(target)


if __name__ == "__main__":
    main()
