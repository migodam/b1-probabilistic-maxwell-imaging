"""Generate a bounded DEVELOPMENT illumination set from disclosed R4 truth.

Default is two predeclared hard high-contrast objects, two angle designs and
four source counts. Only per-source observations are materialized. Fixed-total
energy is a global source scaling, analytically equivalent for the normalized
data objective and normalized dual, so no duplicate inversion is planned.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[1]
R4 = REPO / "research/research_z_alg_r4_diversity"
sys.path.insert(0, str(HERE))
from src.illumination_design import energy_scale, incident_field, source_positions

DEFAULT_CASES = ("dev_high_digit8", "dev_high_ring")


def load_r4_model():
    spec = importlib.util.spec_from_file_location("r4_physics_for_generation", R4 / "src/physics.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.Model


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def generate(out: Path, case_ids: set[str], seed: int) -> dict:
    if out.exists() and any(out.iterdir()):
        raise FileExistsError("Use an empty output directory; preserve prior realization")
    out.mkdir(parents=True, exist_ok=True)
    cases = json.loads((HERE / "data_contract/stage1_cases.json").read_text())["cases"]
    cases = [c for c in cases if c["case_id"] in case_ids]
    if {c["case_id"] for c in cases} != case_ids:
        raise ValueError("Unknown case ID")
    with np.load(R4 / "public_data/operators.npz", allow_pickle=False) as f:
        base = {k: f[k] for k in f.files}
    Model = load_r4_model()
    ops_dir = out / "operators"
    obs_dir = out / "public_data"
    ops_dir.mkdir(); obs_dir.mkdir()
    # Resolve material only in this local generation process. No GT enters public files.
    records = []
    counts = {"LU": 0, "actual_complex_RHS": 0}
    started = time.perf_counter()
    for design in ("uniform", "random"):
        for p in (3, 6, 12, 24):
            positions = source_positions(p, design, base["sources"], seed)
            raw_E0 = incident_field(base["centers"], positions, float(base["k"]))
            assert raw_E0.shape == (576, p)
            # The uniform p=12 construction is the old physical operator.
            if design == "uniform" and p == 12 and not np.allclose(raw_E0, base["E0"], rtol=1e-12, atol=1e-12):
                raise RuntimeError("Disk-pulse formula does not match original E0")
            op_id = f"{design}_p{p}"
            np.savez_compressed(ops_dir / f"{op_id}.npz",
                                G=base["G"], E0=raw_E0, S=base["S"],
                                centers=base["centers"], sources=positions,
                                receivers=base["receivers"], radius=base["radius"],
                                k=base["k"], grid_shape=np.array([24, 24]),
                                material_bounds=np.array([.02, 20.]))
            model = Model({"G": base["G"], "E0": raw_E0, "S": base["S"]})
            for case in cases:
                ident = Path(case["public_path"]).stem
                with np.load(R4 / "private_truth" / f"{ident}.npz", allow_pickle=False) as f:
                    material = np.asarray(f["m"], float)
                state = model.forward(material, label="development_illumination_generation")
                counts["LU"] += 1
                counts["actual_complex_RHS"] += p
                obs_id = f"{case['case_id']}__{op_id}__per_source"
                np.savez_compressed(obs_dir / f"{obs_id}.npz",
                                    Y=state["y"],
                                    operators_path=np.array(f"../operators/{op_id}.npz"),
                                    observation_id=np.array(obs_id))
                records.append({"observation_id": obs_id,
                                "public_path": f"public_data/{obs_id}.npz",
                                "operator_path": f"operators/{op_id}.npz",
                                "operators_path": f"operators/{op_id}.npz",
                                "case_id": case["case_id"], "design": design,
                                "p": p, "energy_protocol": "per_source",
                                "fixed_total_equivalent_scale": energy_scale(p, "fixed_total"),
                                "noise": "clean", "status": "DEVELOPMENT"})
    manifest = {"schema": "fbg-r2-illumination-v1", "status": "DEVELOPMENT_DISCLOSED_GT",
                "source_geometry": "physical disk-pulse positions; 3/6 nested in 12; 24 adds 12 actual sources",
                "receivers": 48, "records": records,
                "forward_solves": counts["LU"], "generation_costcounts": counts,
                "generation_wall_seconds": time.perf_counter() - started,
                "fixed_total_rule": "Analytic global scaling sqrt(12/p) of both E0 and Y; normalized data/dual objective unchanged; no duplicate optimization.",
                "source_mixing_claim": "An invertible mix of fixed sources adds no independent information; whiten mixed residuals using M inverse or covariance.",
                "truth_in_public_data": False}
    write_json(out / "manifest.json", manifest)
    return manifest


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--cases", nargs="+", default=list(DEFAULT_CASES))
    ap.add_argument("--seed", type=int, default=20261030)
    ap.add_argument("--go", action="store_true", help="Explicit developmental generation")
    a = ap.parse_args()
    if not a.go:
        ap.error("Generation requires --go")
    m = generate(a.out.resolve(), set(a.cases), a.seed)
    print(json.dumps({"records": len(m["records"]), "forward_solves": m["forward_solves"],
                      "output": str(a.out.resolve())}))


if __name__ == "__main__":
    main()
