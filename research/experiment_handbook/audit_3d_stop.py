"""Trace a saved R8 3-D Fourier PR run; no solver imports or new numerical run."""
import argparse
import ast
import json
from pathlib import Path
import time
from lookup import frozen_study


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, default=Path(".artifacts/experiment_handbook_replay"))
    args = parser.parse_args()
    wall, cpu = time.perf_counter(), time.process_time()
    root = args.root.resolve()
    r8 = frozen_study(root, 8)
    src = r8 / "src/vector_pr3d.py"
    tree = ast.parse(src.read_text())
    fun = next(x for x in ast.walk(tree) if isinstance(x, ast.FunctionDef) and x.name == "optimize_pr")
    # Verify both the fixed-cycle loop and absence of a coefficient-based early exit.
    loop = next(x for x in ast.walk(fun) if isinstance(x, ast.For) and ast.unparse(x.target) == "it")
    assert ast.unparse(loop.iter) == "range(outer)"
    assert not any(isinstance(x, ast.Break) for x in ast.walk(loop))
    assert not any(isinstance(x, ast.Name) and x.id in ("delta", "delta2", "coefficient_tol") for x in ast.walk(fun))
    folder = r8 / "RESULTS/vector3d/pr240_fourier_coordinates"
    config = json.loads((folder / "PR_FROZEN_CONFIG.json").read_text())
    assert config["outer"] == 240
    runs = {}
    for arm in ("CIE_sourcewise", "Dual_sourcewise"):
        path = folder / "four_cubes" / (arm + "_receipt.json")
        r = json.loads(path.read_text())
        history = [q for q in r["optimizer"] if q["stage"] == "joint"]
        assert [q["cycle"] for q in history] == list(range(240))
        assert [sum(q["mf"] == mf for q in history) for mf in (1, 2)] == [120, 120]
        assert r["counts"]["current_PR_steps"] == 240
        runs[arm] = dict(receipt=str(path.relative_to(root)), completed_current_updates=240,
                         completed_updates_by_stage=[120,120], coefficient_stop_evaluated=False,
                         corrector_stop=r["optimizer"][-1], historical_wall_s=r["wall_s"], historical_cpu_s=r["cpu_s"],
                         historical_physical_LU=r["counts"]["physical_LU"],
                         historical_physical_RHS=r["counts"]["source_RHS"]+r["counts"]["adjoint_RHS"])
    srap_path = root / "research/srap_cie_array_study_20261121/vector/invert.py"
    srap = srap_path.read_text()
    assert "change=(normcols(update)/normcols(g).clamp_min(1e-30)).mean().sqrt().item()" in srap
    assert "change<.001" in srap
    rec = dict(status="FIXED_CYCLE_3D_PATH_NOT_TRUNCATED_BY_MISREAD_THRESHOLD",
               chain=[str((r8/"experiments/vector_pr240_fourier.py").relative_to(root)), str(src.relative_to(root)), "VectorPR3D.optimize_pr", "range(outer=240), no coefficient stop"],
               runs=runs, saved_gamma_change_history=None,
               saved_gamma_change_history_reason="No per-update gamma difference is computed or serialized by this path; material inner optimizer convergence is a separate quantity.",
               primary_paper_equation="https://www.mdpi.com/2313-433X/5/2/27 ; Eq16 sqrt(mean per-source squared relative coefficient changes)<1e-3",
               interpretation="Frozen metadata is incorrect about the paper's squared threshold; this fixed-cycle run is an intentional finite schedule. No evidence predicts the material result of introducing paper termination.",
               SRAP_vector_stop=dict(source=str(srap_path.relative_to(root)), root_statistic=True, threshold=.001, affected_by_metadata_error=False),
               new_inversions=0,new_Maxwell_LU=0,new_physical_RHS=0,truth_arrays_read=0,
               wall_s=time.perf_counter()-wall,cpu_s=time.process_time()-cpu,
               cost_scope="AST/source and saved receipt reads; historical solver costs above are references, not new costs; startup and research overhead uninstrumented")
    out = args.output / "three_dimensional_stop_audit.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(rec,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print(json.dumps({k:rec[k] for k in ("status","wall_s","cpu_s","new_Maxwell_LU")},ensure_ascii=False))


if __name__ == "__main__":
    main()
