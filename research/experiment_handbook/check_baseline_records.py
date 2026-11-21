"""Read frozen tables/source; no GT arrays, physical solves, or optimization."""
import argparse
import csv
import difflib
import json
import math
from pathlib import Path
import statistics
import time
from lookup import frozen_study


def rows(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, default=Path(".artifacts/experiment_handbook_replay"))
    args = parser.parse_args()
    root = args.root.resolve()
    r13_dir, r8_dir = frozen_study(root, 13), frozen_study(root, 8)
    srap_dir = root / "research/srap_cie_array_study_20261121"
    wall, cpu = time.perf_counter(), time.process_time()
    args.output.mkdir(parents=True, exist_ok=True)
    base = r13_dir / "supplement_aplus48/analysis/stack_comparison"
    records = rows(base / "FACTORIAL_SCORES.csv")
    summary = rows(base / "METHOD_VIEW_SUMMARY.csv")
    r13 = {}
    for solver in ("cie", "dual_cie", "g19", "dr", "a_plus_mm_r13", "a_plus_arc_r7"):
        selected = [r for r in records if (r["solver"], r["initialization"], r["view"]) == (solver, "uniform", "primary")]
        assert len(selected) == len({r["case_id"] for r in selected}) == 48
        mean = statistics.mean(float(r["NRMSE"]) for r in selected)
        ref = next(r for r in summary if (r["cohort"], r["solver"], r["view"]) == ("uniform48", solver, "primary"))
        assert abs(mean - float(ref["NRMSE_mean"])) < 1e-14
        r13[solver] = dict(objects=48, mean_NRMSE=mean,
                           mean_external_child_CPU_s=statistics.mean(float(r["external_child_CPU_seconds"]) for r in selected),
                           uniform16_mean_NRMSE=statistics.mean(float(r["NRMSE"]) for r in selected if int(r["case_id"].split("_")[-1]) % 3 == 0),
                           complement32_mean_NRMSE=statistics.mean(float(r["NRMSE"]) for r in selected if int(r["case_id"].split("_")[-1]) % 3 != 0))
    common = rows(base / "COMMON_ACTUAL_CPU_SIX_AND_TWELVE_SUMMARY.csv")
    for solver, rec in r13.items():
        row = next(r for r in common if (r["cohort"], r["solver"], r["view"]) == ("uniform48", solver, "common_actual_CPU_six"))
        assert row["available_cases"] == "48"
        rec["paid_common_CPU_prefix_NRMSE"] = float(row["NRMSE_mean"])
    u11 = {}
    old = rows(srap_dir / "engine_dev/hardware/historical_gain_11.csv")
    for host in ("Mac", "XINAN"):
        for view in ("raw", "light"):
            rs = [r for r in old if (r["host"], r["view"]) == (host, view)]
            assert len(rs) == len({r["case"] for r in rs}) == 11
            gains = [(min(float(r["paper_nie"]), float(r["cie_dual"])) - float(r["u11_routed_ap"])) / min(float(r["paper_nie"]), float(r["cie_dual"])) for r in rs]
            assert max(abs(g - float(r["gain"])) for g, r in zip(gains, rs)) < 1e-14
            u11[f"{host}_{view}"] = dict(objects=11, wins=sum(g > 0 for g in gains), ties=sum(g == 0 for g in gains), losses=sum(g < 0 for g in gains), median_fractional_gain=statistics.median(gains))
    mac = (srap_dir / "engine_dev/hardware/Mac_20261120_runner.py").read_text()
    xinan = (srap_dir / "engine_dev/hardware/XINAN_20261120_runner.py").read_text()
    diff = "".join(difflib.unified_diff(mac.splitlines(True), xinan.splitlines(True), fromfile="Mac_20261120_runner.py", tofile="XINAN_20261120_runner.py"))
    (args.output / "u11_frozen_runner_diff.txt").write_text(diff, encoding="utf-8")
    contract = json.loads((r8_dir / "delegated/paper_baselines/PAPER_CONTRACTS.json").read_text())
    paper_source = (r8_dir / "experiments/paper_benchmark.py").read_text()
    foam_source = (r8_dir / "experiments/fresnel_foam_invert.py").read_text()
    stop = dict(NIE2016_formula=contract["NIE2016"]["stop_formula"],
                NIE2016_root_thresholds=contract["NIE2016"]["stop_delta2D"],
                NIE2016_equivalent_squared_thresholds=[x*x for x in contract["NIE2016"]["stop_delta2D"]],
                CIE2019_frozen_metadata_squared_threshold=contract["CIE2019"]["coefficient_relative_change_squared_stop"],
                CIE2019_primary_Eq16_root_threshold=1e-3,
                CIE2019_primary_equivalent_squared_threshold=1e-6,
                CIE2019_historical_metadata_conflicts_with_primary=True,
                CIE2019_historical_squared_1e_minus3_equivalent_root_threshold=math.sqrt(1e-3),
                primary_reference="https://www.mdpi.com/2313-433X/5/2/27 ; section 3 Eq(16), directly opened 2026-11-21",
                paper_benchmark_explicit_squared_conversion="coefficient_tol=(.01 if j==0 else .001)**2" in paper_source,
                final_Foam_explicit_squared_conversion="coefficient_tol=1e-4 if mf==6 else 1e-6" in foam_source,
                caveat="Global Frobenius and mean sourcewise ratios differ; sqrt conversion alone does not equate their trajectories.")
    assert stop["paper_benchmark_explicit_squared_conversion"] and stop["final_Foam_explicit_squared_conversion"]
    result = dict(status="SAVED_TABLE_AND_SOURCE_REVIEW", r13=r13, u11=u11, stopping=stop,
                  references=[str(p.relative_to(root)) for p in (base / "FACTORIAL_SCORES.csv", base / "METHOD_VIEW_SUMMARY.csv",
                              base / "COMMON_ACTUAL_CPU_SIX_AND_TWELVE_SUMMARY.csv",
                              srap_dir / "engine_dev/hardware/historical_gain_11.csv",
                              r8_dir / "delegated/paper_baselines/PAPER_CONTRACTS.json")],
                  exposure="Existing published results only; no new independent confirmation", truth_arrays_read=0,
                  new_inversions=0, new_Maxwell_LU=0, new_physical_RHS=0,
                  previous_review_passes=[],
                  cost_scope="CSV/source reads and arithmetic; process startup, tool/research/writing overhead uninstrumented",
                  wall_s=time.perf_counter()-wall, cpu_s=time.process_time()-cpu)
    (args.output / "baseline_record_checks.json").write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"r13":r13,"u11":u11,"wall_s":result["wall_s"],"cpu_s":result["cpu_s"]},ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
