"""Final TEST endpoint success available by measured total cost; score tables only.

This is a completion-cost distribution, not an optimizer first-hit trajectory.
No truth, prediction arrays, candidate oracle or physics implementation is read.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

STUDY = Path(__file__).resolve().parents[1]
DENOMINATOR = 108
MAIN_METHODS = (
    ("strict_gn45_anchor", "Strict GN"),
    ("continue_trf", "Continue TRF"),
    ("discrepancy_allscales45_anchor", "All-scales discrepancy"),
    ("cie_dual_cold55", "CIE dual cold"),
    ("v_physics_seed7_cascade2", "Voxel s7"),
    ("j_full_seed7_cascade2", "J-full s7"),
    ("g_weak_seed7_cascade2", "G-weak s7"),
    ("g_weak_seed7_cascade2_const", "G-weak + const s7"),
    ("recurrent_feedback1_seed7_cascade2", "Recurrent feedback1 s7"),
)
COSTS = (
    ("online_total_wall_seconds", "Actual online total wall seconds"),
    ("online_total_LU", "Actual online total LU"),
    ("online_total_RHS", "Actual online total complex RHS"),
)


def flag(value):
    return value is True or str(value).lower() in ("true", "1")


def number(value):
    if value is None or value == "":
        return None
    value = float(value)
    return value if math.isfinite(value) and value >= 0 else None


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def write_csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_scored_table(path):
    """JSON score rows or flattened scorer CSV with its provenance JSON."""
    metadata_path = path.with_suffix(".json") if path.suffix.lower() == ".csv" else path
    metadata = json.loads(metadata_path.read_text())
    if (metadata.get("schema") != "gn-nn-r4-scored-v1" or
            metadata.get("truth_scope") != "NEW_TEST_AFTER_METHOD_AND_PREDICTION_FREEZE"):
        raise ValueError("Use final post-freeze TEST scorer output, not VAL or intermediate tables")
    if path.suffix.lower() == ".csv":
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
    else:
        rows = metadata["rows"]
    if not rows or any(row.get("split") != "test" for row in rows):
        raise ValueError("Require nonempty TEST-only score rows")
    return rows, metadata


def endpoint(row):
    """Use frozen top1 and score flags; never search candidates or recompute GT."""
    if row is None:
        return dict(status="MISSING_SCORE_ROW", present=False, endpoint=False,
                    quality=False, accounted=False, within55=False, eligible=False,
                    costs={key: None for key, _ in COSTS})
    top = row.get("top1") if isinstance(row.get("top1"), dict) else None
    has_endpoint = top is not None or number(row.get("selected_NRMSE")) is not None
    status = row.get("status", "NOT_RECORDED")
    quality = bool(has_endpoint and status != "NOT_RUN" and flag(row.get("joint_success")))
    costs = {key: number(row.get(key)) for key, _ in COSTS}
    accounted = flag(row.get("accounted_within_budget"))
    within55 = flag(row.get("online_total_within55"))
    # Keep both explicit flags and measured totals: no overrun becomes a
    # budget-qualified success even if a stale/incomplete flag says otherwise.
    measured = all(value is not None for value in costs.values())
    caps_met = bool(measured and costs["online_total_wall_seconds"] <= 55
                    and costs["online_total_LU"] <= 600 and costs["online_total_RHS"] <= 36000)
    eligible = bool(quality and accounted and within55 and caps_met)
    return dict(status=status, present=True, endpoint=has_endpoint, quality=quality,
                accounted=accounted, within55=within55, eligible=eligible, costs=costs)


def prepare(rows, registered_runs):
    requested = [run for run, _ in MAIN_METHODS]
    absent = [run for run in requested if run not in registered_runs]
    if absent:
        raise ValueError("Requested main methods are absent from FINAL_RUNS: "+", ".join(absent))
    identities = {0.: {}, .001: {}}
    by_key = {}
    for row in rows:
        noise = number(row.get("noise"))
        if noise not in identities:
            raise ValueError("Expected separate clean=0 and noisy=.001 score strata")
        obs, parent = row["obs_id"], row["parent_id"]
        known = identities[noise].setdefault(parent, obs)
        if known != obs:
            raise ValueError("More than one observation per mother inside a noise stratum")
        key = row["run"], parent, noise
        if key in by_key:
            raise ValueError("Duplicate run/mother/noise score row")
        by_key[key] = row
    if any(len(ids) != DENOMINATOR for ids in identities.values()):
        raise ValueError("Final table must identify all108 planned TEST mothers in each noise stratum")
    if set(identities[0.]) != set(identities[.001]):
        raise ValueError("Clean/noisy mother identities must agree; their results remain separate")
    tables, methodcounts, curves = {}, [], []
    for noise, ids in identities.items():
        for run, label in MAIN_METHODS:
            endpoints = [endpoint(by_key.get((run, parent, noise))) for parent in sorted(ids)]
            tables[run, noise] = endpoints
            counts = dict(run=run, label=label, noise=noise, planned_mothers=DENOMINATOR,
                          score_rows_present=sum(p["present"] for p in endpoints),
                          missing_score_rows=sum(not p["present"] for p in endpoints),
                          NOT_RUN=sum(p["status"] == "NOT_RUN" for p in endpoints),
                          selected_final_endpoints=sum(p["endpoint"] for p in endpoints),
                          no_final_endpoint=sum(not p["endpoint"] for p in endpoints),
                          quality_only_joint_endpoints=sum(p["quality"] for p in endpoints),
                          budget_qualified_joint_endpoints=sum(p["eligible"] for p in endpoints),
                          quality_joint_but_not_budget_qualified=sum(p["quality"] and not p["eligible"] for p in endpoints),
                          quality_joint_method_accounted_false=sum(p["quality"] and not p["accounted"] for p in endpoints),
                          quality_joint_total55_false=sum(p["quality"] and not p["within55"] for p in endpoints),
                          final_quality_rate=sum(p["quality"] for p in endpoints)/DENOMINATOR,
                          final_budget_qualified_rate=sum(p["eligible"] for p in endpoints)/DENOMINATOR)
            for cost, _ in COSTS:
                counts[cost+"_measured"] = sum(p["costs"][cost] is not None for p in endpoints)
                counts[cost+"_unknown"] = sum(p["costs"][cost] is None for p in endpoints)
            methodcounts.append(counts)
        for cost, _ in COSTS:
            observed = [p["costs"][cost] for run, _ in MAIN_METHODS
                        for p in tables[run, noise] if p["costs"][cost] is not None]
            maximum = max(observed, default=0.)
            for run, label in MAIN_METHODS:
                successful_costs = [p["costs"][cost] for p in tables[run, noise] if p["eligible"]]
                # Only observed successful final-endpoint costs create jumps.
                thresholds = sorted({0., maximum, *successful_costs})
                for threshold in thresholds:
                    available = sum(value <= threshold for value in successful_costs)
                    curves.append(dict(run=run, label=label, noise=noise, cost_metric=cost,
                                       threshold_B=threshold, available_budget_qualified_final_joint_success=available,
                                       planned_mothers=DENOMINATOR, success_fraction=available/DENOMINATOR))
    return tables, methodcounts, curves


def draw(curves, tables, destination):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))
    colors = plt.get_cmap("tab10")
    handles = []
    for i, noise in enumerate((0., .001)):
        for j, (cost, xlabel) in enumerate(COSTS):
            ax = axes[i, j]
            for k, (run, label) in enumerate(MAIN_METHODS):
                points = [row for row in curves if row["run"] == run and row["noise"] == noise and row["cost_metric"] == cost]
                line, = ax.step([row["threshold_B"] for row in points], [row["success_fraction"] for row in points],
                                where="post", color=colors(k), lw=1.6,
                                linestyle="--" if k < 4 else "-", label=label)
                if i == 0 and j == 0:
                    handles.append(line)
            ax.set(xlabel=xlabel+" threshold B", ylabel="Successful final top1 available by B /108",
                   title=("Clean" if noise == 0 else "Noisy .001")+" — 108 mothers", ylim=(0, 1.02))
            ax.yaxis.set_major_formatter(PercentFormatter(1.))
            ax.grid(alpha=.2)
            measured = sum(p["costs"][cost] is not None for run, _ in MAIN_METHODS for p in tables[run, noise])
            ax.text(.02, .97, f"Measured endpoint costs: {measured}/{len(MAIN_METHODS)*DENOMINATOR}",
                    va="top", fontsize=8, transform=ax.transAxes)
    fig.suptitle("Final online top1 completion cost — joint success AND both accounting gates\n"
                 "Each curve is an endpoint availability distribution; no inferred first-hit history, no best@K.", fontsize=13)
    fig.legend(handles, [label for _, label in MAIN_METHODS], loc="lower center", bbox_to_anchor=(.5, .065), ncol=3, fontsize=9)
    fig.text(.5, .025, "Missing / NOT_RUN / unsuccessful endpoints remain in each 108-mother denominator. Common GN preparation is already in scored totals.\n"
             "Mixed CPU / MPS / i9 hardware: wall time is deployment-specific. LU/RHS exclude NN matrix products; training and setup are billed separately.",
             ha="center", fontsize=9)
    fig.subplots_adjust(left=.07, right=.985, bottom=.21, top=.88, hspace=.33, wspace=.27)
    fig.savefig(destination, dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=STUDY/"results/TEST_SCORES.json")
    parser.add_argument("--plan", type=Path, default=STUDY/"audit/FINAL_RUNS.json")
    parser.add_argument("--output", type=Path, default=STUDY/"figures/generalization")
    args = parser.parse_args();args.output.mkdir(parents=True, exist_ok=True)
    receipt = args.output/"COMPLETION_COST_STATUS.json"
    if not args.scores.exists():
        write_json(receipt, dict(status="NOT_RUN_NO_FINAL_TEST_SCORES", expected_source=str(args.scores), figures=[], truth_arrays_read=False))
        print("NOT_RUN_NO_FINAL_TEST_SCORES; no scientific figure generated")
        return
    rows, metadata = read_scored_table(args.scores)
    plan = json.loads(args.plan.read_text())
    tables, methodcounts, curves = prepare(rows, plan["runs"])
    name = "completion_success_vs_actual_cost.png"
    write_csv(args.output/"completion_success_vs_actual_cost_CURVES.csv", curves)
    write_csv(args.output/"completion_success_vs_actual_cost_METHODCOUNTS.csv", methodcounts)
    draw(curves, tables, args.output/name)
    write_json(receipt, dict(status="GENERATED_FROM_FINAL_SCALAR_SCORE_TABLE", source=str(args.scores),
        seal_path=metadata.get("seal_path"), plan=str(args.plan), figure=name,
        denominator_per_method_per_noise=DENOMINATOR, methods=[run for run, _ in MAIN_METHODS],
        semantics="fraction of108 whose frozen FINAL online top1 is jointly successful, complete/in method budget, within measured total55 envelope, and has actual total cost <=B",
        total55_envelope=dict(wall_seconds=55, LU=600, complex_RHS=36000),
        gates="joint_success AND accounted_within_budget AND online_total_within55 AND all measured total costs within55/600/36000",
        quality_only_counts="retained separately in METHODCOUNTS.csv; do not appear on budget-qualified curves",
        noise_or_seed_pooling=False, first_hit_trajectory_inferred=False, oracle_bestK_used=False,
        common_GN_cost_added_again=False, budget_cap_substituted_for_actual_cost=False,
        unknown_cost_policy="unknown costs exclude success credit but preserve all108 denominator entries",
        missing_row_policy="registered main-plan method/mother/noise missing rows are explicitly counted as MISSING_SCORE_ROW failures; no completed coverage is invented",
        hardware_scope="mixed CPU/MPS/i9 deployment times; LU/RHS omit NN matrix products; one-time training/setup separately billed",
        truth_arrays_read=False, physics_executed=False))
    print(json.dumps(dict(status="GENERATED_FROM_FINAL_SCALAR_SCORE_TABLE", figure=str(args.output/name))))


if __name__ == "__main__":
    main()
