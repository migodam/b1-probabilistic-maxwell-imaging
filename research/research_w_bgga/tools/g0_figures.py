"""Plot stored G0 geometry qualification results; never assemble or solve a model.

Run through tools/run_job.py so the parent wrapper owns CPU accounting.  The
input RESULTS.json is authoritative.  Missing cases remain NOT_RUN, including
when an old NPZ happens to exist.  This program does not update any gate or
the root report/cost files.
"""
from pathlib import Path
import argparse
import csv
import json
import math
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "research/delegated/results_support"
DEFAULT_RESPONSE_CASES = ("weak_center", "strong_center", "low_loss_high_contrast")
COLORS = {"fixed": "#166e87", "fitted": "#c45842", "missing": "#a2afb8"}
SCOPE_NOTE = ("Geometry prerequisite only. Discrete comparison is not a continuum, "
              "posterior, speed, or resonance certificate.")


def read_json(path):
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def dump_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def relative(path):
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value) if math.isfinite(value) else None


def nonnegative(value, name):
    value = number(value)
    if value is None or value < 0:
        raise ValueError(f"Invalid {name}: a finite nonnegative number is required")
    return value


def save_figure(fig, directory, stem, csv_path):
    paths = []
    for suffix in ("png", "pdf"):
        path = directory / f"{stem}.{suffix}"
        options = {"dpi": 180} if suffix == "png" else {}
        fig.savefig(path, bbox_inches="tight", **options)
        paths.append(relative(path))
    plt.close(fig)
    return {"status": "CREATED", "figures": paths, "data_csv": relative(csv_path)}


def result_table(results_path, specifications, limit):
    source_status = "NOT_RUN_RESULTS_ABSENT"
    records = []
    if results_path.exists():
        records = read_json(results_path)
        if not isinstance(records, list):
            raise ValueError("RESULTS.json must contain the measured case list")
        source_status = "RESULTS_PRESENT" if records else "NOT_RUN_RESULTS_EMPTY"
    expected = {case["id"] for case in specifications}
    measured = {}
    for row in records:
        if not isinstance(row, dict) or row.get("case") not in expected:
            raise ValueError("RESULTS contains an unregistered or invalid case")
        if row["case"] in measured:
            raise ValueError(f"Duplicate result for {row['case']}")
        measured[row["case"]] = row
    rows = []
    for case in specifications:
        source = measured.get(case["id"], {})
        value = number(source.get("fixed_vs_interface_fitted_white"))
        status = ("NOT_RUN" if not source else "MISSING_METRIC" if value is None
                  else "WITHIN_FROZEN_LIMIT" if value <= limit else "OVER_FROZEN_LIMIT")
        backward = source.get("algebraic_backward", {})
        rows.append(dict(case=case["id"], status=status,
                         fixed_vs_interface_fitted_white=value,
                         frozen_independent_whitened_max=limit,
                         response_relative=number(source.get("response_relative")),
                         ndofs=number(source.get("ndofs")),
                         reference_ndofs=number(source.get("reference_ndofs")),
                         public_target_ports=number(source.get("target_ports")),
                         touched_ports=number(source.get("touched_ports")),
                         algebraic_backward=number(backward.get("backward")),
                         polygon_area_error=number(source.get("polygon_area_error")),
                         reported_numerical_qualification=source.get("numerical_qualification"),
                         radius=case["radius"], contrast=case["contrast"],
                         target_loss=case["target_loss"], evidence=source.get("evidence", "NOT_RUN")))
    return source_status, rows, measured


def error_figure(rows, output, limit):
    fields = ["case", "status", "fixed_vs_interface_fitted_white",
              "frozen_independent_whitened_max", "response_relative",
              "algebraic_backward", "reported_numerical_qualification", "evidence"]
    csv_path = output / "G0_FIXED_VS_FITTED_WHITE.csv"
    write_csv(csv_path, rows, fields)
    fig, ax = plt.subplots(figsize=(14, 5.7))
    for index, row in enumerate(rows):
        value = row["fixed_vs_interface_fitted_white"]
        if value is None:
            ax.text(index, limit * .025, row["status"], rotation=90, ha="center", va="bottom",
                    color="#64717a", fontsize=8)
        else:
            ax.bar(index, value, color=COLORS["fixed"] if value <= limit else COLORS["fitted"])
            ax.text(index, value + limit * .018, f"{value:.4g}", ha="center", fontsize=8)
    ax.axhline(limit, color="#9e3b2e", linestyle="--", linewidth=1.4,
               label=f"Frozen independent whitened limit = {limit:g}")
    values = [row["fixed_vs_interface_fitted_white"] for row in rows
              if row["fixed_vs_interface_fitted_white"] is not None]
    ax.set(xticks=range(len(rows)), xticklabels=[row["case"].replace("_", "\n") for row in rows],
           ylabel="Stored whitened response discrepancy (fixed vs interface fitted)",
           ylim=(0, max([limit] + values) * 1.22),
           title="G0 predeclared geometry panel: fixed mesh versus interface-fitted reference")
    ax.tick_params(axis="x", labelsize=8)
    ax.grid(axis="y", alpha=.2)
    ax.legend(loc="upper left", fontsize=9)
    fig.text(.5, .015, SCOPE_NOTE, ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .06, 1, 1))
    return save_figure(fig, output, "G0_FIXED_VS_FITTED_WHITE", csv_path)


def dimension_figure(rows, output):
    csv_path = output / "G0_DOFS_AND_PORTS.csv"
    fields = ["case", "status", "ndofs", "reference_ndofs", "public_target_ports", "touched_ports",
              "radius", "contrast", "target_loss"]
    write_csv(csv_path, rows, fields)
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    for ax, keys, title, ylabel in (
            (axes[0], ("ndofs", "reference_ndofs"), "Actual full-system dimensions", "FE degrees of freedom"),
            (axes[1], ("public_target_ports", "touched_ports"), "Public target ROI and case-touched support",
             "Structural FE coordinates (log scale)")):
        for offset, key, color, label in zip((-.18, .18), keys,
                (COLORS["fixed"], COLORS["fitted"]),
                (("Fixed mesh", "Interface-fitted reference") if ax is axes[0]
                 else ("Entire public target ROI", "Actually touched by this target"))):
            labelled = False
            for index, row in enumerate(rows):
                value = row[key]
                if value is not None:
                    ax.bar(index + offset, value, .35, color=color, label=label if not labelled else None)
                    labelled = True
        ax.set(title=title, ylabel=ylabel)
        ax.grid(axis="y", alpha=.2)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labels, loc="upper right", fontsize=9)
        if ax is axes[1] and any(row[key] is not None and row[key] > 0 for row in rows for key in keys):
            ax.set_yscale("log")
        if not any(row[key] is not None for row in rows for key in keys):
            ax.text(.5, .5, "NOT_RUN: no stored dimensions", transform=ax.transAxes,
                    ha="center", va="center")
    axes[1].set(xticks=range(len(rows)), xticklabels=[row["case"].replace("_", "\n") for row in rows])
    axes[1].tick_params(axis="x", labelsize=8)
    fig.text(.5, .015, "Port counts are structural FE coordinates, not material parameters or evidence of cost savings.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .06, 1, 1))
    return save_figure(fig, output, "G0_DOFS_AND_PORTS", csv_path)


def response_figure(case_id, specification, measured, results_path, output):
    array_path = results_path.parent / f"{case_id}.npz"
    if case_id not in measured:
        return {"case": case_id, "status": "NOT_RUN_CASE_HAS_NO_RESULT", "source": relative(array_path)}
    if not array_path.exists():
        return {"case": case_id, "status": "NOT_RUN_RESPONSE_ARRAY_ABSENT", "source": relative(array_path)}
    with np.load(array_path, allow_pickle=False) as arrays:
        fixed = np.asarray(arrays["fixed_data"])
        fitted = np.asarray(arrays["reference_data"])
    if fixed.ndim != 2 or fixed.shape != fitted.shape or not fixed.size:
        raise ValueError(f"{case_id}: stored responses must have the same nonempty TX-by-RX shape")
    if not (np.all(np.isfinite(fixed)) and np.all(np.isfinite(fitted))):
        raise ValueError(f"{case_id}: response array contains nonfinite values")
    if not (np.iscomplexobj(fixed) and np.iscomplexobj(fitted)):
        raise ValueError(f"{case_id}: complex response arrays are required")
    csv_path = output / f"G0_RESPONSE_{case_id.upper()}.csv"
    data = []
    for name, array in (("fixed", fixed), ("interface_fitted", fitted)):
        for tx in range(array.shape[0]):
            for rx in range(array.shape[1]):
                value = array[tx, rx]
                data.append(dict(case=case_id, response=name, tx_index=tx, rx_index=rx,
                                 real=float(value.real), imag=float(value.imag),
                                 amplitude=float(abs(value)), phase_rad=float(np.angle(value))))
    write_csv(csv_path, data, ["case", "response", "tx_index", "rx_index", "real", "imag", "amplitude", "phase_rad"])
    fig, axes = plt.subplots(fixed.shape[0], 2, figsize=(11, 2.8 * fixed.shape[0] + 1),
                             squeeze=False, sharex=True)
    receivers = np.arange(fixed.shape[1])
    for tx in range(fixed.shape[0]):
        for col, transform, label in ((0, np.abs, "Complex response magnitude |d|"),
                                      (1, np.angle, "Principal phase arg(d), radians")):
            ax = axes[tx, col]
            ax.plot(receivers, transform(fixed[tx]), "o-", markersize=3, linewidth=1.2,
                    color=COLORS["fixed"], label="Fixed mesh")
            ax.plot(receivers, transform(fitted[tx]), "x--", markersize=3, linewidth=1.1,
                    color=COLORS["fitted"], label="Interface-fitted reference")
            ax.set(title=f"TX {tx}: {label}", ylabel="Magnitude" if col == 0 else "Phase / rad")
            ax.grid(alpha=.2)
            ax.legend(fontsize=8)
            if tx == fixed.shape[0] - 1:
                ax.set_xlabel("Stored receiver index (no interpolation)")
    fig.suptitle(f"{case_id}: actual stored complex responses\n"
                 f"radius={specification['radius']:g}, contrast={specification['contrast']:g}, "
                 f"target loss={specification['target_loss']:g}", fontsize=12)
    fig.text(.5, .015, "Single-frequency geometry comparison. Low loss / high contrast is a stress attempt; resonance is not established.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .055, 1, .91))
    record = save_figure(fig, output, f"G0_RESPONSE_{case_id.upper()}", csv_path)
    return dict(record, case=case_id, source=relative(array_path), shape=list(fixed.shape),
                phase_convention="principal np.angle, radians; raw real/imag retained in CSV")


def cost_snapshot(support):
    """FE/RHS: reservation records; CPU: one completion per wrapper job."""
    ledger = ROOT / "PHYSICS_LEDGER.jsonl"
    if not ledger.exists():
        summary = {"status": "NOT_RUN_LEDGER_ABSENT", "source": relative(ledger)}
        dump_json(support / "COST_SNAPSHOT.json", summary)
        return summary
    reservations, jobs, starts, seen, duplicate_lines = [], {}, set(), set(), []
    last_time = None
    for line_number, line in enumerate(ledger.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        record = json.loads(line)
        identity = json.dumps(record, sort_keys=True, allow_nan=False)
        if identity in seen:
            duplicate_lines.append(line_number)
            continue
        seen.add(identity)
        last_time = record.get("time_unix", last_time)
        if record.get("kind") == "reservation":
            for key in ("rhs", "factorizations"):
                value = record.get(key, 0)
                if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                    raise ValueError(f"Ledger line {line_number}: invalid {key} reservation")
            reservations.append(dict(record, ledger_line=line_number))
        if record.get("label") in ("job_start", "job_completion"):
            if not record.get("job") or not record.get("phase"):
                raise ValueError(f"Ledger line {line_number}: wrapper job/phase missing")
            key = (record["phase"], record["job"])
            if record["label"] == "job_start":
                starts.add(key)
            else:
                if key in jobs:
                    raise ValueError(f"Conflicting repeated wrapper completion: {key}")
                nonnegative(record.get("process_cpu_seconds"), "wrapper process CPU")
                nonnegative(record.get("wall_seconds"), "wrapper wall time")
                jobs[key] = dict(record, ledger_line=line_number)
    phases = sorted({row.get("phase", "UNSPECIFIED") for row in reservations} | {key[0] for key in jobs})
    by_phase = {}
    for phase in phases:
        physical = [row for row in reservations if row.get("phase", "UNSPECIFIED") == phase]
        completed = [row for key, row in jobs.items() if key[0] == phase]
        by_phase[phase] = dict(new_full_fe=sum(row.get("factorizations", 0) for row in physical),
                               new_complex_rhs=sum(row.get("rhs", 0) for row in physical),
                               measured_wrapper_process_cpu_seconds=sum(row["process_cpu_seconds"] for row in completed),
                               wrapper_wall_seconds=sum(row["wall_seconds"] for row in completed),
                               completed_wrappers=len(completed))
    baseline = read_json(ROOT / "COST_SUMMARY.json") if (ROOT / "COST_SUMMARY.json").exists() else {}
    summary = dict(status="LEDGER_SNAPSHOT", source=relative(ledger), through_time_unix=last_time,
                   new_full_fe=sum(value["new_full_fe"] for value in by_phase.values()),
                   new_complex_rhs=sum(value["new_complex_rhs"] for value in by_phase.values()),
                   measured_wrapper_process_cpu_seconds=sum(value["measured_wrapper_process_cpu_seconds"] for value in by_phase.values()),
                   wrapper_wall_seconds=sum(value["wrapper_wall_seconds"] for value in by_phase.values()),
                   peak_measured_job_RSS_bytes=max([row.get("children_peak_RSS", 0) for row in jobs.values()] + [0]),
                   by_phase=by_phase, completed_wrappers=len(jobs), reservation_records=len(reservations),
                   identical_duplicate_ledger_lines_ignored=duplicate_lines,
                   unfinished_wrappers=[dict(phase=key[0], job=key[1]) for key in sorted(starts - set(jobs))],
                   authority="FE/RHS from reservation once; CPU/wall from job_completion once; kernel and receipt copies not added",
                   snapshot_note="In-child snapshot excludes unfinished wrappers; the final reporting-wrapper receipt supersedes it.",
                   accounting_gap=baseline.get("accounting_gap", "NOT_MEASURED: unwrapped work is not zero"),
                   historical_bank=baseline.get("historical_bank", {"status": "NOT_READ"}),
                   standalone_deployment_cost=baseline.get("standalone_deployment_cost", "NOT_MEASURED"),
                   physical_speed_claim=False)
    write_csv(support / "WRAPPER_COSTS.csv", list(jobs.values()),
              ["ledger_line", "phase", "job", "exit_code", "process_cpu_seconds", "wall_seconds", "children_peak_RSS"])
    write_csv(support / "PHYSICS_RESERVATIONS.csv", reservations,
              ["ledger_line", "phase", "label", "rhs", "factorizations", "case", "backend", "freq"])
    dump_json(support / "COST_SNAPSHOT.json", summary)
    return summary


def stored_tests():
    directory = ROOT / "results/jobs/implementation-tests-32"
    receipt_path, log_path = directory / "RECEIPT.json", directory / "OUTPUT.log"
    if not receipt_path.exists() or not log_path.exists():
        return {"status": "NOT_RUN_RECEIPT_OR_LOG_ABSENT", "receipt": relative(receipt_path)}
    receipt, log = read_json(receipt_path), log_path.read_text(encoding="utf-8")
    match = re.search(r"^Ran (\d+) tests?\b", log, re.MULTILINE)
    count = int(match[1]) if match else None
    passed = receipt.get("exit_code") == 0 and count == 32 and bool(re.search(r"^OK\s*$", log, re.MULTILINE))
    return dict(status="PASS_STORED_32_TESTS" if passed else "CHECK_STORED_TEST_EVIDENCE",
                count=count, exit_code=receipt.get("exit_code"), receipt=relative(receipt_path),
                log=relative(log_path), new_execution=False,
                scope="Existing algebra/accounting/implementation tests; not new Maxwell or continuum qualification")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "results/g0_qualification/RESULTS.json")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "figures/g0_qualification")
    parser.add_argument("--summary-dir", type=Path, default=SUPPORT)
    parser.add_argument("--response-cases", nargs="+", default=list(DEFAULT_RESPONSE_CASES))
    args = parser.parse_args()
    if not 2 <= len(args.response_cases) <= 3 or len(set(args.response_cases)) != len(args.response_cases):
        parser.error("Choose two or three distinct predeclared response cases")
    specification = read_json(ROOT / "G0_CASES.json")
    cases = specification["qualification"]
    by_id = {case["id"]: case for case in cases}
    if set(args.response_cases) - set(by_id):
        parser.error("Response cases must occur in G0_CASES.json")
    limit = nonnegative(read_json(ROOT / "GATES.json")["G0"]["independent_whitened_max"], "frozen G0 limit")
    if limit <= 0:
        raise ValueError("The frozen G0 limit must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.summary_dir.mkdir(parents=True, exist_ok=True)
    source_status, rows, measured = result_table(args.results, cases, limit)
    figures = [error_figure(rows, args.output_dir, limit), dimension_figure(rows, args.output_dir)]
    responses = [response_figure(case_id, by_id[case_id], measured, args.results, args.output_dir)
                 for case_id in args.response_cases]
    cost = cost_snapshot(args.summary_dir)
    values = [row["fixed_vs_interface_fitted_white"] for row in rows
              if row["fixed_vs_interface_fitted_white"] is not None]
    port_plan_path = args.results.parent / "PORT_PLAN.json"
    plan = read_json(port_plan_path) if port_plan_path.exists() else {"status": "NOT_RUN_PORT_PLAN_ABSENT"}
    mesh = [ROOT / "figures" / f"PUBLIC_MESH_AND_PORTS.{suffix}" for suffix in ("png", "pdf")]
    summary = dict(status=source_status, source=relative(args.results), expected_cases=len(cases),
                   measured_cases=len(measured), metric_cases=len(values),
                   missing_cases=[row["case"] for row in rows if row["status"] in ("NOT_RUN", "MISSING_METRIC")],
                   cases_over_frozen_error_limit=[row["case"] for row in rows if row["status"] == "OVER_FROZEN_LIMIT"],
                   frozen_independent_whitened_max=limit,
                   min_stored_whitened_error=min(values) if values else None,
                   max_stored_whitened_error=max(values) if values else None,
                   figures=figures, responses=responses, port_plan=plan,
                   reused_fixed_mesh_figures=[relative(path) for path in mesh if path.exists()],
                   fixed_mesh_array_source=relative(ROOT / "results/preflight/FIXED_MESH.npz"),
                   implementation_tests=stored_tests(), cost_snapshot=cost,
                   scientific_scope=SCOPE_NOTE,
                   high_sensitivity_status=specification.get("high_sensitivity_status", "NOT_ESTABLISHED"),
                   gates_modified=False, physical_solves_executed=False)
    dump_json(args.summary_dir / "FIGURE_SUMMARY.json", summary)
    print(json.dumps(dict(status=source_status, measured_cases=len(measured),
                          figures_created=sum(len(item["figures"]) for item in figures)
                          + sum(len(item.get("figures", [])) for item in responses),
                          summary=relative(args.summary_dir / "FIGURE_SUMMARY.json"),
                          physical_solves_executed=False)))


if __name__ == "__main__":
    main()
