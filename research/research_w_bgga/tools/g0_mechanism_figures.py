"""Read stored G0 diagnostics and draw factual mechanism figures, without solves.

No physical backend, truth data, posterior reconstruction, or gate decision is
used here.  Each figure requires its own complete input group.  CSVs retain
missing/undefined rows and every plotted number in the stored units.  Run this
program only through the existing serial job wrapper for CPU accounting.
"""
from pathlib import Path
import argparse
import csv
import json
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
TRAINING = ("weak_center", "left_small")
HELD_OUT = ("strong_center", "right_small", "left_large", "right_large", "low_loss_high_contrast")
QUALIFICATION_CASES = ("weak_center", "strong_center", "left_small", "right_small", "left_large", "right_large",
                       "low_loss_high_contrast", "asym_low", "asym_mid", "asym_high", "empty_wall", "two_region_wall")
METHODS = ("full", "woodbury", "schur_from_paid_ports", "krylov", "goal_rom")
METHOD_LABELS = ("Full", "Woodbury", "Schur from paid ports", "Krylov", "Goal ROM")
METHOD_COLORS = ("#37474f", "#247a8a", "#a85c42", "#6b68a6", "#b58c27")
FAMILIES = ("C1", "C2", "C3", "C4")
FIXTURES = tuple(f"{family}_{index}" for family in FAMILIES for index in range(3))
SCOPES = ("pose_only", "pose_intrinsic", "pose_intrinsic_instrument")
FAMILY_COLORS = dict(zip(FAMILIES, ("#247a8a", "#b58c27", "#6b68a6", "#a85c42")))
SETUP_STEPS = (
    "one public fixed P3 FEM assembly", "public target assembler",
    "shared full-2x2 eps=3.75 ell=.125 background assembly",
    "one background sparse LU and two primal RHS",
    "one complete paid public Green-port setup",
    "sixteen paid background receiver adjoint RHS", "fixed joint two-pass goal QR")
MATCHED_NOTE = "12 fixtures; 9 physical worlds. C2/C3 share the same three worlds (matched observation contracts)."


def rel(path):
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def scalar(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        return None
    return value


def vector(value):
    if not isinstance(value, list) or any(scalar(item) is None for item in value):
        return None
    return value


def nested(value, *keys):
    for key in keys:
        if not isinstance(value, dict):
            return None
        value = value.get(key)
    return value


def read_optional(path, kind, issues):
    if not path.exists():
        issues.append(dict(source=rel(path), status="NOT_RUN_FILE_ABSENT"))
        return None
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        issues.append(dict(source=rel(path), status="NOT_RUN_JSON_INCOMPLETE_OR_INVALID"))
        return None
    if not isinstance(value, kind):
        raise ValueError(f"Unexpected JSON structure: {path}")
    return value


def index_rows(rows, keys):
    result = {}
    for row in rows or []:
        if not isinstance(row, dict):
            raise ValueError("Result rows must be JSON objects")
        identity = tuple(row.get(key) for key in keys)
        if identity in result:
            raise ValueError(f"Duplicate stored result: {identity}")
        result[identity] = row
    return result


def csv_file(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, ensure_ascii=False) if isinstance(value, (list, dict))
                             else value for key, value in row.items()})
    return rel(path)


def range_of(values):
    values = [scalar(value) for value in values if scalar(value) is not None]
    return dict(count=len(values), minimum=min(values) if values else None,
                maximum=max(values) if values else None)


def save_plot(fig, output, stem, data_csv, **metadata):
    files = []
    for extension in ("png", "pdf"):
        path = output / f"{stem}.{extension}"
        fig.savefig(path, bbox_inches="tight", **({"dpi": 180} if extension == "png" else {}))
        files.append(rel(path))
    plt.close(fig)
    return dict(name=stem, status="CREATED", figures=files, data_csv=data_csv, **metadata)


def skip_plot(output, stem, reason, missing):
    old = [rel(output / f"{stem}.{extension}") for extension in ("png", "pdf")
           if (output / f"{stem}.{extension}").exists()]
    return dict(name=stem, status="NOT_RUN_INCOMPLETE_INPUT", reason=reason, missing=missing,
                preexisting_files_not_refreshed=old)


def partition(case):
    return "TRAINING" if case in TRAINING else "HELD_OUT"


def operators(rows, plan, setup, core, output):
    indexed = index_rows(rows, ("case", "method"))
    table, missing = [], []
    for case in TRAINING + HELD_OUT:
        for method in METHODS:
            source = indexed.get((case, method), {})
            cpu = scalar(nested(source, "solve_timing", "measured_kernel_cpu_seconds"))
            error = scalar(source.get("full_response_difference_whitened_l2"))
            status = "MEASURED" if cpu is not None and cpu >= 0 and error is not None else "NOT_RUN_OR_METRIC_MISSING"
            if status != "MEASURED":
                missing.append(f"{case}/{method}")
            table.append(dict(case=case, method=method, partition=partition(case), status=status,
                              whitened_response_error=error, query_solve_cpu_seconds=cpu,
                              query_solve_wall_seconds=scalar(nested(source, "solve_timing", "measured_kernel_wall_seconds")),
                              verification_cpu_seconds=scalar(source.get("verification_cpu_seconds")),
                              artifact_cpu_seconds=scalar(source.get("artifact_cpu_seconds")),
                              strict_algebra_status=source.get("strict_algebra_status", "NOT_RUN"),
                              source="results/g0_operators/RESULTS.json"))
    data_csv = csv_file(output / "OPERATOR_QUERY_ERROR_CPU.csv", table,
                       ["case", "method", "partition", "status", "whitened_response_error",
                        "query_solve_cpu_seconds", "query_solve_wall_seconds", "verification_cpu_seconds",
                        "artifact_cpu_seconds", "strict_algebra_status", "source"])
    valid_plan = (plan is not None and tuple(plan.get("training_ids", [])) == TRAINING
                  and tuple(plan.get("held_out_ids", [])) == HELD_OUT)
    if not valid_plan:
        missing.append("frozen training/held-out PLAN")
    stem = "OPERATOR_ERROR_VS_QUERY_CPU"
    if missing:
        query_plot = skip_plot(output, stem, "All seven cases and five methods plus partition plan are required", missing)
    else:
        fig, axes = plt.subplots(1, 7, figsize=(23, 5.1), sharey=True)
        for ax, case in zip(axes, TRAINING + HELD_OUT):
            ax.set_facecolor("#f0f6f1" if case in TRAINING else "#faf4ec")
            for method, label, color, marker in zip(METHODS, METHOD_LABELS, METHOD_COLORS, ("o", "s", "D", "^", "X")):
                row = next(item for item in table if item["case"] == case and item["method"] == method)
                ax.scatter(row["query_solve_cpu_seconds"], row["whitened_response_error"],
                           label=label, color=color, marker=marker, s=46)
            ax.set(xscale="log", title=f"{partition(case)}\n{case.replace('_', chr(10))}",
                   xlabel="Query solve CPU / s")
            ax.set_yscale("symlog", linthresh=1e-11)
            ax.tick_params(labelsize=8)
            ax.grid(alpha=.2)
        axes[0].set_ylabel("Whitened error versus stored full solve (symlog; actual zero retained)")
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=5, bbox_to_anchor=(.5, 1.01), fontsize=10)
        fig.text(.5, .015, "Measured solve kernels only: shared setup, target-update preparation, verification and I/O are separate. Query CPU is not cold end-to-end CPU.",
                 ha="center", fontsize=10)
        fig.tight_layout(rect=(0, .075, 1, .89))
        query_plot = save_plot(fig, output, stem, data_csv, y_display_transform="symlog linthresh=1e-11; data unchanged")
    setup_index = {row.get("step"): row for row in setup or []}
    setup_table = [dict(step=step, status="MEASURED" if step in setup_index else "NOT_RUN",
                        measured_kernel_cpu_seconds=scalar(setup_index.get(step, {}).get("measured_kernel_cpu_seconds")),
                        measured_kernel_wall_seconds=scalar(setup_index.get(step, {}).get("measured_kernel_wall_seconds")))
                   for step in SETUP_STEPS]
    setup_csv = csv_file(output / "OPERATOR_SHARED_SETUP.csv", setup_table,
                         ["step", "status", "measured_kernel_cpu_seconds", "measured_kernel_wall_seconds"])
    setup_missing = [row["step"] for row in setup_table if row["measured_kernel_cpu_seconds"] is None]
    if setup_missing:
        setup_plot = skip_plot(output, "OPERATOR_SHARED_SETUP", "All seven shared setup kernels are required", setup_missing)
    else:
        fig, ax = plt.subplots(figsize=(12, 5.4))
        ax.barh(range(len(setup_table)), [row["measured_kernel_cpu_seconds"] for row in setup_table], color="#247a8a")
        ax.set(yticks=range(len(setup_table)), yticklabels=[row["step"] for row in setup_table],
               xlabel="Recorded setup kernel process CPU / s", title="Paid shared setup, shown separately from query solves")
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=.2)
        fig.text(.5, .015, "Training full solves remain their recorded baseline solves. Setup bars are not a standalone deployment cost or a cold-time total.",
                 ha="center", fontsize=8)
        fig.tight_layout(rect=(0, .055, 1, 1))
        setup_plot = save_plot(fig, output, "OPERATOR_SHARED_SETUP", setup_csv)
    equations = []
    for row in table:
        source = indexed.get((row["case"], row["method"]), {})
        equations.append(dict(row,
            full_equation_absolute=scalar(nested(source, "full_equation_backward", "absolute")),
            full_equation_backward=scalar(nested(source, "full_equation_backward", "backward")),
            full_equation_rhs_relative=scalar(nested(source, "full_equation_backward", "rhs_relative")),
            method_reported_backward=scalar(nested(source, "method_reported", "backward")),
            core_backward=scalar(nested(source, "method_reported", "core_backward", "backward")),
            reduced_equation_backward=scalar(nested(source, "reduced_equation_backward", "backward")),
            field_relative_error=scalar(nested(source, "field_difference_to_full", "relative")),
            response_relative_error=scalar(nested(source, "response_difference_to_full", "relative"))))
    eq_csv = csv_file(output / "OPERATOR_EQUATIONS.csv", equations,
        ["case", "method", "partition", "status", "strict_algebra_status", "full_equation_absolute",
         "full_equation_backward", "full_equation_rhs_relative", "method_reported_backward", "core_backward",
         "reduced_equation_backward", "field_relative_error", "response_relative_error", "source"])
    core_table = []
    for key in ("case", "description", "core", "dimension", "core_one_norm", "inverse_one_norm_estimate",
                "condition_one_norm_estimate_product", "interpretation", "singular_values", "resonance_claim",
                "jitter", "pseudoinverse", "measured_kernel_cpu_seconds", "measured_kernel_wall_seconds"):
        core_table.append(dict(quantity=key, value=core.get(key) if core else None,
                               status="RECORDED" if core and key in core else "NOT_RUN_OR_NOT_RECORDED"))
    core_table.append(dict(quantity="actual_equation_backward", value=nested(core, "actual_equation_backward", "backward"),
                           status="RECORDED" if nested(core, "actual_equation_backward", "backward") is not None else "NOT_RUN_OR_NOT_RECORDED"))
    core_csv = csv_file(output / "CORE_STABILITY.csv", core_table, ["quantity", "value", "status"])
    ranges = {method: dict(white_error=range_of(row["whitened_response_error"] for row in table if row["method"] == method),
                           query_cpu=range_of(row["query_solve_cpu_seconds"] for row in table if row["method"] == method))
              for method in METHODS}
    return [query_plot, setup_plot], [data_csv, setup_csv, eq_csv, core_csv], ranges


def fixture_label(case, indexed):
    parent = indexed.get((case,), {}).get("physical_parent")
    return f"{case}\nP{parent if parent is not None else '?'}{'*' if case[:2] in ('C2', 'C3') else ''}"


def stress_spectra(rows, output):
    indexed = index_rows(rows, ("case",))
    for index in range(3):
        first, second = indexed.get((f"C2_{index}",)), indexed.get((f"C3_{index}",))
        if first and second and first.get("physical_parent") != second.get("physical_parent"):
            raise ValueError("C2/C3 matched-world provenance is inconsistent")
    spectra, statuses, plots = [], [], []
    for scope in SCOPES:
        missing = []
        for case in FIXTURES:
            source = indexed.get((case,), {})
            data = nested(source, "diagnostics", scope, "data_only") or {}
            state = data.get("status", "NOT_RUN")
            ret_state = data.get("generalized_retention_status", "NOT_RUN")
            schur = vector(data.get("wall_data_only_schur_eigenvalues"))
            retention = vector(data.get("generalized_retention_eigenvalues"))
            if not schur or retention is None or (not retention and not ret_state.startswith("UNDEFINED")):
                missing.append(case)
            common = dict(case=case, family=source.get("family", case[:2]), physical_parent=source.get("physical_parent"),
                          matched_C2_C3=case[:2] in ("C2", "C3"), nuisance_scope=scope,
                          recorded_status=state, retention_status=ret_state,
                          retained_wall_dimension=data.get("retained_wall_dimension"), excluded_wall_dimension=data.get("excluded_wall_dimension"))
            statuses.append(dict(common, operator_error=source.get("empirical_operator_error"),
                wall_rank=nested(data, "wall_rank_audit", "rank"),
                wall_rank_lower=nested(data, "wall_rank_audit", "rank_lower"), wall_rank_upper=nested(data, "wall_rank_audit", "rank_upper"),
                nuisance_rank=nested(data, "nuisance_rank_audit", "rank"),
                nuisance_rank_lower=nested(data, "nuisance_rank_audit", "rank_lower"), nuisance_rank_upper=nested(data, "nuisance_rank_audit", "rank_upper")))
            for spectrum, values, status in (("wall_absolute", vector(data.get("wall_absolute_eigenvalues")), state),
                                             ("wall_data_only_schur", schur, state),
                                             ("generalized_retention", retention, ret_state)):
                if values:
                    spectra.extend(dict(common, spectrum=spectrum, eigen_index=index, eigenvalue=value,
                                        value_status=status) for index, value in enumerate(values))
                else:
                    spectra.append(dict(common, spectrum=spectrum, eigen_index=None, eigenvalue=None,
                                        value_status=status if values == [] else "NOT_RUN_OR_ARRAY_MISSING"))
        stem = f"STRESS_SCHUR_AND_RETENTION_{scope.upper()}"
        if missing:
            plots.append(skip_plot(output, stem, "Every registered fixture needs this scope's stored arrays", missing))
            continue
        fig, axes = plt.subplots(1, 2, figsize=(15, 5.8))
        for index, case in enumerate(FIXTURES):
            data = nested(indexed[(case,)], "diagnostics", scope, "data_only")
            for ax, key, status_key in ((axes[0], "wall_data_only_schur_eigenvalues", "status"),
                                        (axes[1], "generalized_retention_eigenvalues", "generalized_retention_status")):
                values = data[key]
                resolved = data[status_key] == "RESOLVED_AT_SUPPLIED_ERROR_RESOLUTION"
                for mode, value in enumerate(values):
                    ax.scatter(index + (mode - (len(values) - 1) / 2) * .16, value,
                               marker=("o", "s", "D")[mode % 3], s=43,
                               edgecolor=FAMILY_COLORS[case[:2]],
                               facecolor=FAMILY_COLORS[case[:2]] if resolved else "none", linewidth=1.2)
                if not values:
                    ax.text(index, .5, "UNDEFINED", rotation=90, ha="center", fontsize=7,
                            transform=ax.get_xaxis_transform())
                excluded = data.get("excluded_wall_dimension", 0)
                if ax is axes[1] and excluded:
                    ax.text(index, .02, f"{excluded} excluded", rotation=90, ha="center", fontsize=7,
                            transform=ax.get_xaxis_transform())
            for ax in axes:
                ax.set(xticks=range(12), xticklabels=[fixture_label(item, indexed) for item in FIXTURES], xlim=(-.6, 11.6))
                ax.tick_params(axis="x", labelsize=7)
                ax.grid(axis="y", alpha=.2)
        axes[0].set_yscale("symlog", linthresh=1e-12)
        axes[0].set(title="Absolute data-only wall Schur spectrum", ylabel="Recorded eigenvalue (signed symlog; no clipping)")
        axes[1].set(title="Retention on retained A response support", ylabel="Recorded generalized retention eigenvalue", ylim=(-.05, 1.08))
        legend = [Line2D([], [], marker="o", ls="", color="#455a64", label="Recorded rank status resolved"),
                  Line2D([], [], marker="o", ls="", color="#455a64", markerfacecolor="none", label="Recorded rank status UNRESOLVED")]
        fig.legend(handles=legend, loc="upper center", ncol=2, bbox_to_anchor=(.5, .96), fontsize=9)
        fig.suptitle(f"{scope}: stored spectra conditional on recorded error-threshold ranks", fontsize=12)
        fig.text(.5, .035, MATCHED_NOTE, ha="center", fontsize=9)
        fig.text(.5, .008, "Excluded A directions have undefined retention; no Schur/projector error envelope is supplied.", ha="center", fontsize=9)
        fig.tight_layout(rect=(0, .085, 1, .86))
        plots.append(save_plot(fig, output, stem, [rel(output / "STRESS_SPECTRA.csv"), rel(output / "STRESS_STATUS.csv")],
                               absolute_spectrum_display="symlog linthresh=1e-12; signed raw values retained"))
    spectra_csv = csv_file(output / "STRESS_SPECTRA.csv", spectra,
        ["case", "family", "physical_parent", "matched_C2_C3", "nuisance_scope", "spectrum", "eigen_index", "eigenvalue",
         "value_status", "recorded_status", "retention_status", "retained_wall_dimension", "excluded_wall_dimension"])
    status_csv = csv_file(output / "STRESS_STATUS.csv", statuses,
        ["case", "family", "physical_parent", "matched_C2_C3", "nuisance_scope", "recorded_status", "retention_status",
         "operator_error", "wall_rank", "wall_rank_lower", "wall_rank_upper", "nuisance_rank", "nuisance_rank_lower",
         "nuisance_rank_upper", "retained_wall_dimension", "excluded_wall_dimension"])
    ranges = {scope: dict(schur=range_of(row["eigenvalue"] for row in spectra if row["nuisance_scope"] == scope and row["spectrum"] == "wall_data_only_schur"),
                          retention=range_of(row["eigenvalue"] for row in spectra if row["nuisance_scope"] == scope and row["spectrum"] == "generalized_retention"),
                          unresolved=[row["case"] for row in statuses if row["nuisance_scope"] == scope and row["recorded_status"] == "UNRESOLVED"],
                          undefined_retention=[row["case"] for row in statuses if row["nuisance_scope"] == scope and row["retention_status"].startswith("UNDEFINED")])
              for scope in SCOPES}
    return indexed, plots, [spectra_csv, status_csv], ranges


def finite_diagnostics(indexed, output):
    variance, pairs, components, checks = [], [], [], []
    for case in FIXTURES:
        source = indexed.get((case,), {})
        common = dict(case=case, family=source.get("family", case[:2]), physical_parent=source.get("physical_parent"),
                      matched_C2_C3=case[:2] in ("C2", "C3"))
        value = scalar(nested(source, "finite_posterior", "wall_variance_fraction_of_prior"))
        variance.append(dict(common, wall_variance_fraction_of_prior=value,
                             status="RECORDED_FINITE15_SCALAR" if value is not None else "NOT_RUN_OR_MISSING",
                             prior_id=nested(source, "finite_posterior", "prior_id"), labels=nested(source, "finite_posterior", "labels")))
        finite = source.get("finite_pair") or {}
        pair = dict(common, measurement=source.get("measurement"),
                    task_measurement="TX0 RX8" if case[:2] == "C3" else "2TX all16RX",
                    current_endpoint_contract="mean_joint endpoints: wall and corrected target pose",
                    task_endpoint_contract="same fixed declared target j0 in both wall endpoints",
                    recorded_task_contract_id=finite.get("task_contract_id"), ratio_status=finite.get("ratio_status", "NOT_RUN"))
        for key in ("current_norm", "task_norm_value", "current_error", "task_error", "current_norm_lower", "current_norm_upper",
                    "task_norm_lower", "task_norm_upper", "task_to_current_ratio", "ratio_lower", "ratio_upper"):
            pair[key] = scalar(finite.get(key))
        pair["status"] = ("RECORDED" if all(pair[key] is not None for key in
            ("current_norm", "task_norm_value", "current_norm_lower", "current_norm_upper", "task_norm_lower", "task_norm_upper"))
            and finite.get("task_norm") == "euclidean" and isinstance(finite.get("task_contract_id"), str)
            and "same fixed declared target j0" in finite["task_contract_id"] else "NOT_RUN_OR_CONTRACT_MISSING")
        pairs.append(pair)
        for name in ("current_difference", "task_difference"):
            values = vector(finite.get(name))
            if values:
                components.extend(dict(common, vector=name, coordinate=index, value=value) for index, value in enumerate(values))
            else:
                components.append(dict(common, vector=name, coordinate=None, value=None))
        duality = source.get("duality") or {}
        dual = scalar(duality.get("wall_JVP_adjoint_relative"))
        checks.append(dict(common, measurement=source.get("measurement"),
            white_mesh_difference=scalar(source.get("white_mesh_difference")),
            empirical_scaled_operator_error=scalar(source.get("empirical_operator_error")),
            wall_JVP_adjoint_relative=dual, duality_status="RECORDED" if dual is not None else duality.get("status", "NOT_RUN_OR_MISSING"),
            adjoint_equation_absolute=scalar(nested(duality, "adjoint_backward", "absolute")),
            adjoint_equation_backward=scalar(nested(duality, "adjoint_backward", "backward")),
            adjoint_equation_rhs_relative=scalar(nested(duality, "adjoint_backward", "rhs_relative")),
            source="results/g0_stress/RESULTS.json"))
    variance_csv = csv_file(output / "FINITE15_WALL_VARIANCE.csv", variance,
        ["case", "family", "physical_parent", "matched_C2_C3", "status", "wall_variance_fraction_of_prior", "prior_id", "labels"])
    pair_fields = ["case", "family", "physical_parent", "matched_C2_C3", "status", "measurement", "task_measurement",
        "current_endpoint_contract", "task_endpoint_contract", "recorded_task_contract_id", "current_norm", "task_norm_value",
        "current_error", "task_error", "current_norm_lower", "current_norm_upper", "task_norm_lower", "task_norm_upper",
        "task_to_current_ratio", "ratio_lower", "ratio_upper", "ratio_status"]
    pairs_csv = csv_file(output / "FINITE_MEAN_JOINT_AND_FIXED_J0.csv", pairs, pair_fields)
    component_csv = csv_file(output / "FINITE_PAIR_COMPONENTS.csv", components,
        ["case", "family", "physical_parent", "matched_C2_C3", "vector", "coordinate", "value"])
    checks_csv = csv_file(output / "STRESS_EQUATION_DUALITY_MESH.csv", checks,
        ["case", "family", "physical_parent", "matched_C2_C3", "measurement", "white_mesh_difference",
         "empirical_scaled_operator_error", "wall_JVP_adjoint_relative", "duality_status", "adjoint_equation_absolute",
         "adjoint_equation_backward", "adjoint_equation_rhs_relative", "source"])
    plots = []
    missing = [row["case"] for row in variance if row["wall_variance_fraction_of_prior"] is None]
    if missing:
        plots.append(skip_plot(output, "FINITE15_WALL_VARIANCE", "All 12 recorded finite15 variance fractions are required", missing))
    else:
        fig, ax = plt.subplots(figsize=(13.8, 5.2))
        values = [row["wall_variance_fraction_of_prior"] for row in variance]
        ax.scatter(range(12), values, c=[FAMILY_COLORS[row["family"]] for row in variance], s=60)
        for index, value in enumerate(values):
            ax.annotate(f"{value:.3g}", (index, value), xytext=(0, 8), textcoords="offset points", ha="center", fontsize=8)
        ax.set_yscale("log" if min(values) > 0 else "symlog", **({} if min(values) > 0 else {"linthresh": 1e-10}))
        ax.axhline(1, color="#708090", ls="--", lw=.8, label="Prior variance fraction = 1")
        ax.set(xticks=range(12), xticklabels=[fixture_label(case, indexed) for case in FIXTURES],
               ylabel="Stored standardized wall variance / finite prior variance",
               title="Actual 15-atom known-support prior: stored wall variance fraction")
        ax.grid(axis="y", alpha=.2)
        ax.legend(fontsize=8)
        fig.text(.5, .015, MATCHED_NOTE, ha="center", fontsize=9)
        fig.tight_layout(rect=(0, .06, 1, 1))
        plots.append(save_plot(fig, output, "FINITE15_WALL_VARIANCE", variance_csv))
    missing = [row["case"] for row in pairs if row["status"] != "RECORDED"]
    if missing:
        plots.append(skip_plot(output, "FINITE_MEAN_JOINT_VS_FIXED_J0", "All finite pair norms and fixed-j0 contracts are required", missing))
    else:
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.7))
        for ax, families, title in ((axes[0], ("C1", "C2", "C4"), "C1/C2/C4: complete 2TX/16RX fixed-j0 task"),
                                     (axes[1], ("C3",), "C3: same single-port fixed-j0 task")):
            for row in pairs:
                if row["family"] not in families:
                    continue
                x, y = row["current_norm"], row["task_norm_value"]
                errors = [x-row["current_norm_lower"], row["current_norm_upper"]-x,
                          y-row["task_norm_lower"], row["task_norm_upper"]-y]
                if min(errors) < 0:
                    raise ValueError(f"Stored finite-pair interval does not contain norm: {row['case']}")
                color = FAMILY_COLORS[row["family"]]
                resolved = row["ratio_status"] == "DEFINED_CONDITIONAL_ON_CALLER_ERRORS"
                ax.errorbar(x, y, xerr=np.array(errors[:2])[:, None], yerr=np.array(errors[2:])[:, None],
                            fmt="o", color=color, mfc=color if resolved else "none", capsize=3)
                ax.annotate(row["case"], (x, y), xytext=(5, 5), textcoords="offset points", fontsize=8)
            local = [row for row in pairs if row["family"] in families]
            if all(row["current_norm_lower"] > 0 and row["task_norm_lower"] > 0 for row in local):
                ax.set(xscale="log", yscale="log")
            ax.set(title=title, xlabel="Current-data distance of mean_joint endpoints", ylabel="Task distance at the same fixed target j0")
            ax.grid(alpha=.2)
        fig.text(.5, .04, "Current endpoints change wall + pose. Task endpoints use one fixed target j0; no joint-target-change task is labeled H(W).", ha="center", fontsize=9)
        fig.text(.5, .012, "C1/C2 current data: 2TX/16RX; C3/C4 current data: TX0 RX8. Error bars are stored empirical pair-error intervals.", ha="center", fontsize=9)
        fig.tight_layout(rect=(0, .10, 1, 1))
        plots.append(save_plot(fig, output, "FINITE_MEAN_JOINT_VS_FIXED_J0", [pairs_csv, component_csv]))
    ranges = dict(finite15_variance_fraction=range_of(row["wall_variance_fraction_of_prior"] for row in variance),
                  current_mean_joint_distance=range_of(row["current_norm"] for row in pairs),
                  fixed_j0_task_distance=range_of(row["task_norm_value"] for row in pairs),
                  wall_jvp_adjoint_relative=range_of(row["wall_JVP_adjoint_relative"] for row in checks),
                  mesh_white_difference=range_of(row["white_mesh_difference"] for row in checks),
                  missing_duality_scalars=[row["case"] for row in checks if row["wall_JVP_adjoint_relative"] is None])
    return plots, [variance_csv, pairs_csv, component_csv, checks_csv], ranges


def c4_future(indexed, output):
    probes, maxima, missing = [], [], []
    for case in ("C4_0", "C4_1", "C4_2"):
        source = indexed.get((case,), {})
        future = source.get("future25") or {}
        distances = vector(future.get("distances"))
        endpoint_errors = vector(future.get("endpoint_numerical_differences"))
        recorded_max, independent = scalar(future.get("max")), scalar(future.get("independent_argmax_distance"))
        complete = (source.get("future25_status") == "EXECUTED_ALL25_AND_INDEPENDENT_ARGMAX"
                    and distances is not None and len(distances) == 25 and endpoint_errors is not None
                    and len(endpoint_errors) == 2 and recorded_max is not None and independent is not None)
        if not complete:
            missing.append(case)
        if distances:
            probes.extend(dict(case=case, probe_index=index, recorded_distance=value,
                               status=source.get("future25_status", "NOT_RUN")) for index, value in enumerate(distances))
        else:
            probes.append(dict(case=case, probe_index=None, recorded_distance=None, status="NOT_RUN_OR_ARRAY_MISSING"))
        error_sum = sum(endpoint_errors) if endpoint_errors is not None and len(endpoint_errors) == 2 else None
        maxima.append(dict(case=case, status=source.get("future25_status", "NOT_RUN"), recorded_max25=recorded_max,
            independent_argmax_distance=independent, argmax=future.get("argmax"), argmax_position=future.get("position"),
            endpoint_error_first=endpoint_errors[0] if endpoint_errors is not None and len(endpoint_errors) == 2 else None,
            endpoint_error_second=endpoint_errors[1] if endpoint_errors is not None and len(endpoint_errors) == 2 else None,
            endpoint_error_sum=error_sum, recorded_empirical_margin_lower=scalar(future.get("empirical_margin_lower")),
            derived_argmax_comparison_upper=recorded_max+error_sum if recorded_max is not None and error_sum is not None else None,
            array_max_matches_record=bool(max(distances) == recorded_max) if distances and recorded_max is not None else None,
            recorded_norm=future.get("norm"), continuum_certificate=future.get("continuum_certificate")))
    probe_csv = csv_file(output / "C4_FUTURE25_DISTANCES.csv", probes,
        ["case", "probe_index", "recorded_distance", "status"])
    max_csv = csv_file(output / "C4_FUTURE25_ARGMAX.csv", maxima,
        ["case", "status", "recorded_max25", "independent_argmax_distance", "argmax", "argmax_position", "endpoint_error_first",
         "endpoint_error_second", "endpoint_error_sum", "recorded_empirical_margin_lower", "derived_argmax_comparison_upper",
         "array_max_matches_record", "recorded_norm", "continuum_certificate"])
    if missing:
        plot = skip_plot(output, "C4_FUTURE25_AND_INDEPENDENT_ARGMAX", "All three complete 25-probe arrays and matched independent argmax checks are required", missing)
    else:
        fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.5))
        for index, row in enumerate(maxima):
            color = ("#247a8a", "#6b68a6", "#a85c42")[index]
            local = [item for item in probes if item["case"] == row["case"]]
            axes[0].plot([item["probe_index"] for item in local], [item["recorded_distance"] for item in local],
                         "o-", color=color, label=row["case"], markersize=3, linewidth=.8)
            axes[1].errorbar(index-.08, row["recorded_max25"], yerr=row["endpoint_error_sum"],
                             fmt="o", color=color, capsize=4)
            axes[1].scatter(index+.08, row["independent_argmax_distance"], marker="x", color=color, s=55)
        axes[0].set(title="All 25 declared fixed probes", xlabel="Stored fixed-grid probe index", ylabel="Whitened complete-response distance")
        axes[0].legend(fontsize=9)
        axes[1].set(title="Recorded max25 and independent check at its argmax", xticks=range(3),
                    xticklabels=[row["case"] for row in maxima], ylabel="Whitened complete-response distance")
        axes[1].legend(handles=[Line2D([], [], marker="o", ls="", color="#455a64", label="Fixed-grid max25"),
                               Line2D([], [], marker="x", ls="", color="#455a64", label="Independent distance at stored argmax")], fontsize=8)
        for ax in axes:
            ax.grid(alpha=.2)
        fig.text(.5, .04, "Bars use the sum of two endpoint numerical differences at the stored argmax; they are not an error envelope for all 25 probes.", ha="center", fontsize=8)
        fig.text(.5, .012, "Every task compares the same declared target probe in both wall endpoints. Recorded empirical margins have no continuum certificate.", ha="center", fontsize=8)
        fig.tight_layout(rect=(0, .10, 1, 1))
        plot = save_plot(fig, output, "C4_FUTURE25_AND_INDEPENDENT_ARGMAX", [probe_csv, max_csv])
    return plot, [probe_csv, max_csv], dict(max25=range_of(row["recorded_max25"] for row in maxima),
        independent_argmax=range_of(row["independent_argmax_distance"] for row in maxima), incomplete=missing)


def mesh_counts(background, qualification, output):
    table = [dict(scope="operator_background", case=None, full_dofs=background.get("full_dofs") if background else None,
                  free_dofs=background.get("free_dofs") if background else None,
                  public_ports=background.get("public_ports") if background else None,
                  reference_ndofs=None, touched_ports=None,
                  status="RECORDED" if background else "NOT_RUN", source="results/g0_operators/BACKGROUND.json")]
    indexed = index_rows(qualification, ("case",))
    for case in QUALIFICATION_CASES:
        source = indexed.get((case,), {})
        table.append(dict(scope="geometry_qualification", case=case, full_dofs=source.get("ndofs"), free_dofs=None,
                          public_ports=source.get("target_ports"), reference_ndofs=source.get("reference_ndofs"),
                          touched_ports=source.get("touched_ports"), status="RECORDED" if source else "NOT_RUN",
                          source="results/g0_qualification/RESULTS.json"))
    return csv_file(output / "MESH_COUNTS.csv", table,
        ["scope", "case", "status", "full_dofs", "free_dofs", "reference_ndofs", "public_ports", "touched_ports", "source"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "figures/g0_mechanisms")
    parser.add_argument("--summary-dir", type=Path, default=ROOT / "research/delegated/results_support")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.summary_dir.mkdir(parents=True, exist_ok=True)
    issues = []
    operator_dir, stress_dir = ROOT / "results/g0_operators", ROOT / "results/g0_stress"
    operator_rows = read_optional(operator_dir / "RESULTS.json", list, issues)
    plan = read_optional(operator_dir / "PLAN.json", dict, issues)
    setup = read_optional(operator_dir / "SETUP.json", list, issues)
    core = read_optional(operator_dir / "CORE_STABILITY.json", dict, issues)
    stress_rows = read_optional(stress_dir / "RESULTS.json", list, issues)
    complete = read_optional(stress_dir / "COMPLETE.json", dict, issues)
    background = read_optional(operator_dir / "BACKGROUND.json", dict, issues)
    qualification = read_optional(ROOT / "results/g0_qualification/RESULTS.json", list, issues)
    figures, csvs, operator_ranges = operators(operator_rows, plan, setup, core, args.output_dir)
    indexed, stress_figures, stress_csvs, spectrum_ranges = stress_spectra(stress_rows, args.output_dir)
    finite_figures, finite_csvs, finite_ranges = finite_diagnostics(indexed, args.output_dir)
    future_figure, future_csvs, future_ranges = c4_future(indexed, args.output_dir)
    mesh_csv = mesh_counts(background, qualification, args.output_dir)
    figures += stress_figures + finite_figures + [future_figure]
    csvs += stress_csvs + finite_csvs + future_csvs + [mesh_csv]
    summary = dict(status="STORED_DIAGNOSTIC_FIGURE_SNAPSHOT", figures=figures, data_csvs=csvs,
        input_issues=issues, operator_rows_present=len(operator_rows or []), stress_rows_present=len(stress_rows or []),
        stress_complete_record=complete, unique_physical_parents_present=len({row.get("physical_parent") for row in stress_rows or []
                                                                          if row.get("physical_parent") is not None}),
        matched_units=MATCHED_NOTE, operator_ranges=operator_ranges, spectrum_ranges=spectrum_ranges,
        finite_ranges=finite_ranges, future25_ranges=future_ranges,
        core_stability_record=core, physical_solves_executed=False, truth_files_read=False,
        gates_modified=False, posterior_or_image_arrays_created=False,
        cpu_scope="Per-method query solve and setup kernel values are displayed, not summed into wrapper accounting.",
        absent_plot_policy="Only CREATED entries are current; any preexisting files listed in skipped entries were not refreshed.")
    destination = args.summary_dir / "MECHANISM_FIGURE_SUMMARY.json"
    destination.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(dict(summary=rel(destination), created_figures=sum(item["status"] == "CREATED" for item in figures),
                          skipped_figures=sum(item["status"] != "CREATED" for item in figures), physical_solves_executed=False)))


if __name__ == "__main__":
    main()
