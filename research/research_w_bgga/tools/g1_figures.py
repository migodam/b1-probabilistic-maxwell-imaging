"""Summarize sealed G1 JSON artifacts without physics, truth, or array access.

This is a descriptive postprocessor, not a gate decision writer.  It imports
neither bgga nor a PDE/learning backend and never opens NPZ files.  Running it
replaces only its three named results tables and files in its own delegated
output directory.  Original replay and normalization-audit artifacts stay read
only.  Missing cases retain the frozen 24-case denominator and NOT_RUN status.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import statistics
import textwrap


ROOT = Path(__file__).resolve().parents[1]
METHODS = tuple(f"M{i}" for i in range(7))
PARENTS = tuple(range(12))
REPEATS = (0, 1)
CASE_DENOMINATOR = 24
PARENT_DENOMINATOR = 12
ATOMS = 2950
WALL_GROUPS = 118
CHECKPOINTS = (16, 64, 256, 1024, 2048, 2949, 2950)
BOOTSTRAP_DRAWS = 2000
BOOTSTRAP_SEED = 27102801
OUT = ROOT / "research/delegated/g1_figures"
SHORT_NAMES = {
    "M0": "M0 exact", "M1": "M1 confidence", "M2": "M2 secant",
    "M3": "M3 task risk", "M4": "M4 interval", "M5": "M5 quotient",
    "M6": "M6 cache control",
}
METRICS = (
    "joint_TV", "CDF_x_sup", "CDF_y_sup", "CDF_max_sup",
    "event_x_ge_0_absolute_error", "event_y_ge_2_1_absolute_error",
    "event_max_absolute_error", "important_joint_mode_max_absolute_error",
    "position_Bayes_regret", "position_Bayes_regret_limit",
    "reference_position_Bayes_risk", "position_Bayes_regret_over_limit",
    "logical_cached_queries", "cached_method_cpu_seconds_excluding_IO",
    "actual_new_full_FE", "actual_new_RHS", "M5_cluster_count",
    "M5_cluster_fraction_of_118",
)
BOOTSTRAP_METRICS = (
    "joint_TV", "CDF_max_sup", "event_max_absolute_error",
    "important_joint_mode_max_absolute_error", "position_Bayes_regret",
    "cached_method_cpu_seconds_excluding_IO",
)
CASE_GATES = (
    "joint_TV_gate", "coordinate_CDF_gate", "event_probability_gate",
    "important_joint_mode_gate", "all_probability_gate", "position_risk_gate",
    "all_probability_and_risk_gate", "original_raw_interval_gate",
)
BOUND_FIELDS = {
    "posterior_atom_lower": "atom_lower_signed_violation",
    "posterior_atom_upper": "atom_upper_signed_violation",
    "unknown_mass_upper": "unknown_mass_upper_signed_violation",
    "log_likelihood_lower": "likelihood_lower_signed_violation",
    "log_likelihood_upper": "likelihood_upper_signed_violation",
}
FIGURE_STEMS = (
    "G1_parent_joint_TV_and_regret", "G1_M4_unknown_mass_queries",
    "G1_M5_cluster_count", "G1_cached_method_CPU",
    "G1_raw_bounds_and_normalization_scope",
)


def number(value):
    """Finite JSON numbers only: do not turn null, booleans, or strings into 0."""
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value) if math.isfinite(value) else None
    return None


def maximum(values):
    values = list(values)
    return max(values) if values and all(v is not None for v in values) else None


def scalar_stats(values):
    values = [float(v) for v in values if number(v) is not None]
    return {
        "n": len(values), "mean": statistics.fmean(values) if values else None,
        "median": statistics.median(values) if values else None,
        "min": min(values) if values else None, "max": max(values) if values else None,
    }


def combine_statuses(statuses):
    statuses = list(statuses)
    if "FAIL" in statuses:
        return "FAIL"
    if statuses and all(s == "PASS" for s in statuses):
        return "PASS"
    if not statuses or any(str(s).startswith("NOT_RUN") for s in statuses):
        return "NOT_RUN"
    return "NOT_DEFINED"


def missing_status(status):
    return "NOT_DEFINED" if "NOT_DEFINED" in str(status) else "NOT_RUN"


def threshold_gate(value, limit, row_state):
    if row_state != "EVALUATED":
        return missing_status(row_state)
    if number(value) is None or number(limit) is None:
        return "NOT_DEFINED"
    return "PASS" if value <= limit else "FAIL"


class Sources:
    def __init__(self):
        self.records = []
        self.problems = []

    def problem(self, label, detail):
        self.problems.append({"source": label, "detail": detail})

    def read(self, relative_path, expected_type):
        path = ROOT / relative_path
        # There is deliberately no generic array loader, source-map resolver,
        # directory traversal, or access to historical banks/GT in this tool.
        if path.suffix != ".json" or not path.resolve().is_relative_to(ROOT.resolve()):
            raise ValueError("only campaign-local, explicitly named JSON is allowed")
        record = {"path": str(path), "status": "NOT_RUN"}
        try:
            value = json.loads(path.read_text())
            if not isinstance(value, expected_type):
                raise ValueError(f"expected {expected_type.__name__}, got {type(value).__name__}")
            record["status"] = "READ"
        except (OSError, ValueError) as exc:
            self.problem(relative_path, f"NOT_RUN: {type(exc).__name__}: {exc}")
            value = expected_type()
        self.records.append(record)
        return value, record["status"] == "READ"


def keyed_rows(rows, keys, sources, label):
    result, duplicates = {}, set()
    for row in rows:
        if not isinstance(row, dict) or any(key not in row for key in keys):
            sources.problem(label, "row lacks its explicit case/method identity")
            continue
        identity = tuple(row[key] for key in keys)
        if identity in result:
            duplicates.add(identity)
        result[identity] = row
    for identity in duplicates:
        result.pop(identity, None)
        sources.problem(label, f"duplicate identity {identity}; retained as NOT_RUN")
    return result


def mode_error(score, threshold):
    """Re-select important *predeclared cells* using the actual GATES cutoff."""
    if number(threshold) is None or not isinstance(score.get("joint_modes"), list):
        return None, None
    relevant = []
    for mode in score["joint_modes"]:
        mass = number(mode.get("reference_mass"))
        if mass is None:
            return None, None
        if mass >= threshold:
            error = number(mode.get("absolute_error"))
            if error is None:
                return None, None
            relevant.append(error)
    # Empty important cells are explicit NOT_DEFINED, never vacuous PASS.
    return (max(relevant) if relevant else None), len(relevant)


def normalized_formula(text):
    return str(text).replace(" ", "").replace("<=", "").replace("lambda0", "")


def source_metric_check(row, score, sources, label):
    checks = {
        "joint_TV": score.get("joint_TV"),
        "CDF_x_sup": (score.get("coordinate_CDF_sup") or {}).get("x"),
        "CDF_y_sup": (score.get("coordinate_CDF_sup") or {}).get("y"),
        "event_x_ge_0_absolute_error": (score.get("event_probability_errors") or {}).get("pose_x>=0", {}).get("absolute_error"),
        "event_y_ge_2_1_absolute_error": (score.get("event_probability_errors") or {}).get("pose_y>=2.1", {}).get("absolute_error"),
        "position_Bayes_regret": score.get("position_Bayes_regret"),
        "position_Bayes_regret_limit": score.get("position_Bayes_regret_limit"),
    }
    mismatched = [key for key, value in checks.items()
                  if number(row.get(key)) is not None and number(value) is not None
                  and row[key] != value]
    if mismatched:
        sources.problem(label, "summary/offline mismatch in " + ", ".join(mismatched))
    return mismatched


def summarize_bound_comparison(parent, repeat, comparison, positives, summaries):
    checkpoint = comparison.get("query_count")
    category_statuses = []
    for category, field in BOUND_FIELDS.items():
        raw = comparison.get(field)
        vector = raw if isinstance(raw, list) else [raw]
        finite = [number(value) for value in vector]
        finite_values = [value for value in finite if value is not None]
        negative_infinities = sum(value == "-Infinity" for value in vector)
        positive_infinities = sum(value == "Infinity" for value in vector)
        unavailable = sum(value is None and raw_value not in ("-Infinity", "Infinity")
            for value, raw_value in zip(finite, vector))
        positive_values = []
        for index, value in enumerate(finite):
            if value is not None and value > 0.0:
                positive_values.append(value)
                positives.append({"parent_index": parent, "repeat_index": repeat,
                    "query_count": checkpoint, "category": category,
                    "component_index": index if isinstance(raw, list) else "SCALAR",
                    "raw_positive_signed_violation": value})
        for index, value in enumerate(vector):
            if value == "Infinity":
                positives.append({"parent_index": parent, "repeat_index": repeat,
                    "query_count": checkpoint, "category": category,
                    "component_index": index if isinstance(raw, list) else "SCALAR",
                    "raw_positive_signed_violation": "Infinity"})
        complete = raw is not None and (not isinstance(raw, list) or bool(raw)) and unavailable == 0
        status = "FAIL" if positive_values or positive_infinities else "PASS" if complete else "NOT_RUN"
        category_statuses.append(status)
        summaries.append({"parent_index": parent, "repeat_index": repeat,
            "query_count": checkpoint, "category": category, "status": status,
            "number_of_components": len(vector) if complete else 0,
            "finite_components": len(finite_values),
            "nonfinite_or_unavailable_components": len(vector) - len(finite_values),
            "negative_infinite_components": negative_infinities,
            "positive_infinite_components": positive_infinities,
            "unavailable_components": unavailable,
            "negative_components": sum(v < 0 for v in finite_values),
            "exact_zero_components": sum(v == 0 for v in finite_values),
            "positive_components": len(positive_values) + positive_infinities,
            "minimum_raw_finite_signed_violation": min(finite_values) if finite_values else None,
            "maximum_raw_finite_signed_violation": max(finite_values) if finite_values else None,
            "minimum_positive_signed_violation": min(positive_values) if positive_values else None,
            "median_positive_signed_violation": statistics.median(positive_values) if positive_values else None,
            "maximum_positive_signed_violation": max(positive_values) if positive_values else None,
            "threshold": "strictly > 0, no roundoff tolerance or clipping",
        })
    observed = combine_statuses(category_statuses)
    flag = comparison.get("any_positive_raw_bound_violation")
    if flag is True:
        return "FAIL"
    if flag is False and observed == "PASS":
        return "PASS"
    return observed if observed == "FAIL" else "NOT_RUN"


def load_cases(sources, gates, plan, score_rows, seal):
    score_map = keyed_rows(score_rows, ("parent_index", "repeat_index", "method_id"), sources, "SCORE_TABLE.json")
    timing_map = keyed_rows(seal.get("method_timings", []), ("parent_index", "repeat_index", "method_id"), sources, "seal.method_timings")
    query_map = keyed_rows(seal.get("logical_queries", []), ("parent_index", "repeat_index", "method_id"), sources, "seal.logical_queries")
    expected = {(p, r, m) for p in PARENTS for r in REPEATS for m in METHODS}
    unexpected = sorted(set(score_map) - expected)
    if unexpected:
        sources.problem("SCORE_TABLE.json", f"unexpected identities excluded from frozen grid: {unexpected}")
    rows, checkpoint_rows, bound_positives, bound_summaries, reference_rows = [], [], [], [], []
    risk_formula = gates.get("position_Bayes_regret")
    saved_risk_formula = (plan.get("scoring") or {}).get("Bayes_squared_position_regret_limit")
    risk_formula_compatible = bool(risk_formula and saved_risk_formula
        and normalized_formula(risk_formula) == normalized_formula(saved_risk_formula))
    if not risk_formula_compatible:
        sources.problem("GATES.G1.position_Bayes_regret", "saved risk convention absent or differs; risk gate NOT_DEFINED")

    for parent in PARENTS:
        for repeat in REPEATS:
            case_name = f"parent_{parent:02d}_repeat_{repeat:02d}"
            case_path = f"results/g1_replay_001/{case_name}"
            offline, offline_ok = sources.read(f"{case_path}/OFFLINE_SCORES.json", list)
            offline_map = keyed_rows(offline, ("method_id",), sources, f"{case_name}/OFFLINE_SCORES.json")
            reference, reference_ok = sources.read(f"{case_path}/OFFLINE_REFERENCE.json", dict)
            reference_rows.append({"parent_index": parent, "repeat_index": repeat,
                "status": "RECORDED" if reference_ok else "NOT_RUN",
                "logZ_exact": number(reference.get("logZ_exact")),
                "p0min": number(reference.get("p0min")),
                "default_unknown_envelope_floor": number(reference.get("p0min_over_Zexact_plus_p0min")),
                "logical_cached_queries": number(reference.get("logical_queries")),
                "actual_new_full_FE": number(reference.get("actual_new_full_FE")),
                "actual_new_RHS": number(reference.get("actual_new_RHS")),
                "computed_after_all_online_outputs_sealed": reference.get("computed_after_all_online_outputs_sealed"),
                "diagnostic_label": reference.get("diagnostic_label", "NOT_RUN"),
                "interpretation": reference.get("interpretation", "NOT_RUN")})
            for method in METHODS:
                identity = (parent, repeat, method)
                summary = score_map.get(identity, {})
                score = offline_map.get((method,), {})
                metadata, metadata_ok = sources.read(f"{case_path}/{method}.json", dict)
                reported_status = summary.get("status", metadata.get("status", "NOT_RUN"))
                all_complete = bool(summary and score and offline_ok and metadata_ok)
                statuses = [reported_status, score.get("status"), metadata.get("status")]
                state = ("EVALUATED" if all_complete and all(s == "COMPLETED" for s in statuses)
                         else "NOT_DEFINED" if any("NOT_DEFINED" in str(s) for s in statuses)
                         else "NOT_RUN")
                cdf = score.get("coordinate_CDF_sup") or {}
                events = score.get("event_probability_errors") or {}
                event_x = summary.get("event_x>=0_absolute_error", events.get("pose_x>=0", {}).get("absolute_error"))
                event_y = summary.get("event_y>=2.1_absolute_error", events.get("pose_y>=2.1", {}).get("absolute_error"))
                mode, important_count = mode_error(score, number(gates.get("important_joint_mode_min_mass")))
                row = {"parent_index": parent, "repeat_index": repeat, "method_id": method,
                    "reported_method_status": reported_status, "comparison_status": state,
                    "summary_source_present": bool(summary), "offline_score_source_present": bool(score),
                    "method_source_present": metadata_ok, "reference_metadata_present": reference_ok,
                    "joint_TV": number(summary.get("joint_TV", score.get("joint_TV"))),
                    "CDF_x_sup": number(summary.get("CDF_x_sup", cdf.get("x"))),
                    "CDF_y_sup": number(summary.get("CDF_y_sup", cdf.get("y"))),
                    "event_x_ge_0_absolute_error": number(event_x),
                    "event_y_ge_2_1_absolute_error": number(event_y),
                    "important_joint_mode_max_absolute_error": mode,
                    "important_reference_joint_cells": important_count,
                    "saved_important_joint_mode_max_absolute_error": number(summary.get("important_joint_mode_max_absolute_error")),
                    "position_Bayes_regret": number(summary.get("position_Bayes_regret", score.get("position_Bayes_regret"))),
                    "position_Bayes_regret_limit": number(summary.get("position_Bayes_regret_limit", score.get("position_Bayes_regret_limit"))),
                    "reference_position_Bayes_risk": number(score.get("reference_position_Bayes_risk")),
                    "logical_cached_queries": number(metadata.get("logical_query_count", summary.get("logical_queries"))),
                    "cached_method_cpu_seconds_excluding_IO": number(metadata.get("deployment_method_process_cpu_seconds_excluding_IO")),
                    "actual_new_full_FE": number(metadata.get("actual_new_full_FE", summary.get("actual_new_full_FE"))),
                    "actual_new_RHS": number(metadata.get("actual_new_RHS", summary.get("actual_new_RHS"))),
                    "ground_truth_loss_status": score.get("ground_truth_loss_status", summary.get("ground_truth_loss_status", "NOT_RUN")),
                    "floating_point_roundoff_certified": metadata.get("floating_point_roundoff_certified", "NOT_RUN"),
                    "label": metadata.get("label", "NOT_RUN"),
                    "M5_cluster_count": None, "M5_recorded_group_count": None,
                    "M5_cluster_fraction_of_118": None,
                }
                mismatch = source_metric_check(row, score, sources, f"{case_name}/{method}")
                if identity in timing_map:
                    sealed_cpu = number(timing_map[identity].get("process_cpu_seconds_excluding_IO"))
                    if row["cached_method_cpu_seconds_excluding_IO"] is None:
                        row["cached_method_cpu_seconds_excluding_IO"] = sealed_cpu
                    elif sealed_cpu is not None and sealed_cpu != row["cached_method_cpu_seconds_excluding_IO"]:
                        mismatch.append("cached_method_cpu_seconds_excluding_IO")
                if identity in query_map:
                    sealed_queries = number(query_map[identity].get("logical_queries"))
                    if sealed_queries is not None and sealed_queries != row["logical_cached_queries"]:
                        mismatch.append("logical_cached_queries")
                if metadata_ok and metadata.get("method_id") != method:
                    mismatch.append("method_id")
                if mismatch:
                    row["comparison_status"] = "NOT_RUN_SOURCE_MISMATCH"
                    sources.problem(f"{case_name}/{method}", "sealed/summary identity or metric mismatch: " + ", ".join(mismatch))
                state = row["comparison_status"]
                row["CDF_max_sup"] = maximum((row["CDF_x_sup"], row["CDF_y_sup"]))
                row["event_max_absolute_error"] = maximum((row["event_x_ge_0_absolute_error"], row["event_y_ge_2_1_absolute_error"]))
                risk_limit = row["position_Bayes_regret_limit"]
                row["position_Bayes_regret_over_limit"] = (
                    row["position_Bayes_regret"] / risk_limit
                    if row["position_Bayes_regret"] is not None and risk_limit is not None and risk_limit > 0 else None)
                row["joint_TV_gate"] = threshold_gate(row["joint_TV"], number(gates.get("joint_TV")), state)
                row["coordinate_CDF_gate"] = threshold_gate(row["CDF_max_sup"], number(gates.get("coordinate_CDF_sup")), state)
                row["event_probability_gate"] = threshold_gate(row["event_max_absolute_error"], number(gates.get("event_probability")), state)
                row["important_joint_mode_gate"] = threshold_gate(mode, number(gates.get("important_joint_mode_error")), state)
                row["all_probability_gate"] = combine_statuses(row[key] for key in CASE_GATES[:4])
                row["position_risk_gate"] = threshold_gate(row["position_Bayes_regret"], risk_limit if risk_formula_compatible else None, state)
                row["all_probability_and_risk_gate"] = combine_statuses((row["all_probability_gate"], row["position_risk_gate"]))
                row["original_raw_interval_gate"] = "NOT_DEFINED"
                if method == "M4":
                    comparisons = score.get("checkpoint_bound_comparisons") or []
                    comparisons_by_count = keyed_rows(comparisons, ("query_count",), sources, f"{case_name}/M4 checkpoints")
                    metadata_by_count = {c.get("query_count"): c for c in metadata.get("checkpoints", [])}
                    checkpoint_states = []
                    for query_count in CHECKPOINTS:
                        comparison = comparisons_by_count.get((query_count,), {})
                        point = comparison.get("queried_conditional_point_estimate_scores") or {}
                        raw_state = summarize_bound_comparison(parent, repeat,
                            {"query_count": query_count, **comparison}, bound_positives, bound_summaries)
                        checkpoint_states.append(raw_state)
                        checkpoint_meta = metadata_by_count.get(query_count, {})
                        upper = number(comparison.get("unknown_mass_upper", checkpoint_meta.get("unknown_mass_upper")))
                        checkpoint_rows.append({"parent_index": parent, "repeat_index": repeat,
                            "query_count": query_count, "status": "RECORDED" if comparison and checkpoint_meta else "NOT_RUN",
                            "original_raw_interval_gate": raw_state,
                            "unknown_mass_upper": upper,
                            "one_minus_unknown_mass_upper": 1.0 - upper if upper is not None else None,
                            "unknown_mass_upper_is_exactly_one_in_saved_float": upper == 1.0 if upper is not None else "NOT_RUN",
                            "actual_reference_unknown_mass": number(comparison.get("actual_reference_unknown_mass")),
                            "unknown_mass_upper_signed_violation": number(comparison.get("unknown_mass_upper_signed_violation")),
                            "maximum_atom_lower_signed_violation": number(comparison.get("maximum_atom_lower_signed_violation")),
                            "maximum_atom_upper_signed_violation": number(comparison.get("maximum_atom_upper_signed_violation")),
                            "queried_conditional_point_joint_TV": number(point.get("joint_TV")),
                            "queried_conditional_point_position_Bayes_regret": number(point.get("position_Bayes_regret")),
                            "interval_label": comparison.get("interval_label", checkpoint_meta.get("interval_label", "NOT_RUN")),
                            "floating_point_roundoff_certified": comparison.get("floating_point_roundoff_certified", "NOT_RUN"),
                            "scope": "original finite numerical envelope; raw floating point violations retained"})
                    row["original_raw_interval_gate"] = combine_statuses(checkpoint_states)
                    row["M4_unknown_mass_upper_at_stop"] = number((metadata.get("details") or {}).get("unknown_mass_upper_at_stop"))
                    row["M4_unknown_mass_stop_threshold"] = number((metadata.get("details") or {}).get("unknown_mass_stop_threshold"))
                if method == "M5":
                    details = metadata.get("details") or {}
                    count = number(details.get("cluster_count"))
                    row["M5_cluster_count"] = count
                    row["M5_recorded_group_count"] = len(details["group_ids"]) if isinstance(details.get("group_ids"), list) else None
                    row["M5_cluster_fraction_of_118"] = count / WALL_GROUPS if count is not None else None
                    row["M5_future_scope"] = details.get("future_scope", "NOT_RUN")
                    row["M5_other_radius_contrast_qualification"] = details.get("G0_other_radius_contrast_task_cluster_qualification", "NOT_RUN")
                rows.append(row)
    return rows, checkpoint_rows, bound_positives, bound_summaries, reference_rows


def parent_table(case_rows):
    indexed = {(r["parent_index"], r["repeat_index"], r["method_id"]): r for r in case_rows}
    parents = []
    for parent in PARENTS:
        for method in METHODS:
            cases = [indexed[(parent, repeat, method)] for repeat in REPEATS]
            row = {"parent_index": parent, "method_id": method,
                "expected_noise_repeats": 2,
                "evaluated_noise_repeats": sum(c["comparison_status"] == "EVALUATED" for c in cases),
                "aggregation": "arithmetic mean of both noise repeats; neither repeat may be dropped"}
            for metric in METRICS:
                values = [number(c.get(metric)) if c["comparison_status"] == "EVALUATED" else None for c in cases]
                complete = all(v is not None for v in values)
                row[f"{metric}_noise_mean"] = statistics.fmean(values) if complete else None
                row[f"{metric}_noise_min"] = min(values) if complete else None
                row[f"{metric}_noise_max"] = max(values) if complete else None
                row[f"{metric}_status"] = "RECORDED" if complete else "NOT_RUN" if any(c["comparison_status"] != "EVALUATED" for c in cases) else "NOT_DEFINED"
            for gate in CASE_GATES:
                row[f"{gate}_both_repeats"] = combine_statuses(c[gate] for c in cases)
                row[f"{gate}_pass_cases_out_of_2"] = sum(c[gate] == "PASS" for c in cases)
            parents.append(row)
    return parents


def method_table(case_rows, parents, plan):
    rows = []
    for method in METHODS:
        cases = [r for r in case_rows if r["method_id"] == method]
        pairs = [r for r in parents if r["method_id"] == method]
        row = {"method_id": method, "definition": (plan.get("methods") or {}).get(method, "NOT_RUN"),
            "case_denominator": CASE_DENOMINATOR, "parent_denominator": PARENT_DENOMINATOR,
            "evaluated_cases": sum(r["comparison_status"] == "EVALUATED" for r in cases),
            "not_run_cases": sum(str(r["comparison_status"]).startswith("NOT_RUN") for r in cases),
            "not_defined_cases": sum(r["comparison_status"] == "NOT_DEFINED" for r in cases),
            "metric_aggregation": "first average the two noise repeats within each parent; then summarize parents",
            "scientific_decision": "PARENT_REVIEW_REQUIRED", "publication_level_pass": "NOT_ASSESSED"}
        for metric in METRICS:
            stats = scalar_stats(r[f"{metric}_noise_mean"] for r in pairs)
            for name, value in stats.items():
                row[f"{metric}_parent_{name}"] = value
            incomplete_run = any(str(r[f"{metric}_status"]).startswith("NOT_RUN") for r in pairs)
            row[f"{metric}_status"] = (
                "RECORDED_ALL_12_PARENTS" if stats["n"] == PARENT_DENOMINATOR else
                "NOT_RUN_INCOMPLETE_PARENT_PAIRS" if stats["n"] or incomplete_run else "NOT_DEFINED")
        for gate in CASE_GATES:
            counts = Counter(r[gate] for r in cases)
            for state in ("PASS", "FAIL", "NOT_DEFINED", "NOT_RUN"):
                row[f"{gate}_{state}_cases_out_of_24"] = counts[state]
            row[f"{gate}_PASS_fraction_denominator_24"] = counts["PASS"] / CASE_DENOMINATOR
            row[f"{gate}_all_24_case_status"] = combine_statuses(r[gate] for r in cases)
        rows.append(row)
    return rows


def paired_bootstrap(parents, sources):
    records = []
    indexed = {(r["parent_index"], r["method_id"]): r for r in parents}
    try:
        import numpy as np
        draws = np.random.default_rng(BOOTSTRAP_SEED).integers(0, PARENT_DENOMINATOR,
            size=(BOOTSTRAP_DRAWS, PARENT_DENOMINATOR))
    except ImportError as exc:
        np, draws = None, None
        sources.problem("numpy bootstrap", f"NOT_RUN: {exc}")
    for method in METHODS:
        for metric in BOOTSTRAP_METRICS:
            pairs = [(indexed[(p, method)][f"{metric}_noise_mean"], indexed[(p, "M0")][f"{metric}_noise_mean"]) for p in PARENTS]
            count = sum(a is not None and b is not None for a, b in pairs)
            row = {"method_id": method, "baseline_method": "M0", "metric": metric,
                "status": "RECORDED_DESCRIPTIVE_ONLY" if count == PARENT_DENOMINATOR and np is not None else "NOT_RUN",
                "paired_parent_count": count, "required_parent_count": PARENT_DENOMINATOR,
                "noise_repeats_averaged_before_bootstrap": True,
                "resampling_unit": "whole parent, same parent draws for every method/metric",
                "draws": BOOTSTRAP_DRAWS, "seed": BOOTSTRAP_SEED,
                "contrast": "method minus M0, in metric units",
                "mean_paired_difference": None, "descriptive_percentile_2_5": None,
                "descriptive_percentile_97_5": None,
                "interpretation": "exposed development cases; descriptive interval only, no publication-level PASS, no G2 speed gate"}
            if row["status"] == "RECORDED_DESCRIPTIVE_ONLY":
                differences = np.asarray([a - b for a, b in pairs], dtype=float)
                samples = differences[draws].mean(axis=1)
                lo, hi = np.percentile(samples, (2.5, 97.5))
                row.update(mean_paired_difference=float(differences.mean()),
                    descriptive_percentile_2_5=float(lo), descriptive_percentile_97_5=float(hi))
            records.append(row)
    return records


def normalization_audit(sources):
    raw, raw_ok = sources.read("results/g1_normalization_repair/RESULTS.json", list)
    decision, decision_ok = sources.read("results/g1_normalization_repair/DECISION.json", dict)
    indexed = keyed_rows(raw, ("parent", "repeat"), sources, "normalization RESULTS.json")
    rows, case_facts = [], []
    for parent in PARENTS:
        for repeat in REPEATS:
            case = indexed.get((parent, repeat), {})
            checkpoints = keyed_rows(case.get("checkpoints", []), ("query_count",), sources, f"normalization parent {parent} repeat {repeat}")
            complete = bool(case)
            for n in CHECKPOINTS:
                check = checkpoints.get((n,), {})
                lower, upper = number(check.get("lower_violations")), number(check.get("upper_violations"))
                evaluated = bool(check and lower is not None and upper is not None and number(check.get("unknown_mass_upper")) is not None)
                complete &= evaluated
                state = "FAIL" if evaluated and (lower > 0 or upper > 0) else "NO_VIOLATION_OBSERVED" if evaluated else "NOT_RUN"
                rows.append({"parent_index": parent, "repeat_index": repeat, "query_count": n,
                    "status": state, "lower_violations": lower, "upper_violations": upper,
                    "unknown_mass_upper": number(check.get("unknown_mass_upper")),
                    "scope": check.get("scope", "NOT_RUN"),
                    "maximum_original_point_error_to_80_digit_normalization": number(case.get("maximum_original_point_error_to_80_digit_normalization")),
                    "original_queries": number(case.get("original_queries")),
                    "new_physics_queries": number(case.get("new_physics_queries")),
                    "measured_certificate_and_independent_audit_cpu_seconds": number(case.get("measured_certificate_and_independent_audit_cpu_seconds")),
                    "cost_scope": "certificate plus independent audit; not an isolated deployment timing"})
            case_facts.append({"parent_index": parent, "repeat_index": repeat,
                "complete": complete,
                "maximum_original_point_error_to_80_digit_normalization": number(case.get("maximum_original_point_error_to_80_digit_normalization")),
                "measured_certificate_and_independent_audit_cpu_seconds": number(case.get("measured_certificate_and_independent_audit_cpu_seconds")),
                "original_queries": number(case.get("original_queries")),
                "new_physics_queries": number(case.get("new_physics_queries"))})
    complete = bool(raw_ok and decision_ok and len(indexed) == CASE_DENOMINATOR
        and all(r["complete"] for r in case_facts)
        and decision.get("status") == "NORMALIZATION_AUDIT_COMPLETE"
        and decision.get("rows") == CASE_DENOMINATOR)
    observed_positive = any(r["status"] == "FAIL" for r in rows)
    return rows, {"status": "COMPLETE_RECORDED_CONDITIONAL_AUDIT" if complete else "NOT_RUN_INCOMPLETE_AUDIT",
        "complete_all_24_cases_and_7_checkpoints": complete,
        "recorded_decision": decision,
        "corrected_checkpoint_rows": len(rows),
        "observed_positive_corrected_bound_violation": observed_positive,
        "case_facts": case_facts,
        "original_point_error_to_80_digit_normalization": scalar_stats(r["maximum_original_point_error_to_80_digit_normalization"] for r in case_facts),
        "certificate_plus_independent_audit_cpu_seconds": scalar_stats(r["measured_certificate_and_independent_audit_cpu_seconds"] for r in case_facts),
        "scope": "directed normalization conditional on saved float64 log-potentials only",
        "not_a_certificate_for": ["Maxwell/FEM forward error", "formation of float64 Gaussian potential", "continuum PDE"],
        "original_raw_violation_gate_unchanged": True}


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            cells = {}
            for key, value in row.items():
                if value is not None:
                    cells[key] = value
                    continue
                not_run = any(str(row.get(name, "")).startswith("NOT_RUN")
                    for name in ("comparison_status", "status"))
                for metric in METRICS:
                    if key == metric or key.startswith(metric + "_"):
                        not_run |= str(row.get(metric + "_status", "")).startswith("NOT_RUN")
                cells[key] = "NOT_RUN" if not_run else "NOT_DEFINED"
            writer.writerow(cells)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def plotting_context(sources):
    # Keep matplotlib's generated font/cache files inside this tool's owned area.
    os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".matplotlib"))
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
            "axes.spines.top": False, "axes.spines.right": False,
            "savefig.facecolor": "white", "pdf.fonttype": 42})
        return np, plt
    except ImportError as exc:
        sources.problem("figure dependencies", f"NOT_RUN: {exc}; no installation attempted")
        return None, None


def save_figure(fig, stem, plt, manifest):
    files = []
    for extension in ("png", "pdf"):
        path = OUT / f"{stem}.{extension}"
        fig.savefig(path, dpi=200 if extension == "png" else None, bbox_inches="tight")
        files.append(str(path))
    plt.close(fig)
    manifest.append({"stem": stem, "status": "GENERATED_REQUIRES_VISUAL_REVIEW", "paths": files,
        "scientific_acceptance": "NOT_ASSESSED"})


def plot_parent_scores(parents, gates, np, plt, manifest):
    fig, axes = plt.subplots(2, 1, figsize=(14.5, 8.5), sharex=True)
    colors = plt.get_cmap("tab10").colors
    markers = ("o", "s", "^", "D", "v", "P", "X")
    for j, method in enumerate(METHODS):
        pairs = [r for r in parents if r["method_id"] == method]
        x = np.asarray(PARENTS) + 1 + (j - 3) * 0.065
        for axis, metric in zip(axes, ("joint_TV", "position_Bayes_regret")):
            y = np.asarray([r[f"{metric}_noise_mean"] if r[f"{metric}_noise_mean"] is not None else np.nan for r in pairs])
            low = np.asarray([r[f"{metric}_noise_min"] if r[f"{metric}_noise_min"] is not None else np.nan for r in pairs])
            high = np.asarray([r[f"{metric}_noise_max"] if r[f"{metric}_noise_max"] is not None else np.nan for r in pairs])
            axis.errorbar(x, y, yerr=(y - low, high - y), fmt=markers[j], ms=5,
                color=colors[j], alpha=.9, capsize=2, linewidth=.8, label=SHORT_NAMES[method])
    tv_limit = number(gates.get("joint_TV"))
    if tv_limit is not None:
        axes[0].axhline(tv_limit, color="black", linestyle="--", lw=1, label="case TV threshold")
    baseline = [r for r in parents if r["method_id"] == "M0"]
    limits = [r["position_Bayes_regret_limit_noise_mean"] for r in baseline]
    axes[1].plot(np.arange(1, 13), [v if v is not None else np.nan for v in limits],
        color="black", linestyle="--", lw=1, label="mean of recorded case regret limits")
    axes[0].set_yscale("symlog", linthresh=1e-5)
    axes[1].set_yscale("symlog", linthresh=1e-9)
    axes[0].set_ylabel("Joint TV (symlog; exact zero shown)")
    axes[1].set_ylabel("Squared position Bayes regret (symlog)")
    axes[1].set_xlabel("Parent (two noise repeats per parent)")
    axes[1].set_xticks(range(1, 13))
    for axis in axes:
        axis.grid(axis="y", alpha=.2)
        axis.set_xlim(.5, 12.5)
    axes[0].legend(ncol=4, loc="upper left", fontsize=9)
    axes[1].legend(loc="upper left", fontsize=9)
    fig.suptitle("G1 exposed development replay: all seven methods, by parent", fontsize=15)
    fig.text(.01, .015, "Markers: mean of both noise repeats; whiskers: their full range. Horizontal offsets reveal overlapping exact-zero methods.\nGates are applied to each original case before aggregation; parent averages do not establish a case or publication-level PASS.", fontsize=9)
    fig.tight_layout(rect=(0, .07, 1, .95))
    save_figure(fig, FIGURE_STEMS[0], plt, manifest)


def plot_m4(case_rows, checkpoints, gates, np, plt, manifest):
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    colors = plt.get_cmap("tab20").colors
    positive_gaps = [r["one_minus_unknown_mass_upper"] for r in checkpoints
        if r["query_count"] <= 2949 and r["one_minus_unknown_mass_upper"] is not None and r["one_minus_unknown_mass_upper"] > 0]
    display_floor = min(positive_gaps) / 10 if positive_gaps else 1e-17
    zero_count = 0
    for p in PARENTS:
        for repeat in REPEATS:
            rows = [r for r in checkpoints if r["parent_index"] == p and r["repeat_index"] == repeat]
            upper = [r["unknown_mass_upper"] if r["unknown_mass_upper"] is not None else np.nan for r in rows]
            style = "-" if repeat == 0 else ":"
            axes[0].plot(range(len(CHECKPOINTS)), upper, style, color=colors[p], alpha=.55, marker=".", lw=.8)
            before = rows[:-1]
            gaps = [r["one_minus_unknown_mass_upper"] for r in before]
            axes[1].plot(range(6), [g if g is not None and g > 0 else np.nan for g in gaps],
                style, color=colors[p], alpha=.65, marker=".", lw=.8)
            for j, gap in enumerate(gaps):
                if gap == 0:
                    axes[1].scatter(j, display_floor, marker="v", color=colors[p], alpha=.5, s=17)
                    zero_count += 1
    axes[0].set_xticks(range(7), [str(n) for n in CHECKPOINTS], rotation=45)
    axes[0].set_ylim(-.04, 1.04)
    # The replay's stopping rule is recorded separately from the TV metric gate.
    stop_thresholds = {r.get("M4_unknown_mass_stop_threshold") for r in case_rows
        if r["method_id"] == "M4" and r.get("M4_unknown_mass_stop_threshold") is not None}
    for threshold in sorted(stop_thresholds):
        axes[0].axhline(threshold, color="black", linestyle="--", lw=.8,
            label=f"recorded M4 stop threshold {threshold:g}")
    if stop_thresholds:
        axes[0].legend(fontsize=8, loc="lower left")
    axes[0].set_title("Saved unknown-mass upper U")
    axes[0].set_ylabel("U (original floating point envelope)")
    axes[0].set_xlabel("Logical cached queries at checkpoint")
    axes[1].set_yscale("log")
    axes[1].set_ylim(bottom=display_floor / 2)
    axes[1].set_xticks(range(6), [str(n) for n in CHECKPOINTS[:-1]], rotation=45)
    axes[1].set_title("Near one, before final enumeration")
    axes[1].set_ylabel("1 - U (log scale; smaller means U nearer 1)")
    axes[1].set_xlabel("Logical cached queries (through 2949)")
    m4 = [r for r in case_rows if r["method_id"] == "M4"]
    for repeat, marker in ((0, "o"), (1, "x")):
        subset = [r for r in m4 if r["repeat_index"] == repeat]
        axes[2].scatter([r["parent_index"] + 1 + (.08 if repeat else -.08) for r in subset],
            [r["logical_cached_queries"] if r["logical_cached_queries"] is not None else np.nan for r in subset],
            marker=marker, label=f"noise {repeat}")
    axes[2].axhline(ATOMS, color="black", linestyle="--", lw=1, label="2950 finite atoms")
    axes[2].set_ylim(0, ATOMS * 1.07)
    axes[2].set_xticks(range(1, 13))
    axes[2].set_title("Recorded M4 stop query count")
    axes[2].set_xlabel("Parent")
    axes[2].set_ylabel("Logical cached queries, each case")
    axes[2].legend(fontsize=8)
    for axis in axes:
        axis.grid(alpha=.15)
    fig.suptitle("M4 default [0,1] unknown envelope and query burden", fontsize=15)
    complete_traces = sum(all(r["status"] == "RECORDED" for r in checkpoints
        if r["parent_index"] == p and r["repeat_index"] == repeat)
        for p in PARENTS for repeat in REPEATS)
    fig.text(.01, .01, f"Complete recorded case traces: {complete_traces}/24; missing checkpoints remain NOT_RUN. Saved U=1 is marked at the display floor ({zero_count} entries); CSV values are unchanged.\nQueries are cached logical accesses; actual new FE/RHS are recorded separately. Original bound violations remain in the raw audit.", fontsize=9)
    fig.tight_layout(rect=(0, .09, 1, .94))
    save_figure(fig, FIGURE_STEMS[1], plt, manifest)


def plot_m5(case_rows, np, plt, manifest):
    fig, axis = plt.subplots(figsize=(12, 5))
    for repeat, marker in ((0, "o"), (1, "x")):
        rows = [r for r in case_rows if r["method_id"] == "M5" and r["repeat_index"] == repeat]
        axis.scatter([r["parent_index"] + 1 + (.09 if repeat else -.09) for r in rows],
            [r["M5_cluster_count"] if r["M5_cluster_count"] is not None else np.nan for r in rows],
            s=40, marker=marker, label=f"noise {repeat}")
    axis.axhline(WALL_GROUPS, linestyle="--", color="black", lw=1, label="118 original wall groups")
    axis.set(xticks=range(1, 13), xlabel="Parent", ylabel="M5 wall cluster count / 118 original groups",
        ylim=(0, WALL_GROUPS * 1.09), title="M5 complete-link quotient size on the exposed finite prior")
    axis.grid(axis="y", alpha=.2)
    axis.legend(loc="lower right")
    fig.text(.01, .01, "Scope: the recorded 25 prior pose nodes with fixed target intrinsic parameters and gamma.\nCluster construction and all table queries belong to M5's deployment work; this figure makes no physical acceleration claim.", fontsize=9)
    fig.tight_layout(rect=(0, .09, 1, 1))
    save_figure(fig, FIGURE_STEMS[2], plt, manifest)


def plot_cpu(case_rows, parents, np, plt, manifest):
    fig, axis = plt.subplots(figsize=(12, 5.5))
    colors = plt.get_cmap("tab10").colors
    any_value = False
    for j, method in enumerate(METHODS):
        values = [r["cached_method_cpu_seconds_excluding_IO_noise_mean"] for r in parents if r["method_id"] == method]
        finite = [(i, value * 1000) for i, value in enumerate(values) if value is not None and value > 0]
        if not finite:
            axis.text(j, .02, "NOT_RUN", rotation=90, transform=axis.get_xaxis_transform(), ha="center")
            continue
        any_value = True
        y = [value for _, value in finite]
        axis.scatter([j + (i - 5.5) * .023 for i, _ in finite], y, color=colors[j], alpha=.7, s=28)
        axis.plot([j - .22, j + .22], [statistics.fmean(y)] * 2, color=colors[j], lw=3)
    if any_value:
        axis.set_yscale("log")
    axis.set_xticks(range(7), [SHORT_NAMES[m] for m in METHODS])
    axis.set_ylabel("Method process CPU ms per case (parent noise mean; log scale)")
    axis.set_title("Cached-array method CPU, excluding I/O", fontsize=15)
    axis.grid(axis="y", alpha=.2)
    zero_physics = len(case_rows) == 168 and all(r["actual_new_full_FE"] == 0 and r["actual_new_RHS"] == 0 for r in case_rows)
    physics_text = "Recorded new FE/RHS = 0 is a replay property, not a physical speedup." if zero_physics else "New FE/RHS records incomplete or nonzero: inspect factual summary."
    fig.text(.01, .01, "Each dot is a complete parent mean after averaging both noise repeats (up to 12); horizontal marks average complete parents.\nArrays were already loaded: NOT cold end-to-end, NOT PDE-speed evidence. " + physics_text + "\nReplay receipt includes separate online I/O and offline reference work; shared G0 and historical bank setup remain paid costs.", fontsize=9)
    fig.tight_layout(rect=(0, .13, 1, 1))
    save_figure(fig, FIGURE_STEMS[3], plt, manifest)


def plot_bounds(positives, checkpoint_rows, repair_rows, repair, np, plt, manifest):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), gridspec_kw={"width_ratios": (1.05, 1.05, 1.2)})
    categories = ("posterior_atom_lower", "posterior_atom_upper", "unknown_mass_upper")
    samples = {name: [math.log10(r["raw_positive_signed_violation"]) for r in positives
        if r["category"] == name and number(r["raw_positive_signed_violation"]) is not None] for name in categories}
    pooled = [value for values in samples.values() for value in values]
    if pooled:
        lower, upper = min(pooled), max(pooled)
        bins = np.linspace(lower - .1, upper + .1, 45 if upper > lower else 3)
        for name, values in samples.items():
            if values:
                axes[0].hist(values, bins=bins, histtype="step", lw=1.5, label=f"{name} (n={len(values)})")
        axes[0].set_yscale("log")
        axes[0].legend(fontsize=8, loc="upper left")
    else:
        axes[0].text(.5, .5, "No positive posterior-bound violation observed\nor original audit NOT_RUN", ha="center", transform=axes[0].transAxes)
    axes[0].set_title("Original raw positive violations")
    axes[0].set_xlabel("log10(raw positive probability-bound violation)")
    axes[0].set_ylabel("Component count (log scale)")
    raw_counts, corrected_counts = [], []
    for n in CHECKPOINTS:
        raw_counts.append(sum(r["query_count"] == n and r["original_raw_interval_gate"] == "FAIL" for r in checkpoint_rows))
        corrected_counts.append(sum(r["query_count"] == n and r["status"] == "FAIL" for r in repair_rows))
    x = np.arange(len(CHECKPOINTS))
    axes[1].bar(x - .17, raw_counts, width=.34, label="original raw / F64 reference", color="#b75346")
    if repair["complete_all_24_cases_and_7_checkpoints"]:
        axes[1].bar(x + .17, corrected_counts, width=.34, label="directed / 80-digit normalization", color="#277f74")
        # Draw visible zero-count markers without altering the reported counts.
        axes[1].scatter(x + .17, corrected_counts, marker="_", color="#277f74", s=100)
    else:
        axes[1].text(.5, .9, "Corrected audit: NOT_RUN / incomplete", transform=axes[1].transAxes, ha="center", fontsize=9)
    axes[1].set(xticks=x, xticklabels=[str(n) for n in CHECKPOINTS], ylim=(-.6, 25),
        title="Cases with any positive violation / 24", xlabel="Logical cached query checkpoint", ylabel="Case count")
    axes[1].tick_params(axis="x", rotation=45)
    axes[1].legend(fontsize=8, loc="upper left")
    axes[2].axis("off")
    point = repair["original_point_error_to_80_digit_normalization"]
    point_range = f"{point['min']:.3g} to {point['max']:.3g}" if point["n"] else "NOT_RUN"
    ll_positive = sum(r["category"].startswith("log_likelihood") for r in positives)
    scope = ("Normalization audit scope\n\n"
        f"Status: {repair['status']}\n"
        "Directed bounds are conditional on saved float64 log-potentials.\n\n"
        "Original point error versus 80-digit normalization:\n" + point_range + "\n"
        "(a normalization diagnostic, not a forward-error bound).\n\n"
        f"Positive raw log-likelihood components: {ll_positive}.\n"
        "No positive raw violation is rounded or clipped away.\n\n"
        "Does not certify Maxwell/FEM forward error, formation of the float64 Gaussian potential, or the continuum PDE.\n\n"
        "Online use adds certificate cost. Original raw results and their FAIL status stay unchanged.")
    axes[2].text(0, 1, "\n".join(textwrap.fill(line, 53) if line else "" for line in scope.splitlines()),
        va="top", fontsize=9, transform=axes[2].transAxes)
    fig.suptitle("M4 original floating point audit and conditional directed normalization", fontsize=15)
    fig.text(.01, .01, "Counts compare distinct recorded numerical audits: original float64 reference versus independent 80-digit normalization.\nThese descriptive plots neither replace the original raw gate nor establish Maxwell certification or a publication-level PASS.", fontsize=9)
    fig.tight_layout(rect=(0, .08, 1, .94))
    save_figure(fig, FIGURE_STEMS[4], plt, manifest)


def make_figures(case_rows, parents, checkpoints, positives, repair_rows, repair, gates, sources, skip):
    manifest = []
    np, plt = (None, None) if skip else plotting_context(sources)
    if plt is None:
        return [{"stem": stem, "status": "NOT_RUN_BY_OPTION" if skip else "NOT_RUN_DEPENDENCY_UNAVAILABLE",
            "paths": [], "scientific_acceptance": "NOT_ASSESSED"} for stem in FIGURE_STEMS]
    jobs = (
        (FIGURE_STEMS[0], lambda: plot_parent_scores(parents, gates, np, plt, manifest)),
        (FIGURE_STEMS[1], lambda: plot_m4(case_rows, checkpoints, gates, np, plt, manifest)),
        (FIGURE_STEMS[2], lambda: plot_m5(case_rows, np, plt, manifest)),
        (FIGURE_STEMS[3], lambda: plot_cpu(case_rows, parents, np, plt, manifest)),
        (FIGURE_STEMS[4], lambda: plot_bounds(positives, checkpoints, repair_rows, repair, np, plt, manifest)),
    )
    for stem, job in jobs:
        try:
            job()
        except Exception as exc:
            plt.close("all")
            sources.problem(stem, f"NOT_RUN: {type(exc).__name__}: {exc}")
            manifest.append({"stem": stem, "status": "NOT_RUN_FIGURE_ERROR", "paths": [], "error": str(exc),
                "scientific_acceptance": "NOT_ASSESSED"})
    return manifest


def factual_summary(case_rows, parents, methods, checkpoints, positives, summaries, references,
                    repair, receipt, seal, plan, gates, sources, figures, bootstrap):
    by_category = {}
    for category in BOUND_FIELDS:
        values = [r["raw_positive_signed_violation"] for r in positives if r["category"] == category]
        by_category[category] = scalar_stats(values)
    before_final = [r for r in checkpoints if r["query_count"] <= 2949]
    m4 = [r for r in case_rows if r["method_id"] == "M4"]
    m5 = [r for r in case_rows if r["method_id"] == "M5"]
    all_original_queries = [r["logical_cached_queries"] for r in case_rows]
    return {
        "schema": "BGGA_G1_FACTUAL_POSTPROCESS_V1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "FACTUAL_SUMMARY_REQUIRES_PARENT_SCIENTIFIC_REVIEW",
        "scope": "saved exposed finite-prior replay JSON only; no new physics, training, GT, or NPZ access",
        "data_status": plan.get("data_status", "NOT_RUN"),
        "expected_design": {"parents": 12, "noise_repeats_per_parent": 2, "cases": 24, "methods": list(METHODS)},
        "case_gate_denominator": CASE_DENOMINATOR,
        "mean_and_median_aggregation": "two noise repeats averaged within each parent before method summaries; missing repeat never dropped",
        "partial_numeric_summary_policy": "available complete-parent descriptive values retain their explicit n; all pass fractions retain denominator 24; paired bootstrap requires all 12 parents",
        "gates_source": str(ROOT / "GATES.json"), "G1_thresholds_verbatim": gates,
        "risk_gate_source": "recorded per-case regret limit, after saved replay formula compatibility with actual GATES; no new lambda0 is inferred",
        "gate_interpretation": "per-case engineering threshold comparisons on exposed development data only; original exact-bound violation uses strict >0 with no tolerance",
        "important_joint_modes": "actual GATES reference-mass cutoff applied to predeclared joint partition cells; not continuous mode discovery",
        "case_method_rows": len(case_rows), "parent_method_rows": len(parents),
        "evaluated_case_method_rows": sum(r["comparison_status"] == "EVALUATED" for r in case_rows),
        "method_facts": methods,
        "paired_parent_bootstrap": {"draws": BOOTSTRAP_DRAWS, "seed": BOOTSTRAP_SEED,
            "same_parent_draws_for_all_methods": True, "baseline": "M0", "interval": "descriptive percentile [2.5%,97.5%] of paired parent mean differences",
            "status_counts": dict(Counter(r["status"] for r in bootstrap)),
            "publication_level_pass": "NOT_ASSESSED", "G2_speed_gate": "NOT_RUN"},
        "accounting": {
            "recorded_run_receipt": receipt, "recorded_seal_case_count": seal.get("case_count", "NOT_RUN"),
            "all_online_outputs_written_before_reference": seal.get("all_method_outputs_written", "NOT_RUN"),
            "reference_input_to_online_methods": seal.get("reference_input_to_online_methods", "NOT_RUN"),
            "all_168_case_methods_record_2950_logical_cached_queries": len(all_original_queries) == 168 and all(q == ATOMS for q in all_original_queries),
            "observed_logical_cached_queries": scalar_stats(all_original_queries),
            "actual_new_full_FE_from_case_method_records": scalar_stats(r["actual_new_full_FE"] for r in case_rows),
            "actual_new_RHS_from_case_method_records": scalar_stats(r["actual_new_RHS"] for r in case_rows),
            "new_PDE_speed_evidence": False, "cold_end_to_end_timing": "NOT_RUN",
            "method_CPU_scope": "cached arrays already loaded; method process CPU excludes I/O; no PDE-speed interpretation",
            "historical_and_G0_setup": "receipt setup snapshots retained separately; not recharged per method and not treated as free",
            "normalization_audit_cost": "additional certificate plus independent-audit CPU, not isolated online certificate CPU; not added to original timing or treated as free",
        },
        "M4": {"stop_logical_cached_queries": scalar_stats(r["logical_cached_queries"] for r in m4),
            "recorded_unknown_mass_stop_threshold": scalar_stats(r.get("M4_unknown_mass_stop_threshold") for r in m4),
            "pre_final_checkpoint_upper_U": scalar_stats(r["unknown_mass_upper"] for r in before_final),
            "pre_final_saved_U_exactly_one_entries": sum(r["unknown_mass_upper"] == 1 for r in before_final),
            "pre_final_complement_1_minus_U": scalar_stats(r["one_minus_unknown_mass_upper"] for r in before_final),
            "original_raw_interval_case_status_counts": dict(Counter(r["original_raw_interval_gate"] for r in m4)),
            "raw_positive_violation_components_by_category": by_category,
            "raw_checkpoint_rows": len(checkpoints),
            "raw_signed_comparison_rows": len(summaries),
            "original_raw_values_preserved": True},
        "M5": {"cluster_count": scalar_stats(r["M5_cluster_count"] for r in m5),
            "cluster_fraction_of_frozen_118_groups": scalar_stats(r["M5_cluster_fraction_of_118"] for r in m5),
            "recorded_group_count": scalar_stats(r["M5_recorded_group_count"] for r in m5),
            "future_scope": sorted({r.get("M5_future_scope", "NOT_RUN") for r in m5}),
            "other_radius_contrast_qualification": sorted({r.get("M5_other_radius_contrast_qualification", "NOT_RUN") for r in m5})},
        "normalization_repair": repair,
        "offline_reference_metadata": references,
        "figures": figures, "figure_visual_review": "NOT_RUN_BY_POSTPROCESSOR",
        "input_inventory": sources.records, "input_or_processing_problems": sources.problems,
        "source_files_modified": False, "gates_modified": False,
        "scientific_decision": "PARENT_REVIEW_REQUIRED", "publication_level_pass": "NOT_ASSESSED",
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-figures", action="store_true", help="write tables but explicitly mark figures NOT_RUN")
    args = parser.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    sources = Sources()
    gates_doc, _ = sources.read("GATES.json", dict)
    gates = gates_doc.get("G1", {})
    plan, _ = sources.read("results/g1_replay_001/PLAN.json", dict)
    score_rows, _ = sources.read("results/g1_replay_001/SCORE_TABLE.json", list)
    receipt, _ = sources.read("results/g1_replay_001/RUN_RECEIPT.json", dict)
    seal, _ = sources.read("results/g1_replay_001/ALL_ONLINE_OUTPUTS_SEALED.json", dict)
    cases, checkpoints, positives, bounds, references = load_cases(sources, gates, plan, score_rows, seal)
    parents = parent_table(cases)
    methods = method_table(cases, parents, plan)
    bootstrap = paired_bootstrap(parents, sources)
    repair_rows, repair = normalization_audit(sources)
    write_csv(ROOT / "results/G1_METHOD_TABLE.csv", methods)
    write_csv(ROOT / "results/G1_PARENT_TABLE.csv", parents)
    tables = {
        "G1_CASE_TABLE.csv": cases, "G1_PAIRED_PARENT_BOOTSTRAP.csv": bootstrap,
        "G1_M4_CHECKPOINTS.csv": checkpoints, "G1_RAW_BOUND_POSITIVES.csv": positives,
        "G1_RAW_BOUND_SUMMARY.csv": bounds, "G1_NORMALIZATION_AUDIT.csv": repair_rows,
        "G1_OFFLINE_REFERENCE_METADATA.csv": references,
    }
    for name, rows in tables.items():
        write_csv(OUT / name, rows)
    figures = make_figures(cases, parents, checkpoints, positives, repair_rows, repair, gates, sources, args.skip_figures)
    write_json(OUT / "FIGURE_MANIFEST.json", figures)
    summary = factual_summary(cases, parents, methods, checkpoints, positives, bounds, references,
        repair, receipt, seal, plan, gates, sources, figures, bootstrap)
    write_json(ROOT / "results/G1_FACTUAL_SUMMARY.json", summary)
    print(json.dumps({"status": summary["status"], "method_rows": len(methods), "parent_rows": len(parents),
        "case_rows": len(cases), "figure_statuses": dict(Counter(r["status"] for r in figures)),
        "input_or_processing_problem_count": len(sources.problems), "output": str(OUT),
        "scientific_decision": "PARENT_REVIEW_REQUIRED"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
