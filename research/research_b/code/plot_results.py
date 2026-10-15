"""Plot saved Research B results only.  No solver, training or inference calls.

Usage from the repository root:
  ../../.venv_nn/bin/python -B research/research_b/code/plot_results.py --stage E0
  ... --stage E1                 # after the saved campaign is complete

E0 is an exposed seven-case panel, not a blind population experiment.  Mode 1
and Mode 2 refer to the selection pool.  EIG/risk-driven action choices use the
declared finite proxy reference even in the Mode 1 selection-pool panel.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import time

_WALL_START = time.perf_counter()
_CPU_START = time.process_time()
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[1]
POLICIES = ["no_measurement", "none", "random", "predetermined_second_frequency", "fixed2.4",
            "fixed_1.2", "fixed_1.6", "fixed_2", "fixed_2.4", "fixed_2.8", "fixed_3.2",
            "uniform_one_block", "uniform", "raw", "classical_Fisher", "joint_projected_Fisher",
            "event_EIG", "expected_material_risk", "robust", "conflict", "conflict_diagnostic",
            "calibration_only", "calibrate_then_risk", "geometry_only", "external_geometry",
            "risk_cost_stop", "robust_abstain", "two_step_bundle"]
POLICY_LABELS = {"no_measurement": "No new data", "none": "No new data", "random": "Random",
    "predetermined_second_frequency": "Fixed kR 2.4", "fixed2.4": "Fixed kR 2.4",
    "uniform_one_block": "Uniform (1 block)", "uniform": "Uniform",
    "raw": "Raw separation", "classical_Fisher": "Fisher",
    "joint_projected_Fisher": "Joint Fisher", "event_EIG": "Event EIG*",
    "expected_material_risk": "Material risk*", "robust": "Robust set",
    "conflict": "Conflict diagnostic", "conflict_diagnostic": "Conflict diagnostic",
    "calibration_only": "Calibration only", "calibrate_then_risk": "Calibrate + risk*",
    "geometry_only": "Geometry only", "external_geometry": "External geometry",
    "risk_cost_stop": "Risk cost stop*", "robust_abstain": "Robust abstain", "two_step_bundle": "Two-step risk*",
    **{f"fixed_{k}": "Fixed kR " + k for k in ("1.2", "1.6", "2", "2.4", "2.8", "3.2")}}
ACTION_LABELS = {"no_measurement": "No new data", "repeat_1.6": "Repeat 1.6",
    "frequency_2.4": "Full 2.4", "frequency_2.4_pol0": "2.4 / pol0",
    "frequency_2.4_front": "2.4 / front", "reference_calibration": "Calibration",
    "external_geometry": "External geometry"}
CASE_LABELS = {"solid": "Solid", "small_cavity": "Small cavity", "large_cavity": "Large cavity",
    "annular_cavity": "Annular cavity", "combined": "Combined", "prior_boundary": "Prior boundary",
    "gain_negative_control": "Gain control"}


def ordered(values, preferred):
    return [x for x in preferred if x in values] + sorted(set(values) - set(preferred))


def label(value, mapping):
    return mapping.get(value, str(value).replace("frequency_", "kR ").replace("_", " "))


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def save(fig, out, stem, note):
    fig.text(.01, .012, note, fontsize=8, va="bottom", color="#364152")
    fig.tight_layout(rect=(0, .075, 1, .96))
    for suffix in ("png", "pdf"):
        fig.savefig(out / f"{stem}.{suffix}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def heatmap(ax, matrix, rows, cols, *, cmap="viridis", fmt=".3f", vmin=None, vmax=None, norm=None):
    matrix = np.asarray(matrix, dtype=float)
    image = ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax, norm=norm)
    ax.set_xticks(range(len(cols)), cols, rotation=35, ha="right")
    ax.set_yticks(range(len(rows)), rows)
    if matrix.size <= 140:
        for i in range(len(rows)):
            for j in range(len(cols)):
                if np.isfinite(matrix[i, j]):
                    rgba = image.cmap(image.norm(matrix[i, j]))
                    lum = .21 * rgba[0] + .72 * rgba[1] + .07 * rgba[2]
                    ax.text(j, i, format(matrix[i, j], fmt), ha="center", va="center",
                            fontsize=8, color="white" if lum < .5 else "#141b23")
    return image


def ranking_data(stage, results):
    if stage == "E0":
        comparisons = read(results / "E0_CRITERION_COMPARISON.json")
        risk = read(results / "E0_ACTION_RANKINGS.json")
        picks = [(case, values["selected_actions"]) for case, values in comparisons.items()]
        rows = [{"case": r["case"], "action": r["action"], "risk_reduction": r["risk_reduction"]} for r in risk]
    else:
        source = read(results / f"{stage}_RANKINGS.json")
        source = [r for r in source if r.get("step", 0) == 0]
        unique = {}
        for r in source:
            key = (r["parent"], r["noise"])
            if key not in unique or (not unique[key].get("utilities") and r.get("utilities")):
                unique[key] = r
        source = list(unique.values())
        picks = [(str(r["parent"]), r["picks"]) for r in source]
        rows = [{"case": str(r["parent"]), "action": action, "risk_reduction": value["risk_reduction"]}
                for r in source for action, value in r.get("utilities", {}).items()]
    return picks, rows


def plot_rankings(stage, results, out):
    picks, rows = ranking_data(stage, results)
    criteria = ordered({p for _, values in picks for p in values}, POLICIES)
    actions = ordered({a for _, values in picks for a in values.values()}, list(ACTION_LABELS))
    case_votes = defaultdict(list)
    for case, values in picks:
        for criterion, action in values.items():
            case_votes[case, criterion].append(action)
    votes = np.zeros((len(criteria), len(actions)))
    for (_, criterion), choices in case_votes.items():
        for action in choices:
            votes[criteria.index(criterion), actions.index(action)] += 1 / len(choices)
    cases = list(dict.fromkeys(r["case"] for r in rows))
    risk_actions = list(dict.fromkeys(r["action"] for r in rows))
    values = defaultdict(list)
    for r in rows:
        values[r["case"], r["action"]].append(r["risk_reduction"])
    matrix = np.array([[np.mean(values[c, a]) if values[c, a] else np.nan for a in risk_actions] for c in cases])
    fig, axes = plt.subplots(1, 2, figsize=(13, max(6.3, len(cases) * .23)))
    image = heatmap(axes[0], votes, [label(x, POLICY_LABELS) for x in criteria],
                    [label(x, ACTION_LABELS) for x in actions], fmt=".1f", cmap="Blues", vmin=0)
    axes[0].set_title("(a) Initial action choices")
    fig.colorbar(image, ax=axes[0], shrink=.65, label="Case/parent-equivalent choices")
    limit = max(float(np.nanmax(np.abs(matrix))), 1e-12)
    image = heatmap(axes[1], matrix, [label(x, CASE_LABELS) for x in cases],
                    [label(x, ACTION_LABELS) for x in risk_actions], cmap="RdBu",
                    norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit))
    axes[1].set_title("(b) Expected material-risk reduction")
    fig.colorbar(image, ax=axes[1], shrink=.65, label="Fixed-scale material MAE reduction")
    fig.suptitle(f"{stage}: action ranking from saved numerical predictions", fontsize=14)
    note = ("E0: seven exposed development cases; descriptive only, no population or blind-test claim.\n"
            if stage == "E0" else ("E2 reuses E1 parents: exposed closed-loop followup. Noise repeats averaged within parent; initial ranking only.\n"
                                   if stage == "E2" else "Noise replicates are averaged within parent; initial ranking only. Saved-results snapshot.\n"))
    save(fig, out, f"{stage}_action_risk_matrix",
         note + "* Reference-driven acquisition. Risk values use a finite continuous-proxy-derived reference, not an exact full-Maxwell posterior.")
    return {"ranking_cases_or_parents": len(cases), "criteria": criteria, "risk_actions": risk_actions}


def policy_data(stage, results):
    if stage == "E0":
        first = read(results / "E0_CANDIDATE_ONLY_REALIZED.json")
        second = read(results / "E0_POLICY_REALIZED.json")
        a = {(r["case"], r["policy"]): r for r in first}
        b = {(r["case"], r["policy"]): r for r in second}
        records = []
        for key in sorted(a.keys() & b.keys()):
            r, s = a[key], b[key]
            records.append({"parent": key[0], "noise": 0, "policy": key[1],
                "blocks": s["measurement_blocks"], "mode1": r["candidate_only_selected_MAE"],
                "mode2": s["selected_material_MAE"], "initial1": r["initial_candidate_only_MAE"],
                "initial2": s["initial_selected_MAE"]})
        return records
    data = read(results / f"{stage}_POLICY_RESULTS.json")
    return [{"parent": str(r["parent"]), "noise": r["noise"], "policy": r["policy"],
             "blocks": r["measurement_blocks"], "mode1": r["mode1_candidate_only"]["selected_MAE"],
             "mode2": r["mode2_reference_selection"]["selected_MAE"]} for r in data]


def parent_average(records, key):
    values = defaultdict(list)
    for record in records:
        values[record["parent"]].append(record[key])
    return {parent: float(np.mean(v)) for parent, v in values.items()}


def one_block_panel(records):
    chosen = {}
    for r in records:
        if r["blocks"] > 1:
            continue
        key = (r["parent"], r["noise"], r["policy"])
        if key not in chosen or r["blocks"] >= chosen[key]["blocks"]:
            chosen[key] = r
    return list(chosen.values())


def plot_mae(stage, records, out):
    panel = one_block_panel(records)
    policies = ordered({r["policy"] for r in panel}, POLICIES)
    by_policy = {p: [r for r in panel if r["policy"] == p] for p in policies}
    common = set.intersection(*({r["parent"] for r in v} for v in by_policy.values()))
    cases = sorted(common, key=lambda c: (0, int(c)) if str(c).isdigit() else (1, c))
    fig, axes = plt.subplots(1, 2, figsize=(max(15, len(policies) * .8), max(5.8, len(cases) * .22)))
    matrices = []
    for mode in ("mode1", "mode2"):
        parent_values = {p: parent_average(rows, mode) for p, rows in by_policy.items()}
        matrices.append(np.array([[parent_values[p][c] for p in policies] for c in cases]))
    vmax = max(float(np.max(m)) for m in matrices)
    titles = ["(a) Mode 1: unweighted candidate selection pool", "(b) Mode 2: finite reference decision pool"]
    for ax, matrix, title in zip(axes, matrices, titles):
        image = heatmap(ax, matrix, [label(x, CASE_LABELS) for x in cases],
                        [label(p, POLICY_LABELS) for p in policies], cmap="magma_r", vmin=0, vmax=vmax)
        ax.set_title(title, fontsize=11)
        fig.colorbar(image, ax=ax, shrink=.7, label="Selected material MAE (lower is better)")
    fig.suptitle(f"{stage}: selected material error at at most one new measurement block", fontsize=14)
    exposure = "Seven exposed cases; no blind-population inference." if stage == "E0" else f"{len(cases)} common saved parents; noise replicates averaged within parent."
    if stage == "E2":
        exposure += " Same E1 parents: exposed adaptive followup."
    pool_note = ("Fixed pools: selection alone cannot improve generator Best-of-K."
                 if stage != "E2" else "Finite catalogue refresh is separate from selection; no unlimited material-discovery guarantee.")
    save(fig, out, f"{stage}_selected_MAE_modes",
         exposure + " * Risk/EIG action selection uses the proxy reference in either panel.\n"
         + pool_note + " No general unknown-topology claim.")
    return {"MAE_common_parents": cases, "MAE_policies": policies, "MAE_block_limit": 1}


def plot_cost(stage, records, out):
    policies = ordered({r["policy"] for r in records}, POLICIES)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, len(policies)))
    max_blocks = int(max(r["blocks"] for r in records))
    curves = {}
    for i, policy in enumerate(policies):
        curves[policy] = {}
        rows = [r for r in records if r["policy"] == policy]
        if stage == "E0":
            blocks = sorted({0, *[r["blocks"] for r in rows]})
        else:
            blocks = list(range(max_blocks + 1))
        for ax, mode, initial in zip(axes, ("mode1", "mode2"), ("initial1", "initial2")):
            xs, ys = [], []
            for block in blocks:
                if stage == "E0":
                    subset = [r for r in rows if r["blocks"] == block]
                    xs.append(block)
                    if block == 0 and not subset:
                        ys.append(float(np.mean([r[initial] for r in rows])))
                    else:
                        ys.append(float(np.mean(list(parent_average(subset, mode).values()))))
                else:
                    last = {}
                    for r in rows:
                        if r["blocks"] <= block:
                            key = (r["parent"], r["noise"])
                            if key not in last or r["blocks"] >= last[key]["blocks"]:
                                last[key] = r
                    subset = list(last.values())
                    xs.append(float(np.mean(list(parent_average(subset, "blocks").values()))))
                    ys.append(float(np.mean(list(parent_average(subset, mode).values()))))
            ax.plot(xs, ys, "o-", ms=4, lw=1.5, color=colors[i], label=label(policy, POLICY_LABELS))
            curves[policy][mode] = {"mean_actual_blocks": xs, "mean_parent_MAE": ys}
    for ax, mode in zip(axes, ("Mode 1 candidate pool", "Mode 2 reference pool")):
        ax.set(title=mode, xlabel="Mean additional blocks actually used", ylabel="Mean selected material MAE")
        ax.set_xticks(sorted({0, *[r["blocks"] for r in records]}))
        ax.grid(alpha=.2)
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8, frameon=False)
    fig.suptitle(f"{stage}: material error versus measurement-block cost", fontsize=14)
    save(fig, out, f"{stage}_measurement_cost_curve",
         ("E0 exposed-case mean only; zero new Maxwell evaluations. " if stage == "E0" else
          (("E2 reuses E1 parents (exposed followup). " if stage == "E2" else "") +
           "Noise repeats averaged within parent. Stopped episodes keep their last result; x is mean paid blocks. ")) +
         "Simulation block units, not hardware seconds/energy.\n"
         "Full prediction/setup/calibration/computation costs are separate in RESOURCE_SUMMARY.json; * denotes reference-driven acquisition.")
    return {"cost_curve_maximum_block_cap": max_blocks,
            "cost_curve_stop_carry": "At each block budget, keep last saved record per parent/noise/policy; stopped episodes stay included. Plot actual mean paid blocks, never pretend a stopped episode bought more data.",
            "cost_curves": curves}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["E0", "E1", "E2"], default="E0")
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    parser.add_argument("--out-dir", type=Path, default=ROOT / "figures")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "savefig.facecolor": "white"})
    metadata = plot_rankings(args.stage, args.results, args.out_dir)
    records = policy_data(args.stage, args.results)
    metadata.update(plot_mae(args.stage, records, args.out_dir))
    metadata.update(plot_cost(args.stage, records, args.out_dir))
    metadata.update(stage=args.stage, source_results=str(args.results.resolve()),
                    output_figures=str(args.out_dir.resolve()), new_maxwell_attempts=0,
                    new_observations_simulated=0, new_training_runs=0,
                    script_wall_seconds=time.perf_counter() - _WALL_START,
                    script_cpu_seconds=time.process_time() - _CPU_START,
                    scope="Saved-results visualization only; E0 exposed cases and finite proxy reference limitations are printed on figures")
    (args.out_dir / f"{args.stage}_FIGURE_RECEIPT.json").write_text(json.dumps(metadata, indent=2))
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
