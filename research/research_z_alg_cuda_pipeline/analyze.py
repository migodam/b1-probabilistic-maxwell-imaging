"""Summarize actual execution receipts and draw performance/queue figures."""
from collections import Counter
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent
RESULTS, FIGURES = ROOT / "results", ROOT / "figures"


def main():
    FIGURES.mkdir(exist_ok=True)
    throughput = json.loads((RESULTS / "throughput_v1.json").read_text())
    backend = json.loads((RESULTS / "backend_choice.json").read_text())
    svd = json.loads((RESULTS / "svd_choice.json").read_text())
    ledger, runs = [], []
    for folder in sorted(RESULTS.glob("endtoend_*")):
        summary = json.loads((folder / "runner_summary.json").read_text())
        total = Counter()
        receipts = []
        for path in sorted(folder.glob("public_*.json")):
            receipt = json.loads(path.read_text())
            z = receipt.get("result", {})
            costs = z.get("counts", {})
            total.update(costs)
            row = {"run": folder.name, "lanes": summary["lanes"], "task_id": receipt["task_id"],
                   "status": receipt["status"], "method": receipt["task"]["method"],
                   "wall_s": receipt["worker_wall_s"], "optimizer_evaluations": receipt["optimizer_evaluations"],
                   "LU": costs.get("LU", 0), "source_RHS": costs.get("source_RHS", 0),
                   "adjoint_RHS": costs.get("adjoint_RHS", 0), "receiver_RHS": costs.get("receiver_RHS", 0),
                   "FullJ": costs.get("FullJ", 0), "Gram": costs.get("Gram", 0),
                   "host_returns": costs.get("host_result_returns", 0),
                   "relative_residual": z.get("relative_residual"),
                   "state_relative_residual": z.get("state_relative_residual"),
                   "truth_access": "NONE_PUBLIC_DATA_RUNTIME_ONLY"}
            ledger.append(row)
            receipts.append(row)
        handoffs = []
        ends = {r["task_id"]: r["finished_seconds"] for r in summary["receipts"]}
        for slot in range(summary["lanes"]):
            starts = sorted((r for r in summary["starts"] if r["slot"] == slot), key=lambda r: r["seconds"])
            handoffs += [current["seconds"] - ends[previous["task_id"]] for previous, current in zip(starts, starts[1:])]
        runs.append({"run": folder.name, "lanes": summary["lanes"], "wall_s": summary["wall_s"],
                     "shared_setup_s": summary["shared_init_wall_s"], "status_counts": summary["status_counts"],
                     "optimizer_evaluations": summary["cumulative_optimizer_evaluations"],
                     "counts": dict(total), "shared_setup_counts": summary["shared_setup_counts"],
                     "handoff_median_ms": float(np.median(handoffs) * 1000) if handoffs else None,
                     "handoff_max_ms": max(handoffs) * 1000 if handoffs else None,
                     "svd_backend": summary.get("svd_backend", "cpu"), "tasks": receipts})
    with (ROOT / "COST_LEDGER.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(ledger[0]))
        writer.writeheader(); writer.writerows(ledger)
    first = next(r for r in runs if r["run"] == "endtoend_1")
    for r in runs:
        r["speedup_vs_one_lane"] = first["wall_s"] / r["wall_s"]
        r["wall_saving_vs_one_lane"] = 1 - r["wall_s"] / first["wall_s"]
    output = {"scope": "ENGINEERING_RUNTIME_ONLY; already-public data; no new blind GT recovery experiment",
              "hardware": throughput["hardware"], "dtype": "complex128", "physical_dimension": 576,
              "sources": 12, "receivers": 48, "throughput": throughput, "backend_selection": backend,
              "svd_selection": svd, "endtoend_runs": runs,
              "occupancy": "NOT_MEASURED; installed torch profiler supports CPU activity only; nsys absent",
              "scientific_algorithm_claim": "NONE_NEW; physics/adjoints/GlobalDual/TRF/L-BFGS-B inherited",
              "synchronization_claim": "No explicit device-wide synchronization in the new hot path; necessary per-result lane events remain",
              "timing_scope": "Executor run wall includes shared GPU setup, observations, all tasks and writer; Python/torch import and SSH transport are outside this clock"}
    (ROOT / "RESULTS_SUMMARY.json").write_text(json.dumps(output, indent=2))
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.4), layout="constrained")
    for axis, obj, title in zip(axs[:2], ["data", "dual"], ["Data objective + adjoint", "GlobalDual + adjoint"]):
        rows = [r for r in throughput["rows"] if r.get("objective") == obj and
                ((r["kind"] == "legacy_processes" and r["concurrency"] == 1) or r["kind"] == "resident_streams")]
        labels = ["Old 1 worker" if r["kind"] == "legacy_processes" else f'{r["concurrency"]} stream(s)' for r in rows]
        bars = axis.bar(labels, [r["evaluations_per_second"] for r in rows], color=["#888888", "#6a9bbf", "#3d7caa", "#125886"])
        axis.bar_label(bars, fmt="%.1f", padding=3)
        axis.set(title=title, ylabel="Evaluations / second")
        axis.tick_params(axis="x", rotation=20)
    rows = [r for r in throughput["rows"] if r["kind"] == "batched_candidate_score"]
    bars = axs[2].bar([str(r["concurrency"]) for r in rows], [r["evaluations_per_second"] for r in rows], color="#44836b")
    axs[2].bar_label(bars, fmt="%.0f", padding=3)
    axs[2].set(title="Full nonlinear candidate scoring", xlabel="Batch size", ylabel="Candidates / second")
    fig.suptitle("RTX 4060: fixed-work Maxwell throughput, complex128")
    fig.text(.01, -.04, "Shared CUDA setup reported separately. Forward-only scoring and objective+gradient are different workloads. Single timing run.", fontsize=9)
    fig.savefig(FIGURES / "cuda_throughput.png", dpi=180, bbox_inches="tight"); plt.close(fig)
    fig, axs = plt.subplots(2, 1, figsize=(12, 7), gridspec_kw={"height_ratios": [1, 1.5]}, layout="constrained")
    labels = [f'{r["lanes"]} lane(s)' + (" + hybrid SVD" if "hybrid" in r["run"] else "") for r in runs]
    bars = axs[0].barh(labels, [r["wall_s"] for r in runs], color=["#83919d", "#659cbd", "#317ca8", "#3d866c"][:len(runs)])
    axs[0].bar_label(bars, fmt="%.2f s", padding=3)
    axs[0].set(xlabel="Shared setup + 12 complete optimizer tasks + output (seconds)", title="Equal task allocations; no GT or accuracy ranking")
    r = next(r for r in runs if r["run"] == "endtoend_4")
    summary = json.loads((RESULTS / r["run"] / "runner_summary.json").read_text())
    colors = {"native_trf": "#4c78a8", "reciprocal_trf": "#f58518", "data_lbfgs": "#54a24b", "globaldual": "#b279a2"}
    starts = {x["task_id"]: x for x in summary["starts"]}
    for row in summary["receipts"]:
        start = starts[row["task_id"]]
        method = next(m for m in colors if row["task_id"].endswith(m))
        axs[1].broken_barh([(start["seconds"], row["finished_seconds"] - start["seconds"])], (row["slot"] - .34, .68), facecolors=colors[method])
    axs[1].set(yticks=range(4), yticklabels=["Lane 0", "Lane 1", "Lane 2", "Lane 3"], xlabel="Executor elapsed seconds", title="Continuous task handoff: CPU + CUDA intervals, not an SM occupancy trace")
    axs[1].legend(handles=[Patch(color=c, label=m) for m, c in colors.items()], ncols=4, fontsize=8, loc="upper center")
    fig.savefig(FIGURES / "end_to_end_queue.png", dpi=180); plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 4.5), layout="constrained")
    for name, marker in [("scipy_gesdd_cpu_1thread", "o"), ("cusolver_gesvd", "s"), ("cusolver_gesvdj", "^")]:
        rows = [r for r in svd["rows"] if r["backend"] == name]
        ax.plot([r["stage"] ** 2 for r in rows], [r["wall_s"] for r in rows], marker=marker, label=name)
    ax.set(xlabel="Material parameter dimension", ylabel="SVD + required host returns (seconds)", title="Existing mature operators: CUDA is worthwhile only for the large SVD")
    ax.legend(fontsize=8)
    fig.savefig(FIGURES / "svd_operator_selection.png", dpi=180); plt.close(fig)
    print(json.dumps({"runs": [{k: r[k] for k in ["run", "wall_s", "optimizer_evaluations", "speedup_vs_one_lane", "handoff_median_ms", "handoff_max_ms"]} for r in runs]}, indent=2))


if __name__ == "__main__":
    main()
