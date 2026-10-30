"""Offline Stage 1 receipt analysis and figures; no solver or new physics.

Reads frozen DEVELOPMENT receipts, already-disclosed R4 material truth, and the
parent's scored summary. Never opens future fresh/private FBG truth.
"""
from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
import numpy as np

HERE = Path(__file__).resolve().parents[1]
R4 = HERE.parent / "research_z_alg_r4_diversity"
RAW = HERE / "results/stage1"
SCORES = HERE / "results/stage1_scores.json"
CASES = HERE / "data_contract/stage1_cases.json"
OUT = HERE / "results/STAGE1_ANALYSIS.json"
FIG = HERE / "figures"
CASE_ORDER = ("dev_low_curved_bar", "dev_moderate_ring", "dev_high_digit8",
              "dev_high_ring", "dev_high_curved_bar")
METHODS_IMAGE = ("native_trf", "reciprocal_trf", "gd_all", "feedback_beta:0.2",
                 "feedback_diag", "feedback_selected")
IMAGE_LABELS = ("Native", "Reciprocal", "GD all", "β = 0.2", "Fitted diag", "Selected")
TRAJECTORY_METHODS = ("gd_all", "feedback_beta:0.2", "feedback_diag", "feedback_selected")
COLORS = {"gd_all": "#3274a1", "feedback_beta:0.2": "#e8871e",
          "feedback_diag": "#777777", "feedback_selected": "#9467bd"}
CASE_LABELS = {"dev_low_curved_bar": "Low curved bar (A=0.4)",
               "dev_moderate_ring": "Moderate ring (A=2.5)",
               "dev_high_digit8": "High digit 8 (A=12)",
               "dev_high_ring": "High ring (A=12)",
               "dev_high_curved_bar": "High curved bar (A=12)"}


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def iou(a: np.ndarray, b: np.ndarray) -> float:
    union = np.count_nonzero(a | b)
    return float(np.count_nonzero(a & b)/union) if union else 1.


def otsu(x: np.ndarray) -> float:
    x = np.asarray(x, float).ravel()
    low, high = float(x.min()), float(x.max())
    if high <= low:
        return high
    hist, edges = np.histogram(x, bins=256, range=(low, high))
    centers = (edges[:-1]+edges[1:])/2
    weights = np.cumsum(hist)
    moments = np.cumsum(hist*centers)
    n, total = len(x), moments[-1]
    den = weights*(n-weights)
    scores = np.divide((total*weights-n*moments)**2, den,
                       out=np.zeros(256), where=den > 0)
    return float(centers[int(np.argmax(scores))])


def material_metrics(m: list, truth: np.ndarray, mask: np.ndarray, A: float) -> dict:
    arr = np.asarray(m, float).ravel()
    if arr.size != 576 or not np.isfinite(arr).all():
        return {"finite_endpoint": False, "NRMSE": None, "IoU_Otsu": None}
    rmse = float(np.sqrt(np.mean((arr-truth)**2)))
    return {"finite_endpoint": True, "NRMSE": rmse/A,
            "IoU_Otsu": iou(arr > otsu(arr), mask)}


def method_key(receipt: dict) -> str:
    return receipt["task"]["method"]


def compact_phase(p: dict, truth: np.ndarray, mask: np.ndarray, A: float) -> dict:
    metrics = material_metrics(p.get("m"), truth, mask, A)
    counts = p.get("counts") or {}
    return {"stage": p.get("stage"), "phase": p.get("phase"),
            "relative_residual": p.get("relative_residual"),
            "ranking_objective": p.get("ranking_objective"),
            **metrics, "LU": counts.get("LU"),
            "source_RHS": counts.get("source_RHS", 0),
            "adjoint_RHS": counts.get("adjoint_RHS", 0),
            "tangent_RHS": counts.get("tangent_RHS", 0),
            "receiver_RHS": counts.get("receiver_RHS", 0),
            "optimizer_wall_s": p.get("optimizer_wall_s"),
            "wrapper_status": p.get("wrapper_status"),
            "wrapper_stop_reason": p.get("wrapper_stop_reason")}


def summarize() -> tuple[dict, dict, dict]:
    cases = {c["case_id"]: c for c in read(CASES)["cases"]}
    score_rows = {(r["case_id"], r["method"]): r for r in read(SCORES)["rows"]}
    raw = {}
    truth = {}
    for case_id, c in cases.items():
        ident = Path(c["public_path"]).stem
        with np.load(R4 / "private_truth" / f"{ident}.npz", allow_pickle=False) as f:
            truth[case_id] = (np.asarray(f["m"], float).ravel(),
                              np.asarray(f["mask"], bool).ravel(), float(c["amplitude"]))
    rows = []
    for path in sorted(RAW.glob("*.json")):
        if path.name == "SUMMARY.json":
            continue
        receipt = read(path)
        task, result = receipt["task"], receipt.get("result") or {}
        case_id, method = task["case_id"], method_key(receipt)
        raw[(case_id, method)] = receipt
        scored = score_rows.get((case_id, method), {})
        t, mask, A = truth[case_id]
        counts = result.get("counts") or {}
        full_rhs = sum(int(counts.get(k, 0)) for k in
                       ("source_RHS", "adjoint_RHS", "tangent_RHS", "receiver_RHS"))
        phase_rows = [compact_phase(p, t, mask, A) for p in result.get("phases", [])]
        feedback = []
        for f in result.get("feedback_records", []):
            feedback.append({"stage": f.get("stage"), "kind": f.get("kind"),
                             "tangent_energy_ratio": f.get("tangent_energy_ratio"),
                             "fit_wall_s": f.get("fit_wall_s"),
                             "directions": f.get("directions"),
                             "safety": f.get("safety"), "whole_box": f.get("whole_box")})
        rows.append({"case_id": case_id, "method": method,
                     "receipt_path": str(path.relative_to(HERE)),
                     "transport_status": receipt.get("status"),
                     "wrapper_status": result.get("wrapper_status"),
                     "stop_reason": result.get("wrapper_stop_reason"),
                     "wall_seconds": receipt.get("wall_s"),
                     "task_wall_s": result.get("task_wall_s"),
                     "LU": counts.get("LU"), "actual_complex_RHS": full_rhs,
                     "source_RHS": counts.get("source_RHS", 0),
                     "adjoint_RHS": counts.get("adjoint_RHS", 0),
                     "tangent_RHS": counts.get("tangent_RHS", 0),
                     "receiver_RHS": counts.get("receiver_RHS", 0),
                     "feedback_fits": counts.get("feedback_fits", 0),
                     "Gram": counts.get("Gram", 0), "FullJ": counts.get("FullJ", 0),
                     "NRMSE": scored.get("NRMSE"), "IoU_Otsu": scored.get("IoU_Otsu"),
                     "relative_residual": result.get("relative_residual"),
                     "data_pass": scored.get("data_pass"),
                     "joint_pass": scored.get("joint_pass"),
                     "phases": phase_rows, "feedback_records": feedback})
    if len(rows) != 64 or len(score_rows) != 64:
        raise RuntimeError(f"Expected 64 raw and scored receipts, got {len(rows)} / {len(score_rows)}")
    group = defaultdict(list)
    for r in rows:
        group[r["method"]].append(r)
    by_method = [{"method": method, "cases_planned": len(group[method]),
                  "joint_pass": sum(r["joint_pass"] is True for r in group[method]),
                  "joint_fail": sum(r["joint_pass"] is False for r in group[method]),
                  "wall_seconds_sum_observed": sum(r["wall_seconds"] or 0 for r in group[method]),
                  "LU_sum_observed": sum(r["LU"] or 0 for r in group[method]),
                  "actual_complex_RHS_sum_observed": sum(r["actual_complex_RHS"] or 0 for r in group[method])}
                 for method in sorted(group)]
    indexed = {(r["case_id"], r["method"]): r for r in rows}
    escape = []
    for case_id in ("dev_high_digit8", "dev_high_curved_bar", "dev_high_ring"):
        beta = indexed[(case_id, "feedback_beta:0.2")]
        gd = indexed[(case_id, "gd_all")]
        diag = indexed[(case_id, "feedback_diag")]
        stage_pairs = []
        for stage in (8, 12, 24):
            for phase in ("proposal", "corrector"):
                def find(row):
                    return next((p for p in row["phases"] if p["stage"] == stage and p["phase"] == phase), None)
                b, g, d = find(beta), find(gd), find(diag)
                stage_pairs.append({"stage": stage, "phase": phase,
                                    "beta_relative_residual": b["relative_residual"] if b else None,
                                    "gd_all_relative_residual": g["relative_residual"] if g else None,
                                    "fitted_diag_relative_residual": d["relative_residual"] if d else None,
                                    "beta_NRMSE_OFFLINE": b["NRMSE"] if b else None,
                                    "gd_all_NRMSE_OFFLINE": g["NRMSE"] if g else None})
        escape.append({"case_id": case_id, "stage_comparison": stage_pairs,
                       "beta_joint_pass": beta["joint_pass"],
                       "gd_all_joint_pass": gd["joint_pass"],
                       "fitted_diag_joint_pass": diag["joint_pass"]})
    result = {"schema": "fbg-r2-stage1-analysis-v1",
              "status": "DEVELOPMENT_DISCLOSED_R4",
              "raw_receipts": 64, "scored_receipts": 64,
              "source_receipts": "results/stage1/*.json",
              "source_archive_frozen": "audit/STAGE1_SOURCE.zip",
              "source_version_scope": "Analyze archived Stage 1 source/receipts only; later local code-v2 budget fixes are not backported into these runs.",
              "score_source": "results/stage1_scores.json",
              "truth_source": "research_z_alg_r4_diversity/private_truth; previously disclosed DEVELOPMENT only",
              "no_new_physics": True,
              "cost_scope": "45 s wall-only task cap; Tfit/feedback counted in task cost; no matched LU/RHS cap. Overlapping task wall is not summed as elapsed batch wall. RHS is source+adjoint+tangent+receiver full-system RHS; Gram excluded.",
              "rows": rows, "by_method": by_method, "high_case_phase_comparison": escape}
    return result, raw, truth


def reconstruction_figure(raw: dict, truth: dict, scores: dict) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(5, 7, figsize=(24, 16), constrained_layout=True)
    norm = PowerNorm(gamma=.46, vmin=.02, vmax=20.)
    last = None
    for i, case_id in enumerate(CASE_ORDER):
        t, _, A = truth[case_id]
        last = axes[i, 0].imshow(t.reshape(24, 24), cmap="magma", norm=norm,
                                 origin="upper", interpolation="nearest")
        axes[i, 0].set_title(f"Truth  ·  A={A:g}", fontsize=9)
        axes[i, 0].set_ylabel(CASE_LABELS[case_id], fontsize=10)
        for j, (method, label) in enumerate(zip(METHODS_IMAGE, IMAGE_LABELS), start=1):
            ax = axes[i, j]
            receipt = raw.get((case_id, method))
            result = receipt.get("result") if receipt else None
            scored = scores.get((case_id, method))
            m = result.get("m") if result else None
            if m is not None and np.asarray(m).size == 576 and np.isfinite(np.asarray(m, float)).all():
                ax.imshow(np.asarray(m, float).reshape(24, 24), cmap="magma", norm=norm,
                          origin="upper", interpolation="nearest")
                n = scored.get("NRMSE") if scored else None
                residual = result.get("relative_residual")
                passed = scored.get("joint_pass") if scored else False
                ax.set_title(f"{label}  ·  {'PASS' if passed else 'FAIL'}\n"
                             f"N={n:.3f}  R={residual:.1e}", fontsize=8.5,
                             color="#157a45" if passed else "#8b2727")
            else:
                ax.set_facecolor("#ededed")
                ax.text(.5, .5, "NO FINITE\nENDPOINT", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"{label}  ·  NOT_RUN/FAILED", fontsize=8.5)
        for ax in axes[i]:
            ax.set_xticks([]); ax.set_yticks([])
    cbar = fig.colorbar(last, ax=axes, location="right", shrink=.72, pad=.01)
    cbar.set_label("Material m · shared [0.02, 20], power γ=0.46")
    fig.suptitle("Stage 1 DEVELOPMENT reconstructions · full nonlinear 24×24 material",
                 fontsize=16, fontweight="bold")
    fig.savefig(FIG / "stage1_reconstruction.png", dpi=160)
    plt.close(fig)


def trajectories_figure(raw: dict, analysis: dict) -> None:
    indexed = {(r["case_id"], r["method"]): r for r in analysis["rows"]}
    fig, axes = plt.subplots(5, 2, figsize=(16, 18), constrained_layout=True)
    phase_x = np.arange(6)
    phase_labels = ["8P", "8C", "12P", "12C", "24P", "24C"]
    method_abbrev = {"native_trf": "N", "reciprocal_trf": "R", "data_lbfgs": "LB",
                     "gd_all": "GD", "gd_late": "GL", "fixedgram": "FG",
                     "feedback_beta:0.2": "β", "feedback_diag": "D", "feedback_block": "B",
                     "feedback_selected": "S", "feedback_self": "Self", "feedback_shuffle": "Sh",
                     "feedback_beta:0.05": "β05", "feedback_beta:1.0": "β1"}
    for i, case_id in enumerate(CASE_ORDER):
        ax, cost = axes[i]
        for method in TRAJECTORY_METHODS:
            row = indexed.get((case_id, method))
            if not row:
                continue
            phase_map = {(p["stage"], p["phase"]): p for p in row["phases"]}
            vals = [phase_map.get((s, p), {}).get("relative_residual", np.nan)
                    for s in (8, 12, 24) for p in ("proposal", "corrector")]
            ax.plot(phase_x, vals, marker="o", ms=4, lw=1.8,
                    label=method_abbrev[method], color=COLORS[method])
        ax.axhline(.001, color="#2c8c61", ls="--", lw=1, alpha=.7)
        ax.set_yscale("log"); ax.set_ylim(1e-5, 2)
        ax.set_xticks(phase_x, phase_labels)
        ax.grid(alpha=.2, which="both")
        ax.set_title(CASE_LABELS[case_id] + " · phase-end physical residual", fontsize=11)
        ax.set_ylabel("Relative data residual")
        if i == 0:
            ax.legend(ncol=4, fontsize=8, loc="upper right")
        case_rows = [r for r in analysis["rows"] if r["case_id"] == case_id]
        for row in case_rows:
            method = row["method"]
            wall, lu = row["wall_seconds"], row["LU"]
            color = COLORS.get(method, "#9da8b2")
            marker = "o" if row["joint_pass"] else "x"
            cost.scatter(wall, lu, marker=marker, color=color, s=52, alpha=.9)
            cost.annotate(method_abbrev.get(method, method), (wall, lu),
                          xytext=(3, 2), textcoords="offset points", fontsize=7)
        cost.axvline(45, color="#b14b48", ls="--", lw=1)
        cost.set_xlim(0, 49)
        cost.set_ylim(0, max(1500, max((r["LU"] or 0 for r in case_rows), default=0)*1.1))
        cost.set_title(CASE_LABELS[case_id] + " · measured cost (○ pass, × fail)", fontsize=11)
        cost.set_xlabel("Observed task wall (s)")
        cost.set_ylabel("Full VIE LU")
        cost.grid(alpha=.2)
    fig.suptitle("Stage 1 DEVELOPMENT · finite paths and measured task costs\n"
                 "P=proposal, C=corrector; wall-only 45 s cap includes feedback fit",
                 fontsize=15, fontweight="bold")
    fig.savefig(FIG / "stage1_trajectories.png", dpi=160)
    plt.close(fig)


def main() -> None:
    analysis, raw, truth = summarize()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(analysis, indent=2, allow_nan=False) + "\n")
    scores = {(r["case_id"], r["method"]): r for r in read(SCORES)["rows"]}
    reconstruction_figure(raw, truth, scores)
    trajectories_figure(raw, analysis)
    print(json.dumps({"rows": len(analysis["rows"]),
                      "joint_pass": sum(r["joint_pass"] is True for r in analysis["rows"]),
                      "figures": [str(FIG / "stage1_reconstruction.png"),
                                  str(FIG / "stage1_trajectories.png")],
                      "analysis": str(OUT)}))


if __name__ == "__main__":
    main()
