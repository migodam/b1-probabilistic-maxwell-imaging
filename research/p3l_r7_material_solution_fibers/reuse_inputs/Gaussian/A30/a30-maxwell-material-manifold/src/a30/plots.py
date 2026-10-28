"""Scientific figures from saved A30 E0/E1 artifacts only.

No forward model, training, gate decision, label generator or locked loader is
imported. Invalid/missing capacity fits remain visible and are not averaged into
complete campaign statistics. All figures export their plotted raw values.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
import csv
import json
import math
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm


FAMILIES = ("shell", "asymmetric_piecewise", "smooth", "mixed_scale")
COMPLEX = frozenset(("shell", "asymmetric_piecewise", "mixed_scale"))
PARENTS = tuple(range(250065, 250081))
HULL_METHODS = ("full_hard", "small_hard", "full_hull", "small_hull", "full_L_mean", "small_L_mean")
CAPACITY_METHODS = ("atlas", "level_set", "sparse", "single_chart")
BUDGETS = (16, 32, 64)
FIXED_CASES = tuple(f"p{p}_weak_a1_snr20" for p in (250065, 250069, 250073, 250077))
LABELS = {"N3": "N3", "full_hard": "Full hard", "small_hard": "Small hard",
          "full_hull": "Full convex", "small_hull": "Small convex",
          "full_L_mean": "Full L mean", "small_L_mean": "Small L mean",
          "atlas": "Topology atlas", "level_set": "Level set",
          "sparse": "Sparse multiscale", "single_chart": "Single chart"}
FAMILY_LABELS = {"shell": "Shell", "asymmetric_piecewise": "Asymmetric",
                 "smooth": "Smooth", "mixed_scale": "Mixed scale", "complex": "Complex (12 parents)"}
COLORS = {"full_hard": "#475569", "small_hard": "#94a3b8", "full_hull": "#087f8c",
          "small_hull": "#55b2ba", "full_L_mean": "#b45309", "small_L_mean": "#e4a453",
          "atlas": "#087f8c", "level_set": "#7c3aed", "sparse": "#c2410c", "single_chart": "#2563eb"}


def _safe(root: Path, relative: str | Path, *, exists=False) -> Path:
    relative = Path(relative)
    if relative.is_absolute():
        path = relative
    else:
        path = root / relative
    if any(p.lower() == "locked" or p.lower().startswith(("locked_id", "locked_ood")) for p in path.parts):
        raise ValueError("LOCKED_PATH_REJECTED")
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError("PATH_ESCAPES_A30_ROOT")
    if exists and not path.is_file():
        raise FileNotFoundError(path)
    return path


def _read_json(root, relative):
    return json.loads(_safe(root, relative, exists=True).read_text(encoding="utf-8"))


def _read_csv(root, relative, *, required=True):
    path = _safe(root, relative)
    if not path.is_file():
        if required:
            raise FileNotFoundError(path)
        return []
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _number(value):
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return value if math.isfinite(value) else None


def _valid(row):
    flag = str(row.get("valid", "false")).strip().lower()
    status = str(row.get("status", "")).strip().upper()
    return flag in {"true", "1", "yes"} and status not in {
        "FAILED", "INVALID", "NOT_RUN", "BUDGET_STOP", "MISSING", "FAIL"}


def _family(parent):
    if int(parent) not in PARENTS:
        raise ValueError("PARENT_OUTSIDE_EXPOSED_WHITELIST")
    return FAMILIES[(int(parent) - 250065) // 4]


def _group_parents(group):
    return [p for p in PARENTS if group == "complex" and _family(p) in COMPLEX or _family(p) == group]


def _ci(values):
    values = np.asarray(values, dtype=float)
    if not len(values) or not np.isfinite(values).all():
        return None, None, None
    rng = np.random.default_rng(291011)
    draws = values[rng.integers(0, len(values), (2000, len(values)))].mean(axis=1)
    lo, hi = np.quantile(draws, [.025, .975])
    return float(values.mean()), float(lo), float(hi)


def _write_csv(path, rows):
    keys = list(dict.fromkeys(key for row in rows for key in row))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _note_panel(ax, title, message):
    ax.set_axis_off()
    ax.text(.5, .66, title, ha="center", va="center", fontsize=18, fontweight="bold", transform=ax.transAxes)
    ax.text(.5, .40, message, ha="center", va="center", fontsize=11, color="#475569", transform=ax.transAxes)


def _style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11,
                         "axes.labelsize": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
                         "svg.fonttype": "none", "axes.titleweight": "semibold", "xtick.labelsize": 9,
                         "ytick.labelsize": 9, "legend.fontsize": 9})


@dataclass
class FigureWriter:
    root: Path
    output: Path
    metadata_output: Path
    mapping: list[dict] = field(default_factory=list)

    def save(self, stem, figure, rows, *, caption_zh, caption_en, sources, status="COMPLETE", notes=None):
        self.output.mkdir(parents=True, exist_ok=True)
        raw = self.output / "raw_plot_data" / (stem + ".csv")
        figure.savefig(self.output / (stem + ".png"), dpi=200, bbox_inches="tight")
        figure.savefig(self.output / (stem + ".svg"), bbox_inches="tight")
        plt.close(figure)
        _write_csv(raw, rows or [{"status": status}])
        self.mapping.append({"figure": stem, "status": status,
                             "png": str((self.output / (stem + ".png")).relative_to(self.root)),
                             "svg": str((self.output / (stem + ".svg")).relative_to(self.root)),
                             "raw_plot_csv": str(raw.relative_to(self.root)),
                             "sources": sources, "caption_zh": caption_zh, "caption_en": caption_en,
                             "notes": notes or [], "Maxwell_calls": 0, "training_updates": 0, "locked_reads": 0})

    def finish(self, timing):
        self.metadata_output.mkdir(parents=True, exist_ok=True)
        (self.metadata_output / "FIGURE_SOURCE_MAPPING.json").write_text(
            json.dumps({"figures": self.mapping, "timing": timing,
                        "timing_is_nested_in_job_wrapper": True, "SHA256_checks": "NOT_RUN_BY_USER_POLICY"},
                       indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _load_hull(root, *, details=False):
    base = "results/a30/hull/"
    summary = _read_json(root, base + "HULL_AUDIT_SUMMARY.json")
    rows = _read_csv(root, base + "HULL_MATERIAL_AND_STRUCTURE.csv")
    cases = _read_csv(root, base + "A29_1_CASES.csv") if details else []
    weights = _read_json(root, base + "A29_1_WEIGHTS.json") if details else []
    with np.load(_safe(root, base + "HULL_RECONSTRUCTIONS.npz", exists=True), allow_pickle=False) as archive:
        images = {key: archive[key].copy() for key in archive.files}
    expected = {f"p{p}_{b}_a1_snr20" for p in PARENTS for b in ("weak", "ring")}
    if set(map(str, images["case_id"])) != expected or len(images["case_id"]) != 32:
        raise ValueError("COMPLETE_32_HULL_RECONSTRUCTIONS_REQUIRED")
    for name in ("truth", "background", "N3", *HULL_METHODS):
        if images[name].shape != (32, 576) or not np.isfinite(images[name]).all():
            raise ValueError("HULL_IMAGE_LAYOUT_OR_FINITE_CHECK_FAILED:" + name)
    if len(rows) != 32 * 9 or any(name not in summary for name in HULL_METHODS):
        raise ValueError("COMPLETE_HULL_METRICS_REQUIRED")
    return summary, rows, cases, weights, images


def plot_hull_family(writer, summary):
    groups = ("complex", *FAMILIES)
    fig, axes = plt.subplots(1, 5, figsize=(15.7, 4.8), sharey=True)
    rows = []
    extrema = [0.]
    for ax, group in zip(axes, groups):
        for i, method in enumerate(HULL_METHODS):
            s = summary[method]["complex"] if group == "complex" else summary[method]["family"][group]
            point, lo, hi = 100 * s["point"], 100 * s["ci95"][0], 100 * s["ci95"][1]
            ax.errorbar(point, i, xerr=[[point - lo], [hi - point]], fmt="o", color=COLORS[method], capsize=3, ms=5)
            extrema.extend((lo, hi))
            rows.append({"group": group, "method": method, "parent_count": s["parent_count"],
                         "gain_fraction": s["point"], "gain_percent": point,
                         "ci95_low_percent": lo, "ci95_high_percent": hi,
                         "bootstrap_samples": s["bootstrap_samples"], "bootstrap_seed": s["seed"],
                         "statistical_unit": s["statistical_unit"]})
        ax.axvline(0, color="#94a3b8", lw=.9)
        ax.grid(axis="x", color="#e2e8f0", lw=.7)
        ax.set_title(FAMILY_LABELS[group])
        ax.set_xlabel("Parent gain vs N3 (%)")
    padding = max(3., .08 * (max(extrema) - min(extrema)))
    for ax in axes:
        ax.set_xlim(min(extrema) - padding, max(extrema) + padding)
        ax.set_ylim(len(HULL_METHODS) - .5, -.5)
    axes[0].set_yticks(range(len(HULL_METHODS)), [LABELS[x] for x in HULL_METHODS])
    fig.suptitle("Saved material-bank capacity and L-weighted estimators", y=1.02, fontsize=14)
    fig.tight_layout()
    writer.save("fig01_hull_family_gains", fig, rows,
                caption_zh="全/小 hard、convex hull 与 L mean 相对 N3 的 parent 配对材料改善。先在每个 parent 内平均 weak/ring 两背景误差，再计算相对改善；点与区间直接引用保存的 2000 次 parent bootstrap。complex 为 12 parents，smooth 单列。所有数据来自已暴露 development；convex 是 offline 容量诊断，L mean 是有限且 data-dependent bank 的 weighted estimator。",
                caption_en="Parent-paired material gain against N3, after averaging the two backgrounds within each exposed parent. Points and exploratory 95% intervals are copied from the saved 2,000-parent-bootstrap summary. Convex results are offline capacity diagnostics; L means are data-dependent-bank weighted estimators, not calibrated posteriors.",
                sources=["results/a30/hull/HULL_AUDIT_SUMMARY.json"],
                notes=["All five panels use one common horizontal scale.", "No gate is assigned by the figure generator."])


def plot_weights_ess(writer, cases, weights):
    by_case = {row["case_id"]: row for row in cases}
    records = []
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.7))
    ess = defaultdict(dict)
    distributions = defaultdict(list)
    for row in weights:
        case = row["case_id"]
        meta = by_case[case]
        pid = int(meta["parent_id"])
        for pool in ("full", "small"):
            value = np.asarray(row[pool + "_oracle"], dtype=float)
            if not np.isfinite(value).all() or value.min() < -1e-12 or abs(float(value.sum()) - 1) > 1e-9:
                raise ValueError("INVALID_SAVED_HULL_WEIGHTS:" + case)
            ordered = np.sort(value)[::-1]
            distributions[pool].append(ordered)
            oracle_ess = float(1 / np.sum(value * value))
            linear_ess = float(meta[pool + "_L_ESS"])
            ess[(pid, meta["background_id"])][pool + "_hull"] = oracle_ess
            ess[(pid, meta["background_id"])][pool + "_L_mean"] = linear_ess
            order = np.argsort(-value, kind="stable")
            ranks = np.empty(len(value), dtype=int)
            ranks[order] = np.arange(1, len(value) + 1)
            actual_indices = list(range(len(value))) if pool == "full" else row["small_indices"]
            for j, weight in enumerate(value):
                records.append({"record_type": "saved_oracle_weight", "case_id": case, "parent_id": pid,
                                "family": meta["family"], "background_id": meta["background_id"], "pool": pool,
                                "candidate_index_in_full_pool": int(actual_indices[j]), "weight_rank": int(ranks[j]),
                                "weight": float(weight), "oracle_ESS": oracle_ess, "L_mean_ESS": linear_ess})
    for pool, method in (("full", "full_hull"), ("small", "small_hull")):
        maxn = max(map(len, distributions[pool]))
        mean = np.mean([np.pad(v, (0, maxn - len(v))) for v in distributions[pool]], axis=0)
        axes[0].step(np.arange(1, maxn + 1), np.cumsum(mean), where="mid", color=COLORS[method], label=LABELS[method], lw=2)
        for rank, cumulative in enumerate(np.cumsum(mean), 1):
            records.append({"record_type": "mean_rank_cumulative_weight", "pool": pool,
                            "weight_rank": rank, "cumulative_weight": float(cumulative), "case_count": 32})
    axes[0].set(xlabel="Descending weight rank", ylabel="Mean cumulative oracle weight", ylim=(0, 1.03),
                title="Offline convex-weight concentration")
    axes[0].legend(frameon=False)
    for method in ("full_hull", "small_hull", "full_L_mean", "small_L_mean"):
        values = []
        for pid in PARENTS:
            value = float(np.mean([ess[(pid, bg)][method] for bg in ("weak", "ring")]))
            values.append(value)
            records.append({"record_type": "parent_ESS", "parent_id": pid, "family": _family(pid),
                            "method": method, "ESS": value, "conditions": 2})
        axes[1].plot(range(16), values, "o-", color=COLORS[method], ms=3.5, lw=1.2, label=LABELS[method])
    axes[1].set_xticks(range(16), [str(p)[-2:] for p in PARENTS])
    axes[1].set(xlabel="Parent ID suffix (2500xx)", ylabel="ESS, averaged over weak/ring", title="Weight ESS on the saved bank")
    axes[1].legend(frameon=False, ncol=2, fontsize=8)
    for ax in axes:
        ax.grid(color="#e2e8f0", lw=.7)
    fig.tight_layout()
    writer.save("fig02_weights_and_ess", fig, records,
                caption_zh="左图为保存的 offline simplex 权重按降序排列后的平均累积质量；不足最大 bank 长度的尾部权重为零。右图为每个 parent 两背景平均的 ESS=1/Σw²：convex ESS 由原保存权重计算，L mean ESS 直接引用原 case 表。ESS 只描述该 bank 权重集中程度，不证明 posterior 校准。",
                caption_en="Saved offline simplex weights are sorted to show cumulative concentration. ESS equals 1/sum(w^2): oracle ESS is computed from the saved hull weights, while L-estimator ESS is copied from the saved case table. This is a descriptive weight statistic, not a posterior-calibration result.",
                sources=["results/a30/hull/A29_1_WEIGHTS.json", "results/a30/hull/A29_1_CASES.csv"])


def plot_structure_blur(writer, hull_rows, summary):
    baseline = {row["case_id"]: row for row in hull_rows if row["method"] == "N3"}
    metrics = (("edge_f1", "edge", "Edge F1 difference", True),
               ("support_iou", "iou", "Support IoU difference", True),
               ("component_abs_error_25", "component_MAE_increment", "Component MAE increment", False),
               ("hole_abs_error_25", "hole_MAE_increment", "Hole MAE increment", False),
               ("intermediate_pixel_fraction", None, "Intermediate-pixel fraction difference", None),
               ("magnitude_TV", None, "Magnitude TV difference", None))
    fig, axes = plt.subplots(2, 3, figsize=(13.7, 8.3))
    records = []
    for ax, (metric, summary_key, title, direction) in zip(axes.ravel(), metrics):
        for i, method in enumerate(HULL_METHODS):
            rows = [row for row in hull_rows if row["method"] == method and row["family"] in COMPLEX]
            parent_differences = []
            for pid in _group_parents("complex"):
                pair = [row for row in rows if int(row["parent_id"]) == pid]
                if len(pair) != 2 or {r["background_id"] for r in pair} != {"weak", "ring"}:
                    raise ValueError("STRUCTURE_PARENT_PAIR_INCOMPLETE")
                candidate = float(np.mean([float(r[metric]) for r in pair]))
                reference = float(np.mean([float(baseline[r["case_id"]][metric]) for r in pair]))
                difference = candidate - reference
                parent_differences.append(difference)
                row = {"record_type": "parent_difference", "metric": metric, "method": method,
                       "parent_id": pid, "family": _family(pid), "candidate_mean": candidate,
                       "N3_mean": reference, "difference": difference, "conditions": 2}
                if metric in {"component_abs_error_25", "hole_abs_error_25"}:
                    prefix = "components" if metric.startswith("component") else "holes"
                    row.update(predicted_count_mean=float(np.mean([float(r[prefix + "_25"]) for r in pair])),
                               truth_count_mean=float(np.mean([float(r["truth_" + prefix + "_25"]) for r in pair])))
                records.append(row)
            if summary_key:
                saved = summary[method][summary_key]
                mean, lo, hi = saved["point"], *saved["ci95"]
            else:
                mean, lo, hi = _ci(parent_differences)
            ax.errorbar(mean, i, xerr=[[max(0., mean - lo)], [max(0., hi - mean)]], fmt="o",
                        color=COLORS[method], capsize=3, ms=4.5)
            records.append({"record_type": "complex_summary", "metric": metric, "method": method,
                            "mean_difference": mean, "ci95_low": lo, "ci95_high": hi,
                            "parent_count": 12, "bootstrap_samples": 2000, "bootstrap_seed": 291011,
                            "positive_difference_meaning": "higher_score" if direction is True else
                            "higher_count_error" if direction is False else "descriptive_only"})
        ax.axvline(0, color="#94a3b8", lw=.9)
        ax.grid(axis="x", color="#e2e8f0", lw=.7)
        ax.set_yticks(range(6), [LABELS[x] for x in HULL_METHODS])
        ax.set_ylim(5.5, -.5)
        ax.set_title(title)
        ax.set_xlabel("Candidate minus N3")
    fig.suptitle("Complex parents: saved structure and mixing diagnostics", fontsize=14, y=1.01)
    fig.tight_layout()
    writer.save("fig03_structure_and_blur", fig, records,
                caption_zh="complex 的 parent 配对结构与混合诊断。edge F1/IoU 的正差表示分数升高；component/hole MAE 的正差表示计数误差增加。计数采用 0.25 阈值、4-connected foreground 与 8-connected background，不能视为连续 topology 证书。intermediate-pixel fraction 和 magnitude TV 是描述性混合指标；不预设其增减等于好坏。前四组区间引用保存汇总，后两组按相同 2000 次 parent bootstrap 计算。",
                caption_en="Complex-parent paired differences against N3. Positive edge/IoU differences raise scores; positive component/hole MAE differences increase count errors. Digital topology uses the saved 0.25 threshold and 4-foreground/8-background connectivity. Intermediate-pixel fraction and magnitude TV are descriptive mixing diagnostics, with no assumed favorable direction.",
                sources=["results/a30/hull/HULL_MATERIAL_AND_STRUCTURE.csv", "results/a30/hull/HULL_AUDIT_SUMMARY.json"])


def _capacity_by_key(rows):
    result = {}
    for row in rows:
        if row.get("method") not in CAPACITY_METHODS:
            raise ValueError("UNKNOWN_CAPACITY_METHOD:" + str(row.get("method")))
        budget = int(row["budget"])
        pid = int(row["parent_id"])
        if budget not in BUDGETS or row.get("family") != _family(pid) or row.get("background_id") not in {"weak", "ring"}:
            raise ValueError("CAPACITY_METADATA_OUTSIDE_FROZEN_CASES")
        expected = f"p{pid}_{row['background_id']}_a1_snr20"
        if row["case_id"] != expected:
            raise ValueError("CAPACITY_CASE_ID_MISMATCH")
        key = (row["case_id"], row["method"], budget)
        if key in result:
            raise ValueError("DUPLICATE_CAPACITY_ROW:" + str(key))
        result[key] = row
    return result


def _capacity_parent(by_key, pid, method, budget):
    pair = [by_key.get((f"p{pid}_{bg}_a1_snr20", method, budget)) for bg in ("weak", "ring")]
    if any(row is None for row in pair):
        return {"parent_id": pid, "family": _family(pid), "method": method, "budget": budget, "status": "MISSING"}
    if any(not _valid(row) or _number(row.get("delta_nrmse")) is None for row in pair):
        return {"parent_id": pid, "family": _family(pid), "method": method, "budget": budget, "status": "INVALID"}
    baseline = [_number(row.get("N3_nrmse")) for row in pair]
    if any(value is None or value <= 0 for value in baseline):
        return {"parent_id": pid, "family": _family(pid), "method": method, "budget": budget, "status": "INVALID_BASELINE"}
    b = float(np.mean(baseline))
    value = float(np.mean([float(row["delta_nrmse"]) for row in pair]))
    return {"parent_id": pid, "family": _family(pid), "method": method, "budget": budget,
            "status": "VALID", "N3_parent_nrmse": b, "candidate_parent_nrmse": value,
            "gain_fraction": (b - value) / b, "gain_percent": 100 * (b - value) / b}


def plot_capacity(writer, capacity_rows, hull_summary):
    sources = ["results/a30/capacity/CAPACITY_CASES.csv", "results/a30/hull/HULL_AUDIT_SUMMARY.json"]
    if not capacity_rows:
        fig, ax = plt.subplots(figsize=(9.5, 3.8))
        _note_panel(ax, "E0 continuous-representation capacity: NOT_RUN",
                    "No saved CAPACITY_CASES rows are available.\nNo capacity curve is inferred from the hull results.")
        writer.save("fig04_capacity_parameter_budget", fig, [{"status": "NOT_RUN", "reason": "saved_capacity_rows_absent"}],
                    caption_zh="材料表示容量尚无实际保存结果，因此该图明确保留 NOT_RUN。16/32/64 参数预算的曲线只会在真实结果存在时生成。",
                    caption_en="No saved continuous-representation capacity rows exist. Parameter-budget curves remain NOT_RUN.",
                    sources=sources, status="NOT_RUN")
        return
    by_key = _capacity_by_key(capacity_rows)
    fig, axes = plt.subplots(1, 5, figsize=(16.0, 4.8), sharey=True)
    raw = [{"record_type": "source_case", **row} for row in capacity_rows]
    complete = True
    for ax, group in zip(axes, ("complex", *FAMILIES)):
        expected = _group_parents(group)
        for method in CAPACITY_METHODS:
            means, lows, highs = [], [], []
            for budget in BUDGETS:
                parents = [_capacity_parent(by_key, pid, method, budget) for pid in expected]
                valid = all(row["status"] == "VALID" for row in parents)
                if valid:
                    mean, lo, hi = _ci([r["candidate_parent_nrmse"] for r in parents])
                else:
                    mean = lo = hi = None
                    complete = False
                means.append(np.nan if mean is None else mean)
                lows.append(np.nan if lo is None else lo)
                highs.append(np.nan if hi is None else hi)
                raw.append({"record_type": "group_summary", "group": group, "method": method,
                            "budget": budget, "status": "COMPLETE" if valid else "INCOMPLETE",
                            "expected_parents": len(expected), "valid_parents": sum(r["status"] == "VALID" for r in parents),
                            "missing_parents": sum(r["status"] == "MISSING" for r in parents),
                            "invalid_parents": sum(r["status"].startswith("INVALID") for r in parents),
                            "mean_parent_nrmse": mean, "ci95_low": lo, "ci95_high": hi,
                            "bootstrap_samples": 2000, "bootstrap_seed": 291011})
            means, lows, highs = map(np.asarray, (means, lows, highs))
            ax.plot(BUDGETS, means, "o-", color=COLORS[method], lw=1.7, ms=4, label=LABELS[method])
            ax.fill_between(BUDGETS, lows, highs, color=COLORS[method], alpha=.11, linewidth=0)
        n3 = hull_summary["N3"]["complex"] if group == "complex" else hull_summary["N3"]["family"][group]
        baseline = float(np.mean([r["baseline"] for r in n3["parents"]]))
        ax.axhline(baseline, color="#475569", ls="--", lw=1, label="N3")
        raw.append({"record_type": "N3_reference", "group": group, "mean_parent_nrmse": baseline})
        ax.set_title(FAMILY_LABELS[group])
        ax.set_xticks(BUDGETS)
        ax.set_xlabel("Real parameter budget")
        ax.grid(color="#e2e8f0", lw=.7)
        missing = [r for r in raw if r.get("record_type") == "group_summary" and r.get("group") == group and r["status"] != "COMPLETE"]
        if missing:
            ax.text(.02, .98, f"{len(missing)} incomplete points\n(see raw table)", ha="left", va="top", fontsize=8,
                    transform=ax.transAxes, bbox={"facecolor": "white", "alpha": .8, "edgecolor": "none"})
    axes[0].set_ylabel("Achieved parent-mean delta-chi NRMSE")
    axes[-1].legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.suptitle("Saved offline representation fits at matched real parameter budgets", fontsize=14, y=1.02)
    fig.tight_layout()
    writer.save("fig04_capacity_parameter_budget", fig, raw,
                caption_zh="四种材料表示在共同 16/32/64 个实参数预算下的 achieved offline 材料误差。先平均每个 parent 的 weak/ring 两背景，再等权平均 parent；区间为 2000 次探索性 parent bootstrap。存在缺失或 invalid parent 的组不绘制完整统计点，原表保留所有失败行。非凸拟合的 achieved error 不是全局容量下界，曲线不提供测量可部署性证明。",
                caption_en="Achieved offline material-fit error for four representations at matched real parameter budgets. Conditions are averaged inside each parent before equal-parent aggregation. Missing or invalid parent pairs suppress the corresponding aggregate point. Achieved nonconvex fits are not global capacity certificates or deployable recovery results.",
                sources=sources, status="COMPLETE" if complete else "INCOMPLETE")


def plot_material_images(writer, hull_images, hull_rows, capacity_rows, capacity_images):
    ids = {str(case): i for i, case in enumerate(hull_images["case_id"])}
    metric = {(r["case_id"], r["method"]): r for r in hull_rows}
    cap_keys = _capacity_by_key(capacity_rows) if capacity_rows else {}
    cap_ids = {str(case): i for i, case in enumerate(capacity_images.get("case_id", []))}
    columns = ("truth", "N3", "full_hard", "full_hull", "full_L_mean", "best64")
    col_labels = ("Truth", "N3", "Full hard oracle", "Full convex oracle", "Full L mean", "Best-64 offline oracle")
    fig, axes = plt.subplots(4, 6, figsize=(14.7, 10.4))
    norm = Normalize(0., 1.2)
    raw, selection = [], []
    mappable = None
    complete = True
    for row_index, case in enumerate(FIXED_CASES):
        i = ids[case]
        bg = hull_images["background"][i]
        images = {method: hull_images[method][i] for method in columns[:-1]}
        best_method = None
        options = []
        for method in CAPACITY_METHODS:
            fit = cap_keys.get((case, method, 64))
            key = method + "_64"
            if fit is not None and _valid(fit) and _number(fit.get("delta_nrmse")) is not None and case in cap_ids and key in capacity_images:
                image = capacity_images[key][cap_ids[case]]
                if image.shape == (576,) and np.isfinite(image).all():
                    options.append((float(fit["delta_nrmse"]), CAPACITY_METHODS.index(method), method, image))
        if options:
            error, _, best_method, image = min(options, key=lambda x: (x[0], x[1]))
            images["best64"] = image
            selection.append({"case_id": case, "method": best_method, "budget": 64, "delta_nrmse": error,
                              "selection_role": "OFFLINE_PER_CASE_ORACLE; stable method order resolves ties"})
        else:
            complete = False
            selection.append({"case_id": case, "status": "NOT_RUN_OR_INVALID", "selection_role": "OFFLINE_PER_CASE_ORACLE"})
        for column_index, method in enumerate(columns):
            ax = axes[row_index, column_index]
            image = images.get(method)
            if image is None:
                ax.set_facecolor("#f1f5f9")
                ax.text(.5, .5, "NOT_RUN / INVALID", ha="center", va="center", fontsize=9, transform=ax.transAxes)
                raw.append({"record_type": "image_status", "case_id": case, "display_method": method, "status": "NOT_RUN_OR_INVALID"})
            else:
                delta = np.asarray(image) - bg
                magnitude = np.abs(delta).reshape(24, 24)
                mappable = ax.imshow(magnitude, origin="lower", cmap="viridis", norm=norm, interpolation="nearest")
                displayed_method = best_method if method == "best64" else method
                for pixel_row in range(24):
                    for pixel_col in range(24):
                        value = delta.reshape(24, 24)[pixel_row, pixel_col]
                        raw.append({"record_type": "material_pixel", "case_id": case, "family": FAMILIES[row_index],
                                    "display_method": method, "source_method": displayed_method,
                                    "pixel_row": pixel_row, "pixel_column": pixel_col,
                                    "delta_chi_real": float(value.real), "delta_chi_imag": float(value.imag),
                                    "delta_chi_magnitude": float(abs(value)), "colorbar_low": 0., "colorbar_high": 1.2,
                                    "clipped_above_colorbar": bool(abs(value) > 1.2)})
                if method != "truth":
                    error = min(options, key=lambda x: (x[0], x[1]))[0] if method == "best64" else float(metric[(case, method)]["delta_nrmse"])
                    subtitle = f"NRMSE {error:.3f}"
                    if method == "best64":
                        subtitle += "\n" + LABELS[best_method]
                    ax.text(.5, -.04, subtitle, ha="center", va="top", fontsize=8, transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
            if row_index == 0:
                ax.set_title(col_labels[column_index], fontsize=9)
            if column_index == 0:
                ax.set_ylabel(FAMILY_LABELS[FAMILIES[row_index]] + "\n" + case, fontsize=8)
    raw.extend({"record_type": "best64_selection", **row} for row in selection)
    fig.subplots_adjust(left=.13, right=.93, top=.91, bottom=.05, wspace=.09, hspace=.31)
    if mappable is not None:
        color_ax = fig.add_axes([.945, .19, .013, .57])
        fig.colorbar(mappable, cax=color_ax, label="Actual |delta chi| (fixed 0–1.2)")
    fig.suptitle("Fixed weak-background cases: saved material increments", fontsize=14, y=.97)
    writer.save("fig05_fixed_material_images", fig, raw,
                caption_zh="每个 family 预先固定的第一个 weak 条件：250065/250069/250073/250077。图中为实际 |χ−χ_b|，所有列共用 0–1.2 色标且不逐图归一化；原表保留实部、虚部、模与超色标标志。best-64 为这些固定 case 上四种合法表示的 offline 误差最小者，明确为 ORACLE；invalid/NaN 不作有效图像。观察图像不替代复杂 parent 总体统计或 topology 证书。",
                caption_en="The first weak-background case of each family is fixed before inspecting error. Images show actual |chi minus background| on one shared 0–1.2 scale; exported pixels retain real/imaginary components and clipping flags. Best-64 is an explicitly offline per-case oracle among valid saved representations. NaN/invalid fits are not rendered as valid materials.",
                sources=["results/a30/hull/HULL_RECONSTRUCTIONS.npz", "results/a30/hull/HULL_MATERIAL_AND_STRUCTURE.csv",
                         "results/a30/capacity/CAPACITY_CASES.csv", "results/a30/capacity/CAPACITY_IMAGES.npz"],
                status="COMPLETE" if complete else "CAPACITY_NOT_RUN_OR_INVALID",
                notes=["Pixel row/column coordinates; no physical coordinate or unit is invented.",
                       "Values above 1.2 saturate visually and remain unchanged in the raw table."])


def plot_parent_gains(writer, summary, capacity_rows):
    methods = list(HULL_METHODS) + [method + "_64" for method in CAPACITY_METHODS]
    data = np.full((16, len(methods)), np.nan)
    raw = []
    cap_keys = _capacity_by_key(capacity_rows) if capacity_rows else {}
    status = {}
    for j, method in enumerate(HULL_METHODS):
        for family in FAMILIES:
            for row in summary[method]["family"][family]["parents"]:
                i = PARENTS.index(int(row["parent_id"]))
                data[i, j] = 100 * float(row["gain"])
                status[(i, j)] = "VALID"
                raw.append({"parent_id": row["parent_id"], "family": family, "method": method,
                            "budget": "saved_bank", "status": "VALID", "baseline_parent_error": row["baseline"],
                            "candidate_parent_error": row["candidate"], "gain_fraction": row["gain"],
                            "gain_percent": data[i, j], "statistical_unit": "EXPOSED_PARENT_EXPLORATORY"})
    for offset, method in enumerate(CAPACITY_METHODS):
        j = len(HULL_METHODS) + offset
        for i, pid in enumerate(PARENTS):
            row = _capacity_parent(cap_keys, pid, method, 64)
            status[(i, j)] = row["status"]
            if row["status"] == "VALID":
                data[i, j] = row["gain_percent"]
            raw.append({**row, "statistical_unit": "EXPOSED_PARENT_EXPLORATORY"})
    finite = data[np.isfinite(data)]
    limit = max(50., math.ceil(float(np.max(np.abs(finite))) / 10.) * 10.)
    cmap = plt.get_cmap("RdBu").copy()
    cmap.set_bad("#e2e8f0")
    fig, ax = plt.subplots(figsize=(13.5, 8.6))
    image = ax.imshow(np.ma.masked_invalid(data), cmap=cmap, norm=TwoSlopeNorm(0., -limit, limit), aspect="auto")
    for i in range(16):
        for j in range(len(methods)):
            value = data[i, j]
            text = f"{value:+.1f}" if np.isfinite(value) else ("INV" if status[(i, j)].startswith("INVALID") else "NR")
            color = "white" if np.isfinite(value) and abs(value) > .58 * limit else "#0f172a"
            ax.text(j, i, text, ha="center", va="center", color=color, fontsize=8)
    labels = [LABELS[x].replace(" ", "\n", 1) for x in HULL_METHODS] + [LABELS[x].replace(" ", "\n", 1) + "\n64" for x in CAPACITY_METHODS]
    ax.set_xticks(range(len(methods)), labels)
    ax.set_yticks(range(16), [str(p) + "  " + FAMILY_LABELS[_family(p)] for p in PARENTS])
    for line in (3.5, 7.5, 11.5):
        ax.axhline(line, color="white", lw=2)
    ax.axvline(5.5, color="white", lw=3)
    ax.set_title("Every exposed parent: two-background material gain against N3", pad=14)
    fig.colorbar(image, ax=ax, fraction=.032, pad=.025, label="Parent gain (%)")
    fig.tight_layout()
    writer.save("fig06_every_parent_gain", fig, raw,
                caption_zh="全部 16 个 exposed parent 的材料相对改善，按 family 原顺序列出。每个单元先平均 weak/ring 误差，再计算 (N3−candidate)/N3。左侧引用保存 hull 汇总；右侧为 budget64 表示的实际保存结果。缺失标 NR、invalid 标 INV，均不替换为零。smooth 单列保留，不进入 12 complex parents 的主统计。",
                caption_en="All 16 exposed parents, with both backgrounds averaged before relative gain is calculated. Hull cells are copied from saved summaries; representation cells use saved budget-64 fits. NR and INV denote missing and invalid evidence, not zero gain. Smooth parents remain distinct from the complex primary group.",
                sources=["results/a30/hull/HULL_AUDIT_SUMMARY.json", "results/a30/capacity/CAPACITY_CASES.csv"])


def plot_e1_certificates(writer):
    candidates = ("results/a30/e1/CERTIFICATE_AUDIT.csv", "results/a30/compensation/CERTIFICATE_AUDIT.csv",
                  "results/a30/CERTIFICATE_AUDIT.csv")
    present = [source for source in candidates if _safe(writer.root, source).is_file()]
    if len(present) > 1:
        raise ValueError("AMBIGUOUS_E1_CERTIFICATE_SOURCES")
    rows = _read_csv(writer.root, present[0]) if present else []
    aliases = {"lower": ("certificate_lower", "lower", "lower_bound", "d_lower"),
               "upper": ("certificate_upper", "upper", "upper_bound", "d_upper"),
               "actual": ("actual_endpoint_distance", "actual_distance", "actual_d", "Maxwell_distance", "d_actual"),
               "linear": ("d_lin", "bounded_d_lin")}
    parsed, invalid = [], []
    for index, row in enumerate(rows):
        values = {name: next((_number(row[key]) for key in keys if key in row and _number(row[key]) is not None), None)
                  for name, keys in aliases.items()}
        status = str(row.get("status", "")).upper()
        if values["lower"] is None or values["upper"] is None or values["actual"] is None or status in {"NOT_RUN", "MISSING", "FAILED"}:
            invalid.append({"row_index": index, "reason": "actual_endpoint_or_interval_unavailable", **row})
            continue
        if values["lower"] < 0 or values["upper"] < values["lower"]:
            raise ValueError("INVALID_SAVED_E1_INTERVAL")
        parsed.append({"row_index": index, "pair_id": row.get("pair_id", row.get("case_id", str(index))), **values,
                       "source_status": status, "interpretation": "SAVED_REAL_ENDPOINT_AND_CERTIFICATE_VALUES"})
    if not parsed:
        fig, ax = plt.subplots(figsize=(9.5, 3.8))
        state = "NOT_RUN" if not rows else "NO_PLOTTABLE_REAL_RESULTS"
        _note_panel(ax, "E1 actual Maxwell certificate evidence: " + state,
                    "No saved actual endpoint distance and interval can be plotted.\nThe E0 hull certificate does not supply this E1 evidence.")
        writer.save("fig07_e1_certificate_status", fig,
                    [{"status": state, "reason": "actual_E1_endpoint_certificate_rows_absent"}, *invalid],
                    caption_zh="E1 只在保存了实际 Maxwell endpoint 距离与对应证书区间时绘制。当前无这类可绘制结果，图中明确记录状态；E0 的 simplex 数值证书不能填充 E1 物理结果。",
                    caption_en="E1 is plotted only when saved actual Maxwell endpoint distances and corresponding intervals exist. The E0 simplex certificate cannot substitute for these physical results.",
                    sources=present or list(candidates), status=state)
        return
    fig, ax = plt.subplots(figsize=(11.5, 4.8))
    x = np.arange(len(parsed))
    actual = np.asarray([r["actual"] for r in parsed])
    lower, upper = np.asarray([r["lower"] for r in parsed]), np.asarray([r["upper"] for r in parsed])
    ax.vlines(x, lower, upper, color="#087f8c", lw=2, label="Saved certificate interval")
    ax.scatter(x, actual, color="#b45309", s=23, label="Actual saved endpoint distance", zorder=3)
    linear = [r["linear"] for r in parsed]
    if all(value is not None for value in linear):
        ax.scatter(x, linear, color="#475569", marker="x", s=28, label="Saved bounded linear distance")
    ax.set_xticks(x, [str(r["pair_id"]) for r in parsed], rotation=90, fontsize=7)
    ax.set(xlabel="Saved pair ID", ylabel="Whitened response distance", title="Actual E1 endpoint distances and saved intervals")
    ax.grid(axis="y", color="#e2e8f0", lw=.7)
    ax.legend(frameon=False)
    fig.tight_layout()
    writer.save("fig07_e1_certificate_status", fig, [*parsed, *invalid],
                caption_zh="原保存 E1 certificate 区间与实际 Maxwell endpoint 距离逐对对照；仅绘制具有实际 endpoint 数值的行，缺失行在原表保留。图本身不改变证书前提、合法性或 parent gate 裁决。",
                caption_en="Saved E1 certificate intervals and actual Maxwell endpoint distances, pair by pair. Missing rows are retained in the raw table. This display does not judge certificate assumptions or assign a scientific gate.",
                sources=present, status="COMPLETE" if not invalid else "INCOMPLETE")


def generate(root: Path, *, output: Path | None = None, metadata_output: Path | None = None,
             include_saved_hull_details=False):
    start, cpu = time.perf_counter(), time.process_time()
    root = root.resolve()
    output = _safe(root, output or "results/a30/figures")
    metadata_output = _safe(root, metadata_output or "research/delegated/a30_plots")
    _style()
    writer = FigureWriter(root, output, metadata_output)
    summary, hull_rows, cases, weights, hull_images = _load_hull(root, details=include_saved_hull_details)
    capacity_rows = _read_csv(root, "results/a30/capacity/CAPACITY_CASES.csv", required=False)
    capacity_images = {}
    path = _safe(root, "results/a30/capacity/CAPACITY_IMAGES.npz")
    if path.is_file():
        with np.load(path, allow_pickle=False) as archive:
            capacity_images = {key: archive[key].copy() for key in archive.files}
        if len(capacity_images.get("case_id", [])) != 32 or set(map(str, capacity_images["case_id"])) != set(map(str, hull_images["case_id"])):
            raise ValueError("CAPACITY_IMAGE_CASE_LAYOUT_MISMATCH")
        for key, value in capacity_images.items():
            if key != "case_id" and value.shape != (32, 576):
                raise ValueError("CAPACITY_IMAGE_MATERIAL_LAYOUT_MISMATCH:" + key)
    if include_saved_hull_details:
        plot_hull_family(writer, summary)
        plot_weights_ess(writer, cases, weights)
        plot_structure_blur(writer, hull_rows, summary)
    plot_capacity(writer, capacity_rows, summary)
    plot_material_images(writer, hull_images, hull_rows, capacity_rows, capacity_images)
    plot_parent_gains(writer, summary, capacity_rows)
    timing = {"wall_seconds": time.perf_counter() - start, "process_cpu_seconds": time.process_time() - cpu,
              "GPU_seconds": 0, "Maxwell_calls": 0, "training_updates": 0, "locked_reads": 0,
              "charge_method": "nested evidence; use the inclusive a30_job wrapper receipt once"}
    writer.finish(timing)
    return {"figure_count": len(writer.mapping), "output": str(output), "metadata_output": str(metadata_output),
            "statuses": {row["figure"]: row["status"] for row in writer.mapping}, "timing": timing}
