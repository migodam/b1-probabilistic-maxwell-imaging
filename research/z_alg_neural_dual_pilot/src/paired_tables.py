"""Build the final Chinese result tables from completed Stage 2 receipts."""

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
from statistics import median


PILOT = Path(__file__).resolve().parents[1]
RESULTS = PILOT / "results"
DYNAMIC = ("gradient1", "gradient5", "neural", "hybrid", "no_feature")
FIXED = ("global_real", "fixed_16_0", "fixed_16_2", "fixed_portfolio")
ARM_NAMES = {
    "gradient1": "1 秒梯度", "gradient5": "5 秒梯度", "neural": "神经提案",
    "hybrid": "神经＋梯度", "no_feature": "固定特征",
    "global_real": "全局实部", "fixed_16_0": "固定 16 区 0",
    "fixed_16_2": "固定 16 区 2", "fixed_portfolio": "三个固定权重择优",
}
STUDY_NAMES = {"primary": "主比较", "supplement": "低噪声补充"}
COHORT_NAMES = {"Fresh6": "原生 6 件", "SRAP22": "SRAP 22 件"}
DECISIONS = ("PASS", "ZERO", "BELOW", "UNKNOWN")
EXPECTED = {"primary": (980, 392), "supplement": (440, 176)}
NUMBERS = (
    "noise", "u", "radius", "lower", "upper", "coarse_d", "proposal_wall", "proposal_cpu",
    "proof_wall", "proof_cpu", "charged_wall", "charged_cpu", "recovery_wall", "recovery_cpu",
)
INTEGERS = ("side", "passed", "coarse_nu", "nu", "nu_lower", "recovery_attempt")


def rows_from_csv(path):
    with path.open(newline="", encoding="utf-8") as f:
        records = list(csv.DictReader(f))
    for row in records:
        for key in NUMBERS:
            row[key] = float(row[key]) if row[key] else None
        for key in INTEGERS:
            row[key] = int(row[key]) if row[key] else None
    return records


def read_final(root):
    summary = json.loads((root / "summary.json").read_text())
    rows = rows_from_csv(root / "cases.csv")
    if summary["dynamic_count"] != sum(a for a, _ in EXPECTED.values()):
        raise ValueError("Wait for all 1,420 final dynamic receipts before generating tables.")
    if len({r["id"] for r in rows}) != len(rows):
        raise ValueError("Final case identifiers are duplicated.")
    for study, (expected_dynamic, expected_fixed) in EXPECTED.items():
        dynamic = [r for r in rows if r["study"] == study and r["arm"] in DYNAMIC]
        fixed = [r for r in rows if r["study"] == study and r["arm"] in FIXED]
        if (len(dynamic), len(fixed)) != (expected_dynamic, expected_fixed):
            raise ValueError(f"Final row counts disagree in {study}.")
        recovery = summary["costs"][study]["recovery_workers"]
        complete = {oid for worker in recovery if worker["finished"] for oid in worker["complete_objects"]}
        if complete != {r["oid"] for r in dynamic}:
            raise ValueError(f"Final recovery completion record is missing in {study}.")
    groups = defaultdict(Counter)
    keys = ("cohort", "noise", "u", "side", "arm")
    for row in rows:
        groups[tuple(row[k] for k in keys)][row["decision"]] += 1
        if row["passed"] != int(row["lower"] > row["radius"]):
            raise ValueError(f"Strict threshold disagrees with the published decision in {row['id']}.")
    existing = {tuple(s[k] for k in keys): Counter({d: s[d] for d in DECISIONS}) for s in summary["strata"]}
    if groups != existing:
        raise ValueError("The final summary and final cases table disagree on the strata.")
    return summary, rows


def med(values):
    values = list(values)
    return median(values) if values else None


def simple(value, digits=3):
    return "—" if value is None else f"{value:.{digits}f}"


def table(headers, rows):
    yield "| " + " | ".join(headers) + " |"
    yield "| " + " | ".join("---" for _ in headers) + " |"
    for row in rows:
        yield "| " + " | ".join(map(str, row)) + " |"
    yield ""


def case_counts(rows):
    count = Counter(r["decision"] for r in rows)
    return "/".join(str(count[k]) for k in DECISIONS)


def coarse_counts(rows):
    measured = [r for r in rows if r["coarse_d"] is not None]
    coarse = [r for r in measured if r["coarse_d"] > r["radius"]]
    return (
        len(coarse), sum(r["decision"] == "PASS" for r in coarse),
        sum(r["decision"] in ("ZERO", "BELOW") for r in coarse),
        sum(r["decision"] == "UNKNOWN" for r in coarse),
        sum(r["decision"] == "PASS" and r["coarse_d"] <= r["radius"] for r in measured),
    )


def context(row):
    return row["study"], row["cohort"], row["noise"], row["u"]


def context_label(key):
    study, cohort, noise, cap = key
    return STUDY_NAMES[study], COHORT_NAMES[cohort], f"{noise:g}%", f"{cap:g}"


def collection(rows):
    by_context = defaultdict(list)
    for r in rows:
        by_context[context(r)].append(r)
    return [(k, by_context[k]) for k in sorted(by_context)]


def dynamic_table(groups, *, transfer=False):
    for key, group in groups:
        for side in (4, 8):
            subset = [r for r in group if r["side"] == side and r["arm"] in DYNAMIC]
            n = len({r["oid"] for r in subset})
            cells = [
                [r for r in subset if r["arm"] == arm] for arm in DYNAMIC
            ]
            if any(len(cell) != n for cell in cells):
                raise ValueError("Dynamic stratum has incomplete object/arm matching.")
            counts = ("/".join(map(str, coarse_counts(cell))) if transfer else case_counts(cell) for cell in cells)
            yield (*context_label(key), side**2, n, *counts)


def fixed_table(groups):
    for key, group in groups:
        subset = [r for r in group if r["arm"] in FIXED]
        n = len({r["oid"] for r in subset})
        cells = [[r for r in subset if r["arm"] == arm] for arm in FIXED]
        if any(len(cell) != n for cell in cells):
            raise ValueError("Fixed stratum has missing object/arm matches.")
        yield (*context_label(key), n, *(case_counts(cell) for cell in cells))


def per_arm(rows):
    for study in EXPECTED:
        for side, arms in ((4, DYNAMIC), (8, DYNAMIC), (0, FIXED)):
            for arm in arms:
                subset = [
                    r for r in rows if r["study"] == study and r["arm"] == arm
                    and (side == 0 or r["side"] == side)
                ]
                yield study, side, arm, subset


def arm_totals(rows):
    for study, side, arm, subset in per_arm(rows):
        ids = {r["oid"] for r in subset}
        passes = {r["oid"] for r in subset if r["decision"] == "PASS"}
        unknown = {r["oid"] for r in subset if r["decision"] == "UNKNOWN"}
        yield (
            STUDY_NAMES[study], side**2 if side else "固定", ARM_NAMES[arm], len(subset),
            f"{len(passes)}/{len(ids)}", len(unknown), case_counts(subset),
            sum(r["recovery_attempt"] for r in subset),
            sum((r["nu_lower"] or 0) >= 10 or (r["nu"] or 0) >= 10 for r in subset),
        )


def arm_costs(rows):
    for study, side, arm, subset in per_arm(rows):
        yield (
            STUDY_NAMES[study], side**2 if side else "固定", ARM_NAMES[arm], len(subset),
            simple(med(r["proposal_wall"] for r in subset)),
            simple(med(r["proof_wall"] for r in subset)),
            simple(med(r["charged_wall"] for r in subset)),
            simple(med(r["charged_cpu"] for r in subset)),
            simple(med(r["charged_wall"] for r in subset if r["decision"] == "PASS")),
            simple(sum(r["charged_wall"] for r in subset)),
            simple(sum(r["charged_cpu"] for r in subset)),
        )


def section_times(summary, rows):
    for study, (expected, _) in EXPECTED.items():
        receipts = summary["costs"][study]
        dynamic = [r for r in rows if r["study"] == study and r["arm"] in DYNAMIC]
        # Hybrid includes the neural proposal shared with the zero-step arm.
        # The 5-second search already includes the 1-second prefix cost.
        unique = [r for r in dynamic if r["arm"] in ("gradient5", "hybrid", "no_feature")]
        recovered = [r for r in dynamic if r["recovery_attempt"]]
        worker_setup = receipts["proposal_setup"]
        native_worker = receipts["native_workers"]
        recovery_worker = receipts["recovery_workers"]
        label = STUDY_NAMES[study]
        times = [
            ("去重后提议时间", len(unique),
             sum(r["proposal_wall"] for r in unique), sum(r["proposal_cpu"] for r in unique)),
            ("补记的 32² 诊断", receipts["diagnostic_coarse_evals"],
             receipts["diagnostic_wall"], receipts["diagnostic_cpu"]),
            ("提议进程初始化", len(worker_setup),
             sum(r["wall"] for r in worker_setup), sum(r["cpu"] for r in worker_setup)),
            ("首轮验证读入端口及导入", len(native_worker),
             sum(r["read_holdout_and_import_setup"]["wall"] for r in native_worker),
             sum(r["read_holdout_and_import_setup"]["cpu"] for r in native_worker)),
            ("96² 限制上筛准备", len(native_worker),
             sum(r["restriction_setup"]["wall_s"] for r in native_worker),
             sum(r["restriction_setup"]["cpu_s"] for r in native_worker)),
            ("96² 首轮验证器加载几何", len(native_worker),
             sum(r["native_setup_and_cache"]["setup"]["load_fixed_geometry"]["wall_seconds"] for r in native_worker),
             sum(r["native_setup_and_cache"]["setup"]["load_fixed_geometry"]["cpu_seconds"] for r in native_worker)),
            ("96² 首轮逐案实际证明", expected,
             sum(r["proof_wall"] - r["recovery_wall"] for r in dynamic),
             sum(r["proof_cpu"] - r["recovery_cpu"] for r in dynamic)),
            ("补证进程加载几何", len(recovery_worker),
             sum(r["setup"]["setup"]["load_fixed_geometry"]["wall_seconds"] for r in recovery_worker),
             sum(r["setup"]["setup"]["load_fixed_geometry"]["cpu_seconds"] for r in recovery_worker)),
            ("补证新增逐案时间", len(recovered),
             sum(r["recovery_wall"] for r in recovered), sum(r["recovery_cpu"] for r in recovered)),
            ("固定对照完整进程", 1,
             summary["fixed_setup"][study]["whole_replay_wall"],
             summary["fixed_setup"][study]["whole_replay_cpu"]),
        ]
        for stage, count, wall, cpu in times:
            yield label, stage, count, simple(wall), simple(cpu)


def fixed_setup_rows(summary):
    for study in EXPECTED:
        for r in summary["fixed_setup"][study]["fixed_profile_setup"]:
            yield STUDY_NAMES[study], f"{r['cap']:g}", ARM_NAMES[r["arm"]], simple(r["wall_seconds"]), simple(r["cpu_seconds"])


def matched_rows(rows):
    index = {(r["study"], r["oid"], r["noise"], r["u"], r["side"], r["arm"]): r for r in rows}
    for study in EXPECTED:
        for side in (4, 8):
            for arm in ("neural", "hybrid"):
                pairs = defaultdict(list)
                for r in rows:
                    if not (r["study"] == study and r["side"] == side and r["arm"] == arm and r["passed"]):
                        continue
                    joint = (study, r["oid"], r["noise"], r["u"], side)
                    controls = [index[(*joint, tested)] for tested in ("gradient1", "gradient5")]
                    costs = [a["charged_wall"] for a in controls if a["passed"]]
                    if costs:
                        pairs[r["oid"]].append(min(costs) / r["charged_wall"])
                reduced = {name: median(g) for name, g in sorted(pairs.items())}
                yield (
                    STUDY_NAMES[study], side**2, ARM_NAMES[arm],
                    len(reduced), sum(map(len, pairs.values())), simple(med(reduced.values()), 2),
                    "、".join(k for k, v in reduced.items() if v >= 5.) or "—",
                    "、".join(f"{k} {v:.2f}×" for k, v in reduced.items()) or "—",
                )


def passing_ids(groups):
    for key, group in groups:
        by_status = defaultdict(list)
        for side, arms in ((4, DYNAMIC), (8, DYNAMIC), (0, FIXED)):
            for arm in arms:
                chosen = tuple(sorted({
                    r["oid"] for r in group if r["arm"] == arm and r["passed"]
                    and (side == 0 or r["side"] == side)
                }))
                if chosen:
                    label = (f"{side**2} 区 " if side else "") + ARM_NAMES[arm]
                    by_status[chosen].append(label)
        for ids, arms in sorted(by_status.items()):
            yield (*context_label(key), "、".join(arms), "、".join(f"`{name}`" for name in ids))


def render(summary, rows):
    groups = collection(rows)
    yield "# 封存比较的数值结果表"
    yield ""
    yield ("主比较包含 28 件材料的 40 个读数、98 个读数与上限组合，补充比较在同一批 22 件 SRAP 上加 44 个读数。"
           "两个区数、五个动态方法共得到 1,420 份完整的原生回执。一个物体的噪声档和上限属于配对条件。各表的对象数按独立材料计。")
    yield ""
    yield "## 严格判定分层"
    yield ""
    yield ("每格依次为 **过线／严格零／已证未过线／未定**。严格过线要求 96² 下端严格大于冻结半径。"
           "严格零要求上端为零。已证未过线要求正上端不超过半径。")
    yield ""
    yield from table(("研究", "材料批次", "噪声", "上限", "区数", "每臂物体数",
                      *(ARM_NAMES[a] for a in DYNAMIC)), dynamic_table(groups))
    yield "### 固定控制"
    yield ""
    yield "固定择优复用前三列的实际证书。其评估依赖共享固定几何准备，费用另列。"
    yield ""
    yield from table(("研究", "材料批次", "噪声", "上限", "每臂物体数",
                      *(ARM_NAMES[a] for a in FIXED)), fixed_table(groups))
    yield "## 粗网过线与原生复核"
    yield ""
    yield ("每格依次为 **粗网数值过线／其中原生过线／其中原生严格零或已证未过线／其中原生未定／粗网未过而原生过线**。"
           "粗网分数只参与诊断。")
    yield ""
    yield from table(("研究", "材料批次", "噪声", "上限", "区数", "每臂物体数",
                      *(ARM_NAMES[a] for a in DYNAMIC)), dynamic_table(groups, transfer=True))
    yield "## 各臂总计"
    yield ""
    yield ("过线与未定的物体数分别按至少一个读数或上限满足条件计，同一物体可以兼有两种结果。"
           "严格状态负方向至少 10 的计数只使用原生回执里的已证信息。")
    yield ""
    yield from table(("研究", "区数", "方法", "配对行", "过线物体", "未定物体",
                      "过／零／未过／未定", "补证行", "已证至少 10 个状态负方向"), arm_totals(rows))
    yield "## 费用与独立配对"
    yield ""
    yield ("逐方法表使用独立臂的计费。最早 1 秒和完整 5 秒是同一次梯度轨迹的两个读点，神经初值也由直接提案与精修共享。"
           "这些逐方法中位数用于配对比较，实际全程计费表只记一次完整梯度搜索与一次神经初值。单位均为秒。")
    yield ""
    yield from table(("研究", "区数", "方法", "配对行", "提议墙钟中位", "实际证书墙钟中位",
                      "独立计费墙钟中位", "独立计费 CPU 中位", "过线行独立墙钟中位",
                      "独立计费墙钟和", "独立计费 CPU 和"), arm_costs(rows))
    yield "### 实付分项"
    yield ""
    yield ("逐案实际证明包含失败和缓存后实际支出的时间。补证行的第一轮支出归入首轮，新增部分只计一次。"
           "各进程并行，墙钟列是各项付费墙钟相加，单个实验的历时见原进程回执。"
           "固定对照完整进程包含其共享准备，以下准备表单独展示其组成。")
    yield ""
    yield from table(("研究", "工序", "记录或进程数", "墙钟秒", "CPU 秒"), section_times(summary, rows))
    yield from table(("研究", "上限", "固定对照的共享准备", "墙钟秒", "CPU 秒"), fixed_setup_rows(summary))
    yield ("训练和冷验证采用不同的计时边界，另见训练、几何和独立冷启动回执。"
           " [训练记录](results/neural_dual_training.json) · [固定特征训练记录](results/no_feature_training.json)"
           " · [状态压缩计时](results/geometry_inertia.json) · [独立冷验证](results/cold_latency.json)")
    yield ""
    yield "### 同证书成功条件下的配对速度"
    yield ""
    yield ("每行先对两边均严格过线的同物体、噪声、上限和分区取最快梯度臂与所列方法的独立付费墙钟之比。"
           "先按物体取中位数，再汇总。没有成对过线时用横线。")
    yield ""
    yield from table(("研究", "区数", "方法", "配对物体", "配对条件", "物体中位比",
                      "物体中位达 5×", "逐物体中位比"), matched_rows(rows))
    yield "## 过线对象编号"
    yield ""
    yield "同一行中列出的各方法在该条件有完全相同的过线对象集合。全零或未过线的臂见分层计数。"
    yield ""
    yield from table(("研究", "材料批次", "噪声", "上限", "方法", "过线对象"), passing_ids(groups))
    yield ("原生半径、区间端点、方法编号和实际费用以逐案表为准。"
           " [逐案 CSV](results/cases.csv) · [汇总 JSON](results/summary.json) · [配对费用 CSV](results/paired_timing.csv)")


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--results", type=Path, default=RESULTS)
    cli.add_argument("--output", type=Path, default=PILOT / "EXPERIMENTAL_RESULTS.md")
    options = cli.parse_args()
    summary, rows = read_final(options.results)
    staged = options.output.with_suffix(".tmp")
    staged.write_text("\n".join(render(summary, rows)) + "\n", encoding="utf-8")
    staged.replace(options.output)


if __name__ == "__main__":
    main()
