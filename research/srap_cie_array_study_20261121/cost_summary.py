"""Account unique trajectories once, plus the two independent correction branches."""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"


def main():
    rows = []
    groups = defaultdict(list)
    for folder in ("predictions", "line_cap", "arc_cap", "lazy_fresnel", "legacy_mac"):
        for path in sorted((OUT / folder).glob("*.npz")):
            with np.load(path) as data:
                info = json.loads(str(data["cost_json"]))
                records = json.loads(str(data["rows_json"]))
            first = records[0]
            raw_wall, raw_cpu, counts = info["raw"]["raw_spent"]
            counts = Counter(counts)
            wall = raw_wall
            cpu = raw_cpu
            for light in info["light"]:
                counts.update(light["counts"])
                wall += light["additional_wall"]
                cpu += light["additional_cpu"]
            cohort = first["cohort"] if folder == "predictions" else folder
            row = dict(cohort=cohort, host=first["host"], case=first["case"], method=first["method"],
                       raw_wall_s=raw_wall, light_additional_wall_s=wall-raw_wall, paid_wall_s=wall,
                       paid_cpu_s=cpu, max_endpoint_overrun_s=max(r["overrun"] for r in records),
                       unscored_prefixes=sum(r["phi"] is None for r in records),
                       status=first["proposal_status"].splitlines()[0],
                       common_dense_LU=counts["LU"], complex_RHS=counts["complex_RHS"],
                       PR_steps=counts["PR_CG_steps"], AP_cycles=counts["AP_cycles"],
                       FFT_CG_steps=counts["fft_pcg_iterations"], full_counts=json.dumps(dict(counts)))
            rows.append(row)
            groups[cohort, first["host"]].append(row)
    with (OUT / "cost_per_trajectory.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    text = ["# 实际计算费用", "",
            "单臂实际墙钟合计包含原始轨迹与两条额外的 15/30 秒轻校正分支。两个 raw 前缀共享轨迹，仅计一次。状态构建、FFT、背景谱、物理审计、未入选状态和截止后才返回的计算均计入。各并行臂墙钟之和计作计算使用量。", "",
            "| 数据组 | 主机 | 轨迹 | 累计臂墙钟 s | 累计单进程 CPU s | 无可用材料的前缀 | 超过预算返回最大 s |",
            "|---|---|---:|---:|---:|---:|---:|"]
    for (cohort,host), jobs in sorted(groups.items()):
        text.append(f"|{cohort}|{host}|{len(jobs)}|{sum(r['paid_wall_s'] for r in jobs):.1f}|{sum(r['paid_cpu_s'] for r in jobs):.1f}|{sum(r['unscored_prefixes'] for r in jobs)}|{max(r['max_endpoint_overrun_s'] for r in jobs):.2f}|")
    totals = Counter(r["status"] for r in rows)
    text += ["",f"共 {len(rows)} 条二维轨迹，累计 {sum(r['paid_cpu_s'] for r in rows)/3600:.3f} CPU 小时。原始求解停止状态为 {dict(totals)}。共同 corrector 的独立状态、截止前完成候选的可用时间见逐例表。所有原始求解臂只有各自 15、30 秒前已完成审计的候选能参加相应预算的比较。单次不可打断的操作可能越过截止，迟到结果无法入选。", "",
             "每个任务读取观测包、导入运行库及最后序列化状态发生在臂计时以外；冷启动物理算子及分解发生在臂内。并发运行会影响可达的阶段，因此同一对象的两臂在同机顺序测量。不同硬件的胜率只在各机配对后计入。", "",
             "[逐轨迹完整计数与超时](results/cost_per_trajectory.csv) · [同工作量跨主机测量](engine_dev/HARDWARE.md) · [二维测量生成计数](data/generation_cost.json) · [补集生成计数](data_followup/generation_cost.json) · [三维 CUDA 计数](vector/REPORT_ZH.md)", "",
             "新 24 件二维材料共享 744 个细网前向右端。生成费用和线性方程残差逐对象存档；实测无需额外目标生成。矩阵自由 H 求解的 FFT 与合法 Fourier 回注的 FFT 分开记数。二维 32² 的原 Maxwell 审计仍用稠密材料相关 LU，60² 使用既有大网格物理正演。不同网格的求解计数应结合实际时间阅读。"]
    for file in (ROOT/"data_followup/generation_cost.json", OUT/"mac_sampled_memory.json"):
        if not file.exists(): continue
        info = json.loads(file.read_text())
        if "max_worker_bytes" in info:
            text += ["",f"Mac 工作进程的5秒间隔采样中，最大单进程驻留内存为 {info['max_worker_bytes']/2**20:.1f} MiB；短于采样间隔的峰值未知。"]
    vec = ROOT/"vector/predictions"
    receipts = [json.loads(p.read_text()) for p in sorted(vec.glob("*/SEALED.json"))]
    if len(receipts)==48:
        text += ["",f"三维 48 个正式冷启动臂合计 {sum(z['total_s'] for z in receipts):.1f} 臂墙钟秒，包含两条 corrector。中途纠正共面圆环为矩形面阵的费用及三网场误差另见三维报告。二维线程与三维作业同时运行；GPU 内核运行时间没有通过设备 profiler 单列。"]
    text += ["", "评分耗时与开发开销分别列在各 cohort 的执行记录、同状态诊断与三维报告。部分早期开发诊断仅记录关键操作耗时，未记录启动到退出的总 CPU，整项研究的 CPU 总额高于二维反演小计。[公开开发诊断](engine_dev/REPORT_ZH.md)"]
    (ROOT/"COST_LEDGER.md").write_text("\n".join(text)+"\n")


if __name__ == "__main__":
    main()
