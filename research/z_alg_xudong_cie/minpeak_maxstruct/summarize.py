#!/usr/bin/env python3
"""Seal the online choices, replay scalar controls and produce the paired report."""

from datetime import datetime, timezone
from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from threadpoolctl import threadpool_limits

from minpeak import HERE, RESEARCH, TRANSFER, clock, load, read, timing, write
from core.certificate import certify
from core.inertia import Prepared


def checked_scalar_controls(entries, G, S, B, Y, radii, lower):
    prior = HERE / "matched_scalar_controls.json"
    if prior.is_file():
        return read(prior)
    started = clock()
    rows = {}
    for i, entry in enumerate(entries):
        case = entry["case"]
        result = lower["cases"][case]
        passing = [trial for trial in result["all_numeric_trials"]
                   if trial["stage"] == "mixed_global_seed"
                   and (trial["lower_bound_numeric"] or 0.) > radii[i]]
        passing.sort(key=lambda r: (-r["cap"], -r["lower_bound_numeric"]))
        tried = []
        bound = result["historical"]["best"]["cap"]
        if passing:
            row = passing[0]
            actual = result["selected"]
            if actual is not None and row == actual["proposal"]:
                proof, reused = actual["proof"], True
            else:
                proof = certify(Prepared(G, S, B, Y[i]), row["cap"], row["sector"][0], row["optical"][0],
                                radii={"original_CN99": radii[i]}, precision=128, threads=1)
                reused = False
            tried.append(dict(cap=row["cap"], mu=row["sector"][0], lam=row["optical"][0],
                              reused_selected_proof=reused, proof=proof))
            if proof.get("strictly_excludes", {}).get("original_CN99"):
                bound = row["cap"]
        rows[case] = dict(archival_cap=result["historical"]["best"]["cap"], final_scalar_cap=bound,
                          top_numerical_candidate_count=len(passing), checked=tried)
    record = dict(scope="Highest passing numerical scalar trial in the same per-case timing trace.",
                  matching="Same 24² geometry, caps, observations, residual criterion and stored radii.",
                  cases=rows, **timing(started))
    write(prior, record)
    return record


def plot_brackets(rows):
    fig, ax = plt.subplots(figsize=(8.25, 4.25), layout="constrained")
    for i, row in enumerate(rows):
        if row["u_feas"] is not None:
            ax.plot([row["u_cert"], row["u_feas"]], [i, i], lw=5, color="#3274a1", alpha=.85,
                    solid_capstyle="round", label="Certified interval for minimum peak" if i == 0 else None)
            if row["truth_peak"] >= row["u_feas"]:
                ax.plot([row["u_feas"], row["truth_peak"]], [i, i], lw=1.7, color="#9e490b", ls="--",
                        label="Verified reduction from generating peak" if i == 0 else None)
            ax.plot(row["u_feas"], i, "o", color="#153f59", ms=5)
        ax.plot(row["u_archival"], i, marker="|", color="#48525d", ms=13,
                label="Inherited lower certificate" if i == 0 else None)
        ax.plot(row["truth_peak"], i, marker="D", color="#9e490b", ms=5)
    ax.set_yticks(range(len(rows)), [r["case"][-4:] for r in rows])
    ax.invert_yaxis()
    ax.set(xlabel="Peak dielectric contrast", ylabel="R14 material",
           title="Minimum peak on each original CN99 noise ball")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=7))
    ax.grid(axis="x", color="#dfe5eb")
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.15), ncol=2, frameon=False, fontsize=8)
    fig.savefig(HERE / "minpeak_brackets.png", dpi=170)
    fig.savefig(HERE / "minpeak_brackets.pdf")
    plt.close(fig)


def plot_materials(rows):
    fig, axes = plt.subplots(2, 4, figsize=(9.3, 5.3), layout="constrained")
    plotted = []
    for ax, row in zip(axes.flat, rows):
        if row["u_feas"] is None:
            ax.text(.5, .5, "Unresolved", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            continue
        with np.load(HERE / "online" / f"{row['case']}_feasible.npz", allow_pickle=False) as f:
            m = f["material"]
        image = ax.imshow(m.reshape((24, 24)), vmin=0., vmax=max(r["u_feas"] or 0. for r in rows),
                          origin="lower", cmap="viridis", interpolation="nearest")
        plotted.append(ax)
        ax.set_title(f"{row['case'][-4:]}  ·  peak {row['u_feas']:.2f}", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
    axes.flat[-1].set_axis_off()
    if plotted:
        fig.colorbar(image, ax=plotted, shrink=.81, label="Dielectric contrast")
    fig.suptitle("Strictly feasible lower-peak materials")
    fig.savefig(HERE / "minpeak_materials.png", dpi=170)
    fig.savefig(HERE / "minpeak_materials.pdf")
    plt.close(fig)


def report(rows, costs):
    lines = [
        "# 七件 R14 对象的最小峰值：证书下界与可行上界",
        "",
        "对于每份原始带噪端口数据，我们寻找在同一误差球内能达到的最低介电峰值。"
        "严格证书排除峰值不超过下端的任意材料，完整全波验证则给出达到上端的一张可行材料图。"
        "两端之间是仍未确定的最小峰值范围。",
        "",
        "| 对象 | 原 CN99 半径 | 历史严格下界 | 同批标量下界 | 64 块搜索后下界 | 全波可行峰值 | 全波残差/半径 | 生成峰值 | 已证实的降峰幅度 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        upper = r["u_feas"]
        u_text = "未确定" if upper is None else f"{upper:.5f}"
        ratio_text = "—" if upper is None else f"{r['residual_over_radius_upper']:.7f}"
        reduction = r["ambiguity_peak_reduction_lower"]
        red_text = "—" if reduction is None else f"{reduction:.5f}"
        lines.append(f"| {r['case'][-4:]} | {r['radius']:.8f} | {r['u_archival']:.2f} | "
                     f"{r['u_scalar']:.2f} | {r['u_cert']:.2f} | {u_text} | "
                     f"{ratio_text} | {r['truth_peak']:.5f} | {red_text} |")
    complete = sum(r["u_feas"] is not None for r in rows)
    ambiguous = sum((r["ambiguity_peak_reduction_lower"] or 0.) > 0 for r in rows)
    extension = sum(r["u_cert"] > r["u_archival"] for r in rows)
    block = sum(r["u_cert"] > r["u_scalar"] for r in rows)
    bracketed = [r for r in rows if r["u_feas"] is not None]
    wide = max(bracketed, key=lambda r: r["bracket_width"], default=None)
    narrow = min(bracketed, key=lambda r: r["bracket_width"], default=None)
    contrast = (f"{wide['case'][-4:]} 的区间最宽，为 ({wide['u_cert']:.2f},{wide['u_feas']:.5f}]。"
                f"{narrow['case'][-4:]} 的区间最窄，为 ({narrow['u_cert']:.2f},{narrow['u_feas']:.5f}]。"
                if wide else "")
    coverage = [r["truth_residual_over_radius"][1] for r in rows]
    lines.extend([
        "",
        f"**{complete}/7 获得严格可行上界，{ambiguous}/7 获得可验证的降峰。** "
        f"在线搜索较已发表下界提高 {extension}/7 件，高于同批被检验标量方案的有 {block}/7 件。"
        "同批标量对照来自本轮记录的候选，覆盖范围由搜索轨迹明示。",
        contrast,
        "",
        "![最小峰值证书区间](minpeak_brackets.png)",
        "",
        "![同色标的严格可行材料](minpeak_materials.png)",
        "",
        "## 完整物理口径",
        "",
        "七件都是 R14 已公开的 24→24、12 发 48 收全环观测，使用原盘脉冲离散标量 TM 算子。"
        "观测包的 G、B、S、Y 和各自半径均与原数组逐项相同。"
        "各源共用 48 个已测接收位置，一次使用全部 576 个复端口。半径沿用原 1% 复 RMS、CN99 噪声合同，"
        "本匹配网格没有额外模型失配补偿。原材料先验为 [0.02,20]。"
        "排除下界覆盖 [0,u]，输出见证落在原先验内，因此区间对两种下限均成立。",
        "",
        "上界搜索只读取合法的历史重建及原始读数。按全波残差和材料峰值选历史初值，"
        "用精确切向的带盒 Gauss–Newton 降低峰值，逐次用全波端口接受下降。"
        "最后对保存候选独立执行区间求解，验证整个原离散状态方程与端口残差。"
        "某个峰值下的局部迭代未进噪声球，只意味着本次找不到该上界。",
        "",
        "下界复用已发表的同算子证书。64 块按每组 3×3 原生单元划分，"
        "在存档标量权重种子上分别更新实部扇区与光学权。每件的在线择权在约一秒内结束。"
        "入选后重新用 128 位区间核对惯性与下界；区间复核独立计时。"
        "本轮线上计算未读取原材料。封存后才读取同一观测的原区间计分标量，"
        f"原生成材料的全波残差区间上端为各自半径的 {100 * min(coverage):.2f}% 至 "
        f"{100 * max(coverage):.2f}%。 "
        "[原存档的离线区间评分](../../z_alg_optical_inertia/validation/max_struct/offline_scalar_scores.json)",
        "",
        "表中最后一列比较已验证的可行材料与生成峰值。如果该差为正，同一批数据至少允许这一幅度的峰值降低。"
        "上下界之间仍包含证书松弛和搜索未解决的部分，无法用这两端分别定量归因，也无法推出精确的最小峰值。",
        "",
        "## 费用",
        "",
        f"- Apple M4 上，新上界搜索及严格审计共 {costs['upper_wall_s']:.2f} 秒墙钟，"
        f"{costs['upper_cpu_s']:.2f} 秒 CPU。在线求解使用 {costs['physical_factorizations']} 次"
        f"原生全波分解、{costs['source_rhs']} 个源右端和 {costs['receiver_rhs']} 个接收右端。"
        f"独立区间审计另用 {costs['strict_inverse_rhs']} 个近似逆右端、"
        f"{costs['strict_source_rhs']} 个源右端、{costs['strict_receiver_rhs']} 个接收右端。",
        f"- 64 块数值择权合计 {costs['gradient_wall_s']:.2f} 秒墙钟，单件最高 "
        f"{costs['gradient_wall_max_s']:.3f} 秒。同一观测的严格证书证明与标量复核分别计费，"
        f"各为 {costs['lower_proof_wall_s']:.2f} 秒及 {costs['scalar_control_wall_s']:.2f} 秒墙钟。",
        f"- 供筛选的 {costs['historical_cost_rows']} 个已计费历史端点原共付 "
        f"{costs['historical_receipt_cpu_s']:.2f} 秒 CPU，"
        f"其中七个最终选作初值的端点为 {costs['historical_selected_seeds_cpu_s']:.2f} 秒。"
        "历史费用与本机回放费用分别记账。初期对两种归档回执及 JSON 序列化的修复见本地失败日志。",
        "",
        "[机器可读区间](summary.json) · [逐例 CSV](brackets.csv) · [在线封存记录](online_seal.json) · "
        "[上界证据和可行数组](online/) · [下界证书](lower_certificates.json) · "
        "[标量复核](matched_scalar_controls.json) · [复现入口](REPRODUCE.md)",
    ])
    (HERE / "REPORT_ZH.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    entries, G, S, B, Y, radii, _, _ = load(TRANSFER / "public_inputs.npz")
    upper = read(HERE / "online/online.json")
    lower = read(HERE / "lower_certificates.json")
    case_names = [e["case"] for e in entries]
    if list(upper["observations"]) != case_names or list(lower["cases"]) != case_names:
        raise ValueError("The seven online cases must all finish before scoring.")
    replay = checked_scalar_controls(entries, G, S, B, Y, radii, lower)
    seal = HERE / "online_seal.json"
    if not seal.is_file():
        write(seal, dict(cases=case_names, sealed_utc=datetime.now(timezone.utc).isoformat(),
                         artifact_inputs=["online/online.json", "lower_certificates.json",
                                          "matched_scalar_controls.json"],
                         truth_access_during_online="NONE"))
    # Independently verified interval scalars already stored against these exact G, B, S, Y and radii.
    offline = read(TRANSFER / "offline_scalar_scores.json")
    scored = []
    for i, name in enumerate(case_names):
        lo = lower["cases"][name]
        hi = upper["observations"][name]
        ref = offline["rows"][name]
        radius = float(radii[i])
        if ref["public_CN99_radius"] != radius or ref["status"] != "inside":
            raise ValueError(f"Reference feasibility mismatch for {name}.")
        audit = hi.get("accepted")
        chosen = lo["selected"]
        numerical_certificate = (chosen["proof"]["certified_frobenius_residual_lower"] if chosen
                                  else lo["historical"]["best"]["lower"])
        ceiling = audit["peak"] if audit else None
        if audit and not lo["final_cap"] < ceiling:
            raise ArithmeticError(f"Crossed strict bounds for {name}.")
        field = dict(case=name, radius=radius, u_archival=lo["historical"]["best"]["cap"],
                     u_scalar=replay["cases"][name]["final_scalar_cap"], u_cert=lo["final_cap"],
                     u_feas=ceiling, certified_lower_field_distance=numerical_certificate,
                     residual_over_radius_upper=float(np.nextafter(audit["full_wave_distance_hi"]/radius,
                                                                   np.inf)) if audit else None,
                     full_wave_residual_interval=[audit["full_wave_distance_lo"], audit["full_wave_distance_hi"]]
                     if audit else None, minimum_contrast=audit["minimum"] if audit else None,
                     truth_peak=ref["truth_peak"], truth_peak_interval=ref["peak_interval"],
                     truth_residual_over_radius=ref["true_residual_over_radius"],
                     ambiguity_peak_reduction_lower=float(np.nextafter(ref["peak_interval"][0] - ceiling,
                                                                         -np.inf)) if audit else None,
                     bracket_width=float(ceiling - lo["final_cap"]) if audit else None,
                     online_outcome=hi["status"], chosen_lower_source=chosen["proposal"]["stage"]
                     if chosen else "inherited")
        scored.append(field)
    upper_rows = list(upper["observations"].values())
    audits = [x for row in upper_rows for x in row.get("final_audit_trials", [])]
    lower_rows = list(lower["cases"].values())
    archival_cpu = {}
    for row in upper_rows:
        for endpoint in row["history_endpoints"]:
            if endpoint["status"] == "numeric_replay":
                path = RESEARCH / endpoint["source"]
                receipt = read(path.parent / "RECEIPT.json")
                archival_cpu[endpoint["method"]] = receipt.get("whole_process_CPU_seconds",
                                                               receipt.get("CPU_seconds"))
    available_cpu = [value for value in archival_cpu.values() if value is not None]
    seed_cpu = [archival_cpu.get(r["warm_start"]["method"]) or 0. for r in upper_rows]
    costs = dict(upper_wall_s=upper["wall_seconds"], upper_cpu_s=upper["cpu_seconds"],
                 physical_factorizations=sum(r["operations"].get("physical_factorizations", 0) for r in upper_rows),
                 source_rhs=sum(r["operations"].get("physical_source_rhs", 0) for r in upper_rows),
                 receiver_rhs=sum(r["operations"].get("physical_receiver_rhs", 0) for r in upper_rows),
                 strict_inverse_rhs=sum(r["inverse_rhs"] for r in audits),
                 strict_source_rhs=sum(r["source_rhs"] for r in audits),
                 strict_receiver_rhs=sum(r["receiver_rhs"] for r in audits),
                 gradient_wall_s=sum(r["online_multiplier_search"]["wall_seconds"] for r in lower_rows),
                 gradient_wall_max_s=max(r["online_multiplier_search"]["wall_seconds"] for r in lower_rows),
                 lower_proof_wall_s=sum(p["proof"]["cost"]["wall_seconds"] for r in lower_rows
                                        for p in r["verified_proposals"]),
                 scalar_control_wall_s=replay["wall_seconds"],
                 historical_cost_rows=len(available_cpu), historical_receipt_cpu_s=sum(available_cpu),
                 historical_selected_seeds_cpu_s=sum(seed_cpu),
                 lower_total_wall_s=lower["wall_seconds"], lower_total_cpu_s=lower["cpu_seconds"])
    write(HERE / "summary.json", dict(scope="R14 seven paired native synthetic observations",
                                      source_model=upper["shared_geometry"], rows=scored, costs=costs))
    with (HERE / "brackets.csv").open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=list(scored[0]))
        writer.writeheader()
        writer.writerows(scored)
    report(scored, costs)
    plot_brackets(scored)
    plot_materials(scored)


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
