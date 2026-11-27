"""Build the matched-observation decision table from sealed campaign receipts."""
from pathlib import Path
import json

from ports import dump
import numpy as np

HERE=Path(__file__).resolve().parent


def read(p):
    return json.loads(p.read_text())


def compact():
    config=read(HERE/"campaign.json")
    cases=[]
    for item in config["objects"]:
        case=item["case"]
        p=HERE/"loop"/case
        design=read(p/"design_frozen.json")
        result=read(p/"attack_results.json")
        assert result["stop"].startswith("equal bounded")
        init=read(p/"initial_receipt.json")
        acq=read(p/"acquisition_receipt.json")
        truth=read(p/"postseal_score.json")
        strict=read(p/"strict_attacks.json")
        assert acq["frozen_design_mtime_ns"]<acq["observation_mtime_ns"]
        arms={}
        for name,row in result["results"].items():
            source=row.get("alias_of",name)
            fit=result["results"][source]
            classes={f:max(r["best"]["residual_over_each_radius"]) for f,r in fit["families"].items()}
            span=fit["numeric_inner_ROI_range"]
            arms[name]=dict(cut=design["actions"][name]["cut"],tasks=classes,
                feasible_family_count=sum(res<1 for res in classes.values()),
                numeric_inner_span=span,numeric_inner_width=span[-1]-span[0] if len(span)==2 else None,
                strict_inner_width=strict["arms"].get(source,{}).get("ROI_width_lower"),
                safe_outer_width=4.5,registered_rho=config["rho"],
                predicted_candidate_D_over_2r=design["actions"][name]["predicted_worst_pair_D_over_2r"],
                median_old_candidate_new_residual=float(np.median(list(fit["fixed_old_material_replay_new_over_radius"].values()))),
                offline_truth_old_residual=truth["sealed_truth_residual_old_over_r"],
                offline_truth_new_residual=truth["sealed_truth_residual_new_over_r"][name],
                observation_id=acq["actions"][name]["observation_id"],
                available_W_sum=acq["actions"][name]["source_energy_sum_W_times_pulse_duration"],
                accepted_W_sum=acq["actions"][name]["real_driven_accepted_W_sum"],
                role_changes=acq["actions"][name]["port_role_changes_from_old_A"])
        n=sum(a["feasible_family_count"]>0 for n,a in arms.items() if n not in ("A_repeat","A_reverse"))
        best_name,best_arm=min(((n,a) for n,a in arms.items() if n not in ("A_repeat","A_reverse")),
                   key=lambda entry:entry[1]["tasks"]["unrestricted"])
        record=dict(case=case,initial_observation_id=init["observation_id"],initial_feasible_candidate_pairs=not design["candidate_pair_fallback"],
            ROI_cells=design["ROI_cells"],patch_cells=design["connected_patch_cells"],
            ROI_area_m2=design["ROI_area_m2"],patch_area_m2=design["patch_area_m2"],arms=arms,
            switched_arms_with_any_feasible_task_candidate=n,
            best_switched_unrestricted_method=best_name,best_switched_unrestricted_residual=best_arm["tasks"]["unrestricted"],
            design_cpu_s=design["cpu_s"],attack_cpu_s=result["cpu_s"],strict_cpu_s=strict["cpu_s"])
        cont=p/"continuation_receipts.json"
        if cont.exists():
            supplement=read(cont)
            if "stop" in supplement:
                record["supplement"]=dict(cpu_s=supplement["cpu_s"], arms={
                    k:{f:max(q["final"]) for f,q in arm["families"].items()} for k,arm in supplement["arms"].items() if "families" in arm})
        cases.append(record)
    p=HERE/"loop"/"fresh_bg_b"
    followup={}
    for name in ("dual_cie","long_gn"):
        src=p/f"{name}_verified.json"
        if src.exists():
            followup[name]=read(src)
    budget=p/"budget_witness_verified.json"
    if budget.exists():
        followup["budget_matching"]=read(budget)
    truth={r["case"]:read(HERE/"loop"/r["case"]/"postseal_truth_verified.json") for r in cases}
    dump(HERE/"active_compact.json",dict(active=cases,solver_followup=followup,offline_truth_verification=truth,
           decisive_outcome="Both B and the compensation-directed cut retain strictly verified 0.5 mean ambiguity on fresh_bg_b within the preregistered exterior area-weighted L2 cap 0.35. Baseline-budget failures on the remaining rotated arms are inconclusive. Full-class outer widths retain 4.5.",
          calibration="fixed nominal exact stored operator",frequency_MHz=400,grid_side=32))
    return cases


def md(cases):
    lines=["## 主动试验结果","",
           "用每臂自己的 A 与追加读数联合复拟合。在四件材料上，重复 A 与 A 整体反转都保留宽的有限面积均值歧义。表中每个宽度来自两个通过 Arb 独立残差审计的材料。旋转布局列给出不限制 ROI 时两个观测球的最大残差倍数，数值低于 1 才有严格审核的资格。","",
           "| 材料 | A 复读严格均值宽度 | A 反向严格均值宽度 | B 无限制最大残差 | T 判别 | 补偿引导 | 同预算旋转布局最小 unrestricted 残差 |",
           "|---|---:|---:|---:|---:|---:|---:|"]
    for r in cases:
        a=r["arms"]
        width=lambda key: f"{a[key]['strict_inner_width']:.3f}" if a[key]["strict_inner_width"] is not None else "未严格审计"
        f=lambda key: a[key]["tasks"]["unrestricted"]
        lines.append(f"| {r['case']} | {width('A_repeat')} | {width('A_reverse')} | {f('B_interleaved'):.3f} | {f('T_discrimination'):.3f} | {f('compensation'):.3f} | {r['best_switched_unrestricted_residual']:.3f}（{r['best_switched_unrestricted_method']}） |")
    counts=sum(r["switched_arms_with_any_feasible_task_candidate"] for r in cases)
    rotated=sum(sum(n not in ("A_repeat","A_reverse") for n in r["arms"]) for r in cases)
    lines+=["",f"{rotated} 个逐材料旋转布局中，{counts} 个在首轮同预算五个材料家族攻击中找到数值相容端点。独立离线复算封存生成材料，四件的旧 A 和所有新布局的严格残差上界均不超过 0.677r。同一数据有物理解。旋转布局下较大的已见残差反映求解支出和初始化，首轮结果不能用于全类排除。全部外宽仍为 4.5。","",
            f"四件材料的标准反演合计 {sum(r['attack_cpu_s'] for r in cases):.1f} 秒 CPU，A 数据下的设计合计 {sum(r['design_cpu_s'] for r in cases):.1f} 秒 CPU。所有历史、动作及逐臂外区变化、真值材料误差和八个未读取终端的仅离线复算比分开保存。","",
            "[分析数据](active_compact.json) · [逐布局及家族表](active_summary.csv) · [首轮残差及见证宽度图](figures/closed_loop_residuals_and_spans.png)",""]
    lines+=supplement_md()
    lines+=["## 96² 原生与更细格点复算","",
            "冻结布局只使用七件 32² 历史见证。两件原生 96² 材料的 A、B 固定对分离及用 32² 数据训练的精确 maximin 布局如下。全部为固定对数值结果，复用已经独立在 96² 求出的旧材料。","",
            "| 物体 | A | B | 32² 训练 maximin | 事后 96² 固定对最高值 |",
            "|---|---:|---:|---:|---:|"]
    native=read(HERE/"native"/"report.json")
    for c in native["native"]:
        f=c["frozen_action_ratios"]
        lines.append(f"| {c['case']} | {f['A']:.3f} | {f['B']:.3f} | {f['exact_train_maximin']:.3f} | {c['best_allcuts_ratio']:.3f} |")
    lines+=["","n96_c01 的两份材料还以分片常数细分到 192²，由重新构造的裸空间算子与目标加载电路独立复算。两份细格点响应相对于原 96² A 带噪观测的残差为 0.692r 和 0.728r。这一对材料的跨格差异已测；目标加载的全材料类模型误差上界未取得。此处未做 96² 主动复拟合。",
            "[原生和更细网格数据与支出](native/report.json)",""]
    return "\n".join(lines)


def supplement_md():
    p=HERE/"loop"/"fresh_bg_b"
    out=read(p/"dual_cie_verified.json")
    long_fit=read(p/"long_gn_verified.json")
    budget=read(p/"budget_witness_verified.json")
    lines=["## 对同一读数加深求解后，B 和补偿引导仍然留有歧义","",
       "在独立连续背景 fresh_bg_b 上，对原 B 和补偿引导开展同数据加强试验。只读每布局自己的 A 与新 Y。补测为 0。先按每组接收端口独立改写 Dual-CIE 的电流辅助空间，保留发布的三个 Fourier 阶段、β、像素读出及站稳停止规则，再从辅助结果和原全波候选作更长的物理 GN。两次 Dual-CIE 加 GN 约各 51 秒 CPU。严格审计确认 B 和补偿引导的零 ROI、连通斑块、低均值、高均值及 unrestricted 五个家族各有相容端点。","",
       "| 数据合同 | 零区材料残差上界，旧 / 新 | 高均值材料残差上界，旧 / 新 | 最大已验证 ROI 均值跨度 | 该对外区面积加权 L2 |",
       "|---|---|---|---:|---:|"]
    for arm,row in out["arms"].items():
        name="交错 B" if arm=="B_interleaved" else "补偿引导"
        zero,high=(row["families"][k] for k in ("zero","mean_high"))
        ratio=lambda q:f"{q['original_A']['upper_over_r']:.3f} / {q['new_layout']['upper_over_r']:.3f}"
        lines.append(f"| {name} | {ratio(zero)} | {ratio(high)} | {row['verified_inner_width']:.3f} | {row['extreme_verified_pair']['exterior_L2_area']:.3f} |")
    lines+=["","增加到 220 步并校正归一化权重的单纯物理 GN 在补偿引导臂也找回零区和高均值端点，严格宽度为 2.168。B 的同预算纯 GN 在 unrestricted 对照仍高于噪声半径。CIE 起点使 B 和补偿臂都达到双球内，首轮的全波失败属于数值求解失败。","",
        "以各臂重新求出的零区材料为中心，再用预登记的较紧外区预算 0.35 检验 38 格、0.1484375 平方米的 ROI。两种布局均找到 ROI 常量和形状自由两类备选。各备选的 ROI 均值相对各自中心严格升高至少 0.5000009999998。材料在 [0,4.5]，外区面积加权 L2 严格上界均小于 0.349999965001。以下列出目标均值处的双球残差上界。","",
        "| 数据合同 | 固定常量 ROI，旧 / 新 | 形状自由均值，旧 / 新 |",
        "|---|---|---|"]
    for arm,row in budget["layouts"].items():
        a,b=(row["families"][key][-1] for key in ("fixed_constant_ROI","shape_free_ROI_mean"))
        name="交错 B" if arm=="B_interleaved" else "补偿引导"
        ratio=lambda q:f"{q['old_residual_upper_over_r']:.4f} / {q['new_residual_upper_over_r']:.4f}"
        lines.append(f"| {name} | {ratio(a)} | {ratio(b)} |")
    lines+=["","因此这两个布局对 fresh_bg_b 的同一预登记任务均未迫使所需外区变化超过 0.35。界限对各自从数据构造的零区参考点成立。严格回放验证了有限个材料端点和中间任务幅度，连续路径和最小外区代价仍待证明。其余对象及布局的宽类保持未决。","",
       "[双球和预算回执](loop/fresh_bg_b/budget_witness_verified.json) · [Dual-CIE 逐点严格回执](loop/fresh_bg_b/dual_cie_verified.json) · [纯 GN 严格回执](loop/fresh_bg_b/long_gn_verified.json) · [同色标材料图](figures/fresh_bg_b_budget_witnesses.png)",""]
    return lines


if __name__=="__main__":
    cases=compact()
    report=HERE/"REPORT_ZH.md"
    txt=report.read_text()
    prior,heading=txt.split("## 复现",1)
    if "## 主动试验结果" in prior:
        prior=prior.split("## 主动试验结果",1)[0]
    report.write_text(prior+md(cases)+"\n## 复现"+heading)
