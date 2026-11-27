"""Build the Chinese readout from certified receipts."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name):
    p = HERE/name
    return json.loads(p.read_text()) if p.exists() else {}


def digits(x):
    return f"{x:.6g}"


def main():
    summary = load("summary.json")
    b = load("B_independent_trajectory.json")
    held = load("B_independent_held.json")
    ab = load("joint_original_AB.json")
    lines = [
        "# 原生带噪观测允许的有限面积连续材料变化",
        "",
        "两套原始布局 A 都允许把一整块原本归零的区域连续增亮。严格核验覆盖整段材料折线及其中每个未采样的材料配置。另一项固定点证书给出精确等复读数的连续曲线，证明的幅度小得多。原 B 的同一份带噪数据也允许零区和大斑块两个重新拟合的物理解。",
        "",
        "在两张固定的原生标称校准网络上，采用 400 MHz 二维标量 TM、24 个 50 Ω 端口、独立实对比度 0 到 4.5 和原封存的复 Frobenius 球。A 或 B 各使用 8 次照明、64 个复读数，每次 0.125 W，源列可用功率共 1 W。源驻留时间未提供。原测量指保存在发布数据中的、已扣除空场并带噪的复端口观测。两套原生网格为不同对象。",
        "",
        "## 原布局 A 的严格连续区间",
        "",
        "| 原生观测 | 任务 | 面积 m² | 任务均值上达 | 全登记区均值上达 | 严格段 | 最大段残差／半径 | 全程外区 L¹ / L² / L∞ 预算上界 |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    challenge = []
    for case, label, method, title in [
        ("n32_c01","_broad","ball_min_margin","固定连接斑块，每格同幅"),
        ("n32_c01","_whole_region","ball_shape_free","全登记区均值，形状自由"),
        ("n96_c01","_broad","ball_min_margin","固定连接斑块，每格同幅"),
        ("n96_c01","_ball","ball_shape_free","斑块均值，形状自由"),
    ]:
        paths = [p for p in summary.get("paths",[]) if p["case"]==case and
            p["archive"].endswith(label+".json") and p["method"]==method]
        if not paths:
            continue
        row = paths[0]
        costs = " / ".join(digits(row["complete_path_exterior_budget_upper"][k])
                            for k in ("L1_m2","L2_m","Linf"))
        lines.append("| " + " | ".join([case,title,digits(row["area_target_m2"]),
            digits(row["mean_task_lower"]),digits(row["mean_region_lower"]),str(row["segments"]),
            digits(row["residual_upper_over_r"]),costs]) + " |")
        challenge.append(row)
    lines += [
        "",
        "表中每个任务区间从零开始。固定斑块均匀增亮，剩余登记区始终固定为零。形状自由的分支允许任务区内每格独立变化，只将区内均值投影回指定值。外区成本相对于该分支自己的真空区锚点，L¹ 单位为 m²，L² 为 m，L∞ 为无量纲对比度。成本对比覆盖每条折线上所有材料。",
        "",
        "原观测的完整编号为 `bmtm26_n32_c01_f400_a08_c48` 和 `bmtm26_n96_c01_f400_a08_c48`。每条路径的最紧段中心另以原网格全波直接求解，完整复数读数和独立原生回放都随材料节点保存。先前的单像素路径分别把一个像素升至 2.25 和 4.5，对全区均值的贡献只有约 0.01308 和 0.00721。固定面积试验分别包含 43 和 156 像素，占原登记区四分之一。",
        "",
    ]
    for p in challenge:
        required = p.get("remaining_to_declared_targets", p.get("three_count_limits_to_goal"))
        t = p["method_clocks"]
        cost = p["reduced_anchor_setup_clocks"]
        lines.append(f"- `{p['case']}` 的 `{p['archive']}` 分支支付了计算折线 {digits(t['wall_s'])} 秒墙钟、{digits(t['cpu_s'])} 秒 CPU，核验原生支撑预处理 {digits(cost['wall_s'])} 秒墙钟、{digits(cost['cpu_s'])} 秒 CPU。该形状任务均值离 2 尚有 {digits(required['task_mean'])}。")
    restart = load("n32_c01_paths_broad.json").get("checkpoint_restart")
    if restart:
        lines.append(f"- 32² 固定斑块先用尽 140 次尝试，在幅度 {digits(restart['start_at_uniform_patch'])} 封存材料。随后从同一份材料独立重启，另用了 {restart['attempts']} 次尝试，增加 {restart['strict_new_segments']} 条严格段。")
    lines += ["",
        "同一 256 外区像素支持的余量加权与最小 L² 对照均保留，包含每次失败的端点、材料余量、实雅可比谱、端点逆缺陷和原物理交叉项。加大支持用旧数据拟合图的差异辅助选取外区像素。每个候选仍用原始观测独立验真。用尽步进预算或降至最小步长只能确定本次找到的区间。",
        "",
        "[原生结果与费用](summary.json) · [路径带图](certified_paths.png) · [复算与全段证明](REPRODUCE.md) · [预算条带图](phase_bands.png) · [条带原数](phase_bands.json)",
        "",
        "## 精确等读数的固定点曲线",
        "",
    ]
    for case, rec in summary.get("exact_constant_field",{}).items():
        lines.append(f"- `{case}` 的固定连接斑块在精确校准存储系数下保持复读数 F(x₀) 不变，严格幅度为 0 到 {digits(rec['amplitude'])}。收缩上界 {digits(rec['contraction_upper'])}，外区材料余量占用上界 {digits(rec['exterior_margin_fraction_upper'])}。 [{case} 证书]({case}_exact_chart.json)")
    for case in ("n32_c01","n96_c01"):
        shape = load(f"{case}_exact_shape.json")
        if shape and shape["strict"]:
            lines.append(f"- `{case}` 允许任务区形状自由。数据雅可比选出的形状方向在各任务格介于 {digits(shape['direction_on_task_min'])} 和 {digits(shape['direction_on_task_max'])}，整区实际变亮。{shape['target_cells']} 格任务区均值严格连续到 {digits(shape['patch_mean_interval'][1])}，所有原复读数等于锚点。 [均值任务等场证书]({case}_exact_shape.json)")
    lines += ["",
        "降阶计算使用原生锚点、全波支撑反馈和真实端口。端点电流、逆、接收伴随沿每段分别插值，保留材料变化与场变化的二阶交叉项。接收修正后的三次 Bernstein 系数保留复数正负方向，最后才取范数并加向外误差。精确等场管单独包围雅可比提升缺陷、切向缺陷、方向二阶反馈和横向余量。",
        "",
        "## 原 B 和 A+B 重拟合",
        "",
    ]
    if b:
        a = b["independent_endpoint_audits"]
        bv,bh = (a[k]["B"]["residual_upper_over_r"] for k in ("vacuum","high_patch"))
        av,ah = (a[k]["A"]["residual_lower_over_r"] for k in ("vacuum","high_patch"))
        path = b["continuous"]
        cb = b["endpoint_pair_exterior_cost"]
        lines += [
            f"原 B 的同一份封存观测同时允许零区和 43 个连接斑块像素逐格不低于 2。两个严格残差上界分别为 {digits(bv)} 和 {digits(bh)} 倍 B 半径，全区均值为 0 和 {digits(b['endpoint_pair_region_means'][1])}。相对两端点的外区变化 L¹、L²、L∞ 分别为 {digits(cb['L1_m2'])}、{digits(cb['L2_m'])}、{digits(cb['Linf'])}。旧研究的 B 真空搜索停在约 4.51 倍半径，这一较长的独立重拟合找出了球内材料。",
            "",
            f"同一对材料复查 A 的残差下界分别为 {digits(av)} 和 {digits(ah)} 倍 A 半径。B 的另一条固定斑块折线严格覆盖均匀振幅 0 到 {digits(path['certified_patch_mean_interval'][1])}，包含 {path['strict_segments']} 段，最大段残差上界 {digits(path['accepted_ball_upper_over_r'])} 倍 B 半径。",
            "",
        ]
        if held:
            vals = [digits(held["scores"][k]["residual_lower_over_r"]) for k in ("vacuum","high_patch")]
            lines += [f"先封存 B 的两个材料，再用未参与求解的原 A 发射列与另外 8 个接收端口评分。两份材料的留出残差下界为 {vals[0]} 和 {vals[1]} 倍留出半径。这个独立留出结果只评分这里实际找到的两份材料。 [材料与留出复波形](B_refit_shapes.png) · [留出核验](B_independent_held.json)", ""]
    if ab:
        for name, label in (("vacuum","真空区域"),("high_patch","固定斑块最低 2")):
            vv = ab["selected"].get(name)
            if not vv:
                continue
            s = ("严格可行" if vv["jointly_verified"] else "未找到联合可行解")
            lines.append(f"- 原 A+B 采集的{label}分支{s}，所选候选两球残差上界中的较大者为 {digits(vv['maximum_A_B_residual_upper_over_r'])}。 [独立复核与全部尝试]({ab['sources']['vacuum' if name=='vacuum' else 'high']})")
        if ab.get("heldout"):
            h = ab["heldout"].get("high_patch")
            if h:
                lines.append(f"- 联合高斑块见证封存后在独立接收弧复查，留出残差下界为 {digits(h['residual_lower_over_r'])} 倍半径。 [复算收据](joint_original_AB.json)")
        lines += ["", "A+B 使用原来各自封存的两份观测，支付 16 发、128 个原始复读数和总计 2 W 源列可用功率。两球分别验收。联合真空分支判定按当前数据与计算范围给出，未覆盖材料全类。 [联合收据](joint_original_AB.json)"]
    (HERE/"REPORT_ZH.md").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    main()
