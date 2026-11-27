"""Format the completed local search, with raw-input and exact core checks."""
from pathlib import Path
from fractions import Fraction
import json
import numpy as np

HERE=Path(__file__).resolve().parent


def display(x):
    return f"{np.ceil(x*1e7)/1e7:.7f}"


def main():
    d=json.loads((HERE/"results.json").read_text())
    lines=[]
    certificates={}
    for key,z in d["outcomes"].items():
        with np.load(HERE/f"{key}_candidate.npz",allow_pickle=False) as q:
            x=q["material_opposite_class"]
            y=q["material_public_same_class"]
            roi=q["roi_indices"]
            assert np.all((x>=0)&(x<=4.5)) and sum(map(Fraction,x[roi]))==Fraction(z["core_opposite_class"])*len(roi)
            assert sum(map(Fraction,y[roi]))==Fraction(z["core_same_class"])*len(roi)
            assert q["observed_old_400"].shape==(8,8) and q["observed_added_scalar_460"].shape==(1,)
        for k in ("opposite","same"):
            record=z[f"independent_strict_{k}"]
            residuals=[record[str(b)]["residual_upper_over_r"] for b in (400,460)]
            status=("两球内" if max(residuals)<1 else "该候选未通过，材料类未决")
            lines.append(f"| {key} | {'高' if z['initial_sealed_decision']=='high' else '低'} |"
                         f" {int(z[f'core_{k}_class'])} | {display(residuals[0])} | {display(residuals[1])} | {status} |")
        certificates[key]=z["certified_possible_mean_diameter_lower"]
    positive=[k for k,d in certificates.items() if d is not None]
    if positive:
        msg="同一组旧矩阵与追加单复值允许精确核心均值 0 和 2。该结论在 "+ "、".join(positive)+" 上严格通过，宽材料盒的可行均值直径至少为 2。"
    else:
        msg="在两件已采样品上，针对反向核心的有限搜索均未取得同时相容的材料。宽材料盒下的反向核心仍未决。"
    verified=HERE/"validation.json"
    note=""
    if verified.is_file():
        j=json.loads(verified.read_text())["point_residual_intervals_over_r"]
        a=j["unit_415"]["opposite_400"][0]
        b=j["unit_415"]["opposite_460"][1]
        note=f"unit_415 的反向候选已进入追加单标量球，残差上界为 {display(b)}r。该点的旧矩阵残差下界为 {np.floor(a*1e7)/1e7:.7f}r，仍超过旧球半径。 [点态双侧区间](validation.json)"
    text=rf"""# 96² 已采单标量的宽先验再攻击

{msg}

只使用 n96/controls.npz 中原 400 MHz 总波矩阵与冻结权重，以及 n96/new_readings.npz 中每件样品实际获得的一个 460 MHz 复数。未读取任何离线样品或未采端口。旧读数为 8×8，各列幅度缩放 \(1/\sqrt 8\)。460 MHz 使用原 fluence_loop.readout，即开放端口的总传递乘固定 12 位收发权重，再乘冻结的幅度 0.375。归一化可用能量上界为 0.140625。两频率分开使用绝对半径 \(r=0.0016152302843858903\sqrt{{\mathrm W}}\)。

两频率使用同一份实材料，宽盒为 \([0,4.5]\)，核心为 164 格、\(0.0711805556\,\mathrm{{m^2}}\)。在已判规格的反向固定核心上，从对应公开规格开始，做一次带盒投影的局部 Gauss–Newton 搜索。搜索覆盖公开支撑与当下导数选出的额外外部格。400 MHz 使用 128 个实数方程，460 MHz 的 8×8 导数严格投影为**一个复数对应的两个实数方程**。最终每一候选按独立裸 VIE 与 25 端电路的开路公式重算。

| 样品 | 封存的校准规格判定 | 固定核心均值 | 400 MHz 残差严格上界 / \(r\) | 460 MHz 残差严格上界 / \(r\) | 状态 |
|---|---|---:|---:|---:|---|
"""+"\n".join(lines)+rf"""

{note}

两样品的搜索和独立验证合计 {d["cpu_s"]:.3f} 秒单线程 CPU，墙钟 {d["wall_s"]:.3f} 秒。逐样品保留数值试探的 LU、导数右端与独立电路复算。未过门槛的点保留为候选，整个连续材料类保持未决。这个试验仅解释本次已采的单标量。两校准邻域原有的单标量证书仍按原合同解释。

[完整收据](results.json) · [unit_415 材料与实际波](unit_415_candidate.npz) · [unit_862 材料与实际波](unit_862_candidate.npz) · [复现入口](run.py)
"""
    (HERE/"REPORT_ZH.md").write_text(text)


if __name__=="__main__":
    main()
