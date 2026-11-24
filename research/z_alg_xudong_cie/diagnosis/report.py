#!/usr/bin/env python3
"""Render the measured diagnostic trajectories and common-scale material maps."""

import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from diagnose import CASES, HERE, Twin, read_npz
from score import read, trace, TRAJECTORIES

LABELS = ("Austria", "Slotted block", "SRAP 004", "Fresnel Twin")
COLORS = ("#1d659a", "#b94d23", "#28745a", "#764293")
PLOTPOINT = (16, 64)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(HERE/f"{name}.{ext}", dpi=145, bbox_inches="tight")
    plt.close(fig)


def trajectories(summary):
    fig, axs = plt.subplots(4, 3, figsize=(11.9, 10.7), layout="constrained")
    for name, label, (a,b,c) in zip(CASES, LABELS, axs):
        obj = summary[name]
        fold = TRAJECTORIES/name
        _, tr = trace(fold)
        x = np.array([r["step"] for r in tr])
        stages = [r["step"] for i,r in enumerate(tr) if i and r["stage"] != tr[i-1]["stage"]]
        for ax in (a,b,c):
            for q in stages:
                ax.axvline(q, color="#b1bac1", lw=.8, ls=":")
            ax.spines[["top","right"]].set_visible(False)
            ax.grid(alpha=.15)
            ax.set_axisbelow(True)
        for key, title, color in (("optical_diag_pair_normalized_rms","Diagonal",COLORS[0]),
                                  ("optical_offdiag_pair_normalized_rms","Off-diagonal",COLORS[1])):
            a.semilogy(x, [r[key] for r in tr], color=color, label=title)
        a.set_ylabel(label)
        b.semilogy(x,[r["aux_data_relative"] for r in tr],color=COLORS[0],label="Auxiliary")
        audits = obj["train_audits"]
        b.semilogy([r["step"] for r in audits],[r["fit_relative"] for r in audits],
                   "o-",ms=2.7,color=COLORS[1],label="Physical, train")
        if name == Twin:
            held = obj["held_source_audits"]
            b.semilogy([int(k) for k in held],[v["relative"] for v in held.values()],
                       "o--",ms=2.3,color=COLORS[2],label="Physical, held transmitters")
            b.legend(frameon=False, fontsize=7)
        next_y = np.maximum([r["sector_positive_Fro"] for r in tr], 1.e-13)
        c.semilogy(x,next_y,color=COLORS[2],label="Global positive part")
        c.semilogy(x,[r["sector_pixel_positive_trace"] for r in tr],
                   color=COLORS[3],label="Sum of pixel positive parts")
    for ax, label in zip(axs[0],("Optical entry-normalized RMS",
                                 "Relative complex port residual",
                                 "Material sector violation")):
        ax.set_title(label,fontsize=11)
        ax.legend(frameon=False,fontsize=7)
    for ax in axs[-1]:
        ax.set_xlabel("Actual Dual-CIE updates")
    axs[0,2].set_ylim(bottom=1.e-13)
    save(fig,"dual_cie_trajectories")


def maps(summary):
    fig, axs = plt.subplots(4,5,figsize=(11.5,8.5),layout="constrained")
    for name,label,axes in zip(CASES,LABELS,axs):
        here = TRAJECTORIES/name
        receipt = read(here/"receipt.json")
        selects = {row["step"]: row["path"] for row in receipt["material_snapshots"]}
        values = (16,64,min(r["step"] for r in receipt["material_snapshots"] if r["stage"]==2),
                  receipt["final_step"])
        for ax,step in zip(axes[:4],values):
            raw = read_npz(here/selects[step])
            n = int(np.sqrt(raw["m"].size))
            im = ax.imshow(raw["m"].reshape(n,n),origin="lower",vmin=0.,vmax=6.,cmap="viridis")
            ax.set_title(f"Step {step} · peak {raw['m'].max():.2g}",fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
        axes[0].set_ylabel(label)
        if name == Twin:
            axes[4].set_axis_off()
            axes[4].text(.5,.5,"Held-source\nscoring only",ha="center",va="center",
                         transform=axes[4].transAxes,fontsize=8)
        else:
            truth = read_npz(here/"offline_reference.npz")["m"]
            axes[4].imshow(truth.reshape(n,n),origin="lower",vmin=0.,vmax=6.,cmap="viridis")
            axes[4].set_title("Coarse reference",fontsize=8)
            axes[4].set_xticks([])
            axes[4].set_yticks([])
    fig.colorbar(im,ax=axs.ravel().tolist(),shrink=.75,label="Dielectric contrast (display saturates at 6)")
    save(fig,"dual_cie_materials")


def power_maps(summary):
    fig,axs = plt.subplots(4,3,figsize=(8.1,9),layout="constrained")
    for name,label,axes in zip(CASES,LABELS,axs):
        here = TRAJECTORIES/name
        raw = read_npz(here/"step_0064.npz")
        n = int(np.sqrt(raw["m"].size))
        for ax,key,cmap in zip(axes,("current_energy","optical_pixel_trace","optical_pixel_offdiag_fro"),
                               ("magma","coolwarm","magma")):
            data = raw[key].reshape(n,n)
            ceiling = float(np.max(abs(data)))
            low = -ceiling if cmap=="coolwarm" else 0.
            shown = ax.imshow(data,origin="lower",cmap=cmap,vmin=low,vmax=ceiling)
            fig.colorbar(shown,ax=ax,shrink=.56)
            ax.set_xticks([])
            ax.set_yticks([])
        axes[0].set_ylabel(label)
    for ax,caption in zip(axs[0],("Current energy","Signed optical trace","Off-diagonal optical norm")):
        ax.set_title(caption,fontsize=9)
    fig.suptitle("Step 64, unscaled per-cell values")
    save(fig,"dual_cie_power_maps")


def brief(score):
    a = score["cases"]
    output = [
        "# 四条原 Dual-CIE 轨迹的电流诊断",
        "",
        "三个合成实例的辅助电流可以把端口残差维持在 1.0%–5.0%，而由同一步材料完整求解的"
        "终点端口残差为 53.6%–84.4%。Fresnel Twin 的辅助残差为 5.71%，实际重算训练端口为"
        "25.19%，未用的九个入射源为 26.30%。这两个残差由各自的实测有效通道计算。",
        "",
        "| 对象 | 真实更新数 | 跨源/对角原始矩阵范数中位比 | 每元素成对归一 RMS 中位比 | 物理训练残差 | 未用入射残差 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name, label in zip(CASES,LABELS):
        q = a[name]
        held = (f"{q['quality_final']['physical_held_source_rel']:.5f}" if name == Twin else "—")
        output.append(f"| {label} | {q['final']['step']} | {q['median_raw_offdiag_to_diag']:.3f} | "
                      f"{q['median_pairwise_normalized_offdiag_to_diag']:.3f} | "
                      f"{q['quality_final']['physical_train_rel']:.5f} | {held} |")
    output += [
        "",
        "合成场的非对角元素共有 90 项，对角有 10 项。把平方和摊到每个元素，"
        "并按照两列电流与总场的能量归一，合成非对角强度与对角接近，实测中则明显更低。"
        "原始范数反映全矩阵中累计违背的强度；两个比值分别回答不同问题。",
        "",
        "![逐步电流与实际端口](dual_cie_trajectories.png)",
        "",
        "| 合成对象 | 终点全图相对 L2 | 目标内相对 L2 | 背景 RMSE | 半峰 IoU |",
        "|---|---:|---:|---:|---:|",
    ]
    for name,label in zip(CASES[:3],LABELS[:3]):
        q = a[name]["quality_final"]
        output.append(f"| {label} | {q['relative_material_L2']:.3f} | {q['inside_relative_L2']:.3f} | "
                      f"{q['background_RMSE']:.3f} | {q['half_peak_IoU']:.3f} |")
    output += [
        "",
        "非对角成对归一 RMS 对全图材料误差的 Spearman 相关依次为 "
        + "、".join(f"{a[name]['correlations']['all']['optical_offdiag_pair_normalized_rms']['relative_material_L2']['rho']:+.3f}"
                   for name in CASES[:3]) + "。它与实际训练端口残差的相关则为 "
        + "、".join(f"{a[name]['correlations']['all']['optical_offdiag_pair_normalized_rms']['physical_train_rel']['rho']:+.3f}"
                   for name in CASES) + "。各轨迹共审计了 18、15、16、14 个非初值步骤。"
        "同一物体内的更新彼此相关，材料与端口结果的方向差异限制了这些统计量的单独解释。",
        "",
        "![共用色标的材料图](dual_cie_materials.png)",
        "",
        "![逐格原始功率](dual_cie_power_maps.png)",
        "",
        "## 模型、复核和费用",
        "",
        "Austria 对应 150→32 交叉网格观测，开槽和 004 对应不同的 96→32 带噪读数。"
        "所有合成正演审计使用各自原反演网格与同一端口。Fresnel 在 60² 离散模型上使用各源原有效接收索引，"
        "另以原保存的九个保留入射源和 2,169 个有限复数读数评分。生成参考仅在所有线上轨迹封存后打开。",
        "",
        "测量用成熟源逐束更新，保持各自 DPCS、Fourier 模式和材料读出。"
        "每步从不改写的生成器读取实时 J、E、材料和阶段，记录带符号的逐源光学对角、Gram 光谱、"
        "局部和全局扇区的正部、辅助电流端口残差与本构缺陷。原合同的 b 取倒数得到公式中的 β；"
        "当前 PR 用首阶段的 b，材料读出使用各阶段的 b。终点的归一光学矩阵在四案分别有 10、10、10、9 个非零方向。"
        "完整物理电流的归一非对角光学违背最大为合成 3.23e−16、实测 4.26e−9。",
        "",
        "实际更新耗时与指标计算、快照保存和全波审计分开。四案成熟引擎墙钟分别为 "
        + "、".join(f"{a[name]['optimizer_wall']:.2f}" for name in CASES) + " 秒。"
        "每案的 receipt.json 给出模型装载、本机总费用、复右端与 Green 作用次数。"
        "所有三段均走完。光学矩阵的每个元素和空间功率可用快照直接复算。"
        "早期历史迭代交叉核对见回执；跨机末段轨迹仍可能积累小的数值差异。",
        "",
        "[逐案机器可读评分](scores.json) · [在线封存](online_seal.json) · "
        "[轨迹、快照和费用](trajectories/) · [复算入口](REPRODUCE.md)",
    ]
    (HERE/"REPORT_ZH.md").write_text("\n".join(output)+"\n")


def main():
    scores = read(HERE/"scores.json")
    trajectories(scores["cases"])
    maps(scores["cases"])
    power_maps(scores["cases"])
    brief(scores)


if __name__ == "__main__":
    main()
