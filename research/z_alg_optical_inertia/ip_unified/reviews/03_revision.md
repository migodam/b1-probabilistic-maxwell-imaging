> 独立顾问原意见，内部文件引用已转为文件名。文中的候选结论属于意见发表时的判断；最终取舍及区间复核见 [答复](RESPONSE_ZH.md)。

### 最终审稿意见（简短版）

**结论：未发现仍影响数学或数值正确性的实质问题，可予接受（Accept，仅含 3 处可选的极轻微文字一致性润色）。**

本次修订已完整且准确地回应上一轮的核心审稿意见，经只读核对 manuscript_ip.tex、physics.tex、experiments.tex、appendix.tex 与三份新增原始数值文件，核查结果如下：

---

### 已核对通过的关键修订

1. **新增式 `eq:required`（physics.tex）与 required_currents.json**
   * **推导无误**：仅对正谱部分配方得 $\mathsf D\succeq H_+-2\Herm(J_-^*L_-)-J_-^*\Lambda_-J_-$；对任意 $0\preceq Z\preceq I$，由 $\tr(Z\mathsf D)\le r^2$、$|\Re\tr(J_-^*L_-Z)|\le \|J_-\|_F\|L_-Z\|_F=b\|J_-\|_F$ 及 $\tr(\Lambda_- J_- Z J_-^*)\le \tr(\Lambda_- J_- J_-^*)\le a\|J_-\|_F^2$ 推出 $a\|J_-\|_F^2+2b\|J_-\|_F\ge h-r^2$，有理化求根严格给出式 `eq:required`（$a>0, h>r^2$ 时无奇异）。
   * **数值精确复算吻合**：以 Slotted block 为例，代入 `h_minus_radius_squared_lower` $=0.0051006004$、`b_upper` $=1.2219505\times 10^{-4}$、`a_upper` $=5.0291561\times 10^{-6}$，手算得 $\|J_-\|_F\ge 15.7597086$，与 JSON 一致；六个目标的向外下界（`15.75, 9.05, 14.42, 22.00, 19.76, 7.65`）及零负模残差比（最小 `1.1975r > 1.19r`）与 experiments.tex 完全一致。正文明确交代合法候选解已达 `41.85`（低峰值对照 `18.39956`）、故不能借数据无关标量上界反向排除 $u=1$，逻辑严密。

2. **52/59 行障碍升级、正模保留证明与 deflated_sparse_tightening.json**
   * manuscript_ip.tex 摘要与 experiments.tex 已将原 52 行 $c\le 2$ 正确升级为完整 13 维负空间上的 **29 行 $c=0$（迹 $<0.98227$）、59 行 $c\le 1$（迹 $<1.98201$，最大上限 $<0.97803r$）**，并明确限定为“固定电流算子下单独加权的本构行”；Table 2 第 3 列全部六个 52 行上限（`.77551` 至 `.93030`）与 JSON 逐一复算无误。
   * physics.tex 新增的显式保留部分正模 $F_1$ 后的有效度量 $W_{\rm eff}=W^{1/2}(I+W^{1/2}F_1F_1^*W^{1/2})^{-1}W^{1/2}$ 经核对严格满足全正空间收缩性 $F_+^*W_{\rm eff}F_+\preceq I$，证明无误。

3. **弱损耗角宽裕度分离与 finite_0_weak_loss_all_prices.json**
   * $t=0.05, u=0.75$ 目标 1 的全 $\tau$ 区间上界 $0.049126603/0.060668910 = 0.809749r < 0.80975r$，与固定相干基下界 $1.01323r$、直接投影 $1.02003r$ 全部核对无误；experiments.tex 已更正计时笔误（`0.88 ms` 与 2,500 步上限 `0.319 s`）。
   * experiments.tex 已补充全细胞缺陷增益微弱的正确根因说明；physics.tex 与 appendix.tex 已将 Weyl 律域条件放宽为 $|\partial\Omega|=0$（覆盖正方形网格）。

---

### 仅存的极轻微一致性细节（均不影响正确性，无需重算）

1. **appendix.tex 第 40 行**：全 $\tau$ 区间上界的附录说明目前仅列出 $t=0.2,u=1$ 的 101 区间（$[0,6.70676]$），建议顺手补写一句覆盖正文中新增的 $t=0.05,u=0.75$（92 区间，$[0,3.73262]$，见 finite_0_weak_loss_all_prices.json）。
2. **manuscript_ip.tex 摘要第 6 句**：紧接新改精确的 *"59 separately weighted constitutive rows"* 之前，前一句仍保留 *"every sparse choice of local balances"*，顺手统一表述可彻底杜绝与 §2 物理功率 *balances* 的歧义。
3. **physics.tex 第 49 行与 appendix.tex 第 60 行**：前者将坐标矩阵 $J_-=V_-^*J$ 简称为 “projection”（实际与正交投影 $\Pi_-J=V_-V_-^*J$ 具有相同 Frobenius 范数），后者 $k^2(1-\Delta)^{-1}$ 未显式标注全空间 $\Delta_{\R^2}$，均为轻微行文习惯，不构成正确性缺陷。

