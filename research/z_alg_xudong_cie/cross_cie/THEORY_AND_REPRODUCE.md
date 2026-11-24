# 在原 Fourier 电流空间内修正光学约束

设总场 $E=B+GJ$。把每个源的两列 $J,E$ 同除以该源入射列的范数，所得矩阵为 $J_w,E_w$。在区域投影 $P$ 内，损失使用

$$
Q_o=\operatorname{Im}_H(J_w^*PE_w),\qquad
Q_r=J_w^*PJ_w-\operatorname{Re}_H(J_w^*PUE_w),
$$

其中 $U$ 是逐点材料上限。将 $Q_o$ 除以当前区域的总场能量。将 $Q_r$ 除以 $(\|PJ_w\|_F+\|PUE_w\|_F)^2$。沿单次线搜索，两项分母冻结。光学可只保留原源对角、保留全块，或在有损符号下只保留正部。扇区使用正部。

若采用实内积梯度并令损失为 $\|Q_o\|_F^2/2$，其未归一梯度为

$$
\nabla_J=-i\big(PE_wQ_o-G^*(PJ_wQ_o)\big)W.
$$

$W=\operatorname{diag}(1/\|B_{:s}\|_2)$。将 $Q_o$ 换成被保留的投影即可得到对应损失的梯度。扇区正部的梯度为

$$
\big[(2PJ_w-PUE_w)Q_{r,+}-G^*(UPJ_wQ_{r,+})\big]W.
$$

所有区域的共轭传播右端先求和，最后只施加一次 Green 伴随。当前的原 CIE 参数 $p=b/(1+bm)$ 给出尺度 $\|p(J-mE)W\|_F^2/\|(J-mE)W\|_F^2$。光学和扇区统一乘以这个尺度，再与原源码的端口及状态项一起做单实步长搜索。修正方向投影到相同的 Fourier 系数与 DPCS 余空间。原材料在这一个线搜索中固定，随后由原读出更新。每次候选线上的物理二次型用系数直接求得，向前裁剪电流步，零步始终可用。

## 逐格光学

对一格的加权行向量，定义 $a=\|j\|^2$、$b=\|e\|^2$、$c=j^*e$，其中 $c$ 表示行元素的共轭内积。无损损失是 $[ab-\operatorname{Re}(c^2)]/[2(b+\delta)]$。$\delta$ 在一次修正内固定为当前平均行总场能量的 $10^{-12}$ 倍，并设置浮点防零保护。$b>0$ 时，消去任意实标量的最小误差 $r_\mathrm{re}^2=a-(\operatorname{Re}c)^2/b$。逐格损失在 $\tfrac12r_\mathrm{re}^2 b/(b+\delta)$ 与 $r_\mathrm{re}^2 b/(b+\delta)$ 之间。这解释了为什么它直接约束共用实材料，但未包含材料盒的额外限制。

被动测量臂采用正特征值平方除以 $b+\delta$。唯一可能为正的特征值为 $[\operatorname{Im}c+\sqrt{ab-(\operatorname{Re}c)^2}]/2$ 的正部。两条梯度都由标量式得到，避免构造每格源矩阵或额外稠密分解。

## 停止、开发和复现

`cross_engine.iterate(model, Y, case, method, m0, upper, ...)` 每步返回原 `(material, step, stage, event)`。`upper` 支持标量或逐点向量，并在每次原材料读出传入约束。`model.clock` 计入修正及额外 FFT。对照 `reset` 每 16 步清空原 PR 历史。冻结的修正臂也执行相同规则，并在原阶段拟提前结束时检查光学步的系数相对位移平方。两种位移须同时通过阶段阈值，每阶段原 500 次上限及外层共同墙钟仍适用。固定设置见 [FROZEN.json](FROZEN.json)。

从本目录使用 SciPy、NumPy 和 threadpoolctl。独立数据的统一裁决不由本目录开发脚本执行。

```sh
python -B check.py
python -B develop.py --arm dual --arm diag16 --arm full_ot16 --arm full16 --arm regional4_4 --arm regional16_1 --arm pixel025 --arm pixel1
python -B develop.py --arm reset --arm diag16_guard --arm full16_guard --arm regional16_1_guard --arm pixel025_guard --summary development_refinement_scores.json
python -B develop.py --wall 30 --case fresnel_twin_8GHz_half180_tx --arm dual --arm reset --arm diag16_guard --arm full16_guard --arm pixel025_guard --summary development_twin_scores.json
python -B develop.py --wall 15 --arm pixel025_accept_only --arm dual_start065 --summary portfolio_synthetic_branches.json
python -B develop.py --wall 24 --case fresnel_twin_8GHz_half180_tx --arm pixel025_accept_only --arm dual_start065 --summary portfolio_twin_branches.json
python -B portfolio.py
```

[实内积梯度与原算法校验](checks.json) · [基本开发组](development_scores.json) · [匹配重启及退出检验](development_refinement_scores.json) · [Twin 的保留入射端口](development_twin_scores.json)
