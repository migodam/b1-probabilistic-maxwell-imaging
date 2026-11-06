# R10：分母中的 Green、Two-fold 与可检验的材料更新

本文记录实际实现所用公式。`DERIVED` 表示有限矩阵推导，不是连续 Maxwell 新定理；实验结论另列。时间约定、有限离散 G/E/S 与 R9 一致，实 contrast 为 m=εr−1，固定 μ。所有投影从公开背景算子构造。

## 1. A+ 现有分母究竟包含什么

固定辅助电流 J，令 Q=Einc+GJ，e=J−D_m Q。A/A+ 定义

\[
d_i=\beta-G_{ii},\quad a_i=1+d_i m_i,\quad z_i=m_i/a_i,
\quad u=\beta D_a^{-1}e,
\quad R_A=(I+D_zG_{\mathrm{off,block}})u.
\]

这里 Goff,block 是互不重叠的 2×2 空间块内非对角项，**不是全域 G 去掉对角线**。A+ 的 Möbius/MM 读出降低同一个 A 材料目标的求解费用；并没有把所有多次散射写进这个分母。完整正演/最终验收仍使用原 G，不受这一局部近似替换。

自项 Gii 来自原物理离散，不能因恢复失败就称它“写错”。当前局限是只在分母吸收自作用、局部非对角作用只传播一次。这是待检验的辅助 metric 选择，而不是 Maxwell 方程本身出错。

## 2. 为什么不能把 Two-fold 奇异值直接塞入逐点分母

取背景投影 P=UU*。一般 D_mP≠PD_m；非均匀材料会把包含和遗漏的子空间重新耦合。G 的奇异值又不是其复特征值，更不能代替有相位的 Gii。逐模式写成 m/(1+(β−σ)m) 不等价于原材料算子，通常也没有原物理解释。

本轮实现真正的矩阵分母：选 K=B U*，定义

\[
A_K(m)=I+\beta D_m-D_mK,\qquad R_K(m,J)=\beta A_K(m)^{-1}(J-D_mQ).
\]

所有 N 个实材料自由度保留。取 β>0，只要 A_K 可逆，R_K=0 当且仅当 e=0；但非零 residual 区域的 objective 和梯度改变。这个零集等价性既不保证相同驻点，也不保证恢复正确材料。

候选：

- Two-fold：U 为接收强右子空间与其正交补中的 G 强右子空间之和，阈值均 σ/σmax≥0.1。
- 等秩 Green 对照：只取 G 的前 r 个右奇异向量，r 与 Two-fold 相同。
- 简化 Green 强子空间：只取阈值0.1对应的 G 强右向量。

主要 K=PGP，即 B=U(U*GU)；固定状态另有 K=GP，即 B=GU 的控制。投影仅定义辅助 residual，不把 K 宣称为新的物理 Green、互易材料或被动 Maxwell 求解器。

历史二维公开算子得到 G 强秩19、接收强秩9、Two-fold秩28；这与反演实际 noise/fixed DPCS rank 是两个不同选择。实际 FFT/APCS 空间完全不变。

## 3. 低成本求逆与精确导数（DERIVED）

令 a=1+βm、z=m/a、F=D_zB、T=I−U*F。Woodbury 给出

\[
A_K^{-1}b=D_a^{-1}b+F T^{-1}U^*D_a^{-1}b.
\]

只分解 r×r 复矩阵。不得把该小分解计成免费的操作，也不得与 N×N 原 Maxwell LU 等价计费。每次材料 value/gradient 计一小 LU、正/伴随小系统 RHS 和低秩乘法。材料先验与原算法一致。

记 w=A_K^{-1}e，固定源 metric Ω、ρ>0，

\[
f_K(m)=\tfrac12\rho\beta^2\operatorname{Re}\langle w,w\Omega\rangle+R(m).
\]

对实方向 h，

\[
dw[h]=-A_K^{-1}D_h\{Q+(\beta I-K)w\}.
\]

令 v=A_K^{-*}(wΩ)，则

\[
(\nabla f_K)_i=-\rho\beta^2\operatorname{Re}\sum_s\overline{v_{is}}
[Q+(\beta I-K)w]_{is}+\partial_iR.
\]

这包括分母依赖材料的导数。不能把 A_K(m) 冻住后声称求了同一目标的精确梯度。实现用8次最多的 L-BFGS-B 材料步，并与旧点、普通 CIE 读出种子比较该候选的完整辅助目标；最终整条轨迹仍只由原 Maxwell Φ 选择。

原 A+ 的四颜色独立条件更新依赖互不重叠局部块；这里 K 非局部，因此没有沿用这个独立性假设。Möbius、Woodbury 和投影本身都是已有数学工具，未宣称新恒等式。

## 4. 小 Green 尾并不自动意味着分母近似准确（DERIVED）

若 A_G、A_K 均可逆，则

\[
R_K-R_G=-A_K^{-1}D_m(G-K)R_G.
\]

故误差取决于遗漏 Green 作用、材料乘法和接收前 residual/resolvent 的共同作用。仅比较 ||G−K|| 或 G 奇异值能量不足以证明材料更新更好。高 contrast 增大 D_m；非均匀材料与 P 的不交换性会造成子空间耦合。

这预测两个可检验反例：同样的静态强子空间，对不同材料状态可能产生不同的物理步方向；辅助目标下降也可能伴随 Φ 上升。固定状态试验已记录两类情况，完整结果见 RESULTS/denominator_fixed_state。

## 5. 这不是全局凸化（DERIVED）

即使只有一个实自由度，合法正分母也不能使材料目标普遍凸。取 J=0、Q=1、a=1+m，则

\[
f(m)=m^2/(1+m)^2,\qquad f''(m)=2(1-2m)/(1+m)^4,
\]

在 m>1/2 时为负。此例只否定“分式归一化保证凸性”，不作为具体 Maxwell 成像反例。

同理，换一个对数据凸的 loss（例如某种交叉熵）一般不能消掉非线性 forward 的曲率项：复合目标的 Hessian 除 J*Hloss J 外，仍有 loss 梯度乘 forward 二阶导数。必须根据实际噪声和数据意义选 loss，不能靠名称承诺全局凸化。

## 6. 能量分解与因果边界

本轮能量为离散 Euclidean/Frobenius 平方范数，不是电磁储能。材料误差投到 G/S 右基是坐标诊断，不是材料 Jacobian 可观测性。辅助 J、材料导出的物理 j(m)、两者差 g 分开；eaux=ephys+g 必须保留交叉项。

FFT/APCS 的 affine 距离下界仅约束该辅助电流族，原物理 corrector 不在该族中。Two-fold 顺序投影、数值交集、投影后张成空间与实际 FFT/APCS 空间分别记账。百分比必须连同绝对带能量报告。

同一 current 算法不代表两条轨迹实际 J 相同。固定状态反事实隔离一步材料读出；完整同-current反演隔离整个材料读出规则；等秩 Green 控制检验 Two-fold 的接收信息是否必要。不得将三个层次混为一个因果证明。

## 7. 原始自项保留控制（DERIVED / 已运行DEV）

为排除“新方法只是丢掉原自项”的解释，另设

\[
K_s=PGP+\operatorname{diag}\{G-PGP\}.
\]

它保留每个原始 Gii，同时引入投影后的非局部耦合。此时 Woodbury 对角基为 a=1+[β−diag(G−PGP)]m，通常复数；伴随求逆必须除以 conj(a)，材料导数仍用同一个 K_s。

两个完整历史DEV反演保留原设置：b20d最终NRMSE .04042；40ac .05396。后者corrector前 .16018，完整corrector后达到上述数值。该控制反对“原Gii错误导致失败”的简单解释，但不单独证明是哪组非局部耦合改变了搜索盆地。它未作为第三候选加入新确认集；新确认只检验事先冻结的Two-fold与更简单的Green强子空间两版。
