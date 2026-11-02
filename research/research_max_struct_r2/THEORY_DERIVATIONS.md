# R2：辅助状态的度量、表示误差与可观测商

状态：以下有限矩阵命题为 PROVED_UNDER_ASSUMPTIONS；经典工具的 CIE 专门化不主张首次。数值收益必须由封存后评分另判。连续 Maxwell 同场族沿用 R1 的适定性前提。R1、Atlas、FBG-R3 的恒等式不重复做精度试验。

## 1. 消去 current 究竟做了什么

固定材料 m，令 L=I−diag(m)G 可逆，j*=L⁻¹diag(m)E，H=P(m)L，P_i=β/(1+βm_i)。为对应本轮总能量标量归一化，明确记 Sbar=S/||Y||_F、ybar=y/||Y||_F、Hbar=P(m)L/||E||_F，要求 ||Y||_F、||E||_F 非零；单源使用同一固定标量。以下本节中的 S、y、H 均为这些归一化后的算子和数据的简写。本段写一个源，多个共同度量源逐列相加。

考虑 Ψ(m,j)=½||S j−y||²+½||H(j−j*)||²。写 z=H(j−j*)，r=Sj*−y，K=SH⁻¹，则内层为 ½||r+Kz||²+½||z||²。正规方程给 z*=−(I+K*K)⁻¹K*r。代回并用 Woodbury，得到

    min_j Ψ = ½ r* (I+KK*)⁻¹ r.

这是 profile 目标，**不是边缘似然**；概率积分还需要相应归一化/行列式项和定义良好的潜变量先验。接近大接收 resolvent 的材料，其物理残差受到更强压低。因此精确 current 优化可以更准确地最小化一个与最终物理目标不同的目标；更彻底的内层并不必然更接近真实材料。这一形式已有 [van Leeuwen–Yang，arXiv2403.15292v1，§2.4 Corollary 1 / 式(10)](https://arxiv.org/html/2403.15292v1) 的 variable metric reduced formulation；更早的精确状态消元与 Schur 导数见 [van Leeuwen–Herrmann，Inverse Problems 32,015007，§4](https://doi.org/10.1088/0266-5611/32/1/015007)。不能作为新数学首创。

若 κ=||K||₂，则 ||r||²/[2(1+κ²)] ≤ Ψ_profile ≤ ||r||²/2。没有材料域上统一的 κ，辅助小值不能提供统一物理残差上界。R3 已有任意辅助 current 的精确残差修正式；这里把它与消去后的度量联系起来，不再重新证明旧式。

一个决定性代数反例：一维 K=κ、r=κ 时，profile=κ²/[2(1+κ²)]≤½，而物理残差κ任意大。这否定无 resolvent 假设的普遍推论，但这个抽象标量反例本身不是一套可辐射 Maxwell 材料。实际物理效果由本轮真实 Hankel–Green 的端点诊断给出。

## 2. Fourier/APCS 缺口与内层未收敛是两个可分开的量

固定 m、共同目标，A=S*S+H*H 正定。J_full 为无约束 current 极小点，允许的 current 是仿射族 J_D+U a。由于 full 解的一阶项为零，配方得到

    Ψ(J_D+Ua)=Ψ(J_full)+½||J_D+Ua−J_full||²_A.

所以 affine 最小值与 full 最小值的差恰为到该仿射族的 A 距离平方的一半；当前 iterate 与 affine 最小值之差是尚未优化掉的部分。嵌套增加 Fourier 空间只能使 profile 值不增，**不能保证材料图更好**。A 中包含传播与接收，普通 current 的 L² SVD 尾不是此缺口。固定 DPCS 偏移也必须保留，不可用不同 affine anchor 作同容量比较。这是加权正交投影/Pythagoras 的应用；[Osborne，ETNA 28 (2007)，pp2–3，式(1.3)–(1.5)](https://etna.ricam.oeaw.ac.at/vol.28.2007-2008/pp1-15.dir/pp1-15.pdf) 已有变量投影依据。

当不同源的 Ωd、Ωs 不成比例，不能假装一个 Cholesky 为全部源服务：正常方程是 S*S J Ωd+H*H J Ωs=B。可用广义特征分解同时合同对角化两个 SPD 度量，随后每个广义源独立求解不同的加权系统。R2 的 full-profile 控制明确选总能量标量归一化，才可对全部源复用一次 factor；它与原生 source-normalized CIE 的辅助目标不同，必须连同额外容量和费用一起标出。

## 3. 新实现需要的导数，而非冻结公式重测

R2 full-current 控制使用 d=SJ−Y、Q=E+GJ、e=J−mQ、w=β/(1+βm)，内部目标 ||d||²/(2||Y||²)+||we||²/(2||E||²)+R(m)。在 current 精确极小点，包络定理消除 dJ/dm；实材料导数为

    g_i = Re Σ_s conj(w_i e_is) [−w_i² e_is−w_i Q_is] / ||E||² + ∂_i R.

−w²e 项在非物理 current 不可省略。新接线只做一次真实4×4模型上的方向差分，计入三次 current 分解；不扩成精度扫参。对 current 的不精确解，包络导数不再精确，需内层误差控制；相关 partial-minimization 与不精确内层分析见 [Aravkin–Drusvyatskiy–van Leeuwen，arXiv1606.02395v3](https://arxiv.org/abs/1606.02395)。本轮没有隐藏这种近似。

## 4. 等变性是完整迭代的性质，不是性能定理

源酉换基 E→EU、Y→YU。Ω→U*ΩU 后目标不变；global Frobenius 梯度、PR 内积、线搜索、材料共享 readout、初始化与 stage zero-padding 都相应运输。只要停止也使用 ||ΔA||²_F/||A||²_F，整个有限轨迹和停止事件等变。原实现即使 sourcewise_cg=False，逐列相对变化均值也不等变。R2 修补的是完整 full-metric/停止合同，不是发明 global CG。

本轮成像强制同样1500步，故停止规则不会成为胜负混杂因素；full-metric 接线另在小模型核对。原始列归一化、sourcewise PR、global PR 都可各有有利轨迹。等变性本身不推出更大 basin、较少工作或正确 GT。

## 5. 新增照明作用于候选商，而不是优化器的奇迹

R1 的旧源同场关系定义合法候选集合中的等价类。旧源任意线性组合与任意新增接收对该类仍恒定；换源基不能打破它。对完整候选 transfer H_a,H_b，额外源c的距离 d_ab(c)=||W(H_a−H_b)c||。单pair、功率 c*Rc≤1 的最大值是 ||W(H_a−H_b)R⁻¹/²||₂。这是标准奇异向量设计。

**有限候选判别命题。** 已知真值属于有限库，噪声半径ε。若所选采集对所有pair有 d_ab>2ε，则任何一次观测至多有一个合格候选。证明：两个候选都在观测ε球内则由三角不等式 pair≤2ε，矛盾。若真值只保证在某候选的预测κ邻域内，则保留候选的半径必须扩为ε+κ；同样充分条件改成 pair>2(ε+κ)。没有覆盖界κ，候选内分离不能保证覆盖真值。

R2实际三维实验把先验明确收紧为 ε=I+t b(x)zzᵀ、t∈[0,3]，b是已知平滑紧支撑。两个沿z的横向旧源对整个区间严格无信息。新源按公开端点完整 transfer 选择；独立生成进程提供一份真实新增的有限DDA观测；拟合进程只读该观测，求一个t。该实验为一个公开族内实例提供新增观测后的相容标量估计；不证明连续族全局可辨识性、拟合唯一性或噪声下的恢复保证。它**不能**证明自由体素张量成像超过CIE，也不保证端点最优等于整个连续族的minimax最优。

已知轴的族还有精确消元：只有 j_z 非零，令 D=diag(t b)，则 j_z=(I−D G_zz)⁻¹ D Einc_z，外场仍由完整 S_z 外传。这是受限本构的精确块化，不是弱散射；本轮为复用同一已资格后端，仍计费完整3N DDA。将此缩减称为通用Maxwell标量化是错误的。

## 6. 当前论文价值判断

变量投影、数据度量、global CG、已知轴消元、候选分离三角界均有标准来源。潜在论文价值来自**被真实证据支持的组合解释或新失效现象**，不是公式改名。若二维盲测未有材料收益，应明确算法贡献不成立；即使胜出也只支持这三个独立对象，不制造统计显著性。

原有 CIE/TSOM 与不同 β 的公开原始方案见 [Zhong–Xu，J.Imaging 5,27 (2019)，§2.2–2.4](https://doi.org/10.3390/jimaging5020027)；测量域 state/data 失配的相邻已知路线见 [Sun 等，CC-CSI (2017)，§3](https://doi.org/10.1109/TAP.2017.2673758)。全文资格及检索失败记录见 [PRIOR_ART](literature/PRIOR_ART.md)。
