> This is the original frozen-stage derivation record. For the subsequent full central-data coverage proof, critical-point certificate and constrained source correction, read amendment_codex/THEORY_CERTIFICATES.md. Original UNKNOWN statements describe the first implementation, not the final central-data certificate.

# R5：核查后的数学与物理合同

以下都是离散二维 scalar-TM VIE 的陈述。波数记作 κ，避免与独立未知材料数 k 混淆。n 是完整电磁状态数；k 是坐标选择矩阵 P 的列数；s 是旧源数；D 是复商长度；L 是可实施源位置数。本轮 Gate1 n=k=4/6、s=2/3、名义 D=6/20，L=s+12。不是二维标量证据到三维向量 Maxwell 的推广。

## 1. 固定宿主消元（历史结构，非 R5 新理论）

J=D_m(Einc+GJ)，A(m)=I−D_mG，Y=SJ。令 m=m0+Pz，P 为不同的单位坐标列；A0 可逆，J0=A0⁻¹D_m0 Einc，E0=Einc+GJ0。相减得到

A0(J−J0)=P D_z[PᵀE0+PᵀG(J−J0)]。

设 U=A0⁻¹P、H=PᵀGU、B=PᵀE0、R=SU，便有

J−J0=U(I−D_zH)⁻¹D_zB，
Y−Y0=R(I−D_zH)⁻¹D_zB。

行列式引理给出 det A(m)=det A0 det(I−D_zH)。该等式保留完整离散多次散射，不是 Born 近似。任意宏观孔洞改变很多行，不能仅因其只有一个 contrast 参数就把变化当作 rank-one。`debug.json` 的非零宿主/四坐标检查相对误差 3.70e−16；这只是实现核查，不是 Gate2。

## 2. 对偶图与 Q6/Q20 的真实含义

单接收行下定义 w=R(I−D_zH)⁻¹D_z、v=R+wH，则 w_i=z_i v_i。旧数据条件 wB=t。B 满列秩且 k=2s 时，取 w0=t B†、NB=0、rank N=s，得到 w=w0+cN。

实材料要求 w_i conjugate(v_i)−conjugate(w_i)v_i=0。复化时 c 与 d=conjugate(c) 先独立，形成 2s 个双一次方程。一般的多齐次交数是 binomial(2s,s)，并非每例都有 D 个有限、互异、实、box 内、无节点/无极点的材料。历史代码用 maximal minors、leading block 消元及乘法矩阵求浮点端点；伴随共轭性、实性、box、node、pole、原正演验收。

本轮没有证明该特定系统的所有投影无穷远分支、节点、重根或 box 外边界均被排除。20/20 数值端点不是20个严格认证复根，也不是材料 box 全根覆盖。Q6/Q20 的系数及解析 Jacobian 在四/六材料测试中核查，`formula_tests.json` 的多项式误差 <1.6e−16、差分导数误差 <1.1e−11。

## 3. 新同频照明与不可消歧空间

对于相同宿主、接收器、频率、支持，真实新源的宿主场限制为 b，响应增量

h_b(c)=(w0+cN)b，Δh_b=(c_i−c_j)Nb。

证明仅使用固定材料下的源线性性与上节参数化。若 b=Ba，则 Nb=0，h_b=ta，对整条严格旧数据等响应纤维成立，无需纤维有限或已被算法发现。它并非 Z-ALG 独占能力；任何正确响应模型均能预测这一不可能性。

历史双根独立验证：最大对偶/完整正演差 1.29e−15，旧源组合差 3.47e−18；新源最大差 4.42e−3。浮点误差数值是观测，不是精确等式的替代证明。数值 NB≈0 与理论 NB=0 分开。

## 4. 有限集的一源分离与功率

实际源字典 Bplus 与系数 a 给 b=Bplus a。每对差异行 d_ij=(c_i−c_j)N Bplus。若所有 d_ij≠0，且硬件可行系数有非空开放区域（或归一化前的可行锥有开放区域），有限个超平面的并不能充满该区域，故一般相干组合可分离有限根集。限定离散菜单时，这个存在性不再自动保证菜单内有合适源；非零也不保证超过误差。

本轮实际源是二维理想 impressed line current：E_l(r)=(i/4)H0^(1)(κ|r−r_l|)。所有源在活动 cell 外，源位置固定。其自由空间远场相位为 exp(−iκ rhat·r_l)。归一化辐射能量

Q_ij=(1/2π)∫ exp(iκ rhat·(r_i−r_j))dθ=J0(κ|r_i−r_j|)，P(a)=a*Qa。

用2048角度积分核查 Q，误差2.45e−15。所有测量 a*Qa=1、持续时间1，且 ||a||2≤4。后者限制极端抵消/超指向电流；并不是完整的有损、互耦、加载天线功率模型。功率仅为校准后的发射自由场辐射能量/单位长度。实验是物理线源模型的数值实施，未做硬件实验。

单对 unconstrained 最优方向 a∝Q†d* 是标准 Cauchy–Schwarz/广义 Rayleigh 结论。本轮许多该方向因 feed norm>4 被拒绝；不暗中放宽功率限制。max-min 是有限可行菜单的枚举最大值，不是求出了连续全局最优照明。Fisher 比较用单解实 Jacobian，最大化旧+新 Gram 的正则化 logdet，也是标准局部设计。

## 5. 直接补基、稳定性和节点

满行秩 Bfull 时，w=tfull Bfull†；方阵时等价右乘逆。若 v_i=R_i+(wH)_i≠0，z_i=w_i/v_i。这里没有代数求根。噪声下先报告原始复 readout、虚部和 real-box projection，再报告独立联合 TRF64，不能把投影后值冒充原始精确解。

两材料的差满足 δw(I−HD_z′)=v(z)D_δz。若整个声明域具有 min_i|v_i|≥ν>0、||I−HD_z′||≤C、σmin(Bfull)≥β>0，则

||δz||2≤C/(νβ)||δtfull||2。

这是标准线性代数的条件稳定性界。本轮没有给整个 box 的 ν 正下界，因此不能据此宣称全域稳定恢复。T3/T4 的 Bfull 条件数约2.68e6/4.06e5，真实揭示补基也可能噪声不稳定；联合非线性 TRF 仍须公平比较。

若 w_i=v_i=0，保持同一 w 而改变 z_i 仍满足构成式；只要 resolvent 保持可逆，所有同频/同接收源均看不到该变化。`formula_tests.json` 用 G=0、S_i=0 的合成节点反例检查处理逻辑；不是声称六个真实 Hankel cell 已发现此精确节点。连续 d 维实纤维每个复标量观测至多去掉两个局部实自由度，ceil(d/2) 只是必要下界，不能当作全局唯一性定理。

## 6. 局部 Krawczyk 证书与序列化审计

F(m)=pack(S A(m)⁻¹D_mE−y) 是实方阵方程。解析导数

∂Y_s/∂m_i=(S A⁻¹)_i (E+GJ)_{i,s}，

注意转置求解是 Aᵀ，不是误用复共轭伴随。取点预条件 C，区间框 X，计算

K=x−CF(x)+(I−C J(X))(X−x)。

区间 Gauss–Jordan 每个 pivot 排除0；若 K⊂int X 且 ||I−CJ(X)||∞<1，则存在且唯一实根在 X。运算使用 mpmath.iv 64→128→256 bit，C/G/E/S/y 是精确二进制64常数。证书没有涵盖 Hankel 函数评价误差、真实积分核误差或连续 PDE 离散化误差。

初始 JSON 框端点从高精度舍入到 double。为防止把舍入报告误当原证书域，`qualify_predictions.py` 将报告框每端向内缩一 ULP，严格检查原 K 的向外包围仍落在这个内框；内框包含于实际认证框，原 contraction 界仍适用。这给出可精确重建的 binary64 安全框。正式使用 `qualification.json` 中 `certified_safe_binary64_box`；早期 qualification 快照保留但不作为最终证书。

对这些安全框，区间完整字典响应再乘精确的 source 权重，得到局部根的响应包围。T0 两根框不交且选定源分离下界为 clean 0.0040682、noisy-central 0.0056466。这些是中心数据的微小根框，不是含初始观测误差的根管！

## 7. 噪声必须移动 w0；中心根不等于可行集

w(t,c)=tB†+cN，δw0=δt B†。冻结 w0 后只移动 c 会遗漏初始观测误差。公式检查见 debug 的 noise_center_error=4.26e−17。

本轮数值局部根管使用原非线性旧+新数据拟合，在冻结的局部半径框内找 componentwise |residual|≤1.2e−4 的 witness。成功只是局部找到相容点；失败只是局部搜索未找到，不能排除该框，也不能排除 box 其他区域。它不等价于严格 interval noise tube。

全噪声可行集 F={m∈box: |F_old(m)−y_old|≤ε+δ} 的覆盖需要对候选框补集做严格排除。本轮没有该证明，因此全域返回 UNKNOWN。中心双根即使都认证并高度分离，也不能据此下 CERTIFIED_UNIQUE 的噪声全域结论。固定半径失配的真实模型误差超过预定 δ 于4/6目标，噪声+失配保障的前提实际失效。

## 8. 非原创部分与尚未成立的命题

宿主消元、对偶源线性性、旧源张成空间不可能性、有限超平面分离、Gram最优匹配、Fisher设计、Krawczyk、满秩补基及上述稳定性界，均是既有数学工具或直接推论。本轮贡献是可复核 qualification/反例与受限实验链，不把这些工具重新命名为新理论。

未成立：box 全根认证、噪声可行集覆盖、独立于普通 reciprocal/homotopy 的测量优势、实际天线/复杂已知宿主/未知支持/3D/PEC 的效果。不得通过“已运行并归档”将其改写为成功。
