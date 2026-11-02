# MAX-STRUCT-R1：主线程数学审查

本文件给出本轮实际使用的证明，不以再次计算 Atlas 恒等式到机器精度作为新证据。状态 `PROVED_UNDER_ASSUMPTIONS` 指下列精确数学模型与假设；`DERIVED` 指在其上得到的条件推论；物理离散和浮点回执另行标 `EMPIRICAL`。基础线性代数、磁场形式、shorted operator、单调性与源编码不主张原创。

## 1. M：固定的是运动学空间，不是材料算子

物理直觉：固定 μ 时，无磁荷约束不随 ε 改变；但恢复 E 需要 H 的 curl，所以 H 看起来平滑不足以证明表示更有效。

从 curl E=iωμH 与 curl H=J_ext−iωεE 得

\[
E=(i/\omega)\eta(\operatorname{curl}H-J_{ext}),\quad
\operatorname{curl}(\eta\operatorname{curl}H)-\omega^2\mu H
=\operatorname{curl}(\eta J_{ext}),\quad \operatorname{div}(\mu H)=0.
\]

仅当 μ 为常数标量时最后一式简化为 div H=0。外加电流附近 ε 已知才可固定右端。观测 E 是 H 的固定线性读出，要求接收器处在所保存全域状态中，或已经完成固定外域边界消元；ROI 的任意 H 通过体极化外传仍有
\(P_c=(i/\omega)(I-\varepsilon_b\eta)\operatorname{curl}H\)，是 η/H 双线性。没有消除材料非线性或共振。

有限兼容 Hilbert 空间中，P、Q 正交互补，ε 实 SPD，aI≤ε≤bI。令 B=QεQ 限制于 ran Q，L=P−QB⁻¹QεP。对 u∈ran P，Lu 是满足 PE=u 的 ε 能量最小延拓：对 Qv 变分有 QεLu=0。因此
\[
a\|Lu\|^2\le\langle Lu,\varepsilon Lu\rangle
\le\langle u,\varepsilon u\rangle\le b\|u\|^2.
\]
得到 \(\|L\|\le\sqrt{b/a}\)。块消元给 S=PεP−PεQB⁻¹QεP=L*εL；由 ε 的块逆，S⁻¹=Pε⁻¹P（ran P 上）。这是经典 Schur/shorted 结构。对边/面离散须使用正确 mass 内积与 adjoint；不对 PML、负介电、ENZ 或点偶极采样直接套用。

对复可逆 ε，代入 η=ε⁻¹ 可得
\(\operatorname{Im}_H\eta=-\varepsilon^{-*}(\operatorname{Im}_H\varepsilon)\varepsilon^{-1}\)。因此被动 ε 对应非正 η 虚部。先验保持 R(η⁻¹)，不能偷偷换为 η 平滑。

**资格：PROVED_UNDER_ASSUMPTIONS。** 本轮没有建立所需全域兼容 H(curl) 数值合同；M 表示和成像比较 NOT_RUN。DDA 的场中心 H 残差不得替代该资格。

## 2. T：有限端口排除的条件证明

原文实际参数 q=1−ε_r⁻¹。正 q、非磁、实无损、D 补集连通及原文其余条件下，Albicker–Griesmaier Thm6.1 对 0<α≤q_min 给出 B⊂D ⇒ αT_B≤fin Re F_q；T_B 使用 curl 磁 Herglotz 场，F 使用带 1/(4πr) 定义的磁远场。一般非零背景比较需要散射算子修正，不能直接作 Re(F₂−F₁)。详见 [原文逐项审计](research/delegated/theory/THEORY_QUALIFICATION.md)。

设有独立已知且对合法材料统一成立的 n₋(A_B)≤k，k<p。压缩 A_p=C*A_B C 的任何严格负子空间，经 C 映射仍为 A_B 的负子空间（负子空间上 C 必单射），故 n₋(A_p)≤k。若 \(\|\widehat A_p-A_p\|\le\delta\)，升序 Weyl 界给
\[
\lambda_{k+1}(\widehat A_p)\ge\lambda_{k+1}(A_p)-\delta\ge-\delta.
\]
逆否得到 λ_(k+1)<−δ ⇒ B⊄D。结论是“不全包含”，不是 B∩D=∅。未触发只能 NOT_EXCLUDED。若从完整算子误差压缩，界乘 ||C||²；若直接证明压缩矩阵误差，则不必恢复被压缩掉的无限维尾部。

**资格：PROVED_UNDER_ASSUMPTIONS；物理前提 INCONCLUSIVE。** 原文的有限维例外空间包含材料、外部 Calderón 传播/迹和紧算子谱；本轮未得到统一可用 k，更未得到全误差 δ。任意有限矩阵都有有限负方向，不能据此认证。没有使用 GT 反推 k，也没有以酉性缺陷冒充全误差界。闭壳与非零有损 host 分别不满足相应条件；没有推导出闭壳外包络定理。

## 3. F：有限幅度同场材料族

物理直觉：材料只通过 ε 对实际存在电场的作用进入 Maxwell。某张量方向从未被照明激发，则无论接收器多密，均看不到沿该方向的本构变化。

固定 μ、源、频率。若锚点总场 U₀ 和电流 j₀满足 j₀=X₀U₀、U₀=B+Gj₀，且 VU₀=0，那么 X_t=X₀+tV 满足 X_tU₀=j₀。由新材料 forward 唯一性，U_t=U₀，j_t=j₀，所有接收 S 的 Sj 相同。在连续模型中，Δε E_s=0 令两条 Maxwell 方程及出射场保持；还需新 ε 的本构合法性与正演唯一适定性。这是有限变化，不只是 Jacobian 核。

对于复互易张量，令 U=PΣQ*，r=rank_C U。取 P=[P_r,W]，V=\bar P ZP* 是保持 Frobenius 范数的酉合同；Vᵀ=V 当且仅当 Zᵀ=Z。VU=0 等价 Z 的前 r 列为零，对称性又令前 r 行为零，所以
\[
V=\bar W K W^*,\quad K^T=K,\quad
\dim_{\mathbb C}\ker\mathcal A_U=(3-r)(4-r)/2.
\]
实对称 V 要求同时湮灭 Re U、Im U，故对 [Re U,Im U] 使用实 SVD，其实核维数同式。已知轴对角时，第 i 个对角值只在整行 U_i,:=0 时自由；块常时叠加块中所有位置全部源，不能取中心点的核。数值近零阈值只能用于浮点秩诊断，不能把弱方向宣布为严格核。

对称矩阵的 Frobenius 正交基是 E_ii 和 (E_ij+E_ji)/√2。在 U 的左奇异坐标里，后者右乘 Σ 的范数平方为 (σ_i²+σ_j²)/2，且六个输出正交。因此六项谱为
\[
\sigma_1,\sigma_2,\sigma_3,
\sqrt{(\sigma_1^2+\sigma_2^2)/2},
\sqrt{(\sigma_1^2+\sigma_3^2)/2},
\sqrt{(\sigma_2^2+\sigma_3^2)/2}.
\]
最弱 dyad 为 \(\bar p_3p_3^*\)，不是一般为 Hermitian 的 p₃p₃*。实对称类以实堆叠奇异值替代；对角类为三行场范数。这个谱没有接收器，不是完整 Maxwell Jacobian 谱。

### 连续解析控制

取 d=e_z，两个横向源 e_x exp(ikz)、e_y exp(ikz)。任意光滑紧支撑 b≥0 下，ε=I+b e_ze_zᵀ 与 ε'=I+3b e_ze_zᵀ 都使 εE=E。因此真空入射 E,H 已是两个材料的完整出射解，散射严格为零；正定性和通常的散射唯一性假设给同场结论。本轮用球内 smooth bump、支撑半径 .1λ，外部有限接收作数值控制。这是同方向多极化不能激发轴向材料的具体反例，不是未知目标恢复。

### 非平凡严格被动锚点与光滑方向

连续两源可取 q=E₁×E₂（双线性叉乘），V=αb qqᵀ，则 qᵀE_s=0，不涉及可能退化的 UᵀU 逆。只在所需内域假设场正则，b 光滑紧支撑即可避免界面问题。无需除以 ||q||；q=0 时方向自然消失。本轮离散实现仅用全局有界尺度固定最大 cell 算子范数 .12，基准 Im_H ε=.25I，故更新损耗下界至少 .13；实部的原有正定裕量更大。

**资格：PROVED_UNDER_ASSUMPTIONS；具体被动锚点数值为 EMPIRICAL_FINITE_DDA。** 每个网格重构核所得材料不完全相同；固定粗网格 q 经线性插值/外推再乘解析 bump 的细网格试验必须单独报告，不能把两个网格各自机器精度同场当作同一连续材料的收敛证明。

## 4. 有限噪声下的一条合法弱射线

写 L₀=I−X₀G，R₀=L₀⁻¹，L_t=(I−tVGR₀)L₀。完整有限差分为
\[
F_t-F_0=t S R_0(I-tVGR_0)^{-1}VU_0.
\]
在 0≤t≤T，若 a=||VGR₀||_F 且 Ta<1，则
\[
\|W S L_t^{-1}\|_2\le
K=\frac{\|W S R_0\|_F}{1-Ta},\quad
\|W(F_t-F_0)\|_F\le tK\|VU_0\|_F.
\]
Frobenius 是算子范数的上界，所以这是区间界而非只在若干 t 采样最大值。实现通过已有 field LU 使用 R₀=I+X₀(I−GX₀)⁻¹G，并收费全部 3N 个额外 RHS。浮点 norm 未作 outward rounding，本轮不称严格机器证书。

若 ||V||mat=1，t_blind=min(T,δ/(K||VU₀||)) 是该射线的构造性歧义下界。实际 y 的锚点残差 r₀ 已占用噪声球半径 ρ，必须改用 δ=ρ−r₀≥0；三角不等式才保证另一候选也解释 y。两个候选间距≤2ρ 仅说明存在相交噪声球，不证明当前 y 同时合格。

本轮弱源振幅预设 .02，三源复秩达到3，因而方向不再精确同场。实测界和实际 y 残差只是有限噪声歧义，不能被包装成严格零核。合法半径由 max_i||tV_i||≤.10 给出，基准损耗 .25 保证整段严格被动；只声明这一条射线，不扩为材料椭球。

## 5. 新照明、接收与信息

所有原场相同意味着任意新增接收器仍相同；旧源线性组合只取已有解的线性组合，所以无新物理信息。真正新源的有限端口分离由候选 transfer 决定，局部核被破坏不足以保证超过噪声。

在已定义的外部源系数功率度量 R>0 下，A=W(H_a−H_b)，约束 c*Rc≤1 的最大差为 σ_max(A R⁻¹/²)。本轮使用理想横向平面波通道的系数 R=I；不是实天线口功率标定，也没有执行真实追加采集。C 臂确实重新求解独立方向/极化对两个材料的响应；最优相干组合仅是候选内预测。此处不能声称 CIE 获得了新数据后的成像优势。

## 6. D：局部 block 模型与源基诊断

精确根、实化坐标、共享材料无节点、内点无 prior 下，δe=Bδa−Cδm；逐 cell 权重使 W_β C=CΔ_β，且 CᵀC 为 cell 对角。因而 (CᵀW²C)⁻¹CᵀW²B=(CᵀC)⁻¹CᵀB，材料 β 在根处一阶抵消。

精确 current/material 两块 sweep 的 M=H_mm⁻¹H_amᵀH_aa⁻¹H_am 与 H_mm⁻¹/²H_amᵀH_aa⁻¹H_amH_mm⁻¹/² 相似。联合 LS Hessian 半正定，其 Schur 补保证 0≤该矩阵≤I，故谱在[0,1]。非零 residual 的权重导数、prior、活动面与单 PR 更新均不能偷套该谱；本轮单 PR 在根的方向差商甚至不构成线性 Fréchet 导数。

完整源酉变换 E'=EU、Y'=YU，右侧源 metric 必须 Ω'=U*ΩU。则 trace(R'Ω'R'*)=trace(RΩR*)。完整矩阵梯度运输且全局 PR 内积、线搜索、历史重置一致时等变；逐列 PR/步长属于依赖源坐标的分解，即便耦合目标正确也可产生不同有限步轨迹。这是现有 source encoding 的机制审计，不是新的观测或新的成像算法。

若连续映射作用于紧正向不变集合且每个非固定点有连续 V 严格下降，则沿轨迹 V 收敛；任一极限点 z 的连续下降差 V(z)−V(Tz)=0，故必固定，排除非平凡周期。双 β 无共同下降函数时此证明不适用，但不能反推 chaos。历史 wall 截断、短时非正规放大和缺失 PR 历史均不提供周期证据。
