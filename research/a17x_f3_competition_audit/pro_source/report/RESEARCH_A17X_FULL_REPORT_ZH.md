# Research A17-X：Maxwell 反馈原生材料歧义探索

## 有限激发坐标、双重反馈修饰的 current exchange 与多源材料可实现性

**日期：2026-10-17。性质：原文审计、独立数学推导、定向原创性审查及下一轮实验协议。**

本报告没有运行新的真实 Maxwell 网格实验、GPU 训练或 P1 实验。随包程序只做小型复矩阵代数核验；不能把它当作三维成像、多解发现或后验采样性能证据。P1 的最终结果保持 UNKNOWN。

## 执行摘要

**科学判断：存在值得研究的有限材料传输结构，但 A17 目前还不是已经成立的“材料多解生成原理”。** 本次找到的最具体路线不是材料 Jacobian SVD 后接 NN，而是：

$$
\text{局部材料变化}
\longrightarrow \text{等效材料激发 }T
\longrightarrow \text{线性接收条件 }OT
\quad+\quad
\text{非线性共享材料相容条件 }T=D(B+HT).
$$

这个表示把观测歧义和材料可实现性分开。在某些明确条件下，等效激发可以沿直线变化，材料却沿一个**精确有理函数族**变化，全部高阶 Maxwell 反馈已经包含在内，不需要逐步 GN 修正。多源条件并不自动满足；它们可能把这些方向全部消灭。

本次最值得保留的三个结果是：

1. **有限双侧 exchange 恒等式与剩余误差界。** 原 A17 的 current-space 修改同时更新内部反馈和接收映射，有限材料响应之差具有“左侧有限反馈修饰—原 current Schur 核—右侧有限激发修饰”的结构。它不是把冻结 Jacobian 差直接冒充有限非线性差。
2. **反馈行相容的精确材料图册。** 一个可计算的线性联合核，给出一类多源共享、receiver-dark、可精确积分的激发方向。该类只是充分子集，不等于完整材料弱空间。联合核非零时可构造有限有理图册；为空时不能推断没有其他弯曲路径。
3. **显式激发坐标中的密度变化。** 单标量照明、场不为零且材料本构可逆时，可以把材料和等效激发互相转换，并计算实 Jacobian 体积因子。这为正确 noisy posterior 提案提供了一个不依赖 many-to-one corrector 的接口。多源和三维各向同性矢量材料不能直接套用无约束单源公式。

**原创性审计改变了定位。** DCTMC（Levinson–Markel，2016）已经把非线性散射写成线性测量的 T-matrix 与局部材料约束；CSI、非辐射电流反演和 2025 年 current-to-permittivity 网络也已有非常接近的思想。因此，“用电流代替材料”“材料由电流解析读出”“把接收条件变线性”本身都不够新。可争取的贡献是：**多源精确有限传输条件 + A17 特定双侧更新 + 未解析缺陷控制 + 在同成本下优于精确局部 T-matrix 和普通 goal-oriented ROM 的真实收益**。

**执行裁决：解析机制与小实验 CONDITIONAL GO；大规模 learned Explorer / Diffusion / Flow NO-GO；强 IEEE TAP 论文尚未成立。** 缺口不是更多网络模块，而是合法材料中的非平凡多源见证、便宜且不空泛的误差控制，以及相对于最强非 A17 对照的总成本优势。

---

# 0. 实际读取范围、证据等级与缺失项

## 0.1 本次实际读取

原始仓库方面，直接读取 A17 的 README、closure_final/START_HERE.md、PAPER_TAP_FINAL.tex、CORE_THEORY.tex、CONDITIONAL_EXCHANGE_THEORY.tex；直接读取最新 Research F 指定 commit 的 RESEARCH_REPORT_ZH.md 核心结果、LXZ 指定 commit 的完整报告、A26 正式报告主要内容、A30 指定版本 E0/E1 决策，以及 C1 START_HERE。具体版本与入口见文末来源清单。[P1–P6]

Library 中取得并读到 Research F 完整中文报告；取得 Research L 正式完整 ZIP 并重点核对其原 A17/局部几何定理，而不再沿用早先“只有代数附件”的缺失状态；取得 Z ZIP 并读取正式报告与理论说明；读取 X 正式报告的主要执行摘要、证据与 target 边界；四命题文件已定位，但未完成每一段的独立复核。检索到此前 inverse-SLAM 的广义特征值与 sin²θ 推导。它们是此前研究材料，不是外部同行评审证据。[H1–H5]

跨会话检索只返回了部分历史上下文，**没有获得用户列出的每个聊天 Session 全文**。因此不声称全文重读全部会话。原始 F 二阶 Maxwell 向量缺失仍按其最新报告处理；本次没有偷偷用 LU 重新生成。P1 没有读到确认结果，不作任何成功推定。

外部文献方面，对 DCTMC 理论原论文的方程、数据相容/局部性迭代及成本部分进行了实际核对；读取 Deep Null Space Learning、e-GLASS、SOM-Net、DINO 几何 MCMC 的可解析方法部分，以及 2025 current-to-permittivity 网络 v1 的方法全文相关部分。非线性 data-consistent networks、1994 非辐射电流研究与 2026 inverse-Born neural operators 尚有全文访问限制，明确只使用取得的摘要/出版页面，不声称逐条排除了其定理。

PDF 截图接口对本轮若干论文返回错误；报告没有使用无法视觉核验的图表数值。数学比较依赖可解析公式和文字，不把图像猜读当证据。

## 0.2 最新已确认事实意味着什么

最新 F 数值报告在两个 development 条件上得到合法、noise-compatible 的 physics-only 重建，但两例来自同一 parent，且限定已知 smooth 家族及若干几何条件。最终残差约 22.2396 / 21.5121，低于原阈值 23.7854；bound-scaled 最优性约 0.04527 / 0.002001，均未过事前 0.001 门槛，故 F2 四臂传输没有启动。不能把“有两张合法图”写成“找到两种解释同一份数据的材料”。[P2]

两个旧 TRAIN smooth 点去除相同材料的权重缩放 gauge 后，有效材料维数与响应数值秩均为 7。它们支持局部有限噪声宽度的研究，不提供非零 exact tangent kernel 的证据。最新审计还说明，大 decoder 二阶响应中存在可消除的坐标加速度，不能全部归因于 Maxwell 高阶反馈。[P2]

A26 未建立显式 OM/feedback 在复杂材料上的独立增量；A30 Sparse-64 的约 35.69% 是 truth 拟合的离线容量，不是在线成像；LXZ 新增显著结构 world 为零，且概率比较是在径向 proxy target 上进行。它们共同排除“已有材料歧义与物理加速都已证实，只差训练网络”的叙事。[P3–P5]

## 0.3 本文命题的等级

下文区分：**原 A17 已证明**、**已有数学直接应用**、**本次在明确假设下推导**、**待实验假设**。本次推导正确，并不等于相关命题在全部文献中从未出现；尤其 T-matrix、Schur、隐函数、谱扰动与 MH 不会被重新命名为原创基础理论。

---

# 1. A17 原论文真正证明了什么

## 1.1 物理导数与多源共享材料

固定几何、频率、源与材料状态，对照明 p：

$$
L j_p=X e_p^{\rm inc},\quad L=I-XG,
\qquad L\,\delta j_p=B_p s,
\quad B_p s=DX[Q_Ts]E_p^{\rm tot}.
$$

Q_T 把共同实材料坐标 s 展开到材料空间，并按指定材料质量矩阵归一化。于是

$$K_p=L^{-1}B_p,\qquad J_p=SK_p.$$

所有照明共享 s，电流分别存在。材料法方程在联合实化后的空间中计算。实际 A17 的 DDA 使用材料相关极化率，不能把其 X、DX 偷换成二维线性 X=diag χ。[P1]

## 1.2 指定 Galerkin 增广

令 V、C 为正交 current bases，V* C=0。假设 V*LV 可逆，定义

$$
R_V=V(V^*LV)^{-1}V^*,\quad K_V=R_VB,
$$
$$
\bar C=(I-R_VL)C,\quad Y_C=S\bar C,
\quad N_C=C^*(B-LK_V),\quad \Gamma_C=C^*L\bar C.
$$

若 Γ_C 可逆，原文精确证明

$$
K_{V+C}-K_V=\bar C\Gamma_C^{-1}N_C,
\qquad J_{V+C}-J_V=Y_C\Gamma_C^{-1}N_C.
\tag{1.1}
$$

消去 retained coordinates 即得。Y_C 同时包含新增电流的直接辐射和经 retained currents 返回的辐射；N_C 包含材料对新增空间的直接及间接注入；Γ_C 包含新增与 retained 空间之间的反馈闭环。

注意：S C=0 不推出 S bar C=0。原始 receiver-dark current atom 经保留空间反馈后可能可见。反过来，强内部响应也可能因接收几何而几乎不可见。[P1]

## 1.3 Paired GN closure 的精确含义

固定相同 B、r、Λ≻0 和线性 prior 项 ℓ，令

$$
M_V=(J_V^*J_V+\Lambda)^{-1},
\quad s_V=-M_V(J_V^*r+\ell).
$$

原文证明

$$
\begin{aligned}
s_{V+C}-s_V={}&-M_VJ_V^*Y_C(\Gamma_C^{-1}N_Cs_{V+C})\\
&-M_VN_C^*\Gamma_C^{-*}Y_C^*(r+J_{V+C}s_{V+C})
\end{aligned}
$$

从而

$$
s_{V+C}-s_V\in\operatorname{Range}M_V[J_V^*Y_C,N_C^*].
\tag{1.2}
$$

这是“指定 current 模型改变所造成的完整冻结材料步差”的支撑空间，不是合法材料歧义空间，不是 posterior tangent，也不是完整非线性更新误差证书。复 current factors 必须联合实化；多源情况下一个复 current atom 不等于一个实材料维度。[P1]

同维交换 [V,A]→[V,C] 使用两个相对于相同 core 的 Schur 差。若共同 core 奇异而端点可逆，应直接算端点，不能强行使用不存在的 core inverse。材料活动约束变化时，原 unconstrained quadratic 定理不能代替新的 KKT 分析。[P1]

## 1.4 三条不能跨越的边界

第一，(1.1) 只比较指定两个 Galerkin 空间。完整误差是

$$J_F-J_V=(J_{V+C}-J_V)+(J_F-J_{V+C}).$$

第二，r=0 且 ℓ=0 时，各模型 GN 步都为零，Jacobian kernel 却可以完全不同。**GN 更新准确并不说明弱空间准确。**

第三，原 Schur inverse 本来就包含保留闭环中的任意次反馈。A17 的“冻结二次”指材料优化目标，不意味着原 Maxwell 物理只保留二阶散射。

---

# 2. 弱方向的来源，以及 A17 能做的诊断

## 2.1 首先限制到真实材料先验

设材料由 x=g(q) 产生。去掉同材料 gauge 后，在材料质量坐标中取合法 tangent frame Q。正确导数是

$$A=W\,\mathcal R\{SL^{-1}B_\chi Q\},$$

其中 W、实化与源堆叠遵守固定数据合同。若处于先验边界，应使用可行切锥及二阶可行条件，不是无限延伸一个线性 tangent space。

box 被动性只是合法性的一部分。固定相位、相关长度、材料混合、峰值归一化、硬边界或频散约束也属于先验；不能从 ambient 像素零空间推断原 prior 自由度。用于几何的 gauge quotient 也不能自动被赋予均匀概率密度。[H1,H2]

## 2.2 注入—内部状态—接收三个量

对单位合法材料速度 v，可分别记录

$$i(v)=\|B_\chi Qv\|,\quad c(v)=\|L^{-1}B_\chi Qv\|,\quad o(v)=\|Av\|.$$

在固定且声明过的度量下，小 i 表示材料几乎不激发电流；大 i、小 c 表示内部方程对该激发的抑制；大 c、小 o 表示强内部变化却弱辐射到已测通道。解释为具体路径相消时，还应计算 direct 与 dressed 通道的复内积，不能仅凭一个范数比宣布发生共振或抵消。

全 Maxwell Jacobian 已经包含这些反馈；A17 的作用是**解释和有选择地计算它们**，不是给同一数据添加新的信息。原 prior 不允许的方向应先拒绝，不能归类为物理暗方向。

## 2.3 对“不同照明相互补偿”的一个重要修正

若各照明分别测量、分别堆叠，且采用相应独立噪声白化，则

$$\|Av\|^2=\sum_p\|A_pv\|^2.$$

一个照明的正响应不会和另一个照明的负响应在这个范数中消掉。Exact-null 要求每个照明都不变。只有设备实际测量相干叠加，或模型中存在被合法处理的 nuisance 混淆时，才有相应跨源抵消的解释。共享材料必须同时满足所有源，是强约束，而不是可自由配对的补偿参数。

## 2.4 真弱、false-weak 与未知

记 A_n=A_V+D_C 为增广结果，E=A_F-A_n。若对待测方向有经验证的 ε(v)≥||Ev||，则

$$
\boxed{
\max(0,\|A_nv\|-\varepsilon(v))
\le\|A_Fv\|\le\|A_nv\|+\varepsilon(v).
}
\tag{2.1}
$$

上界低于噪声/任务尺度才可称可靠弱；旧 ROM 弱而新下界超过尺度，才能确认 false-weak；剩余情况必须标 UNKNOWN。若要认证整个弱空间，必须界定 E 在整个该空间上的算子范数，几个随机 probes 不等于确定性认证。

反例：L=I₃，V=e₁，材料激发 B=e₂，S=e₂*。ROM 给 J_V=0，full 给 J_F=1。加入 e₃ 的 exchange 仍为零，不能据此认证这个材料方向不可观测。缺失 e₂ 是未解析 defect。

式 (2.1) 与类似思想已在 L/F 中提出。本次不将其作为新的核心贡献；真正的新接口在第 5 节的有限反馈版本。[H1,H2]---

# 3. 从冻结材料导数到完整有限反馈

## 3.1 精确恒等式与它包含的物理

以下固定 G、S 与源。设 χ₁=χ₀+h，ΔX=X(χ₁)−X(χ₀)，L₀=I−X₀G，R₀=L₀⁻¹。把全部源的总场写成矩阵 E₀。两端状态方程相减得到

$$
(L_0-\Delta XG)\Delta j=\Delta X E_0,
$$

因此只要两端状态算子可逆，

$$
\boxed{\Delta j=(I-R_0\Delta XG)^{-1}R_0\Delta X E_0.}
\tag{3.1}
$$

R₀ 已包含原材料中的多重散射。相对于强背景的一阶材料变化，并不等于真空 Born 一阶近似。反过来，精确算出 Δj 也不意味着存在非零 h 使 SΔj=0。

令 K=R₀ΔXG、b=R₀ΔXE₀，截取 M 项时有恒等式

$$
\Delta j-\sum_{r=0}^{M-1}K^rb=(I-K)^{-1}K^M b.
$$

若 ||K||≤q<1，余项上界为 q^M||b||/(1−q)。有限维 ρ(K)<1 足以保证 Neumann 收敛；||K||<1 只是方便的充分条件。I−K 可逆却不保证级数收敛，例如 K=2I。局部有理求解只要求相应逆存在，并不要求某个 Neumann 范数小于 1。

这是已有 resolvent/Born 结构。本文推进的是它与 A17 current-space 变化、材料可实现性和可计算误差的耦合，不声称发现了 resolvent identity。[H1,L1,L10]

## 3.2 高阶导数：本构与场变化都必须保留

在当前材料处记 R=L⁻¹，E=e+Gj，X_u=DX[u]。对实材料坐标取方向导数：

$$Dj[u]=R X_uE,$$
$$D^2j[u,v]=R\{X_{uv}E+X_uG Dj[v]+X_vG Dj[u]\},$$
$$
\begin{aligned}
D^3j[u,v,w]=R\{&X_{uvw}E
+X_{uv}GDj[w]+X_{uw}GDj[v]+X_{vw}GDj[u]\\
&+X_uGD^2j[v,w]+X_vGD^2j[u,w]+X_wGD^2j[u,v]\}.
\end{aligned}
\tag{3.2}
$$

一般阶用带标签集合 A 表示：

$$D_Aj=R\left[X_AE+\sum_{\varnothing\ne B\subsetneq A}X_BG D_{A\setminus B}j\right].$$

当 X 线性、所有方向相同为 v 时，

$$D^nj[v,\ldots,v]=n!(RX_vG)^{n-1}RX_vE.$$

混合方向的乘积不能交换顺序。材料 decoder g 还贡献链式项，例如

$$D^2(F\circ g)[v,v]=D^2F[Dgv,Dgv]+J D^2g[v,v].$$

沿材料路径 χ(t)，对固定 ambient 方向 v，

$$
\frac{d}{dt}(J_{\chi(t)}v)
=S R X_{\dot\chi}G R X_vE
+S R X_{\dot\chi,v}E
+S R X_vG R X_{\dot\chi}E.
\tag{3.3}
$$

若 v 自身来自移动 prior frame，还要加 J dot v。这个公式明确显示传播 resolvent、材料注入和内部场同时改变；不能只更新 L⁻¹ 而冻结不该冻结的 injection。

## 3.3 没有额外“材料历史自由度”

静态材料路径上的场可以写成变系数 ODE，但精确场仍由端点材料唯一确定。若 A(χ)=I−GX(χ) 可逆，则

$$E(\chi_1)=A(\chi_1)^{-1}A(\chi_0)E(\chi_0).$$

闭合合法材料路径恢复同一场。Dyson/Magnus 的有序写法不产生新的物理记忆，除非另引入真实时变、滞回或多值本构。算法的路径选择、prior 分支和数值误差可以依赖路径；这是另一回事。[H1]

---

# 4. 主表示：材料端口上的有限等效激发

## 4.1 一个必须明确的额外假设

假设本轮材料变化仅发生在已声明的 s 个物理端口中：

$$\Delta X=U D U^*,\qquad U^*U=I_s.$$

U 可以选择若干 voxel 的标量自由度，或每个三维 voxel 的三个极化分量。标量各向同性材料 D 为对角；三维各向同性为 blockdiag(dᵢI₃)；一般张量为合法的局部块。

这是假设材料算子变化具有局部/低秩支持，**不是**假设 latent 维数低。一个单参数的全域 contrast 缩放也可使 ΔX 满秩。若 U 需要覆盖整个大网格，后面的 s×s 反馈系统可能毫无成本优势。

非线性 DDA 中 dᵢ 是极化率增量，不必等于 χ 增量。每次读出材料后必须使用实际本构逆和总材料约束，不能要求增量自身被动。

## 4.2 精确端口模型

定义

$$
P=R_0U,\qquad H=U^*GP,\qquad O=SP,\qquad B=U^*E_0.
\tag{4.1}
$$

这里 P 是端口到完整电流的传播，不是 source 数；后文照明数记为 n_s。定义等效激发矩阵 T∈C^{s×n_s}：

$$
\boxed{
T=D(B+HT),\quad \Delta j=PT,\quad \Delta F=OT.
}
\tag{4.2}
$$

若 I−DH 可逆，令

$$Q_H(D)=(I-DH)^{-1}D,$$

则

$$T=Q_H(D)B,\qquad \Delta F=OQ_H(D)B.\tag{4.3}$$

**证明。** 将 Δj=R₀U D U*(E₀+GΔj) 中的端口材料激发定义为 T，即得。反向由 (4.2) 构造的 Δj 满足完整差分状态方程。矩阵行列式引理给出

$$\det L_1=\det L_0\det(I-DH),$$

所以端口反馈可逆性与完整新状态可逆性在该精确支持假设下等价。

这里有两个不同的困难：OT≈0 控制接收变化；T=D(B+HT) 控制这些电流是否能由同一个局部材料实现。材料 prior 是第三个独立条件。

## 4.3 为什么这个表示有用，又为什么它本身不新

冻结材料 Jacobian 在 D=0 附近只有 O D B；有限反馈还需要 H。相同 O、B 而不同 H 的系统可以拥有相同一阶响应，却有不同有限材料补偿。因此 frozen SVD 不包含全部有限传输信息；但 refreshed full Maxwell 方法当然可以间接获得它。

另一方面，(4.3) 就是背景修饰后的局部 T-matrix。DCTMC 已使用 A(I−VΓ)⁻¹VB=Φ，并研究线性数据相容与局部材料的关系。因此单独推出 (4.2) 或 (4.3) 不构成 A17 原创性。本文保留这层表示，是因为它让有限 exchange 与多源可实现性能够严格连接，而不是因为它没有先例。[L1]

---

# 5. 命题 A：双重反馈修饰的有限 A17 exchange

这一节是本次从原 A17 到有限材料传输的主要可计算接口。

## 5.1 条件与定义

固定完整 anchor χ₀、E₀、L₀、U。旧 current 空间 V 和增广 C 满足原 A17 的 core/Schur 可逆条件。定义

$$P_V=R_VU,\quad H_V=U^*GP_V,\quad O_V=SP_V,$$
$$\bar C=(I-R_VL_0)C,\quad \Gamma_C=C^*L_0\bar C,$$
$$N_U=C^*(U-L_0P_V).$$

原交换恒等式现在用于端口 RHS U：

$$E_C:=P_{V+C}-P_V=\bar C\Gamma_C^{-1}N_U.\tag{5.1}$$

因此内部传播与接收映射同时改变：

$$H_{V+C}-H_V=U^*GE_C,\qquad O_{V+C}-O_V=SE_C.$$

B=U*E₀ 保持相同。对同一个有限材料 D，假设 I−DH_V、I−DH_{V+C} 都可逆，记 Q_V=(I−DH_V)⁻¹D，Q_n=(I−DH_{V+C})⁻¹D。

## 5.2 结论

两个有限增量预测之间精确满足

$$
\boxed{
\Delta F_{V+C}(D)-\Delta F_V(D)
=
[\,S+O_VQ_VU^*G\,]\bar C\Gamma_C^{-1}N_UQ_nB.
}
\tag{5.2}
$$

左端是指定两个 reduced finite-increment 模型之间的差，不是未经剩余项的 full−ROM 误差。

物理上，原 receiver return 不再仅是 S bar C，而是经过有限材料反馈再次修饰；原端口注入 N_U 也要作用于有限激发 Q_nB。**有限变化会同时修饰 exchange 的左右两端。**

## 5.3 完整证明

由 resolvent 恒等式，

$$Q_n-Q_V=Q_V(H_{V+C}-H_V)Q_n.$$

因而

$$
\begin{aligned}
O_nQ_nB-O_VQ_VB
&=(O_n-O_V)Q_nB+O_V(Q_n-Q_V)B\\
&=[S+O_VQ_VU^*G]E_CQ_nB.
\end{aligned}
$$

代入 (5.1) 即得。没有使用 Neumann 收敛，也没有截断材料阶数。

当 D=tD₁ 时，Q_n=tD₁+O(t²)，故

$$\Delta F_n-\Delta F_V=tS\bar C\Gamma_C^{-1}N_UD_1B+O(t^2),$$

其一阶项恰是对应材料注入 U D₁B 的原 A17 exchange。

另外，P_VQ_VB 恰是**差分状态方程** L₁k=UDB 在 V 上的 Galerkin 解。它不等于声称“完整新状态 j₀+k 在 V 中”；j₀ 是固定 full anchor，可以不属于 V。若改用 ROM anchor，总场 B 的误差必须另计。

## 5.4 反例、边界与创新等级

完整 L₀ 可逆不保证 V*L₀V 可逆。例如 L₀=[[0,1],[1,0]]，V=e₁，full 稳定而 reduced core 为零。类似地，full 新材料稳定不排除 reduced I−DH_V 的伪极点。

(5.2) 是 Schur/resolvent 的组合结果，不应命名成一种新的普适散射定理。它的潜在增量是：利用**同一组 A17 current factors**更新有限材料反馈、接收误差和后续 action 选择，而不反复构造完整材料 Jacobian。真正的速度和稳定性优势仍需实验。

## 5.5 命题 A 的第二部分：未解析缺陷的结构化证书

定义端口残差

$$R_U=U-L_0P_V.$$

则

$$P-P_V=R_0R_U.$$

对 k_V=P_VQ_VB，完整有限差分方程残差精确分解成

$$
\boxed{\rho_D=UDB-L_1k_V=R_UQ_VB.}
\tag{5.3}
$$

**证明。** L₀P_V=U−R_U，代入 L₁=L₀−UDU*G；使用 (I−DH_V)Q_V=D，U 相关项相消，只剩 R_UQ_VB。

因此 full 输出误差为

$$\Delta F_F-\Delta F_V=SL_1^{-1}R_UQ_VB.\tag{5.4}$$

(5.3) 揭示应按 Explorer 当前有限提案产生的 RHS Q_VB 选择 current modes，而非仅按旧 GN residual 选择。

若已取得经验证的

$$\|H-H_V\|\le\varepsilon_H,\quad \|O-O_V\|\le\varepsilon_O,$$

并满足

$$\vartheta=\|Q_V\|\varepsilon_H<1,$$

则完整端口反馈也可逆，且

$$
\boxed{
\|\Delta F_F-\Delta F_V\|_F
\le
\frac{\big(\varepsilon_O+\|O_VQ_V\|\varepsilon_H\big)
\|Q_VB\|_F}{1-\|Q_V\|\varepsilon_H}.
}
\tag{5.5}
$$

**证明。** 记 ΔH=H−H_V、ΔO=O−O_V。因

$$I-DH=(I-DH_V)(I-Q_V\Delta H),$$

Banach 扰动引理给 I−Q_VΔH 可逆，并有

$$Q_HB=(I-Q_V\Delta H)^{-1}Q_VB.$$

再由

$$\Delta F_F-\Delta F_V=(\Delta O+O_VQ_V\Delta H)Q_HB,$$

取谱范数与 Frobenius 范数得 (5.5)。证毕。

一个保守的可计算入口是已验证 α₀≤σ_min(L₀)：

$$\varepsilon_P\le\|R_U\|/\alpha_0,
\quad\varepsilon_H\le\|U^*G\|\varepsilon_P,
\quad\varepsilon_O\le\|S\|\varepsilon_P.$$

更精细的 primal/dual 输出界可以替代全局范数。若 B 或 anchor prediction 有误差，分别加入相应传播项与基准输出误差。白化数据时 S 及多源堆叠的范数必须按实际噪声合同计算。

**这个证书的价值与限制。** 它不要求 ||R₀ΔXG||<1，也不必每步重新估计完整 L₁ 的最小奇异值；代价是需要可信的 anchor/端口误差界。经验残差小、Γ 条件数好或 held-out 平均误差小，都不能代替这些界。θ≥1 只表示该充分证书失效，不表示完整材料无解；此时应刷新、缩步或 full fallback。

## 5.6 同端口 anchor 的精确更新

接受材料 D 后，仍使用相同 U，可精确更新

$$
P'=P(I-DH)^{-1},\quad H'=H(I-DH)^{-1},\quad O'=O(I-DH)^{-1},
$$
$$B'=(I-HD)^{-1}B.\tag{5.6}$$

从新 anchor 再做 −D 可以恢复旧端口算子。此处 D 是本构算子增量，非线性 α(χ) 情形不能把 −D 与任意材料参数线性反向步混淆。

这提供同一局部支持上多步探索的复用机制，但属于 Woodbury/端点 resolvent 代数，不是物理路径记忆。换端口需要新传播列；长期浮点更新需要残差审计与定期重建。---

# 6. 命题 B：不是任意 nullspace——可精确积分的暗激发子空间

## 6.1 一条精确有限材料路径

先考虑标量对角 D。取材料激发系数 a，令

$$T_1=D_a B,\qquad D_a=\operatorname{diag}(a).$$

要求所有源的接收响应为零：OT₁=0。进一步假设，对于每个活动单元 i，存在同一个复数 γᵢ，使

$$
(HT_1)_{i,:}=\gamma_iB_{i,:}.
\tag{6.1}
$$

即有限反馈产生的多源场改变，在该单元与原多源场共线。定义

$$
\boxed{d_i(t)=\frac{t a_i}{1+t\gamma_i}.}
\tag{6.2}
$$

只要分母、I−D(t)H 可逆，并且本构逆得到的材料始终满足原先验，便有

$$T(t)=tT_1,\qquad \Delta j(t)=tPT_1,\qquad \Delta F(t)=0.$$

**证明。** 对每个单元，

$$D(t)[B+tHT_1]=D(t)(I+tD_\gamma)B=tD_aB=tT_1.$$

所以 T=tT₁ 精确满足完整非线性材料反馈方程；唯一状态解条件使其等于真实有限响应。OT₁=0 给同数据。证毕。

这不是三阶近似，而是一整条有限有理材料路径。材料沿曲线变化，内部电流沿直线变化。单元材料二、三阶项由 γᵢ 决定，不需要逐步 GN 去猜。

## 6.2 更强结论：条件本身可以组成一个线性联合核

假设所有使用的 Bᵢ,: 非零。令行投影

$$\Pi_i=\frac{B_{i,:}^*B_{i,:}}{\|B_{i,:}\|^2}.$$

定义两个对 a 线性的算子

$$\mathcal O_B a=\operatorname{vec}(O D_aB),$$
$$[\mathcal C_Ba]_{i,:}=(H D_aB)_{i,:}(I-\Pi_i).$$

以及

$$\mathcal E_B=\ker\mathcal O_B\cap\ker\mathcal C_B.$$

在 E_B 上定义线性的行系数

$$\gamma_i(a)=\frac{(H D_aB)_{i,:}B_{i,:}^*}{\|B_{i,:}\|^2}.$$

则对足够小的 a∈E_B，

$$
\boxed{d_i(a)=\frac{a_i}{1+\gamma_i(a)}}
\tag{6.3}
$$

精确给出一族同数据材料，等效激发为 D_aB。若 E_B 复维数为 r，且材料在原点附近允许这些复自由度，则该映射局部形成实维数 2r 的嵌入图册。

**证明。** 前半部分直接使用 (6.1)–(6.2)，把 t 吸收进 a。Ψ(a)=d(a) 在零点的导数是恒等嵌入，因为 γ(a) 对 a 线性，分母为 1+O(||a||)。故 DΨ(0) 在 E_B 上满秩，得到局部嵌入。小邻域中 I−DH 持续可逆；材料合法性必须独立成立。若材料只允许实参数，应在对应实化子空间中重算维数。

**算法意义。** 不必先构造完整网格材料 J 再做全部 SVD。可以对端口级的 O_B 与 C_B 做 matrix-free 联合约束求解，直接筛选“接收端暗、且反馈保留共享材料方向”的激发。这里没有排斥任何矩阵分解工具：小型 SVD、广义特征值或约束 Krylov 都可用；物理增量在附加的 finite-feedback 相容结构，而不是禁止使用 SVD。

## 6.3 它不是完整 ambiguity 空间

E_B 是“可以让等效激发保持直线”的充分子集，比 ker O_B 更严格。一般合法 fiber 在电流域也会弯曲，T(t)=tT₁+t²T₂+…，不需要满足 C_Ba=0。因此 E_B={0} 不能证明材料局部唯一，也不能否定普通非线性 continuation。

单源标量情况下，每个非零 Bᵢ 只有一个复分量，共线条件自动成立；多源情况下不是。三维各向同性单源时，每个单元已有三个极化分量，共线条件也不自动成立。**单源标量的漂亮构造不能无条件升级为多源三维 Maxwell。**

## 6.4 命题 B 的有限误差版本

对一般 a，将反馈分解为

$$HT_1=D_\gamma B+R,\qquad R_{i,:}\perp B_{i,:}.$$

仍用 (6.2) 提出材料。则真实响应精确满足

$$
\boxed{
\Delta F(t)=tOT_1+t\,OQ_H(D(t))R.
}
\tag{6.4}
$$

**证明。** 以 tT₁ 代入反馈方程的残差为

$$tT_1-D(t)(B+tHT_1)=-tD(t)R.$$

真实激发和该候选的差为 t(I−D(t)H)⁻¹D(t)R=tQ_H(D(t))R。左乘 O 即得。

若 OT₁=0，第一项消失，第二项为 O(t²)，但 (6.4) 本身是精确公式。它区分“初始接收端暗”与“反馈导致多源比例不再相容”。这可作为步长或方向评分依据，而不必用一个 frozen singular value 代表全部有限可靠性。

R≠0 不意味着不存在其他补偿：R 的后续输出可能抵消，也可能通过 T₂ 等弯曲电流分量修复。此式只诊断这条特定有理提案。

## 6.5 先验反例：固定相位会立即破坏这条材料族

设先验只允许 d=c u，其中 u 为实数。即使初始 a=c 合法，

$$d(t)/c=t/(1+t\gamma)$$

在 Im γ≠0 时通常不是实数。因此提案从二阶开始就离开固定相位先验。原 P0 的相位/峰值/形状约束不能靠 “Maxwell exact” 免除。

对 lossless、二元材料、相关随机场或高斯形状先验也有同类问题。要么证明该图册保持其支持，要么在原先验中额外求解相容条件；不能偷偷把 prior 放宽为任意复像素再宣布成功。

若 Bᵢ,:=0，除法图册失效，应保留材料变量而非加 ε。特别地，若改变支持内的原总场对所有源都为零，则 T=0 可精确满足任意小 D 的状态方程，产生“材料变化但电流完全不变”的另一类退化歧义。它与非零暗电流图册不同。

## 6.6 原创性判断

本节的线性联合核、显式有理图册和精确 misalignment 公式是本次独立推导的可检验结果。基础来自 contrast-source/T-matrix 代数；本报告没有证据证明这些特定表达从未被以等价形式研究。可将其作为论文候选命题，但不能现在宣称首次发现非辐射材料变化。它们值得推进，是因为条件可直接计算、失败可明确解释，并且能和第 5 节的有限 A17 exchange 相接。

---

# 7. 命题 C：多源 current–material realizability 的完整条件

## 7.1 任意电流差分，不局限于端口近似

令各照明电流差为 k_p。同数据要求 Sk_p=0；完整状态差分要求

$$L_0k=\Delta X(E_0+Gk).$$

在标量对角材料下，记

$$A_i=(L_0k)_{i,:},\qquad E_i=(E_0+Gk)_{i,:}.$$

当 Eᵢ≠0，存在同一个复材料增量 dᵢ 解释全部照明，当且仅当

$$A_i=d_iE_i,$$

等价于全部二阶 minors 为零：

$$A_{i,p}E_{i,q}-A_{i,q}E_{i,p}=0\quad\text{对所有 }p,q.$$

此时

$$d_i=A_iE_i^*/\|E_i\|^2.$$

若 Eᵢ=0，则必须 Aᵢ=0，且材料不能由这个局部比值确定。上述条件之后仍要检查原 prior、总被动性和新 L 的可逆性。反向从任何合法同数据材料出发，其完整电流必满足这些条件；故这是带可逆性与先验条件的精确对应。[H1]

在端口表示中，将 A 换成 T，E 换成 B+HT，条件相同。其 minors 对 T 至多二次，提供无需 material SVD 的有限约束问题；变量规模可能是 s×n_s，而不是 s。

## 7.2 近似相容不是合法材料，但有清楚残差

用上述最小二乘比值重建 d，记

$$\rho=L_0k-\Delta X(E_0+Gk).$$

真实电流差 k_true 满足

$$k_{\rm true}-k=-L_1^{-1}\rho.$$

所以若 Sk=0，真实接收变化是

$$\Delta F=-SL_1^{-1}\rho.\tag{7.1}$$

这允许用完整/带证书的传播判断一个近似 dark-current 提案是否真实可行。仅报告 row-collinearity loss 很小不够：场接近零、resolvent 大或白化敏感时，微小本构残差仍可能产生显著数据变化。

## 7.3 三维矢量与张量材料

三维各向同性材料在每个 voxel 需要同一个 dᵢ 乘全部三个极化分量和全部照明。把 3×n_s 的局部块展平成一个向量，仍使用“共线”条件。不能逐极化分量各恢复一个 d，却称为标量介质。

对 unrestricted tensor 增量 ΔAᵢ∈C^{3×3}，局部条件是

$$A_i=\Delta A_iE_i,\qquad A_i,E_i\in\mathbb C^{3\times n_s}.$$

它可解当且仅当

$$A_i=A_iE_i^\dagger E_i,$$

等价于 ker Eᵢ⊆ker Aᵢ。全部解为

$$\Delta A_i=A_iE_i^\dagger+Z_i(I-E_iE_i^\dagger).$$

若 rank Eᵢ=3，解唯一。推导来自线性矩阵方程：右乘 E† 构造一个解，再加在 E 的列空间上为零的任意矩阵。

但任意矩阵解未必 reciprocal、被动或属于指定材料族。必须对**总张量**施加对称性、Hermitian imaginary part 半正定以及实际频散约束。多频时每个频率不能独立挑选一个最方便的材料增量；共同物理参数必须同时解释所有频率。

## 7.4 与此前 sin²θ 推导的准确联系

设 K 将合法材料速度映射到完整多源电流速度，D=ker S 为接收暗电流空间。去掉 ker K 后，以 current 能量为分母，广义特征问题

$$K^*P_{D^\perp}Kv=\rho\,K^*Kv$$

的非平凡特征值对应 Range(K) 与 D 的 principal-angle sin²θ；当维数不匹配时，未匹配的正交方向补上特征值 1。等价地，单元最小二乘相容残差满足

$$\min_d\|A_i-dE_i\|^2=\|A_i\|^2\sin^2\theta_i.$$

这延续此前 inverse-SLAM 中“混淆空间夹角”的解释，不是新的主角定理。[H5]

必须注意，真正 likelihood 使用 S*W*WS，不是 P_{D⊥}；除非在特定等距规范下，它的广义特征值不等于裸 sin²θ。先验或 nuisance 正则化也会改变这种几何。角度可以解释强电流为何难见，但不能代替有噪声量纲的 likelihood 宽度。

## 7.5 为什么电流空间没有凭空创造材料自由度

线性化的多源可实现电流满足 L₀k∈Range(B_prior)。因此

$$\{k:Sk=0,\ L_0k\in\operatorname{Range}B_{\rm prior}\}
=L_0^{-1}B_{\rm prior}(\ker J_{\rm prior}),$$

其中 range 必须按共同实材料域理解，并正确处理 B_prior 的同材料 gauge/零注入核。电流域只是重新表示同一个局部问题，不会把原先验七维满列秩自动变成高维材料 fiber。

---

# 8. 非线性近邻解何时存在，何时不能连续探索

## 8.1 精确 fiber：常秩与法向补偿

令 f=F∘g 是合法 prior 坐标中的白化 forward。若在邻域恒秩 r，参数维数 d，则正规 fiber 维数为 d−r。满列秩时局部只有孤立精确解；满行秩时可取弱方向 V_w 张成 kernel，法向 V_s 使 Df V_s 可逆，从隐函数定理得到

$$q(z)=q_0+V_wz+V_sa(z),\quad f(q(z))=f(q_0).$$

a(0)=0、Da(0)=0，二阶补偿满足

$$Df\,V_sD^2a[u,v]=-D^2f[V_wu,V_wv].$$

若 Df 只有恒秩 r<m，需要合法选择独立输出坐标；没有恒秩依据时，不能删掉数据余核约束。隐函数定理是标准数学，本报告不重新将它包装为 A17 定理。[H1,H2]

## 8.2 二阶障碍与 prior 外在曲率

设 w 是真正 prior tangent-null 速度，Q 是材料正交 prior frame。沿 prior 路径的不可由坐标重选消掉的二阶项为

$$C_{\rm int}=D^2F[w,w]+J\,\mathrm{II}_{\rm prior}(w,w).$$

若存在精确同数据曲线，必须

$$P_{\operatorname{Range}(JQ)^\perp}C_{\rm int}=0.$$

该条件是二阶必要条件，不是有限曲线存在的充分条件。f(x,z)=(x,z²) 在原点有 z-null 但只有孤立精确解；换成 (x,z³)，二阶条件通过，三阶仍阻碍。

F 已完成“坐标加速度—prior 第二基本形式—Maxwell Hessian”的理论分离；本次使用它避免把有理图册失败全部归因于物理高阶项。[H1,P2]

## 8.3 有限噪声不能冒充 exact fiber

实际任务通常应研究

$$\Phi(q)=\tfrac12\|f(q)-y\|^2,$$

在 prior 支持内的近最优 likelihood 区域或明确 noise tube。对于近似弱速度，第一阶漂移不为零，不能直接套用 exact-null 二阶障碍去证明“不存在路径”。

合理的有限补偿问题是：尽量接近提议材料，同时满足预定 Φ_cap、原 prior 与最小材料进度。这个约束集合不是等概率区域；Gaussian posterior 也不以某个 95% noise ball 为硬支持。

## 8.4 固定弱/强划分的有效范围

固定 V_w、V_s 可在一个小的良态 chart 内使用；其尺度受 normal 最小奇异值 β、Jacobian Lipschitz 常数、prior 曲率、边界距离与 resolvent 条件控制。通常 β 小、曲率大时可行 chart 变小。

跨越谱簇、active constraint、局部场零点、端口支持边界或反馈极点时，应刷新 chart 或改用分片/全局提案。Pseudo-arclength 可以处理某些投影坐标的 turning point，但不能让孤立解产生路径，也不能把不连通的 prior 世界连接起来。

端点材料接近、预测相同也不保证可连续连接。两端可能属于不同局部解支；有限噪声是否连接它们取决于中间路径的 likelihood 与 prior，而不是所使用的反馈阶数。

## 8.5 一个应当修正的旧实验假设：anchor 不必是优化驻点

构造通过 χ₀ 的等预测 fiber，只要求 χ₀ 合法及 F(χ₀) 是所研究的预测；在 noise tube 内探索，只要求其满足声明的可行性；posterior sampler 也可从任意正密度合法点初始化。它们在数学上都不要求 χ₀ 已达到 MAP/GN 最优性门槛。

因此，最新 F1 未达旧最优性门槛不等于“该点不能做任何局部传输研究”。但旧 F2 的 NOT_RUN 记录必须保留，不能事后宣布过关。新 A17-X 可以事前登记一个独立的 feasible-anchor transport pilot，明确它不声称 anchor 是 MAP、不继承旧 F2 的成功资格，也不降低原 prior、数据与材料进度要求。若复用原七维 P0 点，它仍只是局部几何控制，不是高维复杂多解成像证据。

---

# 9. A17 动态控制探索几何：能证明什么，不能证明什么

## 9.1 Fisher 缺陷与弱投影

在同一实材料坐标下，令 A_n=A_V+D，rank D≤r。则

$$\Delta H=A_V^*D+D^*A_V+D^*D,$$

其 range 包含于 Range[A_V*Y,N*]，rank 至多 2r。多源复 factors 必须先合并实化，不能把 r 偷换成单个复 current basis 的列数。

ΔH 不必半正定。例如标量 A_V=1、D=−0.9，新 Fisher 从 1 降为 0.01。Galerkin enrichment 可以修正高估，也可以修正低估；它不是“增加真实测量”。因此不存在无条件的信息单调增加定理。

若 full A 与 A_n 的算子误差≤ε，则

$$\|H_F-H_n\|\le2\|A_n\|\varepsilon+\varepsilon^2=\eta.$$

在同维、匹配的弱谱簇与剩余谱间隙 g>0、η<g/2 的条件下，标准 Davis–Kahan 型界给

$$\|P_{w,F}-P_{w,n}\|\le2\eta/g.$$

没有谱隙时，小矩阵误差可以造成接近 90° 的弱子空间变化。低秩 Fisher 缺陷也不意味着弱投影之差一定同样低秩。[H1,H2]

## 9.2 材料移动后需要额外的误差项

current enrichment 控制同一 χ 的表示误差；材料移动改变真正 A(χ)。若对齐后的 prior frame 中

$$\|A(\chi')-A(\chi)\|\le L_A\,d_M(\chi',\chi),$$

则总扰动包括表示误差、材料物理变化和 frame 变化。L_A 可由 (3.3)、X 导数、内部场、resolvent 及 prior 导数界构造；它可能很大且昂贵。

A17 单个冻结 exchange 不会自动给出 L_A。只有 uniform 可逆性、材料本构正则、prior chart 正则和弱谱簇分离同时成立，才能用这些项控制弱投影随材料移动的变化。

不需要硬谱簇时，可使用已有的软弱滤波 λ(H+λI)⁻¹并做低秩更新；它避免 cutoff 跳变，但仍有 λ 控制的偏差与条件数，不会成为新的物理歧义定理。[H1]

## 9.3 新 action 应以有限提案为条件

对 Explorer 已经提出的 D，相关的当前误差不是完整 J 的平均误差，而是

$$R_UQ_V(D)B$$

经过接收端有限传播后的影响。候选 C 的精确已解析改变量是 (5.2)。这建议三个可测试的 action 目标：

- 对即将探索的材料方向，减少 false-weak 判别区间；
- 对有限 D，缩小 (5.5) 或更精细 dual bound；
- 对固定材料提案的 likelihood 接受/拒绝，消除决策不确定性。

一个具体评分是

$$\mathrm{score}(C\mid D)=\frac{[b_V(D)-b_{V+C}(D)]_+}{\mathrm{paidCost}(C)},$$

其中 b 是真实有效误差上界；若只有估计量，必须标为 heuristic。可以另乘已冻结的任务价值，但不能用 GT 提供在线方向。

## 9.4 没有免费的收益定理

所有候选的 basis 构造、Green 作用、正交化、Schur core、dual 更新及失败试探都要计费。仅看到 ||(5.2)|| 大不意味着误差减少；它也可能是 current 模型在两个不同方向上偏离 full。

一般非对称 Maxwell Galerkin 没有单调残差或子模收益保证。可在特殊 coercive 范数/特定近似条件下研究 greedy 误差下降，但不能无条件移植 SPD 最佳逼近结论。这里最重要的实验问题是：**按有限提案加权的 action 能否以更少总成本缩小真实决策误差，而不是是否能画出漂亮的 eigenvectors。**---

# 10. 主算法与不能省略的非 A17 对照

## 10.1 建议主方案：反馈可实现激发传输

暂称 Feedback-Realizable Source Transport，FRST。这是本文便于讨论的标签，不表示已经存在经过验证的新算法或原创优先权。

输入是合法材料先验及其可用密度/采样机制、完整 Maxwell 接口、观测与噪声合同、材料距离及总物理预算。起点来自任意合格 physics reconstruction；不必 N3，也不必 GN。若要 exact fiber，起点必须满足所定义的等数据目标；若做 noisy posterior，不必先投影到人为 noise ball。

**第一步：确定可探索的材料支持。** 从原 prior 允许的局部变化、物理相关区域或一个冻结的端口字典选 U；不得由测试 GT 选支持。原 prior 若根本不允许这种局部变化，应标不适用，而不是把它松成任意像素。

**第二步：建立端口接口。** 用精确 R₀U 得 P/O/H/B，或用 A17 current ROM 建立它们及端口残差 R_U。前者是强基线，不是只能给新方法享受的 oracle。端口 setup、缓存与验证全部计入。

**第三步：提出材料变化。** 优先测试第 6 节联合核及其有理图册。若只近似相容，用 (6.4) 直接评估有限响应，或在端口 T、D 上解受约束问题。目标要保留明确材料进度，例如最大化任务差异或与 archive 的分离，同时限制 likelihood 和原 prior；不允许通过回到原点轻易过关。

**第四步：选择 current action。** 按当前 D 下的 R_UQ_VB、finite bound 或 likelihood 判定不确定性，选择是否增广 C、换 chart、缩步或 full fallback。只需更改内部求解精度时，固定材料不动，这样不会改变 proposal law。

**第五步：完整验证与接收。** 所有报告为成功的材料端点，用独立完整 Maxwell 后端核对状态、观测、多源一致性和原 prior。若宣称整条路径可行，还要检查中间点或给区间界；端点通过不保证途中通过。

**第六步：更新与停止。** 同 U 内使用 (5.6) 或重新定 anchor；跨支持支付新 setup。长期位移低于任务尺度、prior 合法率极低、current 表示成本过高或 bounds 无效时停止该路线。

## 10.2 有理闭式条件不满足时，不要强行补一个大网络

可采用小型端口约束问题：

$$
\min_{T,D}\ \|T-T_{\rm proposed}\|_M^2
\quad\text{s.t.}\quad
T=D(B+HT),\quad\Phi(D)\le\Phi_{\rm cap},\quad D\in\mathcal D_{\rm prior},
$$

再加明确的材料进度或任务坐标。这里仍可能需要 SQP、Newton 或其他非线性 corrector，但计算使用有限端口反馈，不要求每步构造完整材料 SVD。

它与 CSI、DCTMC 和标准 nonlinear constrained optimization 有明显重叠。只有 A17 action/certificate 带来独立收益，才称为 A17-centered 方法；否则应诚实归类为已有物理表示上的一种实现。

## 10.3 四个严肃对照

**B0：Full-Jacobian + prior-native nonlinear corrector。** 使用同一材料先验、进度约束与精度；允许 factorization、adjoint、Krylov 复用。不能故意用每个方向重新分解的低效版本。

**B1：精确局部 T-matrix / Woodbury + 同一个 FRST 材料提案器。** 这条完全删除 A17，但保留有限物理机制。它是证明 A17 增量最重要的对照。

**B2：普通 goal-oriented primal/dual ROM 或 residual-driven Krylov。** 使用相同 current budget、同样的有限材料 RHS 与验证精度；用来判断收益是否只是一般 ROM 工程，而非 A17 的双侧 exchange 结构。

**B3：原 prior 坐标的训练自由 constrained exploration / local MH / tempered SMC。** 在低维 P0 上它们很可能更简单更快。概率比较必须固定同一 target；heuristic corrector 只能进入 candidate-discovery 组。

DCTMC 或 CSI 的忠实小型实现可作为 current-domain 补充对照，不应把“换一个变量名字”的方法当作从零开始的新 baseline。[L1–L3]

---

# 11. 命题 D：显式激发坐标与正确后验传输

## 11.1 单源标量材料的双向坐标

考虑固定 U、单标量照明，b∈Cˢ 且各分量非零。令 t 为端口等效激发，定义

$$e(t)=b+Ht,\qquad d_i(t)=t_i/e_i(t).\tag{11.1}$$

在 eᵢ(t)≠0 且 I−D(t)H 可逆的域上，

$$t=(I-DH)^{-1}Db$$

是逆映射。因此这不是先采样后投影的 many-to-one corrector，而是一个明确的非线性变量变换。同数据条件在该变量中为 O t=0。

**证明。** d_i e_i=t_i 直接给 t=D(b+Ht)；逆矩阵存在则 t 唯一。反向从 d 计算 t 后，非零 e 确保比值恢复 d。证毕。

材料映射的复导数为

$$D_td=\operatorname{diag}(e^{-1})(I-DH).$$

所以在完整 2s 维实化复材料空间中，

$$
\boxed{
|\det_{\mathbb R}D_td|
=|\det_{\mathbb C}(I-DH)|^2\prod_i|e_i|^{-2}.
}
\tag{11.2}
$$

这是 holomorphic 变量变换的实 Jacobian 公式。若 d 是非线性极化率增量，还需乘上实际本构逆的实 Jacobian；非 holomorphic 或低维约束不能直接套用平方模公式。

## 11.2 正确 noisy posterior

设材料端口外固定，p_D 是对应的条件 prior 密度，并在合法域内以独立复材料坐标为参考测度。则激发坐标中的 target 为

$$
\pi_T(t\mid y)\propto
p_D(d(t))\,|\det_{\mathbb R}D_td|
\exp\{-\tfrac12\|W\mathcal R[F_0+Ot-y]\|^2\}.
\tag{11.3}
$$

采用有明确正反向密度的 MH proposal 即可保持该目标。无需把每个 noisy 样本强行压回 O t=0。若只在固定暗子空间移动，还需要其他更新覆盖亮坐标，或明确只采样相应条件分布；单个非遍历暗核不等于全 posterior sampler。

这个坐标表述使原先隐藏在“物理 corrector”中的体积变化显式化，但仍不是 MH 的新不变性定理。它能否减少总物理成本取决于端口维数、域的条件数和 prior 在该坐标中是否更易探索。

## 11.3 多源怎么做，哪些东西不能省

可以从已有照明选一个固定 reference 或相干线性组合 b=Bc，用 (11.1) 给材料参数化，再对**全部原始照明**计算 T_p=Q_H(D)b_p 和原完整 likelihood。

这只让 reference 响应线性化，不意味着所有照明都满足 O t_p=0。不能把 reference 当作一份新增独立测量，也不能重复计算其 likelihood。随包小矩阵检查中，reference 的输出变化约 4.4×10⁻¹⁶，而其他源变化约 0.028，直接展示了这一危险。

三维各向同性的 arbitrary t 不一定对应共同标量 d；必须施加第 7 节块共线条件，或回到原生材料参数。固定相位、低维形状等先验在 t 中一般变成约束流形；此时使用完整 2s 维体积因子会改变 target。最安全的实现是保留原 prior coordinates，仅让 finite feedback 提供 proposal 的均值、预条件或求解精度。

## 11.4 A17 近似坐标仍可形成正确 proposal

用冻结的 H_V 构造 t↔d 的近似物理坐标，只要该映射确实可逆，其 Jacobian 就按**实际使用的 H_V**计算。再用完整 F(d) 计算 target，MH 可以纠正 proposal 偏差。不能在密度公式里用 H_full，却在实际生成时用 H_V。

这说明 posterior 正确性不要求 surrogate 完美；但错误的 surrogate 会影响接受率和混合。若在线改变 anchor/basis，必须能计算对应反向 proposal，或在正式采样阶段冻结这些选择。任意历史依赖的自适应不能仅靠“每个瞬时核看起来正确”保证整体正确。

## 11.5 一般原 prior 坐标的正确做法

对 q∈Rᵈ，明确的 target 为

$$\pi(q\mid y)\propto p_0(q)e^{-\Phi(q)}.$$

可以用 finite-compensation 结果设置 proposal mean，但最终取一个满秩 Gaussian 或其他显式密度：

$$q'\sim N(q+\mu(q),\Sigma(q)),\qquad\Sigma(q)\succ0,$$

并使用完整

$$\alpha=1\wedge\frac{p_0(q')e^{-\Phi(q')}q(q\mid q')}{p_0(q)e^{-\Phi(q)}q(q'\mid q)}.$$

反向 mean/covariance、log determinant、prior legality 的成本不能忽略。Gaussian latent 中 pCN 等 prior-reversible proposal 可以只需 likelihood ratio。

**没有显式 prior 密度，不一定完全不能采样。** 若存在严格定义的 prior law、原生 latent 密度或一个精确 prior-reversible 核，仍可正确运行。例如 independent prior proposals 的 prior 因子可在 MH 比率中相消。真正不够的是“有一个看起来合理的材料候选袋，却没有定义其生成概率”，这无法被声明为 Bayesian prior。

## 11.6 多对一 corrector、coarea 与全局概率

先加高斯噪声，再优化投影到相同数据集合，通常改变 proposal law；投影可能产生低维奇异分布或点质量。只保存随机 seed 也不会自动恢复可逆性。Augmented-state MH 需要真正一一映射的扩展状态、维数配平、逆操作与正确 Jacobian；constrained HMC/RATTLE 则需要完整约束几何和相应测度。[L9]

小噪声正规 fiber 的 Hausdorff 密度含 coarea 因子 p₀/√det(Df Df*)，而不是均匀材料弧长。这个结论已在 F/L/X 讨论，不是本次新增。[H1–H4]

候选图像在连续分布中通常各自概率为零。应报告明确材料邻域、结构/任务事件的 posterior mass。保留 MCMC 拒绝重复；candidate archive 可以去重，但不能用去重频率代替 posterior。局部区域内的归一化概率不等于该区域在全球 posterior 中的总质量。

## 11.7 有限反馈证书如何降低 full likelihood 成本

若 ||F−F_hat||≤δ，白化 residual r_hat=F_hat−y，则

$$|\Phi-\hat\Phi|\le\|\hat r\|\delta+\tfrac12\delta^2.$$

(5.5) 可以为这个区间提供 δ。用同一个 uniform 随机数，若接受/拒绝阈值与整个 log-ratio 区间分离，就可得到与准确 MH 相同的判定；否则继续 current enrichment 或 full solve。

没有确定性误差上界时，可用标准 delayed acceptance：第一阶段用冻结 surrogate 筛选，第二阶段用 full target 修正。DINO 已明确实现这种 derivative-informed geometric MCMC 与两阶段纠偏，所以这部分不构成 A17 独占贡献。[L7]

---

# 12. NN 是否需要，以及应该学什么

## 12.1 当前答案：不需要先上 NN

第 6 节闭式图册存在时，解析公式可能比学习更便宜、更稳。图册不存在时，网络也不能凭训练损失创造 prior 支持或消除多源相容障碍。因此先建立真实材料位移与成本收益，再决定是否学习。

最优先的学习对象是**action / step selection**：下一步用哪个已定义的材料提案、是否刷新 basis、加入哪组 current modes、当前误差界是否值得继续收紧。网络不改变材料 proposal 时，概率语义也更容易保持。

第二选择是反复出现、具有摊销价值的 finite compensation map，为小型 constrained solver 提供初猜。第三才是显式密度的局部随机转移或 global proposal。完整材料 GT 回归并非本任务的主对象。

## 12.2 无 GT 的 teacher 与训练目标

Teacher 来自固定原 prior、相同观测合同下的高精度有限 continuation 或端口约束求解，保存真实材料位移、T、row residual、current actions、所需精度与失败原因。GT 只用于隔离的最终风险评分。

可使用以下目标：有限 state residual、dual-weighted output residual、共享材料相容、目标材料进度、有限 secant 预测误差、action 总成本与失败率。No-op 的物理残差很小，因此必须有非零 progress teacher、最小材料距离或有效探索奖励；不能只训练 ||F(χ')−F(χ)||²。

学习 singular vectors 的坐标符号与近重根旋转不稳定；更适合学习 basis-invariant residual actions、投影算子作用、物理材料位移或明确 gauge 下的等变坐标。

## 12.3 稳定性不是一句“物理约束网络”

一个可证明的部署域至少需要：本构和 prior chart 光滑、相关 resolvent 有界、用到的局部场不接近零、有限 feedback pole 有安全距离、所选弱谱簇有间隙或采用软滤波，以及 approximation error 有界。

若在同一合法 chart 内，理想传输场 V 为 L-Lipschitz，学习场误差≤ε，且两条轨迹均不离开该域，则标准 Grönwall 给

$$\|x_\theta(t)-x(t)\|\le e^{Lt}\|x_\theta(0)-x(0)\|+\varepsilon(e^{Lt}-1)/L.$$

它不保证跨 prior 边界、极点或 rank-change 的泛化。可把网络限定为 proposal/初猜，并保留 certificate/full fallback；不能用万能逼近定理代替可验证的物理控制。

## 12.4 何时停止学习路线

若解析端口模型已经便宜、teacher 构建成本高、部署样本量小或跨材料域不稳定，就保留无 NN 方法。正确 break-even 是

$$N_{\rm deploy}>\frac{C_{\rm teacher}+C_{\rm train}}{C_{\rm baseline/valid}-C_{\rm learned/valid}},$$

分母必须为正，并在同样的有效位移、精度或概率误差标准下比较。cached decoder latency 不能替代完整成本。

---

# 13. 原创性审计：最接近的先例与真正剩下的空间

## 13.1 DCTMC：本次最关键的近邻

Levinson–Markel 的 DCTMC 原文给出 A(I−VΓ)⁻¹VB=Φ，并把 T=(I−VΓ)⁻¹V 作为线性数据相容对象；再通过材料 V 的局部/对角结构约束进行迭代。其计算也显式考虑 T↔V 的矩阵逆成本。因此，“数据在电流/T-matrix 中线性、材料约束非线性”已经是明确先例。[L1]

本次可能保留的差别是：固定合法 anchor 附近的**严格共享材料有限图册**、反馈行相容及其精确残差、A17 的双侧 current exchange 和有剩余项的误差控制。不能仅把 V/T 改名 χ/current 就宣布新原理。

## 13.2 非辐射电流与 CSI

Habashy–Oristaglio–de Hoop 1994 的摘要已明确讨论在初始辐射电流之后联合调整非辐射电流与材料；本次未取得完整证明正文，因此它是高度相关先例，不作未读细节推断。[L14]

CSI 通过数据方程与材料状态方程联合优化 contrast source 和 contrast；SOM/TSOM 则分离由测量较可靠确定的电流成分与需要进一步约束的成分。多源共享材料、状态方程与 regularization 并非新增概念。[L2,L3]

本研究的不同任务是有限多解发现，而非单一材料重建；但任务改名本身仍不足，需要本报告显式图册/误差控制和对应收益。

## 13.3 SOM-Net 与 current-to-permittivity 神经反演

SOM-Net 使用 deterministic current、BP 初值与物理展开，并对电流、材料、散射约束训练。它并不是只把测量喂给黑箱 CNN。因此不能把任何 current/material 交替或解析材料更新都说成新物理网络。[L3]

2025 年 arXiv:2506.21349 v1 将方法分为 current estimator 与 permittivity solver，按 transmitter 解析恢复并聚合材料估计，并使用 current、permittivity、scattering losses。本文不是第一个从电流解析读出材料。真正不同之处只能是严格多源相容的**有限材料集合与转移**，而不是监督单图预测。[L11]

该预印本后续 PDF 的标题有所变化；本报告具体方法比较使用明确的 v1，不混用不同版本的性能或理论主张。

## 13.4 Null-space / nonlinear data-consistent networks

Deep Null Space Learning 在线性 forward 的 null-space 中构造数据一致的 learned correction 并研究正则化收敛。非线性 data-consistent networks 又把经典重建与非线性一致映射结合。[L4,L5]

因此“重建后沿弱空间生成再修正”本身无法排除重合。本文有限激发图册的必要比较点是：它不是仅对 frozen material kernel 投影，而是给有限反馈与共享材料的显式关系。非线性网络论文的完整证明访问仍不充分，不能在未读条件上宣称已经严格超越。

## 13.5 e-GLASS：先验 latent 中的歧义探索已有直接先例

e-GLASS 在 GAN latent 中使用测量和感知几何选择扰动方向，并探索差异较大、测量变化较小的候选。其 forward 是线性测量，但 decoder 非线性；它不是 calibrated posterior sampler。[L6]

本研究若退回 prior Jacobian 弱方向加网络步长，没有显著新意。需要显示 Maxwell 有限反馈确实改变了可达范围、合法率、方向可靠性或总成本。

## 13.6 DILI、DINO、manifold lifting 与 global MCMC

DILI 已按 prior-preconditioned likelihood geometry 构造采样；DINO 已用 forward/derivative 学习、几何 proposal、reduced coefficient 操作和 delayed acceptance 降低昂贵 PDE posterior 成本。[L7,L8]

Manifold lifting 将噪声提升为约束流形并处理正确测度；局部 MCMC 与 normalizing-flow 全局提案的组合也已有先例。[L9,L12]

因此“J/H + NN + MH”不构成本文主张。A17 必须减少**有限 Maxwell proposal 的求解/判定成本**，而不是重复既有的不变性证明。

## 13.7 Born / resolvent / 2026 文献边界

电磁 Born 与 inverse Born 收敛、Padé 重求和已有专门研究。2026 年 inverse-Born neural operators 的可检索摘要又覆盖了基于高阶逆算子展开的学习思路。本次未取得其完整论文，不声称排除了全部重叠。[L10,L13]

“高阶反馈可以用 NN 学”不是安全的原创定位。可审计的有限多源材料传输、current exchange 误差和成本才是本次更具体的候选。

## 13.8 删除 A17，算法究竟失去什么？

若保留精确 P=R₀U，删除 A17 后，第 4、6、7、11 节的物理传输与概率变量变换仍成立。这是必须承认的事实：**A17 不是这些物理材料世界存在的原因。**

删除 A17 可能失去的是：用同一组 current factors 同时更新 O/H、按有限提案重用双侧 Schur 信息，以及在受控剩余误差下选择性改进内部求解。可是普通 primal/dual ROM、Krylov 或精确端口求解可能同样甚至更便宜。

所以本项目只有在这层计算和可靠性增量被证明后，才适合以 A17 为论文核心。若比较只损失少量工程细节或完全无变化，应放弃“A17 原生生成算法”的主张，而不是继续包装。---

# 14. 全部计算成本与关键失败模式

## 14.1 成本变量

设完整 current 维数 N，端口维数 s，source 数 n_s，receiver 数 m，current basis 维数 r。dense 后端的新材料分解通常 O(N³)，已有分解每 RHS 为 O(N²)；迭代后端要改用 Green matvec、预条件构造和每 RHS 实际迭代数，不能直接套用 dense 计价。

精确端口 setup 需要 s 个 R₀U RHS，以及 H/O 构造、P 的存储或可复用作用。每个有限提案的核心成本约为 s³+s²n_s+m s n_s；需要完整电流时再付 N s n_s。若 s 很小且同端口提案很多，可能比每个 χ 重新分解便宜；若 s≈N，基本没有这一优势。

A17 压缩可避免精确求全部端口响应，但每次加入 C 都要支付候选构造、L₀C/GC、正交化、Schur solve、H/O 和误差界更新。其实际成本由后端与已有缓存决定，不能仅报告小 Γ 的求逆时间。

## 14.2 Full material Jacobian 也可以便宜

同一个冻结 L 对所有源共享时，接收伴随 L^{-*}S* 可以跨源复用。以原 8-source/32-receiver 标量平台为例，完整 Jacobian 的接收侧构造可能只需 32 个伴随 RHS，而不是 8×32；加上 8 个 source states。低维 d=7 的 JVP 路线也只需约 8d 个 RHS，取更合适实现。[P3]

显式 512×1152 实 Jacobian 约占 4.5 MiB。这个规模下，“完整 J 太大所以只能 A17”并不成立。A17 必须在完整生成流程、多个状态或更大三维问题中证明成本收益。

## 14.3 探索成本应按成功的任务量归一化

建议报告

$$C_{\rm valid\ displacement}=C_{\rm total}/N_{\rm nontrivial,legal,validated},$$

以及达到固定材料覆盖、事件误差或 ESS 的成本曲线。分母为零时必须报告失败，不能转成高通过率的小幅 no-op。

总账包含 anchor、full states、source/adjoint RHS、basis、full J/H、finite feedback、失败和拒绝、reference、teacher 和训练。独立 full endpoint 审计也必须出现在研究总账中；可以另列部署算法必要成本与独立科学审计成本，但两者都不能消失。

如果每个候选最终都必须 full solve，有限端口模型仍可能减少中间 corrector 和失败提案成本；不能因此把最终验证算成零。反过来，精确端口公式在其假设内本来就是 full discrete Maxwell 的代数重用，部署时未必需要重复 LU；实验阶段仍用独立 solver 审计其实现。

## 14.4 最容易导致失败的机制

1. **先验不闭合。** 有理补偿离开固定相位或形状 prior，这是材料问题，不是 current basis 不够。
2. **多源过约束。** 初始 weak directions 存在，但 straight-source 联合核为空，或 finite row residual 很大。
3. **低 latent、高算子秩。** 全域材料变化让 s 接近 N，Woodbury 并不便宜。
4. **伪弱与伪极点。** ROM 在真实可见方向漏响应，或者 reduced feedback 奇异而 full 不奇异。
5. **证书过松。** 非正规 resolvent 的全局范数让 bound 失效，即使实际输出误差很小。
6. **field zero。** 材料比值病态，显式坐标 Jacobian 大；不能靠分母加 ε 隐藏。
7. **先验有宽度但没有重要多模式。** 极微小 noise-compatible 材料变化不等于结构世界发现。
8. **概率 law 缺失或变化。** 改 prior、profile nuisance、删重复、corrector 多对一，都可能改变目标。
9. **baseline 太弱。** 不允许 full-J 复用，或者漏掉精确局部 T-matrix，可能制造虚假 A17 优势。
10. **学习不摊销。** teacher/训练超出任何部署收益，或 learned policy 仅记住少数 anchor。

---

# 15. P1 未知时的条件决策

P1 的具体成功率、维度、prior law、最优性与物理预算尚未确认，本报告均不填造结果。

| 到达的证据 | 正确科学对象 | 下一步 | 不允许的结论 |
|---|---|---|---|
| A：丰富的 prior-supported 连续方向，且有真实有限竞争材料 | 原 prior 中的有限局部探索 | 比较联合核、misalignment、普通 corrector、精确端口与 A17 端口；检查最大合法位移与总成本 | “因为高维所以后验已经多模态” |
| B：大部分方向仍强可辨识，局部宽度小 | 有限噪声 posterior 宽度、负对照或主动测量 | 保留可靠单点物理反演；报告任务尺度上的不确定性，不强行生成结构 | “所有材料在全球唯一” |
| C：有分离的高支持结构，但局部路径不通 | 非局部 candidate/global proposal + 准确验证 | 全局模型或多起点负责跳跃，A17 仅尝试加速有限提案的 likelihood / 初始化 | “多阶反馈自动跨越 prior 不连通性” |
| D：有结构候选，但没有可信生成 law / density / reference | 候选发现与覆盖研究 | 先冻结 prior law 和噪声；仅做 candidate-only 验证，概率保持未知 | “候选 softmax 即 posterior mass” |

D 与 A/B/C 可以同时发生。缺少 density evaluator 但有精确 prior law 或合法 latent 并非绝对障碍；应检查可否使用 prior-reversible sampling，而不是先假造密度。

结果到达后只需填入 target manifest、合法坐标、rank/有限位移、prior law 和成本，再选择分支。第 5–7 节的条件代数不需要随 P1 重做。若实际 prior 不满足局部端口支持或有理图册闭合，就明确判不适用。

---

# 16. 下一轮最小、可证伪的实验

## 16.1 第一阶段：不新增物理的证据与代码审计

运行随包小矩阵程序；读取原缓存，不自动重做已付费物理。建立 target_id、prior_id、noise_id、coordinate_id、source packing、材料质量度量和 hash manifest。区分 known TRAIN geometry control、physics reconstruction、P1 新目标与 radial proxy。

检查旧图中哪些仅 box 合法、哪些真在原 prior；哪些只是 noise-compatible、哪些达到原最优性门槛。F2 原门槛未通过的状态不能被改名后绕过门槛。P0 七维 smooth 只作几何控制。

## 16.2 建议新预算，而非本次已运行预算

建议下一轮硬上限：**128 个新的完整材料 states、4096 个 complex RHS、GPU 训练 0、一般三维新 campaign 0**。既有费用不混入这份新额度，也不覆盖旧失败记录。

建议分配：E0/E3 最多 24 states/768 RHS；E1/E2 最多 80 states/2560 RHS；异常与独立审计预留 24 states/768 RHS。所有失败计入。若实际后端成本不适合该配置，先根据不读取测试质量的 microbenchmark 缩小整个配对矩阵，不按结果为某一臂加预算。

E4 概率 reference 和 E5 训练不在这份额度内，未经前置证据不得启动。

## 16.3 E0：current exchange 对材料弱几何

在同一材料和同一合法坐标，比较 full Jacobian、旧 ROM、A17 增广；同时保存 ΔH 的两侧低秩因子和 full−new 的剩余项。

指标包括弱投影误差、false-weak 数量、方向响应上下界覆盖、Fisher defect、finite secant 误差和全部构造成本。Full J 只作为 evaluator 或所有方法同享的基线，不允许免费给 A17 selector 读取。

若谱簇没有明确间隙，不比较单个 singular vector 的身份；改比整个子空间或软滤波作用。多个随机 probes 的效果不得标 CERTIFIED。

**反证目标：** A17 改善 frozen GN gap 却不改善 weak classification，或者同预算 full-J/Krylov 更好，均应直接报告。

## 16.4 E1：有限材料补偿，增加一个不可缺少的对照

在合格的未知物体 physics anchor 上，预注册 2 anchors×1 个不读 GT 的方向×2 个材料半径，组成 4 个配对块。每臂每块最多 3 个中间新状态加 1 个独立端点审计，避免不平衡地追加迭代。

五个主臂：

- frozen material-J weak/null step；
- refreshed full-J + standard nonlinear corrector；
- 普通 prior/material-native predictor–corrector；
- **精确局部 T-matrix + 与新方法相同的 finite proposal**；
- A17-guided finite transport。

每个配对块使用相同目标材料 RMS 距离、相同 prior 支持与 likelihood 要求。匹配初始材料距离，而不是相同 latent 参数范数。无法实现该距离记 FAIL，不静默缩到 no-op 后算成功。

比较完整预测误差、状态残差、prior/passivity、材料位移、结构/任务差异、重复率、最大可达范围与费用。所有成功候选 independent full-Maxwell 验证。未达到合格 anchor 门槛时，该块 NOT_RUN，不用 known truth anchor 替换后宣称未知成像成功。

**反证目标：** A17 未胜过 exact-port 或普通 full-J corrector，则拒绝其有限传输增量主张。

## 16.5 E2：针对探索的 current selection

在 E1 同一组材料提案上，用同成本比较 forward-residual、GN-update、weak-classification 与 finite-response/likelihood-decision selectors。可以复用相同 full evaluator 缓存，不额外创造有利场景。

把候选 current modes 的生成和所有被评估但未选中的 actions 都计入。核心指标是同预算下真实 finite error、可靠接受决策数及有效材料位移，而不是新增 basis 的数量。

**反证目标：** 如果新的评分只是更贵地算出同一个结果，或者所有 bounds 都空泛，就停止 adaptive enrichment 路线。

## 16.6 E3：多源可实现性

先做完整可验证小例子：单源标量 dark chart、两独立源的人工共线正控制、一般多源不相容反例、固定相位先验反例、场为零与 reduced 伪极点反例。

再在小型真实 TM 后端中检查原 prior 内是否存在非零联合核或可接受的 finite misalignment。三维矢量仅复用已存在可验证资产；没有相关资产则标 NOT_RUN，不用标量结果代称三维。

报告 joint-kernel 数值阈值与 dimension、各源完整数据变化、材料本构恢复、多源 minors、真实 source rank 与条件数。人工设计 B/H/O 的正例只支持代数，不构成实际天线配置中的歧义证据。

## 16.7 E4：概率 pilot 的双门槛

只有已经发现非平凡合法材料变化，且 prior/noise law 明确、reference 有可信误差说明时才运行。先用解析小 target 验证 proposal law，然后申请独立物理 reference 预算。

比较 local MH、full-Maxwell geometric proposal、explicit source-chart proposal、A17 proposal/筛选；C1 global density 仅在同 target 及密度合同成立时加入。普通 deterministic corrector 留在 candidate-only 表，不冒充采样臂。

报告 event mass/CDF 误差、重复状态、mode/region coverage、burn-in、ESS、跨区域 hitting 与总成本。不能用单个 ESS 替代所有稀有事件或跨模式性能。

## 16.8 E5：学习的启动门槛

多个独立合法 anchor 上存在重复的有限传输机制，解析版已胜过强基线，而且预计部署量足以摊销 teacher 成本时，才训练小型 action/compensation 网络。第一版参数量可限制在约 5 万以下作为工程约束，不把它当理论最优值。

否则停止训练建议。没有发现可用机制时，Diffusion/Flow 不是补救措施。

## 16.9 总表必须分别报告的五类结果

**Physics：** 完整状态/数据残差、真实有限 secant、多源一致性、总被动性、prior validity。

**Exploration：** 非平凡材料位移、有效覆盖、结构或任务事件、去重与重复率、可达半径；不能只报图数。

**Probability：** prior 支持与 law、相对 log likelihood、event mass、target 正确性、有限链误差、未覆盖质量。

**Risk：** physics-only Top-1、原 N3 Top-1、learned completion、oracle Best-of-K、可部署选择策略。Best-of-K 不冒充实际选择器。

**Cost：** 全部 source/adjoint、states、factorization、basis、J/H、feedback、失败、reference、teacher/train 和 wall time。共享资产成本规则在方法运行前冻结。

---

# 17. 直接给 Codex 的交接

完整可复制提示词在 `protocol/CODEX_NEXT_ROUND_ZH.md`。其任务是先审计、先运行合成检查、根据 P1 的实际状态分支，再在冻结预算内执行 E0/E3 和合格的 E1/E2；不是照单全跑 E0–E5。

应使用独立工作目录和新结果文件，不覆盖 A17 原论文、旧 F/LXZ 决策或 P1 运行。发布/提交仓库不是本报告默认动作。

---

# 18. 十个核心问题的最终回答

**Q1：A17 能否定义比普通 material SVD 更有意义的探索几何？**

能提供更可解释、对计算缺陷敏感的 injection–feedback–receiver 分解；进一步的有限端口 H 与 row-realizability 包含 frozen J 不具有的有限信息。但 full J 本来就包含一阶反馈，A17 不增加可观测信息。真正有希望的是第 5–6 节 finite structure，而不是把弱 singular vectors 改名。

**Q2：怎样区分真弱、ROM 假弱和 prior 非法？**

先限制到实际 prior；再用方向/子空间的 full−ROM 剩余响应界分类。只看 J_ROM v≈0 或一次 exchange 很小不够。缺少尾项界必须 UNKNOWN。Prior 不允许的方向直接排除。

**Q3：current exchange 能否控制材料移动时弱子空间变化？**

能控制固定状态下的已解析表示改变；再加未解析误差、材料 J 的路径变化界、prior-frame 变化和谱隙，才能控制移动后的弱投影。没有这些条件，答案是否定的。

**Q4：完整有限反馈能否形成非 GN 传输？**

可以。在反馈行相容的暗激发子空间中，dᵢ=aᵢ/(1+γᵢ(a)) 给精确有限图册。单源标量还具有直接 t↔d 坐标。一般多源、三维与窄 prior 下这些充分条件可能失败，届时仍需小型非线性约束求解或停止。

**Q5：能否不显式构造完整材料 SVD？**

可以使用端口等效激发、joint-kernel 作用、多源 minors、块矩阵相容和 matrix-free constrained solve。它不保证更低维：多源电流变量可能比材料更多，必须以压缩和真实成本判断。

**Q6：需要 NN 吗？**

目前不需要。优先解析图册、有限端口求解与误差控制。若确有可重复摊销收益，先学 current action/步长，后学 finite compensation 或显式密度 proposal；不回到 GT 单图补全。

**Q7：什么时候能连续探索，什么时候需要全局提案？**

在合法 prior 的正规非零 fiber、或具有任务尺度宽度的有限噪声区域内可以局部探索。局部满列秩、prior 边界、二阶障碍或分离支持可能阻断。是否全局多解要另证；global proposal 负责跳跃，A17 不自动搭桥。

**Q8：怎样结合正确 posterior sampling？**

使用实际 prior/noise law，选择显式变量变换 Jacobian、正反向 proposal density、prior-reversible kernel 或正确 augmented-state 方法；用 full likelihood、有效区间判决或标准 delayed acceptance 保持 target。不可把投影后的候选直接称 posterior。

**Q9：相比 full-J + standard corrector + learned proposal，新增什么？**

候选新增的是精确有理有限图册、可实现性误差分解，以及用同一 A17 factors 更新内部 H/接收 O 并控制有限响应的接口。但 exact local T-matrix 与普通 dual ROM 是必须比较的强替代。若它们以同成本达到相同效果，A17 的独立贡献未建立。

**Q10：能否成为强 IEEE TAP 论文？**

存在清楚的条件性理论主线，但现在不能据此认定论文已经充分。最关键缺口是：真实多源、原 prior 内的重要有限材料歧义；稳定便宜的残余证书；相对 exact-port/full-J/goal-ROM 的总成本增益；以及明确三维矢量或其他充分有说服力的电磁物理证据。

---

# 19. 论文定位、停止条件与转向

建议候选定位为：**Finite Material Transport and Task-Controlled Current Exchange in Nonlinear Electromagnetic Inverse Scattering**。它首先是物理算子与探索算法论文，不必带 NN，也不必一开始宣称完整 posterior generation。

可形成的贡献组合是：有限双侧 exchange 与结构化尾项界；多源可实现暗激发的有限材料图册及失效条件；以有限探索任务驱动的 current enrichment；同预算真实电磁验证。标准 Schur、隐函数、Davis–Kahan、coarea、MH 放在工具位置，不作为原创标题。

**停止 A17-centered 主线的条件：** 精确端口更简单更快；所有合法材料都需要近满秩端口；误差界长期空泛；新的 selector 不优于普通 dual residual；或有理族只存在于 off-prior 材料。此时可以保留无 A17 的有限材料探索，或转向可信不确定性/主动测量，不再加网络掩盖失败。

**停止多解学习的条件：** 原 prior 在任务尺度上近似可辨识；只得到 no-op 或 tiny perturbations；没有可信 prior law；全局不同结构仅是优化器伪影。保留 N3/physics 单点估计不代表失败，它们解决不同任务。

A17 原论文不需要被这条尚未完成的生成路线“补救”。可以在讨论中补充：原定理控制 current representation 对冻结材料步的影响，并不自动控制 weak geometry 或有限材料多解。没有必要把原论文改写成已经证明 posterior generation 的文章。

---

# 20. 本次实际完成的轻量核验

程序：`checks/check_finite_source_transport.py`，固定随机种子 20261017，依赖 NumPy/SciPy。所有矩阵是合成复矩阵，不是真实 Green 离散，也没有实际天线采集或材料重建实验。

| 核验 | 实际结果 |
|---|---:|
| A17 端口传播 exchange 相对误差 | 5.92×10⁻¹⁶ |
| 完整有限端口响应相对误差 | 1.48×10⁻¹⁴ |
| 双重修饰有限 exchange 相对误差 | 1.31×10⁻¹⁵ |
| reduced 端口与新差分方程 Galerkin 一致性 | 5.75×10⁻¹⁶ |
| 有限 state residual 因式分解相对误差 | 2.44×10⁻¹⁶ |
| 同端口 anchor 正向/反向恢复相对误差 | 约 1.13×10⁻¹⁶ |
| 单源标量 source chart 逆映射相对误差 | 1.62×10⁻¹⁶ |
| source chart 实 Jacobian 与有限差分相对误差 | 9.45×10⁻¹³ |
| source chart log-volume 公式误差 | 1.33×10⁻¹⁵ |
| 两独立源人工共线例的有理路径相对误差 | 2.79×10⁻¹⁶ |
| 联合核图册实际复维数 | 1（人工构造例） |
| 联合核图册相对误差 | 1.80×10⁻¹⁶ |
| 多源 misalignment 精确输出公式相对误差 | 2.59×10⁻¹⁴ |
| generalized sin²θ 一致性误差 | <10⁻¹²（精确值见 JSON） |

一个条件界数值检查中 θ≈0.007518，真实输出误差≈0.018840，上界≈0.106597，覆盖但不紧。这里 α₀ 由浮点 SVD 得到，属于公式实现检查，**不是区间算术认证的物理证书**。

反例同时核对：ROM 响应为 0 而 full 为 1；full 矩阵稳定而 Galerkin core 为 0；K=2 时 resolvent 存在但 Neumann 发散；只保持一个源输出时其他源发生变化；初始实方向经复反馈有理补偿后产生约 −0.03846 的虚增量，从而离开实/固定相位约束。

这些结果支持代数实现，不支持“已在 P1 找到新世界”“已实现正确概率采样”“已降低 Maxwell 总成本”或“已完成强 TAP 实验”。

---

# 21. 来源、版本与访问状态

## 项目原始证据

[P1] `migodam/a17-maxwell-certified-exchange`。直接读取 README 与 `closure_final/START_HERE.md`、`closure_final/manuscript/CORE_THEORY.tex`、`CONDITIONAL_EXCHANGE_THEORY.tex`、`PAPER_TAP_FINAL.tex`。原文定位以仓库路径为准；未在本报告中声明一个未经独立确认的固定 blob。原文是指定冻结 Galerkin/current 模型的定理，不含一般 deterministic full-output certificate。

[P2] `migodam/b1-probabilistic-maxwell-imaging`，commit `8b620b315b2b4b212d0f1e4c143ced49f22f02e0`，`research/research_f_intrinsic_transport/RESEARCH_REPORT_ZH.md`。已读取实际 F0/F1/F2 状态和主要数值；不声称重新执行。

[P3] `migodam/a26-om-feedback-prior-completion`，`A26_RESEARCH_REPORT_ZH.md`。直接读取主要实验/物理合同与负结果。

[P4] `migodam/a30-maxwell-observable-material-manifold`，commit `717d3d26a290248bd2b9aeb55f3d03694b407b84`，`A30_E0_E1_DECISION.md`。直接读取 oracle 容量、有限补偿及执行偏差。

[P5] `migodam/b1-probabilistic-maxwell-imaging`，commit `dde0237eaa0d8c5507fedb07bb704c66a817028d`，`research/research_lxz_20261016/RESEARCH_REPORT_ZH.md`。完整读取。不同目标、proxy/full、known TRAIN/未知重建的区别保留。

[P6] `migodam/c1-physics-state-generative-imaging`，`START_HERE.md`。直接读取入口和 Joint P2 限制；没有独立重读该仓库全部原始链。

## 先前研究资料

[H1] `Research_F_Full_Report_ZH.md`，Library 正式生成报告，1804 行，全文读取。已有有限 feedback、高阶导数、prior 内在曲率、多源 minors、概率与成本边界；这些不计为本次新发现。

[H2] `Research_L_Full_Report.zip` 中 `RESEARCH_L_FULL_REPORT_ZH.md` / PDF。取得正式全文文件，重点复核第 3–6 节与 A17 geometry 理论。不同于早期只包含脚本的 L algebra package。

[H3] `Research_Z_Physics_Context_Latent_Generative_Maxwell_Imaging.zip` 中正式报告及 `theory/THEORY_NOTES_ZH.md`。已读取；其架构和实验设计不转写成运行成功。

[H4] `RESEARCH_X_REPORT_ZH.md` 及四命题资料。重点读取 A17 fixed-material likelihood refinement、target/measure 和多解证据边界。对整段历史会话没有完整访问。

[H5] `TriSpace TSOM SLAM Research1.md` / Obsidian 版本中的广义特征值与 sin²θ 段落。作为先前推导使用；加权 prior/nuisance 情形不无条件沿用裸夹角结论。

## 外部原论文

[L1] H. W. Levinson, V. A. Markel. *Solution of the nonlinear inverse scattering problem by T-matrix completion. I. Theory*. Physical Review E 94, 043317 (2016), DOI 10.1103/PhysRevE.94.043317; arXiv:1401.3319。实际核对正文 Eq. (5)、(9)–(11)、迭代及计算部分。

[L2] P. M. van den Berg, R. E. Kleinman. *A contrast source inversion method*. Inverse Problems 13, 1607–1620 (1997), DOI 10.1088/0266-5611/13/6/013。核对官方出版信息、方法说明及后续原论文中的相应方程；没有声称完整重新审读全部 CSI 文献。

[L3] *SOM-Net: An Unrolling Neural Network for Inverse Scattering Problems*. arXiv:2209.03567 (2022)。读取原文 deterministic current、材料解析更新、训练损失与复杂度部分。TSOM 参考为 Zhong–Chen, Inverse Problems 25, 085003 (2009)，后续 FFT twofold SOM 为 IEEE TAP 59, 914–927 (2011)；对后二者本次未逐式重新审读。

[L4] J. Schwab, S. Antholzer, M. Haltmeier. *Deep Null Space Learning for Inverse Problems: Convergence Analysis and Rates*. Inverse Problems 35, 025008 (2019), arXiv:1806.06137。原文方法/理论可解析部分已读。

[L5] Y. E. Boink, M. Haltmeier, S. Holman, J. Schwab. *Data-consistent neural networks for solving nonlinear inverse problems*. Inverse Problems and Imaging 17(1), 203–229 (2023), DOI 10.3934/ipi.2022037，arXiv:2003.11253。官方摘要/图注与方法说明已读，完整证明访问未闭合。

[L6] A. Montanaro, D. Valsesia, E. Magli. *Exploring the solution space of linear inverse problems with GAN latent geometry*. ICIP 2022, arXiv:2207.00460。实际核对 Algorithm 1 与测量/感知几何。

[L7] L. Cao, T. O’Leary-Roseberry, O. Ghattas. *Derivative-Informed Neural Operator Acceleration of Geometric MCMC for Infinite-Dimensional Bayesian Inverse Problems*. JMLR 26 (2025)，论文编号 24-0745。实际核对 reduced proposal、Eq. (36)–(37) 与 Algorithm 1；不复述未视觉审核的性能表。

[L8] T. Cui, K. J. H. Law, Y. M. Marzouk. *Dimension-independent likelihood-informed MCMC*. Journal of Computational Physics 304, 109–137 (2016), arXiv:1411.3688。核对全文可解析的 prior-preconditioned LIS 方法。

[L9] K. X. Au, M. M. Graham, A. H. Thiery. *Manifold lifting: scaling Markov chain Monte Carlo to the vanishing noise regime*. JRSS B 85(3), 757–782 (2023), DOI 10.1093/jrsssb/qkad023，arXiv:2003.03950。本次核对官方/作者摘要与先前 F 的已读推导；没有逐条重新证明其全文。

[L10] K. Kilgore, S. Moskow, J. C. Schotland. *Convergence of the Born and inverse Born series for electromagnetic scattering*. Applicable Analysis 96(10), 1737–1748 (2017), DOI 10.1080/00036811.2017.1292349；另见 *Padé approximants of the Born series of electromagnetic scattering by a diffraction grating*, Physical Review A 109, 033522 (2024)。本次为出版摘要与已知结构核对，本文高阶公式自行推导。

[L11] Y. Cheng et al. *Generalizable Neural Electromagnetic Inverse Scattering*. arXiv:2506.21349v1 (2025)。直接读取 v1 HTML 的 current estimator、permittivity solver、跨源聚合与训练目标。后续 PDF 标题为 *Electromagnetic Inverse Scattering from a Single Transmitter*；不混用版本结果。

[L12] M. Gabrié, G. M. Rotskoff, E. Vanden-Eijnden. *Adaptive Monte Carlo augmented with normalizing flows*. PNAS 119(10), e2109420119 (2022), arXiv:2105.12603。本次 PDF 访问失败，比较限于取得的作者摘要与此前正式研究的引用，不声称核对全部证明。

[L13] J. C. Schotland, A. Titi, J.-N. Wang. *Inverse Born series based neural operators*. arXiv:2608.18262 (2026)。检索到摘要，全文访问失败；它是需要后续仔细核对的近邻，不据摘要宣称完整原创性排除。

[L14] T. M. Habashy, M. L. Oristaglio, A. T. de Hoop. *Simultaneous nonlinear reconstruction of two-dimensional permittivity and conductivity*. Radio Science 29(4), 1101–1118 (1994), DOI 10.1029/93RS03448。官方摘要明确包含非辐射电流与材料联合迭代；全文未取得。

**最终结论：本次找到了可证明、可失败、可实现的有限材料传输接口，但没有证实它已优于现有强方法。下一步最有价值的是把“原 prior 内存在真实有限传输”和“A17 以更低成本可靠计算它”分开做实验。前者失败，停止生成；后者失败，停止把 A17 放在论文核心。**