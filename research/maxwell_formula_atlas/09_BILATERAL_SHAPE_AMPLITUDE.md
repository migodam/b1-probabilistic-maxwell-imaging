# 09 · 双边传递、形状—幅度分离与按任务控制的物理计算

本章沿着“固定的两端—变化的中间材料”展开，不把标准恒等式当成新定理。R10 的实测问题是：这些公式能否减少完整非线性 Maxwell 的重复费用，或者生成更有用的材料世界？公式成立与算法胜出分别报告。

统一尺寸：ROI 阶数 $q=256$，接收端 $n_r=32$，源端 $n_s=8$。固定宿主消元后的 $C\in\mathbb C^{q\times q}$、$O\in\mathbb C^{n_r\times q}$、$E_R\in\mathbb C^{q\times n_s}$、$F_h$ 均固定；$D$ 是相对宿主的材料增量对角矩阵。时间因子、白化和原先验沿用 F-BASE-001—003。本章只对冻结离散系统作结论。

## F-BI-001 · 双边非线性传递

由 $h=D(E_R+Ch)$ 消去电流得到

$$\mathcal T(D)=(I-DC)^{-1}D=D(I-CD)^{-1},\qquad F(D)=F_h+O\mathcal T(D)E_R.$$

推导使用 $(I-AB)^{-1}A=A(I-BA)^{-1}$，只要求对应 resolvent 可逆，不要求 $D$ 可逆。接收端和输入端的维数小，不意味着中间状态在所有材料下都可以无损缩小。**资格：标准离散 T-operator 恒等式。** 实现：R9 精确 ROI solver；R10 `src/transfer_campaign.py`。

## F-BI-002 · 有限材料差的双边条件

$$\mathcal T(D_2)-\mathcal T(D_1)=(I-D_2C)^{-1}(D_2-D_1)(I-CD_1)^{-1}.$$

在等式左乘 $I-D_2C$、右乘 $I-CD_1$ 后，两侧均化为 $D_2-D_1$，给出直接证明。因此等观测要求

$$O(I-D_2C)^{-1}\Delta D(I-CD_1)^{-1}E_R=0.$$

这是有限非线性响应等价的必要且充分条件，但未知 $D_2$ 仍在左侧 resolvent 中；它没有自动变成固定线性核约束。材料必须仍属于原 decoder 支持。**计算含义：** 允许构造有限响应和误差分析；不能据此宣布找到 exact fiber。

## F-BI-003 · 固定形状的精确加权实现

原 P1 材料增量 $D=\delta A(z)$，其中 $A=\operatorname{diag}(a_i)$，$0\le a_i\le1$，$\delta\in[.8,1.8]+i[.03,.18]$。令

$$P=A^{1/2},\quad B=PCP,\quad L=OP,\quad R=PE_R.$$

则

$$F(z,\delta)-F_h=\delta L(I-\delta B)^{-1}R.$$

证明：$\mathcal T(\delta P^2)=\delta P(I-\delta PCP)^{-1}P$，可直接左右乘原方程验证，完全不需要 $P^{-1}$。所以零占据率 cell 不造成除零。**资格：精确离散恒等式，保留所有多次散射。** 形状改变后 $B,L,R$ 也改变；仅复幅度查询可以共用同一分解。

## F-BI-004 · 传递矩、可达与可观测状态

在 $|\delta|\rho(B)<1$ 的收敛域内，

$$F-F_h=\sum_{j\ge0}\delta^{j+1}LB^jR.$$

系数 $LB^jR$ 是输入输出矩。有限维实现的可达空间为 $\mathcal R=\operatorname{span}\{R,BR,\ldots,B^{q-1}R\}$；不可观测空间为 $\mathcal N=\cap_{j=0}^{q-1}\ker(LB^j)$。最小实现阶数是 $\dim\mathcal R-\dim(\mathcal R\cap\mathcal N)$，也等于充分长块 Hankel 矩阵的秩。这里的交集是理论定义，算法可用 QR/SVD 的投影实现，不要求显式集合求交。

解析延拓可将同一有理响应用于级数收敛域外的非极点位置；不能把级数截断在该域外当成有效近似。R10 的 64 阶结果只是被测试的压缩阶数，未认证最小阶数。**来源：标准有限维实现与 Krylov 理论；最近材料模态先例见 R10 查新表。**

## F-BI-005 · 固定宿主的双边材料 Gram

宿主处一阶变化为 $dF=O\,dD\,E_R$，第 $i$ 个材料列 $J_{(r,s),i}=O_{ri}(E_R)_{is}$，故

$$J^*J=(O^*O)\circ\overline{E_RE_R^*}.$$

证明逐项求和：$\sum_{r,s}\overline{O_{ri}E_{is}}O_{rj}E_{js}=(O^*O)_{ij}\sum_s\overline{E_{is}}E_{js}$。这比仅 receiver SVD 同时保留了输入照明，但仍是宿主一阶信息。若 $D_x=\partial d/\partial x\in\mathbb C^{q\times66}$，真实高斯 latent 上的白化 Gram 为

$$H_x=\frac2{\sigma^2}\operatorname{Re}\{D_x^*(J^*J)D_x\}.$$

**条件：** 参数是实坐标，材料按实 Hilbert 度量；复 Gram 不得直接当成实 Hessian。R10 的谱更贴近 Full-J 并未转化为更好的 GT 覆盖。

## F-BI-006 · 完整非线性灵敏度的同型结构

令 $L_D=O(I-DC)^{-1}$，$R_D=(I-CD)^{-1}E_R$。由 F-BI-002 的极限有

$$dF=L_D\,dD\,R_D,\qquad J_D^*J_D=(L_D^*L_D)\circ\overline{R_DR_D^*}.$$

固定两端宿主 Gram 与当前 Full-J 的差别恰在两侧内部反馈。它解释何时 Born 双边坐标失配，但不允许免费得到状态依赖的两侧场。**资格：标准 distorted-Born 导数结构。**

## F-BI-007 · 共同状态压缩的条件性障碍

考虑允许所有独立对角 $D$ 的材料族。若输入矩阵 $E_R$ 每行非零，则 $\operatorname{span}_{D}\{DE_R\}=\mathbb C^q$：用只保留第 $i$ 个对角元的矩阵即可产生第 $i$ 个坐标方向。因而不存在一个小于 $q$ 的固定状态子空间包含所有这些注入。

这只是充分反例：它限制所有材料通用的精确状态压缩，不否定固定形状的压缩，也不直接适用于受限的 P1 Gaussian64 decoder 像。更小的输入输出商仍要另审可观测性，不能由注入满张成单独排除。

## F-BI-008 · 交换子不是端口误差证书

$$[A,C]_{ij}=(a_i-a_j)C_{ij},\qquad \|[A,C]\|_F^2=\sum_{ij}|a_i-a_j|^2|C_{ij}|^2.$$

它衡量交换性；$[A,C]=0$ 仅说明 $C$ 不混合 $A$ 的不同特征值子空间，不保证共同对角化或低秩（例如 $A=I$、$C$ 为 Jordan 块）。大交换子可能主要落在端口不敏感方向。实际观测误差还需要 $L,R$、resolvent 稳定性和被舍弃状态的响应。**资格：恒等式加反推限制，不是新的普适误差界。**

## F-BI-009 · 固定形状 Schur 与摊销

复杂非 Hermitian $B=QTQ^*$，$Q$ 酉、$T$ 上三角。令 $\widehat L=LQ,\widehat R=Q^*R$，则

$$F-F_h=\delta\widehat L(I-\delta T)^{-1}\widehat R.$$

每节点只作三角求解。若 setup 为 $s$、直接节点费用 $c_d$、复用节点费用 $c_r<c_d$，同一计时边界下的盈亏点为 $M>s/(c_d-c_r)$。不可拿 warm 三角解与 cold 完整求解直接作全系统加速比。R10 精确256阶 Schur 的实测盈亏点21—58节点；压缩后的盈亏点另测。

## F-BI-010 · 复幅度解析导数

在 resolvent 非奇异处，因 $B$ 与其 resolvent 交换，

$$\frac{\partial F}{\partial\delta}=L(I-\delta B)^{-2}R,\qquad \frac{\partial^2F}{\partial\delta^2}=2LB(I-\delta B)^{-3}R.$$

这来自 $\partial[\delta(I-\delta B)^{-1}]=(I-\delta B)^{-2}$。对 $\delta=\alpha+i\beta$，实 Jacobian 两列为 $F_\delta,iF_\delta$。likelihood 包含复共轭绝对值，因而不是全纯函数；优化和概率积分仍在两个实属性维度上进行。

## F-BI-011 · 正确属性条件质量与 profile 的区别

原 Uniform2 经线性映射后密度为 $1/.15$，所以

$$\bar L(z)=\frac1{.15}\int_{.8}^{1.8}\int_{.03}^{.18}\exp[-\Phi(z,\alpha+i\beta)]\,d\beta\,d\alpha.$$

它不等于 $\exp[-\min_\delta\Phi]$。如果可交换导数和积分，

$$\nabla_z\log\bar L=-\mathbb E[\nabla_z\Phi\mid z,y].$$

保留 Gaussian64 先验后才得到形状 posterior。一个由 $y$ 选出的固定形状条件积分，不是 Full66 posterior。R10 对3个固定形状建立数值参考；没有为全部 $z$ 积分。

## F-BI-012 · 条件化能减少什么不确定性

对材料任务 $\varphi$，

$$\mathbb E[\varphi\mid y]=\mathbb E_{z\mid y}\mathbb E[\varphi\mid z,y],$$
$$\operatorname{Var}(\varphi\mid y)=\mathbb E_{z\mid y}\operatorname{Var}(\varphi\mid z,y)+\operatorname{Var}_{z\mid y}\mathbb E[\varphi\mid z,y].$$

精确处理属性只解决第一层条件计算；缺失形状分支仍影响外层质量与均值。Rao–Blackwell 化降低给定正确外层采样器的积分方差，不保证随机优化生成的形状具有正确 posterior 分布。

## F-BI-013 · 随机目标改变的条件协方差

在同一局部线性化，$H=J^TJ+I$，随机数据目标 $\eta\sim N(0,I)$、随机 prior 目标 $\xi\sim N(0,\Sigma_\xi)$ 独立，则 GN 更新中的随机部分满足

$$\operatorname{Cov}(\Delta x)=H^{-1}(J^TJ+\Sigma_\xi)H^{-1}.$$

完整随机目标 $\Sigma_\xi=I$ 时化为 $H^{-1}$。去掉两属性目标后 $\Sigma_\xi=I-P_{amp}$，差为 $-H^{-1}P_{amp}H^{-1}$，虽至多秩2却通过耦合影响全部坐标。修改目标协方差不等于保留原 posterior 的采样器。R10 的 bilateral trace-normalized covariance 是提案启发式，非新增物理信息。

## F-BI-014 · 非正规矩阵的损耗稳定性界

令 $K_\delta=\delta^{-1}I-B$、$\operatorname{Im}_H B=(B-B^*)/(2i)$，设 $\beta_B=\lambda_{min}(\operatorname{Im}_H B)$，则

$$-\operatorname{Im}_H K_\delta\succeq b_\delta I,\quad b_\delta=\frac{\operatorname{Im}\delta}{|\delta|^2}+\beta_B.$$

若 $b_\delta>0$，由 $|v^*K_\delta v|\ge b_\delta\|v\|^2$ 和 Cauchy–Schwarz 得

$$\sigma_{min}(K_\delta)\ge b_\delta,\qquad \|K_\delta^{-1}\|_2\le b_\delta^{-1}.$$

不要求 $B$ 正规，也没有假设 $\operatorname{Im}_HB\succeq0$。R10 的最小特征值实际略负，材料损耗仍使 $b_\delta$ 正。float64 特征值仅给数值资格；没有 outward rounding 就不能称严格机器证书。

## F-BI-015 · 全状态 residual 到 likelihood 误差

对近似解 $\tilde u$，$r=R-K_\delta\tilde u$，真实输出误差为 $LK_\delta^{-1}r$，所以

$$e:=\|W_\sigma(F-\tilde F)\|\le\frac{\sqrt2\|L\|_2\|r\|_F}{\sigma b_\delta}.$$

代码可用较保守但便宜的 $\|L\|_F$ 代替 $\|L\|_2$。若 $\tilde r_y=W_\sigma(\tilde F-y)$，

$$|\Phi-\tilde\Phi|\le\|\tilde r_y\|e+\tfrac12e^2.$$

这是从同一物理 residual 到任务判断的桥梁。R10 冻结门为 $e\le.01$ 且势误差上界 $\le.1$；不合格才 exact fallback。门只覆盖被评价的节点，不凭有限节点推断全先验支持。

## F-BI-016 · 可达 Krylov 压缩与不能省略的检查

取正交 $V_r$ 张成 $\operatorname{span}\{R,BR,\ldots,B^{p-1}R\}$ 的有效列，$H_r=V_r^*BV_r$，

$$\tilde F-F_h=LV_r(\delta^{-1}I-H_r)^{-1}V_r^*R.$$

若前 $p$ 个 Krylov 块完全包含在空间中，Galerkin 传递匹配相应低阶输入矩。这个局部展开性质不能独自保证所有合法复幅度的准确度，故实际算法仍计算完整 $r$，用 F-BI-014—015 决定是否复用或回退。压缩可减少内部解阶数，不减少原材料自由度；一次固定形状的 basis 不跨任意形状通用。

## F-BI-017 · 势误差到概率误差的条件

若同一可归一化先验支持上，$\operatorname{osc}(\Phi-\tilde\Phi)\le\Delta$，则标准 likelihood ratio 界给

$$\|\pi-\tilde\pi\|_{TV}\le\tanh(\Delta/4).$$

若只有 $|\Phi-\tilde\Phi|\le\varepsilon$，可用 $\Delta\le2\varepsilon$。有界实任务误差至多 $(\sup\varphi-\inf\varphi)TV$。必须同时处理积分截断、求积和遗漏支持；R10 节点级误差门尚不满足这个全支持前提。因此没有把节点通过率当作 posterior 认证。

## F-BI-018 · 被动性、互易性与有限观测的边界

对 $e^{-i\omega t}$ 标量 TM、实 $\mu$，$H=(\partial_yu,-\partial_xu)/(i\omega\mu)$，

$$\langle S\rangle=\frac{\operatorname{Im}(\bar u\nabla u)}{2\omega\mu},\quad \nabla\cdot\langle S\rangle=-\frac\omega2\operatorname{Im}\varepsilon|u|^2-\frac12\operatorname{Re}(\bar uJ).$$

若离散 $C^T=C$、$D$ 对角，$\mathcal T(D)^T=\mathcal T(D)$；这不要求矩形带标签 $OTE_R$ 自身对称。32个有限 receiver 不是完整边界功率积分。原 prior 已要求被动性，不能重复把它作为独立新观测。能量型稳定性用于误差控制有计算价值，但没有在本轮新增材料信息。

## 从公式到本轮决策

保留：F-BI-003/009/014—016 组合的固定形状响应复用，检查 setup 摊销与完整 residual。停止扩大：仅凭 F-BI-005 的宿主双边弱方向主张 GT 覆盖改善。保留为独立数值任务：F-BI-011 的属性积分。未启动：Full66 SDE、通用最小状态、跨所有形状的精确公共压缩。

完整推导、先例与实测来源：[R10 理论](../p3l_r10_bilateral_structure/BILATERAL_OPERATOR_THEORY.md)、[能量与误差](../p3l_r10_bilateral_structure/PHYSICS_ENERGY_AND_RECIPROCITY.md)、[10形状物理复用](../p3l_r10_bilateral_structure/REDUCED_RESPONSE_FRESH10.md)、[原创性审查](../p3l_r10_bilateral_structure/PRIOR_ART_AND_NOVELTY.md)。
