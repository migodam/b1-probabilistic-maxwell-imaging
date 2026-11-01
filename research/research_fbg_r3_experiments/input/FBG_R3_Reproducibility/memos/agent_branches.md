# FBG-R3-T 独立研究 C/E：物理可实现电流、有限分支与可控制的非物理通道

日期：2026-10-31。作者角色：并行新变量／多分支研究线。本文为整合报告的工作稿，不是新全波成像实验。

## 0. 独立结论

最值得继续的结构不是再命名一种 $T$-matrix，也不是宣称反馈参数变化自动枚举物理解。更有内容的结构是：**把共享实材料允许的电流流形，与数据允许的电流仿射空间明确分开；有控制地释放两者之间的法向自由度，再把候选校正回真实 Maxwell。** A17 可以决定哪些 current modes 对材料步重要、释放后损失多少材料曲率；Z-ALG 可以提供小块精确多根提案。两种作用必须分别证明。

本线形成两个可直接使用的条件结果：

1. **正定反馈度量无法移动孤立精确根。** 在固定观测、满列秩根附近，只改变材料相关正定 residual metric 的 continuation 所跟踪的精确根恒定。多条 beta 驻点曲线不等于多个材料解分支。
2. **两个同数据物理解之间存在显式的 relaxed-current 通道，并且本构缺陷可精确量化。** 线性 current 和 material 插值的缺陷为 $t(1-t)\Delta X\Delta Q$；profile 掉实材料后，缺陷由两端多源场的实夹角控制。这个结果定量解释了允许暂时非物理的 joint-current 方法为何可能绕开纯材料下降的障碍，但它不寻找未知另一端点，也不保证优化器经过该通道。

新颖性边界：A17-X 原完整报告第4、7节已给出有限端口、共享材料 minors、row least squares、缺陷传播和 current/material 对应。它们不是本轮新发现。本文新的推导重点是物理 current 流形的法向、finite chord 的显式障碍，以及这两点怎样限定分支搜索算法。一般双线性插值、隐函数定理、CSI、WRI、Schur 和代数根枚举均有先例，不能由本稿直接宣布优先权。

## 1. 已实际核对的来源与历史结果

主仓库固定版本：`6659aa49dff2535996d40962791e9a5a6e974adc`。已实际取得并阅读：

- `research/research_z_alg_r3_digits/RESEARCH_REPORT_ZH.md`。
- `research/research_z_alg_r4_basin/RESEARCH_REPORT_ZH.md` 与 `MATHEMATICAL_THEORY.md`。
- `research/research_z_alg_r4_diversity/RESEARCH_REPORT_ZH.md`。
- `research/research_z_alg_r2/RESEARCH_REPORT_ZH.md`、`MATHEMATICAL_THEORY.md` 与 `theory/general_quotient.md`。
- `research/p3l_r10_bilateral_structure/RESEARCH_REPORT_ZH.md`、`BILATERAL_OPERATOR_THEORY.md` 与 `followup_shape_reuse/SHORT_MECHANISM.md`。
- A17 并行研究员取得的正式 `A17X_ORIGINAL_FULL_REPORT.md`，重点复核第4–8节及原有算法和查新边界。

上述原件保存在 `fbg_r3_work/sources/` 及 `fbg_r3_work/a17_sources/`。本线未取得 FP-T R1 正式全文，相关原始来源状态为 **SOURCE UNAVAILABLE（本线）**；下述局部收缩公式独立推导，不冒称重读了该会话。主研究员可按其取得的全文更新来源状态。

必须保留的历史结果：

| 历史工作 | 已成立的结果 | 不能推出的结论 |
|---|---|---|
| R3-DIGITS | 普通自写 GN 清楚恢复4、6，Gram增加9；成熟TRF追加仍有高对比失败 | 九数字普遍恢复、所有Full-J差、多个archive点即多个根 |
| R4-BASIN | 精确有限有理搜索能去到较低数据损失区；新7残差约0.7625降到0.01542 | 成像改善；该例材料RMSE反从约5.05增到10.89，新三目标联合门均失败 |
| R4-DIVERSITY | 低中对比普通成熟Full-J可正常成像；共360项返回 | GlobalDual整体精度或成本占优 |
| Z-ALG R2 | 条件4材料Q6和6材料Q20；GENERAL6固定64起点漏掉一个合法根 | 任意576材料多项式规模全根枚举、monodromy群已求出 |
| R10 | 固定形状rank64，10形状×32节点320/320通过，含setup的32查询工作负载约51.82%节时 | 未知形状反演全过程同样加速 |
| R10追加 | 旧加权状态基跨形状0/192通过，实际回退后约2.249倍直接LU时间；20/24形状新支持激活 | 所有固定TSOM/current bases都无效；另4个无支持激活的形状也失败，支持激活不是唯一解释 |

原文链接统一形如：
`https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/6659aa49dff2535996d40962791e9a5a6e974adc/research/<path>`。

外部一手资料实际打开：

- Levinson–Markel, DCTMC theory, [作者全文](https://arxiv.org/pdf/1401.3319), Phys. Rev. E 94, 043317 (2016), DOI 10.1103/PhysRevE.94.043317。正文允许中间非局部材料，往返于数据相容 $T$ 与近对角材料；一般非线性全局收敛并未由其线性极限分析自动给出。web refs `turn16view1`, `turn18view4`。
- van Leeuwen–Yang, [作者全文](https://arxiv.org/html/2403.15292v2), [作者机构正式发表记录](https://ir.cwi.nl/pub/34830/), Inverse Problems 41(2), 2025, DOI 10.1088/1361-6420/adab89。约束松弛消元后得到参数相关 Gram 度量，直接限制 GlobalDual／新weighted residual首创主张。web refs `turn16view0`, `turn18view2`, `turn16search2`。
- van Leeuwen–Herrmann, A penalty method for PDE-constrained optimization in inverse problems, DOI 10.1088/0266-5611/32/1/015007；作者预印本 https://arxiv.org/abs/1504.02249。此处核对出版/作者摘要，未声称全文逐行读完。

## 2. 固定支持的精确有限电流变量：有用但已存在

令 $X_0=\operatorname{diag}(m_0)$，$L_0=I-X_0G$ 可逆，$j_0=L_0^{-1}X_0E$，$Q_0=E+Gj_0$。$U_A\in\mathbb R^{n\times k}$ 必须是实际活动单元的坐标列；一般稠密 decoder 的 $k$ 个系数不意味着rank-$k$材料改动。

定义

$$V=L_0^{-1}U_A,\qquad C=U_A^*GV,\qquad E_A=U_A^*Q_0,\qquad O=SV.$$

对于 $\Delta X=U_ADU_A^*$、$D=\operatorname{diag}(d)$，有

$$\eta=D(E_A+C\eta),\qquad j_D=j_0+V\eta,$$
$$\boxed{F(m_0+U_Ad)=F_0+O(I-DC)^{-1}DE_A.}$$

证明：两物理方程相减为 $L_0(j_D-j_0)=\Delta X[Q_0+G(j_D-j_0)]$；令括号的活动单元本构激发为 $\eta$，得到以上方程。矩阵行列式引理给出

$$\det L_D=\det L_0\det(I-DC),$$

故在该固定支持问题上，小resolvent与完整resolvent的可逆性等价。

该公式就是A17-X局部T-matrix、Z-ALG和R10已有的精确结构。本轮只把它作为有限branch与法向控制的计算载体。

### 2.1 共享实材料的profile与有限缺陷

任意 $\eta\in\mathbb C^{k\times p}$ 不一定对应合法共享材料。记 $q=E_A+C\eta$。在每行 $q_i\ne0$ 时，无正则实材料最小二乘为

$$d_i(\eta)=\operatorname{clip}_{[l_i,u_i]}\frac{\Re\langle q_i,\eta_i\rangle}{\|q_i\|^2},\qquad e=\eta-Dq.$$

这里 $[l_i,u_i]=[m_{\min}-m_{0i},m_{\max}-m_{0i}]$；要求总材料合法，不要求增量本身为正。若含耦合TV/形状先验，逐行clip不再是整个带先验内层问题的精确解。

对该读出 $D$，只要 $I-DC$ 可逆，

$$\boxed{F(D)=F_0+O\eta-O(I-DC)^{-1}e.}$$

证明：精确物理 $\eta_D$ 满足 $(I-DC)\eta_D=DE_A$，而候选满足 $(I-DC)\eta=DE_A+e$，相减即得。于是完整白化数据的误差桥为

$$\|W_y\operatorname{vec}(F(D)-Y)\|\le\|W_y\operatorname{vec}(F_0+O\eta-Y)\|+\|W_y\operatorname{vec}[O(I-DC)^{-1}e]\|.$$

这不是本轮新公式：A17-X第7节已经给出全电流版。本轮固定支持的价值是所有修正都只需 $k$ 阶求解；但当小核已精确时，直接评价真实物理响应本来同样便宜，不能把这个界包装为无条件计算节省。

### 2.2 物理current流形的法向

在一个 $e=0$、$q_i\ne0$、box不活动的物理点，定义实正交行投影

$$P_{q_i}^{\perp,\mathbb R}v=v-q_i\frac{\Re\langle q_i,v\rangle}{\|q_i\|^2}.$$

令 $A=I-DC$。profile缺陷导数为

$$\boxed{De[\dot\eta]=P_q^{\perp,\mathbb R}A\dot\eta.}$$

证明：$\dot e=A\dot\eta-\operatorname{diag}(\dot d)q$。最小二乘的一阶最优条件在物理点给 $\dot d_i=\Re\langle q_i,(A\dot\eta)_i\rangle/\|q_i\|^2$。代回得式。

因为 $A$ 可逆，每行投影实秩为 $2p-1$，故缺陷导数总实秩为 $k(2p-1)$。常秩定理给出：物理电流集合在该点附近是实维 $k$ 的嵌入流形，切空间为

$$\boxed{\dot\eta=A^{-1}\operatorname{diag}(\dot d)q,\quad\dot d\in\mathbb R^k.}$$

这把可物理实现的材料切向与额外 $2kp-k$ 个自由current法向区分开。后者不是新的物理材料自由度；它们是CSI/WRI暂时松弛物理时使用的辅助方向。

数据切向在该流形上是

$$J_A\dot d=W_y\operatorname{vec}_{\mathbb R}\{OA^{-1}\operatorname{diag}(\dot d)q\}.$$

若 $J_A$ 在邻域恒秩 $r_A$，同数据物理fiber局部维数为 $k-r_A$。若满列秩，局部没有可连续探索的非零exact fiber；不能因为ambient current空间很大就创造多解。

节点、活动box、复材料及张量材料要另论：复独立标量材料的物理维数是 $2k$；三维各向同性需把三极化与所有照明同放一行；任意张量材料不适用这里的标量row投影。

## 3. 定理一：反馈度量continuation不移动孤立精确物理解

令 $r(m)\in\mathbb R^N$ 为固定观测的完整白化残差。设

$$f_\lambda(m)=\tfrac12r(m)^TW_\lambda(m)r(m),\qquad W_\lambda(m)\succ0,$$

其中 $r,W$ 二阶光滑，$m_*$ 是box内部精确根，$r(m_*)=0$，$J_*=Dr(m_*)$满列秩，且 $W$ 在所讨论的参数区间保持正定。

**结论。** 对每个 $\lambda$，$m_*$都是同一非退化局部极小值；任何经 $m_*$ 的局部光滑stationary continuation分支恒为 $m(\lambda)=m_*$。

**证明。**

$$\nabla_m f_\lambda(m_*)=0,$$
$$\nabla_m^2 f_\lambda(m_*)=J_*^TW_\lambda(m_*)J_*\succ0.$$

所有含 $D_mW$ 的梯度或Hessian附加项都至少含一个残差因子，在根处消失。因此 $\partial_\lambda\nabla_m f_\lambda(m_*)=0$。隐函数定理给stationary branch局部唯一，而恒定分支已存在；或直接有

$$\frac{dm}{d\lambda}=-H^{-1}\partial_\lambda\nabla f=0.$$

证毕。

**含义。** 在固定观测下，beta或反馈图册变化能改变非零残差目标的驻点、吸引域、鞍点及box活动面，但不会凭空把一个孤立精确根搬成另一个。把多个beta初值的不同末端称为“stationary branches”可以；要称为物理解分支，必须分别通过原物理数据门，并证明它们的材料区别。局部极小端点不拟合数据时，不属于等测量材料fiber。

若有正则 $\mathcal R(m)$，精确GT不一定是驻点；定理不能直接用于正则化极小点。若Gram奇异或采用截断伪逆，正定条件失败，可能改变零集合，这属于另一个问题，不能当作成功枚举物理解。满列秩不成立时，正则fiber的切向单位模态仍存在，但改变metric并未定义沿fiber的探索动力学。

## 4. 定理二：两个物理解之间的有限current通道及精确障碍

此节是**已知两个端点时的存在性与几何定理**，不是寻找第二端点的算法保证。

设两个合法实对角材料 $X_a,X_b$ 对同一入射场 $E$ 产生电流 $j_a,j_b$，总场 $Q_a=E+Gj_a$、$Q_b=E+Gj_b$，且

$$j_a=X_aQ_a,\quad j_b=X_bQ_b,\quad Sj_a=Sj_b=Y.$$

定义

$$j_t=(1-t)j_a+t j_b,\quad Q_t=(1-t)Q_a+tQ_b=E+Gj_t,$$
$$X_t=(1-t)X_a+tX_b,\quad0\le t\le1.$$

### 4.1 不作profile的精确通道

因为观测对current线性，$Sj_t=Y$ 对全路径成立；因为材料box凸，$m_t$ 始终合法。原Maxwell本构缺陷为

$$\boxed{e_t=j_t-X_tQ_t=t(1-t)(X_b-X_a)(Q_b-Q_a).}$$

**证明。** 展开 $X_tQ_t$ 的四项，端点项用 $j_a=X_aQ_a,j_b=X_bQ_b$ 消去，剩下 $t(1-t)[X_aQ_a+X_bQ_b-X_aQ_b-X_bQ_a]$。即为所示乘积。

对不含材料正则的标准lifted目标

$$\mathcal J_\rho(m,j)=\tfrac12\|Sj-Y\|^2+\tfrac\rho2\|j-X(E+Gj)\|_F^2,$$

有

$$\boxed{\max_{t\in[0,1]}\mathcal J_\rho(m_t,j_t)=\frac\rho{32}\|(X_b-X_a)(Q_b-Q_a)\|_F^2.}$$

所以当 $\rho\downarrow0$，两个同数据物理解在增广变量空间中存在任意小的这个目标障碍。这个路径不需要求物理逆，因此中间材料resolvent即使奇异，lifted路径仍代数有定义；但物理corrector和材料端点验收必须单独排除奇异。若有材料先验，其路径最大正则代价必须另加，不能靠 $\rho\to0$ 消去。

该结论解释的是一个真实的搜索空间扩展机制。它不等于目标全局凸化，不保证所有驻点消失，不保证优化器找到通道，也不提供未知第二端点。单纯可逆材料重参数化无法获得这种“离开物理current流形”的路径。

### 4.2 精确profile后的多源场角障碍

对每个cell，令端点row为 $q_a,q_b\in\mathbb C^p$，材料为实数 $m_a,m_b$；$j_t=(1-t)m_aq_a+t m_bq_b$、$q_t=(1-t)q_a+tq_b$。若 $q_t\ne0$，定义

$$\widehat m_t=\frac{\Re\langle q_t,j_t\rangle}{\|q_t\|^2},\qquad\widehat e_t=j_t-\widehat m_tq_t.$$

则

$$\boxed{\|\widehat e_t\|^2=\frac{t^2(1-t)^2(m_b-m_a)^2}{\|q_t\|^2}\left[\|q_a\|^2\|q_b\|^2-(\Re\langle q_a,q_b\rangle)^2\right].}$$

**证明。** 实最小二乘给

$$\|\widehat e_t\|^2=\|j_t\|^2-\frac{(\Re\langle q_t,j_t\rangle)^2}{\|q_t\|^2}.$$

令 $a=\|q_a\|^2,b=\|q_b\|^2,c=\Re\langle q_a,q_b\rangle$，将两条仿射row代入并展开，分子为 $t^2(1-t)^2(m_b-m_a)^2(ab-c^2)$。实Cauchy–Schwarz保证 $ab-c^2\ge0$。

亦可用复Lagrange恒等式分解：

$$\|j\|^2\|q\|^2-(\Re\langle q,j\rangle)^2=\sum_{s<t}|j_sq_t-j_tq_s|^2+(\Im\langle q,j\rangle)^2.$$

前一项是跨illumination比例不一致，后一项是实材料相位不一致。复材料时后一项不需要，但实材料不能删掉。

若 $c\ge0$，则

$$\widehat m_t=w_a(t)m_a+w_b(t)m_b,$$
$$w_a=\frac{(1-t)[(1-t)a+tc]}{\|q_t\|^2},\quad w_b=\frac{t[(1-t)c+tb]}{\|q_t\|^2},$$

满足 $w_a,w_b\ge0,w_a+w_b=1$，故profile材料仍在两端材料之间，box合法。若 $c<0$，这个充分保证失效，必须实际clip；clip后的缺陷不小于无约束profile缺陷，不能仍使用等式冒充box内最优值。

更精确地，若 $m_c=\operatorname{clip}(\widehat m_t)$，实LS正交性给

$$\|j_t-m_cq_t\|^2=\|\widehat e_t\|^2+\|q_t\|^2(m_c-\widehat m_t)^2.$$

因此box裁剪的额外障碍可以显式加回。若 $q_t=0$，ratio及场角除法式无定义，真实最小缺陷直接为 $\|j_t\|$；不能把分子同时为零当作物理相容。无约束复材料可删去相位项，实材料不可以。上述范数在unitary源基变换下不变；一般可逆源重组需连同源噪声／功率内积一起变换。

若两端row范数同为 $q$ 且实夹角余弦 $\cos\theta=c/q^2\ge0$，路径最大profile缺陷简化为

$$\max_t\|\widehat e_t\|^2=\frac{(m_b-m_a)^2q^2(1-\cos\theta)}8,$$

最大值在 $t=1/2$。对多源场，决定通道障碍的是幅值、材料差及**实多源场夹角**，不是单独减小current tangent范数。

### 4.3 与A17-X的关系

A17-X的精确straight-current图册要求各cell反馈场变化与原多源row相容，实材料还要求正确实相位。上式显示同类约束的有限几何：若材料发生变化而两端row不实共线，直线current一般离开物理流形；若满足相容，则profile缺陷为零。

这不是宣布所有材料解可由原A17-X有理图册连接。恰相反，它明确指出什么时候必须暂时释放非物理法向，或者使用弯曲物理路径。

## 5. 一个完整两单元代数例：两孤立物理解与通道

取

$$G=\begin{pmatrix}-.06+.02i&.065+.007i\\.065+.007i&-.15+.02i\end{pmatrix},\quad E=(1,.7)^T,\quad S=(1,1.3),$$

材料box为 $[.02,20]^2$。$G^T=G$，其Hermitian虚部正定；但**没有声称该矩阵正好来自用户的圆盘Hankel几何**。它是精确复系数代数VIE，用来检查逻辑和公式。

令真实材料 $(2,7)$，响应

$$y=5.792448701397235+0.5010146107905695i.$$

写 $a=G_{11},d=G_{22},b=G_{12}$，$A_0=S_1E_1,B_0=S_2E_2$，$C_0=b(S_1E_2+S_2E_1)-dA_0-aB_0$。由二维逆公式，数据方程等价于

$$A m_1+B m_2+C m_1m_2-y=0,$$

其中 $A=A_0+ya,B=B_0+yd,C=C_0-y(ad-b^2)$。消去

$$m_2=\frac{y-Am_1}{B+Cm_1},$$

其为实数的条件是实二次多项式

$$\Im(-A\overline C)m_1^2+\Im(y\overline C-A\overline B)m_1+\Im(y\overline B)=0.$$

本例系数约为 $(-.0242418611053,.158556850286,-.220146256151)$。两根读出

$$m^{(a)}=(2,7),\qquad m^{(b)}\approx(4.5406220091,2.0214837292).$$

两者原响应残差约 $10^{-16}$，最小物理resolvent奇异值约1.046和1.080，实Jacobian最小奇异值分别约0.00958和0.01989。因此两者均为局部孤立根；beta正定metric continuation从其中一根不能移动到另一根。

中点线性材料/current缺陷范数为约0.42086；实材料profile中点材料约 $(3.13932,3.80423)$，缺陷降为约0.01519，且在端点box间。精确弦线式相对双精度误差约 $10^{-15}$，profile场角式约 $10^{-13}$。

加一条真正不同入射场 $E_2=(.2,1)^T$ 后，两材料预测分别约 $5.39402+.47505i$ 与 $3.24013+.23685i$，此额外数据可区分它们。将原唯一照明乘非零标量不产生同样新信息。

数值文件：`branches_checks.py` 与 `branches_checks.json`。另以9维随机复对称模型检查小块有限缺陷身份，误差 $3.74\times10^{-16}$；检查profile缺陷法向导数，中央差分相对误差 $5.12\times10^{-11}$。这些均为代数检验，不是新成像性能证据。

## 6. 三种真正不同的算法六元组

### C/E-1：受控法向释放的joint-current反演（优先）

**变量：** 全材料 $m$ 与current系数 $a$，$j=j_D+U_ca$；保留完整材料自由度。

**约束：** $m$在原box／原先验；状态约束以有限正罚项允许暂时违反；最后使用完整物理求解验收。

**目标：**

$$\min_{m,a}\tfrac12\|W_y(S(j_D+U_ca)-Y)\|^2+\tfrac\rho2\|W_c[j-X(E+Gj)]\|^2+\mathcal R(m).$$

**更新：** 固定m时a是线性最小二乘；固定a时实材料box是行LS，含TV则解明确的带box二次/凸子问题。也可真variable projection消去a；必须包括内层解的正确梯度。A17另提供物理tangent误差与joint自由current的曲率损失检查，避免把两个Schur对象混淆。

**子空间：** 从原FFT-TSOM APCS开始；新增atoms先按物理tangent／normal分解，比较A17任务相关增广与同成本residual-Krylov。不能仅按current能量排名。

**Globalization：** 在物理残差和数据残差双尺度上调整rho；遇停滞可短暂释放选定normal modes。分别保留少量不同material端点，再同预算完整物理corrector；最终按原数据+相同先验排序。

**归因：** 整体仍属CSI/WRI近邻。独立贡献的候选是有任务条件的normal释放规则、profile曲率保留门、finite bridge角障碍与真实成像增量。仅交替LS、ALM或rho continuation不够新。

### C/E-2：固定小支持的代数多分支提案（独立次优）

**变量：** 活动材料增量 $d\in\mathbb R^k$ 或等效激发 $\eta\in\mathbb C^{k\times p}$；outside材料在本次提案阶段冻结，之后全材料corrector释放。

**约束：** $\eta=D(E_A+C\eta)$、真实box、$I-DC$可逆。

**目标：** 提案阶段解选定测量的精确或噪声容忍兼容性；保留多个合法根，随后用全部测量的原目标corrector。原正则必须同样提供给普通baseline。

**更新：** $k=4,6$时，从已有测量选一个接收组合及 $s=k/2$个照明，使对应source kernel满秩；调用已成立的Q6/Q20条件代数核。其余测量用于筛查/排序；带噪时只把方程根当proposal，不能称原全数据exact根。

**子空间：** 真正坐标支持 $U_A$。选择依据只使用公共位置、当前材料、原观测、A17双侧材料步相关性、残差与预定空间字典。不能把GT支持或错误图像的幅度方向固定为唯一候选。

**Globalization：** 按少数空间不同的支持枚举，每块保留所有或预定限额的合法分离候选，再做共同全材料corrector；明确保留未经校正而未选的分支记录。支持替换让各轮可以改变全部材料，不等于一次冻结所有outside已知。

**核心条件包含性。** 若某个真实完整数据根与当前材料只在该支持内不同，且选定通道的Q6/Q20先决条件成立、根为简单非节点非极点并全部枚举，则它一定在原始proposal集合中。证明：该根满足选定测量与实材料兼容式，条件商代数完整列出该有限系统全部根，合法性检查不会删除它。若候选限额截断、浮点root过滤或outside错误，这个包含性保证不再成立。

**成本：** 参考完整LU已付费时，$V=L_0^{-1}U_A$需 $k$个完整RHS；形成 $GV,O$与核也需计费。一个proposal只需 $k$阶LU和 $p$个小RHS，原始数据评分全部可由精确小核计算。改变支持／锚点通常需更新核。Q6/Q20代数装配、特征求解、失败回退和最终full corrector不是免费。普通block-TRF和多起点必须共享同一个小核。

**风险：** 四/六像素支持可能太小，未包含有意义的shape变化；多源/多接收全数据使局部exact方程往往不相容；投影测量会制造许多假branch；Q维数随k组合增长。不是继R4之后“更多uniform rays”的替代命名，应先用16/36像素真实TM小网格测能否形成空间结构变化。

### C/E-3：真正的多分支shape–property profile（限定场景）

**变量：** 只有在原材料族确为 $D=\delta A(z)$时使用z与复属性delta。

**物理：** R10精确有理响应 $F(z,\delta)-F_h=\delta L_z(I-\delta B_z)^{-1}R_z$。

**目标：** 对每个z保留属性目标的若干局部极小分支 $\delta_b(z)$，而非一次GN的2×2线性消元。分支目标 $\phi_b(z)=\mathcal L(z,\delta_b(z))$。

**更新：** 有条件内层优化，外层用envelope gradient；若内层Hessian可逆，

$$\nabla_z^2\phi_b=\mathcal L_{zz}-\mathcal L_{z\delta}\mathcal L_{\delta\delta}^{-1}\mathcal L_{\delta z}.$$

内层分支相交时全局lower envelope不光滑，应保留分支身份，不把其导数当单个光滑函数。

**子空间与Globalization：** 固定shape多属性查询才复用R10核；每次shape改变后做原residual门控，不能冻结旧rank64。外层多个结构proposal来自与普通baseline相同的先验。

**归因与优先级：** variable projection和R10响应都已有。可能价值在多属性分支不提前丢弃及同shape复用；但对本轮576独立实材料会引入额外结构假设，故优先级低于前两条。

## 7. 分支术语与continuation必须分清

| 对象 | 数学定义 | 本轮可以怎么用 |
|---|---|---|
| 精确物理解 | $F(m)=Y$且材料合法 | 分别完整验收；孤立根可多个 |
| 正则fiber | 常秩映射的等值流形 | IFT局部切向/法向校正；不能由有限噪声近似代替 |
| 噪声相容区域 | 白化残差处于固定阈值内 | 可以厚、弯曲或连通；不自动等于posterior mass |
| 非零残差驻点分支 | $\nabla f_\lambda(m)=0$ | beta/metric continuation主要跟踪这一对象 |
| 多起点候选 | 不同初始化后的端点 | 没有分支方程时只叫候选；archive迭代点不自动是分支 |
| 精确群作用 | 对所有允许材料有 $F(gm)=F(m)$ | 必须证明作用域、观测标签与box保持；几何置换不能任意套到非对称采集 |
| Monodromy | 非分歧复代数覆盖中沿闭参数路延拓引起的根置换 | Q6/Q20可提出研究，但未确定其群；路径通常经过复材料，不保证物理可行 |

对stationary continuation，若Hessian $H$非奇异，$m'(\lambda)=-H^{-1}\partial_\lambda\nabla f$。此时局部只有一条光滑分支；“每步把beta加倍并复制状态”不是数学branch splitting。Hessian奇异、box活动变化或Gram rank变化时需要分别分类。Pseudo-arclength可越过某些fold，但也可能跟到鞍点；它不是全局成像保证。

## 8. 局部corrector与transverse contraction的准确角色

设完整实白化forward为f，在常秩正则根附近使用

$$T_y(m)=m-J(m)^\dagger(f(m)-y).$$

在精确根 $m_*$，

$$DT_y(m_*)=I-J_*^\dagger J_*.$$

在可观测normal方向导数为0，在fiber切向为I。满列秩时是普通局部GN/Newton的零导数；它不证明存在多解fiber，也不把未知远端候选带入收缩邻域。

对于LM $T=m-(J^TJ+\lambda I)^{-1}J^Tr$，在根处normal特征值为 $\lambda/(\sigma_i^2+\lambda)$，null方向仍为1。这个公式只在光滑内点／固定活动面、准确Jacobian与相同残差下成立。非零残差必须考虑残差Hessian；有正则时固定点也改变。

因此分工应明确：有限proposal负责进入新的候选区域；物理corrector负责局部收敛；验收负责区分真正数据解与较低但仍不合格的极小点；先验／可辨识性决定多个合格物理解如何解释，而不是由corrector挑出GT。

## 9. 可直接执行的最小实验

### E0-C：先验证“机制不同”，不做576大跑

1. 跑本稿两单元解析例：正定beta metric从一个精确根出发应不移动；代数二次求根应得到两根；current chord数据全程一致、物理defect与定理一致。
2. 在原Hankel scalar-TM的4/6单元真实配置上，复用Z-ALG原多根对象，检查多源profile角障碍。不能把本稿抽象G当真实TM证据。
3. 对比真实CIE、标准CSI、material-profile CSI、受控normal释放，统一current空间、物理成本和初值。普通CSI已能同样穿过低障碍通道时，不将通道存在性当A17独立增量。
4. 所有初始化和proposal不读取另一真实端点；已知端点只用于离线检查定理。若测试算法用到了第二端点，则只算路径构造演示。

### E1-C：小块空间branch的决定性实验

16或36全未知像素、12源×若干接收；k=4/6支持字典预声明。比较：同成本block-TRF、多起点block-TRF、Q6/Q20选定数据proposal、随机支持、A17支持选择；全部使用同一精确核与全材料corrector。目标先覆盖新空间形状而不是旧错误形状幅度。记录真实原数据合格数、分离材料候选数、空间变化、best@K与在线top1、全RHS/小RHS/LU/wall；不以GT选支持。

### E2-C：normal释放与current增广分开测

固定同一FFT-TSOM起点，分别只增广物理tangent表示、只增加自由normal current、同时增广并设曲率保留门。匹配总current维数和成本，比较材料步偏差、profile Hessian特征、最终原数据达标、材料重建；若只看辅助残差下降，则不足以证明材料收益。

### 进入576的判据

必须至少出现一个相对强基线可重复的独立机制：例如同current预算下A17 normal门减少错误材料更新，或同物理预算下空间小块代数proposal增加真正合格材料／改善在线成像。单纯有限辅助目标更低、已知端点间存在通道、Q6复现旧小模型根，均不够支持576新大跑。

## 10. 给主研究者的建议

建议将正文主线定位为“材料相关current自由度的受控释放”，以原A17物理tangent选择和新增joint-current Schur曲率门作为可计算核心；finite chord／场角障碍说明为什么允许非物理current可能产生不同于纯材料metric的搜索机制。多分支提案保持为独立算法对照，优先复用已有Q6/Q20而不从零发明群胚语言。

最不建议的三件事：把A17-X已有row ratio再包装为新定理；继续扫uniform/full-image amplitude rays后只报告更低残差；用beta doubling产生的不同错误端点冒充等测量材料branch。

上述建议不是停止反馈方向。CIE仍是成熟强基线和可用state chart；但新数学应落在**current物理流形法向、释放后的材料曲率、有限分支通道与原数据恢复**，而不是形式上更复杂的 $\exp/\tanh$ 材料权重。
