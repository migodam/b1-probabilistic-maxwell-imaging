# 最终分支算法与第7/8节严谨性审阅

本文件供主稿整合，未修改 `report_body.md`。2026-10-31。

## 1. 审稿结论

主稿“正定度量不能移动孤立精确根”“finite-current弦线障碍”“profile实材料场角公式”“常秩corrector”四处数学一致。必要边界已经基本保留：固定观测、正定权重、box内点、共同实材料、非零row以及噪声相容并非严格商关系。

需保持的七条限定：

1. 精确根不动定理只针对不含额外材料正则的正定weighted residual。它不否认反馈改变非零残差驻点与吸引域。
2. finite-chord定理的两个同数据物理端点是前提。该定理证明通道存在和障碍值，不负责生成未知另一端点。
3. 通道中间的联合变量不满足本构，因此一般不是原物理forward；物理corrector与最终评分必须重算或使用严格等价有限核。
4. real-profile式必须使用 $(\Re\langle q_a,q_b\rangle)^2$，不可用模平方替代；否则漏掉实材料相位相容。
5. 发生clip时加上 $\|q\|^2(m_{\rm clip}-\widehat m)^2$；row为零时直接取 $\|j\|$，不使用比值。
6. Q6/Q20的全部根属于明确的固定支持、选定测量的有限系统。outside错误、所选通道不足以及数值退化都可以使它与完整成像目标不一致。
7. 相同响应不能作为材料archive唯一去重标准；不同材料正是多解任务需要保留的对象。反过来，大材料差但未通过原数据门的端点，只是候选，不是等响应物理解。

排版建议：主稿signed-gain段末有一个单独 `$` 接换行后 `\mathcal G_C>0`，建议并回一个完整数学定界；本文未修改。

## 2. 第二算法：真实空间小支持的精确有限响应分支提案

建议工作名称：**Finite-Support Branch Proposals（FSBP）**。这是描述性名称，不暗示已经取得论文级优势。它是完整材料反演的proposal层，与主算法“FFT-TSOM上的受控辅助电流释放”独立对照。

### 2.1 六元组

| 项目 | 可执行定义 |
|---|---|
| Optimized variables | 提案阶段为实际活动单元 $d\in\mathbb R^k$，$k=4$或6；等价current变量 $\eta\in\mathbb C^{k\times p}$。全材料corrector阶段恢复全部 $m\in\mathbb R^n$。 |
| Physical constraints | $m_0+U_Ad$ 满足原box与先验，$\eta=D(E_A+C\eta)$，$I-DC$可逆。共享材料用于全部源，不为每源单独拟合材料。 |
| Objective | 代数提案解选定原观测通道的兼容方程；全部原数据与共同正则用于候选评价；最后同预算最小化原物理 $\Phi(m)$，原数据门单独判断。 |
| Update | 支持内调用Q6/Q20条件商核，输出多个实合法候选；有限核计算全部源和接收的响应；从选出的材料候选运行共同的全材料TRF/LM/已冻结corrector。 |
| Subspace | $U_A$为真正坐标列。公共2×2和2×3空间块是初始字典，可加预先声明的分离cell组合。不能用“6个稠密形状系数”冒称rank-6材料更新。 |
| Globalization | 在空间不同的小块、有限代数根和全材料corrector之间搜索；不要求每个原始proposal比incumbent损失更低，但保存原incumbent。材料多样性和固定预算控制候选数量。 |

### 2.2 唯一使用的精确有限核

参考 $m_0$ 必须只由观测与公共初始化得到。令

$$L_0=I-X_0G,\quad j_0=L_0^{-1}X_0E,\quad Q_0=E+Gj_0,\quad F_0=Sj_0.$$

对真实坐标支持 $A$（$|A|=k$），

$$V_A=L_0^{-1}U_A,\quad C_A=U_A^*GV_A,\quad E_A=U_A^*Q_0,\quad O_A=SV_A.$$

对于 $D=\operatorname{diag}(d)$，

$$\eta_D=(I-DC_A)^{-1}DE_A,$$
$$\boxed{F_A(d)=F_0+O_A\eta_D.}$$

此式是完整离散非线性Maxwell在“outside暂时固定”条件下的精确消元；不是Born。$E_A$保留全部 $p$ 个原始源，$O_A$保留全部接收器。后面用于代数枚举的通道缩减**只发生在提案方程**，不改变这个all-data核。

### 2.3 k=4/6为何可执行：已有Z-ALG条件结构

取 $k=2s$，要求原观测已有源数 $p\ge s$；选一个已有接收通道 $\ell$ 和 $s$ 个已有照明索引 $I_s\subseteq\{1,\ldots,p\}$。$p\ge s$只保证有足够列可选，不保证后面的子核rank条件。若$p<s$，当前Q6/Q20配置不能运行，不能以虚构新照明补足。置

$$o=O_A[\ell,:],\quad q=E_A[:,I_s],\quad t=(Y-F_0)[\ell,I_s].$$

代数问题是

$$t=o(I-DC_A)^{-1}Dq.$$

定义接收对偶行

$$w=o(I-DC_A)^{-1}D,\quad v=o+wC_A,$$
$$wq=t,\qquad w_i=d_iv_i.$$

若 $\operatorname{rank}q=s$，用 $w=w_0+cN$ 消去测量，其中 $c\in\mathbb C^s$、$Nq=0$。实材料相容给

$$w_i(c)\overline{v_i}(z)-v_i(c)\overline{w_i}(z)=0,\quad i=1,\ldots,2s.$$

复代数求根中 $z$ 独立，最后必须检查 $z=\overline c$。将系统写成 $A(c)(1,z)^T=0$，取全部最大子式。

已成立的Z-ALG条件定理：若最高次系数块 $T$ 精确可逆，配对结构给有限唯一lift，商维数为 $\binom{2s}{s}$；若scheme约化且线性分离式有效，小乘法矩阵可恢复全部简单兼容候选。于是

| 活动实材料 | 原接收×原源 | 最高次块 | 商／特征矩阵 |
|---:|---:|---:|---:|
| 4 | 1×2 | 4×4 | 6×6 |
| 6 | 1×3 | 15×15 | 20×20 |

这里复用的是已经推导和运行的Q6/Q20，不宣称本轮重新发现该代数。更大k的商维数组合增长，不能外推为576未知量的可行全根求解器。

**必须逐项检查：** 源子核rank、最高次块条件、共轭lift、材料实性、恢复节点 $v_i$、原box、原材料先验以及 $I-DC_A$ 可逆性。$T$退化只表示当前消元图册失效，不能说物理无解；重根、rank不足和有限精度无法判定的情况应记录 `UNRESOLVED` 并付费fallback或跳过。

### 2.4 支持怎样选择、何时替换

第一版应使用容易复核的公共规则，避免一开始再引入未经辨别的A17 selector。

1. **公共空间字典。** 所有2×2连续块作为k=4候选、2×3及3×2连续块作为k=6候选；边界规则和stride在看测试观测前固定。若增加分离支持，规定“两个／三个公共相邻cell对”的组合规则与数量上限，不从GT形状提取。
2. **数据相关优先级。** 在同一anchor已付的全数据伴随梯度 $g=\nabla\Phi$ 上，为支持计算 $s(A)=\|g_A\|_2/\sqrt{k}$，或使用已付物理corrector的材料步分量。此分数只用于提案次序，不声称近似全局价值。
3. **保留非局部探索配额。** 每轮支持配额中，至少一份来自固定循环的“尚未尝试公共块”。若投影梯度接近零，不能停止所有支持尝试；原目标错误局部极小正是非局部提案需要处理的情况。
4. **A17独立消融。** 若已有同锚点、同残差下A17双侧步差 $\Delta p_C$，可用其在块内的幅度对公共字典重排；该计算全部收费，并与梯度排序、随机排序及循环排序等预算比较。它是算法候选，不是已证最优支持规则。
5. **支持替换。** 一块的所有预定通道／候选配额用尽后就切换字典下一块；不因GT分数或图像“看着像数字”而追加预算。anchor改变后核一般重建；只有代数上正确并计费的更新才能复用。不能沿用旧 $C_A,E_A,O_A$ 当作新anchor真实核。

这一规则保持全过程可修改所有材料，只是每次提案冻结outside。它没有证明小块足以突破全图形状盆地；这正是16/36像素最小实验要回答的问题。

### 2.5 通道选择不引入新信息

对 $E_A$ 做带列主元QR，选 $s$ 个线性独立原源列；若数值源rank不足，尝试预先允许的其它源子集或记录失败。接收通道可按 $\|O_A[\ell,:]\|_2$ 排序，以固定tie-break选首个；每块最多尝试预先固定的少数通道，避免无限挑选直到成功。

也可以用已有源／接收的已知线性组合，但最小实现先用实际通道索引，避免混合后噪声度量的额外混淆。这些都是现有数据的投影，**不是新增照明信息**。若新增真实频率或源，必须作为另一观测合同和费用单独比较。

### 2.6 三层验收不可混为一谈

**层1：代数与物理合法性。** 实材料、box、先验、无节点误除、有限核resolvent可逆。通过者是一个可完整正演的材料，不代表拟合全部观测。

**层2：全部数据评价。** 对每个合法d，用全部源 $E_A$ 求 $\eta_D$，得到全部接收的 $F_A(d)$，计算原白化残差与共同正则。未通过原数据门的点可以作为有标签的待校正proposal，不能称为等响应根。也不能因为其原始分数不是最小就全部丢弃：R4已显示raw线低谷排序与corrected排序可能不同。

**层3：全材料校正与最终资格。** 对进入配额的材料种子运行同一个完整n材料corrector，不再冻结outside。最终原物理数据门决定其是合格解还是仅优化候选；材料NRMSE/IoU只在封存后用GT离线评分，不用于选择或停止。

候选限额应预先规定。一个可复核选择器是：先按归一化材料距离聚类，每簇保留原物理分数最低者；预算中保留一份给最低原分数，一份给与incumbent材料相距较大的不同簇。空间多样性指标、配额和剩余预算规则须事前冻结。此选择可能漏掉有价值根，故一旦截断，不能再声称proposal集合完备。

### 2.7 伪代码

```text
INPUT:
    Fixed G, E, S, all measured Y, original noise metric
    Material box and common regularizer R
    Public initializer, support dictionary D4/D6
    Fixed budget B; support / channel / seed quotas
    Frozen Q6/Q20 solver; common full-material corrector

m0, paid_state = COMMON_ANCHOR(public_initializer, Y, B_anchor)
incumbent = m0
qualified_archive = []
proposal_archive = []

while paid_budget < B:
    Evaluate or reuse same-anchor L0 factorization, j0, Q0, F0
    Obtain paid original-data gradient / last physical material step
    Select next spatial support A using frozen queue and quota
    if no support remains: stop

    k = |A| in {4, 6}; s = k/2
    if number_of_existing_sources < s: record UNSUPPORTED; skip support
    Solve L0 V_A = U_A                            # k full RHS
    Form C_A = U_A^* G V_A, E_A = U_A^* Q0, O_A = S V_A
    Select s original source columns by rank-revealing QR
    Select frozen number of original receive channels

    for each selected receive channel ell:
        q = E_A[:, source_indices]
        o = O_A[ell, :]
        t = (Y-F0)[ell, source_indices]
        candidates, status = Q6_or_Q20(C_A, q, o, t)
        if chart fails:
            record UNRESOLVED; execute only prepaid fallback
        for each candidate d:
            check realness, conjugate lift, nodes, original box/prior
            check I-diag(d) C_A nonsingular
            if any check fails: retain rejection reason and paid cost
            else:
                eta = solve(I-diag(d) C_A, diag(d) E_A)  # all p sources
                F_candidate = F0 + O_A eta               # all r receivers
                phi = ORIGINAL_SCORE(F_candidate, Y, R(m0+U_A d))
                store material seed, all-data residual, phi, provenance

    seeds = FROZEN_SCORE_AND_DIVERSITY_SELECTOR(valid_candidates, quota)
    for seed in seeds:
        m_end = FULL_MATERIAL_CORRECTOR(seed, all Y, per-seed paid quota)
        verify full original Maxwell / exact equivalent response
        if original data gate passes:
            add to qualified_archive, retaining distinct materials
        else:
            retain as nonqualified optimization candidate
        update incumbent only by common original score

    choose next anchor by frozen corrected-score/archive policy
    refresh kernel after anchor/support change; never reuse stale state

OUTPUT:
    incumbent, separately labelled qualified material set,
    proposal/rejection archive, all paid costs, unresolved / truncated items
```

没有根、没有合格根、当前图册失效、预算没运行完是四个不同状态。程序必须分别输出，不把它们合成“无解”。

### 2.8 能证明什么与不能证明什么

**条件包含性。** 若某个完整数据物理解与当前anchor只在A内不同；选定测量的Q6/Q20先决条件成立；相应根简单、非节点非极点；枚举器在精确代数下返回全部根且未截断，那么该物理解必在proposal集合中。因为它满足选定测量和共享实材料兼容式，条件商枚举完整覆盖这些根，正确合法性检查不会删除它。

这不是对未知576材料GT的覆盖保证：GT通常并非只在4/6个像素上不同；真噪声下选定方程根也并非全数据精确解；数值容差和候选截断会破坏精确包含性。所有结论应以实际预算下的合格根／成像收益检验。

### 2.9 真实成本与强对照

已有同锚点完整LU时，每块新增 $k$ 个完整RHS，并支付 $GV_A$、端口装配和QR。每个有限候选需要k阶LU、全部p个小RHS及输出乘法；完整数据评分在精确核下不需要再做n阶LU，但最终独立核验或全材料corrector有其真实费用。改变anchor一般触发新的n阶LU和基点场求解。

必须比较：

- 同一支持、同一精确核的成熟block-TRF，多公开起点；
- 相同总预算的普通全材料多起点；
- Q6/Q20与专门双线性／行列式求解器；
- 同k、同预算的循环／随机／梯度／A17支持选择。

不能让Q6/Q20免费使用有限核，却让对照每次重做完整LU；也不能用更多种子带来的效果声称新几何优势。root枚举、系数装配、失败、fallback与最终corrector全部入账。

## 3. Continuation与branch机制表

| 操作 | 改变的对象 | 是否增加真实信息 | 合法的理论表述 | 主要风险 |
|---|---|---|---|---|
| 固定beta多起点 | 初始化及辅助目标参数 | 否 | 预算内不同吸引域的搜索 | 机会数增加冒充几何增量 |
| beta／反馈图册continuation | 非零残差目标的度量与材料读出 | 否 | regular stationary branch在H可逆时由IFT延拓 | 孤立精确根不动；pole/Gram病态/box事件不同 |
| current-subspace continuation | 联合变量中可用current自由度 | 否 | 新模式降低局部剖面目标；材料曲率可能downdate | current替材料拟合数据而不恢复物理 |
| 物理current ROM增广 | 对同一真实Jacobian／finite响应的近似精度 | 否 | A17有符号Schur exchange、误差与材料步差 | 与自由current释放混淆 |
| 材料8→12→24 | 材料参数空间 | 否 | 低维先验到全材料自由度的逐步释放 | coarse解约束错误shape；不代表物理网格变小 |
| PDE penalty rho continuation | 允许的本构违约程度 | 否 | finite-chord障碍与rho线性缩放；需要恢复物理 | rho过小数据拟合但材料任意，过大恢复原难度 |
| 选定源／接收→全部原数据 | 已有观测在目标中的约束集合 | 否 | observation-subset homotopy或proposal构造 | 子集根未必拟合全部数据，不能称新信息 |
| 新频率／新照明采集 | 真实测量合同 | 是，需按噪声关系分析 | 新数据可增加可辨识性 | 新观测收益混入算法收益；额外采集未收费 |
| Stationary pseudo-arclength | $(m,\lambda)$驻点曲线的参数化 | 否 | 某些fold附近继续跟踪 | 可能进入鞍点；不保证到真实材料 |
| 精确fiber局部continuation | 固定观测常秩等值流形 | 否 | 法向corrector＋切向移动 | rank change、prior边界、孤立根不支持连续fiber |
| Q6/Q20多根proposal | 固定支持兼容代数的不同根 | 否 | 明确条件下有限根枚举 | 节点／极点／实性过滤及组合复杂度 |
| Monodromy | 明确复代数覆盖上闭参数路引起的根置换 | 否 | 避开判别集的代数延拓 | 复路径不保持实box；已有Q6/Q20未求全群 |

## 4. 两单元例的精确有理审计

原JSON未改，图仍只使用该已有模型。正文应将观测定义为 $y:=F(2,7)$，小数仅作近似显示。

矩阵和端口均可视为精确有理复数：

$$G=\begin{pmatrix}-3/50+i/50&13/200+7i/1000\\13/200+7i/1000&-3/20+i/50\end{pmatrix},$$
$$E=(1,7/10)^T,\qquad S=(1,13/10).$$

直接有理运算给

$$\boxed{y=\frac{1826362941900+157970241200i}{315300667481}.}$$

消去第二材料后的实二次方程，除非零常数因子，可精确写成

$$\boxed{(m_1-2)(764347498747m_1-3470613075400)=0.}$$

因此判别式与

$$\boxed{(1941918077906)^2>0}$$

成正比例；两根简单。完整第二材料回代得到

$$\boxed{m^{(b)}=\left(\frac{3470613075400}{764347498747},\frac{394249570000}{195029801279}\right).}$$

本轮使用Python标准库Fraction精确核验：二次式在 $m_1=2$ 精确为零、第二根回代虚部精确为零、两材料分量严格在 $(.02,20)$、恢复分母及物理det精确非零。无须依赖浮点Jacobian奇异值来证明两根是局部孤立。分母不为零意味着二次消元在两根附近等价；简单根给局部孤立性。

### 全正材料域无物理极点

Hermitian实部

$$H=\frac{G+G^*}{2}=\begin{pmatrix}-.06&.065\\.065&-.15\end{pmatrix}$$

首主子式为负，行列式 $.004775=191/40000>0$，故 $H\prec0$。任意严格正对角X上，

$$L=I-XG=X^{1/2}\widetilde L X^{-1/2},\quad\widetilde L=I-X^{1/2}GX^{1/2},$$
$$\frac{\widetilde L+\widetilde L^*}{2}=I-X^{1/2}HX^{1/2}\succ I.$$

因此 $\widetilde L$ 严格accretive且可逆，L也可逆；整个正box没有物理极点。相似变换非unitary，不能直接推断 $\sigma_{\min}(L)\ge1$；该不等式只直接适用于 $\widetilde L$。若要全box粗界，可由相似变换得到 $\sigma_{\min}(L)\ge\sqrt{m_{\min}/m_{\max}}$。图中两条具体路径测到 $\min\sigma(L)\ge1.046$，属于数值结果。

### 图的准确证据范围

| 量 | 线性材料 | Profile实材料 |
|---|---:|---:|
| 原current的state defect最大值 | .4208596 | .0151926 |
| 材料重算物理后的相对data residual最大值 | .0859374 | .00367337 |

current chord本身相对data residual不超过 $3.06\times10^{-16}$，与上表右列完全是不同对象。profile全程在端点box中。该模型为复对称、虚部正定的代数VIE，**不宣称是实际Hankel圆盘采集网格、没有新增真实Maxwell成像benchmark，也没有证明未知端点搜索成功。**

## 5. CC-CSI补充核验，供总算法表使用

已实际取得并阅读作者全文 https://arxiv.org/pdf/1906.10864 的Eq4–22及Algorithm1，web ref `turn36view0`。变量是共享contrast $\chi$ 与各源contrast current $j_p$；$Q=E^{inc}+G_Dj$。

在本报告符号 $R=Sj-Y,e=j-XQ$ 下，原文三残差为 $\rho=-R,\gamma=-e,\xi=Y-SXQ=-(R-Se)$，目标为

$$\eta_S\sum_p\|R_p\|^2+\eta_D\sum_p\|e_p\|^2+\eta_S\sum_p\|R_p-Se_p\|^2,$$

其中 $\eta_S^{-1}=\sum_p\|Y_p\|^2$，$\eta_D^{-1}=\sum_p\|XE_p^{inc}\|^2$。current与材料使用交替Polak–Ribiere共轭梯度，材料线搜索还处理归一化，用Brent法；不是逐行闭式LS的简单重复。

第三项用背景current-to-receiver算子S把state与data缺陷交叉起来；它不同于完整材料resolvent传播 $SL(m)^{-1}e$ 的精确物理输出证书。正文数学/算法已核验，实验图逐幅未核验，不从论文图猜新数值。DOI 10.1109/TAP.2017.2673758。
