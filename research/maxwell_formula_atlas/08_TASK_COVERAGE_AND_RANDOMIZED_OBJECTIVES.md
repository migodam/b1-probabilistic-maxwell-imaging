# 08 · 径向—角向覆盖、随机化目标与任务方向

本章把 `POST_E5_THEORY_DECISION.md` 中F-POST-001—003正式纳入持续维护的公式集。保持原编号，不把标准推导重新命名为新定理。F-POST-001是材料空间中的精确几何关系；引用E5的解释属于 `POSTHOC_EXPLANATORY_ONLY`。F-POST-002与003是已知线性代数/随机优化构件的推导，其作为新候选生成机制的状态均为 **`NOT_IMPLEMENTED / NOT_RUN`**。

E5实际实施与评分仍以 [最终运行协议](../config/possibility_volume_protocol.json) 和 [冻结统计协议](../config/possibility_statistics_frozen.json) 为准。添加公式不改变候选、预算、主统计门或已有结果。本章新增物理RHS、训练、数值测试均为0。

## F-POST-001 · 材料邻域覆盖的精确径向—角向分解

### 空间、假设与推导

原材料 $m\in\mathbb C^{576}$，采用全网格复材料RMS距离。令

$$z(m)=\frac1{24}\begin{bmatrix}\operatorname{Re}m\\\operatorname{Im}m\end{bmatrix}\in\mathbb R^{1152},\qquad d_M(m,m')=\|z(m)-z(m')\|_2.$$

固定anchor $m_a$、真实材料 $m_*$ 和一个候选 $m_j$。定义 $g=z(m_*)-z(m_a)$、$u=z(m_j)-z(m_a)$，以及 $R=\|g\|$、$\rho=\|u\|$。若 $R\rho>0$，令 $\cos\theta=g^Tu/(R\rho)$。展开平方即得

$$\boxed{d_M(m_*,m_j)^2=\|g-u\|^2=R^2+\rho^2-2R\rho\cos\theta.}$$

这在材料嵌入空间中精确成立，不要求decoder线性、弱散射、固定Jacobian或GT方向各向同性。若 $R=0$，误差就是 $\rho$；若 $\rho=0$，候选就是anchor，误差就是 $R$，两种情形均无需定义角度。

当 $R\rho>0$ 时，覆盖 $d_M(m_*,m_j)\le\varepsilon$ 等价于

$$\cos\theta\ge c(R,\rho,\varepsilon):=\frac{R^2+\rho^2-\varepsilon^2}{2R\rho}.$$

由于 $\cos\theta\le1$，必要条件为 $|R-\rho|\le\varepsilon$。当 $c>1$ 时任意方向都不可能命中；当 $c\le-1$ 时任意方向都可命中；$c=1$时只有严格同向的退化角帽。**径向可达只是必要条件，不能代替角向覆盖。**

与anchor比较，候选改善要求

$$d_M(m_*,m_j)^2<R^2\quad\Longleftrightarrow\quad \rho<2R\cos\theta\qquad(\rho>0).$$

因此沿一个固定正方向，误差关于长度为二次函数 $R^2+\rho^2-2R\rho\cos\theta$，其在 $\rho\ge0$ 上的最优长度为

$$\rho_* =\max(0,R\cos\theta).$$

若允许沿同一直线取signed长度，则最优signed系数是 $R\cos\theta$；E5的weak-line系数允许signed变化，但其材料方向也随decoder/corrector变化，不能把这个一维oracle直接当作当前算法。$R$与$\theta$都依赖未知GT，以上最优长度是解释工具，不是部署规则。

### 线性子空间与有限K的条件解释

仅在另行指定的固定线性材料空间 $\mathcal S$ 内，分解 $g=g_\parallel+g_\perp$，其中 $g_\parallel=P_\mathcal Sg$。设候选增量为 $u=\rho v$、$v\in\mathcal S$、$\|v\|=1$。由正交性，

$$\|g-\rho v\|^2=\|g_\perp\|^2+\|g_\parallel\|^2+\rho^2-2\rho\|g_\parallel\|\cos\theta_\parallel.$$

若 $\|g_\perp\|>\varepsilon$，无限增加该空间内的点也不能覆盖GT。若 $\rho\|g_\parallel\|>0$，覆盖对应一个角帽，其阈值为

$$c=\frac{\|g_\perp\|^2+\|g_\parallel\|^2+\rho^2-\varepsilon^2}{2\rho\|g_\parallel\|}.$$

若进一步假设 $d=\dim\mathcal S\ge2$，提案方向独立且均匀分布于 $S^{d-1}$，那么对固定GT方向，$t=\cos\theta_\parallel$ 的密度为

$$f_d(t)=\frac{\Gamma(d/2)}{\sqrt\pi\Gamma((d-1)/2)}(1-t^2)^{(d-3)/2},\qquad -1<t<1.$$

这是单位球面按第一坐标切片得到的标准密度；其归一化由Beta积分给出。于是对 $-1<c<1$，单次命中概率与独立K次至少命中一次的概率分别为

$$p_d(c)=\int_c^1 f_d(t)dt,\qquad P_K=1-[1-p_d(c)]^K\le Kp_d(c).$$

这个概率解释假设的是**提案方向均匀**，不需要GT方向随机或均匀。不同半径对应不同 $p_k$ 时，独立性下为 $1-\prod_k(1-p_k)$；没有独立性时最多保留适用的union bound，不能继续使用乘积。antithetic配对、动态frame、非线性decoder、corrector以及资格筛选都可能破坏均匀性或独立性。$d=1$只有正负两个方向，需按离散分布单独处理。

**不把此条件模型冒充E5概率。** 不能直接代入d=8或66推断本轮命中率，也不能以“高维”作为全部失败原因。真实候选不保证留在初始线性空间；当前corrector也会改变方向。

### E5的posthoc含义与边界

既有结果显示：B的新增主ε=.075命中来自.05半径，.10/.15没有新增主命中；保存同一批路径全部visited噪声合格点后仍为6/32。较大半径不是自动收益，与上述精确几何关系一致；没有据此声称已经逐候选计算了GT角度。

实际候选仿射空间的GT投影是另一个事后容量诊断。若 $U$的列为实际候选相对anchor的增量，则

$$\|(I-UU^\dagger)g\|\le\min_j\|g-Ue_j\|.$$

这是因为每个实际候选均位于该仿射空间。B有19/32个parent的这一oracle剩余RMS仍高于.075，另7个parent的仿射oracle可进入.075但有限候选未命中；前者揭示实际方向张成的遗漏，后者揭示点/系数放置的遗漏。数值秩依据现有诊断阈值；仿射投影可能使用负系数、违反decoder或数据门，并且使用GT，不能作为新成像结果。事实来源：[posthoc汇总](../results/possibility_volume/posthoc_geometry/SUMMARY.json)、[parent诊断](../results/possibility_volume/posthoc_geometry/PARENT_DIAGNOSTICS.csv)、[封存后理论裁决](../POST_E5_THEORY_DECISION.md)。

**资格：** 精确材料几何身份 + 明示假设的球面概率解释；E5连接为 `POSTHOC_EXPLANATORY_ONLY`。平方展开、正交投影与球面角帽公式均为标准数学，不主张原创定理。

## F-POST-002 · 随机化data/prior目标产生完整局部方向

### 原目标与局部模型

采用Gaussianized实latent $x\in\mathbb R^{66}$，原prior为标准Gaussian。白化数据residual $r(x)\in\mathbb R^{512}$，原目标为

$$\Psi(x)=\tfrac12\|r(x)\|^2+\tfrac12\|x\|^2.$$

在一个固定anchor $x_a$，记 $r_a=r(x_a)$、$J=\partial r/\partial x|_{x_a}\in\mathbb R^{512\times66}$。固定这些量并作局部线性化 $r(x_a+p)\approx r_a+Jp$。令

$$H=J^TJ+I_{66}\succ0,\qquad g=J^Tr_a+x_a.$$

取彼此独立的 $\eta\sim\mathcal N(0,I_{512})$ 与 $\xi\sim\mathcal N(0,I_{66})$，分别随机化白化data目标与prior中心。局部随机目标为

$$q_{\eta,\xi}(p)=\tfrac12\|r_a+Jp-\eta\|^2+\tfrac12\|x_a+p-\xi\|^2.$$

对p求导并置零，得到唯一最优位移

$$\boxed{p=H^{-1}(J^T\eta+\xi-g).}$$

均值为 $-H^{-1}g$；只有anchor恰为该局部目标驻点、$g=0$时才是零均值。这个项不能因anchor被称作MAP而未经检查省略。

独立性给出

$$\operatorname{Cov}(J^T\eta+\xi)=J^TJ+I=H,$$

因而

$$\boxed{\operatorname{Cov}(p)=H^{-1}HH^{-1}=H^{-1}.}$$

它在全部66个latent方向上正定，并按数据曲率缩放，而不预先把weak8之外的方向置零。只有随机化起点、却把同一严格凸目标优化到底，仍会返回同一点；这里改变的是优化目标的随机线性项。

### 一般线性Gaussian问题中的精确分布资格

若完整forward本身线性，$y=Ax+\epsilon_{obs}$、$\epsilon_{obs}\sim\mathcal N(0,\Gamma)$，prior为 $\mathcal N(m,C_0)$，且 $\Gamma,C_0\succ0$，独立抽 $\eta_y\sim\mathcal N(0,\Gamma)$、$\xi_x\sim\mathcal N(0,C_0)$，解

$$x_* =\arg\min_x\left[\tfrac12\|Ax-y-\eta_y\|_{\Gamma^{-1}}^2+\tfrac12\|x-m-\xi_x\|_{C_0^{-1}}^2\right].$$

令 $H_0=A^T\Gamma^{-1}A+C_0^{-1}$。由正规方程，

$$x_*=H_0^{-1}\{A^T\Gamma^{-1}(y+\eta_y)+C_0^{-1}(m+\xi_x)\}.$$

它是Gaussian，均值为 $H_0^{-1}(A^T\Gamma^{-1}y+C_0^{-1}m)$，协方差为 $H_0^{-1}$，与该线性Gaussian posterior一致。这项精确资格来自完整线性模型和精确优化，不自动延伸到非线性Maxwell、截断域或受限corrector。

### 实现位置、成本与失败条件

固定anchor局部采样可以复用H的一次factorization；但若沿每个随机目标继续非线性优化，r、J、H都随状态变化，需要支付真实物理与导数费用。固定host/ROI算子仍可共享，各候选的变化材料系统不能免费冻结。使用近似解、局部极小点、硬几何约束、拒绝或去重，均会改变输出分布。

作为可能性云方法时，候选最终必须对**原观测y**检查数据资格；不能用其随机化后的观测目标放行。它可能扩展当前局部探索，但不保证数据合格、高GT覆盖或正确posterior。随机化目标、其有限非线性优化与过滤分布的关系，需要独立分析与实验。

**当前状态：`NOT_IMPLEMENTED / NOT_RUN`。** E5未运行此机制，不把局部Gaussian公式当作已发生的Full66采样。提出它是为了给下一次机制裁决一个成熟强基线，而不是替当前weak8失败补造结果，也不授权追加实验。

### 先例与原创性边界

优化产生随机候选的直接先例为 [Bardsley等，RTO，2014](https://doi.org/10.1137/140964023) 与 [Wang、Bui-Thanh、Ghattas，rMAP，2018](https://doi.org/10.1137/16M1060625)。rMAP作者明确说明特定随机data/prior中心下与RML的关系；RTO作者明确区分线性直接posterior样本与非线性需密度/接受校正的情形。本轮新核验范围是出版社/作者摘要，不冒称逐式审核其全部非线性算法。上述局部和线性Gaussian公式为独立标准推导，不把它等同于完整非线性RTO实现，也不主张原创。证据层级与检索429记录见 `PAPER_CORE_AND_ORIGINALITY.md`。

## F-POST-003 · 在曲率预算下选择任务敏感方向

### 假设与推导

令任务 $T(m)\in\mathbb R^s$，$W\succeq0$为预先确定的任务单位/重要性metric，$m=\mathcal G(x)$。假设复材料decoder作为实映射以及任务在anchor附近可微。定义

$$B=\left.\frac{\partial[W^{1/2}T(\mathcal G(x))]}{\partial x}\right|_{x_a}\in\mathbb R^{s\times66},\qquad H=J^TJ+\lambda I\succ0,\quad\lambda>0.$$

当λ=1且J按原噪声白化时，H是当前Gaussian prior对应的GN曲率；其他λ是另行声明的正则尺度，不能静默当作原prior。给定 $b\ge0$，考虑局部方向问题

$$\max_{p^THp\le b^2}\ \|Bp\|^2.$$

令 $u=H^{1/2}p$。约束变成 $\|u\|\le b$，目标变成

$$u^T\underbrace{H^{-1/2}B^TBH^{-1/2}}_{C_T\succeq0}u.$$

由对称矩阵Rayleigh商，最优值为

$$\boxed{b^2\lambda_{max}(C_T).}$$

若最大特征向量 $v_T$ 取单位范数，则一个最优位移为 $p_*=bH^{-1/2}v_T$。将 $C_Tv_T=\mu v_T$ 换回原坐标，得到等价广义特征问题

$$\boxed{B^TBp=\mu Hp.}$$

当B为零时最优值为零，不能因此声称任务在有限非线性区域恒定；当最大特征值重复时，最大特征子空间内方向不唯一，不应给单向量附加无根据的物理身份。

### 固定弱空间中的受限版本

若仅允许 $p=Qa$，其中 $Q\in\mathbb R^{66\times k}$满列秩，则 $Q^THQ\succ0$，问题变为

$$\max_{a^TQ^THQa\le b^2}a^TQ^TB^TBQa,$$

对应 $Q^TB^TBQ\,a=\mu Q^THQ\,a$。Q不必Euclidean正交，但必须保留真实 $Q^THQ$，不能因为用了材料归一化基便把该metric替换为单位阵。限制在weak8内只比较其中任务方向，不能解决其空间外的覆盖遗漏。

### “预算”与实际物理资格不是一回事

$p^THp$是明确的局部二次曲率预算。若anchor梯度 $g\ne0$，原二次目标变化为 $g^Tp+\tfrac12p^THp$，而不是仅有二次项；即使g=0，有限步的非线性余项仍需处理。因此该约束不能自动解释为真实likelihood变化上限，也不等于材料半径约束或绝对噪声资格。

任务重要性与数据不可辨识性是不同指标：某方向数据很弱但任务几乎不变，另一个方向可能在相同曲率代价下产生明显任务变化。该广义谱给出这种区别的局部量化。W必须由任务尺度/容差事先确定，不得看到GT后调整。

H、B通常随状态变化；若任务仅是材料decoder的函数，构造B主要需要任务/decoder导数，但构造当前J仍要付物理成本。hard threshold在阈值处不可微，不可直接套用；soft occupancy可以求导，但局部导数不保证硬事件正确。任务数小使 $\operatorname{rank}(B^TB)\le s$，这不意味着剩余材料方向可删除，它们仍可能影响GT邻域、其他任务或非线性行为。

最大task spread也不保证围绕GT展开；它可能仅把错误候选拉得更远。因此该方向准则必须与真实物理资格、集合范围、有限候选预算及独立GT/task覆盖一起裁决。F-VOL-015说明为什么task dispersion与near-GT权重还会影响最终加权动作。

### 资格与已知背景

**当前状态：`NOT_IMPLEMENTED / NOT_RUN`。** E5计算了任务表现，没有实施这里的任务导向proposal或槽位分配。Rayleigh商/广义特征值是标准数学；目标导向降维也不是本轮新概念。已核验的 [Cui等2014 LIS](https://arxiv.org/html/1403.4680v2) 与 [Zahm等 certified reduction](https://arxiv.org/html/1807.03712v4)构成prior/likelihood降维背景，但不能把它们说成已经逐式提出本章具体任务目标，也不能因只核验这几篇就声称本构造没有先例。

若未来把任务重复性、原物理误差界与这个局部方向准则组合成独立算法，需要补做最接近的目标导向降维/实验设计先例核验，并与成熟随机目标方法和当前精确ROI强基线比较。这里仅保留推导与成立条件，不进行新颖性认证，不启动新的物理实验。

## 使用与状态规则

F-POST-001说明现有材料点为何可能“大而不覆盖”；F-POST-002说明完整局部方向可以通过已知随机目标机制被激发；F-POST-003说明任务方向与单纯数据弱方向不同。三个公式分别属于几何解释、未运行的随机目标机制和未运行的任务导向机制，不应合并成一个已经验证的新算法。

实际E5结果只能由封存后的证据支持；局部Gaussian、球面独立均匀提案、任务可微性与曲率预算等条件必须随公式一起引用。新增本章不改变任何历史门或已冻结配置。
