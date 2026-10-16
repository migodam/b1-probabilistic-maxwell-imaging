# Research F — 非线性 Maxwell 材料解集、有限反馈与可学习探索

## 从“剩余误差在弱方向”到“哪些材料世界可达、可采样、值得学习”

**研究日期：2026 年 10 月 16 日。**  
**交付性质：证据审查、独立数学推导、算法设计与下一轮实验规范。**  
本报告没有启动新的 Maxwell 网格仿真、GPU 训练、Mie 计算或历史实验重跑。随附计算仅为六维随机矩阵与两单元代数模型的轻量检查，不是新的电磁成像证据。

---

## 执行摘要：科学裁决

**应当继续研究材料探索，但不应把“训练一个沿 N3 的 nullspace 行走的网络”作为目前的主路线。** 更稳健的研究对象是：

\[
\boxed{\text{合法材料先验上的有限噪声后验，以及其中受物理控制的有限材料传输。}}
\]

精确等数据 fiber 是这个对象在特定常秩、小噪声条件下的几何极限，不是所有实例都具有的先验事实。Explorer 学到的是高效移动的计算机制；它不会凭空恢复观测没有提供的信息，也不会因为能生成很多图就成为正确后验采样器。

这次研究有六项关键结论。

**第一，弱差量不等于局部自由度。** 若正向映射二阶导数有界，GT 与重建足够近且数据相同，则 GT 差量的强方向分量可以被二阶项控制，而不必严格为零。反过来，大部分误差能量落在冻结 Jacobian 的弱空间，既不证明 prior 内存在非平凡 fiber，也不证明 GT 在同一连通分支。最新诊断中的 N3 98.09%–99.06% 弱能量只是一个特定阈值、完整复材料空间中的离线几何观察。[R1]

**第二，当前两个合法 P0 smooth chart 更符合“局部可辨识、有限噪声有宽度”，而非“已有精确局部解流形”。** 去除相同材料的权重缩放 gauge 后，两个七维材料切空间的响应矩阵数值满列秩。若该性质在邻域中成立，精确 fiber 局部是孤立点。它不排除其他地方存在离散竞争解释，也不排除有限噪声后验有值得采样的范围。[R1,R2]

**第三，完整反馈确实能精确描述有限材料变化，但不能创造补偿解。** 最新诊断把 N3→GT 的一阶响应误差约 8.26% 降至三阶约 0.0203%，有理形式达到数值精度；然而 GT 差量仍不是等数据方向。有限 resolvent 的主要价值是准确验证、隐式传输、曲率计算与受控近似，不是证明“高阶一定能绕到 GT”。[R1]

**第四，本次最值得优先验证的新判断是：先前所称 decoder 曲率，尚未分离坐标加速度与真实先验第二基本形式。** 对组合映射 \(F\circ g\)，\(J D^2g\) 的大小依赖参数化。一个完全平坦的材料平面，也能通过弯曲坐标产生任意大的这一项。真正不能被坐标重选消除的是材料先验的外在曲率，以及其与 Maxwell Hessian 组合后落入数据余核的部分。由此可以设计一个不支付新正向求解的缓存审计，先判断此前巨大曲率比中究竟有多少是坐标问题。

**第五，N3 是应当保留的单图决策基线，不是 Explorer 必须修补的主体。** 实际 N3 使用十二通道确定性特征、五层残差 CNN、统一材料 box，以及归一化复材料 MSE 加 0.05 SSIM 损失；SSIM 对实虚部使用不同固定动态范围。它不是理想条件均值的已验证实现，更不是后验。尤其，固定复相位先验的精确条件均值应保留该相位；实际 N3 的相位偏离不能仅由“对多个合法世界取平均”解释。[R3–R5]

**第六，A17 的潜在增量应该转向“有限材料传输的双侧误差控制”，而不是继续追逐冻结 GN gap。** 这里确有一个清楚的数学接口：有限电流残差、接收端伴随、材料 realizability 和可信步长。是否优于普通完整求解或 Krylov，尚无证据。已有 A17 证明的是指定 current-space 改变如何影响冻结材料更新；它不自动控制弱子空间、fiber 曲率、后验密度或采样混合。[R6]

最终推荐：

> **先完成 intrinsic-geometry audit 与合法的独立 physics anchor，再实现材料度量下的有限噪声探索器。只有它能稳定走出有意义的合法材料距离，才学习有限补偿映射。概率任务优先保留明确密度的 MH/pCN；只有发现独立高支持区域，才引入全局生成跳跃。**

本报告给出了六组有明确假设的命题、两个 Maxwell-form 反例/恒等式、一种不依赖 N3 的主算法，以及本质不同的 current-space 备选。这里的命题不等于六个原创论文定理：常秩、coarea、MH、谱扰动与 Grönwall 是已有数学；本次贡献候选是它们在实际 Maxwell/prior 结构中的组合、诊断方式和可计算接口。

---

## 0. 证据边界与阅读范围

### 0.1 实际读到的内容

本次直接读取了五个必读仓库的有关入口和关键报告；重点核对了固定 commit 的 finite-feedback 诊断、METRICS.json 中的核心记录、LOCAL_PATH_AND_CURVATURE、REPRODUCE、THEORY_QUESTIONS_FOR_PRO，LXZ 的完整中文报告，A26 decoder 与训练损失/归一化源码、配置和研究报告，以及 A17 正式稿的 CORE_THEORY.tex。[R1–R7]

A30 的主要容量/补偿负结果也已核对。其 Sparse-64 的 35.69% 是离线 oracle 表达容量，不是在线反演收益；不能拿来证明 posterior mode 已存在。[R7]

跨会话检索返回了 Dyson、材料优先概率反演及相关研究的部分历史上下文，但**不是用户列出的每个 Session 的全文**。独立 L、Z、X 原始完整 ZIP 未在本次逐文件展开；LXZ 对它们的归纳属于仓库二手整合证据，而不是本报告假装完成了对 41 份原始材料的全文复核。

### 0.2 没有完成的验证

公开诊断报告提供了原始 NPZ 的发布位置，但本次未下载并逐端点重算那 45 份二进制数组。因此，下文把这些数值称为“仓库报告/机器可读记录的结果”，**不称本次独立 full-Maxwell 复现实验**。原始端点重放留给下一轮只读 Codex 审计。

文献方面，已取得并阅读 Montanaro 等人的 e-GLASS、Deep Null Space Learning、Manifold Lifting 和 DINO 的可解析全文有关部分。Boink 等人的非线性数据一致网络取得了原论文出版信息、摘要与方法描述，但未成功完整取得其证明正文。2026 年 8 月 inverse-Born neural operator 新文目前只有摘要级核对，不作全文新颖性排除。部分 PDF 截图接口失败；本报告不依赖未经视觉核对的图表读数。

### 0.3 证据不能混用

A17 的三维 vector-Maxwell/DDA、A26/A30/LXZ 的二维 scalar-TM/VIE、Joint P2 的径向 vector-Mie，以及 C1+MH 的四维 proxy posterior 是不同目标与数值层级。它们不能互相替代完整成本或概率保证。

同样，P0 在不同早期原型中曾有不同语境。本报告的“原 P0”特指 finite-feedback/LXZ 审计所使用的原材料生成先验；新的 P1 会明确另立统计目标。

---

## 1. 七种集合，以及真正要研究的对象

令材料采用已实化、质量加权的坐标 \(x\in\mathbb R^n\)，数据采用固定噪声白化坐标。已知背景、源、接收器、频率及固定 gain 均视作已给条件。正向映射记为 \(F(x)\)。

| 对象 | 定义示意 | 不能由它推出什么 |
|---|---|---|
| 精确等数据 fiber | \(\mathcal F_{y_0}=\{x:F(x)=y_0\}\) | 不自动满足原材料先验 |
| 有限噪声可行 tube | \(\mathcal T_\tau=\{x:\|F(x)-y\|\le\tau\}\) | 不是后验可信集，更不是等概率集 |
| 近最优 likelihood 集 | \(\Phi(x)\le\inf\Phi+\Delta\) | 未加入 prior，也不保证连通 |
| prior-supported 子集 | 与 \(\operatorname{supp}P\) 相交 | support 内也可能是极低密度 |
| posterior 高密度区 | 由 \(\pi(x\mid y)\) 和参考测度定义 | 高密度不等于大概率质量 |
| 数值优化 basin | 某算法收敛到同一驻点的初值集 | 与物理 fiber 的连通分支不是同一对象 |
| 局部材料连通分支 | 指定合法集合内路径连通的部分 | 不含其他不可连续抵达的解释 |

这里 \(\Phi(x)=\frac12\|F(x)-y\|^2\)。高维“密度大”和“质量大”应严格区分：某个极窄高峰可能比一个较宽的较低峰拥有更少质量。

**重要修正：** 如果 \(F(x_0)\neq F(x_*)\)，GT 就不属于穿过 \(x_0\) 的精确 fiber，无论两幅图多像。实际有噪声时，正确对象通常是有限噪声 posterior，而不是强迫所有样本满足 \(F(x)=F(x_0)\)。

最新诊断的一个案例中，GT 自己的白化残差为 24.090925，超过预先冻结的 95% noise 阈值 23.785449。这可以是正常噪声尾部，并不自动意味着模型错了。若把该球作为硬支持，GT 会被定义性地排除；Gaussian posterior 则仍给它非零密度。不能事后改 sigma 或扩大阈值来把它“救回”。[R1]

### 1.1 材料先验是集合，也是测度

设合法材料由 \(x=g(z)\) 表示。需要分别记录：

\[
\mathcal P=\operatorname{Image}g,\qquad z\sim p_Z.
\]

\(\mathcal P\) 决定哪些材料可出现；\(p_Z\) 决定它们有多受支持。改变 decoder、soft rasterization、潜变量维数、材料相位或物理 mixture 解释，可能改变集合；去 gauge 后假设新坐标均匀，即使集合不变，也会改变测度。

因此“同一个 prior”至少要同时保证支持、生成映射、参数分布、已知条件和离散化规则一致。

---

## 2. Q1：剩余误差何时可以解释为弱自由度？

### 2.1 有限差量与切向量是不同对象

令 \(h=x_*-x_0\)，\(J_0=DF(x_0)\)。精确恒等式是

\[
F(x_*)-F(x_0)
=\int_0^1 DF(x_0+t h)\,h\,dt.
\tag{2.1}
\]

即使两端数据相同，也只能得到路径平均 Jacobian 作用为零，不能推出 \(J_0h=0\)。

**命题 1A：有限差量的强方向界。**  
假设连接两点的线段位于 \(F\) 的 \(C^2\) 区域，且
\(\|D^2F(x)\|\le M_2\)。令 \(P_s\) 投影到 \(J_0\) 奇异值至少为 \(s>0\) 的右奇异子空间，则

\[
\boxed{
\|P_s h\|
\le
\frac{\|F(x_*)-F(x_0)\|+\frac12M_2\|h\|^2}{s}.
}
\tag{2.2}
\]

证明：Taylor 余项满足 \(\|R_2(h)\|\le M_2\|h\|^2/2\)，而
\(s\|P_s h\|\le\|J_0 h\|\)。代入
\(J_0h=F(x_*)-F(x_0)-R_2(h)\) 即得。这里的范数及阈值必须使用声明过的材料/数据度量。

因此，若 \(\Delta F=0\)，强分量至多为二阶量；它不必等于零。最简单的例子是

\[
F(x,z)=x+z^2,\qquad (x,z)=(-t^2,t).
\]

沿 fiber 走，冻结弱方向是 \(z\)，有限位移却必须带上 \(x=-t^2\) 的强方向补偿。把任何强分量都认作“尚未恢复可观测信息”是不正确的。

若 \(h\) 不小、\(M_2\) 大或谱阈值 \(s\) 很小，式 (2.2) 可以没有实际约束力。它不是 GT 可达证书。

### 2.2 先验切空间才是正确的局部对象

假设 \(g\) 在去 gauge 后为 \(d\) 维光滑 immersion。应研究

\[
A_z=D(F\circ g)(z)=DF(g(z))Dg(z),
\tag{2.3}
\]

而不是直接把整个复像素空间 \(\ker J\) 都解释为合法材料自由度。

若 \(\operatorname{rank}A_z=r\) 在邻域内恒定，则标准常秩定理给出局部 fiber 维数 \(d-r\)。这是已有数学结果。尤其，当 \(r=d\) 时，局部精确解是孤立点；观测维数大于 \(d\) 并不妨碍这一结论。

最新两份 TRAIN smooth chart，八个生成参数中有一个相同材料 gauge，七个独立材料方向均有非零响应奇异值。该数值事实支持“当前点附近没有被观察到的精确先验内 null”，但单点浮点 SVD 不是整个邻域常秩的严格证书。[R1,R2]

### 2.3 “敏感信息已恢复”需要一个可检验定义

单说 full-GN 收敛、材料 MAE 小或残差低，都不等价于恢复所有可靠可观测信息。一个可操作定义至少包括：正确先验支持、正确完整 likelihood、驻点/KKT 误差、可观测方向的局部曲率与置信尺度，以及正向模型误差。

如果用强坐标 \(q(x)\) 表示局部观测可确定的商空间信息，那么“强信息准确”应定义为 \(q(x_0)\) 与相容真值集合的 \(q\) 接近，而非相对某个冻结像素 SVD 的坐标恰好一致。这样的商坐标通常本身是非线性的。

现有 physics_R_base 只有一次 ROM 更新，不能拿来否定“理想物理反演已经恢复敏感信息”的思想实验。该前提在最新诊断中没有完整实现。[R1]

---

## 3. Q2：GT 是否在同一局部 fiber？一个 Maxwell-form 反例

### 3.1 足够条件与不能保证的部分

足够条件可以明确写出：\(x_0,x_*\) 都在同一个正规 prior chart 中，\(F(x_0)=F(x_*)\)，并且都落在该 chart 的单个隐函数邻域内。在这个邻域中，fiber 是图像，二者之间存在合法局部路径。

其中每个条件都重要。若仅知道两端材料接近、数据接近，既不能推出相同精确 fiber，也不能推出同一局部分支。连通性是一个有限区域性质，不是某一点 Jacobian 的性质。

### 3.2 被动、互易两单元代数模型

考虑具有 Maxwell VIE 形式的两单元系统：

\[
G=
\begin{pmatrix}0&\gamma\\ \gamma&0\end{pmatrix},
\quad e=\binom10,\quad
X(q)=\operatorname{diag}(a+q,a-q),
\quad S=(0,1).
\]

固定 \(a\) 满足 \(\operatorname{Re}a>|q|\)、\(\operatorname{Im}a>0\)；两单元损耗为正。\(G\) 为对称耦合，假设所用区间分母不为零。直接解得

\[
F(q)=j_2(q)=
\frac{\gamma(a^2-q^2)}
{1-\gamma^2(a^2-q^2)}.
\tag{3.1}
\]

取合法先验为“总 contrast 固定、两个单元之间重新分配”的一维材料线。则

\[
F(q)=F(-q).
\]

对非零且足够小的 \(q\)，\(q\) 与 \(-q\) 是两个分离的局部解；在这条 prior 线上，保持同一数据不能连续穿过 \(0\)。两幅材料的距离为 \(2\sqrt2|q|\)，可以任意小。

所以：

\[
\boxed{\text{材料很近 + 数据完全相同，也不保证 prior 内 fiber 局部连通。}}
\]

这是一个可解析的 Maxwell-form 离散代数反例，不是已验证的连续三维实验。源向量是抽象可控激励；本报告不声称它已经对应某一既定远场测量装置。它足以反驳只依赖给定算子形式和“材料近”的普遍命题。

有限噪声时，若 \(F(0)\) 与观测间的差异落入所选容差，两个区域可能由 noise tube 连接；若不落入，则仍分离。连接与否取决于容差和 prior，而非“用了几阶反馈”。

### 3.3 边界、离散几何与 rank change

材料 box 的内部可用切空间描述；在边界只能使用可行切锥。被动性例如 \(\operatorname{Im}\chi\ge0\)，其边界上的某些双向切向运动不合法。对真实频率扫描，还要保证 dispersion 的物理相容性，不能在每个频率随意独立改变材料来制造歧义。

Hard rasterization 常导致参数小变而像素完全不变，随后突然跳格。这种导数为零不是一族不同材料的精确不可观测性。相同材料的重复坐标、Gaussian 权重整体缩放和 decoder 非单射也必须先排除。

Rank change、fold、分支交汇会破坏单 chart 的正规结构。Pseudo-arclength continuation 可穿过某些普通参数化的 turning point，却不能穿过不存在的合法分支，也不能让一个孤立点自动长出 fiber。

---

## 4. 完整 Maxwell 有限反馈：从精确恒等式到高阶导数

### 4.1 基本假设与精确有限响应

固定 \(G,S,e\)，允许 \(X(x)\) 是线性或非线性本构/离散极化率映射。设

\[
L_0=I-X_0G,\quad R_0=L_0^{-1},\quad
j_0=R_0X_0e,\quad E_0=e+Gj_0.
\]

对于有限变化 \(\Delta X=X(x_0+h)-X(x_0)\)：

\[
L_h(j_0+\Delta j)=(X_0+\Delta X)e,
\quad L_h=L_0-\Delta XG.
\]

减去背景方程：

\[
L_h\Delta j=\Delta X E_0.
\tag{4.1}
\]

故只要 \(L_0,L_h\) 可逆，

\[
\boxed{
\Delta j
=(I-R_0\Delta XG)^{-1}R_0\Delta X E_0.
}
\tag{4.2}
\]

记 \(K=R_0\Delta XG,\ b=R_0\Delta X E_0\)，则
\(\Delta j=(I-K)^{-1}b\)。观测变化为 \(S\Delta j\)，最后按预定顺序实化和白化。多源分别有 \(E_{0,p}\)，共用同一个 \(\Delta X\)。

**\(R_0\) 已经包含原材料背景中的完整多重散射。** 在它上面只保留 \(b\) 是“相对于当前背景的一阶材料变化”，不等于忽略了背景内部的全部多次散射。

### 4.2 Neumann 级数及其真正适用范围

\[
\Delta j=\sum_{k=0}^{\infty}K^k b
\]

在例如 \(\|K\|<1\) 时成立。截取 \(p\) 项的精确余项恒等式是

\[
\Delta j-\sum_{k=0}^{p-1}K^kb
=(I-K)^{-1}K^p b.
\tag{4.3}
\]

若 \(q=\|K\|<1\)，则

\[
\|\mathrm{remainder}\|
\le \frac{q^p}{1-q}\|b\|.
\tag{4.4}
\]

有限维下 \(\rho(K)<1\) 足以保证算子 Neumann 级数收敛；\(\|K\|\ge1\) 只表示这个特定范数充分条件失效。非正规矩阵可以有很大的范数和较小的谱半径。反过来，\(I-K\) 可逆也不意味着级数收敛，例如 \(K=2I\)。

最新诊断中的 1.129 与 1.032 是所用范数**上界**，不是已测实际范数越过 1；二十阶收敛不与之矛盾。[R1]

Padé 或 Schur 可以在某些多项式展开失效的范围内提供良好有理近似，但可能有伪极点，不能无条件稳定。2024 年已有用 Padé 处理发散 vectorial electromagnetic Born series 的光栅工作；“用有理近似重求和电磁多重散射”不是新的通用概念。[L8]

### 4.3 一阶、二阶、三阶与一般阶导数

以下 \(R=L(x)^{-1}\)、\(E=e+Gj(x)\)，所有方向导数对实材料坐标取，复数只是状态表示。记 \(X_u=DX[u]\)。

\[
Dj[u]=R X_u E.
\tag{4.5}
\]

若 \(X\) 非线性，二阶为

\[
D^2j[u,v]
=R X_{uv}E+
R X_uG\,Dj[v]+R X_vG\,Dj[u].
\tag{4.6}
\]

在线性 \(X\) 下第一项消失，从而

\[
D^2j[u,v]
=R X_u G R X_vE+R X_v G R X_uE.
\tag{4.7}
\]

三阶的非线性本构形式为

\[
\begin{aligned}
D^3j[u,v,w]=R\{&
X_{uvw}E+
X_{uv}G Dj[w]+X_{uw}G Dj[v]+X_{vw}G Dj[u]\\
&+X_uG D^2j[v,w]+X_vG D^2j[u,w]+X_wG D^2j[u,v]\}.
\end{aligned}
\tag{4.8}
\]

在线性 \(X\) 下，它展开成六个有序乘积；不能把不同方向的 \(X_uG\) 交换顺序。

一般情况下，把方向标签组成集合 \(A\)，\(X_B\) 表示标签集合 \(B\) 对应的高阶导数，则

\[
D_Aj
=R\left[X_AE+
\sum_{\varnothing\neq B\subsetneq A}X_BG D_{A\setminus B}j\right].
\tag{4.9}
\]

同一方向 \(v\) 重复 \(n\) 次时：

\[
j^{(n)}
=R\left[X^{(n)}E+
\sum_{k=1}^{n-1}\binom nk X^{(k)}G j^{(n-k)}\right].
\tag{4.10}
\]

若 \(X\) 线性，

\[
D^n j[v,\ldots,v]
=n!(R X_vG)^{n-1}R X_vE.
\tag{4.11}
\]

所有观测导数在这些表达式左侧乘 \(S\)。导数阶数升高不意味着必须显式储存高阶张量：给定方向的递推只需逐阶解线性系统；已有分解可重复使用。相同方向的二、三阶响应各增加一组源 RHS，而不是显式形成全部 \(n!\) 项。

A17 的三维 DDA 使用材料相关极化率。不能把二维 \(X=\operatorname{diag}\chi\) 的零高阶本构导数直接搬过去；需要保留实际 \(DX,D^2X,D^3X\)。[R6]

### 4.4 材料 decoder 的链式法则不能省略

对 \(f=F\circ g\)：

\[
D^2f[v,v]=D^2F[Dgv,Dgv]+J\,D^2g[v,v].
\tag{4.12}
\]

三阶为

\[
D^3f[v,v,v]
=D^3F[w,w,w]+3D^2F[w,D^2g[v,v]]+J D^3g[v,v,v],
\quad w=Dgv.
\tag{4.13}
\]

因此，“给 Maxwell 加到三阶”并不等于“给 prior 参数路径加到三阶”。有限参数变化应首先形成实际材料 secant，再决定是否需要高阶物理近似。最新 weak chart 的巨大线性失真恰好说明了这一点。[R1,R2]

### 4.5 正向反馈与 inverse Born 不同

这里从 \(\Delta X\) 计算 \(\Delta F\) 的是正向有限反馈。Inverse Born series 则试图从数据扰动表示反向材料映射；它需要另选逆/正则化分支，不能自动表示一整条非唯一 fiber。

电磁正向及 inverse Born 的收敛条件已有专门研究。[L7] 2026 年 8 月出现了 inverse-Born-based neural operators 的新预印本；摘要描述的是学习非线性逆映射的 operator expansion，并没有据此建立材料多解动力学或 posterior kernel。本次未取得其全文，因此这里只把它列为必须继续核对的直接新颖性风险。[L9]

---

## 5. Dyson、Magnus 与静态场：不存在凭空产生的路径记忆

静态 Neumann 是固定算子 \(K\) 的几何级数。时间/路径有序 Dyson 是变系数微分方程的解，两者不能混称。

令

\[
A(t)=I-GX(x(t)),\qquad E(t)=A(t)^{-1}e.
\]

假设整个路径上可逆，且 \(x(t)\) 光滑、源 \(e\) 固定。微分得到

\[
\dot E(t)=\mathcal H(t)E(t),\qquad
\mathcal H(t)=A(t)^{-1}G\dot X(t)=G R(t)\dot X(t).
\tag{5.1}
\]

其形式解确实可以写成时间有序 Dyson 展开，或者 Magnus 指数展开；例如 Magnus 第二项包含不同时刻 \(\mathcal H\) 的交换子。它们属于变系数线性 ODE 的已有方法。[L6]

但这里存在更强的精确表达：

\[
\boxed{
U(t,s)=A(t)^{-1}A(s),\qquad E(t)=U(t,s)E(s).
}
\tag{5.2}
\]

证明只需微分 \(A(t)^{-1}A(s)\)，并验证 \(U(s,s)=I\)。因此任意两条合法路径只要端点相同，得到的完整静态场相同；闭合路径的精确传播为 \(I\)。

等价地，令材料坐标连接矩阵
\(\mathcal H_i=A^{-1}G\,\partial_iX\)，则在这一约定下

\[
\partial_i\mathcal H_j-\partial_j\mathcal H_i-
[\mathcal H_i,\mathcal H_j]=0.
\tag{5.3}
\]

这不是宣布“所有材料几何都平坦”。它只说这个由可逆静态场方程诱导的场传输是纯坐标变换。Fiber 的第二基本形式、先验曲率、切空间转动以及探索算法的路径依赖仍然可以很大。

**研究含义：** 不应因为看到非交换 Dyson 项，就推断存在额外的物理历史自由度或自然跨模式机制。若材料有滞回、时变本构、真实动态激励或求解不唯一，则模型已不同，本结论不再直接适用。

---

## 6. T1/T2：二阶障碍、先验内在几何与有限可延拓性

### 6.1 二阶条件的严格含义

设 \(f:\mathbb R^d\to\mathbb R^m\) 为 \(C^2\)，路径
\(z(t)=z_0+t v+\frac12t^2a+o(t^2)\) 满足 \(f(z(t))=f(z_0)\)。令 \(A=Df(z_0)\)，则

\[
Av=0,\qquad Aa=-D^2f[v,v].
\tag{6.1}
\]

令 \(Q_{\rm out}=I-AA^\dagger\)，即 \(\operatorname{Range}A\) 的正交补投影，则

\[
\boxed{Q_{\rm out}D^2f[v,v]=0}
\tag{6.2}
\]

是二阶必要条件，也是该二阶线性补偿方程可解的充要条件。它只是**二阶 jet 的可解性**，不单独保证真实非零曲线存在。

例如 \(f(x,z)=(x,z^2)\)，在原点 \(z\) 是 Jacobian null 方向，但二阶项落在余核，精确 fiber 只有原点。若换成 \((x,z^3)\)，二阶条件通过，三阶仍失败。

三阶必须满足

\[
A b=-D^3f[v,v,v]-3D^2f[v,a],
\tag{6.3}
\]

且右端的余核投影为零。更高阶继续递推。在一般 \(C^\infty\) 情形，所有有限阶条件甚至全部 Taylor 系数都通过，也不能排除 flat-function 障碍；例如 \(e^{-1/z^2}\) 的光滑延拓。这个 flat 反例不是解析 Maxwell 本体的例子，但可以警示任意光滑 decoder。解析情形也不能把“检查到第三阶”当作充分存在性证明。

### 6.2 命题 2：把坐标曲率与真实 prior 曲率分开

以下在质量加权后的欧氏材料空间中推导。令 \(g\) 为 \(C^2\) immersion，

\[
Q_z=Dg(z),\quad H_z=Q_z^\top Q_z,\quad
P_{\mathcal P}=Q_zH_z^{-1}Q_z^\top.
\]

对参数速度 \(v\)，材料速度 \(w=Q_zv\)。定义

\[
\Gamma_z(v,v)=H_z^{-1}Q_z^\top D^2g[v,v],
\]

\[
\mathrm{II}_{\mathcal P}(w,w)
=(I-P_{\mathcal P})D^2g[v,v].
\tag{6.4}
\]

于是精确分解

\[
D^2g[v,v]=Q_z\Gamma_z(v,v)+
\mathrm{II}_{\mathcal P}(w,w).
\tag{6.5}
\]

第一项在 prior 切空间内，是坐标加速度；第二项是 prior 在材料空间中的第二基本形式。

选择参数加速度 \(a=-\Gamma_z(v,v)\)，即在该点沿诱导材料度量取 geodesic 二阶运动，材料加速度便只剩 \(\mathrm{II}_{\mathcal P}\)。沿这条路径的二阶数据响应为

\[
\boxed{
C_{\rm int}(w)=D^2F[w,w]+J\,\mathrm{II}_{\mathcal P}(w,w).
}
\tag{6.6}
\]

这说明应把此前的两项升级成三个报告项：

\[
\text{Maxwell feedback}
+\text{真实 prior 外在曲率}
+\text{纯坐标加速度}.
\]

而不是只比较 \(D^2F\) 与 \(JD^2g\) 的范数。

更重要的是，对 \(A=JQ_z\) 的余核，

\[
Q_{\rm out}D^2(F\circ g)[v,v]
=Q_{\rm out}C_{\rm int}(w),
\tag{6.7}
\]

因为 \(Q_{\rm out}JQ_z\Gamma=0\)。因此二阶障碍在正确处理速度后具有坐标不变性，巨大 raw decoder 项却不一定具有这种意义。

**证明要点：** 式 (6.5) 是材料切空间的正交分解；链式法则代入 \(z''=-\Gamma\) 得 (6.6)；再利用余核投影消去 \(A\Gamma\) 即得 (6.7)。非欧氏但固定材料质量矩阵时，把转置内积改成相应质量内积即可。

### 6.3 一个直接反驳“decoder 大曲率必然难学”的例子

取

\[
g(u,v)=(u+c v^2,v),\qquad F(x_1,x_2)=x_1.
\]

先验 image 是整个平面，真实 \(\mathrm{II}_{\mathcal P}=0\)，但原点的
\(J D^2g[(0,1),(0,1)]=2c\)，可以任意大。

固定坐标直线 \((u,v)=(0,t)\) 造成数据漂移 \(ct^2\)；改成
\((u,v)=(-ct^2,t)\)，材料路径就是 \((0,t)\)，数据严格不变。

所以，最新诊断的 decoder/Maxwell 比值 38361 和 42.4，**不能单独证明真实 prior 流形具有同量级的不可消除曲率**。这不否认已有诊断，而是指出下一步应从其缓存中进一步分解哪些弯曲来自 Gaussian-CDF 坐标、峰值归一化和真实材料 image。[R1,R2]

### 6.4 Fiber 在材料空间中的第二基本形式

在 prior 切空间取正交 frame \(Q\)，令 \(A=JQ\)。对 \(w\in T_x\mathcal P\cap\ker J\)，若二阶兼容条件成立，则选择最小 prior 法向补偿得到

\[
\boxed{
\mathrm{II}_{\mathcal F}(w,w)=
\mathrm{II}_{\mathcal P}(w,w)
-
Q A^\dagger
\left(D^2F[w,w]+J\mathrm{II}_{\mathcal P}(w,w)\right).
}
\tag{6.8}
\]

这里第一项垂直于 prior，第二项位于 prior 内部、垂直于 fiber 切空间。正规 fiber 存在时，这正是材料空间中的 fiber 第二基本形式。若余核兼容失败，式右侧只能给最小二阶残差修正，不能被称为真实 fiber 曲率。

**这比把 Maxwell Hessian 单独叫“流形曲率”更准确。** 真正决定可延拓性的对象是物理反馈、材料先验外在几何和可观测 normal 逆共同组成的表达式。

### 6.5 命题 3：非零余核曲率给出有限可行步长障碍

假设一族 \(C^3\) prior 合法路径在 \(z_0\) 具有相同非零切向 \(v\in\ker A\)，其加速度和三阶余项有统一界。设

\[
\kappa=\|Q_{\rm out}D^2f[v,v]\|>0,
\]

且投影余项范数不超过 \(C_3|t|^3\)。任何这样的路径均满足

\[
\|f(z(t))-f(z_0)\|
\ge \frac{\kappa}{2}t^2-C_3|t|^3.
\]

当 \(|t|\le\kappa/(4C_3)\) 时，

\[
\|f(z(t))-f(z_0)\|\ge\frac{\kappa}{4}t^2.
\tag{6.9}
\]

因此局部数据容差 \(\epsilon\) 限制

\[
|t|\le 2\sqrt{\epsilon/\kappa},
\tag{6.10}
\]

但该结论只适用于声明的局部半径和统一路径正则性范围，不能当作所有远距离补偿路径的全局上界。

对近似弱向量，取最小二阶 normal correction，得到实用上界

\[
\|\Delta f(t)\|
\le |t|\|Av\|
+\frac12t^2\|Q_{\rm out}C_{\rm int}\|
+C_3|t|^3+\epsilon_{\rm solve}.
\tag{6.11}
\]

这给出一阶漂移、不可补偿二阶弯曲、高阶误差与 solver 误差的明确分工。


## 7. T3：什么时候局部探索可以稳定学、稳定走？

### 7.1 命题 4A：一个带明确半径的有限补偿 chart

考虑 \(f:\mathbb R^d\to\mathbb R^m\)，在 \(z_0\) 满行秩。取正交基
\(T\) 张成 \(\ker Df(z_0)\)，\(N\) 张成其正交补，令

\[
B=Df(z_0)N,\qquad \sigma_{\min}(B)=\beta>0.
\]

假设 \(Df\) 在半径 \(\rho\) 内为 \(L\)-Lipschitz。对

\[
r\le \min\{\rho/2,\beta/(4L)\}
\]

以及 \(\|s\|\le r\)，存在唯一的小 normal correction \(a(s)\)，使

\[
f(z_0+Ts+Na(s))=f(z_0),\qquad
\|a(s)\|\le \frac{2L}{\beta}\|s\|^2.
\tag{7.1}
\]

“唯一”只指这里指定的局部 normal 球，不是全局唯一材料解。\(L=0\) 时取相应线性极限。

**证明。** 定义

\[
\mathcal T_s(a)=-B^{-1}
[f(z_0+Ts+Na)-f(z_0)-Ba].
\]

在 \(\|a\|\le2L\|s\|^2/\beta\) 的球上，Taylor 余项使该映射保持球内，导数范数至多 \(2Lr/\beta\le1/2\)。压缩映射定理给出唯一解。该论证是经典隐函数定理的定量版本，不应单独包装成新的 Maxwell 定理。

若 \(f\) 只恒秩 \(r<m\)，可在已证明常秩的局部选择独立输出坐标；若没有常秩依据，不能随意删去余核输出后宣布得到精确 fiber。

这里可学习的有限对象其实非常明确：**\(a(s)\)，即给定合法切向位移所需的 normal compensation。** 它不需要知道 GT。若 \(\beta\) 很小或 \(L\) 很大，则 chart 半径小、补偿导数大；“学习固定 anchor 的全局函数”因此可能不合适，需要重心刷新与 atlas。

### 7.2 谱稳定性：exact null 与人为 weak cutoff 不一样

对于相同秩的两个矩阵 \(A,\widetilde A\)，若 \(A\) 的最小非零奇异值至少为 \(\beta\)，且

\[
\|\widetilde A-A\|\le\epsilon_A<\beta/2,
\]

则其 exact-null 投影满足标准子空间扰动界

\[
\|P_{\ker A}-P_{\ker\widetilde A}\|
\le 2\epsilon_A/\beta.
\tag{7.2}
\]

可以直接用
\(\|\widetilde A^\dagger\widetilde A P_{\ker A}\|
\le \|\widetilde A^\dagger\|\epsilon_A\)
及对称方向证明。

若投影的是“奇异值低于人工阈值的 weak cluster”，则需要的是该 cluster 与其余谱之间的分离。对 Gram 矩阵，误差为

\[
\|\widetilde A^\top\widetilde A-A^\top A\|
\le (\|A\|+\|\widetilde A\|)\epsilon_A.
\tag{7.3}
\]

投影稳定性受此误差与 cluster eigen-gap 的比控制。小的 JVP 平均误差或良好的 GN 更新，并不足以控制它。

反例：

\[
A_\epsilon=(\epsilon,0),\qquad
\widetilde A_\epsilon=(0,\epsilon).
\]

矩阵差趋于零，null 方向却正交。失效原因正是没有统一谱隙。对近重根的单个奇异向量做监督还存在任意旋转/定号问题；更合理的是学习子空间投影、坐标无关的材料向量，或明确 gauge 后的 frame-equivariant 表示。

### 7.3 一个不需要硬谱隙的实用替代

在 prior 的材料正交坐标中，定义软弱滤波器

\[
W_\lambda(A)=\lambda(A^\top A+\lambda I)^{-1},
\qquad \lambda>0.
\tag{7.4}
\]

这是标准 Tikhonov 型谱滤波，不是新定理。它避免了 cutoff 处的跳变，且有明确的泄漏界：

\[
\boxed{\|A W_\lambda(A)\|\le\sqrt{\lambda}/2.}
\tag{7.5}
\]

证明：对每个奇异值 \(\sigma\)，滤波后的响应为
\(\sigma\lambda/(\sigma^2+\lambda)\)，其最大值是 \(\sqrt\lambda/2\)。

由 resolvent identity，

\[
\|W_\lambda(A)-W_\lambda(B)\|
\le \frac{\|A\|+\|B\|}{\lambda}\|A-B\|.
\tag{7.6}
\]

它不需要人为 weak cluster 的谱隙，但 \(\lambda\to0\) 时仍会变得敏感。这是一种清楚的偏差—稳定性交换：不承诺 exact-null 运动，换取受控有限数据漂移。

令材料速度

\[
\dot x=Q_x W_\lambda(A_x)b_\theta(x,y,\xi),
\tag{7.7}
\]

其中 \(Q_x\) 是 prior 切空间正交 frame，\(\|b_\theta\|\le B_0\)，并确保运动始终留在合法 prior 内部。则连续时间的观测漂移有

\[
\|F(x(T))-F(x(0))\|
\le T\,B_0\sqrt\lambda/2.
\tag{7.8}
\]

这里每一时刻使用当前 \(A_x=J_xQ_x\)。冻结初始 \(A\) 并不能获得同样保证；离散实现还要加入重traction、模型误差和积分误差。

### 7.4 命题 4B：学习场的有限时间稳定性

在一个固定正规 chart 域内，设理想 prior 合法切向场 \(V\) 为 \(L_V\)-Lipschitz，学习/近似场 \(\widehat V\) 的材料范数误差统一不超过 \(\epsilon_V\)。两条解在时间 \(T\) 前都没有离开该域、撞到边界或穿越 rank-change，初始误差为 \(e_0\)。则

\[
\|\widehat x(t)-x(t)\|
\le e^{L_Vt}e_0+
\frac{\epsilon_V}{L_V}(e^{L_Vt}-1).
\tag{7.9}
\]

当 \(L_V=0\) 时取 \(e_0+t\epsilon_V\)。这是 Grönwall 的直接应用。对精确 null 投影构造的场，\(\epsilon_V\) 中必须包含
\(O(B_0\epsilon_A/\beta)\)；\(L_V\) 也受到 prior frame 变化、Jacobian 变化和谱间隙影响。

若另有

\[
\|DF(x)\widehat V(x)\|\le\eta
\]

及每步重traction/求解误差记录，则连续数据漂移至多 \(T\eta\)，离散项另计。使用二阶补偿 predictor 可以在 \(C^3\) 条件下把单步数据余项从 \(O(h^2)\) 降为 \(O(h^3)\)，但这不是无需 corrector 的全局稳定保证。

**结论：** 可学习性需要一个范围明确的、条件良好的局部问题。它不要求流形平坦；要求曲率、谱/软滤波条件数、材料参数化及模型误差在部署域内可控。万能逼近定理不提供这些条件。

### 7.5 高曲率不等于没有泛化，低曲率也不等于全局安全

曲率大可能只要求更小步长、重新定中心或多个 chart。真正困难的是：两条分支很接近、prior support 不光滑、局部逆不稳定、场约束突然失效，或训练未覆盖该区域。

Reach 还取决于不同区域之间的距离；不能只拿 \(1/\|\mathrm{II}\|\) 当作全局 reach。两条平行、各自零曲率的分支距离很近，整体 reach 仍可很小。

两个映射 \(F_\pm(x,z)=x\pm c z^2\) 在原点的 Jacobian 完全一样，有限 normal compensation 却相反。仅向网络提供冻结 Jacobian 不足以识别有限传输。提供当前材料、内部场或有限反馈动作，可以补充**计算表示**；但它们仍来自同一观测/假设状态，不是新增测量信息。

---

## 8. 一个本质不同的备选：把精确材料 fiber 写成电流兼容集合

用户的目标不必永远从 material tangent 出发。在线性、各向同性 \(X=\operatorname{diag}\chi\) 下，有一个精确的有限状态表示值得单独研究。

### 8.1 命题 5：多源电流变化的精确材料 realizability

令 \(k_p=j'_p-j_{0,p}\)，把所有源排成矩阵 \(k\)。相同数据需要

\[
Sk_p=0\quad\text{对所有 }p.
\tag{8.1}
\]

状态方程等价于

\[
L_0k=\operatorname{diag}(h)(E_0+Gk).
\tag{8.2}
\]

记每个材料单元的多源行向量为

\[
a_i(k)=(L_0k)_{i,:},\qquad
b_i(k)=(E_0+Gk)_{i,:}.
\]

若 \(b_i\ne0\)，存在单一材料增量 \(h_i\) 同时解释所有源，当且仅当

\[
a_i=h_i b_i.
\tag{8.3}
\]

这等价于所有二阶 minors 为零：

\[
a_{i,p}b_{i,q}-a_{i,q}b_{i,p}=0
\quad\text{对所有源对 }p,q.
\tag{8.4}
\]

此时材料由

\[
h_i=\frac{a_i b_i^*}{\|b_i\|^2}
\tag{8.5}
\]

恢复。式 (8.4) 对 \(k\) 至多是二次多项式，而 (8.5) 是有理映射。

**证明。** 式 (8.2) 是有限差分状态方程的直接重排；行向量共线的必要充分条件给出 (8.4)，在 \(b_i\ne0\) 时相应投影系数唯一。再检查新材料在 prior 中及 \(L(\chi_0+h)\) 可逆，即可从这些代数约束还原一个真正的 forward solution。反向从任意合法同数据材料出发，其完整状态自然满足这些等式。

若 \(b_i=0\)，必须有 \(a_i=0\)；该处材料不由此局部除法确定，应更换 chart 或保留未消去变量，不能加一个 epsilon 后假装数学等价。

### 8.2 它为何重要，又为何可能很差？

这给出“receiver-null 电流 + 多源共同材料兼容 + prior”的有限等数据描述。所有多重散射已经在方程内，不是二阶 Born 截断。因此，有时可以在电流域做约束 continuation，而不直接积分材料 nullspace。

但**任意 receiver-dark 电流绝不自动对应合法材料**。跨源共线、被动性、固定相位、空间先验和 \(L\) 可逆性都是独立约束；它们可能把表面上巨大的暗电流空间压成孤立点。各向异性 tensor 材料应改成相应矩阵匹配问题；DDA 极化率需先恢复极化率增量，再合法反演本构映射。

这与 contrast source inversion 的双变量状态/材料思想有明确祖先关系。CSI 早已利用数据方程与本构状态方程联合优化；本报告不声称首次发现 bilinear lifting。[L10] 这里值得研究的不同点是：**用于精确同数据解集的多源代数消元、prior 合法性、局部分支与有限传输**，而不是再次提出一种通用 CSI 重建器。

计算上，电流维数是材料维数乘源数，可能比材料 chart 昂贵得多。只有压缩后的 dark-current 坐标保留了多源兼容性，并真正降低有限解搜索成本，这条备选才值得推进。

---

## 9. T6：A17 可以贡献什么，不能贡献什么？

### 9.1 原 A17 的正确位置

**符号区分：** 本报告的 \(R_0=L_0^{-1}\) 是完整背景 resolvent；A17 原文中的旧 Galerkin 算子在这里应另记 \(R_V=V(V^*LV)^{-1}V^*\)，不能混同。

SOM/TSOM 的电流分解也不能直接等同于材料 \(DF\) 的 kernel。材料方向须先经过当前总场诱导的 injection，再经内部 current response 和接收映射；这正是 A17 所强调的任务联系。[R6]

A17 在固定材料状态、固定注入与共同二次目标上，对指定 current-space 改变给出

\[
\Delta J=Y_C\Gamma_C^{-1}N_C.
\tag{9.1}
\]

它把材料激发、内部反馈和接收返回联系起来；由此得到 paired material update-space。它告诉我们：**某种 current 表示改变如何传递到材料优化任务**，而不是“这些电流方向就是解流形自由度”。[R6]

最直接的反例是 \(r=0,\ell=0\)：所有正则 GN 步都为零，即使两个 Jacobian 的 kernel 完全不同。因而对这一个残差的 update fidelity 完美，仍可以对探索几何毫无控制。

A17 的 Schur inverse 也不是“Maxwell 只做到二阶”。它对所保留/增广的 current block 消元可包含无限次线性反馈；“冻结二次”说的是材料优化目标，而不是整个内部多次散射只保留两次。这一区别不要求把原论文推翻。

### 9.2 有限材料变化的双侧接口

对有限变化写

\[
L_h k=b_h,\qquad b_h=\Delta X E_0.
\]

给定近似电流变化 \(\widehat k\)，其完整方程残差

\[
\rho=b_h-L_h\widehat k.
\]

若 \(Z\) 满足 \(L_h^*Z=S^*\)，则

\[
\boxed{\Delta F-S\widehat k=Z^*\rho.}
\tag{9.2}
\]

这是精确 goal-oriented residual identity。若只有近似伴随 \(\widehat Z\)，并记
\(\eta=S^*-L_h^*\widehat Z\)，则

\[
\Delta F-S\widehat k
=\widehat Z^*\rho+\eta^*L_h^{-1}\rho.
\tag{9.3}
\]

若已验证 \(\sigma_{\min}(L_h)\ge\alpha_h>0\)，

\[
\left\|\Delta F-S\widehat k-\widehat Z^*\rho\right\|
\le\|\eta\|\|\rho\|/\alpha_h.
\tag{9.4}
\]

它把两侧截断误差、有限材料改变和目标观测直接连接。A17 型 Schur 或 primal/dual basis 可以用来近似这些因素，但上述恒等式本身属于既有残差/伴随代数，不是自动的原创性。

关键难点是 \(\alpha_h\) 和完整残差是否便宜可得。仅仅观察到 reduced core 条件数良好，不能代替完整 \(L_h\) 的下界。可使用例如

\[
\sigma_{\min}(L_h)\ge\sigma_{\min}(L_0)-\|\Delta XG\|
\]

这样的局部充分界；失效时需要其他可靠界或 full fallback。界太保守也可能完全没有筛选价值。

### 9.3 如何连接后验，而不是只连接 forward error？

若完整白化预测误差有界 \(\|F-\widehat F\|\le\epsilon_F\)，令
\(\widehat r=\widehat F-y\)，则

\[
|\Phi-\widehat\Phi|
\le\|\widehat r\|\epsilon_F+\tfrac12\epsilon_F^2.
\tag{9.5}
\]

这可以用于 likelihood 区间、lazy MH 判定或 delayed acceptance。若区间无法确定接受/拒绝，就调用完整模型。没有被证明的误差上界时，不能把经验 surrogate error 当作 deterministic certificate。

几何控制还需要另外的 \(\epsilon_A\) 与二阶误差。对正常可逆的 normal block，曲率误差的典型尺度是

\[
O(\epsilon_2/\beta)+O(M_2\epsilon_A/\beta^2),
\tag{9.6}
\]

另加 prior 几何误差及切向对齐误差。这里隐含同秩、受控方向和 uniform bounds；不能由一个 GN gap 数字推出它。

### 9.4 什么才算 A17 增量成立？

需要同时满足三个条件：相同有限材料/likelihood 精度、相同目标分布或探索覆盖指标、计入离线准备与全部 full-state/JVP/VJP/RHS/拒绝成本后仍更好。

给定一方向，完整 Hessian-vector 二阶响应可以通过复用分解和少量 RHS 获得；小型 dense TM solver 已很便宜。此时 A17 维护候选空间、验证与更新的成本可能更高。DINO 一类 derivative-informed surrogate 已直接针对昂贵 PDE 后验几何加速，不能用一个低标准普通 MCMC 作为唯一对照。[L4]

**本次裁决：A17 有值得检验的接口，没有已经成立的新探索收益。** 不建议仅因本次高阶理论修改原稿的主结论；可增加一段明确的适用边界，并把 finite-transport 作为独立后续问题。

---

## 10. N3 的真实数学对象

### 10.1 源码核查结论

实际网络并不是抽象的 \(f_\theta(y)\)。它通过冻结 encoder 构造十二通道特征，并以一个物理 base 加残差：

\[
\widehat m
=
\operatorname{Box}\!\left(
m_{\rm base}
+s_{\rm TRAIN}\sqrt{\Delta V}\,
N_\theta(c_y)
\right).
\tag{10.1}
\]

五个卷积层、宽度 32、dilation \(1,2,4,2,1\)、GELU、零初始化输出 head，共 31,810 参数。输出是两通道实/虚完整材料增量，不限制在 weak space。Box 是
\(\operatorname{Re}\chi\in[-0.5,4]\)、\(\operatorname{Im}\chi\in[0,0.5]\)。[R3]

通道按 TRAIN-only pooled 均值/标准差归一化，输出两部分共享一个 TRAIN 材料增量 RMS。实际损失是

\[
\mathcal L
=
\mathbb E_{\rm minibatch}
\left[
\frac{\|\widehat{\Delta\chi}-\Delta\chi\|_{\rm pixel}^2}
{s_{\rm TRAIN}^2}
+0.05(1-\operatorname{SSIM}_\Delta)
\right].
\tag{10.2}
\]

SSIM 分别对实虚两通道做 \(7\times7\) 局部计算，动态范围使用 4.5 与 0.5，随后平均。不是只对材料幅值计算的 SSIM。训练图中没有 full-Maxwell、measurement loss 或在线反馈求解。[R4,R5]

训练基于 64 个 TRAIN parents、多背景/幅度/噪声条件；开发集用于 checkpoint 选择。它是有限数据、有限网络类、有限优化的经验风险最小化器，不是一个可以直接套用理想 Bayes 极限定理的对象。[R5]

### 10.2 条件均值结论到底要求什么？

若去掉 SSIM，采用固定正定权重的纯平方损失，允许所有可测预测器并达到总体最优，则

\[
f^*(C)=\mathbb E[\chi\mid C],
\tag{10.3}
\]

其中 \(C\) 是实际提供给网络的全部条件，包含特征、base 及已知背景等。

只有 \(C\) 对材料相对于 \(Y\) 充分，或
\(\mathbb E[\chi\mid C]=\mathbb E[\chi\mid Y]\) 恰好成立，才能把它写成给定原始数据的条件均值。确定性特征维数多不证明充分性；维数少也不单独证明丢失了目标相关信息。

N3 的 SSIM、box、模型容量、weight decay、有限数据与优化都使理想结论不能直接用作它的精确解释。

### 10.3 四个假说的证据等级

| 假说 | 可以说什么 | 尚不能说什么 |
|---|---|---|
| H1 条件平均 | 纯 MSE 总体极限支持这种机制 | 实际 N3 就等于正确 posterior mean |
| H2 局部 prior 投影 | 残差 CNN 可学习类似“向训练结构拉回”的行为 | 它是正交投影、MAP proximal 或保持 Maxwell fiber |
| H3 有限反演更新 | base+residual 的形式容许近似更新算子 | 更新具有 Newton/GN、下降、收敛或迭代一致性 |
| H4 分布依赖回归器 | 这是结构与训练目标直接支持的事实 | 它没有学到任何物理有用规律 |

区分 H1/H2/H3 需要额外实验：对可计算 posterior 的小问题比较条件均值；重复应用网络检查 idempotence/迭代行为；在 held-out 分布中检查相位、数据误差与 sensitivity-aligned response。只看一张好看的图不能区分。

### 10.4 为什么条件均值不一定数据相容？

即使 noiseless posterior 的每个样本都满足 \(F(\chi)=y\)，非线性通常导致

\[
F(\mathbb E[\chi\mid y])\ne
\mathbb E[F(\chi)\mid y]=y.
\]

在局部小协方差条件下，设 \(m=\mathbb E[\chi\mid y]\)，则

\[
F(m)-y\approx
-\frac12D^2F(m):\operatorname{Cov}(\chi\mid y),
\tag{10.4}
\]

并有更高阶余项。Gaussian noisy posterior 下，连
\(\mathbb E[F(\chi)\mid y]=y\) 也不是一般恒等式。

均值留在 prior support 的凸包，而不是必然留在 support 本身。二元材料、孔洞数、固定峰值归一化等非凸条件可能被平均破坏；凸 box 与非负材料分数则可保留。

### 10.5 固定相位提供了一个比“像均值”更强的反证

如果正确条件下的材料始终为

\[
\chi=\chi_b+c\rho,\qquad \rho\ \text{为实值、非负},
\]

且 \(\chi_b,c\) 已知，那么精确条件均值必为

\[
\mathbb E[\chi\mid C]
=\chi_b+c\,\mathbb E[\rho\mid C].
\tag{10.5}
\]

它必然保留固定增量相位。峰值 normalization 不一定保留，但固定相位必须保留。

因此实际 N3 的 phase violation 说明：“合法多世界平均”不是全部解释。SSIM 对实虚部分的非等变处理、函数类误差、有限样本、优化、输出约束与实际条件信息都可能贡献偏差；本次不能为它们分配精确百分比。只因为 N3 MAE 较低，也不能推出其偏差全来自真正的不可辨识性。[R1,R4]

### 10.6 一个精确的风险分解

在平方损失下，令 \(C=C(Y)\)，包含全部固定已知条件，并记 \(m_Y=\mathbb E[\chi\mid Y]\)、\(m_C=\mathbb E[\chi\mid C]\)，则

\[
\begin{aligned}
\mathbb E\|\chi-f(C)\|^2
={}&\mathbb E\|\chi-m_Y\|^2\\
&+\mathbb E\|m_Y-m_C\|^2\\
&+\mathbb E\|m_C-f(C)\|^2.
\end{aligned}
\tag{10.6}
\]

依次对应原观测下不可约 Bayes 风险、条件特征损失，以及实际预测器相对该条件均值的偏离。后者还包含目标函数偏移、容量、训练与优化等效应，不能无证据继续唯一细分。

**N3 的地位：保留 Top-1 基线与可选 warm start；不把“改善 N3”定义为 Explorer 的本质目标。**

---

## 11. T4：有限噪声 posterior、fiber 上的测度与正确随机核

### 11.1 首选统计目标

在原生 prior 坐标 \(z\) 中，

\[
\pi_\sigma(z\mid y)
\propto
p_Z(z)
\exp\left[-\frac{1}{2\sigma^2}
\|\Sigma_0^{-1/2}(F(g(z))-y)\|^2\right].
\tag{11.1}
\]

若已经把 \(\sigma\) 吸收进白化映射，则删去外面的 \(\sigma^{-2}\)，不能重复白化。所有讨论必须固定这一约定。

Bayesian posterior 不要求初始化先落进某个 95% residual 球。初始化应该合法且计算上有用，但好 sampler 可以从低密度点开始；代价是 burn-in 和混合。相反，精确 fiber continuation 必须从满足其等式定义的合法点出发。

### 11.2 正规小噪声极限的 coarea 因子

设 \(f:\mathbb R^d\to\mathbb R^r\) 在目标 fiber 附近满行秩，\(p_Z\) 连续且可积，fiber 正规，所需尾部/局部化条件保证极限归一化成立。则对 noiseless 数据 \(y_0\)，小噪声极限在 fiber 上的密度相对于其 Hausdorff 测度为

\[
\boxed{
d\pi_0(z\mid y_0)
\propto
\frac{p_Z(z)}
{\sqrt{\det(Df(z)Df(z)^\top)}}\,
d\mathcal H^{d-r}(z).
}
\tag{11.2}
\]

这是 coarea 的标准后果：normal 方向 likelihood tube 的体积并不处处相同。不能只在 fiber 上按 \(p_Z\) 或按均匀 arclength 采样，就声称获得小噪声 posterior。[L3]

若转到材料流形体积，且 \(g\) 单射、\(H=Dg^\top Dg\)，prior 的材料体积密度是 \(p_Z/\sqrt{\det H}\)，再除以限制在 prior 正交切空间上的 forward normal Jacobian。选择哪一种坐标都可以，但不能混用密度与参考测度。

若原始参数有 gauge，pushforward prior 要对同材料 preimage 做积分。简单固定一个 gauge 并宣布剩余七参数均匀，通常改变 posterior。Gaussian 权重归一化中的共同尺度尤其需要这样处理。[R2]

一个更简单的实施选择是：几何计算去掉 gauge，但概率采样暂时保留原始八参数及其已知 Gaussian-CDF prior，让相同材料的 preimage 通过原目标自然积分。只有推导了正确商空间密度，才把 sampler 真正降到七维。

若 rank 改变、存在不同维数 strata，或者数据维数高于有效 rank，不能不加说明地把 determinant 换成 pseudodeterminant 就完成概率理论。应明确局部独立观测坐标、目标空间测度和噪声投影；否则保留有限噪声式 (11.1) 最安全。

### 11.3 没有局部 fiber 也可以有全局多解 posterior

若 \(d\le m\)，多个分离点 \(z_k\) 都满足 \(f(z_k)=y_0\)，且每个 \(A_k=Df(z_k)\) 满列秩，在适当隔离与尾部条件下，小噪声 posterior 可趋近离散混合，其局部权重正比于

\[
p_Z(z_k)\,[\det(A_k^\top A_k)]^{-1/2}.
\tag{11.3}
\]

因此，当前七维 chart 满列秩不等于全局无多解。它只说明每个正规解附近没有非零 exact-fiber 自由度。要找其他分支，应做 global discovery，而不是无限延长一个不存在的局部 manifold。

### 11.4 Manifold lifting：可信的备选，而不是新造的概率概念

引入噪声变量 \(\eta\)，约束

\[
C(z,\eta)=f(z)+\sigma\eta-y=0.
\]

在 \(\sigma>0\) 时，该约束 Jacobian 含 \(\sigma I\)，因此满行秩。正确 lifted density 相对于约束流形的 Hausdorff 测度是

\[
\pi_{\rm lift}(z,\eta)
\propto
\frac{p_Z(z)\varphi(\eta)}
{\sqrt{\det(DC\,DC^\top)}}.
\tag{11.4}
\]

其材料边缘就是有限噪声 posterior。这是 Au、Graham、Thiery 已经研究的方法；他们使用包含体积项的 constrained HMC，并关注小噪声极限下的混合。[L3]

对当前 \(m\) 很大、prior dimension 很小且 full-column identifiable 的情形，lifting 不会神奇消除所有病态：当 \(\sigma\to0\)，约束 Gram 的某些方向仍可能退化；其计算维度也增加。应与原生低维 MH/SMC 比较，不预设 lifting 胜出。

### 11.5 Explorer 的随机性不是 posterior 正确性

若希望在正规 fiber 上采样密度 \(\rho\)，一种清楚的连续生成元是

\[
\mathcal L\phi
=
\rho^{-1}\operatorname{div}_{\mathcal F}
(\rho\nabla_{\mathcal F}\phi).
\tag{11.5}
\]

它揭示了几何体积与目标密度的作用。简单“投影高斯噪声 + 梯度”若忽略 Itô 曲率漂移、度量导数或边界条件，通常不是这个生成元。可用 Stratonovich 或 intrinsic 写法正确构造，但离散后仍需误差分析/校正。

有限噪声下最直接的可验证方案是：在有明确密度的原生 prior 坐标中使用 MH。给定提案 \(q(z'\mid z)\)，

\[
\alpha(z,z')
=
1\wedge
\frac{p_Z(z')e^{-\Phi(z')}q(z\mid z')}
{p_Z(z)e^{-\Phi(z)}q(z'\mid z)}.
\tag{11.6}
\]

必须保留拒绝后的重复状态。不能把 chain 中去重后的漂亮图按相同权重当作 posterior。

**特别危险的操作：** 先从高斯 proposal 采样，再做一个 many-to-one physics corrector，最后只用 likelihood ratio 接受。Corrector 改变了 proposal law；没有逆映射、Jacobian 或其他正确 augmented-state 推导，原 MH 比率不再适用。

### 11.6 局部概率与全局概率

如果为了计算方便只在固定区域 \(D\) 中运行正确采样器，得到的是

\[
\pi_D(\cdot\mid y)=\pi(\cdot\mid y,D),
\]

不是整个 posterior。除非另外估计 \(\pi(D\mid y)\)，不能把这个区域内部的归一化权重当作全局材料概率。

单张连续材料图通常具有零概率质量。应报告明确事件、材料邻域或结构/功能事件的 posterior mass，并保留未覆盖质量。Top-K 代表图不应重新归一化为 100%。


## 12. T5：GT 可达、覆盖概率、Bayes 风险与“学到了什么信息”

### 12.1 可达路径存在，不等于算法会找到它

若 GT 位于同一正规 fiber 的连通局部区域，可存在合法路径；但寻找它需要足够的方向覆盖、局部 chart 半径、数值稳定性和预算。一个固定向量场可能只沿某条闭合轨道运动，永远覆盖不到同一分支的大部分区域。

随机扩散若在连通紧致流形上具有合适的非退化切向噪声、边界条件和正密度，才可能具有相应遍历性质；有限预算的访问概率仍取决于 bottleneck 和混合。仅知道“连通”没有数值效率保证。

GN 的目标是降低一个标量误差。到达等数据解时，数据梯度通常消失；继续运行同一个 GN 并不会主动沿 null 方向探索。若另加 prior penalty，算法会倾向同一个 MAP 点，而不是在 posterior 上走动。需要显式 proposal/progress coordinate、随机机制或非局部跳跃。

### 12.2 命题 6A：不可辨识风险与样本风险

若两个合法材料 \(x_1,x_2\) 对观测诱导完全相同的分布，先验各为 \(1/2\)，令 \(D=\|x_1-x_2\|\)，则任何只使用该观测的点估计器的平均平方风险至少为

\[
D^2/4.
\tag{12.1}
\]

证明：条件后验始终是等权两点分布，平方损失最优点为中点，风险就是两点方差。观测随机化也不能降低该 Bayes 风险。

若两种观测分布为 \(P_1,P_2\)，一般两点检验论证给出

\[
\inf_{\widehat x}\frac12\sum_{i=1}^2
\mathbb E_i\|\widehat x-x_i\|^2
\ge \frac{D^2}{8}(1-\operatorname{TV}(P_1,P_2)).
\tag{12.2}
\]

对于同协方差高斯观测，等先验最优判别错误率是
\(\Phi_{\rm N}(-\delta/2)\)，\(\delta\) 为两均值间白化距离。这里 \(\Phi_{\rm N}\) 是标准正态 CDF，与 likelihood cost 不是同一个符号。

若 GT 自 posterior 产生，而独立 posterior sample \(\widetilde X\) 在给定 \(Y\) 后与 GT 条件独立，则

\[
\mathbb E[\|\widetilde X-X_*\|^2\mid Y]
=2\operatorname{tr}\operatorname{Cov}(X_*\mid Y).
\tag{12.3}
\]

条件均值的对应风险只有
\(\operatorname{tr}\operatorname{Cov}(X_*\mid Y)\)。因此，**正确 sampler 不应被要求必然击败接近条件均值的 N3 的 Top-1 MSE**。它解决的是表示不确定性、事件概率、覆盖与后续决策，不是违反 Bayes 风险下界。

对于 MAE，条件中位数是相应一维分量损失的 Bayes 决策；不能把 MSE 的精确二倍关系照搬成 MAE 结论。

### 12.3 GT 邻域命中与 Best-of-K

给定 \(y\) 与某个固定真实 \(x_*\)，设

\[
p_\epsilon(y,x_*)=
\pi(\|X-x_*\|\le\epsilon\mid y).
\]

若得到 \(K\) 个独立 posterior samples，则

\[
P(\text{至少一次命中}\mid y,x_*)
=1-(1-p_\epsilon)^K.
\tag{12.4}
\]

Best-of-K 的 oracle 风险可由这一分布积分计算。例如平方误差下，设
\(L_K=\min_k\|X_k-x_*\|^2\)，则

\[
\mathbb E[L_K\mid y,x_*]
=\int_0^\infty
[1-\pi(B_{\sqrt t}(x_*)\mid y)]^K\,dt,
\tag{12.5}
\]

在可积条件下成立。在线算法并不知道哪个候选距离 GT 最近，因此 Best-of-K 不能冒充在线 Top-1 选择性能。

对相关 MCMC 链，不能简单把 \(K\) 换成某个 ESS 就把 (12.4) 当精确命中概率。ESS 往往针对具体统计量；rare event hitting、跨 mode 和均值估计可能具有不同时间尺度。

一个较强但明确的充分条件是：独立 proposal \(q\) 满足
\(\pi(z)/q(z)\le M\)。其 MH kernel 有 minorization

\[
K(z,\cdot)\ge M^{-1}\pi(\cdot),
\]

因此任何初始状态下，连续 \(K\) 次都未进入事件 \(B\) 的概率不超过

\[
(1-\pi(B)/M)^K.
\tag{12.6}
\]

证明来自每一步至少 \(\pi(B)/M\) 的进入概率。现实中 \(M\) 可能很大、甚至无法可靠估计；这是一项条件保证，不是现有 C1 的已验证界。

### 12.4 额外测量什么时候真正提供新信息？

对旧数据 indistinguishable 的两个材料，新增测量要使其条件观测分布不同。固定高斯噪声下，新增均值差带来的 KL 为

\[
\frac12\|\Sigma_{\rm new}^{-1/2}
(F_{\rm new}(x_1)-F_{\rm new}(x_2))\|^2.
\tag{12.7}
\]

局部上，应让 \(DF_{\rm new}\) 对旧合法 tangent null 方向有非零作用。新增频率、极化、照明、接收位置或近场信息是否有效，需要具体计算，不能凭“更多数据”泛泛保证。若 gain/geometry 等 nuisance 也未知，要在相同 nuisance 模型下比较；不能允许一组 free gains 随意吸收新信息后又按已知 gain 报告辨识收益。

### 12.5 一份 GT 监督能否学到整个 posterior？

“每个连续 \(y\) 只有一份 GT”不意味着条件生成在原理上不可能。若训练样本来自正确联合分布，在足够数据和合适函数类下，条件最大似然总体目标能识别 \(p(x\mid y)\)。

但有限训练集本身并不唯一确定每个 \(y\) 附近的多解分布；需要可信的 prior、forward/noise model、跨样本共享结构或额外 posterior/physics 训练。仅用材料 MSE 会把分布压成一个 Bayes 决策，而不是学到全部概率。

另外，**posterior 与探索动力学不是同一个可识别对象**。同一个 invariant distribution 可以对应无数不同速度、不同转移、不同 nonreversible circulation 的 kernel。Maxwell 数据规定约束与 likelihood，不规定某条唯一“真实材料探索时间轨迹”。我们学习的是受计算成本约束的算法动力学，而非实验中可直接观察的物理演化定律。

---

## 13. 主算法、神经网络选项与本质不同的备选

### 13.1 最终主选：材料度量下的有限噪声传输探索

暂称 **Material-Native Finite Transport，MNFT**，只是本报告的算法标签，不是已经成立的新方法名称或原创性声明。

必要输入为：完整 forward/adjoint 接口、合法 prior generator 与其导数、原生 prior 密度、正确噪声模型、已知背景/采集信息、材料度量及有限计算预算。N3 不是必要输入。

**Anchor。** 从原生 prior 的均值/中性初始化或预先冻结的少量 TRAIN-only 初值出发，在 prior 中求 MAP 或 likelihood 受约束的物理重建：

\[
z_0\approx
\arg\min_{z\in\mathcal Z}
\left\{\Phi(g(z))-\log p_Z(z)\right\}.
\tag{13.1}
\]

可用 trust-region GN、LM、精确 Newton 或其他合适求解器，不预设 GN 最优。记录实际驻点误差、合法性和完整 residual。N3 可仅作额外 warm start；若投影回 prior，应把投影代价和改变后的 initialization 明确记录，不能把 N3 本体称为合法 prior anchor。

**局部材料几何。** 用 \(H=Dg^\top Dg\) 和去 gauge frame 形成按材料长度归一化的 \(A=JQ\)。可选择真实 tangent-null、稳定 weak cluster，或更稳健的软弱速度

\[
w=QW_\lambda(A)b,\qquad \|w\|=1.
\tag{13.2}
\]

若归一化前 \(\|w\|\) 太小，应停止或改 proposal；归一化之后式 (7.5) 的速度界必须按实际 \(b/\|QW_\lambda b\|\) 重算，不能继续使用原 \(B_0\)。工程上也可不归一化，并直接由材料距离控制时间步。

方向 \(b\) 可来自固定随机种子、任务梯度或 archive repulsion；不要让 GT 进入选择器。Repulsion 是探索策略，不是 posterior score。

**Predictor。** 对目标材料弧长 \(s\) 用 material-native retraction：

\[
z_p=\arg\min_{z'\ {\rm local}}
\frac12\|g(z')-(g(z)+s w)\|^2,
\tag{13.3}
\]

或用显式 geodesic 二阶 predictor \(z+s v-\frac12s^2\Gamma(v,v)\)。只有在需要时才加入式 (6.8) 的 fiber normal correction。

**Corrector。** 如果已有真正的正规 fiber 且存在非零 tangent，可以求

\[
\begin{cases}
F(g(z'))=F(g(z)),\\
\langle g(z')-g(z),w\rangle=s,
\end{cases}
\tag{13.4}
\]

并在这些约束下最小化 correction norm，限制在当前小 chart 内。这保留了进度坐标，避免退回原点。维数、独立约束数量和可解性必须检查；没有非零 tangent 时不要使用这个精确等式系统。

当前 P0 更适合有限噪声版本：

\[
\min_{z'}\;
\frac12\|g(z')-g(z_p)\|^2
\]

满足

\[
\Phi(g(z'))\le \Phi_{\rm cap},\quad
z'\in\mathcal Z_{\rm legal},\quad
\langle g(z')-g(z),w\rangle\ge s_{\min}.
\tag{13.5}
\]

\(\Phi_{\rm cap}\) 必须预先定义，例如由固定参考最优值加 \(\Delta\) 给出。它是探索域定义，不是自动的 posterior credible region。需要概率时另外运行针对 (11.1) 的正确 kernel。

**有限反馈与步长。** 使用精确 resolvent、完整隐式解或经验证的低阶近似评估端点。按实际材料距离而非 prior 参数范数比较。若 residual 预测误差、prior 边界、谱变化或 correction 过大，则缩步和刷新 chart。Endpoints 合格不代表路径全区间合格；要求路径证书时需额外中点/区间界，成本单独计入。

**停止。** 数据相容位移长期小于预先设定的物理分辨阈值，或成本没有超过简单 prior-coordinate 方法的价值，就停止 local-explorer 分支，不靠多画几张近重复图维持“多解”叙事。

### 13.2 概率主备选：明确密度的 prior-coordinate MH

对于 Gaussian prior 坐标，最简单基线是 pCN：

\[
z'=\sqrt{1-\beta^2}\,z+\beta\xi,\qquad \xi\sim N(0,I),
\]

接受率只需正确 likelihood ratio，因为 proposal 对 Gaussian prior 可逆。若 prior 经 CDF 映射合法化，应保持该映射完整，并保留其统计含义。

更一般可用满秩、可计算密度的局部 Gaussian：

\[
q_\theta(z'\mid z,y)
=N\!\left(z+\mu_\theta(z,y),
\ \varepsilon^2 B_\theta(z,y)B_\theta(z,y)^\top+\delta^2I\right).
\tag{13.6}
\]

它可利用 \(A,H\)、软弱滤波或学到的有限位移设置均值/协方差；MH 必须计算反向密度。固定阶段结束后冻结自适应参数，或使用有证明的 adaptive-MCMC 条件；不能认为每个瞬时 kernel 都不变，就推出任意 online adaptation 自动正确。

这种方法未必利用最高阶几何，但概率语义强、baseline 清楚、容易拆分 proposal 提升与 target 正确性。当前证据中 C1 global-MH 在四维 proxy 上已经胜过所测 local/hybrid；因此“更复杂 local geometry 一定更好”没有支持。[R8]

### 13.3 神经 Explorer：学习有限补偿，不学习 GT 图像

若 E2 证实有值得走的合法局部位移，推荐学习

\[
a_\theta=\mathcal N_\theta(z,c_y,\text{current/field summaries},v,s),
\tag{13.7}
\]

即给定当前状态、意图方向和有限步长所需的补偿。网络必须知道当前状态，不能只有初始 \(c_y\)。输出可包括 normal coefficient、步长可信度或 proposal covariance；不要求大 CNN。

Teacher 可来自无 GT 的高精度 material-native continuation，训练目标例如

\[
\mathcal L_{\rm trans}
=
\|a_\theta-a_{\rm corr}\|_{H_N}^2
+\lambda_F\|F(g(z_\theta'))-F(g(z))\|^2
+\lambda_p
\big(\langle g(z_\theta')-g(z),w\rangle-s\big)^2.
\tag{13.8}
\]

有限噪声版本应把第二项改为明确的 likelihood-domain penalty，而不是所有样本都强迫 exact-fiber。可另加 JVP/Sobolev 监督，以控制对方向和步长的导数；这与 DINO 的 derivative-informed 思想有共同点，不是独占贡献。[L4]

进度项和真实非零 teacher 位移防止“输出 no-op 即物理损失为零”的塌缩。Prior legality 由坐标保证，不能仅靠软 penalty。若 exact teacher chart 尚未建立，就不能声称正在监督学习真实 fiber map。

部署时网络只提供初猜；完整 corrector 验证。如果将它用于 posterior，应把它放入有明确密度的 proposal 或具有可验证逆与体积的 augmented transport 中。直接把“预测+投影”接到普通 MH 上不充分。

**NN 必要性的判据是摊销收益：**

\[
N_{\rm break-even}
>
\frac{C_{\rm training}+C_{\rm teacher}}
{C_{\rm baseline/sample}-C_{\rm learned/sample}},
\tag{13.9}
\]

分母必须为正，且两种方法达到同一材料覆盖/概率误差标准。以缓存 CNN latency 代替总成本是不成立的。

### 13.4 全局生成扩展：负责跨区域提案，不必取代所有局部物理

只有确认存在其他高支持区域而 local kernel 难以抵达时，才训练或复用

\[
q_\theta(z\mid y).
\]

推荐先用明确密度的低维 MDN/flow 作为独立 proposal，防御混合

\[
q_{\rm mix}=(1-\eta)q_\theta+\eta p_Z,\qquad \eta>0,
\tag{13.10}
\]

并通过准确 MH 纠偏。不同 kernel 用固定混合概率组合即可保持共同 posterior。若混合概率依赖当前状态，则需要正确处理整个转移概率；不能套用固定混合的结论。

全局 generator 应输出合法材料坐标或非局部 initialization/mode jump，而不是必须先把 hole count 分类。对于有限噪声、有低密度桥的连续 support，普通连续 flow 可以表达多个 mode；若 prior support 真正不连通，而 latent support 连通且 generator 连续，则不能在不穿过非法区域的同时覆盖两个分支，需要离散潜变量、多 chart 或其他非连续/分片结构。这是 support 问题，不是“所有多峰都需要拓扑分类”。

生成器的 physics context、k-space、current 或 feedback features 不会增加 \(I(X;Y)\)，但可能改善计算条件数、表达效率与迁移。必须以相同观测、容量、训练预算消融，不能把同一个 \(y\) 的多种确定性表示当成多次独立 likelihood。

### 13.5 电流域备选

本质不同的备选是第 8 节的 receiver-null current continuation：在 \(Sk=0\)、多源 minors、prior 与可逆性约束上做有限探索。它不依赖 N3，也不依赖一个固定材料 tangent chart。

该路线只有在 current-space 压缩确实降低成本且不破坏材料 realizability 时才值得推进。否则，原生七维 prior 上的简单物理 MH 或 constrained optimization 更合理。这里不建议把两条路线一开始拼成一个庞大系统。

### 13.6 计算账：有限 resolvent 不是免费算子

设 current 维数为 \(n_j\)，源数为 \(P\)，有效 prior 维数为 \(d\)。Dense 实现中，一个新状态通常需要 \(O(n_j^3)\) 的分解及 \(O(Pn_j^2)\) 的源求解；已有分解上的一个方向 JVP 约为 \(P\) 个 RHS，完整 \(d\) 方向 tangent 约为 \(dP\) 个 RHS。给定一方向的二阶和三阶递推各再增加约 \(P\) 个 RHS；混合 Hessian 可复用两个一阶响应后再解合并 forcing。实际数目必须按后端调用记录，而不只按理论式计数。

完整输出的 VJP/接收伴随可以按源、任务或接收维数组织，成本不能固定假设为一次；只估计一个任务标量与保存整个 \(m\times n\) Jacobian 是不同工作。对三维迭代后端，应另记每 RHS 的迭代数、Green matvec、预条件构建与内存，而不是把 dense 复杂度直接套用。

特别注意：

\[
\boxed{\text{低维材料 latent}\ \not\Rightarrow\
\text{低秩 }\Delta XG.}
\]

一个只有单参数的全域 contrast 缩放，也可以改变满秩对角 \(X\)。只有局部少单元改变、可验证低秩结构或准确压缩时，Woodbury/Schur 才可能把有限更新变成小系统。把 \(d=7\) 当作 rank-7 Maxwell operator update 是错误的。

一个网络步的完整成本包括当前材料的 physics context、proposal/decoder、corrector 的所有失败尝试、J/metric 刷新和最终 likelihood。若用式 (13.8) 的真实 physics loss 训练，还必须计入训练中对新预测材料的 F/J 调用。较省成本的版本只用已保存 teacher 的位移与导数监督，把真实 F 留给预先冻结的验证样本；这改变了训练约束强度，应单独标注。


---

## 14. 现有失败意味着什么？如何设计 P1 而不偷换目标？

原 P0 的结果没有给出“全局材料已经唯一可辨识”的证明，也没有给出“只是搜索器太弱”的证明。当前证据支持几项并存的解释。

**Anchor 问题有直接证据。** N3 和一次 physics-R base 没有满足原 P0 的完整资格，不能用于证明该先验内的局部解族。

**坐标与有限步问题有直接证据。** 合法 TRAIN chart 的 finite decoding 明显不同于初始切向；但本次推导表明 raw decoder 曲率还需拆成 intrinsic 与 coordinate 部分。

**原 prior 内局部可辨识有数值支持，但只是局部。** 七维有效 tangent 满列秩意味着这里可能只有小噪声宽度。Prior-Constrained 未找到竞争 pair，是方法/预算下的负结果，不是全局唯一性定理。

**高阶物理补偿是真机制，但未证明足够大且合法。** 精确 resolvent 能准确计算一条 off-prior oracle chord，并不意味着该 chord 可以被合法化或能通向新的高 posterior mass。

**数据不足以支撑复杂多模式学习，是可信限制而非唯一归因。** 64 个基础 TRAIN parents、狭窄生成规则与少量公开校准点，不足以独立验证广义复杂二维 material posterior 的学习性能。[R1,R5,R8]

### 14.1 合理的 P1 应从材料物理定义

一种候选是固定频率/已知频散下的被动混合材料场：

\[
\chi(r,\omega)=\chi_b(r,\omega)+
\sum_{a=1}^{A} f_a(r)c_a(\omega),
\qquad f_a\ge0,\quad\sum_a f_a\le1.
\tag{14.1}
\]

其中 \(c_a(\omega)\) 为预先定义的被动材料响应，\(f_a\) 由带空间相关长度与制造/组织尺度限制的随机场产生。可以用多尺度 Gaussian field、稀疏小波或 level-set/soft mixture 实现，而不是预定义孔洞数量。

它的物理意义是局部混合或占据率不确定性；若研究对象实际是锐利二元介质，应另外定义边界模型，不能把 arbitrary soft fractions 当成同一二元材料 prior。

P1 应同时允许局部连续变化和区域形态改变，并具有明确 prior density 或可精确采样的生成机制。相关长度、contrast range、材料类别和 noise/acquisition 必须在看测试观测之前冻结。不能事后放宽损耗、增加 free gains 或独立实虚像素来制造“成功”。

### 14.2 一个合格 ambiguity benchmark 需要什么？

先在冻结 prior/采集/噪声下寻找可验证的竞争区域，并估计其质量或至少给可信局部下界；它们必须材料上有意义，而非仅为极微小像素噪声。再冻结为 discovery set，另设未用于挑选的 evaluation set。

小噪声下既可能出现连续 fiber，也可能是多个分离的局部可辨识解；两类都值得研究，但应分别报告。若 P1 的 posterior 仍基本单峰，就把问题转向可信不确定性宽度或 active sensing，而不是继续强迫生成多种结构。

---

## 15. 五篇最接近的核心先例与原创性边界

下述比较只使用取得的原论文材料；未取得证明全文的条目明确标注。比较结论中“可迁移到 Maxwell”指数学接口可容纳非线性 \(F\)，不等于论文已经做了 Maxwell 实验。

### 15.1 Montanaro–Valsesia–Magli：e-GLASS，ICIP 2022

它从 GAN latent 的一个反演点出发，用测量敏感矩阵和感知差异指标挑选有限扰动方向。forward measurement 是线性的，但 decoder 非线性；探索的是候选差异，不是校准 posterior。它没有用完整 Maxwell finite response 解决有限步补偿。**本任务与它最直接重合的是“生成先验内的弱方向探索”；差别应体现在 intrinsic prior curvature、有限补偿与概率纠偏，而不是声称首次探索 latent nullspace。**[L1]

### 15.2 Boink–Haltmeier–Holman–Schwab：非线性数据一致网络，IPI 2023

这是 classical reconstruction 后接 data-consistent learned map 的理论框架，关注噪声趋零的正则化与稳定性；它不是 posterior sampler。其非线性 forward 框架比线性 null-space network 更接近本题。“学习映射后再投影回相同数据”已经有先例。本次未取得完整证明正文，因此不把未核对的条件当作已排除的理论差异。[L2]

### 15.3 Au–Graham–Thiery：Manifold Lifting，JRSSB 2023

它把参数和噪声一起提升到约束流形，使用 coarea volume term 与 constrained HMC；目标是正确 posterior，而不是增加图像数量。非线性几何由约束积分、投影和可逆性检查处理。**“在弯曲可行流形上正确采样”已有成熟先例；Maxwell 的增量必须来自具体算子计算、材料 prior 与完整成本。**[L3]

### 15.4 Cao–O’Leary-Roseberry–Ghattas：DINO 几何 MCMC，JMLR 2025

它用 forward 值与参数 Jacobian 联合训练 neural operator，再以 surrogate 几何提案和 delayed acceptance 加速昂贵 PDE posterior sampling。原实验是 diffusion–reaction 与 hyperelastic material，并非本任务 Maxwell。它直接占据“导数学习 + 后验几何 + 正确 MH + 摊销成本”这条主线。新的工作要比较有限 Maxwell transport 是否带来同成本增量，而不是仅换一个名字重新组合这些组件。[L4]

### 15.5 Gabrié–Rotskoff–Vanden-Eijnden：NF 增强自适应 Monte Carlo，PNAS 2022

它结合局部 MCMC 与 flow 非局部 proposal，并利用链样本更新生成模型。其理论与实验明确区分 proposal 分布和目标分布，且研究的跨 basin 设置包含粗略已知 mode 位置。它支持 global/local 组合，但不替我们证明未知复杂 Maxwell mode 可以自动发现。本次核对到原论文摘要与公开方法说明，未完整复核全部理论证明。[L5]

### 15.6 其他不可忽略的边界

Deep Null Space Learning 对线性算子提供结构化数据一致性与正则化分析；它学习一个解选择/修正，而不是所有解的概率。[L11] DILI 已利用 prior-preconditioned likelihood geometry 做 MCMC，不能把“在 likelihood-informed 子空间单独移动”当作新概念。[L12]

电磁 Born/inverse Born 收敛、Padé 多重散射以及 CSI 的双变量结构均有先例。[L7,L8,L10] 2026 年 inverse-Born neural operator 新文又进一步压缩了“把高阶级数做成 NN”这一表层创新空间。[L9]

因此，这个项目较可信的贡献候选不是“Dyson + NN + MH”，而是：

\[
\boxed{
\text{先验商几何}
+\text{Maxwell 有限反馈}
+\text{多源材料兼容}
+\text{受控有限传输}
+\text{可验证的总成本收益}.
}
\]

---

## 16. 最小决定性实验：不再重复支付已知负结果

### 16.1 阶段安排

| 阶段 | 核心区分问题 | 当前执行范围与通过条件 |
|---|---|---|
| E0 几何审计 | 弱 GT chord 是否真是合法局部自由度？ | 先只读旧 J/预测/材料；把 oracle chord、prior tangent、anchor 合法性分开 |
| E1 intrinsic audit | 大 decoder 曲率是坐标效应还是 prior 真弯曲？ | 用旧 Dg、J、decoder 导数与已有二阶响应做三项分解；不得偷偷重新解旧物理 |
| E2 局部探索 | 材料重traction是否带来更大同支持位移？ | 新、独立合法 physics anchor；同材料距离、同 full 成本比较 |
| E3 学习传输 | 是否值得摊销一个非平凡有限补偿 map？ | 只有 E2 有实质合法范围和可学习重复结构才训练小网络 |
| E4 概率 | 正确 kernel 是否混合并估计对事件质量？ | 需要独立可信 reference；正确 invariant law 与有限链误差分开 |
| E5 全局 | 是否存在 local 无法访问的其他高支持区？ | 先有竞争区域证据，再做全局提案；不预设 topology 类 |

### 16.2 最应立刻做的零新增 F/J 审计

对两份已有合法 TRAIN chart，使用保存的 material tangent 与 Jacobian：

\[
d_2=D^2g[v,v],\quad
d_{2,\parallel}=P_{\mathcal P}d_2,\quad
d_{2,\perp}=(I-P_{\mathcal P})d_2.
\]

逐项计算

\[
C_M=D^2F[w,w],\quad
C_{\rm coord}=Jd_{2,\parallel},\quad
C_{\rm prior}=Jd_{2,\perp},
\]

并检查其向量和、方向、余核投影，而不是只比较范数比。

使用 SVD 去 gauge 后构造材料正交 frame；参数先验不因此自动降成均匀七维。对原 decoder 的不同差分步长做收敛检查，并报告 active max cell、CDF 饱和、H 条件数与数值 rank。若缓存缺少 \(C_M\)，标注缺失；不要悄悄调用 \(R_0\) 再声称“零物理”。已有显式 J 的矩阵乘法与 decoder-only 求导可以做，新的 current linear solve 必须计 RHS。

这个实验直接区分：是否值得先换 chart，而不是换 Maxwell 高阶模型。如果不可消除的余核项很小而 coordinate 项很大，材料原生重traction应成为首选；若余核项很大，则局部步长上限是真障碍，训练 tangent NN 不会使它消失。

### 16.3 独立 anchor 与 E2 的新预算提案

这是一份给下一轮 Codex 的**新预算建议，不是本次已经执行的计算**：

- 总上限：64 个新完整材料 states、2048 个 complex RHS、GPU 训练 0、非径向 3D 0；所有失败与拒绝计入。
- 先从两个预先固定的 exposed smooth development 条件验证条件 smooth-P0 内的 physics-only anchor。该子问题的 prior 条件必须公开，不称一般未知结构反演。
- Anchor 最多占用 24 states；其余留给同材料距离的两方向/小步 continuation 对照。若无法形成平衡比较，保留 NOT_RUN，不从剩余额度临时挑赢家。
- 背景初始化、prior mean 或 TRAIN-only 初值可以使用；GT 只进入独立离线评估器。N3 warm start 作为另列基线，不污染“physics-only”定义。
- 比较 frozen tangent、原参数 curvature predictor、material-native retraction + refreshed corrector；高阶物理阶数的旧 E1 不重复跑。
- 阶段总预算先冻结，各臂同时报告自己的必要 setup/refresh 成本与可共享缓存成本，不能只给新方法享受共享优惠。

非平凡位移阈值可预先设为 TRAIN 材料增量 RMS 中位数的 5%，再加上与实际研究任务一致的多尺度结构差异标准；这是建议的工程 gate，不是理论常数。若研究目标是局部数值不确定性而非结构 mode，位移与事件尺度应相应重定义，并公开说明。

### 16.4 E3–E5 的门槛

E3 只在多个独立 anchor 的 material-native 轨迹超过预注册有意义尺度、prior 合法、完整 likelihood 合格后启动。训练的 teacher 数据成本必须计入总账。建议第一版不超过约 5 万参数；不是因为这个数有理论最优性，而是避免用容量掩盖不存在的有限传输。

E4 不能用 64-state 预算假装建立一般七维 posterior 的高精度 reference。可先用解析代数问题验证 kernel，再申请独立 reference 预算；二维条件切片也可用，但必须明确它不是完整七维 posterior。既有四维 C1/MH chains 可做只读概率一致性检查，不重复付费制造“新增结果”。

E5 必须有至少两个材料上分离且概率质量有依据的高支持区域。多张 event 代表图不等于分离 mode；同一连续 ridge 上不同事件也可以有显著质量。

### 16.5 评价指标必须分开

必须分别报告：N3 Top-1、physics-only Top-1、oracle Best-of-K、局部材料覆盖、posterior event/CDF 误差、完整 Maxwell residual 与总成本。

新增诊断至少包括：prior 合法率、独立材料去重数量、相同材料距离上的数据漂移、路径中间点检查覆盖率、gauge/active-boundary/rank 状态、每有效位移成本、每有效事件估计成本。不要用 prior 参数步长相同来声称材料距离公平。

N3/physics 的误差使用所有预定案例；不能只在通过 noise 球的少数例子上报告效果。已暴露 development 不称独立 blind test。新的 P1 不与原 P0 成功率直接合并。

---

## 17. 论文边界、No-Go 与十个核心问题的最终答案

### 17.1 值得成为理论论文的部分

一种清楚的理论主线是：

**Maxwell 有限材料 fiber 的 prior-aware 几何与阻碍。**

核心不是重命名常秩定理，而是建立并验证式 (6.6)–(6.8) 的物理/先验/坐标分离、多源电流兼容的有限代数接口，以及能判别有限补偿可达范围的可计算界。还需要真实 vector-Maxwell 或充分明确的电磁适用范围与物理 witness。

若只得到几个普通矩阵恒等式、两个 scalar-TM 例子和未量化的速度优势，尚不足以认定强 IEEE TAP 成果。理论正确与投稿竞争力是两个判断。

另一条可独立成文的方向是：

**Learning finite material transport for posterior exploration。**

它需要多个有可信 posterior/ambiguity 的复杂材料问题、训练/测试分离、与 DINO、普通 geometric MCMC、global density + MH 的强对照、完整摊销成本和 calibration。它更偏生成式计算成像，不必把 A17 放成所有贡献的起点。

原 A17 仍保留其 current-representation/material-update 的定位；不需要用未完成的新探索系统来替它补一个并不存在的主定理缺口。

### 17.2 明确 No-Go 条件

若去 gauge 后 prior 内统一满列秩、有限噪声宽度低于任务尺度，停止 exact-local-fiber 学习。若 intrinsic audit 证明主要是坐标问题，先换 chart/retraction，不训练“补高阶 Maxwell”的网络。

若完整有限 feedback 虽准确却没有带来更大合法位移或更低成本，停止将 Dyson/Padé/Schur 当作算法卖点。若 prior 放宽才出现解族，作为 P1 新目标报告，不宣称解决了 P0。

若候选只是在极低 prior 密度区制造相容性，或者没有可信事件质量，不能称合理 posterior worlds。若 proposal density/投影 Jacobian 不可得而又没有另一个不变性推导，不报告精确 posterior sampling。

若 learned explorer 没有越过完整 teacher+训练成本的 break-even，保留 physics-only。若没有竞争高支持区域，暂不训练大规模 global generator；更可靠的替代是单点估计加局部不确定性、或者主动测量。

### 17.3 Q1–Q10

**Q1：何时剩余误差对应 local fiber？**  
两端必须属于同一合法正规 prior 区域，满足相同数据，且 prior-restricted Jacobian 在邻域恒秩并有非零 kernel。有限强分量通常是曲率补偿，不必为零。仅有低 MAE、低 residual 或大 weak energy 不充分。

**Q2：GT 可能在其中吗？**  
可能，但没有一般保证。同一 IFT graph 邻域给局部可达；noise tube、相同精确 fiber、分离全局解是不同对象。两单元反例表明材料任意接近也可处于不同局部解点。

**Q3：Frozen Jacobian 与完整 feedback 差在哪里？**  
前者只给一阶响应；后者准确给有限 secant、曲率、高阶相互作用和失效尺度。可通过复用 current resolvent 的方向递推高效计算，但是否值得由实际成本决定。

**Q4：怎样构造运动？**  
先在真实 prior 的材料度量中构造方向，使用 intrinsic/geodesic retraction、有限反馈 corrector 与显式进度约束，动态刷新。没有精确 fiber 时改走有限 likelihood 域，不伪造 nullspace。

**Q5：如何正确给概率？**  
固定 prior/noise/target，用明确密度 MH/pCN 或正确 constrained kernel；小噪声 fiber 需要 coarea/volume 因子。局部条件概率与全局质量分开，图像代表与事件概率分开。

**Q6：Explorer 相对 N3 解决什么新问题？**  
它解决分布覆盖、合理备选、事件不确定性与后续信息采集，不创造新的观测信息，也不保证降低 Top-1 MSE。N3 应作为单图风险基线保留。

**Q7：高曲率会让学习无用吗？**  
不必。需先除坐标伪曲率，然后检查 uniform 曲率、谱/软滤波稳定性、prior 边界和模型误差。满足这些条件时有有限时间稳定界；不满足时需要小步、atlas、global jump，甚至停止学习。

**Q8：A17/Dyson/Schur 的独特增量在哪？**  
候选增量在 Maxwell-specific 多源材料 realizability、有限 primal/dual feedback error 与复用计算，而非仅把已有 Jacobian/MCMC 重新命名。本次没有证明同成本优势。

**Q9：不依赖 N3 选什么 anchor？**  
首选 prior 内的独立 full-physics MAP/likelihood reconstruction，求解器按问题选。N3 在已验证能降低总初始化成本时保留 warm start；它本身不自动合法。

**Q10：到底学习 manifold、field、kernel 还是 posterior？**  
几何由 forward+prior 决定，能算就不必整张学；有限补偿 map 可用于快速移动；kernel 再加目标密度与不变性；全局 generator 用于跨区域。它们不能无条件互换。最值得优先学习的是有可验证物理 corrector 的有限传输 proposal，而不是声称学到唯一“真实材料动力学”。

---

## 18. 本次轻量数学检查与结果边界

随附 `checks/algebra_checks.py` 固定随机种子 20261016，仅调用 NumPy/SciPy 在小矩阵上检查。结果包括：

| 检查 | 数值结果 |
|---|---:|
| 有限 resolvent 与直接差分状态解相对误差 | \(2.04\times10^{-15}\) |
| 多源兼容式恢复已知材料增量相对误差 | \(2.50\times10^{-16}\) |
| 路径有序场 ODE 与端点 \(A(1)^{-1}A(0)\) 相对误差 | \(1.63\times10^{-13}\) |
| 双侧有限输出残差恒等式相对误差 | \(1.34\times10^{-14}\) |
| 非正规 \(K\)：范数 / 谱半径 | \(5.013 / 0.3\) |
| 上例四十项 Neumann 相对误差 | \(3.94\times10^{-16}\) |
| 两单元 \(q\) 与 \(-q\) 的数据差 | 数值 0 |
| 平坦 prior 的 raw decoder 二阶项 | 2000，而 intrinsic 项为 0 |

有限差分一、二、三阶相对误差分别约
\(2.30\times10^{-10},7.97\times10^{-8},1.01\times10^{-3}\)；
三阶差分含明显减法消去误差，不能拿它宣称第三阶导数已达到机器精度。精确递推公式的推导与该有限精度核验分别保留。

这些检查支持代数实现一致性，不证明 continuum-Maxwell 离散误差、复杂材料 posterior 正确性、神经网络泛化或论文实验收益。

---

## 参考与访问标识

以下 URL 用于复现定位；本报告所有实验数字按其原目标与证据层级引用。完整访问状态另见 `sources/SOURCE_AUDIT.md`。

### 项目原始材料

**[R1]** `migodam/b1-probabilistic-maxwell-imaging`，commit  
`c76733676a7746627dfbf710c5f31b6375ae4335`，  
`research/finite_feedback_fiber_diagnostic_20261016/` 下的
`START_HERE.md`、`DIAGNOSTIC_REPORT_ZH.md`、`METRICS.json`（核心前段）、
`REPRODUCE.md`、`THEORY_QUESTIONS_FOR_PRO.md`。  
https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/c76733676a7746627dfbf710c5f31b6375ae4335/research/finite_feedback_fiber_diagnostic_20261016

**[R2]** 同上，`LOCAL_PATH_AND_CURVATURE.md`。  
https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/c76733676a7746627dfbf710c5f31b6375ae4335/research/finite_feedback_fiber_diagnostic_20261016/LOCAL_PATH_AND_CURVATURE.md

**[R3]** `migodam/a26-om-feedback-prior-completion`，`src/a26/decoder.py`；本次关联源码快照 main tree 为 `2603f3bd16c48e48587fedc4363975779def17f5`。  
https://github.com/migodam/a26-om-feedback-prior-completion/blob/2603f3bd16c48e48587fedc4363975779def17f5/src/a26/decoder.py

**[R4]** 同仓库同 commit，`src/a26/training.py`，本次直接读取含 310–558 行的归一化、SSIM 与 delta_loss。  
https://github.com/migodam/a26-om-feedback-prior-completion/blob/2603f3bd16c48e48587fedc4363975779def17f5/src/a26/training.py

**[R5]** 同仓库，`A26_RESEARCH_REPORT_ZH.md`、`configs/a26/pilot.json`、`REPRODUCING_A26.md`。  
https://github.com/migodam/a26-om-feedback-prior-completion

**[R6]** `migodam/a17-maxwell-certified-exchange`，读取时默认分支的
`closure_final/START_HERE.md`、`closure_final/manuscript/CORE_THEORY.tex`、
`PAPER_TAP_FINAL.tex`。理论 blob：`c31c5c6790820b153d2b6a738df94fda0c12165d`。  
https://github.com/migodam/a17-maxwell-certified-exchange/blob/main/closure_final/manuscript/CORE_THEORY.tex

**[R7]** `migodam/a30-maxwell-observable-material-manifold`，
commit `717d3d26a290248bd2b9aeb55f3d03694b407b84`，`A30_E0_E1_DECISION.md`。  
https://github.com/migodam/a30-maxwell-observable-material-manifold/blob/717d3d26a290248bd2b9aeb55f3d03694b407b84/A30_E0_E1_DECISION.md

**[R8]** LXZ，commit `dde0237eaa0d8c5507fedb07bb704c66a817028d`，  
`research/research_lxz_20261016/RESEARCH_REPORT_ZH.md`。  
https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/dde0237eaa0d8c5507fedb07bb704c66a817028d/research/research_lxz_20261016/RESEARCH_REPORT_ZH.md

**[R9]** C1 与 B1 入口及 Joint P2 状态。  
https://github.com/migodam/c1-physics-state-generative-imaging/blob/main/START_HERE.md  
https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/START_HERE.md

### 外部原论文

**[L1]** A. Montanaro, D. Valsesia, E. Magli. *Exploring the solution space of linear inverse problems with GAN latent geometry*. ICIP 2022. arXiv:2207.00460.  
https://arxiv.org/abs/2207.00460

**[L2]** Y. E. Boink, M. Haltmeier, S. Holman, J. Schwab. *Data-consistent neural networks for solving nonlinear inverse problems*. Inverse Problems and Imaging 17(1), 203–229, 2023; early online 2022. DOI:10.3934/ipi.2022037.  
https://www.aimsciences.org/article/doi/10.3934/ipi.2022037

**[L3]** K. X. Au, M. M. Graham, A. H. Thiery. *Manifold lifting: scaling Markov chain Monte Carlo to the vanishing noise regime*. JRSS Series B 85(3), 757–782, 2023. DOI:10.1093/jrsssb/qkad023.  
https://academic.oup.com/jrsssb/article-pdf/85/3/757/50859790/qkad023.pdf

**[L4]** L. Cao, T. O’Leary-Roseberry, O. Ghattas. *Derivative-Informed Neural Operator Acceleration of Geometric MCMC for Infinite-Dimensional Bayesian Inverse Problems*. JMLR 26, 2025.  
https://www.jmlr.org/papers/volume26/24-0745/24-0745.pdf

**[L5]** M. Gabrié, G. M. Rotskoff, E. Vanden-Eijnden. *Adaptive Monte Carlo augmented with normalizing flows*. PNAS 119(10), e2109420119, 2022. DOI:10.1073/pnas.2109420119.  
https://arxiv.org/abs/2105.12603

**[L6]** S. Blanes, F. Casas, J. A. Oteo, J. Ros. *The Magnus expansion and some of its applications*. Physics Reports 470, 151–238, 2009. arXiv:0810.5488.  
https://arxiv.org/abs/0810.5488

**[L7]** K. Kilgore, S. Moskow, J. C. Schotland. *Convergence of the Born and inverse Born series for electromagnetic scattering*. Applicable Analysis 96(10), 1737–1748, 2017. DOI:10.1080/00036811.2017.1292349.  
https://researchdiscovery.drexel.edu/esploro/outputs/journalArticle/Convergence-of-the-Born-and-inverse/991019167764504721

**[L8]** *Padé approximants of the Born series of electromagnetic scattering by a diffraction grating*. Physical Review A 109, 033522, 2024. DOI:10.1103/PhysRevA.109.033522. 本次核对出版摘要，不据此声称全文逐式核验。  
https://journals.aps.org/pra/abstract/10.1103/PhysRevA.109.033522

**[L9]** J. C. Schotland, A. Titi, J.-N. Wang. *Inverse Born series based neural operators*. arXiv:2608.18262, submitted 18 August 2026. 本次摘要级核对，未取得全文。  
https://arxiv.org/abs/2608.18262

**[L10]** P. M. van den Berg, R. E. Kleinman. *A contrast source inversion method*. Inverse Problems 13, 1607–1620, 1997. DOI:10.1088/0266-5611/13/6/013.  
https://doi.org/10.1088/0266-5611/13/6/013

**[L11]** J. Schwab, S. Antholzer, M. Haltmeier. *Deep Null Space Learning for Inverse Problems: Convergence Analysis and Rates*. Inverse Problems 35, 025008, 2019. arXiv:1806.06137.  
https://arxiv.org/abs/1806.06137

**[L12]** T. Cui, K. J. H. Law, Y. M. Marzouk. *Dimension-independent likelihood-informed MCMC*. Journal of Computational Physics 304, 109–137, 2016. arXiv:1411.3688.  
https://arxiv.org/abs/1411.3688

---

**最终建议：** 下一轮最有信息价值的动作不是训练更大的网络，而是利用已有缓存，计算“可消除的坐标弯曲”与“不可消除的 prior/Maxwell 二阶障碍”各占多少。随后以独立合法物理重建为起点，把有限噪声材料传输做成一个能失败、能停、能计费、能验证概率边界的算法。若这一步成功，再学它；若失败，就把结论转化为可信的不确定性与选测问题，而不是继续承诺从相同数据中找出唯一真值。
