# FBG-R3-T：Agent A+F 研究备忘录

日期：2026-10-31。范围：成熟 CIE-I / CSI / FFT-TSOM、最近邻、固定电流材料读出的缺陷几何，以及受控 APCS 释放。本文为整合用中间稿，不是已经完成的成像性能报告。

## 1. 核心判断

最值得保留的新连接，不是另选一个标量分式，而是区分**物理材料切向电流**与**自由辅助电流**。两者在当前 CIE 叙事中常被混用，但它们离开材料—电流可实现流形的阶次不同。

对固定候选电流、固定正权重、实材料的理想逐点材料 LS，CIE 读出与普通 CSI 读出的差严格等于一个带符号明确的可实现性缺陷平方项。在正确分支上，该差非负。围绕物理解，任意线性电流路径一般产生二阶差；沿 Maxwell 物理切向直线，则首次可能出现四阶差；沿精确有限物理路径，差恒为零。这是本次推导、并经 root 与 A17 agent 独立交叉核验的数学结果，尚不声称文献中的首次发现。

它有三个直接后果。

1. “modified contrast 的导数小”不能单独证明材料读出更稳定：在物理点，逆变换恰抵消这一导数缩小，材料读出的一阶导数与 β 无关。
2. CIE 的有意义作用位于离物理流形的联合搜索、状态残差度量、电流自由度约束与分阶段路径；不能归结成材料坐标变化。
3. 新方法可以把 A17 的物理切向能力与辅助电流的法向自由度分开设计，后者通过 profile curvature 的损失量控制释放。但该控制是有条件的局部保证，不是全局凸化或 basin escape 定理。

成熟 CIE 已含 DPCS/APCS、FFT-TSOM、小到大的 Fourier 子空间、β continuation 和两 β 更新。2026 LUCIE 已直接研究“CIE 电流恢复后线性恢复材料”。2021 VE-NIE-CSI 已结合 CIE/NIE 与虚拟实验、电流正则。DCTMC、WRI 和 variable projection 覆盖了多种 lifting / 消元想法。因此，这些组合或公式本身不能作为新颖性中心。

## 2. 统一记号与材料可实现流形

采用离散标量形式；3D 各向同性介质可将每个体素的三个场分量和全部源堆成一行，以下逐行结论照用：

\[
J\in\mathbb C^{n\times p},\quad E\in\mathbb C^{n\times p},\quad
Q(J)=E+GJ,\quad X=\operatorname{diag}(m),\quad Y=SJ.
\]

物理状态约束、当前 Maxwell 算子为

\[
e(m,J)=J-XQ(J)=0,\qquad L(m)=I-XG.
\]

实材料可实现流形为

\[
\mathcal M=\{(m,J):m\in\mathbb R^n,\ J=X(E+GJ)\}.
\]

当 L(m) 可逆时，J(m)=L(m)^{-1}XE。可实现性不只是一个正则先验：它要求所有照明源在每一体素共享同一实材料值。任意数据相容电流一般不满足这个共享关系。

以下内积包含任何**预先固定的正源权重**。若材料实，必须用 \((u,v)_\mathbb R=\Re(u^*v)\)，把复源行视作实 Hilbert 空间。投影 \(\Pi_q^{\perp,\mathbb R}\) 表示去掉 q 的**实线性张成**，不是去掉复线性张成。两者不可混用；单源时，后者会错误地消除相位不一致。

## 3. 成熟 CIE 的完整 formulation 与机制分解

### 3.1 状态改写与联合度量

令 β>0，\(R_\beta=\beta X(I+\beta X)^{-1}\)。统一到上述记号，CIE 状态残差为

\[
r_\beta(m,J)=\beta J-R_\beta\{E+(G+\beta I)J\}.
\]

直接展开得到精确恒等式

\[
\boxed{r_\beta(m,J)=D_\beta(m)e(m,J),\qquad
D_\beta(m)=\beta(I+\beta X)^{-1}.}
\]

因此，物理零点集合不变，而联合材料—电流空间上的残差度量变化。若每一步已精确消去 J，仅留下相同的数据误差，那么 rβ 恒为零；此时 m 与 R 的可逆坐标替换本身不消除局部极小。反之，CSI/CIE 允许 J 暂时偏离物理状态，改变该残差度量会改变联合优化路径。

对 m≥0，Rβ∈[0,1)，但这不等于整个场迭代算子收缩。令 w=βJ，状态方程是

\[
w=R_\beta E+R_\beta(I+G/\beta)w.
\]

大 β 减小显式 G/β 项，却保留趋近单位的局部回路 Rβw。一维代数例子中，m,g>0 的迭代乘子为 \(m(\beta+g)/(1+\beta m)\)，只要 gm>1 就对所有 β 大于 1。这个例子否定“有界 modified contrast 自动蕴含统一 Banach contraction”，不声称任意这样的 g 都对应连续 Maxwell 自作用核。

还应注意：在 2019 的未整体重标定状态目标中，β 趋零使 Dβ≈βI，等于减小状态残差相对数据残差的权重；不可把 β=0 的整个目标直接称为普通 CSI。后文的 \(m_{\rm CSI}\) 是普通材料 LS 读出，亦是 β 正向趋零的读出极限。

### 3.2 六元组 A0：成熟单 β CIE-I/FFT-TSOM

令接收算子 S 的截断 SVD 为 \(U_L\Sigma_LV_L^*\)。定义

\[
J_D=V_L\Sigma_L^{-1}U_L^*Y,\qquad
U_A=(I-V_LV_L^*)F_{\rm low},\qquad J=J_D+U_AA.
\]

F_low 是选中的低频 Fourier 列；U_A 一般并不正交，可数值正交化。DPCS 固定已经稳定测得的电流分量，APCS 位于其正交补，包含弱可见、不可见但可能影响内部多次散射的电流。第二重 Fourier 限制才是稳定这些自由度的重要正则。

| 六元组分量 | 具体内容 |
|---|---|
| Optimized Variables | 修改后的材料 R 与每个源的 APCS Fourier 系数 A |
| Physical Constraints | 多源共享材料；物理状态通过 CIE 残差弱实施；DPCS 固定；并非每个迭代都解精确 Maxwell |
| Objective | \(F_\beta(R,A)=\sum_s\{\|Sj_s-y_s\|^2/\|y_s\|^2+\|\beta j_s-R(q_s+\beta j_s)\|^2/\|E_s\|^2\}\) |
| Update | 每次一轮 APCS Polak–Ribière CG，随后材料 R 的逐点 LS；不是每块全局精确求解 |
| Subspace | 外部 SVD 的 DPCS/APCS 分解，加内部低频 Fourier 截断；随阶段放宽 |
| Globalization | 先大 β、小子空间；后小 β、大子空间，并沿用上一阶段结果作为初始值 |

对应的复材料 LS，固定 j、q 后逐体素为

\[
R_i=\frac{\sum_s w_s\overline{(q_{is}+\beta j_{is})}\,\beta j_{is}}
{\sum_s w_s|q_{is}+\beta j_{is}|^2}.
\]

2016 Eq.14 是复 R 的式子；本文第 4 节实 R 定理是其状态方程所诱导的实材料理想块，不能冒称完全复刻论文的全算法。2019 的数值方法同样包含复介电率例子。2016 原文已允许 β 随位置变化，故空间自适应 β 的概念也非新颖性依据。

### 3.3 两 β 版本的严格限度

2019 在后续阶段允许用 β_J 更新电流、用 β_R 更新材料，且 β_J>β_R。两次更新使用不同的状态残差度量，不能直接套用“同一个固定目标的 block coordinate descent”下降证明。相同物理零点集合不推出共同 Lyapunov 函数。

具体复现须写清：在两块之间先恢复同一物理 m，再将其映射到对应 β 的 R；不能把同一个数值 R 在不同 β 下解释成同一材料。即使如此，两目标也仍不同。本次未从文献取得一个覆盖该两 β 方案的共同单调性定理。

### 3.4 已按正文确认的 schedules

下表是原文配置事实，不应把数值直接移植到不同 Green 归一化。

| 原文实例 | Fourier mask 参数 M_F | β / 两 β | 观察 |
|---|---:|---:|---|
| 2016 Austria，εr=3.5 | 5→7→10 | 6→2→0.5 | 另组 5→3→1 也给相近结果 |
| 2016 Austria，εr=4 | 5→8→10 | 8→2→0.5 | 出现内部伪目标，原文承认该设置达到限制 |
| 2019 例 1、例 2 | 6→10 | 6→1 | 两阶段基线 |
| 2019 例 1第二阶段 | 10 | β_J=6，β_R=1 | 加速 |
| 2019 例 2第二阶段 | 10 | β_J=6，β_R=1；再试 3,1 | 前者过早停滞；后者得到较好收敛 |
| 2019 例 3，改进配置 | 4→8→10 | 6→1→0.5 | 清除两阶段结果中的中央伪影 |

这里 M_F 是论文定义的频域 mask 尺寸，并非“总计 M_F 个电流模式”；3D mask 取八个低频角区。2019 使用 60 个偶极照明、60 个接收位置的三个场分量，合成与反演分别 60³、30³，10% 高斯噪声。停止条件是 APCS 系数相对变化小于 10⁻³。因此 100 次未成熟辅助交替不能代表成熟方法的性能上限。

来源：2016 作者全文与 2019 全文。genui{"citation":{"refs":["turn19view0","turn1view0"]}}

## 4. 核心命题：CIE 材料读出的有限缺陷偏置

### 4.1 实材料精确公式

固定一行 j、q=E+Gj，以及固定 β>0、固定正源权重。记

\[
a=\|j\|^2,\qquad b=\Re\langle q,j\rangle,\qquad c=\|q\|^2>0.
\]

普通 CSI 的实材料理想 LS 为 \(m_0=b/c\)，缺陷 \(e_0=j-m_0q\)。CIE 的实 R 理想 LS 为

\[
R_\beta^*=\frac{\beta(b+\beta a)}{c+2\beta b+\beta^2a}.
\]

只要分母有效、逆变换有定义，恢复材料为

\[
\boxed{m_\beta=\frac{R_\beta^*}{\beta(1-R_\beta^*)}
=\frac{b+\beta a}{c+\beta b}.}
\]

从而

\[
\boxed{
m_\beta-m_0
=\frac{\beta(ac-b^2)}{c(c+\beta b)}
=\frac{\beta\|e_0\|^2}{c+\beta b}.}
\tag{C1}
\]

证明只用展开、实 LS 正交性 \(\Re\langle q,e_0\rangle=0\) 和 \(\|e_0\|^2=a-b^2/c\)。在 c+βb>0 的分支，差非负。这给出**CIE 相对 CSI 理想读出的上偏置**，不等于相对真材料上偏置，也不等于重建更差或更好。偏置量精确由多源共享实材料不可实现缺陷决定。

条件边界：β、当前电流、权重必须固定；不含激活的材料 box、空间正则、TV 或非精确 R 步；这些都会改变直接公式。若分母趋零，公式揭示病态而不提供稳定性。该命题也不是任意两个 CIE/CSI 迭代轨迹的比较，因为它们通常拥有不同的电流。

### 4.2 复材料版本及一源退化

若 R、m 为复，改用 b=\(\langle q,j\rangle\)；普通复 LS 为 m0=b/c，则

\[
m_\beta=\frac{b+\beta a}{c+\beta\bar b},\qquad
m_\beta-m_0=\frac{\beta(ac-|b|^2)}{c(c+\beta\bar b)}.
\tag{C2}
\]

无实数序关系。单一复源/单标量场分量时 ac=|b|²，对任意非零 q 的候选电流都有 mβ=m0=j/q；β 对固定电流读出没有作用。实材料即使只有一个源，仍可能存在相位不相容，不能得到该退化结论。

### 4.3 一阶 β 抵消

在物理 anchor \(j_*=m_*q_*\)，对任意电流方向 V，q 的变化是 GV。令 L*=I−X*G，则

\[
\boxed{Dm_\beta[V]_i=
\frac{\Re\langle q_{*i},(L_*V)_i\rangle}{\|q_{*i}\|^2}.}
\tag{C3}
\]

它与 β 无关。可由 C1 看出：e0 为一阶小量，所以两读出差至少为二阶；也可直接微分得到。这里只需 q_i 非零、\(1+\beta m_{*i}\ne0\)，无需 G 正定或正规。

对于复材料，去掉 Re，结论同样成立。不能由此否定 FBG-R2 的 reduced Gram 目标：该目标改变精确物理状态上的数据度量，并非本节固定候选电流的材料读出。

## 5. 与 A17 的 0/2/4 阶次桥接

以下均为本次推导，A17 agent 与 root 已独立核验。它是机制辨识命题，不是“重建误差四阶收敛”的算法定理。

### 5.1 任意自由辅助电流：一般二阶

在物理 anchor，取直线 \(J(t)=J_*+tV\)，\(Q(t)=Q_*+tGV\)。逐行实材料 CSI 缺陷满足

\[
e_0(t)=t\Pi_{q_*}^{\perp,\mathbb R}LV+O(t^2).
\]

故

\[
m_\beta(t)-m_0(t)
=\frac{\beta t^2\|\Pi_{q_*}^{\perp,\mathbb R}LV\|^2}
{\|q_*\|^2(1+\beta m_*)}+O(t^3).
\tag{C4}
\]

这清楚指出 β 调节的局部差异出现在哪个法向通道，而不是所有电流方向都被相同“减非线性”。

### 5.2 Maxwell 物理切向直线：四阶

给定实材料方向 h，物理切向电流

\[
T_h=L_*^{-1}\operatorname{diag}(h)Q_*.
\]

令 V=T_h。由于 LV=hQ*，沿这条**线性切向外推**有精确恒等式

\[
j_i(t)=(m_{*i}+th_i)q_i(t)-t^2h_i(GV)_i.
\]

因此

\[
e_{0i}(t)=-t^2h_i\Pi_{q_i(t)}^{\perp,\mathbb R}(GV)_i,
\]

及

\[
\boxed{
m_{\beta,i}(t)-m_{0,i}(t)
=\frac{\beta t^4h_i^2
\|\Pi_{q_i(t)}^{\perp,\mathbb R}(GV)_i\|^2}
{c_i(t)+\beta b_i(t)}.}
\tag{C5}
\]

若分母在邻域远离零，则为 O(t⁴)。四阶的是两种理想读出的差；它们相对真实材料曲线 m*+th 的误差仍通常为二阶，不能混淆。

### 5.3 精确有限物理路径：恒为零

若 J(t)=J(m*+th) 是精确 Maxwell 解，则 j_i(t)=(m_i*+th_i)q_i(t)，缺陷为零，两个理想读出都精确恢复 m*+th。A17-X 已有有限响应公式可描述这类路径；本稿不重新主张其有限 transport 公式为新贡献。

三类路径的区别：

| 电流路径 | 实材料法向缺陷 | CIE−CSI 理想读出差 |
|---|---:|---:|
| 任意辅助电流直线 | 通常 O(t) | 通常 O(t²) |
| 物理材料切向直线 | O(t²) | O(t⁴) |
| 精确物理材料路径 | 0 | 0 |

由此，A17 的物理切向增广与 CSI 自由辅助电流增广不能解释成同一种几何操作。一个改善物理响应近似，一个允许暂时违反物理实现关系。

## 6. 候选 F：有可观测曲率条件的 APCS 释放

这是建议进一步评估的完整 formulation，而不是已经证明超过成熟 CIE 的成像算法。其意义是把“当前少模式，稍后多模式”的经验调度提升为可审查的选择条件；是否具有独立论文价值取决于此条件能否在具体 Maxwell 可观测子空间上给出实质性能与成本收益。

### 6.1 六元组

| 分量 | formulation |
|---|---|
| Optimized Variables | 实材料 m 与辅助电流系数 A；每阶段 J=J_D+U_AA。Rβ 作为状态残差坐标，非唯一主变量 |
| Physical Constraints | 多源共享实材料；DPCS 约束；\(J-X(E+GJ)\) 的软约束；阶段输出重新求解或校正为物理解 |
| Objective | 一个固定阶段的 \(\tfrac12\|W_d(SJ-Y)\|_F^2+\tfrac\rho2\|W_sD_\beta(m)[J-X(E+GJ)]\|_F^2+\mathcal R(m)\)。β、ρ、权重在该阶段固定 |
| Update | 对当前子空间作材料—电流 Gauss–Newton 或 variable projection；材料 box 用明确约束；子空间释放须通过下述 profile curvature 条件；阶段切换可单独 continuation |
| Subspace | 保留 Fourier APCS；候选来自 A17 材料方向电流、Maxwell/Krylov 与当前状态残差响应，经 DPCS 正交化。物理配对方向与自由 nuisance 电流要分别记账 |
| Globalization | 子空间 trust/释放条件、有限原数据校正、同预算多条候选路径。中间 lifted 状态可暂时非物理；最终或保留的分支须给出原数据误差与物理缺陷 |

若在 off-manifold 点使用 m 而非 R 变量，真实 Jacobian 包含 \(D_\beta'(m)e\)；不能只保留 Dβ 对状态残差的乘法。物理 anchor 上 e=0，该项才消失。

候选构造要避免一个常见错误：将 T_h 作为任意自由电流后，不能仍把它当成与材料 h 一一配对的物理 tangent。配对步 \((h,T_h)\) 的线性状态缺陷为零；单独开放 T_h 的独立系数属于 nuisance 扩张，会改变材料可辨识曲率。

### 6.2 可证明的 controlled-release 命题

本节结构由 A17 agent 给出，本 agent 核查成熟 CIE/TSOM 全文中未发现同样的显式 profile-curvature safeguard。Schur 消元本身是标准工具，不能单独作为新颖性声明。

在任意一次固定实化、加权局部线性化中，考虑

\[
\min_{h,a}\ \frac12\|r+Bh+Aa\|^2+\frac12h^T\Lambda h.
\]

其中 h 为实材料方向，a 为既有自由电流系数。设新电流列为 C，\(P_A=AA^\dagger\)，\(Z=(I-P_A)C\)，取正交基 Q_Z 张成 ran Z。定义

\[
U=Q_Z^T(I-P_A)B,\qquad v=Q_Z^T(I-P_A)r.
\]

消元新旧电流后，材料二次模型满足

\[
H_{new}=H_{old}-U^TU,\qquad g_{new}=g_{old}-U^Tv.
\tag{F1}
\]

自由辅助电流扩张总能弱化材料 profile Hessian；数据/状态残差下降并不意味着材料信息增加。在可识别材料子空间上，若 H_old 正定且

\[
\eta=\|UH_{old}^{-1/2}\|_2^2\le\bar\eta<1,
\]

则

\[
H_{new}\succeq(1-\bar\eta)H_{old}.
\tag{F2}
\]

对该二次模型的无约束步 p_old、p_new，有精确关系与稳定性界

\[
p_{new}-p_{old}=H_{new}^{-1}U^T(v+Up_{old}),
\]

\[
\|p_{new}-p_{old}\|_{H_{old}}
\le\frac{\sqrt\eta}{1-\eta}\|v+Up_{old}\|.
\tag{F3}
\]

该条件保护的是选定材料子空间中的局部曲率。若仅靠很大 ridge 使 H_old 正定，得到的只是正则模型稳定性，不能声称原测量的可观测信息得到保护。活动 box、非二次正则、强非线性有限步需另外处理。

### 6.3 为什么它还不是全局成功定理

F2 也可能阻止一次本来有利的跨 basin 松弛；保留局部曲率与改善远程搜索存在真实取舍。所以下一步必须比较：同维数但盲目 Fourier 释放、同维数但无约束 A17 模式释放，以及使用 F2 的释放。不能只展示新方法每步局部条件数较小。

最有可能形成独立贡献的版本是：将 C1/C4/C5 给出的材料可实现性法向通道，与 F1–F3 给出的材料曲率损失通道联系起来，构造一个不依赖 GT 的模式选择及有限步可靠性条件。仅“增加 A17 modes”“加入一个 penalty”或“自适应 β”不足以构成该贡献。

## 7. 最小判别实验与当前核验状态

### 7.1 可以立即复现的代数判别

固定一个非奇异离散 Maxwell 型系统和同一物理 anchor，取多源 Q。比较三条路径：任意 V、V=T_h、精确 J(m+th)。对同一候选 J(t)，分别计算普通实 CSI LS 与实 R CIE LS，不运行两条不同优化轨迹。

观测量为 C1 的等式误差、\(\|m_\beta-m_0\|\) 的 log-log 斜率及单体素法向投影。预期斜率为 2、4、零值；用多个 β 检查系数，而非只找最佳 β。小 t 下应从缺陷平方计算差，避免两个相近读出直接相减的消减误差。

A17 agent 已把此检查纳入 `/workspace/scratch/44631410e194/fbg_r3_work/check_a17_math.py` 与相邻 JSON，采用小型复对称 G 和多源系统。其最后一轮步长减半观测阶次为 2.00041 与 4.01631，有限恒等式绝对误差约 10⁻¹⁶。它是有限维代数核验，**不是 24×24 Maxwell 成像 benchmark，也不是超过 CIE 的性能证据**。本 agent 另独立验证过 C1、C3：在 4×4 复对称 G、3 源、β∈{0.2,1,4} 的测试中，有限公式和缺陷恒等式的最大绝对误差小于 4×10⁻¹⁴；一阶近似误差随扰动幅度线性缩小。

### 7.2 真正区分新候选与标量 CIE 的最小成像实验

保留成熟 CIE 的 β schedules、初始化、DPCS 阈值、优化器、停止规则与目标误差归一化。只改变 APCS 的释放规则，至少比较：

1. 成熟 FFT-TSOM 的阶段性 Fourier 释放；
2. 同增量维数的 A17/Maxwell 候选模式，无曲率保护；
3. 同增量维数、受 F2 控制的候选模式释放。

在 user 既有失败高对比目标和成功目标各取一个预注册实例即可作首个判别，但不能据此声称跨结构泛化。β 与所有选择阈值不得利用待测目标 GT。记录原始数据误差、非物理缺陷、材料误差（仅评估）、\(\lambda_{min}\) 或可识别子空间的相对曲率损失、模式数、G/G* 应用次数、全 Maxwell 求解与线性化成本。最终必须对材料重新作原始 Maxwell forward；辅助电流的好数据拟合不能代替该验证。

若 3 只降低 surrogate residual，却没有改善有限物理校正后的 basin 成功率或成本，则 controlled release 仅是稳定性诊断工具，而未达到“超越成熟 CIE”的目标。若 2 较 1 改善但 3 无额外收益，则贡献在模式生成而非此 safeguard。若 3 主要作用是降低维数，要用相同实际维数的 Fourier 对照排除一般容量正则解释。

## 8. 主要失败机制

| 机制 | 可检验征兆 | 解释边界 |
|---|---|---|
| R 接近 1 的逆变换放大 | \(\partial m/\partial R=1/[\beta(1-R)^2]\) 大 | 不能由 R 图像平稳推出材料图像平稳 |
| off-manifold 电流的材料偏置 | C1 缺陷大、理想 CIE 与 CSI 读出分离 | 与材料真误差的方向没有普遍等价 |
| APCS 过度释放 | 联合残差下降、材料 profile 曲率损失大 | 电流可“代替”材料解释数据 |
| 子空间过小 | 强目标的高频内部电流无法表示 | 2016 已直接讨论该冲突，非新发现 |
| 大 β 电流步停滞 | 两 β 下 APCS 相对变化很小却结果差 | 2019 例 2已有实例；停止量不是正确性证书 |
| 当前源弱激发某体素 | c_i=||q_i||² 很小 | 所有逐点 readout 均病态；β 不创造信息 |
| 权重或参数在证明内变化 | 冻结公式无法预测真实一步 | 对固定块的命题不能当全算法定理 |
| 物理模型近共振 | L 条件数大，finite path 接近极点 | CIE/A17 都需明确非奇异域和有限步控制 |

## 9. 最接近的 prior art 与不可越过的新颖性边界

### 9.1 LUCIE：2026 年最直接的材料恢复近邻

Bevacqua、Isernia、Crocco 的 *A Low-Nonlinearity Approach to Inverse Scattering via an Unconventional Exploitation of the Contraction Integral Equation*，DOI **10.1109/LAWP.2026.3660123**。

它已采用两步：先以大 β 的 CIE 联合优化估计电流与近似支撑；再由 \(\widehat Q=E+G\widehat J\) 构造 \(m\mapsto S\operatorname{diag}(m)\widehat Q\) 的线性数据方程，用 Tikhonov/TSVD 恢复材料，可结合所估支撑。第二步不是单纯逐行 J/Q。文中约束 R 为实，实验可加正性；其数值例的高 β 为 10。电流恢复需要大量迭代，不能只拿第二步耗时称整法“几秒”。该工作不提供 user 的 χ≈20、完整 576 变量、成本匹配优越性证明。

作者接受稿全文已取得并阅读；数学 OCR 有损，未取得本地可渲染 PDF。确切 online 日期没有核实：机构库及稿件年份为 2026，Crossref published 仅 year=2026，created=2026-02-02；机构页 2026-01-01 不能当成已核验的正式 online 日期。

机构页：[https://iris.unirc.it/handle/20.500.12318/167028](https://iris.unirc.it/handle/20.500.12318/167028)。

作者全文：[ResearchGate 全文](https://www.researchgate.net/publication/400394332_A_Low-Nonlinearity_Approach_to_Inverse_Scattering_via_an_Unconventional_Exploitation_of_the_Contraction_Integral_Equation)。genui{"citation":{"refs":["turn12view1","turn14view3","turn14view4","turn17view2"]}}

### 9.2 2021 smart rewritings 与 Virtual Experiments

正式题目为 *Quantitative Non-Linear Inverse Scattering: A Wealth of Possibilities Through Smart Rewritings of the Basic Equations*，DOI 10.1109/OJAP.2021.3063248。它不仅回顾 CS-EB、NIE/CIE、Y0，还提出 Y0-NIE-CSI 与 VE-NIE-CSI；后者组合状态改写、虚拟照明重组与集中/圆对称电流正则。因此“反馈 + illumination 重组 + 电流结构”的概念早已有直接先例。虚拟源重组不增加测量信息，但改变正则或目标后仍可能改善算法；不能用信息不增论把这些优化效果全部排除。genui{"citation":{"refs":["turn30view0","turn31view0","turn31view1"]}}

### 9.3 DCTMC：数据相容 T 与局部材料交替

Levinson–Markel 2016 DCTMC I：\(ATB=\Phi\)，\(T=(I-V\Gamma)^{-1}V\)，\(V=T(I+\Gamma T)^{-1}\)。方法允许中间非局部 V，交替处理 T 数据相容性与 V 的对角/局部性，之后再映回 T 并覆盖数据确定块。强制对角步骤可含带距离权重的非对角汇总，不仅是取 diag。

因此，先完成非局部响应、再恢复局部材料不是新想法。局部性投影可破坏数据相容，数据投影可破坏局部性；强非线性情形没有一般收敛定理。完整 T 的稠密储存和反演也不能忽略。genui{"citation":{"refs":["turn14view1","turn21view2","turn22view1"]}}

### 9.4 WRI、variable projection 与 residual Gram

WRI 的典型 penalty 为

\[
\min_{m,u}\ \tfrac12\|Pu-d\|^2+\tfrac\rho2\|A(m)u-b\|^2.
\]

对 u 精确 variable projection 已是成熟结构。2021 WRI* 的作者全文直接给出数据空间 dual、\((\lambda^2I+FF^*)\bar y=r\) 与 relaxed state 的关系。因此“消去状态得到小型 residual Gram”“数据空间 dual”“先估计场再材料更新”均不是新颖性中心。近似消元若忽略其 m 依赖，则不能继承精确 variable projection 的梯度公式。genui{"citation":{"ref":"turn27view1"}}

van Leeuwen–Yang 2025 *An analysis of constraint-relaxation in PDE-based inverse problems*，DOI 10.1088/1361-6420/adab89，给出残差 Gram 加权形式

\[
J_\rho(c)=\tfrac12\sum_j e_j^*(I+\rho^{-1}G(c))^{-1}e_j.
\]

其凸性结果需要完整状态 Hilbert 空间恰由源/测量 Riesz 表示函数张成；配合固定内积与对材料仿射的 A(c)，极限目标为二次凸函数。Remark 6 明确不保证唯一性，所选小空间也未必逼近真实 PDE。因此不能把该凸性条件直接搬到有限照明、完整 576 材料变量的问题。genui{"citation":{"refs":["turn19view2","turn21view0","turn22view0"]}}

## 10. 来源获取审计

“正文可用”表示本次确实读取正文；“摘要/元数据”不算全文。已被网站临时阻断的路径不会被冒称成功下载。

| 工作 | 正文状态 | primary / 作者 URL | 本次可用 web ref 或文件 |
|---|---|---|---|
| Zhong 等 2016 CIE | **FULL TEXT READ**：12 页作者 PDF，web 可读；本地下载 403 | [作者 PDF](https://www.researchgate.net/profile/Marc-Lambert-2/publication/296475773_A_New_Integral_Equation_Method_to_Solve_Highly_Nonlinear_Inverse_Scattering_Problems/links/57385f7f08aea45ee83ea955/A-New-Integral-Equation-Method-to-Solve-Highly-Nonlinear-Inverse-Scattering-Problems.pdf)；DOI 10.1109/TAP.2016.2535492 | 本文 A0 citation；原 web PDF ref turn19view0，后续 turn20view0–2 |
| Zhong–Xu 2019 CIE | **FULL TEXT READ**：全文 XML、数学 MathML 已核 | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8320915/)；[EuropePMC XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8320915/fullTextXML) | `cie2019.xml`、`cie2019_body.txt`；turn1view0 |
| LUCIE 2026 | **AUTHOR FULL TEXT READ**；精确 online 日未核；本地 PDF unavailable | 上述机构页、RG 全文；[Crossref](https://api.crossref.org/works/10.1109/LAWP.2026.3660123) | `lucie_crossref.json`；turn12view1、turn14view3/4、turn17view2 |
| Bevacqua–Isernia 2021 | **FULL TEXT READ relevant theory**：作者接受稿 HTML；机构 PDF 最初可读、后 403 | [作者全文](https://www.researchgate.net/publication/349740905_Quantitative_Non-Linear_Inverse_Scattering_A_Wealth_of_Possibilities_Through_Smart_Rewritings_of_the_Basic_Equations)；[arXiv](https://arxiv.org/abs/2012.15252) | turn30view0、turn31view0/1 |
| DCTMC I 2016 | **FULL TEXT READ relevant theory/algorithm** | [arXiv 1401.3319](https://arxiv.org/abs/1401.3319)；[HTML](https://arxiv.org/html/1401.3319v3)；DOI 10.1103/PhysRevE.94.043317 | `dctmc.html`；turn14view1、turn21view2/22view1 |
| DCTMC II | 摘要；**SOURCE UNAVAILABLE: full text not inspected** | [arXiv 1505.06777](https://arxiv.org/abs/1505.06777) | turn11academia39 |
| WRI* 2021 | **FULL TEXT READ formulation sections** | [作者 HTML](https://slim.gatech.edu/Publications/Public/Journals/Geophysics/2021/rizzuti2020dfw/rizzuti2020dfw.html) | turn27view1 |
| van Leeuwen–Herrmann 2016 penalty | 机构 PDF 可被检索；本次仅核 abstract/conclusion，**不标全文已读** | [机构 PDF](https://dspace.library.uu.nl/bitstream/handle/1874/329651/Leeuwen_Herrmann_2016_A_penalty_method_for_PDE_constrained_optimization_in_inverse_problems.pdf?sequence=1)；DOI 10.1088/0266-5611/32/1/015007 | turn29search6、turn29search8 |
| van Leeuwen–Yang 2025 | **FULL TEXT ACQUIRED**；主要定理及条件已读，未逐字审核所有数值附录 | [CWI](https://ir.cwi.nl/pub/34830/)；[PDF](https://ir.cwi.nl/pub/34830/34830.pdf) | `constraint2025.pdf/.txt`；turn19view2、turn21view0、turn22view0 |
| CSI 1997 | 官方摘要；材料 LS另由 2016全文核验；**SOURCE UNAVAILABLE: original full text not inspected** | DOI [10.1088/0266-5611/13/6/013](https://doi.org/10.1088/0266-5611/13/6/013) | turn24search0 |
| FFT-TSOM 2011 | 原文未独立取得；**SOURCE UNAVAILABLE**；其实现从 2016/2019正文复原 | DOI [10.1109/TAP.2010.2103027](https://doi.org/10.1109/TAP.2010.2103027) | 2019 bibliography + Sec.2 |
| FBE-CIE 2020 | 摘要/书目信息；**SOURCE UNAVAILABLE** | DOI [10.1109/TMTT.2020.2978165](https://doi.org/10.1109/TMTT.2020.2978165) | turn23search3/4 |
| Multiresolution CIE 2020 | 机构摘要，PDF restricted；**SOURCE UNAVAILABLE** | [机构页](https://iris.unitn.it/handle/11572/260153) | turn24search3 |
| Proximal ADMM CIE 2023/24 | 摘要；**SOURCE UNAVAILABLE** | DOI [10.1109/TMTT.2023.3298203](https://doi.org/10.1109/TMTT.2023.3298203) | 只核 metadata/abstract，不归纳其具体收敛定理 |
| LRWI 2021 | 由 branches agent 提供近邻；本 agent 未读全文 | DOI [10.1109/TGRS.2020.3037856](https://doi.org/10.1109/TGRS.2020.3037856)；[arXiv](https://arxiv.org/abs/2006.16492) | rank-2 moment lifting + PDE relaxation，完整核验见其报告 |
| CC-CSI 2017 | 由 branches agent 提供近邻；本 agent 未读全文 | DOI [10.1109/TAP.2017.2673758](https://doi.org/10.1109/TAP.2017.2673758)；[arXiv](https://arxiv.org/abs/1906.10864) | 数据/状态交叉残差，完整核验见其报告 |

注意两个时间陷阱：DCTMC arXiv HTML 的自动生成页面日期不可替代 v3 的 2016 历史；LUCIE 机构库年度占位日期不可当正式 online 日。

## 11. 交给主研究者的建议

主理论可按三个层次书写：已知 CIE 联合残差改写；本次实材料缺陷读出的有限恒等式及 0/2/4 阶次桥接；由此产生的物理切向与自由辅助电流分工、以及有条件的受控释放。把“算法已经超过成熟 CIE”“拥有全局收敛性”“不存在假解”留为未证，不要让一个漂亮局部命题替代它们。

最可辨认的独立贡献是**normal-defect 阶次与材料 profile-curvature 损失的联合设计原则及条件**；成立范围足够精确，也给出可被推翻的实验。若做不到这一连接的实际选择价值，则应将其作为 CIE/A17 机制论文或诊断结果，而不是宣布新的全局反演体系。
