# B1 下一阶段研究报告
## 统一材料后验、可辨识结构与可信多假设成像

研究日期：2026-10-13  
固定审计版本：`migodam/b1-probabilistic-maxwell-imaging@db5aacc9b7f1cc46118a216657df9af8d3d9d1f7`  
本轮执行：原始文档与关键代码审读、定向文献全文核对、数学推导、零新增物理调用的小矩阵验证。  
本轮未执行：新的 Maxwell campaign、GPU 训练、原仓库完整复现、真实仪器校准。

> **研究裁决：推进“统一材料空间中的 nuisance-aware、event-aware Maxwell posterior inference”。不要把下一阶段定义成三个模板增加一个半径，也不要现在训练大型生成网络。**
>
> 最有价值的科学问题是：**哪些材料结构差异在允许材料常数、位置和校准误差变化后仍留下可检测的 Maxwell 信号；为了可靠比较这些结构事件，哪些近似物理计算必须纠正？**
>
> 未知半径是检验这些问题的低成本工具，不是最终架构。拓扑标签应主要用于输出解释和评价，不应成为预先规定的生成专家列表。

---

# 0. 阅读导航与证据等级

本报告分为：1 研究裁决；2 B1 审计；3 输出对象与统一先验；4 理论包；5 nuisance 与物理约束；6 文献和新颖性；7 算法竞争；8 推荐架构；9 最小实验与预算；10 Codex 执行；11 论文定位；12 最终路线图。独立执行 prompt 见 `CODEX_PROMPT_ZH.md`，可复核代数见 `verify_theory.py` 和 `checks.json`。

证据标签含义：

- **Known**：文献或标准概率/线性代数结论，不作为原创性申报。
- **Proven here**：本报告给出自足推导；这表示推导成立，不表示文献首次发现。
- **Conditional**：定理依赖明确假设，尤其 Maxwell 正则性、稳定性或余项控制；这些假设尚未在复杂几何上认证。
- **Numerically checked**：小矩阵或已保存数字核验，不等于新的物理证据。
- **Proposed / Unresolved**：可执行研究方案或仍缺证明/实验的主张。

文献编号 `[R#]` 与仓库编号 `[B#]` 在末尾给出来源和阅读深度。原仓库的自述结果与本轮独立算术、理论推导明确分开。某些大型 JSON 仅核对了相关案例和字段；没有声称逐个数组完整重算。A17、A30直接核对了关键入口/裁决，A12、A23–A29、A31的全部底层材料未在本轮重新审计。

# 1. Executive research decision

## 1.1 B1 已经证明什么

B1 已经实现一个统计目标明确、物理模型真实但支持很窄的演示：给定 16 个复散射通道，在三个固定同心分层球模板及连续复介电常数先验内，计算事件 evidence、事件后验质量、代表图、条件材料不确定性和预测诊断。[B1–B7]

它没有证明任意三维几何可恢复，没有证明生成器能够自动发现未注册拓扑，没有建立全局数值概率证书，更没有证明现实仪器下的概率校准。竞争案例不是“网络生成了三种未知结构”，而是对三个已注册结构做了连续参数积分。

## 1.2 最大瓶颈的排序

第一是 **material support 与可解释事件的定义**，第二是 **非线性结构变化与 nuisance 的可辨识性**，第三才是 **可靠后验探索和积分成本**。在实际非同心三维问题出现之前，尚无证据表明大型神经生成器是关键瓶颈。

已有物理账本记录 384 次调用的 kernel 包围计时总和约 2.27 秒。因此，在此解析 Mie 台架上急于训练网络替代前向，不符合目前成本结构。这并非一般三维 VIE/FEM 也便宜，而是当前实验不能证明它们的加速需求。[B4]

## 1.3 是否马上做 E1

**可以继续现有 E1，但改变其科学问题；不要并行重开一个同类 campaign。** 先读取正在运行的 B1 实验、授权和结果，将本报告的连续参数误差检查、nuisance 消融和有限割线诊断接入现有工作。未知半径实验即使成功，也只通过“连续几何概率参考”门槛，不自动通过论文新颖性或开放世界成像门槛。

首个决定性问题应是：在真实 full-Mie 中，允许 ε 与增益等变化后，孔洞响应剩余多少？该剩余能否预测连续几何 posterior、结构事件 Bayes factor 和增加频率后的变化？如果只能说明“噪声越大越不确定”，贡献不够。

## 1.4 明确选择与停止

**Primary route**：统一材料先验 + tempered population inference + nuisance-aware finite-secant proposals + 可检查的事件概率误差。先无 NN，后按实测瓶颈引入 learned proposal。

**高风险备选**：不优先学习生成器，而是做 event-targeted active sensing，选择最能区分竞争材料事件的新频率/极化/入射方向。

**现在停止的主张**：三个模板概率就是开放世界概率；局部 Fisher 满秩就是全局拓扑唯一；高 ESS 就没有漏模式；全波 residual 小就概率正确；A17 低秩恒等式自动带来全局非线性后验降维；参数少的 Mie 网络加速等价于一般三维成像加速。

B1 负责概率目标、物理参考、理论和校准；C1 以后负责摊销 proposal/transport。两者共享统一材料表示和事件接口，不按孔洞数分设网络。

# 2. B1 critical audit

## 2.1 物理模型与实现

直接审读 `src/b1/physics.py` 和 `configs/protocol.json`。模型使用 `exp(-iωt)`，非磁、各向同性、正实部且非负虚部的介电材料，外部真空，外半径 1、kR=1.6、默认 nmax=10。每个多极阶分别计算电/TM 与磁/TE 系数，匹配切向边界条件；主方法不是 Born 数据生成。[B3,B6,B7]

8 个散射角、两个偏振/投影产生 16 个复振幅，实化为 32 维。这里的“三维”是矢量分层球 Maxwell，不是任意三维网格目标、各向异性材料或三维接收阵列实验。

主求解器做边界线性系统匹配，另有 admittance recursion；两者独立装配但共享部分特殊函数，不应称完全无共同数值依赖。已有 treams 对照是另一实现的有用证据；不过保存的 `FULL_FORWARD_CHECKS.json` 中四个独立记录只有输入、差异和环境，没有 oracle 原始系数/数据。后续 `refine.py` 已增加这些字段，**代码已修复不等于历史数组已补齐**。[B8,B9]

重要的可复用接口是 `Forward.run(..., radii=..., eps_layers=..., k=...)`：可变同心层边界不需要重写电磁内核。偏心空腔不能沿用同心公式假装完成。

## 2.2 概率目标与实现

当前 H0 为实心；H1 为半径 .50 的中心空腔；H2 为半径 .24 的中心空腔加 .46–.60 的真空球壳。H2 不是两个任意位置的小球孔洞。各介质区域共享 ε′∈[2.2,2.8]、ε″∈[.04,.16]，独立均匀 proper prior，模板 prior 各 1/3。[B3,B6]

`inference.py` 的 Gauss–Legendre 节点权重是在连续均匀 prior 上的积分权重，未把 prior 偷换成离散候选库。分支 evidence 用 log-sum-exp，之后乘模板 prior 归一化。实 Gaussian 噪声尺度固定，不随候选响应重标定。对实增益 g~N(1,s_g²) 的边缘化保留候选依赖的 covariance determinant。此次静态审读与小矩阵核验没有发现这些核心公式的归一化错误。[B7]

Top-K 保留第三事件剩余质量；`unmodeled_geometry_mass=None` 是正确语义。`outside_adopted_prior_support_mass=0` 只描述这个先验排除了外部结构，不表示观测证明外部结构不可能。

PPC 采用后验平均的 χ² discrepancy tail。它是模型内 discrepancy 检查，不是经证明均匀分布的 frequentist p-value，也不是先验支持完整性的证书。

## 2.3 关键结果应如何重读

| 项目 | 已存证据 | 可以成立的结论 | 不能成立的结论 |
|---|---|---|---|
| 竞争案例 | H1 50.6418%，H2 40.4265%，H0 8.9317% | 三个注册事件存在竞争 | 自动发现任意材料世界 |
| 信息不足 | 概率接近等 prior | 该实验数据对这三个事件区分弱 | 所有观测配置都不能分辨 |
| 未建模增益 | 错误 H2 获约92.50% | 校准失配可导致错误集中 | 使用 full Maxwell 就能避免模型错配 |
| 增益边缘化 | H1/H2 约49.73%/50.27% | 正确 nuisance 模型恢复合理歧义 | 真结构已唯一恢复 |
| 模板外半径 .68 | 注册 H1 仍获约83.53% | support 缺失不能靠积分补救 | PPC 总能发现所有 support 缺失 |
| Top-2 覆盖 | 86.5278%，报告质量86.6522% | 模型内探索性整体一致 | 已完成真实世界校准 |

这些概率来自经验细化参考。原生四类 reference failure 必须保留；高阶 polynomial-forward likelihood 积分及有限点交叉检查不是 uniform full-Maxwell certificate。[B4,B5,B8]

竞争案例设置的相对实噪声尺度为 .65。它展示高噪声下模型事件竞争，不单独证明低噪声、复杂几何下存在彼此远离的非线性材料模式。下一阶段必须区分连续 ridge、真正分离 modes、prior 所制造的模式以及数值伪峰。

## 2.4 本轮新增审计发现：拓扑质量稳定会掩盖材料方差失真

`distinguishable` 案例的 H0、Re(ε) 保存了以下数值：[B10,B11]

| 量 | 原生 order7 | 经验细化 |
|---|---:|---:|
| H0 概率 | 0.9999999971 | 0.9999999952 |
| 条件均值 | 2.50002435 | 2.51225976 |
| 条件方差 | 3.2290906×10⁻⁶ | 7.3325252×10⁻⁴ |
| 条件标准差 | 0.00179697 | 0.02707864 |

后者方差约为前者 **227.08 倍**，标准差约 **15.07 倍**，而模板概率几乎不变。这是对已保存数字的新算术核验，不意味着细化结果已经是认证真值。

结论是严格的：**只检查 ΔP(H) 和 ΔlogBF，不足以验证整个 posterior。** 离散积分恰好落在一个狭窄峰附近时，模板质量可以饱和而条件矩严重失真。原先 prior-mass-bin quantile 近似又是另一种数值对象，必须独立检查 CDF、均值、方差和分位数收敛。

### 更直接的一致性反例：原生分位数与方差不可能同时精确

同一H0条件后验还报告q₀.₀₅=2.44357576、q₀.₉₅=2.55644933，均值m=2.50002435。如果它们来自同一个精确分布，左右各至少5%的质量给出必要条件

\[
\operatorname{Var}(\epsilon')\ge
0.05(m-q_{0.05})^2+0.05(q_{0.95}-m)^2
=3.1851109\times10^{-4}.
\]

但已报方差只有3.2290906×10⁻⁶，比这一必要下界小约98.64倍。[B7,B10] 所以这组数不能同时是同一个连续posterior的精确summaries。

来源可以从代码解释：covariance来自离散quadrature节点质量，quantile来自把质量铺到prior-mass bins的另一重构。两套不同近似被放进同一个summary，因此出现**summary coherence failure**；这不是evidence归一化公式错误，也不直接否定已验证的竞争案例拓扑质量。

E0应对原生与细化输出都加入这一必要条件检查，并在同一个明确的连续近似target上计算CDF与moments，再分别做积分收敛。小幅Monte Carlo/数值误差需允许容差，但这里约99倍的下界违反不能当作浮点波动。修复可以利用既有paid polynomial-forward bank，不必先花新的Maxwell预算。
## 2.5 对 sampler 与校准的判断

同一暴露 panel 上，Gauss5 名义 cold 成本75状态，两个容差均通过720/720；固定 importance 是27次 proposal setup+48次采样，同成本下两个门槛仅124/720。它是来自单个 y0 的冻结 proposal bank，不是逐数据自适应 SMC、MCMC 或所有 importance 方法的失败。[B4,B8]

720 数据来自24个材料 parent×30次噪声，独立物体数不是720。Top-2 coverage−reported mass 的 parent-bootstrap 区间约[−.04583,.03644]；即使点估计很接近，也仍允许不小误差。原报告保留部分可靠性分箱偏差，并没有通过全面 calibration gate。[B8,B12]

## 2.6 上游连续性，不偷换证据

A17 的正确价值是解释 current modification 如何影响材料 curvature 和 residual pullback。最终入口报告三个完整 nonlinear pair 的 complex-material error 分别下降10.84%、0.22%、7.16%，但 child-wall ratio为2.68、1.74、1.46；不是部署加速。其主要更强证据仍是冻结 material-step fidelity，而不是 posterior sampling。[B13]

A30 固定版本的 Sparse-64 离线容量改善35.69%，真实 L-weighted mean却恶化3.32%；合法补偿找到的最小 full-F 距离仍8.84，未达到原门槛。更重要的是该 A30 平台为24×24二维 scalar TM，不能并入 B1 三维 Maxwell 成功证据。[B14]

因此不能预设“我们已经证明补偿使复杂材料不可辨”。下面的 nuisance-aware 理论要接受可能被新物理实验否定的检验，不能把 A30 的失败改写成支持。

# 3. 从固定模板转向统一材料概率空间

## 3.1 最基本的对象

设 χ 表示物理材料场，ν 为 nuisance，a 为采集配置，先验为 proper 概率测度 μ₀。定义

\[
\Pi^y(d\chi,d\nu)
=Z(y)^{-1} L(y\mid\chi,\nu,a)\,\mu_0(d\chi,d\nu).
\]

这个定义不需要先出现“0孔、1孔、2孔”的分类器。可以在统一 latent field u 上定义 μ₀，再通过 χ=T(u) 推送到材料空间。拓扑只是生成材料之后计算的一个可测量特征。

对一般区域材料，可使用多相 level-set、稀疏多尺度系数或 implicit material field；表面和连续本构参数分开但共同推断。Bayesian level-set 本身是既有方法，不是 B1 的新发明。[R2]

统一不代表没有先验。体素分辨率、level-set 长度尺度、wavelet 尾分布、最大对比度等都定义了支持。有限 latent decoder 尤其不能凭“连续生成”就声称支持所有材料。

## 3.2 每张图的概率应该是什么

定义互不相交的物理事件 A₁,…,A_J，连同剩余集合构成材料空间划分。对每个事件给出

\[
p_j=\Pi^y(A_j),\qquad
\widehat\chi_j\in A_j,\qquad
\Pi^y(\cdot\mid A_j).
\]

图是事件的代表，数字是事件质量，第三项是事件内部的不确定性。单个连续样本不携带可与 LLM 离散 token 概率直接类比的非零点质量。若先验含原子，则某些精确事件可以有原子质量，需要另说。

建议每个 entry 的接口至少包括：事件定义与空间分辨率、代表材料、质量及其数值/Monte Carlo 误差、条件几何与材料区间、prior mass、posterior/prior odds、adequacy 标记和未探索状态。

应区别两个“其他”：

- 已定义模型内、未展示或未探索的事件：有数学后验质量，可能尚未被可靠估计。
- 模型根本没有定义的材料世界：不能填一个凭空的百分比，应标 `unquantified`。

## 3.3 事件如何定义而不重新锁死拓扑

推荐先选一个与应用有关的有限空间尺度 ℓ，以及物理材料距离，例如平滑到 ℓ 后的复介电场误差、区域占据差异、边界距离和材料常数差异的组合。用该距离选 medoid/代表图，再构造带有确定 tie-breaking 的 Voronoi 事件。也可以直接用空腔体积范围、空腔质心区域、区域平均介电常数等可测事件。

这些事件是输出层的表达，不必是 generator 的输入分支。不能只按 measurement-space 距离聚类：那会恰好把“图像很不同但电磁观测相近”的竞争世界合并，丢掉研究对象。应同时报告物理材料距离和观测距离。

聚类依赖尺度与规则。数据自适应划分也可定义合法条件于 y 的事件，但校准时必须重复整个“采样→聚类→选Top-K”过程，不能事后固定最好看的划分再报覆盖率。

## 3.4 “有孔与无孔”不必靠离散标签实现

在单一半径连续先验中，r=0通常质量为0；若希望真正的 solid event有正质量，可使用 spike-and-slab，或者更自然地采用统一 level-set field：整个指定域内u>0是一个可能具有正 prior mass 的函数事件，u跨零时自然出现空腔。这样没有外置孔洞数量专家。

uniform radius与uniform cavity volume不是同一个物理 prior。若 v=(r/R)³ 且v均匀，则r的密度是3r²/R³。坐标不变的要求是：变坐标时推送同一个 prior measure；不是在每种坐标下重新宣称“均匀”然后要求Bayes factor不变。

材料 grammar 亦有多路径问题：多个token序列可能表示同一个物体，需规范化排序或正确累加所有路径的概率。不能把几何重复编码偷偷变成更大的物理先验质量。

# 4. Theory package

## 4.1 T-A：不存在只由连续概率密度唯一决定的“天然Top-K图像”

**状态：Known / Proven here；方法学限制，不作为 Maxwell 原创新定理。**

一维严格正连续密度 p(x) 可通过 t=F(x) 变换为均匀密度。原坐标中的峰，在该坐标中消失。一般变坐标时密度还乘 Jacobian，因此“最高density”“density basin”“density level-set”不是无需指定基准测度的物理对象。

更强地，在标准 nonatomic 概率空间上，不存在既非平凡、又对所有保概率可测同构都不变的有限事件划分：若某一块A与补集都有正质量，可取两边等质量子集并交换；概率测度没变，A却变了。因而不额外给材料度量、分辨率或任务，不可能唯一选出若干物理世界。

这不妨碍事件质量坐标不变。把物理事件A与先验一同推送，Π(A)不变。正确要求应是**事件定义公开、物理意义稳定、概率积分可验证**，而不是幻想后验自动给出唯一的世界分组。

对固定事件划分，选质量最大的K项最小化未覆盖事件概率。若目标是几何保真，应改用

\[
\min_{\widehat\chi_1,\ldots,\widehat\chi_K}
\mathbb E_{\Pi^y}\big[\min_j\ell(\chi,\widehat\chi_j)\big],
\]

并说明诱导的事件划分；这与单点MAP、Top-K density不是同一任务。

## 4.2 T-B：有限非线性材料变化的 exact Maxwell secant

**状态：Proven here（离散模型恒等式）；属于既有 resolvent/反馈代数，不声称新发明。**

在明确的电流离散模型中，设

\[
j_X=X(e^{inc}+Gj_X),\qquad L_X=I-XG,
\qquad e_X=e^{inc}+Gj_X,
\]

且L_X、L_{X+h}可逆。测量F_X=Sj_X，S固定。则

\[
\boxed{F_{X+h}-F_X=S L_{X+h}^{-1}h e_X.}
\tag{1}
\]

证明：新旧方程相减，得到L_{X+h}(j_{X+h}−j_X)=h e_X，左乘S L_{X+h}^{−1}即可。没有线性化h；这是完整有限变化。

若h=UDU*只修改q个局部自由度，记R_X=L_X^{−1}，则

\[
\boxed{
F_{X+h}-F_X
=S R_X U D
\big(I-U^*G R_XUD\big)^{-1}U^*e_X.
}
\tag{2}
\]

这来自push-through/Woodbury恒等式。它包含背景内场注入、背景传播、局部反馈与接收返回，能解释为何直接辐射弱的结构改变仍可产生重要信号。

**算法价值有条件。** 背景R_XU和接收返回要真实付费；q大、全局形变、近共振、许多新位置、频率改变时，缓存或低秩优势可能消失。连续介质中的“局部区域”不是自动有限rank，离散低秩必须伴随分辨率收敛。DDA polarizability与χ的关系若非线性，应将X定义为真正进入离散状态方程的材料算子，不能直接偷换为ε。

本轮随机复矩阵验证：式(1)相对误差2.61×10⁻¹⁵，式(2)为2.64×10⁻¹⁵。它们只核验代数，不验证任何三维网格离散的物理正确性。

## 4.3 T-C：两世界和复合事件的可辨识性必须分开

**状态：Known Gaussian decision identity + Proven here coupling consequence；后者是一般概率工具，不独占 Maxwell。**

先实化并按共同噪声covariance白化。对两个固定完整世界w₀,w₁，有

\[
d^2=\|F(w_1)-F(w_0)\|_{\Sigma^{-1}}^2,
\quad D_{KL}(P_1\|P_0)=d^2/2,
\quad e^*=\Phi_N(-d/2)
\tag{3}
\]

其中e*为等prior二元最小误判率，Φ_N为标准正态CDF。它不是两个拓扑大类的Bayes error。

对于材料事件A、B，令μ_A、μ_B为完整prior条件于事件的分布，包含材料与nuisance；其测量分布是Gaussian混合P_A、P_B。对任意coupling γ∈Γ(μ_A,μ_B)，利用TV的凸性与Gaussian公式：

\[
\operatorname{TV}(P_A,P_B)
\le\int\left[2\Phi_N(d(w,v)/2)-1\right]d\gamma(w,v),
\]

故

\[
\boxed{
 e^*(A,B)\ge
 \sup_{\gamma\in\Gamma(\mu_A,\mu_B)}
 \int\Phi_N(-d(w,v)/2)\,d\gamma(w,v).
}
\tag{4}
\]

证明：用同一coupling表达两边混合密度，三角不等式将混合TV界为成对TV的平均，再代入e*=(1−TV)/2。上式对等event prior成立；其他prior应使用加权检验风险，不能照抄等prior公式。

**意义**：一个很相似、但prior几乎为零的材料对，不能证明两个事件整体不可区分。式(4)要求将相似的完整世界以足够prior质量配对，补上单个最小secant的缺陷。它也明确区分“局部补偿存在”与“复合Bayes factor真正不利于分辨”。

当增益已解析边缘化后，候选covariance不同，式(3)的共同covariance公式不能直接用。但可回到包含增益的完整世界，保持测量噪声Σ共同，再使用式(4)。离散样本上的transport只是近似；其经验最优值不是未经误差控制的总体界。

**可能的Maxwell贡献**：从式(1)–(2)构造物理约束的低成本pair cost、补偿coupling和失效边界；不能把一般Gaussian/transport引理单独报为创新。

## 4.4 T-D：nuisance-profiled 小空腔可辨识系数

**状态：Conditional / Proven here；小体积Maxwell展开本身已有文献，当前对象需验证余项与适用条件。[R7]**

假设半径a的小空腔位于光滑背景内部、远离接口和其他小物体；材料对比度在固定有界范围、频率固定、相关Maxwell resolvent稳定。假设存在展开

\[
F(a,\nu_0)=F(0,\nu_0)+a^3 b+o(a^3).
\tag{5}
\]

ν包含允许补偿的本构/几何/校准变量。设

\[
A=\Sigma^{-1/2}D_\nu F(0,\nu_0),
\quad c=\Sigma^{-1/2}b,
\quad c_\perp=(I-P_A)c.
\]

若nuisance均值流形在ν₀附近正则、参数在约束内部，并只允许局部δν=O(a³)的补偿，则

\[
\boxed{
\min_{\delta\nu=O(a^3)}
\|\Sigma^{-1/2}[F(a,\nu_0)-F(0,\nu_0+\delta\nu)]\|
=a^3\|c_\perp\|+o(a^3).
}
\tag{6}
\]

严格地说，O(a³)集合应选为足以包含最小范数线性补偿的固定有界缩放邻域；不是允许隐含常数任意随a变化。证明：代入(5)及nuisance Taylor展开，令δν=a³v，leading term是a³ min_v||c−Av||，最小二乘投影给式(6)。

这个公式的直觉是：**可检测的不是孔洞引起的全部响应，而是不能被允许nuisance模仿的那一部分。** c的大小与c_perp的大小可能完全不同。

如nuisance受边界约束，应将线性子空间投影替换为切锥上的最小化。若Jacobian退化或全局有远处补偿解，局部公式不涵盖它们。νprior很紧时，profile允许的补偿可能受到prior惩罚，所以式(6)不是marginal evidence公式。

若c_perp=0，leading order完全被吸收；只能结论剩余是o(a³)。**不能没有额外Maxwell推导就宣布普遍a⁵定律。** 本轮玩具函数F(a,θ)=(θ+a³,a⁵)显示raw slope≈3、补偿后≈5，但它只是反例，不能替代Maxwell证据。

在指定两点比较、固定coefficient、独立同质重复观测和噪声σ条件下，d≈C a³√N/σ。要求二元误判率≤α，给出

\[
a\gtrsim
\left[\frac{2\Phi_N^{-1}(1-\alpha)\,\sigma}{C\sqrt N}\right]^{1/3}.
\tag{7}
\]

这里C应是适当nuisance处理后的物理系数；相关通道不能按√N任意累加。式(7)不是统一“λ/某常数”的分辨率定理，也不是复合结构后验质量阈值。近共振虽可能增强内场，也可能放大模型误差并破坏uniform remainder；高对比度随a变化的情形需另行分析。

**下一项关键理论工作**不是再写投影引理，而是给实际 Maxwell 小空腔 b 与余项可核验估计，证明哪些ε/gain/采集配置会让c_perp消失或显著减小，并用全波检验。

## 4.5 T-E：从物理近似误差到材料事件 odds

**状态：Known change-of-measure + Proven here bounds；算法价值取决于能否得到非空误差控制。**

使用同一个proper prior。设近似负log likelihood为Φ_a，full为Φ_f=Φ_a+δ。则

\[
\frac{d\Pi_f}{d\Pi_a}
=\frac{e^{-\delta}}{\mathbb E_{\Pi_a}e^{-\delta}}.
\]

对正质量事件A、B：

\[
\boxed{
\log\frac{\Pi_f(A)}{\Pi_f(B)}
-\log\frac{\Pi_a(A)}{\Pi_a(B)}
=\log\mathbb E_{\Pi_a}[e^{-\delta}\mid A]
-\log\mathbb E_{\Pi_a}[e^{-\delta}\mid B].
}
\tag{8}
\]

证明：分别将事件积分写成Π_a(A)乘条件期望，取比值，归一化常数抵消。

式(8)说明 **影响odds的是不同材料事件的likelihood tilt差异，不是单独的最大forward误差**。在所有世界中相同的δ会完全抵消。

令Π_t∝exp(−Φ_a−tδ)μ₀，有精确导数

\[
\frac{d}{dt}\log\frac{\Pi_t(A)}{\Pi_t(B)}
=-\mathbb E_{\Pi_t}[\delta\mid A]
+\mathbb E_{\Pi_t}[\delta\mid B].
\tag{9}
\]

这是event-aware full-solve分配的合理目标；但“按估计差异分配计算”仍是待验证策略，不因该恒等式就成为最优调度定理。

### 事件质量区间

设A_j构成完整划分，近似未归一化质量q_j=∫_{A_j}exp(−Φ_a)dμ₀，且整个A_j上δ∈[l_j,u_j]。定义m_j^-=q_jexp(−u_j)、m_j^+=q_jexp(−l_j)。则

\[
\boxed{
\frac{m_j^-}{m_j^-+\sum_{k\ne j}m_k^+}
\le\Pi_f(A_j)\le
\frac{m_j^+}{m_j^++\sum_{k\ne j}m_k^-}.
}
\tag{10}
\]

若q_j本身只有积分上下界，可再分别代入相应端点。公式的证明仅用分子单调增、其他分母质量单调减。

如果只探索了部分区域，未探索U必须有其未归一化质量上界才能认证；不能直接从式(10)中删除。标准Gaussianlikelihood上界乘prior(U)给出的界通常很松，特别在evidence很小时可能完全无用。因此本报告**不建议重新把全域uniform certificate作为启动E1的硬前提**。先做exact-corrected inference及经验误差；严格界作为条件通过后的研究支线。

### 从forward error到δ

白化残差r_a=y_w−F_a,w，白化误差e=F_f,w−F_a,w，则

\[
\delta=-r_a^Te+\tfrac12\|e\|^2,
\qquad |\delta|\le\|r_a\|\,\|e\|+\tfrac12\|e\|^2.
\tag{11}
\]

带候选依赖covariance的边缘化likelihood还需要logdet和逆covariance差异，不能只用式(11)。全域oscillation有界时，已知TV≤tanh(oscδ/4)；有限验证点上的最大值不能直接当全域oscillation。

## 4.6 T-F：Maxwell residual与receiver-dual的误差出口

**状态：Known goal-oriented residual identity；可用于式(11)，不是新的误差认证理论。**

离散系统Lj=b，近似解j_a，residual r=b−Lj_a。令L*Z=S*，则

\[
F-F_a=Z^*r.
\tag{12}
\]

对近似dual Z_a，若存在可信β≤σ_min(L)，则

\[
\|(F-F_a)-Z_a^*r\|
\le\beta^{-1}\|S^*-L^*Z_a\|\,\|r\|.
\tag{13}
\]

证明：j−j_a=L⁻¹r，而Z−Z_a=L^{-*}(S*−L*Z_a)。矩阵范数取相容operator norm；必要时可用更保守的Frobenius norm。

这里β往往是瓶颈。passivity本身不自动给出近共振/高对比度下有用的β。精确dual、缓存构建、多源多频的成本都必须计入。若没有β或可验证局部替代，输出只能标“经验误差估计”，不能宣称certified。

本轮仅随机复矩阵测试式(12)，相对误差约4.10×10⁻¹³。

## 4.7 T-G：漏模式、ESS与未知世界的no-go

**状态：Known / Proven here。**

若S是已探索区域，真实后验漏质量η=Π(Sᶜ)，则

\[
\operatorname{TV}(\Pi,\Pi(\cdot\mid S))=\eta.
\tag{14}
\]

证明：在S上密度变为1/Π(S)，Sᶜ上归零，按TV定义积分。对0≤loss≤L_max的决策，使用截断后验选择的Bayes action，其相对真实Bayes action regret≤2L_maxη。

但是η不知道时，有限样本ESS不能估计它。反例：proper uniform proposal下，一个宽10⁻⁴的区间likelihood为10⁴，其余为1，则该小区间后验质量约0.500025。1024个样本完全漏掉它的概率约0.902664；一旦漏掉，所有已见权重相同，ESS=1024，报告的该事件质量仍为0。本轮固定随机种子实际出现了这个情况。

因此“proposal有全支持”和“样本权重看起来健康”都是必要但不充分的实践条件。反过来，**若样本确实来自正确后验且具有有效收敛保证，事件频率当然是合法质量估计**；不能笼统禁止频率。应禁止的是把未经修正/校准的生成频率直接叫物理后验。

### 无法保证识别所有模型外世界

若模型外w_out与模型内w_in产生同一测量分布，任何只基于y的adequacy/abstention规则在两者下都同分布；故不可能同时对内低误报、对外必检出。更一般地，检测power与false-alarm之差至多是两分布TV。

因此必须声明有限尺度、材料contrast和采集限制。任何有限SNR下对任意小孔、任意薄层均保证拓扑正确的主张，在没有最小特征限制时都不成立。PPC可以发现某些失配，但不可能解决这个no-go。

## 4.8 T-H：后验降维不能等同于删除弱材料方向

**状态：Known likelihood-informed reduction；A17 reweighting部分为既有局部身份。**

如果L(y|χ)=g(Tχ)，则条件于Tχ，后验在剩余自由度上的分布等于prior条件分布。一般情形中，对投影P的自然ridge likelihood是

\[
\bar L(Pu)=\mathbb E_{\mu_0}[L(u)\mid Pu],
\tag{15}
\]

不是L(Pu)，不是把弱方向置零，也不是exp(E[logL|Pu])。这与已有认证降维理论一致。[R3]

A17的ΔJ=YΓ⁻¹N可让冻结affine模型改变的tilt依赖双侧签名(Ns,YᵀJ₀s)，因此保留给定签名的条件分布。这仅是 **model-change-only sufficiency**，不是整个nonlinear posterior低维，也不是跨所有几何分支的共同充分坐标。[B3]

局部Jacobian rank、Fisher spectrum和少量latent坐标，不能保证全局mode coverage。一条一维非线性响应曲线也可在32维数据空间具有32维线性span。

若统一材料先验有相关尾部，SMC/MCMC必须保留它们的正确条件采样；projection后冻结weak directions会重复A30固定chart的问题。

### A17到nonlinear posterior的具体余项

为避免仅口头复用A17，给出同一白化、同一anchor下的明确公式。令r=F₀−y，两个局部模型为F₀+J₀s与F₀+(J₀+D)s，D=YΓ⁻¹N。其负log likelihood差是

\[
\delta_{aff}(s)=(r+J_0s)^TDs+\tfrac12\|Ds\|^2
=(Y^Tr+Y^TJ_0s)^T\Gamma^{-1}Ns
+\tfrac12\|Y\Gamma^{-1}Ns\|^2.
\tag{16}
\]

所以reweighting依赖Ns与YᵀJ₀s，加上固定r。若使用不同anchor，还必须保留constant offset；若跨事件使用不同local chart，事件身份和各归一化常数不能省略。

真实非线性模型写成F_i(s)=F₀+J_i s+R_i(s)，记r_i=r+J_i s。则

\[
\delta_{full}-\delta_{aff}
=r_1^TR_1-r_0^TR_0
+\tfrac12(\|R_1\|^2-\|R_0\|^2).
\tag{17}
\]

若在明确chart内二阶导有界，使||R_i(s)||≤κ_i||s||²/2，则

\[
|\delta_{full}-\delta_{aff}|
\le \tfrac12(\kappa_1\|r_1\|+\kappa_0\|r_0\|)\|s\|^2
+\tfrac18(\kappa_1^2+\kappa_0^2)\|s\|^4.
\tag{18}
\]

式(16)来自二次型展开，式(17)是精确分解，式(18)是条件余项界。它们可接入式(8)–(10)，但κ_i、chart外质量和多拓扑转换的边界仍必须被控制。近共振或大几何变化下κ可能很大，证明存在不等于界有用。B1目前两参数rank2，不具降维收益；此推导没有推翻其native A17 NOT_RUN结论。
## 4.9 候选T1–T15的裁决映射

| 原任务候选 | 本报告处理 | 当前状态 |
|---|---|---|
| T1,T3 有限尺度/补偿 | 式(5)–(7)，要求Maxwell系数与余项验证 | Conditional；主理论优先 |
| T2 非线性secant | 式(1)–(4) | 离散恒等式/概率后果成立；物理适用需验证 |
| T4,T9,T10 概率稳定/多保真/反馈 | 式(8)–(13) | 已推导；可用certificate与加速未建立 |
| T5 mode coverage | 式(14)、高ESS漏峰反例 | 成立；η可估性未解决 |
| T6 后验坐标 | 式(15)、局部签名边界 | 既有理论；不支持固定小chart |
| T7 跨拓扑proposal | 统一field MH或正确RJMCMC | 正确性机制已知；Maxwell效率Proposed |
| T8 层级因子化 | future evidence链式分解 | 已知；不构成Transformer必要性 |
| T11 主动测量 | I(A;Y_new|y,a)或Bayes risk | 已知准则；Maxwell目标与算法待检验 |
| T12 scaling/no-go | 低维积分成本、q更新成本、漏峰/支持界 | 部分界成立，无普遍加速定理 |
| T13 prior dominance | posterior/prior odds及mixture区分 | 已知正确对象；需实际敏感性图 |
| T14 Top-K partition | 第4.1节 | 唯一天然划分一般不可能 |
| T15 open-world adequacy | 第4.7节 | 不能统一保证识别观测等价的模型外目标 |

# 5. Nuisance、先验与物理约束

## 5.1 边缘化与profile不能混用

Bayes factor需要积分nuisance，不是各分支把gain/ε拟合到最优后比较最低residual。最大likelihood缺少prior volume和Occam代价。式(6)是局部几何诊断，不能替代该积分。

对线性Gaussian nuisance，设实化模型y=f(χ)+B(χ)ν+n，ν~N(m,Λ)，n~N(0,Σ)，則

\[
y\mid\chi\sim N(f(\chi)+B(\chi)m,\,\Sigma+B(\chi)\Lambda B(\chi)^T).
\]

这可解析边缘化，但B依赖候选时必须保留logdet。复gain可用两个实参数和相应2列实化矩阵表示。正幅度/log-normal、有限区间phase等并不自动服从该Gaussian公式。

## 5.2 对称性与校准

未知gain、phase、source/receiver误差可能吸收材料差异。nuisance维数越大并不必然改善推断；它可能使本来可辨的结构变得不可辨。需要有物理依据的prior尺度、参考target或附加校准测量，而不是“放开所有gain”后继续宣称材料信息不变。

多频测量共享仪器误差时，不应把每个频率当独立重抽的gain。如果把真实共享误差错当独立噪声，evidence会重复计信息。新增频率还应配套合理dispersion model，不可让每频ε完全自由而同时宣传causality约束。

## 5.3 哪些物理约束适合作prior

当前exp(-iωt)约定下，被动正介电材料取Imε≥0。正定/半正定张量、合理的材料范围及有物理来源的色散模型可限制support。互易性首先取决于材料与系统假设；Maxwell boundary conditions通常由forward强制，不应把重复写一项residual当成新增先验信息。

Kramers–Kronig是跨频约束，有限窄带数据不足以无假设地确定带外贡献。可选proper Debye/Lorentz参数先验，但必须明确频带、loss及参数范围。几何最小尺度、区域平滑性和稀有配置的prior，是统计选择，不应包装成电磁定律。

## 5.4 何时abstain

建议独立输出四种状态：`NUMERICAL_UNRESOLVED`、`MODEL_ADEQUACY_WARNING`、`DATA_AMBIGUOUS`、`PRIOR_SENSITIVE`。它们不等价：数值没收敛不能称物理不确定；prior主导不代表计算错；多个高质量模式不代表模型失配。

给出一张“证据强度面板”：prior mass、posterior mass、Bayes factor、prior敏感性范围、PPC discrepancy和数值精度。要量化unknown model probability，必须真的定义该model class的prior和likelihood；否则保持未量化。

# 6. Literature and originality map

本轮是有针对性的研究核对，不是保证覆盖所有截至2026-10-13论文的系统综述。关键区分如下。

| 工作 | 核对内容及物理范围 | 已覆盖的思想 | B1不应重复申报的内容 |
|---|---|---|---|
| Carpio–Iakunin–Stadler，Bayesian topological priors [R1] | 全文方法；2D Helmholtz/偏振成像，形状参数、MCMC、未知对象数model selection | 拓扑信息辅助prior、连续几何不确定性、多模式 | “Bayes+未知孔洞/对象数+多张图”本身 |
| Iglesias–Lu–Stuart，Bayesian level-set [R2] | 方法/定理框架；非Maxwell示例 | 统一隐式函数表达接口、Bayesian posterior与MCMC | “level-set不用预先列拓扑”本身 |
| Zahm等，Certified dimension reduction [R3] | 条件期望ridge与理论正文 | posterior-preserving低维近似及误差界 | “LIS+prior尾部+后验降维”一般原理 |
| Deep Injective Prior [R4] | 正文V及实验；2D介电成像 | latent prior、posterior variational approximation、样本/UQ | “已有人只做MAP、从未做posterior samples”的说法是错的 |
| Back-Projection Diffusion [R5] | posterior测试小节；2D wideband scattering | 物理结构conditioning、conditional sampling、multimodal诊断 | “物理backprojection+diffusion采样”本身 |
| Physics-guided conditional diffusion，2025 [R6] | 方法及选图步骤；2D microwave | 100候选、forward scattered-field MSE选择 | “Diffusion生成很多图+Maxwell选最低residual” |
| DiffusionEMIS，TGRS2024 [R8] | 官方代码README/出版信息；全文未完整取得 | 3D EM diffusion、point-cloud/electrical parameter相关先例 | 不能声称“首个3D EM diffusion”；概率校准细节仍Unresolved |
| Del Moral–Doucet–Jasra，SMC [R9]；Christen–Fox，DA [R10] | 方法和正确性机制 | tempered population、normalizer估计、近似筛选加精确修正 | SMC/RJ/DA的普通组合 |
| Modrák等，2025 SBC [R11]；Säilynoja等，2026 posterior SBC [R12] | 文章方法/结果 | data-dependent test quantities及条件于观测的计算检查 | 单独增加SBC并非新理论 |

几个需要准确表述的细节：

**Carpio等并未简单把同一份数据重复使用两次。** 全文第3节前明确将数据分为一部分生成topological prior、另一部分形成posterior likelihood。B1使用measurement-conditioned proposal时可以直接保留原prior并用正确q比修正，不必照抄其prior设计。[R1]

**Deep Injective Prior确实进行posterior建模。** 它在latent中用围绕MAP的对角Gaussian variational family，并固定均值为z_MAP、优化尺度；材料空间的非线性pushforward可复杂，但该近似本身不是分离latent modes完整质量的独立认证。[R4]

**Back-Projection Diffusion确实检查多模式与概率相关指标。** 其posterior小节明确没有ground-truth posterior，用多样本data misfit等诊断；不能将其说成完全不做UQ，也不能把这类结果当成已经解决结构事件mass calibration。[R5]

**新的SBC检查也要选对量。** 仅对材料参数做rank直方图，可能漏掉“算法忽略数据，永远返回prior”等错误。下一阶段同时检查joint log likelihood、几何–材料补偿函数、事件指标以及连续参数；对离散/tied量做随机化处理。[R11,R12]

## 6.1 最可能存在的研究空缺

不是“未知半径+SMC+Top-K”。更有希望的是以下联合问题：

> 在真实三维矢量Maxwell、连续几何和合理nuisance条件下，建立结构事件的可辨识性系数与非线性有限secant联系，并将接收相关物理误差转化为对这些事件概率/决策的可检查控制，使统一材料后验探索在相同概率精度下真正减少物理成本。

其中各基础构件多数已有先例。论文必须在Maxwell系数、非空有效条件、跨几何传播或计算收益上作实质推进。当前检索没有证明上述整体方向无人做过；提交前仍需针对最终得到的具体定理与算法做定向novelty核对。

# 7. Competing algorithm families

## 7.1 A：Adaptive integration / response surrogate

适合3–4个连续参数的校准台，先用精确Mie或低误差response surrogate，再自适应积分likelihood/CDF。复杂度tensor约q^d；sparse grid在有足够混合光滑性时可减轻，但拓扑事件边界、排序约束、共振和窄峰会破坏优势。不存在不看posterior结构就能给出的统一“超过6维必失败”门槛。

若用surrogate，必须分别核验forward误差、likelihood积分误差和连续CDF误差。Gauss多次细化很稳定仍不自动给uniform certificate。它是参考方法，不是最终开放世界生成架构。

## 7.2 B：统一latent的tempered SMC / parallel tempering

设u为统一材料变量，γ_β(u)=μ₀(u)L(u)^β，β从0升到1。初始化按proper prior，增量权重L^(Δβ)，必要时resample，然后用对当前γ_β正确不变的MCMC move。归一化常数由增量权重平均的连乘估计；标准条件下有相应理论，但有限样本normalized ratio及adaptive实现仍需单独分析。[R9]

level-set在固定latent维数下即可改变孔洞拓扑，因此不必为了跨拓扑硬上RJMCMC。针对Gaussianprior的pCN move可保持prior；跨大几何需要多尺度/global move和tempering，不能只沿局部GN椭球走。

主要风险是低prior质量但重要mode根本没有粒子、temperature跃迁过大、resampling祖先退化和local mixing。多seed一致性、ESS、acceptance都只是诊断；需要小维独立reference、mode stress cases和已知failure对照。

## 7.3 C：Maxwell-informed finite edits + delayed acceptance

proposal可以利用式(1)–(2)的finite edit响应，而不只用Jacobian预测大拓扑变化。局部形状/区域编辑以接收residual、nuisance投影及多尺度prior引导，同时保留global prior moves。

若proposal密度q可计算，标准MH用完整prior、likelihood与反向/正向q比；真正变维时使用RJ维度匹配、birth/death概率和Jacobian。选择fixed-D统一field时不需要人为维度跳跃。

Delayed acceptance先用近似targetπ_a筛选：

\[
\alpha_1=\min\{1,\pi_a(u')q(u|u')/[\pi_a(u)q(u'|u)]\},
\]

再对通过者做full evaluation：

\[
\alpha_2=\min\{1,\pi_f(u')\pi_a(u)/[\pi_f(u)\pi_a(u')]\}.
\]

在support与proposal条件满足、按正确方式组合时，保持full target。该正确性是已有DA机制，不是B1发明。[R10]

所谓exact correction表示对声明的数值full solver target正确，仍有continuum discretization和现实模型误差。在线不断改surrogate/metric需冻结epoch、使用适当adaptive理论或明确有限适应；不能默认所有state-dependent更新都保留target。

**采用条件**：实际节省full solves且最终事件/CDF误差不恶化。**放弃条件**：RbU/adjoint准备比直接solve更贵、stage2拒绝率高、proposal使mode丢失、效果只在已知相同几何缓存下出现。

## 7.4 D：learned proposal / transport / autoregression

C1优先学习proposal，而不是宣称网络softmax等于物理posterior。Autoregressive discrete-continuous模型可显式计算q；normalizing flow可提供latent density；所有与prior的变换需保持同一measure。含λμ₀的defensive mixture保证support，只保证不为零，不保证有限预算命中rare modes。

diffusion/flow-matching产生样本，不一定自动提供便宜且精确的endpoint density。没有q就不能机械使用L/q；应采用有证明的path-space correction、可计算密度transport或声明approximate sampling并验证。训练时把physics residual加到loss也不是Bayesian正确性的替代。

自回归条件概率严格可由future evidence定义：Z(s)=Σ_t p₀(t|s)Z(st)，p(t|s,y)=p₀(t|s)Z(st)/Z(s)。困难在子树积分与漏质量，不在是否调用Transformer。beam search截断必须保留剩余质量不确定性，不能把保留下来的路径重新归一化冒充完整posterior。

启动NN前计算break-even：若训练/数据setup成本为C_train，传统每例成本C_base，learned+correction每例C_new，则摊销需N(C_base−C_new)>C_train且精度相同。若C_new≥C_base，无论生成图片多漂亮都没有计算理由。

## 7.5 E：主动测量作为备选

对事件A的后验不确定性，可选测量a最大化I(A;Y_new|y,a)，或最小化posterior decision risk。它不同于最大化所有参数Fisher determinant，后者可能主要学习gain等nuisance。

如果一个合理新频率能显著区分两个竞争世界，先测比训练大生成器更有实际意义。比较时计入测量代价、共享校准、噪声相关性和新频率dispersion先验；不能只挑一对MAP的最大距离当整个posterior的EIG。

# 8. Recommended architecture

## 8.1 第一版，无NN

1. **统一材料prior**：多尺度implicit/level-set field和区域本构参数；无孔洞数输入标签；保留适当prior tails。
2. **Population posterior engine**：tempered SMC，global prior move、几何move和物理引导move混合。
3. **物理proposal模块**：局部small edit先用Jacobian，较大edit使用finite secant/局部反馈；所有proposal保持正确q或正确不变核。
4. **精度模块**：低保真筛选、高保真校正；先经验验证，能够获得有用界时才加入事件odds区间。
5. **事件输出模块**：在声明材料尺度上选择代表与事件，报质量、分支内不确定性、未探索/模型外状态。
6. **可靠性模块**：prior sensitivity、PPC、likelihood-aware SBC、continuum solver核验和独立test场景。

这是同一posterior架构，球、孔洞、多个区域只是它可表达的材料结果。概率参考的校准台可以更简单，但不能用校准台的表示支持冒充最终开放模型。

## 8.2 A17保留什么

保留物理的材料注入–反馈–接收分解、精确finite-change响应和receiver-dual误差出口。只有成本实测通过时才保留online反馈低秩更新。停止强制current selection、固定OPM chart和“native A17必然加速”的要求。

## 8.3 C1的条件升级

当某类proposal的full-physics准备很贵但跨对象重复、或传统population在已核验的multimodal benchmark上明显低效，才让C1学习该proposal/transport。网络可输入measurements、current/field、adjoint残差、几何和frequency state；它输出统一材料编辑或latent proposal，不直接输出未经验证的event confidence。

# 9. Minimal decisive experiment package

## 9.1 E0：零新增物理的审计与接入

读取当前B1运行状态，确认已有E1进度及授权；不启动重复任务、不修改旧384账本、不打开旧locked truth。

复核现有数组上的evidence、gain logdet、事件总质量、CDF与moments；新增本报告的variance比检查。核对代码修复与历史artifact缺口，明确不能由新代码推定旧数据已存在。

补齐单位测试：变坐标保持同一prior时odds不变；radius-prior与volume-prior改变时结果确实可变；pure-prior inference负控制；高ESS漏峰；事件partition稳定；representation/unknown mass语义。现有测试要先检查是否触发physical调用，零预算阶段不可盲跑`audit/pilot/refine`。

输出：`E0_AUDIT.md`、`RUN_STATE.json`、`AUTHORIZATION.json`、新阶段ledger策略、可复用资产清单及未满足的source/reference项。

## 9.2 E1-M：可辨识机制，而非仅未知半径

### 物理与场景

优先接入正在运行的variable-radius实验；复用已声明3D vector-Mie。固定外边界，比较central cavity radius变化、复ε变化、实gain已知/未知与多频观测。kR=1.6保留原接点；另选小kR和较高kR作机制对照，先核验nmax收敛，不能把小体积/低频asymptotics扩展为任意频率定律。

最小未知量为r、ε′、ε″，gain可解析边缘化。事件采用预注册有限尺度radius/void-volume区间；若存在solid atom明确给其properprior。如果现有E1已做该模型，不再从零重复。

### 必须比较的量

raw finite secant、只允许ε补偿后的真实full secant、再允许gain后的secant；局部投影预测；完整marginal evidence与posterior radius–ε相关性；增加频率前后的变化。

**不要把profile minimizer直接当后验。** 同时报告补偿所需参数是否靠近prior边界、prior mass是否很小。A30未发现强歧义，故不以“必须找到补偿导致不可辨”作为筛选真值的理由。

### 对照与失效条件

用不同噪声尺度、半径大小、material contrast和near-resonant例检验；每类保留失败。若c_perp预测没有超出raw distance的解释力，或在真实Maxwell中未进入其asymptotic范围，降级或终止“补偿辨识新理论”主张。

若r posterior是一条连续ridge或单峰，不得人为切成若干“发现的modes”。这仍可用于验证合法事件质量，但不能声称多峰生成成功。

## 9.3 E1-U：最小统一生成表示，不按拓扑分类

作为低成本表示验证，可用径向latent field

\[
u_z(s)=1+z_c B_c(s)+z_sB_s(s),\quad s=r/R,
\]

其中B_c为中心附近光滑compact bump，B_s为内部球壳附近compact bump，外边界邻域保持u=1。z_c,z_s用proper连续Gaussian或明确有界prior。u<0表示vacuum，u>0表示共享ε介质；root求解得到同心层radii，调用已有Mie接口。

这会在同一个连续latent模型内出现solid、central cavity、annular cavity及组合，而不是先选择H再采样。之后才计算void topology作评价。两个geometry latent+复ε共4个实参数，仍可能用自适应积分与独立quadrature交叉参考。

**限制必须写明**：这只是统一架构的径向校准切片，中心和壳位置仍受basis限制，不是任意三维结构生成。若prior自然产生非常薄层或接近合并的roots，必须预先定义最小尺度规则/退化层处理及其诱导prior，不能数值失败后删掉困难样本而不报告。

## 9.4 E2-Q：非同心三维solver资格，条件通过才推断

最低有意义升级是偏心空腔或非同心多区域。优先评估nested vector-spherical-wave boundary matching/translation方法；它保留球面界面的解析优势，但内嵌cavity不等于普通同背景多球cluster。不能未经核验就假设treams的标准cluster API已直接支持该嵌套物理。[R13]

备选是真正3D VIE/DDA或BEM/FEM独立参考。需要检查偏心→同心极限、vacuum/homogeneous极限、multipole/mesh收敛、passive energy、multiple illumination与独立实现一致性。

先做solver资格，不同时启动高维posterior和大型training。资格失败则不做E2 posterior，保留失败与成本。

## 9.5 E2-P：后验算法和物理贡献

在solver资格、prior支持及参考可计算性通过后，比较：generic tempered SMC/pCN；几何metric proposal；Maxwell finite-secant proposal；Maxwell proposal+full correction。不能把额外adjoint/RbU成本隐藏在proposal中。

至少选择一个独立参考可处理的低维非同心切片，其余较高维只称reference-limited study。多个seed一致并不构成高维真posterior；若模式质量无法独立验证，禁止宣称概率加速。

reference重点不是更多漂亮图，而是每个有实质mass的事件、关键CDF与rare-mode stress case。展示conditional材料sample或medoid，并对代表图做forward核验；条件均值即使合法也未必是高likelihood典型样本。

## 9.6 统一评价协议

| 方面 | 主指标与检查 |
|---|---|
| 概率 | event mass绝对误差、重要event的log odds、normalizer误差、Monte Carlo误差与solver误差分列 |
| 连续不确定性 | 关键CDF sup差、分位数、均值与SD收敛，尤其狭窄峰 |
| 探索 | 独立种子、known-mode stress、prior-only负控制、mode遗漏而非仅ESS |
| 成像 | 未知geometry、within-event uncertainty、finite-scale物理材料距离；不要只报PSNR |
| 校准 | 独立material parents、parent-cluster bootstrap、likelihood-aware SBC、Brier/log scores、PPC和mismatch |
| 成本 | cold/warm、setup/cache、full/adjoint、多RHS、失败/重复调用、wall、CPU/GPU memory |
| 可推广性 | held-out geometry families、support边界、gain/phase/noise/dispersion失配、不同resolution |

建议的**研究容差，不是已证明界**：主要事件mass误差≤.01；mass≥10⁻³的事件log-odds误差≤.05；关键normalized CDF sup差≤.01；主要事件均值/SD在refinement中按预注册尺度稳定。极小概率事件单独报上界/绝对误差，避免用无意义relative error驱动无限计算。若资源不足，不放宽到产生正结果，而是标reference unresolved。

独立parent初期可先64个，每个少量noise；这是探索，不足以宣称全面现实校准。稀有几何需单独stress设计，并与prior-predictive加权评价分开，不把人为oversampling频率当真实prior。

## 9.7 staged预算建议

**下表是下一阶段建议，不是本轮已执行或对旧账本的扩容。** 每个state按材料×频率计费；额外incidence/多RHS计算同时记录，不能只报一个粗state count。先以少量资格calls测solver实际成本，时间或计数上限谁先到就停。

| 阶段 | 建议上限 | 分配与权限 |
|---|---:|---|
| E0 | 0新增物理 | 已有source/数组与小矩阵代数；可立即做 |
| P1：E1-M及必要E1-U pilot | 4,000 Mie材料–频率状态 | 400机制/物理资格，2,000自适应reference和算法，400fresh truth/parents，400独立/端点核验，800失败与追加refinement；总上限不因类别未用完而抬高 |
| P1参考不足时 | 暂不追加 | 先报告何处reference失效，再决定是否值得单独追加；不能自动转训练 |
| E2-Q | 128非同心qualification states | 独立新授权；CPU/GPU和时间上限按首次benchmark确定 |
| E2-P | 建议7,000 states | 2,000探索/算法、4,000概率reference、1,000失败和复核；仅资格通过且预测可支付时另行授权 |

P1默认CPU，GPU0；单次任务建议设置2小时host wall保护上限，远低于该上限也可因科学门槛失败提前停。2小时是允许消耗的上限，不是运行时预测。新kR、更多层、困难prior会改变单次成本，不将历史2.27秒线性外推成承诺。

E2可用RTX4060 8GB，但只有在矩阵/场storage实际可容纳时；不得为使用GPU而强行选择更复杂solver。large NN training仍不在以上授权范围。

## 9.8 Go / No-Go

G0：prior、measure、likelihood、numerical uncertainty、ledger语义无未解决核心错误。失败只修复，不做算法优越性结论。

G1：统一径向模型确实自然跨事件，并有独立连续posterior参考。通过只说明架构切片成立；不等于一般3D。

G2：nuisance-aware full-secant指标对预注册实验的odds/歧义变化具有raw distance没有的解释价值，且失败条件吻合。否则终止将该机制作为主贡献，不能把σ调大到所有case模糊来“验证”。

G3：在同等事件/CDF精度、相同支持与全部成本下，物理proposal有一致收益。建议预注册至少20% inclusive wall或full-solve节约作为工程推进阈值；该20%是研究决策阈值，不是理论常数。若只在Mie上看不到加速，不必否定一般solver，但也不能据此自动启动昂贵E2。

G4：非同心solver与支持资格通过。否则不许声称一般3D，不启动C1大型训练。

# 10. Autonomous Codex implementation plan

建议新建隔离目录 `research/b1_next/` 或新branch，原结果只读。先读取源入口、当前运行记录和本包，再写执行plan；不覆盖正在运行的E1，不与其重复统计同一物理状态。

模块建议：`material_prior.py`、`event_partition.py`、`finite_secant.py`、`nuisance.py`、`population_inference.py`、`reference.py`、`probability_diagnostics.py`、`resource_ledger.py`。这些名称是建议，不要求为清单创建空壳。

每次迭代记录：科学假设、唯一主要改变、比较基线、已付费与新增资产、结果、失败原因、下一步裁决。预算内允许改proposal、temperature、积分策略；改变prior/support/noise model或主评价规则时新建protocol版本，不能把不同posterior拼成同一个实验。

粒子、proposal setup、truth generation、独立solver、失败、multifrequency和adjoint都进入账本。cold/warm成本分离；固定几何跨观测复用只在warm栏体现。所有claim同时有source路径、状态标签和计算范围。

详细可复制prompt另见 `CODEX_PROMPT_ZH.md`。默认只自动批准E0；后续P1的数值预算在prompt配置块中显式设置。未授权时完成可做的零物理工作并输出决策，不反复询问或暗中重置预算。

# 11. Paper positioning

建议工作题目：

**Nuisance-Aware Material-Event Inference for Three-Dimensional Maxwell Inverse Scattering**

只有取得真正有效的概率误差界后，才在标题使用Certified；只保持full discrete target的MH correction，不等于continuum或现实校准认证。

三个贡献候选：

1. **物理理论**：小空腔/区域改变的full-Maxwell secant，在材料和校准补偿后如何留下可辨识结构信号；给出系数、余项范围及实测失败边界，而非仅写投影矩阵。
2. **概率计算**：接收相关物理误差如何影响material-event odds，并形成实用的验证/修正与计算分配机制；严格证书若无效则诚实保留为经验方法。
3. **统一成像**：无需预设拓扑生成标签，在非同心三维连续几何中输出多个代表、可靠事件mass及内部uncertainty，并对prior/support/nuisance和多模式覆盖做独立检查。

审稿人最可能拒绝的原因：Mie台架太简单；经典Bayes/SMC/DA重新命名；所有模式来自手工事件列表；所谓概率来自未校准生成器；对同一approximate solver做self-consistency验证；成本不含setup/adjoint；复杂场景只展示图而无可靠mass reference；论文不清楚哪些条件真正Maxwell-specific。

达到强论文标准需要非平凡物理结论及至少一种一般基线没有的可靠能力。只提高PSNR、生成更多样本、把未知半径算出来，都不足以单独建立上述贡献。

若主机制失败：优先转向 **主动结构区分与不足以判断的成像**，研究有限测量下哪些additional measurements真正解决材料歧义。若物理proposal不省成本但概率参考有效，保留为benchmark/工具，不持续将其包装成加速论文。若统一prior严重漏support，优先扩大材料表示；不能继续优化一个排除真值的sampler。

# 12. Final ranked roadmap

**第一优先**：接入正在运行的B1，完成E0；补上continuous uncertainty误差门槛和已有证据缺口。

**下一项定理**：式(6)的Maxwell系数/余项与实际采集配置，连到式(4)的复合事件区分，而不是只研究local Fisher spectrum。

**下一项决定性实验**：同一variable-radius/full-Mie台架上比较raw response、ε补偿、ε+gain补偿以及多频后验；不预设一定出现多个modes。

**下一项算法实验**：统一implicit材料prior下的tempered inference，先与独立参考比event mass与CDF，再看物理proposal的全成本收益。

**高风险备选**：event-aware active measurement；只有通用proposal确实成为瓶颈才进入C1 learned transport。

**终止路线**：固定孔洞数专家作为最终架构；未纠正diffusion频率作为posterior；在二维变量上强推A17降维；再以一个未经认证的uniform norm界作为所有研究的启动前提；以低residual单赢家替代多世界概率。

最终判断不是“下一步上SMC还是Diffusion”，而是：**先把何种材料差异可由数据决定、哪些概率数字可信这两件事接起来。统一生成是表达方式，Maxwell结构辨识与事件概率控制才是最值得争取的科学贡献。**

---

# 附录A：本轮数值验证与局限

运行 `python verify_theory.py`，依赖NumPy/SciPy，固定seed20261013。脚本不导入B1solver，不访问locked数据，不触发Maxwell调用。

| 核验 | 本轮结果 | 解释 |
|---|---:|---|
| Gaussian gain解析边缘化 vs dense covariance | 最大logL差8.19×10⁻¹¹ | 浮点代数一致 |
| finite material secant | 相对误差2.61×10⁻¹⁵ | 离散恒等式 |
| 局部feedback更新 | 相对误差2.64×10⁻¹⁵ | 非Maxwell随机矩阵 |
| receiver-dual residual identity | 相对误差4.10×10⁻¹³ | 未验证真实solver β |
| event reweighting | 最大差5.55×10⁻¹⁷ | 条件期望公式 |
| log-odds导数 | 差4.75×10⁻¹² | finite-difference核对 |
| event bounds | 所有真实toy masses落入区间 | 界很宽，不能宣称实用证书 |
| raw/profiled toy slope | 3.00005 / 5.00000 | 非Maxwell反例 |
| rare-mode | 真mass .500025；ESS1024但estimate0 | full support+高ESS不足 |
| 复合Gaussian coupling | 数值Bayes error 0.258543；下界0.219930 | 混合分布验证，非Maxwell |
| A17 affine tilt / nonlinear remainder | 误差2.36×10⁻¹⁶ / 3.65×10⁻¹⁶ | 显式有界余项核验 |
| B1原生/细化variance比 | 227.0771 | 已保存摘要的算术，不是新reference |

首次写JSON时遇到NumPy boolean序列化错误，已在脚本转换为原生类型后重跑；该错误不涉及物理或理论计算。最终所有asserted checks通过。

文献PDF解析文本可访问，但数次web截图返回Internal Error；本报告没有据此作无法核验的图形测量或视觉比较。部分关键论文使用可访问HTML正文核对。DiffusionEMIS仅核对官方代码说明与出版信息，不宣称已审完其全文。未完整下载仓库二进制bank，因此原项目报告的所有数值不等于本轮全部独立复现。

# 附录B：源清单

## B.1 固定B1仓库来源

下列[B1]–[B12]均相对于：
`https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/db5aacc9b7f1cc46118a216657df9af8d3d9d1f7`

[B1] `START_HERE.md`；完整入口审读。  
[B2] `RESEARCH_REPORT_ZH.md`；主体与尾部核对。  
[B3] `THEORY.md`；公式与声明审读。  
[B4] `RESULTS.md`；结果及成本表审读。  
[B5] `DECISION.md`、`GATE_DECISION.json`、`RESEARCH_LOG.md`；科学状态、历史失败和预算。  
[B6] `configs/protocol.json`；prior、noise、物理与容差。  
[B7] `src/b1/physics.py`、`src/b1/inference.py`；直接代码审读。  
[B8] `research/delegated/probability_review/REVIEW_NATIVE.md`；关键review段落，不声称全部附属JSON重新计算。  
[B9] `src/b1/refine.py`、`results/refinement/FULL_FORWARD_CHECKS.json`；新保存代码与历史oracle字段缺口直接对照。  
[B10] `results/PRIMARY_SUMMARY.json`；相关case及字段核对。  
[B11] `results/refinement/SUMMARY.json`；对应case及reference字段核对。  
[B12] `README.md`；入口、命令、预算保护和复现范围。

[B13] `https://github.com/migodam/a17-maxwell-certified-exchange/blob/main/closure_final/START_HERE.md`；本轮读取main所展示的closure入口，不声称固定A17完整物理重放。另读根README了解efficiency followup范围。  
[B14] `https://github.com/migodam/a30-maxwell-observable-material-manifold/blob/717d3d26a290248bd2b9aeb55f3d03694b407b84/A30_E0_E1_DECISION.md`；固定裁决关键正文。

## B.2 定向文献

[R1] Carpio, Iakunin, Stadler. *Bayesian approach to inverse scattering with topological priors*. arXiv:2003.09318v2. 正文§2–3、模型和方法核对；2D范围、data split及model selection。`https://arxiv.org/abs/2003.09318`

[R2] Iglesias, Lu, Stuart. *A Bayesian level set method for geometric inverse problems*. Interfaces and Free Boundaries 18(2),181–217,2016. DOI:10.4171/IFB/362. `https://arxiv.org/abs/1504.00313`

[R3] Zahm, Cui, Law, Spantini, Marzouk. *Certified dimension reduction in nonlinear Bayesian inverse problems*. Mathematics of Computation 91,1789–1835,2022. DOI:10.1090/mcom/3737. 条件期望ridge理论正文核对。`https://arxiv.org/abs/1807.03712`

[R4] *Deep Injective Prior for Inverse Scattering*. arXiv:2301.03092，取得的PDF版本正文§V–VI及局限核对。`https://arxiv.org/abs/2301.03092`

[R5] *Back-Projection Diffusion: Solving the Wideband Inverse Scattering Problem with Diffusion Models*. arXiv:2408.02866v1；尤其§5.9 posterior evaluation。`https://arxiv.org/abs/2408.02866`

[R6] *Physics-guided Conditional Diffusion Networks for Microwave Image Reconstruction*. arXiv:2510.25729v1；正文latent diffusion和forward-MSE candidate selection核对。`https://arxiv.org/abs/2510.25729`

[R7] Ammari, Vogelius, Volkov. *Asymptotic formulas for perturbations in the electromagnetic fields due to the presence of inhomogeneities of small diameter II. The full Maxwell equations*. J. Math. Pures Appl.80(8),769–814,2001. DOI:10.1016/S0021-7824(01)01217-X. 原作者机构出版记录及相关原始研究核对；本文没有据未全文核对的假设直接宣布uniform remainder。相关：Ammari, Moskow, Vogelius, ESAIM COCV9,49–66,2003, DOI:10.1051/cocv:2002071。

[R8] Bi, Chen, Li. *DiffusionEMIS: Diffusion Model for 3-D Electromagnetic Inverse Scattering*. IEEE TGRS62,1–12,2024. DOI:10.1109/TGRS.2024.3349681. 官方代码：`https://github.com/Xueting-Bi/diffusionEMIS`。全文概率机制未完整核对。

[R9] Del Moral, Doucet, Jasra. *Sequential Monte Carlo samplers*. JRSS B68(3),411–436,2006. DOI:10.1111/j.1467-9868.2006.00553.x。

[R10] Christen, Fox. *Markov Chain Monte Carlo Using an Approximation*. JCGS14(4),795–810,2005. DOI:10.1198/106186005X76983。

[R11] Modrák et al. *Simulation-Based Calibration Checking for Bayesian Computation: The Choice of Test Quantities Shapes Sensitivity*. Bayesian Analysis20(2),461–488,2025. DOI:10.1214/23-BA1404。出版页/方法结论核对。

[R12] Säilynoja et al. *Posterior SBC: simulation-based calibration checking conditional on data*. Statistics and Computing36,78,2026. DOI:10.1007/s11222-026-10825-9；2026-01-31正式发表。正文的likelihood test quantities与posterior-local检查核对，不将其当现实模型正确性的证明。

[R13] treams官方T-matrix文档与Beutel等CPC2024方法背景。`https://tfp-photonics.github.io/treams/tmatrix.html`。层球/常规cluster能力可查；没有确认其普通cluster调用可直接完成嵌套偏心cavity。

## B.3 未完成的验证，不应隐去

尚未完整逐源重审A12、A23–A29与A31；尚未检查所有原始bank和所有topology review附件；未对所有相关2026生成式论文完成全文系统综述；未验证真实非同心Maxwell的式(5)余项与式(13)β；未证明event-aware调度最优；未获得开放支持的完整posterior mass reference；未建立主路线首创性或IEEE TAP接收资格。

这些限制不改变可立即执行的路线选择，但限定本报告能作出的科学承诺。
