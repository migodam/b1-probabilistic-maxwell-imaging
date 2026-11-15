# 双度量与有限交替动力学：最邻近文献证据

2026-11-15。这里只提取先例和适用条件，最终新颖性、公式与实验解释由主线程裁决。两条ScholarQA检索受429限制、没有候选返回；因此本报告是固定种子原文核查，不能支持“未找到即不存在”。完整回执与原文见SEARCH_LOG及fulltext目录。

## 可直接核实的四项区别

1. **消去辅助状态后出现材料相关残差度量已有直接PDE先例。** van Leeuwen–Yang的§2.3–2.4给representer与reduced formulation；§3.3.3 Eq23采用随c变化Gram，Eq24明确改为数据估计的固定Gram。作者指出后者只需对data residual求导。两种目标在文中被明确分开；不能把逐迭代数值冻结可变Gram说成原可变目标的完整导数实现。这里未核其结果与当前Maxwell算子的逐式同构。[原文](https://arxiv.org/html/2403.15292v2)

2. **测量域耦合状态缺陷与数据缺陷已有直接CSI先例。** CC-CSI的Eq11–12增加cross-correlated residual；Eq14–20和Alg1使用顺序current/material PR方向及一维步长。Eq20在material line search中保留随试探材料变化的归一化分母；Alg1 step7又在current更新前刷新该量。必须区分目标、方向公式、line search和冻结时机，不能只看到两个块目标名称不同就判定非势，也不能把该实现称为完整精确block最佳响应。Eq15/18显式依赖旧gradient/direction；其实际一步映射含memory。[原文](https://arxiv.org/html/1906.10864v1)

3. **多损失的势/旋转分解是梯度场先例。** Balduzzi等§2定义simultaneous gradient和无穷小动力学，Lemma1分解其Jacobian为对称/反对称部分。它提供分析多目标耦合的语言；并未在该文证明当前顺序精确best-response、有限PR、stage/restart及box程序的谱或稳定性。把它当成非势直觉的背景引用可以；把其SGA保证直接移植给实际Maxwell迭代，没有本文证据。[原文](https://proceedings.mlr.press/v80/balduzzi18a/balduzzi18a.pdf)

4. **有限内解不必等于没有收敛理论，但必须核具体假设。** Attouch–Bolte–Svaiter §2.3 H1/H2/H3分别要求共同函数充分下降、相对subgradient误差、连续聚点条件，Theorem2.9再用KL性质。§6 Alg4/Theorem6.2覆盖受控误差的proximal Gauss-Seidel，度量特征值有统一正界。Remark6.1明确零proximal项的经典Gauss-Seidel不在该框架内。故“每块各自目标下降”“只做固定PR次数”“两个metric都是SPD”均不是这篇定理的完整验收条件。[作者原文](https://bolte.perso.math.cnrs.fr/MPA.pdf)

## 对当前对照命名的证据限制

以下是从原文条件得出的**核对项目**，不是新定理：

- 分别记录目标中的材料相关空间metric、源metric、归一化权重，以及各自更新/冻结的时机；“物理模型相同”不足以说明优化目标相同。
- 区分完整共同目标、冻结metric代理目标、预条件方向、非匹配block目标；这四个标签应由实际公式和程序确定。
- 固定次数PR与exact best response比较时，保留gradient/direction、sourcewise/global系数、line-search、stage/restart及box的状态范围。CC-CSI原文已显示类似memory，不能把它当当前项目独有结构。
- Hager–Zhang的§1 Eq1.2–1.3/Table1.1进一步给出非线性CG/PR递推，说明新方向依赖前方向与gradient差。[作者PDF](https://people.clas.ufl.edu/hager/files/cg_survey.pdf) 这只是算法身份背景，不是当前程序任何收敛证明。

## 已核实书目

- Tristan van Leeuwen, Yunan Yang. *An analysis of constraint-relaxation in PDE-based inverse problems*. Inverse Problems41(2),025009,2025. DOI[10.1088/1361-6420/adab89](https://doi.org/10.1088/1361-6420/adab89). arXiv2024；S2year2024不应误写期刊年。原文HTML生成日期显示2026，投稿元数据仍为v2/2024-11-29。
- Shilong Sun, Bert Jan Kooij, Tian Jin, Alexander G. Yarovoy. *Cross-Correlated Contrast Source Inversion*. IEEE TAP65(5),2592–2603,2017. DOI[10.1109/TAP.2017.2673758](https://doi.org/10.1109/TAP.2017.2673758). arXiv2019。
- David Balduzzi, Sebastien Racaniere, James Martens, Jakob Foerster, Karl Tuyls, Thore Graepel. *The Mechanics of n-Player Differentiable Games*. ICML2018,PMLR80:354–363. [官方书目](https://proceedings.mlr.press/v80/balduzzi18a.html)。
- Hedy Attouch, Jérôme Bolte, Benar Fux Svaiter. *Convergence of descent methods for semi-algebraic and tame problems: proximal algorithms, forward–backward splitting, and regularized Gauss–Seidel methods*. Mathematical Programming137(1–2),91–129,2013；online2011. DOI[10.1007/s10107-011-0484-9](https://doi.org/10.1007/s10107-011-0484-9)。
- William W. Hager, Hongchao Zhang. *A survey of nonlinear conjugate gradient methods*. Pacific Journal of Optimization2,35–58,2006. [作者出版目录](https://people.clas.ufl.edu/hager/gradient-methods-and-software/)。

本轮无实验、无GT、无新公式有效性判断、无全局新颖性裁决。未检索到的文献不作为负证据；未从这些论文中建立current/material四组合对照的恢复预测。更完整综述应在检索恢复后补直接alternating-best-response、nonlinear Gauss-Seidel局部谱与variable projection原始文献，并逐式核对现有项目历史。
