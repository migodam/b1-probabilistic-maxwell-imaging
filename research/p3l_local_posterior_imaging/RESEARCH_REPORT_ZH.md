# P3-L：非线性 Maxwell 材料概率成像的实际实验与裁决

**最终裁决：REDIRECT。当前“弱密度 × strong Gaussian”成像路线 NO-GO；没有建立合格 FULL posterior，也没有得到可信 Bayesian 成像改善。**

128个事前冻结的新对象全部进入统计。近似材料均值的complex MSE为**0.0615021**，probit坐标点估计为**0.0110309**：MSE增加**457.5%**，配对差的95%区间**[0.0427873,0.0582865]**完全为正。原P1、噪声和Maxwell模型均未改变；没有通过放宽先验或筛选漂亮案例挽救结论。

这不反驳真正posterior mean的Bayes最优性。当前失败的是具体近似和有限采样；没有可靠FULL参考，就不能证明正确Bayes材料平均无益，也不能宣称问题已达到信息论极限。

## 1. 问题、物理与概率目标

本轮要回答：在物理重建已恢复敏感信息后，把剩余不确定性正确平均，是否能比单点重建得到更低材料风险？强方向的概率厚度与弱方向的分布必须共同保留，不能把优化后的中心面当完整后验。

沿用P1的64维标准Gaussian场坐标和2维独立uniform属性，经原C²饱和decoder形成24×24被动材料、256-cell ROI。8个照明共享同一材料，32个complex receivers；gain固定1。host不是空气，拓扑不是generator标签。Full仅指原生二维scalar-TM/VIE的完整离散求解，不是三维vector-Maxwell或连续介质认证。

噪声sigma=0.00119211798514519固定，Re/Im各方差sigma²/2，不随对象重标定。Phi=||F-y||²/sigma²=.5||white(F-y)||²；pi(q|y)正比exp(-Phi-||z||²/2)。每个source按Re32后Im32 packing。u=normalCDF(eta)使原Gaussian64×uniform2精确变成Gaussian66测度，保留正确Jacobian；这不允许在原uniform坐标里随意添加Gaussian Hessian。

noise-ball是反演资格诊断，不是Gaussian posterior硬支持。数学上合格的采样核、经验混合、物理相容与Bayes风险是四个不同结论。

## 2. 对象、隔离与基线

固定800个独立parent：32个development、512个joint-supervision train、128个validation和128个新locked。development内24/8用于几何与方法开发；网络的128validation是另外的独立集合。仅依赖y的anchor特征选择16个深参考条件，实际深入3个，其余未完成深采样。

数据生成、在线反演和离线评分分进程。所有locked模型、原始parent ID、预算、失败fallback和主比较在观测生成前冻结；预观察的报告合同修复另有记录。128份在线结果和Full-IS全部关闭后，独立evaluator才读取GT。后来这128个对象已经暴露，不能供下一轮调参后继续称新held-out。

32个development基础解全部进入原noise-ball，27个达到冻结KKT标准。128个新对象的native96有79个、probit64有115个达到各自驻点标准，所有失败与未驻点对象均保留。这里不是全局MAP证明。

probit坐标点的MSE比native点低11.71%，配对相对改善95%区间[10.57%,12.85%]。这是坐标与额外优化后的点估计收益，不是Bayesian平均或新增观测信息。MAP不是重参数化不变量，不能把两个点估计器混叫同一MAP。

## 3. Full posterior reference实际做到了什么

| 对象 / 方法 | 每链 production | 材料函数最大 Rhat | 最小 bulk ESS | 均值 RMS MCSE | 资格 |
| --- | ---: | ---: | ---: | ---: | --- |
| [DEV_000/mala](A/references/DEV_000/mala/REFERENCE_STATS.json) | 4096 | 3.20089 | 4.45452 | 0.0556548 | UNQUALIFIED |
| [DEV_028/mala](A/references/DEV_028/mala/REFERENCE_STATS.json) | 4096 | 5.39208 | 4.17595 | 0.105155 | UNQUALIFIED |
| [DEV_030/ess_mode](A/references/DEV_030/ess_mode/REFERENCE_STATS.json) | 4096 | UNASSESSED | UNASSESSED | 0.0114013 | UNQUALIFIED |
| [DEV_030/hmc_pilot](A/references/DEV_030/hmc_pilot/REFERENCE_STATS.json) | 1024 | 1.26968 | 11.4748 | 0.0101516 | UNQUALIFIED |
| [DEV_030/mala_mode](A/references/DEV_030/mala_mode/REFERENCE_STATS.json) | 4096 | 1.92985 | 5.55516 | 0.0182821 | UNQUALIFIED |
| [DEV_030/mala_retracted](A/references/DEV_030/mala_retracted/REFERENCE_STATS.json) | 4096 | 2.29 | 5.00098 | 0.0239059 | UNQUALIFIED |
| [DEV_030/pcn_mode](A/references/DEV_030/pcn_mode/REFERENCE_STATS.json) | 4096 | UNASSESSED | UNASSESSED | 0.010337 | UNQUALIFIED |
| [DEV_030/spine_hmc_pilot](B/references/DEV_030/spine_hmc_pilot/REFERENCE_STATS.json) | 1024 | 1.6054 | 6.85148 | 0.0108156 | UNQUALIFIED |
| [DEV_030/mala_decoder](C/references/DEV_030/mala_decoder/REFERENCE_STATS.json) | 4096 | 4.18571 | 4.26597 | 0.0230438 | UNQUALIFIED |
| [DEV_030/surrogate_force_linear](C/references/DEV_030/surrogate_force_linear/REFERENCE_STATS.json) | 4096 | UNASSESSED | UNASSESSED | 0.0166675 | UNQUALIFIED |
| [DEV_030/surrogate_force_material](C/references/DEV_030/surrogate_force_material/REFERENCE_STATS.json) | 4096 | 1.63156 | 6.66187 | 0.016892 | UNQUALIFIED |

全部11组均未达到原参考标准。共44条链、155,648个保留production状态，包含所有拒绝重复；这个数量不等于有效样本数。原生Elliptical Slice提供独立机制核验，也未合格。SMC的128粒子只剩1个祖先，不能作为teacher。

这说明目前的主要瓶颈仍是sampling/reference，而不是已经证明不存在材料不确定性。较稳定的Phi轨迹、较高接受率或无数值异常都不能替代材料函数混合。

## 4. 条件积分、材料矩与近似误差

对Gaussian64采用完整正交frame，属性显式保留。FULL弱marginal必须积分strong方向的概率体积；冻结是SLICE，优化是PROFILE，有限区域条件化是LOCAL-C，它们不是同一目标。

40组强条件实验覆盖4个对象、弱维数2/4/8/16/24及两个外层坐标。无一达到严格strong-gradient驻点标准；128点Full重要性核验中38组ESS<2，27组最大权重>.99。10组的局部线性梯度logZ修正上界<.01，仍无法消除非二次范围上的权重退化，因此不能把失败全部归因于中心未优化。

实现了原饱和多项式decoder在声明Gaussian条件下的解析一、二阶材料矩，经过独立quadrature/MC与tail/小方差检查。它准确积分supplied Gaussian，只消除内层材料统计费用；并不证明该Gaussian等于Maxwell posterior。

32个development对象中，固定属性Gaussian材料均值MSE .0206962，高于native点 .0138212。联合属性的probit Gaussian修订也运行GH5/9/17，GH9到17最大均值RMS变化2.8963e-6；联合材料均值MSE仍为 .0185700，高于相应坐标点 .0120131，四个Full校正ESS均约1。积分稳定与参考分布正确是两回事。

独立三维SLICE的阶数细化已收敛，但没有积分其余63个自由度。不能用这项成功给FULL posterior背书。

## 5. 独立128对象的主要成像结果

| 冻结估计器 | 128-parent complex MSE |
| --- | ---: |
| 16-state context 点估计 | 0.0128417288 |
| 独立 native96 点估计 | 0.0124938122 |
| probit64 坐标点估计 | 0.0110309173 |
| 弱密度 × strong Gaussian 材料均值 | 0.0615021168 |
| strong conditional-mode plugin 平均 | 0.0565169523 |
| 零观测 prior 材料均值 | 0.1754206811 |
| 有限 Full-SNIS，未合格 | 0.0240009114 |

主比较为材料均值减probit点：MSE差+0.0504712，95%区间[0.0427873,0.0582865]。conditional-mode plugin的差+0.0454860，[0.0376429,0.0537659]；有限Full-SNIS的差+0.0129700，[0.0106119,0.0156354]。三者均未改善主风险。

所有128对象、7个估计器均有输出，无raw失败或fallback，无truth缺失；不是删除困难对象后的统计。bootstrap以独立parent为单位、2000次、seed2026102299，比较总体平均风险之比，未平均每对象的百分比。独立审计重算所有均值与区间，数值差为零。

128个Full-SNIS的ESS最小/中位数均1、最大1.99706。权重公式正确没有使128-draw近似自动可靠，有限SNIS偏差仍在。不能将这列称合格posterior mean。

图像使用冻结前四个对象，共享色标。近似平均在小结构周围出现弥散激活；这不是新发现的高概率材料世界。[实部图](D/locked_evaluation/figures/fixed_first4_epsilon_real.png)、[损耗图](D/locked_evaluation/figures/fixed_first4_epsilon_imag.png)和[逐列说明](IMAGE_READING_GUIDE.md)可直接核验。

## 6. 学习密度、覆盖与失败

CNN+MLP、MLP diagonal/full Gaussian、MDN和小Flow各3个seed，实际完成15次joint-supervised训练。CNN/full-MLP约467k参数、规模差<.05%。全部在validation patience停止，最佳checkpoint均update0；训练NLL下降而validation恶化。独立审计没有找到已确认的batch、缩放或恢复错误，不能把初始化被选中称作训练收益。

随后运行48个简单Gaussian正则化对照，validation选择PCA32/ridge1000/shrink.5，NLL14.8138，优于无条件Gaussian15.1007。它是冻结128评估的实际密度，不是最后临时选择的神经模型。

新对象的weak8/probit2 proper NLL均值14.3094，parent SE .2091，分量90%coverage在85.94%–92.97%。这只是弱坐标的joint-simulation预测诊断，十个相关分量不是1280个独立对象，也不是材料credible coverage或FULL校准。接近名义的单分量coverage可以与严重错误的材料平均同时出现。

在4个development对象上，条件与无条件密度各256步Full-MH均0接受。归一化和稍好的低维NLL没有修复当前strong Gaussian参考的失配。没有不合格teacher distillation；RL、bandit或policy学习没有启动，因为缺少可信material-error奖励参考。

## 7. 科学修订：保留decoder的代理力

Gaussian平均失败后，我们没有继续增大网络，而是检验更具体的物理机制：用相同Full端点预算，比较latent-linear代理力与保留完整非线性decoder的冻结material-response代理力。每条内部轨迹不调用Full梯度；4次真实Jacobian setup照常计费。production代理、质量矩阵和步长冻结，端点使用真正Full Hamiltonian接受判定，目标不变。

两臂各4链，每链1024warmup+4096production，每臂164,384 RHS。材料代理把付费端点平均绝对Phi预测误差由8.27591降至.545205，白化响应RMS误差由.114001降至.028495，说明decoder保留确实改善了代理响应。

但材料均值MCSE没有改善：linear .0166675，material .0168920；包含audit的wall由245.671s升至393.749s。材料代理Rhat1.63156、最小bulk ESS6.66，仍未合格。linear有2个常数未评估分量，保留null。更大步长、更准确forward或更高接受率均没有转化成参考资格和成像增益，因此当前成本收益NO-GO。

普通Full-HMC pilot的warmup/长度不同，只作为实际成本背景，不冒称单独替换力的配对消融。通用代理力+真实Hamiltonian MH已有先例，本轮不将其正确性写成原创定理。[完整对照](C/SURROGATE_COMPARISON.md)。

## 8. 合法材料投影与几何

补做32个暴露DEV的合法reference投影，一次固定anchor起点、256decoder迭代上限、零物理/GT。全部32份原prior合法，实际2.089s。随后所有32份Full验证共256 RHS，仅1份进入noise-ball；独立离线MSE .0208401，仍高于原点 .0138212，未修复Gaussian平均问题。它是wrong/unqualified参考均值的局部投影，非全局restricted Bayes action；mean动作本身也不必自动满足noise-ball。

两个真实P1 anchor的投影场Lie-bracket有限差分已经运行，物理材料度量下normal leakage约.002064/.002739；冻结Maxwell、只改变decoder约.002021/.002680。该方向对、阈值、谱隙和步长均保存。它既不是一般Frobenius定理，也不证明存在exact material fiber。

## 9. Maxwell复用与A17

16个固定提案比较Full、exact T-matrix、GMRES与goal-ROM。相同current space的A17/Galerkin差≤2.02e-15；并无额外精度。rank16条件浮点界在5/16次可提前判断，11次fallback、观察到决定不一致0。总960 RHS包含576 RHS inverse setup与各128 RHS的T-matrix/Full/GMRES。

median Full .01303s，T-matrix .01098s，GMRES .009225s，但stability构造 .01476s、bounds+ROM约 .02128s，cold inverse另 .07649s。尚无端到端加速。这些界依赖IEEE/BLAS浮点模型，不是定向区间认证；未声称certified early rejection或A17独立优势。

## 10. 费用、内存恢复与封账

| 阶段 | 实付 complex RHS | state attempts | LU | 受控 wall 小时 | 剩余额度已封账 |
| --- | ---: | ---: | ---: | ---: | ---: |
| A | 5,542,136 | 514,404 | 514,404 | 1.8512 | 257,864 |
| B | 1,132,136 | 85,796 | 85,796 | 0.4037 | 1,367,864 |
| C | 659,584 | 61,552 | 61,536 | 0.3332 | 840,416 |
| D | 1,040,616 | 45,576 | 45,576 | 1.2165 | 459,384 |
| E | 0 | 0 | 0 | 0.0000 | 700,000 |

合计 **8,374,472 complex RHS、707,328 state attempts、707,312 LU**。受控wall **3.8046小时**，已记录physics wall **2.0024小时**，GPU训练0；new Mie/3D均0。没有外部付费服务。剩余 **3,625,528 RHS已封账**，不转移、不自动延续。

追加式账本覆盖失败、拒绝、forward、adjoint、J、LU和恢复。旧背景G/S及校准数组已付费共享，读取/审计另记CPU；此次增量费用不能冒称从零建立solver的cold部署成本。controller等待、封包和网络发布时间另列，不当作物理加速。

内存溢出来自transfer对旧LU的强引用与并发context worker。修复为有界4项缓存、单一物理worker，并保留中断前109,488 RHS和缺失尾部wall保守上界。3份已预约attempt缺少终点数组的历史缺口保留，未重跑掩盖。最后代理runner的重复metadata关键字在Book前报错，0物理；原stderr和源码保留，最小parser修复后沿用同一paid job ID与原授权，未清账。

## 11. 对核心问题的回答与下一步

**Probability：** 尚未更准确、可靠地计算原FULL posterior。11组参考不合格、SMC祖先与IS权重退化，不允许可信teacher、FULL credible coverage或escape mass声明。

**Imaging：** 当前近似概率成像明确失败，128对象主风险恶化；只有坐标点估计的11.71% MSE改善有独立统计支持，不能归因Bayesian平均。

**Computation：** 解析Gaussian材料矩保留为实用工具。保留decoder可降低代理误差，但未降低同RHS下的材料MCSE；A17/ROM没有总成本增量。现在没有等统计精度下更便宜的概率成像结论。

**具体裁决：** 停止把weak-density×single-Gaussian strong参考作为本P1的posterior成像器；停止当前训练设置与TM24 A17整合路线，不扩大NN或RL。继续价值只在可靠FULL参考、非Gaussian强条件与材料风险导向的推断机制；下一阶段需独立授权，不能花本次余额继续找正例。

**论文：** 当前不足以支撑强IEEE TAP算法论文。LIS、Rao–Blackwell、Laplace、transport、代理力HMC与goal-ROM有成熟先例；ScholarQA发现查询429和全文缺口均保留，不能据此宣称原创。缺少的决定性实验是：同一Full P1目标下可信独立reference、多个新对象的≥5%材料MSE改善及配对区间排除0，并证明Maxwell专用机制在完整成本上优于强传统方法。

本轮结论限于具体近似、采样和计算机制，不是Bayesian逆散射不可能，也不是原P1没有重要不确定性的证明。

[完整参考审计](POSTERIOR_REFERENCE_AUDIT.md) · [数学合同](THEORY_AND_ALGORITHM.md) · [风险解释](THEORY_RISK_INTERPRETATION.md) · [迭代](ITERATION_LOG.md) · [失败](FAILURE_ANALYSIS.md) · [复现](REPRODUCTION.md) · [账本](RESOURCE_LEDGER.json) · [原始风险](D/LOCKED_IMAGE_RISK.json)
