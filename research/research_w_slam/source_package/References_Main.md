# 参考文献：外部研究

以下 S 编号与主报告一致。书目核查截至 2026 年 10 月 24 日（UTC）；期刊卷年、在线发表日期和预印本日期分别标注，未核实的月日不补写。**F** 表示取得原始全文并检查相关模型或方法章节，并不表示逐句复核全部证明；**A** 表示核查原始摘要和出版元数据；**M** 表示作者或机构原始出版记录。与正文结论直接有关的访问限制列在相应条目下。项目内部证据另用 W/U 编号。

## [S1] iSAM2

Michael Kaess, Hordur Johannsson, Richard Roberts, Viorela Ila, John J. Leonard, Frank Dellaert. “iSAM2: Incremental Smoothing and Mapping Using the Bayes Tree.” *The International Journal of Robotics Research*, **31**(2):216–235, 2012.

DOI：[10.1177/0278364911430419](https://doi.org/10.1177/0278364911430419)。[作者原始 PDF](https://people.csail.mit.edu/kaess/pub/Kaess12ijrr.pdf)。

核查 **F**：重点检查选择性重线性化、部分状态更新和复杂度条件。其计算收益依赖图结构，不能直接推到消场后的全局 Maxwell 材料因子。

## [S2] Information-guided gating / selective partial optimization

Reza Arablouei. “Efficient Incremental SLAM via Information-Guided Gating and Selective Partial Optimization.” *Robotics*, **15**(5):87, **2026 年 4 月 27 日**。

DOI：[10.3390/robotics15050087](https://doi.org/10.3390/robotics15050087)。[期刊原始页](https://www.mdpi.com/2218-6581/15/5/87)。对应预印本标题为 “Efficient Incremental SLAM via Information-Guided and Selective Optimization”，[arXiv:2601.08110](https://arxiv.org/abs/2601.08110)，**2026 年 1 月 13 日**；[预印本全文](https://arxiv.org/html/2601.08110)。

核查 **F（预印本）/A（正式版）**。本报告对纯信息门控的反例是独立数学结果；预印本中信息量与 inexact-GN 残差之间的论证仍需额外核验，不能据此否定未经逐式审阅的期刊终版全部结论。

## [S3] NF-iSAM

Qiangqiang Huang, Can Pu, Kasra Khosoussi, David M. Rosen, Dehann Fourie, Jonathan P. How, John J. Leonard. “Incremental Non-Gaussian Inference for SLAM Using Normalizing Flows.” *IEEE Transactions on Robotics*, **39**(2):1458–1475, **2023 年 4 月**；在线发表年份为 2022。

DOI：[10.1109/TRO.2022.3216498](https://doi.org/10.1109/TRO.2022.3216498)。[IEEE 原始页](https://ieeexplore.ieee.org/document/9946412)。[arXiv:2110.00876v2 全文](https://arxiv.org/html/2110.00876v2)，修订于 **2022 年 7 月 2 日**。[IEEE 当期目录](https://www.ieee-ras.org/images/publications/t-ro/TRO.39.2.Table_of_Contents.pdf)。

核查 **F**：已研究非 Gaussian 联合后验、连续/离散不确定性及增量更新；其图结构带来的计算收益不自动适用于 dense Maxwell 因子。

## [S4] 多地图特征的 multipath BP-SLAM

Xuhong Li, Xuesong Cai, Erik Leitinger, Fredrik Tufvesson. “A Belief Propagation Algorithm for Multipath-based SLAM with Multiple Map Features: A mmWave MIMO Application.” **arXiv:2403.10095v1，2024 年 3 月 15 日**。

原始来源：[arXiv 全文](https://arxiv.org/html/2403.10095v1)；[作者机构出版记录](https://tugraz.elsevierpure.com/en/publications/a-belief-propagation-algorithm-for-multipath-based-slam-with-mult/)。

核查 **F**：采用路径参数和不同地图特征模型，并有 mmWave MIMO 实测；不是以体介电材料为未知量的复电场 Maxwell 反演。此处引用明确固定到已审核预印本版本。

## [S5] Radio-based SLAM 概率占据栅格

Xuhong Li, Erik Leitinger, Fredrik Tufvesson, Florian Meyer. “Probabilistic Occupancy Grid for Radio-Based SLAM.” **arXiv:2603.03559v1，2026 年 3 月 3 日**。

原始来源：[arXiv 全文](https://arxiv.org/html/2603.03559v1)。

核查 **F**：forward radio sensor model 显式考虑 cell dependencies，支持地图扩展和细化。物理模型是二维水平镜面路径、最高双反射、反射系数材料描述；不应误写为独立 voxel 模型，也不应等同于任意体材料 full-Maxwell posterior。

## [S6] Forward-model occupancy mapping

Sebastian Thrun. “Learning Occupancy Grid Maps With Forward Sensor Models.” *Autonomous Robots*, **15**:111–127, 2003.

DOI：[10.1023/A:1025584807625](https://doi.org/10.1023/A:1025584807625)。[作者原始稿](https://robots.stanford.edu/papers/thrun.occ-journal.pdf)。

核查 **F**：其核心包括处理格点间的统计依赖。不能把“只存独立 voxel 边缘”归为该文的已有结果或固有限制。

## [S7] 分层环境逆散射不确定性

Carolina Abugattas, Ana Carpio, Elena Cebrián, Gerardo Oleaga. “Quantifying uncertainty in inverse scattering problems set in layered environments.” *Applied Mathematics and Computation*, **500**:129453, **2025 年 9 月**。

DOI：[10.1016/j.amc.2025.129453](https://doi.org/10.1016/j.amc.2025.129453)。[出版社原始页](https://www.sciencedirect.com/science/article/pii/S0096300325001808)。[arXiv:2504.05776v1 全文](https://arxiv.org/html/2504.05776v1)。

核查 **F**：已知层状背景中的 inclusion 几何/材料推断，采用标量纵向 P-wave 模型，结合 Laplace 与 MCMC。它不是完整 Maxwell 未知墙–目标联合概率地图。

## [S8] Bayesian time-lapse FWI with HMC

P. D. S. de Lima, M. S. Ferreira, G. Corso, J. M. de Araújo. “Bayesian time-lapse full waveform inversion using Hamiltonian Monte Carlo.” *Geophysical Prospecting*, **72**(9):3381–3398, 2024；在线发表 **2024 年 9 月 8 日**。

DOI：[10.1111/1365-2478.13604](https://doi.org/10.1111/1365-2478.13604)。[出版社原始页](https://onlinelibrary.wiley.com/doi/full/10.1111/1365-2478.13604)。[arXiv:2311.02999v2 全文](https://arxiv.org/html/2311.02999v2)。

核查 **F**：实现将 baseline posterior 的 Gaussian 近似用于 monitor prior，比较序列与平行方案。物理材料变化与静态地图估计修正必须区别。

## [S9] 墙体和目标参数的 DNN 联合估计

Fardin Ghorbani, Hossein Soleimani. “Simultaneous Estimation of Wall and Object Parameters in TWR Using Deep Neural Network.” *International Journal of Antennas and Propagation*, **2022**:7810213, 10 pages，**2022 年 12 月 28 日**。

DOI：[10.1155/2022/7810213](https://doi.org/10.1155/2022/7810213)。[出版社原始页](https://onlinelibrary.wiley.com/doi/10.1155/2022/7810213)。对应 [arXiv:2111.04568](https://arxiv.org/abs/2111.04568) 的提交日期为 **2021 年 10 月 30 日**；[原始预印本 PDF](https://arxiv.org/pdf/2111.04568)。

核查 **F（2021 预印本）/A（正式版）**：本报告对模型维数和先验范围的描述以已检查的预印本为依据，未将其有限测试精度解释为 calibrated posterior 或未知拓扑普适成功；正式版未逐节复核。

## [S10] IEEE TAP：delay-only 墙参数估计

Pavel Protiva, Jan Mrkvica, Jan Macháč. “Estimation of Wall Parameters From Time-Delay-Only Through-Wall Radar Measurements.” *IEEE Transactions on Antennas and Propagation*, **59**(11):4268–4278, **2011 年 11 月**。

DOI：[10.1109/TAP.2011.2164206](https://doi.org/10.1109/TAP.2011.2164206)。[作者 PDF 地址](https://elmag.fel.cvut.cz/sites/default/files/users/machac/files/ap_2011-2.pdf)。[作者机构出版目录](https://fel.cvut.cz/en/faculty/people/874-jan-machac/publications)。

核查 **A/M**：已核原始摘要索引和出版记录，本次未能获取原始 PDF 全文。可确认其墙厚/介电率估计主题；不能据此假定它也完成目标 posterior 或复杂体材料联合推断。

## [S11] IEEE TIM：分布式 TWI 墙–目标联合估计

Yu Yao, Jiahui Chen, Shisheng Guo, Xiaojian Hao, Nian Li, Guolong Cui, Lingjiang Kong. “Target Localization and Wall Parameters Estimation via Distributed Through-Wall Imaging Radar.” *IEEE Transactions on Instrumentation and Measurement*, **74**:1–11, article **6500611**, 2025；在线 DOI 对应年份为 2024。

DOI：[10.1109/TIM.2024.3502875](https://doi.org/10.1109/TIM.2024.3502875)。[IEEE 原始页](https://ieeexplore.ieee.org/document/10758832)。

核查 **A**：原始摘要明确同时进行墙参数和目标定位，并利用分布式观测与目标一致性。全文获取受限，不对未检查的后验表示、迭代和回访细节作否定性断言。

## [S12] Full-wave 墙传播与 SAR 目标反射率估计

Francis M. Watson, Daniel Andre, W. R. B. Lionheart. “Resolving Full-Wave Through-Wall Transmission Effects in Multi-Static Synthetic Aperture Radar.” *Inverse Problems*, **40**(8):085009, **2024 年 8 月**；在线发表 **2024 年 7 月 9 日**。

DOI：[10.1088/1361-6420/ad5b83](https://doi.org/10.1088/1361-6420/ad5b83)。[原始预印本全文](https://arxiv.org/html/2403.10354v1)。[作者机构正式版记录](https://dspace.lib.cranfield.ac.uk/items/279aa6cc-1df0-458f-9d5d-a9fb1de82a09)。

核查 **F**：墙背景 full-wave BEM/Green、POD 与 variable projection；目标反射率侧采用 Born。墙体的全阶传播不等于目标–墙任意次数的完整双侧反馈。

## [S13] 可回访浅层的 stepped layer stripping

Nan Huai, Zhaofa Zeng, Jing Li, Yingwei Yan, Qi Lu. “Model-based layer stripping FWI with a stepped inversion sequence for GPR data.” *Geophysical Journal International*, **218**(2):1032–1043, **2019 年 8 月**；在线发表 **2019 年 5 月 8 日**。

DOI：[10.1093/gji/ggz210](https://doi.org/10.1093/gji/ggz210)。[出版社原始全文](https://academic.oup.com/gji/article/218/2/1032/5487078)。

核查 **F**：§2.5.3 的阶梯频率顺序明确重新访问此前层。该文是“layer stripping 不一定永久冻结”的具体电磁 FWI 先例。

## [S14] 任务导向后验低秩近似

Alessio Spantini, Tiangang Cui, Karen Willcox, Luis Tenorio, Youssef Marzouk. “Goal-Oriented Optimal Approximations of Bayesian Linear Inverse Problems.” *SIAM Journal on Scientific Computing*, **39**(5):S167–S196, 2017；在线发表 **2017 年 10 月 26 日**。

DOI：[10.1137/16M1082123](https://doi.org/10.1137/16M1082123)。[出版社原始页](https://epubs.siam.org/doi/10.1137/16M1082123)。[原始预印本 PDF](https://arxiv.org/pdf/1607.01881)。

核查 **F**：线性 Gaussian 框架下 QoI 后验协方差与均值的最优低秩近似、广义谱和加权 Bayes 风险。一般任务谱压缩不能作为本研究新发现。

## [S15] Goal-oriented optimal experimental design

Ahmed Attia, Alen Alexanderian, Arvind K. Saibaba. “Goal-Oriented Optimal Design of Experiments for Large-Scale Bayesian Linear Inverse Problems.” *Inverse Problems*, **34**(9):095009, **2018 年 9 月**。

DOI：[10.1088/1361-6420/aad210](https://doi.org/10.1088/1361-6420/aad210)。[arXiv:1802.06517v2](https://arxiv.org/abs/1802.06517v2)，修订于 **2018 年 6 月 11 日**；[原始全文 PDF](https://arxiv.org/pdf/1802.06517)。

核查 **F**：检查 §2–3 的线性前向、Gaussian prior/noise、线性 QoI，以及 A-GOODE/D-GOODE。结果不能无条件替换为非线性多模态 Maxwell 任意损失的 Bayes 风险设计。

## [S16] Risk-aware goal-oriented Bayesian OED

John D. Jakeman, Rebekah White, Bart van Bloemen Waanders, Drew P. Kouri, Alen Alexanderian. “Risk-Aware Goal-Oriented Bayesian Optimal Experimental Design.” **arXiv:2609.08023v1，2026 年 9 月 7 日**。

原始来源：[摘要与元数据](https://arxiv.org/abs/2609.08023)；[全文](https://arxiv.org/html/2609.08023v1)。

核查 **F（框架相关章节）**。该文是预印本；任务域、预测不确定性和未来数据的多层风险设计已有相邻工作，不能作为 Maxwell 系统的已验证性能。

## [S17] PDE 正演近似与 Bayesian posterior 误差

S. L. Cotter, M. Dashti, A. M. Stuart. “Approximation of Bayesian Inverse Problems for PDEs.” *SIAM Journal on Numerical Analysis*, **48**(1):322–345, 2010.

DOI：[10.1137/090770734](https://doi.org/10.1137/090770734)。[作者机构正式稿](https://eprints.maths.manchester.ac.uk/2211/1/090770734.pdf)。[arXiv:0909.2126 原始稿](https://arxiv.org/pdf/0909.2126)，首次提交 **2009 年 9 月 11 日**。[作者机构出版记录](https://research.manchester.ac.uk/en/publications/approximation-of-bayesian-inverse-problems-for-pdes/)。

核查 **F**：重点 Theorem 2.4、Corollary 2.5、Theorem 2.6；将具有适当可积性与稳定性控制的前向/势函数误差转化为 Hellinger 后验误差。所需统一条件不是任意 Maxwell 材料 prior 自动满足的性质。

## [S18] Sequential Monte Carlo for PDE inverse problems

Nikolas Kantas, Alexandros Beskos, Ajay Jasra. “Sequential Monte Carlo Methods for High-Dimensional Inverse Problems: A Case Study for the Navier–Stokes Equations.” *SIAM/ASA Journal on Uncertainty Quantification*, **2**(1):464–489, **2014 年 9 月**。

DOI：[10.1137/130930364](https://doi.org/10.1137/130930364)。[原始预印本 PDF](https://arxiv.org/pdf/1307.6127)。[作者机构出版记录](https://discovery.ucl.ac.uk/id/eprint/1425838/)。

核查 **F（问题与算法框架）**：序列 PDE 观测的高维 posterior 更新、SMC/MCMC/tempering 已有成熟先例；其相对成本结论不能直接转用为 W-SLAM 的加速证据。

## [S19] Flow matching posterior sampling with correction

Jan Blechschmidt, Oliver G. Ernst, Moritz Poguntke, Björn Sprungk. “Flow Matching for Fast Posterior Sampling in Bayesian Inverse Problems.” **arXiv:2610.02377v1，2026 年 10 月 1 日**。

原始来源：[摘要与元数据](https://arxiv.org/abs/2610.02377)；[全文 PDF](https://arxiv.org/pdf/2610.02377)。

核查 **F（引言与方法框架）**。近似 posterior 的 TV/KL 精度估计与 Metropolization 已有邻近工作；渐近校正不代表有限样本模式覆盖，也不意味着训练和校正成本为零。

## [S20] Dual-weighted residual 与任务误差控制

Roland Becker, Rolf Rannacher. “An optimal control approach to a posteriori error estimation in finite element methods.” *Acta Numerica*, **10**:1–102, **2001 年 5 月**；出版社在线日期 **2003 年 1 月 9 日**。

DOI：[10.1017/S0962492901000010](https://doi.org/10.1017/S0962492901000010)。[出版社原始页](https://www.cambridge.org/core/journals/acta-numerica/article/abs/an-optimal-control-approach-to-a-posteriori-error-estimation-in-finite-element-methods/5C67A03F528C6FA69F37A97DF5C3BE19)。

核查 **A**：原始摘要已明确局部残差、近似伴随权重、任务量误差和自适应计算。一般 dual-weighted residual 是继承技术；本研究需另外证明其在联合 posterior、有限 Maxwell 更新与可撤销冻结中的特定误差桥梁。

## [S21] 多 stand-off 观测破解未知墙歧义

Genyuan Wang, Moeness G. Amin. “Imaging Through Unknown Walls Using Different Standoff Distances.” *IEEE Transactions on Signal Processing*, **54**(10):4015–4025, **2006 年 10 月**。

DOI：[10.1109/TSP.2006.879325](https://doi.org/10.1109/TSP.2006.879325)。[作者原始 PDF](https://www1.villanova.edu/content/dam/villanova/engineering/research/BestPaper2010IEEETransactions.pdf)。

核查 **F**：多次改变 stand-off 并利用目标一致性降低墙–目标歧义已有直接物理先例。新的贡献必须落在更明确的 full-wave、共享状态和非线性条件上。

## [S22] 非线性 Bayesian 降维证书

Olivier Zahm, Tiangang Cui, Kody Law, Alessio Spantini, Youssef Marzouk. “Certified dimension reduction in nonlinear Bayesian inverse problems.” *Mathematics of Computation*, **91**:1789–1835, 2022.

DOI：[10.1090/mcom/3737](https://doi.org/10.1090/mcom/3737)。[原始预印本 PDF](https://arxiv.org/pdf/1807.03712)。

核查 **F**：KL 误差控制与梯度矩阵谱尾依赖相应分析假设；不自动覆盖离散拓扑模式、任意先验或共振附近的 Maxwell 参数族。

## [S23] Factor graphs for robot perception

Frank Dellaert, Michael Kaess. “Factor Graphs for Robot Perception.” *Foundations and Trends in Robotics*, **6**(1–2):1–139, **2017 年 8 月**。

DOI：[10.1561/2300000043](https://doi.org/10.1561/2300000043)。[作者原始出版页](https://www.cs.cmu.edu/~kaess/pub/Dellaert17fnt.html)。[作者原始书稿](https://dellaert.github.io/21S-8803MM/Readings/Dellaert17fnt.pdf)。

核查 **A/M**，定位原始书稿但未逐章审计。用于标准因子图和平滑框架的定位，不据此宣称 Maxwell 因子具有相同的稀疏复杂度。

## [S24] Information-theoretic active SLAM

Thomas Kollar, Nicholas Roy. “Efficient Optimization of Information-Theoretic Exploration in SLAM.” *Proceedings of the AAAI Conference on Artificial Intelligence*, 2008.

原始来源：[AAAI 出版页](https://www.aiinternational.org/Library/AAAI/2008/aaai08-217.php)；[作者团队出版目录](https://robustrobotics.github.io/publications.html)。

核查 **A/M**：已把信息目标、地图协方差与重访结合。此处不补写未核实的 DOI 或页码，也不将其几何 SLAM 信息准则直接等同于 Maxwell 任务 Bayes 风险。

