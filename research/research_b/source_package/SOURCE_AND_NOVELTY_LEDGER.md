# 来源、阅读层级与新颖性账本
检索与核对日期：2026-10-14。引文为来源，不构成本项目的性能证明。只使用论文、作者/机构页面、官方技术文档及用户研究仓库；检索摘要与实际全文阅读分开。

## 一、项目原始资料

| ID | 文件/路径 | 本次读取与适用范围 |
|---|---|---|
| P2-MAIN | [Joint P2 RESEARCH_REPORT_ZH](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/joint-p2-20261014/research/joint_p2/RESEARCH_REPORT_ZH.md) | 分段读取；径向 vector-Mie 共同后验、实验结论与科学边界 |
| P2-AMB | [AMBIGUITY_AND_COMPENSATION](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/joint-p2-20261014/research/joint_p2/AMBIGUITY_AND_COMPENSATION.md) | 已读；补偿见证、非全局性质、局部机制负结果 |
| P2-DOWN | [DOWNSTREAM_DECISION_TEST](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/joint-p2-20261014/research/joint_p2/DOWNSTREAM_DECISION_TEST.md) | 已读；独立几何 context、参考粒子、代表点不足 |
| P2-REF | [POSTERIOR_REFERENCE](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/joint-p2-20261014/research/joint_p2/POSTERIOR_REFERENCE.md) | 已读；proper shared gain、连续数值参考、full/proxy correction 与误差范围 |
| P2-COVER | [MODE_COVERAGE](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/joint-p2-20261014/research/joint_p2/MODE_COVERAGE.md) | 已读；事件与模式区别、remainder、漏事件与 ESS 反例 |
| P2-HANDOFF | [THEORY_HANDOFF](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/joint-p2-20261014/research/joint_p2/THEORY_HANDOFF.md) | 已读；事件内 likelihood 平均、TV/regret、latent/material KL 分解 |
| P2-STATUS | [THEORY_STATUS](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/joint-p2-20261014/research/joint_p2/THEORY_STATUS.md) | 已读；A17/A30 模块未接入 P2、未取得一般概率认证 |
| P2-GATE | [GATE_DECISION](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/joint-p2-20261014/research/joint_p2/GATE_DECISION.json) | 已读；3459 attempts、剩余 541 不属于 B 的预算 |
| C1-ENTRY | [C1 START_HERE](https://github.com/migodam/c1-physics-state-generative-imaging/blob/main/START_HERE.md) | 已读；当前 P2 与历史 P1 导航 |
| C1-P1 | [C1 P1 报告](https://github.com/migodam/c1-physics-state-generative-imaging/blob/main/p1/C1_P1_RESEARCH_REPORT_ZH.md) | 阅读主要报告段落；二维连续 latent 与 13-cell vector-DDA，不是一般形状模型 |
| C1-P2 | [C1 P2 报告](https://github.com/migodam/c1-physics-state-generative-imaging/blob/main/p2/C1_P2_RESEARCH_REPORT_ZH.md) | 阅读主要报告段落；最终事件质量与 conditional CDF、单频 C-target 限制 |
| A17-BASE | [A17 原研究报告](https://github.com/migodam/a17-maxwell-certified-exchange/blob/main/A17_FINAL_RESEARCH_REPORT_ZH.md) | 已读；原始局部 current exchange 与反馈材料机制 |
| A17-CLOSE | [A17 closure_final](https://github.com/migodam/a17-maxwell-certified-exchange/blob/main/closure_final/START_HERE.md) | 已读；覆盖原报告后续噪声与非线性事实，不能沿用过时“未运行”结论 |
| A17-ENTRY | [A17 START_HERE](https://github.com/migodam/a17-maxwell-certified-exchange/blob/main/START_HERE.md) | 已读；效率改进相对原交换路线，receiver-only 仍更快 |
| A23 | [A23 中文报告](https://github.com/migodam/a23-fast-material-imaging/blob/main/A23_RESEARCH_REPORT_ZH.md) | 阅读物理、公式、协议与失败部分；没有同质量加速证明 |
| A24 | [A24 中文报告](https://github.com/migodam/a24-source-anchored-om-feedback/blob/main/A24_RESEARCH_REPORT_ZH.md) | 阅读报告前部与裁决；低秩反馈不足、没有 3D/NN 扩张 |
| A24-ALG | [RELATIVE_FEEDBACK_OM_ARCHITECTURE](https://github.com/migodam/a24-source-anchored-om-feedback/blob/main/theory/RELATIVE_FEEDBACK_OM_ARCHITECTURE.md) | 阅读精确相减、反馈压缩、复杂度、伴随及适用条件 |
| A30 | [A30 E0/E1 DECISION](https://github.com/migodam/a30-maxwell-observable-material-manifold/blob/717d3d26a290248bd2b9aeb55f3d03694b407b84/A30_E0_E1_DECISION.md) | 阅读主要结果与补偿失败；24² scalar TM，不是一般 3D 证据 |

P2 读取时树 commit：`7ed4faef9bc3dd214d4de689c2203cc383f7495f`。其余 main 文件可能继续变化；Codex 应记录自己实际读取的 revision，不把文件名中的历史日期当不可变身份。本包不附带原始资产拷贝，也没有从标题推测任何本地路径。

## 二、最接近的文献

### R01 — 复合模型区分已有经典先例
A. C. Atkinson and V. V. Fedorov, “The design of experiments for discriminating between two rival models,” *Biometrika* 62(1), 57–70, 1975.
DOI: https://doi.org/10.1093/biomet/62.1.57
另见多模型版本：https://doi.org/10.1093/biomet/62.2.289

读取层级：出版社摘要/元数据；没有声称本次完整阅读付费原文。原文已明确顺序、非顺序及非线性参数模型区分。现代 T-optimal/model-discrimination 文献已覆盖对 rival parameters 再拟合这一基本思想。因此“测量后允许重新补偿再区分”不能单独称为新数学理论。

### R02 — EM inverse scattering 的 adaptive sensing 不是空白
A. E. Fouda and F. L. Teixeira, “Bayesian compressive sensing for ultrawideband inverse scattering in random media,” *Inverse Problems* 30, 114017, **2014**.
DOI: https://doi.org/10.1088/0266-5611/30/11/114017
作者稿：https://arxiv.org/abs/1401.1092

读取层级：公开 PDF 文本相关章节，特别是 adaptive sensing；不是只看标题。数据库有 preprint 日期/版本日期，本文使用 journal year 2014。框架包含 Born 线性化、Bayesian compressive sensing、以 posterior covariance/信息选择新增传感器、distorted-Born 高对比扩展与 time reversal。不能说旧方法都只做固定采集。

### R03 — 已有 full nonlinear EM VIE + Bayesian design
A. I. Sandhu, B. M. Dia, O. Dorn, P. Soupios, “Bayesian Experimental Design for Efficient Sensor Placement in Two-Dimensional Electromagnetic Imaging,” *IEEE Access* 11, 65649–65662, 2023.
DOI: https://doi.org/10.1109/ACCESS.2023.3288997
机构页：https://pure.kfupm.edu.sa/en/publications/bayesian-experimental-design-for-efficient-sensor-placement-in-tw/

读取层级：Undermind 全文阅读，核对 §II-A、II-C、III。2D scalar TM、full nonlinear volume-integral forward，GMRES；不是纯 Born。实际 cylinder radius/材料固定、推断中心位置。EIG 采用单峰 MAP 附近 Laplace，比较传感器/频率配置；没有建立材料多分支、共享校准、完整顺序闭环或 inclusive runtime 优势。本次没有复制其表格数据。

### R04 — 非线性 goal-oriented EIG 已存在
S. Zhong, W. Shen, T. A. Catanach, X. Huan, “Goal-Oriented Bayesian Optimal Experimental Design for Nonlinear Models using Markov Chain Monte Carlo.”
DOI: https://doi.org/10.1137/24M1649344
2024 preprint；此处不把 preprint 年份自动当 journal 发表年份。

读取层级：Undermind 全文。Eqs.10–14 为 QoI EIG，Eq.17 和 Algorithm1 为 nested MC/MCMC/KDE。它不是直接 Bayes 最优决策损失。§3.2.1 内层 MCMC 使用生成合成观测的参数作初始化，作者指出不可辨识/多峰时局限；B 的在线模块不得读取评价 truth 来复制这种加速。案例是非 EM 的 source inversion。

### R05 — nuisance-aware BED 已存在
A. Bartuska, L. Espath, R. Tempone, “Small-noise approximation for Bayesian optimal experimental design with nuisance uncertainty,” *CMAME* 399, 115320, **2022**.
DOI: https://doi.org/10.1016/j.cma.2022.115320
作者稿：https://arxiv.org/abs/2112.06794
机构页：https://nottingham-repository.worktribe.com/output/9400253

读取层级：作者摘要/期刊元数据。针对 nuisance 引入的内层积分做 small-noise 近似，可结合 Laplace；包括复合材料 EIT。不能把“对 nuisance 边缘化后选测”称为 B 独有。

### R06 — prior/weight ambiguity 的 robust EIG 已存在
J. Go and T. Isaac, “Robust Expected Information Gain for Optimal Bayesian Experimental Design Using Ambiguity Sets,” 2022 preprint.
https://arxiv.org/abs/2205.09914

读取层级：作者摘要与网页内容。基于 KL 邻域的 prior ambiguity set，构造 EIG 的 affine relaxation 和稳定化近似。本报告不沿用未全文核对的具体最优性条件。

### R07 — risk-aware goal-oriented design（2026 新近先例）
J. Jakeman, R. White, B. G. van Bloemen Waanders, D. P. Kouri, A. Alexanderian, “Risk-Aware Goal-Oriented Bayesian Optimal Experimental Design,” 2026-09-07 preprint.
https://arxiv.org/abs/2609.08023

读取层级：Undermind 全文。§3.2 Eq.13 组合 posterior deviation、QoI-domain risk 与 dataset risk；Gaussian weighted likelihood，continuous simplex designs；advection-diffusion 案例。没有 Maxwell 分支/共享仪器补偿的实证。并非所有 risk criterion 都等价于本报告的 posterior optimal material decision loss。

### R08 — 直接 decision-focused BOED（2026 先例）
J. Go, X.-N. Qian, B.-J. Yoon, “Goal-driven Bayesian Optimal Experimental Design for Robust Decision-Making Under Model Uncertainty,” 2026-05-25 preprint.
https://arxiv.org/abs/2605.26093

读取层级：公开 HTML 相关方法/理论。摊销 variational posterior + differentiable convex decision layer，目标直接对接决策；理论讨论无关参数方向的梯度。案例 source localization、epidemic、pharmacokinetic，不是 Maxwell compensation。直接优化下游风险也不是 B 首创。

### R09 — offline-online、低秩与目标 EIG 误差已有先例
K. Wu, P. Chen, O. Ghattas, “An Offline-Online Decomposition Method for Efficient Linear Bayesian Goal-Oriented Optimal Experimental Design: Application to Optimal Sensor Placement,” *SIAM J. Sci. Comput.*, 2023.
DOI: https://doi.org/10.1137/21M1466542

读取层级：出版社详细摘要、元数据。线性 Bayesian GOOED，两个 Hessian-related 算子的低秩近似、offline-online decomposition、目标 EIG 误差分析、swapping greedy。不能只用“离线 low rank + online no PDE”申报新颖性；Maxwell 非线性候选依赖和稳定常数才需检验。

### R10 — BAE model-error / nuisance-aware OED
K. Koval and R. Nicholson, “Non-intrusive optimal experimental design for large-scale nonlinear Bayesian inverse problems using a Bayesian approximation error approach,” *Journal of Scientific Computing*, journal DOI 2025.
DOI: https://doi.org/10.1007/s10915-025-03008-7

读取层级：Undermind 全文。Theorem1 与 Corollaries1/3 在完整 parameter-error cross covariance 的 BAE 近似中给出 linearization invariance，实际统计需 MC 收敛。允许零线性 surrogate，做 marginal A-optimal nuisance design；Darcy 与 tsunami 案例。**这不是 nonlinear full posterior action ranking 的 uniform certificate。**

### R11 — PDE nonlinear robust design
A. Chowdhary, A. Attia, A. Alexanderian, “Robust Optimal Experimental Design of Infinite-Dimensional Bayesian Nonlinear Inverse Problems.”
DOI: https://doi.org/10.1137/24M1693921

读取层级：文献摘要/元数据；preprint 与 journal 年份不混写。max-min、uncertain prior/model/noise、可扩展 EIG 近似均有相关框架。

### R12 — 观测相关性不能丢掉
A. Attia and E. Constantinescu, “Optimal Experimental Design for Inverse Problems in the Presence of Observation Correlations,” 2022.
DOI: https://doi.org/10.1137/21M1418666

读取层级：出版社元数据/摘要。用于说明 observation correlation-aware OED 已有基础；本报告的联合高斯条件公式为标准推导。

### R13 — 顺序假设检验及停止
M. Naghshvar and T. Javidi, “Active sequential hypothesis testing,” *Annals of Statistics*, 2013.
DOI: https://doi.org/10.1214/13-AOS1144
https://arxiv.org/abs/1203.4626

读取层级：作者摘要/元数据。主动选择行动、序贯停止、误判代价和 dynamic programming 已有理论；不把 stopping rule 当新概念。

### R14 — 异质成本顺序设计
G. Vershinin, A. Cohen, O. Gurewitz, “Active Sequential Hypothesis Testing with Non-Homogeneous Costs,” 2025-09-15 preprint.
https://arxiv.org/abs/2509.11632

读取层级：作者摘要。异质成本下比值及 policy 处理已存在；本报告不使用其未核对的具体渐近定理。

### R15 — BOED 稳定性与非 KL 准则（2026）
D. Wu, L. Liang, H. Yang, “Beyond Expected Information Gain: Stable Bayesian Optimal Experimental Design with Integral Probability Metrics and Plug-and-Play Extensions,” 2026-04-23 preprint.
https://arxiv.org/abs/2604.21849
https://arxiv.org/html/2604.21849v1

读取层级：HTML 相关 theory/assumptions。bounded-kernel MMD、Wasserstein、energy distance 在不同 tail/smoothness 条件下给 utility perturbation/stability；部分高维 plug-and-play 实验明确超出理论覆盖。它支持研究稳定效用，但不自动给 B 的材料 MAE 与 Maxwell 误差常数。

### R16 — EIG 与 SNR/Fisher 界的极新近先例
M. Doumbouya, A. Vidard, O. Zahm, “Certified Bayesian optimal experimental design: from Expected Information Gain to Signal-to-Noise Ratio,” 2026-09-26 preprint.
https://arxiv.org/abs/2609.32282

读取层级：**本次仅检索到作者摘要；PDF/HTML 获取失败。** 摘要报告基于 dimensional logarithmic Sobolev inequalities 的上下界、covariance/Fisher 结构及 incremental sequential bounds。不能基于摘要声称其理论适用于多峰 Maxwell 或反过来声称 B 首次有 certified utility。正式写论文前必须全文核对。

### R17 — control-oriented 设计
M. Madhavan 等, “A Control-Oriented Approach to Optimal Sensor Placement,” *SIAM/ASA J. Uncertainty Quantification* 14(2),711–741, 2026-06-18.
DOI: https://doi.org/10.1137/25M1733537

读取层级：出版社摘要/元数据。线性 PDE Bayesian inverse + control-oriented sensor placement；已有将设计对接 downstream control objective 的方法。

### R18 — 色散模型官方技术说明
COMSOL, “Fitting the Debye Dispersion Model to Experimental Data.”
https://www.comsol.com/support/learning-center/article/fitting-the-debye-dispersion-model-to-experimental-data-78171

读取层级：官方技术文档。用于因果/频域材料模型背景，不把软件文档视为研究新颖性来源。本报告明确采用自己的 \(e^{-i\omega t}\) 约定，正负虚部须随约定核对。

## 三、检索范围与限制

关键词覆盖 Maxwell/electromagnetic inverse scattering、adaptive sensing、sensor placement、Bayesian OED、goal-/decision-/risk-oriented design、nuisance uncertainty、model discrimination/T-optimality、sequential hypothesis testing、robust EIG、ROM/offline-online design、surrogate utility error 与 passivity/dispersion。使用 web 和 Undermind 的语义检索及所列 primary sources。

本次不是全学科系统综述，也没有穷尽所有 IEEE TAP/adaptive microwave imaging 文献。四篇关键论文通过全文阅读工具核对，其余层级如上；网页 PDF 截图服务返回错误，未据此宣称成功进行截图视觉审查，也没有抄录未核对图表的数值。B 的潜在创新被定位为待实验验证的组合机制，而非“全球首次”。

## 四、不得写入论文的未证实主张

“最优新频率是 2.4”；“所有角度对球都无用”；“两种极化总比一种获得不同结构信息”；“ROM 已证明加速”；“高 ESS 证明没有漏模式”；“EIG 减少保证每个样本 MAE 减少”；“归一化 proposal 等于 posterior”；“best-found 补偿距离大证明不存在更近解”；“全部三维 Maxwell 材料都可由有限菜单区分”；“本文首次提出任务风险/复合假设/稳健信息设计”。
