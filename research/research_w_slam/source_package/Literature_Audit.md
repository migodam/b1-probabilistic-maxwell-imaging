# 文献与原创性审计附录

**Research W-SLAM — 从子空间可观测性到可回访的概率 Maxwell 地图**

检索截止：2026 年 10 月 24 日（UTC）。本附录为主报告提供文献证据、方法边界及原创性比较；不构成对一般未知材料类的算法认证。附录使用独立的 L01–L27 来源编号，主报告的 S 编号另行编排。

## 1. 审计结论

最有力的原创路线不是“把 SLAM 引入 Maxwell”，也不是“概率地图可以回访”。这些能力分别已由因子图平滑、iSAM2、非高斯增量 SLAM、占据地图、RF-SLAM 和序列 Bayesian PDE 反演提供。当前真正需要证明的是：**在消去场后非局部耦合、墙–目标有限幅度补偿和 PDE 成本占主导的条件下，任务误差证书能否决定哪些计算可以延迟，以及哪些新观测能消除跨观测一致的物理补偿。**

尤其需要补入两条需要补充的强邻近工作：

1. Arablouei 2026 已提出 information-guided gating + selective partial optimization（IGG/SPO），用信息量决定局部或全局增量 SLAM 更新，保留全部观测。新方法不能只对比 iSAM2。
2. Watson–Andre–Lionheart 2024 已联合未知墙体 nuisance 参数与目标反射率，使用墙体 full-wave BEM/POD 和 variable projection。其目标侧采用 Born 近似，因而“全波墙体补偿”不是本研究的新概念，但双侧所有阶相互散射、后验分支及任务更新证书仍有明确物理差异。

严谨性要求：以下“未解决”均指本文实际建模范围或本次已核查章节，并非声称整个学界不存在该结果。全文可访问不等于逐句复核所有证明。访问等级 F=原始全文，已检查相关方法/模型章节；A=原始摘要/出版社元数据，技术细节仍待全文；M=作者/机构原始出版目录或原始书目，仅核实出处。预印本和期刊版本分开，不能将 v1 的证明问题直接套给终版。

## 2. Originality matrix（覆盖任务要求的 14 类）

|已有方法 / 主要来源|已知结果|本次核查可确认的边界|W-SLAM 可能超越之处|必须做的实验或理论比较|
|---|---|---|---|---|
|Incremental smoothing / iSAM2 [L01]|Bayes tree、fluid/selective relinearization、局部变量更新、循环闭合后历史修正|复杂度收益依赖条件独立与局部稀疏；没有 Maxwell PDE 端到端成本定理|任务损失证书控制延迟更新，且统计状态始终保留|同一后验目标下 iSAM-like 阈值重线性化、周期 batch、证书策略；实际 PDE 次数和误差|
|信息门控增量 SLAM，IGG/SPO 2026 [L02]|logdet 信息门控、受影响集合内多次 GN、保留全部观测|信息/Hessian 增量本身不依赖 innovation；不能单独推出任务修正界|以实际 stationarity/dual residual + nonlinear remainder 取代纯信息阈值|加入大 innovation / 小 Fisher 增量反例；测任务误差与计算量|
|Factor graphs / map–pose covariance [L03]|联合图推断、Schur elimination、协方差与稀疏增量计算|PDE 网格稀疏不蕴含消场后材料信息矩阵稀疏|可验证的低秩全局消息近似及误差/成本界|材料 dense Schur 与 latent-field saddle system 对照；不能假称 Bayes tree 保持普通 SLAM 复杂度|
|Active SLAM / belief-space exploration [L04]|用地图不确定性规划观测轨迹，重复访问旧区域可以有价值|经典 SLAM 信息目标与电磁任务损失不同|在合法频率/孔径合同下，真 Maxwell 任务 Bayes 风险收益|风险 OED、D-opt、angle-only、random/acquisition-cost matched|
|相关 occupancy mapping [L05]|forward sensor model 已保留格点依赖；独立 cell posterior 是近似|该工作用已知位姿与几何距离传感器，非介电材料|占据、连续被动材料及墙–目标模式相关性的联合物理后验|独立边缘 vs correlated Gaussian vs mixture vs reference posterior|
|Multipath BP-SLAM 2024 [L06]|多类地图特征、虚拟锚与点散射体切换、联合位置姿态、真实 mmWave 数据|SAGE 提取路径参数的几何/射线模型，非复杂场材料 Maxwell likelihood|不依赖路径提取的复场联合材料–目标反演|公平说明模型不同；不能宣称首次 RF SLAM 或首次未知墙联合定位|
|Radio occupancy SLAM 2026 [L07]|层次化反射表面–占据栅格、相关格点、forward radio sensor model、先验地图修正|二维水平镜面路径至双反射，墙为反射系数；不是体介电率/损耗全波推断|全局多次散射材料图、非线性补偿与任务风险证书|至少承认其 probability mapping/revisit 已成立；做射线近似失效区间测试|
|TWR DNN joint wall/object [L08]|墙厚/介电率与一个或两个目标位置/材料联合预测|二维 TMz FDTD 数据；目标数预定；网络回归非校准后验|未知合法拓扑、多模态概率、measurement-only 锚点|目标数/形状先验外泛化、posterior calibration，而非只比图像|
|Unknown-wall multiple-stand-off TWI [L09]|不同 stand-off 观测可通过一致目标位置约束破解墙参数歧义|简化墙参数/传播与成像模型|所有阶 full-wave 的跨观测补偿一致性条件、有限幅度可实现见证|同单次距离、不同 stand-off、真实重新优化各分支；避免把矩阵堆叠当新发现|
|TAP delay-only wall estimation [L10]|多 Tx/Rx 间距延迟估计墙厚与介电率|低维墙参数；本次未能获取全文，不能臆断所有假设|空间变化复杂材料及联合任务风险|作为传统墙参数估计基线；不制造“从未联合估计”叙述|
|TIM distributed wall/target estimation [L11]|分布式 TWI 同时估计墙参数和目标位置；目标一致性用于墙补偿|仅核查原始摘要；全文细节待核验|宽材料类 posterior、多次散射、风险证书，需具体对比|获取全文后实现/匹配其假设基线；不能用缺少摘要描述来断言无回访|
|Indoor joint tracking/wall estimation [L12]|多径鬼影识别、目标跟踪与未知墙估计、实测室内雷达|当前为原始摘要级证据|体材料、运动目标–墙–pose–仪器相关性和 full-wave likelihood|不同运动先验与真实墙变化的识别实验|
|Full-wave TWI wall nuisance + Born target [L13]|任意形状近似已知墙体的 full-wave BEM 修正 Green 函数，POD，variable projection|目标 Born、反射率单次作用；墙所有阶不等于目标–墙全部双侧反馈|完整耦合 Maxwell + posterior branches + 计算收益证书|必备 ordinary reduced-model / variable projection 基线|
|Known-wall nonlinear microwave TWI [L14]|非均匀有限墙作为已知散射体，LM 非线性目标反演|墙已知，不是未知墙 posterior|未知墙–目标联合，但成熟 FWI 依然是强基线|known-wall oracle 只能诊断，不得作为可实现主方法|
|Bayesian inverse scattering in layered media [L15]|复杂 inclusion 的材料/几何 UQ，Laplace 与 MCMC，存在附加模式|已知层状背景，标量纵波力学模型；不是 Maxwell 未知墙|完整 Maxwell/宽墙材料图/跨时间相关后验|小维参考 posterior 与网格/数值误差验证应继承|
|Bayesian time-lapse FWI [L16]|baseline posterior 作为 monitor prior、HMC、时间相关信息复用|实现用 baseline Gaussian 近似；物理变化与静态地图修正须分开|静态 W + 动态 X/pose 的联合 smoothing，保留历史相关与多模式|parallel、sequential、batch joint；避免 prior+旧 likelihood 双重计数|
|Layer stripping / revisits [L17]|空间 window、频率 continuation、step sequence 可以重新更新浅层|“layer stripping”不是永久冻结的同义词|任务选择而非预定深度/频率顺序，并给风险/成本证书|hard freeze 仅一条人为 ablation；必须对比 iterative/stepped stripping|
|Principal-angle information geometry [L18]|CRB 与敏感性子空间主角度的联系已有经典文献|广义谱恒等式不能单独 claim new；有 prior 需改定义|跨时共享状态约束 + Maxwell 可实现性 + 有限非线性辨识|线性 toy 仅定理单元测试；需 full-wave 实例和重新优化|
|Goal-oriented posterior spectral approximation [L19]|任务量的最优低秩后验协方差/均值近似、广义谱结构|Gaussian linear setting，有具体误差度量|非线性多模态传播任务商空间及可撤销计算证书|同概率精度下同成本 comparison；不能只展示“task eigenspace”|
|Certified nonlinear dimension reduction [L20]|在 log-Sobolev 等假设下，梯度二阶矩谱尾可界 likelihood ridge 近似的 KL|证书条件与梯度期望估计有成本；不自动适配离散拓扑模式|Maxwell 专用、跨分支任务风险证书/后验反馈|报告 certificate overhead 与失效；不能把近似 KL 当 exact posterior|
|Risk-aware / goal-oriented OED [L21]|后验任务量、风险、未来数据多层期望的实验设计已有框架|新预印本不等于成熟已复现方案；非 Maxwell|Maxwell 模式重新优化后仍可分的采集设计；采集与计算分开|exact lowdim EVSI reference + nuisance marginalization + nested new data|
|Sequential infinite-dimensional Bayesian PDE OED [L22]|序列实验设计、MDP、参数/状态降维、可微 surrogate、策略优化|其 Laplace/eigen 信息奖励不等于本任务 Bayes loss|主动观测联合模式风险、独立全波校正与成本透明|强 greedy / myopic / direct posterior 基线足够；无须先引入 RL|
|Incremental nonGaussian SLAM / NF-iSAM [L23]|Bayes-tree 非高斯联合密度、continuous/discrete uncertainty、normalizing flows|计算收益依赖图结构；近似 posterior 需评估|非局部 PDE 下可校正 posterior 与风险成本证书|nonGaussian smoother 或小问题 direct posterior；mixture 本身非新|
|Physics generative inverse scattering [L24]|score prior + 材料/电流表示 + Lippmann–Schwinger data/state loss|ScoreField 是单实例优化，不应当作 calibrated sampler|合法多分支测度、corrected proposal、在线物理回访|普通 proposal+MH、unlearned priors、相同 Maxwell budget|
|Flow posterior sampler + exact correction [L25]|条件 flow matching、可算密度诊断、independence-MH 渐近校正|训练覆盖与误差测量有成本；非 Maxwell online mapping|面向特定 Maxwell 模式和任务的测度/成本收益|无训练 direct MCMC/SMC、MH-corrected vs uncorrected proposal|
|Goal-oriented dual-weighted residual [L26]|对偶加权局部 residual 估计目标量误差并驱动自适应计算，早已成熟|FEM 离散误差证书不自动等于 posterior approximation/branch risk 证书|有限 Maxwell inverse update + 联合 posterior + 延迟计算的特定可算界|必须把标准 DWR 步骤标成继承；只claim新增的统计/物理耦合部分|
|Sequential Bayesian PDE inverse inference [L27]|SMC 按新增观测推进高维 PDE posterior，tempering/MCMC adaptation/parallelization|不是“SLAM 首次带来 PDE 概率时序更新”；并未声称任意 full-wave 成本获益|实际宽材料模式覆盖、可回访计算预算和真实性能比较|direct SMC 是强基线，source-progressive 更名不构成进步|

## 3. 不能再使用的弱原创性表述

- “首次将 SLAM 引入电磁/射频地图”“首次同时估计墙和目标”：已有明确反例 [L06–L14]。
- “先前 layer stripping 都冻结浅层”：Huai et al. 的 stepped sequence 明确在深层开始后返回浅层 [L17]。
- “我们第一次让历史地图再次优化”：平滑与 loop closure 已有 [L01–L03]。
- “我们用信息量决定局部还是全局更新”：IGG/SPO 已直接提出 [L02]。
- “我们不独立存每个 voxel，而保存关联”：forward occupancy map 和 radio occupancy 已有 [L05,L07]。
- “Bayesian + mixture/flow + factor graph 是新方法”：NF-iSAM 已覆盖相当范围 [L23]。
- “Schur/主角度/广义特征值推出概率置信度”：主角度是响应几何；先验、绝对幅值、任务、非线性、模式质量均未由角度决定。[L18–L20]
- “只要更宽带固定 witness 距离变大就证明全局唯一”：必须允许其他材料/目标/仪器分支重新优化；这不是已有统计定理能补上的缺口。

## 4. 值得保留为新命题的具体对象

### 4.1 跨观测补偿一致性，而不是简单堆叠

多 stand-off 破解未知墙歧义已有物理先例 [L09]。新的命题必须明确：共享静止目标、动态状态限制、共享 pose gauge、跨时仪器误差分别约束哪些 compensation coefficients；新增数据能够让局部共享 null directions 消失，但 finite Maxwell alternative 可能仍存在。证明应附可计算的 lower singular bound、非线性余项和合法完整 Maxwell 实例，而不是把一条投影恒等式包装成定理创新。

### 4.2 任务误差 + 计算误差的可撤销更新证书

iSAM2 的 β/α 阈值是状态变化和局部解传播阈值，不是指定任务 Bayes risk 证书。[L01] IGG/SPO 的信息门控已更接近，但仅看 Hessian/logdet 仍看不到 innovation。因此差异应落在**真实 residual/dual residual 参与、跨相关块任务误差被界住、有限 Maxwell nonlinear correction 被验证、证书本身便宜**。若证书比直接更新昂贵，则只能成为理论结果，不能声称计算优势。

### 4.3 多模态传播等价的合法压缩

Gaussian goal-oriented low-rank approximation 已很成熟 [L19]。新命题可以研究多模态 joint map–target belief 的任务传播等价类：给定未来采集集合与损失函数，在预测分布距离可控时压缩分支，并控制后续决策风险。必须保留模式质量，合并不能只是 residual softmax 或材料图像相近。未来新观测变化时旧 equivalence 可能失效，故压缩需存可恢复 provenance。

## 5. Information-only gating 的可证明局限（本次推导）

这是对一类判据的反例，不是对 [L02] 全文或终版结论的否定。

令旧局部目标为

$$Q_{\mathrm{old}}(x)=\frac12x^{\mathsf T}H_{\mathrm{old}}x,\qquad
H_{\mathrm{old}}=\begin{bmatrix}1&a\\a&1\end{bmatrix},\qquad 0<|a|<1,$$

旧最优为 $x=0$。新增因子为 $\tfrac12(\varepsilon x_1-b)^2$。Hessian 增量是 $\varepsilon^2 e_1e_1^{\mathsf T}$，信息门控量

$$\Delta\log\det=\log\!\left(1+\frac{\varepsilon^2}{1-a^2}\right)$$

与观测创新 $b$ 无关。只更新 $x_1$、冻结 $x_2=0$ 得

$$x_1^{\mathrm{local}}=\frac{\varepsilon b}{1+\varepsilon^2},$$

而完整最优解为

$$x_1^*=\frac{\varepsilon b}{1+\varepsilon^2-a^2},\qquad
x_2^*=-\frac{a\varepsilon b}{1+\varepsilon^2-a^2}.$$

取任务 $h(x)=x_2$。在任意固定 $\varepsilon>0$、$a\ne0$ 下，局部冻结的任务误差平方随 $b^2$ 无界增长，即使 $\Delta\log\det$ 任意小。故“不大的 Fisher 增量 ⇒ 不需任务回访”在一般线性 Gaussian 模型中不成立。需要限制观测创新，或直接计算、界住残差与交叉响应。该反例可纳入 E1/E2，无需 PDE 求解。

原始 arXiv2601.08110 版本 §4 约 Eq.(25) 用信息门控量引入 inexact-GN 残差 forcing 参数；本次未找到由 logdet 单独推出该 forcing inequality 的论证。可据此指出这是需额外条件的桥梁；不能据此声称期刊终版全部收敛结论错误，终版证明尚待逐式审阅。

## 6. 各必读论文的精确范围与易错点

### L01 — iSAM2

Kaess, Johannsson, Roberts, Ila, Leonard, Dellaert. *iSAM2: Incremental Smoothing and Mapping Using the Bayes Tree*. IJRR 31(2):216–235, 2012. DOI: 10.1177/0278364911430419.

原始 PDF：https://people.csail.mit.edu/kaess/pub/Kaess12ijrr.pdf

访问 F；重点 §4.1–4.3。fluid relinearization 用增量超过 β 触发，partial state update 遇变化小于 α 停止子树传播；一般 loop closure 仍可能最坏 O(n³)，有界局部结构才支持常数级叙述。

### L02 — 2026 information gating / selective optimization

Reza Arablouei. *Efficient Incremental SLAM via Information-Guided Gating and Selective Partial Optimization*. Robotics 15(5):87, 2026. DOI: 10.3390/robotics15050087.

期刊：https://www.mdpi.com/2218-6581/15/5/87

预印本（标题较短）：https://arxiv.org/html/2601.08110

访问 F（arXiv），期刊元数据/摘要已核。Jan 2026 arXiv；publisher版本2026。用 logdet gate 选择 global/local GN，动态图受影响集合，保留全部 measurements。收敛声称需要独立逐式核实，见上述反例。

### L03 — Factor graphs

Frank Dellaert, Michael Kaess. *Factor Graphs for Robot Perception*. Foundations and Trends in Robotics 6(1–2):1–139, 2017.

作者原始页面：https://www.cs.cmu.edu/~kaess/pub/Dellaert17fnt.html

原始书稿：https://dellaert.github.io/21S-8803MM/Readings/Dellaert17fnt.pdf

访问 M/原始 PDF 已定位，未逐章核查。用于标准框架定位，技术模型与稀疏性应由本报告自己明确。

### L04 — Active SLAM

Thomas Kollar, Nicholas Roy. *Efficient Optimization of Information-Theoretic Exploration in SLAM*. AAAI, 2008.

原始出版页：https://www.aiinternational.org/Library/AAAI/2008/aaai08-217.php

作者出版目录：https://robustrobotics.github.io/publications.html

访问 A/M（原始出版页访问不稳定）。已把信息目标、地图协方差与重访联系起来，不能将 active revisit 单独当原创。

### L05 — Correlated occupancy

Sebastian Thrun. *Learning Occupancy Grid Maps With Forward Sensor Models*. Autonomous Robots 15:111–127, 2003. DOI: 10.1023/A:1025584807625.

原始作者稿：https://robots.stanford.edu/papers/thrun.occ-journal.pdf

访问 F，检查独立性讨论。已指出给定完整 map 的观测条件独立，不能替换为给定单 cell 的独立；EM forward model 保留格点依赖。

### L06 — Li 2024 multi-feature MP-SLAM

Xuhong Li, Xuesong Cai, Erik Leitinger, Fredrik Tufvesson. *A Belief Propagation Algorithm for Multipath-based SLAM with Multiple Map Features: A mmWave MIMO Application*. arXiv:2403.10095, 2024.

原始全文：https://arxiv.org/html/2403.10095v1

机构记录：https://tugraz.elsevierpure.com/en/publications/a-belief-propagation-algorithm-for-multipath-based-slam-with-mult/

访问 F。特征可在虚拟锚/点散射模型之间作 Markov 切换；联合位置姿态，采用 SAGE 估计的路径距离/AOA/AOD/幅度，并有室外 mmWave MIMO 测量。二维同步模型，不是频域复电场的体材料逆散射。

### L07 — Li 2026 radio occupancy

Xuhong Li, Erik Leitinger, Fredrik Tufvesson, Florian Meyer. *Probabilistic Occupancy Grid for Radio-Based SLAM*. arXiv:2603.03559v1, 3 March 2026.

原始全文：https://arxiv.org/html/2603.03559v1

访问 F，重点 §II/III。层次表面特征与 occupancy grid，forward sensor model 显式考虑 cell dependencies；二维水平镜面路径最高双反射，材料用 reflection factors，已允许 prior map extension/refinement。不能误写成独立 occupancy 或只点状地图。

### L08 — DNN TWR joint estimation

Fardin Ghorbani, Hossein Soleimani. *Simultaneous Estimation of Wall and Object Parameters in TWR Using Deep Neural Network*. arXiv:2111.04568 (2021); International Journal of Antennas and Propagation, article 7810213 (2022). DOI: 10.1155/2022/7810213.

原始预印本：https://arxiv.org/pdf/2111.04568

正式出版页：https://onlinelibrary.wiley.com/doi/10.1155/2022/7810213

访问 F。2D TMz FDTD，均匀非磁墙厚/介电率；已知目标数一或二，输出中心/介电率等5或8参数。文中精度百分数依赖其有限测试合同，不等于 calibrated posterior 或未知拓扑普适成功。

### L09 — Unknown walls, multiple standoff

Genyuan Wang, Moeness G. Amin. *Imaging Through Unknown Walls Using Different Standoff Distances*. IEEE TSP 54(10):4015–4025, 2006. DOI: 10.1109/TSP.2006.879325.

原始作者 PDF：https://www1.villanova.edu/content/dam/villanova/engineering/research/BestPaper2010IEEETransactions.pdf

访问 F。不同距离造成的目标位移曲线交点/共同墙参数聚焦已有，属于多观测 wall-target ambiguity 被打破的直接前例。

### L10 — IEEE TAP wall estimation

Pavel Protiva, Jan Mrkvica, Jan Macháč. *Estimation of Wall Parameters From Time-Delay-Only Through-Wall Radar Measurements*. IEEE TAP 59(11):4268–4278, 2011. DOI: 10.1109/TAP.2011.2164206.

原始作者 PDF（全文获取未成功）：https://elmag.fel.cvut.cz/sites/default/files/users/machac/files/ap_2011-2.pdf

作者出版目录：https://fel.cvut.cz/en/faculty/people/874-jan-machac/publications

访问 A/M：已核对原始 PDF 的摘要索引，但原始 PDF 的全文获取未成功。modified common midpoint / 多间距 delay-only 墙厚介电率方法；不应从当前访问程度推断更细的未解决项。

### L11 — IEEE TIM distributed joint estimation

Yu Yao, Jiahui Chen, Shisheng Guo, Xiaojian Hao, Nian Li, Guolong Cui, Lingjiang Kong. *Target Localization and Wall Parameters Estimation via Distributed Through-Wall Imaging Radar*. IEEE TIM 74:1–11, article 6500611, 2025; online DOI year2024. DOI: 10.1109/TIM.2024.3502875.

IEEE：https://ieeexplore.ieee.org/document/10758832

访问 A：IEEE 原始摘要可检索，全文获取受限。明确 simultaneous wall/target estimation 和 target position consistency，用分布式 TWI 融合；不能把未读全文的概率/迭代细节说死。

### L12 — Indoor joint target tracking / walls

Feng et al. *Multipath Ghost Recognition and Joint Target Tracking With Wall Estimation for Indoor MIMO Radar*. IEEE Transactions on Radar Systems 2:154–164, 2024. DOI: 10.1109/TRS.2024.3354509.

IEEE：https://ieeexplore.ieee.org/document/10400524

访问 A；作者全名单本轮未核实，这里不补写未经核实的姓名。原始摘要已涉及未知室内墙边界、跟踪融合和实测。

### L13 — Full-wave walls + Born targets

Francis M. Watson, Daniel Andre, William Robert B. Lionheart. *Resolving Full-Wave Through-Wall Transmission Effects in Multi-Static Synthetic Aperture Radar*. Inverse Problems, 2024. DOI: 10.1088/1361-6420/ad5b83.

原始全文：https://arxiv.org/html/2403.10354v1

机构正式版记录：https://research.manchester.ac.uk/en/publications/resolving-full-wave-through-wall-transmission-effects-in-multi-st-2/

访问 F，特别 §3.2 Eq17/18。F(V;m)=F0(m)+F1(V;m)：墙背景 full-wave Green，目标 Born；能包括入射/散射途中多次墙反射，但未包括任意目标再次交互的完整联合多次散射。POD interpolation + variable projection，重要 reduced backend 基线。

### L14 — Nonlinear TWI, known wall

Yanqing Chu, Kuiwen Xu, Yu Zhong, Xiuzhu Ye, Tianyi Zhou, Xudong Chen, Gaofeng Wang. *Fast Microwave Through Wall Imaging Method With Inhomogeneous Background Based on Levenberg–Marquardt Algorithm*. IEEE TMTT 67(3):1138–1147, 2019. DOI: 10.1109/TMTT.2018.2878703.

作者机构摘要：https://pure.bit.edu.cn/en/publications/fast-microwave-through-wall-imaging-method-with-inhomogeneous-bac

IEEE：https://ieeexplore.ieee.org/document/8537959

访问 A/M。原始机构摘要明说非均匀有限墙视为 known scatterer；不能将其描述为未知墙联合后验。

### L15 — Layered-environment inverse scattering UQ

Carolina Abugattas, Ana Carpio, Elena Cebrián, Gerardo Oleaga. *Quantifying uncertainty in inverse scattering problems set in layered environments*. Applied Mathematics and Computation 500:129453, 2025. DOI: 10.1016/j.amc.2025.129453.

原始全文：https://arxiv.org/html/2504.05776v1

出版页：https://www.sciencedirect.com/science/article/pii/S0096300325001808

访问 F。§2 是标量纵向 P-wave，背景层已知，未知 inclusion 形状和材料；自适应 FE/AD LM-Fletcher、Laplace、affine-invariant MCMC。不能被本报告称为“全 Maxwell 未知墙 Bayesian 同类成功”。模式、网格与频率讨论可作方法启发。

### L16 — HMC time-lapse FWI

P. D. S. de Lima, M. S. Ferreira, G. Corso, J. M. de Araújo. *Bayesian time-lapse full waveform inversion using Hamiltonian Monte Carlo*. Geophysical Prospecting 72(9):3381–3398, 2024. DOI: 10.1111/1365-2478.13604.

正式版：https://onlinelibrary.wiley.com/doi/full/10.1111/1365-2478.13604

原始稿：https://arxiv.org/html/2311.02999v2

访问 F；online8Sep2024。baseline posterior → monitor prior；实现中采用其 Gaussian 近似，讨论 cross-covariance 与 parallel/sequential 比较。后续物理变化不是“同一静态墙先估错再改正”的同一事件。

补充强候选但本次未做全文审计：*Bayesian variational time-lapse full waveform inversion*, GJI237(3):1624,2024，https://academic.oup.com/gji/article/237/3/1624/7641425 。可作为延伸阅读，但不能仅凭题名推导其贡献。

### L17 — Layer stripping that revisits

Nan Huai, Zhaofa Zeng, Jing Li, Yingwei Yan, Qi Lu. *Model-based layer stripping FWI with a stepped inversion sequence for GPR data*. GJI218(2):1032–1043,2019. DOI: 10.1093/gji/ggz210.

原始全文：https://academic.oup.com/gji/article/218/2/1032/5487078

访问 F，§2.5.2/2.5.3。矩形窗口版本每步只更新当前层，Hanning window 重叠平滑；stepped sequence 明确 L1(f1,f2)→L2(f1)→L1(f3)→L2(f2)→L3(f1)…，**此前层会被回访**。其最终按确定序列结束不同于永久锁死起始层。这里是2D EM FWI真实反例。

其他暂未全文审计：Bian/Yu2011，*Layer-stripping full waveform inversion with damped seismic reflection data* DOI10.1007/s12583-011-0177-6，https://journal.hep.com.cn/jes/EN/10.1007/s12583-011-0177-6；*Localized adaptive waveform inversion: theory and numerical verification* GJI2023 DOI10.1093/gji/ggac496，https://doi.org/10.1093/gji/ggac496。只能说传统方法含 data weighting/continuation，不能默认硬冻结。

### L18 — Principal-angle CRB geometry

Louis L. Scharf, L. Todd McWhorter. *Geometry of the Cramer-Rao bound*. Signal Processing31(3):301–311,1993.

原始出版页：https://www.sciencedirect.com/science/article/pii/016516849390088R

访问 A；摘要直接把敏感性子空间 principal angles 与 CRB 几何联系起来。本次不据未读全文给出其全部定理版本。

### L19 — Goal-oriented optimal posterior approximation

Alessio Spantini, Tiangang Cui, Karen Willcox, Luis Tenorio, Youssef Marzouk. *Goal-Oriented Optimal Approximations of Bayesian Linear Inverse Problems*. SIAM J Sci Comput39(5),2017. DOI: 10.1137/16M1082123.

原始全文：https://arxiv.org/pdf/1607.01881

访问 F。任务 QoI 最优低秩协方差和后验均值近似，广义 eigenpencil 平衡 prior/noise/forward/task；用于限定 task-weighted spectral compression 的已有范围。

### L20 — Certified nonlinear Bayesian dimension reduction

Olivier Zahm, Tiangang Cui, Kody Law, Alessio Spantini, Youssef Marzouk. *Certified dimension reduction in nonlinear Bayesian inverse problems*. Mathematics of Computation91:1789–1835,2022. DOI: 10.1090/mcom/3737.

原始全文：https://arxiv.org/pdf/1807.03712

访问 F。log-Sobolev 条件、likelihood-informed gradient moment 的谱尾 KL 界与采样估计，不自动覆盖任意被动 Maxwell 拓扑 prior。

### L21 — Risk-aware goal-oriented design (2026)

John D. Jakeman, Rebekah White, Bart van Bloemen Waanders, Drew P. Kouri, Alen Alexanderian. *Risk-Aware Goal-Oriented Bayesian Optimal Experimental Design*. arXiv:2609.08023v1,7Sep2026.

原始全文：https://arxiv.org/html/2609.08023v1

访问 F（重点框架）；完整作者名已从 arXiv metadata 核实。风险有后验 QoI、任务域、未来数据三层；advection-diffusion 例子，不是 Maxwell。

补充原始结果：Shijie Zhong, Wanggang Shen, Tommie Catanach, Xun Huan. *Goal-Oriented Bayesian Optimal Experimental Design for Nonlinear Models Using Markov Chain Monte Carlo*, SIAM/ASA JUQ14(1):19–47,2026 DOI10.1137/24M1649344；原始全文 https://arxiv.org/html/2403.18072v2 ，已检查框架：QoI EIG 的 nested MC + posterior MCMC + predictive-density KDE，设计空间 Bayesian optimization。属于任务 EIG 而非任意指定损失的 exact EVSI；原始预印本2024，期刊2026。

### L22 — Sequential Bayesian PDE OED (2026)

Kaichen Shen, Peng Chen. *Sequential Bayesian Optimal Experimental Design in Infinite Dimensions via Policy Gradient Reinforcement Learning*. arXiv:2601.05868,9Jan2026.

原始全文：https://arxiv.org/html/2601.05868v1

访问 A（已核原始摘要与可检索的方法段，未逐节审核全文）；finite-horizon MDP、降维、derivative-informed LANO、Laplace D-opt reward。不要将论文性能倍数当作本系统性能证据。

相关2025 JCP论文 *Sequential infinite-dimensional Bayesian optimal experimental design with derivative-informed latent attention neural operator*, JCP532:113976，DOI10.1016/j.jcp.2025.113976；https://arxiv.org/abs/2409.09141 ，原稿 https://par.nsf.gov/servlets/purl/10628357 。完整作者名单尚未核实，保留为邻近文献线索。

### L23 — NF-iSAM nonGaussian incremental inference

Qiangqiang Huang, Can Pu, Kasra Khosoussi, David M. Rosen, Dehann Fourie, Jonathan P. How, John J. Leonard. *Incremental Non-Gaussian Inference for SLAM Using Normalizing Flows*. IEEE T-RO39(2):1458–1475,2023（online2022）. DOI10.1109/TRO.2022.3216498. arXiv:2110.00876v2.

原始全文：https://arxiv.org/html/2110.00876v2

IEEE：https://ieeexplore.ieee.org/document/9946412

访问 F。已近似增量完整联合非高斯 posterior，利用 Bayes tree conditional、triangular normalizing flow，处理连续/离散不确定性。可作为“多模式不是新”证据，但不能假定 dense Maxwell 复制其复杂度。卷页/DOI 已由作者机构目录与 IEEE April2023 目录核实：https://about.uq.edu.au/experts-publication/47889/4748 ，https://www.ieee-ras.org/images/publications/t-ro/TRO.39.2.Table_of_Contents.pdf 。

补充：Olson/Agarwal *Inference on Networks of Mixtures for Robust Robot Mapping*, IJRR2013 DOI10.1177/0278364913479413，https://april.eecs.umich.edu/media/pdfs/olson2013ijrr.pdf 。Max-mixture optimization 不能等同 calibrated sum-mixture posterior。Huang/Papalia/Leonard RA-L2022 nested sampling 可作小 SLAM posterior reference，作者目录 https://qq-huang.github.io/ 。

### L24 — ScoreField inverse scattering (2026)

Wenhan Guo, Yuan Gao, Yu Sun. *ScoreField: Neural Inverse Scattering with Score-Based Generative Priors*. arXiv:2608.02937v1,3Aug2026.

原始全文：https://arxiv.org/html/2608.02937v1

访问 F，完整作者名已从原始全文核实。材料 contrast 与 induced current INR，Lippmann–Schwinger data/state consistency，加 pretrained score prior 的实例优化；不应宣称其已输出 calibrated posterior samples。

### L25 — Flow posterior sampling with MH (2026)

Jan Blechschmidt, Oliver G. Ernst, Moritz Poguntke, Björn Sprungk. *Flow Matching for Fast Posterior Sampling in Bayesian Inverse Problems*. arXiv:2610.02377v1,1Oct2026.

原始 PDF：https://arxiv.org/pdf/2610.02377

摘要：https://arxiv.org/abs/2610.02377

访问 F，重点引言/框架。conditional FM、可算密度和精度指标，independence-MH 校正；包括 EIT。训练 proposal 与正确 MH 校正本身已有，不是 W-SLAM 独立原创。完整作者名经 arXiv metadata 核实。

### L26 — Classical dual-weighted residual

Roland Becker, Rolf Rannacher. *An optimal control approach to a posteriori error estimation in finite element methods*. Acta Numerica10:1–102,2001（online2003）. DOI10.1017/S0962492901000010.

原始出版页：https://www.cambridge.org/core/journals/acta-numerica/article/abs/an-optimal-control-approach-to-a-posteriori-error-estimation-in-finite-element-methods/5C67A03F528C6FA69F37A97DF5C3BE19

访问 A：publisher 摘要明确 local residual × adjoint sensitivity 用于 task quantity error 与 adaptive mesh。故 goal-oriented dual residual 本身应标为继承标准技术，新增贡献需是 Maxwell inverse posterior/冻结分支中的特定可算证书。

### L27 — Sequential Bayesian PDE inference

Nikolas Kantas, Alexandros Beskos, Ajay Jasra. *Sequential Monte Carlo Methods for High-Dimensional Inverse Problems: A Case Study for the Navier–Stokes Equations*. SIAM/ASA JUQ2(1):464–489,2014. DOI10.1137/130930364.

原始全文：https://arxiv.org/pdf/1307.6127

作者机构出版记录：https://discovery.ucl.ac.uk/id/eprint/1425838/

访问 F（重点问题/算法框架）；PDE 初值场从序列离散观测的 Bayesian inference，SMC + adaptation/tempering/MCMC。正当 likelihood 逐次更新不是 SLAM 专属新思想，本项目必须直接比较 direct SMC 的精度与成本。

## 7. 对主论文选择的建议

**优先 Direction A，附一个可验证且小而硬的 B 结果。** A 必须加 finite-amplitude Maxwell 与共享状态耦合，不能只重述用户已有广义谱；B 必须输出实际任务误差证书，并在成本记录中胜过 ordinary incrementals/partial updates。C 作为后续系统化工作，当前 measurement-only 数据相容锚点未通过、宽材料 mode discovery 未证实，做完整 C 风险最高。

最低可接受新贡献组合：

1. “跨观测共享补偿一致性 + 非线性局部余项”的辨识命题，明确 gauge/prior/动态限制，full-wave 可实现例。
2. “innovation-aware goal residual + Schur coupling”的 freeze/reopen 任务误差界，证书不需要完整 batch 解；说明失效域。
3. 相同概率精度和相同合法材料合同下的真实成本优势；若与 direct posterior/成熟 batch 或增量法无优势，则不把 revisit policy 作为主贡献。

审计停止条件：若算法仅依赖 trace/sum sin²、gradient magnitude 或 logdet 来选区域，不能证明任务风险/近似误差桥梁；若只在预定两墙上胜过硬选单分支；若均值预测成功但 mode probability 不校准；若只比 hardfreeze；若低秩/局部证书总成本高于直接全局更新——均应停止当前主张，转向可辨识性机制论文或负结果报告。

