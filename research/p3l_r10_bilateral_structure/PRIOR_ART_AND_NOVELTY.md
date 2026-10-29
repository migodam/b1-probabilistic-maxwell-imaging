# R10 最接近先例与有界新颖性判断

检索日期：2026-10-29 UTC / 2026-10-30 Asia/Singapore。使用 `scholarqa-research` 的 Evidence QA 与引用核验流程；这是有明确问题的有限先例审查，不是系统综述，也不是全局新颖性保证。最终科学判断由主线程作出。本文件只保存短摘述与证据位置，不重分发第三方全文。

## 1. 直接判断

固定接收 SVD、两阶段电流子空间、T-matrix 模态分析、随机优化产生多个反演样本、prior/likelihood 降维、任务导向降维、Schur 消去和最小实现均已有明确先例。R10 不能把其中任一构件或标准恒等式作为新理论贡献。

**更直接的形状—材料先例已补核验：** Forestiere–Miano (2016) 明确给出固定几何下与介电常数无关的模式及介电常数的有理响应系数；Chen–Bergman–Sivan 的 GENOME (2019) 采用 eigenpermittivity 而非 eigenfrequency。它们比一般 T-matrix SVD 更接近 B8 的科学思想，必须作为主要冲突文献，而不能只放背景参考。当前以复 Schur 避免非正规本征向量问题、并做原噪声尺度残差/费用裁决，可形成明确实现选择，但不是首次分离几何和材料。

已完成证据给出分开的结论。固定形状、每形状 32 次幅度查询时，64 维标准 Krylov/Schur 核通过 320/320 完整残差双门，冷核比 full Schur 节省 63.76%、比直接 LU 节省 51.82%。但 Fresh10 K8 的完整 RML/shape/bilateral 主尺度覆盖分别为 3/10、3/10、2/10，双边 target 没有新增优势证据。这里可保留的是具体残差门、实现与成本结果以及明确负对照，不是新的材料谱原理或已验证的覆盖机制。整体新颖性仍为 `NOT_ESTABLISHED`；最终论文价值由主线程判断。

## 2. 检索范围、成功与缺口

研究问题分成五项：固定双边传播与已有 SOM/TSOM；形状幅度条件族和实现理论；随机目标候选与 rMAP/RTO；prior/任务相关方向；T 算子的能量与频率边界。

ScholarQA 发现查询为：`twofold subspace optimization inverse scattering`、`T matrix singular value decomposition inverse scattering geometry material`、`randomize then optimize randomized maximum a posteriori inverse problems`，每项 paper-limit=7、snippet-limit=5。六个发现端点均返回 HTTP 429，candidate IDs 为零；记录完整保留于 [scholarqa_collect.json](theory/scholarqa_collect.json)。这不是“没有相关文献”。不继续高频重试，使用已知 DOI 的 batch verification 与出版社/作者原始页面填补缺口。

指定七篇 DOI 已全部解析 7/7，见 [scholarqa_verified.json](theory/scholarqa_verified.json)。扩展检索解析见 [extended](theory/scholarqa_extended_verified.json) 与 [additional](theory/scholarqa_additional_verified.json)。一个未核实的 goal-oriented DOI 候选 `10.1137/16M1071266` 返回 unresolved，已排除；最终正确 DOI 为 `10.1137/16M1082123`，由出版社与 batch record 双重确认。

针对明确出现的 shape/amplitude 证据缺口，追加 `electromagnetic scattering eigenpermittivity modes material geometry separation`、`geometry material separation scattering eigenpermittivity modes spectral representation`、`Bergman spectral representation composite dielectric geometry material` 三项定向主来源搜索。其四篇身份全部解析，见 [eigenpermittivity verification](theory/scholarqa_eigenpermittivity_verified.json)。GENOME 正式年份为 2019，S2 的 2017 对应预印本。

Semantic Scholar 多篇 year 对应预印本而非正式出版：rMAP 是 **2018**（预印本 2016），Le Ru 等是 **2013**（预印本 2012），Spantini 等 goal-oriented 是 **2017**（预印本 2016），Zhang 等 matrix oscillators 是 **2023**（预印本 2022），Golub–Pereyra 是 **1973**（记录 1972）。下表采用正式出版社年份。

IOP/IEEE 部分原文页面不可访问，TSOM/FFT-TSOM 原篇本轮为身份 C 级；实质机制支持来自同作者 2012 PIERS 摘要以及 2010 作者综述的可检索摘要，不能冒称逐式审核原篇。未进行完整 citation-network 扩张、文献全集筛查、SPECTER 重排或完整 Scideator 40-shot 形式新颖性分类；本文件不输出“全球首次”。

## 3. 最接近工作证据账本

证据等级：A=当前获取的可识别全文支持；B=当前摘要直接支持；C=仅身份；D=题名/关系推断（不支持实质结论）。所有“R10 边界”列是本研究的解释，不是作者原话。

| ID 与正式文献 | 可核验的实质内容；证据等级 | 与 R10 的冲突/边界 |
|---|---|---|
| P01 Zhong & Chen (2009), *Twofold subspace-based optimization method for solving inverse scattering problems*. [DOI](https://doi.org/10.1088/0266-5611/25/8/085003), S2 `70fd11af9932f84523334119c3ca9a183aa3c03b` | 原篇身份 C。作者后续 [PIERS 2012 摘要](https://www.piers.org/pierspublications/PIERS2012KularLumpurAbstracts.pdf), p37，为 B：内部 current-to-field 谱进一步缩小待优化电流空间。 | 固定内部/外部传播子空间并非新概念。R10 应强调合法材料 target 和集合覆盖，不能把 TSOM 简化为两个空间的交集。 |
| P02 Zhong & Chen (2011), *An FFT Twofold Subspace-Based Optimization Method for Solving Electromagnetic Inverse Scattering Problems*. [DOI](https://doi.org/10.1109/TAP.2010.2103027), S2 `e5eff7c7a6a23de2571cded519dd3f38585cd9a7` | 原篇身份 C；同作者 PIERS 摘要 B：Fourier 近似避免部分 SVD，FFT 实现减少构造开销。作者 [发表目录](https://www.ece.nus.edu.sg/stfpage/elechenx/Publications.htm) 确认年份、页码 914–927。 | 固定 Fourier/子空间复用已有先例。R10 不能仅以“fixed basis is cheaper”宣称方法新颖。 |
| P03 Chen, Zhong & Agarwal (2010), *Subspace Methods for Solving Electromagnetic Inverse Scattering Problems*. [DOI](https://doi.org/10.4310/MAA.2010.v17.n4.a6), S2 `cf2adb06a3d58279724fc2f09b14fe68b6d6a072` | 出版社 PDF 可检索摘要 B：综述小散射体 MUSIC/强度恢复与扩展散射体 SOM。全文下载本轮受限；R9 逐式阅读只作为历史背景，不升级当前阅读等级。 | 不把位置/强度分开、固定 current 表示或多阶段逆散射当新发明。 |
| P04 Bardsley, Solonen, Haario & Laine (2014), *Randomize-Then-Optimize: A Method for Sampling from Posterior Distributions in Nonlinear Inverse Problems*. [Publisher](https://epubs.siam.org/doi/10.1137/140964023), S2 `b7868bd26be1d7cd6ca039900983bb900c569d2e` | B：随机优化在线性情形可直接给 posterior 样本；非线性需处理提案密度并可用 MH/重要性校正。 | 随机优化与多候选已有直接先例。R10 的未校正、过滤、有限预算云只能作可能性集合，不直接是准确 posterior。 |
| P05 Wang, Bui-Thanh & Ghattas (2018), *A Randomized Maximum A Posteriori Method for Posterior Sampling of High Dimensional Nonlinear Bayesian Inverse Problems*. [Publisher](https://epubs.siam.org/doi/abs/10.1137/16M1060625), S2 `fc845bcfa838e2c0c8c03ed257a9c5b3ec09ea75` | B：特定 data/prior 均值随机化对应 RML；逐样本求 PDE 约束优化，并讨论近似 Metropolization 降偏。 | R10 完整 RML 是必要的成熟强基线；shape/bilateral target 协方差是其候选变体，不能改名为新采样定律。 |
| P06 Cui, Martin, Marzouk, Solonen & Spantini (2014), *Likelihood-informed dimension reduction for nonlinear inverse problems*. [Author full text](https://arxiv.org/html/1403.4680v2), [DOI](https://doi.org/10.1088/0266-5611/30/11/114015), S2 `8b4dd26286ce50efb6a7137e895e86d8d6a1a084` | A/B：比较 likelihood 与 Gaussian prior 的影响，构造 LIS 并在其补空间保留 prior。 | “只看数据弱方向”与 prior-aware 子空间不同；R10 不能将 prior/likelihood 广义谱称为新概念，也不能把低 rank 直接等同材料覆盖。 |
| P07 Spantini, Cui, Willcox, Tenorio & Marzouk (2017), *Goal-Oriented Optimal Approximations of Bayesian Linear Inverse Problems*. [Publisher](https://epubs.siam.org/doi/10.1137/16M1082123), [Author text](https://arxiv.org/html/1607.01881v2), S2 `137f439e51982c7219aaa40d82e7e63372989c48` | A/B：在线性 Gaussian 问题中，面向 QoI 的低秩近似同时考虑 prior、forward、noise 和任务。 | “按任务分配方向”已有紧邻先例。R10 需以不同的非线性材料集合机制与实证说明差别；当前实际候选臂尚非完整任务导向算法。 |
| P08 Suryadharma, Fruhnert, Rockstuhl & Fernandez-Corbaton (2017), *Singular-value decomposition for electromagnetic-scattering analysis*. [Publisher](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.95.053834), S2 `fe3060d7b3cdbe064804024033e6ad2189e5a4c2` | B：用 T-matrix SVD 分解散射响应，区分正交 singular modes 与非正交 eigenmodes。 | T-matrix SVD 不是 R10 新构件；其模式针对给定散射体，不自动是所有材料共享的低维空间。 |
| P09 Le Ru, Somerville & Auguié (2013), *Radiative correction in approximate treatments of electromagnetic scattering by point and body scatterers*. [Publisher](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.87.012504), [Author PDF](https://www.wgtn.ac.nz/scps/research/research-groups/raman-lab/eric-le-ru/publications/2013-PRA-Radiative-correction-for-electromagnetic-scattering.pdf), S2 `23f12cf5324617ebd4c85c78fa9ea62dc909c192` | A/B：K/reactance 表示组织能量守恒及辐射校正，可推广到较高多极阶。 | 标准能量约束不是新理论；当前有限近场 OTE 需独立推导归一化，不能直接视为完整 power-normalized T 或 S。 |
| P10 Molesky, Chao, Jin & Rodriguez (2020), *Global T operator bounds on electromagnetic scattering: Upper bounds on far-field cross sections*. [Publisher full text](https://doi.org/10.1103/PhysRevResearch.2.033172), S2 `5d9d0a2427489d1cb0a51eb1e174d7802d3d00e7` | A：利用净实功率与无功功率约束进行散射设计上界；其结构放松明确区分必要约束与真实几何可实现性。 | 能量满足不等于存在某个原 decoder 材料，正如 receiver-null 电流不等于合法材料。此工作研究设计界，非本轮材料候选覆盖。 |
| P11 Zhang, Monticone & Miller (2023), *All electromagnetic scattering bodies are matrix-valued oscillators*. [DOI](https://doi.org/10.1038/s41467-023-43221-2), [Author PDF](https://millergroup.yale.edu/sites/default/files/files/2023_11-Zhang.pdf), S2 `e027cb5357a2d766cd0943648d59372ccef2f140` | A：频率响应表示内嵌 passivity/causality，全文及补充材料说明通道与因果条件。 | R10 的 \(\delta\) 是固定频率下材料幅度，不是 \(\omega\)；不能从 \(\delta\)-resolvent 导出宽频 sum rule、因果色散或频谱界。 |
| P12 Golub & Pereyra (1973), *The Differentiation of Pseudo-Inverses and Nonlinear Least Squares Problems Whose Variables Separate*. [Publisher](https://epubs.siam.org/doi/abs/10.1137/0710036), S2 `5d38e06821c7f88d27b5b28039383c10a5c00fe3` | 身份 C；作者 1976 [chapter summary](https://www.sciencedirect.com/science/article/pii/B9780125142502500123) 为 B：可分离线性/非线性参数的 variable projection。 | 完整多重散射对幅度有理非线性，不是直接套线性 VarPro。当前 block-GN Schur 只是标准二次步求解；完整非线性 profile `NOT_RUN`。 |
| P13 Wang, Sun, Dai, Zhang & Su (2024), *Cross-Correlated Subspace-Based Optimization Method for Solving Electromagnetic Inverse Scattering Problems*. [DOI](https://doi.org/10.1109/TAP.2024.3450328), S2 `c0bc431a56b37f41d6f9fc27a9a7335a4aa6fd30` | B（S2 所取原论文摘要）；IEEE [issue contents](https://engage.ieee.org/rs/756-GPH-899/images/T-AP_Nov2024-Electromagnetics.pdf?version=0) 核验身份：结合 CC-CSI 与 SOM，报告质量/复杂度取舍。未读全文。 | SOM 家族持续发展，不能只对比最早 TSOM 就宣布全面领先。该摘要尚不支持其已实现 R10 的 finite-set coverage 目标。 |
| P14 Forestiere & Miano (2016), *Material-independent modes for electromagnetic scattering*. [Publisher](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.201406), [author preprint](https://arxiv.org/abs/1611.04501), S2 `a6e63d4dfd0828a9499c60767dccc4c96c7bb8ff` | B：任意形状均匀介质物体的全波响应使用与介电常数无关的模式，系数为介电常数的有理函数。 | **直接冲突：固定形状、变化材料幅度的模式/有理分离并非 R10 首创。** 本轮 soft occupancy、宿主与端口合同不同，但不足以自动构成新物理。 |
| P15 Chen, Bergman & Sivan (2019), *Generalizing Normal Mode Expansion of Electromagnetic Green’s Tensor to Open Systems*. [Publisher](https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.11.044018), S2 `9a1546efe58f0c5953dcff167abfa77ad0b9775d` | B：以介电常数作为本征值，展开开放系统 Green 张量，并允许源位于散射体内外。 | **直接冲突：eigenpermittivity 与固定几何复用已有成熟先例。** 本轮 Schur 实现和数值 residual guard 的收益要独立测量，不宣称新模态原理。 |
| P16 Bergman, Chen & Farhi (2020), *Scattering electromagnetic eigenstates of a two-constituent composite and their exploitation for calculating a physical field*. [Publisher](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.102.063508), S2 `841fc92908b3005b45d17ce222104998d44992f6` | B：固定宿主和微结构的双正交电磁本征态，可复用于由入射场或源电流产生的场计算。 | 可复用材料/结构谱理论不限于最初的静电近似。R10 的普遍“几何固定即可复用”不能作为独立新颖性。 |

补充排除：Ye 等 (2012), [*Electromagnetic imaging of separable obstacle problem*](https://pubmed.ncbi.nlm.nih.gov/22330461/), DOI `10.1364/OE.20.002206`，B 级摘要表明“separable”指未知散射体与已知障碍物空间不重叠，并非本轮 shape/amplitude 参数分离。不能仅凭标题把它当同一机制。

更早的材料谱背景为 Bergman (1979), [*Dielectric constant of a two-component granular composite: A practical scheme for calculating the pole spectrum*](https://doi.org/10.1103/PhysRevB.19.2359)，S2 `f81fd8233c742844c674b36e78e6762360b7a9fb`。B 级摘要已说明对给定微结构，把有效介电常数视为组成材料介电参数的极点/留数问题；其自伴算子背景不可直接套为当前开放、非正规 B 的 Hermitian 性。

最小实现史料：[Ho & Kalman (1966), NASA record](https://ntrs.nasa.gov/citations/19670049337)，C 级身份；原始全文未获取。R10 的 reachable/observable quotient、Hankel rank 与相似实现均给出标准独立推导，不将文献身份误当物理实验验证。

## 4. 主张—冲突表

| 拟议主张 | 最近冲突 | 本轮可保留的更窄表述 | 尚需证据 |
|---|---|---|---|
| 首次固定收发双边结构 | Born/VIE、SOM/TSOM、T-matrix 均有成熟结构 | 给出该离散合同下可复用的明确算子与依赖关系 | 不作为新颖性主张 |
| 首次利用谱分开可见/不可见 | P01–P03、P08、经典 realization | 区分 current、材料 tangent、传递族最小实现 | 固定 all-D 或 decoder 内压缩需另证 |
| 首次随机优化形成概率体积 | P04–P05 | 在原物理噪声门下比较有限材料集合及完整费用 | 独立测试覆盖、范围/重复性、强 RML 对照 |
| 首次 prior 或任务导向方向 | P06–P07 | 固定宿主双边 prior pullback target covariance 的具体可复现变体 | 方向 ablation、不同任务与真实覆盖 |
| 形状幅度解耦使非线性已消失 | P12、P14–P16；标准 resolvent、Schur | 固定形状多幅度计算可复用；联合 GN 内幅度块可标准消去 | per-shape setup、摊销、内层分支、原支持 |
| transfer 多样性证明材料不同模式 | P08；经典 moment/realization | 有限 counterfactual transfer sketch 与材料距离联合分配 | 同池对照及对真实材料覆盖的增益 |
| 能量约束给物理合法材料证书 | P09–P11 | 必要能量条件可筛错；共同构成/decoder 仍需满足 | 正确通道/归一化与可实现性 |
| \(\delta\) 的极点是宽频共振定律 | P11 研究真实频率的因果被动性 | 这里只是固定频率上的材料幅度有理结构 | 多频物理、色散模型与新推导均 `NOT_RUN` |

## 5. 三机制的相对新颖性与可证伪点

**M1 双边 prior/task 形状提案。** 最接近的机制组合是 TSOM 的固定传播组织与 LIS/goal-oriented prior–data–task 比较。单纯把两者拼接不够。实际 R10 的 covariance eigenclip/trace 规则是启发式；三种新臂用相同 solver/caps 比较，但实付 RHS 略不同。Fresh10 双边对完整 RML 的覆盖差 −0.1，98.75% 配对区间 [−0.4,0]；均匀动作 MSE 差 +0.0000368256，区间 [−0.000234365,+0.000361607]。当前裁决为 `NO_ADDED_BENEFIT_SUPPORTED`，非全部双边机制的普遍不可能性。总体新颖性：`NOT_ESTABLISHED`。

**M2 分离随机反演。** 最近先例是 rMAP/RML、variable projection 与 P14–P16 的固定几何材料谱。实际 block-GN 对同一二次模型代数等价于 joint GN，所以不能把其写成新求解原理。完整非线性 amplitude profile 仍为 `NOT_RUN`，仅形状 target 在 Fresh10 与完整 RML 同为 3/10 覆盖。独立幅度计算核取得成本收益，但未集成进候选反演，不能代称“分离反演已经更快”。当前为 `STANDARD_COMPONENTS / SCOPED_KERNEL_BENEFIT / INVERSION_BENEFIT_UNPROVEN`。

**M3 transfer-aware 多样性。** 标准 moment/realization 只提供表示。对同样物理合格的 finite cloud，基于 transfer 与材料联合距离分配预算具有可检验性，但可能删去重要响应等价材料。须比较同池随机/材料最远点策略。当前 `PROPOSED / NOT_RUN`；尚不足以宣称唯一或首次。

## 6. 论文叙事应如何受到证据约束

最终叙事采用三个分别裁决的结果：A，固定形状重复幅度计算有实测费用收益；B，新 target 未显示材料覆盖优势；C，三个条件数值参考通过任务门，但 Full66 posterior 为 `UNQUALIFIED_NOT_ESTIMATED`。R9 的合法宽云覆盖不足仍是动机，不能在 R10 中用计算核速度取代这一未解决的科学目标。

64 维由八个 DEV 形状的 32/64/96/128 阶扫描选定，再在十个另组锚点形状上验证；这建立经验有效阶，不建立严格最小实现阶。320 节点核的材料真值在其他支线已曝光，当前仅相对 DEV 压缩 heldout，不能写新的盲法材料测试。只有 M32 是该核实际完整查询工作量，其余 M 是实测均值外推。标准 Krylov/Galerkin、Schur 与 eigenpermittivity 思想均不重新命名为新物理。

Fresh10 的 N10 配对统计不足以声称高概率覆盖；三种新臂 matched solver/caps，历史三臂有更小 cap 和不同约束，不能将六臂称为费用匹配。已有证据不支持原创覆盖机制成功，也没有否定所有可能的结构提案。全文证据链接为 [候选评分](results/candidates_fresh10/STATISTICS.json)、[独立统计与审计](INDEPENDENT_QA.md)、[另组形状核](REDUCED_RESPONSE_FRESH10.md)。有限检索不认证全球新颖性。
