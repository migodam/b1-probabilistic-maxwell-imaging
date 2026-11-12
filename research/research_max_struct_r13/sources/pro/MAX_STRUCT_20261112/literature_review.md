# MAX-STRUCT 文献与原创性审查

审查日期：2026-11-11（UTC）。本文件是并行文献研究的独立记录；不代替对仓库实验的审计。它区分读过全文、读过作者或出版社摘要、仅追溯到书目信息。检索覆盖了 CIE 的近期后续工作，不能等同于穷尽全部先行文献。

## 1. 对研究方向影响最大的结论

**最值得警惕的重合，不是“别人也做过 DR”，而是已有研究分别完成了你可能打算作为核心创新的三个步骤：扩展场空间缓解非凸性、从数据估计内部场后线性反演材料、以及在可观测空间外加入先验。** 强论文必须说明新结构把这些步骤连接后增加了什么可计算、可检验的能力。

1. CIE 2019 已有两块更新不同 β、逐轮减小 β 并增加 Fourier 子空间、全三维矢量实验。[R05] 因此普通双 β、continuation、扩大 Fourier rank 和“3D 可用”不能单独构成新意。
2. WRI/IR-WRI 已经将物理约束暂时放松，允许非物理场，利用双线性交替更新和乘子记忆改善错误 basin。[R08–R11] “放开内部场+共享材料+记忆”本身是成熟工具组合；DR 与 ADMM 的关系必须讲清楚。
3. CC-CSI 早已把状态残差和观测残差联系到接收域。[R12] 仅把本构残差经传播算子送到接收器再惩罚，必须与它比较。
4. 数据驱动 ROM、Lippmann–Schwinger–Lanczos（LSL）直接从边界数据构造内部场，再使用双场乘积形式的 Lippmann–Schwinger 恒等式反演材料。[R21–R23] 因此双侧传播、有限差分恒等式、内部场先于材料，都不是未经探索的新方向。
5. **2026 年 Bevacqua–Isernia–Crocco 的 CIE 两阶段方法尤其接近**：大 β 恢复电流与 support，再重建 total field，再解正则化线性材料问题。[R07] 把 G19/DR 输出电流拿来做材料线性读出，不能声称首次提出。
6. 对固定线性观测保持数据一致并学习其零空间的框架已有严格理论。[R20] Maxwell 的真正难处是材料到电流的非线性兼容性，以及不同状态上的观测纤维如何弯曲；把 Jacobian 小奇异方向称为“可辨识商空间”尚不够。

**可保留的潜在原创性空间**：对受限辅助表示导致的材料搜索障碍，给出与完整 resolvent、当前照明和接收器共同有关的定量证书；用该证书指导低成本表示补全/反射搜索；对噪声球内分支另外提供可辨识性诊断和材料选择规则。这仍只是值得验证的定位，不是已证明的新定理或已有算法胜利。

## 2. 25 篇核心文献：贡献、重合与边界

每一项只摘取本任务所需内容；“本研究含义”是本次分析，不是原作者的结论。

### A. Maxwell / CIE / CSI

**R01 — van den Berg & Kleinman (1997), “A contrast source inversion method”, Inverse Problems 13, 1607–1620.**

- DOI：[10.1088/0266-5611/13/6/013](https://doi.org/10.1088/0266-5611/13/6/013)。访问：出版社摘要；CSI 实现机制同时由 R02 全文支持。
- 贡献：把材料与各照明下 contrast source 联合作为未知量，允许简单加入先验；基础对照。
- 本研究含义：共享材料 (J_ell=mQ_ell)、交替 field/current/material 更新、避免每次精确正演，均有明确前史。创新必须在约束几何、受限表示、数据与材料指标的区别上，而不是变量命名。

**R02 — van den Berg & Abubakar (2001), “Contrast Source Inversion Method: State of Art”, PIER 34, 189–218.**

- DOI：[10.2528/PIER01061103](https://doi.org/10.2528/PIER01061103)。[作者论文全文](https://www.jpier.org/ac_api/download.php?id=0106113)。访问：出版社页面及全文抽取。
- 贡献：讨论 CSI 的局部极小值，采用乘法正则化和边缘保留的加权 TV；目标包含数据与对象方程残差。
- 本研究含义：不要把“低 residual 并不等于正确材料”作为发现本身；应提供以前没有的可计算分类与预测。MR-CSI 是比裸 CSI 更公平的传统对照。

**R03 — Zhong & Chen (2011), “An FFT twofold subspace-based optimization method for solving electromagnetic inverse scattering problems”, IEEE TAP 59, 914–927.**

- DOI：[10.1109/TAP.2010.2103027](https://doi.org/10.1109/TAP.2010.2103027)。[作者出版目录](https://www.ece.nus.edu.sg/stfpage/elechenx/Publications.htm)。访问：作者目录和后续 CIE 全文对其分解的明确复述，未获得此论文全文。
- 贡献边界：观测算子的确定电流与模糊电流分解，以及第二折 Fourier 子空间正则化。
- 本研究含义：DPCS/APCS/Fourier 子空间不是实现细节，而是方法所容许的状态集合。固定受限空间可降低方差，也可排除恢复真实材料所需要的内部电流方向；需以原算法相同空间与完整空间作因子实验。

**R04 — Zhong, Lambert, Lesselier & Chen (2016), “A New Integral Equation Method to Solve Highly Nonlinear Inverse Scattering Problems”, IEEE TAP 64, 1788–1799.**

- DOI：[10.1109/TAP.2016.2535492](https://doi.org/10.1109/TAP.2016.2535492)。访问更新：本轮后续成功读取[共同作者 Marc Lambert 上传的原始全文](https://www.researchgate.net/publication/296475773_A_New_Integral_Equation_Method_to_Solve_Highly_Nonlinear_Inverse_Scattering_Problems)，实际审读§II的LS→NIE推导、§III相关解释与§IV合成/实测实验部分；不声称逐字审过全部论文。页面标明作者上传日期2016-05-15。本条取代较早的原文访问缺口。
- 原式(3)：\(\beta I=R(\beta I+E^{inc}+GI)\)，\(R=\beta\chi/(1+\beta\chi)\)，是等价变换。原文已讨论Austria环内假目标、Fourier限制与噪声的冲突，以及观测不能充分区分的材料“plateau”。
- 本研究含义：辅助表示造成的恢复限制及低残差错误材料并非本轮首次发现。新贡献须是精确定量、可预测性和受控算法干预；原始二维实验可按其设置复现，不能再以原文不可读为由替代成小网格演示。本轮通过上述作者上传全文完成访问。

**R05 — Zhong & Xu (2019), “Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems”, Journal of Imaging 5, 27.**

- DOI：[10.3390/jimaging5020027](https://doi.org/10.3390/jimaging5020027)。[开放全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8320915/)。访问：完整理论、算法步骤与实验设置。
- 核心式：式(3)为 \(\beta J=R(\beta J+E^{inc}+G_DJ)\)，\(R=\beta\chi/(1+\beta\chi)\)；式(15)是 field/current 两项目标。第一折确定电流来自接收算子 SVD，第二折为 Fourier APCS。第2.4节明确采用 current/material 不同 β。
- 实验：60 个电偶极子、三圈布置、接收三分量；60³ 正演、30³ 反演，10% 白噪声。
- 本研究含义：CIE 的变量变换不改变物理前向映射的可辨识性。应分析有限残差时的度量、坐标和子空间如何改变搜索。

**R06 — Bevacqua & Isernia (2024), “Support Reconstruction of Dielectric and Metallic Targets via the Contraction Integral Equation”, IEEE TAP 72(3), 2643–2653.**

- DOI：[10.1109/TAP.2024.3353314](https://doi.org/10.1109/TAP.2024.3353314)。[作者机构记录](https://iris.unirc.it/handle/20.500.12318/144708)。访问：机构摘要与元数据，未读全文。
- 贡献：强目标或合适超参数使 CIE auxiliary contrast 接近目标内1、目标外0，利用 CSI 恢复 support；包含非凸介电/金属目标、实测、多频与色散。
- 本研究含义：大 β 的材料压缩既可帮助 support，又会损害幅值分辨；若新方法主张 support/contrast 解耦，应直接比较。

**R07 — Bevacqua, Isernia & Crocco (2026), “A Low-Nonlinearity Approach to Inverse Scattering via an Unconventional Exploitation of the Contraction Integral Equation”, IEEE AWPL.**

- DOI：[10.1109/LAWP.2026.3660123](https://doi.org/10.1109/LAWP.2026.3660123)。[机构摘要](https://iris.unirc.it/handle/20.500.12318/167028)。访问：检索返回完整机构摘要；直接打开被403阻挡。
- 贡献：大 β CIE 先估 contrast sources/support；由电流得到 total field；然后把数据与未知 contrast 建成线性关系，解带 support 的正则化逆问题；合成和实测验证。
- 本研究含义：这是对“稳定电流先恢复、材料后读出”最直接的近期威胁。新的场读出必须提供额外证书、误差界或跨材料分支访问能力。

### B. 扩展场、变量消元与残差耦合

**R08 — van Leeuwen & Herrmann (2013), “Mitigating local minima in full-waveform inversion by expanding the search space”, GJI 195, 661–667.**

- DOI：[10.1093/gji/ggt258](https://doi.org/10.1093/gji/ggt258)。[作者实验室的 WRI 说明](https://slim.gatech.edu/research/full-waveform-inversion)。访问：出版社摘要、作者页面；出版社全文重定向无法读取。
- 贡献：放松波动方程，先估满足数据与物理折中的全场，再在固定场下更新材料。
- 本研究含义：非物理内部状态能够扩大搜索空间是重要基线，不能视为 DR 独有机制。必须比较完整场 AP/WRI/ADMM 与 DR，才能归因反射记忆。

**R09 — van Leeuwen & Herrmann (2016), “A penalty method for PDE-constrained optimization in inverse problems”, Inverse Problems 32, 015007.**

- DOI：[10.1088/0266-5611/32/1/015007](https://doi.org/10.1088/0266-5611/32/1/015007)。访问：出版社摘要、作者机构记录、R10–R11 对公式的完整复述；原PDF下载403。
- 贡献：二次罚函数与变量投影消去扩展场，得到 reduced penalty objective。
- 本研究含义：\(\min_u\|Pu-d\|^2+\lambda\|A(m)u-b\|^2\) 的 Schur complement、场消元和诱导数据度量已成熟。重新写成 VIE 不足以形成基本理论创新。

**R10 — Aghamiry, Gholami & Operto (2019), “Improving full-waveform inversion by wavefield reconstruction with the alternating direction method of multipliers”, Geophysics.**

- DOI：[10.1190/geo2018-0093.1](https://doi.org/10.1190/geo2018-0093.1)。[作者全文](https://arxiv.org/pdf/1809.00891)。访问：全文。
- 贡献：IR-WRI 用增广拉格朗日/ADMM 累积数据和源残差，避免无限增大 penalty 导致的病态；波动方程双线性使每块是线性问题。
- 本研究含义：“反射记忆”与“残差累积乘子记忆”不是同一更新，但应在相同完整表示下比较。若只胜过无记忆 AP，难以证明优于成熟扩展场方法。

**R11 — Rizzuti, Louboutin, Wang & Herrmann (2021), “A dual formulation of wavefield reconstruction inversion for large-scale seismic inversion”, Geophysics 86.**

- DOI：[10.1190/geo2020-0743.1](https://doi.org/10.1190/geo2020-0743.1)。[作者全文](https://slim.gatech.edu/Publications/Public/Journals/Geophysics/2021/rizzuti2020dfw/rizzuti2020dfw.html)。访问：全文。
- 贡献：数据噪声球作为约束、物理残差作为目标的 dual WRI；使用 source weighting；用常规时间域传播降低扩展场的存储/求解成本。
- 边界：论文讨论三维可扩展性，但主文明确未做三维实验。
- 本研究含义：噪声球约束不自动产生新的 inverse selection 原则。需说明你的材料先验与 source residual 选择不同在哪里。

**R12 — Sun, Kooij, Jin & Yarovoy (2017), “Cross-Correlated Contrast Source Inversion”, IEEE TAP 65(5), 2592–2603.**

- DOI：[10.1109/TAP.2017.2673758](https://doi.org/10.1109/TAP.2017.2673758)。[作者机构摘要](https://research.tudelft.nl/en/publications/cross-correlated-contrast-source-inversion/)，[代码](https://github.com/TUDsun/CC-CSI)，[预印本](https://arxiv.org/abs/1906.10864)。访问：机构摘要，预印本全文读取被禁用。
- 贡献：在测量域耦合 state error 和 data error 的 cross-correlated cost。
- 本研究含义：观测加权状态误差本身已有先例；若构造 \(S R_m e\) 或其近似，必须明确与 CC-CSI 的算子、导数与优化块有什么实质差别。

### C. 投影、反射与局部几何

**R13 — Hesse & Luke (2013), “Nonconvex Notions of Regularity and Convergence of Fundamental Algorithms for Feasibility Problems”, SIAM J. Optimization.**

- DOI：[10.1137/120902653](https://doi.org/10.1137/120902653)。[作者全文](https://arxiv.org/pdf/1212.3349)。访问：全文。
- 贡献：非凸集合的局部 regularity 与 coercivity 给出 AP 和部分 DR 情形的局部线性收敛。
- 本研究含义：在 Maxwell 上验证/使用这些条件有价值，但把一般定理改写成 \(\mathcal A,\mathcal C\) 不构成新的全局 basin 理论；可行错误材料同样可能是稳定交点。

**R14 — Bauschke, Bello Cruz, Nghia, Phan & Wang (2014), “The rate of linear convergence of the Douglas–Rachford algorithm for subspaces is the cosine of the Friedrichs angle”, J. Approximation Theory 185, 63–79.**

- DOI：[10.1016/j.jat.2014.06.002](https://doi.org/10.1016/j.jat.2014.06.002)。[全文](https://arxiv.org/pdf/1309.4709)。访问：全文。
- 贡献：对子空间 DR 得到锐的 \(\cos\theta_F\) 因子，并与 AP 比较；迭代本身和投影 shadow 的极限要区分。
- 本研究含义：同样局部切向下 AP 的每周期渐近因子可为 \(\cos^2\theta\)，所以 DR 的经验逃逸能力不能解释成“局部收敛总更快”。应记录有限距离状态离开哪些稳定/活动区域。

**R15 — Liang, Fadili & Peyré (2017), “Local Convergence Properties of Douglas–Rachford and ADMM”, JOTA.**

- DOI：[10.1007/s10957-017-1061-z](https://doi.org/10.1007/s10957-017-1061-z)。[作者全文](https://arxiv.org/abs/1606.02116)。访问：摘要与作者公开 PDF 片段。
- 贡献：在凸、部分光滑及相应非退化条件下有限时间识别流形并进入局部线性收敛；局部 polyhedral 情况速率与切向 Friedrichs angle 有关。
- 本研究含义：不可直接把此凸函数定理用于非凸 Maxwell 本构集合；活动边界、材料盒约束和不唯一投影必须单独处理。

**R16 — Elser (2003), “Phase retrieval by iterated projections”, JOSA A 20, 40–55.**

- DOI：[10.1364/JOSAA.20.000040](https://doi.org/10.1364/JOSAA.20.000040)。[作者预印本](https://arxiv.org/abs/math/0111080)。访问：作者摘要、期刊摘要。
- 贡献：difference map 把两种投影组成带参数的搜索映射，应用于非凸相位恢复。
- 本研究含义：非物理 iterates、投影 shadow 与反射式绕开停滞都是已知思想。Maxwell 新内容应该是约束构造、材料共享导致的几何和对搜索成败的可证预测。

### D. 可观测性、数据一致与选择

**R17 — Darbas, Heleine & Lohrengel (2019), “Sensitivity analysis for 3D Maxwell’s equations and its use in the resolution of an inverse medium problem at fixed frequency”.**

- DOI：[10.1080/17415977.2019.1588896](https://doi.org/10.1080/17415977.2019.1588896)。[作者预印本](https://arxiv.org/abs/1911.02409)。访问：出版社理论介绍、作者摘要。
- 贡献：对介电常数与电导率进行场灵敏度分析，并用边界变化定位/估计小幅扰动，包含三维数值例。
- 本研究含义：\(DF(m)h\) 的精确 Maxwell 推导属于基础。创新需要超出小扰动敏感度：有限距离分支、受限表示误差、对算法的可计算作用。

**R18 — Bucci & Isernia (1997), “Electromagnetic inverse scattering: Retrievable information and measurement strategies”, Radio Science 32, 2123–2137.**

- DOI：[10.1029/97RS01826](https://doi.org/10.1029/97RS01826)。访问：出版社摘要。
- 贡献：由散射场近似带限结构分析固定频率有限尺寸目标可恢复信息维度，并研究单/多照明测量策略。
- 本研究含义：离散模型 full rank 不等于噪声下所有像素可辨识；“自由度”要注明分辨阈值与 admissible model class，不能从一个 Jacobian 数值秩推出全局唯一性。

**R19 — Bouchet, Rotter & Mosk (2021), “Maximum information states for coherent scattering measurements”, Nature Physics 17, 564–568.**

- DOI：[10.1038/s41567-020-01137-4](https://doi.org/10.1038/s41567-020-01137-4)。[作者预印本](https://arxiv.org/abs/2002.10388)。访问：作者摘要；PDF超出网页读取大小限制。
- 贡献：最优入射态来自散射矩阵参数导数构成的 Hermitian Fisher information 算子，实验验证精度改善。
- 本研究含义：强传播/透射通道不必是强参数敏感通道。新工作可把材料向量、接收 aperture、noise whitening 与分支辨别联合起来；若只是优化 \(DF^*DF\)，属于信息设计的标准化应用。

**R20 — Schwab, Antholzer & Haltmeier (2019), “Deep Null Space Learning for Inverse Problems: Convergence Analysis and Rates”, Inverse Problems 35, 025008.**

- DOI：[10.1088/1361-6420/aaf14a](https://doi.org/10.1088/1361-6420/aaf14a)。[全文](https://arxiv.org/pdf/1806.06137)。访问：全文。
- 贡献：在线性 \(A\) 下 \(x\mapsto x+(I-A^+A)N(x)\) 严格保持数据一致；与传统正则化结合得到收敛与速率。
- 本研究含义：材料非线性 \(F\) 的 Jacobian 核只提供一阶一致，有限步必须 retraction 或重新检查正演。电流零空间 \(\ker S\) 严格保持电流数据，却未保证共享材料本构可实现。

### E. 数据驱动内部场、有理结构及实测基准

**R21 — Borcea, Druskin, Mamonov, Zaslavsky & Zimmerling (2020), “Reduced Order Model Approach to Inverse Scattering”, SIAM J. Imaging Sciences 13(2), 685–723.**

- DOI：[10.1137/19M1296355](https://doi.org/10.1137/19M1296355)。[作者稿](https://arxiv.org/abs/1910.13014)，[全文](https://websites.umich.edu/~borcea/Publications/ROM_INVERSION.pdf)。访问：作者摘要和全文片段。
- 贡献：从时间分辨的阵列 scattering matrix 构造波传播 ROM，利用 Galerkin 结构定量反演；与 Born 数据变换和仅 support 成像区分。
- 本研究含义：先提取稳定的可观测算子再恢复材料有强前史。可考虑作独立备选，但单频有限孔径数据通常不足以直接构造相同 ROM。

**R22 — Druskin, Moskow & Zaslavsky (2021), “Lippmann–Schwinger–Lanczos algorithm for inverse scattering problems”, Inverse Problems 37, 075003.**

- DOI：[10.1088/1361-6420/abfca4](https://doi.org/10.1088/1361-6420/abfca4)。[全文](https://arxiv.org/pdf/2101.12317)。访问：全文。
- 核心式：式(2.16–2.19)用 resolvent 恒等式将两个材料的数据差写成材料差与两侧内部场的积分；Lanczos ROM 给出未知内部场近似，使材料求解近似线性。
- 边界：背景/真值的正交化内部基近似一致是文中 Conjecture 2.2 的关键，不应误称全面证明；最初模型是 Schrödinger 型。
- 本研究含义：有限 secant 恒等式不是新理论；新意只能来自 Maxwell 约束下的证书/选择/算法和可验证的误差控制。

**R23 — Borcea, Liu & Zimmerling (2024), “Electromagnetic inverse wave scattering in anisotropic media via reduced order modeling”, JCP 515, 113272.**

- DOI：[10.1016/j.jcp.2024.113272](https://doi.org/10.1016/j.jcp.2024.113272)。[全文](https://arxiv.org/html/2403.03844v1)。访问：全文。
- 贡献：从瞬态 Maxwell 数据构造保留因果传播结构的 ROM，估计内部波并用于成像和介电张量反演。
- 边界：无损非磁各向异性模型；理论针对3D，但数值为 cylindrical/orthotropic 约化二维。各向异性边界数据存在 boundary-fixing diffeomorphism 非唯一性。
- 本研究含义：适合多频/瞬态分支。不能把文中较好目标形状当作所有高对比、有损、稀疏孔径问题的保证。

**R24 — Chen, Bergman & Sivan (2019), “Generalizing Normal Mode Expansion of Electromagnetic Green’s Tensor to Open Systems”, Physical Review Applied 11, 044018.**

- DOI：[10.1103/PhysRevApplied.11.044018](https://doi.org/10.1103/PhysRevApplied.11.044018)。[作者全文早期版](https://arxiv.org/pdf/1711.00335)。访问：全文早期版。
- 贡献：以介电常数为本征值，固定频率、几何和背景的 eigenpermittivity mode 展开开放系统 Green tensor；通过 LS 方程处理内外源。
- 本研究含义：极点—留数有价值但结构依赖参数化。固定 support、单材料标量变化可以复用模态；任意逐像素材料变化一般不共享固定极点，更不必是单一 Möbius 圆。

**R25 — Geffrin & Sabouroux (2009), “Continuing with the Fresnel database: experimental setup and improvements in 3D scattering measurements”, Inverse Problems 25, 024001.**

- DOI：[10.1088/0266-5611/25/2/024001](https://doi.org/10.1088/0266-5611/25/2/024001)。[官方数据库](https://www.fresnel.fr/3Ddatabase/)，[参数](https://www.fresnel.fr/3Ddatabase/parameters.php)。访问：论文摘要与官方实验设置。
- 贡献：3D目标、弱散射测量、漂移补偿、五个公开目标包括当时的 mystery target。
- 本研究含义：公共数据已经暴露，不能称 blind target；可做真正独立的标定/测量检验。官方说明入射为远端抛物面天线，可近似 plane wave；频率3–8 GHz，低频需注意信噪比。不要把所有 Fresnel 套件都按相同入射建模。

## 3. 补充近期风险与可用资源

- Du et al., *Fast Physics-Driven Untrained Network for Highly Nonlinear Inverse Scattering Problems*, [arXiv:2602.13805](https://arxiv.org/html/2602.13805v1)，完整 HTML 已查。其组合为 Fourier 电流参数化+CIE+无训练网络+材料补偿+bridge penalty。它可作速度/表征对照，但作者称 CIE 物理电流更低频的表述需要独立验证；等价状态方程变换不会改变同一物理材料的精确电流。用定量原始电流误差而非复述该主张。
- Zhong et al., *A Multiresolution Contraction Integral Equation Method…*, [10.1109/TMTT.2019.2956939](https://doi.org/10.1109/TMTT.2019.2956939)，机构摘要已查。多分辨率 CIE 不是空白。
- *Wavelet-Based Subspace Regularization…*, [10.3390/electronics9111760](https://doi.org/10.3390/electronics9111760)，出版社摘要已查。换 Fourier 为 wavelet 有前例。
- Borcea et al., *Waveform inversion via reduced order modeling*, [arXiv:2202.01824](https://arxiv.org/abs/2202.01824)，Geophysics DOI [10.1190/geo2022-0070.1](https://doi.org/10.1190/geo2022-0070.1)。完整阵列响应与 collocated source/receivers 是关键；可由插值/reciprocity 近似，但不能默认为无代价。
- *Regularized Reduced Order Lippman–Schwinger–Lanczos Method…*, [arXiv:2311.16367](https://arxiv.org/abs/2311.16367)，机构摘要已查，明确研究 resolvent ROM 的正则化。研究极点学习应先比较该路线。
- CC-CSI 开源：[TUDsun/CC-CSI](https://github.com/TUDsun/CC-CSI)。WRI dual 频域演示：[slimgroup/Software.rizzuti2019SEGadf](https://github.com/slimgroup/Software.rizzuti2019SEGadf)。需要许可证、单位、噪声模型与停止准则审计之后再纳入统一对照。

## 4. 对新数学主张的逐项审查

### 4.1 哪些推导应作为引理，而不是论文主定理

令 \(R_m=(I-D_mG)^{-1}\)，\(Q_m=E^{inc}+GJ_m\)。有限维、逆存在时：

\[
DF(m)h=S R_m D_h Q_m,
\qquad F(m_2)-F(m_1)=S R_{m_2}D_{m_2-m_1}Q_{m_1}.
\]

这些是标准 resolvent/对比介质恒等式，第二式与 R22 的双场乘积结构同源。它们非常适合统一本任务，但不能把简短代数证明等同于发现新 Maxwell 结构。

相同地，\(\ker DF(m)\) 是一阶不可观测材料方向；\(DF^*\Sigma^{-1}DF\) 是 Gauss–Newton/Fisher 信息；Schur complement 消去场得到有效材料度量；线性投影反射速率由 principal angle 控制。这些都是研究工具。

### 4.2 什么可能形成真正有内容的 Maxwell 命题

以下是本次审查提出的研究目标，不是文献空白认证：

1. **受限表示的障碍界。** 固定 DPCS/Fourier 电流集合 \(\mathcal V\) 下，真实材料或某个正确材料方向至少需要多大的本构/物理残差；界能否用 \((I-P_{\mathcal V})R_mD_hQ_m\)、接收 whitened sensitivity 与 active material tangent 计算，而不需要GT。
2. **量纲一致的几何关系。** 把场—电流约束切向的广义角度与可辨识材料信息联系起来，显式消除纯场自由度、材料尺度和与数据无关的公共切向；得到可算 generalized eigenproblem，而不是任意欧氏角。
3. **有限步证书。** 将线性化不可观测方向与有限 secant 数据差分开，导出材料沿某方向可走多远且不离开噪声球的界，显示在何种照明/共振条件下 Jacobian 结论失效。
4. **可验证选择。** 在多个噪声球内候选间，明确每条先验或新测量选择的材料信息，给出不可由单组数据消除的 ambiguity，而非把corrector下降当作GT证书。

这些目标必须对既有 ring/shell/bridge 正反例给出彼此不同的预测，并在新对象上成立；若所有指标仅与 residual 单调等价，贡献仍不足。

### 4.3 五个最危险的审稿反例

1. 两个材料都落在同一噪声球。任何固定数据的 deterministic solver 都无法保证选择GT；几何再漂亮也不能制造测量信息。
2. 改变 \(Q\) 与 \(J\) 的单位/权重，principal angle 变化而物理信息不变。需要物理/噪声归一化或广义度量。
3. 将 G19 电流空间补全后收益完全追上 DR。则先前“reflection优势”主要来自表示；应如实接受并把论文转成表示障碍与自适应补全。
4. 完整场 ADMM/IR-WRI 同成本复制 DR。则贡献可能是 Maxwell 实现与实测，而不是新的搜索动力学理论。
5. 新 prior 使开发形状更好，却把真实 bridge 或连续材料抹掉。选择原则必须报告错误先验反例、候选多样性与 holdout measurement，而不只NRMSE。

## 5. 期刊定位：事实与研究判断分开

官方范围：TAP 发表电磁波与材料相互作用、天线传播/感知的原创贡献；SIIMS 强调成像的数学和计算方法并要求正式分析与有竞争力的实用结果；Inverse Problems 强调对广泛逆问题读者有意义的显著原创新科学。来源分别为：[TAP scope](https://ieeeaps.org/ieee-tap/about-the-transactions/scope-of-articles-in-the-transactions)、[SIIMS editorial policy](https://epubs.siam.org/siims/editorial-policy)、[Inverse Problems scope](https://publishingsupport.iopscience.iop.org/journals/inverse-problems/about-inverse-problems/)。以下是判断，不是期刊承诺。

| 目标 | 最有说服力的主贡献 | 决定性证据 | 主要风险 |
|---|---|---|---|
| 强 TAP | 与材料恢复相关的传播/接收反馈、明确打破辅助表示障碍且计算合理 | 原 CIE 2D/3D重现，真实矢量3D与Fresnel标定，损耗/孔径变化，成本公平 | 仅为CIE/CSI的一种新权重；实测增益来自额外先验或标定 |
| 强 SIIMS | 统一且可计算的材料可观测几何，导出表示自适应/选择算法 | 新命题+反例，表示×动力学因子实验，独立预测验证，充分对照 | 推导等价于GN/变量投影；理论只局部而结论宣称全局 |
| 强 Inverse Problems | 从 finite secant、quotient 与 admissible class 得到新的稳定性/非唯一性或信息界 | 严格条件、可达到/不可达到界、噪声/失配理论与少量有判别力数值 | 只是离散矩阵恒等式；没有真正超出经典resolvent/正则化理论 |

建议先组织为一篇聚焦论文：**“表示障碍与材料可观测几何如何控制非线性 Maxwell 可恢复路径”**。若最后只有工程改进，则走 TAP 叙事并诚实弱化新理论；若确实获得独立的稳定性/非唯一性定理，再拆理论文。不要先把一组实验切成两篇高度重复论文。

建议章节：问题与已有正反例 → Maxwell精确增量/受限表示 → 定量障碍与几何证书 → 由证书导出的最简算法 → 因子机制实验 → 独立对象/3D/实测 → 非辨识与失配边界。

## 6. 给母任务的可直接引用来源

核心全文 URL：

- CIE2019：https://pmc.ncbi.nlm.nih.gov/articles/PMC8320915/
- IR-WRI：https://arxiv.org/pdf/1809.00891
- WRI dual：https://slim.gatech.edu/Publications/Public/Journals/Geophysics/2021/rizzuti2020dfw/rizzuti2020dfw.html
- DR angle：https://arxiv.org/pdf/1309.4709
- 非凸DR：https://arxiv.org/pdf/1212.3349
- Null-space learning：https://arxiv.org/pdf/1806.06137
- LSL：https://arxiv.org/pdf/2101.12317
- EM ROM：https://arxiv.org/html/2403.03844v1
- GENOME：https://arxiv.org/pdf/1711.00335
- CIE2026直接竞争：https://iris.unirc.it/handle/20.500.12318/167028


## 7. 定向追加：DPCS 解锁、精确 profile 与有限 DR 分支事件

本节是在阅读仓库既有 Pro 报告第3节和参考文献后，针对本轮新候选做的额外审查；不是把旧报告书目再包装成新发现。

### 7.1 最直接先行工作：NFFT-SOM 已经专门纠正固定 DPCS 的噪声误差

**R26 — Wei, Chen, Zhao & Chen (2016), “Two FFT Subspace-Based Optimization Methods for Electrical Impedance Tomography”, PIER 157, 111–120.** DOI [10.2528/PIER16082302](https://doi.org/10.2528/PIER16082302)，[出版社全文](https://www.jpier.org/ac_api/download.php?id=16082302)。已读全文主要推导和实验讨论。NFFT-SOM 将传统 SOM 的小奇异值子空间换成完整 Fourier 空间；LF-SOM 又改变确定电流的构造。它分析了高噪声和截断维数敏感性。因此“DPCS 不是不可改动的真值”在 SOM 社群内已有清楚前史。

**R27 — Wei & Chen (2018), “New FFT Subspace-based Optimization Methods for Inverse-Scattering Problems”, IEEE APS/URSI conference contribution.** [原始两页论文](https://usncursi.org/archive/aps-ursi/2018/papers/1927.pdf)。已读全文；DOI 尚未核实，不能虚构。其式(3)明确为

\[
J_p=J^s_p+F\alpha_p,\qquad \alpha_p\in\mathbb C^M,
\]

这里 F 是**完整** Fourier 基，而不是先投影到 dominant singular vectors 的正交补。作者把这点直接解释成可纠正 deterministic current 中的噪声误差；使用 Austria 和 30% 加性白噪声进行比较。它仍最小化加权数据/状态残差，不是显式噪声球约束。该论文是新方法必须正面对比的最近邻，而非可略去的两页旁支。

**R28 — Zhang, Ma, Xu & Zhong (2020), “Wavelet-Based Subspace Regularization for Solving Highly Nonlinear Inverse Scattering Problems with Contraction Integral Equation”, Electronics 9,1760.** DOI [10.3390/electronics9111760](https://doi.org/10.3390/electronics9111760)，[全文](https://www.mdpi.com/2079-9292/9/11/1760)。已读§3相关公式：式(14)直接写 APCS 为 \(W^{-1}\gamma\)，没有显式 \((I-V_DV_D^*)\) 投影；逐轮加入 wavelet 系数。不能假定所有 TSOM 后续都硬锁 DPCS，复现实验应核对具体实现。

对照之下，**R05 CIE2019 式(10)** 明确为

\[
J^a=(I-V_DV_D^*)\operatorname{IDFT}\gamma.
\]

所以仓库 G19 的固定 DPCS 与 Fourier 双重限制有原论文依据，不是仓库无意的编码错误；但“去掉它”也不是首次提出的算法思想。

由此建议把新贡献限定为：识别固定 DPCS 和 Fourier 截断分别在**材料 profile**上产生何种可计算偏置，区分合理降噪正则化与伤害真实材料分支的表示障碍，并用受噪声预算约束的低维修正实现可预测的收益。必须包含四个同动力学对照：原空间；完整模糊空间但保持 DPCS；允许 DPCS 校正但保持低 Fourier；NFFT/全 current 空间。若解锁收益由旧 NFFT 完全复制，应把它归为机制解释与再设计，不能宣称新原理。

### 7.2 精确 profile 恒等式的数学地位

**R29 — van Leeuwen & Yang (2025), “An Analysis of Constraint-Relaxation in PDE-Based Inverse Problems”, Inverse Problems 41,025009.** DOI [10.1088/1361-6420/adab89](https://doi.org/10.1088/1361-6420/adab89)，[全文](https://arxiv.org/html/2403.15292v2)。已读 Theorem1、Corollary1、极限解释及讨论。其 representer theorem 将松弛源写成接收伴随场的有限组合；消元后给出 \(\tfrac12 e^*(I+\rho^{-1}G(m))^{-1}e\)，G 是材料相关 Gram 矩阵。作者明确指出实际大问题的度量近似和噪声影响仍待研究。有限非共址阵列也在 Remark1–2 中处理，只是每个照明可能使用不同 Gram。

**R30 — Golub & Pereyra (1973), “The Differentiation of Pseudo-Inverses and Nonlinear Least Squares Problems Whose Variables Separate”, SIAM J. Numer. Anal.10,413–432.** DOI [10.1137/0710036](https://doi.org/10.1137/0710036)。已查出版社摘要。其内容是消去线性变量、对常秩伪逆/投影求导，不能把新的 current 消元作为首次变量投影。

**R31 — Golub & Pereyra (2003), “Separable nonlinear least squares: the variable projection method and its applications”, Inverse Problems19,R1–R26.** DOI [10.1088/0266-5611/19/2/201](https://doi.org/10.1088/0266-5611/19/2/201)。已查出版社摘要及元数据，未逐页读全文。

设固定 m 后的 current 目标是严格凸二次型，H(m) 正定，完整最小解为 j★，受限空间是 \(j_D+\operatorname{range}(V)\)。则

\[
f_V(m)=f_{\mathrm{full}}(m)+\frac12\operatorname{dist}_{H(m)}^2
\big(j^\star(m),j_D+\operatorname{range}(V)\big)
\]

直接来自配方与 H-正交投影。嵌套空间的差又是加权 Pythagoras。**它可以成为关键解释性引理，不能单独充当高影响力新定理。** 若 H 仅半正定，需另给可达性、零空间/商空间定义和最小解选择，不能直接引用正定结论。

更有内容的目标是证明/验证材料排序翻转：物理较好的 m₁ 可能因表示惩罚较大而有 \(f_V(m_1)>f_V(m_2)\)，从而将“辅助投影损失”提升为实际算法偏好错误材料的因果机制。还应区分计算 exact profile 的昂贵离线诊断与实际算法只使用的廉价误差估计，不能把完整解的代价隐藏在证书里。

对于

\[
J=J_D+V_Db+B_Fa,
\]

把 data residual 平方分成 \(\|\Sigma_Db\|^2+\|SB_Fa-Y_\perp\|^2\) 需要 \(U_D^*SB_F=0\)，例如 \(B_F=(I-V_DV_D^*)F_{\rm trunc}\)。若直接用未正交化 Fourier 基，就必须保留交叉项。显式噪声预算和低秩校正相对 NFFT2018 有设计差别，但噪声球/约束正则化本身也不是新数学。

### 7.3 DR 的法向漂移已有前史；材料锥的精确首事件可以更具体

**R32 — Elser, Rankenburg & Thibault (2007), “Searching with iterated maps”, PNAS104,418–423.** DOI [10.1073/pnas.0606359104](https://doi.org/10.1073/pnas.0606359104)。检索获得正文“The Difference Map”完整段落；PMC直接访问此次被验证码阻挡。该段已清楚解释 AP 在非零最近点对停滞，而 difference map 在两条错开的正交线模型上执行 \(z\mapsto z+\beta c\mapsto z+2\beta c\)，沿近解间的位移离开。故“反射累计法向记忆所以能够逃离 AP 陷阱”不是新定性机制。

**R33 — Elser (2017), “Matrix product constraints by projection methods”, J. Global Optim.68,329–355.** DOI [10.1007/s10898-016-0466-9](https://doi.org/10.1007/s10898-016-0466-9)，[全文预印本](https://arxiv.org/pdf/1601.01003)。已读矩阵积投影及§4。它在 proximal-pair trap 中明确将 ADMM 变量解释为 discrepancy accumulator，逐次重新居中投影并逃离陷阱，再给出 RRR 与 ADMM 的关系。这比只引用一般DR收敛定理更贴近本任务的动力学叙事。

**R34 — Dizon, Hogan & Lindstrom (2022), “Circumcentering Reflection Methods for Nonconvex Feasibility Problems”, Set-Valued Var. Anal.30,943–973.** DOI [10.1007/s11228-021-00626-9](https://doi.org/10.1007/s11228-021-00626-9)，[全文](https://arxiv.org/pdf/1910.04384)。已读问题设定、主要保证及误差讨论。通过反射点外心加速，在若干非凸原型下有局部超线性/二次性质，也指出数值误差能破坏理论速度。它不是沿材料锥法向到首分支事件，但“几何外推加速DR”不能泛称新思想。

**R35 — Themelis, Stella & Patrinos (2022), “Douglas–Rachford splitting and ADMM for nonconvex optimization: accelerated and Newton-type linesearch algorithms”, Comput. Optim. Appl.82,395–440.** DOI [10.1007/s10589-022-00366-y](https://doi.org/10.1007/s10589-022-00366-y)，[全文](https://arxiv.org/pdf/2005.10230)。已读 AssumptionI 和算法框架。它使用 Douglas–Rachford envelope 与 quasi-Newton/Anderson方向加速；主要非凸DR假设包括一项具有 Lipschitz 梯度，不能直接当作两个非凸indicator投影问题的全局保证。

本轮候选的新意应精确写成：在 \(c=P_Ca,\ a=P_Ac,\ d=a-c\) 的非零 AP 最近点对处，利用 A 凸与材料共享锥 C 的特殊结构，给出保持 \(P_C(c+td)=c\) 的**有限法向区间**，计算第一个投影分支改变事件；它预测何时 shadow 必须改变，而不预测哪一改变更接近 GT。

对无边界实 rank-one 锥，阈值 \(\|v\|/\|w\|\) 本质是正交 rank-one 分量的奇异值交叉，属于 Eckart–Young 低秩投影的直接推论。真正可能提升内容的是共享复照明实化后的具体几何、材料盒约束的端点事件、按体素取首事件、近似 AP 对的扰动界、以及事件外推与若干普通 DR 步严格等价的条件。

若想形成论文主结果，需要三个额外层次：

1. **稳健性：**实验状态不是精确最近点对；给出 projector 唯一性裕度与 AP误差如何控制预测事件区间。接近奇异值相等点，投影不连续/多值必须显式处理。
2. **可预测性：**在独立对象上事先预测 DR 停台长度/首材料改变位置，与实测轨迹比较；不能在看到轨迹之后调阈值。
3. **算法价值：**事件跳跃必须减少真实CPU并保持成功率，比较同成本普通DR、惯性/Anderson与简单倍增外推；不能只演示公式在解析玩具模型上正确。

事件改变既可能进入好分支，也可能进入错误分支；若多个材料本来都数据可行，该理论不提供 GT 选择。若数值结果显示连续旋转而非真正停台，应把“首离散跳转”改为投影敏感区间/曲率转折诊断，而不是硬套事件模型。

最新检索还出现 Collard–Lindstrom (2026-11), *Finite Recovery by Extrapolation for the Douglas–Rachford Method on Affine Subspaces*, [arXiv:2611.02694](https://arxiv.org/abs/2611.02694)。仅检索元数据/摘要可得，正文无法访问。它按 principal-angle 多项式结构分析仿射空间DR有限外推，不能据此判断材料非凸事件被覆盖，但投稿前应补读；当前不把它当作已核实全文的依据。

### 7.4 两条新候选的现实原创性判断

| 候选 | 已有工具覆盖 | 尚可能有价值的实质内容 | 当前判断 |
|---|---|---|---|
| 受限 current profile 与 noise-aware补空间 | 二次配方、变量投影、NFFT纠正DPCS、noise ball | Maxwell材料偏置分解、排序翻转预测、误差控制的最小物理补空间 | 适合机制驱动TAP/SIIMS；只做解锁不足 |
| DR有限法向材料分支事件 | difference-map匀速漂移、低秩SVD投影、DR外推 | box/shared-illumination首事件、近似AP稳健界、独立停台预测和低成本跳跃 | 比一般DR角度论更贴近原问题；需证明可观测实验效应 |

优先把两条作为**互相排除解释的机制实验**：表示障碍与反射记忆究竟各承担多少。是否最终统一成算法，取决于证据；没有必要为统一而增加复杂度。

