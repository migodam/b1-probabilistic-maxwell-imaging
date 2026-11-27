# 来源、先例与本轮新增结果

本轮先核对用户提供的 22 页报告、附包中的推导/搜索/认证代码和数值记录，并读取了《Maxwell逆散射结构改革》和上一轮《高对比度 Maxwell 逆散射研究裁决》的研究回答。三格、五格、反馈尾部的原结论均以 Pro 为提出方。用户独立 PDF 与 ZIP 中报告正文和版面相同，独立 PDF 多一项包装元数据。本文中「新增」指相对这些可查的工作；全球优先权需要进一步的专门检索。

|已有工作|已建立的内容|本轮增量|
|---|---|---|
|Z-ALG potential、R16、光学惯性|接收零空间、材料与电流不定性、共享材料读出、材料切向的电流表示上限、局部光学实性/扇区与峰值排除。R16 在 96 个旧条件上比较过最佳暗电流与普通方向。|固定外积曲率的**泛型秩证明、唯一常数左核时的完整源对谱、真实近零光学权的信号强度**。三者分别回答有多少非线性可消、余下多少约束、是否高于误差尺度。 [材料工作](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/z_alg_secant_tasks) · [R16](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/research_z_alg_r16)|
|本次 GPT Pro|固定 $\Gamma_N$、子式分解、固定 PSD 全域凸化分类、五格线性恢复、一个特殊全局可恢复条件，以及带完整反馈尾部的充分条件和反例。原报告已预警只剩互易性的危险。|在实际 24²、32² 算子上定量；五格真实材料盒内与三格材料边界各得一个**区间认证的直接材料错误极小**；实性、部分投影与随机同秩消融。|
|R7 与此前 HC-CIE|支撑受限的完整有理响应，含二格互耦端口交叉原子；公共单 Möbius 因子在内部二阶驻点上的有限下降障碍；两阶反馈存在严格有限逃逸例。|物理割线外积局域恒等式、同一内部/端口系数与可见性、对停滞材料的**全部两格必要排除**。原物理驻点障碍不自动覆盖含 prior、活动盒和有限步未收敛的历史端点。 [R7 推导](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/research_max_struct_r7/THEORY_DERIVATIONS.md)|
|MAX-STRUCT R6、R14、R15、A17、FBG、SRAP|已有 Fourier/DPCS、β 延拓、A/A+、G19/S19、受限反馈、全场 DR 和源相关阵列研究。SRAP 的半环收益、线阵失效、真实测量留出及三维面阵失败均已有公开记录。|沿同一冻结端点的精确联合响应、受限负曲率及有限合法步；新测量的是额外救援与支撑可达性。原有算法历史胜负沿用各自的曝光史、任务和预算。 [SRAP 完整正负结果](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/srap_cie_array_study_20261121/) · [A17/全基线裁决](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/top_venue_killtest_20261120)|
|SOM、NIE、CIE|直接识别强可测电流成分，随后优化剩余变量；积分方程变换、源共享与辅助优化有成熟文献。|外积消去保留的信息维数与可修复格数是另一问题。完整成像性能以真基线端点和实际新增预算为依据。 [2016 NIE](https://doi.org/10.1109/TAP.2016.2535492) · [2019 CIE](https://www.mdpi.com/2313-433X/5/2/27)|
|强散射稀疏成像及 DCTMC|Foldy–Lax 模型的有效源、多照明支撑定位及第二步材料幅度读出已有先例；T 矩阵数据相容与材料局域性也已有系统研究。|已知背景上实现这些思想是公平对照。两点交互是否能在单格端口张成空间之外量出、多少错误端点连放宽后的两点都不可达，是更窄的数学问题。 [Chai、Moscoso、Papanicolaou](https://math.stanford.edu/~papanico/pubftp/msiSubmitMiguel.pdf) · [Levinson、Markel](https://arxiv.org/abs/1401.3319)|

光滑有界集上的 Carleman 加权凸化有不同的 PDE、额外数据和变量假设。本文的障碍涉及有限矩阵上固定权重的**指定外积平方目标**；它不能否定其他凸化方法。Woodbury、正交投影、实化 Hessian、Eckart–Young 下界、两侧 MUSIC、信号子空间定位与条件最小二乘均按已有工具处理。对有限读数下的必要排除，数据需要真实支持该结论，不能用计算更快或重参数化来代替。

**独立批评。** Neon 在独立 CLI 工作区读取原 Pro 推导、算子、原始观测与既往光学结果，在 C/Accelerate 做小模型和已知开缝支撑试算。首次输出因时限停在工作进度，留下的 8 格已知支撑是 oracle 诊断。第二次直接核对本轮证明与两份端点数组，给出了完整评议。我们接受其两项限制：源依赖掩码必须使用对应的测量图；放宽到复对称矩阵只是必要条件。它提出在残差定位的边带中保持精确加权剖面。没有非 oracle 含真支撑的证明，现阶段不足以把这条建议写成恢复结果。仅在 $r=d\le s$ 情况下，常数互易约束的缺额取 $d(d+1)/2$；正文使用包含 $r,s,d$ 的完整秩公式。[理论及反例](THEORY.md)

查新覆盖了附包的相关工作、历史研究的实际推导与上述作者/出版方原文。自动学术检索遇到索引限流，文献判断采用可读取的原论文；没有声称已排尽所有等价代数结果。
