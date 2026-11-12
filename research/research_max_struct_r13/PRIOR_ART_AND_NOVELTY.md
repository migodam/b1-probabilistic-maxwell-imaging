# R13：先例、新增内容与不能声称的创新

最终新颖性判断由父线程作出。本文件不是穷尽查新证明。原ZIP的35项参考文献保留各自阅读深度，7项元数据独立核验；若干ScholarQA/Semantic Scholar请求429，原始失败回执没有改写。完整旧工作对应表在 [PRIOR_ART](delegated/evidence/PRIOR_ART.md) 与 [CLAIM_LEDGER](delegated/evidence/CLAIM_LEDGER.csv)。

## 与已有数学的区别

DR/RRR、非物理中间状态、法向漂移、反射记忆和非凸可行性搜索都是前人方法。[Elser的RRR分析](https://arxiv.org/html/1601.01003v1) 已讨论其动力学与停滞结构，不能将“离开AP停滞”本身作为新的Maxwell定理。

[Aragón Artacho–Borwein–Tam的矩阵补全工作](https://carmamaths.org/jon/DR_MatrixCompletion.pdf) §4.1 已将非凸rank约束放入DR。故每cell双场的rank-one表示和SVD投影不是本轮首创。这里额外指定实材料角度区间、跨源/三分量共享材料，并研究normal ray上的第一次材料分支竞争；其闭式阈值适用于更一般的共享实双线性约束，并非Maxwell专属。

[Hesse–Luke的局部正则性分析](https://doi.org/10.1137/120902653) 和 [Dinh–Jansen–Luke的非一致非凸松弛DR](https://arxiv.org/html/2502.12285v1) 给出各自假设下的局部定量结果。本轮没有验证足以把这些结果外推到所有材料竞争点的正则性；相等奇异值附近的投影病态尤其不能忽略。

[Collard–Lindstrom的有限外推](https://arxiv.org/html/2611.02694v1) 研究仿射子空间DR与Krylov/外推空间。本轮有界材料锥不是两个仿射子空间，因此其有限恢复结果不能直接作为事件材料恢复保证。反过来，本轮也不能把沿不变影子的整数快进宣称为全新的外推思想。

## 与Maxwell逆问题算法的区别

|先例/既有方法|已经具有的结构|R13实际比较/增量|边界|
|---|---|---|---|
|CSI、NIE/CIE、FFT-TSOM|辅助current、传播与本构分离、材料变换、Fourier/APCS、双β|保留真实成熟实现，新增同审计时钟、完整current表示消融|不将block优化或变量改名算新机制|
|NFFT-SOM|不依赖截断Fourier的完整current修正与DPCS释放|已声明Eq3–6/joint-L-BFGS适配控制|原短文未完全指定CG/停止；不是原作者成熟实现复现|
|DCTMC|数据相容T矩阵与物理材料结构之间的迭代|双场状态按真实有限照明保持共享cell材料，避免把它写成同一个T矩阵算法|可行性分裂思想已有，不能仅以更换变量声称原创|
|WRI/约束松弛、EM-ADMM、PADMM-CIE|场提升、非物理中间状态、约束残差与对偶记忆|比较具体有限材料事件、原Maxwell审计和恢复能力|没有建立DR普遍优于这些方法的理论或全规模基线证据|
|项目A/A+/G19|从CIE发展出的材料残差传播、Möbius读出及完整Green分母导数|研究表示profile、原生迭代完成度、反射与合法重启|G19是本项目算法；“full derivative”不等于“full current”|
|R8 NN初值|轻量物理反投影网络提供学习先验|冻结旧权重，公平比较NN→DR/G19/direct|无新网络创新；训练先验与训练费用不免费|
|R12共享材料DR|固定背景传播投影与共享实材料精确投影|有限事件与normal-memory控制、独立新对象、条件救援、corrector分解|R13不是首次提出共享材料DR|

NIE/CIE原文、DCTMC、WRI及ADMM对应原始来源和读取状态见上述claim ledger。CIE2019原文是 [Zhong–Xu, J. Imaging 5(2),27](https://doi.org/10.3390/jimaging5020027)。比较中 `global_dual` 指 global-PR Dual-CIE，`physical_globaldual` 才是既有总场source-Gram的物理GlobalDual（β=0）；不能混用名称填补基线。

父线程补读了 [DCTMC Part I](https://arxiv.org/pdf/1401.3319) 的迭代与计算步骤，以及 [Part II](https://arxiv.org/pdf/1505.06777) 的§VI对照。允许非局部中间相互作用、保持数据相容再逼近局部材料，明显早于本项目。其T↔V和force-diagonalization循环不能直接改名为本轮正交PC/PA；Part II的三维实验是**标量波、clean数据**，也不能作为已完成的三维矢量Maxwell比较。该文已展示物理拟合继续改善而材料误差变差，故这一现象本身不是R13原创。R13没有完成其有资格的同数据算法复现，不能宣称优于DCTMC。

另核对到2026年 [FBE-CIE＋L2/3＋ADMM](https://doi.org/10.1016/j.apm.2025.116666)。本轮取得出版方摘要，全文403；它说明CIE＋稀疏正则/ADMM已有直接先例。没有审读完整算法或运行该baseline，也不把其不同prior的文献结果与本轮表格混为同条件比较。未检索到某个词组不构成不存在前人的证明。

## 本轮可保留的贡献单元

1. 继承并独立审查有界实材料首次竞争定律，区分精确快进和近似对新动作；给出真实Maxwell中的可证伪预测，而不是重复机器精度恒等式。
2. 由Gram扰动得到给定状态球下的分支偏好界，并实际记录其不可判范围；不伪造未来轨迹误差界。
3. 推导近似AP对的接收伴随修正，明确真正控制压力解释误差的是 resolvent 作用后的pair defect。它是标准KKT/伴随结构的具体组合，贡献需由其诊断价值支撑。
4. 构造同材料影子、不同完整状态与未来物理候选的实际反例，否定“跨同一个事件就有相同恢复能力”。这不是物理多解或chaos定理。
5. 事件方向与等范数随机、简单normal倍增、A-shadow重置和物理记忆重置的严格消融；以及与受限/完整current、NN、共同corrector分离的实验。实际恢复与费用的增量只能由最终封存结果支持。

## 发表主张的取舍

不要求每个对象全胜才有算法论文，但需要一个可说明、可重复、费用合理的适用区域。单个已暴露B正例、精确事件预测或较低辅助目标都不够支撑强TAP算法主张。真正三维DDA、实测留出预测和跨网格失败必须同时放入证据边界。

若独立恢复增量不足，合理稿件应围绕“有限分支事件、扩展状态记忆与表示/选解混杂的机制研究”，不能把它包装为已经全局凸化或普遍击败CIE的新算法。若有清晰条件增量，则论文应围绕那个条件陈述，而不是用总体平均隐藏失败。连续Maxwell唯一性/稳定性与因果材料色散结论本轮没有建立，不能借有限矩阵实验替代Inverse Problems级连续理论。
