# R9先例与主张边界

当前工作支持机制辨析、可复现反例与基线修复；没有证据把它命名为普遍优于成熟CIE的新凸化算法。论文价值需要来自可推广且可反驳的认识，不能仅来自运行数量、精细一维读出或工程加速。

|主张|真正新增的证据|已有部分/不能主张的部分|当前资格|
|---|---|---|---|
|A-material改变材料路径|在保留同global current时，公开弯条有很强的材料恢复正例；多种新对象与三维也有负例|空间metric、block优化、预条件思想已有；不同current更新不能算材料metric增量|EMPIRICAL，适用区域未稳定界定|
|A+减少A材料读出成本|Möbius条件结构与接受规则的实际实现，和完整反演费用的比较|可逆坐标不消除原Φ的局部极小；MM/条件最小化本身已有|DERIVED＋EMPIRICAL，无普遍恢复保证|
|噪声适配DPCS修复|公开σ选连续奇异前缀，使新六对象所有普通CIE及A+均2/6→4/6|截断逆与discrepancy原则已有；收益必须同时给基线，不是A+独有新理论|DERIVED＋EMPIRICAL，小型确认|
|低罚权完整current初始化|公开40ac壳从错误终点转为恢复；删去后续profile阶段仍成功、费用更低|CSI/WRI/variable projection与罚权延拓已有；新确认仍2/6|EMPIRICAL，互补开发正例|
|换元/单调loss不自动凸化|用局部极小的序关系与驻点Hessian合同变换解释A+的限制|基本优化事实，不是原创定理|PROVED_UNDER_ASSUMPTIONS|
|真三维与标准离散修正|统一current后29项实材料矢量反演；Mie及粗细网格反例；标准单元积分固定诊断|不是作者同离散复现；IGT、自项积分和界面平均有成熟先例|EMPIRICAL，连续高对比精度仍有缺口|
|resident CUDA执行|完整RHS及完整更新的实际速度、残差与费用|Torch/cuBLAS/cuFFT和GMRES成熟；工程速度不增加测量信息、不证明材料更准|EMPIRICAL，限定工作量|

## 直接文献及其真正作用

- [Zhong等，NIE，2016](https://doi.org/10.1109/TAP.2016.2535492)与[Zhong–Xu，vector CIE，2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC8320915/)是实际基线来源。原论文的强散射比较必须连同采集、噪声、未知量、Fourier模态、初始化、停止及LSIE对照一起读，不能仅模仿Austria或四立方体形状。
- [van Leeuwen–Herrmann，PDE罚法](https://doi.org/10.1088/0266-5611/32/1/015007)是完整current松弛/延拓的直接近邻。R9采用它的思想作机制控制，不据一个壳正例声称发明了WRI。
- [van Leeuwen–Yang，约束松弛分析](https://doi.org/10.1088/1361-6420/adab89)的特殊有限空间二次目标，需要整个状态空间由匹配source/measurement的Riesz表示张成。读取版本为[arXiv v1](https://arxiv.org/html/2403.15292v1#S2.SS7)，不是最终期刊全文；版本间和v1公式存在的差异均在证据账本标出。当前完整内部状态与外部矩形端口不满足其特殊资格。
- “交叉熵”可能指随机群体搜索、KL数据loss、熵正则运输或相对熵先验，四者的信息与计算合同不同。R9没有用对运输计划的凸性替代对Maxwell材料的凸性，也没有实施大规模随机搜索。已核验文献、读全文/仅摘要的区别和检索429失败都在[定向证据账本](delegated/nonconvex_prior_art/PRIOR_ART_ZH.md)。
- [DDA收敛分析](https://arxiv.org/abs/0704.0033)、[矩形偶极/IGT](https://www.sciencedirect.com/science/article/abs/pii/S0022407315000382)与[高折射率filtered coupled dipoles](https://doi.org/10.1103/PhysRevE.82.036703)说明传播离散本来就需要单独分析。R9的标准立方体积分不是新Maxwell机制。
- [Farjadpour等的界面平滑](https://math.mit.edu/~stevenj/papers/FarjadpourRo06-SPIE.pdf)和[高介电CDA的weighted discretization](https://www.sciencedirect.com/science/article/abs/pii/S0030401898006452)为法向/切向有效材料提供近先例。局部静态平均适配到当前脉冲配点后，不能自动继承原FDTD收敛定理，也不能称完整WD-FCD复现。

## 目前能写什么，不能写什么

能够准确组织的稿件是：**高对比CIE中噪声硬约束、辅助材料路径与物理模型误差的机制研究**。它需要把历史弯条正例、噪声修复对所有基线的收益、壳正例的不泛化以及三维资格限制一并展示。

不能据现有证据写成“全局凸化Maxwell反演”“A+普遍击败CIE”“新增网络解决所有盆地问题”或“有限DDA低残差等于连续三维恢复”。新六例不足以给总体非劣性/成功率显著性结论。原论文对象的完整结果和同成本曲线还在独立分支更新，不能从期望结论反选测试对象。

下一项有价值的数学问题是如何在不锁死噪声的条件下保留可观测current方向，以及低罚权阶段何时能产生不同于普通CIE的有效材料方向。这需要一个可预测的条件和新的同成本确认，而不是把现有失败端点换坐标后重新命名。
