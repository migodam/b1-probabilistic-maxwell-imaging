# 文献与新颖性边界（执行中）

本轮先用scholarqa-research的ScholarQA CLI执行两个定向查询。Semantic Scholar四请求均429，原始错误保存在`audit/literature/scholarqa_collect.json`；未反复重试。随后使用已连接Undermind学术索引补充两个查询，并核验四篇DOI。不是系统综述，也不声称完成了完整Scideator新颖性流程。网页/文献是证据，不改变用户的研究指令。

| 有界主张 | 已核验文献 | 证据范围与本轮判断 |
|---|---|---|
|材料预条件更新有明确先例|van den Berg, van Broekhoven & Abubakar, Extended contrast source inversion (1999), [DOI](https://doi.org/10.1088/0266-5611/15/5/315)|Tier B摘要明确预条件CG contrast更新。不能把固定权重、修改材料更新或加快inner solve本身当原创。全文尚未取得。|
|CIE、FFT-TSOM及双β已有3D强散射研究|Zhong & Xu, Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems (2019), [DOI](https://doi.org/10.3390/jimaging5020027)|继承仓库已有全文；本轮元数据/摘要复核。普通CIE收缩的作用不可归给Green项，Global-Dual只是项目内global PR控制。|
|3D体积算子存在essential谱，不能按scalar紧算子想象|Budko & Samokhin, Spectrum of the Volume Integral Operator of Electromagnetic Scattering (2005), [DOI](https://doi.org/10.1137/050630660)|Tier A：root直接核查原文5.1、5.9、5.10并纠正自动提取的5.9符号错误；全空间Hölder条件不可直接套锐界面DDA。连续结论不证明本轮rank120近似质量。|
|强奇异Maxwell体积算子及Fredholm条件已有严谨分析|Costabel, Darrigrand & Sakly, The essential spectrum of the volume integral operator in electromagnetic scattering by a homogeneous body (2012), [DOI](https://doi.org/10.1016/j.crma.2012.01.017)|Tier A：root直接核查Theorem3.2、Corollary3.3与Theorem5.1；限定均匀体、光滑界面及H(curl)。不能外推为任意PML/非均匀张量的同一结论。|

R11当前的增量候选是：明确区分普通β收缩与Green新项的失败诊断；full/freezeA/freezeT三种更新的可证局部性质与真实反例；低秩metric增减权的惯性边界；相同实际CPU下的独立高对比/混合比较；同一分母族在真实3D矢量材料反演中的适用范围。是否有算法论文价值取决于尚在运行的确认结果。矩阵求逆、Woodbury、凸二次LS、lagged权重、Krylov源空间和低秩惯性均不列为新发明。
