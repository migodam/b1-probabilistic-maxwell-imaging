# R7 聚焦先例证据核验

2026-10-28。结论限于本次取得的证据：任意阶散射导数、current 子空间分解、非均匀背景更新、prior-relative 广义谱以及多照明共享材料均有明确先例。R7 的区别必须落在具体补偿后响应、有限尺度预测与用途上；本审计没有决定其最终原创性，也没有证明检索之外不存在更接近工作。

## 范围、失败与复核

按 `scholarqa-research/SKILL.md`、retrieval/reporting 完整流程执行一次focused collect：三条query，每条paper/snippet各limit5。六个操作全部HTTP429，返回候选0；原结果保存在 `SCHOLARQA_COLLECT.json`，没有换检索词无限重试。随后按已获授权，利用R7T已有DOI/arXiv身份检查primary文本。批量verify成功5/5，`unresolved_ids=[]`；2015补充primary成功1/1。SOM/TSOM的定向abstract读取仍为null，publisher-elided，OA状态CLOSED，原文缺口保留。

完整逐claim证据等级、canonical S2 ID、短原文、章节/公式定位、URL及限定见 [CLAIM_LEDGER.tsv](CLAIM_LEDGER.tsv) / [JSON](CLAIM_LEDGER.json)。[检索合同](SEARCH_PLAN.md)、[费用与检索收据](RETRIEVAL_RECEIPT.json)、[primary访问记录](WEB_PRIMARY_READBACK.json)均在本目录。没有新增Maxwell、训练或R7 GT读取；没有新worker。R7T的26条文献仅作为种子目录，不能写成26篇均获本轮原文复核。

## 五篇指定文献与直接证据

|文献身份|取得证据|能支持的先例边界|仍缺什么|
|---|---|---|---|
|Arhab, Anagnostou, Joelson (2018), Wave Motion83:67–79, DOI10.1016/j.wavemoti.2018.07.008|**A**：作者preprint§§2.3–2.4、eqs14–27、递推图2及有限收敛数值例|将扰动材料×场作为新的等效source，逐次得到任意阶field导数并展开数据；raw current/field jets不是新概念|这不等于R7补偿后首个有效阶、共同U冲突预测或同成本优势。[原稿](https://arxiv.org/html/1803.00002v1)|
|Xudong Chen (2010), IEEE TGRS48(1):42–49, DOI10.1109/TGRS.2009.2025122|**C原篇**：S2身份+作者目录；**A后续primary归述**：Xu2015 eqs3–5，ref26|SOM把current-to-data GS的可辨/歧义部分分开，不需要把整个current都当未知优化|原2010全文和abstract未取得，不能写成原篇全文核验。[作者目录](https://www.ece.nus.edu.sg/stfpage/elechenx/Publications.htm)|
|Yu Zhong & Xudong Chen (2009), Inverse Problems25:085003, DOI10.1088/0266-5611/25/8/085003|**C原篇**；**A后续primary归述**：Xu2015 eqs6–7，ref29|TSOM再用内部GD的谱限制歧义current搜索空间；不只是对外部GS做一次截断|原IOP全文restricted，abstract被elide；不得伪称原文直接证明材料fiber。[DOI](https://doi.org/10.1088/0266-5611/25/8/085003)|
|Xiuzhu Ye & Xudong Chen (2017), IEEE TAP65(12):7224–7232, DOI10.1109/TAP.2017.2766658|**B**：机构存档原abstract|S-DBIM每步更新非均匀背景Green函数，并线性恢复deterministic current以改善field估计；这是必须保留的非线性强对照|只核了abstract，未复制其数值优势，也未核完整优化约束。[机构记录](https://pure.bit.edu.cn/en/publications/subspace-based-distorted-born-iterative-method-for-solving-invers/)|
|Tiangang Cui, James Martin, Youssef Marzouk, Antti Solonen, Alessio Spantini (2014), Inverse Problems30:114015, DOI10.1088/0266-5611/30/11/114015|**A**：作者稿§3.2 eqs7–10|local LIS比较GN数据敏感性与prior precision；global LIS跨posterior平均prior-preconditioned GN谱|R7的材料RMS pullback与prior metric不同；这种差别本身不自动构成新理论或算法收益。[原稿](https://arxiv.org/pdf/1403.4680)|

**书目更正：** R7T `source/09_sources.tex:55–57` 用了arXiv题名“diffracted field”搭配正式journal DOI。S2 DOI规范记录和publisher检索结果的正式题名为“High-order functional derivatives of the **scattered** field according to the permittivity-contrast function”；arXiv1803.00002仍为“**diffracted** field”。属于同篇版本题名差异，应同时标明而不是默默当两个不同论文。历史输入未改写。[正式publisher身份](https://www.sciencedirect.com/science/article/abs/pii/S0165212518301161)

## primary补充：current子空间与共同材料

为弥补SOM/TSOM原文未开放，读了作者公开的 Xu, Zhong, Song, Chen & Ran (2015), *Multiplicative-Regularized FFT Twofold Subspace-Based Optimization Method for Inverse Scattering Problems*, IEEE TGRS53(2):841–850, DOI10.1109/TGRS.2014.2329032。该作者proof在eqs1–4为每个incidence使用不同current、同一chi，eqs5–7明列GS/GD两类子空间，并把其来源分别归于Chen2010和Zhong2009。原PDF有proof/终稿重复块，定位用印刷页码与公式号。[作者公开PDF](https://renchengsong.github.io/papers/TGRS2329032.pdf)

这提供**同材料多照明是标准setup**的直接证据，但未在本次有限阅读中取得与R7“共享补偿坐标的加权冲突量、最早阶、有限工作区间预测”完全相同的具体结果。这里的“未取得”不是“不存在”。没有把一般joint inversion、两类current谱或高阶导数单独改名为R7新理论。

## 对R7材料广义可观测性的有限比较

以下是分析者对已取得公式的结构比较，不是文献已经证明R7命题：

- SOM/TSOM的原生未知分量是source-specific current；GS描述向receiver的映射，GD描述内部传播。current-space不可见不自动给所有源共享的合法material变化。
- LIS的谱以prior为基准衡量参数对数据的局部影响，并通过posterior支持的平均生成全局空间。R7若以材料RMS去除decoder坐标条件数，须明确实材料切空间、rank变化和gauge；不能仅因换metric就作原创性结论。
- Arhab提供high-order field导数及source递推的直接竞争；Ye提供动态背景更新的传统算法竞争。R7新增主张需要相对这些强普通方法核验实际有限响应或用途，不能只展示Taylor小步匹配。

本次429使邻域发现受限；未追引用图、未核完所有26文献、未作2026近期全面更新，也未检得针对共同补偿冲突的穷尽先例。建议顶层将这些结果作为claim边界与强比较依据；最终新颖性和科学裁决仍留在顶层。
