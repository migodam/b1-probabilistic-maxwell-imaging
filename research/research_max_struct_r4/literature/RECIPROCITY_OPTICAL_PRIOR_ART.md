# R4：互易性与光学定理约束的实现先例

2026-11-03。只核对最近实现先例；保留原有 `PRIOR_ART.md`。没有运行测试、读取 GT、提交或发布。父线程提供的 `EᵀJ` 对称性与矩阵光学定理仍由父线程审查适用性、符号、离散资格及算法。

## 直接结论

**已经找到“互易性/守恒派生约束改善 lifted nonlinear inversion”的明确实现先例：DCTMC Part II。** 这比泛泛的散射定理引用接近得多。它每轮投影单体素 polarizability 以满足被动/无吸收条件，并利用互易数据副本约束实验 T-matrix；算法在迭代中仍允许非局部 interaction matrix。不能把“在物理方程完全满足以前加入这些物理约束”整体说成首次提出。([Levinson–Markel, 2016](https://doi.org/10.1103/PhysRevE.94.043318))

本次限定搜索**没有核实到**直接给 CIE/CSI/SOM 的辅助电流矩阵施加 `EᵀJ=(EᵀJ)ᵀ` 与 `Im_H(E*J)=J*ΓJ` 的逐式同一实现。这个空缺为 **UNRESOLVED**，不是“没有先例”，更不是原创结论。

## 最近证据与差异

|来源|实读证据与具体操作|与父线程候选的区别|
|---|---|---|
|**R1** Howard W. Levinson & Vadim A. Markel (2016), *Solution of the nonlinear inverse scattering problem by T-matrix completion. II. Simulations*. [正式DOI](https://doi.org/10.1103/PhysRevE.94.043318), [作者预印本v3](https://arxiv.org/abs/1505.06777v3)|**A，完整作者稿。** §III.B，PDF p4，Eq26–28，`local_only/dctmcII2016.txt:216–264`：被动性给逆 polarizability 虚部不等式；无吸收时变成等式；将相应投影用于每次 DCTMC 迭代。§V.B，PDF pp12–13 / Fig8，`txt:737–813`：加入 source/receiver 互换所得 reciprocal 数据，不进行新物理测量，使实验 T 更接近对称；直接比较迭代行为。|能直接排除宽泛“互易性或能量约束从未用于 lifted 逆散射”表述。它的守恒约束是**单体素材料/polarizability**投影，互易性通过**外部数据增广与 T-matrix completion**实现；所读代码描述不是对有限源辅助 current 的集体 Gram 等式/PSD 约束。约束位置、变量大小、材料依赖和成本应逐项比较。|
|**R2** Jing Tu & Edwin A. Marengo (2016), *Generalized likelihood ratio test change detection with optical theorem constraint*. [出版社](https://opg.optica.org/josaa/abstract.cfm?uri=josaa-33-11-2225), [DOI](https://doi.org/10.1364/JOSAA.33.002225)|**B，出版社摘要。** 把已知 passive/lossless/active 类型对应的光学定理约束并入 GLRT 的非线性规划，使用背景场信息；适用于一般 probing fields/backgrounds，并给检测例子。Crossref 已核作者、题名、2016。|这是**有约束散射数据估计/检测**的直接算法先例，不只是定理。没有取得完整公式页，不能声称与父线程矩阵式等价。它不重建共享逐点材料，也不是 CSI current 更新；不得把检测性能当材料恢复证据。|
|**R3** Gaku Umezu, Yoshihiro Yamauchi & Shouhei Kidera (2022), *Contrast Source Inversion Enhanced Confocal Imaging for Highly Heterogeneous Breast Media in Microwave Mammography*. [DOI](https://doi.org/10.1109/JERM.2022.3196504), [作者机构全文](https://www.ems.cei.uec.ac.jp/parts/article/2022/Contrast%20Source%20Inversion%20Enhanced%20Confocal%20Imaging%20for%20Highly%20Heterogeneous%20Breast%20Media%20in%20Microwave%20Mammography.pdf)|**A，全文。** IEEE pp495–496，§II，尤其 CSI cost Eq5 与其后的 Green-function 修改；`local_only/umezu2022.txt:169–232`。CSI 先重建 total field，再借互易性产生 Tx-to-cell / cell-to-Rx Green 函数用于 confocal 成像。|标题/正文确实同时出现 CSI 与 reciprocity，但**所读算法没有把集体互易条件加入 CSI current 子问题**。应视为相邻运用，不能错误宣称父线程方案已经被该文逐式实现。|
|**R4** Levinson & Markel (2016), *…T-matrix completion. I. Theory*. [DOI](https://doi.org/10.1103/PhysRevE.94.043317)|**A，前轮已读全文，今回重读。** §III Eq11 后及 Remark5，`local_only/dctmc2016.txt:265–285,372–390`：真实 T 在实空间对称，但截断得到的实验 T 通常不对称；迭代让 T 数据兼容、V 逐渐接近局部。|解释 R1 为什么互易“冗余数据”对 relaxed 变量仍能改变算法。该文当时只提一般对称性，并未给出父线程固定 E 上的 current 集体约束。不能由此断言这些条件充分刻画合法逐点材料。|

## 未升级为实质证据的线索

- Ottaviani & Pierotti, *Reconstruction of scattering data by the optical theorem*, DOI [10.1109/ULTSYM.1989.67122](https://doi.org/10.1109/ULTSYM.1989.67122)：Crossref 身份可核，发现阶段的第三方全文索引描述有限孔径散射数据重建约束；IEEE 页面仅返回机器人验证。故本轮仅 **C**，不拿它支撑精确算法主张。
- Taniguchi & Takayanagi (2018), *Inverse Scattering Problem and Transformations among Potentials*, DOI [10.7566/JPSCP.23.013010](https://doi.org/10.7566/JPSCP.23.013010)：Crossref 身份可核，原PDF403；第三方发现摘要涉及 generalized optical theorem 下的 half-on-shell T 变换。仅 **C**，不是本轮已核的 CSI/current 实现。
- 关于源/接收互换影响 CSI 输出的会议摘要、一般 optical theorem、统计散射和界限理论，没有升级成所问集体约束实现的证据。
- 一条未核实的 point-inhomogeneity 候选没有完成身份解析，已排除；不引用失败猜测的 arXiv 地址。

## 给主线程的具体对照项

以下仅是文献差异记录，非数学裁决：

1. 将候选写成 **固定照明子空间上的辅助 current 集体约束**，明确它不同于 R1 的逐点 polarizability 投影或 reciprocal 数据副本；不要使用“首次把互易/守恒用于逆散射”的宽泛说法。
2. 与“只补齐互易外部数据”及“只施加实材料/被动材料条件”区分，才能看出改善来自新增的 relaxed-current 约束还是已知处理。
3. 记录 off-diagonal source-pair 条件是否有额外作用；仅逐源功率标量与整个矩阵条件不应混称。
4. 必要的物理兼容条件不等于完整点态 constitutive law、可辨识性或 GT 证书。候选是否保持非平凡 lifting，以及加入条件后剩余自由度，留给父线程公式分析。

## 检索与收据

按要求只进行3条 web发现查询：`"contrast source inversion" "reciprocity"`；`"inverse scattering" "optical theorem" "constraint"`；`"contrast source" "energy conservation" inversion`。另沿**已识别的 DCTMC Part I→Part II**核对直接算法，不展开大范围引用网络。

ScholarQA单查询 paper search 仍HTTP429，exit1，零候选；失败保存于 `receipts/scholarqa_reciprocity_collect.json`。另保存 focused batch verify 的实际结果；独立身份由 Crossref/作者 arXiv/出版社核对。原始检索页面与第三方全文只在忽略的 `local_only/`；公开收据是元数据、查询和定位，不分发论文全文。

详细收据：`receipts/RECIPROCITY_SEARCH_LOG.json`、`receipts/reciprocity_downloads.json`、`receipts/reciprocity_pdf_locations.json`、相关 `crossref_*.json`。没有全球或二元新颖性标签。
