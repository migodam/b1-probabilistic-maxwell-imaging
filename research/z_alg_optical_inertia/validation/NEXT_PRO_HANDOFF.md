# 给 ChatGPT 的交接：峰值下界已跨材料验证

[22 页英文修订稿](manuscript_ip.pdf) · [中文报告](RESEARCH_REPORT_ZH.md) · [一键下载完整研究包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-optical-inertia-validation-v1/Z_ALG_Optical_Inertia_Validation.zip) · [平铺 LaTeX 源包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-optical-inertia-validation-v1/Z_ALG_Optical_Inertia_LaTeX.zip)

论文目标是用不完整的多次照明散射读数，排除“所有位置的介电对比度都低于给定阈值”，无需先反演一幅图。研究对象为完整非线性二维标量 TM。固定已知阵列、背景和误差半径。一个结果能证明峰值下界，并不能同时给出位置、图像误差或初始化方式。

## 投稿时可以说到哪里

1. 主定理沿用正谱中扣除状态负方向数的残差下界。既有六件材料的原生 96² 证明、两条全局单源约束的严格分离、四块与十六块空间单源的强弱变化均见英文稿。新的数值验证保留这套定理和旧几何，检查能否换材料。
2. 在选定的全部七份 R14 24² 同阵列观测上，直接复用 R10 环的 $u=5,\mu=10^{-3},\lambda=0$。七件都落在各自原 CN99 球内，其中 5 件严格排除上限 5。6 件真峰值超过 5 的物体中漏了 0009（真峰 16.25）。旧权重 $\nu=4$，任何固定相干单束取同权时下界为零。搜索过的权重中仍为 0/7 单束成功；连续全族的可行性见证未通过，不能把旧稿的连续分离定理统计扩成五件。 [原权重与逐件区间](max_struct/REPORT_ZH.md) · [连续权重的未决项](max_struct_single/README.md)
3. 在 SRAP 另 22 件原 10% 带噪的全环观测上，原生 96² 对上限 0.5 是 22/22 联合、20/22 同权最佳相干单束。0.75 的五种原权重为 0/22，其严格证书上界最高只达半径的 16.81%。004/014 的 5 个历史 SRAP-CIE 终点作为像素均匀的介质图都被排除；其原粗网投影本身也有接近 100% 的残差。77 个 MAX-STRUCT 历史终点都高于各自认证下界。新的独立证据是排除整个低峰值材料盒。 [区间](srap/README_ZH.md) · [原始历史终点对照](reconstructions/README_ZH.md)
4. Manitoba EIGOR 另有四件实际介质，7 个频率、288 个成对复数端口。定义 $R=\max_j\|Y_j-Y_0\|_F$，四个 $Y_j$ 来自两次空场与两个 PEC 柱交叉校准。在允许被动损耗的扇区合同下，主分析 $d_0-R>0$ 的数目为 20/28、5/28、3/28、1/28（上限 0.25/0.5/0.75/1）。在 0.5 留下正余量的三件是大尼龙柱、双尼龙柱、聚乙烯 E。对 E 的 3.5 GHz 记录，余量是主场的 1.358%；另一件 4 GHz 双柱余量只有 0.0525%。互易平均后逐项通过集合相同。四种校准处理的凸包有严格保证；真实介质所在的物理三维误差集合未被覆盖。 [全部频点](transfer/README_ZH.md) · [University of Manitoba 数据来源](https://home.cc.umanitoba.ca/~lovetrij/EMILab/mwt_data.html)
5. 对有损材料，证的是 $|\chi|^2/\operatorname{Re}\chi$ 的峰值。另给 $|\operatorname{Im}\chi|\le\kappa\operatorname{Re}\chi$，才能按 $1+\kappa^2$ 把它转回实部峰值。公开材料的典型损耗值不等于所有频率与位置的独立上界。 [推导与反例](transfer/README_ZH.md)

## 留给作者核稿的三个判断

- 理论定位宜放在**共享材料的跨照明不定惯性 + 同数据单束可行性见证**。实部/光学守恒、区域约束、正负惯性及校准误差的三角不等式各有既有数学基础。扩大结果说明整条方法具有迁移能力，没有增加一个新的电磁定理。
- 真实暗室最实在的贡献，是把每个材料结论尚可容忍的三维及天线误差量化出来。若要把特定数字写作实物峰值，需有与加载介质对应的独立误差上界。已有档案的四种校准和互易性只能约束它们观测到的误差部分。
- 文章可以进入作者核稿。当前尚未向期刊投稿。请保留 0009 的漏检、0.75 在全部 22 件的失败、4 块和 16 块强单束的区别，以及 Manitoba 在 5.5 GHz 四件均未留下正余量的结果。 [前一版论文与附录](../submission/index.html) · [复现](REPRODUCE.md)
