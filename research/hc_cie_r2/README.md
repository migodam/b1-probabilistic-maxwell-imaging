# HC-CIE-R2　有限观测下的结构消元与材料变化

[阅读中文图文报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/hc_cie_r2/) · [下一轮研究交接](CHAT_HANDOFF.md) · [逐案索引](RESULTS.md) · [代码、数据和运行方法](REPRODUCE.md)

576 格、12 个照明、48 个接收器时，完整曲率投影至多约束 6,336 个暗电流复坐标中的 66 个。另一种保留端口信息的做法，可以直接排除在指定错误重建上任选两个格修复的全部可能。在已知背景的两点变化监测中，二维和三维实算测出了这项结构对阵列与噪声的依赖。一般未知整图的同端点救援未增加材料恢复数。[研究报告原文](RESEARCH_REPORT_ZH.md)

|需要查什么|入口|
|---|---|
|代数证明、实材料条件、可见性与辐射|[理论](THEORY.md) · [与先例的区别](NOVELTY_AUDIT.md)|
|决定性小模型与成像规模上的检验|[严格局部反例](TOY_AND_CERTIFICATES.md) · [固定谱](SPECTRAL_RESULTS.md) · [实性方程](REAL_WEIGHTS.md) · [前端](REAL_FRONTEND.md)|
|错误端点、负曲率与实际恢复|[两格排除与变化恢复](PAIR_RESULTS.md) · [同端点干预](EVIDENCE_IMAGING.md) · [局部动力学](LANDSCAPE.md)|
|跨网格、实物、三维|[网格与噪声](MODEL_GAP.md) · [Austria 与 Fresnel](CLASSICAL_AND_MEASURED.md) · [双极化三维](VECTOR_CHANGE.md)|
|实验如何选择与计费|[冻结合同](EXPERIMENT_PROTOCOL.md) · [算法](ALGORITHM.md) · [成本](COST.md) · [研究经过](RESEARCH_LOG.md)|
|哪些方向被否定、下一项检验|[负结果](NEGATIVE_RESULTS.md) · [下一项研究](NEXT_RESEARCH.md) · [英文研究短稿](PAPER_DRAFT.md)|

本轮源码、数据表与图保存在本目录；历史输入和逐候选数组可作为同版附件下载。原 Institut Fresnel 测量按原作者渠道取得。[研究附件](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/hc-cie-r2-20261128)
