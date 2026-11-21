# 交接 · 已知非均匀宿主中的材料折叠

请阅读中文报告、文献对照与机器证书。讨论时分别跟踪折叠处正演的可逆性、整个材料盒沿指定观测线的 **2 → 1 → 0** 根数，以及同一个噪声球中的材料距离和路径。后续比较沿用固定的原模型和误差合同。

[中文报告与图表](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/z_alg_host_fold/) · [命题级文献对照](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/z-alg-host-fold-20261121/research/z_alg_host_fold/PRECEDENTS.md) · [固定版本研究包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-host-fold-20261121/Z_ALG_host_fold_2026-11-21.zip)

## 读数和材料

这是一套固定的二维标量 TM 圆盘脉冲模型。36 个状态单元中，30 个组成已知非均匀宿主，6 个已知位置上的实介电对比度未知。六个材料独立落在 $[0.5,30]$。三源一收返回三条复读数。宿主消元精确保留该离散模型的多次散射，观测在同一个几何、宿主与频率下比较。[模型、坐标与参数](results/host_model.json) · [冻结算子与原始读数](results/host_operators.npz) · [解析离散算子的区间包围](results/host_operators_interval.json)

## 可使用的证据

| 问题 | 本研究结论 |
|---|---|
| 折叠处状态方程会失去唯一性吗？ | 完整 36 阶正演在折叠认证区间的最小奇异值严格大于 **0.26124**。全部噪声可行材料的完整正演最小奇异值都大于 **0.25260**。另证整个合法材料盒的正演可逆。 |
| 材料到测量的局部映射发生什么？ | 指定材料线段上认证一个单零方向的普通 fold。左右零方向的二阶配对离开零，回读所需场量具有严格非零余量。 |
| 两个固定读数各有多少合法材料？ | 在整个 $[0.5,30]^6$ 中，$y_-$ 恰有 **2** 份，$y_+$ 恰有 **0** 份。两套 $20$ 个复根均经过认证并按物理合法性分类。用同样系数跑 HomotopyContinuation.jl，两端分别找到 20/20 个，得到相同分类。 |

[英文论文稿](manuscript.pdf) · [全部噪声可行材料的正演余量](results/host_noise_margin.json) · [折叠与状态奇异值](results/fold_certificate.json) · [全盒正演可逆性](results/host_state_certificate.json) · [全盒逐观测根数](results/point_homotopy_certificate.json) · [独立求根记录](results/point_fibers.json)

## 分支与带噪范围

<!-- HOST_FOLD_HANDOFF_SCOPE -->
全盒根数在两端之间也经验证。所给读数线长约 $2\times10^{-6}$。对线上的每一个读数，临界读数前恰有两根，临界读数恰有一根，临界读数后恰无根。55 段经认证的局部区间和较强直接法外包的 222 个材料盒形成 1,306 个已证交集，覆盖全盒在该线上仍可能有解的区域。结论是固定几何中指定一参数截线的完整根数。[全盒截线证书](results/normal_section.json) · [验证代码](src/normal_section.py)

以固定中心取三条复读数整体的欧氏噪声球，半径为 $5\cdot2^{-22}\approx1.19209\times10^{-6}$。两份 $y_-$ 材料在同一个球的逆像内，RMS 距离约为 **0.22922**。任何给同一个带噪读数返回一份材料的估计器，对这两份候选的最坏误差至少为 **0.11461 RMS**。另证两者在球内可沿连续材料路径连接。这个连接证书覆盖这两个候选之间的路径。[距离、连接及严格内点](results/noisy_diagnostics.json)

同一全盒、相同保留体积分数 $10^{-6}$ 时，中心化直接法与二次图覆盖分别用 **0.1237 秒**和 **0.0298 秒**，共享宿主预处理单列。达到共同最大细分深度时，直接法留下 **222** 个盒，二次图留下 **345** 个盒。一个正体积小盒通过了二次图的十五项逐项必要检验，完整的联合噪声检验却严格排除整盒。这说明更快的外包可以保留一些实际不相容的材料。[覆盖成本与外包前沿](results/noisy_cover.json) · [严格排除的小盒](results/noisy_diagnostics.json)
<!-- /HOST_FOLD_HANDOFF_SCOPE -->

[覆盖记录](results/noisy_cover.json) · [覆盖盒数组](results/noisy_cover_cells.npz)

## 继续讨论的边界

现有文献已给出全阶传递矩阵、根枚举、常规 fold 认证、材料比值回读和有界误差集合反演。本研究继承 R6 的实二次兼容推导，并用成熟的同伦、区间包围与逐盒排除核查已知非均匀宿主中的具体模型。上述证据条件限定在六个未知实材料、指定三条复读数以及固定有限离散。讨论扩展到换几何、未知宿主、更多材料或连续 Maxwell 方程时，先说明哪些物理假设和区间证据需要重建。[原论文公式和逐项比较](PRECEDENTS.md)
