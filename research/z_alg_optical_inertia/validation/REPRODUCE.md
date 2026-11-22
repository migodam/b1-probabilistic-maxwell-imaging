# 复算扩大验证

[完整科研包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-optical-inertia-validation-v1/Z_ALG_Optical_Inertia_Validation.zip) · [英文 LaTeX 源包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-optical-inertia-validation-v1/Z_ALG_Optical_Inertia_LaTeX.zip)

使用 Python、NumPy、SciPy、python-flint、threadpoolctl 和 Matplotlib。完整科研包中先前论文的几何算子、区间计算和源文件按原路径放置。下列入口位于 research/z_alg_optical_inertia/validation。机器数据记录实跑墙钟和 CPU 时间。

| 实验 | 可复算的内容 |
|---|---|
| SRAP，全环 22 件 | 从 srap 目录执行 python -B online.py --public-bundle public_observations.npz。可从独立原 SRAP 归档复算真材料覆盖；本包内保留原次逐件评分。[数据来源与复算说明](srap/README_ZH.md) |
| MAX-STRUCT，R14 七件 | 从 max_struct 依次运行 scan.py、verify.py 和 summarize.py。原合成真材料受项目封存约定，公开文件包含可独立复算在线证书的全部复场及原次离线覆盖标量。[复算说明](max_struct/REPRODUCE.md) |
| 旧重建对照 | 从已发布的两套历史结果读取 253 个终点的峰值和原残差。本包公开逐条数值与五个均匀细化的核对结果。[来源与命令](reconstructions/README_ZH.md) |
| Manitoba 实测 | 从大学取得第 1 号层析原始 ZIP，保持原文件名。physical.py 先生成校准后的本地数据，scan.py 做有限网格搜索及选中乘子的区间核对，summarize.py 出表。[下载、参数及复算命令](manitoba/REPORT_ZH.md) · [官方数据](https://home.cc.umanitoba.ca/~lovetrij/EMILab/mwt_data.html) |
| 四种校准和互易性 | 从 transfer 运行 calibration_hull.py，参数 --raw 和 --prepared 可指向原 ZIP 和 physical.py 的本地输出。运行 manitoba/reciprocal_replay.py 回放互易投影，再从 transfer 执行 margins.py 和 plot.py。各选定通道的复读数只在本地处理，公开严格余量。[模型假设与逐频表](transfer/README_ZH.md) |
| R14 连续单束对照的跟进 | 从 max_struct_single 运行 global_witness.py。保留七个找不到符号见证的案例。[适用范围](max_struct_single/README.md) |

英文修订稿沿用 IOP 期刊类。解开单层 LaTeX 源包，在其中用 tectonic -X compile manuscript_ip.tex 或 latexmk -pdf manuscript_ip.tex。公开报告沿用同一 B1 论文主页，原有六件实验和冻结版 PDF 可单独阅读。 [旧稿](../submission/index.html)
