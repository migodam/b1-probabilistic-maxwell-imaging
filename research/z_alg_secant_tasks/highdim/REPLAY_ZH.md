# 高维证据的独立复算

公开包含五份合成旧观测、实际候选材料、完整离散物理矩阵及全部求解入口。`public_vie_inputs.npz` 为 2.6 MiB，其中的两个 `G_*` 是逐元素无损的完整 Green 矩阵，五份观测各有自己的 `E0` 入射场和 `S` 接收矩阵。每份观测用 `<ID>__G_key` 指向其 Green 矩阵。`results/<ID>.npz` 给出原带噪 `Y_archive`、区域任务向量及模型候选。等预算实验另有两个 `__comparison_arrays.npz`。这些输入包没有历史逐像素材料标签。

在解压后的项目包根目录运行，可用安装了 NumPy、SciPy 的 Python 3.10 或更新版本。

```sh
python highdim/replay.py --output replay.json
```

这一步只读取 `highdim/` 自己的文件。程序以复数双精度逐端求解 `A(m)j=diag(m)E0`，输出 45 对材料的实际任务差、正负端距原观测和中心场、端点之间场距、噪声球半径与合格情况。它也复算两例等预算比较的每个保存的候选与可行跨度。脚本比较复算场和封存场，并验证有明确区间判定的噪声球合格情况。新输出是双精度重放，原算例的数值误差区间见 `results/<ID>.json`。

以下两步会从保存的公开历史中心重新计算所有方向、两侧全波端点、数值区间和配对 GN。`recomputed/` 用于保留封存的旧候选，便于比对。

```sh
python highdim/run.py --output recomputed
python highdim/compare.py --output recomputed
```

独立包会自动使用保存的物理输入。在同时保有历史工作树的环境，设置 `HIGHDIM_USE_FROZEN=1` 可指定使用同一包内输入。历史目录位置若偏离默认布局，可用 `HIGHDIM_WORKSPACE` 指向包含 `research/delegated/max_struct_r1` 和 `research/srap_cie_array_study_20261121` 的工作区。`archives.py` 用于从历史工作树导出包内全物理矩阵，独立复算无需运行它。

历史研究的原始生成、评分和完整成本账本见下列入口。重做本轮中方向、GN 与数据球计算无需重新生成旧观测。`score.py` 用于有权限的工作区内真值评分，单独需要历史 R10 `sealed_truth/twod/`、R8 `data/paper_benchmarks_private/` 及 SRAP `data/truth/` 中的参考文件。公开包提供原已报告的离线评分和半径，不打包这些参考数组。

- R10 环形案例及困难起点追溯见 R12 补充运行包。 [MAX-STRUCT R12 配对补充实验](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r12-paired-20261112)
- Austria 原 150² 观测见 R8 发布的纸面基准附件。 [MAX-STRUCT R8 原观测与基准](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r8-20261106)
- SRAP 三份历史 96² 观测及配对端点从二维阵列发布包获取。 [SRAP-CIE 原观测、数组与说明](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/srap-cie-array-study-20261121)
