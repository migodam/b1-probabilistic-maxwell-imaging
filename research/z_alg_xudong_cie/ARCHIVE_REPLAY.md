# Track B 归档与复跑

复现包按仓库根路径展开，包含 Track B 源码、报告、算子、合成观测、预测、计时和相关历史依赖。先取同一 B1 研究仓库，在仓库根展开发布附件 `Z_ALG_Xudong_CIE_Research.zip`。[项目仓库](https://github.com/migodam/b1-probabilistic-maxwell-imaging) · [本轮发布](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/z-alg-xudong-cie-20261124)

```sh
unzip /path/to/Z_ALG_Xudong_CIE_Research.zip -d .
python -m pip install numpy scipy threadpoolctl python-flint matplotlib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
study=research/z_alg_xudong_cie
```

## 包中各类输入

| 用途 | 目录及边界 |
|---|---|
| B3 合成训练与原生复核 | `.artifacts/z_alg_xudong_cie/b3_holdout/public/` 含原 32²、96² 传播频谱和端口、48 份成对训练观测及半径。离线参考材料和隐藏观测放在同任务的 `offline/`，训练加载器只读 `public/`。 |
| 网格可表示性 | `holdout/projection/` 含全部真材块平均及其原生复场、固定算例恢复到噪声球内的 32² 材料数组与区间复核；这些数据只在评分后打开。 |
| 追加逐格相位检验 | `mechanism/acquisition/public/` 与同级 `offline/` 保存四件全新材料、两种孔径和 2% 噪声的读数及原生参考；24 组输出在 `mechanism/predictions/`。 |
| 支持域剪枝对照 | `.artifacts/z_alg_xudong_cie/site_control/public/` 与 `offline/` 保存对应的冻结训练与评价数组。上限、分支、端口掩码和试验计时见 `pruning/`。 |
| B2 开发 | 包含三件合成的原训练数组、全部逐步指标与电流快照。实测部分保留作者产生的材料图、标量结果、端口来源与掩码定义。实测原读数、校准后的复数据及可反推复数据的电流/场快照应在本地复建。 |
| B1 原生与面积 | 六件原生材料的匹配全环读数、几何、公开区间证书，以及用于 R10 面积证书的原保存算子和读数。历史出厂参考的覆盖和面积可直接查阅本包所存区间计分标量。 |
| B1 七件 MAX-STRUCT | 包含对应 24² 公开算子、七份重复合同核对所需的原观测、候选重建材料和精简历史计费回执。已有真值评分以区间标量保留。 |
| 历史 Python 依赖 | 原 SRAP 的 `engine.py`、`legacy_la.py` 及实际调用的 R2/R6/R8/R10 物理与读出代码保留在原相对路径。Python 依赖从常用发布源安装。 |

机器可读文件清单 `research/z_alg_xudong_cie/archive_dependencies.json` 随 ZIP 生成。通过结果、未通过结果和真实计算费用保存在对应目录。使用历史档案时沿用其作者署名和许可。 [光学惯性验证包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/z-alg-optical-inertia-validation-v1) · [固定 Green 原公开算子](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/z-alg-secants-20261122) · [SRAP 归档](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/srap-cie-array-study-20261121) · [R8](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r8-20261106) · [R14](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r14-dualmetric-20261115)

## 从官方档案复建暗室输入

Fresnel 从 Geffrin、Sabouroux 与 Eyraud 的 2005 论文补充材料获取原 `supp.zip`，或将其中 `FoamTwinDielTM.exp`、`FoamDielIntTM.exp` 与 `FoamDielExtTM.exp` 解到同一目录。IOP 页面可能要求在浏览器完成正常下载流程。只需这三个 TM 文件。复建采用原空场拟合和发射留出掩码。 [Fresnel 采集论文](https://doi.org/10.1088/0266-5611/21/6/S09) · [IOP 补充数据](https://iopscience.iop.org/article/10.1088/0266-5611/21/6/S09/data) · [Institut Fresnel 数据说明](https://www.fresnel.fr/3Ddatabase/)

Manitoba 从 University of Manitoba 官方档案取得第 1 号原始 ZIP。输入第一次空场和大金属柱进行端口标定，第二次空场和小柱仅给出单独校验。复建不改变训练、留出或有效通道。 [University of Manitoba 数据页](https://home.cc.umanitoba.ca/~lovetrij/EMILab/mwt_data.html) · [采集论文](https://doi.org/10.1109/MAP.2011.6138442)

```sh
mkdir -p .artifacts/z_alg_validation/manitoba
curl -fL 'https://home.cc.umanitoba.ca/~lovetrij/EMILab/research/EIGOR/Tomographic_Dataset_1.zip' \
  -o .artifacts/z_alg_validation/manitoba/Tomographic_Dataset_1.zip
python -B "$study/rehydrate_chambers.py" --fresnel /path/to/supp.zip \
  --manitoba .artifacts/z_alg_validation/manitoba/Tomographic_Dataset_1.zip
```

脚本将 Fresnel 训练与留出视图展开到历史 `.runs` 路径，Manitoba 的复数缓存留在本地 `.artifacts`。只复跑一间暗室时，单独传相应输入。脚本优先复用已有输入。 [本轮实测拆分和评分](chambers/README.md)

## 实际程序入口

只用合成输入即可执行单观测验证。以下命令写到该例的新目录。

```sh
python -B "$study/benchmark/run.py" --case b3m_0001__pw180_10x20__n02 \
  --method dual --output "$study/benchmark/replay_example" --workers 1
```

全部实验顺序、重新出图方法及离线评分见 [B 复现说明](REPRODUCE.md)。B1 可分别从 [原生 96²](minpeak_native/REPRODUCE_ZH.md)、[七件 R14](minpeak_maxstruct/REPRODUCE.md) 和 [面积证书](area/README.md) 开始。B2 可从 [电流诊断](diagnosis/REPRODUCE.md) 与 [光学修正的公式及开发对照](cross_cie/THEORY_AND_REPRODUCE.md) 开始。剪枝可从[支撑区条件](pruning/THEOREM.md)和[弱端口互易对照](pruning/RECIPROCITY.md)开始。复对称压缩基与见证数组随 ZIP 分发。归档阅读页是静态 HTML，生成阶段不需要 Chrome。
