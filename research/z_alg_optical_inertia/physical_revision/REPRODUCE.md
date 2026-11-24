# 复算与数据获取

本轮证据包解压后直接保留 `research/` 层级。以下命令在解压目录根执行。合成观测、原生 \(96^2\) Fourier 核、六件参考材料、预先选择的乘子、源代码和结果均在包内。完整复算需要 Python 3.13、NumPy、SciPy、python-flint、threadpoolctl 和 Matplotlib。英文图使用 Matplotlib 自带字体，中文电尺寸和反演图另需系统安装 CJK 字体，例如 Noto Sans CJK。用一份解压副本运行会写入结果与图件。

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
rev=research/z_alg_optical_inertia/physical_revision
~~~

## 合成数据与图

| 研究任务 | 已有输入和输出 | 复算入口 |
|---|---|---|
| 阵列扫描 | 96 份复场和半径，四档发射数与两档接收数，480 条逐例结果，有符号源谱与区间 | [逐条命令](array/ARRAY_SCAN_ZH.md#能量噪声与复算) |
| 负状态方向和暗电流 | 同一原生 96² 参考材料与旧全环读数，原生四模、FFT-TSOM 子空间、坐标数组 | [取数与独立核对](modes/REPORT_ZH.md) |
| Born 模型范围 | 全环原始强散射与三个弱对照，两类 Born 限制和固定权重 | [主结果](born/REPORT.md) · [弱对照及数组](born/controls/RESULTS.md) |
| 反演影响 | 原阵列 18 份观测，八个有资格的拟合比较，另一个电紧凑物体及四个噪声球 | [冻结规则](inversion/PROTOCOL.md) · [结果与代价](inversion/RESULTS_ZH.md) |
| 电尺寸、角向采样和暗室 | 频率与阵列统计、已发表对照、派生最短距离和经验校准量 | [图表解释](electrical/FIGURE_NOTES.md) · [暗室复算](#暗室数据的本地复算) |

只重画已有图件时，可直接运行下列命令。原生阵列、模式和反演的数值结果已经随包提供。

~~~sh
python -B "$rev/array/summarize.py"
python -B "$rev/modes/plot.py"
python -B "$rev/inversion/plot.py"
python -B "$rev/electrical/plot_electrical.py"
~~~

下列命令重放各原生模型。阵列分支先使用包内固定数据，也可按其复算说明先从六件参考材料重新生成新阵列读数。附加合成反演用正常入口，离线评估只在预测产生后读参考图和留出端口。

~~~sh
python -B "$rev/modes/analyze.py"
python -B "$rev/born/run.py"
python -B "$rev/born/extinction.py"
python -B "$rev/born/controls/replay.py"
python -B "$rev/electrical/angular.py"
python -B "$rev/inversion/gate.py"
python -B "$rev/inversion/fit.py" --out "$rev/inversion/replayed_fits"
python -B "$rev/inversion/score_original.py" --input "$rev/inversion/replayed_fits"
python -B "$rev/inversion/compact.py"
python -B "$rev/inversion/compact_followup.py"
~~~

## 来源与历史入口

- 六件材料的原始对象、公开观测与评分真值位于 `research/z_alg_optical_inertia/synthetic/{public,evaluation}/`。新阵列、负状态模态和原件反演引用的是同一批材料。本包带齐合成输入与原生 Green 乘子。重做原始数据生成另见先前交接中的依赖和数据生成入口。Annular inclusion、Austria 的固定仿真算子位于 `research/z_alg_secant_tasks/highdim/public_vie_inputs.npz`。 [最初光学惯性证据包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-optical-inertia-20261122/Z_ALG_Optical_Inertia_2026-11-22.zip)
- 原始全环噪声曲线、外向距离、4×4 权重和校准容差出自既有发布。需要的 冻结 FFT 算子、原始 Y、区间工具与权重依据保存在 `results/theory/`、`submission/curve/`、`submission/single/`、`core/` 和 `validation/srap/online.py`。Fresnel 适配器的公开 R8 源码保留在 `research/delegated/max_struct_r1/r8_worktree/`，早期 FFT-TSOM 基函数在 R14 冻结源代码路径。 [原全环噪声曲线和权重出处](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-optical-inertia-submission-20261122/Z_ALG_Optical_Inertia_Submission_2026-11-22.zip)
- 跨实验暗室表保留原先发表的校准误差标尺、选定的残差界与资料出处。旧 4 区/16 区单源对照完整材料继续沿用第二版投稿页。 [实测校准与历史结果](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-optical-inertia-validation-v1/Z_ALG_Optical_Inertia_Validation.zip) · [原单源和分区证据](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/z_alg_optical_inertia/ip_submission_v2/)

## 暗室数据的本地复算

已发布的标量表、图和结果用于离线重画。完整实测复算需要从采集方获取复读数，采用原数据使用条件，并在本地生成校准依赖。两种实测均比较经验误差尺度，完整结果仍受采集和模型误差合同制约。

1. **University of Manitoba.** 将官方 EIGOR `Tomographic_Dataset_1.zip` 放到 `.artifacts/z_alg_validation/manitoba/`。以包内的 `validation/manitoba/physical.py` 重做原校准并输出本地 `prepared_private.npz`。 [采集方数据页](https://home.cc.umanitoba.ca/~lovetrij/EMILab/mwt_data.html) · [采集论文](https://doi.org/10.1109/MAP.2011.6138442)
2. **Institut Fresnel.** 从 2005 年 S09 论文的 `supp.zip` 提取 `FoamDielIntTM.exp`、`FoamDielExtTM.exp` 和 `FoamTwinDielTM.exp`。保存至 `research/delegated/a2_literature/data/2005_iop_21_6_S09/`。以下命令经随包提供的原 R8 加载器生成本地 `mask` 和 `raw_scattered` 文件。 [采集论文](https://doi.org/10.1088/0266-5611/21/6/S09) · [出版社补充数据](https://beta.iopscience.iop.org/article/10.1088/0266-5611/21/6/S09/data)

~~~sh
study=research/z_alg_optical_inertia
python -B "$study/validation/manitoba/physical.py"
python -B - <<'PY'
import numpy as np
from research.z_alg_optical_inertia.fresnel import physical
for frequency in physical.FREQUENCIES:
    for name in physical.NAMES:
        rows, measured, _ = physical.load_file(physical.RAW / (name + ".exp"), frequency)
        transmitters = np.unique(rows[:, 0]).astype(int) - 1
        field = measured[:, transmitters]
        np.savez_compressed(
            physical.RESULTS / f"{name}_{frequency}GHz_ports.npz",
            mask=np.isfinite(field), raw_scattered=field,
        )
PY
python -B "$rev/electrical/scale.py"
python -B "$rev/born/ports/replay.py"
~~~

公开 ZIP 的逐项路径在 [归档入口](package_release.py)。实测原始文件、本地校准的端口缓存，以及从实测空场拟合得到的复权重留在本地。可以读取包内的 `born/ports/results.json` 重查经验标尺比较。
