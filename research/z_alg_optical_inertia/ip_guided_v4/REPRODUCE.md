# 数值证据和复现

将 **Z_ALG_Guided_Evidence.zip** 解压到空目录。下列命令从解压根运行。包内保留本轮原复场、材料候选、严格端点、绘图数据、结构化失败记录及费用。所有本轮仿真和 A13 推理所需的历史输入都按原相对路径提供。暗室原复场通过原采集方取得。重新计算会写入解压副本的对应目录，保留另一个副本可随时比较原回执。论文 PDF 和匿名稿件源码为独立附件。

~~~sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install numpy 'scipy>=1.15' python-flint threadpoolctl matplotlib
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
study=research/z_alg_optical_inertia
guided="$study/ip_guided_v4"
pool="$guided/pool"
routing="$guided/routing"
~~~

## 先从存档重算证据

| 问题 | 数据和回执 |
|---|---|
| 非交换背景、增广矩阵和奇异情形 | [两单元区间形式检查](theory/check_forms.json)包括背景场、接收、带边同余、无损光学和有损旋转。另有一单元奇异块的精确有理数对照。 |
| 同原算子的局部盒排除、盘体宽度及两份材料的非唯一性 | [逐例配对表](boxes/results/paired_summary.csv)，[原生局部盒](boxes/results/verified.json)，[盘体物理宽度](boxes/results/disk_radius.json)，[GENERAL6_00](boxes/results/ambiguity.json)和 [离散输入及回放数组](boxes/results/ambiguity.npz)。 |
| 一个源阈值能否筛掉材料盒 | [四件开发](threshold/development.json)，[八件独立确认](threshold/confirmation/offline.json)，同级 public、scoring、online.json 和 predictions.npz 保留在线/离线阶段的独立输入输出。 |
| 两个分支方法用相同停机准则 | [配对说明](threshold/same_stop.md) · [八件各次运行及费用](threshold/same_stop.json)。两路共享确认观测、盒、精修预算和首个经区间核对的可行点停机准则。 |
| 背景格也可变的材料族 | [两处局部假设与六件原生读数](families/results.json)、[逐级背景上限](families/nested.json)和[同位置低峰对照](families/same_site.json)。两份 host_patch 原读数和两份低峰新读数同时保留原复端口、半径与离线参考。 |
| 有界精修、排序和分支裁剪的真实费用 | [24² 全链费用](pool/structured_results/timed_pipelines.json)，[保留原生 A13 料盒证据](pool/flow_results/native_method_paths.json)，同级 JSON/NPZ 记录生成、优化、严格端点、未决盒和离线留出。 |
| 同一合成测量的纯 CIE、共同校正与直接全波 | [逐件配对表](pool/cie_pair/per_object.csv)、[纯 CIE](pool/cie_pair/results/online.json)、[CIE 后统一物理校正](pool/cie_pair/corrected_results/online.json)、[直接全波](pool/cie_pair/direct_results/online.json)。每组配有预测数组、离线端口和用时。另保留[直接精修追加预算诊断](pool/cie_pair/direct_results/diagnostic.json)。 |
| Born 和匹配全波，含校准实测 | [全部配对臂](routing/paired_arms.csv)，[主对象 28 件](routing/principal_28.csv)，[端点、迭代和停止](routing/summary.json)，[实测 37 份](routing/measured_port_cases.csv)。online/offline 和 CIE 各子目录分别保留训练复场、候选、留出端口和全部预算。 |
| 数据计数和复现依赖 | [对象及观测目录](benchmark/objects.csv)。包内 benchmark/dependency_manifest.json 逐文件说明归档输入、截取的 NPZ 键、A13 行和实测数据来源。 |

以下命令直接重读已封存的证据，无需先生成新材料或重新求候选。A13 的独立区间核查含 1024 维精确二进制 Green 系统，用时可高于重画图和合并标量表。逐审计费用见对应 JSON。

~~~sh
python -B "$guided/theory/check_forms.py"
python -B "$guided/boxes/ambiguity.py" --replay-npz "$guided/boxes/results/ambiguity.npz"
python -B "$guided/boxes/disk_radius.py"
python -B "$guided/boxes/paired_summary.py"
python -B "$guided/boxes/plot.py"
python -B "$guided/threshold/confirm.py" score
python -B "$guided/families/nested.py"
python -B "$guided/families/same_site.py"
python -B "$pool/replay.py" --support-only
python -B "$pool/flow_verify.py" all
python -B "$pool/flow_shift_strict.py"
python -B "$pool/flow_neumann_ideal.py"
python -B "$pool/flow_ambiguity.py"
python -B "$pool/flow_arms.py"
python -B "$pool/plot_pool.py"
python -B "$pool/cie_pair/summarize.py"
python -B "$routing/summarize.py"
python -B "$routing/plot.py"
python -B "$guided/benchmark/build.py"
~~~

## 从观测开始重新计算

### 局部盒和背景材料族

原生 96² 使用 v3 保存的 192² Fourier 符号、原发收矩阵和原复读数。24² 环体读取旧 R10 存档的对应观测，32² 独立样本读取同格的原 Green 和发收算子。已知低峰盘仅用于独立离线见证。GENERAL6_00 可重建 R2 几何或直接读取本轮证明使用的精确二进制数组。[局部盒冻结](boxes/FREEZE.md) · [非真空背景材料族](families/DESIGN.md)

~~~sh
python -B "$guided/boxes/run.py"
python -B "$guided/boxes/verify.py" --spectral --dark
for arm in broad disk host; do python -B "$guided/boxes/followup.py" "$arm"; done
python -B "$guided/boxes/disk_radius.py"
python -B "$guided/boxes/ambiguity.py"
python -B "$guided/families/run.py"
python -B "$guided/families/nested.py"
python -B "$guided/families/same_site.py"
~~~

### 合成材料的排序和前瞻阈值

在新输出根可完整重做合成实验生成、在线精修和离线评分。阈值扫描与八件确认读取存档的原合成实验几何。已封存的八个确认观测直接用于在线和离线复算。若要重做生成，先另存解压副本的 threshold/confirmation 目录，再在空的对应路径运行 confirm.py generate。

~~~sh
export GUIDED_POOL_OUTPUT_ROOT="$PWD/replay_pool"
python -B "$pool/structured_pool.py" generate
python -B "$pool/structured_pool.py" online
python -B "$pool/structured_pool.py" offline
python -B "$pool/replay.py"
python -B "$pool/replay.py" --offline
python -B "$pool/replay.py" --support-only
unset GUIDED_POOL_OUTPUT_ROOT
(cd "$guided/threshold" && python -B run.py)
python -B "$guided/threshold/confirm.py" online
python -B "$guided/threshold/confirm.py" score
~~~

对已曝光的八件确认对象重做同停机准则对照，使用 confirmation/public 内的冻结观测和同级 online.json、predictions.npz。默认每件十次，平衡两条路径的执行顺序。逐盒独立区间复核和离线评分包含在同一个回执中。

~~~sh
python -B "$guided/threshold/same_stop.py" --repetitions 10
~~~

### A13 原生推理和盒内续算

另在兼容环境装 PyTorch。原运行使用 CPU 推理，保留 seed 11 的 flow 和 probe 两个完整推理检查点、latent_stats、训练说明和其 32² 物理实现。A13 的 y.npy、chi.npy 保持原索引，三个本轮对象位于 0、37、111，其余行置空为 NaN。本轮在线推理只读三个对象的 16 个拟合端口，已封存结果中的另 16 个端口和对应材料用于离线评分。

~~~sh
python -m pip install torch
export GUIDED_POOL_OUTPUT_ROOT="$PWD/replay_pool"
python -B "$pool/flow_audit.py" online
python -B "$pool/flow_audit.py" offline
python -B "$pool/flow_continue.py"
python -B "$pool/flow_continue.py" --offline
python -B "$pool/flow_verify.py" all
python -B "$pool/flow_screen.py"
python -B "$pool/flow_shift_strict.py"
python -B "$pool/flow_neumann_ideal.py"
python -B "$pool/flow_ambiguity.py"
python -B "$pool/flow_arms.py"
unset GUIDED_POOL_OUTPUT_ROOT
~~~

### 匹配的 Born、原生全波和 CIE

合成分支可直接重跑。离线评分读已公开的 optfresh 材料和 SRAP 再生成对照、R14 路径下的七件配对合成材料、无损/有损几何和原复接收端口。原合成评分数组在 online 产物落盘后读取。详见 [原配对命令和误差合同](routing/REPORT_ZH.md#复现)。

~~~sh
for suite in native ring loss weak mild_mid; do
  python -B "$routing/routing_run.py" --suite "$suite"
done
python -B "$routing/cie.py" --branch principal
python -B "$routing/cie.py" --branch principal --box4 --out "$routing/cie_principal_common_box"
python -B "$routing/cie.py" --branch diagnostic
python -B "$routing/cie.py" --branch diagnostic --diagnostic-extra
python -B "$routing/cie.py" --branch principal --box4 --steps 96 --full-schedule --objects optfresh_01 optfresh_05 --out "$routing/cie_principal_longer"
python -B "$routing/early_exit.py"
python -B "$routing/early_exit.py" --halfwidth .005
python -B "$routing/score.py" --suite native --with-cie
for suite in ring loss weak mild_mid diagnostic; do
  python -B "$routing/score.py" --suite "$suite"
done
python -B "$routing/summarize.py"
~~~

### 同一 24² 测量的 CIE 配对

以下命令重读四件开发和八件确认的冻结观测、候选条件及原 SRAP/CIE 算子，分别执行原 CIE、加共同物理校正的 CIE 和独立直接全波校正。在线预测先保存，真值及留出端口只在离线评分时用到。重算可用空的新目录。直接精修追加预算诊断读取同一个已曝光失败对象，脚本写回解压副本中保存诊断的目录。 [合同和逐组耗时](pool/cie_pair/REPORT_ZH.md)

~~~sh
pair="$pool/cie_pair"
python -B "$pair/paired.py" online --out "$PWD/replay_cie_only"
python -B "$pair/paired.py" score --out "$PWD/replay_cie_only"
python -B "$pair/paired.py" online --physical-correction --out "$PWD/replay_cie_corrected"
python -B "$pair/paired.py" score --out "$PWD/replay_cie_corrected"
python -B "$pair/direct.py" online --out "$PWD/replay_cie_direct"
python -B "$pair/direct.py" score --out "$PWD/replay_cie_direct"
python -B "$pair/direct_diagnostic.py"
~~~

### 两所暗室的原数据入口

包中可直接读取本轮实测拟合的材料、预测复端口、三种端口掩码、误差表和校准时钟。获取原数据后可重新拟合。数据文件及由其生成的原接收和空场缓存保存在本地。

- University of Manitoba。将 Tomographic_Dataset_1.zip 放在 .artifacts/z_alg_validation/manitoba/，运行 validation/manitoba/physical.py 生成本地 prepared_private.npz。 [原数据入口](https://home.cc.umanitoba.ca/~lovetrij/EMILab/mwt_data.html) · [采集论文](https://doi.org/10.1109/MAP.2011.6138442)
- Institut Fresnel。从 2005 年 S09 补充文件提取 FoamDielIntTM.exp、FoamDielExtTM.exp、FoamTwinDielTM.exp，放在 research/delegated/a2_literature/data/2005_iop_21_6_S09/。包中保留原 R8 加载器。 [采集论文](https://doi.org/10.1088/0266-5611/21/6/S09) · [出版社数据入口](https://beta.iopscience.iop.org/article/10.1088/0266-5611/21/6/S09/data)

原数据到位后，运行下列命令。Fresnel 的两个留出掩码按原几何代码中的规则生成，缓存带原值，仅留本地。

~~~sh
python -B "$study/validation/manitoba/physical.py"
python -B - <<'PY'
import numpy as np
from research.z_alg_optical_inertia.fresnel import physical

for frequency in physical.FREQUENCIES:
    for name in physical.NAMES:
        rows, measured, _ = physical.load_file(physical.RAW / (name + '.exp'), frequency)
        tx = np.unique(rows[:, 0]).astype(int) - 1
        field = measured[:, tx]
        mask = np.isfinite(field)
        source_degrees = np.arange(len(tx)) * (360 // len(tx))
        relative_angle = (np.arange(360)[:, None] - source_degrees) % 360
        held_tx = mask & (tx % 4 == 0)[None, :]
        held_rx = mask & ~held_tx & ((relative_angle - 60) % 5 == 0)
        np.savez_compressed(physical.RESULTS / f'{name}_{frequency}GHz_ports.npz',
                            mask=mask, raw_scattered=field,
                            nominal_held_tx=held_tx, nominal_held_rx=held_rx)
PY
python -B "$routing/measured.py"
~~~

## 历史材料和本包的边界

主要旧数据已通过原项目发布。本包抽取实际被本轮脚本访问的文件和数组。历史 u100_best.npz 内仅提供本轮 CIE 诊断读取的 negative_vec、negative_eig、center；完整历史谱证据仍沿用原发布。A13 权重来自工作区本地项目资产。两份 FLOW 模型及第 0、37、111 条已评分材料已随本包提供。

- 原生区间工具、有损合同和暗模数组。 [统一论文证据包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/optical-inertia-unified)
- 原生全环阵列、Born 权重、校准标量记录。 [v3 物理证据包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/optical-inertia-ip-v3)
- R10 环形观测、R14 七件原读数、24² 紧凑域与 SRAP 的历史验证。 [验证发布](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/z-alg-optical-inertia-validation-v1)
- 原生全波与 CIE 求解和历史 SRAP/R8 算子来源。 [CIE 研究与复现发布](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/z-alg-xudong-cie-20261124)
- GENERAL6 六单元合同、Q20 原解和独立路径记录。 [R2 原报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r2/)

重建发布 ZIP 的入口为 benchmark/dependencies.py，使用完整研究工作区已有归档输入。该命令用于从完整源工作区重新打包。用于计算复现的命令见前文。

~~~sh
python -B "$guided/benchmark/dependencies.py" --output Z_ALG_Guided_Evidence.zip
~~~
