# R12 复现入口与已完成产物（初稿）

本文件是实际脚本/路径的导航，不触发新实验。完成结果保持不动；在新的代码 checkout 和空输出目录中复现，不能对当前已封存目录直接重跑。固定 HERE 输出的脚本会拒绝既有 runs，NN `prepare` 还会重写该副本的协议/初值，因此必须使用新副本。没有 SHA 检查要求。

## 1. 已验证的本地环境

```sh
cd "/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r12_worktree"
export R12_PY="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export MKL_NUM_THREADS=1
```

使用该实际 Python 路径，不解析或替换它的链接。本轮使用 CPU、NumPy/SciPy/PyTorch；未安装新依赖或 GPU 服务。命令中的相对路径均以 worktree 根为起点。新 checkout 必须保留这里列出的真实 sibling 数据/继承代码依赖，或显式将它们迁移到对应新位置；不得编造历史代码路径。

共享物理入口为 `research/research_fbg_r3_experiments/src/physics.py::MaxwellModel`、`research/research_max_struct_r8/src/vector_pr3d.py::VectorPR3D`、`research/research_max_struct_r8/src/vector_cie3d.py::VectorCIE3D` 和 `research/research_max_struct_r1/src/vector_backend.py::TensorDDA`。本轮 core/adapter 不替换这些 backend。

## 2. 双场 core 与已完成窄检查

源代码和 API：

- [src/field_current.py](src/field_current.py)：`DenseFieldCurrentProjectors`、`PhysicalAudit`、`run_projection`、完整 audit 存储。
- [src/README.md](src/README.md)：真实参数、每次 callback 的增量 cost 语义、布局与 budget 语义。
- [src/wiring_check_report.json](src/wiring_check_report.json)：已完成 scalar AP/DR 6步集成、两 cell 三分量共享材料、mask/whitener、不可达球、伪逆极限和零余量预算检查。
- [src/projection_grouping_benchmark.json](src/projection_grouping_benchmark.json)：已完成的单个固定 A workload v1/v2 吞吐记录。

在新 checkout 中，窄检查的实际入口为：

```sh
"$R12_PY" research/research_max_struct_r12/src/wiring_check.py
"$R12_PY" research/research_max_struct_r12/src/benchmark_grouping.py
```

这不是必须重复的资格门槛；现有结果已留存。benchmark 不运行反演、不读取 toy truth 标签，只取原 public Y 与已有 operators 构造固定投影 workload。不要把2.40×单 workload 结果当成重建速度承诺。

`PhysicalAudit(objective,J,Q,metrics,costs)` 的 objective 必须是比较中的同一原 data+prior；costs 是本次 callback 的增量，不是累计计数。core 接受 scalar `(N,P)` 或 cell-major vector `(3N,P)`，没有 truth 输入。CPU budget 从 runner 初始化审计开始；setup 另测，调用者须先扣除 setup，边界内原子工作可以有小幅 overrun。

## 3. 冻结 R8 NN 四对象接口

真实输入来源：

- 当前 worktree 的 `research/research_max_struct_r8/RESULTS/training/physics/weights.pt`，原 physics checkpoint epoch86；没有 retraining。
- sibling R8 worktree：`../r8_worktree/research/research_max_struct_r8/` 下的 `configs/warm_primary_manifest.json`、其指定的 `data_public/operators/*.npz`、`data_public/warm_test/` 中四个准确 case ID 的 `.npz`，以及原 public 初始化文件。
- 原 `src/warm_features.py::FeatureOperator` 和 `src/warm_network.py::load_predictor`，十通道定义不变。

四个准确 ID 是 `r8_test_0001,0004,0022,0027`。0001为6源，其余12源；不要将它们当作 R10/R11 的同编号对象。准备期重推断与历史 public NN 初始化的相对差为0，见 [PREPROCESS_RECEIPT.json](delegated/neural/PREPROCESS_RECEIPT.json)。

新 checkout 的实际顺序：

```sh
"$R12_PY" research/research_max_struct_r12/delegated/neural/run_nn_interface.py prepare
"$R12_PY" research/research_max_struct_r12/delegated/neural/run_nn_interface.py run
# 仅在整组 ALL16_PREDICTIONS_AND4_INITIALIZATIONS_SEALED 后：
"$R12_PY" research/research_max_struct_r12/delegated/neural/score_nn_interface.py
```

`prepare` 在推断/新结果前冻结四 ID 和四臂；`run` 是一条单线程本地有限队列，依次做 uniform-direct、NN-direct、NN-DR、NN-CIE。DR为3000步/audit50/grouped v2；CIE为原 R8 mf3/6/10/12 与500/500/500/1500调度，保留 `initialize_physical_current=False`。四臂共享原最多3000迭代/120CPU秒 corrector，cap仅属于 corrector，proposal另外计费。

生成文件均位于 [delegated/neural/](delegated/neural/)：`FROZEN_PROTOCOL.json`、`PREPROCESS_RECEIPT.json`、`public/`、`runs/`、`ALL_PREDICTIONS_SEALED.json`、`OFFLINE_SCORES.json/.csv`、`OFFLINE_CPU_PREFIXES.json`、`OFFLINE_DR_AUDITS.json`、`COST_SUMMARY.json`、`materials.png`、`objective_cpu.png` 和 `RESULTS_ZH.md`。所有16个最终预测封存后才打开 sibling R8 `sealed_truth/warm_test/`；scorer 不调用新 Maxwell solve。

这些是四个暴露对象的开发诊断，不是独立24确认。材料误差采用原 R8 amplitude-normalized RMSE 与 prediction-only Otsu IoU；不要与3D的 relative L2混为一个指标。

## 4. 高对比度 true3D：两个独立几何加一个配对

新 checkout 的实际顺序：

```sh
"$R12_PY" research/research_max_struct_r12/delegated/vector_high/generate_high.py
"$R12_PY" research/research_max_struct_r12/delegated/vector_high/run_high.py run
# 仅在 ALL12_PREDICTIONS_AND_PREFIX_CHOICES_SEALED 后：
"$R12_PY" research/research_max_struct_r12/delegated/vector_high/score_high.py
```

`generate_high.py` 先冻结 seed2026111213、两独立几何及所有预算，再生成 dense dyadic DDA 观察。连续分区是带空腔随机长方体并集/材料优先分区；6³/8³材料是解析体积平均。`random_high8_crossgrid8` 与 `random_high8` 是相同几何的配对，不增加独立样本数。公开 solver 文件只有 Y、sigma、grid、bounds；材料与 clean 数据在独立 `sealed_truth/`。

`model_adapter.py` 重用原 VectorPR3D 的 physical/current/material/prepare，只给本补充的优化循环显式 `[0,8]` bounds、启用既有 exact value cache 并补充实际 action 计数；共享 R8 `[0,4]` 默认未改。逆解为216材料/648 field、60 dipoles/180电接收通道；fine forward为8³/1536 field。四臂 CIE/Dual/DR/AP 共享初始chi=.25及原quadratic prior。基线120cycles/inner12，DR/AP3000/audit50/grouped v2/320CPU cap（含setup）；共同 corrector至多1000迭代/60CPU秒。

已完成输出在 [delegated/vector_high/](delegated/vector_high/)：协议、generation receipt、12个 `RESULTS/{case}/{method}/`、完整投影 archives、`ALL12_PREDICTIONS_SEALED.json`、前缀 seal、`OFFLINE_SCORES.json/.csv`、`OFFLINE_PREFIX_SCORES.json`、`COST_LEDGER.json`、`MATERIAL_SLICES.png`、`PHI_CPU.png` 和 `RESULTS_ZH.md`。

高index cell相位约1.885（6³）/1.414（8³）；这是 coarse finite-DDA 挑战，不是 continuum-accuracy 资格。mixed 连续材料含2/8，但6³真实体积平均最大值为7.80691，已如实保留，没有按结果重画几何。cross-grid clean response差异为26.1567%，不能用1%噪声 alone 解释。

保存状态分解的实际入口如下；它不运行新反演，不做 Maxwell forward/adjoint，也不重新读材料 truth：

```sh
"$R12_PY" research/research_max_struct_r12/delegated/vector_high/decompose_saved.py
```

该脚本只求原廉价 quadratic R，使用已有 offline metrics，并对 AP/DR 已保存物理 current 做固定 S readout。产物是 `RAW_CORRECTOR_DECOMPOSITION.json`、`RAW_CORRECTOR_DECOMPOSITION_ZH.md`、`RAW_CPU_PREFIX_COMPARISON.json`、`DECOMPOSITION_RECEIPT.json`。CIE/Dual raw exact physical current未保存；精确全局残差标缺，但 source-normalized norm下界已能认证其球外。选中 corrector evaluation 可恢复，selected accepted iteration未记录；不补造。

## 5. 其他已存在 R12 入口的路由

|工作|真实入口/记录|当前复现边界|
|---|---|---|
|Pro四toy重放|`experiments/pro_reproduction/consensus_reflection_toy.py`、`REPLAY_RECEIPT.json`|原自含 toy 生成/评分代码与新的 GT-free core 分开；不宣称 bitwise 跨BLAS复现|
|同toy成熟基线|`delegated/baseline/toy_baselines.py`、`INVENTORY_AND_REPRODUCE.md`|两套 explicit Fourier/box 适配均保留，不能择优隐藏|
|R11真实坏终点|`experiments/basin_escape.py`、`BASIN_ESCAPE_REPORT_ZH.md`、`audit/BASIN_BUDGET_ACCOUNTING.md`|开发非确认；v1与实际远程成本记录单列|
|原低对比vector补充|`delegated/vector/run_vector.py`、`continue_physical.py`、`REPORT.md`|chi≤4；原120corrector与另行共同 continuation不混为同预算|
|实测/外部锚点|`delegated/external*/REPORT.md`|各自真实mask、模型/预算/数据访问限制以原report为准，不迁移synthetic operator|

这些条目只是路径导航，未在本整理轮重跑。

## 6. 独立24：已完成执行、封存和评分

当前权威合同分别为：

- [ROOT_CORE_FREEZE.json](delegated/confirmation/ROOT_CORE_FREEZE.json)：primary/core在24生成前冻结；8方法、3000 DR/audit50/grouped v2、240/144/96 CPU与30/60/120/240前缀。
- [CONFIRM24_GENERATOR_FROZEN.json](delegated/confirmation/CONFIRM24_GENERATOR_FROZEN.json)：24独立材料、4家族、high/mixed与8/18分布；public与sealed truth分开。
- [SECONDARY_SELECTION_FREEZE.json](delegated/confirmation/SECONDARY_SELECTION_FREEZE.json)：24数据已生成、确认逆解启动0、确认GT评分0时注册 noise-feasible minimum-same-R secondary。此注册晚于数据生成，不伪称生成前冻结。
- [CURRENT_EXECUTION_HANDOFF.md](delegated/confirmation/CURRENT_EXECUTION_HANDOFF.md)、运行时 manifests、修复与 transport receipts：实际队列、失败尝试、现有资源分配和费用。

确认入口是 `delegated/confirmation/generate_observations.py`、`run_baseline.py`、`finite_fifo.py`、`seal_public_cohort.py` 及 `experiments/run_confirmation_search.py`。这批队列已结束。复现时使用新的输出目录和对应runtime manifest，保留原失败与有效输出；旧README/handoff中被明确标为historical的命令不代表当前有活跃调度。

Secondary只筛选已有、合法、完整物理、cap内、CN99球内审计候选，再最小化同一 R，按Phi破同分；无候选明确 unavailable。它不改 primary或数值轨迹，新增检查/复制/写出费用计入相同预算。192确认与48开发均已封存评分；见CONFIRMATION_RESULTS_ZH.md和RESULTS/offline_scoring_v1。此文件是运行入口，不替代科学审查。


### 仅重现确认评分与图像（不启动逆解）

恢复Git内compact predictions、sealed_truth和对应公开seals后，在R12目录运行下例；指定新输出路径，避免覆盖原评分：

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python delegated/confirmation/score_confirmation.py --output /tmp/r12-reproduced-scores --figures /tmp/r12-reproduced-figures
```

这是原评分入口的已实现参数；本交付没有为生成说明而重跑科学。scorer保留10000次按对象bootstrap和主/次NA语义。额外保存的GT噪声资格、rawDR抽取与坏端点pooloracle各有独立script/receipt；它们不改变冻结算法。大状态恢复路径和来源见ARTIFACTS.md及release内容索引，完整历史轨迹不必载入即可读取compact最终材料与评分。

## GitHub release传输分片的恢复

大逻辑包上传发生HTTP408后，发布改用128MiB运输分片。它们不是新实验或新数据；按PUBLICATION_PARTS_MANIFEST.json中顺序连接便还原原验收逻辑包，不能单独unzip某个part。所有目录条目/字节索引仍保留。

只取需要的逻辑包，可用本目录`restore_release_assets.py`，必须显式指定输出目录（建议外置盘）：

```sh
python3 research/research_max_struct_r12/restore_release_assets.py \
  --asset max-struct-r12-remote-01.zip.zst \
  --output-dir /path/to/external-disk/r12-archives --decode
```

`--decode`需要已有zstd命令；得到ZIP后按ARTIFACTS的目录映射恢复。补充包本身为.zip，不需zstd。脚本只下载/连接，不运行科学。首次归档中的文稿可能早于最终修订；以Git仓库最终报告为准，原科学数组/源快照不被覆盖。
