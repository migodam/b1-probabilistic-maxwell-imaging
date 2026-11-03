# R5 复现入口与冻结证据保护

本文件来自实际源码接口检查，本次没有执行以下物理命令，也未新增tests或hash检查。只读成本重算结果见 `COST_LEDGER.md`、`RESULTS/COST_LEDGER.json`，原始成本仍保留在各SUMMARY/RECEIPT。

## 环境和源代码

工作目录是本仓库根目录，包含 `research/research_max_struct_r5/` 以及引用的R2/R4、FBG-R3、Z-ALG-R4源与公共数据。Python依赖从真实imports确认：NumPy、SciPy、threadpoolctl；绘图另需matplotlib。CPU complex128/float64；各数值runner设BLAS/OpenMP单线程。本轮不用PyTorch/CUDA、MPS、远端或服务。

历史回执记录的解释器是 `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python`。这是历史执行路径，不保证异机存在；系统python在库存阶段缺NumPy，故不可假定系统python已具备依赖。本轮没有安装依赖或验证新环境。可以在已有兼容环境或独立虚拟环境中提供这些依赖；没有据空缺版本记录补造精确版本锁。

精确执行来源按批次归档：`audit/trajectory_replay_source_snapshot/`、`audit/factorial_source_snapshot/`、`audit/factorial_extension_source_snapshot/`、`audit/handoff_source_snapshot/`。源码与输入路径是复现依据，不需要hash验证。

## 输出路径是固定的，先准备隔离副本

replay仅有`--run-id`；handoff仅有`--handoff-id`；诊断只有互斥`--run-id`/`--all-ready`/`--legacy-stage-states`。它们没有`--output`参数，ROOT由脚本自身位置决定。不能把本例伪装成可传一个新output参数。

在另一全新目录中准备所需仓库结构、冻结源码/config与历史public/reference输入；保留R2/R4历史RESULTS和reference文件，但不要把原R5已生成的RESULTS、scoring seals/receipts当成可覆写输出带入新运行。新的R5 `RESULTS`由runner建立；若需源代码快照，应在隔离副本中恢复相应批次源码，不改原始冻结快照。以下命令只能在这种**全新隔离副本**中运行，绝不能删除原交付目录来腾位置。

`run_trajectory_replay.py`、`run_factorial.py`、`run_stage_handoff.py`使用`mkdir(exist_ok=False)`拒绝现有run目录。`analyze_trajectories.py`会跳过已完成回执并复用已有checkpoint JSON，故在原冻结目录运行不会构成新的独立重放；必须隔离。离线评分也有固定路径和新建文件保护，且会读GT，须待全部新预测完成封存后才执行。

## 8个辅助轨迹replay

从`configs/trajectory_replay_frozen.json:runs`读取全部8个匿名run_id；每个执行一次，至多2个进程、每进程单线程。例如：

```bash
python research/research_max_struct_r5/experiments/run_trajectory_replay.py --run-id obs_b20d5a55677d4eca__global
python research/research_max_struct_r5/experiments/run_trajectory_replay.py --run-id obs_b20d5a55677d4eca__onehop
```

其余case为`obs_42cb233e4f224bb3`、`obs_63abafa750fc483c`、`obs_100bd65498644322`，各取`__global`、`__onehop`。不改冻结日程、初始值、box/prior、checkpoint位置，不增加随机seed或corrector。保存的21 checkpoint/3 stage snapshots和auxiliary candidate材料是后续诊断输入。历史replay成本只有auxiliary及既有initial/stage物理验证，不包含后续诊断。

## 174个固定状态诊断

在隔离副本的8 replay全部完成后：

```bash
python research/research_max_struct_r5/experiments/analyze_trajectories.py --all-ready
python research/research_max_struct_r5/experiments/analyze_trajectories.py --legacy-stage-states
```

第一条处理168 checkpoint；第二条使用历史已付费状态做6个新增诊断。也可用`--run-id <完整run_id>`处理单条轨迹，不能与另外两选项同用。诊断会付新的forward/adjoint/Jacobian费用，不是零物理评分。原首次失败保留在冻结证据中；复现不应人为再次制造失败或把其费用从历史账本抹除。

`ZERO_PHYSICS_CORRECTIONS_RECEIPT.json`描述历史后处理，当前并无单独通用CLI可承诺复现该整理步骤；冻结修订后的原始JSON/NPZ是其证据。勿把诊断runner重新执行称为零物理后处理。

## 2 factorial、3 extension、2 handoff

```bash
python research/research_max_struct_r5/experiments/run_factorial.py --arm current_only
python research/research_max_struct_r5/experiments/run_factorial.py --arm material_only
python research/research_max_struct_r5/experiments/run_factorial.py --arm material_only --config factorial_extension_frozen.json --case obs_42cb233e4f224bb3
python research/research_max_struct_r5/experiments/run_factorial.py --arm material_only --config factorial_extension_frozen.json --case obs_63abafa750fc483c
python research/research_max_struct_r5/experiments/run_factorial.py --arm material_only --config factorial_extension_frozen.json --case obs_100bd65498644322
python research/research_max_struct_r5/experiments/run_stage_handoff.py --handoff-id global_stage0_to_material_only
python research/research_max_struct_r5/experiments/run_stage_handoff.py --handoff-id material_only_stage0_to_global
```

factorial的`--config`接受绝对路径，或解析为R5 configs目录内相对路径；默认`factorial_frozen.json`。`--case`只能选所给合同内case；extension只允许material_only。factorial虽然支持config内output_subdir，数值合同仍须冻结；不要假设其他runner也支持该配置项。handoff读取固定`handoff_frozen.json`并需要上面新生成的global/material_only stage_0与其native coefficients，先运行相应源臂。corrector迭代上限1000属于冻结停止规则，保留其stop，不增加预算填满或换初值。

## 封存后的零物理评分

```bash
python research/research_max_struct_r5/experiments/offline_score_stage_i.py
python research/research_max_struct_r5/experiments/offline_score_factorial_extension.py
```

无CLI参数；ROOT相对脚本，写固定RESULTS/audit文件，故同样只在新副本运行。它们会先建立/检查seal，再读已披露开发对象的m/mask，仅用于离线评分；不能将任何GT方向或指标回送在线选择。handoff专用评分入口现已提供，见下方补充。评分回执新增LU/RHS=0，与诊断runner有本质成本区别。

## 回放可比性的限制

冻结数组与回执是实际历史证据。新机器BLAS、浮点与计时不同可改变PR/优化器轨迹；固定迭代/阶段记录与实际stop应原样保存，不承诺逐位或同耗时重现。source snapshot、公开文件路径、合同和全部失败记录一起保留；不以再次求解覆盖原始交付。网页/绘图和文献检索若另有产物，仅消费冻结结果与既有资料，独立记录耗时，不混成新增物理。

## 最终补充入口：已知材料vector pair、保存数组分析与图

这些真实入口均无CLI参数，ROOT相对脚本，固定输出；同样只在新隔离副本中运行，勿覆盖原冻结结果。

```bash
python research/research_max_struct_r5/experiments/vector_pair_mechanism.py
python research/research_max_struct_r5/experiments/analyze_block_modes.py
python research/research_max_struct_r5/experiments/offline_score_handoff.py
python research/research_max_struct_r5/experiments/summarize_mechanisms.py
python research/research_max_struct_r5/experiments/make_figures.py
```

`vector_pair_mechanism.py`读取固定`configs/vector_pair_frozen.json`，用R1 TensorDDA，写`RESULTS/vector_pair/`，`mkdir(exist_ok=False)`拒绝已有目录。它实际新增1次6×6 LU、2复RHS，外部near API同时计算E/H但不另解线性系统。已知两cell各向同性材料、固定照明与两beta的机制例，不含反演、GT评分、任意三维成像或新的独立Mie验证。

`analyze_block_modes.py`读取168个保存checkpoint NPZ与冻结public G矩阵，做4×4 block矩阵分析，不调用forward或优化器；输出`RESULTS/block_modes/`及`audit/BLOCK_MODES_RECEIPT.json`。源码允许已有目录并会写结果，故必须隔离，不能把存在目录的保护误写为全脚本通用特性。

`offline_score_handoff.py`在两handoff预测完成后封存并独立评分，读取开发对象m/mask，写`RESULTS/OFFLINE_HANDOFF_SCORES.*`和相应seal/receipt；新增物理0。`summarize_mechanisms.py`读取已经封存的离线评分和诊断JSON，写`RESULTS/MECHANISM_SUMMARY.json`及`audit/MECHANISM_SUMMARY_RECEIPT.json`，不重新求解。二者数据依赖决定它们位于数值运行之后。

`make_figures.py`读取保存诊断、预测、离线指标及已披露truth进行离线展示，写`FIGURES/`和`audit/FIGURE_RECEIPT.json`，不增加物理，需matplotlib。实际figure receipt没有wall/cpu计时字段，成本保持null，不能以生成时间戳推导耗时。

历史数字1–9参考表见`RESULTS/HISTORICAL_DIGITS_REFERENCE.md/json`：仅摘录FBG-R3成熟E3已发布评分18行（single/dual CIE），没有新runner或重评分。R5各方法九数字均NOT_RUN，历史45秒是nominal预算，不能当误差/求解难度上限证书；相关历史LU/RHS不计R5新成本。

## 补充图与发布导出

`python research/research_max_struct_r5/experiments/plot_vector_pair.py` 只读已保存的vector-pair结果，绘制图07；无新增物理。`python research/research_max_struct_r5/experiments/export_report.py --archive-dir /ABSOLUTE/NEW/RELEASE_DIRECTORY` 仅渲染现有中文报告及打包证据，不运行科学实验。
