# Classic2D实际执行命令

以下在仓库根目录执行。科学合同是 `research/research_max_struct_r9/configs/classic_extension_main.json`（两例×八臂）；Windows迁移记录为 `classic_extension_xinan.json`，CPU前缀为 `classic_extension_paid_budgets.json`。不读取或上传GT。远端已有 `D:/python/python.exe`、NumPy2.2.6、SciPy1.15.3、psutil7.2.2；本地分析Python还需Matplotlib。

原执行root为 `D:/AI/MAX_STRUCT_R9_CLASSIC_20261106`。`prepare/run`只对空的运行目标做一次；已有QUEUE_STARTED或terminal输出时程序拒绝覆盖。**新clone也可能含已发布RESULTS；独立checkout本身不能保证输出为空。** 本页先保留历史原命令，后给出不覆盖旧RESULTS的全新输出命令。连接账户和key由本机私有环境变量提供，不把tailnet地址或私钥写入公开仓库。

## 全新复现：独立本地输出，不依赖旧Mac子任务

原程序已有历史迁移guard：历史 `classic_extension_xinan` 评分要求两份Mac遗留子任务terminal。它只适用于原研究记录，不是新cohort的科学条件。下面的 `--reproduction` 入口**只要求新16项全部terminal/封存**，不读取旧Mac控制、旧Mac预测或旧CUDA费用。它仍在全部预测封存后才读GT；公开Phi的已付CPU前缀选择也先于GT。

先确认原远端任务已结束、取得四P核资源时隙。选一个从未运行过的远端root；本例的新本地输出为 `reproduction/classic`。若该路径已有数据，请换成另一个新路径，fetch会在写入前检查并拒绝覆盖任何已有文件，不删除参考结果。

```sh
export R9_PYTHON=python
export R9_REMOTE_ROOT=D:/AI/MAX_STRUCT_R9_CLASSIC_REPRODUCTION_20261106
export R9_REPRO_OUTPUT=reproduction/classic
: "${R9_SSH_HOST:?Set the authorized user@host in your local environment}"
: "${R9_SSH_KEY:?Set the local SSH identity file path}"

"$R9_PYTHON" research/research_max_struct_r9/delegated/classic_extension/remote_transport.py prepare \
  --host "$R9_SSH_HOST" --identity "$R9_SSH_KEY" --remote-root "$R9_REMOTE_ROOT" \
  --local-output-root "$R9_REPRO_OUTPUT"

"$R9_PYTHON" research/research_max_struct_r9/delegated/classic_extension/remote_transport.py run \
  --host "$R9_SSH_HOST" --identity "$R9_SSH_KEY" --remote-root "$R9_REMOTE_ROOT" \
  --local-output-root "$R9_REPRO_OUTPUT"

# In another terminal, using the same exported variables:
"$R9_PYTHON" research/research_max_struct_r9/delegated/classic_extension/remote_transport.py status \
  --host "$R9_SSH_HOST" --identity "$R9_SSH_KEY" --remote-root "$R9_REMOTE_ROOT" \
  --local-output-root "$R9_REPRO_OUTPUT"

# Fetch once all16 are terminal. It maps the remote cohort into the new local root.
"$R9_PYTHON" research/research_max_struct_r9/delegated/classic_extension/remote_transport.py fetch \
  --host "$R9_SSH_HOST" --identity "$R9_SSH_KEY" --remote-root "$R9_REMOTE_ROOT" \
  --local-output-root "$R9_REPRO_OUTPUT"

"$R9_PYTHON" research/research_max_struct_r9/experiments/classic_extension.py score \
  --reproduction --output-root "$R9_REPRO_OUTPUT"
"$R9_PYTHON" research/research_max_struct_r9/experiments/classic_extension.py figures \
  --reproduction --output-root "$R9_REPRO_OUTPUT"
"$R9_PYTHON" research/research_max_struct_r9/experiments/classic_extension.py ledger \
  --reproduction --output-root "$R9_REPRO_OUTPUT"
"$R9_PYTHON" research/research_max_struct_r9/experiments/classic_extension.py report \
  --reproduction --output-root "$R9_REPRO_OUTPUT"
```

新文件全部在 `reproduction/classic`：原始16预测、`transport/`、`offline_truth/`、OFFLINE_SCORES、`FIGURES/`、PHASE2_COST_LEDGER及CLASSIC_BENCHMARKS。原 `research/research_max_struct_r9/RESULTS/`、原FIGURES、原sealed_truth与根CLASSIC_BENCHMARKS均不写入。新账本只计本次16项和本次transport/analysis，不重复计旧Mac9、旧CUDA单点或phase1。prepare包只需要代码、config和data_public，不需要旧Mac结果或GT。

## 历史原执行命令（含迁移专用guard）

以下是原run记录，**不是在含旧结果的新clone中直接复跑的入口**；新cohort请使用上节完整命令。

```sh
# Activate the intended local analysis environment first.
export R9_PYTHON=python
export R9_REMOTE_ROOT=D:/AI/MAX_STRUCT_R9_CLASSIC_20261106
: "${R9_SSH_HOST:?Set the authorized user@host in your local environment}"
: "${R9_SSH_KEY:?Set the local SSH identity file path}"

"$R9_PYTHON" research/research_max_struct_r9/delegated/classic_extension/remote_transport.py prepare \
  --host "$R9_SSH_HOST" --identity "$R9_SSH_KEY" --remote-root "$R9_REMOTE_ROOT"

# Keep this foreground SSH transport connected until all16 workers finish.
"$R9_PYTHON" research/research_max_struct_r9/delegated/classic_extension/remote_transport.py run \
  --host "$R9_SSH_HOST" --identity "$R9_SSH_KEY" --remote-root "$R9_REMOTE_ROOT"
```

可在另一个本地终端只读查询；全部terminal后再fetch一次：

```sh
"$R9_PYTHON" research/research_max_struct_r9/delegated/classic_extension/remote_transport.py status \
  --host "$R9_SSH_HOST" --identity "$R9_SSH_KEY" --remote-root "$R9_REMOTE_ROOT"

"$R9_PYTHON" research/research_max_struct_r9/delegated/classic_extension/remote_transport.py fetch \
  --host "$R9_SSH_HOST" --identity "$R9_SSH_KEY" --remote-root "$R9_REMOTE_ROOT"
```

Windows FIFO源码是 `remote_fifo.py`。它只开4个worker，分别固定P核logical0/2/4/6，在NumPy导入前设置亲和性；BLAS各1线程、CUDA不可见，不动其他进程。输出固定在仓库相对目录 `research/research_max_struct_r9/RESULTS/classic_extension_xinan/`。外部worker CPU为Windows GetProcessTimes user+kernel，不能称Mac wait4。

原连接内实际等价的Windows入口为（只在新的空run root中执行一次）：

```powershell
Set-Location D:/AI/MAX_STRUCT_R9_CLASSIC_20261106
& D:/python/python.exe -B research/research_max_struct_r9/delegated/classic_extension/remote_fifo.py
```

本地取回后，**保留以下output-root参数**。score内部先检查Windows16项及迁移时留下的Mac两子任务全部terminal；随后先保存公开Phi的已付CPU前缀选择，最后才构造公开Austria几何真值并评分。

```sh
"$R9_PYTHON" research/research_max_struct_r9/experiments/classic_extension.py score \
  --output-root research/research_max_struct_r9/RESULTS/classic_extension_xinan

"$R9_PYTHON" research/research_max_struct_r9/experiments/classic_extension.py figures \
  --output-root research/research_max_struct_r9/RESULTS/classic_extension_xinan

"$R9_PYTHON" research/research_max_struct_r9/experiments/classic_extension.py ledger \
  --output-root research/research_max_struct_r9/RESULTS/classic_extension_xinan

"$R9_PYTHON" research/research_max_struct_r9/experiments/classic_extension.py report \
  --output-root research/research_max_struct_r9/RESULTS/classic_extension_xinan
```

figures只首次运行一次并检查实际导出的全幅材料图。ledger生成增量canonical，report生成最终材料／80→300／CPU前缀表，并将自身分析CPU下界加入canonical、设置 `finalized_for_parent_compiler=true`。最终文件为 `RESULTS/classic_extension_xinan/PHASE2_COST_LEDGER.json`；`by_platform.Mac.known_CPU_seconds` 与 `Windows.known_CPU_seconds` 已各自去重，二者和就是 `actual_new_known_CPU_seconds`。根编译器使用此canonical后，不再另加Mac9、Win16、CUDA单点或transport原回执。

原113次phase1保持不变。Mac9已完成项及迁移后7项NOT_RUN分别记录在 `MAC_PARTIAL_EXECUTION.json`；它们不和Windows16混在同一费用／效果对比表。公开Y原始生成费为R8历史复用，未重新生成或重复收费。

GPU单点回调是另一个已授权、已结束的工程测量，不是16项队列的一部分。只有GPU独占时隙和host P12被明确授权时，才可用以下原执行入口；已有benchmark terminal文件也拒绝覆盖：

```sh
"$R9_PYTHON" research/research_max_struct_r9/delegated/classic_extension/run_cuda_once.py \
  --host "$R9_SSH_HOST" --identity "$R9_SSH_KEY"
```

该单点只比较已封存m80的完整value+gradient，无优化或GT。原结果与分项费用已在 `CUDA_CALLBACK_BENCHMARK.md`，不需为复核材料指标再运行此GPU任务。
