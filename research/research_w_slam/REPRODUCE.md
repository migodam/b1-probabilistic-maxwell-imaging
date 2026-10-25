# 保存证据复核与新运行说明

本说明是命令和来源文档，没有执行安装、求解、重画或新试验。先阅读`PHYSICS_PROTOCOL.md`与`PRIOR_AND_LIKELIHOOD.md`。保存数据读回、图表重画与新增物理运行是不同动作；已有receipt不证明新机器数值等价。

## 可移植环境与依赖

`requirements.txt`记录NumPy/SciPy/Matplotlib、scikit-fem12.0.2、Triangle20250106的版本；Python3.13来自原环境记录。新平台需自行满足这些版本并验证数值结果。示例只创建本地可移植环境，不依赖某个私人解释器绝对路径：

```sh
python3.13 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

旧实际job的解释器路径仅是receipt中的历史来源，不是复现必需路径。`tools/run_job.py`给子进程设置`PYTHONPATH=ROOT/src:ROOT/vendor:ROOT.parent/research-w-continuation/.deps`并固定BLAS线程1。当前代码的`wslam.fem_backend`已内置薄适配；在独立环境中安装skfem/Triangle后，不必依赖外部`.deps`。但若相邻旧`.deps`真实存在，它可能优先于环境安装，须在新运行目录审阅该路径，不能默默混用依赖。

## 默认：已有数组和原始记录

从以下保存文件读取即可复核已有指标，不需要构造Book或FEM：

- `results/e0/qualification/RESULTS.json`、`results/e0/refinement01_executable/RESULTS.json`、`results/e0/gate.json`。
- `results/e0/four_blocks/RESULTS.json`、`BLOCKS.npz`及`contrast_*.npz`。
- `results/e1/dynamic/RESULTS.json`、原始`pose*_f*.npz`、信息与补偿端点NPZ；额外秩审核在`results/e1/rank_audit/`。
- `results/jobs/<原标签>/RECEIPT.json`、`OUTPUT.log`和共享`results/budget_ledger.json`。
- `figures/actual_e0_e1/`中的PNG/PDF/CSV、双语图注和schema缺口审计。

图表helper`tools/wslam_figures.py`只读取保存数据，不求解PDE，但绘图CPU仍应记账。它会重写固定`figures/actual_e0_e1`输出，须在新展示目录或经审阅的独立执行副本运行，不能误称只读。原独立绘图未测得精确CPU，由父线程保守计180CPU秒，保留其上界属性，不再次收费。

## 新授权工作目录中的命令形式

下列示例**只适用于另行准备、审阅并获授权的新执行根目录**。不要在已有付费结果目录直接运行。每个物理步骤经串行wrapper预约；阶段限额不因改job标签重置。新目录应使用新的campaign ledger和明确输入清单，原账本、失败和sealed证据保持原样；这不是创建旧结果备份或清空历史的步骤。

```sh
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export WSLAM_LEDGER="$PWD/results/budget_ledger.json"
PYTHON="$PWD/.venv/bin/python"

"$PYTHON" tools/run_job.py E0E1 new-qualification "$PYTHON" tools/e0_qualify.py
"$PYTHON" tools/run_job.py E0E1 new-refinement "$PYTHON" tools/e0_refine.py
"$PYTHON" tools/run_job.py E0E1 new-e0-decision "$PYTHON" tools/e0_decide.py
"$PYTHON" tools/run_job.py E0E1 new-e1 "$PYTHON" tools/e1_dynamic.py
"$PYTHON" tools/run_job.py E0E1 new-e1-rank "$PYTHON" tools/e1_rank_audit.py
"$PYTHON" tools/run_job.py A0 new-a0-preparation "$PYTHON" tools/a0_prepare.py
"$PYTHON" tools/run_job.py A0 new-a0-fitting "$PYTHON" tools/a0_run.py --root "$PWD"
"$PYTHON" tools/run_job.py A0 new-a0-scoring "$PYTHON" tools/a0_score.py --root "$PWD"
```

这些是现有入口的调用形式，**不是完整自动重现队列或预先放行**。修复plan、补充合同、输入和父线程门判断必须先就位。`e0_refine.py`会复用同一新campaign里已完成qualification的独立均值；当前可执行参考是`.17/384/grid3`，不是已被资源守卫拒绝的`.13/512/grid3`。`e0_decide.py`还包含新的外孔径PDE，不能当作纯JSON读回命令。四块端口准备也收费，且必须纳入总E0E1上限。

有限实验是条件阶段。确认其补充合同、阶段预算和前置条件后，入口形式为：

```sh
"$PYTHON" tools/run_job.py E2E3 new-finite-bank "$PYTHON" tools/finite_bank.py --output "$PWD/results/e2e3/finite_bank"
"$PYTHON" tools/run_job.py E2E3 new-revisit "$PYTHON" tools/e2_revisit.py
"$PYTHON" tools/run_job.py E0E1 new-saved-array-figures "$PYTHON" tools/wslam_figures.py
```

有限bank唯一显式恢复入口是`--resume`，仅在其checkpoint、精确完成前缀和共享账本一致且恢复获授权时使用。存在不确定已付费row时程序应停止，不能自动重试。`e2_revisit.py`仍从固定ROOT结果路径读取；把bank写入任意其他`--output`并不会自动重定向所有后续入口。

## 固定输出与可变CLI的限制

多数driver在导入后即准备固定ROOT输出，并采用`mkdir(exist_ok=False)`拒绝覆盖。改变wrapper标签不会改变这些科学输出目录。`a0_run/a0_score --root`只覆盖各自支持的根路径，不是全项目重定向；`finite_bank --output`也只控制bank自身。源码以后可能继续演进，因此历史receipt中的命令、源提交、冻结补充协议和实际输出是历史复核依据，当前CLI不是过去执行行为的证明。

不删除旧结果、旧账本、失败receipt来取得一次“干净重跑”；不额外执行哈希/校验和步骤。此目录也不重授第三方许可。FEM源来自B1研究仓库固定提交的薄适配，详见`THEORY_AND_SOURCE_AUDIT.md`与上游`LICENSE_STATUS.md`；公开可读不自动授予新复用许可，Triangle保留其自身许可条件。
