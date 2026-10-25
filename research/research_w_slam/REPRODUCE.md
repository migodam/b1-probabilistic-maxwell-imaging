# 保存证据复核与新运行说明

本说明是命令和来源文档，没有执行安装、求解、重画或新试验。先阅读`PHYSICS_PROTOCOL.md`与`PRIOR_AND_LIKELIHOOD.md`。保存数据读回、图表重画与新增物理运行是不同动作；已有receipt不证明新机器数值等价。

## 可移植环境与依赖

`requirements.txt`记录NumPy/SciPy/Matplotlib、scikit-fem12.0.2、Triangle20250106的版本；Python3.13来自原环境记录。新平台需自行满足这些版本并验证数值结果。示例只创建本地可移植环境，不依赖某个私人解释器绝对路径：

```sh
python3.13 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

旧实际job的解释器路径仅是receipt中的历史来源，不是复现必需路径。开发工作目录的`tools/run_job.py`仍给子进程设置`PYTHONPATH=ROOT/src:ROOT/vendor:ROOT.parent/research-w-continuation/.deps`并固定BLAS线程1，保留正在执行job的原来源。`tools/prepare_publication.py`只在**公开副本**中移除最后这个旧相邻依赖路径，并将该差异写入`PUBLICATION_FILES.json`；不会修改活跃工作目录的wrapper。当前`wslam.fem_backend`已内置薄适配，独立环境安装requirements即可提供skfem/Triangle；`vendor`路径不存在时不构成必需依赖。若使用开发源而非公开副本，相邻旧`.deps`可能优先于环境包，不能将混合环境误写成独立新安装。

当前wrapper使用Unix的`fcntl`、`resource`与`preexec_fn`。本轮未验证Windows入口，不能声称原生Windows可直接运行。原研究包PDF重建是可选工具，另需Pandoc、XeLaTeX与CJK字体；这些不是NumPy/FEM缓存读回的依赖。

## 公开源码与Release的目录拼接

公开stage位于现有仓库的`research/research_w_slam/`。源码、文档、图表、合成输入JSON、结果JSON/CSV、原账本和job receipt/log在Git小包中；NPZ/NPY及未进小包的原始图片另列于stage内`RELEASE_LAYOUT.json`。该清单只描述保存文件路径与字节数，不创建、上传或验证Release。

数组归档成员须采用**相对于stage根的路径**，不多加一层目录。例如归档中的`inputs/online/COARSE_BANK.npz`与`results/e2e3/finite_bank/online/FINITE_GRID.npz`应解包到`research/research_w_slam/`下面。`results/e2e3/inference/POSTERIORS.npz`、`inputs/offline/CONTINUOUS_CLEAN_MEANS.npz`及各阶段端点/复场亦保持原相对路径。只下载源码小包不能完成缓存统计复算；缺失数组不应通过启动新付费物理求解来悄悄补齐。

`--partial`发布仍是运行中快照。Release清单需在相关job终止并完成原始receipt后重新生成；其存在不证明数组集完整。历史JSON中保留的绝对输入/预测/解释器路径是来源信息，不是可移植的运行路径。已有`a0_analysis`会对缺失的旧`results/...`路径在当前root下重定位；其他读回应优先按本说明的stage相对路径，不能依赖旧机器的挂载点。详细源码/路径审计见`research/delegated/publication/PORTABILITY_AUDIT.md`。

## 默认：已有数组和原始记录

从以下保存文件读取即可复核已有指标，不需要构造Book或FEM：

- `results/e0/qualification/RESULTS.json`、`results/e0/refinement01_executable/RESULTS.json`、`results/e0/gate.json`。
- `results/e0/four_blocks/RESULTS.json`、`BLOCKS.npz`及`contrast_*.npz`。
- `results/e1/dynamic/RESULTS.json`、原始`pose*_f*.npz`、信息与补偿端点NPZ；额外秩审核在`results/e1/rank_audit/`。
- `results/jobs/<原标签>/RECEIPT.json`、`OUTPUT.log`和共享`results/budget_ledger.json`。
- `figures/actual_e0_e1/`中的PNG/PDF/CSV、双语图注和schema缺口审计。

图表helper`tools/wslam_figures.py`只读取保存数据，不求解PDE，但绘图CPU仍应记账。它会重写固定`figures/actual_e0_e1`输出，须在新展示目录或独立执行副本运行，不能误称只读。原独立绘图未测得精确CPU，由父线程保守计180CPU秒，保留其上界属性，不再次收费。

## 独立复现工作目录中的命令形式

下列示例适用于保留原始证据的独立复现目录。不要在已有付费结果目录直接启动全量运行。每个物理步骤经串行wrapper预约；阶段限额不因改job标签重置。新目录应使用新的campaign ledger和明确输入清单，原账本、失败和sealed证据保持原样；这不是创建旧结果备份或清空历史的步骤。

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

这些是现有入口的调用形式。冻结补充合同、输入、阶段预算和各入口的真实前置条件必须先就位；完整运行消耗相应物理预算。`e0_refine.py`会复用同一新campaign里已完成qualification的独立均值；当前可执行参考是`.17/384/grid3`，不是已被资源守卫拒绝的`.13/512/grid3`。`e0_decide.py`还包含新的外孔径PDE，不能当作纯JSON读回命令。四块端口准备也收费，且必须纳入总E0E1上限。

有限实验是条件阶段。以下顺序补全了终态费用finalizer、独立bank验证与E2–E5入口；每个科学输出目录须不存在，cached-statistics另述：

```sh
"$PYTHON" tools/run_job.py E2E3 finite-bank-001 "$PYTHON" tools/finite_bank.py --output "$PWD/results/e2e3/finite_bank"
# 等待上一个wrapper完全终止；不能只看到bank内COMPLETE就并发运行finalizer。
"$PYTHON" tools/run_job.py E2E3 bank-cost-final-001 "$PYTHON" tools/finalize_bank_cost.py --root "$PWD"
"$PYTHON" tools/run_job.py E2E3 bank-validation-001 "$PYTHON" tools/finite_bank_validate.py
"$PYTHON" tools/run_job.py E2E3 revisit-001 "$PYTHON" tools/e2_revisit.py
"$PYTHON" tools/run_job.py E2E3 auxiliary-controls-001 "$PYTHON" tools/e3_negative_controls.py --bank "$PWD/results/e2e3/finite_bank"
"$PYTHON" tools/run_job.py E4 approximation-001 "$PYTHON" tools/e4_compare.py --root "$PWD"
"$PYTHON" tools/run_job.py E5 sensing-001 "$PYTHON" tools/e5_sensing.py --bank "$PWD/results/e2e3/finite_bank" --reference-gate "$PWD/results/e2e3/inference/gate.json"
```

有限bank唯一显式恢复入口是`--resume`，仅在其checkpoint、精确完成前缀和共享账本一致时使用。存在不确定已付费row时程序应停止，不能自动重试。`finalize_bank_cost.py`的wrapper标签须避开`finite-bank-`前缀；它需要所有该前缀job终止，包含失败001的真实CPU，写`online/JOB_COST_FINAL.json`，不修改原在途snapshot。`finite_bank_validate.py`另付12个固定索引的独立4TX验证（48RHS/12FE）；它不是读取旧PASS。E2要求该独立检查PASS，E5另要求有限参考gate中的`reliable_actual_reference=true`；失败时保留依赖阶段NOT_RUN，不调整科学条件。E4/E5使用finalizer终态费用，不能把bank内在途CPU写成完整设置成本。`e2_revisit.py`仍从固定ROOT结果路径读取；把bank写入任意其他`--output`并不会自动重定向所有后续入口。


## 缓存统计复算：不新增FEM

在已经拼接全部数组的stage根执行以下命令。它们仅读取保存数组并写新统计/图表目录；CPU由wrapper记录，RHS/LU为零。不要把这些命令与上面的付费bank验证混淆。`a0_analysis`保留原始有shape错误的分析，复算写独立目录；统计复算不重新选科学门：

```sh
PYTHON="$PWD/.venv/bin/python"
export PYTHONDONTWRITEBYTECODE=1
export WSLAM_LEDGER="$PWD/results/budget_ledger.json"
"$PYTHON" tools/run_job.py A0 cached-a0-analysis-001 "$PYTHON" tools/a0_analysis.py --root "$PWD" --output "$PWD/results/a0/analysis_recomputed" --clean-mean-convention includes-gamma
"$PYTHON" tools/run_job.py E0E1 cached-summary-001 "$PYTHON" tools/summarize_pilot.py --root "$PWD" --output-json "$PWD/results/PILOT_SUMMARY_RECOMPUTED.json" --tables "$PWD/results/pilot_summary_tables_recomputed"
"$PYTHON" tools/run_job.py E0E1 cached-figures-001 "$PYTHON" tools/pilot_figures.py --root "$PWD" --output "$PWD/figures/pilot_recomputed"
"$PYTHON" tools/run_job.py E0E1 evidence-export-001 "$PYTHON" tools/export_evidence.py --root "$PWD"
```

上述label不会重置共享上限；若原campaign剩余CPU不足，应在独立复现campaign中记录本次费用，保留原账本与receipt，不能清空旧账本。`analysis_recomputed`、`pilot_summary_tables_recomputed`等输出目录须不存在。汇总脚本按父世界先平均两个噪声repeat；24噪声条件不是24独立父样本。直接读取JSON/CSV或数组本身无需新增物理调用。

## 数组Release归档

`tools/package_release.py`以实际保存NPZ/NPY为准，结合`RELEASE_LAYOUT.json`（若有）和有界科学目录枚举，保留失败/raw数组与完整场快照；不载入数组、不产生哈希、不删除源文件。输出`release_assets/`由根`.gitignore`排除，归档成员均是stage相对路径。完整模式先检查bank/inference终态，不要求科学PASS；若前置科学失败使依赖阶段明确NOT_RUN/blocked，可提供真实终态状态清单：

```sh
"$PYTHON" tools/run_job.py E0E1 release-package-001 "$PYTHON" tools/package_release.py --root "$PWD" --mode final
# 仅用于已有显式失败/NOT_RUN终态声明的发布：
"$PYTHON" tools/run_job.py E0E1 release-failure-package-001 "$PYTHON" tools/package_release.py --root "$PWD" --mode final --allow-not-run --status-manifest "$PWD/GATE_SUMMARY.json"
```

归档清单`release_assets/ASSET_MANIFEST.json`记录每个成员路径、源字节数、归档成员大小、压缩大小、总archive bytes及缺失/陈旧layout项，不声称科学完整。失败发布仍只打包真实文件，不补造缺失数组。运行中仅用`--mode partial`并保留其snapshot标签。先解包源码stage，再在仓库根将每个数组包解到同一stage：

```sh
unzip research-w-slam-E0E1.zip -d research/research_w_slam
unzip research-w-slam-A0.zip -d research/research_w_slam
unzip research-w-slam-controlled_E2E3.zip -d research/research_w_slam
unzip research-w-slam-E4E5.zip -d research/research_w_slam
unzip research-w-slam-inputs.zip -d research/research_w_slam
# 如清单包含辅助证据包，同样解包到上述stage根。
```

只执行清单实际存在的包；没有某阶段归档不等于该阶段已通过。归档创建与网络上传是分开的动作。

## 固定输出与可变CLI的限制

多数driver在导入后即准备固定ROOT输出，并采用`mkdir(exist_ok=False)`拒绝覆盖。改变wrapper标签不会改变这些科学输出目录。`a0_run/a0_score --root`只覆盖各自支持的根路径，不是全项目重定向；`finite_bank --output`也只控制bank自身。源码以后可能继续演进，因此历史receipt中的命令、源提交、冻结补充协议和实际输出是历史复核依据，当前CLI不是过去执行行为的证明。

不删除旧结果、旧账本、失败receipt来取得一次“干净重跑”；不额外执行哈希/校验和步骤。此目录也不重授第三方许可。FEM源来自B1研究仓库固定提交的薄适配，详见`THEORY_AND_SOURCE_AUDIT.md`与上游`LICENSE_STATUS.md`；公开可读不自动授予新复用许可，Triangle保留其自身许可条件。

最终图表补充：默认pilot_figures之后运行 `tools/refine_plot_labels.py`，只修正两张图的坐标、统一概率色标和方法编号，不修改已保存数组与CSV。E2图采用原两副本posterior平均，是重复实验描述，不是平均数据的posterior。最终独立端点审核为 `tools/posterior_endpoint_audit.py`，须在E2完整权重保存之后执行；最多96RHS/24FE，共用E2E3原预算，既不重加权posterior也不新增概率门。

本轮公开数组入口：https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-slam-20261025 。六个ZIP均按stage相对路径解包；无需执行新Maxwell即可读取已有物理均值、Jacobian、联合权重及复场快照。归档成员与实际字节数见RELEASE_ASSETS.json，没有哈希步骤。
