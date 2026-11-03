# R5 成本口径与来源

本次只读汇总，不运行正演/伴随/Jacobian，也没有运行新tests或hash验证。输出仅为文档和JSON，不改冻结回执。

## 唯一加总层

- 8 replay：`RESULTS/TRAJECTORY_REPLAY_SUMMARY.json:totals_sum_of_run_costs`，与8个run `physical_model_cost`逐一求和一致。
- 174诊断：唯一计入 `RESULTS/trajectory_diagnostics/COST_TOTALS.json:total_new_diagnostic_physical_model_cost`、`total_timed_wall_seconds`、`total_timed_cpu_seconds`。174=168新checkpoint+6历史固定态；历史态重新诊断是本轮新增开销。
- 2 factorial：`RESULTS/FACTORIAL_SUMMARY.json:runs[].physical_model_cost` 与外层wall/cpu。
- 3 extension：`RESULTS/FACTORIAL_EXTENSION_SUMMARY.json:runs[].physical_model_cost` 与外层wall/cpu。
- 2 handoff：`RESULTS/HANDOFF_SUMMARY.json:runs[].physical_model_cost` 与外层wall/cpu。

后3组与对应目录逐run RECEIPT中的LU/actual_complex_RHS核对一致。没有把汇总和明细再次相加。`OFFLINE_FACTORIAL_EARLY_SUMMARY.json`是离线指标，不是新数值运行，不重复计费。

## 嵌套和计量单位

`corrector.cost`已经在该run `physical_model_cost`内，auxiliary initial/stage validations也在其中。callback/observer/Green/setup写盘等已在外层计时内。replay内层auxiliary timer与callback wall不可加到process wall；corrector wall不能加到run wall。`actual_complex_RHS`已经包括source/adjoint/receiver/tangent等完整物理RHS，不能再把分项加一次。`full_factorization_dimension`是维数累加，不是LU次数。

small4真实Gram求解、block RHS、scalar divisions、Green apply columns、material inner evaluations是不同单位；保留在原始SUMMARY/RECEIPT，不能冒充完整Maxwell LU/RHS。可直接读`trajectory_diagnostics/SUMMARY.json:auxiliary_operations`与`legacy_fixed_states.auxiliary_operations`；其中material_inner_evaluations包含源码已有3次direct exact-state-gradient请求，另有明确计数，不把其当额外物理solve。

历史G/E/S导入无需重建，其历史建立费用不转记本轮。实际load/setup耗时仍在外层任务计时，`cached_operator_imports`只作为导入计数。历史参考端点、历史auxiliary/corrector费用和handoff被复用的host Phi不重新收费。新回执确实产生的诊断full forward/adjoint仍照实收费。

## 首次失败

`trajectory_diagnostics/obs_b20d5a55677d4eca__global/RECEIPT.json`保留`prior_attempts`、`failed_attempts`与IndexError堆栈。异常在`frozen_ls`对一维a作二维索引时发生。该次已付1LU/24complex RHS、wall .028636291157454252s、CPU .02110899999999999s；COST_TOTALS已经包含一次。异常发生前辅助操作计数没有返回，无法恢复，保持UNKNOWN/null，不用零填补，也不重复计数嵌套的同一次失败。

## 时间、平台和未知费用

最终总账sum_task_wall=290.6700045838952s，sum_task_CPU=283.911895s；这是记录内的任务计时总和，不是并行批次历时。replay允许至多2进程、每进程单线程，任务有重叠。诊断timers在imports后开始；其他runner scope也保留各自最终receipt写入等排除项。未记录的启动、导入、最终写入、整理及人类等待时间不可补造。

独立评分与零物理后处理按自己的receipt另列；FIGURE_RECEIPT已存在但没有wall_seconds/cpu_seconds键，因此图耗时仍null；文献/文档整理无独立计时仍列null。其新增物理为0，不等于耗时0。美元和LLM token无可靠计费源，均null。全部本地CPU、GPU任务/worker/seconds=0、remote任务/seconds=0；不从“未使用GPU”推断CPU计算免费。

本次核对仅为解析JSON和求和，不是新物理或算法测试；不从成本数字推导科学结论。


## 最终追加，不回算历史

`RESULTS/vector_pair/RECEIPT.json`新增1LU/2complex RHS、wall .2739342921413481s、CPU .137581s，已独立追加到总账；operator_setup、solve_receipt与external_near为嵌套细项不再加。它为已知材料固定两cell机制例，imaging_status=NOT_AN_IMAGING_EXPERIMENT。

最终物理合计7717LU/183686complex RHS。独立非物理追加：BLOCK_MODES_RECEIPT wall .5344545003026724/CPU .302718s；MECHANISM_SUMMARY_RECEIPT wall37.69185700034723/CPU .050568s；OFFLINE_HANDOFF_SCORING_RECEIPT wall .01060175010934472/CPU .010134s。计入各自真实计时，不从wall/CPU差异推断额外计算。FIGURE_RECEIPT仅有generated_utc和零新物理计数，无耗时键，列null。

历史九数字参考来自已发布`research/research_fbg_r3_experiments/results/e3_all_scores.json:rows`，18行完整摘录single/dual CIE的online selected指标、停止与实际费用；无重评分、无GT数组读取、无新物理。历史费用只在参考附录展示，不计本轮成本。
