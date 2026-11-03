# R5 成本账本

只汇总实际回执，未运行正演或新测试；未作科学解释。机器可读源为 `RESULTS/COST_LEDGER.json`。历史已付费算子和历史九数字参考成本不重新收费。

| 新工作分组 | 完成数量 | sum_task_wall/s | sum_task_CPU/s | 新LU | 新complex RHS |
|---|---:|---:|---:|---:|---:|
| trajectory_replay | 8 runs | 115.943460 | 113.480139 | 32 | 384 |
| trajectory_diagnostics | 174 states + 1 failed attempt | 9.685399 | 9.346890 | 439 | 9672 |
| factorial | 2 runs | 43.562694 | 42.599467 | 2062 | 49416 |
| factorial_extension | 3 runs | 90.157787 | 88.085800 | 3121 | 74784 |
| handoff | 2 runs | 31.046731 | 30.262018 | 2062 | 49428 |
| vector_pair | 1 runs | 0.273934 | 0.137581 | 1 | 2 |
| 已记录合计 | — | 290.670005 | 283.911895 | 7717 | 183686 |

`vector_pair`为已知材料两cell机制例，1LU/2RHS计入物理成本，**不是三维成像**。174诊断包含168 checkpoint与6历史固定态新增诊断；首次失败1LU/24RHS已计一次。

sum_task_wall不是并行历时；batch elapsed未知null。auxiliary、observer、Green、corrector及callback计时嵌套在run外层，不能再加一次；physical_model_cost已含corrector。全部本地CPU，GPU任务/worker/seconds=0，remote任务/seconds=0；美元费用/LLM tokens未知null。

独立非物理工作：

| 工作 | wall/s | CPU/s | 新LU/RHS |
|---|---:|---:|---:|
| OFFLINE_SCORING_RECEIPT | 0.27090995805338025 | 0.18098599999999998 | 0/0 |
| OFFLINE_FACTORIAL_EXTENSION_SCORING_RECEIPT | 0.021885375026613474 | 0.014455999999999997 | 0/0 |
| ZERO_PHYSICS_CORRECTIONS_RECEIPT | 3.3207458341494203 | 3.033463 | 0/0 |
| literature | UNKNOWN/null | UNKNOWN/null | 0/0 |
| cost_accounting_and_documentation | UNKNOWN/null | UNKNOWN/null | 0/0 |
| BLOCK_MODES_RECEIPT | 0.5344545003026724 | 0.302718 | 0/0 |
| MECHANISM_SUMMARY_RECEIPT | 37.69185700034723 | 0.050568 | 0/0 |
| FIGURE_RECEIPT | UNKNOWN/null | UNKNOWN/null | 0/0 |
| OFFLINE_HANDOFF_SCORING_RECEIPT | 0.01060175010934472 | 0.010134000000000004 | 0/0 |

已记录非物理timers小计wall=41.850454s、CPU=3.592325s，非完整非物理成本。FIGURE_RECEIPT已存在，但无wall_seconds/cpu_seconds键，因此不能补造计时。文献/整理费用也未知。具体口径见 `audit/COST_ACCOUNTING_NOTES.md`。
