# R4 开发评分与成本独立审查

范围：仅历史 R2 公开开发案例 `obs_b20d5a55677d4eca`。分类为 **POST_SEAL_PUBLIC_DEVELOPMENT**；GT 已公开，不是盲测。科学解释及最终结论由父线程负责。

## 封存与评分

先保存原九端点结果、评分代码和既有 seal 至 `results/scoring_history/previous_nine/`，再新增四个已封存预测。`audit/DEVELOPMENT_PREDICTIONS_SEALED_V3.json` 登记全部 13 个预测的路径、大小、mtime 和 receipt 状态，写入后才读取公开开发 GT。未读取任何 confirmation GT。

原九行评分逐字段保持一致：`True`。新 seal 时间 `2026-11-03T04:57:20.755924+00:00`，本轮首次开发 GT 读取时间 `2026-11-03T04:57:20.756364+00:00`。

沿用 R2 冻结评分函数 `score_blind.py::{otsu,iou,prior}`：RMSE/12 ≤ 0.15；预测自身 256-bin Otsu IoU ≥ 0.7；原观测 99% CN 白化残差门；公开材料边界及原数值容差。没有新增门，也没有重新选择预测。

| 端点 | RMSE/12 | Otsu IoU | 相对残差 | 白化残差 | 原联合门 |
|---|---:|---:|---:|---:|---|
| A_block_absorption | 0.042653872 | 0.928571429 | 0.000967014434 | 23.208405 | PASS |
| global_common_corrector | 0.389182533 | 0.422680412 | 0.107820545567 | 2587.699647 | FAIL |
| self_absorption | 0.389417644 | 0.422680412 | 0.107820651645 | 2587.702193 | FAIL |
| block_absorption_cached | 0.041691092 | 0.928571429 | 0.000969536827 | 23.268943 | PASS |
| native_extended | 0.557673743 | 0.255813953 | 0.097706558936 | 2344.963353 | FAIL |
| onehop_absorption | 0.041672481 | 0.939759036 | 0.000967887255 | 23.229353 | PASS |

原 99% 白化半径为 25.165422196633。从保存的预测复核物理目标与 receipt 的最大差为 5.55e-17，相对残差差为 0；无需 forward。

## 同 prior、物理目标及停止规则

六个比较端点共用冻结标量 TM 模型、同一公开 Y、host=0.1、材料范围 [0.02,20]、历史 smooth=TV=1e-4 prior。原物理目标为 Phi=0.5||Ypred−Y||²/||Y||²+R(m)。R(m)=1e-4 mean(d²)+1e-4 mean(sqrt(d²+0.01²)−0.01)，d 为相邻材料差/19.5。runner 将同一个 `model.reg.evaluate(m)[:2]` 传入辅助读出；未用 GT 选端点。辅助变量中的状态残差目标随方法不同，不能因此称辅助目标相同。

证据：`research_fbg_r3_experiments/src/physics.py::{MaxwellModel.physical_score,data_value_gradient}`；`research_z_alg_r4_diversity/src/physics.py::Regularizer`；各端点保存的配置、runner 源码、AUXILIARY 和 CORRECTOR receipt。

A 直解、self、cached 和 onehop 使用 mf3/6/10 各 500 次全局 PR 更新；源指标、DPCS 和 beta 日程保持原设置。物料读出最多 60 次、ftol=1e-12、gtol=1e-8。self 仅保留对角吸收；cached 以逆缓存/矩阵乘代替小块重复求解，浮点轨迹不同；onehop 保留一次跨格残差交换。native_extended 保留原 sourcewise CIE，额外执行 mf12 的 6500 次更新，合计 8000 次，其辅助计算预算不同。

global_common_corrector 复用历史 mf3/6/10 全局 CIE 辅助端点；原历史终点也作为 incumbent 保留。各方法都用原 Phi 排序。所有六个端点调用同一 `cooperative_paths.refine`：L-BFGS-B，目标乘 1e6，maxiter=1000、maxfun=4000、ftol=1e-12、缩放后 gtol=0.01（原目标 1e-8）、maxls=30。全部实际在 1000 次迭代上限停止。此记录保留预算终止状态，不认证收敛或局部最优。

| 端点 | 辅助更新 | 内层物料评估 | corrector评估 | corrector停止 |
|---|---:|---:|---:|---|
| block_absorption | 1500 | 61440 | 1026 | MAXITER 1000 |
| self_absorption | 1500 | 70051 | 1048 | MAXITER 1000 |
| block_absorption_cached | 1500 | 58599 | 1019 | MAXITER 1000 |
| native_extended | 8000 | 40000 | 1011 | MAXITER 1000 |
| onehop_absorption | 1500 | 57669 | 1025 | MAXITER 1000 |
| global_common_corrector | 1500（历史复用） | 7500 | 1011 | MAXITER 1000 |

辅助各阶段均记录 MAXITER；内层读出的具体停止计数保存在 `DEVELOPMENT_COST_AUDIT.json`，不将达到外层日程称为收敛。

## 实际成本

| 端点 | 新物理 LU | 新物理 complex RHS | 进程 wall s | 进程 CPU s | 辅助 wall s |
|---|---:|---:|---:|---:|---:|
| block_absorption | 1033 | 24708 | 59.423727 | 56.562409 | 43.246323 |
| self_absorption | 1055 | 25236 | 37.285989 | 34.663160 | 26.061381 |
| block_absorption_cached | 1026 | 24540 | 46.374740 | 43.600083 | 35.428253 |
| native_extended | 1015 | 24324 | 73.594180 | 68.743973 | 59.806890 |
| onehop_absorption | 1032 | 24684 | 42.165716 | 39.800771 | 31.142698 |
| global_common_corrector | 1013 | 24288 | 10.607888 | 10.154329 | 5.181024（历史） |

global 的新 wall 仅含本轮加载、排序、公共 corrector 与保存；另有历史辅助 5.181085 s、3 LU/36 RHS，不重复计入新成本。它还保留旧全球终点；完整历史全球运行本身为 689 LU/16452 RHS、wall 9.220951 s，见 global receipt，不能将其全部算作新的辅助成本。其余完整新运行包含 imports/公共数据读取、辅助执行、物理排序、corrector 与产物写入；不含未重跑的历史算子组装，最终 receipt 自身写入除外。

| 端点 | 局部工作（与物理 LU/RHS 分开） |
|---|---|
| block_absorption | 4x4 分解 18342720；求解 RHS 261872640 |
| self_absorption | 标量复除法 1166529024 |
| block_absorption_cached | 4x4 分解 8366544；求解 RHS 33466176；逆矩阵应用 RHS 252054144 |
| native_extended | source-mode-iterations 48408000；该代理量不是 FLOPs 或等时指标 |
| onehop_absorption | 标量复除法 1064386944；onehop 块矩阵应用 RHS 248840064 |
| global_common_corrector | 复用辅助端点；当前不重跑局部操作 |

direct 与 cached 的物理设置、辅助日程和 corrector 上限相同，但求解顺序改变了轨迹、内层次数与 corrector 评估数。这里保留每次实际 wall 及工作量，未做重复计时，也未归因全部 wall 差。不同方法的局部工作、物料内层次数和 native 额外容量不同，因此不是等算力比较。

六个运行 receipt 均记录 failed_attempts=[]。既有启动失败位于原审计记录；本次评分未启动求解器、未重跑检查，评分脚本本身新增 0 LU / 0 RHS，无 hash 检查，无确认集 GT 读取。费用/货币成本未提供，保留为未知。

本审查输出是评分、协议一致性与成本事实，不替代父线程对新颖性、机制、成像收益或外推范围的判断。
