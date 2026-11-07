# Final online endpoint completion-cost plotting

交付状态：**PREPARED / SCIENTIFIC_PLOT_NOT_RUN**。本任务只准备代码和说明，没有打开 TEST score、truth 或预测数组，没有运行绘图、物理计算、新测试或数值核验。

父线程完成最终 seal 与 TEST scoring 后，在 study 目录执行：

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" experiments/plot_completion_costs.py
```

默认输入为 `results/TEST_SCORES.json` 与 `audit/FINAL_RUNS.json`。也可使用同名 JSON provenance 仍在的 scorer CSV：`--scores results/TEST_SCORES.csv`。无需读取任何 NPZ、truth 或候选材料。缺少最终分数时，只输出 `COMPLETION_COST_STATUS.json: NOT_RUN_NO_FINAL_TEST_SCORES`，不画假数据。

## 曲线的精确定义

clean 与 noisy 两行分开，每个方法每个噪声层的分母固定为108个 TEST 母体；不合并噪声，不合并模型 seed。三列横轴 B 分别为 scorer 已记录的实际 `online_total_wall_seconds`、`online_total_LU`、`online_total_RHS`。曲线值为：

> 108个母体中，冻结的最终 online top1 已满足 joint_success，且 `accounted_within_budget` 与 `online_total_within55` 均为真、三项实际总成本全部已记录且在55秒/600LU/36000RHS总包络内，并且本列的实耗≤B的比例。

这描述**最终端点何时/以多少实耗可用**，不是优化过程中的首次命中时间，不从 candidate/best@K 或中间轨迹推断更早成功。每个合格最终端点只在它的已记录实耗处增加一次；B是累计预算阈值，不是虚构的评估时刻。预算上限不会被代填成实耗。总成本直接取 scorer，已含适用方法的共同 GN 准备；本脚本不再加一次。

未成功、无最终 top1、NOT_RUN 和实际缺失 score 行都保留分母；注册主方法缺失行在 methodcounts 中单列 `MISSING_SCORE_ROW` 数量，不伪称这些母体已经完成CIE或其他传统方法。缺失成本不补零，不给预算合格成功信用。仅质量合格却超预算或会计不完整的最终端点，只计入单独的 quality-only endpoint 数，不进入主曲线。保留两个预算门同时约束，避免 post 部分超限但总计暂未超55秒的端点被误作其原方法预算内成功。

## 固定主方法与输出

九个方法：Strict GN、Continue TRF、all-scales discrepancy、CIE dual cold、Voxel seed7、J-full seed7、G-weak seed7、G-weak+const seed7、recurrent feedback1 seed7。只用 FINAL_RUNS 已注册的这些精确 run 名；没有按 TEST 表现再选方法。

- `figures/generalization/completion_success_vs_actual_cost.png`：clean/noisy × wall/LU/RHS六格图。
- `completion_success_vs_actual_cost_CURVES.csv`：紧凑阶梯曲线阈值、成功数与固定108分母。
- `completion_success_vs_actual_cost_METHODCOUNTS.csv`：每方法/噪声的真实score覆盖、NOT_RUN、无最终端点、quality-only数、预算合格数、质量合格但非预算合格数、各项成本已记录/未知数。
- `COMPLETION_COST_STATUS.json`：实际输入、方法、精确定义与禁止推断范围。

混合 CPU/MPS/i9 下，wall曲线是本部署的实际端点费用，不能解释为硬件无关的算法速度。LU/RHS不包含NN固定矩阵乘法和卷积；一次性训练、setup及其成本仍由独立账单列出，不能据LU/RHS曲线声称所有方法等计算。此图不提供未知中断成本的完整总项目费用。

当前只有实现交付；实际图像生成、视觉验收与科学解释均留待父线程使用最终封存分数执行。
