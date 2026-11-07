# Recurrent-depth VAL plotting handoff

状态：**绘图入口已实现，未执行诊断或绘图，未等待输入文件。** 本任务只读了 `experiments/recurrent_iterations.py` 与预声明 `audit/RECURRENT_ITERATION_DIAGNOSTIC.json`，并做静态语法检查；没有调用 physics、网络推理、训练、SVD 或 GT 文件。

实际输入文件名已从诊断代码确认：

```text
results/recurrent_iterations/<checkpoint.parent.name>/trajectory.csv
results/recurrent_iterations/<checkpoint.parent.name>/RECEIPT.json
```

新绘图脚本 `experiments/plot_recurrent_iterations.py` 只读这两类结果以及固定的六案例声明。它不读取 checkpoint 权重、stage NPZ、观测、anchor 或 truth 数组，也不会启动诊断。三个目录的 CSV/收据缺一时，它立即返回 `NOT_READY`（退出码2）、打印缺失路径，不生成输出、不启动计算、不进入等待循环。

## 父线程在实际诊断完成后的命令

从仓库根目录执行：

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 "/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" research/research_gn_nn_preservation_r4/experiments/plot_recurrent_iterations.py
```

默认目录映射为：

| 图中模式 | 结果子目录 |
|---|---|
| refreshed-q | recurrent_refreshed_q_seed7_cascade2 |
| quadratic feedback | recurrent_feedback1_seed7_cascade2 |
| hinge feedback | recurrent_feedback_hinge_seed7_cascade2 |

诊断代码按 checkpoint 的父目录命名输出。如果实际父目录不同，明确传入全部三个映射，不按结果好坏选择目录：

```sh
python research/research_gn_nn_preservation_r4/experiments/plot_recurrent_iterations.py --runs refreshed_q=REFRESHED_DIRECTORY feedback1=QUADRATIC_DIRECTORY feedback_hinge=HINGE_DIRECTORY
```

`--input-root` 可指向另一个实际诊断结果根目录。模式名来自显式目录映射，图表 JSON 同时保留收据的 checkpoint 路径，绘图脚本不会通过加载模型来猜测模式。

## 一张科学图的固定布局

输出：`figures/recurrent/recurrent_depth_val.png`。

图为六行 × 六列；每一行是一个预声明观察，不取均值、不筛选 best case，每个子图都保留三模式：

1. digit2、parent val_p0028，clean 与 noise=.001；
2. digit8、parent val_p0036，clean 与 noise=.001；
3. islands、parent val_p0044，clean 与 noise=.001。

这合计是**三个母体的六个观测**，不是六个独立母体；全部 A=.2。完整 obs_id 保留于端点 CSV/JSON，顺序来自固定声明。

每一列直接使用已保存字段，不重新计算任何范数或物理量：

| CSV 字段 | 图中含义 | 原诊断采用的归一化 |
|---|---|---|
| material_RMS_step | 每步材料更新 RMS | 已存 √mean((m_t−m_{t−1})²)，单位为材料 contrast RMS |
| current_relative_step | 每步网络辅助电流更新 | 已存 ‖j_t−j_{t−1}‖/‖j0‖；不是实际 Maxwell current 的相邻步差 |
| constitutive_relative_defect | 网络状态 constitutive defect | 已存 ‖j_t−m_t q_t‖/‖j0‖ |
| true_FullVIE_relative_data | 当前材料的实际 Maxwell 数据残差 | 已存 ‖F(m_t)−Y‖/‖Y‖ |
| physical_current_gap_relative | 实际 Maxwell 电流与网络辅助电流之间的差 | 已存 ‖J(m_t)−j_t‖/‖j0‖ |
| physical_internal_field_gap_relative | 实际 Maxwell 内场与网络场之间的差 | 已存 ‖Q(m_t)−q_t‖/‖Q0‖ |

这里 `true_FullVIE_relative_data` 中的 true 是“真正执行原始 Maxwell forward”之意，参照的是公共观测 Y，**不是隐藏 GT 电流/材料误差**。所有物理计算属于父线程运行的既有诊断；本绘图脚本的数值求解调用为零。收据原有 wall、NN inference 与 physical_costcounts 被原样收集到端点 JSON，不用预期调用次数代替实测记录。

注意本诊断的 constitutive defect 分母是固定 anchor ‖j0‖，与其它 VAL 收据中分母为 ‖j_aux‖ 的同名诊断可能不同；不能跨表直接比较数值大小。

每个面板在训练深度2处画竖虚线，并对2之后区域着色；图标题明确：**VAL public-input recurrent-depth diagnostic，NOT independent TEST；3–16步超出训练深度。** 没有新的成功阈值或 PASS 线。正值使用 log 显示；若实际记录精确零或负值，则采用仅供显示的 symmetric-log 轴（线性范围由实际最小正值确定），不把显示参数当作科学阈值。

## 实际迭代2与16端点输出

绘图成功后另写：

```text
results/recurrent_iterations/recurrent_depth_endpoints.csv
results/recurrent_iterations/recurrent_depth_endpoints.json
```

CSV 为18个模式×观测组合、每组合六指标，共108行；各行包含实际 iteration2、iteration16、带符号 delta16_minus2、原始字符串及端点可用性标签。JSON 为18个组合的紧凑嵌套结构，保存六指标、各自缺失迭代、完整案例身份和输入收据来源。

只做端点标量相减，不做公式验证、相邻步恒等式检查、收敛拟合、任意 PASS/FAIL 阈值或 best-case 排名。iteration2或16缺失时端点保留 `NOT_RECORDED`；非有限原值保留 `NONFINITE_RECORDED` 和原始字符串，JSON 不用非标准 NaN/Infinity 冒充数值。图中缺失点留断线并标记，非有限点显式标记，不把不良案例删掉后重画。意外的额外观测、非VAL行或重复迭代会明确拒绝，避免混入其它结果。

## 解释边界

曲线趋平、更新变小、端点残差下降或16步保持有限，只是这六个预声明 VAL 观察、这些已训练参数在超训练深度上的经验轨迹。它们**不保证固定点存在、唯一性、收缩、GT收敛、任意迭代稳定性或新 TEST 泛化**。图和 JSON 均保留这条边界。最终解释由父线程完成。
