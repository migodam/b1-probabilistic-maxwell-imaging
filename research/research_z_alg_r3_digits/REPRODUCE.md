# 代码、输入和复现边界

从本目录运行。完整可下载包与 Git 仓库都包含公开物理算子、匿名观测、冻结源码、冻结预测和评分后披露的材料；`private_truth` 与生成 seed 不发布。数值重求解使用完整 576 阶 VIE，机器负载改变墙钟停止位置，故不保证逐位重现原预测。原始冻结预测是报告的实际证据。

## 1. 环境

Python 依赖为 NumPy、SciPy、Pillow、matplotlib、threadpoolctl，见 `requirements.txt`。不需要 CUDA、GPU、神经网络或字体下载。已生成数据的反演与评分重放都不需要字体文件。新生成字形时使用本机字体及代码中的备选列表。

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

CLI 自动把 BLAS/OpenMP/Accelerate 工作线程设为 1，单次方法无需另启并行。原实验为两个独立单线程 worker；进程 wall 与 solver 数值 wall 分别记录。

## 2. 从一幅匿名观测重新反演

下面三条命令只读公开算子与观测，不读任何材料或数字标签。输出写入新目录，程序拒绝覆盖旧文件。

```bash
python run_solver.py --public public_data/holdout/obs_20ca8c325456.npz --operators public_data/operators.npz --method fullj_native --config configs/main_methods.json --output /tmp/r3_digits_native.json
python run_solver.py --public public_data/holdout/obs_20ca8c325456.npz --operators public_data/operators.npz --method fullj_reciprocal --config configs/main_methods.json --output /tmp/r3_digits_reciprocal.json
python run_solver.py --public public_data/holdout/obs_20ca8c325456.npz --operators public_data/operators.npz --method globaldual --config configs/main_methods.json --output /tmp/r3_digits_globaldual.json
```

主合同是每方法每图 **90 秒 solver 数值预算**，不是 CLI 默认的 DEV 75 秒。必须传入 `configs/main_methods.json`。读盘、Python 启动及写 JSON 位于数值计时之外，原批次另存完整外层进程时间。公共粗层在每个方法内重新计算并付费；不是使用知道 GT 的预计算图像。

九个输入全部在 `public_data/holdout`。可依次运行以上入口，保持冻结数据、物理算子、配置一致。生成或反演全部九幅无需数字标签；披露后的排列只用于结果展示。

### 计算预算补足与普通伴随控制

原主协议在目标生成前冻结，但其中 GN 有1000次迭代上限。观察到一些普通臂未花满90秒后，研究者在**尚未揭晓任何主 GT** 时另冻结以下两个补充，保留原始27份预测。它们不能冒称与原协议一样在数据生成前已冻结。

```bash
python budget_fill.py --original results/holdout/fullj_native/obs_20ca8c325456.json --public public_data/holdout/obs_20ca8c325456.npz --operators public_data/operators.npz --output /tmp/r3_digits_native_remaining_budget.json
python data_adjoint_control.py --public public_data/holdout/obs_20ca8c325456.npz --operators public_data/operators.npz --output /tmp/r3_digits_ordinary_adjoint.json
```

补足只使用 `90−原数值耗时`，不是再给90秒。普通伴随控制从相同公共粗层重新付费，使用原数据目标的 L-BFGS＋Full-J，以检查改目标与改优化器的混淆；它自身总预算为90秒。追加模型准备与输入装载、GN迭代上限及继续时重启的细小差异，均在独立审查中保留。

## 3. 不增加物理计算的评分复现

```bash
python replay_scores.py --split holdout --output /tmp/r3_digits_score_replay.json
python replay_scores.py --split dev_v2 --output /tmp/r3_digits_dev_replay.json
python replay_scores.py --split budget_filled --output /tmp/r3_digits_budget_filled_replay.json
python replay_scores.py --split data_adjoint --output /tmp/r3_digits_adjoint_replay.json
```

此命令只读 `truth_after_freeze` 中无 seed 的材料/掩模和 `experiments/frozen_predictions`。它重新计算 RMSE、MAE、前景/背景 RMS、SSIM、IoU 和 oracle best@K，并检查在线 top1 与冻结物理分数的最小候选相同。新增 Maxwell LU/RHS 为 0；这是**已披露真值的事后重放，不是新的盲测**。

`experiments/score_predictions.py` 是原始受控盲测评分器，读取当时私有 mapping 并先要求预测快照存在；公开包请使用上述独立重放入口，不需要恢复或伪造私有 seed。

## 4. 审核真实盲法及完整物理

按时间核对 `audit/HOLDOUT_FREEZE.json`、`audit/EXPERIMENT_GENERATION_HOLDOUT.json`、`audit/EXPERIMENT_BATCH_HOLDOUT.json`、`audit/EXPERIMENT_PREDICTIONS_FREEZE_HOLDOUT.json` 和 `audit/EXPERIMENT_VALIDATION_HOLDOUT.json`。算法/超参数先冻结，再生成九个新目标；27 份预测先快照，再打开 GT。完整预测、失败、wall、CPU、LU 与 RHS 不因截图好坏筛选。

独立验收解总场系统 $(I-GD_m)Q=E_0$，优化解电流系统 $(I-D_mG)J=D_mE_0$；两者是相同离散 VIE 的不同组装，不是跨离散或连续 PDE 验证。验收还用已有 LU 的 LAPACK 条件估计，不做完整逆矩阵/SVD。

## 5. 报告和归档

`aggregate_costs.py` 重建已知成本账本；缺失的旧 DEV 失败费用保留 `UNKNOWN`，不补造数字。`package_release.py --build` 仅打包现有文件，不增加物理解；它排除私有 GT、seed、缓存和压缩包自身。

网页构建器 `build_site.cjs` 只渲染 Markdown，保留公式。冻结的源码和预测目录不会被渲染器改写。科研图由冻结评分结果生成，GT/重建共享 `[0.5,20]` 色阶。
