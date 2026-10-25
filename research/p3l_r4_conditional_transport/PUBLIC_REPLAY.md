# 换一块磁盘后，如何重放 R4 已保存的风险分数

这是一条**离线、标签已暴露**的统计复现入口：只读取已保存的材料风险图像和 GT `chi`，不加载网络、不训练、不调用 Maxwell，也不重新选择方法。它不检查旧 manifest 中的绝对路径，因此无需恢复原 `/Volumes/...` 挂载点。

工具：[tools/replay_r4_public_risk.py](../../tools/replay_r4_public_risk.py)。只需 Python 和 NumPy；不需要 PyTorch、SciPy 或 A30 后端。该工具与本说明属于封包后追加的便携入口，请从最新仓库/补充材料取得；已有大包、原数组、manifest 和闭合科研 ledger 均未修改。

## 运行

完整证据包解压后的相关布局是：

```text
NEW_DISK/extracted/
├── workspace/
│   └── research/p3l_r4_conditional_transport/
│       ├── inputs/offline_only/R4L1_*.npz, R4L2_*.npz
│       └── results/fresh/{locked,confirmation}/
│           ├── BOUND_RISK.json
│           └── risk45/
│               ├── online/<method_id>.npz
│               └── scores/RESULTS.json
└── dependencies/...
```

```bash
python /path/to/repository/tools/replay_r4_public_risk.py \
  --evidence-root /NEW_DISK/extracted \
  --out /NEW_DISK/replay_run_001 \
  --offline-exposed-data
```

`--evidence-root`也可直接指向`extracted/workspace`、正常仓库根目录，或`research/p3l_r4_conditional_transport`目录。工具只尝试这几个固定相对布局，**不遍历166万个文件**。`--out`必须是新目录；缺少显式`--offline-exposed-data`确认则拒绝执行。

最小输入是上面列出的两个`BOUND_RISK.json`、两个闭合`RESULTS.json`、90份方法导出，以及768份GT文件。GT文件只读取`chi`；方法导出只读取`parent`、`prediction_raw`、`prediction_clipped`、`fallback_anchor`。不读取raw currents、模型权重、latent truth、CFM samples或原物理cache。完整分组包中，这些材料主要来自`compact_core`、`raw_neural_and_scores`和`new_R4_inputs`；统计重放不需要体积最大的`raw_physics_and_arrays`流。

## 输出与含义

- `METHODS.csv`：全部45方法在512/256对象上的raw/clipped全576-cell complex MSE、MAE及fallback数。
- `FAMILIES.csv`与两个`*_family_parent_losses.npz`：先在每个parent内平均seed **loss**，不是先平均图像。
- `COMPARISONS.csv`：原固定Spline对R2、strong CNN、Diag fusion和permuted fusion的比较，4000次父样本paired bootstrap，seed20261025099、98.333333%区间。
- `RESULTS.json`和`RECEIPT.json`：与紧凑闭合结果对账、读取文件数、实际wall/CPU/RSS、明确0 RHS/0训练/0推断/0新假设。

统计单位仍是512和256个独立parent，三个训练seed不会把样本量扩大为1536/768。主族仍按原validation规则固定为Spline，不因重放结果而换方法或改确认集规则。支持性比较保留原描述范围，不新增familywise或posterior声明。数值对账的`1e-12`仅用于独立浮点算式差异，不改变任何科学阈值或历史失败判定。

## 已实际验证一次

在现有闭合数组上完成一次运行：读取862个所需文件，失败/fallback为0；均值最大绝对差`6.94e-18`，比较/区间最大差`2.15e-16`。实际wall **1.321 s**、CPU **0.543 s**、RSS **84.4 MB**。没有Book，也没有向闭合科研ledger追加记录。

实际回执：[AUDIT/public_replay/run_once_20261025_v1/RECEIPT.json](AUDIT/public_replay/run_once_20261025_v1/RECEIPT.json)。这是保存动作的分数复现，不是重新完成训练、采集或盲测；两组标签现均为EXPOSED。最终科学解释仍以主报告及主线程审核为准。

## 如要另行开展新物理工作

这不属于本入口。证据包另含`workspace/src/p1_material`、`workspace/src/p3r2`、P1几何、冻结R2 baseline以及`dependencies/a30`下的A30源码与完整嵌套vendor。A30依赖其vendor的相对位置，不能只复制一个`physics.py`。

新的、另行授权的物理campaign应显式配置后端位置：R4 `core.Target`支持`B1_P3L_A30_ROOT`；部分历史驱动直接使用相邻`A30/a30-maxwell-material-manifold`布局，须保留该布局或写明新的路径适配。**不要静默重写闭合manifest的旧路径、解除原campaign关闭标记，或把本统计重放当成新物理资格。** 旧绝对路径保留作历史provenance；本工具完全不解引用它们。
