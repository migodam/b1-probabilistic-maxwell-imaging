# 复现与证据读取

本目录是独立 P3-L campaign，不修改旧 P1/P2/F3 及其账本。首先读 START_HERE.md、概率合同与最终科学裁决，再查看具体作业的 PROTOCOL、SOURCE、CHAIN 和 CLOSED 回执。网页可达、代码正常退出或归一化密度都不等于完整 posterior 已被验证。

## 环境与依赖

使用 Apple Silicon CPU、Python、NumPy、SciPy、Matplotlib；神经模型另用 PyTorch。本轮实际训练为 CPU，GPU 使用为零。源码和冻结环境记录保存在 `AUDIT/` 及每个作业的 SOURCE 中。

Full solver 复用现有 A30 的 `FrozenMaxwellAdapter`，不是新的 solver。其版本及历史背景数组来源见 SOURCE_INVENTORY 和原始 provenance 记录。按该记录准备公开 A30 checkout，通过 `B1_P3L_A30_ROOT` 指定；不要把它指向修改过物理目标的另一个实现。P1 decoder 源码保存在同一 B1 的 `src/p1_material/`。

## 原始证据

Git 保存报告、配置、图、源码、诊断摘要与小型成本快照。完整 chains、逐次尝试/端点、模型权重、背景/current/LU 缓存及追加式 `RESOURCE_LEDGER.jsonl` 放同一个 B1 repository 的 P3-L evidence Release。分卷包使用路径、文件大小、分片偏移及 Git revision；没有 SHA-256 流程。

恢复工具是 `tools/package_p3l_evidence.py --restore NEW_DIRECTORY --bundle-dir DOWNLOADED_BUNDLE`。目标目录必须不存在；它重组完整大文件，逐项检查路径、大小和偏移，不运行物理。不要直接改写原 campaign 的数组或账本。

## 运行入口与含义

所有入口位于 B1 根目录的 `tools/`。`run_p3l_job.py` 执行基础反演/完整目标链，`run_p3l_hmc.py` 和 `run_p3l_spine.py` 是已保存的独立 HMC 对照；warmup、production、RNG 与真实成本以各作业冻结协议为准。`run_p3l_surrogate_sequence.py` 只执行已登记的两臂代理力试验。

`run_p3l_joint_gaussian.py` 计算声明的联合 Gaussian reference，不能当 Full posterior；`run_p3l_legal_reference_projection.py` 只投影其材料均值，不能自动赋予 Maxwell 相容性。Full validation 与独立 DEV evaluator 分进程。

`run_p3l_locked_imaging.py`、`run_p3l_locked_importance.py` 和 `evaluate_p3l_locked_imaging.py` 分别是冻结的在线估计、Full importance 校正与独立离线评分。配置和代码版本在任何 locked 观测生成之前保存；原始 128 parent 一个都不删除。不要根据真值重新选 checkpoint、链、半径或估计器。

`summarize_p3l_campaign.py` 与 `publish_p3l_progress.py` 只读取保存的摘要及小型轨迹，不运行 Maxwell 或读取 locked truth。它们不能自动授予 posterior/reference 资格；最终状态只转述主线程 CAMPAIGN_CLOSURE。

## 费用与重放

每个 Full 状态立即保存 q、材料、散射 prediction 和 backward residual；truth、anchor、检查点另存 current/field/LU。大规模链不永久保存每个中间 LU，缓存的有效性和恢复费用都有独立回执。拒绝样本保留重复状态。

原封账物理作业拒绝覆盖或自动重跑。读取、统计和出图可以从保存数组复现；新的 Full 重算仍需真实新增授权，不能删除 ledger、改 run ID 或新建 checkout 来伪装旧费用未发生。本轮未耗尽额度也不自动延续到下一个研究阶段。

原始中断、数值失败、未合格参考、update0 神经检查点及未执行模块都属于证据。请保留它们，并同时区分数学正确性、数值近似、经验混合诊断和独立材料风险。
