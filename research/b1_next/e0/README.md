# 复现B1 E0

在本隔离worktree根目录安装原requirements。主入口：`python -m src.b1 e0`。已完成阶段禁止原地覆盖；独立复现需用新output目录并复制AUTHORIZATION.json（新增物理授权仍0）。不会导入/调用物理solver或写旧账本。

已付费观测成像：`python -m src.b1 infer --model continuous --observation results/cases/gain_marginalized/posterior.npz --output research/b1_next/e0/results/new_observation_output`。输出目录必须不存在。新refined/continuous入口使用同一连续density的moments与直接CDF；full-grid保留legacy节点诊断。NPZ必须有gain_sd，JSON可明确选择known/marginalized gain。

实际执行程序与原始数组见src/b1/continuous.py、e0.py、e0_checks.py及本目录results。报告程序为tools/e0_report.py，仅读取已完成结果。禁止在旧384目录再次运行audit/pilot/refine，不删除或降低账本。


独立零物理重放示例（从worktree根执行）：

```bash
mkdir -p research/b1_next/e0_reproduction
cp research/b1_next/e0/AUTHORIZATION.json research/b1_next/e0_reproduction/
python -m src.b1 e0 --output research/b1_next/e0_reproduction
```

该重放仍不会增加或重置物理计数，CPU/polynomial成本另记；不要覆盖本次原始结果。
