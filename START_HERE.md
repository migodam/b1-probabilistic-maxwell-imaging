# 最新研究：B1 × C1 Joint P2

[中文报告](research/joint_p2/RESEARCH_REPORT_ZH.md) · [交互Top-K图库](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/joint_p2/) · [科学裁决](research/joint_p2/GATE_DECISION.md) · [完整证据下载](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/joint-p2-20261014)

同一连续材料prior保留实心/中心/环形空腔解释，提供事件质量、条件不确定性、full-Mie补偿与独立传感器更新。新增3459/4000次物理attempt；C1未达到完整连续posterior目标。下面原有内容保留为历史研究。

# 历史入口：B1 E0

已完成零新增物理的连续posterior摘要修复：[E0入口](research/b1_next/e0/START_HERE.md)、[独立中文报告](research/b1_next/e0/RESEARCH_REPORT_ZH.md)、[资源与原始数组](research/b1_next/e0/results/SUMMARY.json)。旧报告/图/384次账本保留；新阶段公开更新到同一仓库。

以下为历史固定模板baseline入口。

# B1：多假设 Maxwell 概率成像

这是一个实际运行的3D矢量Mie概率成像研究原型。输入16个复散射通道，在三个固定同心材料模板和连续复介电常数先验下，计算各拓扑事件的后验质量、代表性实/虚部材料切面与不确定性。

- [完整中文研究报告](RESEARCH_REPORT_ZH.md)：背景、数学目标、两轮实验、失败归因、结论。
- [在线重建图与概率面板](https://migodam.github.io/b1-probabilistic-maxwell-imaging/)。
- [离线面板](dashboard.html)；[主要重建图](figures/fig01_competitive_representatives.png)。
- [理论](THEORY.md)、[结果](RESULTS.md)、[科研裁决](DECISION.md)、[研究日志](RESEARCH_LOG.md)。
- [复现方法与可执行入口](README.md)；[原始概率表](results/CASE_PROBABILITIES.csv)。
- [逐次物理资源账本](results/PHYSICS_LEDGER.jsonl)、[运行账本](results/JOB_LEDGER.jsonl)。

竞争案例：中央空腔50.64%，中央空腔＋环形空腔40.43%，实心剩余质量8.93%。Top-2不重新归一化。百分比是连续材料拓扑分支的近似后验质量，单个精确材料图像没有非零质量。模板外几何的质量尚未估计。

结论：受限Top-K原型已实现；一般几何恢复、严格全域概率认证、真实仪器统计校准与强TAP原创贡献尚未建立。原始参考失败、错误重放、历史负结果与访问缺口均保留。384/384次新物理评估全部记录；没有大型NN训练或旧locked数据访问。

## 给网页版GPT的阅读顺序

先读本入口和中文报告，再读THEORY.md与RESULTS.md；涉及概率可靠性时查阅DECISION.md、results/refinement/SUMMARY.json及research/delegated/probability_review/REVIEW_NATIVE.md。不能把公开发布状态、软件测试或development重复观察当作独立科学确认。
