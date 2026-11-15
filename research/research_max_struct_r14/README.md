# MAX-STRUCT R14 — Dual-Metric Maxwell Recovery Dynamics

本轮研究原A+非局部材料反馈、两块metric不一致与实际恢复之间的关系。完成576次全程反演及完整state干预，不覆盖R4–R13。

| 入口 | 内容 |
|---|---|
| [中文报告](RESEARCH_REPORT_ZH.md) | 最终判断、真实恢复、边界与论文资格 |
| [简明交接](CHAT_HANDOFF.md) | Chat快速入口 |
| [Pro研究资料](PRO_HANDOFF.md) | 精确公式、完整state与开放问题 |
| [理论推导](THEORY_DERIVATIONS.md) | 理想与真实程序的区别 |
| [轨迹分析](TRAJECTORY_ANALYSIS.md) | 因果干预、失败、非根力 |
| [实现/基线](BASELINE_AUDIT.md) | 真正的M0/M1、β、material导数、预算偏差 |
| [新24逐对象](analysis/independent24_v1/report.html) | 11方法、全部正负结果 |
| [同成本附注](analysis/pairwise_cost_independent24_v1/report.html) | 成对CPU与LU/RHS已付费前缀 |
| [复现](REPRODUCE.md) | 真实入口、数据恢复、封存和评分 |
| [成本](COST_LEDGER.md) / [失败](FAILURES.md) | 真实费用与所有偏差 |
| [英文论文稿](PAPER_DRAFT.md) | 主张由结果决定，非投稿完成声明 |

A+在独立24中为12/24，CIE与Dual各9；一致非局部11独有一个原栈全部失败的新高contrast恢复，DR在该组混合材料更强。简单晚期一致化未新增恢复。完整source/public、原始全字段候选/checkpoint、离线truth和Pro包在最终DOWNLOADS索引；material-only review不替代原档。

原用户任务与Pro ZIP读取范围分别保存；附带文档的命令与建议只作为研究来源，不自动当作新的用户指令。NN、实测或连续物理等未完成项目有明确边界，未用空白模板冒充执行。
