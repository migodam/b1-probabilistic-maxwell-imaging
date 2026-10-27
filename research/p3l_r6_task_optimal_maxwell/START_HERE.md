# Research P3L-R6：任务最优Maxwell材料计算

**本轮科学实验已完成并封账。** 6,860,984/8,000,000 RHS；512个新parent的两档噪声确认；33个新训练阶段；所有主要正负结果、强化控制、计时漂移与失败保留。

先读[中文完整报告](RESEARCH_REPORT_ZH.md)与[逐路线判定](GO_NO_GO.md)。原冻结Full基线范围内，候选在两噪声批次节时34.6%/32.6%并通过1%风险NI；更强Full-q低噪NI仍HOLD，跨窗口20%下界未全部确认。完整66维posterior未认证，TAP结论REDIRECT。

- [新噪声独立确认](STAGE_D_LOCKED_GENERALIZATION.md)
- [数学、条件与反例](THEORY_AND_COUNTEREXAMPLES.md)
- [A：分层概率](STAGE_A_STRATIFIED_PROBABILITY.md)、[补偿](STAGE_A_MAXWELL_COMPENSATION.md)、[任务积分](STAGE_A_TASK_INTEGRATION.md)
- [B：功能学习](STAGE_B_FUNCTIONAL_DISTILLATION.md)、[seed分歧](STAGE_B_SEED_VARIANCE_AND_ACTION_RISK.md)、[任务响应矩](STAGE_B_MAXWELL_CONTEXT_ABLATION.md)
- [C：锚点费用](STAGE_C_PHYSICS_ANCHOR_COST.md)、[停止与预算](STAGE_C_RISK_AWARE_STOPPING.md)、[风险—成本](RISK_COST_PARETO.md)
- [机器结果](RESULTS_SUMMARY.json)、[资源](RESOURCE_SUMMARY.json)、[闭账](CLOSURE.json)、[失败分析](FAILURE_ANALYSIS.md)
- [统计合同](SHARED_STATISTICAL_CONTRACT.md)、[历史证据](PRIOR_EVIDENCE_AND_FAILURES.md)、[复现](REPRODUCE.md)
- [GPT Pro交接](GPT_PRO_HANDOFF/START_HERE.md)

本轮没有重跑旧R5-B512，也没有改写其至少80%收益保留率HOLD。新风险图像不被称为posterior样本或Maxwell合法代表。原始精度数组以独立无损证据包交付，报告和关键图可直接阅读。
