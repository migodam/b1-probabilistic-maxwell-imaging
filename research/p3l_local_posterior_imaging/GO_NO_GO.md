# 最终裁决：REDIRECT

科学裁决与发布状态分开。本轮计算已结束，未用余额封账。

| 路线 | 决定 | 直接证据 |
| --- | --- | --- |
| weak-density × strong Gaussian成像 | NO-GO | 128新parent MSE .0615021 vs .0110309，差95%区间完全为正 |
| 固定/联合GN Gaussian材料平均 | NO-GO当前近似 | DEV32均恶化，Full权重退化 |
| 合法reference投影 | NO-GO修复效果 | DEV32 MSE .0208401 vs .0138212 |
| FULL概率参考 | REDIRECT | 11组四链0合格；SMC祖先1；128 IS均ESS<2 |
| 解析材料矩 | 保留工具 | supplied Gaussian内数值核验通过，不宣称Full准确 |
| 当前CNN/MLP/MDN/Flow | NO-GO当前训练设置 | 15次validation均选择update0 |
| decoder material-response代理力 | NO-GO当前统计成本增量 | forward误差下降但MCSE不改善，wall更高 |
| A17 TM24加速 | NO-GO当前增量 | 同空间同精度，setup/bounds/fallback无总wall优势 |
| posterior distillation / RL | NOT_RUN | 没有合格teacher/统计奖励参考 |
| 强IEEE TAP算法论文 | NO-GO当前证据 | 没有可信Full概率与Bayesian材料风险增量 |

值得保留的是概率合同、原始失败、解析矩和可复现benchmark。下一步应先解决Full参考与非Gaussian条件，而不是增加网络容量。若不能建立参考，则应终止当前局部Gaussian概率成像机制，不宣称整个Bayesian逆散射不可能。

[主线程关闭记录](CAMPAIGN_CLOSURE.json) · [完整报告](RESEARCH_REPORT_ZH.md) · [独立统计复核](AUDIT/locked_results_review/REVIEW.md)
