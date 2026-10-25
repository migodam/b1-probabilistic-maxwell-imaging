# Research P3L-R5：A/B 实验结果

**A：局部概率矩有进步，完整覆盖与成本门未通过，REDIRECT。B：统计压缩通过开发门，独立确认未运行，CONDITIONAL GO。没有新的 Full-66 posterior 或 Maxwell 专属风险贡献主张。**

- **实际Maxwell**：3个二维子模型、6个parent、18份配对观测；32,705状态attempts、291,784 RHS，预算内封账。
- **实际训练**：42个网络fit、24个最终选定RRR模型（另含网格与OOF重复拟合）；B新增Maxwell费用0，训练约460秒MPS。
- **最有价值的概率发现**：实际decoder存在正体积rank0宿主层；接缝Full-Maxwell响应斜率约2.997，不能用固定正秩chart覆盖。
- **最有价值的风险结果**：214,272参数统计压缩器在已有VAL256保留约92.96%的fusion增量；CPU单线程、batch64、每parent的额外provider耗时约原K32的2.54%，完整附加链约54.57%。共同物理锚点成本仍在。
- **明确缺口**：新256+256同配置确认受32,768 RHS硬预算限制，NOT_RUN；严格CPU/MPS逐像素重放容差有失败，已保留。

[完整中文报告](RESEARCH_REPORT_ZH.md) · [网页](report.html) · [最终裁决](GO_NO_GO.md) · [冻结协议](FROZEN_PROTOCOL.md) · [复现说明](REPRODUCE.md) · [机器结果](RESULTS_AND_METRICS.json)

建议阅读：

1. [A：完整实验与失败](stage_a_response_coordinates/RESEARCH_REPORT_ZH.md)及[18观测逐行表](stage_a_response_coordinates/OBSERVATION_SUMMARY.csv)。
2. [B：矩、训练、成本与边界](stage_b_material_risk/REPORT_ZH.md)及[真实VAL图库](stage_b_material_risk/GALLERY.html)。
3. [数学与数值核验](THEORY_AND_NUMERICAL_CHECKS.md)、[独立审计](audit/FINAL_SOURCE_QA.json)和[资源账本](RESOURCE_LEDGER.json)。
4. [下一轮Pro交接](GPT_PRO_HANDOFF/START_HERE.md)：具体新现象、未解决条件与三个理论问题。

旧R4测试数据没有被重新称作本轮盲测；开发集收益不等于新对象泛化。风险图、概率样本、合法代表继续分开。所有NOT_RUN和旧失败保留。
