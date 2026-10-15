# Research B：从多个材料世界到有价值的新测量

[在线中文报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_b/) · [完整报告](RESEARCH_REPORT_ZH.md) · [全部原始证据ZIP](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-b-20261015)

**实际完成E0、E1、E2和四例真实A包交接。** 新Maxwell 3284/4096，GPU0。复杂选测在32parent上未稳定超过固定3.2；加入独立corrector才明显改善四例A的候选支持。新观测history参考仍UNRESOLVED，不宣称full posterior或一般3D成像。

建议按以下顺序阅读：

1. [中文研究报告](RESEARCH_REPORT_ZH.md)：首屏八个问题、真实结果和解释。
2. [科学判断](GATE_DECISION.md)、[下一阶段](NEXT_STAGE_PROPOSAL.md)、[失败](FAILURES.md)。
3. [实际ranking](ACTION_RANKING_RESULTS.md)、[材料MAE](MATERIAL_RISK_AND_MAE.md)、[补偿机制](SHARED_COMPENSATION_TEST.md)、[漏候选](CANDIDATE_COVERAGE_STRESS.md)。
4. [费用](RESOURCE_SUMMARY.md)、[复现](REPRODUCIBILITY.md)、[源资产](SOURCE_MANIFEST.json)、[协议差异](PROTOCOL_RECONCILIATION.md)。
5. [A→B数据交接](AcquisitionRecommendation/README.md)、`code/ab_adapter.py`及tests；`results/REMAINING_WORLDS.json`保留多个候选，权重none。

结果入口：`results/SUMMARY_METRICS.csv`、`PARENT_METRICS.csv`、`PAIRED_COMPARISONS.csv`；E3另有`results/E3/SUMMARY_METRICS.csv`及`PARENT_METRICS.csv`。大JSON/数组/forward caches随证据ZIP提供，仓库小索引链接不等于缺失实验。

历史P2、A、B1/C1配置/账本不改。源Research Pack原指令从属于本次用户任务；原schema保留。没有SHA检查、额外发布前实验或审批流程。
