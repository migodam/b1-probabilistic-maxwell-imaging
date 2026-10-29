# 从这里开始：P3L-R9-S

本轮已完成，包含32个新parent的完整Maxwell候选云和持续维护的公式集。

- [中文主报告](RESEARCH_REPORT_ZH.md)：主要结果、十个问题、范围与成本。
- [可搜索公式总集](../maxwell_formula_atlas/index.html)：122条推导，固定/变化量、优化含义、源码来源。
- [最终裁决](GO_NO_GO.md)：候选能力、GT覆盖、费用、条件概率分别判断。
- [32parent对照CSV](POSSIBILITY_VOLUME_COMPARISON.csv)与[逐parent评分](POSSIBILITY_PARENT_SCORES.csv)。
- [复现说明](REPRODUCE.md)及[完整证据包](EVIDENCE_PACKAGE.md)。
- [GPT Pro交接](GPT_PRO_HANDOFF/START_HERE.md)。

固定方案354个合格非anchor候选、动态362、随机103；主容差GT覆盖6/32、5/32、5/32。固定略省wall但未达20%；exact ROI核含setup节时58.65%。没有把候选权重称为准确posterior，也没有把尚不够高的覆盖称作高概率GT保证。
