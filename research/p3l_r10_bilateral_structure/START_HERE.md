# R10-BI：从双边结构到实际费用与材料覆盖判断

**已完成本轮实验。** 保留固定形状的物理响应复用；停止扩大没有覆盖收益的双边提案。强TAP贡献仍HOLD。

- [完整中文报告](RESEARCH_REPORT_ZH.md)：真实物理、候选、条件任务、失败与十问。
- [科学决策](GO_NO_GO.md) / [机器可读结果](EXPERIMENT_RESULTS.json)。
- [公式集新增18条](../maxwell_formula_atlas/09_BILATERAL_SHAPE_AMPLITUDE.md)：推导、成立条件、计算后果。
- [新10对象材料比较](CANDIDATE_FRESH10_RESULTS.md) / [配对统计](STATISTICAL_ANALYSIS.md)。
- [10形状响应核](REDUCED_RESPONSE_FRESH10.md) / [836节点条件任务](PROBABILITY_WORKLOAD_REPLAY.md)。
- [独立QA](INDEPENDENT_QA.md) / [失败分析](FAILURE_ANALYSIS.md) / [复现](REPRODUCE.md)。
- [中文论文叙事](PAPER_NARRATIVE_V2_ZH.md) / [英文稿](PAPER_NARRATIVE_V2_EN.md) / [查新](PRIOR_ART_AND_NOVELTY.md)。
- [下一轮交接](GPT_PRO_HANDOFF/START_HERE.md)。

固定形状32次查询，rank64双门相对直接256阶LU冷形状费用降低51.82%；一个开发条件积分工作量节时59.22%。这些不是全系统成像加速。Fresh10上Full66/Shape/Bilateral覆盖分别3/10、3/10、2/10，原锚点2/10，无新覆盖优越性结论。

全部新账目：306096 RHS256 +16928小系统RHS，9273 LU256，0新RHS576；历史264 RHS576 setup另列。完整精度数组、付费失败、源码与必要历史输入无损保留在同一B1 Release。最终公开版本与下载入口见 PUBLICATION.md。
