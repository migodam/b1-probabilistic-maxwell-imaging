# 原材料先验内的竞争解释：本轮已结束

**结论：NO_GO_FOR_GENERATOR_COMPARISON。** E0 与 E1 已实际执行；0 组竞争 pair，E2 不启动。不是 GLOBAL_UNIQUE，也不是原先验不存在多解的证明。

- [中文研究报告](RESEARCH_REPORT_ZH.md)
- [图片与精确数据入口](gallery.html)
- [阶段裁决](E0_E1_E2_DECISION.md)
- [资源](RESOURCE_SUMMARY.json) / [失败](FAILURE_LEDGER.jsonl) / [执行偏差](EXECUTION_DEVIATIONS.json)
- [冻结协议](PROTOCOL_FROZEN.json) / [精确指标](METRICS.json)

主要发现：历史 T/G 的 noise-ball 端点没有原 A25/A30 的固定增量相位及生成器支持；noise-ball 通过也不等于点似然接近。原 TRAIN64 的 8,064 个缓存比较及本轮 45 states / 680 RHS 未得到原先验内的合格竞争见证。N3 仍保留强 Top-1，未改历史结果。所有 posterior_mass=null。

仅本目录为新阶段。T/G/A30 历史预算、truth、reference、裁决保持原样。所有旧 development 条件均为历史已暴露，未读取 locked truth。

物理代码：`src/prior_constrained/`。新求解器代码未从头开发，复用原 A30 frozen adapter。已付费数组在 `offline_only/`，只能用于复核本轮离线见证构造，不得给未来 online 搜索泄漏端点、参数、clean means。

RUN_CLOSED.json 已阻止继续求解。本轮不会在后台继续消耗剩余资源。
