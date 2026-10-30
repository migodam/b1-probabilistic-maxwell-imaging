# R12-BF 阅读入口

先读 [中文研究报告](RESEARCH_REPORT_ZH.md)，再看 [九数字实际图](DIGITS_123456789_RECONSTRUCTION.png)、[GO/NO-GO](GO_NO_GO.md) 与 [费用](RESOURCE_SUMMARY.json)。

本轮实际完成四个暴露P1开发对象的八臂比较、九个冻结N数字的六臂有界材料探索，保存全部失败和候选。主结论是：精确双边分解组织了原Full-J信息，但分层刷新没有改善主材料覆盖或数字双孔恢复。停止递归扩大，不以未认证posterior为失败借口。

- [双边公式](BILATERAL_MODE_DERIVATION.md)；[nested SVD冗余](NESTED_SVD_REDUNDANCY_AUDIT.md)
- [固定点谱](FIXED_POINT_SPECTRUM_AUDIT.md)；[一阶横向校正](TRANSVERSE_CONTINUATION.md)
- [分层算法](HIERARCHICAL_ALGORITHM.md)；[视觉阅读](DIGITS_VISUAL_REVIEW_ZH.md)
- [材料覆盖CSV](MATERIAL_COVERAGE_COMPARISON.csv)；[逐层CSV](DEPTH_ABLATION.csv)
- [全部原始账本](COST_LEDGER.jsonl)；[数值汇总](RESULTS_SUMMARY.json)
- [复现与资产](REPRODUCE.md)；[下一轮交接](GPT_PRO_HANDOFF/START_HERE.md)

新费用163,994 RHS、6,974 LU，均在300,000/12,000上限内。没有新的独立确认、H3、Full66后验或精确材料纤维结论。R11的超预算失败仍属于R11，没有覆盖。

完整数组以同一B1仓库的 [Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr12-bilateral-fixedpoint-20261030) 分片提供；阅读包不能冒充全部原始状态。
