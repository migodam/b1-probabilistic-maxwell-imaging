# MAX-STRUCT-R6：A路径扩大对照与FFT-TSOM接入

本轮回应“扩大实验数据、比较优化后的A与CIE/Dual-CIE、检查能否建立在FFT-TSOM之上”。保留R1–R5证据，在独立工作树中完成21材料×4臂主比较与2材料×2臂FFT工程复跑。

- [中文完整报告](RESEARCH_REPORT_ZH.md)
- [Chat简要入口](CHAT_HANDOFF.md)
- [数学与组合边界](THEORY_DERIVATIONS.md)
- [实际数据范围](BENCHMARKS.md)、[数据合同](DATA_CONTRACT.md)
- [成熟基线审计](BASELINE_AUDIT.md)、[NIE审计](NIE_BASELINE_AUDIT.md)
- [FFT实现](FFT_TSOM_INTEGRATION.md)
- [固定实验协议](BENCHMARK_PROTOCOL.md)、[研究过程](RESEARCH_LOOP.md)
- [费用](COST_LEDGER.md)、[失败](FAILURES.md)、[复现](REPRODUCE.md)

`src/`和`experiments/`为真实入口；`RESULTS/`保存全部终点、历史与封存后评分，`FIGURES/`为汇总图；`sealed_truth/`在首次反演完成前由独立端保管，发布后已属于公开数据，不能再次称盲测。没有测试GT在线选根。
