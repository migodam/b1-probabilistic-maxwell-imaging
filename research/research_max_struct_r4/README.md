# MAX-STRUCT-R4：CIE公式机制的推导、验证与边界

本轮是冻结R3之后的独立续轮，保留R1/R2/R3原始证据。完成用户引用session的A/B/C验证、结构修正、新鲜确认、独立评分、跨平台/网格检查和费用审查。

**结论：A在一个公开困难开发例实现材料恢复，但六对象确认与成熟CIE同为4/6，没有新成功率或等成本优势；B/C已测版本未突破。尚未得到成熟、可推广的优越新算法。** 正面机制、反例和失败后修正都保留，不把运行完成当科学通过。

- [中文完整报告](RESEARCH_REPORT_ZH.md)
- [原候选A/B/C逐项汇报](ABC_VALIDATION_ZH.md)
- [数学推导与假设](THEORY_DERIVATIONS.md)、[理论状态](THEORY_STATUS.md)
- [六对象逐例评分](results/CONFIRMATION_SCORES.md)、[独立科学审查](audit/FINAL_SCIENTIFIC_REVIEW.md)
- [研究循环](RESEARCH_LOOP.md)、[数据合同](DATA_CONTRACT.md)、[基线审计](BASELINE_AUDIT.md)
- [复现](REPRODUCE.md)、[总费用](COST_LEDGER.md)、[失败记录](FAILURES.md)
- [Chat读取入口](CHAT_HANDOFF.md)、[下一轮交接](CODEX_TO_PRO_HANDOFF.md)

物理范围：真实有限2D scalar-TM Hankel VIE，24×24主网格和一个32×32开发诊断；没有把它改名为3D。无NN、无β扫参、无测试GT在线选根。六对象确认在算法冻结后生成，GT直到18预测封存后才开放；后续修正明确为已暴露开发版本，不改变原分母。

复用现有成熟CIE/Fourier-TSOM和原物理corrector。来源基版本为R3发布 `1bd07e49`，仅作普通git来源记录，不做SHA核验。执行源码快照、配置、原始数组与收据均保留。

网页：[中文报告与图像](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r4/)。

Chat原始入口：[CHAT_HANDOFF.md](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/CHAT_HANDOFF.md)。
