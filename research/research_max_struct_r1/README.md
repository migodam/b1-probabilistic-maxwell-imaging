# Research MAX-STRUCT-R1

2026-11-02 完成并按后续授权公开的最小机制研究。**没有新的成像优势主张。**

[中文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r1/) · [Chat 可读取的原始报告](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r1/RESEARCH_REPORT_ZH.md) · [完整复现包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r1-20261102) · [发布与封存边界](PUBLICATION.md)

- [中文完整报告](RESEARCH_REPORT_ZH.md) · [HTML 阅读版](REPORT.html)
- [主线程数学推导](THEORY_DERIVATIONS.md) · [理论/实验状态](THEORY_STATUS.md)
- [数据物理合同](DATA_CONTRACT.md) · [成熟基线审计](BASELINE_AUDIT.md)
- [复现说明](REPRODUCE.md) · [成本](COST_LEDGER.md) · [失败和未运行项](FAILURES.md)
- [后续交接](CODEX_TO_PRO_HANDOFF.md) · [独立审查](audit/FINAL_INDEPENDENT_REVIEW.md)

F：实 SPD 解析同场反例、严格被动复对称有限 DDA 分支、源/接收三臂、有限噪声射线；公开开发证据，不是盲恢复。

M：缺兼容 H(curl)/弱残差合同，表示比较 NOT_RUN。T：未获得可用统一 k 和完整 δ，非空证书 INCONCLUSIVE。新鲜成像及 ROM NOT_RUN。

主要数据：`results/families/SUMMARY.json`、`results/families/RANK_CLASSES.json`、`results/backend_qualification/QUALIFICATION.json`；动力学：`research/delegated/audit/D1.json`、`D3.json`。标准运行/不运行回执：`audit/receipts/`。所有“PASS”应结合对应范围，运行结束不是科学通过。

已有基线公开出处：[冻结主仓库](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/7317601b43b914525c21f09e5a60020552ea677b)、[A17 仓库](https://github.com/migodam/a17-maxwell-certified-exchange)。这些是历史来源；本轮入口列于页首。
