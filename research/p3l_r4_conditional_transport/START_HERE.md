# Research P3L-R4：从材料风险改善到完整条件概率的实测边界

本轮真实训练、Maxwell 计算与两批新对象评价已完成。**监督成像 GO；完整后验仍未合格；固定核的数值接受决定加速 GO（限定范围）；动态非完整探索没有独立增量。**

512/256 个全新 parent 上，主要风险估计的全网格 complex MSE 为 **0.006806 / 0.007378**；冻结 R2 CNN 为 **0.007832 / 0.008323**，分别改善 **13.10% / 11.36%**。与本轮重新训练、数据更多的强 CNN 相比仍改善 **5.36% / 4.90%**。这些是监督风险动作，不是经过认证的 posterior mean。

[完整中文报告](RESEARCH_REPORT_ZH.md) · [网页](report.html) · [图库](gallery.html) · [科研裁决](GO_NO_GO.md) · [机器指标](RESULTS_AND_METRICS.json) · [复现说明](REPRODUCE.md) · [证据包入口](EVIDENCE_PACKAGE.md)

建议按以下顺序阅读：

1. [材料风险结果](MATERIAL_RISK_RESULTS.md)：为什么很差的全局生成均值仍能帮助局部 CNN，以及普通 Gaussian 控制为何重要。
2. [完整条件概率](CONDITIONAL_FLOW_RESULTS.md)和[校准审计](PROBABILITY_CALIBRATION_AUDIT.md)：真正的相对概率进步，以及尚未解决的严重 Maxwell 失配。
3. [物理计算复用](MAXWELL_REUSE_RESULTS.md)：相同 Full-MH 路径上的实际节省，以及缓存、近似和数值范围。
4. [几何实验](NONHOLONOMIC_GEOMETRY_RESULTS.md)：有限材料候选确实存在，四段路径却没有胜过直接提案。
5. [理论与方法](THEORY_AND_ALGORITHMS.md)、[完整失败记录](FAILURE_ANALYSIS.md)、[下一轮 Pro 交接](GPT_PRO_HANDOFF/START_HERE.md)。

历史 P1/P2/P3-L/R2 不改写。新对象在冻结方法和两批预测完成后才进入独立 GT 评分；现在两批标签均为 EXPOSED，不能再称未来盲测。原始失败、重复状态、拒绝、过时中期日志和费用仍完整保留。
