# Track B 交接摘要

已完成全部8类表示的 development、三 seed head/蒸馏修订、诚实 nested2fold18baselinefit 与成本复核。没有新物理或新 truth。

- 原 strongCNN / K32 / deterministic frozen / deterministic matched raw MSE：0.00746135 / 0.00714550 / 0.00704106 / 0.00703415。
- 最终214272参数 compact raw-y provider：raw MSE 0.00716774，保留原增量 **92.96%**；clipped保留 **91.72%**。
- 最终 canonical batch64额外provider费用 **2.54%**，完整provider+features+head **54.57%**。原provider压缩门在VAL上PROVISIONAL_PASS；完整pipeline≤1/3强门FAIL。
- 四类RRR在honest nested OOF仍变差；没有物理独立增量证据。oldfullTRAIN几何是y-only条件frame，最终halfTRAIN residual/fullTRAIN baseline差异已披露。
- 新物理RHS **0**；新256+256同anchor独立确认 **RESOURCE_LIMITED / NOT_RUN**（相同历史512context约1.935M RHS，独立cap32768）。旧已暴露R4 TEST未用于选择。
- 全部GPU训练fit 460.0s，whole-job保守预留 506.3s。一次初始5.77GB MC materialization失败保留，后续streaming/训练峰值<3GiB，未把首次资源门失败抹掉。
- 两个prospectivefinalists与三seedcheckpoint在SELECTION_FROZEN.json冻结：deterministic_matched、distilled_w128_n32768。**没有继续训练/物理计划**。

完整报告：[REPORT_ZH.md](REPORT_ZH.md)；可视化：[GALLERY.html](GALLERY.html)；逐模型与门槛：[RESULTS.json](RESULTS.json)；推断入口：[inference.py](inference.py)；[COST_LEDGER.json](COST_LEDGER.json)。统计区间来自选择用VAL，仅exploratory，不是独立确认。

复现门补充：CPU/MPS原3e-6门FAIL保留（max3.945e-6，meanriskdiff≤9.31e-10）；独立CPU1两次新Action全量重复逐数组相同PASS。retention98.333%bootstrap区间仍宽，92.96%仅VAL点估计。最终canonicalbootstrap、教师摊销、固定VAL0..3实际material图、phasepercell与RRRrank成本曲线均已补齐。
