# 原168任务评分与能量分析：最终交付入口

**COMPLETE。** 168有效端点、24对象全图、0缺失/0评分错；material scorer先写公开seal再按root授权读取GT，之后energy诊断顺序完成。只新增24个OFFLINE_SCORING GT-current LU/288 RHS，没有补算端点或阶段current，没有改方法/预算/阈值。Fresh Twofold实际r=30的图标签已更正；历史DEV28未改。

- [主结果、成本预算、完整负例与有限噪声解释](FRESH_TWOD_REVIEW.md)
- [0001/0004/0022的native末stage、公开选择起点、corrector及current解释](THREE_CASE_INTERPRETATION.md)
- 全对象原表：`SCORES_FINAL.csv`、`OBJECT_WIDE_COMPARISON.csv`、`PAIRED_BOOTSTRAP_SUMMARY.csv`；全部103个material FAIL（含内部GlobalDual）在`ALL_MATERIAL_FAILURES.csv`。
- 预算：`BUDGET_SCORES.csv`保留2016单元，`BUDGET_SUMMARY.csv`为定点汇总；内部clock、physical LU、保守外部包络严格分列，不将内部曲线称equal-total-CPU。
- 能量：`delegated/twod_energy/scored/`内的全tau表、aux/physical/gap/cross及96份oracle；20个rank/mf空间cache共享给七臂，没有按GT选tau。
- 有限噪声双材料表：`CIE_GT_CN99_FEASIBLE_PAIRS.csv`；posthoc selector脚本/表：`posthoc_portfolio.py`、`POSTHOC_PORTFOLIO_OBJECTS.csv`。
- 费用入口：`POSTSEAL_ANALYSIS_COSTS.json`；原inverse完整child成本仍在`COST_ACCOUNTING.json`及原twod总账。

主要结论：CIE/Dual/A+global/A+source/Twofold30/G19成功数为11/10/8/10/8/8；新候选没有总体收益，Twofold/G19完整child CPU为CIE的1.2208/1.1899倍。0001是明显正例，0004/0022是明显反例；0000/0021主要是.15附近约.0038–.0061的NRMSE门变化，不能将所有退化夸大成大失败。

CIE22/24拟合noise但仅11/24材料成功；另外11个data-pass/material-fail对象中，GT与估计都合法且同在实际y的单case CN99球。这是有限噪声两合法材料实例，非精确多解、信息论下界或24联合99%；soft prior仍可偏好其中一个。

最低原Phi的CIE+G19纯posthoc诊断平均NRMSE降至.149599，但success为10/24，丢掉0000/0021并新增0001；两臂总成本2.1899×CIE，未选臂照付。它不是新确认算法或数学创新。

全24三页材料、energy绝对/比例、同case oracle、endpoint cost及clock明确的预算图均已目检，路径在`FIGURES/twod_comparison/`。Root之后的right3及conditional4未混入本冻结比较；需另行封存、授权和收费。本次现有分析到此完成。
