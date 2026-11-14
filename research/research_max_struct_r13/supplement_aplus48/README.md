# R13 A+48与统一技术栈补充（已完成）

[网页](index.html) · [中文报告](RESEARCH_REPORT_ZH.md) · [Chat交接](CHAT_HANDOFF.md) · [全部图](index.html#galleries) · [原始附件](ARTIFACTS.md)

原A+按原配置补齐相同48对象；另完成48份R7最终读出控制、64份冻结NN接口、40份native状态重放。新增184运行，使用已有公开材料，不增加独立样本、不称新盲测。原R4–R13结果保留。

全48 primary平均NRMSE：CIE .19454、Dual .17985、G19 .16758、DR .17672、A+-MM .15210、A+-Arc .15214。同实际CPU前缀：Arc .15282、MM .15746；A+对CIE/Dual的平均材料增量保留。Arc最终几乎同MM而平均CPU少约5.1%。

NN最终主要帮助G19、对DR有小收益、对两版A+几乎无新增收益。正确性、材料误差、raw/corrected、费用和接口分别报告；没有GT在线选择。

## 继续阅读

- [轨迹和失败](TRAJECTORY_ANALYSIS.md)：完整状态与物理误差、current梯度变化。
- [Theory2素材](THEORY2_MATERIALS.md)：公式、假设、已否定的简单解释与后续问题。
- [历史比较](HISTORICAL_COMPARABILITY.md)：为何不能把跨轮排名当同一比赛。
- [成本](COST_LEDGER.md)、[失败记录](FAILURES.md)、[复现](REPRODUCE.md)、[英文论文补充](PAPER_ADDENDUM.md)。

数据与状态经公开验收及父审封存后评分。大数组与完整候选池通过附件下载，紧凑表格和论文图随仓库发布。私有连接脚本及本机D:路径别名不发布。任务执行状态见TASK_RESUME，不能替代科学证据。
