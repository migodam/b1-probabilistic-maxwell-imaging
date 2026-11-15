# 下载、恢复与使用边界

实际附件名、字节数和匿名访问结果统一见`DOWNLOADS.json`，不猜测尚不存在的URL。网页evidence入口列全部可用附件。

## 快速阅读和后续数学研究

先读中文RESEARCH_REPORT_ZH、THEORY_DERIVATIONS和PRO_HANDOFF。Pro有限ZIP包含原完整checkpoint、force数组和紧凑表，不包含GT、聊天缓存或第三方全文；成员清单与START_HERE说明哪些计算可直接重放。它不是全部576运行的替代档案。

## 完整数据档

- `dualmetric_pilot12_v1_part*.zip`及其ARCHIVE_RECEIPT：首12个曝光材料的48次运行。
- `r14_raw_dualmetric_remaining36_v1_part*.zip`、`r14_raw_self_only16_v1_part*.zip`、`r14_raw_stage_control16_v1_part*.zip`与`terminal3_PACKAGE_PLAN.json`：其余曝光36、self控制与stage控制。
- `r14_raw_independent24_v1_part*.zip`与`independent24_PACKAGE_PLAN.json`：新24×11方法的全部264原始运行。`independent24_v1_ANALYSIS_FULL.zip`仅为对应派生分析。
- `r14_local_fullstates_shared_readout4_v1_part*.zip`、`r14_local_fullstates_vector4_v1_part*.zip`及local_fullstates计划：Apple统一readout16和vector24原始字段。所有attempt保留，成功shared明确是attempt_02。
- `r14_sources_configs_public_part*.zip`和`r14_analysis_diagnostics_part*.zip`配最终PACKAGE_PLAN：本轮源/配置/公开观测、全部本地诊断和离线候选注释。大数值表放Release，不为缩小网页而删除。
- offline-truth独立ZIP保存已用于封存后评分的合成GT；只读原件归档，不属于online输入，不能混入Pro或公开观测输入。

每个组的所有part需齐备，不能将同名PACKAGE_PLAN相互覆盖。先给不同manifest单独目录，再按对应RESTORE说明恢复到新的外盘目录。文本是无损ZIP_DEFLATED；NPZ等二进制原字节ZIP_STORED。普通ZIP CRC/成员尺寸检查，不使用SHA。所有原始文件保留。

`src/package_evidence.py restore`适用于其PACKAGE_PLAN schema，pilot使用其ARCHIVE_RECEIPT原成员路径。它们重建R14 study目录；不在恢复后的冻结RESULTS上重新运行优化器。公开Git clone只含审阅表、源码和图片，大NPZ必须另外恢复。历史物理/算法依赖按REPRODUCE指定的真实旧Git版本取得。

不同位置的旧绝对路径seal不自动有效；保留旧seal，按REPRODUCE生成当前路径的未授权恢复记录，核对后另存授权。不要用虚拟D盘alias或全局替换路径绕过此边界。阅读已有表图无需重新评分或真值。

## 成本和科学判断

出版下载核对只证明附件可访问；它不重新验证数学、物理或材料恢复。完整solver费用、离线分析和传输费用分开，失败也保留。无lossy数组压缩、无GT在线选根、无修改历史R4–R13结论。

## 本轮真实发布状态

DOWNLOADS记录60个匿名可读附件，总计11,808,228,940bytes。小Pro ZIP为R14_PRO_PACKET_v1.zip；独立离线GT为R14_SYNTHETIC_OFFLINE_TRUTH_v1.zip。source+diagnostics共用根PACKAGE_PLAN.json及RESTORE.md；terminal3、independent24、local_fullstates各用各自前缀manifest，不混淆。源码ZIP中的HTML为有时间戳的出版前预览，可从随包Markdown重新build；最终公开网页及最后的费用/访问回执以Git主分支为准。
