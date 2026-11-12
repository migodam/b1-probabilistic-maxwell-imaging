# R13 原始附件与实际发布清单

完整公开入口由 [DOWNLOADS.json](DOWNLOADS.json) 的实际上传与公网size-readback记录生成。报告与图可直接在仓库/网站阅读；大数组通过现有项目Release提供，不另建仓库。

## 附件组织

|附件|内容与用途|
|---|---|
|`confirmation_core48_full_results.zip`（按顺序分片）|48×9的材料、候选池、物理状态、成本和阶段记录|
|`confirmation_mechanism16_full_results.zip`|固定16的A+/NN/AP/Anderson/在线机制等176份|
|`confirmation_online_controls16_full_results.zip`|48份额外公开三臂控制|
|`confirmation_rescue48_full_results.zip`|9×5救援、全48触发规则/旁路出处|
|`confirmation_local_reference16_full_results.zip`|M4充分原生日程96份参考|
|`confirmation_local_trf16_full_results.zip`|16新对象+1旧B的Full-J TRF|
|`physics_v1_recovered_full_results_v2.zip`|规范60份物理结果：18原Austria+42独立恢复|
|`physics_v1_full_results.zip`|第一次物理尝试：18完成和42导入失败，保留失败证据|
|`physics_recovery_v1_full_results.zip`|原42份恢复来源，与规范包重叠，不重复计费/计样本|
|development/representation完整归档|已暴露对象、反例与完整current消融；不混成独立确认|
|`runtime_*.zip`|原执行冻结代码/公开在线输入与精确版本|
|`max-struct-r13-postseal-inputs-sources.zip`|当前研究源码、Pro来源、公开输入、封存GT、来源及政策索引|
|`max-struct-r13-postseal-analysis-tables.zip`|全量离线评分、付费前缀、oracle等大型表和验证数组|

确切名称/大小/parts以DOWNLOADS.json为准；上表是内容索引，不冒充所有资产名的穷举。NN权重沿用仓库R8已跟踪文件，没有本轮新训练。

## 读取与复现

- [复现命令和依赖](REPRODUCE.md)
- [紧凑逐对象CSV](analysis/review_tables/README.md)
- [全部真实图集](GALLERY.html)
- [费用](COST_LEDGER.md) / [失败](FAILURES.md)

下载分片只按公布字节与顺序连接，不运行SHA256校验。恢复工具要求新/空输出位置；旧数组保留。不同物理尝试分别解压，不能把失败包覆盖到恢复后的规范结果。

## 科学与隐私边界

release只含本项目研究来源及自身结果，排除私有Chat原文缓存、SSH/密钥/临时签名URL、private operational目录、伪造集成夹具、环境缓存和外部完整论文。原始Fresnel文件不二次分发，适配输入及出处/取得方式保留。公开GT已经属于封存后的离线材料；不可拿它回灌在线候选选择。

材料池与代表性全物理场可用，不代表所有native PR梯度/方向与每一步辅助current都已保存。完整状态保存范围见`delegated/full_current/AUXILIARY_STATE_COVERAGE.md`。发布后的新费用回执/说明可能晚于ZIP快照，以仓库最终文档与费用清单为准。
