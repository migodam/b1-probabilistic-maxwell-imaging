# 原始附件与下载

六个附件已发布到[同仓库Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r13-aplus48-20261114)，总计1215382781字节。逐个通过匿名公网请求读取HTTP响应和ZIP头；没有计算或验证SHA，也不把4字节读回称为完整文件完整性证明。

|附件|字节|内容/还原位置|
|---|---:|---|
|[aplus48_dynamics8_full_results.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/research-max-struct-r13-aplus48-20261114/aplus48_dynamics8_full_results.zip)|404129586|40记录副本/完整native状态：解压到supplement_aplus48/dynamics_run。|
|[aplus48_extension32_full_results.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/research-max-struct-r13-aplus48-20261114/aplus48_extension32_full_results.zip)|125854107|原A+缺失32：解压到supplement_aplus48，保留RESULTS/queues。|
|[aplus_stack112_execution_sources.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/research-max-struct-r13-aplus48-20261114/aplus_stack112_execution_sources.zip)|74364|112实际source/runtime/config（20文件）；配合原C2公开输入与legacy依赖。|
|[aplus_stack112_full_results.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/research-max-struct-r13-aplus48-20261114/aplus_stack112_full_results.zip)|458815258|48固定Arc＋64NN接口：解压到supplement_aplus48/stack_run。|
|[max-struct-r13-aplus48-analysis.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/research-max-struct-r13-aplus48-20261114/max-struct-r13-aplus48-analysis.zip)|183412845|全评分、索引、gradient banks、原生标量表、图、独立审查：解压到supplement_aplus48。|
|[max-struct-r13-aplus48-public-inputs-and-sources-v2.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/research-max-struct-r13-aplus48-20261114/max-struct-r13-aplus48-public-inputs-and-sources-v2.zip)|43096621|原48/NN16公开输入及本补充源码、配置和审计：按research/...布局解压到仓库。最新网页/文档以Git为准。|

旧原16、四基线48、NN/G19/DR/direct及原权重/真值从[原R13附件](../ARTIFACTS.md)复用，不混成新增消费。原物理/单位/噪声合同不变，复现见[REPRODUCE](REPRODUCE.md)。

本机D:硬链接别名、私有SSH控制脚本与传输日志不发布。记录副本不能替代主预测；所有GT评分为封存后离线。来源commit与Release目标在GitHub可审查，正文的后续出版说明不会反向改变冻结数组。
