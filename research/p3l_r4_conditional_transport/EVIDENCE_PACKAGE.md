# 阅读包、完整原始证据与复现

本轮继续使用原B1仓库，不新建项目。公开状态、最后读回时间与远端名称/大小核对，以最终`PUBLICATION_RECEIPT.json`为准；打包完成、上传完成、公开可访问和科学资格是四个不同事实。

## 先读研究结果

[中文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p3l_r4_conditional_transport/) · [实际数组图库](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p3l_r4_conditional_transport/gallery.html) · [匿名Markdown](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/p3l_r4_conditional_transport/RESEARCH_REPORT_ZH.md)

[既定Release入口](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-p3lr4-20261025)包含两种用途不同的材料：

- `Research_P3LR4_Complete_Results_and_Handoff.zip`：适合下一轮GPT Pro阅读的中文报告、指标、关键图表、R4及P1核心源码和数学/失败交接。它不包含全部大型权重与物理缓存，不以“阅读包”代替原始证据。
- 完整原始证据：**120分片，175,313,076,233字节**。含实际模型、输入、生成样本、完整材料与prediction、current/field/J、路径、拒绝重复、失败和原始账本。不是只打包成功图。

## 原始封存与后补文件的边界

完整分卷以提交`60066394392fe04fb6a5b30371ffa34c978ade9d`及科研账本sequence **1522869**为来源边界。实际选入字节及mtime另有逐文件清单，提交号不被当作全部未跟踪科学数组的内容证明。打包COMPLETE并释放锁后，没有再修改原始分卷。

`Research_P3LR4_Postpack_Supplement.zip`单独保存后来的发布说明、训练解释、标准TV/KL下界、可移植风险重放工具与交付重试记录。它更新文档/入口，不改模型、真值、原始评分或封账科学ledger。补充包内说明原文件与后补文件的版本，不能把后来生成的内容冒称归档当时已有。

最终上传及公开读回回执在Release小附件和本仓库后续交付提交中；它们可以晚于阅读包/补充包产生，不用反复改写已封存大包以追求自引用。

## 下载与恢复

每个分片不超过1.5GB。六类是同一恢复树的不同tar流，不是六份可替代的数据。四位分片编号按字典顺序拼接后解压，不能逐片解压。

|流|分片数|字节数|
|---|---:|---:|
|compact_core.tar.gz|1|1,323,321,353|
|raw_physics_and_arrays.tar|113|169,457,408,000|
|raw_neural_and_scores.tar|3|3,071,047,680|
|reused_inputs.tar|1|227,850,240|
|derived_training_inputs.tar|1|1,208,289,280|
|new_R4_inputs.tar|1|25,159,680|

精确文件名/大小见Release附件`R4_FULL_EVIDENCE_PARTS.json`与[本地可读清单](AUDIT/delivery_layout/R4_FULL_EVIDENCE_PARTS.tsv)。完整命令和依赖布局见[恢复说明](AUDIT/delivery_layout/RESTORE_LAYOUT_REVIEW.md)。请解到全新目录，保留原始分卷；建议至少约400GB可用空间以容纳下载、约166万文件和文件系统开销。流式解包不额外保存一份巨大拼接tar。

无需下载最大的物理流就能先做[材料风险离线重放](PUBLIC_REPLAY.md)：核心、neural及新对象输入三类包含所需862个文件。工具仅需Python+NumPy，已实际复核512/256对象、全部45方法，均值差至多6.94e−18，比较差至多2.15e−16。该重放不运行网络或Maxwell，不是重新盲测。

## 依赖与完整性范围

`R4_DEPENDENCIES.json`区分固定A30来源、实际P1/R2运行模块、历史baseline和系统依赖；`R4_PACK_RECEIPT.json`保存完整打包回执。Python/PyTorch/NumPy/OS本身不打包，实际版本见[RUNTIME_VERSIONS.json](RUNTIME_VERSIONS.json)。旧manifest保留历史绝对路径；新的路径映射应明确记录，不静默改写旧身份。

使用精确文件名、字节数、mtime和远端uploaded状态核对。未做SHA检查，也不把元数据核对称为密码学内容证明、连续PDE认证或新的科学复现。所有原始本地文件与失败记录继续保留。

## 费用

科研封账为10,719,752 complex RHS、337个关闭作业；[科研费用](RESOURCE_SUMMARY.md)不因出版而重置。打包、网络上传、重试、恢复说明与后补统计解释另列`DELIVERY_LEDGER.jsonl`及交付回执，仍受本轮总资源上限约束。网络等待不混入算法性能表，公开成功也不改变FULL后验未合格等科学裁决。
