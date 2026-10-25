# R5源码、阅读包与原始证据

本轮在现有B1仓库交付，不新建研究仓库。静态网页、Git代码和Release运输不改变科学门槛。

- [原始证据Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr5-parallel-ab-20261026)
- [紧凑Pro阅读包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/p3lr5-parallel-ab-20261026/Research_P3LR5_Codex_Results_and_Pro_Handoff.zip)
- [机器可读归档清单](EVIDENCE_ARCHIVES.json)
- [精确依赖路径与尺寸](audit/EVIDENCE_INPUT_REFERENCES.json)
- [源码与冻结配置](SOURCE_INVENTORY.json)

## 五个独立无损归档

|文件|内容|
|---|---|
|[r5-track-a.tar.zst](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/p3lr5-parallel-ab-20261026/r5-track-a.tar.zst)|A全部实际观测、积分节点/权重、材料/预测、局部RTO及失败版本|
|[r5-track-b.tar.zst](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/p3lr5-parallel-ab-20261026/r5-track-b.tar.zst)|B实际模型、teacher缓存、42次网络fit、24个最终选定RRR模型（网格与OOF重复拟合另计）、VAL预测和统计|
|[r5-shared-physics.tar.zst](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/p3lr5-parallel-ab-20261026/r5-shared-physics.tar.zst)|共享物理适配器、累计预约/费用原始日志、qualification、closure、未执行确认种子|
|[r5-parent-probes.tar.zst](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/p3lr5-parallel-ab-20261026/r5-parent-probes.tar.zst)|父线程曲率界/逆动作与真实Full-Maxwell接缝测试|
|[r5-essential-reused-inputs.tar.zst](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/p3lr5-parallel-ab-20261026/r5-essential-reused-inputs.tar.zst)|复现必要的已付费历史背景、选定冻结模型、训练/验证与context，以及所需源代码；不是新物理数据|

合计压缩字节为806,494,073，约0.806GB；按同名JSON清单列出每个原始成员。归档逐成员完整字节比对通过，没有SHA-256检查、精度转换或删除原科学源文件。复用输入包只包含必要文件，不是175GB R4大包副本。

## 恢复

归档成员路径均相对于`Gaussian/`。先恢复到新的空目录核对清单，保持B1/A30 sibling布局，再按[复现说明](REPRODUCE.md)配置。不要对历史工作树盲目解包覆盖。原始字节保留并不要求inode/ctime相同。

阅读包只包括报告、紧凑指标、关键图、代码和交接，不含完整训练大数组/checkpoint；完整内容按上表分卷获取。模型推理可以复现，但付费Maxwell阶段已关闭，不因下载包自动获得新计算授权。

文件存在、网页可达或上传成功，不等于posterior正确、独立确认完成或IEEE TAP贡献成立。科学资格以[主报告](RESEARCH_REPORT_ZH.md)与[裁决](GO_NO_GO.md)为准。
