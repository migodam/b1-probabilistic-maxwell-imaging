# Research T：Maxwell物理、基础反演与补偿成核

**已实际执行并封口。8/8条件有合格默认基础材料；40个固定结构提案产生32个合格补偿端点记录。数据相容成功，材料精度没有改善。**

|轮次|基础规则／主要变化|raw通过|补偿后通过|独立parent有成功|
|---|---|---:|---:|---:|
|1|最低残差base＋普通全材料补偿|0/16|15/16|4/4|
|2|同提案＋非core正材料相保护|0/16|10/16|4/4|
|3|最小N3改动的合格base；固定两个编辑|0/8|7/8|4/4|


跨轮按相同观测、固定材料距离与相位门槛构建archive，48个保存材料端点（16基础变体＋32补偿结果）全部保留。每份ring观测9／10／12／9个；每份weak观测2个。数量不是独立连通解分支、后验模式或完整覆盖认证。

- [中文研究报告](RESEARCH_REPORT_ZH.md) · [成像与散射图库](offline/index.html) · [裁决](GATE_DECISION.md)
- [G物理交接](T_PHYSICS_HANDOFF/README.md) · [默认base与48份真实材料索引](T_PHYSICS_HANDOFF/CANDIDATE_INDEX.json)
- [逐端点指标](PER_EDIT_METRICS.csv) · [离线成像误差](offline/OFFLINE_SCORING_REPORT_ZH.md) · [完整成本](RESOURCE_SUMMARY.md)
- [原始数组与可复现证据Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-t-20261015)。下载解压至仓库根目录，恢复相同路径；读取结果无需重跑。

**必须同时读的负结果：**最低残差base的平均complex MAE为0.147656，原N3为0.041756；改为最小材料改动选择后为0.050757。补偿候选的MAE约0.49–0.51，Best-of-4包含基础图时改善为0。通过噪声球不能代替成像准确、posterior质量或统计不可区分。

新授权消费192/192 full-state、5832/6144 complex RHS；受控求解jobs共213.213秒。GPU0，新3D0，旧余额0。剩余312RHS不能越过已耗尽的state预算。历史A30/A/B、Joint P2和A17冻结稿保持原样。

范围为历史暴露的四parent×weak/ring，24×24二维scalar TM；不是独立盲验证、连续网格证书或一般三维Maxwell成像。所有概率为null/candidate_only。公开状态与科研结论分别记录。
