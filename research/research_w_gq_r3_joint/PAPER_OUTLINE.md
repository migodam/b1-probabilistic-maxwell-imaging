# 论文定位：先证明成像增量，再决定是否立项

本轮不把“联合 GN + Schur + SVD + 全波回退”命名为新的 Maxwell 理论。A17/A17-X 已有的材料注入、双侧反馈、paired 更新后果、有限提案条件 action 与计费边界均作为前置成果；[归属表](research/delegated/sources/A17X_REUSE_CROSSWALK.md)逐项列明。

## 可以支撑正文的实际内容

1. **问题与证据身份。** 静态、受限六宏格未知墙体与连续椭圆未知目标，三个频率、四源十六接收。主指标为目标位置、形状和材料，不以墙体精确恢复作为要求。
2. **公式而非命名。** 从 Born 双 Green 积分、固定对比源 Green、投影式 TSOM，推到实际 FEM 材料载荷及墙 nuisance 消元。明确本实现只是 SOM-inspired 载荷 chart。
3. **可运行方法。** 测量初始化、四起点连续有界联合优化、完整 Maxwell 接受与测量规则选支。公开普通 Full-GN 强基线和全部失败支。
4. **决定性对照。** 相同未知材料坐标、prior、初始化与上限下比较 Full-GN、fixed chart、adaptive chart；另列先墙后目标、已知墙 oracle 和 prior-only。fine-grid结果与盲测评分分开。
5. **为何固定谱不够。** 分开参考谱尾、当前材料激发与传播变化；用已经付费的同状态 Jacobian 检查目标 Schur 后果，不重复制造机器精度代数实验。
6. **成本与限制。** rank16 chart 的传播RHS并不天然少于16接收伴随。若额外审核/回退不能带来同质量净收益，保留普通求解器，删除加速主张。

## 暂不能作为论文原创点

标准 Schur 消元、current SVD、resolvent identity、普通多起点FWI、全波下降接受、条件最小二乘材料读取；以及A17-X已经推导过的材料激发、反馈与更新后果框架。F35/F36的目标 nuisance 专门化来自标准包络求导，本身不具有已确认的新颖性。

## 真正仍需决定的贡献

只有在独立材料/几何范围内，付费的目标切向相关模式选择比普通完整/加权伴随求解器获得可重复的目标质量或端到端成本优势，才有理由形成新的算法论文。少量盲测成功可以作为可运行证明例，不能替代一般准确性或论文级确认。

本轮最终科学裁决与实际数字以 [REPORT_ZH.md](REPORT_ZH.md) 和 `GO_NO_GO.json` 为准。若无独立增量，本文件是后续论文的证据缺口清单，而不是宣称已经完成 IEEE TAP 论文。

## 现象先例与系统贡献不是同一个判断

用户独立提出的方向应准确保留。2016年TSOM穿墙实测已有目标恢复而墙体失真的现象先例，这只限制“首次现象”的主张，不自动覆盖材料约束下的目标分辨率、有限非线性补偿、适用域条件或新的计算算法。更严谨、更成体系可以是贡献，但须写成具体可检验命题及新增证据，不能仅凭篇幅或独立想到自动判定原创。现有实验尚未支持当前adaptive实现超越普通Full-GN；不等于整个目标充分传播路线没有理论价值。

由独立Astra撰写的完整读者叙事、摘要和审稿反驳见 [PAPER_NARRATIVE.md](research/delegated/narrative/PAPER_NARRATIVE.md) 与 [SCIENTIFIC_COUNTERARGUMENTS.md](research/delegated/narrative/SCIENTIFIC_COUNTERARGUMENTS.md)。最近邻原文与数据兼容性见 [文献对照](research/delegated/literature/LITERATURE_AND_DATA_ZH.md)。
