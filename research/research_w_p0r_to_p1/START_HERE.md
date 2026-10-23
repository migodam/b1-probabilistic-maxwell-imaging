# Research W：墙–目标补偿、概率分支与递进成像

本轮已从P0-R执行到受控P1-B。总体裁决为 **REDIRECT**：保留真实二维全波联合补偿见证和可核对的概率基准；当前一般未知墙体在线成像、递进调度及A17独立优势均未建立。最终费用与封账状态以 [CLOSED.json](CLOSED.json) 和 [资源总结](RESOURCE_SUMMARY.md) 为准。

一个具体结果是：同质量目标移动 **0.1λ₀**，由墙体实介电与损耗重分配补偿，两个完整joint均值的白化距离为 **0.798716**。两墙calibration距离却为 **3.560512**，因此这是joint补偿，不是calibration近等价对。独立粗细网格差均低于0.1。

必须同时保留一个范围限制：原 `minimum_feature=0.1` 没有完整数学定义，本轮没有追认原宽几何先验PASS。Oracle构造也不是在线发现。两墙材料均满足明确固定支撑家族的材料范围、被动性、面积、查询域及跨频规则。

随后另行冻结proper两分支prior，再生成独立观测。P1-A/P1-B各16个独立参数parents、每parent两份噪声，仅含两种固定墙几何。reference求积通过经验细化审核。固定observation000保留约 **54.235% / 45.765%** 的两个事件及条件不确定性；家族外质量为null。

在全部新P1-B观测上，多假设均值比actual-Full验证的连续点控制降低约 **37.03%** 位置平方损失，配对95%区间排除0；但不优于本批零观测prior均值。Progressive相对Direct SMC没有独立风险优势，且一项ensemble CDF误差超过0.02。全部失败留在总体统计中。

- [中文网页报告与实际材料图库](report.html)
- [完整中文报告（Markdown）](RESEARCH_REPORT_ZH.md)
- [最终科学裁决](SCIENTIFIC_DECISION.md) / [各阶段Gate](GATE_SUMMARY.json)
- [可交给下游重新加权的概率材料包](MaterialHypothesisPackage/README.md) / [EXAMPLE.json](MaterialHypothesisPackage/EXAMPLE.json)
- [全部费用与未用额度](RESOURCE_SUMMARY.md) / [原始追加账本](PHYSICS_LEDGER.jsonl)
- [复现：纯数组回放与需新授权的物理执行分开](REPRODUCE.md)
- [完整响应、材料、粒子及明确全场子集的证据Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-p0r-to-p1-20261024)

先读前两项即可理解物理问题、正负结果和十二问。P0-R审计、B2存在性、P0-C辨识性、P1-A概率正确性与P1-B调度贡献各自保存，不能把其中一项成功当作其余阶段PASS。历史目录、旧失败与已封账预算均保持原样。
