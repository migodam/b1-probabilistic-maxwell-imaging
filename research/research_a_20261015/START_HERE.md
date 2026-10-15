# Research A：发现、验证并交接多个材料解释

**已经完成真实实验：球对称 Mie 平台可以保留多种相容结构，但没有建立“历史反馈生成器”的额外优势；A30 困难二维实验没有新增结构收益。** 本轮应保留结构分层搜索、候选交接和新增测量研究，停止当前反馈/ROM/复杂学习模型的优越性主张。

|第一屏问题|本轮实际回答|
|---|---|
|找到了几个物理不同的世界？|32份不同测量的静态有限编辑 archive 合计193个材料世界，每份测量3–10个，平均6.031；不是同一测量下193个。按冻结结构尺度选出的代表平均5.0个。|
|比普通 multi-start 多几个？|相同每对象12次新 full endpoint 下，multi-start 平均2.562个，静态方法净多3.469个/对象；32个 archive 条目总数净多111。它是计数差，不是证明111个均为对照永远找不到的分支。|
|Best-of-4 是否改善？|离线归一化材料 MAE：0.136339→0.112949，改善17.1558%；配对差95%区间[-0.035784,-0.012455]。|
|在线选一张是否改善？|没有。按最小 full residual 选图的 MAE=0.152732，基础图=0.151770。Best-of-4只用于离线评分。|
|新解如何产生？|预先覆盖不同径向结构区间，保持编辑区域，用有界ε/损耗和共享gain补偿，再用原 full-Mie验证。独立同区间重启能复现端点，增益主要归于结构分层搜索。|
|adaptive 是否更好？|未证明。Best4比静态方法坏0.7436%，平均搜索wall更高；旧C1 MDN proposal也没有胜过multi-start。|
|困难二维条件是否通过？|没有。8个历史暴露条件中5个无合法解；余3个的Best4各方法相同。没有≥6格结构变化的新相容族。|
|TSOM/ROM/feedback有不可忽视增量？|未建立。本轮没有公平端到端TSOM/ROM增量对照；同结构区间的独立重启几乎复现有限编辑端点。|
|下一轮继续什么、停止什么？|保留结构分层baseline和B测量交接；停止当前历史自适应、feedback/ROM优越性升级与E3。|
|强算法/TAP论文是否成立？|本轮不足。Best4改善成立，但预注册coverage/成本联合条件未满足，且归因对照不支持新反馈机制；E2是明确负结果。|

阅读顺序：[中文研究报告](RESEARCH_REPORT_ZH.md) → [方法比较](METHOD_COMPARISON.md) → [消融与归因](ALGORITHM_ABLATIONS.md) → [失败记录](FAILURE_CASES.md) → [三轮研究迭代](ITERATION_LOG.md)。完整裁决、理论边界和下一阶段授权建议由主线程写入本目录对应文件。

[Top-K材料图库](index.html) 展示基础图、在线选择、离线评分用真值、不同候选及散射拟合；球对称中心切面只是显示。精确三维球层材料保存在 [MaterialHypothesisPackage](MaterialHypothesisPackage/) 的径向区间复ε数组中。32份真实交接包均为 **candidate_only**；概率和未覆盖质量不可量化，不使用 residual softmax。

可核查案例见 `confirm_00`：中心孔半径约0.0983R的材料与接近无中心孔的材料，材料L2距离0.04088、白化预测距离0.04067，均满足同一观测阈值。它们的ε接近2.2+0.16i的支持边界，gain约1.276；必须保留这一边界压力，不能把它称为高后验概率。

E1新预算与历史P2独立，不使用旧余额。E2使用328/384个full-state attempts、9214/9216个复LU RHS；只余2 RHS，不能再算一个8源状态。最终费用以本目录资源汇总和append-only账本为准。GPU训练、新大型bank和非同心三维campaign均未运行。

发布说明：沿用现有B1公开网页；固定入口见文末，不新建仓库。


公开入口：[可浏览材料图库](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_a_20261015/) · [本次证据Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-a-20261015)。所有历史阶段在同一B1仓库保留。
