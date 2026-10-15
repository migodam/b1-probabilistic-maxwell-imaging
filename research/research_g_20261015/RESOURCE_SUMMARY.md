# Research G 独立成本

|费用|实际|硬上限|
|---|---:|---:|
|新full-state attempts|291|384|
|complex full-solver RHS|6936|12288|
|受控数值wall秒|76.305076|10800|
|CPU秒|72.441919|单列，不与wall互换|
|GPU|0|0|

三个独立base验证=3states/24primalRHS；72轨迹各raw+两步+独立验证=288states。共291次LU；2328primal RHS+4608receiver-transfer adjoint RHS=6936。没有隐藏JVP/VJP/newtruth solve。已保存background数组只读setup，wall仍计入，不称免费物理加速。

E0新state/RHS0。每个方法每轮36states/864RHS，公共setup/base只记一次；per-method随机顺序、原始wall及全部prefix见CSV。完整job wall包含准备/优化/序列化/独立full验证，nested ACTION不再相加。

总数值时间包含E0的独立回执、首次失败offline路径桥接、最终独立评分与75图、侧车导出及witness分析；子回执和rollup不重复计费。交互阅读/写代码/文档/上传时间不伪作数值CPU。当前物理solver异常0，数据不相容54端点、QP未收敛24条另列，所有尝试已付费。

剩93states/5352RHS，不用来救结果。T192states、旧A30/A/B/C1/P2账本不动，无预算转账。source/COST_LEDGER+FAILURE_LEDGER为追加记录，RESOURCE_SUMMARY为最新快照。没有SHA检查或发布前重跑实验。
