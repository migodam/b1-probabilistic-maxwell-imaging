# R13 表格聚合附件

来源：development_v1, development_v2, representation_v1。所有数值来自已授权封存后评分表；本脚本未打开材料/真值NPZ或调用求解器。

保留全部估计器与缺失分母，按物理域、host、网格及采集条件分层。样本单位为材料对象，重复运行/初值/策略不能当作新增独立对象。配对变化与bootstrap是描述性统计，不新增成功门。

条件策略使用全48原G19分母：未触发保留原primary，触发使用实际rescue估计器；缺失触发结果保持不可用。每种部署策略收费包括全48parent加实际触发rescue；竞争策略总费用不能相加冒充实验消费。Apple本地参考与i9分开。

原始表：ALL_ESTIMATOR_ROWS、STRATIFIED_SUMMARY、PAIRED_OBJECT_CHANGES、PAIRED_OBJECT_BOOTSTRAP及SOURCE_RUN_COSTS。图只绘实际评分/费用记录，没有重构或模拟材料图。最终科学解释、负结果讨论与发布由父线程完成。
