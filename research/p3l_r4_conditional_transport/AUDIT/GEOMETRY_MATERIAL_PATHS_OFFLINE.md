# 已封账几何材料路径：隔离离线评分

先验证两次在线物理作业 COMPLETE/no_GT 和不同 PID，再在本独立 BookIII 分析进程读取原前8DEV真值。已评分 7688 个已付费端点、1791 个完成的有限路径；缺失 0 项。没有 Full、LU、RHS 或训练。

每条路径列出原锚点、终点和已采样路径点到真值的材料 RMS=||chi−chiGT||/24，以及实际有限材料移动与匹配标记。最接近真值的点仅为离线事后 oracle 描述，不反向挑选候选、修正在线图像或改变教师资格。Heun 四个校正端点与四个 predictor 付费场探针分开，predictor 的较近距离不冒充物理路径终点；Euler 仍保留已知 O(h²) 自身场缺陷，不能当成精确交换子或理论反证。

COMPLETED_FINITE_PATH_SCORES.csv 列原路径所有中间 Full-J/最终 Full 成本。原三档FD共享480RHS，Heun小步门禁原1728RHS，原在线上下文实际费用也分列；setup-inclusive列是原共享研究取得成本加路径，不是宣称每方法都必需同样冷成本，各路径共享项不能相加。ALL_PAID_ENDPOINT_MATERIAL_SCORES.csv 保留包括 FD 和未完成末段在内的全部端点，未把部分路径当完成。

GROUPED_DESCRIPTIVE.csv 按原算法/维数/秩/半径汇总，不以重复配置增加独立父对象数；不同实际移动不能假称同材料距离比较。输入都是已暴露开发对象，没有后验参考或连续路径证书。完整机器证据位于 results/offline_only/geometry_material_paths_v1/RESULTS.json；严禁把这里的真值误差或 oracle 索引作为在线特征。

附表为 full66/秩4/请求RMS=.02 的完整八对象切片，描述性展示，全部配置仍保留在机器表。共同锚点平均GT-RMS=.1157463942。

|研究/臂|实际末端RMS|末端GT-RMS|最近已采样GT-RMS（含锚点）|路径RHS|
|---|---:|---:|---:|---:|
|dynamic_geometry_first/circuit_forward|0.050984|0.129439|0.112707|128|
|dynamic_geometry_first/circuit_reverse|0.061469|0.128305|0.113666|128|
|dynamic_geometry_first/direct_bracket_forward|0.020000|0.117463|0.114189|8|
|dynamic_geometry_first/direct_bracket_reverse|0.020000|0.117420|0.114386|8|
|dynamic_geometry_first/dynamic_weak|0.020000|0.116507|0.114296|48|
|dynamic_geometry_first/frozen_weak|0.020000|0.116356|0.114326|8|
|dynamic_geometry_first/ordinary_field_curvature|0.020000|0.120581|0.115384|8|
|dynamic_geometry_first/ordinary_field_secant|0.020000|0.118037|0.114696|8|
|dynamic_geometry_heun_r1/Heun_forward|0.023066|0.118472|0.112964|288|
|dynamic_geometry_heun_r1/Heun_reverse|0.022874|0.116448|0.113186|288|
|dynamic_geometry_heun_r1/direct_forward|0.020000|0.117465|0.114178|8|
|dynamic_geometry_heun_r1/direct_reverse|0.020000|0.117429|0.114392|8|

最近距离的候选集合含原锚点，所以它逐对象不大于锚点距离是定义保证，不能由此推导算法收益。不同方法支付点数不同，事后取最近点也有不同oracle机会；均不是在线可部署选择。

后续距离门禁复核见 GEOMETRY_MATCHED_DISTANCE_GATE_REVIEW.md。上表未通过匹配距离的平均GT误差只描述原始失败结果，严禁据此给算法排优劣；Heun的.0229/.0231与直接法.0200不是同距离比较。新的固定非GT门禁先保留全部计划尝试/失败分母，再给出双方均满足距离条件的原同父对象配对，原始数组、门限和1%教师flags全部未改。
