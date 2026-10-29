# 固定子空间的证据与合同

R7 的八个对象在本轮全部视为已暴露开发对象。新对象只属于 `R9V_FRESH_000`—`031`，其预测与GT评分另行封存。历史R7结论、门和账本没有改写。

## 必须区分三种固定

1. 接收算子 $S$ 与背景内部Green $K$ 固定。其SVD可以一次分解。
2. 材料/latent分析坐标 $Q$ 固定。系数仍可变化；decoder与完整forward仍然非线性。
3. 局部材料敏感矩阵 $J(x)$ 固定。它一般随材料、内部总场、decoder导数改变，不能由前两项推出。

旧Frozen-J corrector属于第三种；本轮fixed-subspace + dynamic Full-J属于第二种。报告与代码保留不同名称。

## 实际稳定性

八条R7已保存轨迹在RMS约0.02处的弱空间最大主角度分别为9.79、6.34、6.53、52.75、20.18、5.83、11.24、6.88度。有效材料rank均为66。谱簇projector、材料metric、谱隙与单根奇异向量的符号/交换分别处理。跨parent汇总基只作迁移控制，不冒称同轨迹稳定性。

角度大不自动证明固定方案失败；角度小也不保证有限材料半径上的corrector可行。最终以实际材料、完整Maxwell、任务与费用为准。

## 已完成开发比较

原E2中固定弱幅度同时固定材料半径，可能把strong correction空间压缩为零。其14/16个B/C相同端点和失败原样保留，不作固定探索不可能的证明。E3_REV1允许弱射线幅度随decoder径向retraction调整，仍保持固定探索方向；B/C在两个半径均8/8合格，修复说明约束几何而非公式误差是原障碍。

修复后固定/动态的费用优势不稳定，未达到20%门。32个新对象的E5因此研究新的、完整的候选云输出，不把旧开发结果称独立确认。

## 数学与源码入口

- [材料几何与corrector推导](formula_collection/03_MATERIAL_GEOMETRY_AND_CORRECTION.md)
- [候选体积推导与实际引擎定义](formula_collection/07_MATERIAL_POSSIBILITY_VOLUME.md)
- [冻结候选协议](config/possibility_volume_protocol.json)
- [实际稳定性表](STABILITY_METRICS.csv)

E5固定/动态corrector均限制为strong complement加一条weak line；随机控制使用完整latent切向，因此是较强的普通corrector类。三臂比较是整个算法策略比较，不声称只改变一个提案方向的严格因果隔离。
