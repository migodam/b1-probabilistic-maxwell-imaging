# 有限材料传输

**NOT_RUN：F1在预分配24 states内没有达到冻结的完整anchor驻点资格。**

A Frozen Tangent、B Updated Tangent、C Material-Native Retraction、D Intrinsic Curvature + Finite Feedback 四臂均未执行；无曲线、成功率、材料位移或加速结果可以比较。已冻结的目标RMS为0.00527554628，来自原TRAIN增量RMS中位数的5%；未因失败放宽。结构检查也已固定在configs/EXECUTION_OVERLAY.json，但本轮没有使用它宣称成功。

F0完成的是两个旧已知TRAIN材料点的曲率拆分；F1完成的是同一个development parent的两个条件反演。两者不是同一anchor，也不拼成一条新传输路径。优化器12步历史是反演轨迹，不是相容材料探索样本，更不是posterior chain。

原来的Neumann/resolvent精度仍是历史证据，未重算或当作本轮探索收益。未执行任何NN训练、Mie、3D或current-space生成分支。正式传输teacher pairs数量=0。
