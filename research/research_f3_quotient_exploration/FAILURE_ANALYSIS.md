# 失败、恢复与未测内容

1. 原F失败不改写。本轮新商坐标标准实际通过，旧八维Gate2仍FAIL。
2. 首attempt到C_native第二距离时q7≈1.156e−17，规范w0浮点值.39999999999999997，小于原严格支撑.4；validator在fullstate前正确拒绝。失败在decoder-only重放复现。修复只把canonical权重舍入回数学闭区间，并用2L*min(t,1−t)稳定计算density边界距离。未放宽材料、noise或概率目标。
3. 首attempt8states/224RHS与closure保留，一次恢复在同ledger继续，所有已付费state/J直接载入。没有重算旧F1.5，也没有新campaign额度。零物理诊断脚本曾有一次SCALE import错误，修正后未触发物理求解；不计为物理state。
4. B第二距离两例均不达标、A-ring第二距离不达标；原始实际RMS保留，不用更短距离冒称成功。
5. ring无competitive任务端点；weak的三端点彼此很近，0拓扑变化。D没有同成本增量。
6. 9个midpoint均通过宽noise-ball但不competitive；仅否定指定q线段，不证明无其他路径或分离模式。此检查根据已有端点追加，在求解前冻结，属于后续诊断，不称盲测。
7. 8步MH仅核验实现/发现尺度，无混合、事件概率、CDF、KL或credible coverage认证。被拒绝proposal不是样本。
8. 只有同parent的两个已暴露背景，额外已知家族/global geometry；不声称一般2D或3D成像、跨parent泛化或N3同信息量优越性。
9. 独立A17 Schur消融、普通非线性完整corrector竞赛、P1、3D、Mie、NN、严格后验reference均NOT_RUN。
