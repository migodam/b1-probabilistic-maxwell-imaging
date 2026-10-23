# 概率正确性、数值准确度与校准边界

这里可以支持“在明示受控家族内，分支与目标CDF具有经过检查的数值准确度”，不能支持“原任意未知墙体posterior已经校准”。有限元响应精度、求积误差、Monte Carlo误差、统计coverage和模型外support是不同对象。

## 已经核验

1. 原proper prior的分支质量与连续volume均保留；full quadrature weights之和1。
2. gamma在完整likelihood内共享并边缘化一次；分源同化最后各likelihood因子只使用一次。
3. Full5/9不同物理节点集与直接截断CDF交叉核对，32观测全部通过0.02事件/CDF门槛。
4. 18-node response插值与强直接三阶求积各32/32通过同一门槛。
5. SMC独立审计核验了有限目标的prior-reversible全局提案、对称局部索引提案及正确MH比；拒绝重复状态、重采样、祖先信息全部保留。

独立实现审计发现并修复了零权重粒子槽掩盖漏分支的诊断错误，以及允许空完成标记通过的离线GT保护错误。修复前反例、修复后纯数组核验均保存；P1-B推断在修复后执行。它们不是新增Maxwell成果。

## 不能推出的结论

5个CDF检查点的数值收敛，不是完整CDF的uniform认证，更不自动认证任意95%区间coverage。P1-B因此单列有限原子测度的equal-tail区间，保持无bin展宽，并在离线评分中报告经验coverage；只16个parent，统计能力有限。

两个branch在新prior中是明确事件，未宣称它们是原连续高维材料posterior的两个自然模式。Outside-family mass为null，不能从两个事件在受控模型内质量和1推出原未知世界覆盖100%。

posterior证据是观测密度，绝对logZ可为正；它不是“观测事件概率”。分支odds由proper prior和完整边缘likelihood得出，未使用residual softmax。概率误差和模型外未知量不会被候选图的美观替代。
