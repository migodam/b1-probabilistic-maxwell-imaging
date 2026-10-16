# 失败与限制

1. **历史缓存缺向量。** 两旧TRAIN anchor的Maxwell directional Hessian只存了norm，q/bmax未持久化。F0把cosine、三项总向量、精确intrinsic/cokernel范数记MISSING，用旧norm条件下的三角不等式区间表达可知范围。未调用LU补齐。
2. **坐标与先验不能混为一谈。** weak原decoder二阶材料能量只有0.003286%在prior法向；ring为14.7715%。数据映射后两项可以抵消，不能简单相加范数。两例真实prior余核响应均非零；它们的方向是近似弱方向，并非exact null，所以不能援引exact-fiber二阶必要条件宣布无路径。
3. **F1驻点资格未达。** 两例已合法、噪声相容，但固定24-state分配未达0.001最优性。优化器和边界处理不足的可能性不能被当成物理不存在性。
4. **权重比例边界。** 两Gaussian权重比接近原支持最小.4/.9；共同缩放gauge不能改变比例。忽略边界的GN会预测不可行的大步；零物理的有界线性子问题给出的可行下降已很小。该结果不改变原冻结Gate。
5. **样本范围有限。** F0两公开TRAIN数学点；F1同一公开development parent的两背景，已知global geometry。不是两独立物体，更不是盲测。
6. **F2未运行。** 无非平凡传输、候选多样性、posterior校准、跨分支收益或learned Explorer证据。也没有对一般材料P1做排除。
7. **软件/物理层次。** 本轮是2D scalar-TM/VIE；A17 3D vector-DDA的非线性极化率必须另外处理。没有continuum离散误差认证或一般3D结论。

原始负结果、旧账本与历史报告保持不变；所有新付费端点已保存，没有LXZ首批缺数组问题。
