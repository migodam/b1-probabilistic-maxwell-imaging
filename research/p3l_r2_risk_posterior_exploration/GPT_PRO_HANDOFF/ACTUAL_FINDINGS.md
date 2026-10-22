# 可用于理论推理的实际现象

1. 首批256fresh，固定局部CNN风险降低25.27%，调整CI[22.34%, 27.99%]；不是后验reference结论。第二批256fresh的21模型拆分控制结果单独保留，未换第一批winner。
2. P1灵敏度图精确因式分解为decoder导数×固定空间包络。完整灵敏度＋回投影 相对 仅残差回投影图在这批确认对象上风险较低。 完整灵敏度＋回投影 相对 仅训练谱灵敏度图的调整后区间跨零，独立增量没有建立。 仅训练谱灵敏度图 相对 仅随机谱 decoder 灵敏度图在这批确认对象上风险较低。 这要求区分Maxwell frame、decoder几何与卷积空间归纳偏好。
3. 约99.10%风险收益位于基础点过渡相；修正13.03%能量在prior线性切空间外，但仅约5.83%收益来自该部分。不能简单说NN只在nullspace填充或曲率解释全部风险。
4. 66维Student、sourceSMC、动态MALA、learned10×prior56都没得到合格FULL reference；projection NLL有改进但全维校正仍失败。新Proposal的真实ΔΦ与latent线性预测差数千倍，decoder/Maxwell余项存在抵消。
5. 相同field2 SLICE的直接CDF积分与GH节点CDF差约0.28，已有pCN与直接CDF差约0.016；不能将坏quadrature误差全归于sampler。Jensen gap约7.49e-4，严格logZ refinement仍未过。
6. rank4动态弱分布在两个真实P1 anchor上表现为非对合，normal leakage约0.93/0.82，差分和反向回路有数值对照；不是toy反例，也不是全P1定理。
7. 1934/2048个所选轨迹端点在固定任务尺度和数据判据下合格，但冻结弧长平均风险未改善。数据可行移动、材料风险动作与正确后验积分是不同对象。

所有数学身份、真实Maxwell结果、探索性离线归因和未执行猜想分开；本包不提供FULL posterior probabilities，不提供原P1之外的mass。
