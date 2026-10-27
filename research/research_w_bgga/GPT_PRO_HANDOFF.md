# 给下一轮GPT Pro：哪些机制值得保留，哪些应删除

先读RESEARCH_W_BGGA_REPORT_ZH、GATE_DECISION及COST_SUMMARY，不把准备页当最终结果。

1. 保留：公共target ROI固定P3离散；材料依赖b_X的正确增量式；完整墙–目标反馈；full/Schur/resolvent一致。12案例独立贴体白化差≤.02033。
2. 保留问题：材料不确定性与指定传播任务不确定性不等价。C2高相对保留仍98%prior方差；C4_1单端口d=.729、future25 max3.068且自己的argmax经独立复核，但主观测d约2.77，目标不移动。
3. 删除普遍单点安全主张：最佳MAP-wall仍11/24 jointTV失败；目标行动可好而joint/cross-covariance错误。TV≥1−P(W=g)只约束Dirac墙替代，不能用它否定精确Green重参数化。
4. 删除当前full-port加速：2257ports、1.026GB images、setup2257RHS，query .572–.871CPU秒，direct .340秒。普通ROM .011秒但held-out误差至6.517。当前无同质量突破。
5. 默认[0,1]区间策略没有省查询：24case全查2950。current/future complete-link没有合并，始终118簇。一次高精度归一化只修浮点包络，不能替代物理误差界。
6. G1为旧bank重放，不是新的冷部署。G2/未知支撑/新确认未跑，普通强smoother/FWI没有被击败。旧A0未达仍未达。

最有价值的新问题：能否在保持原joint与随机残差材料的情况下，仅计算当前任务确实需要的Green actions，并得到比直接full和ordinary goal-ROM更便宜的likelihood/task包络？应先给可核对的小成本模型；若必须先查全表或全后验，判经济性NO-GO。

不建议增加NN、Flow、A17模块或先跑30parent确认。未证明新物理信息、posterior calibration或TAP级算法原创。
