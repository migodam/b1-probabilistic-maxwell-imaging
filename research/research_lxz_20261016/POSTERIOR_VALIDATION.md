# 概率正确性与有限链误差

有限状态详细平衡核验：selector修正不变性误差2.78e-17，不修正偏差0.042；global正确qratio误差5.55e-17，漏qratio偏差0.06。高ESS漏一半mass反例保留。未重跑旧合成检查。

独立审阅确认当前kernel数学target未发现新偏差。几何FD索引修复前仍是冻结prior-preserving kernel，但不能解释成Jacobian；旧两臂131072次全保留，再修复重跑131072次。采样初始点为64个prior中最佳likelihood，SMC从全部64点开始，早期差异包含初始化差异。

| 方法 | event MAE | 材料CDF grid误差 | recall* | 4096次平均arm秒 |
|---|---:|---:|---:|---:|
| local | 0.012039 | 0.071142 | 1.0000 | 0.4998 |
| geometry | 0.026158 | 0.178926 | 0.9479 | 0.5615 |
| geometry_nuisance | 0.011720 | 0.075228 | 1.0000 | 0.5531 |
| global | 0.004364 | 0.028285 | 1.0000 | 0.5079 |
| hybrid | 0.007243 | 0.042557 | 1.0000 | 0.5231 |
| tempered | 0.016145 | 0.103910 | 0.9896 | 0.0127 |
| prior | 0.011206 | 0.074878 | 1.0000 | 0.4307 |

*固定参考事件质量≥.05、经验估计≥.005；这是本轮false-exclusion诊断阈值，不同于旧P2的.01 discovery。CDF仅保存grid最大差。SMC为64个加权粒子，MCMC为保留拒绝重复的链平均，ESS不可直接互换。


事件、mean、SD、CDF、quantiles均从同一chain/weighted empirical measure求取，条件参考来自原同一continuous target积分，见POSTERIOR_FUNCTIONALS.json。reference数值误差另保留；不将empirical CDF叫continuous supnorm。constant-event ESS为null，不能删掉后声称全mode混合。peak-probability excess只是signed peak差，不等同校准或错误自信发生率。

teacher/bank是旧物理+后验额外监督，C1不是免费或同信息量优势。有限4seed链不是iid exactposterior，SMC finiteN也可能漏区域。不声明posterior calibration、fullMie正确性、continuum或模型外coverage。

实际proxy事件包在MaterialHypothesisPackage/PROXY_*.json；展示mass不重归一化，candidate点概率null。代表图的22次额外response不叫full验证，已计入独立ledger。
