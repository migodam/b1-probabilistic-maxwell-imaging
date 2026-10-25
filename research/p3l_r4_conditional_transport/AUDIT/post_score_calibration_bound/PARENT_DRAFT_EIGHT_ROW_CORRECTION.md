# 附加诊断：不依赖MCMC参考的材料后验误差下界

这是在全部模型和两批评分封存后，从原先已冻结的PIT>0.975事件推导的post-hoc数学解释。没有训练、重选、GT初始化或新增Maxwell计算。它使用标准数据处理与二项区间，不宣称新Maxwell定理。

设真实联合law为 P(dy,dc)=p_Y(dy)π_G(dc|y)，模型law为 Q(dy,dc)=p_Y(dy)Q_G(dc|y)。c是复材料场，允许plateau原子和退化支撑；不需要像素Lebesgue密度。定义可测事件 E={GammaCDF_256(Phi(c;y))>0.975}。在原离散P1、已验证噪声/packing及精确目标下，P(E)=a=0.025；数值/物理模型失配会破坏此前提，所以结论不外推到连续PDE或真实实验设备。

记q=Q(E)。相同p_Y给出
E_Y TV(π_G,Q_G)=TV(P,Q)≥|q−a|。
二分事件的KL链式分解也给出
E_Y KL(π_G||Q_G)=KL(P||Q)≥d_Bern(a||q)。
若绝对连续性不成立，KL可能无穷大，下界仍成立。

每个全新parent只有一个事前固定样本index0，跨parent独立生成，因此k/n估计q，不选择低残差样本。在冻结模型、理想独立随机流和原joint law前提下，以单侧Clopper–Pearson下端 L=Beta^{-1}(0.05/8;k,n−k+1) 为八项固定模型×cohort比较作Bonferroni校正。若L>a，则
平均条件材料TV ≥ L−a；
平均条件材料forward-KL ≥ a log(a/L)+(1−a)log((1−a)/(1−L))，
同时置信度至少95%。这里的置信度是有限独立parent Monte Carlo频率学意义，不是Bayesian可信度、posterior校准或浮点形式认证。

原rawSpline固定seed1：L1有510/512，L2有255/256。对应L为0.982548/0.972303，给出平均条件材料TV下界0.957548/0.947303、forward-KL下界3.83065/3.38056nat。几何CFM L1/L2的TV下界为0.957548/0.955370，KL下界3.83065/3.71603nat。CFM不需要已合格log-density即可使用这个sample-law事件界。

这个界与Spline相对Gaussian改善5.568/6.171nat完全相容：相对KL变好，不代表绝对law接近；少数投影的coverage看似合适，也不能推翻如此强的response-coupling失配。TV下界是对随机Y平均的条件误差，不能升级为每一观测都至少有同一误差。

当前八项比较与阈值是原预冻结PPC面板，不由新图挑阈值；本下界本身是在查看结果后推导，必须标POST_HOC_INTERPRETATION。原生成模型与两批主要比较均不据此更改。已有R2/R4proper-score收益继续保留；只加强“完整后验未合格”的可解释证据。

