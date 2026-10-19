# 理论：哪些成立，哪些不是新定理

## 材料尺度与prior测度

令M=[Re Dg;Im Dg]/24，M=Q S V^T，T=V S^{-1}，则MT=Q。K=J_white T。加入T^T H_prior T描述原Gaussian局部代价；它不是新的材料volume density。Dg的数值核、uniform property边界和可行切锥仍须处理。

可逆线性坐标变换q=C a时，M、J、H分别推送为MC、JC、C^T H C；对相同物理材料速度，Q W Q^T不依赖该参数标签。实测一处缓存相对差5.27e-14。非线性重参数化须另处理坐标加速度，这里没有宣称一般hessian tensor不变。

## soft滤波的连续性

对PSD H，W(H)=(I+s²H)^{-1}。精确恒等式W(B)-W(A)=-s² W(B)(B-A)W(A)，故norm差≤s² norm(B-A)。它解释soft filter避免hard截断在谱交叉处跳变；这是标准resolvent分析，不是Maxwell新定理，也不说明一定更便宜。

## 有限数据竞争的条件

实际材料增量h=g(q+d)-g(q)，完整白化变化为J_chi h+e。令a=J_chi h，则

ΔPhi = r^T a + ||a||²/2 + (r+a)^T e + ||e||²/2。

若有严格||e||≤epsilon，线性预测的Phi误差不超过||r+a||epsilon+epsilon²/2。这里只测secant与full端点，没有uniform bound，因此刷新是经验规则，不能标certified early acceptance或路径证书。decoder余项h-Dg d应先与Maxwell余项分离。

## 先验边界与精确fiber

原P1局部数值满秩不提供非零exact-fiber tangent；finite-noise不确定性无需exact null。合法端点不推出整个路径合法/相容，GT差量弱能量不推出可达性，corrector低residual不推出新结构或posterior质量。

## 与已有研究的关系

[Cui et al., 2014, LIS](https://arxiv.org/abs/1403.4680)已用prior与likelihood相对影响构造低维推断。这里增加材料任务metric，不把谱分解本身写作原创。

[Cui, Law, Marzouk, 2016, DILI](https://arxiv.org/abs/1411.3688)已研究operator-weighted、likelihood-informed MCMC；本轮candidate corrector不是DILI或合法posterior transport的替代。

[Girolami and Calderhead, manifold MCMC](https://doi.org/10.1111/j.1467-9868.2010.00765.x)已有位置几何采样。本轮没有证明新的不变测度结果。

A17 Schur/current-space恒等式保留。相同空间精度相同；新增价值只能来自选择、复用或误差控制。本轮该费用收益为负。
