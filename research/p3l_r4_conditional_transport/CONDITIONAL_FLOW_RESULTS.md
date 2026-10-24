# 完整66维条件概率：模型、正面结果与资格边界

## 实际实现

Full-66 rational-quadratic Spline（8层、约0.683M）、同容量diagonal Gaussian（约0.683M）、较小Gaussian、固定16/50三角分解（约0.549M）、raw-y CFM及匹配容量CFM均真实训练。保留三seed、4k/16k/32k数据规模、learning-rate与容量检查。模型在原Gaussianized latent上有完整支持；uniform属性由probit变换正确处理。

三角法为 q(a|y)q(b|a,y)。固定TRAIN Maxwell frame 与随机正交frame享有相同数据、训练容量和概率体积。另运行真实付费anchor-conditioned none/GN/random坐标对照，满秩GN矩阵的logdet完整保留；GN预条件在当前拟合下恶化NLL，未强制作为主法。

有界几何CFM保持所有66维，以网络坐标 J_r=J_x T 构造 W=I+3Vdiag(s²/(1+s²))Vᵀ。W仅依赖y和其在线context且正定，不硬投影；随机控制保留同一谱。三seed×geometry/random实际六次拟合，固定3,000updates。原始坐标的度量是T^{-T}WT^{-1}，不能误报原Euclidean界。

## 两批独立评分

Spline NLL77.509/77.278，matchedGaussian83.077/83.449；物理三角NLL76.592/76.259，随机三角77.950/77.860。Spline相对差5.568/6.171nat，物理三角相对差1.358/1.602nat，parent区间均支持改善。

几何CFM energy score0.130661/0.138822，对随机控制0.165310/0.174586改善20.96%/20.49%。它仍远弱于好的材料风险动作；第一批模型样本均值MSE约0.03568。显式Spline样本均值约0.02278，不能用最接近GT的样本替代均值或在线决策。

这些proper scores可提供真实相对conditional-distribution证据，不能被“没有MCMC reference”抹掉。它们又不能给出绝对π距离或概率认证。[标准评分解释](AUDIT/JOINT_SIMULATION_PROPER_SCORE_INTERPRETATION.md)逐项说明。

## 实际替代机制与失败

三类physics-trained mixed bidirectionalKL与成对supervised控制真实运行。physics目标改善平均response却显著损害joint NLL和覆盖，未过复制gate，因此停止扩大。所用source-adjoint8与forward8、validation Full查询及失败均记费；不是只写理论式。

CFM完整散度审核采用CPU双精度DOP853、全66维trace、两预定tolerance，895秒预算后保留部分完成结果。某对象NLL0.136nat的容差差超过0.05；离散32步与更精细积分也不一致。CFM数值logq不用于exact-MH，未完成求解为NOT_RUN。样本energy score仍是合法的模型样本评分。

真正FULL资格见[概率审计](PROBABILITY_CALIBRATION_AUDIT.md)：本轮未通过。相对评分进步与绝对失配同时成立。
