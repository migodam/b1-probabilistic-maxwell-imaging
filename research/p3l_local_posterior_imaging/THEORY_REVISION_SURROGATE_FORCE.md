# 从 Gaussian 材料平均失败到保留 decoder 的代理力

这是本轮实际负结果引出的新机制试验；不是已完成的新算法成绩。固定 P1、noise、prior 和已登记的独立风险评价均保持不变。

## 失败首先发生在哪里

解析材料矩准确积分了声明的 Gaussian 参考，但该参考没有通过真实 likelihood 校正。加入两个属性的不确定性及其与空间场的相关性后，32 个 development 对象的平均材料 MSE 仍从新坐标 mode 的 0.0120131 升到 0.0185700，四个 importance 诊断 ESS 均约 1。数值 quadrature 已稳定，不能再仅靠增加 Gaussian 积分精度修复此项偏差。

饱和 decoder 的一个具体风险是：当前导数为零或很小，不代表有限 prior 扰动不能越过饱和边界。只用 latent Jacobian 的局部二次模型，可能对这些有限材料变化赋予错误厚度。此解释须接受实际 full-response 对照，不能从局部导数直接推定完整 posterior。

## 两种相同付费权利的代理

在已付费的材料点上，保存 receiver transfer 与全部照明的总场。由此形成材料响应矩阵 C，满足该点的材料一阶响应为 C·delta_chi。普通 latent 代理冻结白化 Jx，使用 f0+Jx·delta_x。

材料代理使用 Fhat(x)=F0+C·(g(Tx)-chi0)。它完整保留原 decoder 与属性 probit 变换，只冻结 Maxwell 材料响应。它没有计算精确有限 feedback，也不声称有限响应误差已认证。两个代理使用相同的已付费背景、Full-J setup 和缓存权利。

材料代理的 log-density 力为 `-2 Re[Dg^H C^H(Fhat-y)]/sigma²`，再乘正确的 probit 链式导数并加 `-x`。它只做矩阵作用和 decoder 计算，不在内部轨迹调用 Full solver。

## 为什么允许使用近似力

固定 SPD 的 inverse mass C_m，动量服从 N(0,C_m^-1)。每个固定力 kick 和位置 drift 都是行列式为 1 的剪切映射；对固定步数的对称 leapfrog，加动量翻转得到可逆的体积保持提案。端点用真实 `U_full=Phi(Tx)+||x||²/2` 和真实动能做 Hamiltonian MH，而不是用代理能量代替。

因此，固定代理力可以不同于真实梯度，接受修正仍针对声明的完整数值 posterior。这里要求同一次正反轨迹使用同一代理、质量矩阵和步数；生产阶段不能随起点重新建立代理后仍套用这个简单接受率。新试验只在丢弃的 warmup 中刷新，生产期冻结。非有限轨迹、舍入问题和混合失败都必须单独报告。

对于确实来自代理势能的力，端点能量误差可分成两项：代理 integrator 误差，以及 Full–proxy 势能差在两个端点的变化。接受率低时，这提供了区分积分不准与 Maxwell 代理漂移的诊断，不能仅凭接受率评价混合。

## 已有工作与本轮可检验增量

近似力加真实 Hamiltonian 接受修正已有成熟先例，不能写成本轮原创定理。Li、Holbrook、Shahbaba、Baldi 的 [Neural Network Gradient Hamiltonian Monte Carlo](https://arxiv.org/abs/1711.05307) 已研究此类正确性；作者提供的 [正式论文](https://andrewjholbrook.github.io/papers/Li.pdf) 给出了相应机制。本轮使用解析的 Maxwell 状态响应，没有训练该文的梯度网络。

本轮要检验的是更窄的问题：在相同新 Full 端点预算下，保留非线性材料 decoder 的冻结物理响应，能否比冻结 latent 一阶代理和普通 Full-gradient HMC 提供更多材料函数 ESS、更小 MCSE，并在总 wall 计入廉价力的全部矩阵作用后仍有收益。仅仅减少 RHS、不改善统计精度，或者只有接受率变高，都不能构成成功。

新试验固定为已暴露 DEV_030 的机制对照。即使该对象出现正结果，也不代表已完成独立泛化、可信 Bayes 成像或强 IEEE TAP 贡献。128 个新对象的已冻结估计器评价不因该结果事后改变。
