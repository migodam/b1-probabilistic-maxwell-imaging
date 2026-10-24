# 由真实响应诊断触发的条件密度修订

本记录先于新增物理训练及新 locked 数据生成。原始 supervised Full-66 模型的负结果保留。

## 诊断与待证伪假设

对16个固定 development 对象、每模型64份预定样本，完整 Maxwell 响应存在两种误差：响应均值偏离观测，以及样本响应离散过大。恒等分解 `mean Phi = Phi(mean response) + response dispersion` 已用已付费数组逐项核验；不能仅把整个问题称为“方差过宽”。另一个精确望远镜分解显示 decoder secant 与 latent 一阶项具有强抵消，不能用各项范数百分比声称因果归属。

新的可证伪假设是：在保持完整66维可归一化密度及 joint supervision 的同时，用真实 Maxwell likelihood 对当前生成样本训练，可以减少上述误差；其代价也可能是覆盖收缩及 proper NLL 恶化。因此同时保留纯监督继续训练和完整概率诊断。

## 数学目标

令真实联合生成 law 为 `p(y) pi(x|y)`，显式模型为 `q_theta(x|y)`。训练目标为

`E_joint[-log q_theta(X|Y)] + lambda E_Y E_q[Phi(X,Y) + ||X||^2/2 + log q_theta(X|Y)]`。

在相应积分存在时，除去不依赖 theta 的常数，这是正向 KL 与反向 KL 的加权和。它没有改变原 P1 posterior 的定义；有限容量、有限训练和随训练改变的 lambda 都不赋予模型精确 posterior 身份。反向项通过真实 Full forward 与 source-adjoint VJP 求导；样本始终由完整可逆66维模型产生，不使用确定性 GN corrector 后冒用原密度。

两个训练臂从同一 checkpoint 开始，使用相同 inclusive minibatch 顺序及固定验证噪声。纯监督臂不填充无用 RHS，仅支付相同验证查询；比较真实成本而非强行耗尽同一预算。选择标准和后续复制 gate 必须在首个 fit 前冻结。

## 先例与创新性边界

- [Karumuri 与 Bilionis，2024](https://arxiv.org/abs/2305.20004) 已研究 Gaussian 与 conditional Flow 的 amortized variational inverse inference。本轮使用反向 ELBO 不是新原理。
- [Siahkoohi 等，2021](https://arxiv.org/abs/2101.03709) 已提出以预训练 conditional normalizing flow 初始化高保真 KL 训练。预训练后物理微调不能单独作为原创贡献。
- [LazyDINO，JMLR 2026](https://jmlr.org/papers/v27/25-0858.html) 使用离线导数信息代理与在线 lazy transport variational inference。R4 的局部 value/JVP 代理只是对照原型，并不是完整 LazyDINO 复现；不能以该原型失败宣称超过 LazyDINO。

以上三个主来源的官方摘要和元数据已在本轮核对。摘要范围没有核验的具体 theorem 或 benchmark 不作引用主张。是否存在 Maxwell 特有增量仍需独立实验和最终科学审查。

## 数值桥接修正

首次实际64 RHS 预检的跨设备梯度比较 FAIL 原样保留。它比较了 Float64 anchor 与 Float32 舍入后的不同材料输入。对已持久化的同一 Float32 输入，其完整 Maxwell 真梯度与 MPS 回传相对差为 `2.4005685983571236e-8`。两项原 CPU 方向有限差分均满足原阈值。修正使用已付费 query6 原始数组，新增 states/RHS 均为0；没有放宽阈值。详见相邻两个预检 JSON。

状态：实现与物理梯度桥接合格；新条件密度的训练效果及 posterior 资格尚未建立。
