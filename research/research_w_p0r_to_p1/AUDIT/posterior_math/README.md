# 概率代数合同：仅有限 quadrature

`src/wr/posterior_math.py` 只接受调用者提供的数组，不读真值、实验结果或数据文件，不生成物理观测，不调用 PDE、Book、RHS 或 LU。当前 P1 gate 为 **NOT_EVALUATED**。`CHECK_RESULTS.json` 是确定性小型合成数组的代数检查，不是 Maxwell、posterior 或 P1 结果。

## 共享 nuisance 与完整绝对似然

默认数组末三轴为 `(frequency,source,receiver)=(5,4,12)`，所以复观测数 `M=240`，实正交分量数为480。每个实分量的标准差固定为 `σ=0.01`；不能将它当成每个复样本模长的标准差。f、y 的前导轴可以广播，支持单个 y 或一批 y。

全频率、全部源和接收器共享同一个 `γ=a exp(iφ)`；`a∼U[0.99,1.01]`、`φ∼U[−0.01,0.01]` 独立。没有按频率或按源自由拟合，也没有 profile likelihood。

定义 `F=Σ|f|²`、`C=Σ f conj(y)`、`Y=Σ|y|²`，则

`Φ(a,φ)=[a²F−2a Re(exp(iφ)C)+Y]/(2σ²)`，

`log p(y|f,a,φ)=−2M log(σ√(2π))−Φ(a,φ)`。

完整高斯常数在默认形状下为约 `1769.391193336041`。不除以通道数或其平方根，不能把这个常数丢掉再报告绝对 evidence。

`log_marginal_likelihood(f,y,n=11|21)` 对 nuisance 使用二维 Gauss–Legendre 积分。两个均匀先验的映射 Jacobian 与密度抵消，最终节点权重为 `w_i w_j/4`，总和为1；不能额外乘 rectangle 面积，也不能积分后挑最大节点。返回：

- `log_normalizing_constant`：上述完整高斯常数；`logZ`：包含该常数的绝对边缘对数密度。
- `posterior_nuisance_moments`：幅度、相位和复 γ 的均值、方差与协方差，均按后验 quadrature 节点权重计算。
- 可选 `posterior_node_weights`、节点与积分权重；以及浮点消去诊断、实际累加类型和机器精度。

`source_blocks=True` 分成四个源块，也可显式给出源索引分区。实现先对块累加 F、C、Y，然后只构造一次完整 likelihood、积分一次共享 γ。重复或遗漏源会拒绝。逐块分别 marginalize 再相乘会改变模型，因为它隐含四个独立 nuisance；这不是合法加速。

为稳定计算，先减去最小节点残差再计算 logsumexp，后验节点权重另行归一化。统计量消去仍可能造成精度损失，`longdouble` 在本机是 float64 等效精度；机器精度尺度不是认证误差上界。代码仅将64eps尺度内的微小负残差截成0，更大的负值报错；截断计数会返回，不隐藏。

## 有限参数网格后验

`finite_grid_posterior(log_likelihood,prior_density,cell_weights)` 使用

`q_i ∝ exp(log_likelihood_i) × prior_density_i × cell_weights_i`。

它不会把网格节点当作等概率样本。原 quadrature 的先验总质量可能不等于1，因此同时暴露 `prior_quadrature_mass`、原始 quadrature evidence 和归一化后的有限先验 evidence；不偷偷重标绝对 evidence。

`finite_grid_summary` 计算分支权重、布尔事件概率以及目标位置均值、方差与协方差。`weighted_cdf` 返回有限节点原子分布的阶梯 CDF，明确不是连续 CDF，也不是 node-bin 连续后验概率。`weighted_crps` 是该有限原子分布的 CRPS；`brier_score` 对输入的事件概率和0/1结果评分。它们都没有把有限 quadrature 自动升级成连续后验结论。

## 检查与保留的误差

`check_math.py` 独立逐节点计算复残差，并使用独立的 logsumexp 实现核对统计量公式；没有调用生产函数来构造参考 likelihood。

- 高斯常数、零场时 nuisance 的均匀先验矩、很大的共同残差下的节点权重均通过代数核验。
- 普通小数组的绝对 logZ 与直接复残差相差约 `2.30e-12`；11阶与21阶相差约 `6.25e-13`。
- 四源块与完整计算相差约 `2.22e-12`；错误的逐源独立 nuisance 积分改变 logZ 约 `0.08199`，这一负对照保留。
- 批量 f/y 广播输出形状 `(2,3)`，与独立残差参考相差约 `3.41e-12`。
- 非均匀 prior density、非均匀体积、branch/event weights、位置矩、CDF、CRPS/Brier 与独立显式求和一致；极小似然使用 logsumexp 避免下溢。

**11/21阶不是通用收敛保证。** 刻意构造的尖锐小数组 likelihood 给出 `logZ11≈154.5658591`、`logZ21≈170.6755696`，差约 **16.1097105**；各自仍与同阶直接复残差求和相容（最大差约 `2.45e-9`）。这说明代数实现正确也不能替代 quadrature 收敛检查。该项记录为 `DISCREPANCY_RETAINED`，没有扩大物理预算、提高噪声、profile nuisance 或宣布 P1 合格。

运行已有环境中的纯代数检查：

```sh
PYTHONPATH=Gaussian/B1/b1-w-p0r-to-p1/research/research_w_p0r_to_p1/src Gaussian/.venv_nn/bin/python -B Gaussian/B1/b1-w-p0r-to-p1/research/research_w_p0r_to_p1/AUDIT/posterior_math/check_math.py
```

检查只写本目录 `CHECK_RESULTS.json`。没有采样、真实数据生成或真值访问；合成数组由固定公式构造。物理与 P1 资格仍须主线程另行验证。
