# 风险目标、概率目标与几何机制

本文件区分数学性质与待检验算法。统计目标始终是原 P1；所有 FULL 指原 24×24 二维 TM/VIE 的完整离散模型。

## 材料 MSE 监督不需要 MCMC teacher

设材料在复数 Hilbert 空间中，内积为全网格平均实内积。令 μ(y)=E[χ|y]。对任意只依赖观测的估计 T(y)，条件期望的正交性给出

`R(T) = E||χ−μ(y)||² + E||T(y)−μ(y)||²`。

因此 joint simulation 下的实际材料 MSE 是合法的条件均值学习目标；需要独立 parent 的泛化评价，不要求先求出每份观测的 posterior。有限训练数据、模型容量和优化不保证得到 μ，也不保证可信区间或完整概率密度。原 P3-L 的失败 Gaussian 平均不构成此恒等式的反例。

若 q0、J 和残差都是 y 及已知 solver 的确定性函数，它们不增加新的观测信息。它们可能改进归纳偏好和有限样本/有限算力下的学习效率；独立增量必须通过同输入信息、容量与完整成本的消融检验。不能把多做一次 Jacobian 称作取得了新传感信息。

## 材料风险与弱方向不是同一排序

`D=Dg(q0)` 给出材料速度，`J` 是正确白化的 Maxwell 导数。任务尺度下可考虑 `(||Jv||²+||v_z||²)/||Dv||²_material`。分母核或极小奇异值必须单独处理，不能把 latent 饱和当真实材料自由度。该广义谱是局部候选生成的工具，不是 posterior 的新参考测度；不插入任意 Gram determinant 改变原 prior。

在局部线性 Gaussian 近似中，材料风险贡献涉及 `D C D*`，而非单独 J 的最小奇异值。真实问题的 C 未知；本轮分别比较观测几何、监督剩余材料风险与普通不分割估计，不能把一个 GN covariance 自动认作真实 C。

有限移动还包含 decoder secant 余项和 Maxwell response 余项。`Dg δq` 很小不排除有限 δq 跨越饱和界后产生大材料变化。轨迹记录实际 RMS `.02/.04`、signed progress、原 prior cost、完整 Φ 与谱旋转，不用参数距离代替材料运动。

## 概率合同与新 non-Gaussian 提案

原 q=(z64,u2)，z 标准 Gaussian，u 在单位方框均匀；`Φ=||F−y||²/σ²`，σ固定 .00119211798514519，实虚噪声方差各 σ²/2。在 `u=NormalCDF(η)` 下 x=(z,η) 的原 prior 恰为标准 Gaussian66。

本轮多中心 Student-t 非线性三角 proposal 采用 `η=cη+s tη` 与 `z=cz+L tz+B(NormalCDF(η)−NormalCDF(cη))`。其 Jacobian 为块三角，绝对 determinant 是 `det(L) prod(s)`；非线性 shear 不改变 determinant。每个完整支撑的 mixture density 可明确求值。独立 MH 使用真实 FULL likelihood 及两端 proposal ratio；数学不变性不保证有限时间 mixing。

另一条路线按八个 illumination 分块渐进同化；每个阶段保留累计 likelihood、增量权重、祖先、重采样和合法 mutation。终点必须真实到达全部八源原目标，未完成时只报告当时的 tempered/partial target。

固定 z 后对原 u 方框积分属于 **SLICE**，其 du 测度与 FULL marginal 不同。它可以独立验证 non-Gaussian kernel、normalization 与材料矩，不会使 66维参考自动合格。

## 轨迹与风险策略

零宽轨迹没有连续多维 posterior 的普通体积。几何轨迹本轮先明确作为 candidate-only 或监督风险策略的动作环境；其访问频率不作 posterior mass。上下文 bandit 直接拟合允许的 TRAIN GT 风险改善；有限时域 model-based policy 同时拟合下一步在线特征变化，规划未来风险回报。二者均在未读取测试 GT 时冻结，并须实际运行独立在线 rollout。稀疏行为覆盖和未来局部几何的近似不构成离策略价值证明。

## 先例及证据范围

- [Cui 等，2014，LIS](https://doi.org/10.1088/0266-5611/30/11/114015)：prior-relative 子空间和 Rao–Blackwell 已有成熟先例，不能宣称首次强弱分解。
- [Adler 等，Task adapted reconstruction](https://doi.org/10.1088/1361-6420/ac28ec)：监督风险与条件均值成像已有理论基础。Semantic Scholar 的 arXiv 年份为2018，期刊发表信息另按 DOI 记录，不混作同一日期。
- [Greenberg 等，2019，APT](https://arxiv.org/abs/1905.07488)：联合模拟下的归一化 density learning 与 proposal correction 有明确先例。弱坐标预测正确不等于 FULL posterior 正确。

三条记录已通过 ScholarQA 批量核验；本次发现查询部分429，保留原始返回。上述定位不是穷尽新颖性审查。独立 Maxwell 算法贡献取决于真实风险/统计精度及总成本的结果。
