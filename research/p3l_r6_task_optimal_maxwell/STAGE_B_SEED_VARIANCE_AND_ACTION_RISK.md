# Stage B：有限seed分歧与action风险

固定三个模型之间的分歧，可以精确分解已有预测的平方风险，但不能据此估计总体训练方差。历史N32768学生的mean-seed风险为.00716774，三seed平均图像风险为.00693198，有限预测分歧为.000235754，恒等式残差约6e-10。所有评价均用已有VAL256，三模型不是从未知训练分布获得的大样本方差估计。[有限模型结果](track_b/SEED_COUPLING.json)逐项保存mean-seed与ensemble两个口径，没有互换分母。

进一步交换provider_i/head_j形成完整3×3表。原基准的九组合预测分歧可正交分为provider约4.01195e-5、head约1.91752e-4、interaction约9.97294e-6，九组合平均风险.00716842、grand prediction风险.00692658。这是平衡有限表的代数分解，不是population bias/variance。原基准off-diagonal减diagonal风险约1.025e-6，而direct_gt约2.460e-4、GT+action约3.437e-4，显示这些已训练组合存在配对依赖；不能把跨seed交换的损失当成独立复现实验。

固定teacher action a_T、学生动作偏移d时，有限评价集上风险差可写成2平均〈a_T−GT,d〉+平均||d||²。实际[有限task风险恒等式](track_b/FINITE_TASK_RISK_IDENTITIES.json)残差小于约1.43e-9；原学生三个seed的cross term均为负。因此“学生特征误差更小”与“最终材料风险更低”不是同一句话。但GT展开没有提供真实条件均值或Bayes风险，不能由cross term直接声称teacher总体bias已被测准。

[主理论文件](THEORY_AND_COUNTEREXAMPLES.md)的teacher-bias说明保持原意：固定正teacher-action权重，在无限容量平方损失理想化下会混入teacher action偏差；这不是对当前有限神经优化器的数值保证。本轮GT+action的matched风险.00790106，不能用这个单点倒推所有teacher混合目标都无效。新四臂只有1024个监督action上下文，历史基准32768个provider对象；有限容量、选择、训练样本与表示误差均不能混成一个“纯bias”结论。

所有九组合预测、逐parent风险、有限分解与通道控制保留在[seed文件](track_b/SEED_COUPLING.json)、[diagnostics目录](track_b/diagnostics/)和[主B报告](track_b/REPORT_ZH.md)。这里没有新训练或物理调用，没有读取新独立truth，也不把历史VAL2048称作新locked；最终总体解释留给主线程。
