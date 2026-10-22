# FULL posterior 参考资格审计：进行中的检查点

本文件依据已完成的 REFERENCE_STATS.json 摘录；不读取链大数组、GT 或追加账本，不重新计算物理。状态为 **RUNNING / FULL_REFERENCE_UNQUALIFIED**，最终科学审查由主线程负责。

## 为什么已经运行的链仍不能当 teacher

四条链各自看起来稳定，或者有较高接受率，都不足以说明它们覆盖了同一目标。当前冻结门槛要求全部相关 Rhat<1.01、材料函数 bulk 和 tail ESS 均≥1000、完整网格材料均值 RMS MCSE≤0.001。还要检查 ROI 未评估常数分量，并在同机制门槛后做独立机制核验。下面保留原诊断中的缺失值；不能把无法评估的常数轨迹填作无限 ESS。

| 已完成对象 / 机制 | 每链 production | 材料函数最大 Rhat | 最小 bulk ESS | 最小 tail ESS | 材料均值 RMS MCSE | 当前资格 |
| --- | --- | --- | --- | --- | --- | --- |
| [DEV_000/mala](A/references/DEV_000/mala/REFERENCE_STATS.json) | 4096 | 3.20089 | 4.45452 | 1 | 0.0556548 | UNQUALIFIED |
| [DEV_028/mala](A/references/DEV_028/mala/REFERENCE_STATS.json) | 4096 | 5.39208 | 4.17595 | UNASSESSED | 0.105155 | UNQUALIFIED |
| [DEV_030/hmc_pilot](A/references/DEV_030/hmc_pilot/REFERENCE_STATS.json) | 1024 | 1.26968 | 11.4748 | 23.4906 | 0.0101516 | UNQUALIFIED |
| [DEV_030/mala_mode](A/references/DEV_030/mala_mode/REFERENCE_STATS.json) | 4096 | 1.92985 | 5.55516 | 8.28549 | 0.0182821 | UNQUALIFIED |
| [DEV_030/mala_retracted](A/references/DEV_030/mala_retracted/REFERENCE_STATS.json) | 4096 | 2.29 | 5.00098 | 6.39396 | 0.0239059 | UNQUALIFIED |
| [DEV_030/pcn_mode](A/references/DEV_030/pcn_mode/REFERENCE_STATS.json) | 4096 | UNASSESSED | UNASSESSED | UNASSESSED | 0.010337 | UNQUALIFIED |
| [DEV_030/mala_decoder](C/references/DEV_030/mala_decoder/REFERENCE_STATS.json) | 4096 | 4.18571 | 4.26597 | 4.97157 | 0.0230438 | UNQUALIFIED |

HMC pilot 有四条链，每链剔除 512 warmup、保留 1024 production；其余上表实验每链剔除 1024 warmup、保留 4096 production。HMC 的材料 Rhat 约 1.270、bulk ESS 约 11.47、材料均值 RMS MCSE 约 0.01015，仍未达到冻结门槛。它是独立机制的负诊断证据，不能因为优于某些 MALA 指标而升级为合格 reference。

已运行的 tempered SMC 也不能仅靠 beta=1 或权重 ESS 获得全局参考资格。当前链与粒子结果只支持诊断，不用于 posterior distillation、可信完整后验区间或已验证 posterior mean 风险声明。

## 不能合并成同一个目标的证据

[GN departure 审计](B/GN_DEPARTURE_AUDIT.json) 的 40 个强条件参考中，38 个重要性 ESS<2，最大约 5.269。10 组的局部线性梯度 logZ 修正上界小于 0.01，仍不能修复这些参考在整个抽样范围内的非二次误差。Gaussian 解析材料矩准确表达的是给定 Gaussian reference 的矩。

[三维切片](B/slice3/RESULTS.json) 固定其他 61 个 Gaussian 坐标与两个属性，并条件于预先冻结的有限 box。最后两阶材料均值 RMS 差为 8.45585e-10、logZ 差为 6.64293e-8；这只是该条件切片的有限阶数收敛证据，不是原 FULL marginal、box 外质量或经过认证的 quadrature 误差。

[材料风险](B/GAUSSIAN_IMAGE_RISK.json) 保留全部 32 个 development parent。GN Gaussian 材料均值相对优化 anchor 的平均 complex MSE 差为 +0.00687494，配对 bootstrap 95% 区间 [0.00296121, 0.0120149]。这是具体近似的实际恶化，不能推出正确 FULL posterior mean 无益。GT 只在此前独立离线评分中读取，本审计只读取已保存指标；数据已暴露，不是 locked 泛化。

## 允许怎样解释当前结果

物理单点基础图、SLICE 数值参考、Gaussian 参考矩、joint-simulation supervision 和 FULL posterior reference 分别保留其含义。归一化网络和合法 MH 公式都不能替代实际材料函数混合诊断。后续方法的最终比较还需真实 full likelihood 校正、风险、setup/refresh/fallback 与完整费用。尚无可信 FULL posterior teacher 或最终新颖性结论。

[概率合同](THEORY_AND_ALGORITHM.md) · [中文记录](RESEARCH_REPORT_ZH.md) · [资源快照](RESOURCE_LEDGER.json) · [内存恢复](AUDIT/MEMORY_RECOVERY_REVIEW.md)
