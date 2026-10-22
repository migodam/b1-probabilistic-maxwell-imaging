# 完整 posterior 参考：最终资格审计

**11组四链，0组合格。** 原门槛保持 Rhat<1.01、材料函数 bulk/tail ESS≥1000、完整网格材料均值RMS MCSE≤.001，并需独立机制核验。以下 UNASSESSED 不替换为无限 ESS。

| 对象 / 方法 | 每链 production | 材料函数最大 Rhat | 最小 bulk ESS | 均值 RMS MCSE | 资格 |
| --- | ---: | ---: | ---: | ---: | --- |
| [DEV_000/mala](A/references/DEV_000/mala/REFERENCE_STATS.json) | 4096 | 3.20089 | 4.45452 | 0.0556548 | UNQUALIFIED |
| [DEV_028/mala](A/references/DEV_028/mala/REFERENCE_STATS.json) | 4096 | 5.39208 | 4.17595 | 0.105155 | UNQUALIFIED |
| [DEV_030/ess_mode](A/references/DEV_030/ess_mode/REFERENCE_STATS.json) | 4096 | UNASSESSED | UNASSESSED | 0.0114013 | UNQUALIFIED |
| [DEV_030/hmc_pilot](A/references/DEV_030/hmc_pilot/REFERENCE_STATS.json) | 1024 | 1.26968 | 11.4748 | 0.0101516 | UNQUALIFIED |
| [DEV_030/mala_mode](A/references/DEV_030/mala_mode/REFERENCE_STATS.json) | 4096 | 1.92985 | 5.55516 | 0.0182821 | UNQUALIFIED |
| [DEV_030/mala_retracted](A/references/DEV_030/mala_retracted/REFERENCE_STATS.json) | 4096 | 2.29 | 5.00098 | 0.0239059 | UNQUALIFIED |
| [DEV_030/pcn_mode](A/references/DEV_030/pcn_mode/REFERENCE_STATS.json) | 4096 | UNASSESSED | UNASSESSED | 0.010337 | UNQUALIFIED |
| [DEV_030/spine_hmc_pilot](B/references/DEV_030/spine_hmc_pilot/REFERENCE_STATS.json) | 1024 | 1.6054 | 6.85148 | 0.0108156 | UNQUALIFIED |
| [DEV_030/mala_decoder](C/references/DEV_030/mala_decoder/REFERENCE_STATS.json) | 4096 | 4.18571 | 4.26597 | 0.0230438 | UNQUALIFIED |
| [DEV_030/surrogate_force_linear](C/references/DEV_030/surrogate_force_linear/REFERENCE_STATS.json) | 4096 | UNASSESSED | UNASSESSED | 0.0166675 | UNQUALIFIED |
| [DEV_030/surrogate_force_material](C/references/DEV_030/surrogate_force_material/REFERENCE_STATS.json) | 4096 | 1.63156 | 6.66187 | 0.016892 | UNQUALIFIED |

计划选择16个深参考条件，实际深入3个独立对象，并集中运行不同原理的失败诊断。其余13个没有完成深参考，不能计入已运行。44条链共保存155,648个production状态；它们不是155,648个有效独立样本。

MALA、pCN、Elliptical Slice、普通HMC、非线性spine-HMC和两种代理力的生产目标均为原FULL数值posterior。生产期参数冻结，拒绝保留重复状态，非有限事件单独记录。合法不变核不保证有限时间混合。

独立128粒子SMC到beta=1后仅剩一个祖先，不能做teacher。128个新对象的Full-SNIS均ESS<2，不能以精确权重掩盖有限粒子退化。未进行posterior distillation、材料credible coverage或完整posterior质量声明。

三维SLICE3_BOX最后两阶材料均值RMS差8.45585e-10、logZ差6.64293e-8；它固定其余坐标和有限box，不是FULL marginal，也不估计box外质量。LOCAL-C逃逸质量仍未知。

[强条件重要性失败](B/GN_DEPARTURE_AUDIT.json) · [实际两臂对照](C/SURROGATE_COMPARISON.md) · [冻结128对象的独立审计](AUDIT/locked_results_review/REVIEW.md)
