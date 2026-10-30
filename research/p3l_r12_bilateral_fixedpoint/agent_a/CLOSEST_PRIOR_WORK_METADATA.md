# R12 最近先例身份审查：受限检索，不作全局新颖性判断

2026-10-30 UTC，按 `scholarqa-research` 技能用Semantic Scholar CLI检索。问题聚焦双边distorted-Born对角材料Jacobian、端口对/Gram，以及SOM/TSOM。执行两条初始query（每条至多6条paper、3条snippet）和一条减少负载的paper-only follow-up。初始四个操作中三个遭HTTP429，只有TSOM paper检索成功；follow-up的distorted-Born paper检索也遭429。因此没有获得可用的逐段全文证据，不能完成实质最近工作/反证覆盖。

以下5条由 `scholarqa-cli verify` 批量再次核实，5/5 resolved、unresolved=0。证据等级均为**Tier C身份信息**。即使个别检索记录含abstract，本小审查也不从中推导方法效果、等价性或新颖性。排除了一条神经网络3D重建题目，因与R12冻结局部几何的直接关系不足。

| 已核实身份 | 年份 | 稳定标识 | 本次可支持的内容 |
|---|---:|---|---|
| Zhong、Chen，Twofold subspace-based optimization method for solving inverse scattering problems | 2009 | [DOI10.1088/0266-5611/25/8/085003](https://doi.org/10.1088/0266-5611/25/8/085003) | twofold subspace inverse scattering先例的出版身份 |
| Zhong、Chen，An FFT Twofold Subspace-Based Optimization Method for Solving Electromagnetic Inverse Scattering Problems | 2011 | [DOI10.1109/TAP.2010.2103027](https://doi.org/10.1109/TAP.2010.2103027) | FFT/twofold主题先例的出版身份 |
| Xu、Zhong、Song、Chen、Ran，Multiplicative-Regularized FFT Twofold Subspace-Based Optimization Method for Inverse Scattering Problems | 2015 | [DOI10.1109/TGRS.2014.2329032](https://doi.org/10.1109/TGRS.2014.2329032) | MR-FFT-TSOM身份；年份为2015，DOI含2014不等于出版年 |
| Zhang、Ma、Xu、Zhong，Wavelet-Based Subspace Regularization for Solving Highly Nonlinear Inverse Scattering Problems with Contraction Integral Equation | 2020 | [DOI10.3390/electronics9111760](https://doi.org/10.3390/electronics9111760) | wavelet/subspace regularization主题先例的出版身份 |
| Huang、Hu，Accelerating Uncertainty Quantification for Nonlinear Inverse Scattering Problems With High Contrast Media by Direct Envelope Methods and Krylov Subspace Iterative Integral Equation Solvers | 2021 | [DOI10.1109/TGRS.2020.3035641](https://doi.org/10.1109/TGRS.2020.3035641) | 高对比度逆散射/UQ/Krylov主题先例的出版身份 |

这份身份清点不证明R12算法相对上述文献有增量。R12当前推导中的resolvent导数、Khatri–Rao列、Hadamard Gram、实decoder拉回与对同一Jacobian嵌套SVD的冗余均是线性代数结构；本轮不把它们作为新定理。潜在算法贡献须来自原先验支持内、真实非线性校正后的材料候选、覆盖、稳定性及总费用，最终由主线程审查。

检索证据：`SCHOLARQA_SEARCH.json` 保留partial成功与三个operation_errors；`SCHOLARQA_VERIFIED.json` 保留五条核实记录；`SCHOLARQA_DISTORTED_BORN.json` 保留follow-up的429。空ID的follow-up verify因无候选返回用法错误，不计作成功核实。此次检索未覆盖全部distorted-Born、Jacobian张量分解、实验设计、广义奇异值或prior-aware流形优化文献，不作“首次”或“全局原创”结论。
