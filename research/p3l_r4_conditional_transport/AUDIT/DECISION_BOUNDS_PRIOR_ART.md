# 接受决定驱动求解的先例与本轮主张边界

2026-10-25核对以下作者/期刊或机构主来源，避免把已有原则重新命名为新 Maxwell 理论。

1. [Chanialidis、Evers、Neocleous、Nobile，2014，Retrospective sampling in MCMC with an application to COM-Poisson regression](https://eprints.gla.ac.uk/96718/)。机构正式记录说明利用 likelihood 上下界来约束 MH 接受概率，并在足以接受/拒绝时停止加精。该工作针对离散分布，不是 Maxwell 求解器，但“先决定，再按需要收紧界”已有明确先例。
2. [Sherlock 等，Efficiency of delayed-acceptance random walk Metropolis algorithms](https://arxiv.org/abs/1506.08155)。便宜近似筛选与完整目标校正的 delayed-acceptance 机制亦为已有方法。它通常不要求与一次标准 Full-MH 的每一步随机决定完全相同，应与本轮 common-random-number 重放的目标区分。
3. [Adaptive, delayed-acceptance MCMC for targets with expensive likelihoods](https://arxiv.org/abs/1509.00172)。近似自适应及代价权衡不是本轮原创性来源。

本轮实际贡献候选仅限于：原P1离散TM24中，怎样利用已付费状态逆作用、残差界与普通迭代方法，在完整计入setup、数据持久化和回退后达到相同Full数值决定；并与享有同样R0缓存权限的普通global-GMRES比较。不得以普通Krylov也能实现的收益证明A17独立增量。

本轮的包络使用float64数值范数及回退，并非定向舍入的连续PDE证书。即使所有已测试接受决定、状态路径与势包络均通过Full重放，也只能报告已验证数值一致性，不能将其升级为任意状态下的数学机器精确证书，或完整posterior已经混合的证明。

上述文献核对主要为摘要/正式记录范围；未读核验的具体定理、复杂度常数和数据集不作转述。最接近的Maxwell/PDE特化prior art仍不构成全面穷尽证明。最终TAP原创性裁决由主线程结合真实实验完成。
