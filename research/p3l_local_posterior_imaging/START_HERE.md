# Research P3-L：非线性 Maxwell 材料概率成像

**运行中，2026-10-22。完整 P1 posterior reference 尚未合格，不能宣称已实现可信 Bayesian 成像或概率校准。**

本轮在原 P1 的 66 维 prior、TM24 full-Maxwell、噪声和观测合同下运行。研究目标是降低材料风险；候选数量、低 residual、接受率和单点图像都不代替 posterior 正确性。

已经实际完成：

- 32 个独立 development 对象的 full-physics 基础反演：32 个进入原 noise-ball，27 个达到冻结驻点标准。
- 24 条完整 66 维采样链，每条 1,024 warmup + 4,096 production；已保存拒绝时的重复状态。目前六组四链参考全部未通过混合诊断。
- 另完成四条 HMC pilot，每条 512 warmup + 1,024 production。材料函数最大 Rhat 为 1.26968，最小 bulk ESS 为 11.47，材料均值 RMS MCSE 为 0.01015；比先前 MALA 更好，仍未合格。
- 4 个对象、5 种弱维数、两种外层坐标，共 40 个强条件实验。128 点 full-likelihood 重要性校正中，38 组 ESS 小于 2，不能作 teacher。
- 独立 128 粒子 tempered SMC 到达完整 likelihood，但祖先退化到 1；有限粒子结果仍未合格。
- 同一原物理模型上的三维受限条件积分：最终两阶材料均值 RMS 差约 `8.46e-10`，log-normalizer 差约 `6.64e-8`。这是 **SLICE3_BOX**，不是积分其余 63 维后的 FULL marginal。
- 全部 32 个 development 对象的离线成像对照：固定属性 GN Gaussian 参考的材料均值 MSE 为 `0.0206962`，优化基础图为 `0.0138212`。近似均值反而恶化约 49.7%；这不是对真正 posterior mean 的否定。
- 512 train / 128 validation 的独立观测、在线物理 context 与合法 joint-supervision 标签。128 个 locked parent 尚未读取或生成。

内存中断已单独记账。上游缓存对旧 LU 的强引用导致增长；本轮适配器现只保留 4 个 transfer，完成对象后释放引用，全局锁限制一次一个 worker。640 个在线 context 已恢复完成，实测约 443–448 MiB。中断前的 109,488 RHS 和失败不删除、不重新归零。

正在进行联合属性–空间场 Gaussian 参考的独立修订，以及保留完整 66 维支持的非线性三角补偿坐标。随后运行 joint-supervised Gaussian、MDN、Flow 与 CNN 对照；未合格链绝不用于 posterior distillation。

16 份固定提案的实际求值对照也已完成：普通 Galerkin 与同空间 A17 的差异仅为浮点误差；rank-16 条件误差区间可提前判定 5/16 次接受决定，但 setup、界计算与 fallback 均需计费，尚无总成本加速结论。这里的区间基于浮点误差模型，不是定向区间运算认证。

[中文研究记录](RESEARCH_REPORT_ZH.md) · [概率合同与算法](THEORY_AND_ALGORITHM.md) · [冻结配置](configs/PROTOCOL.json) · [资源快照](RESOURCE_LEDGER.json) · [内存恢复证据](AUDIT/MEMORY_INTERRUPTION_RECOVERY.json)

原始 chains、端点、LU cache 和追加账本保留在独立目录，待形成证据 Release。当前不是最终结论，也未将正在运行的模块标为 PASS。
