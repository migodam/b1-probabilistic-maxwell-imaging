# B1 × C1 Joint P2：同一测量保留不同材料世界

实验阶段编号20261013；公开交付日期2026-10-14。最新主参考为 [ACCEPTED_REFERENCE_AUDIT](results/ACCEPTED_REFERENCE_AUDIT.json)：全部128新观测的event/mean/SD已复核，三例细化后无未解决项；CDF仍限于声明的数值检查范围。

已实际实现并运行：统一连续材料先验、3Dvector-Mie、正确的 nuisance 边缘化、连续概率reference、C1四维density、合法full-likelihoodcorrection，以及可供下游独立context重新加权的多候选package。它是**同心径向材料校准台**，不是一般未知形状三维成像。

最直接的demonstration是预注册`solid`观测：散射后验同时保留实心 **40.77%**、较大中心空腔 **41.09%**、中等中心空腔 **11.63%**。传入独立带噪声几何传感器后，实心变为 **72.97%**，较大空腔仍保留 **9.01%**。这些是互斥材料区域的mass，不是单张continuous图像的点概率；余量、条件不确定性和未知模型外状态都保留。

按阅读目的选择：

- [可浏览Top-K图库](index.html)：材料实部/损耗截面、事件质量、补偿与下游图。
- [中文独立研究报告](RESEARCH_REPORT_ZH.md)：从物理歧义到概率语义、实际结果和十项科学回答。
- [large-cavity hypothesis JSON](results/packages/large_cavity.json) / [solid JSON](results/packages/solid.json)：机器可读事件、材料、mass/error、forwardpredictions、uncertainty和contextupdate。
- [共同后验参考](POSTERIOR_REFERENCE.md) / [理论状态](THEORY_STATUS.md)：测度、nuisance、连续摘要、full/proxy分离与误差界限。
- [方法比较](METHOD_COMPARISON.md) / [覆盖](MODE_COVERAGE.md) / [下游decision](DOWNSTREAM_DECISION_TEST.md)：同target公平比较及完整freshpanel。
- [科研裁决](GATE_DECISION.md) / [资源](RESOURCE_SUMMARY.md) / [失败](FAILURES.md)：真实正负结果、3459/4000attempt与保留的故障。
- [代码与复现](README.md) / [研究日志](RESEARCH_LOG.md) / [文献主张账本](LITERATURE_AND_CLAIM_LEDGER.md)。

C1模型只在共同proxytarget上训练。采用最终接受的事件参考，三seedfull/diagonaldistillation事件质量MAE **.013215/.013028**，差的区间跨零。旧13阶node-expectation latentKL约 **.21572/.27482**保留为历史经验近似，未重算成完整材料KL。最终acceptedCDF条件误差 **.141761/.161768**仍远高于.01目标；fullcorrection不提供finite-sample全mode或continuum认证。

数学交接条件与误差界见 [THEORY_HANDOFF.md](THEORY_HANDOFF.md)：eventmass必须乘事件内contextlikelihood期望；代表图一般不足以提供该期望。TV/regret界是条件数学关系，没有在本实验取得所有context的uniform误差证书。

当前科研方向是继续受限ambiguitybenchmark与下游handoff、优先主动测量；停止通用r⁵leadingcancellation及NN全面优越性主张。大三维非同心campaign、大型GPU训练、旧lockedtruth和A17稿修改均NOT_RUN。旧B1/C1历史结果原样保留，公开更新使用现有两个仓库，不新增repo。
