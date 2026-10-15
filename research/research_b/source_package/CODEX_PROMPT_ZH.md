# 给 Codex 的执行 Prompt：Research B

你负责独立项目 **Research B — Event-Aware Maxwell Active Sensing**。这是 A 的并行项目，不等待 A 新生成器，也不把 B 的预算借自 A/B1/C1 的历史剩余额度。

## 先读
本研究包：
- RESEARCH_B_REPORT_ZH.md
- SOURCE_AND_NOVELTY_LEDGER.md
- configs/EXPERIMENT_PROTOCOL.json
- schemas/MaterialHypothesisPackage.schema.json
- schemas/AcquisitionRecommendation.schema.json
- code/acquisition_core.py 与 code/verify_theory.py
- checks/THEORY_CHECKS.json

实际来源：
- migodam/b1-probabilistic-maxwell-imaging 的 joint-p2-20261014 分支
- research/joint_p2/ 下 RESEARCH_REPORT_ZH、AMBIGUITY_AND_COMPENSATION、DOWNSTREAM_DECISION_TEST、POSTERIOR_REFERENCE、MODE_COVERAGE、THEORY_HANDOFF、THEORY_STATUS、GATE_DECISION
- ACCEPTED_REFERENCE_AUDIT、FRESH_CONTINUOUS_REFERENCE_ACCEPTED 与真实 package/bank 的 README/provenance
- C1 p1/p2 的报告及 continuous_accepted 结果
- A17 closure_final、A23、A24、A30 对应已冻结研究结论。先记录实际 revision，不依赖名称推测当前版本。

## 本次可直接做的工作：仅 E0，Maxwell 新调用上限 0

1. 建独立工作目录，不覆写历史报告、posterior target 或原始数组。
2. 先执行本包 11 项无物理调用检查；检查代码的输入合同是否足够，保留失败与修复，不把通过 unit tests 写成物理 PASS。
3. 查找实际已保存响应、系数和样本；建立 E0_ACTION_INVENTORY。未知缓存 location 必须找到后再读，不编造本地路径。cache miss 要 fail closed，绝不暗中启动 full solver。
4. 实现 P2-to-MaterialHypothesisPackage adapter：32 real channels/frequency 的实际 packing、实噪声 sigma²、proper truncated shared real gain、原 coefficient-uniform prior 与几何零原子、全部事件和未展示 remainder、full/proxy/reference 身份必须保留。
5. 保存 theta/gain 联合粒子，或从已有数据重建条件 gain posterior。新频率的增量 likelihood 必须是 log p(old,new|theta)-log p(old|theta)，不得每频重新积分原始独立 gain prior。
6. 在已有 1.6/2.4 响应上执行 E0。比较 preset2.4、repeat1.6、可合法提取的接收/偏振子块、no-measurement。所有新 noise 都标模拟。已被当前 posterior 使用的历史通道不能当新增数据；retrospective subset 必须先移除并重算 posterior。
7. 用共同 predictive random numbers 比较 raw medoid separation、经典 Fisher、堆叠后共享 nuisance Schur、event EIG、expected material MAE/候选选择风险、鲁棒复合假设。完整样本、medoid-only 与有限条件样本在同 target 上比较。禁止使用 truth-state J、真实事件标签或最有利新噪声来选测。
8. 先冻结材料 MAE 度量、事件定义、参考数值误差与统计单位。不要把 medoid selection risk 当无约束 posterior median MAE；topology/boundary 单列。
9. 校准错误模型、prior perturbation、event threshold perturbation、sphere rotation/polarization span、rare/missing-mode 负控制分别报告。不要从某一有利频率结果推成优化成功。
10. 不能做到的项标 NOT_RUN/NOT_AVAILABLE，保留原因。只有参考数值精度足够时才给概率 utility；proposal score 不足时改为声明集合/设计先验下的鲁棒结果。
11. 输出 E0_DECISION.md、逐 parent/行动原始 utility、数值误差、条件 event update 分布、材料 MAE 与事件风险、ActionRecommendation examples、failure ledger 与总成本。实际 E0 的新 Maxwell attempts 必须为 0。

## E1/E2/E3 的边界
配置中的 E1 4096 是理论建议，不是已授权执行。先完成一个冻结的 E1 plan 和可实施仪器/模拟 menu，成本包括四新频率 bank、独立核对、fresh truth、finalist verification、失败与细化。不要借 Joint P2 剩余 541 次。

只有用户另行明确批准 E1 的独立预算后，才把 execution_authorized 改为 true 并启动对应 campaign；E2/E3 各自也需要独立预算。不要自行大规模训练、上 3D/GPU、增加场景或改变噪声来追求成功。

## 研究主线
优先验证“跨配置共享补偿不相容”是否比代表图 raw separation/Fisher 更能预测材料风险收益；允许校准测量胜出，允许 uniform sweep 胜出，允许所有可付费动作都不值得做。球对称已有双极化仅支持冗余/SNR 控制，不能冒充非球矢量极化辨识。

图和 ROM 先作为可删除模块。ROM 排名区间必须区分真正 bounds、MC-only 区间和有限 empirical audit；只有 inclusive time-to-quality 改善才保留。A17 是机制依据，不是已经证明的本项目加速。

不要强迫输出唯一图像。明确区分 STOP_TASK_SUFFICIENT、STOP_NOT_WORTH_COST、PRESERVE_AMBIGUITY、EXPAND_SUPPORT_FIRST、数值排名未解决和预算用完。

先给一份短实施计划和 cache inventory，然后直接完成零调用 E0 的可做部分；不要用重新写一篇泛泛综述代替实现。遇到不可恢复的资产缺失，保留具体阻塞和已经完成的产物，不编造实验结果。
