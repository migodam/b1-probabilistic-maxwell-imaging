# 给 Codex 的执行任务：MAX-STRUCT 有限材料事件与可访问性

## 目标

从固定仓库提交 `ec3a37da6f3ee90beb61e753067d4f639f4118f5` 和本包出发，落实一个可审计、低额外成本的整数事件proposal，并用真正独立的确认检验它。不要从零重复R4–R12，不要将已有投影、NFFT解锁或变量投影恒等式当新发明。

先读 `report.md`、`theory_audit.md`（含最后的Maxwell桥接与条件审查）、`evidence_audit.md`、`paper_experiment_protocol.md`。后者给出完整schema/资源/原论文复现路径。沿用repo的实际AGENTS与operator catalog；不得覆盖历史结果。

## 立即实施的工作

1. 在独立研究目录建立current adapter与候选审计调度器。G19对照至少包括legacy、从第一步full APCS、NFFT full current；固定相同材料块、完整分母导数、beta和PR约定。soft DPCS与额外ball约束分开登记。保留R8 physical-complement作历史控制。
2. 实现实共享材料事件：内点和活动盒边界闭式t_star，暗cell/无事件、多个并列事件、整数tie。仅沿一个全局scalar放大法向，不逐cell独立跳。精确pair跳步与近似pair proposal使用不同mode。
3. 先用本包冻结的3个AP状态确认实际大跳变时刻/位置与保存材料；不是只做机器精度恒等式。对全部11个历史起点开展因子研究，并保留所有失败。
4. 统一付费物理审计机会，明确mf12中途候选。分别运行“相同proposal工作+审计槽”和“相同总CPU预算”的比较；两者不是同一种公平性。
5. 比较ordinary DR、integer event、原fractional event、AP、reset DR、同范数随机扰动、普通倍增外推及适当加速基线。准备AP状态、检查条件、被拒绝提案和所有物理审计全部计费。
6. 保存lambda来源和源—接收KKT缺陷。不得混用PA(c)输出a'的lambda与旧a的残差来宣称精确身份。优先复用G*/S*与已付数据，仅小诊断子集新增独立adjoint solve；proxy一致不算独立机制验证。
7. 分离搜索、物理repair与prior选择。保留raw/repair/corrected/selected；已在ball内时最小位移repair原样返回。该规则只保证不移动可行材料，不保证GT。
8. 冻结生产式事件触发/候选选择后生成96核心+24分布外的新2D对象。主确认含matched/cross-grid、cold/NN、9方法；按object block保存全部预测，再允许独立scorer读GT。不得用既有暴露对象反向挑阈值后称盲测。
9. 继续原NIE/CIE设置、成熟MR/CC-CSI和field-ADMM，24个真三维对象与Fresnel标定/验证。复用R9既有三维结果及其跨网格误差；不把点DDA matched-grid图当连续Maxwell证明。

## 必须保留的理论边界

- 事件主定理覆盖固定度量、非负有界实各向同性材料；complex loss/tensor需独立投影推导。
- 非零AP gap的receiver-adjoint公式需要精确pair、可逆resolvent、delta>0与Slater；近似时记录缺陷。
- 精确跳步不改变后续轨迹，因此不可能解释恢复改善；近似事件改变轨迹，无GT保证。
- 规定法线射线上阈值正确、实际DR首事件正确、最终球内、材料正确、端到端更快，是五个不同命题。
- 已经数据可行的错误物理材料可以是AP/DR固定点，不能用纯反射保证GT。
- fullcurrent profile罚项的精确Hessian离开兼容零残差根不必PSD；legacyG19也未被证明有单个共同profile目标。

## 输出与论文

输出 `case_manifest, runs, candidate_registry, physical_audits, events, cost_ledger, profile_diagnostics, offline_scores`，字段遵照完整协议。记录GT访问、配置hash、停止原因、墙钟/CPU/峰值内存、物理和辅助算子RHS、setup/训练/拒绝trial。缺失用null。

主图先展示首事件预测及失败，再展示匹配表示/审计控制、corrector材料损伤、独立成本—质量分布、真正三维与实测的物理残差结构。不要只排平均NRMSE排行榜。

负结果应引出误差条件、活动面、可观测性或模型误差分析，不设置任意成功率gate来终止全部研究。若事件只是普通扰动的另一写法，明确削减算法主张并继续完成表示/选择机制研究。主稿以有限材料分支访问为中心；强TAP依赖可靠电磁验证，SIAM JIS依赖理论与独立预测，Inverse Problems需要更深的可辨识/稳定性内容。
