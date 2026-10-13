# 理论状态：哪些关系成立，哪些机制仍有限

实验阶段20261013；公开交付2026-10-14。Fresh主概率参考与范围以 [ACCEPTED_REFERENCE_AUDIT](results/ACCEPTED_REFERENCE_AUDIT.json) 和 [全部fresh细化记录](results/ALL_FRESH_REFERENCE_REFINEMENT.json) 为准；旧13阶与中间结果保留为历史近似。

P2 已把材料 prior、nuisance likelihood、事件质量和下游更新放入同一数学 target。它没有证明一般 Maxwell inverse problem 的 uniqueness，也没有证明 learned generator 的全域 posterior calibration。

| 项目 | 状态与实际证据 | 允许的结论 |
|---|---|---|
| Coefficient prior 到 material pushforward 的混合测度 | 已解析推导，E0 与独立积分差 ≤2.22e−16 | zero-radius 原子与 inactive latent 必须保留；radius-uniform 是另一 prior |
| Proper truncated shared gain marginal | 独立 quadrature/Cholesky loglik 差 3.55e−15 | logdet、prior 截断归一化与两频共享 gain 在实现中保留 |
| 连续 posterior summary coherence | 七个 C development cases 采用同密度的 moments、直接 CDF 和反演 quantiles | 不能再混合节点 moments 与 prior-bin quantiles |
| Fresh numerical reference | 全128event/mean/SD13/17核对，12/39/48升级17并对照21后无未解决项；全128positiveCDFgrid单调缺陷0 | 不作uniform-bin展宽/quantile；CDFrefinement为固定8+升级3，不认证所有阈值/continuum |
| 局部 projected Jacobian | 五个固定半径的数值诊断已保存 | 描述局部 nuisance 可吸收方向；不推出大空腔的非线性等价性 |
| Full nonlinear material secant | 单频 r=.25R 的 ε+gain endpoint distance .90017，双频 2.46794 | 存在一个合法、有界补偿 witness；不是全局最小或严格等价 |
| Posterior event ambiguity | Full-Mie residual reference 保留实心、中心空腔等显著事件 | 连续材料族可以同时有质量；不是给代表点非零概率 |
| Full/proxy event correction | 实际 full solves、known proposal、2000 residual bootstraps | empirical 数值参考；没有 uniform event/CDF error certificate |
| Full-likelihood MDN correction / delayed-acceptance SMC | 合法权重或 MH correction，实际全成本记录 | correction 指向声明 numerical full target；有限 draws/particles 仍可漏事件 |
| Distillation covariance 优勢 | 同 teacher 的 fresh comparison；event MAE 差区间跨零 | full covariance 可改善 latent density，但不能宣称事件质量全面胜出 |
| 下游 independent context update | 明确 Gaussian geometry likelihood；实际64parent paired regret下降.247748 | 仅在声明条件独立下合法；单点baseline非全域MAP，模型外概率仍无法量化 |
| A17 原双侧材料因子、A30 全材料可观测谱 | NOT_RUN：未接入为 P2 算法模块 | 只能作理论背景，不能声称产生已验证计算收益 |
| 非同心三维、真实仪器、多物体、无限维 posterior | NOT_RUN：不在授权范围 | 径向 vector-Mie 证据不能升级为一般开放世界成像 |

核心推导及其范围在 [THEORY.md](THEORY.md)。Bayesian scattering、多解释输出、Gaussian MDN、tempered SMC 和 importance correction 都有既有先例；有限文献核对不足以建立全球原创性，原始 429/全文访问缺口保留于 [文献账本](LITERATURE_AND_CLAIM_LEDGER.md)。

[THEORY_HANDOFF.md](THEORY_HANDOFF.md)明确证明eventconditionalcontextlikelihood期望、TV/regret误差界及latent/materialKL链式分解。它们为已知概率关系，误差界需要另有uniform条件常数；本实验没有因此获得全context或现实probability认证。

科学上最可继续的方向是 event-aware active sensing：在保留 ambiguity 的前提下判断何种额外 acquisition 能减少与实际决策有关的材料歧义。新增 2.4 频率的补偿距离与 posterior event 变化提供了具体起点，但尚未完成 acquisition optimization 或一般预测理论。
