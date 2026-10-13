# Joint P2：有限文献核对与主张账本

2026-10-13。本文件核对当前实现能与既有方法区分到什么程度；最终科学方向、原创性与论文裁决由父 Codex 作出。

**现有证据不支持把“Bayesian 散射、多张重建、full-covariance MDN、SMC 或下游重新加权”单独作为原创贡献。** 更窄、可以由本项目实际数据支撑的内容是：在一个冻结的径向 vector-Mie prior/likelihood 中，定量分开材料与 gain 的补偿，保留显著材料事件的后验质量，并把条件材料不确定性传递给独立情景更新。它目前是一个可核验的受限实验结果，不能升级为全球首次、一般三维唯一性解决方案或 continuum 概率认证。

## 实际检索范围与访问缺口

使用本地 `scholarqa-research` skill 与 ScholarQA CLI 做 evidence-first 核对，没有调用 DeepSeek、其他外部推理模型或 Maxwell solver。原始检索与核验记录在 [`research/delegated/p2_literature`](../delegated/p2_literature/)；该目录是本项目证据归档，不是 Ai2 Scholar QA 的官方代码或发布。

预先限定两条 formulation：

1. `Bayesian inverse scattering nonuniqueness nuisance marginalization`
2. `posterior mixture density distillation sequential Monte Carlo likelihood correction`

一次 `collect` 的 paper/snippet 四项操作全部返回 **HTTP 429**，候选集为空，退出码 1，wall 36.5241 s。没有重试、翻页、引用图递归或额外宽泛检索。随后对已有报告中的 seed 和已识别基础方法做两次 batch verification：书目核验及带 abstract 的记录核验，五项均 resolved。直接查看已识别论文的 arXiv/publisher 内容属于补充原文核对，没有增加检索 formulation。

- Carpio–Iakunin–Stadler 与 Papamakarios–Murray 的 arXiv 原文可读，指定段落为 Tier A。
- Del Moral–Doucet–Jasra 的 abstract 可读，方法支持为 Tier B；publisher full text 重定向后访问失败。Semantic Scholar 返回 year=2002，而 DOI 的 Crossref 元数据明确给期刊发表日期 2006-05-15，卷 68(3)，411–436；本文件按 **2006 期刊记录** 引用，保留 metadata 冲突。
- Christen–Fox 的书目 resolved，但 abstract 被 publisher elide，DOI 访问 403；此处仅 Tier C，不用标题代替算法细节证据。
- Stuart 的书目 resolved；abstract 被 elide，Caltech full-text URL 访问失败；publisher 页面可读但未取得相应论证原文，此处仅 Tier C。

这不是系统综述或全面 novelty search。尤其未建立完整的 Maxwell calibration/nuisance、分层球形状–材料互补偿、posterior distillation 及下游决策文献覆盖；没有搜到某个机制不能说明它不存在。

## 已取得的主要先例

**散射后验不只包含 MAP。** Carpio–Iakunin–Stadler 的二维 Helmholtz/holography 实验展示：MCMC 可以保留与观测相容的不同形状家族，MAP 附近的 Gaussian 近似可能遗漏它们；第 7 节加入未知材料参数，第 8 节对模型参数积分计算 evidence。它直接排除了“多个散射材料解释”本身的新颖性。该论文的二维、形状参数及数据驱动 topological prior 与 P2 径向 vector-Mie 共同 target 不同，不能直接转用其数值成绩。[原文，§6–8](https://arxiv.org/pdf/2003.09318)

**条件材料概率生成和 full covariance MDN 已有先例。** Papamakarios–Murray 使用 full-covariance Gaussian MDN 表示参数条件密度，并在非原始 prior 下模拟时显式保留 prior/proposal density 修正。它研究 likelihood-free inference；P2 有可计算 Maxwell likelihood、独立参考 teacher 和可选 full-likelihood correction。因此它是明确的方法先例，并不是本项目的 teacher 成本、物理精度或泛化保证。[原文，§2.2–2.5 与 Appendix B](https://arxiv.org/pdf/1605.06376)

**温度路径与加权粒子不是新推断机制。** Del Moral–Doucet–Jasra 对一系列仅知未归一化密度的分布进行 sequential Monte Carlo，支持 weighted particles 和 normalizing-constant ratio 的方法来源。abstract 没有支持“有限高 ESS 保证不漏模式”；本项目仍须报告 seed 差异、粒子数、实际 full costs 与事件参考误差。[期刊 DOI](https://doi.org/10.1111/j.1467-9868.2006.00553.x)

## 主张—证据—限制

| 有界主张 | 支持与层级 | P2 中允许的表达 | 不允许的升级 |
|---|---|---|---|
| Bayesian inverse scattering 可以具有多个观测相容的形状后验家族 | Carpio et al.，§6，Tier A | 多解释的概率成像有明确先例 | P2 首次解决散射多解问题 |
| geometry 与未知材料可以联合放入 scattering posterior | Carpio et al.，§7，Tier A | 材料自由度可以与几何共同推断 | 本文已经证明 P2 的 gain 补偿或三维一般结论 |
| evidence 是 prior-weighted likelihood 积分，不能用最佳 residual 替代 | Carpio et al.，§8 Eq.18，Tier A；P2 自身概率定义 | 补偿最小 residual 是辨识诊断，posterior odds 由积分给出 | 小 residual 直接等于高 posterior mass |
| full-covariance MDN 与 prior/proposal 修正已有方法先例 | Papamakarios–Murray，§2/Appendix B，Tier A | P2 适配既有 normalized conditional density | full covariance、MDN 或简单密度变换本身是原创 Maxwell 算法 |
| sequential weighted samples 与 evidence ratio 有 SMC 先例 | Del Moral et al.，abstract，Tier B | generic tempered SMC 是认真基线 | 有限 ESS 或合法 MH 自动证明精确事件质量/全模态覆盖 |
| 原文曾提出 approximation-based MCMC | Christen–Fox，metadata，Tier C | 列为需要补全文的关联文献 | 用当前未取得的原文支持全部 delayed-acceptance 实现细节 |
| 概率 inverse-problem 理论有基础综述 | Stuart，metadata，Tier C | 作为背景书目保留 | 把未核验定理当作 P2 的 continuum error bound |
| posterior-weighted NLL 是 teacher cross entropy，等于 KL 加 teacher entropy | 本项目明确 target 下的概率恒等式；代码 `c1p2/train.py` | 已知概率训练目标，不是新定理；teacher 是额外监督 | distillation 的改善全部归于网络结构或免费 teacher |
| proper nuisance 应边缘化；独立 context 按 likelihood 重新加权 | 本项目明确联合分布与条件独立假设下的 Bayes 恒等式 | 自身数学 target 必须正确声明，未经核验不引外文结论 | 下游使用真值、重复使用 scattering y 或把 profile likelihood 当 evidence |
| full-minus-proxy residual importance 可以修正 event integrals | `src/b1_joint/correction.py` 的定义与已保存真实 full solves | 属于已知 importance/control-variate 组合；误差分层披露 | 本次有限核验点建立了 uniform Maxwell probability bound |

最后三项属于本项目数学/实现核对，而不是从 Tier C 文献推断的新科学发现。对应的新推导及最终有效性由父 Codex 审核。

## 当前数据支持的窄事实

在已保存的 `COMPENSATION.json` 中，固定参考实心与中心半径 .25R 的另一世界，在单频 1.6 下 raw whitened distance 为 **1.5296523**，只调 ε 后为 **1.2952214**，ε+gain 后为 **0.9001700**。后者 ε≈2.46662+0.13001i、gain≈1.05322，处于声明支持内部；增加 2.4 频率后相应 ε+gain 距离为约 **2.468**。这些是可行、best-found full-Mie endpoints，仍是补偿最小值的上界，不是全局最优或 posterior odds。

独立的 `FULL_MIE_CORRECTED_REFERENCE.json` 对预注册 large-cavity observation，给实心事件 **0.1625397**、中心半径 (.10,.20] 且无 annular cavity 事件 **0.0718280**、中心半径 (.20,rmax] 且无 annular cavity 事件 **0.6772260**。报告的是材料区域质量；不能给单张连续图像这些概率，也不能未经 mode analysis 称它们为三个分离 modes。该 empirical residual reference 的 Monte Carlo intervals 不包括所有 solver/global-support 误差。

因此，最可用的研究陈述是“在冻结 prior/noise 的受限 vector-Mie 校准台中，以 full nonlinear secant 与 nuisance-marginalized event reference 共同量化结构歧义”。具体单频/双频差异、reference 数值问题、生成器的 false exclusion/false confidence 和下游 regret 才是需要真实结果判别的内容。**本账本不宣布 IEEE TAP 原创贡献成立**，也不把目前检索缺口解释成新机制不存在或存在。

## 仍不能由本次证据认证的内容

1. 一般非同心三维未知形状、无限维材料支持、真实仪器 calibration 或模型外 posterior mass。
2. 有限点 response error 到所有重要事件 odds/CDF 的 uniform bound。
3. full covariance distillation 在所有 Maxwell 后验上优于 diagonal 或 Bayesian sampling。
4. sampled generator frequency 自动等于 full-Mie posterior；support 保留或高 ESS 自动排除漏 mode。
5. Top-K 代表图可替代事件内部条件样本；context 更新必然优于 MAP。这需要本项目实际下游比较。
6. “Maxwell 补偿 + 概率保留 + downstream handoff”相对全部既有文献的全球原创性。

## 已核验书目

- Ana Carpio, Sergei Iakunin, Georg Stadler (2020), *Bayesian approach to inverse scattering with topological priors*, Inverse Problems 36, 105001. [DOI](https://doi.org/10.1088/1361-6420/abaa30), [arXiv](https://arxiv.org/abs/2003.09318). S2 `cd4deaaee416743d1b4a783f51dc9edee6a82a3e`。
- George Papamakarios, Iain Murray (2016), *Fast ε-free Inference of Simulation Models with Bayesian Conditional Density Estimation*, NeurIPS 2016. [arXiv](https://arxiv.org/abs/1605.06376). S2 `4123ac46635b015227f06370e367c828d6ccffd3`。
- Pierre Del Moral, Arnaud Doucet, Ajay Jasra (2006 journal record), *Sequential Monte Carlo Samplers*, JRSS B 68(3), 411–436. [DOI](https://doi.org/10.1111/j.1467-9868.2006.00553.x). S2 `3ce93c1a772061161948a6728716f6f4d752809b`；S2 year=2002 的冲突保留于 raw record。
- J. Andrés Christen, Colin Fox (2005), *Markov Chain Monte Carlo Using an Approximation*. [DOI](https://doi.org/10.1198/106186005X76983). S2 `7e0f1514ae7bae557116e8b689d1fb2b2c08960c`。仅书目核验，全文缺口保留。
- Andrew M. Stuart (2010), *Inverse problems: A Bayesian perspective*, Acta Numerica. [DOI](https://doi.org/10.1017/S0962492910000061). S2 `a136f1140f2686edf51eb3cf5b44f36c51ef730a`。仅书目/关联背景，未用其未取得全文支持定理。
