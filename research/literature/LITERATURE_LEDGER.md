# B1 文献与论断台账

核对日期：2026-10-13。范围是指定论文的归属、证据与缺口；这不是系统综述或全球创新性认证。主线程 Codex 负责最终数学、创新性和科学解释。本任务没有启动外部模型或 DeepSeek worker，没有运行新的物理计算，也没有修改源码、配置和结果。

## 可以直接使用的结论

| ID | 有界论断 | 证据与 B1 的边界 |
|---|---|---|
| L01 | 拓扑敏感性先验与 Bayesian scattering 模型选择已有来源。 | [Carpio–Iakunin–Stadler 2020](https://arxiv.org/abs/2003.09318)，B：摘要。其示例是二维合成全息问题；不能借此推断 B1 的三维空腔概率已校准。 |
| L02–L03 | 给定投影的条件似然 profile 与 KL 最优性、带条件的误差控制已有来源。 | [Zahm 等 2018 preprint / 2022 journal](https://arxiv.org/abs/1807.03712)，A：v4 §2.1，式(4)–(9)；§2.2，Assumption2.5。这些条件不自动给 B1 的 polynomial 插值提供认证。 |
| L04 | 非高斯低维、稀疏或可分解 transport couplings 已有理论。 | [Spantini–Bigoni–Marzouk 2018](https://jmlr.org/papers/v19/17-747.html)，B：正式 JMLR 摘要。它不证明 B1 的 Maxwell 专属签名或算法收益。 |
| L05 | 防御混合重要性采样是已有方法。 | [Hesterberg 1995](https://doi.org/10.1080/00401706.1995.10484303)，B：被索引的原出版社摘要。B1 的 defensive component 是 proper prior，权重含 prior/proposal 密度修正；有限样本漏模态与归一化偏差仍须保留。 |
| L06 | 多层球的电磁 Mie 计算是已有物理方法。 | [Peña–Pal 2009](https://doi.org/10.1016/j.cpc.2009.07.010)，B：被索引的原出版社摘要与 program summary；scattnlay 实现 Yang 的算法。论文不能替代 B1 实现的独立核验，也不能外推至任意几何。 |
| L07–L08 | SBC 检查指定生成模型内的 Bayesian 计算一致性，不保证模型捕捉现实真值。 | [Talts 等 2018，v2 2020](https://arxiv.org/abs/1804.06788)，A：§4，Algorithm1。其程序使用 posterior rank statistics；B1 当前的 probability bins、scores、coverage 本身不能被称为完成了该 rank-SBC。 |

A 表示取得可识别全文的相关段落，B 表示取得摘要；二者是文本可用性分级，不是论文质量或 B1 成功评级。JSON 中的每条实质文献论断都有 A/B 来源、定位与限制。

## 已核对的书目信息

| 文献 | 作者与年份 | 稳定标识与状态 |
|---|---|---|
| Bayesian approach to inverse scattering with topological priors | Ana Carpio, Sergei Iakunin, Georg Stadler；2020 | DOI10.1088/1361-6420/abaa30；arXiv2003.09318；S2与arXiv均解析。 |
| Certified dimension reduction in nonlinear Bayesian inverse problems | Olivier Zahm, Tiangang Cui, Kody Law, Alessio Spantini, Youssef Marzouk；preprint2018，journal2022 | arXiv1807.03712；DOI10.1090/mcom/3737。S2以arXiv解析、year=2018；[作者HAL条目](https://cv.hal.science/olivier-zahm)核对2022期刊及DOI。原 DOI unresolved 记录保留；AMS页面403。 |
| Inference via Low-Dimensional Couplings | Alessio Spantini, Daniele Bigoni, Youssef Marzouk；2018 | JMLR19(66):1–71，17-747；正式JMLR页面核对。S2不能解析该URL，不等于论文身份未确认。 |
| Weighted Average Importance Sampling and Defensive Mixture Distributions | Tim Hesterberg；1995 | DOI10.1080/00401706.1995.10484303；Technometrics37(2):185–194。S2身份解析；摘要由出版社索引取得。 |
| Safe and Effective Importance Sampling | Art B. Owen, Yi Zhou；2000 | [DOI10.1080/01621459.2000.10473909](https://doi.org/10.1080/01621459.2000.10473909)；JASA95(449):135–143。身份核对，但摘要未取得，保留为C元数据；不据此断言特定控制变量或方差定理。 |
| Validating Bayesian Inference Algorithms with Simulation-Based Calibration | Sean Talts, Michael Betancourt, Daniel Simpson, Aki Vehtari, Andrew Gelman；2018，v2修订2020 | arXiv1804.06788；S2、arXiv及全文核对。 |
| Scattering of electromagnetic radiation by a multilayered sphere | O. Peña, U. Pal；2009 | DOI10.1016/j.cpc.2009.07.010；CPC180(11):2348–2354；S2身份和原出版社索引核对。 |

## B1 当前边界与缺失证据

下面是对已有源码与回执的读取，不是本任务重新运行或独立确认物理结果。

- 算法结构是固定径向几何的显式 topology mixture 与连续复介电常数先验，分支 quadrature、防御混合重要性采样和 polynomial forward 积分；没有 NN 训练。
- `results/refinement/SUMMARY.json` 将改进参考标为经验细化，明确没有 uniform interpolation remainder；原 full-grid 的四个 `REFERENCE_UNRESOLVED` 保留。不能把 fine integration 一致性升级为严格 full-Mie 后验认证。
- `results/importance/RESULT.json` 记录每分支16个样本、合计48次新评价，以及 finite-N normalized-ratio bias 和漏模态风险。其相对 order7 的最大概率差为0.02565，不能写成采样准确性或成本优势已获证明。
- `results/CALIBRATION_RESULTS.json` 的720次噪声实现来自24个独立材料父样本；这是模型内开发检查。posterior-rank SBC、外部物理校准与任意几何支持证据均未建立。
- 源码明确让未建模几何的概率保持未定义。不能把模型内 Top-K 概率解释为所有可能真实材料世界的概率。
- Maxwell 专属 likelihood-ratio signature、topology coupling 的创新性以及理论是否被实验支持，仍待主线程的核心数学与最终科学审查。此台账没有作该判断。
- `results/RESOURCE_SUMMARY.json` 记载384/384次新增物理评价、remaining=0，training=0，GPU=0，locked arrays read=0。文献任务没有增加这些计数。

## 检索与访问限制

原 `scholarqa_evidence.json` 有6个HTTP429错误且无候选论文；没有重跑宽泛检索。本任务仅增加一次7篇的 ScholarQA metadata/abstract batch：6篇解析，Spantini的JMLR URL unresolved，未出现新429。该批次写入 `metadata_abstract_verification.json`，原失败记录未覆盖。

Hesterberg、Owen–Zhou和Peña–Pal的S2摘要由出版社隐藏；Hesterberg与Peña–Pal只使用被索引的原出版社摘要，明确标为B。Owen–Zhou停留在C。Hesterberg的CMU全文请求超时，Stanford报告库失败；AMS403；部分T&F/ScienceDirect/作者PDF直接访问失败。Zahm与Talts全文可取；Spantini取得正式摘要。一次误填的arXiv1805.09294指向无关论文，已排除并纠正为Talts1804.06788。

本任务没有做 citation-network 扩展、近邻全面覆盖或新颖性认证。所有论断、来源URL、错误、缺口和局部回执定位见 `CLAIM_LEDGER.json`。
