# R6论文 prior-art 最小核对与claim ledger

检索日2026-10-27。使用 **scholarqa-research**，primary transport=`scholarqa-cli`；小型证据QA，不是系统综述或全球新颖性认证。数学有效性、论文贡献与最终创新裁决留给root。

## 检索及负证据

本轮四个主题：①KKT补全残差梯度图像feature；②native→probit阶段预算；③任务加权分层后验积分；④conditional gauge/copula造成材料均值偏差。初始3查询+1阶段分配查询；对两个429缺口仅用较宽术语各补1查询、关闭snippet，不再循环重试。

共10个检索operation，其中7个429；有少量paper/abstract成功。12个candidate记录大多不直接相关，最终仅保留3个近邻。2+1条metadata batch verify均成功，无unresolved ID。不可访问的DOI路由错误保留为限制；随后仅通过已识别PMCID访问作者全文。没有使用泛用web搜索替代失败的学术检索，也未重复换MCP transport搜索。

原bundle：`FINAL_SCHOLARQA_BUNDLE.json`、`FINAL_STAGE_ALLOCATION_SCHOLARQA.json`、`FINAL_BROADENED_GAPS_SCHOLARQA.json`；核验：`FINAL_SCHOLARQA_VERIFY.json`、`FINAL_KKT_NEIGHBOR_VERIFY.json`。primary网页仅用于核对已识别论文，不作为另一次搜索。

## Claim ledger

|R6拟讨论的有界主张|论文/标识与证据|tier / relationship / confidence|R6仍未解决的界限|
|---|---|---|---|
|模型导数/优化更新作为学习重建信息已有先例|Manninen、Gröhl、Lucka、Hauptmann，arXiv:2510.27487，S2 `37fcdf951ef28c540f3671f09950fba925d47da6`。作者摘要描述把模型信息送入学习迭代，并比较gradient、Gauss–Newton和Quasi–Newton更新。[作者预印本](https://arxiv.org/abs/2510.27487)|B；contextualizes；对广义feature先例中等，对R6精确闭合匹配低|未在该摘要看到把残差反投影以Gaussian-prior驻点z替换、保持同一旧head或σ比例闭合的精确算法。不能由未看到宣布首次。|
|有预算的先行阶段→常规优化器handoff已有近邻|Ruo-Qian Wang，2026，arXiv:2610.22379，S2 `6dd690f3af54565170e5f1432349d37f394a6727`。作者摘要提出先由LLM给出warm start，再把已累积评估交给常规优化器。[作者预印本](https://arxiv.org/abs/2610.22379)|B；contextualizes；对广义两阶段预算中等，对native/probit精确匹配低|这是LLM/常规优化的工程校准，非本项目两种坐标目标的Maxwell优化。不能当作R6精确算法的重复证明，也不能把分阶段本身当新贡献。该页面是预印本。|
|QoI导向的Bayesian inverse误差控制与停止有明确先例|Kerfriden、Kundu、Claus，2019，DOI10.3390/ma12040642，S2 `3940e1a37a614b500fafd92da8f909630ca5946c`。作者全文引言/方法概述以QoI影响选择离散分辨率，联合控制MCMC sampling与离散误差，并据误差决定停止或refine。[作者全文PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6416561/)|A；qualifies broad novelty；高(广义先例)，低(精确算法等价)|不是R6的饱和平台、源域纤维、task-amplitude遗漏界或分层ratio分配实现。其存在足以要求R6收窄“首个task-oriented inference/stopping”式叙述。|
|条件gauge/copula可改变材料均值而边缘保持相同|两次相关查询都429；无本轮合格外部文献记录。内部E[QD]−E[Q]E[D]=Cov(Q,D)和纤维条件论证仍可由root直接核对。|无tier；SEARCH_GAP；无文献结论|不能称首个发现copula依赖影响pushforward mean；不能把缺少检索结果当作反证或新颖性支持。也不能由真实Cov非零推出任意具体Gaussian/copula方法必然同方向偏差。|

## 书目信息差异及核对

- Manninen等：S2题名为 *Toward Robust…*，作者arXiv v2为 *Towards robust…*；S2 year2025对应预印本首次提交，作者页显示v2修订2026-04-23。期刊DOI元数据为TCI2026。本文引用明确使用[arXiv v2](https://arxiv.org/abs/2510.27487v2)，不把预印本年份与期刊年份混写成冲突结果。
- Kerfriden等：PMC作者/2019/卷12(4):642/DOI与S2 verify一致；直接DOI入口本次工具不可访问，但已识别PMCID的作者全文可读。Tier A限于实际全文支持的QoI导向误差/停止主张。
- Wang：作者arXiv页的题名、单作者、2026提交与S2一致；本轮只读摘要，标B，不把其摘要报告的经验数字作为本项目机制证据。

## 对论文措辞的建议（分析者建议，不是创新裁决）

可清楚陈述R6做了一个具体的**旧权重、已知noise倍率、低成本anchor的KKT特征近似**，并把同端点actual-gradient控制与成本门分开检验。基础stationarity恒等式、单位球投影1-Lipschitz、QoI导向计算、两阶段handoff与Cov乘积偏差不宜单列为通用原创定理。

可能需要root最终评估的贡献只在：特定饱和Maxwell decoder/目标转换如何形成该闭合；独立两噪声实验是否同时支持风险与成本；同端点控制是否支持机制额外增量；分层材料商与任务积分在真实模型中是否给出可复核的额外资格或费用优势。此最小检索没有足够证据裁定精确KKT闭合或native/probit预算组合的全球新旧。
