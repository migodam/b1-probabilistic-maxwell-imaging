# 已有工作与本轮可主张的范围

本轮2026-10-26首先实际调用 `scholarqa-research`。沙箱内DNS失败；获准公开检索后Semantic Scholar返回429，原始失败均保留于audit/SCHOLARQA*.json。没有将检索缺失解释成原创性。以下为R5已识别论文的作者原始页面核对，不称完成全面文献检索。

|已核对的原始来源|可直接支持的已有机制|限制本轮怎样主张|
|---|---|---|
|[Bardsley等：Scalable optimization-based sampling on function space](https://arxiv.org/abs/1903.00870)|随机扰动优化生成近似样本，并以MH或IS恢复原后验；子空间加速与小噪声实验已有|RTO和非线性优化proposal不能改名为新机制；有限Newton分支密度须另核验|
|[Schillings等：Laplace approximation and noise-level robustness](https://arxiv.org/abs/1901.03958)|Laplace为重要性/准蒙特卡洛基准可在文中条件下对集中后验保持噪声鲁棒性|不能预设小噪声必然使成熟Gaussian基线失败；新响应图须与强Laplace比较|
|[Spantini等：Goal-oriented optimal approximations](https://arxiv.org/abs/1607.01881)|线性Gaussian任务后验均值与协方差的低秩风险最优近似|RRR/任务降维与平方风险的条件均值原则属于已有工具|
|[Au等：Manifold lifting](https://arxiv.org/abs/2003.03950)|通过增广流形和受约束HMC研究小噪声采样效率|未实现的完整lifting基线须标NOT_RUN，不能由旧bracket失败推断它失败|

证据层次为作者页面摘要与元数据。对应原理的本轮具体实现和数值验证另见A/B结果；这些摘要本身不验证R5的Maxwell性能。R5其他2026预印本引用保留来源身份，但本轮没有逐一重新核验其全部正文，不将其当作新检索成果。

响应换元行列式、Gaussian截断矩、RRR、条件期望和Neumann稳定性工具不作为原创定理。可能值得继续的贡献限于：实际Maxwell分层响应的可计算概率与失败边界；确定性任务特征在真实风险—成本比较中的净收益；以及真正通过独立评价的物理特征增量。未完成的新确认不能由开发集显著性替代。
