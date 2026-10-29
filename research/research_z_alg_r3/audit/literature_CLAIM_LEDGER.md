# R3 聚焦先例证据账本

检索日：2026-10-29 UTC / 2026-10-30 Asia/Singapore。主问题：R3的接收对偶硬数据约束变量投影，与已有CSI/SOM/WRI和变量投影究竟相差什么？

|ID/主张|来源与稳定标识|证据定位与等级|关系、置信度和限制|
|---|---|---|---|
|L1：contrast-source状态/数据lifting已有历史|van den Berg、Kleinman (1997), A contrast source inversion method, [DOI](https://doi.org/10.1088/0266-5611/13/6/013), S2 15a610c3f005b6979d624e55eb408a1cc40b902b|本轮S2身份验证C；冻结R2出版社摘要B。本轮出版社robots拒绝|历史背景，高；不以此未取全文证明全部零空间公式|
|L2：确定/模糊current子空间与共同材料已有直接结构|Chen、Zhong、Agarwal (2010), Subspace Methods for Solving Electromagnetic Inverse Scattering Problems, [DOI](https://doi.org/10.4310/MAA.2010.v17.n4.a6), S2 cf2adb06a3d58279724fc2f09b14fe68b6d6a072|冻结R2 audit/novelty.md已定位作者全文pp416–418式(28)–(34)为A；本轮出版社PDF首页与S2为C，页内打开/截图超时|直接重合风险，高；完整零空间与截断子空间的数值差别保留|
|L3：线性变量精确消元后优化非线性参数是变量投影|Golub、Pereyra (1973), The Differentiation of Pseudo-Inverses and Nonlinear Least Squares Problems Whose Variables Separate, [DOI](https://doi.org/10.1137/0710036), S2 5d38e06821c7f88d27b5b28039383c10a5c00fe3|SIAM身份C；出版社提交的Crossref原论文摘要B直接给出伪逆消元后的目标；本轮矩阵最小二乘关系另为自身推导|核心机制非新，高；S2年份1972与期刊1973冲突，以SIAM卷年为准|
|L4：场/材料联合搜索与方程放松可改变优化景观|van Leeuwen、Herrmann (2013), Mitigating local minima in full-waveform inversion by expanding the search space, [DOI](https://doi.org/10.1093/gji/ggt258), S2 df16b0fc59c3aaabbb54b4a5aaea8639953516bd|出版社摘要B：数据项+方程罚项，交替场与材料；作者2013全文链接可定位|“放松状态逃离局部极小”概念已有，高；不证明本轮相同目标或保证成功|
|L5：波动反演状态变量投影及包络梯度直接先例|van Leeuwen、Herrmann (2016卷年/2015在线), A penalty method for PDE-constrained optimization in inverse problems, [DOI](https://doi.org/10.1088/0266-5611/32/1/015007), S2 fe5601e99ae3138e6bd7269d9f2ccd183f07fe39|作者全文§4式(18)–(21)、§5.1式(31)，TierA：内层状态LS、消元目标、包络梯度与QR|直接方法族重合，高；本轮硬数据/ΔF并非其标准固定源罚目标的逐式同一|
|L6：参数相关Gram度量及逆Gram极限已有直接先例|van Leeuwen、Yang (2025期刊/2024预印本), An Analysis of Constraint-Relaxation in PDE-Based Inverse Problems, [arXiv](https://arxiv.org/abs/2403.15292), [全文](https://arxiv.org/html/2403.15292v2), DOI 10.1088/1361-6420/adab89, S2 09f55503f0273f288aafe81b3cb4e0553431428a|TierA：§2.4 Corollary1式(10)，§2.5 Lemma2；官方arXiv提交/作者身份C|本轮ψ概念的新颖性被直接削弱，高；矩阵朝向及特定CSI约束未判为逐式同一。Crossref确认2025-01-28在线/2025-02-28卷期；arXiv v2记录2024，HTML额外2026日期异常保留|
|L7：小型波动逆问题多项式化/多个材料解搜索已有先例|Everett (1996), Homotopy, Polynomial Equations, Gross Boundary Data, and Small Helmholtz Systems, [DOI](https://doi.org/10.1006/jcph.1996.0070)；Everett (2026), Multiple solutions to geophysical inverse problems via tracking flights across a deeply incised model landscape, [DOI](https://doi.org/10.1093/gji/ggag432)|冻结R2 audit/novelty.md及retrieval_primary_metadata.json，出版社摘要B；本轮未再次抓取全文|限制A的一般首创措辞，中高；不推出特定Q6公式已被发表|

## 本轮直接推导与外文献的分界

W的转置严格变为PU=d、U=Δ(F+HU)，由R3自己定义直接得到；这是对结构等价的审查推理，不冒充历史论文的逐字定理。ψ的零点关系、梯度与极点反例来自R3理论；本账本只指出机制先例和不合格的首创表述。

## 实际搜索与停止

3个ScholarQA查询保存在literature_scholarqa_collect.json；snippet/paper全部6项429。没有以这些空返回支持“未见先例”。单次6标识批验证全部resolved，见literature_scholarqa_verify.json。随后针对已知论文的一手网页搜索补缺，并发现L6。未做广泛citation crawl、专利审查、SPECTER排序或完整Scideator分类；停止理由为已找到足以否定强概念新颖性的直接先例，继续广搜对当前有限论文定位的价值低。

失败：CSI出版社robots；2013出版社完整页重定向失败但摘要已读取；SOM全文初次返回首页与结构，指定段落/截图超时；L6 DOI页面访问失败，但arXiv全文与S2身份可核验。Crossref书目记录见literature_primary_metadata.json。所有引文仅支持表内范围。
