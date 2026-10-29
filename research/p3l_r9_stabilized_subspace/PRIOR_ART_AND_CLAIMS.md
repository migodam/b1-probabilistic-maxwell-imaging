# R9-S：限定先例与可主张范围

2026-10-29。结论是：**固定/动态子空间、current两层分解、Fourier表示、Gaussian坐标保持、条件likelihood积分、Schur/Woodbury、任务误差停止均已有明确先例或标准数学。R9可能有价值的增量必须是原冻结 Maxwell/prior/task合同下的可复核成本、稳定性或资格改善；名称组合不构成原创性。** 本审计不认证全球新颖性。

## 检索执行与证据等级

使用 scholarqa-research 的 Evidence QA模式，限定两查询：

1. `subspace optimization method twofold subspace optimization FFT electromagnetic inverse scattering`
2. `likelihood informed dimension reduction nonlinear inverse problems posterior approximation task error`

collect的2个snippet与2个paper操作全部429，退出1、零候选。一次batch verify及一次较后保守重试也429，均原样保存。没有把它当作“无先例”。随后只补取R8P/R7/Z已列出之论文的原arXiv/作者/出版社身份与正文；未扩展成泛化网页综述。SOM系列可读正文使用Z本地留存的作者综述；不在本研究包复制其PDF或全文。LIS与Zahm作者全文本次通过arXiv HTML直接读取。

`audit/theory/scholarqa_collect.json`、`scholarqa_verify.json`、`scholarqa_verify_retry.json`保留错误。`doi_metadata_verification.json`保存精确DOI身份补核，身份与内容支持分开。历史R7/Z的成功S2身份核验只作为已存在来源记录，不冒称本次请求成功。新文献发现覆盖不足，不能称最近全部先例或系统综述。

Tier A=实际作者/出版社全文支持；Tier B=摘要支持；Tier C=只核验身份；Tier D推断不能独立支持实质文献结论。以下使用忠实概要，无大段逐字转录。

## 最强重叠

|先例|直接支持的已知构件|本次证据与定位|对R9的具体约束|
|---|---|---|---|
|Chen、Zhong、Agarwal，2010，Subspace Methods for Solving Electromagnetic Inverse Scattering Problems|SOM从接收算子SVD分出确定current与待优化current，并联合满足场/构成关系；模糊部分可含小奇异值的辐射current，不只是严格receiver null。|A：本地原PDF文本§3.1，式28–34，pp416–418；[DOI](https://doi.org/10.4310/MAA.2010.v17.n4.a6)。|不能以“强/弱current分开”或“共享材料合法性”主张新概念；R9 latent/material弱方向也不是同一个current空间。|
|Zhong、Chen，2009 TSOM|用内部current-to-field映射的主方向进一步限制模糊current，使搜索更稳定。|原篇身份C；作者综述A，§3.3式46–47，pp420–421；[原篇DOI](https://doi.org/10.1088/0266-5611/25/8/085003)。|强弱两次分解及对内部物理重要性的筛选已有直接先例。不能把一个普通谱阈值命名为新的概率方向。|
|Zhong、Chen，2011 FFT-TSOM|用离散Fourier基替代昂贵的内部SVD，并用FFT构造current；降低基构造/作用成本是原方法明确目的。|作者综述A，§3.3 pp421–422；原篇身份由[作者目录](https://www.ece.nus.edu.sg/stfpage/elechenx/Publications.htm)核对，IEEE TAP59:914–927。|FFT是表示/算术优势，不自动保留材料prior或后验事件。只做完整换基没有新信息；截断需明确支持变化。|
|Cui等，2014 LIS|按prior精度比较likelihood敏感性；局部GN谱跨posterior平均形成global LIS；近似后验保留complement prior，提供Rao–Blackwell任务估计。|A：作者稿§3.2式7–10，§3.3式11；[全文](https://arxiv.org/html/1403.4680v2)。|R9从少数anchor构造的弱frame不是自动global LIS；其metric若是材料RMS，应与prior metric区分。|
|Zahm等，2022 certified reduction|以prior条件likelihood E₀[L|Pz]构造KL最优ridge近似；基于log-Sobolev条件和posterior下log-likelihood梯度二阶矩给误差控制，并讨论采样近似。|A：作者稿§2.1式4、8–9，§2.2、§3及§4与LIS比较；[全文](https://arxiv.org/html/1807.03712v4)。|条件积分及证书不是新概念。遗漏分支的样本均值、一个点的谱尾不能免费给适用的global二阶矩和误差界。|

**对LIS必须公平且准确。** Cui2014的具体likelihood近似是 L(Pz)（更一般有prior均值移位），同时保留complement prior；它不是将complement posterior冻结成点质量。Zahm明确区分它与 E₀[L|Pz]；后者对KL有最优性。R9比较时应注明究竟实现哪一个，不把“只固定弱变量、删除其prior”的更弱方法标成完整LIS。

目前证据也不支持“所有先例只找一幅图”的叙述；LIS直接面向posterior expectation。R8P指出的certified reduced-basis Bayesian expectation先例也应保留；本次没有重新扩张审核其全部误差理论。

## R9 claim ledger

|ID|有限主张|当前证据/性质|必须保留的反证门|
|---|---|---|---|
|K1|完整固定正交Gaussianized66坐标保持原prior|标准Gaussian旋转；THEORY_AUDIT推导|平移prior不能置零均值；旋转截断盒不能重设轴对齐支持；动态Q不能省Jacobian。|
|K2|固定弱frame与动态frame在同一GN目标下可量化局部限制损失|标准二次投影恒等式，不是新定理|全非线性/硬事件/全域概率不能由局部损失推出。|
|K3|非零host固定幅度响应可精确复用|R8P矩阵恒等式；已读R9实现保留O=S(P+D_hB)|点值浮点一致只属NUMERICAL_PHYSICS，须计264完整RHS setup及全部reduced work。|
|K4|任务中心化包络可比独立分子分母界更紧|R8P证明与同一一维冻结模型数据；给定信息类的robust积分|普通quad仍快得多；最紧不是Maxwell可实现类最紧。|
|C1|固定Q在原任务下能稳定地省总成本|**待R9实测**；B/C当前Full-J同付，Q刷新只增边际CPU|少端点、少覆盖或少验收不构成快；动态也可更少trial胜出。|
|C2|动态刷新只在可观察转动/模型失效时有用|**条件机制候选**；分析可用principal-angle、实际预测defect、局部遗漏下降量|平稳例dynamic不胜、转动例fixed不败都应记录；没有无条件刷新最优定理。|
|C3|真实低维task积分能在旧R6严格数值门下更便宜|**待数值reference与公平总成本**|reference未全域资格则不能报accuracy/speed win；同evaluator两因子对照区分物理与积分收益。|
|X1|Full66 posterior/连续Maxwell/exact physical fiber已认证|**不成立/未建立**|局部图册、有限路径、rank66和近等价候选均不够。|
|X2|fixed-Q、TSOM思想、Schur或conditional expectation首次用于这些问题|**不允许**|上述直接先例已覆盖广义构件；本次检索失败不提供新颖性证据。|

## 只有何种结果值得继续写成独立贡献

首先需要原合同的稳定正面结果：相同材料RMS、同合法性/竞争门、相同corrector与验收预算，固定Q保留全部目标成功且降低全费用；或动态刷新在有实测转动失效时恢复成功而代价可摊销。单个四维frame冻结本身只是一项普通数值策略。

概率上更强的增量是：在R6原prior、真实Full Maxwell、原事件和严容差下，建立完整支持的数值参考或严格证书，并通过同evaluator强普通积分对照验证总费用。物理复用成功而certification过贵是有意义的分开结论；不能由前者暗示后者。

以上是针对当前合同的研究判断和可证伪要求，最终GO/NO-GO由主Codex结合实际结果作出。没有新物理、训练或云调用；没有第三方全文再分发。
