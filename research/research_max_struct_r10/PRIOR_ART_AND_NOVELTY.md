# 已有工作与本轮增量边界

检索日期：2026-11-07（本地）。结论仅针对检索到的资料，不作全球首创保证。历史DEV、24对象主确认、后续四对象确认及曝光后的结构诊断分别列账。

先执行ScholarQA CLI的两种查询：`inverse scattering contraction integral equation subspace Green operator nonlinear inverse`、`twofold subspace optimization inverse scattering Green operator`；四个Semantic Scholar请求均429，完整失败bundle保存在audit/literature。未重复高频重试。随后使用已连接的Undermind学术索引补充三个定向查询，并核验DOI；这不是完成了一次系统综述，也不是完整Scideator流程。2019 CIE全文另复用仓库已有原文摘录；本轮publisher重新读取遭429。

|有限主张|已有证据|级别与限制|对R10的要求|
|---|---|---|---|
|材料块预条件与改变contrast更新已有先例|van den Berg, van Broekhoven & Abubakar, *Extended contrast source inversion* (1999), [DOI](https://doi.org/10.1088/0266-5611/15/5/315)：摘要明确preconditioned CG contrast update、TV及乘法正则化|Tier B，摘要可见，未取得全文|不能把“改材料步/加预条件”本身当创新；R10保留原prior，不能把prior变化算作结构增益|
|Two-fold subspace optimization已有明确先例|Zhong & Chen, *Twofold subspace-based optimization method for solving inverse scattering problems* (2009), [DOI](https://doi.org/10.1088/0266-5611/25/8/085003)|追加取得出版社摘要（Tier B）：利用内部current-to-field算子谱性质，先限制较低维current得到初值，再增大current空间；全文PDF未取得，Undermind亦标无PDF|不宣称发现Two-fold；本轮固定current表示、改变材料分母，与摘要所述current-space continuation角色不同；完整公式是否覆盖该构造仍需全文对照|
|CIE和FFT-TSOM、不同current/material β已用于3D强散射|Zhong & Xu, *Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems* (2019), [DOI](https://doi.org/10.3390/jimaging5020027)|Tier A：仓库`research_fbg_r3_experiments/audit/literature_cie2019_body.txt`全文正文；Tier B：本次摘要；原文强调CIE与FFT-TSOM、relaxed双β|成熟CIE/Dual是主对照；Global-Dual只是本项目global PR内部消融，不是新文献基线|
|多轮Fourier contrast-source展开已用于CIE稳定化和成本降低|Xu, Zhang & Wei, *Fourier Bases-Expansion Contraction Integral Equation for Inversion Highly Nonlinear Inverse Scattering Problem* (2020), [DOI](https://doi.org/10.1109/TMTT.2020.2978165)|Tier B摘要、已核实DOI，未取得全文|R10继承FFT/APCS，不能把这些已有表示收益归给A+；另列actual affine floor与相同current的材料消融|

## R10实际做了什么

R10改变的是固定current材料残差：从局部one-hop形式改为 `β(I+βD_m−D_mK)⁻¹e`，K由背景G的低秩投影构造；全部材料自由度保留。新实现包含分母导数，以小Woodbury系统计算。主CIE current、FFT/APCS、初始化、prior和原Maxwell最终选择继承不变。

这是一种**已实际独立确认、结果有正有负的材料metric候选**。Woodbury、Möbius、投影、预条件和同零集等价性均不是本轮原创定理。24材料中新增一个高对比环的恢复，同时损害桥和开缝壳路径；Green19平均NRMSE .17669对CIE .16287、完整CPU增加约19%。该证据支持条件性路径差异，不支持默认替换CIE。可发表的独立增量仍需更完整的closest-prior全文对照和可预测适用条件。

目前历史DEV的关键实证是：Green19和Two-fold28都能保留b20d恢复，并使40ac从原A+的NRMSE约.626降到约.055；原始Gii保留控制也成功。因Green19已足够，这些数据**不支持Two-fold接收信息独有增益**，也不支持“原Gii有实现错误”或“目标变凸”。

论文可讨论受限场景中可重复的增量，不要求所有测试对象都获胜。但需要完整报告普通场景退化、共同失败、cost和物理离散边界；不能把单一开发对象、更多预算、额外信息或换prior包装成普遍突破。

追加的分母资格推导使用标准辐射Gram和Rellich论证。Bessel积分式来自[NIST DLMF 10.9.1](https://dlmf.nist.gov/10.9.E1)；Rellich原理的经典来源见[Colton, MSRI 47 (2003), Theorem 3.1](https://library.slmath.org/books/Book47/files/colton.pdf)。本轮给出具体盘脉冲自项与非对角项匹配后，双侧投影材料分母在精确解析有限模型中的可逆性推论；不宣称首次发现耗散算子理论，也不把它当浮点区间证书、条件数界或材料唯一性定理。

按用户要求未新增NIE反演，是本轮资源与主对照选择；不据此宣称NIE在所有实现、参数和任务上都被CIE数学支配。
