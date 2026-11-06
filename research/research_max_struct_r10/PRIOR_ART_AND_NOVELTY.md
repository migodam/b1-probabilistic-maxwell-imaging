# 已有工作与本轮增量边界

检索日期：2026-11-07（本地）。结论仅针对检索到的资料，不作全球首创保证。新数据确认尚未完成时，不把历史DEV正例写成已确立的新算法优势。

先执行ScholarQA CLI的两种查询：`inverse scattering contraction integral equation subspace Green operator nonlinear inverse`、`twofold subspace optimization inverse scattering Green operator`；四个Semantic Scholar请求均429，完整失败bundle保存在audit/literature。未重复高频重试。随后使用已连接的Undermind学术索引补充三个定向查询，并核验DOI；这不是完成了一次系统综述，也不是完整Scideator流程。2019 CIE全文另复用仓库已有原文摘录；本轮publisher重新读取遭429。

|有限主张|已有证据|级别与限制|对R10的要求|
|---|---|---|---|
|材料块预条件与改变contrast更新已有先例|van den Berg, van Broekhoven & Abubakar, *Extended contrast source inversion* (1999), [DOI](https://doi.org/10.1088/0266-5611/15/5/315)：摘要明确preconditioned CG contrast update、TV及乘法正则化|Tier B，摘要可见，未取得全文|不能把“改材料步/加预条件”本身当创新；R10保留原prior，不能把prior变化算作结构增益|
|Two-fold subspace optimization已有明确先例|Zhong & Chen, *Twofold subspace-based optimization method for solving inverse scattering problems* (2009), [DOI](https://doi.org/10.1088/0266-5611/25/8/085003)|本次只核实Tier C书目信息，未取得全文；其算子定义另继承项目已冻结的旧Two-fold实现|不宣称发现Two-fold；具体低秩分母是否被该文或后续论文覆盖仍需全文对照|
|CIE和FFT-TSOM、不同current/material β已用于3D强散射|Zhong & Xu, *Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems* (2019), [DOI](https://doi.org/10.3390/jimaging5020027)|Tier A：仓库`research_fbg_r3_experiments/audit/literature_cie2019_body.txt`全文正文；Tier B：本次摘要；原文强调CIE与FFT-TSOM、relaxed双β|成熟CIE/Dual是主对照；Global-Dual只是本项目global PR内部消融，不是新文献基线|
|多轮Fourier contrast-source展开已用于CIE稳定化和成本降低|Xu, Zhang & Wei, *Fourier Bases-Expansion Contraction Integral Equation for Inversion Highly Nonlinear Inverse Scattering Problem* (2020), [DOI](https://doi.org/10.1109/TMTT.2020.2978165)|Tier B摘要、已核实DOI，未取得全文|R10继承FFT/APCS，不能把这些已有表示收益归给A+；另列actual affine floor与相同current的材料消融|

## R10实际做了什么

R10改变的是固定current材料残差：从局部one-hop形式改为 `β(I+βD_m−D_mK)⁻¹e`，K由背景G的低秩投影构造；全部材料自由度保留。新实现包含分母导数，以小Woodbury系统计算。主CIE current、FFT/APCS、初始化、prior和原Maxwell最终选择继承不变。

这是一种**具体的、待独立确认的材料metric候选**。Woodbury、Möbius、投影、预条件和同零集等价性均不是本轮原创定理。低秩矩阵分母是否具有可发表的算法增量，要由匹配current的实验、同成本结果及更完整的closest-prior全文对照判断。

目前历史DEV的关键实证是：Green19和Two-fold28都能保留b20d恢复，并使40ac从原A+的NRMSE约.626降到约.055；原始Gii保留控制也成功。因Green19已足够，这些数据**不支持Two-fold接收信息独有增益**，也不支持“原Gii有实现错误”或“目标变凸”。

论文可讨论受限场景中可重复的增量，不要求所有测试对象都获胜。但需要完整报告普通场景退化、共同失败、cost和物理离散边界；不能把单一开发对象、更多预算、额外信息或换prior包装成普遍突破。
