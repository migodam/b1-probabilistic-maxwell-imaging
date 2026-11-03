# R5：材料相关残差变换与冻结 metric 的文献证据

核对日期：2026-11-03。范围为证据提取，不作原创性裁决、不推导 R5 核心定理；未读取 GT、未运行反演。复用 R4 已核实账本，仅补查本轮问题。A=已定位作者/出版全文；B=摘要；C=身份；R4-A=继承 R4 的已读全文记录，本轮没有重新取得原文。

## 直接相关证据

|文献|级别、定位|可支持的有限陈述|不能推出|
|---|---|---|---|
|Abubakar & van den Berg (2000), *Three-Dimensional Inverse Scattering Applied to Cross-Well Induction Sensors*, TGRS 38(4),1669–1681，DOI [10.1109/36.851966](https://doi.org/10.1109/36.851966)|R4-A 的周围正文；本轮 [TU Delft身份页](https://research.tudelft.nl/en/publications/three-dimensional-inverse-scattering-applied-to-cross-well-induct/) C；[作者机构PDF](https://repository.tudelft.nl/file/File_e636219e-fe7f-4c5c-b870-558ee7dacd45) 本轮timeout|R4记录 p1671 Eq21为逐位置3×3 EBA近似逆；p1673 Eq30把预条件放入CSI罚目标；Eq33后说明随conductivity contrast变化逐轮重算。|**尚不能核实材料梯度是否微分预条件器**。本轮未获得公式图像；不可声称此文“漏掉dP”。空间块耦合与逐点三分量耦合也须分开。|
|van Leeuwen & Yang (2024 v1), *A data-driven approach to PDE-constrained optimization in inverse problems*, [arXiv:2403.15292v1](https://arxiv.org/abs/2403.15292v1)|A：[全文](https://arxiv.org/html/2403.15292v1) §2.4 Cor1/Eq10；§3.3.3 Eq22–24后正文|状态松弛后得到材料相关Gram加权残差；明确比较G(c)与data确定的固定Gtilde，并说明固定data metric无需对Gram求导，只需misfit的adjoint-state导数。|这是直接相邻的变量/固定metric算法结构，不是R5单频、非同置、有损Maxwell的同一算法或恢复保证。Gtilde由数据识别的附加条件不能省略。|
|Sun, Kooij, Jin & Yarovoy (2017), *Cross-Correlated Contrast Source Inversion*, DOI [10.1109/TAP.2017.2673758](https://doi.org/10.1109/TAP.2017.2673758)|A：[作者全文](https://arxiv.org/html/1906.10864) §3.2 Eq11–20|state error另映到测量空间，交替更新current与contrast；current更新固定上轮contrast，contrast更新使用state与cross-correlated两项并给出材料梯度。|不能把它简称为完整当前材料resolvent的同一读出；也不能把block-coordinate中固定contrast误写为材料更新故意丢弃某个dP项。|
|van Leeuwen & Herrmann (2016卷期/2015 online), *A penalty method for PDE-constrained optimization in inverse problems*, DOI [10.1088/0266-5611/32/1/015007](https://doi.org/10.1088/0266-5611/32/1/015007)|R4-A：`research/research_max_struct_r4/literature/PRIOR_ART.md` P8，§4 Eq13/19–22；本轮Crossref核身份|固定材料精确消去线性状态，profile gradient与Schur结构已有。|envelope theorem消除最优状态的隐导数，并不自动允许忽略任意显式材料权重依赖；具体等价仍交父线程逐式审查。|
|Zhong & Xu (2019), *Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems*, DOI [10.3390/jimaging5020027](https://doi.org/10.3390/jimaging5020027)|R4-A：`BLOCK_SELF_CIE_PRIOR_ART.md` §2，Eq3/5/13–15；本轮出版社429、PMC验证码，Crossref身份|CIE电流残差、材料变换R=βχ/(1+βχ)、归一化罚项和交替反演已有。|坐标/残差尺度改变不等于所有有限残差目标与轨迹相同。当前空间块、目标归一化及导数要单独比较。|
|Zhong, Lambert, Lesselier & Chen (2016), NIE，DOI [10.1109/TAP.2016.2535492](https://doi.org/10.1109/TAP.2016.2535492)|C：R4及本轮Crossref身份；未读原始全文|只作为原始NIE文献身份线索。|不能据题名/二手转述声称原文梯度如何处理dP。R4已读综述的公式证据应另标综述来源。|

## 本轮问题的可审查边界

- **已有直接邻近算法：** 材料相关残差权重、从松弛状态消元得到的variable metric，以及固定data metric及其简化求导，均有已读全文证据。因而“离开root后改变轨迹”或“冻结权重省掉导数”单独不能当作检索上的空白。
- **需要父线程数学裁决：** R5 的 P(m) 具体作用于什么残差、消元或变坐标前后固定了什么、哪些项是真正的显式材料依赖、各坐标是否一一映射及normalizer是否变化。本报告不把通用乘积规则当作R5新定理，也不根据文献名称判断等价。
- **保留两种不同冻结：** 由data预先构成且整个求解固定的metric，与每个outer iteration重算、inner solve固定的metric，合同不同。上述van Leeuwen–Yang已读段落直接支持前一种；不能把后者冒认为同一实现。
- **2000关键缺口仍在：** 不能判定其材料更新是exact derivative、冻结近似或利用特殊恒等式。R4全文访问缺口没有被本轮metadata刷新补上。

## 版本与检索收据

ScholarQA CLI作为首选学术transport，精确2查询、4操作全部429，零候选；第一批5篇batch verify全部resolved，补充其余4篇身份记录；不重试429的发现查询。无可用S2 MCP。依父任务授权，转向已经识别的作者/机构/出版社页面；Crossref仅用于身份，绝不替代内容证据。详情见 `receipts/scholarqa_collect.json`、`scholarqa_verify.json`、`scholarqa_verify_additional.json`、`crossref_verified.json`、`WEB_ACCESS.json`、`EVIDENCE_LEDGER.json`。

van Leeuwen–Yang按arXiv v1身份引用（2024-03-22），不采用HTML自动显示的2026日期。其HTML §3.3.3 Eq23/24显示的权重与§2.4逆权重表记不一致，本报告仅采用正文明确的variable/fixed metric与求导差别，**不转录案例Eq23/24为可靠实现公式**。R4已提示同一版本风险。不是穷尽检索；429绝非不存在先例的证据。最终原创性、算法等价与是否值得实验留主线程。

身份冲突必须保留：S2将arXiv:2403.15292归到后续题名 *An analysis of constraint-relaxation in PDE-based inverse problems*（DOI 10.1088/1361-6420/adab89），本报告内容证据仍严格绑定已读2024 v1题名/版本；S2将EISPY2D正式DOI记录年份返回2021，而Crossref及官方仓库给正式刊载2025，故不混用预印本与期刊年份。S2 bibliographic resolved不代表跨版本公式已核验。

## 补查：高contrast衰减/饱和是否已有讨论

父线程补充问题后做一次有限web补查（2个query；不扩大S2 discovery）。新增直接证据为Bevacqua & Isernia，*Smart Rewritings of the Basic Equations for Quantitative Non-Linear Inverse Scattering*，2020预印本；2021正式题名为 *Quantitative Non-Linear Inverse Scattering: A Wealth of Possibilities Through Smart Rewritings of the Basic Equations*，DOI [10.1109/OJAP.2021.3063248](https://doi.org/10.1109/OJAP.2021.3063248)，[作者PDF](https://arxiv.org/pdf/2012.15252)，**A**：§III.B p4 Eq13后明确指出β增大使辅助contrast R趋近support indicator，导致由R恢复实际contrast困难；§III.D pp5–6再次指出大β使state equation趋于identity、internal radiation operator消失。这是参数化饱和/材料读出恶化的强邻近证据。

但上述讨论**不是**已确认的“有限非零constitutive residual下，metric derivative有统一向高contrast压力”的同一梯度分解/定理。本轮没有充分来源确认父线程给出的负项公式或self复杂d版本已经被文献明确讨论，也没有取得normalization-induced blow-up的直接原文结论。状态为 **UNRESOLVED_AS_EXACT_PRECEDENT**，不能称首次发现；饱和病态、罚项衰减与完整耦合轨迹三者不能混为一谈。父线程推导的正确性/新颖性/科学意义留父线程审查。

补查综述的S2 verify返回429（`receipts/scholarqa_verify_review.json`），不把它列为S2已resolved；其DOI/身份由Crossref与原作者PDF核对（`receipts/crossref_review.json`），延续R4已核身份。此失败未重试。
