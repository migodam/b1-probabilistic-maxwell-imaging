# MAX-STRUCT-R2：电流消元、源协变与接收零空间的邻近文献

检索日期：2026-11-02。此文件提供竞争性证据，最终理论与新颖性裁决由主线程负责。采用 scholarqa-research 的证据审计流程：先查指定机制，再核对原文及身份；没有运行科研实验。本次检索有明显限流，不能称穷尽搜索或系统综述。

## 直接结论

“对有限罚二次目标精确消去 current/state”属于变量投影与 reduced penalty 的已有结构；“消元后得到材料相关的数据残差度量”已有直接原文公式；“固定仿射 current 空间里的 profile 大于完整 profile”属于普通二次投影分解。可以形成有用的 **CIE 专门化、反例与成本审计**，不能仅据这些结构宣称新通用定理或新信息。

全局源协变更新与逐源 PR 的差异已在 R1 D3 得到有限诊断。保持完整源度量与全局矩阵内积能去掉对源坐标的算法依赖，但不会增加物理测量，也不自动改善材料重建。接收暗空间与共享材料合法空间不相同；SOM/TSOM 已把接收谱及内部传播/傅里叶约束用于 current 表示。

## 最接近的已核原文

|文献与稳定链接|本轮证据|直接覆盖的机制|不能外推的结论|
|---|---|---|---|
|van Leeuwen & Herrmann，*A penalty method for PDE-constrained optimization in inverse problems*，Inverse Problems32,015007(2016)，[DOI](https://doi.org/10.1088/0266-5611/32/1/015007)，[作者预印本](https://arxiv.org/abs/1504.02249)|A：`receipts/vanleeuwen2015_penalty.txt:339–407`，PDF pp7–8，§4 Eq13,19–22；S2身份成功|固定参数时精确最小化线性状态；envelope gradient；约化 Hessian 是联合 Hessian Schur 补；数据参与修改后的状态估计|论文研究一般 PDE penalty；不是本项目 CIE 双β迭代的共同下降证明，也不是消除所有局部极小值|
|van Leeuwen & Yang，arXiv2403.15292，**v1**题名 *A data-driven approach to PDE-constrained optimization in inverse problems*，[v1全文](https://arxiv.org/html/2403.15292v1)|A：`vanleeuwen_yang2024.txt:248–315`，§2.3 Thm1、§2.4 Cor1 Eq10；S2同ID现题名为 *An analysis of constraint-relaxation in PDE-based inverse problems*|直接给出消元后 `½tr[E*(I+ρ⁻¹G(c))⁻¹E]`，G 是伴随解的 Gram；representer theorem 与小数据空间求解；材料相关 residual metric 已有|v1正文以弱PDE、有限有界观测、可逆状态算子等为条件；其数据估计 Gram 的特定模型资格不能无条件移植到一般非同置有损 Maxwell|
|Aravkin, Drusvyatskiy & van Leeuwen，*Efficient quadratic penalization through the partial minimization technique*，[arXiv1606.02395v3](https://arxiv.org/abs/1606.02395)|A：`aravkin2017_partial.txt:113–170,275–370`；§II，Thm1与 inexact analysis；S2身份成功|partial minimization 的条件性分析；不精确内层求解所需精度及外层梯度误差；强调高罚参数并不意味着约化目标同样恶劣|原假设含凸光滑 data fit、平滑可逆 A(u) 等；不保证任意 CIE prior/box 下全局收敛；不得把精确LS对一PR步的优势当成新算法优势|
|Golub & Pereyra，*Separable nonlinear least squares: the variable projection method and its applications*，Inverse Problems19,R1–R26(2003)，[DOI](https://doi.org/10.1088/0266-5611/19/2/201)|C：S2及Crossref题名/作者/DOI核验；本轮出版社全文未取得|变量投影的经典综述身份；实质推导由下行已读 Osborne 原文支撑|不据未取得全文声称某一具体 CIE 梯度或成本结果已被该文逐式证明|
|M.R. Osborne，*Separable least squares, variable projection, and the Gauss–Newton algorithm*，ETNA28,1–15(2007)，[出版社原文](https://etna.ricam.oeaw.ac.at/vol.28.2007-2008/pp1-15.dir/pp1-15.pdf)|A：`osborne2007.txt:120–160`，PDF pp2–3 Eq1.3–1.5；出版社PDF首页核验题名作者年份|正交 projector 表示、profile残差加线性 nuisance mismatch 的平方分解；标准 VP 与 reduced Gauss–Newton 的关系|不能把一般正交分解重新命名为新 Maxwell 定理；该文的统计大样本结论不直接适用于本轮少通道|
|Zhong & Xu，*Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems*，J.Imaging5,27(2019)，[开放原文](https://doi.org/10.3390/jimaging5020027)|A：已有原始 XML/正文在 `../research_fbg_r3_experiments/audit/literature_cie2019.xml` 与 `_body.txt`，§2.2–2.4 Eq15；S2/Crossref身份成功|已有 CIE-I modified contrast、DPCS/APCS、Fourier TSOM、每轮一PR current step和LS material；已有不同β用于两块的relaxed方案|该原文并未在本轮查到精确 current profile作为新方法；但这只能说明所读§2实现不同，不能据此认证完整文献没有CIE-VP|
|Sun, Kooij, Jin & Yarovoy，*Cross-Correlated Contrast Source Inversion*，IEEE TAP65,2592–2603(2017)，[DOI](https://doi.org/10.1109/TAP.2017.2673758)，[作者全文](https://arxiv.org/abs/1906.10864)|A：`cc_csi2017.txt:214–280`，§3.1–3.2 Eq10–12；S2/Crossref身份成功|明确指出单独 data/state误差可能掩盖测量域的交叉失配；引入第三 cross-correlated误差及其优化|与 R2 的精确 profile/resolvent认证不是同一个公式；只能作为必须面对的相邻已知机制，不能声称它已提供本轮相同保证|
|Zhong & Chen，*Twofold subspace-based optimization method for solving inverse scattering problems*，Inverse Problems25,085003(2009)，[DOI](https://doi.org/10.1088/0266-5611/25/8/085003)|C：本轮Crossref身份核验；原文未取得|TSOM原始工作的存在；具体 current两层表示由已读2019及项目真实实现交叉支撑|不凭题名判断其每个 projector 公式或算法严格等价|
|Jbilou, Messaoudi & Sadok，*Global FOM and GMRES algorithms for matrix equations*，Appl.Num.Math.31,49–63(1999)，[DOI](https://doi.org/10.1016/S0168-9274(98)00094-4)|B：出版社摘要、Crossref身份；未取得全文|处理多右端矩阵方程的 global矩阵Krylov/projection已有直接先例|不是 nonlinear PR-CG同一算法；不能把这篇论文当成 R1 D3全局PR等变证明，后者主线程已给出独立代数论证|

证据等级：A为可定位完整原文，B为出版社/作者摘要，C仅出版身份。论文身份和它是否支撑某个主张分别记录，未读全文不升级为A。

## 与项目既有知识的交叉边界

- R1 `THEORY_DERIVATIONS.md` §6及 `BASELINE_AUDIT.md` D1/D3：精确块根谱、单PR非线性方向差商、全局源基等变已经属于既有项目结果。
- R1保存的 `research/delegated/theory/ATLAS_READ_EVIDENCE.md` 中 Atlas06 F-OPT-003 已完整给出 current变量投影/有限罚边界；F-OPT-004 指出端口消元的变化零空间及成本；F-OPT-009 保留共享材料合法性。Atlas10 F-BI-020–022 已给出 TSOM投影像、含噪强模态修正、current系数最小二乘消元。
- R3 实现 `src/cie_tsom_baseline.py` 保留原始Fourier系数坐标、PR+、非下降重置、充分迭代对照及共同corrector。不得以只跑一小步的该基线来隔离“变量投影”的固有收益。
- 因此最可信的新工作表述是：在一个公开冻结的CIE合同中，量化 current容量/精确消元/源坐标规则如何改变profile及物理验收，并提供不支持成像泛化的反例。是否具有足够论文贡献由主线程判断。

## 检索与核验收据

`receipts/scholarqa_collect.json` 保存三查询：variable projection contrast source inversion；separable nonlinear least squares electromagnetic inverse scattering contrast sources；global conjugate gradient multiple right hand sides source encoding inverse scattering。6个搜索操作全部HTTP429，exit1；没有据空检索结果宣称没有先例。首轮 `scholarqa_verify.json` 也为顶层429。

稍后一次降低范围的 `scholarqa_verify_retry.json` 成功解析 **6/6**：vanLeeuwen–Herrmann、arXiv2403.15292、arXiv1606.02395、CIE2019、CC-CSI2017、Golub–Pereyra2003。`crossref_verification.json` 另核对6篇正式DOI，含TSOM2009/global1999。Osborne由出版社PDF首页核验。vanLeeuwen–Herrmann S2年份2015为online/preprint，与卷期2016不同；Aravkin S2/arXiv首年2016，本轮读v3日期2017；2403.15292的v1/最新题名变化已明确，按同一arXiv身份去重。

检索中额外使用出版社、作者预印本及机构站点填补已识别论文的全文缺口；不把一般搜索排序当学术结论。原PDF与逐行文本供本地审计，主报告只发表自己的概括和定位，不应把第三方完整PDF随仓库公开发布。

本轮未运行多起点、sweep、反演或数值恒等式验证；没有更改R1或历史研究。全文下载日志及页码索引保留，下载/原文阅读时间不冒充物理计算成本。
