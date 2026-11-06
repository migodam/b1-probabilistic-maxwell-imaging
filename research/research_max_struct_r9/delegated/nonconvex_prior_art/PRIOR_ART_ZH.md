# R9：非凸性干预与“交叉熵”的定向证据

2026-11-06。检索到的直接先例支持**改变搜索路径、改变松弛目标、或在特定有限表示下得到凸子问题**；没有证据允许把它们概括为任意高对比 Maxwell 材料反演的全局凸化。A/A+是否优于成熟CIE，仍须由同对象、原物理终点、材料指标和完整成本决定。本次只整理证据，不运行反演、不读取测试GT、不作最终新颖性判断。

ScholarQA Evidence QA为主流程：两种查询的4个snippet/paper操作均429，失败原件保留；随后引用身份核验**8/8成功**，7个期刊DOI另由Crossref核验。按先前授权使用作者/机构/出版社补全文。选8篇、读6篇全文，另2篇仅出版社摘要；不是系统综述。证据、版本、定位、限制见[EVIDENCE_LEDGER.json](EVIDENCE_LEDGER.json)，不把429当成没有先例。

## 已有负例对新建议的约束

已读R2 `src/profile_current.py`及R2/R3报告。R2是固定β=1的完整current精确消元，data/state用总能量标量归一化；有Cholesky内层与完整包络导数。它不是普通LSIE的罚参数延拓。原盲测profile成功1/3、native/global CIE各2/3、direct0/3；弯条全部失败。其约化式是材料相关的残差度量，不能称“消去非线性”。扩大current空间也不保证材料更准。

R3已量化含噪DPCS硬锁在GT处的加权本构代价.0289064，噪声模型期望.0226618、99%上界.0424996；但clean-CIE和强模态释放仍未恢复，11个公开开发控制均未通过原联合材料指标。因此**精确消元、释放强模态、去噪、单纯增加APCS容量不能作为尚未做过的新发明**。一次普通LSIE penalty continuation或同数据homotopy仍可作为有区别的已知控制；是否值得继续由父线程判断。

## 八篇核验文献与适用范围

|证据|原文支持的机制|不能推出什么|
|---|---|---|
|P1 van Leeuwen–Herrmann，2016卷期/2015预印本，[A penalty method for PDE-constrained optimization in inverse problems](https://doi.org/10.1088/0266-5611/32/1/015007)，**A全文**|二次PDE罚法、精确状态消元；§6实验把罚权逐级增加并用上级解初始化；分析目标是给定容差的Lagrangian驻点。|不是任意材料全局最优定理；罚权选择、延拓和增加的总成本有明确限制。R2固定CIE profile只覆盖其一个近邻控制，未覆盖所有延拓。|
|P2 van Leeuwen–Yang，2025期刊，[An analysis of constraint-relaxation in PDE-based inverse problems](https://doi.org/10.1088/1361-6420/adab89)，本次读**A：题名不同的arXiv2403.15292v1**|消元后是inverse-Gram残差度量；§2.7在特殊有限解空间和仿射系数下给出二次材料目标。|一般情况不保证凸性、可辨识性或该小空间逼近真实PDE；见下方原公式交接。不得混称读过最终期刊全文。|
|P3 Sun等，2017，[Cross-Correlated Contrast Source Inversion](https://doi.org/10.1109/TAP.2017.2673758)，[作者全文](https://arxiv.org/abs/1906.10864)，**A**|§3 Eq11–13在data/state误差之外加入测量域cross-correlated误差，耦合本构预测和端口数据；2D TM/TE与CSI、MR-CSI比较。|不是已证明优于A/A+；新误差用背景传播的一步本构预测，不等于当前材料的完整Maxwell响应。无全局凸化证明。|
|P4 Bao–Liu，2003，[Numerical Solution of Inverse Scattering Problems with Multi-experimental Limited Aperture Data](https://doi.org/10.1137/S1064827502409705)，**B摘要**|多实验正则同伦，随后沿波数作recursive linearization；应用合同明确使用多频、有限孔径far-field数据。|摘要中的global-minimum措辞不能脱离原假设；不能把新增低频观测或频率路径当成现有单频Y的无信息变换。本次未读其证明。|
|P5 Rubinstein，1999，[The Cross-Entropy Method for Combinatorial and Continuous Optimization](https://doi.org/10.1023/A:1010091220143)，**B摘要**|KL辅助优化、重要性采样、分布/温度调整用于非凸问题的随机搜索。|凸有界定义域不等于凸目标；摘要的组合Markov-chain收敛不能套到有限样本576像素连续Maxwell。没有检索到本合同下低成本全局恢复证据。|
|P6 Cuturi，2013，[Sinkhorn Distances](https://arxiv.org/abs/1306.0895)，[作者全文](https://www.marcocuturi.net/Papers/cuturi13sinkhorn.pdf)，**A**|熵正则使固定边缘分布下的transport-plan问题严格凸，并可用Sinkhorn缩放高效求解。|这里凸的变量是运输矩阵，不是非线性forward里的材料。该文是histogram/分类，不是Maxwell反演。|
|P7 Qiu，2021期刊/2020预印本，[Analysis of seismic inversion with optimal transportation and softplus encoding](https://doi.org/10.1088/1361-6420/ac1511)，[作者全文](https://arxiv.org/abs/2009.04138)，**A**|softplus＋归一化把实seismogram变成概率；讨论平移/伸缩的渐近凸性，数值验证Marmousi。|不能泛化为所有空间材料的凸性；编码/窗口可破坏平移凸性，大斜率有梯度/溢出问题；单频复场的相位与幅度需另作合格编码。|
|P8 Cruz等，2021，[Adding Prior Information in FWI through Relative Entropy](https://doi.org/10.3390/e23050599)，[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8152768/)，**A**|在L2数据项之外，用参考模型/井信息构造相对熵prior；数值比较不同形式。|是额外先验信息；原文一种平方熵形式有两个低谷且一个模型失败。它不是CE群体搜索，也没有普遍消除局部极小值。|

## 父线程优先阅读：P2条件凸化究竟在哪里

完整原文缓存`primary_texts/metric2024.html`；精确抽取`metric2024_section2_7.txt`及`metric2024_section3_3.txt`。更紧凑的LaTex交接见[METRIC_CONVEXITY_HANDOFF.md](METRIC_CONVEXITY_HANDOFF.md)。[原文§2.7直达](https://arxiv.org/html/2403.15292v1#S2.SS7)。

其一般松弛消元为

`Jρ(c)=½ Σj e_j(c)* (I+ρ⁻¹G(c))⁻¹ e_j(c)`。

§2.6的PDE残差投影空间`W_n(c)`依赖当前材料，不能只因本构对材料仿射就宣布一般目标凸。§2.7改选**整个解空间**`U=span{p_i}`，其中`p_i`为同组source/measurement泛函的实值、足够正则、线性无关Riesz表示；该段沿用前文无噪与well-posed设置。Corollary2的Eq13为

`J0(c)=½ tr[(M−A(c)M⁻¹D)* M⁻¹ (M−A(c)M⁻¹D)]`，

`Mij=<pi,pj>U`，`A(c)ij=overline(Ac(pj,pi))`。若`Ac`对材料系数c仿射，且M不依赖c，这是关于c的二次目标。**M是所选基的内积Gram，不能混作一般状态Gram G(c)**。Remark3进一步给出`A(c)=M D⁻¹ M`线性求解；它需方形、可逆D。这个特定Galerkin结构和直接线性fit不同于普通罚权延拓，但是否保真、是否可识别由父线程审查。

原文Remark2明确保留材料可辨识与物理近似两个缺口。我们的外部矩形48×6/12端口响应、576材料，并不自动构成上述匹配泛函小解空间。§3.3.2的1D Helmholtz数据Gram Eq21还需要频率导数`d′ij(k)`、边界响应`bi(k)=ui(1;k)`、其导数及边界声速；不能说单频Y直接给出所需Gram。

**原文版本问题**：v1总理论Eq10含inverse metric，但§3.3例子的Eq23–24印成无inverse的`I+ρ⁻¹G`；此处保留不一致，不代作者修公式。最终期刊题名/年份与v1不同，已由S2＋Crossref核验；真正采用公式前应核对最终版本。

## “交叉熵”至少要区分四种问题

1. **CE优化器**：改变材料候选的采样分布，仍可用原物理Phi排序。它能尝试跳出当前局部轨迹，但每个全波样本都要收费。有限样本成功不等于全局证书。P5支持这一类，不能用P6的凸性给它背书。
2. **CE/KL数据loss**：例如对归一化概率p、q写`−Σp log q`。即使它对自由q凸，`q=encoding(F(m))`的材料复合目标仍可能非凸。本次未核实某篇原论文证明它能凸化本项目单频复Maxwell；这项是概念辨析，**不是已验证性能结论**。必须说明复数如何编码、如何保留相位与总幅度、噪声似然为何允许换成该loss。
3. **熵正则OT/softmin**：P6的熵作用于运输计划；P7的编码用于波形几何。平滑一个内层问题与凸化材料外层是不同主张。归一化只比较质量分布时，总幅度需要单独保留；不能把丢信息后的低loss称为更准物理反演。
4. **相对熵prior**：P8把参考材料信息加进目标；额外井/参考模型必须计入信息合同。没有合法参考时，不能用GT、支撑或已知形状补上。原文自身有负例，熵不是通用免局部极小值工具。

## 可测试的简单机制：建议与已知结论分开

以下只是基于证据的候选控制，**不授权扩展实验，也不预告会改善**。

- 普通LSIE完整current penalty continuation：区别于R2固定CIE权重profile；固定少量罚权，最终回到原物理目标，和同额外CPU的原corrector配对。先看是否改变可达端点，再判断材料，而非只看松弛损失变小。P1是直接先例。
- 同数据空间homotopy：可以预先定义从可计算起点响应到原Y的路径；必须注明这是改变中间目标的启发式。P4的多频递归线性化不等于该单频公式；本次未确认二者算法同一性。终点仍必须原Y、原forward、原材料自由度。
- CC-CSI作为成熟耦合控制：若机制是“状态误差在端口被放大”，它是直接近邻；不能把普通cross-residual另命名成首创。增加的接收/Green动作和最终完整正演都计费。
- CE群体搜索若采用：分布维数、样本预算与更新规则先冻结；与等成本随机或multistart比较，保留失败样本费用。低维协同方向可以是proposal，但最终恢复必须释放全部材料像素；低维成功不能冒充全空间全局优化。

当前证据足以阻止“变量改写、精确消元、交叉熵或熵平滑自动全局凸化”的宽泛表述；不足以判断哪一种固定干预能修复R9具体失败，也不足以证明A/A+普遍胜过成熟CIE。那些判断留给父线程的公式资格与封存实验。

原始检索/核验/下载回执均在本目录。`primary_texts/`为研究阅读缓存，含第三方全文，不应直接作为可再分发自有成果发布；Utrecht403下载失败文件与成功作者版本分别保留。
