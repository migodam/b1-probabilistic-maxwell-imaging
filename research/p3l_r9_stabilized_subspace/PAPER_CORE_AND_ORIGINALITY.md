# 概率体积成像：论文核心、已有构件与可证伪的新机制

2026-10-30。独立理论综合；不读取E5真值、候选结果或评分，不新增物理调用。E5实施依据为已批准的 [最终运行协议](config/possibility_volume_protocol.json) 与 [最终统计协议](config/possibility_statistics_frozen.json)，DRAFT仅保留设计历史。本文件讨论论文的科学对象和可支持的贡献层级，不预判主实验结果；完成状态与结论须来自正式封存证据。最终科学方向与新颖性判断由主线程负责。

## 1. 论文首先回答哪个问题

给定一份有噪声散射观测，输出一组**原材料模型内、通过真实非线性物理检查、具有可解释规模的可能材料**。核心问题是：在有限计算费用下，这组可能性是否比单个anchor更经常包含GT附近的材料，或覆盖GT对应的任务值？若能做到，它已经具有成像价值，不必先完成准确posterior重建。

这里的“概率”首先指对声明的材料/噪声总体作重复实验时的覆盖频率；“体积”首先指有限合法材料云的邻域并集。它们分别不是单次观测下的可信后验质量，也不是逐像素区间任意拼接成的576维合法材料体积。两类更强主张可以以后独立升级。

不应把论文动机写成“过去只能输出一张CNN图”。多解反演、贝叶斯成像、随机优化与ensemble方法已有长期研究；单点CNN只是可能的部署基线之一，并非整个领域的能力边界。当前真正需要解决的是：**在相同物理模型和费用约束下，如何得到更多有用、彼此有区别、能够覆盖GT/任务的合法可能性。**

一个目前可用但不预置成功的题目是：**“利用固定传播算子探索非线性材料可能性：覆盖、范围与计算费用”**。若E5支持具体优势，再把其明确机制写入标题；若不支持，则保留方法限制与负结果，不能用“概率体积”重命名失败。

## 2. 连贯方法链，而不是四项不相干的工具

**固定物理承载层。** 背景、采集与ROI固定时，把宿主消去，得到

$$L_Dh=DE_R,\quad L_D=I-DC,\quad F=F_h+Oh.$$

$C,E_R,O,F_h$固定，变化材料只进入对角$D$。完整256格ROI保留所有原非均匀Full66材料与多重散射。其作用是让每个新增可能性更便宜地被同一物理模型验证；它没有自动减少材料不确定性，也不赋予候选概率。精确消元和40个reduced RHS的完整Jacobian结构均为标准代数，见F-OP-008—010。

**合法材料探索层。** 提案从decoder的latent空间出发；固定/动态弱方向组织多方向、多半径探索，完整物理Jacobian用于受限corrector，真实材料半径与signed progress保留多样性。B/C释放弱幅度及其符号，但限制步在弱补空间加当前一条weak直线中；最终合同中的历史 `*_RAY` 名称不表示强制正向射线。R允许更大的corrector空间。它们生成的是可验证的材料点，不是自由接收零空间电流，也不是未经校正的Born图。

**覆盖与范围层。** 输出有限云及其材料/任务邻域，报告GT覆盖、oracle容量、云直径、合格率和总费用。离线GT只用于封存后评分；在线候选选择只看真实dataPhi、合法性与几何门。更多候选不应被自动解释成更可靠：同等覆盖下越小的云、同等范围下越高的覆盖、同等效果下越低的费用才构成清晰增益。

**可选权重与概率升级层。** 均匀或loss软权重可以用于排序与任务动作，明确称候选评分或ensemble权重。要称原posterior，必须另外处理proposal、corrector推送、多对一吸引域、拒绝、去重及原prior质量。这个升级不阻塞前三层的科学检验，也不能从前三层的成功中被默认推出。

## 3. 最接近的先例：仅作窄范围核验

本次沿已有ScholarQA流程，仅增加一个查询：`randomize then optimize randomized maximum a posteriori nonlinear inverse problems`。snippet与paper检索均返回429；随后3个DOI的batch verify同样429。错误保存在 `audit/theory/scholarqa_volume_collect.json`、`scholarqa_volume_verify.json`，没有当作“没有先例”。依据主线程明确授权，改用已知论文的出版社/作者机构/arXiv原始页面核验。没有扩张citation graph、训练模型或启动另一轮广泛综述。

下表除SOM已有全文审阅外，新增加的三篇以**原始摘要Tier B**支持具体陈述；并未取得其完整定理条件的逐式审阅资格。DOI与作者、年份已从出版社/作者机构记录核实。此处是有限证据审计，不是完整Scideator分类复现，也不输出全球新颖性标签。

|工作与已核验出处|直接支持的事实|对当前论文的约束|
|---|---|---|
|Chen、Zhong、Agarwal，2010，*Subspace Methods for Solving Electromagnetic Inverse Scattering Problems*，[DOI](https://doi.org/10.4310/MAA.2010.v17.n4.a6)；已读作者全文§3.1、3.3，Tier A|固定接收SVD组织current，进一步按内部作用筛选或采用Fourier表示，并与场/构成关系联合处理。|固定子空间、两层筛选、current合法性不是新概念。当前material/latent候选空间与其current空间不同，但不同变量名不构成创新。|
|Bardsley、Solonen、Haario、Laine，2014，*Randomize-Then-Optimize: A Method for Sampling from Posterior Distributions in Nonlinear Inverse Problems*，[出版社](https://epubs.siam.org/doi/10.1137/140964023)，SIAM JSC36:A1895–A1910；Tier B|通过随机优化产生候选；线性Gaussian情形可直接得到posterior，非线性情形则推导proposal密度，并配合MH或importance sampling。|“随机提案再物理优化得到多幅图”已经有直接先例。候选产生与posterior校正的分层也不是新概念。|
|Wang、Bui-Thanh、Ghattas，2018，*A Randomized Maximum A Posteriori Method for Posterior Sampling of High Dimensional Nonlinear Bayesian Inverse Problems*，[出版社](https://epubs.siam.org/doi/10.1137/16M1060625)，[作者稿](https://arxiv.org/abs/1602.03658)；Tier B|rMAP面向大规模非线性逆问题的近似posterior样本；特定随机数据和prior均值下归入RML，讨论PDE约束优化和近似Metropolization降低偏差。|不应声称先例必须先做到精确posterior才输出多可能性，也不能把当前未经校正的云称为首次放宽这一要求。|
|Iglesias、Law、Stuart，2013，*Ensemble Kalman methods for inverse problems*，[出版社](https://doi.org/10.1088/0266-5611/29/4/045001)，Inverse Problems29:045001；Tier B|使用迭代ensemble Kalman处理逆问题，估计位于初始ensemble张成空间，并与该空间内的最优逼近和最小二乘比较；文中将其作为无导数优化方法。|ensemble与初始固定空间均有先例。一个ensemble存在并不等于已经建立其posterior资格；低维空间可能形成覆盖上限。|

RML的原始历史本文不另作未经全文核验的优先权归属；使用rMAP作者对其关系的明确说明即可。没有比较实现的RTO/rMAP/EKI，E5就不能声称对这些方法有经验优势。

### F-PAPER-001 · 随机优化为何早已能生成多可能性

以下是独立推导的标准线性Gaussian例子，说明随机化**目标**与只随机化**起点**不同。令 $y=Ax+\eta$，$\eta\sim\mathcal N(0,\Gamma)$、$x\sim\mathcal N(m,C_0)$，其中两协方差正定。独立抽 $\epsilon\sim\mathcal N(0,\Gamma)$、$\xi\sim\mathcal N(0,C_0)$，解

$$x_* = \arg\min_x\left\{\tfrac12\|Ax-y-\epsilon\|_{\Gamma^{-1}}^2+\tfrac12\|x-m-\xi\|_{C_0^{-1}}^2\right\}.$$

令 $H=A^T\Gamma^{-1}A+C_0^{-1}$。则

$$x_*=H^{-1}\{A^T\Gamma^{-1}(y+\epsilon)+C_0^{-1}(m+\xi)\},$$
$$\mathbb E x_*=H^{-1}(A^T\Gamma^{-1}y+C_0^{-1}m),\qquad \operatorname{Cov}(x_*)=H^{-1}.$$

均值与协方差正是线性Gaussian posterior。非线性目标的argmin推送一般不是同一个分布；多分支、局部极小值及Jacobian因子均需处理。这与RTO/rMAP摘要中明确区分的资格相符，但上述公式不是声称复写了全部非线性RTO算法。

相反，在固定严格凸目标中只改变初值，精确优化仍返回同一解。当前E5通过材料半径/signed-progress改变可行探索集合来维持多样性，属于另一种候选生成策略。它是否更适合有限费用下的覆盖，是要实证的问题，而不是“输出多幅图”本身已经证明的优势。

## 4. 哪些是标准构件，哪些才可能形成独立贡献

|内容|当前资格|成为可发表贡献还缺什么|
|---|---|---|
|固定host/ROI消元、共享源/receiver adjoint、完整Jacobian|标准精确离散代数；重要的可复用工程基线|正确实现和实测全费用足以形成工程结果，但不能称新散射定理|
|固定/动态弱空间、随机方向、多半径、GN校正|已有数值思想的明确组合|具体机制在相同模型和费用下稳定提高合法覆盖，且有负对照解释为何|
|候选云而非单点输出|长期存在的集合/ensemble目的|不同于已有随机优化方法的清楚目标和优势；不能以名称取代比较|
|coverage、width、best-of-K、真值隔离、全费用|标准统计/评估原则，当前合同可提高证据可信度|强且可复核的结果；更好的审计本身不自动成为新算法|
|任务弱敏感但高任务变化的提案分配|本文F-PAPER-005给出局部推导，**未实施**|超出标准广义Rayleigh quotient的独立机制、最接近先例核验、与当前强基线的裁决|
|带剩余漏覆盖控制的任务驱动停止|目前只有条件概率分解，没有全域漏失上界|可计算而非oracle的遗漏界，并验证其紧度与成本；否则不能升级为覆盖证书|

可能真正独立的机制，必须明确改变“哪些候选被产生、何时停止、如何保证不漏掉重要材料/任务”中的至少一个环节，并产生现有方法无法直接给出的可检验结果。例如：用可计算的任务残余而非局部谱尾分配候选，并在控制集合范围的同时减少实际GT/task漏覆盖。单独把SOM、Schur、弱空间和ensemble相连，不足以宣称这一点已经做到。

**当前允许的贡献形式是条件性的：** 本工作在原冻结非线性材料合同下建立可复用物理评估与合法候选探索，评估有限预算中的覆盖—范围—费用取舍；是否有算法性新机制，由E5结果、机制证据及后续最接近方法比较决定。没有这些证据时，可以是一项诚实的边界研究，不能提前称突破。

## 5. E5究竟能证伪什么

最终冻结E5的32个新parent、固定K16非anchor槽位、四个半径、固定/动态weak8及R随机控制，使下述具体假设面临失败。此处描述已冻结合同，不读取当前执行结果。

|有限假设|支持它所需的现象|直接反证或削弱|E5不能推出的更强结论|
|---|---|---|---|
|候选云提高表示容量|同一parent中合格云相对合格anchor降低oracle距离、提高预定ε覆盖|非anchor大多不合格、重复或离GT更远；增益只在极宽ε出现|部署时一定能选中oracle候选|
|固定框架足够稳定|B在同K与全费用下保持C的合格覆盖/范围表现|C显著恢复B遗漏的有效区域；B常发生几何失败|所有材料/噪声下无需动态更新|
|弱方向有实际覆盖价值|B/C在相近云宽与费用下优于R|R同等或更好；弱云较窄却持续漏GT；优势仅来自更大半径|已优于RTO/RML、SOM或其他ensemble方法|
|候选可支持任务动作|明确的GT不可见动作在预定任务上有收益，或任务区间以可接受宽度覆盖GT任务|oracle提升但均匀/软权重动作恶化；任务range变宽无覆盖增益|这些权重是posterior概率|
|额外物理计算有价值|all-attempt成本下的增益有可接受大小|花费很多候选求解而coverage和任务无改善|固定ROI消元导致了所有arm之间的差异|

全部臂共享固定ROI solver，所以B/C/R差异不能归因于这一公共承载层。物理核的单独效率证据与E5候选策略证据应分开。B/C与R的corrector几何不同，故对照评估的是整体策略，不是仅改变proposal分布的纯机制消融。

32个parent给出有限总体证据；方向、源、半径和计时重复不是额外独立parent。即使32/32命中，一侧95%二项下界仍约0.911，只能按冻结的总体/ε范围表述。严格posterior门不应阻塞这个结论，但也不能用“经验”一词抹去小样本区间。

E5可能证明固定弱策略失败而随机云仍有效。这将支持“可能性体积值得输出”，同时否定“固定弱8是合适生成机制”。也可能全部云失效；那需要保留主张被反证的事实，不能转而只展示少数好看的材料图。

## 6. 覆盖、任务和加权动作之间的可推导关系

### F-PAPER-002 · 材料覆盖如何传给平滑任务

令 $T:\mathcal M\to\mathbb R^p$ 为任务，$W\succeq0$ 给出任务metric，$z(m)=W^{1/2}T(m)$。假设在相应材料域有

$$\|z(m)-z(m')\|\le L_T d_M(m,m').$$

若合格云存在 $m_j$ 满足 $d_M(m_j,M)\le\varepsilon$，则任务点云存在一个点满足 $\|z_j-z(M)\|\le L_T\varepsilon$。这是材料coverage蕴含任务coverage，反向通常不成立：低维任务相同的材料可能相距很远。

原ROI平均实材料 $T(m)=q^{-1}\sum_{i\in ROI}\operatorname{Re}m_i$，$q=256$。由Cauchy–Schwarz和全网格RMS分母24，

$$|T(m)-T(m')|\le\frac{24}{\sqrt{256}}d_M(m,m')=1.5d_M(m,m').$$

平均soft occupancy $q^{-1}\sum_i\operatorname{sigmoid}((\operatorname{Re}m_i-t)/\tau)$ 的导数上界为 $1/(4\tau)$，因此 $L_T\le3/(8\tau)$。在 $\tau=0.02$ 时上界为18.75；对 $\varepsilon=0.075$ 得1.40625，已经超过任务[0,1]范围，说明这一全局界无信息。不能用材料coverage自动宣称该任务的精确覆盖；应直接评估预定task。

硬阈值不全局Lipschitz。若最近材料候选的距离≤ε，令GT在阈值t附近δ带内的ROI像素比例为

$$a_\delta(M)=q^{-1}\#\{i:|\operatorname{Re}M_i-t|\le\delta\}.$$

则两个材料的阈值占用比例差满足

$$|T_{hard}(m_j)-T_{hard}(M)|\le a_\delta(M)+\frac{576\varepsilon^2}{q\delta^2}.$$

证明：GT远离阈值超过δ而发生翻转的每个像素至少贡献δ²材料误差；总平方误差不超过576ε²。这个margin界说明硬事件的难点在哪里，但 $a_\delta(M)$ 依赖GT，只能用于离线解释，不能免费作为在线证书。

### F-PAPER-003 · 覆盖不等于加权均值准确：一个明确不等式

对同一非空合法云，令任务值 $z_j=W^{1/2}T(m_j)$、非负权重 $w_j$ 满足 $\sum_jw_j=1$，定义

$$\bar z=\sum_jw_jz_j,\qquad V_w=\sum_jw_j\|z_j-\bar z\|^2
=\tfrac12\sum_{i,j}w_iw_j\|z_i-z_j\|^2.$$

设材料ε邻域内候选的索引集合为 $G_\varepsilon=\{j:d_M(m_j,M)\le\varepsilon\}$，其总权重为 $a=\sum_{j\in G_\varepsilon}w_j>0$。由F-PAPER-002，近GT组的加权任务均值 $\bar z_G$ 距离真实任务至多 $L_T\varepsilon$。分组方差分解给出

$$V_w\ge a(1-a)\|\bar z_G-\bar z_{G^c}\|^2
=\frac{a}{1-a}\|\bar z_G-\bar z\|^2\quad(0<a<1).$$

结合三角不等式得到

$$\boxed{\|\bar z-z(M)\|\le L_T\varepsilon+\sqrt{\frac{1-a}{a}V_w}.}$$

当 $a=1$ 时第二项为零。该式解释了两个独立条件：云里要有接近GT的候选，而且权重不能把这些候选几乎清零；同时其余候选的任务分散度要可控。只满足第一项的“best-of-K很准”，不足以保证ensemble平均好。

均匀权重下，只知道至少有一个近GT候选时 $a\ge1/K_y$，于是得到保守界 $L_T\varepsilon+\sqrt{(K_y-1)V_w}$。增加K可改善coverage，却使这个最坏情形动作界变差。这不是矛盾：扩大集合容量与得到可靠单一动作是不同目标。

这个条件不等式是标准方差分解与三角不等式的直接结果，**不作为新数学定理主张**。完整稳定公式登记在 [第07章的F-VOL-015](formula_collection/07_MATERIAL_POSSIBILITY_VOLUME.md)。它的用途是要求论文同时报告coverage、任务dispersion和动作风险，并禁止把oracle容量误写为可部署收益。$a$含GT，在线未知；$V_w$在线可算。因此该式本身不是在线可信风险证书，也不改变冻结E5评分合同。

### F-PAPER-004 · 从有限coverage到总体任务风险的条件桥梁

若任务值所在集合的直径至多 $B_T$，所有实际候选权重在所考察规则下有统一下界 $w_{min}>0$，且云非空，则在coverage事件 $H=\{D_K\le\varepsilon\}$ 上，F-PAPER-003中的 $a\ge w_{min}$。对任意 $\gamma>0$，由 $(u+v)^2\le(1+\gamma)u^2+(1+\gamma^{-1})v^2$，

$$\mathbb E\|\bar z-z(M)\|^2
\le(1+\gamma)L_T^2\varepsilon^2
+(1+\gamma^{-1})\frac{1-w_{min}}{w_{min}}\mathbb E[V_w\mathbf1_H]
+B_T^2\Pr(H^c).$$

均匀有限云有 $w_{min}\ge1/K$；指数loss权重虽对有限实数严格为正，却可能极小，不能凭正性声称有有用的统一下界。若qualified云为空，需按冻结规则报告失败或使用明确fallback；该公式不能静默丢掉空云parent。

与F-VOL-013结合，在噪声资格正确时 $\Pr(H^c)\le\alpha+\beta_{K,\varepsilon}$。因此剩余任务风险可由“材料近似误差、加权任务分散、漏覆盖”三项解释。只有最后两项也有可计算/实测上界时，才得到有用结论；没有把posterior误差当作必经中间量。

注意 $\bar z=\sum_jw_jW^{1/2}T(m_j)$ 是任务动作。一般 $T(\sum_jw_jm_j)\ne\sum_jw_jT(m_j)$，且材料均值可能不在decoder像内，二者不能混写。

## 7. 一个可用来判断未来机制的局部推导：NOT_IMPLEMENTED

### F-PAPER-005 · 数据弱方向不必是任务重要方向

在anchor附近，用real latent位移 $p$。令 $J$ 为数据Jacobian，$B=\partial(W^{1/2}T(\mathcal G(x)))/\partial x$ 为任务Jacobian，取 $H=J^TJ+\lambda I\succ0$，其中 $\lambda$ 是明确的prior/正则尺度。局部问题

$$\max_{p^THp\le b^2}\ \|Bp\|^2$$

的最优值与方向由

$$b^2\lambda_{max}(H^{-1/2}B^TBH^{-1/2}),\qquad B^TBv=\mu Hv$$

决定。证明是令 $u=H^{1/2}p$，应用Rayleigh quotient。若要限制到当前材料弱空间 $p=Q_w a$，则把两侧分别替换为 $Q_w^TB^TBQ_w$ 与 $Q_w^THQ_w$。

这说明任务覆盖不应只看最小数据奇异值。某方向可能很难被观测识别，却对所关心任务几乎没有影响；反之，弱而任务变化大的方向更值得花一个候选槽位。任务权重W应由预先声明的单位/容差决定，不能看GT后调节。

广义特征值推导本身是标准目标导向降维思想，不构成新机制。可能的后续算法是，在真实物理相容与云范围约束下，为不同任务方向分配有限候选，并利用已生成候选的任务重复性减少无效槽位。若提出此算法，必须另行核验最接近的目标导向实验设计/降维先例，且与当前B/C/R强基线比较；当前E5没有实现它，也不因此追加实验。

局部任务rank≤任务数，不能据此删掉其余material方向：未被列出的任务、非线性曲率、GT邻域和硬事件仍可能依赖那些方向。高任务spread也不保证围绕GT分散；它可能仅扩大错误区域。最终仍需真实资格、GT/task coverage及width共同裁决。

## 8. 论文如何写才不会把目标与资格混为一谈

建议正文按因果链组织：先说明同一观测对应多种材料的实际任务需求；再给原材料、噪声与输出集合定义；随后推导固定算子复用和合法候选生成；接着给coverage—width—cost与task-action关系；最后用完整E5结果检验假设，并明确哪些优势没有出现。

严格posterior推断可置于独立扩展部分：它回答候选如何获得原测度质量，不是当前合法云能否有用的先决问题。反过来，合法云成功也不能代替posterior资格。近似权重、经验覆盖、条件概率证书应有不同名称与图例。

可写的论文主张应使用条件句或结果驱动措辞，例如：“在所声明的Full66材料/噪声总体与有限费用下，我们评估了利用固定传播算子的合法材料云，检验其是否改善GT邻域与任务覆盖。”待数据完成后，再填入实测效果、区间与失败，不先填“显著”“稳健”或“高概率包含GT”。

如果数据支持，最强的连贯贡献将是：一种明确的、计算可承受的合法材料可能性生成机制，在有限费用和集合范围约束下改善可复核的GT/task覆盖，并清楚区分集合容量与部署动作。若只证明普通精确ROI核更快，则贡献应收缩到算子复用；若候选云改善而weak机制不优，则保留云目的、放弃weak优势。任何一项都不需要伪装成准确posterior，亦不应把标准组合冒称新理论。

## 证据与范围记录

- 新增文献范围：RTO2014、rMAP2018、EKI2013三篇；新检索一个query，ScholarQA collect与verify均429，已保留。出版社/作者机构/arXiv核验用于有限摘要级主张；arXiv不存在的`1602.03658v3`页面曾返回404，随后核对真实版本为v1，没有掩盖访问失败。
- 已有SOM/TSOM/LIS/Zahm证据沿用 `PRIOR_ART_AND_CLAIMS.md`；未另行复制第三方全文。
- 数学：F-PAPER-001—005均注明假设；标准推导不主张原创定理，未来机制未实施。
- 实验：本文件不读E5 GT、不读E5结果、不新增RHS、不运行小型数值测试、不改变冻结合同。
- 新颖性：只能指出与已核验先例的重叠和可能的增量；本文件不认证全球新颖性，也不替主线程作最终科学裁决。
