"""Assemble bounded scientific reports from closed R7 receipts, zero physics."""
from pathlib import Path
import json,time
import numpy as np
from physics import OUT,write


def put(name,text):
    p=OUT/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text.strip()+'\n')


def main():
    res=json.loads((OUT/'audit/results/METRICS.json').read_text())
    rr=json.loads((OUT/'RESOURCE_SUMMARY.json').read_text())
    mat=json.loads((OUT/'analysis/MATERIAL_PHYSICS.json').read_text())
    source=json.loads((OUT/'analysis/SOURCE_RECOVERY_AND_EARLIEST_CONFLICT.json').read_text())
    flags=dict(GO_HIGHER_ORDER='NO_GO',GO_PHYSICAL_FIBER='NOT_ESTABLISHED',
        GO_FINITE_NOISE_AMBIGUITY='SCOPED_GO',FULL66_POSTERIOR='UNQUALIFIED',
        TASK_IDENTIFIABILITY='LOCAL_OR_EMPIRICAL_ONLY',GO_TAP='REDIRECT')
    compact=dict(campaign='P3L-R7',status='COMPLETED_E0_E1_E2',gates=flags,
      independent_decoder_parents=256,development_contexts=4,locked_new_parents=8,
      zero_physics_rank66_counts={'relative1e-6':88,'relative1e-8':163,'relative1e-10':186,'relative1e-12':191},
      material_radii=[.0025,.005,.01,.02],endpoint_statuses_locked=128,
      at_material_RMS_02={'A_frozen':{'competitive':8,'near_pair_le1':6},
        'B_fullJ':{'competitive':8,'near_pair_le1':7},'C_jets':{'competitive':6,'unavailable':2,'near_pair_le1':6},
        'D_random':{'competitive':1,'near_pair_le1':0}},
      sampled_paths={'B_chord':{'qualified':8,'total':8},'C_raw_cubic':{'qualified':6,'unavailable':2,'total':8},'nodes':[0,.25,.5,.75,1],'continuous_certificate':False},
      higher_order_comparison=res['paired_bootstrap'],resource=rr,
      source_groups=dict(A=[0,1,2,3],B=[4,5,6,7],A_nearpair_B_revealed=5,total=8,new_layout=False),
      qualifications=['finite-noise likelihood-competitive does not imply credible mass','pair whitened distance<=1 is separate from DeltaPhi<=2',
        'C fixes first-order velocity; only second/third coefficients projected; no exact profile or cubic-null claim',
        'E1 equal-material witnesses are numerical and stratum-specific; no universal Full66 rescaling quotient',
        'four controlled copula laws are not actual conditional prior or posterior samples'],created=time.time())
    write(OUT/'RESULTS_SUMMARY.json',compact)
    put('START_HERE.md','''# P3L-R7：从参数 gauge 到真实材料弱响应

本轮已完成零物理 E0/E1 与真实 E2，并封账 **37,456 / 50,000 complex RHS**。原 R2–R6 数据与结论未改写。

最清楚的结果是：材料图确实可以在 0.02 RMS 尺度发生变化而仍竞争，但普通冻结/Full-J 正确器已经覆盖主要收益；本轮高阶机制没有胜出。精确 Maxwell material fiber 和 Full66 posterior 都没有建立。

|对象|本轮结论|
|---|---|
|Full66 decoder gauge|分层依赖；不存在普遍一维尺度 gauge。|
|八个新对象的 0.02 RMS 材料候选|Full-J 8/8 likelihood-competitive，其中7/8另满足均值响应 pair 距离≤1。|
|高阶二/三阶预测器|6/8竞争，另2例局部多项式路径不可行；未显示相对Full-J独立收益。|
|路径|Full-J八条五节点均过竞争门；不认证节点间连续纤维。|
|概率与论文|Full66 posterior UNQUALIFIED；强TAP主张 REDIRECT。|

先读[完整中文报告](RESEARCH_REPORT_ZH.md)，再读[逐项裁决](GO_NO_GO.md)、[理论与证明](THEORY_AND_PROOFS.md)、[长期研究背景](RESEARCH_MEMORY_R7.md)。实际数组、封存预测、失败与中央账本见[证据包](EVIDENCE_PACKAGE.md)，运行入口见[复现说明](REPRODUCE.md)。下一轮从[GPT Pro交接](GPT_PRO_HANDOFF/START_HERE.md)进入。
''')
    put('E1_FULL66_GAUGE_AUDIT.md','''# E1：Full66材料商与参数退化

使用原P1 decoder、Gaussianized66、完整576格复材料实Hilbert度量。256个独立prior点在计算前固定，另15个确定性边界/平台控制不计独立样本。新增Maxwell RHS、LU、forward和训练均为0。

|相对截断|1e−6|1e−8|1e−10|1e−12|
|---|---:|---:|---:|---:|
|数值rank66的对象数|88|163|186|191|

191个点的transition场矩阵在machine尺度满64列秩；65个点有数值构造的饱和不活跃方向。transition格数10–164，pure格数0–242。不能从一个cutoff把这些数字升格为精确代数秩或统一商维数。

每列三档中心差分的最大相对误差1.53e−8；独立SciPy gesvd重算全部秩一致，最大奇异值差1.14e−15。rank公式在1e−10处有一个数值判定不一致，保留为conditioning警告。完整逐点结果在E1_RANK_STRATA.json，初版英文记录与原始收据保留在audit。

四个无pure、transition行数分别19/40/10/14的对象092/201/233/251，可用inverse-quintic修正场与共同属性尺度，四档参数步长.001/.01/.1/.5全部合法、材料RMS差≤2.42e−16。它们是数值同材料见证；直线沿一阶null方向通常会弯离材料等值集。另一个无pure对象136在较紧阈值满66，不能作相同推广。

数学上，严格stratum且A_T满行秩、无pure、属性在内点时，可在局部构造等材料曲线；有pure时径向属性变化不能由该pure格的一阶场导数补偿。详见[条件证明](THEORY_AND_PROOFS.md)。因此**局部gauge存在，不等于Full66普遍可以降去一个尺度；单点DGv≈0，也不等于有限exact gauge**。

![真实256点材料谱与分层](figures/E1_FULL66_SPECTRA_AND_STRATA.png)

![有限null方向及接缝控制](figures/E1_FINITE_KERNEL_AND_SEAMS.png)
''')
    put('THEORY_AND_PROOFS.md',r'''# R7：数学对象、证明和未满足条件

以下把标准恒等式、条件理论与实际数值证据分开。复材料/观测均作为实Hilbert空间处理；当前结论只针对原TM24离散模型，不是连续Maxwell PDE认证。

## 1. 材料秩与真正的decoder gauge

写材料增量为δq(z)，δ≠0为复属性，q为实相分数，B=D_zq。属性位于内点且q≠0时，

\[
\operatorname{rank}_{\mathbb R}DG=\operatorname{rank}B+1+\mathbf1\{q\notin\operatorname{range}B\}.
\]

证明：旋转复相位令δ为正实数，场导数张成实子空间range(B)。两属性张成q、iq；iq独立于所有实列，q是否增加一维取决于其是否在range(B)。Gaussianized属性的normal-pdf缩放在有限内点严格为正，不改变精确秩，但会改变数值conditioning。host的q=B=0须另列rank0；pure的B=0、q≠0为rank2。

在严格饱和分层，固定属性、令T为transition格。若A_Tv=0，则沿z+tv所有transition场不变；只要其余格不跨接缝，材料严格不变。若q无pure格，A_T满行秩，给属性正尺度α，令transition目标分数为q_T/α，再用quintic在(-1,1)的光滑逆恢复场，并用A_T的右逆修正z。连续性保证足够小α下plateau不跨界、属性仍合法，因此得到局部同材料曲线。数值SVD与有限点等式只是该条件的诊断，不是外向舍入证明。

反之，DG满列秩时选取一个可逆66输出minor，逆函数定理给局部单射，故无非平凡局部同材料gauge；它不排除远处不连通的相同材料原像。原prior包含不同开放strata，所以不能将“某regular stratum满秩”改成整个prior几乎处处无gauge。

## 2. 材料归一化与物理弱性

令M=[Re DG; Im DG]/24。对保留的有效SVD写M=UΣVᵀ，取T=VΣ⁻¹、Q=MT。K=J_xT是单位真实材料RMS的白化响应。其谱才区分材料变化与参数缩放；ker(DG)本身不是Maxwell kernel。

R7主方向最小化wᵀ(KᵀK+TᵀT)w、||w||=1。这只是方向选择的prior成本意识，不向posterior添加新的material Gram prior。原目标始终Phi+||x||²/2。纯Maxwell最弱方向及多cutoff结果单独报告，不因巨大潜变量代价而偷称可行材料自由度。

对正则材料弧m'(0)非零，若DF在其材料切空间单射，则响应严格首阶为1。R7数值中的非零c1不满足exact-null前提；只能研究有限噪声弱响应，不能调用exact-null的高阶纤维存在结论。

## 3. Maxwell与decoder二阶项分别来自哪里

\[
D^2(F\circ G)[v,v]=D^2F[DGv,DGv]+DF\,D^2G[v,v].
\]

在当前材料度量下，可将D²G分成有效材料切空间投影和法向余量。切向部分依赖坐标选择，法向部分还依赖所声明的stratum/数值cutoff；它们不是全部Maxwell多次散射。分量向量可相消，不能把norm平方百分比当作数据影响的可加归因。

本轮保存完整decoder2、propagation2、b2、m2_tangent/normal等向量；没有重演历史F只保留norm的证据缺口。

## 4. 实际实现的二/三阶current递推

设L0=I−diag(m0)G_D，j0为原场，E0=e_inc+G_Dj0。用不带阶乘的Taylor系数m(h)=m0+Σhⁿm_n、j(h)=j0+Σhⁿj_n，收集各阶得到

\[
L_0j_n=\operatorname{diag}(m_n)E_0+\sum_{k=1}^{n-1}\operatorname{diag}(m_k)G_Dj_{n-k}.
\]

它是已知高阶散射导数的具体实现，每阶真实支付八个源RHS，复用同一LU。C固定第一阶材料速度||m1||M=1，用共同实强坐标U解z_n=−A⁺b_n，A=J_xU，n=2,3。所有补偿系数为实数，未使用复材料坐标的非法复最小二乘。

**本轮C并非严格response-aligned profile。** 固定prior-aware首速度后，Aᵀc1不要求为0，只投影二/三阶系数。原始c1中可由该top8强基投影的范数比例为.27%–17.59%，大部分对象很小，但不能假装其恒为0。严格profile的IFT理论还需处理z1及后续共同低阶系数。当前预测器和普通Full-J的竞争结论只针对所实现的C，不是否定所有高阶理论。

quintic全局只有C²。严格stratum内三阶表达成立；跨seam后的误差只能由实际有限响应检查支持。本轮全部六条可用C原始路径最终跨1–6个seam格，有限预测准确不自动给全区间Taylor定理。

## 5. 共享补偿冲突

在同一低阶材料路径和共同实强坐标下，第a组的本阶项为A_a z+b_a。令H_a=A_aᵀA_a、z_a=−A_a⁺b_a、c_a=b_a+A_az_a，zbar=(ΣH_a)⁻¹ΣH_az_a；在所需满列秩条件下，

\[
\min_z\sum_a\|A_az+b_a\|^2=\sum_a\|c_a\|^2+\sum_a\|z_a-\bar z\|_{H_a}^2.
\]

证明是正交分解A_a(z−z_a)+c_a及加权平方完成。秩亏时直接做stacked实最小二乘，不套用可逆H公式。本式是标准加权LS，不能命名为新Maxwell定理。R7一阶就已有正冲突，故二阶恒等式数值正确不证明高阶指标优于一阶Fisher/Schur。

## 6. 精确核、弱谱与谱稳定性

常秩映射的精确kernel可积；截断bottom-k弱谱分布未必可积。有限噪声near-equivalence也不是传递的等价关系。冻结rank不能当解流形维数。

在**固定材料frame**的自伴矩阵KᵀK上，若分隔谱簇的gap为g、扰动norm为η且η<g/2，可用Davis–Kahan型sin-angle界≤2η/g。动态prior切空间还需运输frame并控制材料投影变化；直接对齐SVD向量不足。实际R7用投影Fro范数、主角度、J变化和多cutoff报告，未取得便宜的全区间谱/余项证书。rank-changing seam附近不外推这些regular条件。

## 7. 材料任务的有限见证与两点界

ROI均值为线性材料函数，|Δmean|≤sqrt(576/256)||Δm||M=1.5||Δm||M。固定right-half ROI的二阶材料任务Q=mean|m−host|²，局部导数可直接计算；有限余项≤(576/n_right)||Δm||M²。阈值格占比/事件不享有普通光滑导数界。

固定两个材料世界后，Gaussian数据均值的白化距离d给KL=d²/2、TV=2Φ(d/2)−1。标准两点检验归约给指定标量任务的两世界最大平方风险下界Δψ²(1−TV)/8。R7保存这些数值，但只作**未来固定候选对的两点统计子实验**条件界，不是原prior Bayes风险，也不是用同一观测搜索后获得的posterior可信区间。

## 8. 同材料纤维与概率

观测仅经材料生成，给定材料后decoder纤维条件核仍由原prior决定。若一个局部同材料构造满足相分数p_i=p0/α_i、属性δ_i=α_iδ0，则独立组合后均值增量为(m0−host)EαE1/α。EαE1/α−1=E[(α−α')²/(2αα')]≥0，是标准Cauchy–Schwarz推论。

四个非均匀Full66受控等权离散纤维验证了此映射效应，并有128 RHS物理控制。等权四点不是实际prior条件核，不给Full66 posterior质量。generic满秩区域不能照搬这个尺度商；需要真正保留其联合推前测度。

## 9. 两种cubic的条件渐近

若材料prior有正原子w0、分支测度a(h)h^(α−1)dh，响应首项c h^p非零，任务首项b h^s，并且evidence和任务分别有支配/尾部控制且原子主导归一化，则h=σ^(1/p)t换元给分支质量O(σ^(α/p))、任务漏项O(σ^((α+s)/p))。

decoder seam的(p,s,α)=(3,3,1)与正则材料物理cubic的(3,1,1)不同；无原子时分支质量可恒为1。R7只在解析toy核对这些公式，独立积分误差4.92e−16，斜率约.333/.666/1.333；没有真正Maxwell三阶null弧见证，不能写成Full66概率率。

## 10. 先例与创新边界

[Arhab等高阶散射导数原稿](https://arxiv.org/abs/1803.00002)、[Cui等LIS](https://arxiv.org/abs/1403.4680)、[Diaconis等流形采样](https://arxiv.org/abs/1206.6913)、[Spantini等目标低秩](https://arxiv.org/abs/1607.01881)已覆盖相关广义工具。ScholarQA本轮collect六操作429保留，五篇身份核验成功；具体全文/摘要证据和题名版本差异见audit/literature。检索缺口不作为原创性依据。
''')
    put('E2_MAXWELL_FIBER_RESULTS.md','''# E2：完整非线性Maxwell有限材料实验

四个DEV为已暴露R6物理中心，八个新locked只由观测从zero进行Gaussianized direct-Probit64重建。原prior、sigma、TM24、8×32复通道不变。E0/E1新增RHS=0；整个R7支付37,456 RHS。两个locked anchor触及严格状态cap，保留非驻点/失败状态；几何研究不以MAP驻点为必要条件。

## 0.02 RMS主结果

成功同时要求原decoder合法、材料距离误差≤2%、signed progress≥一半、原noise-ball、DeltaPhi≤2。合理prior增量≤10为另列启发式，所有成功例也满足。

|路线|竞争成功/8|另满足pair距离≤1/8|说明|
|---|---:|---:|---|
|A Frozen|8|6|冻结J仍很强，不能省略这个简单基线。|
|B Dynamic Full-J|8|7|普通材料有效坐标GN与原decoder回缩。|
|C Maxwell jets2/3+同正确器|6|6|003/004的局部多项式路径不能在latent半径4内到达.02，计失败，不缩小距离补齐。|
|D Matched Random+同正确器|1|0|原先会放大近null坐标的随机实现已在DEV修正；本表为冻结后的强控制。|

B的pair范围约.173–1.121。003虽然DeltaPhi=−4.987且竞争，却不满足pair≤1；两种门不能混同。8个点的pair关系更不构成严格等价关系。

两种小距离.0025/.005所有四臂都8/8成功；.01时随机仅4/8，其他8/8。全部128个方法/距离状态和未选试探均保留。没有GT在线选样本或停机。

## C是否超过B？

两者在.02均有端点的n=6，C−B的Phi均值−.01154，parent配对描述95%区间[−.06022,.02133]。端点RHS相同，C另支付24 RHS setup/方向；四距离摊销后每个端点另6 RHS，冷单端点另24 RHS。不能因C有两例没跑端点而称它更快。

全8成功率差−.25，描述区间[−.625,0]。只有8个parent，不作总体优势/劣势的高精度统计；就冻结的竞争门与完整费用，**未建立高阶独立收益，GO_HIGHER_ORDER=NO_GO**。这个结论不是对所有高阶方法的数学否定。

两个DEV另给普通Full-J最多20态的参考；实际分别272/304 RHS，包括恢复中心。它们仍是budgeted profile，不宣称严格收敛。原始几何候选及本轮C的没有收益不归因于“Maxwell无解”。

## 路径资格

B的8条anchor到.02端点参数弦，在t=0/.25/.5/.75/1的五节点全部过noise-ball和竞争门；C的6条可用原始cubic路径五节点也过门，另2条不可用。这里只能称sampled-path-qualified。C的有限预测相对误差在六个端点约.86%–4.73%，但全部跨seam，且c1非零：不是精确cubic物理纤维。

![八对象真实材料，按固定ID展示](figures/LOCKED_MATERIAL_GALLERY.png)

![成功率与实际费用](audit/results/01_success_cost.png)

![五节点路径](audit/results/03_five_node_paths.png)

完整逐对象表见audit/results/ENDPOINTS.csv、SUMMARY.csv与E2各parent的RESULT.json；原始latent/material/prediction/current/field/J在Release无损证据内。Full66后验、模式质量、连续exact fiber均未认证。
''')
    put('TASK_OBSERVABILITY.md','''# 材料任务与照明条件可辨识性

任务在计算前固定：T1=ROI平均复材料，T2=right-half ROI的mean|chi−host|²，T3=ROI Re(chi)≥.75的格占比及其≥.25事件。T2是材料二阶函数，不是用GT定义的在线risk；T3依赖真实材料而非latent phase标签。

在B的八个.02 RMS端点上，平均实材料变化绝对值最大约.000632，平均虚材料最大约.000211；区域二阶任务变化最大.005014，阈值占比变化最大3/256。八个硬事件均未改变。它们说明这些具体弱材料弧对低阶任务影响小；不证明所有相容材料的任务都唯一，也不等于posterior方差或credible interval。

局部线性任务反演敏感性通过L K_g⁺计算，并分别保存A/B/联合三组。只有所声明的有效材料frame、rank和局部线性条件下可解释；不是Bayesian covariance。标准两点Gaussian KL/TV任务界也保留，明确用于固定候选对的未来两点子实验。

## 第二照明组是否揭示第一组的弱变化？

使用原8源中固定互斥A={0,1,2,3}、B={4,5,6,7}。A专属弱方向的线性B/A响应norm比为1.78–2.83。A-only预算正确器的.02端点，5/8满足pair_A≤1但pair_B>1，给出有限材料见证。所有源仍共享同一材料，保存的是完整八源预测。

但是这些增量主要已经出现在一阶谱：共同strong基的一阶补偿冲突占联合剩余能量约7.24%–16.21%，二阶约.46%–3.35%。二阶冲突恒等式误差约1e−15，不证明高阶指标有额外采集选择价值。没有得到“线性指标相当而高阶冲突独立预测更好”的确认。

加入B也增加数据量；本轮没有新旋转layout、没有等测量数的最优采集设计、没有新增观测降低Bayes风险的posterior计算。kernel交集关系是标准线性代数。

![材料广义谱与照明子组](audit/results/04_source_groups.png)

![共同材料补偿及有限响应](audit/results/05_group_compensation.png)

机器结果：analysis/SOURCE_RECOVERY_AND_EARLIEST_CONFLICT.json、analysis/MATERIAL_PHYSICS.json及audit/results/SOURCE_GROUP_ENDPOINTS.csv。任务资格为LOCAL_OR_EMPIRICAL_ONLY。
''')
    put('JACOBIAN_STABILITY.md','''# Jacobian、材料切向与rank变化

E1的数值rank66数量随cutoff由88变为191，说明不能把冻结latent Jacobian的rank当物理解集维数。真实plateau-null、近null与seam临界必须分开。E2八个重建中心在主cutoff下材料有效rank均66，最小材料响应奇异值约.0308–.0956；选择的prior-aware单位材料方向c1 norm约8.33–16.53，显然并非exact-null。

五节点B路径上，最大latent J相对变化约5.70%，材料切向投影Fro变化可达3.252。单个弱SVD向量可能旋转/换序，因此用完整投影和主角度比较；不直接凭向量符号对齐下结论。实际数据在analysis/MATERIAL_PHYSICS.json和每个NODE文件。

固定frame的谱簇若有正gap且扰动小于gap一半，可使用条件sin-angle界。动态材料像还需控制frame运输和seam，R7没有获得廉价全区间稳定证书。有限node检查不排除node间变化，也不证明非对合控制的新可达性。

材料二阶项的切向能量比例在八个点约45.3%–65.0%，剩余为当前有效frame的法向部分；这不是数据响应归因。decoder与Maxwell二阶数据向量cosine可到−.910，实际相消须保留向量，不能按norm百分比相加。

![材料与物理曲率、实际Jacobian变化](figures/CURVATURE_AND_J_STABILITY.png)
''')
    put('GO_NO_GO.md','''# R7逐项科研裁决

|判定|状态|证据与限制|
|---|---|---|
|GO_HIGHER_ORDER|NO_GO|冻结C未优于普通B，覆盖6/8对8/8，另有setup费用；只针对本实现，不否定所有高阶理论。|
|GO_PHYSICAL_FIBER|NOT_ESTABLISHED|没有连续exact Maxwell fiber证明，c1非零、仅有限节点。|
|GO_FINITE_NOISE_AMBIGUITY|SCOPED_GO|真实材料RMS.02、B8/8竞争且7/8 pair≤1，原prior/物理/噪声内成立。|
|FULL66_POSTERIOR|UNQUALIFIED|本轮没有试图用12个上下文认证全后验。|
|TASK_IDENTIFIABILITY|LOCAL_OR_EMPIRICAL_ONLY|固定任务与照明的局部导数及有限见证，不是credible interval/Bayes风险。|
|GO_TAP|REDIRECT|高阶导数、商测度、LS冲突均有先例；未建立优于强普通方法的Maxwell新机制。|

Full66存在stratum-specific decoder gauge与generic regular点无尺度gauge同时成立；没有普遍一维商降维。四个非均匀受控copula构造能推广R6的映射反例，但未求原Full66条件prior或posterior权重。

最值得下一步的是：控制一阶谱后再检验有限补偿冲突的独立用途，或先求局部gauge上的真实条件prior对材料任务的影响。不要扩大本C Explorer网络，不重启一轮大Flow，不将任务积分理论与本次候选数量合写为Bayesian算法成功。
''')
    put('FAILURE_ANALYSIS.md','''# R7失败、纠错与未运行内容

历史失败只读；新增原始失败版本和账本保留。全部物理动作数值失败为0，但这不代表没有开发错误、拒绝步或算法失败。

1. 初始新worktree引用本地旧main，未开始物理；立即切换明确的R6源177a9ea，记录R7T报告基线d1bdc8与公开namespace打包差异，不混用成一个实验。
2. 零物理GN检查发现null(material-progress)仍包含decoder kernel，可只降低prior而不移动材料；修成有效材料子空间后再投影。修复前原prior反例和验证保留。
3. DEV0随机控制在材料frame投影后lift，放大近null潜变量，造成假不可行。修成native随机位移在真实材料strong约束下的Euclidean投影；原DEV0及付费重跑都计费，不以弱随机baseline制造优势。
4. 独立审查发现原锚点64次计数未覆盖LU恢复，probe verdict未显式纳入有限jet门，距离回缩后未复查trust。全部在新locked生成前纠正。DEV共150步中4步轻微trust超界，仅DEV2/C，原结果保留、整DEV2付费重放；其他146步零物理复算符合相同约束。
5. 原64cap采用预留最终LU恢复槽的严格实现。新003、007触及cap且梯度仍非零；原样保留，不用GT再初始化，不叫充分MAP收敛。它们仍在noise-ball内，足以研究该合法几何中心的有限响应。
6. C在locked003/004的.02档局部多项式路径不可行，记失败。普通A/B可到同距离，故不能把它解释为原材料类不存在该距离的世界。
7. 六条C路径跨C² seam，不能从良好有限预测误差宣称全区间三阶正则或strict fiber。c1非零，未建立物理cubic-null机制。
8. 一次行政性分析用了无NumPy的system Python而失败，改用已有研究venv后完成，0物理费用；不影响任何封存预测。

主locked预测先封存、后单独GT评分；相区来自公开decoder状态，任务定义不用GT。无新网络训练、无posterior概率softmax、无MCMC/SMC资格升级、无新layout、无3D、无实验测量数据。原R5-B保留率HOLD和R6强化控制HOLD仍保持。
''')
    put('PRIOR_EVIDENCE_AND_FAILURES.md','''# 继承证据与本轮边界

历史固定报告、JSON和commit的逐项来源在[长期研究背景](RESEARCH_MEMORY_R7.md)及audit/history。R6已封账6,860,984 RHS不转入R7；Full66 posterior未认证；S3的V=QD完全补偿为同材料gauge；原Full比较下节时PASS不删除，强化Full-q和跨窗口HOLD也不升级。

R5-B的512新对象确认已完成，80%蒸馏收益保留门仍HOLD。R2/R4/R5材料风险正面结果、相对NLL正面结果与绝对posterior失败同时继承。F3的端点/中点资格、A17同current space无独立ROM精度增量，均未改写。

本轮只付费37,456 RHS，不训练新CNN/Flow，不重复历史512对象。理论包为47页报告与TeX、claim ledger及建议脚本，不是完整原始训练资产。真正复用输入来自已核验的本地原worktree；source/asset身份与归档范围分列。
''')
    put('RESEARCH_REPORT_ZH.md','''# P3L-R7：真正材料弱响应存在，高阶方法尚未胜过普通正确器

本轮完成了256个新prior点的零物理材料秩审计、四个已暴露开发上下文、八个全新测量驱动上下文，以及有限路径和独立数值审查。总新增 **37,456 complex RHS、1,325次state/LU**，50,000硬上限内封账。没有训练网络，也没有把有限候选当成posterior。

最重要的答案是：**去掉decoder参数重复后，原P1中仍能得到实际0.02材料RMS的有限噪声竞争材料；但普通冻结/动态Full-J已经做到，高阶二/三阶方法没有独立增量。** 精确材料solution fiber、Full66后验和强TAP新机制仍未建立。

## 1. 为什么必须先分清四个对象？

不同x表示同一材料，是decoder gauge；材料确实改变但一阶散射弱，是材料可观测性；不同材料在完整非线性forward下有限噪声相近，是候选几何；给这些材料分配质量并计算任务，是另一项概率问题。R6的均匀V=QD补偿只涉及第一项，不能直接证明第三项。本轮按E0→E1→E2执行，没有用第一个问题的正面结果替后面背书。

## 2. 固定物理、资产与独立性

保持原Gaussian64×Uniform2、精确probit Gaussian66、quintic平台、host、被动材料、TM24/VIE、256 ROI、8源×32复接收、sigma=.00119211798514519及source-major sqrt2/sigma白化。原R4 prior/target和A30 frozen solver只读复用。输入理论包CRC与BLAKE2b身份、47页PDF及TeX逐节核查；不把精简包当全部原始数据。

E1用固定2026102807；locked八对象用独立SeedSequence([2026102809,i,1])与noise子流2，未按秩或效果重抽。历史59项元数据、38,507个parent ID和初始化流去重审计有明确覆盖限制。DEV只复用R6 CDEV0000–0003；新数据生成与在线程序分进程，全部八对象620份预测/结果先封存，再读GT评分。代码在生成前冻结，后续未因GT改变方法。

## 3. E1：Full66不是统一尺度商，也不是处处无gauge

在相对阈值1e−6/1e−8/1e−10/1e−12下，rank66分别88/163/186/191个。191个transition场矩阵machine尺度满64列，65个存在饱和不活跃方向；这些数值不能当精确代数秩。

四个无pure且transition行可独立调整的点有弯曲同材料数值路径，四档参数步都合法，最大材料RMS残差2.42e−16。另一个无pure点并不满足这种降维。严格局部IFT解释了何种stratum允许尺度补偿；有pure格或材料Jacobian满66时不能普遍删除一个尺度。独立gesvd与三档全列差分完成，数值阈值敏感和一个rank公式不一致例保留。

![真实Full66材料谱和分层](figures/E1_FULL66_SPECTRA_AND_STRATA.png)

## 4. E2方法：强普通对照与真实高阶向量

主方向由材料正交frame上的KᵀK+prior代价选择，材料速度归一为1；不简单取最小latent奇异向量。四臂共用decoder距离匹配、signed progress、原Gaussian prior目标及真实Full-Maxwell评价。

- A：冻结J的合法有限步正确器。
- B：动态Full-J、有效材料坐标GN和原decoder回缩，是主要普通强基线。
- C：同一首速度、共同top8实strong坐标中的二/三阶current递推，加同一动态正确器。
- D：native随机位移投影、材料速度归一与有限材料RMS匹配的物理提案，加同一正确器。

每档最多三个新trial，明确叫budgeted corrector；另两例DEV给最多20态普通参考。每次新state8 RHS、首次receiver-J32 RHS、每阶八源jet8 RHS，所有拒绝和恢复均付费。首次实际探针200 RHS，伴随误差1.95e−15，chain/递推误差约1e−13，小有限材料响应预测误差2.2e−9/1.4e−8。完整decoder与Maxwell二阶向量均保存。

C固定首速度而只补偿高阶系数；其c1并不为零，也不自动满足严格response-aligned profile。这是本实现的明确范围，不能以它的NO-GO否定所有高阶方法。quintic只有C²，跨seam以后依靠实际有限验证，不能调用全局C³定理。

## 5. 八个新对象：真实不同材料成立，高阶优势不成立

|0.02 RMS路线|竞争成功/8|pair≤1且成功/8|
|---|---:|---:|
|A Frozen|8|6|
|B Full-J|8|7|
|C Jets|6|6|
|D Matched Random|1|0|

竞争门同时要求距离±2%、进度≥一半、合法、原noise-ball与DeltaPhi≤2。它不是credible region。B的003端点虽然更能拟合含噪数据，pair距离1.1206，不能叫两个均值的≤1近等价。C的003/004在局部latent半径4内不能完成.02，照记失败，不用较小距离填数。

在两方法都有端点的六对象上，C−B Phi均值−.01154，描述95%配对区间[−.06022,.02133]；C还要24 RHS setup，四档摊销每档另6 RHS。缺少端点导致费用较少不是加速。8 parent只支持机制小样本，不是总体算法冠军。

![全部新对象，固定ID顺序](figures/LOCKED_MATERIAL_GALLERY.png)

![实际成功率和费用](audit/results/01_success_cost.png)

## 6. 有限路径与物理阶数

B八条弦在五个固定t节点均竞争，C六条可用原始cubic路径也通过。它们只取得sampled-path资格，不证明所有节点之间相容，不证明连续exact fiber。C端点预测相对误差约.86%–4.73%，仍跨1–6个decoder接缝。

八个中心的材料有效rank都是66，最小材料响应奇异值约.0308–.0956；所选方向c1 norm约8.33–16.53。材料首速度为1，响应一阶非零，因此没有把decoder三次平坦冒称为Maxwell三阶null。条件渐近仅由明确toy验证，不写成Full66物理率。

![五节点完整响应](audit/results/03_five_node_paths.png)

## 7. 多照明：五个有限揭示见证，主要已是一阶效应

原八源分为A0–3、B4–7；没有新layout。A专属弱方向的B/A线性敏感度比约1.78–2.83。在A-only正确器得到的.02材料中，5/8对A的pair≤1、对B>1，说明单组弱材料变化可被另一组揭示。

共同strong基下的补偿冲突分解数值正确，最大相对defect约1e−15；但一阶已有约7.24%–16.21%冲突，二阶仅约.46%–3.35%。因此未证明高阶冲突优于普通一阶谱/Fisher。联合也增加了数据量，不能作为等预算最优采集定理。

## 8. 材料任务：图像变化不代表每个任务同样不确定

固定任务包括ROI复均值、right-half材料二阶函数及真实材料阈值事件。B的.02材料变化中，平均实/虚材料变化最大约.000632/.000211，区域二阶任务最大.005014，阈值占比最多变3/256，八个硬事件都不变。它们只是这些有限弱弧上的经验可辨识性；没有覆盖所有相容材料，更没有算posterior质量。

局部task dual、Gaussian两点KL/TV下界及其条件保存于analysis。不能把它们称作已获得Bayesian risk或credible interval。GT评分只作封存后几何诊断，默认成像动作仍是基础点，没有GT选候选或用Best-of-K冒充在线重建。

## 9. R6 copula现象能推广到哪里？

四个非均匀Full66同材料构造中，分别交叉四个场与四个属性，全部64个点合法；独立组合制造材料均值偏差/方差，EαE1/α公式误差<3.5e−17。额外128 RHS核验显示同材料端点响应差约1e−14，交叉端点可达约1.035白化L2。

这是原decoder的真实非均匀映射反例，却是受控等权离散纤维，不是原prior条件抽样或Full66 posterior。它说明在有gauge的stratum必须保留依赖；不说明所有66维点都能用S3商降维。

## 10. 科学与工程资格

共37,456 RHS，10,600 forward、26,856 derivative，1,325 state/LU；19个物理job，物理动作失败0。物理job总wall72.21秒、动作wall14.77秒、CPU67.47秒，进程峰值约190 MB。E1另7.18秒、490 MB；编辑、检索及交付时间不伪称被这个物理账本完整profile。没有GPU/网络训练。

三个前置代码问题和一次随机控制问题都在新locked生成前修复；DEV原结果与付费重放保留。新对象003/007触及anchor状态cap的记录不删除；科学负结果与数值失败分开。实际成本不能与R6不同优化配置的1.73秒直接换成加速比。

## 11. 七个最终问题

1. **Full66 gauge？** 在饱和/特定无pure低transition stratum有局部自由度；不存在普遍一维rescaling，数值近零不等于exact。
2. **真正材料弱响应？** 有。.02 RMS的B八例竞争，其中七例另满足pair≤1，非同材料重复。
3. **高阶比Full-J好？** 当前实现没有，GO_HIGHER_ORDER=NO_GO；不扩大其神经Explorer。
4. **另一组照明能恢复？** 五个有限见证；主要由一阶效应解释，不是新layout泛化。
5. **什么任务更稳定？** 这些弧上的均值/区域函数较稳定，阈值有限变化；仅局部/经验，未认证概率任务。
6. **R6能推广什么？** 非均匀特定gauge上的依赖破坏反例成立；通用尺度商、原prior条件权重和Full66 posterior未建立。
7. **强TAP证据够吗？** 尚不足。标准高阶散射、LIS、coarea和LS已有先例；本轮复杂C未胜普通方法，采集收益也未超过一阶解释。

## 12. 下一轮的最小可证伪问题

先控制一阶材料谱和同等测量费用，找补偿冲突不同而普通Full-J仍无法解释的有限材料区间，再讨论高阶采集指标；若找不到，停止独立高阶算法主张。另一条数学路线是求非均匀局部gauge上的**真实原prior条件核**及其材料任务积分，而不是继续等权候选加softmax。未来若使用更少通道创造exact kernel，必须作为新观测合同，不冒称本轮原采集合格。

本轮最终为：有限噪声材料歧义限定GO；exact physical fiber未建立；Full66 posterior未合格；高阶算法NO-GO；强TAP贡献REDIRECT。完整代码、向量、失败、统计及长期背景均可复核。
''')
    put('REPRODUCE.md','''# 可复现入口

本轮已封账，`CLOSURE.json`阻止在同一账本继续求解。**不要删除closure、清空ledger或覆盖raw来复跑。**

## 零物理复核

使用原Gaussian研究venv（NumPy/SciPy/Matplotlib/threadpoolctl），从独立R7 worktree执行：

```bash
PYTHONDONTWRITEBYTECODE=1 ../../.venv_nn/bin/python research/p3l_r7_material_solution_fibers/audit/results/audit_results.py
PYTHONDONTWRITEBYTECODE=1 ../../.venv_nn/bin/python research/p3l_r7_material_solution_fibers/src/theory_controls.py
```

第一条读取无损证据中的固定数组与中央账本，独立重算128状态、五节点、成本及封存身份，不读取GT。第二条仅是解析toy。数学/decoder QA脚本位于audit及src/validate_decoder_audit.py，历史原始输出不要覆盖，复核请写独立输出目录或使用包内复核收据。

## 原始执行顺序

E1冻结256点→decoder QA→E2物理probe→4DEV→保留错误并纠正→LOCKED_FINAL_FREEZE→generate8→parent逐个单worker→seal→单独offline score→注册DEV profile/gauge控制→close。真实命令和日志已保留。生成与online共用原forward，但在线call graph没有offline_truth读取；score是唯一正式GT评价入口。

R7的src/physics.py使用原R4/P1与A30 frozen vendor；原目标配置、geometry、background及可重建源码作为必要复用输入单独无损归档，见ASSET_MANIFEST。新物理复跑必须在独立目录与明确新预算下进行，不能误改本轮封账或历史R6。单物理worker、BLAS1、最多一个live LU、complex128/float64与原packing必须保持。

原始数组不放普通Git blob；Release完整证据恢复后才可用np.load。归档保留原文件bytes与相对路径，BLAKE2b/逐字节readback收据与代码身份随包提供；不做SHA-256附加检查。网页和Markdown只是阅读层，不能替代原始证据。
''')
    put('GPT_PRO_HANDOFF/START_HERE.md','''# R7给下一轮GPT Pro：有限歧义成立，复杂高阶没有胜出

先读[完整报告](../RESEARCH_REPORT_ZH.md)、[数学与证明](../THEORY_AND_PROOFS.md)、[裁决](../GO_NO_GO.md)和[长期记忆](../RESEARCH_MEMORY_R7.md)。本轮已真正运行37,456 RHS，不是计划。

## 可以继续推理的真实现象

1. Full66秩分层真实而强烈：256点中数值rank66随阈值88→191；四个无pure低transition点有弯曲同材料见证，generic regular点不支持统一尺度商。
2. 8新测量驱动中心均能用普通正确器产生.02 RMS竞争材料；7例均值pair≤1，五节点样本也竞争。它不是连续exact fiber或正概率质量证明。
3. 当前jets数值系数正确，六条原始路径有限响应预测误差低，但C仅6/8到.02、B为8/8。C固定首速度，仅投影二/三阶；c1非零。不能继续以网络复杂度救这项NO-GO。
4. 第一组照明的弱材料在另一组有五个有限揭示见证，但一阶谱已经解释主要作用；二阶冲突还没有独立预测价值。
5. 非均匀Full66的四个受控纤维中，独立场/属性组合制造错误材料与实际散射。它不是原prior条件核；下一步真正的概率难题正在这里。

## 必须保留的资格

普通冻结J也很强；C失败不代表物理材料不存在。Full66 posterior UNQUALIFIED，任务只有LOCAL_OR_EMPIRICAL，TAP仍REDIRECT。R6原窗口PASS与强化控制HOLD、旧R5-B80%HOLD均不改。新增网络训练0，当前8parent不够支持概率校准或大型泛化。

建议优先回答[三个理论问题](NEXT_THEORY_QUESTIONS.md)，而非重抄CNN/Flow/MCMC清单。实际向量、局部路径、同材料控制、scalar toy、配对CSV与所有费用见证据包。未取得廉价全域概率证书。
''')
    put('GPT_PRO_HANDOFF/NEXT_THEORY_QUESTIONS.md','''# 三个由R7实证产生的问题

## 1. 非均匀gauge上的真实条件prior

严格无pure、A_T满行秩时局部尺度曲线可构造；有pure或DG满66时不能沿用。能否在这类实际stratum上计算Gaussian64×Uniform2的disintegration核，量化把正确边缘独立化给材料函数造成的偏差？R7四点等权控制只是映射反例，不是答案。需要保留plateau不活跃子空间、非线性inverse-quintic、其余原像和正确体积，不能用S3单变量Q替代Full66。

## 2. 控制一阶谱后的真正高阶用途

R7 c1非零且第一阶采集冲突已达7%–16%，普通冻结/Full-J覆盖主要用途。能否构造相同真实材料距离、相同首阶/同源能量、不同高阶冲突的材料弧，并给非空、超过seam局部尺度的有效域？若不能，当前P1不应继续包装高阶算法。必须同时比较完整二/三阶向量、普通profile和有限端点，不只比较log斜率。

## 3. 材料图不唯一时，怎样真正认证一个任务？

R7的.02材料弧上均值变化很小、阈值占比偶有变化，但没有覆盖整个相容集和prior质量。是否存在在声明材料类中可计算的全域任务振幅/物理误差界，价格低于直接搜索/积分？标准两点下界和局部task dual还不是这样的证书。先固定一个连续任务与一个材料事件，在有可靠reference的小模型验证，不能由8个几何parent宣布Full66 Bayes风险。
''')
    memory=OUT/'RESEARCH_MEMORY_R7.md'
    if '## R7实际执行补记' not in memory.read_text():
        with memory.open('a') as f:f.write('''\n## R7实际执行补记\n\n新增37,456 RHS、1325 state/LU、256个零物理prior点与8个新locked测量驱动中心。E0/E1、主128状态、五节点路径、附加普通profile与非均匀gauge/coplanar控制都实际执行；具体原始来源为本目录RESULTS_SUMMARY.json、RESOURCE_SUMMARY.json、audit/results、analysis与E2封存。B在.02为8/8竞争、7/8 pair≤1；C6/8，未胜B；D1/8。rank66在多cutoff为88/163/186/191；四个局部非均匀尺度曲线支持stratum-specific gauge，而非普遍Full66商。GO_HIGHER_ORDER=NO_GO，EXACT_FIBER=NOT_ESTABLISHED，FINITE_NOISE=SCOPED_GO，FULL66_POSTERIOR=UNQUALIFIED，TAP=REDIRECT。新layout、完整概率体积、neural Explorer均NOT_RUN。所有修正前DEV和费用保留，原R6封账与历史HOLD均不变。\n''')
    print('Wrote final reports from closed receipts; new physical RHS=0')


if __name__=='__main__':main()
