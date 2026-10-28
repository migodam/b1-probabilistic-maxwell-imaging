# R7：数学对象、证明和未满足条件

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
