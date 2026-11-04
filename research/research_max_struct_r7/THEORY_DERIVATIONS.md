# R7 公式、证明与实际算法边界

本轮从Pro ZIP的数学建议实现；附件不是运行回执。下列标准代数不主张首创。状态标签：DERIVED_FINITE_MODEL、PROVED_UNDER_ASSUMPTIONS、EMPIRICAL、NOT_RUN。连续Maxwell物理精度另由验证判断。

## 1. 固定物理和变化变量

时间exp(-iωt)，二维TM，已知μ=1。实contrast m=εr−1，X=diag(m)，G含k²和真实圆盘积分/自项，Einc、S固定。辅助电流J满足

Q=Einc+GJ，r=J−XQ。

原物理解A=I−XG，J*=A⁻¹XEinc，F=SJ*；不要把I−GX与I−XG交换。完整原物理目标Phi=.5||F−Y||²/||Y||²+R(m)。所有辅助算法的候选最终按这个目标在线保留，而非按材料GT或辅助目标。

## 2. A-material 的精确Möbius坐标（DERIVED_FINITE_MODEL）

b_i=β−G_ii，d_i=1+b_i m_i，z_i=m_i/d_i，H_i=Q_i+b_i J_i。

u_i=β(J_i−m_iQ_i)/d_i=β(J_i−z_iH_i)。以每个固定2×2cell块的真实零对角Goff构造

RA=β(I+D_z Goff)(J−D_zH)。

展开只有z_i、z_i z_j且Goff_ii=0，因此固定其他坐标时整个残差（所有受影响行与源）对z_i严格仿射：RA=A_i+z_iB_i。它不是共同线性化，多个坐标仍有交叉乘积。邻居残差不能遗漏。

固定源metric已吸收入A_i/B_i时，q=||B_i||²，c=<B_i,A_i>，条件状态损失=.5q|z|²+Re(conj(z)c)+const。实m的可行z是Möbius圆弧；若v=Im b≠0，Im z=−v|z|²，即|z+i/(2v)|=1/(2|v|)。没有prior时是复点−c/q到合法圆弧的投影。q=0由prior/box决定；非零虚部的圆弧不能当整个复平面或整个圆。

实m下，条件状态损失可写

const+[n1 m+n2 m²]/[1+2Re(b)m+|b|²m²]，n1=Re c，n2=.5q+Re(cb)。

导数分子只有二次。分母在本次正材料盒内无零；不能外推ENZ/负介电/一般复张量。

### 原prior的MM而非更换prior

实际D边差分已除19.5，R=λsm mean[(Dm)²]+λtv mean[sqrt((Dm)²+eps²)−eps]。sqrt关于平方量凹，故其切线给精确相切上界。固定邻居时条件上界=.5h m²−s m+const。把该上界与精确有理状态损失相加，驻点条件乘以|1+bm|⁴后至多五次。

实现用批量5×5companion根、合法实根、端点和旧值比较；浮点根不是严格符号证书。原完整固定电流损失最后再验，劣化时退回旧值/更佳普通CIE seed。原TV/smooth数值不改。

四颜色=(x mod2,y mod2)。每色cell不共享固定2×2块或nearest-neighbor prior边，所以同色条件问题相互独立。固定源数、块大小下残差操作为O(N ns)，还需计求根常数和prior评估。精确条件最小化不推出整个交替算法、原Phi或材料误差单调。

### R7候选的代数简化：状态二次驻点+完整prior-MM接受

五次求根的实际常数可能超过原readout。条件状态损失的驻点只需解

n1+2n2 m+(n2 d1−n1 d2)m²=0，d1=2Re b，d2=|b|²。

新候选集合仅含这至多两个合法实根、两个box端点和当前值。每个候选仍用精确条件状态损失加**原prior的MM二次上界**评价。选择最小上界候选，因旧值在集合中，依然有 f(new)<=Q(new|old)<=Q(old|old)=f(old)。不需要把正则化条件问题解到全局最小。该候选是有原prior接受准则的近似条件优化，不可称精确正则化投影，也不删prior。

四颜色批量执行保留同一独立性条件，不再求5×5companion特征值。最大2sweeps固定，并使用继承material_inner_gtol=1e−8的精确(dP保留)投影梯度停止。这个容差来自原算法，不是开发后调整的成功门。浮点完整loss fallback保留；数学非增结论只在相应原loss/算术前提下成立，物理成像效果仍由全轨迹决定。

### 为什么便宜的内更新会改变材料恢复

精确目标的一次下降、原L-BFGS-B充分内优化、固定次数MM、投影梯度，都是不同的非根更新映射。CIE的PR历史、材料状态和阶段transport是耦合的；即使新方法每次材料块目标不增，也不能推出相同外循环路径。R7完整轨迹用于检验这一点，不能只比较readout微基准。

## 3. 完整有限材料响应（标准resolvent/Atlas既有）

在m0上A0=I−X0G，E0=Einc+GJ0，O=SA0⁻¹，C=GA0⁻¹。若变化ΔX=PDPᵀ，两个完整物理矩阵可逆，则

ΔF=OP (I−D KP)⁻¹ D EP，KP=PᵀCP，EP=PᵀE0。

证明：A0ΔJ=ΔX(E0+GΔJ)。用q个被改变cell的局部激发变量消元，再乘S即可。背景多次散射全部在O/K/E中；没有Born/Taylor截断。q是被改变状态的秩；两块各4cell的两幅度动作是rank8，32cell两幅度是rank32，不是rank2。3D各向同性体素通常贡献3个状态分量。

两个标量cell，D=diag(s,t)：δ=1−s k11−t k22+st(k11k22−k12k21)，

ΔF=(sU+tV+stH)/δ，U=O1E1，V=O2E2，H=−k22U−k11V+k12O1E2+k21O2E1。

对固定复线性白化W，e=W(F0−Y)，取T=[e,WU−k11e,WV−k22e,WH+(k11k22−k12k21)e]。Γab=Re<Ta,Tb>，p=(1,s,t,st)，则原数据损失=pᵀΓp/(2|δ|²)。加原prior的局部边差分增量即可精确排名。四个端口阵列和Gram需要已计费的setup，不能免费。

## 4. 一般q的端口Gram消元（DERIVED_FINITE_MODEL，标准二次型收缩）

本实验W仅为除||Y||的标量。令T=(I−DK)⁻¹D，Aq=OPᴴOP/||Y||²，Bq=EP EPᴴ，Zq=OP*(F0−Y)EPᴴ/||Y||²。则

Phi_data(D)=Phi_data(0)+Re tr(Zqᴴ T)+.5 tr(Tᴴ Aq T Bq)。

证明：展开||e+OP T EP||²，用trace循环恒等式。每次q维候选不需要重访全部端口。一般非可分离W不能随意沿用这个Aq/Bq因式分解；rank2四阵列公式仍须按实际W构造。

实现可带非均匀幅度权重：32个cell分成增加/减少两组，单元变化=公开合法slack×组幅度，组幅度在[0,1]。这是32状态秩、2维搜索。二次型收缩、Woodbury与Schur并非新发现；真正未决点是动作能否带来原方法没到达的好材料。

## 5. 动作选择与失败分析

首轮固定48个cell对（远距残差曲率、邻接对照、随机对照）及至多24个远距2×2块对，9点每轴粗搜索再从三个最优网格点局部精确优化；零点局部优化是同动作族的对照。候选与GT无关，非全盒全局最优证明。

若小动作不产生实质下降，一次结构修正扩大到rank32。使用完整物理Hessian的弱方向、当前材料平移梯度、负数据梯度及同秩随机对照；选增/减各16cell，有限动作仍按原Phi排名。大rank的setup/求解全部计费。没有无限递增Taylor阶数或未声明扫参。

原物理Hessian对实m为GN+residual反馈+原prior真Hessian。GNij=Re[(OᴴO)ij sum_s conj(Eis)Ejs]/||Y||²；令Z=Oᴴe Eᴴ/||Y||²，反馈=Re[Cij conj(Zij)+Cji conj(Zji)]。TV的真二阶权重是eps²/(d²+eps²)^(3/2)，不是IRLS的1/sqrt(d²+eps²)。代码使用前者。

自由face正Hessian只针对内点变量子空间；不包含活动bound的所有可行方向，不是约束局部极小或盆地证书。负特征值但budget未收敛也不能叫算法卡死。有限响应精确与动作族有效是不同问题。

每个top端点用原完整Maxwell验证并接同一corrector；另保留原锚点直接续跑。小有理解在直接续跑没有更好结果时，不记成新的basin能力。rank2/8逐端点预算数组在同一model中累计先前对照的费用，不能当各方法独立matchedCPU曲线；总研究成本与每段corrector自身成本分开列。

## 6. 极点、学习、网格与维度

材料方向t的det(I−tDC)=0、物理频率共振、Möbius坐标分母零、优化迭代特征值是四个不同概念。对可对角化DC，端口留数同时含源激发与接收可见性；极点距离或留数模长不包含与当前残差的相位配合。非正规算子还可有无有限极点的强多项式响应。这里直接稳定小矩阵求解，不学习排序裸极点。

固定四cell二维圆盘积分Green在h→0时Goff=O(h²|log h|)，固定β/正盒/物理域下额外one-hop修正可趋零；不推出非凸算法轨迹一致。保持物理半径会增加邻居数和成本。三维dyadic近场不继承此标度，纵/横极化不同。

3D本轮用独立标注的vector-adapted joint/dual-like/A与direct，不假装已有成熟二维CIE原论文的3D版本。DDA自项、材料链式导数、Mie几何误差和实际反演见BENCHMARK_3D。Fresnel的天线模型与缺失通道资格见FRESNEL_AUDIT。

若F(m1)=F(m2)，任何同数据确定变换都无法区分它们。R7没有以NN补足缺失信息；学习分支是否值得，仅能由候选生成/误差/成本证据判断，未训练不能推断学习普遍无用。
