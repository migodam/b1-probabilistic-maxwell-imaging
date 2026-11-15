# 双度量、有限 PR 与材料恢复：实际采用的数学

状态：主体局部推导先于本轮恢复评分；§7是根据非根程序问题补充的离线诊断分解，不用于修改独立确认算法。有限矩阵、理想局部映射、实际程序、连续 Maxwell 分开。以下线性代数、链式法则和精确线搜索性质不宣称为新发明。

## 1. 同一个度量必须是同一个材料相关函数

设 e=J−D_m Q，Q=Einc+GJ，L_m=I−D_mG。两个残差变换为

\[
P_0(m)=\operatorname{diag}\frac{b_0}{1+b_0m},\qquad
P_1(m)=(I+D_aG_o)D_p,
\]
\[
d_i=1+(b_1-G_{ii})m_i,\quad a_i=m_i/d_i,\quad p_i=b_1/d_i.
\]

G_o只含原不重叠四cell块的非对角项。度量是M_j=P_j*P_j，而不是P_j。b0=原current的阶段beta，b1=原A+material阶段beta；相同下标在两个块中代表完全同一个函数。

每个stage的current目标
\[
f_c=\tfrac12\|SJ-Y\|_{\Omega_d}^2+
\tfrac\rho2\|P_c(m)e\|_{\Omega_s}^2,
\]
material目标为其对应状态项加R(m)。00/11在单个stage共享势函数；01/10一般不共享势。阶段改变b1时不声称全程同一势。原Dual-CIE的两beta scalar不能被误标为00。

固定m时current系数a的梯度和Hessian作用来自
\[
H_c(a)=P_cL_mUa,\quad
H_c^*z=U^*L_m^*P_c^*z.
\]
\[
P_1^*z=\overline p\odot[z+G_o^*(\overline a\odot z)].
\]
这里只是原onehop的合法current接线，不是改成全局Green逆或增加传播阶数。

### 实际参数下没有利用奇异变换制造额外零点

继承R4的解析行和界，实际24²原核、合法m∈[.02,20]上，
\(\|D_aG_o\|_\infty\)统一上界在b1=1、.2分别为.0378433、.168100，分母最小模大于1。3D的6³ CM+RR核、χ∈[0,8]、24分量块，对b1=6、1分别为.0771186、.341525，分母最小模为1。这些界都小于1，故I+DaGo可逆、P1可逆，metric严格正定；P1 e=0等价e=0。

这是既有界在本轮真实参数上的资格，不是新恒等式发现。保守统一最小奇异值下界见`audit/METRIC_QUALIFICATION.json`。它只针对冻结的有限核及box；不证明L_m全域条件数、可辨识性或连续Maxwell精度。首次运行遇到JSON构造重复scope关键字的输出错误，修复命名后完成；0额外forward/LU，失败CPU未精确捕获，保留为unknown。

## 2. material中的度量导数不可省略

记u=P对角部分作用于e，v=G_o u，r=u+a⊙v。固定J/Q时
\[
\partial_{m_i}u_i=-b_1(Q_i+(b_1-G_{ii})J_i)/d_i^2,
\quad \partial_{m_i}a_i=1/d_i^2.
\]
因此对任意实方向h，
\[
dr=(I+D_aG_o)D_{du}h+D_{da\odot h}G_ou.
\]
目标导数是ρ Re⟨rΩs,dr⟩+dR。原R4/R7精确导数、Möbius候选和完整MM先验接受均保留。冻结P1后仅拟合e是另一个算法，本轮四臂没有混入它。

## 3. 理想精确最佳响应的结论与其边界

在精确无噪声物理根、固定可微活动面、无prior、固定current切空间U上，把所有复变量实化。令K=L_*U，Ch=D_hQ_*，D=WSU，源权重与ρ吸入相应块。
\[
H_c=K^TM_cK+D^TD,\quad B_m=C^TM_mC,
\]
\[
T=B_m^{-1}C^TM_mK H_c^{-1}K^TM_c C.
\]
要求Hc、Bm正定；冗余Fourier系数需先在实际独立空间中解释，不能对奇异正规矩阵强行求逆。根处e=0，dP乘e的一阶项消失；这不能用于离根忽略度量导数。

当Mc=Mm=M时，令A=K^TMC、B=C^TMC。联合二次型
\[
\begin{bmatrix}H_c&-A\\-A^T&B\end{bmatrix}\succeq0
\]
给出0≼A^THc^{-1}A≼B。因此T在B内积下自伴、谱在[0,1]。这属于加权块GS/投影的已知结构。Schur补B−A^THc^{-1}A是消元辅助模型曲率，不能直接叫统计物理Fisher。

等号方向要求可用current变化使Kδa=Ch且Dδa=0；在满足表示条件时对应当前端口不能看见的材料切向。严格收缩还需可观测性及足够谱隙；一致不代表容易恢复。Mc≠Mm时上述合同对称化不成立，可有局部扩张；扩张既不保证逃逸坏盆地，也不意味着必然破坏实际成像。

## 4. 非精确块误差给出实际比较应测什么

在同一个线性化模型中，用真实更新定义正规方程残差
\[
r_c=H_c\delta a-K^TM_cCh,\qquad
r_m=B_m\delta m-C^TM_mK\delta a.
\]
直接代入得
\[
\delta m=Th+
B_m^{-1}C^TM_mK H_c^{-1}r_c+B_m^{-1}r_m.
\]
这是精确线性分解，不是恢复误差上界。定义
\[
\gamma_{cm}=\|B_m^{-1/2}C^TM_mK H_c^{-1/2}\|_2,
\]
则同一线性模型中
\[
\|\delta m-Th\|_{B_m}\leq
\gamma_{cm}\|r_c\|_{H_c^{-1}}+\|r_m\|_{B_m^{-1}}.
\]
一致度量时由上一节Schur补得到γcm≤1；非一致时需要实际计算或界，不能默认。平方根都是正定矩阵的对称平方根，norm均在同一声明单位下。

只有两个余项确实比( T−I )h等理想效应小，才能把理想谱作为有限程序主要解释。例如理想扩张只有10^-7量级，并不意味着固定一次PR的10^-2或更大的方向误差可以忽略。这是理论效应可分辨性的要求，不是强制把所有内解做精确的算法建议。离根、box转面、prior、有限扰动另有余项；本轮保存它们而不把误差统一塞进一个“理论不准”的标签。

## 5. 有限精确线搜索 PR 不自动有一个根处Jacobian

即使第一步没有PR历史，固定m的二次current目标一次steepest-descent精确线搜索也给
\[
\delta a^+(h)=\alpha(h)b(h),\quad
b(h)=K^TM_cCh,
\quad \alpha(h)=\frac{b(h)^Tb(h)}{b(h)^TH_cb(h)}.
\]
该方向响应齐次、奇对称，却通常不满足加法性。简单Hc=diag(1,2)、b(h)=h时，两个坐标方向响应分别为(1,0)、(0,1/2)，合方向响应为(2/3,2/3)。故根处未必存在一个统一线性Fréchet导数矩阵。正负有限差分吻合本身不足以证明存在Jacobian，更不能把精确最佳响应T直接当一次PR程序DT。

这不是新的非光滑优化定理，而是本任务必须检查的实现边界。实际材料readout还含root候选切换、box、gtol停止；即便方向差分随epsilon收敛，也必须比较h1+h2的加法性及分支事件。非根且分母非零、分支固定时可计算局部JVP；PR历史运输/重置则是不同完整状态干预。

**本轮有限程序反例（EMPIRICAL，非稳定性定理）。** 同一个8²人工物理根、同完整Fourier/DPCS起态、同β与零prior，取理想01闭环的最大实特征方向（λ=1.002712449）。按精确最佳响应线性图，32轮的该mode系数为λ³²=1.09054859。实际有限PR／两扫描Arc，在±.001两个扰动下32轮的材料扰动范数比为.691965、.692145；每轮重置PR控制为.633511、.633683。没有换方向或搜索振幅，另两种符号沿同一方向。

因此理想root扩张没有在这个实际有限程序的声明尺度/时段中复现，不能据谱半径直接选择反馈。这个结果没有证明实际map严格ε→0的Fréchet导数不存在、正根渐近稳定或任意起点恢复；数值floor、gtol、Arc候选和有限迭代仍必须保留。详细512次native更新、48次新物理审计与所有完整state见`audit/finite_mode_trajectory_v1/`。

### 5.1 一个必须排除的短期解释：只扰动材料会触发场记忆回拉

**DERIVED_UNDER_ASSUMPTIONS，标准松弛模型的应用，不是实际PR的定理。** 记Gc=Hc^{-1}KᵀMcC、Am=Bm^{-1}CᵀMmK、T=AmGc。假设current只做固定松弛
\[
a^+=\theta a+(1-\theta)G_ch,\qquad h^+=A_ma^+,\quad0\leq\theta<1.
\]
若Th=λh，完整状态方向\((G_ch/λ,h)\)的放大率为r=θ+(1−θ)λ。因此λ>1时r仍大于1，单纯这种松弛不会消掉扩张。

但是只扰动材料、保持原current不动，即a0=0、h0=h，会先回拉：
\[
h_k=(1-\theta)\lambda[\theta+(1-\theta)\lambda]^{k-1}h_0.
\]
例如固定θ=1/2（仅解析示例，没有拟合实际PR），λ=1.002712449，32步材料误差仍约为初始的一半，同时完整状态系统的渐近mode仍扩张。故上一段32步收缩不能被包装为推翻Pro渐近不稳定的定理；它首先反驳的是对材料边缘图直接使用T³²的程序预测。

本轮已经补完这一诊断：同时按\((Gc h/λ,h)\)扰动current和材料，沿用原先的同一h、±.001、表示与PR策略，完成16条轨迹、512次实际更新。01的128次material readout均nit=0，投影guard梯度为2.56e−9至2.89e−9，低于原1e−8；32步材料扰动范数比保持1。完整state和费用见`audit/coupled_mode_trajectory_v1/`。这表明声明尺度下的实际程序进入有限guard deadband，既没有复现理想扩张，也没有证明渐近稳定。未更换epsilon或gtol搜索期望现象；固定θ解析模型不能自动套用到实际PR。

## 6. 原物理实现与辅助缺陷的关系

继承R13已建立的缺陷运输，不重跑恒等式：u=J−D_mQ，v=Q−Einc−GJ，
\[
J-J_{phys}=L_m^{-1}(u+D_mv),\quad
r_{phys}=r_{aux}-WSL_m^{-1}(u+D_mv).
\]
因此新metric减少auxiliary loss、材料步更大、或更强局部扩张，均不能代替端口可见缺陷的实际变化。当前物理审计产生这些量需收费，不把它称免费在线预测器。材料已经数据可行而错误时，纯改变metric不新增信息；五个共同球外失败存在合法GT可行点，也不证明唯一性。

## 7. 本轮允许检验的因果结论

四组合保持同globalPR/表示/初始完整状态时，可检验当前实验中哪一种块耦合改变了恢复集合；单独比较历史sourcewise CIE不够归因。00/11在固定stage的真实辅助势应非增，若违背先查实现/容差；这不是原物理Phi或GT误差非增的保证。

01对11的材料更新相同，但current metric不同；01对00的current metric相同，但材料目标、分母与读出不同。self-only控制才能进一步区分P1对角自作用与四cell非局部项；PR reset控制才可识别继承记忆的影响。结果必须与实际成本、已付费物理候选、prior、corrector及失败对象同时解释。

已有变量消元/可变Gram、CC-CSI、非势博弈与受控非精确交替优化都需引用；本文件不把这些结构更名为原创理论。本轮真实归因与新颖性判断见`PRIOR_ART_AND_NOVELTY.md`和主报告。


## 8. 根附近消失、非根状态保留的metric导数力

这是一项用于已有完整state的链式法则诊断，不是新的优化器。固定本轮实际J/Q，令e=J−D_mQ，r=P(m)e，C h=D_hQ。源metric为Ωs，则
\[
f_m(m;J,Q)=\frac\rho2\operatorname{tr}(r^*r\Omega_s)+R(m).
\]
对实h，
\[
dr=-P C h+(DP[m]h)e,
\]
因而
\[
g_m=g_{\rm const}+g_{\rm metric}+\nabla R,
\]
\[
g_{\rm const}=-\rho\,\Re C^*P^*(r\Omega_s),\quad
\langle g_{\rm metric},h\rangle
=\rho\Re\langle r\Omega_s,(DP[m]h)e\rangle_F.
\]
其中P0和P1都依赖材料，不能把P0的metric导数当成零。已有onehop full derivative包含两项；数值诊断直接复用原函数，再减构成力与prior，不写一个漏导数的替代求解器。

若声明box上P连续可微，\(\|DP[m]h\|\le L_P\|h\|\)，则
\[
\|g_{\rm metric}\|\le
\rho L_P\|P\|\|\Omega_s\|\|e\|_F^2.
\]
这解释了为什么精确物理根的理想一阶闭环不包含该项，却不能在离根的大构成缺陷状态省略它。\(g_{\rm const}\)只有O(\|e\|)的上界，没有普遍正下界；比值大可能来自构成力方向抵消，不能直接解释为新增材料信息或更强物理力。

本轮只在原7checkpoint×8干预枝的first completed step计算：m取step0，J/Q取current更新后step1，分别用实际material family、β与Ωs，比较三项对真实Δm的内积和角度。无需新forward、GT、有限差分或改变算法。它量化理想根模型遗漏的非根机制；是否能够预测恢复仍须由后续独立研究回答。


### P0归一化的一个明确符号后果

对正实b与P0=b/(1+bm)，每cell的metric导数分量精确为
\[
(g_{\rm metric})_i
=-\rho\frac{b^3}{(1+bm_i)^3}\|e_i\|_{\Omega_s}^2\le0.
\]
因此在固定非物理J/Q时，仅通过提高材料、减小归一化残差就存在一个同号梯度分量；它不是新的测量信息。完整构成力、metric导数力与prior一起决定真实readout。P1含复self分母和邻居传播，不能未经推导将这个逐cell同号结论直接套用。

本轮56个预定first-step状态，\(\|g_{metric}\|/\|g_{const}\|\)范围.1003–2.0673，各(pair,PR策略)七样本组的中位数约.907–1.105。构成力单独对真实Δm的work为正38/56，metric导数力为正6/56，完整梯度work为正0/56。它说明把真实material读出解释为“冻结metric的构成残差投影”会在这批非根state上误读方向。该样本全来自已预选state，不是恢复预测器，也不是支持修改material目标的新确认结果。
