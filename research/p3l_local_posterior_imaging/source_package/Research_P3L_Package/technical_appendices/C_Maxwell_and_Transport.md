> 技术附录采用独立章节、公式与文献编号；主报告为最终研究判断与实验配额的权威版本。本附录提供扩展推导，不报告新的 Maxwell 实验。

# Research P3-L：Maxwell 有限反馈、概率决策认证与流形概率 PDE

## 技术补充：证据边界与主要判断

本稿已经阅读本轮下载的 P1《MODEL_AND_PRIOR》、P2《RESEARCH_REPORT_ZH》《METHOD_COMPARISON》《RESOURCE_SUMMARY_FINAL》以及 A17-X × F3《RESEARCH_REPORT_ZH》。没有启动新的 Maxwell forward、训练或付费实验。数学推导针对有限维离散模型；从 TM24 外推到一般三维矢量 Maxwell 的地方均另列条件。

最重要的判断是：**物理几何的优先用途，应是控制概率计算的误差与费用，而不是把弱方向轨迹直接解释为概率分布。** 最有希望的物理分支，是把有限 Maxwell 反馈、原始/伴随残差与概率接受决策连起来，使一次已支付的状态求解服务于多次后验提案，同时保留原 Bayesian 目标。它可以服务于简单 pCN、blocked MCMC、MALA、conditional transport，未必需要新的神经网络。

但这种路线也不是一片空白：LIS 与 state-ROM 联合缩减、Stein variational reduced basis、DINO 加速 geometric delayed-acceptance MCMC 均已有直接先例。普通 resolvent、Davis–Kahan、dual-weighted residual、流形 Langevin 本身不能作为新理论。新的论文必须给出 Maxwell 材料模型中特有且实际非空的误差界、刷新准则、预算选择规律，以及相对于优化 Full/Krylov/ROM 的有效样本与材料积分成本优势。

### 既有成本证据必须保留

P2 的共同评价器审计如下；fallback 使用已付费 full predictions 提供的离线误差，因此尚非可部署在线控制器。

|评价器|setup + query RHS|fallback/48|setup 摊销 ms|
|---|---:|---:|---:|
|optimized Full|384|0|5.346|
|exact T|1056|0|9.487|
|普通 Krylov|88|7|23.563|
|goal-ROM|192|4|38.249|
|A17 exchange|296|17|59.490|

A17 摊销费用约为 Full 的 11.1 倍；它没有因为 RHS 数较少就更快。competition audit 的 7 个商坐标对应 575 个非零对角材料改变量，最优 rank-7 截断仍遗失约 58%–64% 的算子谱范数；小 latent 维数不等于小端口秩。相同 current trial/test spaces 下，A17 Schur exchange 与直接 Galerkin 等价；可争取的差异只能来自更好的空间、便宜更新或更有效的误差认证。

## 1. 统一概率约定

P1 使用 circular-complex Gaussian：

\[
\eta\sim\mathcal{CN}(0,\sigma^2 I_{256}),\qquad
\Phi(q)=\frac{\|F(g(q))-y\|_{\mathbb C}^{2}}{\sigma^2}.
\]

把每个 source 的 Re/Im 正确拼接，并用 \(W=\sqrt2 I/\sigma\) 白化后，令 \(f=W\mathcal R F\)、\(y_w=W\mathcal R y\)，则

\[
\Phi(q)=\tfrac12\|f(q)-y_w\|_{\mathbb R^{512}}^2.
\]

下文全部 likelihood 误差式均使用这一实白化定义。若代码使用 complex residual 除以 \(\sigma\)，前面不能再乘 \(1/2\)。任何潜在 gain/nuisance 边缘化都必须另加到原模型；P1 固定 gain=1，不可默默 profile gain。

令 \(p_0(q)\) 为 P1 的 Gaussian64 × uniform2 原始 latent prior。无论提案使用何种材料度量，目标首先是

\[
\pi(q\mid y)=Z_y^{-1}p_0(q)e^{-\Phi(q)}.
\]

一个声明过的局部集合 \(\mathcal C(y)\) 可定义局部目标
\(\pi_\mathcal C=\pi\,1_{\mathcal C}/\pi(\mathcal C)\)。有限噪声球可用作局部研究范围，但 Gaussian likelihood 本来并不把球外概率置零。

## 2. 有限 Maxwell 材料反馈：可复用对象究竟是什么

### 2.1 多 source 共同材料的 current 形式

设 \(X=\mathrm{diag}(\chi)\)，\(G\) 是离散内部 Green 算子，\(S\) 是接收算子。多个 source 的电流组成矩阵 \(j\)，共同满足

\[
L(X)j=XE^{\mathrm{inc}},\qquad L(X)=I-XG,\qquad F=Sj.
\]

在锚点 \(X_0,j_0\)，令 \(E_0=E^{\mathrm{inc}}+Gj_0\)，
\(\Delta X=X_1-X_0\)、\(k=j_1-j_0\)。直接相减得到

\[
\boxed{(L_0-\Delta XG)k=\Delta X E_0.}\tag{P1}
\]

这一式同时保持所有 illuminations 的共同 \(\Delta X\)。如果从 current 轨迹倒推材料，则每个像素必须有一个相同的复材料解释所有 source；逐 source 的 receiver-dark current 不自动可实现。

令 \(R_0=L_0^{-1}\)、\(H_\Delta=R_0\Delta XG\)、
\(v_\Delta=R_0\Delta X E_0\)。只要 \(I-H_\Delta\) 非奇异，

\[
k=(I-H_\Delta)^{-1}v_\Delta,\qquad
\Delta F=S(I-H_\Delta)^{-1}R_0\Delta X E_0.\tag{P2}
\]

这是真实有限变化，而不仅是 Jacobian。低秩 current approximation 近似的是受右端 \(\Delta X E_0\) 激励、又被 \(S\) 读出的响应；它不要求 \(\Delta X\) 低秩。

### 2.2 可计算的充分半径，以及它的限度

在一致算子范数下，若有认证上界 \(\|H_\Delta\|\le \nu<1\)，则

\[
\|L_1^{-1}\|\le \frac{\|R_0\|}{1-\nu},\qquad
\|k-v_\Delta\|\le \frac{\nu}{1-\nu}\|v_\Delta\|,
\]

\[
\|\Delta F-Sv_\Delta\|
\le \frac{\|S H_\Delta\|}{1-\nu}\|v_\Delta\|.\tag{P3}
\]

证明是 Neumann 级数或 resolvent identity。这个标准定理的作用是给出一个允许复用锚点的充分区域，不是证明区域外存在病态。已读 competition audit 的 \(\theta=2\) 反例已经说明：充分界失效时，完整 feedback 仍可能稳定。

认证 \(\nu\) 与 \(\|R_0\|\) 本身有成本。把浮点 SVD、少数随机向量上的最大放大率或历史真实误差直接当确定性上界是不合法的。大型 Maxwell 问题可能需要连续能量范数中的 inf-sup 下界、可验证的预条件残差界或严格的矩阵范数上界；这些界可能很松，应允许自动退回 Full。

### 2.3 Decoder 曲率与物理反馈曲率可以精确分离

对 \(C^2\) 材料路径 \(q(t)\)，\(X(t)=\mathrm{diag}(g(q(t)))\)。微分 current 方程：

\[
L\dot j=\dot X E,\qquad
L\ddot j=\ddot X E+2\dot XG\dot j.\tag{P4}
\]

其中

\[
\ddot X=\mathrm{diag}\!\left(
D^2g(q)[\dot q,\dot q]+Dg(q)\ddot q
\right).
\]

对于直 latent 路径 \(q_0+tv\)，

\[
\boxed{
\frac{d^2F}{dt^2}
=SR\,\mathrm{diag}(D^2g[v,v])E
+2SR\,\mathrm{diag}(Dg\,v)GR\,\mathrm{diag}(Dg\,v)E.
}\tag{P5}
\]

第一项是 decoder 的二阶材料变化，第二项是 Maxwell 再散射反馈。两者可能同向增强或相互抵消；仅查看 \(Dg\) 或仅查看 Maxwell Hessian 都不能确定总曲率。Research L/F 已经包含相关分解，因此这里是统一后的计算接口，不是新的首创恒等式。

给定同一 LU，\(\dot j\) 与 \(\ddot j\) 分别只需相同 \(L\) 的多 RHS 回代；每个指定方向的二阶变化不必重做 LU。是否比有限差分或 JVP/adjoint 更便宜，应由 source 数、receiver 数和当前缓存决定。

### 2.4 比条件 MAP ridge 更温和的 \(C^2\) 中心曲面

选择固定的强 left-data frame \(C_s\)，并写 \(q=q_0+V_wa+V_sb\)。定义

\[
h(a,b)=C_s^\top[f(q_0+V_wa+V_sb)-f(q_0)]=0.\tag{P6}
\]

若 \(D_bh=C_s^\top J_0V_s\) 可逆，隐函数定理给出局部 \(C^2\) 图 \(b=\beta(a)\)。当 frame 来自一阶强弱 SVD 时，\(C_s^\top J_0V_w=0\)，因此

\[
\beta'(0)=0,\qquad
\beta''(0)[u,v]=-(C_s^\top J_0V_s)^{-1}
C_s^\top f''(q_0)[V_wu,V_wv].\tag{P7}
\]

它只保持选定的强观测分量，不能称为 full constant-data fiber。它可用作后验 tube 的中心，再对强方向余量积分。

正则性区别：一般 \(g\in C^2\Rightarrow U\in C^2\)，用 \(U_b=0\) 定义的条件 MAP ridge 只自动是 \(C^1\)，经典二阶曲率需要额外正则性。P1 的具体 smootherstep 具有局部 \(C^{2,1}\) 正则性，故 ridge 可取得 \(C^{1,1}\) 与几乎处处曲率；固定 saturation stratum 内则更光滑。

特别地，tube
\[
q=q_0+V_wa+V_s\{\beta(a)+u\}
\]
具有三角 Jacobian；若 \([V_w,V_s]\) 固定可逆，则全维行列式就是 \(|\det[V_w,V_s]|\)。曲面自身的 \(\sqrt{\det G}\) 不应被重复乘进这个全维 tube 的概率。

## 3. 原始与伴随残差：由状态误差走到概率误差

### 3.1 接收任务的确定性认证

以下用一个统一向量表示全部 source；复数可经实化处理。设

\[
A(q)x=b(q),\quad f(q)=C x,\quad
r_p=b-A\widetilde x.
\]

若 \(\beta(q)>0\) 是 \(\sigma_{\min}(A(q))\) 的可靠下界，
\[
\|x-\widetilde x\|\le \|r_p\|/\beta.
\]

但完整状态误差可能远大于接收误差。令 adjoint block 满足
\(A^\ast Z=C^\ast\)，近似 \(\widetilde Z\) 的 dual residual 为
\[
R_d=C^\ast-A^\ast\widetilde Z.
\]
则
\[
f-C\widetilde x=\widetilde Z^\ast r_p+
R_d^\ast A^{-1}r_p.
\]

因此对 dual-corrected output
\[
\widetilde f_c=C\widetilde x+\widetilde Z^\ast r_p
\]
有
\[
\boxed{\|f-\widetilde f_c\|
\le \frac{\|R_d\|\,\|r_p\|}{\beta}.}\tag{P8}
\]

这是标准 goal-oriented error identity。Maxwell 特点是锚点接收 adjoint、已有 A17 左右 current spaces、强观测 frame 可以复用；代价取决于需要覆盖多少接收任务。对于 32 receivers，多 RHS adjoint 不应逐 latent 方向重复计费。

### 3.2 更便宜的 likelihood 专用 adjoint

若目标只是 \(\Phi\)，不必总构造所有 receiver adjoints。令
\(\widetilde r=C\widetilde x-y_w\)，
\(A^\ast z=C^\ast\widetilde r\)。近似 \(\widetilde z\) 满足
\(r_d=C^\ast\widetilde r-A^\ast\widetilde z\)。恒等式为
\[
\Phi-\widetilde\Phi
=\Re\langle\widetilde z,r_p\rangle+
\Re\langle r_d,A^{-1}r_p\rangle+
\tfrac12\|CA^{-1}r_p\|^2.\tag{P9}
\]

令 \(e_f\) 为可靠 output-error 上界，
\(e_{\rm lin}=\|r_d\|\|r_p\|/\beta\)，可得到不对称区间：
\[
\Phi\in\left[
\widetilde\Phi+\Re\langle\widetilde z,r_p\rangle-e_{\rm lin},\
\widetilde\Phi+\Re\langle\widetilde z,r_p\rangle+e_{\rm lin}
+\tfrac12e_f^2\right].\tag{P10}
\]

这对 MH 决策比“相对 forward error 很小”更直接。已有 Stein variational reduced basis 和 DINO 工作已经研究 potential/gradient 误差；新的增量应是证明该界在所研究 Maxwell acquisition、噪声尺度与 trajectory 上足够紧，并使刷新/接受决策更便宜。

### 3.3 只有 output norm error 时的通用转换

若 \(\|f-\widetilde f\|\le \varepsilon_f\)，
\(\widetilde r=\widetilde f-y_w\)，则
\[
|\Phi-\widetilde\Phi|
\le \|\widetilde r\|\varepsilon_f+\tfrac12\varepsilon_f^2.\tag{P11}
\]
更准确的区间是
\[
\Phi\in
\left[\tfrac12\max(0,\|\widetilde r\|-\varepsilon_f)^2,\
\tfrac12(\|\widetilde r\|+\varepsilon_f)^2\right].\tag{P12}
\]

噪声白化后的总 residual 往往是 \(\sqrt{512}\) 的尺度。因此“forward error 小于单通道噪声”仍可能产生不可忽略的 log-likelihood 误差；不能用一个与测量数无关的相对阈值保证 posterior 正确。

### 3.4 A17 在这里应控制哪一个量

在**固定的候选材料** \(q\) 上，令已选 current space 为 \(V\)，
\[
R_V=V(V^\ast LV)^{-1}V^\ast,\quad x_V=R_Vb,\quad
r_V=b-Lx_V.
\]
新增 current space \(C_e\) 的 A17 factors 为
\[
\bar C_e=(I-R_VL)C_e,\quad
\Gamma_e=C_e^\ast L\bar C_e,\quad
N_e=C_e^\ast r_V,\quad Y_e=C_{\rm obs}\bar C_e.
\]
当相应 Schur block 非奇异时，扩大 Galerkin space 的观测改变量为
\[
\delta f_e=Y_e\Gamma_e^{-1}N_e.\tag{P29}
\]
这里 \(N_e\) 是本次 fixed-material state residual 的注入，不是用原 A17 tangent \(B\) 冒充状态右端。P29 只是 A17 的 Schur 恒等式用于 state solve；和同空间直接 Galerkin 完全相同。

与当前 likelihood residual \(\widetilde r=f_V-y_w\) 结合，
\[
\widetilde\Phi_{V+C_e}-\widetilde\Phi_V
=\Re\langle\widetilde r,\delta f_e\rangle
+\tfrac12\|\delta f_e\|^2.\tag{P30}
\]
这显示：单看 \(\|Y_e\|\)、\(\|N_e\|\) 或 \(\|\delta f_e\|\) 都不足以判断后验决策影响，还需要和当前数据残差的干涉。该 signed change 仍只是两个 ROM 之间的变化；剩余 tail 必须用 P8–P10 认证。

推荐新选择目标是：在一组当前 MH 候选状态对上，选择能以最小新增 RHS/刷新费收窄 P13 区间的 \(C_e\)。这比“尽可能提高全局 Jacobian 保真度”更直接服务概率推断，也与 A17 原来的 observation/injection 双侧思想一致。若 ordinary goal-ROM 用同样候选空间与 Schur cache 得到相同收益，应把结论归于 posterior-directed current-space selection，而不是归于 A17 公式的独占优势。

## 4. 认证式接受决策：计算缓存可以适应，统计目标不必改变

### 4.1 与完整 MH 逐步一致的决策

固定一个可评价的 proposal \(Q(q'|q)\)，其 log MH ratio 为
\[
\ell(q,q')=\log p_0(q')-\log p_0(q)
+\Phi(q)-\Phi(q')
+\log Q(q|q')-\log Q(q'|q).
\]
如果在两个状态有认证区间 \(\Phi(q)\in[L_q,U_q]\)，则令已知的 prior/proposal 部分为 \(c\)：
\[
\ell\in[\ell_-,\ell_+]
=[c+L_q-U_{q'},\ c+U_q-L_{q'}].\tag{P13}
\]

抽一次 \(u\sim U(0,1)\)，令 \(s=\log u\)。因为 \(s<0\)，MH 接受条件等价于 \(s<\ell\)：

1. 若 \(s<\ell_-\)，直接接受；
2. 若 \(s>\ell_+\)，直接拒绝；
3. 其余情况加精 Krylov/ROM、补充 dual space 或调用 Full；始终保留原来的 \(u\)。

**命题：** 若所有区间真实包含相同离散 full model 的 \(\Phi\)，加精区间收敛到真值，proposal 不随加精被重新选择，则该算法与使用同一 proposal 和随机数的 Full-MH 路径逐步一致；边界等号以概率零发生，故终止几乎处处成立。

证明仅是区间包含性与 MH 阈值比较。它不是新 Bayesian 不变性定理，但提供很强的工程性质：仅用于缩紧认证区间的缓存可任意依历史适应，因为它没有改变一次统计接受决策。若缓存同时改变 proposal，则需要另行满足 adaptive MCMC 条件。

这里的“精确”是一个条件数学结论。有限精度 LU 也有误差：若要确定性认证，应把 backward residual、可靠稳定性下界和必要的舍入界包含在最后区间里，而不是把 Full 调用的浮点返回值视为数学真值。若只与规定 tolerance 的 Full reference 比较，应写“与该数值参考逐步一致”，不升级为无限精度认证。P2 当前尚未实现这种在线证书。

### 4.2 认证精度直接控制 fallback 频率

给定 \(q,q'\) 后，\(\log u\) 的密度为 \(e^s1_{s\le0}\le1\)。所以初始区间宽度 \(w=\ell_+-\ell_-\) 满足
\[
\boxed{\Pr(\text{需要加精}\mid q,q')\le \min(1,w).}\tag{P14}
\]

因此可以训练/选择 current basis 来缩小**将被采样的状态对上的 log-ratio 区间宽度**，而不是要求所有 states 的全场误差都相同。P8–P10 将这一宽度进一步连到 paired primal/dual residual、Maxwell 稳定性与当前 noise。

这给出可验证的算法预测：在相同 proposal 与随机数下，认证 evaluator 应与 Full-MH 返回完全相同接受序列，且新增 Full 调用比例应被区间宽度统计解释。若界太松，P14 会诚实地退化为无保证；不能把经验分位数直接叫确定性认证。

### 4.3 与 delayed acceptance 的严格区别

给定固定且正的 surrogate density \(\widetilde\pi\)，标准 delayed acceptance 使用
\[
\alpha_1=1\wedge\frac{\widetilde\pi(q')Q(q|q')}
{\widetilde\pi(q)Q(q'|q)},\qquad
\alpha_2=1\wedge
\frac{\pi(q')\widetilde\pi(q)}{\pi(q)\widetilde\pi(q')}.
\tag{P15}
\]
两阶段核保持 \(\pi\) 不变，但通常比同一 proposal 的 Full-MH 每步拒绝更多；只能靠更便宜的每步费用补偿相关性。认证区间 evaluator 的接受序列则与原 Full-MH 相同。

若 surrogate 在原 posterior 正质量区域人为置零，可能破坏支持和遍历性。若 surrogate 随当前状态/历史变化，不能在两个方向上偷偷用不同函数却沿用 P15。最简单的研究设计是训练/适应后冻结，或明确使用有理论条件的 state-dependent DA；普通“每步重训、再做一次 MH”并不自动正确。Conrad et al. 的渐近精确局部近似有专门 refinement/ergodicity 条件，不等于任意在线 surrogate 都合法。

### 4.4 近似概率认证不能冒充确定性认证

若每个区间只有条件失败概率 \(\delta_t\)，有限 \(T\) 步可通过 union bound 得到路径耦合失败概率不超过 \(\sum_{t\le T}\delta_t\)，前提是概率保证对自适应历史也成立。固定点验证集上的 99% coverage 不能直接升级为所有自适应 MH 查询的 99% 保证。无限长链若使用固定非零错误概率，也不能据此宣称 exact posterior；可选可求和错误预算或 full correction。

## 5. 从 forward 偏差到后验与材料风险

### 5.1 Uniform potential error 的一个方便界

在共同支撑 \(\mathcal C\) 上令
\(\delta(q)=\widetilde\Phi(q)-\Phi(q)\)，并假设
\(\operatorname{osc}_{\mathcal C}\delta
=\sup\delta-\inf\delta\le\omega\)。同一 prior 下归一化后
\[
\boxed{\|\widetilde\pi_\mathcal C-\pi_\mathcal C\|_{\rm TV}
\le \tanh(\omega/4).}\tag{P16}
\]
特别地 \(|\delta|\le\varepsilon\) 给出
\(\mathrm{TV}\le\tanh(\varepsilon/2)\)。

证明要点：写 \(\widetilde\pi=w\pi/\mathbb E_\pi w\)，
\(w=e^{-\delta}\in[a,b]\)、\(b/a\le e^\omega\)。在给定均值时，凸性使最坏情形出现在 \(a,b\) 两端；对端点质量优化，最大 TV 为
\((\sqrt b-\sqrt a)/(\sqrt b+\sqrt a)\)。该界强调 potential 的加性常数无关紧要。它是一般测度换权的初等结果，不主张原创。

PDE Bayesian inversion 的 Hellinger 稳定性和 approximation error 理论已有成熟文献（Stuart 2010；Zahm et al. 2022）。本研究应明确哪个误差是物理 surrogate、哪个是弱坐标截断、哪个是局部集合漏质量。

### 5.2 为什么这个界与成像直接相关

若所有合法材料在材料 Hilbert norm 下直径不超过 \(D_\chi\)，最大耦合给出
\[
\|\mathbb E_{\pi}g-\mathbb E_{\widetilde\pi}g\|
\le D_\chi\,\mathrm{TV}(\pi,\widetilde\pi).
\]
真正 posterior mean \(m\) 与近似 mean \(\widetilde m\) 的平方损失条件风险满足
\[
R_\pi(\widetilde m)-R_\pi(m)
=\|\widetilde m-m\|^2
\le D_\chi^2\mathrm{TV}(\pi,\widetilde\pi)^2.\tag{P17}
\]
这控制的是相对正确 posterior mean 的额外风险，不保证每个真实物体都改善，也不说明相对基础重建一定存在足够 headroom。

若 \(\alpha_y=\pi(\mathcal C^c|y)\)，则
\(\mathrm{TV}(\pi,\pi_\mathcal C)=\alpha_y\)。即使局部 sampler 完全精确，遗漏质量仍构成独立误差。没有全维 reference 或可证明的 tail bound 时，只能报告“局部 posterior”，不能从好看的 likelihood fit 认证 \(\alpha_y\approx0\)。

不推荐以此 uniform 最坏界作为唯一指标：它可能很保守。可以同时报告 certified bound 与 posterior-weighted 经验误差，但两者标签必须分开。

## 6. 几何刷新与计算复用

### 6.1 预测“何时弱方向不再弱”

设广义谱为 \(Hv=\lambda Mv\)，\(M\succeq mI>0\)，
\(K=M^{-1/2}HM^{-1/2}\)。沿路径，
\[
\|\dot K\|
\le m^{-1}\|\dot H\|+m^{-2}\|H\|\|\dot M\|.\tag{P18}
\]

对 \(H=J^\top J+R\)，
\(\|\dot H\|\le2\|J\|\|\dot J\|+\|\dot R\|\)；
对材料 Gram \(M=Dg^\top W_\chi Dg+\epsilon R\)，
\(\dot M\) 来自 \(Dg,D^2g\)。P5 可为 \(\dot J\) 给出方向性估计或邻域上界。

若 \(\|K(q)-K(q_0)\|\le L_K\ell\)，其中 \(\ell\) 是所用坐标/度量的累计路径长度，则 Weyl 给出
\[
|\lambda_i(q)-\lambda_i(q_0)|\le L_K\ell.
\]
若 weak/strong cluster 的初始 gap 为 \(\gamma\)，且
\(L_K\ell\le\gamma/4\)，Davis–Kahan 给出诸如
\[
\|\sin\Theta\|\le 2L_K\ell/\gamma.
\]
这可决定下一次谱 refresh 的最晚时机。它只有充分意义；小 gap 意味着方向标签不稳定，宜跟踪子空间/平滑谱滤波而不是强制逐个特征向量连续。

Research L 已有这一思路。新的重点应是：上界是否可在线取得、是否比直接刷新便宜、是否能正确预测实际失效点。若统计位移很小而 decoder 接近饱和接缝，材料 Gram 条件数变化可能比 Maxwell 曲率更早触发 refresh。

### 6.2 动态低秩提案的支持问题

若在每个状态只做
\(q'=q+V_w(q)\xi\)，其 proposal 支撑是一个 \(d\) 维仿射平面。方向一旋转，通常 \(q\notin q'+\mathrm{Range}(V_w(q'))\)，逆向 proposal 为零。不能只写一个 likelihood ratio 就认为 MH 纠正了它。

可选方案：

- 使用 full-rank soft mobility，在强补空间保留非零方差；
- 在固定且已知的 chart 中提案，并使用正确 Jacobian；
- 使用有双向投影检查和体积修正的 constrained integrator；
- 先冻结一个 full-dimensional coordinate decomposition，分块更新强、弱坐标。

位置依赖 MALA 即使省略连续 SDE 所需的 divergence drift，只要把实际 Gaussian proposal 及其反向密度准确代入 MH，仍可保持目标不变；省略项会影响效率。省略 reverse proposal density 则会破坏正确性。

### 6.3 公平 break-even 判据

令
\[
c_{\rm cert}=C_{\rm setup}/N_q+C_{\rm query}+C_{\rm bound}
+p_{\rm fb}C_{\rm Full}+C_{\rm refresh}/N_q.
\]
若 proposal 与原 Full-MH 相同且认证决策逐步一致，纯 evaluator 胜出的必要比较就是
\(c_{\rm cert}<C_{\rm Full}\)。若核被 DA/RL/geometry 改变，还必须比较
\[
\boxed{c_{\rm alg}\tau_h^{\rm alg}<
C_{\rm Full}\tau_h^{\rm Full},}\tag{P19}
\]
其中 \(\tau_h\) 应针对材料均值、材料区域平均、risk 等任务统计量，而不只针对 log likelihood。

推荐账本同时保存：Full states、LU、source/adjoint/JVP RHS、矩阵作用、current-space rank、谱分解时间、证书时间、refresh、fallback、训练、冷启动与摊销总 wall。RHS 只是其中一个维度。

Full-J 可通过 32 receiver adjoints 与 8 source states 共享 factorization，不应默认每个 66 latent 方向都支付 8 个独立 tangent RHS。dense TM24 的高效 BLAS 与 cache 基线尤其重要；A17 或 goal-ROM 若需大量 Python callback、逐步重建小矩阵，其名义低秩可能没有实际意义。

## 7. Maxwell 能量约束的正确用途

以 \(e^{-i\omega t}\) 为时间约定，被动各向同性电介质有
\(\operatorname{Im}\varepsilon\ge0\)。无源材料体内的耗散功率为
\[
P_{\rm abs}=\frac{\omega}{2}
\int_\Omega
\{E^\ast\operatorname{Im}\varepsilon E+
H^\ast\operatorname{Im}\mu H\}\,dx\ge0.
\]
对有适当辐射条件、背景和完整通量面的问题，Poynting identity 给出输入、输出和耗散平衡。在常规无耗 exterior 的散射设置，
\[
P_{\rm ext}=P_{\rm sca}+P_{\rm abs}.
\]

Miller et al. 2016 的 Eq. (5)–(7)、(23) 清楚区分 extinction、scattering 与 absorption；其有限响应界依赖真实材料损耗，不能用 lossless conservation 替换。

对于 reciprocal constitutive tensors、匹配的 ports/归一化，完整散射算子满足相应互易对称；无耗情形是通量度量下的 unitary，被动有耗情形是 contractive。有限、带噪、未校准的 32 个近场读数不是完整 scattering matrix，不能直接要求它们 Euclidean-unitary。P1 的背景 Green operator 与宿主吸收约定也必须核对；不能把真空 exterior 的 optical theorem 原封不动用于不同背景。

这些条件可以：

- 检查 approximate current 是否违反基本物理；
- 在有严格正耗散下帮助建立稳定性常数；
- 为 surrogate 提供状态、损耗、输出耦合特征；
- 在 material prior 未规定被动性时提供先验支持信息。

但 P1 的 decoder 已保证被动性；同一材料精确求解后能量平衡本来成立，再给通过平衡的材料额外 posterior 权重，相当于无依据地改写目标。满足能量恒等式也不意味着解释了这份观测。互易性同样是物理法则，不是新的独立 noisy dataset。

Maxwell 完整数据唯一性结果也不能推出 8 illuminations × 32 receivers 条件下的数值识别性。反之，当前 Jacobian 满列秩已经足以排除局部正维 exact-data fiber；不需要引用全局连续逆散射唯一性来论证这一点。

## 8. 局部曲面上的概率 PDE：先固定参考测度

### 8.1 度量不等于概率

设局部 surface \(\psi:A\subset\mathbb R^d\to\mathcal M\) 为浸入，正定 metric matrix 为 \(G(a)\)。它可取材料 pullback metric，也可取计算上选择的几何 metric。体积元为
\[
d{\rm vol}_G=\sqrt{|G(a)|}\,da.
\]

如果来自原 latent prior 与强方向积分的目标是
\(p_a(a)\,da\)，那么同一测度相对于该体积元的密度必须是
\[
\rho(a)=p_a(a)/\sqrt{|G(a)|},\qquad
V=-\log\rho=-\log p_a+\tfrac12\log|G|.\tag{P20}
\]

不能把 \(\sqrt{|G|}\) 当成“更物理”的附加 prior。Xifara et al. 2014 正是澄清 Riemannian volume 与 Lebesgue invariant density 的差别。P1 saturation 可能使 material Gram 退化；若没有 immersion，应换 chart、分层处理或在 latent 中保留正定 mobility，而不是对奇异 Gram 求逆。

### 8.2 不变性与 Fokker–Planck

取生成元
\[
\mathcal L f=\Delta_G f+
\langle\nabla_G\log\rho,\nabla_G f\rangle_G.
\tag{P21}
\]
对相对于 \(d{\rm vol}_G\) 的时间密度 \(h_t\)，
\[
\partial_t h_t
=\operatorname{div}_G(h_t\nabla_GV+\nabla_Gh_t)
=\operatorname{div}_G\left(\rho\nabla_G(h_t/\rho)\right).
\tag{P22}
\]
在无边界或零通量边界
\[
n_G\cdot(h_t\nabla_GV+\nabla_G h_t)=0
\]
下，\(\rho\) 为 invariant density。分部积分给出
\[
\frac{d}{dt}{\rm KL}(h_t\|\rho)
=-\int h_t\|\nabla_G\log(h_t/\rho)\|_G^2d{\rm vol}_G\le0.
\tag{P23}
\]
这不自动给出快速收敛；指数 KL 收敛还需要 log-Sobolev 常数等条件，瓶颈或多 chart geometry 可使其很差。

坐标 Itô SDE 若直接以 \(p_a\,da\) 为目标，则是
\[
da^i=
\{G^{ij}\partial_j\log p_a+\partial_jG^{ij}\}\,dt
+\sqrt2\,B^i_k\,dW^k,\quad BB^\top=G^{-1}.
\tag{P24}
\]
若以 \(\rho=e^{-U}/Z\) 相对于 \(d{\rm vol}_G\) 为目标，则漂移写成
\[
-G^{-1}\nabla U+\operatorname{div}G^{-1}
+\tfrac12G^{-1}\nabla\log|G|.
\]
两式经过 P20 是相同动力学，不能混用。

连续动力学不需 \(Z\)，因为只使用 log density 梯度。有限步 Euler 的 stationary law 通常有偏；若要 exact discrete target，可使用 MH correction，或为离散误差给出并验证界。拒绝为零不是 posterior 正确的证明。

### 8.3 Probability flow ODE 的真正未知量

对 \(\partial_t h=-\operatorname{div}_G(hb)+\Delta_Gh\)，相同边际的 continuity equation 速度为
\[
v_t=b-\nabla_G\log h_t.
\tag{P25}
\]
因此 probability-flow ODE 需要**时间依赖的当前密度 score**，不能把它直接替换成最终 posterior score。若 \(b=\nabla\log\rho\) 且把 \(h_t\) 偷换成 \(\rho\)，速度会变为零，根本不能从 prior 搬运到 posterior。

可以研究 annealing path
\[
\rho_t(q)\propto p_0(q)e^{-t\Phi(q)}1_\mathcal C(q).
\]
固定 reference volume 时，
\[
\partial_t\log\rho_t=-\Phi+\mathbb E_{\rho_t}\Phi.
\]
令 \(v_t=\nabla_G\psi_t\)，需解
\[
\operatorname{div}_G(\rho_t\nabla_G\psi_t)
=\rho_t(\Phi-\mathbb E_{\rho_t}\Phi).
\tag{P26}
\]
它是高维 weighted Poisson problem；不需要已知 normalizer 的数值，却需要全局 posterior expectation/particles 与准确 PDE 解。把概率问题改写成 PDE 并未消除统计困难。

### 8.4 多 chart、边界与流形支持

chart \(a\mapsto b\) 转换时
\[
p_b(b)=p_a(a(b))|\det Da/Db|.
\]
intrinsic density \(\rho\) 对合法几何变换不变。覆盖有重叠时可使用 partition of unity \(w_k(q)\ge0,\sum_k w_k(q)=1\)，把 chart index 作为扩展随机变量，目标设为
\(\pi(q,k)=\pi(q)w_k(q)\)。移动 chart 必须带正确权重和坐标 Jacobian。

数据依赖的少数 trajectories 即使都留在高 likelihood 区域，也不决定各自 tube 的质量。精确低维曲面上的分布与全维 posterior 奇异；必须明确它是条件曲面模型，还是有 normal/strong conditional distribution 的 tube approximation。

### 8.5 概率 PDE 方法比较

|方法|是否需已知 posterior \(Z\)|拒绝/误差来源|本问题最主要代价与边界|
|---|---|---|---|
|Euclidean/metric Langevin|不需|连续不变；离散有偏或需 MH|Maxwell gradient，位置 metric/divergence|
|Riemannian MALA/HMC|不需|MH 可修正数值积分；需正确 reverse/volume|动态 J、metric、隐式或约束积分|
|Fokker–Planck 直接网格|需质量归一到1，不需预先给 Z|空间/时间离散误差|2D 可行，4D 已重，8D grid 通常不可行|
|Probability-flow ODE|score 不需 Z|近似 score 与 ODE 数值偏差|需要所有时间的 density score，非只终点梯度|
|Riemannian score diffusion|通常用目标 samples，不直接需 Z|teacher bias、score bias、离散误差|目标 samples/噪声 kernel；边界/退化 geometry|
|Flow matching|用正规化目标 samples 或正确加权近似；不直接需 Z|teacher 与 vector-field approximation|无 teacher 时仍需先解决 inference；复杂曲面 premetric|
|Schrödinger bridge|可用样本/相对密度，无需把 Z 作为输入|桥势、IPF、网络与数值误差|未知时空控制，比 posterior gradient 更难|
|Annealed weighted-Poisson transport|Z 导数以 expectation 替代|particle/Poisson/discretization error|可与 SMC weights 校正，但高维算子求解昂贵|

Riemannian score modelling、general-geometries flow matching、Riemannian Schrödinger bridge 都已有论文。只有把 Maxwell 可认证局部 chart/metric、original prior 与 strong marginalization 接起来，并证明对物理费用或风险有实质改进，才值得作为新的具体方法。

## 9. 两个可执行的物理算法

### 算法 A：相同 MH 核的认证 Maxwell evaluator

~~~
输入：已声明目标 π 或 π_C，固定 proposal Q，已付费 anchor/cache
对每一步：
  1. 从 Q(q'|q) 生成提案；记录正反向 log density。
  2. 只抽一次 u ~ Uniform(0,1)。
  3. 用已有 current/adjoint spaces 计算 Φ(q), Φ(q') 的认证区间。
  4. 若 log u 在 log-ratio 区间外，完成精确 accept/reject。
  5. 否则比较：一个 primal RHS、一个 dual RHS、ROM enrichment、
     Krylov 加精或 Full 的预期区间缩小量/费用。
  6. 选择一项加精，返回步骤 3；到预算边界时调用 Full fallback。
  7. 接受或拒绝后保存缓存；统计目标与本次 proposal 不变。
输出：与指定 Full-MH 耦合一致的 posterior chain、成本账本。
~~~

核心机制：paired residual 决定 likelihood uncertainty，而 MH 随机阈值决定本次需要多高精度。最小决定性实验不是“发现几个新材料”，而是至少数万共同提案上的逐决策一致性、零认证失误、fallback 比例、ESS/总 wall。成本不胜 Full 时应停止 evaluator 分支。

### 算法 B：强观测中心图 + 完整概率 tube

~~~
1. 从 physics-only q0 及其付费 J0 构造固定 Vw, Vs, Cs。
2. 在 low-dimensional a 上用 P6 continuation 构造中心 β(a)。
3. 不丢弃强变量：q(a,u)=q0+Vw a+Vs[β(a)+u]。
4. 用条件 Laplace 作强方向 proposal/importance distribution；
   需要精确性时对 u 做条件 MH，或对近似边缘做纠正。
5. 在 a 上使用简单 pCN/MALA/low-dimensional quadrature；
   evaluator 可使用算法 A。
6. 汇总 E[g(q)|a,y]，而不是只用中心材料 g(q(a,0))。
7. 通过 enlarged charts / full posterior reference 审计漏质量与边界偏差。
~~~

该中心图不是新的 prior，也不是 full fiber。其可利用性在于：曲率仅用 F''，P5 的两项可解释何时强观测补偿需要改变；tube 的概率 Jacobian 明确，强方向并未冻结。论文价值取决于它是否胜过普通 fixed linear decomposition 与条件 Laplace/Gibbs；若中心图构建费用较高、strong conditional 宽度并不小，则没有优势。

### 9.1 算法 B 的正反提案、边界与概率密度

令 \(T=[V_w,V_s]\in\mathbb R^{66\times66}\) 可逆，维数为
\(d+r=66\)，且在一次统计采样 epoch 内固定。令 \(x=(a,u)\)，
\[
\Xi(a,u)=q_0+V_wa+V_s\{\beta(a)+u\}.
\]
从 \(q\) 恢复 \(x\) 的步骤是：先算
\((a,b)=T^{-1}(q-q_0)\)，再令 \(u=b-\beta(a)\)。只要 \(\beta(a)\) 是同一个确定性函数，\(\Xi\) 是一一对应的全维坐标变换；不要求 \(\beta\) 真是条件 MAP，也不要求物理 correction 完美。

\[
D\Xi=T
\begin{pmatrix}I_d&0\\D\beta&I_r\end{pmatrix},\qquad
|\det D\Xi|=|\det T|.
\]

因此相对于 \(da\,du\) 的精确目标为
\[
\bar\pi(a,u)\propto
p_0(\Xi(a,u))e^{-\Phi(\Xi(a,u))}
1_{\mathcal C}(\Xi(a,u))1_A(a)\,|\det T|.\tag{P27}
\]

最后两个指标只有在研究者事先声明这些局部限制时才加入。若 \(T\) 混合了原 Gaussian 与 uniform 坐标，\(p_0(\Xi(a,u))\) 一般不能分解为独立的 \(p(a)p(u)\)。原 uniform2 边界成为
\(0\le E_{\rm prop}\Xi(a,u)\le1\) 的弯曲边界，不能在新的 \(a,u\) 轴上盲目做逐坐标反射。

一个完整且最稳妥的 baseline 是：

1. 固定 full-rank covariance \(\Sigma\succ0\)，在 \(x\) 空间提案
   \(x'=x+\sqrt h\,\Sigma^{1/2}\xi\)，\(\xi\sim N(0,I_{66})\)。
2. 由同一 \(\beta\) 计算 \(q'=\Xi(x')\)。如果超出原 uniform bounds 或已声明局部支撑，直接拒绝；边界拒绝不需要人为截断 Gaussian proposal 的正规化常数。
3. 支撑内的对称 proposal 接受率是
   \(1\wedge[p_0(q')e^{-\Phi(q')}]/[p_0(q)e^{-\Phi(q)}]\)，因为 \(|\det T|\) 常数相消。
4. 可以分块更新 \(a\) 或 \(u\)，但两块都必须在完整 joint target 中更新；只固定 \(u=0\) 将产生另一个目标。

动态 full-rank mobility 也可合法。设
\[
Q_x(x'|x)=N\!\left(
x+h\Sigma(x)\nabla_x\log\bar\pi(x),\ 2h\Sigma(x)\right).
\]
用真正的正反向 Gaussian 密度
\[
\alpha(x,x')=1\wedge
\frac{\bar\pi(x')Q_x(x|x')}{\bar\pi(x)Q_x(x'|x)}
\tag{P28}
\]
即可。这里不必为保持 MH 正确性额外计算 \(\operatorname{div}\Sigma\)，但反向均值、协方差及其 determinant 不能省略。若使用算法 A 认证 \(\Phi\) 而 proposal 还依赖近似 gradient/metric，必须把**实际用于生成 proposal 的近似函数**冻结并用于正反向 \(Q_x\) 计算；否则 Full target acceptance 并不能修复 proposal bookkeeping。

图中心的数值实现也有一个隐蔽风险：对同一 \(a\) 若以不同历史初值求 correction，可能选到不同 \(\beta(a)\)。最初实现应采用冻结的插值 chart、确定性 continuation 规则或明确带 branch index 的 extended state。不能在接受率中默认两次调用对应同一坐标变换。

如果改用随位置变化的 normal coordinates
\(q=\psi(a)+N(a)u\)，则全维 Jacobian 通常不再常数，必须加入
\(|\det D_{a,u}q|\)。在 \(u=0\)、Euclidean orthonormal normal frame 下，该因子含 \(\sqrt{\det(D\psi^\top D\psi)}\)；这不是任意补一个材料 Gram 的理由。vertical tube 的三角结构是推荐它作为首个概率实现的实际原因。

## 10. 建议的决定性实验与拒绝条件

本节仅为未来设计，不执行也不自动授权新付费实验。

1. **缓存认证台架。** 从 P1 新独立观测建立合法 anchor，冻结 3–4 类 proposal，保存所有提案及随机数。用 Full 生成统一审计 reference；在线 controller 不准读取 reference error。比较 P8、P10、Neumann bound 与普通 residual bound 的覆盖率、效用和费用。
2. **局部半径预测。** 沿相同有限材料路径，比较 decoder curvature、feedback curvature、合并 bound 对 J 漂移、谱旋转和 log-likelihood error 的预测能力。应包含 saturation 接缝附近与 Maxwell amplification 大的困难对象。
3. **概率核审计。** 认证 evaluator 必须与 Full-MH 接受序列逐步相同；DA 单独比较，因为其核不同。先冻结 proposal，避免把 proposal adaptation 与 evaluator acceleration 混淆。
4. **物理规模检验。** TM24 先作为失败也能发表的正确性台架；当且仅当成本不等式显示可扩展可能，才增加网格或进入 vector Maxwell。大规模比较必须采用优化 Full/Krylov preconditioner 与共同 cache 权利。
5. **材料风险终点。** 用相同固定样本费用报告 posterior mean 的材料 MSE、MC standard error、可计算统计量 ESS/s。几何方法只减少相关性/算时，不应改变已声明 posterior；若同参考下 risk 差异超出 MC 误差，应先检查偏差。

拒绝条件：certificate 太松导致几乎次次 Full；setup amortization 在合理样本量内无法收回；只减少 RHS 却增加 wall；local restriction 漏质量主导误差；mean headroom 太小；复杂概率 PDE 不比简单 blocked sampler 更准确或更便宜。

## 11. 创新性排序：本分支的独立论文潜力

### 优先 1：Maxwell posterior 决策误差与计算复用

可研究的组合是：有限 feedback + 共同材料约束 + primal/dual current spaces + log-ratio 认证 + 缓存经济学。真正待证明/待测量的内容是 tightness 与 break-even，以及如何按后验决策分布选择空间。P1–P19 的多数数学工具都是标准；把它们连接起来仍必须产生此前没有的 Maxwell 可检验机制，才能支撑强 TAP。

### 优先 2：强观测补偿图与 conditional tube 的 Maxwell 曲率控制

优点是清楚地区分“曲面中心”和“强方向概率厚度”，并在 C² decoder 条件下有可证明的局部结构。与一般 nonlinear coordinate transform、IFT/continuation 重合很高；若无法用 feedback 分解得到明显更廉价的 chart transport，不应独立包装成新几何学。

### 暂缓：完整流形 diffusion / Schrödinger bridge

当前最大困难是正确 local target、strong marginalization、reference 与 cache 费用。复杂概率 PDE 同时要求未知 geometry 和未知时间密度，增加了两个尚未解决的学习问题。可以保留 2D/4D 的概念验证或几何可视化；在 reference 之前训练大网络不合适。

## 12. 原始文献与本轮核验索引

以下为原始论文/作者或出版者来源；摘要用于创新边界，未对所有论文全文逐定理复核。文献使用永久 URL。

1. **Cui, Martin, Marzouk, Solonen, Spantini (2014), Likelihood-informed dimension reduction for nonlinear inverse problems.** [arXiv:1403.4680](https://arxiv.org/abs/1403.4680). 已有 prior/likelihood 分解、LIS complement prior、Rao–Blackwellization；不等于冻结弱/强变量。
2. **Cui, Law, Marzouk (2016), Dimension-independent likelihood-informed MCMC, JCP 304:109–137.** [arXiv:1411.3688](https://arxiv.org/abs/1411.3688), DOI 10.1016/j.jcp.2015.10.008。 operator-weighted proposals 与 function-space reference；单一 anchor weak spectrum 不是 posterior-global LIS。
3. **Cotter, Roberts, Stuart, White (2013), MCMC Methods for Functions: Modifying Old Algorithms to Make Them Faster, Statistical Science 28(3):424–446.** [arXiv:1202.0709](https://arxiv.org/abs/1202.0709), DOI 10.1214/13-STS421。 pCN 的 mesh robustness 已成熟。
4. **Zahm, Cui, Law, Spantini, Marzouk (2022), Certified dimension reduction in nonlinear Bayesian inverse problems, Mathematics of Computation.** [arXiv:1807.03712](https://arxiv.org/abs/1807.03712), DOI 10.1090/mcom/3737。 non-Gaussian prior/likelihood 及 subspace log-Sobolev error certificate，防止把普通弱谱误差重新命名。
5. **Stuart (2010), Inverse problems: A Bayesian perspective, Acta Numerica 19:451–559.** [作者存档](https://authors.library.caltech.edu/records/fk9j9-5tq15), DOI 10.1017/S0962492910000061。 函数空间 Bayesian inverse problems 的 well-posedness/approximation 基础；本文 P16–P17 是有限局部支撑上的自包含 elementary bound。
6. **Lieberman, Willcox, Ghattas (2010), Parameter and State Model Reduction for Large-Scale Statistical Inverse Problems, SISC 32:2523–2542.** [作者 PDF](https://acdl.mit.edu/LiebermanSISC2010.pdf), DOI 10.1137/090775622。 参数与 PDE state 同时降维已长期存在。
7. **Cui, Marzouk, Willcox (2016), Scalable posterior approximations for large-scale Bayesian inverse problems via likelihood-informed parameter and state reduction, JCP 315:363–387.** [arXiv:1510.06053](https://arxiv.org/abs/1510.06053), DOI 10.1016/j.jcp.2016.03.055。 是“weak/strong + ROM + posterior”的直接先例。
8. **Spantini, Cui, Willcox, Tenorio, Marzouk (2017), Goal-Oriented Optimal Approximations of Bayesian Linear Inverse Problems, SISC 39(5).** [arXiv:1607.01881](https://arxiv.org/abs/1607.01881), DOI 10.1137/16M1082123。 已有 task posterior covariance/mean 与 weighted Bayes risk 的最优低秩理论；本任务非线性 decoder 风险不能直接套用。
9. **Christen and Fox (2005), Markov chain Monte Carlo Using an Approximation, JCGS 14(4):795–810.** [作者 PDF](https://physics.otago.ac.nz/data/fox/publications/ApproxMCMCPreprint.pdf). DA、局部 forward approximation 已成熟；P15 应承认来源。
10. **Conrad, Marzouk, Pillai, Smith (2016), Accelerating Asymptotically Exact MCMC for Computationally Intensive Models via Local Approximations, JASA 111(516):1591–1607.** [arXiv:1402.1694](https://arxiv.org/abs/1402.1694). 有专门随机 refinement/ergodicity 机制，不能推广为任意 adaptive surrogate。
11. **Chen and Ghattas, Stein variational reduced basis Bayesian inversion.** [arXiv:2002.10924](https://arxiv.org/abs/2002.10924). goal-oriented model reduction、potential/gradient error、posterior KL 与 transport sample error 已联结。
12. **Cao, O’Leary-Roseberry, Ghattas (2024), Derivative-informed neural operator acceleration of geometric MCMC for infinite-dimensional Bayesian inverse problems.** [arXiv:2403.08220](https://arxiv.org/abs/2403.08220); [全文](https://arxiv.org/html/2403.08220v2). Sections 2.8、4.3、5、Appendix C 直接覆盖 DA、训练 RHS 成本、posterior geometry 与 misfit approximation；因此 neural surrogate + geometry + DA 不是本轮新点。
13. **Girolami and Calderhead (2011), Riemann manifold Langevin and Hamiltonian Monte Carlo methods, JRSSB 73(2):123–214.** [出版者](https://academic.oup.com/jrsssb/article-abstract/73/2/123/7034367), DOI 10.1111/j.1467-9868.2010.00765.x。 位置依赖 metric/MCMC 是成熟方法。
14. **Xifara, Sherlock, Livingstone, Byrne, Girolami (2014), Langevin diffusions and the Metropolis-adjusted Langevin algorithm.** [arXiv:1309.2983](https://arxiv.org/abs/1309.2983). 直接讨论 Riemannian volume 与 Lebesgue invariant density 的差异；P20–P24 是关键参考。
15. **Ma, Chen, Fox (2015), A Complete Recipe for Stochastic Gradient MCMC, NeurIPS.** [arXiv:1506.04696](https://arxiv.org/abs/1506.04696). position-dependent diffusion 的 divergence correction 非 Maxwell 特有。
16. **Lelièvre, Rousset, Stoltz (2012), Langevin dynamics with constraints and computation of free energy differences, Mathematics of Computation.** [arXiv:1006.4914](https://arxiv.org/abs/1006.4914), DOI 10.1090/S0025-5718-2012-02594-4。 约束曲面 canonical measure、Metropolis-corrected integrators/free-energy 先例。
17. **Song et al. (2021), Score-Based Generative Modeling through Stochastic Differential Equations, ICLR.** [arXiv:2011.13456](https://arxiv.org/abs/2011.13456). reverse SDE/probability-flow ODE 所需的是 time-dependent perturbed-density score。
18. **De Bortoli et al. (2022), Riemannian Score-Based Generative Modelling, NeurIPS.** [arXiv:2202.02763](https://arxiv.org/abs/2202.02763). 流形 score diffusion 已有。
19. **Chen and Lipman (2024), Flow Matching on General Geometries, ICLR.** [arXiv:2302.03660](https://arxiv.org/abs/2302.03660); [官方会议页](https://proceedings.iclr.cc/paper_files/paper/2024/hash/d1f9936d3be6997ffffab692977eebe6-Abstract-Conference.html). premetric/geodesic/spectral-distance flow matching 已有，复杂 geometry 的预处理费用需计入。
20. **De Bortoli, Thornton, Heng, Doucet (2021), Diffusion Schrödinger Bridge with Applications to Score-Based Generative Modeling, NeurIPS.** [arXiv:2106.01357](https://arxiv.org/abs/2106.01357). 桥问题与 IPF/score 学习已成熟。
21. **Levy, Hoffman, Sohl-Dickstein (2018), Generalizing Hamiltonian Monte Carlo with Neural Networks, ICLR.** [arXiv:1711.09268](https://arxiv.org/abs/1711.09268). 可逆 learned MCMC 提案不构成本轮独立创新。
22. **Wang, Chen, Kanagawa, Oates (2025), Reinforcement Learning for Adaptive MCMC, AISTATS 258:640–648.** [官方论文页](https://proceedings.mlr.press/v258/wang25b.html). RL 设计 MH kernel 有理论化先例。
23. **Wang, Fisher, Kanagawa, Chen, Oates (2026), Harnessing the Power of Reinforcement Learning for Adaptive MCMC, AISTATS 300:631–639.** [官方论文页](https://proceedings.mlr.press/v300/wang26a.html); [arXiv:2507.00671](https://arxiv.org/abs/2507.00671). 其研究报告 acceptance 与 ESJD 奖励对 RLMH 训练未必提供足够信号，提出 contrastive-divergence 奖励。P3-L 应把物理预算决策当可能增量，而非声称首次 RL+MH。
24. **Roberts and Rosenthal (2007), Coupling and Ergodicity of Adaptive Markov Chain Monte Carlo Algorithms, Journal of Applied Probability 44(2):458–475.** [出版者](https://www.cambridge.org/core/journals/journal-of-applied-probability/article/coupling-and-ergodicity-of-adaptive-markov-chain-monte-carlo-algorithms/91713CBCA261758088FCDFBF8DB43FE0), DOI 10.1239/jap/1183667414。 diminishing adaptation 与 containment 的经典依据。
25. **Miller et al. (2016), Fundamental limits to optical response in absorptive systems, Optics Express 24:3329–3364.** [arXiv:1503.03781](https://arxiv.org/abs/1503.03781); [作者 PDF](https://www.mit.edu/~soljacic/limits_OE.pdf), DOI 10.1364/OE.24.003329。 Eq. (5)–(7) 与 (23) 支持有耗材料的 extinction/scattering/absorption 区分与响应界；无关 posterior 权重。
26. **Ola, Päivärinta, Somersalo (1993), An inverse boundary value problem in electrodynamics, Duke Mathematical Journal 70(3):617–653.** [官方期刊卷](https://projecteuclid.org/journals/duke-mathematical-journal/volume-70/issue-3), DOI 10.1215/S0012-7094-93-07014-7。 完整边界数据类唯一性与有限测量数值满秩是不同命题。
27. **Ola, Päivärinta, Somersalo (2003), Inverse Problems for Time Harmonic Electrodynamics, MSRI Publications 47.** [作者章节 PDF](https://library.slmath.org/books/Book47/files/ola.pdf). 作者总结 Maxwell 材料系数逆边界与散射的条件；可据此进一步查相应原始唯一性定理，不能直接外推到 P1 acquisition。

