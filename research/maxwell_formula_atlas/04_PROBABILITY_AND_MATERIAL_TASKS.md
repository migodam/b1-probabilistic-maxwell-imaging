# 04 概率、原材料先验与任务积分的公式集

本页是可持续增补的推导目录，编号 `F-PROB-xxx` 保持稳定。它整理原 P1/R5/R6 与 R8P/R9 实际使用的表达式、条件、源码位置和当前资格；不新增计算、物理实验或公式级微型测试。下文“精确”默认指给定离散矩阵及明确 prior 下的数学恒等式，不自动包含浮点舍入、Hankel生成、离散化或连续Maxwell误差。Full66始终为 `FULL66_UNQUALIFIED`。

记 \(\varphi,\Phi_N\) 为标准正态密度/CDF；用 \(\mathcal V\) 或 \(v\) 表示复材料幅度，避免与CDF混名。\(G_{\rm dec}\) 为材料decoder，\(K\) 为离散Green算子，\(F\) 为完整8源32接收forward。所有任务积分都针对同一个固定观测 \(y\)，除非另立模型合同。

## F-PROB-001｜原圆对称复噪声、白化与likelihood

**合同。** 256个复测量的噪声独立满足 \(\eta_k\sim\mathcal{CN}(0,\sigma_c^2)\)，即 \(E|\eta_k|^2=\sigma_c^2\)，实、虚部方差各为 \(\sigma_c^2/2\)。原值为 \(\sigma_c=0.00119211798514519\)。令 \(\operatorname{pack}\) 逐源排列Re32、Im32，则

\[
f(z)=\frac{\sqrt2}{\sigma_c}\operatorname{pack}F(G_{\rm dec}(z)),\quad y_w=\frac{\sqrt2}{\sigma_c}\operatorname{pack}y,
\quad \mathcal P(z)=\frac12\|f(z)-y_w\|^2=\frac{\|F(G_{\rm dec}(z))-y\|_F^2}{\sigma_c^2},
\quad L(z)=e^{-\mathcal P(z)}.
\]

**推导。** 每个复密度为 \((\pi\sigma_c^2)^{-1}e^{-|\eta|^2/\sigma_c^2}\)。乘256个测量，负log的材料相关项恰为上式。白化后的512实分量协方差为单位阵。删去常数 \((\pi\sigma_c^2)^{-256}\) 不改变同模型posterior；跨噪声/采集/model evidence比较必须恢复常数。Z-ALG“每实/虚σ=.001”是另一个合同，不能代入本式。

**源码。** [R4 `Target.phi/white`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport/src/p3r4/core.py:236)；[R6 `Target.white/phi`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/shared/physics.py:143)；[R9 `pack`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_transfer.py:44)。

## F-PROB-002｜原Gaussian64×Uniform2与Gaussianized66

原native坐标 \(q=(z_{0:63},u_R,u_I)\)，\(z_{0:63}\sim N(0,I_{64})\)，\(u_R,u_I\stackrel{\rm iid}{\sim}U(0,1)\)。令 \(x_{64}=\Phi_N^{-1}(u_R),x_{65}=\Phi_N^{-1}(u_I)\)，则 **完整原先验** \(x\sim N(0,I_{66})\)：因为 \(du_j=\varphi(x_j)dx_j\)，两项Jacobian正好补成两个Gaussian密度。

原decoder为

\[
s_i=\frac{b_i+A_i x_{0:63}}{\theta},\quad Q_i=H(s_i)1_{i\in ROI},\quad
D(x)=.8+\Phi_N(x_{64})+i\{.03+.15\Phi_N(x_{65})\},\quad
m_i=.1+.04i+D(x)Q_i.
\]

这里 \(\theta=.5\)，offset列 \(A_{i0}=.65\)。\(H\) 是原C² quintic saturation，不能换成sigmoid/clip后仍声称原prior相同。native Uniform边界对应无穷probit坐标；仅初始化时的向内clip不是精确概率变换。

**源码。** [R4 `to_probit/from_probit/probit_derivative`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport/src/p3r4/core.py:294)；[原 `MaterialPrior.decode/derivative`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport/src/p1_material/prior.py:50)；[`saturation`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport/src/p1_material/prior.py:10)。

## F-PROB-003｜材料推前、evidence与同材料条件核

\[
\mu_0=(G_{\rm dec})_\#\pi_0,\quad
Z_y=\int L_y(m)\mu_0(dm)>0,\quad
\pi_M(dm\mid y)=Z_y^{-1}L_y(m)\mu_0(dm).
\]

对任何有界可测材料任务 \(\tau\)，

\[
E[\tau(M)\mid y]=\frac{\int\tau(G_{\rm dec}(x))L_y(G_{\rm dec}(x))\pi_0(dx)}{\int L_y(G_{\rm dec}(x))\pi_0(dx)}.
\]

若 \(\kappa_0(dx\mid m)\) 是原prior对decoder的disintegration，由 \(L_y\) 仅依赖 \(m\) 可得 \(\pi(dx\mid m,y)=\kappa_0(dx\mid m)\)：数据不会在**精确同材料**gauge内重新赋权。相反，不同材料即使同forward，仍按各自材料prior质量分配；finite-noise近邻关系不能当作精确等价关系。

**源码/理论。** [R8P measure章节](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/02_measure.tex:5)；R9实现为[ `summary`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:12)，它计算有限求积的分子/分母，未实现Full66 disintegration。

## F-PROB-004｜平方损失的Bayes动作与固定数据风险

对实任务 \(T=\tau(M)\)，若 \(E[T^2\mid y]<\infty\)，令 \(\mu_y=E[T\mid y]\)，则

\[
R_y(a)=E[(a-T)^2\mid y]=\operatorname{Var}(T\mid y)+(a-\mu_y)^2,
\quad a^*=\mu_y.
\]

**推导。** 展开 \(a-T=(a-\mu_y)+(\mu_y-T)\)，条件交叉项因 \(E[T-\mu_y\mid y]=0\) 消失。复材料向量实化后，对半正定损失矩阵 \(W\)，同理

\[
R_y(a)=\operatorname{tr}(W\Sigma_y)+\|a-\mu_y\|_W^2.
\]

所以更精确计算动作只减少第二项，固定同一个 \(y\) 的posterior variance不因计算轮数下降。若均值区间半宽≤ε，中点动作的excess risk≤ε²。

**源码/理论。** [R8P square loss与停止](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/06_algorithm_cost.tex:35)；R6已实现的有限VAL恒等式是另一个对象，见F-PROB-021与[ `track_b/report.py`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/track_b/report.py:24)。

## F-PROB-005｜有限特征Bayes动作与信息损失

若 \(S=S(Y)\)，总体平方风险最优动作分别是 \(a_Y=E[T\mid Y]\)、\(a_S=E[T\mid S]\)。条件期望正交性给

\[
E\|T-a_S\|^2=E\|T-a_Y\|^2+E\|a_Y-a_S\|^2.
\]

进一步对任何只依赖S的学习动作 \(\hat a(S)\)，

\[
E\|T-\hat a(S)\|^2=E\|T-a_Y\|^2+E\|a_Y-a_S\|^2+E\|a_S-\hat a(S)\|^2.
\]

**条件与边界。** S必须是Y的可测函数、二阶矩有限。有限已暴露VAL上的teacher/student损失，不等于这里未知的population conditional mean。Gaussian四通道可作为任务有用特征，但特征有用不能反推其生成Gaussian是校准posterior。

**来源。** [R5任务理论说明](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/THEORY_AND_NUMERICAL_CHECKS.md:27)；这是解释学习目标的恒等式，R5/R6未计算未知 \(a_Y,a_S\) 的真实总体分解。

## F-PROB-006｜固定正交Full66坐标与anchor平移

设 \(Q=[Q_s,Q_w]\in O(66)\) **固定**。若 \(x=Q_s s+Q_w w\)，原 \(x\sim N(0,I)\) 等价于独立 \(s\sim N(0,I),w\sim N(0,I)\)。若为了局部几何改成

\[
x=x_*+Q_s s+Q_w w,
\]

则

\[
(s,w)=Q^T(x-x_*)\sim N(-Q^Tx_*,I),\qquad
\rho(s,w)\propto\exp\{-\tfrac12\|x_*+Q_s s+Q_w w\|^2\}.
\]

**推导。** 正交变换保留范数和Lebesgue体积；平移改变均值。忽略anchor-induced线性prior项会丢失prior tilt。固定 \(w=w_0\) 的slice是正则条件模型，不能代替对w的prior conditional积分，也不等于Full66材料边缘。Gaussian盒截断时需另保留截断归一化与domain；旋转后的原盒一般不再是独立轴盒。

**来源。** [R8P正交坐标原则](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/02_measure.tex:54)；[R5 anchor中心响应坐标约定](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/THEORY_AND_NUMERICAL_CHECKS.md:5)。固定正交与移动校正frame的区别亦见同目录03章；本式未被R9用来生成Full66 posterior。

## F-PROB-007｜非正交lift与非线性chart的prior

若满秩线性chart \(x=x_*+T\xi\) 且T为方阵可逆，则

\[
\rho_\xi(\xi)=\varphi_{66}(x_*+T\xi)|\det T|,
\quad \text{precision}=T^TT,\quad E\xi=-T^{-1}x_*.
\]

薄矩阵T只描述子空间slice；可与补空间组成完整变换后再条件化。材料lift \(T=V\Sigma^{-1}\) 通常不正交，不能赋予iid标准Gaussian系数。对可逆非线性 \(x=\Psi(s,w)\)，必须使用 \(\rho=\varphi_{66}(\Psi)|\det D\Psi|\)。状态依赖的移动frame \(Q(x)\) 不自动定义一个全域单射变换，不能逐点写det=1就略去全局映射。

**来源。** [R8P chart Jacobian](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/02_measure.tex:58)；实际实现资格只在明确固定chart内成立；一般全域非线性chart在R9 `NOT_RUN`。

## F-PROB-008｜条件积分不是profile：截面体积、Hessian与分支

固定弱坐标w，令含原prior/Jacobian的联合密度为ρ，则

\[
Z(w)=\int e^{-\mathcal P(s,w)}\rho(s,w)ds,
\quad N_\tau(w)=\int\tau(G_{\rm dec}(s,w))e^{-\mathcal P(s,w)}\rho(s,w)ds.
\]

只有对全部w/源层求和再取 \(N/Z\) 才是任务。\(\mathcal P_{\rm prof}(w)=\min_s\mathcal P(s,w)\) 缺少截面宽度、prior tilt和其他极小点。若s维数r，且每个孤立内部峰 \(s_k(w)\) 的likelihood Hessian \(H_k\succ0\)，Laplace主项是

\[
Z(w)\approx\sum_k e^{-\mathcal P(s_k,w)}\rho(s_k,w)(2\pi)^{r/2}|H_k(w)|^{-1/2}.
\]

若对 \(\mathcal P-\log\rho\) 做MAP展开，应改用该势能的Hessian并避免重复乘ρ。边界截断峰、seam、奇异Hessian需另处理；要成为误差证书还须控制局部余项、峰外evidence及遗漏峰，不能仅比较局部GH阶数。

**来源/实现。** [R8P条件积分](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/02_measure.tex:66)；[R9 fresh mode/GN](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:18)与[ `gh`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:27)使用真实likelihood校正，未把Laplace当oracle；[1D局部Laplace](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:45)只标局部近似。

## F-PROB-009｜响应坐标的体积因子

为避免把原噪声再除一次，使用已白化f。固定w及固定输出正交frameU，设 \(h=U^T(f(s,w)-y_w)\)。若 \(s\mapsto h\) 在第k分支上单射可微，且 \(B_k=U^TD_s f\) 可逆，则

\[
ds=|\det B_k|^{-1}dh,
\quad Z(w)=\sum_k\int e^{-\frac12\|h\|^2-\frac12\|P_{U^\perp}(f-y_w)\|^2}\,
\frac{\rho(s_k(h,w),w)}{|\det B_k|}dh.
\]

若使用未白化输出r及 \(h=U^T(r-y)/\sigma\)，则对应 \(\sigma^r/|\det(U^TD_s r)|\)。这里r为实噪声尺度所在合同，不能把σc与每实分量σ混用。多个逆分支必须相加；chart重叠需原域不交分割/partition of unity，不能重复同一prior质量。有限求根成功不是全域单射或覆盖证明。

**来源。** [R5响应坐标原式](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/THEORY_AND_NUMERICAL_CHECKS.md:5)。R9主线使用原source/材料幅度积分，没有新建非线性响应atlas。

## F-PROB-010｜材料quotient/coarea与混合秩参考测度

在一个d维原latent源层 \(S\) 上，若 \(G_{\rm dec}\) 为常秩r并满足coarea适用条件，设 \(J_rG\) 为非零奇异值之积，则相对于材料像的 \(r\) 维Hausdorff测度，

\[
\frac{d\mu_{0,S}}{d\mathcal H^r}(m)
=\int_{S\cap G^{-1}(m)}\frac{p_0(x)}{J_rG(x)}d\mathcal H^{d-r}(x).
\]

**推导。** 对coarea公式选择被积函数 \(p_0/J_rG\)，先沿原像纤维积分，再沿材料像积分。材料图的几何体积 \(\sqrt{\det(Dm^TDm)}\) 本身不是原材料prior；还需原latent密度与所有原像。

不同秩的源层分别推送并以**测度**相加，不能把不同维density作为同单位普通函数相加。最稳妥的通用表示为 \(\mu_0=\sum_jG_\#(1_{S_j}\pi_0)\)，原 \(S_j\) 不重叠。材料像重叠时各原prior贡献相加；同一个原prior区域被多个chart重复积分才是错误。

**来源。** [R8P coarea及混合层](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/02_measure.tex:62)、[同节源层公式](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/02_measure.tex:122)；[R6 `source_nodes/quotient_nodes`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/track_a/stratified.py:71)。

## F-PROB-011｜host原子、pure层与Full66的区别

设 \(H(s)=0\) 于 \(s\le-1\)、\(H(s)=1\) 于 \(s\ge1\)。原Full66中

\[
p_h=\Pr_0\{s_i\le-1\;\forall i\in ROI\},\quad
\mu_0\supset p_h\delta_{m_h},\quad
\Pr(M=m_h\mid y)=p_h L(m_h)/Z_y.
\]

这是一个高维Gaussian多面体概率，本轮没有计算。全部ROI纯相时若两个属性仍活动，材料像是二维连续层而不是单个原子；只有属性也固定的条件slice才把pure平台变为单材料原子。Full66不能套用uniform slice的单变量PH。finite-noise下若forward有限，则L严格正；浮点下溢0不能证明任一正prior平台的posterior质量为0。

**来源/实现。** [原 `MaterialPrior.decode`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport/src/p1_material/prior.py:50)；R6 uniform平台质量在[ `stratified.py`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/track_a/stratified.py:18)；R9下溢修复在[1D atom bounds](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:17)。

## F-PROB-012｜R6复幅度V=QD的正确联合prior与Q⁻² fibre

R6 S3明确采用三个活动Gaussianized坐标独立截断[-4,4]，其余场坐标0。记

\[
Z_0=\Phi_N(4)-\Phi_N(-4),\ A=1/1.3,\ U_- =\Phi_N(-4),\ U_+=\Phi_N(4),
\quad P_H=\frac{\Phi_N(-A)-\Phi_N(-4)}{Z_0},
\]

\[
Q(z)=H(1.3z),\quad R_D\in[R_-,R_+]=[.8+U_-,.8+U_+],\quad
I_D\in[I_-,I_+]=[.03+.15U_-,.03+.15U_+].
\]

由于 \(du=\varphi(x)dx\)，截断probit属性推成这些区间的**归一化Uniform**。在transition \(-A<z<A\)，联合prior为 \(\varphi(z)dz\,dR_DdI_D/(Z_0\Delta R\Delta I)\)。固定Q的变换 \((R_D,I_D)\mapsto(R,I)=Q(R_D,I_D)\) 的实Jacobian是 **Q²**。所以

\[
\rho_{\rm tr}(R,I)=\frac1{\Delta R\Delta I}
\int_{z_-}^{z_+}\frac{\varphi(z)}{Z_0Q(z)^2}dz,
\]

\[
q_- =\max(R/R_+,I/I_+),\quad q_+=\min(R/R_-,I/I_-,1),\quad
z_\pm=\frac{2I^{-1}_{q_\pm}(3,3)-1}{1.3}.
\]

当 \(0<q_-<q_+\) 时上式成立，否则transition密度0。pure层密度为 \(\rho_{\rm pure}=P_H1_{[R_-,R_+]\times[I_-,I_+]}/(\Delta R\Delta I)\)，host质量 \(P_H\delta_0\)。总材料prior为

\[
\mu_V=P_H\delta_0+(\rho_{\rm tr}+\rho_{\rm pure})\,dR\,dI.
\]

**源码。** [R6 `fibre`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/track_a/stratified.py:31)；[R9 `fibre/prior`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_transfer.py:50)。\(Q^{-1}\)只对应一维实乘积例，不能替代这里的 \(Q^{-2}\)。

## F-PROB-013｜V固定后的条件依赖、源层权重与copula

在transition条件于 \(V=R+iI\)，

\[
p(z\mid V,\mathrm{tr})\propto\varphi(z)Q(z)^{-2}1_{z_-<z<z_+},
\quad D=V/Q(z).
\]

因此Q与D不再独立。总source标签的条件概率为 \(\Pr(\mathrm{tr}\mid V)=\rho_{\rm tr}/(\rho_{\rm tr}+\rho_{\rm pure})\)，pure同理；这些同材料贡献不应被MERGE删掉。事件 \(Q>.5\) 的transition任务密度只把积分区间改为 \((\max(z_-,0),z_+)\)，再加pure密度。原source坐标直接积分时没有Q⁻²，因为尚未换成材料Lebesgue坐标。

历史实标量控制 \(Q,D\stackrel{\rm iid}{\sim}U[1,2]\)、条件QD=2，则

\[
p(q\mid QD=2)=\frac1{q\log2},\quad D=2/q,\quad
E(QD)=2,\quad \operatorname{Cov}(Q,D)=2-(1/\log2)^2.
\]

若独立生成正确的两个边缘，反而得到 \(E(QD)=(1/\log2)^2\ne2\)。这是已有条件依赖反例，不是本轮新测试。

**源码。** [R6 fibre phase numerator](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/track_a/stratified.py:43)；[R9 `source_integral`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:50)；[既有 `qd_reproduction`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/experiments/r8p_benchmarks.py:166)。

## F-PROB-014｜条件prior平均与真实conditional likelihood的差

固定s，以原conditional prior记 \(\bar L(s)=E_0[L\mid s]\)。直接代入归一化定义可得

\[
E[\tau\mid s,y]-E_0[\tau\mid s]
=\frac{\operatorname{Cov}_0(\tau,L\mid s)}{\bar L(s)}.
\]

若任务跨度为Δτ，且整个条件域有 \(0<\ell\le L\le u\)，则sharp bounded-tilt界给

\[
|E_L\tau-E_0\tau|\le\Delta\tau\frac{\sqrt u-\sqrt\ell}{\sqrt u+\sqrt\ell}
=\Delta\tau\tanh\{\tfrac14\log(u/\ell)\}.
\]

条件kernel不变需要**likelihood对该变量恒定**，不只是任务恒定。任务与弱变量无关时，弱方向likelihood平均仍可能改变强方向的边缘权重，不能固定成一个代表。

**来源。** [R8P条件均值与tilt界](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/02_measure.tex:98)。R9未以单点残差估计替代整个条件域 \(\ell,u\)。

## F-PROB-015｜Gaussian checkpoint返回原坐标后的均值/协方差

若网络坐标 \(\xi\sim N(m_\xi,\operatorname{diag}v_\xi)\)，历史row-vector约定为 \(x=(\xi\odot s+b)F\)，则

\[
m_x=(m_\xi\odot s+b)F,\quad C_x=F^T\operatorname{diag}(v_\xi\odot s^2)F.
\]

diagonal网络输出不意味着原材料坐标独立。对场 \(S_i=(b_i+A_ix_{0:63})/\theta\)，

\[
E S_i=(b_i+A_i m_z)/\theta,\quad
\operatorname{Cov}(S_i,S_j)=A_iC_{zz}A_j^T/\theta^2,\quad
\operatorname{Cov}(S_i,x_{64:65})=A_iC_{z,u}/\theta.
\]

只有在检查到相应交叉块为0后，Gaussian独立性才允许把属性/场期望相乘。这个Gaussian是近似q；其pushforward矩计算准确不等于q是真posterior。

**源码。** [R5 `raw_gaussian`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:73)；[ `diagonal_image`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:81)显式拒绝非零场–属性cross；[ `gaussian_image`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:92)。

## F-PROB-016｜Gaussian probit属性的一阶、二阶与相关矩

若 \((U,V)\) 联合Gaussian，均值μU、μV，方差vU、vV，covariance cUV，则

\[
E\Phi_N(U)=\Phi_N\!\left(\frac{\mu_U}{\sqrt{1+v_U}}\right),
\quad E\Phi_N(U)^2=\Phi_2(t_U,t_U;v_U/(1+v_U)),
\]

\[
E[\Phi_N(U)\Phi_N(V)]=\Phi_2\!\left(t_U,t_V;\frac{c_{UV}}{\sqrt{(1+v_U)(1+v_V)}}\right),\quad t_U=\mu_U/\sqrt{1+v_U}.
\]

**推导。** 引入与(U,V)独立的标准正态E1,E2，把CDF乘积写为事件 \(E_1\le U,E_2\le V\)，再标准化Gaussian差值。对 \(D=a+b\Phi_N(U)\)，\(ED=a+bE\Phi_N(U)\)，\(ED^2=a^2+2abE\Phi_N(U)+b^2E\Phi_N(U)^2\)。

**源码。** [R5 `bvn_same/bvn`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:17)；[ `correlated_cell`条件probit矩](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:67)。二元CDF数值积分本身仍需精度资格。

## F-PROB-017｜相关Gaussian通过quintic decoder的材料矩

对单像素联合Gaussian \((S,U_R,U_I)\)，记 \(v_S>0\)、\(c=\operatorname{Cov}(U,S)\)，则

\[
U\mid S=s\sim N\left(\mu_U+\frac c{v_S}(s-\mu_S),\ C_{UU}-\frac{cc^T}{v_S}\right).
\]

由F-PROB-016得到 \(d_a(s)=E[D_a\mid S=s]\)、\(d_{ab}(s)=E[D_aD_b\mid S=s]\)，并有

\[
\mu_{\Delta m,a}=\int H(s)d_a(s)\varphi_{\mu_S,v_S}(s)ds,
\quad E[\Delta m_a\Delta m_b]=\int H(s)^2d_{ab}(s)\varphi_{\mu_S,v_S}(s)ds.
\]

减去 \(\mu_a\mu_b\) 得实虚2×2 covariance。vS=0时PSD要求c=0，取S确定值极限即可。只有场–属性独立时才简化为 \(E[HD_a]=E[H]E[D_a]\)。

四通道只含Re/Im均值和各自variance；复Hermitian variance为VarR+VarI，pseudo-covariance为VarR−VarI+2iCov(R,I)，二者不同。跨像素相关还需 \((S_i,S_j,U)\) 的联合积分；逐像素四通道不是完整空间posterior sampler。

**源码。** [R5 `correlated_cell`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:58)；[ `gaussian_image`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:92)。

## F-PROB-018｜quintic的Gaussian截断矩与平台

\[
H(s)=\begin{cases}0&s\le-1,\\
\frac12+\frac{15}{16}s-\frac58s^3+\frac3{16}s^5&-1<s<1,\\1&s\ge1.\end{cases}
\]

对S=μ+σZ，设a=(-1−μ)/σ、b=(1−μ)/σ，\(I_k(a,b)=\int_a^b z^k\varphi(z)dz\)。分部积分给

\[
I_0=\Phi_N(b)-\Phi_N(a),\quad I_1=\varphi(a)-\varphi(b),\quad
I_k=a^{k-1}\varphi(a)-b^{k-1}\varphi(b)+(k-1)I_{k-2}.
\]

将 \(H(\mu+\sigma z)=\sum_{k=0}^5c_kz^k\) 展开，则 \(EH(S)=\Phi_N(-b)+\sum c_kI_k\)；H²用10阶矩，平台tail同样只加一次 \(\Phi_N(-b)\)。σ=0取H(μ)。极窄/极尾区间有相消，应使用稳定CDF差或seam-split积分，不能从符号精确推断浮点实现自动严格。

另有 \(H(-1+\epsilon)=\frac54\epsilon^3-\frac{15}{16}\epsilon^4+\frac3{16}\epsilon^5\)，说明decoder接缝的三次起始；它不证明Full66存在物理精确解纤维。

**源码。** [R5 `truncated_phase_moments/phase_moments`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:35)；[R9 1D Gaussian moments](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:48)。

## F-PROB-019｜材料阈值与256-ROI硬事件的prior矩

规定

\[
\tau_{\rm mean}=\frac1{256}\sum_{i\in ROI}\Re m_i,\quad
E=\{\#\{i\in ROI:\Re m_i\ge.75\}\ge64\}.
\]

固定合法正实属性 \(D_R\)，单像素阈值转成 \(H(s_i)\ge q_*=(.75-h_R)/D_R\)。若0<q*<1，\(s_i\ge 2I_{q_*}^{-1}(3,3)-1\)。沿原Gaussian坐标t，若 \(s_i=a_it+b_i\)，每像素至多一个阈值根，ai<0时不等号反向。排序至多256个根后，每段计数恒定，prior质量用真实conditional Gaussian CDF差相加；ai=0或q*越界需单独判定，原子阈值使用≥。

这给 \(\alpha_j=\pi_0(C_j\cap E)\)，而不是posterior事件概率。属性仍活动时外层还要对属性conditional prior积分；不能把逐像素Bernoulli边缘当独立概率相乘。平滑均值很接近不保证硬事件接近。

**来源/实现。** [R8P event/prior slicing](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/03_certificates.tex:67)；[R9 1D event masses](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:76)；[R9 uniform summary](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:12)。

## F-PROB-020｜已使用RRR：predictable residual不是GT covariance PCA

将任务损失权重吸入 \(U=W^{1/2}(M-m_0(Y))\)，用只依赖Y的居中特征ψ，设 \(B=E[\psi\psi^T]\)、\(C=E[U\psi^T]\)。对于rank≤k的修正A，考虑

\[
\min_{\operatorname{rank}A\le k}E\|U-A\psi\|^2+\lambda\|A\|_F^2.
\]

令 \(B_\lambda=B+\lambda I\succ0\)，配方得 \(\|AB_\lambda^{1/2}-CB_\lambda^{-1/2}\|_F^2+\mathrm{const}\)。截断SVD因此给

\[
A_k=[CB_\lambda^{-1/2}]_kB_\lambda^{-1/2}.
\]

λ=0时训练/总体最优线性目标的收益为前k个奇异值平方和；λ>0该谱收益对应惩罚目标，不可直接当未惩罚验证风险改善。代码使用 \(\lambda=\mathrm{ridge}\cdot\operatorname{tr}(B)/p\)，并保留输入/输出均值作为intercept。GT residual PCA只看 \(E[UU^T]\)，没有可从Y预测的交叉矩C，不能替代RRR。历史嵌套OOF没有收益，这个负结论只约束冻结线性类，不否定一般非线性表示。

**源码。** [R5 `covariance_fit/prediction`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/rrr_controls.py:15)；[nested RRR](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/nested_rrr.py:51)。

## F-PROB-021｜feature/action蒸馏及有限风险恒等式

R5 teacher四通道为 \(T_q(Y)\)，student provider \(g_\theta(Y)\)。归一化feature损失是 \(L_F=E\|S_T^{-1}(g_\theta-T_q)\|^2\)。若冻结head为 \(h(Y,T)\)，记student动作 \(a_S=h(Y,g_\theta)\)、teacher动作 \(a_T=h(Y,T_q)\)，则R6使用

\[
L_{GT}=E\|a_S-M\|^2/s_A,\quad
L_{F+A}=L_F+E\|a_S-a_T\|^2/s_A,\quad
L_{GT+A}=(E\|a_S-M\|^2+E\|a_S-a_T\|^2)/s_A,
\]

其中sA是冻结训练目标component variance标度。对任何同一评价分布，记d=aS−aT、e=aT−M，则

\[
R_S-R_T=E\|d\|^2+2E\langle e,d\rangle,
\quad |R_S-R_T|\le D^2+2\sqrt{R_T}D,\quad D^2=E\|d\|^2.
\]

**推导。** 展开 \(\|e+d\|^2\)，再Cauchy–Schwarz。若head对T具有已证明统一Lipschitz常数L，可再用 \(D\le L\|g-T_q\|_{L^2}\)；本项目没有提供紧全局深网L。teacher不是Bayes真值，故交叉项通常不为0；有限VAL展开只是empirical identity，不能升级为真实posterior风险分解。

**源码。** [R5 width128 feature蒸馏](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/distill_compact_w128.py:19)、[loss](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/distill_compact_w128.py:38)；[R6 action losses](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/track_b/run.py:47)；[R6 empirical expansion](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/track_b/report.py:24)。

## F-PROB-022｜整个cell的响应范数界转likelihood界

若对整个cell C有 \(\|f(x)-f(a)\|\le e_C\)，令 \(r_a=f(a)-y_w\)，三角不等式给

\[
\max(0,\|r_a\|-e_C)\le\|f(x)-y_w\|\le\|r_a\|+e_C,
\]

\[
\ell_C=e^{-(\|r_a\|+e_C)^2/2}\le L(x)\le
u_C=e^{-\max(0,\|r_a\|-e_C)^2/2}.
\]

**关键条件。** eC必须覆盖整个cell、全部测量和正确白化；锚点残差、少量节点误差或中心条件数不能替代它。若 \(\mathcal P\) 用公共shift P*缩放，ℓ、u、Z全部乘 \(e^{P_*}\)，任务比值不变；跨记录恢复evidence时 \(\log Z=\log Z_{\rm scaled}-P_*\)。

**源码/来源。** [R8P cell residual bounds](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/04_maxwell.tex:5)；[R9 `global_audit`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:86)、[1D `cell`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:54)。

## F-PROB-023｜evidence分母、遗漏质量与数值代理的区别

若不重叠cells覆盖全部原prior，\(p_j=\pi_0(C_j)\)，且 \(\ell_j\le L\le u_j\)，则

\[
Z_L=\sum p_j\ell_j\le Z\le\sum p_ju_j=Z_U.
\]

遗漏域U的evidence上界 \(B_U\ge\int_U Ld\pi_0\)，已处理域D的正下界 \(Z_{D,L}>0\)，给

\[
\pi(U\mid y)\le\frac{B_U}{Z_{D,L}+B_U}\le B_U/Z_{D,L}.
\]

只有prior尾质量δ而L≤1时只能给BU≤δ；高SNR下Z极小，δ小未必表示posterior tail小。若U为所有被排除cell的并，整体prior质量≤1，且每个cell满足 \(\mathcal P\ge P_*+\gamma\)，则 \(B_{U,\rm scaled}\le e^{-\gamma}\)；这不需要给每个cell各分配质量1再相加。

另一量 \(Z_L/Z_U\) 是保证纳入的子概率质量，\(1-Z_L/Z_U\) 是未分配预算，不是已发现新模式的质量。R9 denominator为 \(\min Z_{\rm numerical}-2(\max Z_{\rm numerical}-\min Z_{\rm numerical})\) 的有限细化代理，**不是证明的Z下界**；因此相应遗漏posterior值只标numerical proxy，strict字段null。

**来源/源码。** [R8P evidence与遗漏界](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/03_certificates.tex:4)；[R9 whole-support audit](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:86)；[v3 qualification postprocess](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_qualification_v3.py:10)。

## F-PROB-024｜正负任务矩与归一化中心化包络

对有界任务τ，定义原prior正负矩

\[
A_j^+(c)=\int_{C_j}(\tau-c)_+d\pi_0,\quad
A_j^-(c)=\int_{C_j}(c-\tau)_+d\pi_0.
\]

由ℓ≤L≤u分别控制正负部分，得

\[
H_L(c)=\sum_j[\ell_jA_j^+(c)-u_jA_j^-(c)]
\le\int(\tau-c)Ld\pi_0
\le H_U(c)=\sum_j[u_jA_j^+(c)-\ell_jA_j^-(c)].
\]

若ZL>0，\(H_L,H_U\)各自零点cL、cU满足 \(c_L\le E_\pi\tau\le c_U\)。

**推导与sharpness。** 每个可行likelihood的中心化积分是 \(Z_L^{\rm actual}(E\tau-c)\)，对c严格递减，斜率≤−ZL。为了最小化积分，在τ>c处取ℓ、τ<c处取u；最大化相反，故上下零点是给定包络信息类的可达端点。该“最紧”并不保证极值likelihood能由Maxwell产生。

若A矩也有区间，安全下函数取 \(\ell A_{+,L}-uA_{-,U}\)，上函数取 \(uA_{+,U}-\ell A_{-,L}\)。只在符号被严格确定时排除根区间，不能对未控制误差的moment点估计声称严格bisection。

**源码/来源。** [R8P theorem](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/03_certificates.tex:25)；[R9 `probability_envelopes.run`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_envelopes.py:13)的prior矩是GL数值点估计；[1D `H`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:67)使用quintic/Gaussian矩，仍无定向舍入资格。

## F-PROB-025｜硬事件保持共享分母的sharp区间

令 \(\alpha_j=\pi_0(C_j\cap E)\)，\(p_j-\alpha_j\) 是同cell补事件质量。则

\[
P_L=\frac{\sum\ell_j\alpha_j}{\sum\ell_j\alpha_j+\sum u_j(p_j-\alpha_j)},\quad
P_U=\frac{\sum u_j\alpha_j}{\sum u_j\alpha_j+\sum\ell_j(p_j-\alpha_j)}.
\]

**推导。** 概率为事件evidence/(事件+补事件evidence)，分别让事件端取下/上、补事件取上/下即得。若α只有区间且pj精确，下端取αL、上端取αU。区间分子与分母互相包含，不应把N与Z当独立误差任意相减。任务在cell内恒为0/1，只省了任务边界分辨率，层间evidence仍决定概率。

**源码。** [R8P event公式](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/03_certificates.tex:67)；[R9 1D EL/EU/CL/CU](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:76)；[二维numerical envelope](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_envelopes.py:13)。

## F-PROB-026｜物理近似误差与任务中心化，不能漏掉求积项

令 \(\widetilde L\) 为近似物理likelihood，\(\hat N,\hat Z\) 为其数值积分，\(\hat\mu=\hat N/\hat Z\)。在 \(\hat Z>0,Z\ge Z_L>0\) 下，有恒等式

\[
\mu-\hat\mu=\frac{\int(\tau-\hat\mu)L\,d\pi_0}{Z}.
\]

若中心化求积误差满足

\[
\left|\int_D(\tau-\hat\mu)\widetilde Ld\pi_0-(\hat N-\hat\mu\hat Z)\right|\le\epsilon_{Q,c},
\]

则由括号内最后一项为0得到

\[
|\mu-\hat\mu|\le\frac{\epsilon_{Q,c}+\int_D|\tau-\hat\mu||L-\widetilde L|d\pi_0+\int_{D^c}|\tau-\hat\mu|Ld\pi_0}{Z_L}.
\]

这分别包含求积、物理、遗漏域三项。若每cell \(|\log L-\log\widetilde L|\le\delta_j\)，则 \(|L-\widetilde L|\le\widetilde L(e^{\delta_j}-1)\)。若只给统一δ且积分精确，bounded likelihood tilt还给 \(\mathrm{TV}(\pi,\widetilde\pi)\le\tanh(\delta/2)\)；但单点forward误差不能提供统一δ。

**来源。** [R8P对旧task-centered物理误差的审计](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/01_audit.tex:39)、[cell/task envelope实际统一处理](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/03_certificates.tex:25)。R9不填strict总误差证书，原因正是这些项未全部外包。

## F-PROB-027｜任务停止、事件动作与新观测价值

若可信任务区间 \([\mu_L,\mu_U]\) 半宽≤ε，中点均值误差≤ε，平方loss excess risk≤ε²。事件误报/漏报代价为cFP、cFN时，动作1阈值为 \(p\ge c_{FP}/(c_{FP}+c_{FN})\)；整个概率区间在阈值同侧即可停止该动作，但概率值可能仍不精确。

对一般有限动作集合，有可信风险区间时，若 \(R_U(a)-\min_bR_L(b)\le\epsilon\)，该动作regret≤ε。真正新增数据Yd才有

\[
\operatorname{Var}(T\mid y)-E[\operatorname{Var}(T\mid y,Y_d)\mid y]
=\operatorname{Var}(E[T\mid y,Y_d]\mid y)\ge0.
\]

**推导。** 前两式分别来自F-PROB-004和线性事件风险；最后式为条件total variance。重新打包同8源并不新增观测。计算完成只缩小数值动作误差，不会降低固定y的真实posterior variance。

**来源。** [R8P stopping/action](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/06_algorithm_cost.tex:29)、[new data value](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/06_algorithm_cost.tex:79)；[1D numerical停止轨迹](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:79)。

## F-PROB-028｜固定ROI均匀复幅度的精确传递

设field方程 \(A_hE_h=E^{inc}\)，\(A_h=I-KD_h\)，更新 \(D_{\chi(v)}=D_h+vPP^T\)，P选256个ROI。预计算

\[
B=A_h^{-1}KP,\quad C=P^TB,\quad E_R=P^TE_h,\quad
O=S(P+D_hB),\quad F_h=SD_hE_h.
\]

场差满足 \(E-E_h=BvP^TE\)，故 \(P^TE=(I-vC)^{-1}E_R\)，而current差为 \(D_h(E-E_h)+vPP^TE\)，得到

\[
F(v)=F_h+O\,v(I-vC)^{-1}E_R.
\]

它保留全部离散多次散射，不是Taylor/Born截断。必须包含原非零host \(.1+.04i\)；shape固定但latent活动维数小不等于材料更新rank1。改变ROI形状或非均匀幅度时，标量v与C的Schur复用不再直接适用。

**源码。** [R9 `Physics.setup`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_transfer.py:116)、[ `Physics.reduced/direct`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_transfer.py:143)；[R8P exact amplitude theorem](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/04_maxwell.tex:89)。完整物理推导另见02章，此处只固定概率节点使用的接口。

## F-PROB-029｜Schur全cell majorant、导数与每概率节点费用

令 \(C=QTQ^*\)，\(O'=OQ,E'=Q^*E_R\)。则 \(F(v)=F_h+O'v(I-vT)^{-1}E'\)，

\[
F'(v)=O'(I-vT)^{-2}E',\qquad
F(v)-F(a)=(v-a)O(I-vC)^{-1}(I-aC)^{-1}E_R.
\]

对中心a、半径r的圆盘，\(d_i=|1-aT_{ii}|-r|T_{ii}|>0\)，vmax=|a|+r。若b是右端绝对值，逆序递推

\[
x_i=\frac{|b_i|+v_{\max}\sum_{j>i}|T_{ij}|x_j}{d_i}
\]

由逐行三角不等式包住整个圆盘的 \(|(I-vT)^{-1}b|\)。对|E′|应用两次，再乘|O′|并取Frobenius norm，给 \(e_C=\sqrt2\,r\||O'|x\|_F/\sigma_c\)。d≤0需细分/保留未知域；非正规性可能使界很松。普通Schur残差、Q非酉、基列误差与系数rounding均须另传到输出，才能升级strict。

费用按数学对象分开：setup=1个576 LU+8 host源RHS+256基列RHS+1个256 Schur；每v节点=8个256三角RHS+读出；majorant还需两次256逆序绝对值递推。RHS意义不同，reduced工作不能记“0物理费用”。

**源码。** [R9 `majorant_derivative`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_transfer.py:159)、[ `reduced`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_transfer.py:143)、[ `Book.reserve/span`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_transfer.py:85)；[R8P majorant/cost](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/04_maxwell.tex:123)。

## F-PROB-030｜原P1一维正交条件模型：Gaussian尾部变平台原子

本轮追加条件模型固定 \(x_{1:63}=0,x_{64:65}=0\)，native属性均为.5，只保留原 \(x_0=z\sim N(0,1)\) 未截断。于是

\[
D=1.3+.105i,\quad v(z)=D H(1.3z),\quad A=1/1.3,
\]

\[
\mu_v=\Phi_N(-A)\delta_0+\Phi_N(-A)\delta_D+
(v)_\#(1_{-A<z<A}\varphi(z)dz).
\]

任意任务的evidence分子（τ=1给分母）为

\[
N_\tau=\Phi_N(-A)\{\tau(0)L(0)+\tau(D)L(D)\}
+\int_{-A}^{A}\tau(v(z))L(v(z))\varphi(z)dz.
\]

令t=ΦN(z)，有限积分亦等于 \(\int_{\Phi_N(-A)}^{\Phi_N(A)}\tau(v(\Phi_N^{-1}t))L(v(\Phi_N^{-1}t))dt\)，构成不同坐标的普通数值交叉。均值任务为.1+1.3H(1.3z)，ROI硬事件在此均匀slice中恰为z≥0；它是很容易的阈值配置，不能代表一般空间事件。固定观测为已暴露R6 REGULAR，不是原Full66边缘。

**源码。** [R9 1D contract与integrals](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:6)、[CDF coordinate](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:23)、[exact polynomial moments](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:48)。Gaussian无限tail质量被显式纳入，未以[-4,4]改写原prior。

## F-PROB-031｜任务–响应covariance、观测敏感性与源层分解

用 \(r=\sqrt2\operatorname{pack}F\)、\(y_r=\sqrt2\operatorname{pack}y\)，原势能为 \(\|r-y_r\|^2/(2\sigma_c^2)\)。允许微分穿过积分时，商法则给

\[
\nabla_{y_r}E[\tau\mid y]=\sigma_c^{-2}\operatorname{Cov}_{\pi_y}(\tau,r).
\]

对白化观测y_w，等价 \(\nabla_{y_w}E\tau=\operatorname{Cov}(\tau,f)\)。这说明“弱响应变化”与任务相关性也有关；一个单点奇异值并不直接给posterior任务误差。

源层标签J可用total covariance分解

\[
\operatorname{Cov}(\tau,r)=E[\operatorname{Cov}(\tau,r\mid J)]
+\operatorname{Cov}(E[\tau\mid J],E[r\mid J]).
\]

R6因此分别保存within、between、total任务响应covariance。pure/transition材料像有重叠，所以between-source项不能改名between-physical-mode项。

**源码。** [R6 `summarize`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/track_a/stratified.py:125)，其中fw=√2packF，再除σc²得到上述未白化r观测的敏感性口径。

## F-PROB-032｜数值资格、严格证书与公平费用分解

严格任务证书需要F-PROB-023的**真分母下界**，加完整prior/strata/尾域、整个cell的响应/likelihood界、可靠prior任务矩及算术外包。finite refinement差、QUADPACK error estimate、浮点Schur/majorant和正数denominator proxy均不自动满足这些条件。

R6旧门为均值max absolute≤1e-5 **且** ≤.01参考SD；源层/冻结事件绝对≤.002。reference代理使用再紧五倍的2e-6、.002SD、.0004。本轮数值达门不改变strict字段null。近零SD必须同时报告尺度退化。

费用定义

\[
C_{total}=C_{setup}+C_{online\ physics}+C_{prior/tasks}+C_{profile}+C_{enclosure}+C_{reference/qualification}+C_{rejected/replay}.
\]

二因子分别是“相同概率规则下direct/同一reduced物理”和“相同物理下不同概率规则”。同GH节点、prior与任务的direct/reduced对照只检验物理复用；不能归于新概率理论。setup+积分不等于完整cold start，共同reference/验证也不能不记。已存在结果只生成费用CSV，不重跑任何计时。

**源码/证据。** [R9 canonical v3](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/results/probability/final_numerical_comparison_v3.json)、[aggregate costs](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/results/probability/aggregate_costs.json)、[ledger](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/audit/probability/ledger.jsonl)。详见本页配套 `04_PROBABILITY_COST_TWO_FACTOR.csv`；未执行格保持 `NOT_RUN`，无完整cold wall处留空。

## F-PROB-033｜GN适应Gaussian只作求积坐标：GH importance校正

R9先用真实forward残差定位 \(\hat v\in\mathbb R^2\)，取GN矩阵 \(H_{GN}=J^TJ\succ0\)，\(C=H_{GN}^{-1}=LL^T\)。定义辅助密度 \(q(v)=N(\hat v,C)\)，则原材料积分恒等地改写为

\[
Z_{\rm tr+pure}=E_{q}\!\left[\frac{L_y(V)\{\rho_{\rm tr}(V)+\rho_{\rm pure}(V)\}}{q(V)}\right].
\]

若 \((g_i,w_i)\) 为 \(\int e^{-g^2}h(g)dg\) 的Gauss–Hermite规则，节点 \(v_{ij}=\hat v+\sqrt2L(g_i,g_j)^T\)，缩放权重为

\[
a_{ij}=\frac{w_iw_j}{\pi}\,e^{P_*-\mathcal P(v_{ij})}
\frac{\rho_{\rm tr}(v_{ij})+\rho_{\rm pure}(v_{ij})}{q(v_{ij})},
\quad \log q=-\log(2\pi)-\tfrac12\log\det C-\tfrac12\|L^{-1}(v-\hat v)\|^2.
\]

材料mean/variance、source概率与事件取同一组归一化 \(a_{ij}/\sum a\)，host原子另加。**推导。** 先乘除q，再标准化Gaussian；GH的 \(\sqrt2\) 与 \(1/\pi\) 来自二维标准Gaussian。这里GN不是精确势能Hessian也无妨，因为节点上仍用真实L/q修正；但有限阶GH仍有求积误差与域外遗漏，不能由坐标良好推断reference合格。直接把q或profile当posterior会删除真实likelihood非线性与prior fibre。

**源码。** [R9 `locate`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:18)、[ `gh`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:27)；原R6 [ `local_material_reference.main`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/track_a/local_material_reference.py:22)。

## F-PROB-034｜已知交换对称quotient的分支条件概率

既有Z-ALG二单元例有材料 \((x_1,x_2)\in[.5,4]^2\)、非对称prior \(\rho(x_1,x_2)=x_1/C\)，forward在交换下不变。对有序根 \(a>b\)，两个原材料分支是(a,b)、(b,a)，故

\[
\rho_Q(a,b)=\rho(a,b)+\rho(b,a)=(a+b)/C,
\quad \Pr(x_1=a\mid\{a,b\},y)=a/(a+b).
\]

**推导。** 每个交换分支Jacobian绝对值为1，两个分支L相同，所以Bayes条件概率中的L抵消，只剩原prior比。于是

\[
E[x_1\mid\{a,b\},y]=\frac{a^2+b^2}{a+b},\qquad
\Pr(x_1>x_2\mid\{a,b\},y)=a/(a+b).
\]

若保留正确总densityρQ却错误给两分支各1/2，evidence完全相同，事件却恒为1/2。固定无序pair{1,3}上的分支权重为1/4、3/4；有限噪声时还须对整个有序quotient积分，其事件概率不是固定pair的3/4。这个反例说明evidence正确不能单独资格材料任务。

**源码。** [既有 `integrate_ordered_quotient`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/inputs/zalg/Research_Z_ALG_Research_Kit/zalgo_probability/e4_two_cell_posterior.py:148)；固定pair与有限噪声结果引用既有记录，没有本次重算。此对称例不提供一般非对称Maxwell图册。

## F-PROB-035｜有限求积权重、population variance与K32的ddof0

对给定非负有限求积 \(a_n\)，\(p_n=a_n/\sum a_n\)，将材料向量实化为\(m_n\)，数值输出为

\[
\hat\mu=\sum_np_nm_n,\quad
\hat\Sigma=\sum_np_n(m_n-\hat\mu)(m_n-\hat\mu)^T,
\quad \hat P(E)=\sum_np_n1_E(m_n).
\]

这描述有限积分近似的矩，不是一个连续posterior的TV证书：若真实材料posterior在该层非原子，有限节点测度与该连续测度互相奇异，在该层通常TV=1。合适的资格是连续积分误差、任务误差/包络或另外构造的连续近似，而不是把节点权重相近称为posterior TV小。

另对q的一个实材料分量M的K个iid样本，历史K32使用

\[
\hat v_{K,0}=K^{-1}\sum_{k=1}^K(M_k-\bar M)^2,
\quad E\hat v_{K,0}=\frac{K-1}{K}\operatorname{Var}_q(M).
\]

**推导。** 展开 \(\sum(M_k-\bar M)^2=\sum(M_k-\mu)^2-K(\bar M-\mu)^2\)，取期望。确定性decoder积分给population q variance，没有凭空乘31/32；有限K统计差异必须随接口口径保留。二者仍都只是近似q的矩，不认证真实posterior。

**源码。** [R9 `summary`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:12)；[R6 `summarize`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/track_a/stratified.py:125)；[R5 Gaussian矩与ddof边界](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/THEORY_AND_NUMERICAL_CHECKS.md:21)；[R5 `material_moments.py` population定义](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:1)。

## 最小局限账

1. **响应残差界与遗漏prior质量是不同输入。** F-PROB-022给每cell likelihood界；还需原prior质量/任务矩及F-PROB-023的正分母。大量小残差候选或宽域“没有找到新峰”均不能替代这些量。
2. R9是普通浮点全支持审查与数值交叉，不是strict certificate。1D upper likelihood floor1e-300是保留下溢不确定性的数值处理，不能认证Hankel/Schur系数或极小概率的有效位数。
3. 原R6先验截断模型与原1D未截断条件模型分别有自己的normalizer/atoms；二者都不提供Full66 posterior。Gaussian decoder矩、RRR、蒸馏风险收益也不提供该资格。
4. 全部historical负结果、first-run下溢报告、qualification修复、未运行direct矩阵格保留。本文新增完整RHS、reduced likelihood调用、LU、benchmark调用均为0。

## 配套费用CSV的读法

`04_PROBABILITY_COST_TWO_FACTOR.csv`只读取既有canonical v3、1D v2与receipt，未调用physics模块。56行区分：共同一次setup、在线同节点direct/reduced积分、独立reference、两轮全支持审查、数值任务包络、下溢修复重放与未执行格。`aggregation_role=COMPONENT`是组件；`NONADDITIVE_VIEW/CUMULATIVE_TRAJECTORY/JOB_SUBTOTAL/CAMPAIGN_TOTAL`等带NONADDITIVE的行与组件重叠，禁止直接把全CSV相加。

`denominator_rigorous=FALSE`贯穿数值comparison；分母列写清有限细化代理或普通quad值，不能读取为严格Z下界。`full_cold_start_wall=NOT_RECORDED`保留完整冷启动未测的事实；既有setup+integral组件只在自己的列出现。direct source积分与direct task envelope均明确`NOT_RUN`。总账仍是3624完整complex RHS、421个完整LU、1个Schur、1452408个reduced RHS；这次公式与CSV整理新增费用全部为0个forward/likelihood/benchmark调用。
