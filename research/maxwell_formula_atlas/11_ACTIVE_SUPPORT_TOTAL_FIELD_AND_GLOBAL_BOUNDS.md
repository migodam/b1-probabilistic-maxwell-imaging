# 11 · 精确活动支撑、总场投影与全形状损耗界

R11-S3，2026-10-30。接续 F-BI-019—024。以下为离散 Maxwell 的解析推导、标准投影与条件矩阵界；不以机器精度重放证明其新颖性。计算收益与数字图像资格见 [R11 报告](../p3l_r11_stable_structure_digits/RESEARCH_REPORT_ZH.md)。R11 的全阶 LU 自设上限失守，实验总状态为 `FAIL_CAP_EXCEEDED`；这不改写以下恒等式，也不允许把数值报告升级为预算合规 PASS。

统一尺寸：$n=256$，$n_r=32$，$n_s=8$。$C\in\mathbb C^{n\times n}$、$O\in\mathbb C^{n_r\times n}$、$E_R\in\mathbb C^{n\times n_s}$、$F_h\in\mathbb C^{n_r\times n_s}$ 固定；$a\in[0,1]^n$、$A=\operatorname{diag}(a)$、$D=\delta A$ 改变。$\delta$ 在原被动物性矩形内，$\sigma$ 和 source-major 白化沿用 F-BASE-002。$^*$ 表示共轭转置；解析导数转置求解与 Hermitian 内积不可混同。

## F-SS-001 · 严格零支撑的精确电流消元

设 $\mathcal S=\{i:a_i>0\}$，$\mathcal I$ 为其补集。由 $h=D(E_R+Ch)$ 逐行得 $h_\mathcal I=0$，故

$$M_S=I_S-\delta A_SC_{SS},\quad h_S=M_S^{-1}\delta A_SE_{R,S},\quad F=F_h+O_Sh_S.$$

按 $S,I$ 排序，$I-DC$ 是 $\begin{bmatrix}M_S&-D_SC_{SI}\\0&I_I\end{bmatrix}$，因此它与 $M_S$ 的可逆性等价。只删**精确零**；小正数删减是另一项近似。全部源共享一次分解。空支持时 $F=F_h$。资格：标准精确块消元；实现 `agent_a/solver.py:Solver.evaluate`。

## F-SS-002 · 完整端口灵敏度的活动/非活动分块

定义 $T=O(I-DC)^{-1}$ 和 $Z=O_SM_S^{-1}$。对上述块三角矩阵求逆直接得

$$T_{:,S}=Z,\qquad T_{:,I}=O_I+ZD_SC_{SI},\qquad e=E_R+C_{:,S}h_S.$$

只需活动系统的 $n_s$ 个正演和 $n_r$ 个转置右端；不需要另做全 $n$ 阶 LU。即使材料有零格点，仍要保留全部 $n$ 个材料导数列。资格：标准逆矩阵分块；不能把活动电流状态数当材料参数数。

## F-SS-003 · 支撑出生的非零材料导数

由 $(I-DC)\,dh=(dD)e$，得到

$$dF=T(dD)e,\qquad \frac{\partial F_{rs}}{\partial a_i}=\delta T_{ri}e_{is},\qquad \frac{\partial F}{\partial\delta}=T A e.$$

当 $a_i=0$ 时，$h_i=0$，但 $\partial F/\partial a_i$ 一般非零。因此原生占比优化可以用这些列激活新支持。对原 P1 latent 则须再乘 $Da/Dz$；平台内该因子可为零，不能由本式推出 latent 梯度必能出生新笔画。实虚物性导数分别为 $F_\delta$ 与 $iF_\delta$，Gaussianized 属性还须乘正态 CDF 导数。资格：精确一阶导数，不是新测量信息。

## F-SS-004 · 不随支持改变的总场注入

从 $e=E_R+Ch$、$h=\delta Ae$ 得

$$M_e e=E_R,\quad M_e=I-\delta CA,\quad F=F_h+\delta OA e.$$

$E_R$ 始终固定，避免旧加权状态 $u=\delta\sqrt A e$ 的注入 $\sqrt A E_R$ 随形状改变。这里不包含对低维共同不变空间的任何保证。资格：精确状态变换；零 $a_i$ 无需除法。

## F-SS-005 · 固定总场 Galerkin 与真实缺陷

取固定正交基 $V\in\mathbb C^{n\times r}$，例如 $\operatorname{orth}[E_R,CA_0E_R,\ldots]$ 或 $\operatorname{orth}[E_R,CE_R,\ldots]$。每个当前材料解

$$(I_r-\delta V^*CAV)t=V^*E_R,\quad\widetilde e=Vt,\quad\widetilde F=F_h+\delta OAVt,$$

$$r_e=E_R-(I-\delta CA)\widetilde e.$$

固定基不等于固定 $V^*CAV$。每次形状变化仍需更新 $CAV$、$OAV$ 和投影；其成本必须计入。非正规 Galerkin 的稳定性不由基正交自动保证。资格：标准近似；实现 `agent_b/fixed_total_field.py`。

## F-SS-006 · 全部占比形状的 Hermitian 虚部下界

令 $H_C=(C-C^*)/(2i)$、$P=\sqrt A$、$B=PCP$，有 $\operatorname{Im}_H B=PH_CP$。记 $\lambda=\lambda_{\min}(H_C)$，$\beta_*=\min(0,\lambda)$。对任意 $v$，

$$v^*PH_CP v\ge\lambda\|Pv\|^2\ge\beta_*\|v\|^2,$$

最后一步在 $\lambda<0$ 时用 $\|Pv\|\le\|v\|$，在 $\lambda\ge0$ 时直接用非负性。因此 $\lambda_{\min}(\operatorname{Im}_H B)\ge\beta_*$。资格：条件标准矩阵序不等式；不能把负 $\lambda$ 随意抹为零。

## F-SS-007 · 被动物性矩形上的统一强制性常数

对 $\delta=x+iy$，本合同 $0<y_{min}\le y\le y_{max}<x_{min}\le x\le x_{max}$。$y/(x^2+y^2)$ 随 $x$ 递减、随 $y$ 递增，故

$$b_*={y_{min}\over x_{max}^2+y_{min}^2}+\beta_*.$$

当 $b_*>0$，令 $K_\delta=\delta^{-1}I-B$，则 $-\operatorname{Im}_H K_\delta\succeq b_*I$。Cauchy–Schwarz 给

$$b_*\|v\|^2\le|v^*K_\delta v|\le\|v\|\,\|K_\delta v\|\Longrightarrow\|K_\delta^{-1}\|_2\le b_*^{-1}.$$

若 $b_*\le0$，本充分界失效，不能改用经验正数。本合同有限浮点结果 $\lambda=-4.5739299047\times10^{-5}$，$b_*=0.009210948658>0$；不是向外舍入的连续 PDE 证书。

## F-SS-008 · 总场缺陷到端口误差的精确矩阵顺序

由 push-through 恒等式，$P K_\delta^{-1}P=\delta A(I-\delta CA)^{-1}$。因此

$$F-\widetilde F=OP K_\delta^{-1}P r_e.$$

这个式子在 $A$ 奇异时仍成立，不需要 $P^{-1}$。直接展开 Neumann 多项式并由可逆区域延拓，或用 $(I-AB)^{-1}A=A(I-BA)^{-1}$，均可得到。资格：离散精确恒等式，不需要 $C,A$ 对易。

## F-SS-009 · 全形状白化端口误差界

将前两式取范数并用 source-major 等距实化，有

$$\|W(F-\widetilde F)\|_2\le\epsilon_W={\sqrt2\over\sigma b_*}\|OP\|_2\|P r_e\|_F,$$

$$\|OP\|_2\le\|O\|_2.$$

当前 $\|OP\|_2$ 可用 $32\times32$ 行 Gram 的最大本征值求得；全统一 $\|O\|_2$ 更便宜但更保守。资格：在上述耗散假设下的矩阵范数界；实现是 float64 数值条件界。

## F-SS-010 · likelihood 势误差与回退门

$\Phi(F;y)=\frac12\|W(F-y)\|^2$。若 $\|W(F-\widetilde F)\|\le\epsilon_W$，则

$$|\Phi(F;y)-\Phi(\widetilde F;y)|\le\|W(\widetilde F-y)\|\epsilon_W+\tfrac12\epsilon_W^2=\epsilon_\Phi.$$

接受门为 $\epsilon_W\le0.01$ 且 $\epsilon_\Phi\le0.1$；不通过有限增阶后精确回退。实际配对误差只能做事后评估，不能反过来替代在线界。这是误差控制，不是 posterior 认证。

## F-SS-011 · 近似物理导数不自动继承响应误差门

对精确 $T$、$e$，$J_i=\delta T_{:,i}e_{i,:}$。若用近似 $\widetilde T,\widetilde e$，则

$$J_i-\widetilde J_i=\delta[(T-\widetilde T)_{:,i}e_{i,:}+\widetilde T_{:,i}(e-\widetilde e)_{i,:}].$$

只有响应残差门并不控制全部 Jacobian。R11 的 B/C 在需要 Full-J 时使用精确一致的 F/J；只对 forward trial 尝试降阶。因此不能宣称降阶同时省了32个导数右端。资格：精确误差分解；受控降阶导数尚未实现。

## F-SS-012 · 三类成本与图像收益不可互换

对于相同观测、参数化、正则、初值和迭代配置，报告

$$C_{system}=C_{load}+C_{basis}+\sum C_{trial}+\sum C_J+C_{gate}+C_{fallback}+C_{validation}+C_{optimizer}+C_{save}.$$

活动支撑减少的是物理状态阶数，不改变可用观测；固定总场与耗散界也不增加材料信息。故“相同材料结果、较低 $C_{system}$”是求解器收益；“同预算更低材料风险/更好拓扑”需另有实验。R11 已有 72 份真实反演，但 B/C 全部 forward trial 回退，图像没有新增改善；支撑的正面核结果不能掩盖全阶 LU 预算失败。资格：评价定义及可审计计算分解，不是物理定理。

### 本轮原生反演目标与IRLS的具体含义

N模型令 $x=(a,\xi)$，$a\in[0,1]^{256}$，$\delta(\xi)=0.8+\Phi_N(\xi_1)+i[0.03+0.15\Phi_N(\xi_2)]$。$L$为ROI相邻水平/竖直差分，冻结目标为

$$\mathcal J_N=\Phi(F;y)+0.4\sum_j[\sqrt{(La)_j^2+0.05^2}-0.05]+0.01\|a\|^2+0.05\mathbf1^Ta+\tfrac12\|\xi\|^2.$$

$w_j=\sqrt{(La)_j^2+0.05^2}$，TV梯度为 $0.4L^T[(La)/w]$，代码采用正定IRLS majorizer $0.4L^T\operatorname{diag}(1/w)L$，而非TV精确Hessian $0.4L^T\operatorname{diag}(0.05^2/w^3)L$。阻尼GN解线性步后对$a$裁到$[0,1]$并用原目标回溯；未读GT。该新原生正则不等于原P1 Gaussian64 prior。

P模型仍为原Gaussianized66目标 $\mathcal J_P=\Phi(F(G(x));y)+\tfrac12\|x\|^2$。这一点估计目标不是posterior均值，也不等于native Uniform2坐标下的点目标。N与P不能当成纯后端对照。

### 新N模型的尺度gauge与材料任务

当 $a$ 为自由占比而 $\delta$ 未知，对任意使两端均合法的 $c>0$，有

$$a'=a/c,\qquad\delta'=c\delta,\qquad\delta'\operatorname{diag}(a')=\delta A.$$

这是同一材料的参数化gauge，不是新的Maxwell材料世界。占比阈值事件 $\{a_i>0.5\}$ 一般不沿该gauge不变；复材料误差或固定 $\chi_i$ 阈值才是材料任务。这一结论仅针对本轮新N自由占比参数化，不能从中推断原Gaussian64 quintic image具有普遍全局重缩放gauge。
