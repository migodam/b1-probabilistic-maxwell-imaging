# B1 × C1：材料歧义的概率定义与可检验数学

实验阶段20261013；公开交付2026-10-14。Fresh主概率参考与范围以 [ACCEPTED_REFERENCE_AUDIT](results/ACCEPTED_REFERENCE_AUDIT.json) 和 [全部fresh细化记录](results/ALL_FRESH_REFERENCE_REFINEMENT.json) 为准；旧13阶与中间结果保留为历史近似。

本研究问的是：散射数据还容许哪些不同材料世界，以及各材料区域应获得多少概率。它不把低残差、生成器 confidence 或连续图像的点密度直接解释为重建概率。下述数学用于冻结的同心分层球校准台；它没有证明一般三维逆散射的唯一性或无限维后验误差界。

## 1. 同一个材料 prior，而非预设结构专家

设连续系数与材料参数为

\[
x=(z_c,z_s,\epsilon',\epsilon'')\in[-3,1]^2\times[2.2,2.8]\times[.04,.16],
\]

四者独立均匀。径向 implicit field 为

\[
u(s)=1+z_c[1-(s/.45)^2]_+^2+z_s[1-((s-.65)/.18)^2]_+^2.
\]

两个 bump 的支持不重叠；`u<0` 的区间取严格真空，其他区间取同一个复介电常数。系数生成并没有接受 topology label；中心、环形或组合空腔是 decoder 的事后属性。

对任一个 bump，几何尺度 \(g\)（中心半径或环形半宽）在 \(z<-1\) 时满足

\[
g=a\sqrt{1-(-z)^{-1/2}},\qquad z=-[1-(g/a)^2]^{-2}.
\]

原 coefficient prior 的 pushforward 是一个混合测度：零尺度的质量为 0.5，正尺度的密度为

\[
\rho_a(g)=\frac{g}{a^2}[1-(g/a)^2]^{-3},\quad 0<g\le a\sqrt{1-1/\sqrt3}.
\]

导数给 \(|dz/dg|=4g a^{-2}[1-(g/a)^2]^{-3}\)，乘原密度 \(1/4\) 即得上式。连续部分积分为 0.5，不能再对它归一化后丢掉零尺度原子。两个 inactive 区间各质量 0.5，所以实心材料 prior mass 为 0.25。改成 radius-uniform 是修改 prior，不能称为坐标变换。

inactive coefficient 的条件 prior 仍是 `[-1,1]` 上均匀。不同 inactive 系数可以产生完全相同的材料，因此 latent decoder 不是一一映射。网络在 latent 空间的 KL 不是 physical material KL；对确定 decoder，数据处理不等式只给

\[
\mathrm{KL}(D_\#P\Vert D_\#Q)\le\mathrm{KL}(P\Vert Q).
\]

本实现没有直接估计左侧，也不把右侧冒称材料场精度。

## 2. 噪声与 proper nuisance 的完整 target

每频率使用 16 个 complex channels，先 Re16、后 Im16，构成 32 个实通道。每个实通道的标准差为冻结 \(\sigma\)；即 complex noise 的 real/imag 分量各有方差 \(\sigma^2\)。数据模型为

\[
y=gF(x)+e,\qquad e\sim N(0,\Sigma),\qquad g\sim N(1,.15^2)\mid .7\le g\le1.3.
\]

这里 gain 是一个**实数幅度** nuisance；没有额外的未知复相位。两频共享同一 gain、材料与复 ε。频带内 ε 不随频率变化是显式的受限模型假设。

先考虑未截断 gain。它的 marginal covariance 为

\[
C=\Sigma+s_g^2FF^T,
\]

且 \(y\sim N(F,C)\)。令 \(a=\Sigma^{-1/2}F\)、\(r=\Sigma^{-1/2}(y-F)\)，则

\[
\log L_\infty=-\tfrac12\left[r^Tr-\frac{s_g^2(a^Tr)^2}{1+s_g^2a^Ta}+\log|\Sigma|+\log(1+s_g^2a^Ta)+n\log(2\pi)\right].
\]

未截断的条件 gain 为 Gaussian，方差 \(v=(F^T\Sigma^{-1}F+s_g^{-2})^{-1}\)，均值 \(m=v(F^T\Sigma^{-1}y+s_g^{-2})\)。合法的 bounded marginal 为

\[
L_T=L_\infty\frac{\Phi((1.3-m)/\sqrt v)-\Phi((.7-m)/\sqrt v)}{\Phi(2)-\Phi(-2)}.
\]

logdet 和截断归一化都不能省略。对两频通道先拼接，再边缘化一次；分别积分会引入两个 gain，定义另一 target。实现用稳定的 tail log-CDF 差，零物理审计以独立一维积分和 dense Cholesky 交叉核对。

冻结 prior 与上述 likelihood 后，

\[
p(x\mid y)=\frac{L_T(y\mid x)p_0(x)}{Z(y)},\qquad
P(A_j\mid y)=\frac{\int_{D^{-1}(A_j)}L_Tp_0\,dx}{Z(y)}.
\]

十二个事件是中心半径四个 bins 与环形半宽三个 bins 的笛卡尔积。它们互不重叠并覆盖声明支持，但并不是十二个预设生成专家，也不保证每个 bin 是分离 mode。

## 3. 补偿诊断与 evidence 必须分开

固定一个参考世界，比较另一个不同几何世界的 raw whitened secant，以及 ε、ε+gain 可允许补偿后的 secant。对小变化，若 nuisance/material compensation Jacobian 为 \(J_\nu\)，无约束局部 least-squares 的残差为

\[
(I-P_{\Sigma^{-1/2}J_\nu})\Sigma^{-1/2}\Delta F.
\]

该投影解释一阶可吸收方向，但不保证大结构变化、边界 constraints 或 proper-prior 体积。非线性实验因此最终重新调用 full-Mie，对实际补偿 endpoint 验算。优化找到的合法 endpoint 给最小距离的**可行上界**；未证明全局最小。

即使某个 endpoint residual 很低，其 posterior event mass 仍可能很低：积分还取决于可行邻域体积、ε 范围、gain prior 和其他候选。反之，测量不能区分的材料族不应仅因一个 MAP 更好而删除。P2 的 witness 与 posterior mass 分别来自 `COMPENSATION.json` 和独立参考，不能互相替代。

## 4. Full-minus-proxy correction

给定一个连续 proxy response，可先算 \(Z_P(A)\)。从已知、完整支持的事件条件 proposal \(q_A\) 抽样，有

\[
Z_F(A)=Z_P(A)+E_{q_A}\left[(L_F-L_P)\frac{p_0}{q_A}\right].
\]

本实现用 0.95 adaptive-cell proposal 加 0.05 事件条件 proper prior，并在原 latent `dx` 中算准确 proposal density。identity 与控制变量形式是已知 importance 方法的组合，不是新的 Maxwell 定理。估计事件积分再取比值，在有限样本下有 ratio bias；bootstrap 仅量化已访问 residual 的 Monte Carlo 误差。它不包含 cubature 误差、未检测 narrow region 或模型外误差。

网络 correction 使用另一种合法方案：\(q=.9q_\theta+.1p_0\)，full-likelihood SNIS weight 为 \(w\propto L_Fp_0/q\)。正常 density、Jacobian 与 full support 使 target 清楚，但 64 draws 不提供完整 mode discovery 保证。高 ESS 也不能排除遗漏重要区域。

## 5. 条件 density 与 teacher 的信息量

bounded MDN 在 logit 坐标输出 normalized Gaussian mixture，并显式保留变换 Jacobian。posterior-distilled NLL 满足

\[
E_{p(x\mid y)}[-\log q_\theta(x\mid y)]
=\mathrm{KL}(p\Vert q_\theta)+H(p).
\]

teacher target 固定时优化 NLL 对应降低该 latent KL。posterior teacher 与物理 bank 是额外监督，不能把它们当作免费信息；full/diagonal distillation 共用 teacher 的比较才隔离 covariance 表达差别。

## 6. 下游信息的合法传递

若独立 geometry sensor \(c\) 在给定材料后与 scattering measurement 条件独立，则

\[
p(x\mid y,c)\propto p(c\mid x)p(x\mid y).
\]

P2 的声明 sensor 是归一化真空体积加 \(N(0,.01^2)\) 噪声。truth 只生成模拟 sensor，并不直接输入推断。没有新信息的常数 context likelihood 保留原权重。Top-K 代表点的 sensor likelihood 通常不等于该事件内部的条件平均 likelihood，因而代表点更新可能丢掉相关性与 tail；条件样本加 remainder 才能实现有限 particle package 的完整重新加权。

事件mass加代表图一般不是后续Bayesian更新的充分交接。条件期望的严格式、TV/regret误差界与latent/material KL链式分解见[THEORY_HANDOFF.md](THEORY_HANDOFF.md)。这是已知概率关系的推导，不作为原创定理；本实验没有取得其uniform误差界所需的全context常数。

这些结论支持的是正确的概率接口和受限数值 target。它们不提供现实模型 adequacy、无限维 posterior support、一般三维 uniqueness、所有 mode 发现或 IEEE TAP 原创性保证。

代码证据：[material](../../src/b1_joint/material.py)、[likelihood](../../src/b1_joint/likelihood.py)、[reference](../../src/b1_joint/reference.py)、[correction](../../src/b1_joint/correction.py)、[proposal correction](../../src/b1_joint/proposal_correction.py)、[package](../../src/b1_joint/package.py)。文献与访问限制见 [主张账本](LITERATURE_AND_CLAIM_LEDGER.md)。


材料事件context更新的精确期望、TV误差界和latent/material KL链式分解见 [THEORY_HANDOFF.md](THEORY_HANDOFF.md)。这些标准恒等式解释交付信息需求，不被申报为新原创定理。
