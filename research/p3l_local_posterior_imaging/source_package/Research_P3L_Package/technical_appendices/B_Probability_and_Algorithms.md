> 技术附录采用独立章节、公式与文献编号；主报告为最终研究判断与实验配额的权威版本。本附录提供扩展推导，不报告新的 Maxwell 实验。

# P3-L 概率推断、学习与 RL 专题研究稿

撰写日期：2026-10-19。供主报告整合的独立中间稿；本次没有运行新的 Maxwell、训练或付费实验。已实际阅读 `evidence/p1/MODEL_AND_PRIOR.md`、`evidence/audit/RESEARCH_REPORT_ZH.md`、`evidence/p2/RESEARCH_REPORT_ZH.md`、`evidence/p2/THEORY_REFINEMENTS.md`。文献来源列于末章。下文区分标准结论、针对本模型的推导、待证机制。

## 0. 研究判断

最值得先做的不是训练“沿流形给每点一个概率”的网络，而是建立**保留强方向条件不确定性、以材料条件矩为核心统计量的局部 Bayesian benchmark**。弱方向上的一个好概率密度，并不能补回被强方向冻结、局部支持限制或错误 decoder 平均丢掉的材料信息。

建议主攻“材料风险驱动的分块 Bayesian 推断与条件积分”，并行“轻量条件密度的 amortization”和“保持正确目标的计算预算调度”。RL 可以帮助决定如何花 Maxwell 费用；目前没有理由让它负责定义概率。全局模式跳跃不在本专题范围内。

P2 的负结果与本轮兼容：它的训练集仅由 8 个 parent 的 54 个成功探索更新组成，学习的是选择机制产生的更新分布，既没有 GT proper-score 训练，也没有 posterior teacher；所有 posterior probability 仍为 null。因此这些结果反对原来的小样本 learned proposal 配置，不能反对有正确统计目标的 posterior estimator。反过来，P2 的大量有限噪声候选也不是“局部 Bayesian 方法必然有收益”的证据。

## 1. 必须保留的概率模型与坐标

### 1.1 原 prior 不是任意材料体积测度

原始坐标为

\[
q=(z,u),\quad z\in\mathbb R^{64},\quad u\in[0,1]^2,\qquad
p_0(q)=(2\pi)^{-32}e^{-\|z\|^2/2}\mathbf1_{[0,1]^2}(u).
\]

材料参数为 \(A=0.8+u_1\)、\(B=0.03+0.15u_2\)，材料

\[
\chi=g(q)=0.1+0.04i+\phi(z)(A+iB).
\]

已知 ROI 有 256 个 cells；\(\phi\) 是 64 个 Gaussian DCT 系数经 C² 饱和 decoder 得到的 fraction。饱和意味着全局多对一映射，也可能在材料空间出现原子。原后验始终相对于 \(dz\,du\) 定义；不能给原 prior 随手乘上材料 Gram determinant，再声称概率未变。当前 A anchor 的局部 \(Dg\) 数值满秩与全局可能存在多对一、饱和原子并不矛盾。

每个 complex measurement 的噪声为 \(CN(0,\sigma^2)\)，\(\sigma=0.00119211798514519\)，256 complex channels。定义

\[
r(q)=\frac{\sqrt2}{\sigma}\operatorname{stack}_{\mathbb R}(F(g(q))-y)\in\mathbb R^{512},
\quad \Phi(q)=\tfrac12\|r(q)\|^2.
\]

这等于项目采用的 complex 范数 \(\|F-y\|^2/\sigma^2\)。混淆 complex 方差和实部方差会使 posterior 温度差一倍。

### 1.2 全维局部区域和低维曲面必须分开

给定一个只由观测和预先冻结规则定义的局部区域 \(\mathcal C_y\)，令

\[
\pi_C(q\mid y)=Z_C^{-1}p_0(q)e^{-\Phi(q)}\mathbf1_{\mathcal C_y}(q).
\]

\(\mathcal C_y\) 可以是 prior-scaled trust region、固定 chart 覆盖区或给定 basin，但其边界不能由 held-out GT 选择。原 Gaussian likelihood 没有 noise-ball 硬支持；若为了局部研究额外截断，必须称为局部条件后验，并报告截断的证据和未知漏失质量。

设 \(q=\Psi_y(a,b)\) 是 \(\mathcal C_y\) 内的全维可逆 chart，\(a\in\mathbb R^d\)、\(b\in\mathbb R^{66-d}\)。相对于 \(da\,db\)，

\[
\widetilde\pi(a,b\mid y,C)\propto e^{-U(a,b)},
\quad U=\Phi(\Psi)-\log p_0(\Psi)-\log|\det D\Psi|.
\tag{1}
\]

原变量的 property box 变成 chart 中的真实可行域；它往往随 \(a\) 改变。只有仿射 chart 才能忽略一个常数 Jacobian，只有先验正交分块并无额外耦合截断时，才能把先验写成独立的 \(p_0(a)p_0(b)\)。广义材料方向混合 Gaussian 与 uniform 坐标后，不再具有独立 Gaussian prior。

为了第一版代码稳定，建议先只在 64 Gaussian 坐标中取欧氏正交 \(V_w,V_s\)，并把两项 uniform property 单独保留为条件采样块。使用绝对坐标 \(a=V_w^Tz,b=V_s^Tz\)，原 Gaussian 分解仍然独立。若使用 \(z=z_0+V_wa+V_sb\)，对应 prior 的中心是 \(-V_w^Tz_0,-V_s^Tz_0\)，不是零。

可用 \(u_j=\Phi_{\mathrm N}(v_j)\) 将 uniform 精确高斯化，但边界 anchor 会映到无穷，且尾部的坐标压缩容易混淆材料敏感性。P1 A 的损耗位于上界，故不能把此变换未经处理直接用于主方案。原始 pCN + reflection 是更可靠基线。

## 2. 固定、条件优化、边缘化分别计算什么

### 2.1 固定强方向

固定 \(b=b_0\) 后，

\[
\pi_{\mathrm{fix}}(a)\propto e^{-U(a,b_0)}.
\]

它是 chart 中关于 \(b=b_0\) 的 regular conditional distribution，既不是弱边缘后验，也不是完整后验。它能用于受控诊断，但未表示任何强方向随机性。不同 chart 下“固定同一坐标值”可能是不同事件。

### 2.2 优化 ridge

定义 \(b_*(a)=\arg\min_b U(a,b)\)。profile 权重

\[
\pi_{\mathrm{prof}}(a)\propto e^{-U(a,b_*(a))}
\]

描述最优截面上的密度值。它可以被解释为人工坐标 \(c=b-b_*(a)\) 下的 \(c=0\) 条件分布，但这个人工条件不来自 Maxwell 观测，通常不等于 \(\pi_w(a)\)。如果实际 corrector 优化的是纯 residual 或 \(\Phi+\|z\|^2/2\)，而没有 (1) 中的 Jacobian，其 mode 连该 chart 下的 posterior mode 也未必是。

### 2.3 正确弱边缘与 Laplace 修正

\[
L(a)=\int_{B(a)}e^{-U(a,b)}db,
\quad \pi_w(a)=L(a)\big/\int L(\bar a)d\bar a.
\tag{2}
\]

若内部有唯一、非退化条件 mode，\(H_b(a)=U_{bb}(a,b_*)\succ0\)，则局部 Laplace 近似为

\[
\log L(a)\simeq-U(a,b_*)+\frac{s}{2}\log(2\pi)
-\frac12\log\det H_b(a)+\log\tau(a),\quad s=66-d.
\tag{3}
\]

\(\tau(a)\) 是 Gaussian 近似落入实际条件域 \(B(a)\) 的概率。边界远且尾部可忽略时才可取 1。uniform property 接近上界时，半 Gaussian、截断 Gaussian 或低维 property 数值积分通常比忽略边界合理。边界 mode 还可能带非零法向一阶项，不能不加分析套用内部 mode 的半 Gaussian 公式。

\(-\frac12\log\det H_b\) 是条件概率体积：相同 ridge 高度但较宽的强方向条件分布具有更大弱边缘质量。若 \(U(a,b)=a^2/2+b^2/[2s(a)^2]\)，profile 给出 \(e^{-a^2/2}\)，正确边缘却为 \(s(a)e^{-a^2/2}\)。仅在 \(s(a)\) 为常数时二者一致。这个反例不依赖多拓扑或多模式。

Gaussian 联合模型可更明确地比较。若精度矩阵为

\[
H=\begin{pmatrix}H_{aa}&H_{ab}\\H_{ba}&H_{bb}\end{pmatrix},
\]

固定 slice 的精度为 \(H_{aa}\)，弱边缘精度为 Schur complement
\(H_{aa}-H_{ab}H_{bb}^{-1}H_{ba}\)。此时 profile 与边缘仅差常数 \(\det H_{bb}^{-1/2}\)，所以可以有相同弱边缘形状；但强方向不确定性和材料非线性平均仍然丢失。不能用线性 Gaussian 中这项偶然一致性推广到 nonlinear Hessian 随 \(a\) 变化的情况。

### 2.4 C² 正则性限制及可检验 Laplace 误差

P1 decoder 仅保证 C²。因此 \(U\) 一般仅 C²，\(U_b=0\) 的 IFT 只保证 ridge 为 C¹。经典 ridge 曲率与三阶 Laplace 余项需要 C³、分片 C³、或额外的 Hessian Lipschitz/连续模量假设。

不引入隐藏正则性的方法是：在条件 Gaussian 主要质量所在的椭球 \(E_a\) 上验证

\[
\left|U(a,b)-U(a,b_*)-\tfrac12(b-b_*)^TH_b(a)(b-b_*)\right|\le\delta_a.
\]

若还分别控制真实条件质量与 Gaussian 质量在 \(E_a^c\) 的尾部，则有局部积分比的 \(e^{\pm\delta_a}\) 控制和相应尾部加项。没有尾部控制，就只能称为局部诊断。若所有 \(a\) 上得到统一 \(|\log L-\log\widetilde L-c|\le\delta\)，归一化后 density ratio 位于 \([e^{-2\delta},e^{2\delta}]\)，从而能给弱边缘 TV 或 bounded material expectation 误差上界。数值 Hessian 好看并不等于该条件成立。

### 2.5 surface 体积不应与 Laplace 体积混淆

令 \(q_s(a)=\Psi(a,b_*(a))\)，材料曲面诱导 metric

\[
G_s(a)=D(g\circ q_s)(a)^T W D(g\circ q_s)(a).
\]

如果把 (2) 的 \(\pi_w(a)da\) 推到曲面，关于曲面体积 \(d\mathrm{vol}_{G_s}=\sqrt{\det G_s}\,da\) 的密度是

\[
\rho_s(a)=\pi_w(a)/\sqrt{\det G_s(a)}.
\]

这只是同一个概率的参考测度改写。再乘一次 \(\sqrt{\det G_s}\) 会改变概率；它不等于强方向积分中的 \(\det H_b^{-1/2}\)。退化处 \(G_s\) 不可逆或 decoder 多对一时，需分层/多张 chart 与 multiplicity 处理，不能用全局材料 Lebesgue density 替代 latent pushforward。

## 3. 比采样坐标更重要的量：条件材料均值与条件材料协方差

定义

\[
m(a)=\mathbb E[g(\Psi(a,b))\mid a,y,C],
\quad S(a)=\operatorname{Cov}[g(\Psi(a,b))\mid a,y,C].
\]

则严格有

\[
\mu_C=\mathbb E_{\pi_w}m(a),\qquad
\Sigma_C=\mathbb E_{\pi_w}S(a)+\operatorname{Cov}_{\pi_w}m(a).
\tag{4}
\]

**真正有用的弱方向模型应该输出或积分 \(m,S\)，而不只是输出一张强方向修正后的中心材料。** 这是条件期望和全协方差公式的应用，并非新的概率定理，但它改变 P3-L 的算法设计。

即使条件 \(b\) 近似 Gaussian，

\[
m(a)\approx g(q_*)+\tfrac12\sum_{ij}(H_b^{-1})_{ij}
\partial_{b_i b_j}(g\circ\Psi)(a,b_*),
\]

还要加上条件均值与 mode 不同产生的一阶偏移。P1 的饱和 decoder 和均匀参数边界恰好是这些偏移可能显著的机制。应分别测量“弱边缘密度误差”和“强条件材料矩误差”。

对独立 posterior draws，将材料 \(g(q)\) 替为 \(m(a)\) 的 Rao–Blackwell 化能减少 Monte Carlo 方差。对于不同 MCMC 核或任意相关序列，不能无条件声称渐近方差必降；必须实际比较材料统计量的 ESS 和 MCSE。

### 3.1 Bayes 风险恒等式和五项误差账本

对平方材料损失 \(\ell_W(\widehat\chi,\chi)=\|\widehat\chi-\chi\|_W^2\)，设全 posterior 均值为 \(\mu\)。任意估计严格满足

\[
\mathbb E_\pi[\ell_W(\widehat\chi,\chi)]
=\operatorname{tr}(W\Sigma)+\|\widehat\chi-\mu\|_W^2.
\tag{5}
\]

因此，posterior approximation 影响材料风险的量是**材料均值偏差**，不必是整个 latent density 的每一处误差。若局部区域漏失质量为 \(\epsilon=\pi(C^c\mid y)\)，材料集合的 \(W\)-直径为 \(D\)，则

\[
\|\mu_C-\mu\|_W\le\epsilon D,
\quad\text{局部截断导致的额外平方风险}\le\epsilon^2D^2.
\tag{6}
\]

P1 材料有界，故这样的 \(D\) 可由 property bounds、ROI 和任务 norm 保守计算。它不是关于 unbounded latent 的错误有界假设。

若算法用弱分布 \(\nu\)、近似条件矩 \(\widetilde m\) 和有限采样，材料均值总误差可拆为：

\[
\|\widehat\mu-\mu\|_W\le
\underbrace{\epsilon D}_{\text{局部支持}}
+\underbrace{D\,\mathrm{TV}(\nu,\pi_w)}_{\text{弱概率}}
+\underbrace{\mathbb E_\nu\|\widetilde m-m\|_W}_{\text{强条件矩}}
+\underbrace{e_{\rm phys}}_{\text{forward/likelihood偏差}}
+\underbrace{e_{\rm MC}}_{\text{有限采样}}.
\tag{7}
\]

若 \(\nu\) 已包括 physics 近似，则物理误差不能重复计入；实际账本应以一条明确的插入中间目标链定义五项。(7) 与 (5) 给出可实测的失败定位框架，而不是将所有误差归咎 Flow 容量或 weak space 太小。这个 ledger 的概率不等式是标准工具；独立贡献必须来自 Maxwell 对 \(e_{\rm phys}\)、几何刷新与 \(e_{\rm MC}\) 成本的可计算控制。

### 3.2 MMSE、MAP、MAE 不同，输出的物理资格也不同

在线性 forward、线性 material decoder、Gaussian prior 的无边界模型中，收敛 MAP 就是材料 MMSE。此时增加 posterior samples 不能带来真实 Bayes MSE 改进，反而有有限 MC 误差。P3-L 的正面收益必须来自真实非线性、边界/skew、基础估计偏差、任务损失不匹配等机制；“有弱不确定性”本身不是理由。

P1 posterior mean \(\mathbb E[\phi(x)(A+iB)]\) 可能在不同像素具有不同有效实虚比，因而未必属于原来的单一材料相 decoder 像。它保持凸性所保证的被动性，却未必对应某个合法 \(g(q)\)，更不保证 \(F(\mathbb E\chi)=\mathbb E F(\chi)\)。若最终必须输出合法 decoder 图像，正确的受限 Bayes action 是

\[
\widehat q_{\rm action}\in\arg\min_q\|g(q)-\mu\|_W^2,
\]

即把材料均值投影到允许的行动集合。它不是 posterior MAP，投影会付出可量化的额外平方风险。

若 complex MAE 定义为 \(N^{-1}\sum_x|\widehat\chi_x-\chi_x|\)，每像素的 Bayes action 是 \(\mathbb R^2\) 中的几何 median；分别对 Re/Im 取 scalar median 优化的是分量 L1。论文必须先冻结到底使用哪一种 MAE。

### 3.3 弱方向可以只能恢复 prior，而不能发现 GT

若先验因子化且 \(\Phi(a,b)=\Phi(b)\)，则 \(\pi(a\mid y)=p_0(a)\)：没有算法能从这份数据确定那一次 GT 的弱坐标。正确 prior 平均仍可能比一个任意弱坐标点降低平均材料风险，但不可把这种统计收缩写成恢复了隐藏结构。

更一般地，另用 \((s,c)\) 表示 Gaussian prior 的正交分解，c 是拟用 prior 替代的未显式推断块。用精确边缘构造 \(\widetilde\pi(s,c)=\pi_s(s)p_{0,c}(c)\)，有标准 log-Sobolev 类型界

\[
\mathrm{KL}(\pi\|\widetilde\pi)
\le\tfrac12\mathbb E_\pi\|\nabla_c\Phi\|^2.
\tag{8}
\]

这是 certified LIS 理论的特殊情形 [R3]，不是新定理。它说明要保留哪些方向应由 posterior 区域的影响决定；anchor 的小奇异值不能代替右侧的 posterior 平均量。截断区域、uniform box 或旋转后条件支持需要对应 conditional reference 的 log-Sobolev 常数，不能直接抄 (8)。结合 bounded material、Pinsker 与 (5)，可把该 density 误差转成材料均值/额外风险的保守上界。

线性 Gaussian 诊断中，prior-whitened data Hessian 的特征对 \((v_j,\lambda_j)\) 对材料 posterior variance 的贡献为

\[
s_j=\frac{\|Dg\,v_j\|_W^2}{1+\lambda_j}.
\tag{9}
\]

很弱但几乎不动材料的方向可以不值得精细建模；很强但材料放大量大的方向也不一定能冻结。实用分块可有三类：material-relevant weak \(a\)、likelihood-informed strong \(b\)、task-negligible complement \(c\)。对 \(c\) 用 reference 积分或有界近似比固定零更合理。

## 4. 至少 18 类具有实质区别的架构

| 架构 | 数学对象与机制 | 物理/几何主要成本 | 正确性及主要失败 |
|---|---|---|
| 1. 原 prior pCN + property reflection | 全维 \(\pi\) 或明确 \(\pi_C\)，prior-reversible MH | 每 proposal 一次 forward，原值缓存 | 固定 proposal 下 MH 严格；强似然方向导致步长小、慢混合 |
| 2. Elliptical slice + property Gibbs | Gaussian field 沿椭圆条件 slice，property 单独更新 | 每次接受需随机个 forward | 精确条件更新；强约束使 bracket 收缩、forward 数可能高 [R5] |
| 3. 冻结分块 pCN/MALA | weak 用 pCN，strong 用条件 MALA，uniform 独立反射或切片 | weak 一 forward，strong 一梯度+forward | 每块正确 MH 保持联合目标；冻结几何可能差但不产生目标偏差 |
| 4. LIS/DILI | 后验平均的 informed geometry + reference complement | pilot Jacobian/Hessian 与后续低秩更新 | 已有成熟路线，不能改名为新 Maxwell 理论 [R1,R2] |
| 5. 预条件 HMC | 全维或块条件 Hamiltonian trajectory | L 个梯度，可能复用每点 LU | 最终 MH 修正正确；强曲率、边界和梯度成本限制 |
| 6. 位置相关 MALA / RMHMC | metric 随状态变，按正反 normalized proposal 修正 | 双端 metric/logdet，RMHMC 更昂贵 | 正确测度、drift、积分器必须区分 [R6,R7] |
| 7. Laplace-collapsed weak MCMC | 对 (3) 的近似边缘采样，同时积分条件材料矩 | 每 weak 点条件优化、Hbb/logdet | 目标是近似 marginal；边界、skew、条件多峰使其失效 |
| 8. Pseudo-marginal collapsed MCMC | 用非负无偏条件积分估计替换 L(a) | 每 proposal K 个条件 importance forward | 正确扩展目标给精确边缘；权重方差使链黏住 [R8] |
| 9. Tempering / parallel tempering | \(p_0e^{-\beta\Phi}\) 的温度梯与交换 | 多链多温，所有费用计入 | 交换比正确仍需冷链混合；局部单 basin 未必值回成本 |
| 10. SMC tempering | 加权粒子群、温度递进、rejuvenation | N粒子×T温×mutation | 有一致性条件，有限N不是精确 posterior；权重塌缩 [R9] |
| 11. CNN/MLP 全协方差 Gaussian | 直接学习 \(\pi_w(a\mid y)\) 或近似密度 | 离线合成与anchor；在线近乎无 physics | proper NLL 在总体上成立；Gaussian 不足、训练过拟合 |
| 12. MDN | 若干 Gaussian 的归一化 mixture | head 随 Kd² 增长 | 能表达局部 skew/弯曲，component 不是物理模式 |
| 13. Conditional normalizing flow | 可逆 reference-to-weak transport，显式 Jacobian | teacher/simulation训练；评估便宜 | 明确归一化近似；对漏失区域无自带保证 [R10,R11] |
| 14. EBM / ratio estimation | \(p_0e^{-E_\theta}\) 或 joint/product ratio | 训练负相/partition、可能 physics | score不是自动正规化密度；候选集softmax只给离散概率 |
| 15. Variational posterior | 最小 \(KL(r_\theta\|\pi)\) 的 ELBO | 每训练样本forward+adjoint，或surrogate | 不需 Z；reverse KL/窄族可能漏质量，与 supervised NLL 不同 |
| 16. Score/diffusion/flow-matching | 学习整条分布 transport | teacher或score条件期望，多步生成 | 网络/离散化误差不自动有正确目标；可后加IS/MH |
| 17. Foliation / trajectory ensemble | 对横向坐标积分，并沿条件轨迹求矩 | 轨迹查询、tube/atlas交叠密度 | 少量零宽曲线无法代表全维连续 posterior，需横向测度 |
| 18. RL policy + MH/SMC | 学方向、步长、refresh、stop 的计算决策 | 训练rollout + 在线物理；正反policy评估 | 奖励不定义 posterior；完整 proposal/路径权重不可省 |
| 19. Teacher–student 条件矩蒸馏 | 学 \(\pi_w,m(a),S(a)\)，后加校正 | 多对象可信teacher是主要成本 | 不能把探索候选当teacher；教师偏差会被忠实蒸馏 |
| 20. PDE/物理代理 + exact correction | 学或构造 forward，再用 delayed acceptance/残差界 | setup、bound、fallback，全计 | 近似只影响proposal时可保正确目标；未经界/修正会偏置 |

其中 1–10 是本轮概率 reference 主力，11–13 是朴素学习主力，18 应以成本改进为研究目标。16 的模型复杂度不应先于 posterior reference 的质量。

## 5. 主方案一：分块 Bayesian 推断与 Rao–Blackwell 材料成像

### 5.1 先实现一个真正正确、容易审查的全维基线

原坐标 proposals 为

\[
z'=\sqrt{1-\beta^2}z+\beta\xi,\quad\xi\sim N(0,I),
\qquad u'=\operatorname{Reflect}_{[0,1]^2}(u+\varepsilon),
\]

其中 \(\varepsilon\) 的分布固定、对称。两者乘积对原 prior 可逆，因而

\[
\alpha=1\wedge \exp\{\Phi(q)-\Phi(q')\}
\]

即可用于全后验。若目标是 \(\pi_C\)，还要拒绝不在 C 的提案。reflect 不等于 clip：把越界值压到端点会制造原 prior 中没有的原子。

冻结 Gaussian 正交 frame 后，weak pCN、strong MALA 与 property reflection 可以按固定顺序组成一次 sweep。每块以真实条件 posterior 做正确 MH，联合目标保持不变；不要求每块都完全混合，也不要求组合核可逆。强方向不能被最终 deterministic corrector 覆盖，否则会把分布压到低维集合。

```text
输入 y、原 prior/decoder、预先定义 C；只用 y 得到 q0
warmup：选择 Gaussian weak frame，估计各块尺度，随后冻结
初始化至少 4 条分散链；缓存各自完整 forward
for production sweep:
    a' <- 对原 Gaussian 条件 reference 的 pCN proposal
    用 full likelihood 和 C 判据执行 weak MH
    b' <- 冻结预条件 MALA proposal
    用完整正反 Gaussian density 执行 strong MH
    u' <- 2D reflection / slice proposal
    用 full likelihood 执行 property update
    记录当前状态，包括拒绝后保留的重复状态
    按预先确定节奏计算材料条件矩或直接材料样本
检查材料函数的 split-Rhat、bulk/tail ESS、MCSE 和跨链一致性
输出 posterior mean/median、条件矩误差、总费用和未认证项
```

**成本：** pCN 每步通常 1 新状态、8 complex RHS 和 1 新 factorization；gradient/MALA 的实际 forward/adjoint RHS 按实现计账，不能把 gradient 当免费。生产阶段谱冻结使 eigendecomposition 成本主要在 setup。具有同等缓存权利的普通 Full/MALA 是主强基线。

**决定性实验：** 用同一批新对象比较固定 b、优化 b、联合 Bayesian 三种估计；同时报告材料均值差、strong 条件材料 variance 占比和真实 MSE。若 strong 条件方差显著而 ridge 版本区间偏窄，便能证伪“优化即可替代积分”。若联合 posterior mean 与收敛 deterministic estimator 没有超出 MCSE 的风险差异，先保留校准价值，不扩大网络。

### 5.2 局部 Laplace 只用来加速，不暗中替换 exactness

Laplace conditional 可有三种用途，必须命名：

1. **近似成像器：** 对 (3) 做 weak MCMC，结果是 Laplace-local posterior，必须给近似误差；
2. **proposal：** 用 Gaussian conditional 提议 b，然后完整 MH，目标仍是原 posterior；
3. **importance reference：** 用 \(r(b\mid a)\) 计算 \(L(a)\) 的非负无偏估计，进入 pseudo-marginal；目标边缘可严格保持，但成本和方差可能高。

第三种令

\[
\widehat L(a)=\frac1K\sum_{j=1}^K\frac{e^{-U(a,b_j)}}{r(b_j\mid a)},\quad b_j\sim r(\cdot\mid a).
\]

提议 a' 时重新产生一组内层样本，接受比用 \(\widehat L(a')/\widehat L(a)\) 和 weak proposal 比；拒绝时保留当前那组内层样本与当前 \(\widehat L\)。每次把分母重新抽一遍通常不再是正确 pseudo-marginal，不能以“估计无偏”豁免。[R8]

对 66 维问题，减少 outer 维数不代表内层积分免费。d=2 的 weak 网格仍可能每格需要 64 维条件积分；如果其他坐标固定，只能叫二维条件切片 reference。Stage I 应用可解析 Gaussian 模型校验测度，再用 full 66 维长链的边缘作数值 reference，或明确生成一个低维 conditional-submodel toy。三者不能混报。

## 6. P1 特化的重要构造：条件 Gaussian 下精确积分 C² decoder

这是比“更大的 Flow”更值得优先检验的结构。两项 uniform properties 始终留在 outer 坐标 \(x=(a,u)\)，强条件只对 Gaussian field 的剩余 \(b\) 近似：

\[
r(b\mid x)=N(m_b(x),C_b(x)),\quad z=V_wa+V_sb.
\]

它可以来自条件 Laplace、Gaussian teacher fit 或 linearized posterior。**这不是声称真实 Maxwell 条件后验 Gaussian；下面“精确”仅针对给定 Gaussian reference 的 decoder 积分。**

### 6.1 每像素只需要一维截断正态矩

因为 P1 的 pre-decoder field 是线性的，令第 i 像素的 mode 行向量为 \(k_i^T\)，则

\[
f_i\mid x\sim N(\mu_i,s_i^2),\quad
\mu_i=k_i^T(V_wa+V_sm_b),\quad
s_i^2=k_i^TV_sC_bV_s^Tk_i.
\]

在 ROI 内，P1 fraction 是

\[
\phi(f)=\begin{cases}0&f\le-1/2,\\
p(f)=\tfrac12+\tfrac{15}{8}f-5f^3+6f^5&|f|<1/2,\\
1&f\ge1/2.\end{cases}
\tag{10}
\]

这是项目 \(h(f/0.5)\) 的精确多项式展开。设 \(\alpha=(-1/2-\mu)/s\)、\(\beta=(1/2-\mu)/s\)，标准正态 pdf/cdf 记为 \(\varphi,\Phi_N\)。未归一化的截断标准矩为

\[
T_0=\Phi_N(\beta)-\Phi_N(\alpha),\quad T_1=\varphi(\alpha)-\varphi(\beta),
\]
\[
T_k=\alpha^{k-1}\varphi(\alpha)-\beta^{k-1}\varphi(\beta)+(k-1)T_{k-2}.
\tag{11}
\]

原始截断矩

\[
M_n=\mathbb E[f^n\mathbf1_{|f|<1/2}]
=\sum_{k=0}^n{n\choose k}\mu^{n-k}s^kT_k.
\]

若 \(c=(1/2,15/8,0,-5,0,6)\)、\(d=c*c\) 为多项式卷积系数，则

\[
\mathbb E_r\phi=1-\Phi_N(\beta)+\sum_{n=0}^{5}c_nM_n,
\quad
\mathbb E_r\phi^2=1-\Phi_N(\beta)+\sum_{n=0}^{10}d_nM_n.
\tag{12}
\]

\(s=0\) 直接使用 decoder 值。极端尾部必须用稳定的 cdf/survival difference、对数尾概率和消去误差控制；不能让矩公式的浮点负 variance 被误解为理论失败。

### 6.2 条件材料均值、MSE trace 与 Jensen 误差

记 \(c_u=A+iB\)，则

\[
m_{r,i}(x)=\chi_b+c_u\mathbb E_r\phi_i,
\quad
v_{r,i}(x)=|c_u|^2\{\mathbb E_r\phi_i^2-(\mathbb E_r\phi_i)^2\}.
\tag{13}
\]

全材料 RMS 对应的条件 MSE trace 为 \(24^{-2}\sum_i v_{r,i}\)，因此无需内层 decoder Monte Carlo。完整空间 covariance 则涉及二元 Gaussian 的截断多项式积分；(12) 不自动给出像素间 correlation、连通事件 probability 或 simultaneous credible region。

\(m_r-g(V_wa+V_sm_b,u)\) 是可直接算出的 decoder Jensen discrepancy。它在 pure saturation/transition boundary 附近可能很大，而且不会被一个 latent Gaussian 的均值图表示出来。可检验机制是：在足够准确的 conditional Gaussian reference 下，先把这部分 cheap decoder integration 做对，再比较是否真的需要更复杂的 weak density。

Hessian Cholesky 后可用 \(s_i^2=\|L_b^{-1}(V_s^Tk_i)\|^2\) 计算所有像素方差，额外代价约 \(O(256s^2)\)，moment recurrence 每像素仅到 10 阶；这部分没有 Maxwell RHS。Hessian 与 conditional mode 的成本仍必须计入，不能因 (12) 免费就称整个算法免费。

条件 fraction 的 marginal CDF 也可由单调 p 的反函数和一维 Gaussian CDF 得到，并含 \(\phi=0,1\) 的原子。对 outer samples 取 mixture 后可计算 pointwise intervals；这些区间只对相应 Gaussian-conditional approximation 成立。

### 6.3 一个完全不需要物理计算的精确 prior 对照

原 \(f_i\) 是零均值 Gaussian，(10) 满足 \(\phi(-f)=1-\phi(f)\)，所以 ROI 内 \(E_{p_0}\phi_i=1/2\)。properties 与 field 独立，\(EA=1.3, EB=0.105\)，故原 prior 的材料均值精确为

\[
\mu_{0,i}=\begin{cases}0.75+0.0925i&i\in ROI,\\0.1+0.04i&i\notin ROI.\end{cases}
\tag{14}
\]

在这个特别对称且 property 独立的 prior 中，\(g(Eq)\) 恰好也等于 (14)；这是已证明的先验对称性结果，不能推广为 \(g(E[q\mid y])=E[g(q)\mid y]\)。不需要花 Monte Carlo 估计 prior mean。prior pointwise variance 同样可用 (12) 和
\(E(A^2+B^2)=1.786233\ldots\) 计算。

跨 joint prior predictive 的平方风险严格满足

\[
R(\mu_0)-R(\mu(Y))=\mathbb E\|\mu(Y)-\mu_0\|_W^2.
\tag{15}
\]

这给出“数据实际为材料均值增加了多少信息”的 exact 基准；单靠宽 posterior、候选多、甚至标准 rank-SBC 通过都不能替代它。评价必须包括 (14)、physics-only anchor、optimized MAP、posterior mean、受限 Bayes action。

### 6.4 无偏积分的 control variate，以及 pseudo-marginal 的正确材料均值

设 outer x 的未归一化目标为 \(L(x)=\int\gamma(x,b)db\)，取 normalized Gaussian reference \(r_x\)，令
\(w=\gamma(x,b)/r_x(b)\)、\(h(b)=g(\Psi(x,b))\)、\(m_r=E_{r_x}h\)。由 (12) 已知 \(m_r\)。用独立训练/cross-fitting 的标量或矩阵 \(C_x\) 构造

\[
\widehat N_C(x)=\frac1K\sum_j\{w_jh_j-C_x(h_j-m_r)\},
\quad \widehat L(x)=\frac1K\sum_jw_j.
\tag{16}
\]

在 reference 抽样下 \(E\widehat N_C=\int\gamma h\)、\(E\widehat L=L\)，但一般
\(E[\widehat N_C/\widehat L]\ne N/L\)。若用同一小批数据拟合 C 而不作 correction，连 numerator 无偏性也可能失去；应独立 pilot 拟合或 cross-fit。

然而正确 pseudo-marginal 扩展目标为

\[
\widetilde\pi(x,B)=\frac{\widehat L(x,B)\prod_{j=1}^Kr_x(b_j)}{Z},
\quad Z=\int L(x)dx.
\]

若对相关 outer x 有 \(\widehat L(x,B)>0\) reference-almost surely，并满足支持与可积性条件，则在它的 stationary 分布下，精确有

\[
\mathbb E_{\widetilde\pi}\frac{\widehat N_C(x,B)}{\widehat L(x,B)}
=\frac1Z\int E_{B\sim r_x^K}\widehat N_C(x,B)\,dx
=\mathbb E_\pi h.
\tag{17}
\]

严格正性下，证明由分母与扩展目标中的 size-bias 因子抵消，再用 (16) 的零均值 control variate 得到。\(\widehat N_C\) 可以有负分量；likelihood estimator 本身仍须非负，不能给它加可能变负的普通 control variate 后还声称 pseudo-marginal 正确。初始化需正权重，拒绝时保留当前整组内层样本和权重。

**只有非负性时必须另审零权重事件。** 原始 pseudo-marginal 的 outer 边缘仍正确，但 (17) 不是自动成立。扩展目标在 \(\widehat L=0\) 处质量为零，实际 cancellation 只能得到

\[
\mathbb E_{\widetilde\pi}\frac{\widehat N_C}{\widehat L}
=\frac1Z\int E_r[\widehat N_C\mathbf1_{\{\widehat L>0\}}]dx.
\]

因此需要额外证明 \(\widehat N_C=0\) 在零权重事件上成立，或至少每个相关 x 的 \(E_r[\widehat N_C\mathbf1_{\{\widehat L=0\}}]=0\)。不能只在零目标质量处任意定义比值而忽略 reference 积分中的缺项。LOCAL-C 配未截断 Gaussian reference 时，K 个条件点可能全部位于 C 外，CV numerator 一般不随全部权重一起归零。

反例：取 K=1、\(r=N(0,1)\)、\(\gamma(b)=r(b)\mathbf1_{\{b>0\}}\)、\(h(b)=b\)、\(m_r=0,C=1\)。则 \(\widehat L=\mathbf1_{\{b>0\}}\)、\(\widehat N_C=b(\mathbf1_{\{b>0\}}-1)\)。在 PM 扩展目标内 b>0，CV ratio 恒为 0，而正确条件均值是 \(\sqrt{2/\pi}\)。FULL 正 likelihood 且无强坐标硬截断可满足严格正性；一般 LOCAL-C 不能直接沿用此保证。若不能证明零事件条件，应停用该 CV 比值，恢复 \(\widehat N=K^{-1}\sum_jw_jh_j\)：这个普通 numerator 在全零权重事件上自动为零，保留正确 PM 材料矩身份。

此构造有一个明确增量：利用 P1 decoder 的闭式材料矩降低“材料 numerator”的方差，同时用真实 Maxwell 权重修正 Gaussian conditional approximation。它不需要声称 conditional posterior Gaussian；仍需检查 likelihood estimator 的 log-variance、链黏滞、K 值和总体 RHS。若 Maxwell 使 reference 极差，(16) 不能拯救 \(\widehat L\) 的权重崩塌。

使用 control variate 的 \(\widehat N_C/\widehat L\) 是 moment observable，不是一个 posterior material sample。有限样本下它可能超出材料取值范围；不能把这些值当物理候选、用其频数构造 credible regions。需要 posterior samples 时，应在对应正确扩展目标中按 importance 权重抽选内层真实样本，或使用独立正确 mutation；moment 与 sample 输出要分开。

**创新边界：** 截断 Gaussian 矩、Rao–Blackwell、importance control variates 和 pseudo-marginal 都是标准工具；Tallis 的截断正态矩研究早于本任务数十年，Casella–Robert 等已有采样 Rao–Blackwell 化 [R19,R20]。值得研究的是 P1 饱和材料 law 让条件材料 moments 可精确处理、与 full-Maxwell correction 结合后的风险/费用增量。若只在此人工 prior 有效，论文应明确 law-specific，而不能写成一般三维 Maxwell theorem。

## 7. 主方案二：真正可训练的 CNN + MLP Gaussian/MDN baseline

### 7.1 在一个新观测只有一份 GT 时，proper 监督仍成立

从原 joint law 独立生成 \(q\sim p_0,y\sim p(y\mid q)\)，定义只由 y 决定的 context \(c(y)\) 和 weak label \(a_y(q)\)。总体 log-score 为

\[
\mathcal L(\theta)=E[-\log r_\theta(a_y(q)\mid c(y))]
=H(A\mid C)+E_C KL(\pi(A\mid C)\|r_\theta(A\mid C)).
\tag{18}
\]

若 context 包含完整 y，且其余输入都是 y 的 deterministic function，理想 optimum 就是所定义 weak coordinate 的真实条件分布。每份 y 只有一个 GT 并不构成理论障碍；有限数据下高容量网络记忆每个 y 并把 variance 缩成零则是实际严重风险。[R10,R12]

若仅保留 compressed summary，模型学习的是 \(\pi(A\mid C)\)，除非证明 summary 充分，否则不能称其等于 \(\pi(A\mid y)\)。CNN 图像、J 和 anchor 都不会在已给定完整 y 后增加新的统计信息，但能极大改善有限容量网络的表示和样本效率。CNN 是否有价值必须通过 matched-parameter MLP ablation 判断。

### 7.2 一个约 0.5M 参数的确定配置

第一版用固定 train-derived Gaussian frame，d=16，避免每个对象独立特征向量引入 gauge 混乱。所有 object-level train/validation/test 划分在计算 anchor 前冻结。

| 输入/模块 | 形状与层 | 参数数（带bias、不计无参数激活） |
|---|---|---:|
| 材料 CNN | Reχ0、Imχ0、φ0、ROI：4×24×24；3×3 stride2 conv 4→16→32→64；flatten 576→128 | 97,584 |
| 数据 MLP | y real-stack512 与 whitened residual512 拼接：1024→256→128 | 295,296 |
| 物理 context MLP | q0(66)、冻结谱(66)、projected gradient(66)、material variance proxy(66)：264→128→64 | 42,176 |
| 融合层 | 128+128+64=320→128 | 41,088 |
| Full Gaussian head | d+d(d+1)/2=152 outputs | 19,608 |
| 合计 | Gaussian \(\mu\) 与 Cholesky covariance | **495,752** |
| 可替代 MDN head | K=4；K[1+d+d(d+1)/2]=612 outputs | 78,948 |
| MDN合计 | mixture weights softmax、4均值、4Cholesky | **555,092** |

Gaussian diagonal 经 softplus 加一个固定、数值合理的 floor；floor 和正则强度仅用 validation 选。MDN 的 component 数不解释为物理结构类别。为了公平，MLP-only 对照应替换 CNN 而维持相近参数总量。

训练标签是原 prior GT 的 weak 投影，或经过可信 full posterior teacher 得到的 weighted samples；二者分开实验。不能把 P2 的成功更新轨迹频率当成目标密度。直接监督建议至少数千独立 prior parents，多份独立 noise realization 可以复用同一 truth forward，但每份 y 的 anchor/J 费用仍需支付。同一 parent 的所有 noise replicas、轨迹和 samples 必须留在同一数据划分。

### 7.3 动态谱坐标如何保持表示一致性

符号翻转、谱交叉、近重根旋转会令 a 标签跳变。最大分量取正的 sign convention 只能消除部分歧义，不能消除近重根簇内旋转。推荐两级处理：

1. 第一版冻结 frame；报告在当前 y 上是否仍弱，作为纯学习概率的隔离实验；
2. 动态版在 ambient 64 Gaussian 坐标输出 \(m_\theta,D_\theta,U_\theta\)，其中 \(D\succ0\)、\(U\in\mathbb R^{64\times8}\)，然后

\[
\mu_a=V_w^Tm_\theta,
\quad\Sigma_a=V_w^T(D_\theta+U_\theta U_\theta^T)V_w.
\tag{19}
\]

该 head 需要 64+64+64×8=640 outputs，总参数约 558,704。对 \(V_w\to V_wO\)，均值/covariance 自动按 \(O^T\mu,O^T\Sigma O\) 变换，表示同一个 latent material distribution。几何特征本身也应用 projectors、谱簇不变量或完整 ambient 对象，不能同时喂入任意有序 eigenvectors 然后宣称 gauge-invariant。广义非正交 frame 需额外 coordinate Jacobian；用其 span 的正交 QR 概率坐标能简化第一版。

### 7.4 density、score 与 physics loss 的不同任务

连续 density 的值不是点 probability，可大于 1；单点质量仍为 0。候选 scores 的 softmax 只在给定有限候选集上正规化。若候选从 proposal r 抽出，正确 importance weight 是 \(\gamma(q)/r(q)\)，若候选来自 quadrature，则还须乘 quadrature cell volume。重复候选和人为 farthest-first 去重都会改变采样设计。

Gaussian/MDN/Flow 有显式正规化；EBM 的 \(r_\theta\propto e^{-E_\theta}\) 有未知 partition，需要训练负相、MCMC、数值积分或合法 density-ratio 方法。不能靠最后一层 sigmoid 产生正确 ambient density。

从 joint GT 的 NLL 得到一个 posterior 网络后，再额外加 \(\lambda E_{r_\theta}\Phi\) 通常改变 optimum、过度收缩；把 posterior 网络作为 prior 又乘一次 likelihood 也会重复使用数据。若选择 variational 目标，应明确用

\[
KL(r_\theta\|\pi)=E_{r_\theta}[\log r_\theta-\log p_0+\Phi]+\log Z,
\]

未知 \(Z\) 不影响优化，但它与 (18) 的 forward KL 是不同训练任务。

### 7.5 训练选择和 local-family misspecification

只由 y 决定的 anchor-success 筛选，在给定 y 后不会进一步改变 \(\pi(q\mid y)\)，但会限制训练/评价的 y 范围。若成功标签使用 GT、oracle MAE 或“GT 真的落在此 surface 内”，则训练 law 被选择因子 \(P(S\mid q,y)\) 改变，必须单列或纠偏。

若对 GT 做 nearest-surface projection，再把所得 a 当目标，网络学习的是该投影映射的 pushforward law，未必等于原 chart 的 weak marginal。即使 weak label 能拟合很好，被丢掉的 strong residual 仍可能使材料重建偏差很大。必须同时记录 GT 到 surface 的距离、原 weak 投影能量、实际 credible support 和是否触及 chart 边界。

### 7.6 teacher distillation 的数据预算和使用方式

teacher 应以 independent observation 为统计单位。每个 teacher 必须保留链数、warmup、ESS/MCSE、失败与 bias audit，不得把每条相关链的几千点当成几千个独立训练 parent。建议先有 128–256 个 development reference observations 比较 Bayesian 方法；只有在至少数百个多样 teacher targets 通过诊断、且已积累至少 2,048 个独立 train parents 后才比较以上网络。更大规模可取 4,096 train、512 validation、512 final test parents，但这是拟议预算，不是已授权执行。

蒸馏可学习 \(r_\theta(a\mid y)\) 与条件材料矩 \(m_\theta(a,y),S_\theta(a,y)\)。教师有有限 ESS 时，应按每个 observation 近似等权，并记录有效样本，而不让混合很慢却保存点很多的对象支配 loss。教师的 Laplace 偏差不会因为学生 held-out NLL 好而消失。

最稳健的上线形式是学生给 normalized full-dimensional proposal（弱 density × 强 conditional density × property density），再做完整 importance weighting 或 independence MH。只有 weak rθ 而没有 strong proposal/conditional integral，不能直接计算全维 importance weight。

## 8. 主方案三/并行方向：预算感知 trajectory policy，概率由校正核保证

### 8.1 先区分三种“风险下降”

对同一份固定 y，在 prior 与 likelihood 不变时，真实 posterior covariance 和最小 Bayes risk \(\operatorname{tr}W\Sigma\) 已经确定。增加 samples 不能减少这份数据固有的材料不确定性，只能减少我们对 posterior moments 的计算误差。只有增加新测量才可能改变真正的 posterior uncertainty。RL 若奖励输出 covariance 越小越好，很容易学会错误收缩。

可以严格区分：

1. **物理不确定性：** 真实 \(\Sigma\)，不能由多付计算费降低；
2. **推断误差：** 近似 density、未混合、Monte Carlo error，可由计算降低；
3. **决策风险：** \(\|\widehat\mu-\mu\|^2\) 代表估计器距离 Bayes action 的损失，可通过更好推断降低。

RL 应主要优化第 2、3 项的单位成本改善。

### 8.2 十类奖励及其局限

| 奖励 | 合理用途 | 不会自动保证的性质/失败 |
|---|---|---|
| 高 likelihood | warmup、候选修复 | 容易去 MLE；不是 posterior occupancy |
| 高 log prior + log likelihood | 接近 posterior mode | 只优化点密度，忽略体积、不能代替 sampling |
| 输出 posterior variance 下降 | 仅在与可信 reference 比较“估计误差”时合理 | 直接奖励会 mode collapse |
| 材料 Bayes estimator 误差下降 | teacher 的 \(\|\widehat\mu_T-\mu_{ref}\|^2\) | teacher 偏差和有限样本噪声 |
| 轨迹可达/稳定 | 学步长、几何 refresh | 曲线全程低残差不等于 probability 正确 |
| 未探索区域覆盖 | 防局部停滞，探索 warmup | 可能过度访问低概率尾部 |
| material ESJD | 相比 latent ESJD 不易钻 decoder gauge | 跳跃大不保证独立样本/关键函数混合 |
| 每 RHS 的 MCSE 改善 | 直接优化 expensive inference 费用 | 短轨迹 MCSE 不稳定，必须用 teacher/离线评估校正 |
| terminal GT MSE | 学近似 Bayes decision | 最优决策器不等于 posterior sampler |
| 综合概率/曲率/成本奖励 | 学一个 proposal/control policy | 所有额外偏好需只改 proposal，或说明目标已变 |

PPO/SAC/offline RL 都可以用，但这里 action space 不大，contextual bandit、Bayesian optimization、固定 schedule 或简单阈值控制可能更省训练成本。建议先给每条物理轨迹可比较的 terminal reward：

\[
R_T=-\sum_{k\in\mathcal H}\omega_k
\frac{(\widehat E_T h_k-E_{ref}h_k)^2}{\operatorname{Var}_{ref}(h_k)+\varepsilon}
-\lambda\,\mathrm{Cost}_T,
\tag{20}
\]

其中 \(\mathcal H\) 包括材料均值的 PCA coordinates、整体损耗/介电常数、边界区域 fraction 与 nonlinear material test quantities；不要只用 latent coordinates。reference unavailable 时可用多链 batch MCSE、cross-chain disagreement 和 jump proxies，但只能称启发式 reward。

2025 的 RLMH 已把 MH 设计表为 RL 并研究 ergodicity；AISTATS 2026 后续进一步发现 acceptance-rate 和 ESJD reward 在其设置中训练信号不足，提出 contrastive-divergence reward。[R17,R18] 因此“PPO/SAC+MH”或“trajectory reward”本身不足为新贡献。P3-L 有价值的差异应是 Maxwell 求解成本、decoder secant/response curvature、Jacobian refresh 和材料函数 MCSE 的联合调度。

### 8.3 固定 policy 仍需正反 proposal 概率

设策略给动作 \(j\) 的概率 \(w_j(q)\)，动作 j 产生 normalized proposal \(Q_j(q,dq')\)。边缘 proposal 是

\[
Q_\theta(q,q')=\sum_jw_j(q)Q_j(q,q').
\]

正确 MH 用

\[
\alpha(q,q')=1\wedge
\frac{\pi(q')Q_\theta(q',q)}{\pi(q)Q_\theta(q,q')}.
\tag{21}
\]

可避免计算整个 mixture，改用 auxiliary action j 并在接受比中包含

\[
\frac{\pi(q')w_j(q')Q_j(q',q)}{\pi(q)w_j(q)Q_j(q,q')}.
\tag{22}
\]

只说“策略从几种各自正确的 MCMC 核里挑一种”不够。反例：状态 \(\{0,1\}\)，目标均匀；恒等核和 swap 核都保持均匀。若在 0 选恒等、在 1 选 swap，两状态下一步都到 0，目标立刻被破坏。state-dependent mixing weights 必须进校正或构造受证明约束的 lifted kernel。

history-dependent policy 或不断变化的 geometry cache 还带来 adaptive MCMC 问题。第一版最稳妥的实现是只在 warmup 学 policy/atlas，production 冻结参数与 lookup 规则；若 production 仍学习，则需 diminishing adaptation 与 containment 等适用证明，不能仅因每个冻结核有 invariant distribution 就宣称非齐次链正确。

### 8.4 deterministic corrector 与整条 trajectory 的合法性

把全维随机提案经过 \(b\leftarrow b_*(a)\) 的 many-to-one corrector 后，输出落在低维曲面；它对完整连续 posterior 的 proposal 可能是奇异的，简单乘一个 endpoint likelihood MH 比不能修复。可选路径：

1. 只在明确的 collapsed target 上使用 corrector，并处理 L(a)；
2. corrector 后加 normalized full-rank conditional noise，且能算完整 proposal density；
3. 构造可逆、可计算 Jacobian 的增广 trajectory map，用整条正反路径密度；
4. 用 SMC 在路径空间携带正确 importance weights。

若 policy 的动作是“提议多长 trajectory”，反向路径的动作选择、随机动量、方向变量、步长密度和非体积保持 Jacobian 都属于 proposal。L2HMC、A-NICE-MC 已展示 learned proposals 可与正确校正结合 [R15,R16]；P3-L 不应省略这些技术条件。

```text
训练期（独立 development observations）:
    用可信 teacher 得到材料函数参考值
    rollout policy：动作是核类别/步长/刷新请求/预算分配
    每次 proposal 仍执行 (21) 或 (22) 的正确校正
    terminal reward 用材料矩误差与全部物理费用
    学完后冻结 policy 与 cost model
测试期:
    用全新 y 初始化多链
    policy 根据当前状态给动作概率，生成 tractable proposal
    正反动作概率、proposal密度和full目标一起执行MH
    按迭代记录所有状态，不按访问频率再softmax
    报告 material ESS/second、总 risk、fallback 与训练摊销
```

wall-time stopping 与 state-dependent 计算时间可能影响“随机时刻观察到的当前状态”的分布；统计量按 Markov iteration 的状态序列计算，不以实际停留秒数额外加权。基于当前估计“看起来够好了”停止会引入有限样本选择效应；应报告 stop rule 并在独立 benchmark 上测量。

### 8.5 SMC 是 trajectory policy 的另一种合法接口

在固定 prior 下使用 \(\gamma_t(q)=p_0(q)e^{-\beta_t\Phi(q)}\)。标准 tempering incremental weight 为

\[
w_t^{(i)}\leftarrow w_{t-1}^{(i)}
\exp[-(\beta_t-\beta_{t-1})\Phi(q_{t-1}^{(i)})].
\]

随后进行 \(\pi_t\)-invariant mutation。若 policy 产生不保持 \(\pi_t\) 的非平衡 transition，需完整的 forward/backward kernel SMC weight，而不是继续沿用上式。SMC 可以学习温度、move 数、proposal、branch 资源分配；finite-particle 结果有 MC 误差，importance ESS 不能证明未漏掉大质量区域。[R9]

## 9. Trajectory、trajectory ensemble 与 surface 的概率比较

### 9.1 零宽曲线的限制是测度障碍

对全维连续 posterior，有限条 C¹ 曲线的并集质量为 0。一个仅支持这些曲线的分布与原 posterior 的 TV 为 1，正向 KL 为无穷；即使图片看起来丰富，也不能把它称完整 posterior 近似。

沿一条曲线定义 \(q(t)\) 的 arc-length density 是另一种受限概率。如果某 density \(\rho\) 是相对于弧长定义，则 t 密度为 \(\rho(q(t))\|q'(t)\|\)。用 \(e^{-U(q(t))}\) 在任意 t 上归一化，结果随 parameterization 改变，除非明确定义它的参考测度或补上适当 Jacobian。

### 9.2 少量轨迹可以准确积分某些统计量，但需横向概率

有用的组织方式是 foliation：以横向坐标 c 索引轨迹 \(q=\Psi(t,c)\)，

\[
E_\pi h=\int\int h(\Psi(t,c))\,
\pi(\Psi(t,c))|\det D\Psi(t,c)|\,dt\,dc.
\]

对 t 做高质量 quadrature、对 c 做 Monte Carlo/SMC 或 sparse quadrature，可以构成合法的混合积分器。一条轨迹的“权重”来自其横向 cell/tube 的积分，不来自沿线测了多少个点。若弱 posterior 是 d 维，则每条局部 trajectory 至少需要 d−1 维横向测度；把 strong directions 再积分掉的问题仍需单独处理。

少量曲线即使无法在 density 意义上逼近，也可能在材料均值上足够。若 \(g\) 在 local region 对某 metric 是 L-Lipschitz，\(|E_\pi g-E_\nu g|\le L W_1(\pi,\nu)\)。这属于目标统计量/transport 距离近似，而不是 density 或 credible set 完整性。若实际材料分布的 quantization dimension 较低，它可能有价值；这个维数必须从 held-out reference 验证，不能由弱 Jacobian 维数直接推断。

### 9.3 atlas overlap 的两个正确处理方式

方式一：预定 partition of unity \(\omega_k(q)\ge0\)、\(\sum_k\omega_k=1\)，定义 chart 分量 \(\gamma_k=\omega_k\gamma\)。分量权重是 \(Z_k/Z\)，需积分或相对 normalizer estimation。

方式二：各 chart/tube 有 normalized proposal \(r_k(q)\)，混合 \(r=\sum_k\alpha_k r_k\)，使用 \(\gamma(q)/r(q)\) 的 importance weights。分母自动处理重叠，不能改成仅除以当前采样 chart 的 density 后仍解释为同一 mixture 权重。

| 表示 | 能有效回答的问题 | 主要偏差 | 建议定位 |
|---|---|---|
| 单 trajectory | continuation、局部一维条件统计、Jacobian刷新测试 | 横向质量完全缺失、GT未必附近 | geometry/calculation primitive |
| 多 trajectory + 合法 tubes/foliation | 各区域矩积分、预算调度、非均匀 quadrature | transverse coverage、overlap/normalizer | 值得条件推进 |
| d维 weak marginal + strong conditional | 局部材料 posterior mean/covariance、credible intervals | C截断、Laplace、有限采样 | 当前最稳健主表示 |
| 完整66维 local posterior | 联合校准、检查weak restriction | full physics成本、mixing | 必要reference与fallback |

## 10. 概率验证不能再只看 SBC 或候选数

### 10.1 标准 rank-SBC 的零物理盲区

设 GT 从原 prior 产生，算法完全忽略 y，永远输出独立 prior samples。对只依赖 q 的标量 rank statistic，GT 与这些 samples 在 prior-predictive 平均下仍可交换，rank 可以是 uniform。因此单独的原始 parameter rank-SBC 不能证明算法使用了 Maxwell 数据。

Talts 等建立了 SBC 工作流 [R13]；Modrák 等的原始论文明确研究 test quantities 对敏感性的影响，并给出包括 posterior=prior 在内的盲区 [R14]。本轮至少包含：

1. 原 latent/property 与材料函数的 rank-SBC；
2. 数据依赖的 test quantities，例如 \(\Phi(q;y)\)、residual 的指定投影；
3. 相对 (14) 的真实 material MSE 改善，以及相对 prior 的 posterior contraction；
4. posterior predictive checks 和独立/保留 measurement 的预测；
5. 按 anchor 成功、谱、boundary proximity、decoder saturation 分层的 conditional coverage；
6. 对主要 material moments 的多链 MCSE、rank-normalized split-Rhat、bulk/tail ESS。

校准与 sharpness 必须一起看。高 residual predictive fit 本身也不能证明 latent posterior 正确，尤其大量参数可补偿时。新增 measurement holdout 必须从 acquisition design 预先冻结，不能看到误差后挑接收点。

### 10.2 局部可信区间不能冒充全后验覆盖

若 \(D_y\subset C_y\) 且 \(\pi_C(D_y\mid y)=1-\alpha\)，在 full posterior 下

\[
\pi(D_y\mid y)=(1-\epsilon_y)(1-\alpha),
\quad\epsilon_y=\pi(C_y^c\mid y).
\]

不知道 \(\epsilon_y\) 就不能给全 posterior 的 95% 声明。实际 GT 不在 local support 不等于其概率为零；应报告 full-cohort failure 和 support miss，保留 full-dimensional escalation 条件。

### 10.3 proper GT density 与材料风险有关系但不单调等价

在 bounded material 条件下，若 student 对正确 target 的 forward KL 为 \(\delta\)，其材料均值偏差受 \(D\sqrt{\delta/2}\) 控制，额外平方风险受 \(D^2\delta/2\) 控制。这说明总体 proper NLL 有真实统计意义。

但某模型 NLL 改善而 MAE/MSE变差完全可能：density 可在材料不重要的方向改善，在材料重要尾部恶化；MAE 与平方损失最优行动不同；MC samples 数不足也影响估计器。实际选择应同时看 proper score、material risk、coverage、成本，而非把单个 GT 密度提高视为成功。

## 11. 三条最值得推进的方案及决策条件

### 第一优先：material-aware blocked posterior + P1 条件矩积分

**机制：** strong 条件不确定性被保留，Gaussian reference 的 decoder moments 按 (12) 做精确 cheap integration，weak 概率由全维/分块 MCMC 或经验证 collapsed 方法建立。真正可能改善的是非线性材料平均和 Bayesian action，不是候选去重数量。

**最小决定性比较：** 原 prior exact mean；optimized deterministic；full joint sampler；ridge+profile；ridge+Laplace volume；Laplace+analytic decoder moments；exact-corrected importance/PM moments。每种比较配对相同对象并给 MCSE，报告 strong 条件 variance、Jensen discrepancy、relative error 与真实 MSE。

**扩大条件：** reference 不同方法一致、held-out MSE 增量超出 MCSE 与对象级抽样误差，或在 MSE无增量时能给显著更可靠的不确定性。若 profile/Laplace/真实边缘接近且 posterior mean≈MAP，应停止复杂生成分支，把论文转向计算成本/校准边界。

### 第二优先：conditional density + conditional moments 的轻量蒸馏

**机制：** 直接以 joint GT proper NLL 或统计合法 teacher 学 \(\pi_w\)，以 cheap analytic/reference moments 学 m,S；优先 full Gaussian/MDN，再决定 flow。用原始测量保持任务定义，physics context 提供归纳偏置。

**决定性比较：** matched CNN/MLP；diagonal/full Gaussian/MDN；固定/ambient-equivariant geometry；有无 actual-y 输入；teacher原始/解析conditional moments；无校正/相同预算MH校正。零物理 (14) 和仅先验/仅anchor网络必须在主图内。

**扩大条件：** 在足够多全新 parents 上改善 proper score、material risk 和 cold/amortized成本，并且 exact correction 没有把加速全部吃掉。P2 的 54 successful updates 不满足数据条件。

### 第三优先（条件并行）：geometry/physics 成本感知 policy + 正确核

**机制：** 学何时使用已有几何、何时刷新、选择步长/核、分配轨迹预算，目标是 material moments 的 ESS per total second，概率由 (21)/(22) 或合法SMC保证。

**决定性比较：** frozen schedule、简单curvature threshold、contextual bandit、RL；相同proposal family与缓存；训练wall/physics一起摊销。若 RL 仅把proposal多样性提高，却不改善 material ESS/费用或 risk，应停。

**论文界限：** generic learned MCMC、RL+MH、LIS、Rao–Blackwell、Laplace determinant 都已有成熟文献。最可能形成独立贡献的是 Maxwell-specific geometry 对 refresh/likelihood decision 的可靠控制，以及 P1 decoder 条件积分对材料 Bayes action 的量化作用。后者若只在特定 law 成立，应作为可解释结构化 benchmark 与算法模块，而不是强行宣称普遍理论。

## 12. 已核验原始文献与可引用链接

下面每项均由实际检索或打开得到；文献采用正式论文链接。每篇的概括刻意控制为短说明；本文的数学扩展和反例是独立推导，不能作为该文献原结论引用。

| 编号 | 原始论文与链接 | 直接相关内容 |
|---|---|---|
| R1 | Cui, Martin, Marzouk, Solonen, Spantini (2014), *Likelihood-informed dimension reduction for nonlinear inverse problems*, Inverse Problems 30,114015. https://arxiv.org/abs/1403.4680 ; DOI 10.1088/0266-5611/30/11/114015 | prior/likelihood 相对谱结构、LIS |
| R2 | Cui, Law, Marzouk (2016), *Dimension-independent likelihood-informed MCMC*, JCP 304:109–137. https://arxiv.org/abs/1411.3688 | informed subspace与Gaussian reference complement |
| R3 | Zahm, Cui, Law, Spantini, Marzouk (2022), *Certified dimension reduction in nonlinear Bayesian inverse problems*. https://arxiv.org/abs/1807.03712 ; DOI 10.1090/mcom/3737 | nonlinear dimension reduction 的非渐近误差与log-Sobolev工具 |
| R4 | Cotter, Roberts, Stuart, White (2013), *MCMC Methods for Functions: Modifying Old Algorithms to Make Them Faster*, Statistical Science 28:424–446. https://arxiv.org/abs/1202.0709 | Gaussian-reference preserving function-space MCMC |
| R5 | Murray, Adams, MacKay (2010), *Elliptical Slice Sampling*. https://arxiv.org/abs/1001.0175 | Gaussian prior下elliptical slice |
| R6 | Girolami, Calderhead (2011), *Riemann Manifold Langevin and Hamiltonian Monte Carlo Methods*, JRSS B73:123–214. https://doi.org/10.1111/j.1467-9868.2010.00765.x | 几何MALA/HMC成熟先例 |
| R7 | Xifara, Sherlock, Livingstone, Byrne, Girolami (2014), *Langevin diffusions and the Metropolis-adjusted Langevin algorithm*, Stat Prob Letters91:14–19. https://arxiv.org/abs/1309.2983 | 澄清 manifold 与Lebesgue invariant density |
| R8 | Andrieu, Roberts (2009), *The pseudo-marginal approach for efficient Monte Carlo computations*, Ann Stat37:697–725. https://arxiv.org/abs/0903.5480 | 非负无偏估计与扩展目标 |
| R9 | Del Moral, Doucet, Jasra (2006), *Sequential Monte Carlo Samplers*, JRSS B68:411–436. https://doi.org/10.1111/j.1467-9868.2006.00553.x | SMC路径目标、normalizer ratio、rejuvenation |
| R10 | Papamakarios, Murray (2016), *Fast ε-free Inference of Simulation Models with Bayesian Conditional Density Estimation*, NeurIPS. https://arxiv.org/abs/1605.06376 | 条件density/MDN式simulation-based inference |
| R11 | Rezende, Mohamed (2015), *Variational Inference with Normalizing Flows*, ICML32:1530–1538. https://arxiv.org/abs/1505.05770 | 正规化可逆流的VI目标 |
| R12 | Gneiting, Raftery (2007), *Strictly Proper Scoring Rules, Prediction, and Estimation*, JASA102:359–378. https://doi.org/10.1198/016214506000001437 | proper score、calibration与sharpness |
| R13 | Talts, Betancourt, Simpson, Vehtari, Gelman (2018), *Validating Bayesian Inference Algorithms with Simulation-Based Calibration*. https://arxiv.org/abs/1804.06788 | SBC基本工作流 |
| R14 | Modrák et al., *Simulation-Based Calibration Checking for Bayesian Computation: The Choice of Test Quantities Shapes Sensitivity*, Bayesian Analysis, DOI10.1214/23-BA1404. https://arxiv.org/abs/2211.02383 | parameter-only SBC盲区、data-dependent quantities |
| R15 | Levy, Hoffman, Sohl-Dickstein (2018), *Generalizing Hamiltonian Monte Carlo with Neural Networks*, ICLR. https://arxiv.org/abs/1711.09268 | learned HMC、保持校正与训练mixing目标 |
| R16 | Song, Zhao, Ermon (2017), *A-NICE-MC: Adversarial Training for MCMC*, NeurIPS. https://arxiv.org/abs/1706.07561 | learned可逆proposal的既有先例 |
| R17 | Wang, Chen, Kanagawa, Oates (2025), *Reinforcement Learning for Adaptive MCMC*, AISTATS. https://arxiv.org/abs/2405.13574 ; https://proceedings.mlr.press/v258/wang25b.html | RL用于MH策略与ergodicity条件 |
| R18 | Wang, Fisher, Kanagawa, Chen, Oates (2026), *Harnessing the Power of Reinforcement Learning for Adaptive MCMC*, AISTATS/PMLR300:631–639. https://proceedings.mlr.press/v300/wang26a.html | acceptance/ESJD训练信号局限，contrastive-divergence reward |
| R19 | Tallis (1961), *The Moment Generating Function of the Truncated Multi-Normal Distribution*, JRSS B23:223–229. https://doi.org/10.1111/j.2517-6161.1961.tb00408.x | 截断Gaussian矩的经典先例 |
| R20 | Casella, Robert (1996), *Rao-Blackwellisation of sampling schemes*, Biometrika83:81–94. https://doi.org/10.1093/biomet/83.1.81 | 采样中的条件平均与variance reduction |

本稿没有声称已经完成真实 posterior reference、新训练或大型实验。对 posterior 数值正确性的建议是未来 protocol；已存在 P2 结果不能被追认为这些新方法的 validation。

### 12.1 本稿完成的零 Maxwell 数学核验

用五组一维 Gaussian 均值/方差，对 (12) 的闭式 recurrence 与独立数值积分作比较，最大绝对差为 \(6.1\times10^{-13}\)；只核验标量 decoder moment 公式，不是 Maxwell 实验或 posterior 近似验证。网络三个配置的层级参数计数也已独立算术核对。除此之外没有生成新的物理结果。


## 实验配额说明

最终实验规模、训练时限及RHS cap以主报告第15章和EXPERIMENT_PROTOCOL.yaml为准。
