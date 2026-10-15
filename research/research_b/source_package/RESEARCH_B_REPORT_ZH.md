# Research B：事件感知 Maxwell 主动测量
## 从带补偿的材料多解，到值得付费的新观测

**研究日期：2026-10-14。状态：理论、文献审查与实验设计；未运行新的 Maxwell 实验。**

本报告与 Research A 并行。A 负责从基础重建生成物理可行的不同材料解释；B 负责判断哪些解释仍值得区分、何种新增观测最有价值、何时应停止。B 不要求 A 提供归一化的完整后验。当前可用基线是 Joint P2 的声明径向材料模型与参考样本。

本次实际完成：指定仓库关键文件核对、最接近文献检索、数学命题与反例、算法和 E0–E3 协议、两个 JSON Schema、无 Maxwell 调用的概率/代数代码及 11 项检查。**Joint P2 缓存上的 E0 数值排名尚未执行，E1/E2/E3 未执行，因此没有选出新的最优频率，也没有测得主动测量的 MAE 或成本优势。** 不把下述建议写成已验证结果。

---

# 0. 研究裁决

**建议推进 B，优先级高于继续把原始生成器打磨成“完整正确 posterior”；但论文主线不应是“Maxwell + EIG”。**

最值得检验的窄主张是：

> 对当前数据允许的材料分支，主动测量应破坏它们跨观测共享的补偿关系，而不是最大化某一张图的局部灵敏度；这种破坏必须落实为材料或下游决策风险下降，并足以支付测量、校准与计算费用。

建议两条主算法与一条条件性工程支线：

1. **参考后验驱动的补偿感知任务风险设计**：完整条件样本、共享 nuisance、材料 MAE/决策损失为主目标，event EIG 为辅助比较。
2. **候选集合驱动的鲁棒复合假设设计**：没有可靠概率时，做声明集合内的最坏情形/权重敏感性设计，配合漏候选挑战；不伪造真实 EIG。
3. **带效用区间的多保真排序**：用 A17/A23 结构减少筛选费用，只对竞争测量做全波复核；若完整成本不回本，去掉 ROM。

强 TAP 论文需要超过径向球台架：至少还应有一个真正破坏球对称性的矢量 Maxwell 测试族，以及共享校准/跨频失效的负控制。径向台架很适合建立可核查的理论反例和 E0/E1，但不能独自支撑一般开放三维成像的结论。

---

# 1. 已有证据与不能跨越的边界

## 1.1 Joint P2：物理见证、事件参考和下游收益分别成立

主来源为 `migodam/b1-probabilistic-maxwell-imaging` 的 `joint-p2-20261014` 分支，读取时 tree commit 为 `7ed4faef9bc3dd214d4de689c2203cc383f7495f`。文件及来源编号见 `SOURCE_AND_NOVELTY_LEDGER.md`。

同心球台架的材料系数为
\[
x=(z_c,z_s,\epsilon',\epsilon''),\quad
(z_c,z_s)\in[-3,1]^2,\quad
\epsilon'\in[2.2,2.8],\quad\epsilon''\in[.04,.16].
\]
系数先验独立均匀，几何解码含零空腔原子；换成半径均匀会换先验，而不是换积分坐标。观测为 \(kR=1.6\) 下八角度、两个偏振的 16 个复通道，实化为 32 维。实部/虚部各自噪声方差为 \(\sigma^2\)，\(\sigma\) 由固定参考响应冻结；不能混入 A23/A30 使用的“复方差为 \(\sigma^2\)”约定。

gain 为一个跨频共享的、正规截断实高斯变量：
\[
\nu\sim \mathcal N(1,.15^2)\ \text{truncated to }[.7,1.3].
\]
当前模型不是任意复增益或相位漂移模型，也没有验证真实宽带色散。

| 结构 | 已有测量 | 最佳已找到的 ε+gain 补偿白化距离 |
|---|---|---:|
| 中心空腔 \(r=.08R\) | \(kR=1.6\) | 0.035593 |
| 同上 | \(1.6+2.4\) | 0.107504 |
| 中心空腔 \(r=.25R\) | \(1.6\) | 0.900170 |
| 同上 | \(1.6+2.4\) | 2.467944 |
| 环形空腔 \(h=.08R\) | \(1.6\) | 4.450536 |
| 同上 | \(1.6+2.4\) | 8.022080 |

这些是**可行端点的 full-Mie 见证，不是全局补偿最小值**。第二频率 2.4 是预先指定的。小孔仍难分、环形结构相对容易分，均应保留，而不是只呈现大孔正例。[P2-AMB]

大空腔 case 中，单频 full-Mie residual-corrected reference 的 solid / medium central / large central 事件质量约为 .162540 / .071828 / .677226。单频到双频的 .677377→.966140 属于连续 proxy 的信息消融，不是全部 D 条件的 full-Mie 后验认证。[P2-REF,P2-GATE]

下游独立几何传感器采用
\[
c=\text{normalized void volume}+\eta,\qquad \eta\sim N(0,.01^2).
\]
完整带权参考后验的平均事件决策 regret .00007225，明显低于简单 Top-K 代表点的 .063901。此结果来自 **B1 reference，不是 C1 原始生成器**；它验证条件分布交接价值，不证明主动 EM 测量策略已经成立。[P2-DOWN]

## 1.2 C1：粗事件准确不等于测量设计可靠

C1-P2 full/diagonal 蒸馏的事件质量 MAE 约 .013215/.013028，但重要事件内部 conditional CDF 误差约 .141761/.161768。未来测量的似然可以对事件内部连续参数非常敏感，因此仅事件质量较准仍可能排错 acquisition。[C1-P2,P2-COVER]

高 ESS、合法重要性修正、所有粗事件至少命中一次，都不能保证稀有模式、事件内部尾部或模型外材料被覆盖。有限原子样本与连续参考之间的 TV 通常甚至是 1；不能把样本数量增加当成 TV 误差认证。

## 1.3 A17/A23/A24/A30 对 B 的实际价值

A17 当前 closure 证明的是指定条件下的局部材料更新保真度；其三组完整非线性对照的 child-wall 比约 2.68、1.74、1.46，不是加速。双侧材料注入与接收反馈可以指导 B 的近似算子，但不自动成为新观测的效用。[A17-CLOSE]

A23 全材料少次反馈的 H1/H2/H3 均失败；A24 rank16 压缩不足、共同质量区未达到。相对 Maxwell 恒等式保留，已证实的同质量加速不存在。[A23,A24]

A30 Sparse-64 的 35.69% 是 truth-assisted 离线表示容量；实际 L-weighted mean 退化，局部补偿测试也没有找到协议要求的近等价对。应分开“候选里有好图”“数据能选中好图”“测量真的不足”三件事。[A30]

---

# 2. 统一模型：材料、几何、共享仪器误差与行动

令 \(\theta\) 表示完整物理材料世界：材料参数、几何以及需要时的色散参数。几何可嵌入共同空间的材料场 \(\chi_\theta(r,\omega)\)，避免把移动网格的差分误当材料导数。令 \(\nu\) 为持久仪器 nuisance，\(\xi_a\) 为按明确动力学变化的短时误差。

历史为
\[
D_t=\{(a_s,y_s):s\le t\},\quad
\pi_t(d\theta,d\nu)=p(d\theta,d\nu\mid D_t).
\]
行动 \(a\) 包含频率、源、极化、接收配置、驻留时间、校准配置和成本。新观测模型为
\[
Y_a=f_a(\theta,\nu,\xi_a)+\varepsilon_a.
\]
适应性选择行动本身不额外乘一个与材料有关的“行动似然”：策略若只依赖已经条件化的历史，则条件化该历史后，行动选择不会增加独立物理证据。

## 2.1 正确的共享 nuisance 更新

独立热噪声只意味着**给定持久 nuisance 后**独立：
\[
p(z\mid\theta,D_t,a)
=\int p(z\mid\theta,\nu,a)\,\pi_t(d\nu\mid\theta).
\]
若已有材料边缘似然，等价实现是
\[
\boxed{\ell_a(\theta;z)
=\frac{p(D_t,z\mid\theta,a)}{p(D_t\mid\theta)}.}
\]
随后
\[
\pi_{t+1}(d\theta)\propto \ell_a(\theta;z)\pi_t(d\theta).
\]

**禁止**每个频率分别对同一个原始 gain prior 积分，再把这些边缘似然相乘；那是另一个允许 gain 每频重新抽样的模型。

两种合法实现：保存 \((\theta,\nu)\) 联合加权粒子；或对每个材料计算上述边缘似然比。只有材料样本、没有条件 gain 信息时，必须重建 \(\pi_t(\nu\mid\theta)\)，不能把 nuisance 直接忘掉。[P2-REF；R05]

若热噪声也有跨块相关性，则应先用联合高斯条件分布：
\[
m_{a|0}=f_a+C_{a0}C_{00}^{-1}(y_0-f_0),\quad
C_{a|0}=C_{aa}-C_{a0}C_{00}^{-1}C_{0a},
\]
再处理 nuisance。白化须针对联合或条件协方差，不能凭通道数量宣称增加了同等数量独立信息。[R12]

## 2.2 物理线性化与精确有限差异

对固定频率/离散化，
\[
j_s=A_\chi(e_s^{inc}+G_\omega j_s),\quad
L_\chi=I-A_\chi G_\omega,\quad
F_{\omega,s,r}=S_{\omega,r}j_s.
\]
\(A_\chi=\mathrm{diag}(\alpha(\chi))\) 是正确的本构/极化率算子，不应无条件用裸 \(\chi\) 代替。材料 Jacobian 为
\[
J_a=W_a\,\mathcal R\!\left[
 S_aL_\chi^{-1}B_{\chi,s}
\right],\qquad
B_{\chi,s}\delta\chi=
\mathrm{diag}(\alpha'(\chi)\delta\chi)E_{\chi,s}.
\]
\(\mathcal R\) 表示按合同实化并施加材料质量坐标。跨照明共享同一个材料变量。

对两个材料状态，令 \(E_2=e^{inc}+Gj_2\)。直接相减给出**精确有限 secant 恒等式**：
\[
\boxed{j_1-j_2=L_1^{-1}(A_1-A_2)E_2,}
\quad
\boxed{F_a(\chi_1)-F_a(\chi_2)
=S_aL_1^{-1}(A_1-A_2)E_2.}
\]
证明：把 \(j_1-j_2\) 代入 \(L_1\) 并使用两条状态方程即可。要求共同嵌入、该行动的相同 Green/接收映射以及 \(L_1\) 可逆。它不是 Born 截断，也不保证便宜；它说明大分支差异受**材料差、本构、总场、内部反馈、接收耦合**共同控制。

---

# 3. 最关键设计规律：先合并测量，再消去补偿

## 3.1 局部可补偿方向必须共享同一组参数变化

在共同局部状态，白化材料与 nuisance Jacobian 记为 \(B_0,N_0\) 和 \(B_a,N_a\)。当前局部材料不可辨识方向满足
\[
\mathcal N_0=\{h:\exists v,\ B_0h+N_0v=0\}.
\]
新旧联合不可辨识方向是
\[
\boxed{
\mathcal N_{0,a}
=\{h:\exists\text{同一个 }v,\ 
B_0h+N_0v=0,\ B_ah+N_av=0\}.
}
\]
所以
\[
\mathcal N_{0,a}\subseteq \mathcal N_0\cap\mathcal N_a,
\]
且可严格包含：右侧允许两个不同补偿，左侧要求一个物理仪器解释全部数据。

最小反例：
\[
y_0=\chi+\nu+\epsilon_0,\qquad z=\nu+\epsilon_a.
\]
第一条里 \(\chi\) 被 \(\nu\) 吸收；第二条单独对 \(\chi\) 没有直接灵敏度。但联合可估计 \(\chi=y_0-z\)。两条分别投影的材料 Fisher 都为 0，联合后为 \(1/2\)（单位独立噪声）。**不直接照射未知物体的 reference calibration 也可能具有很高材料价值。**

这是可直接用于 Maxwell gain/reference 的机制，不是声称一般线性代数为新理论。

## 3.2 正确的 nuisance-projected Fisher

将观测联合白化后一次性堆叠 \(B,N\)，对局部高斯 nuisance prior 精度 \(\Lambda_\nu\)：
\[
\mathcal I_{\rm eff}
=B^\top B-B^\top N(N^\top N+\Lambda_\nu)^{-1}N^\top B.
\]
无正规 prior 时使用相应可识别子空间上的伪逆；有限边界应使用可行切锥或实际边缘化，不能用无约束投影当完整 Bayes evidence。

**共享 nuisance 的 Schur complement 不能按行动分别算完再相加。** 同一个 prior 也只能计一次。若材料与 nuisance 当前后验存在交叉精度，应对完整联合精度更新后再取 Schur complement，而不是使用独立 prior 近似。

增加正确建模的数据与 A17“更换近似 Jacobian”不同：前者是统计实验增广；后者的曲率差含交叉项，可能不定。不能拿 A17 的负增益说明正确 Bayes acquisition 的期望信息会变负。

## 3.3 可检验的补偿冲突分解

对固定材料分支方向，局部补偿损失可写成
\[
q_b(v)=c_b+(v-v_b)^\top H_b(v-v_b),\qquad H_b\succ0.
\]
则
\[
\boxed{
\min_v[q_0(v)+q_a(v)]
=c_0+c_a+
(v_0-v_a)^\top(H_0^{-1}+H_a^{-1})^{-1}(v_0-v_a).
}
\]
证明：联合最小点
\(v_*=(H_0+H_a)^{-1}(H_0v_0+H_av_a)\)，代回并完成平方。

最后一项是**共享补偿不相容代价**。即使新频率本身的 raw contrast 不大，若它要求与旧数据完全不同的 gain/材料修正，也能产生联合区分能力。

此为标准二次优化恒等式的 Maxwell 解释。对有边界、非线性和多分支问题，仅保留
\[
\min_v(q_0+q_a)\ge\min_vq_0+\min_vq_a;
\]
在紧致集合且极小点存在时，等号要求存在共同极小点。局部 Hessian 分解不能跨越任意非线性分支。

还需区分：**每个假设世界可以有不同的 gain 值**；共享只要求“同一个世界跨所有数据块的 gain 一致”。把两个不同世界强制使用同一个未知 gain，反而会夸大可辨识性。

---

# 4. 代表图分离、完整事件分离与误差下界

## 4.1 完整事件应通过预测分布比较

设声明事件 \(E_j\) 构成全部已建模支持的分割，
\[
p_j=\pi_t(E_j),\qquad
q_j^a(z)=\int p(z\mid\theta,\nu,a)\,
 \pi_t(d\theta,d\nu\mid E_j).
\]
事件更新：
\[
p_j(z)=\frac{p_jq_j^a(z)}{\sum_kp_kq_k^a(z)}.
\]
内层平均同时包含事件内材料、形状、色散与 gain 分布。代表图的单点 likelihood 一般不能代替它。[P2-DOWN,P2-COVER]

多类 0–1 最优事前错误率为
\[
R_{01}^{+}(a)=1-\int\max_j p_jq_j^a(z)\,dz.
\]
两个等先验事件为
\[
P_e(a)=\tfrac12[1-\operatorname{TV}(q_1^a,q_2^a)].
\]
这明确比较的是**预测混合分布重叠**，不是两幅 medoid 的距离。

若事件内各样本具有相同已知高斯协方差，则代表预测的 TV 误差可由
\[
\operatorname{TV}(q_j^a,N(\mu_{\rm rep},\Sigma))
\le E_j\!\left[2\Phi(\|\Sigma^{-1/2}(\mu-\mu_{\rm rep})\|/2)-1\right]
\le E_j\|\Sigma^{-1/2}(\mu-\mu_{\rm rep})\|/\sqrt{2\pi}
\]
控制。这给出了“代表图什么时候足够”的可测条件；条件 nuisance 导致不同协方差时，不能直接使用此式。

## 4.2 最难补偿对：适合鲁棒下界，不等于 Bayes 平均

令 \(\mathcal C_j(D_t)\) 是仍符合历史数据与物理约束的事件内状态集。集合设计可考虑
\[
d_{jk}^2(a)=
\inf_{\eta_j\in\mathcal C_j,\eta_k\in\mathcal C_k}
\|\mu_{0,a}(\eta_j)-\mu_{0,a}(\eta_k)\|^2,
\]
其中 \(\eta=(\theta,\nu)\)，联合白化和共享 nuisance 约束必须正确。新行动单独的距离与历史联合距离都可报告，但选择须以条件实验的实际风险为准。

“找到一个很近的对”给出合法歧义见证；“多起点优化找到的最近对仍很远”不证明所有对都远。非凸优化得到的是这个 infimum 的**上界**，不是下界。

固定两点、高斯同协方差白化距离 \(d\) 的等先验错误率为
\[
P_e=\Phi(-d/2).
\]
代入已有端点，仅作噪声尺度解释：
\(d=.90017\Rightarrow P_e=.326325\)，
\(d=2.46794\Rightarrow P_e=.108607\)，
\(d=.107504\Rightarrow P_e=.478566\)。
这不是 P2 复合事件的 Bayes error；两种频率配置的优化端点也未必相同。

对材料距离为 \(\Delta\) 的两点，若损失为该材料度量，
\[
\inf_{\widehat\chi}\frac12
[E_1L(\widehat\chi,\chi_1)+E_2L(\widehat\chi,\chi_2)]
\ge \frac{\Delta}{2}(1-\operatorname{TV}(P_1,P_2))
=\Delta P_e.
\]
证明用三角不等式 \(L(d,\chi_1)+L(d,\chi_2)\ge\Delta\)，在两密度公共部分积分。它也给出包含此二点的最坏情形风险下界。但要据此宣布“全部可用测量都无法解决”，必须控制整套允许行动及有限预算/顺序历史，不能只看一个频率。

## 4.3 Fisher 什么时候能预测事件区分？

充分的局部工作条件是：两个事件的主要条件质量集中在同一个正则局部参数域；white response 的二阶余项相对噪声和事件分离裕量足够小；nuisance 共享结构正确且其局部解稳定；事件边界不切过高质量 ridge；没有近共振/分岔造成的大曲率；先验质量没有主导性失配。

例如 \(F\) 的二阶导数在半径 \(\rho\) 内满足 white 范数界 \(M\)，则一阶 secant 误差至多 \(M\|h\|^2/2\)。只有该余项远小于候选行动的预测距离差，Fisher 排名才有可解释的精度。

反例：\(y=\theta^2+\epsilon\)，事件为 \(\theta>0\) 与 \(\theta<0\)，对称条件分布。在 \(\theta=\pm1\) 局部 Fisher 非零且可很大，但两个事件的预测分布完全相同，符号 EIG 为零。局部导数刻画每条分支内部，不会发现远端镜像分支。

---

# 5. 到底优化什么？

## 5.1 主目标：预期材料/决策风险下降

指定决策空间 \(\mathcal D\) 与物理损失 \(L\)：
\[
R_L(\pi)=\inf_{d\in\mathcal D} E_\pi L(d,\theta),
\quad
U_L(a)=R_L(\pi_t)-E_{Z_a|D_t}R_L(\pi_{t+1}).
\]
**第一主终点建议为固定物理尺度下的体积加权复材料分量 MAE**：
\[
L_{\rm MAE}=
\frac{1}{V s_{\rm ref}}\int_\Omega
\big(|\Re(\widehat\chi-\chi)|+
|\Im(\widehat\chi-\chi)|\big)\,dr.
\]
\(s_{\rm ref}\) 在评价前冻结，不依赖 truth。损耗 MAE、支撑、边界、空腔/拓扑事件单列；不同物理单位不能未经声明直接相加。

无约束逐像素分量 MAE 的 Bayes 决策为各分量后验中位数；若损失是复平面模长，则每像素是二维几何中位数。平方误差对应后验均值；MAP 对应密度最大化；最大事件质量对应事件 0–1 损失。它们不是同一决策。

若只允许从 A 的 K 张物理候选中选一张，则应使用
\[
R_K(\pi)=\min_{k\le K}\sum_iw_iL(\chi_k,\chi_i).
\]
这是受限候选选择风险，不能称为所有材料场上的最优 MAE。无约束中位图还可能不满足材料拓扑、全局色散或指定物理 family；应把“图像 MAE 决策”和“物理可行候选输出”分开报告。

## 5.2 event EIG：适合事件，但没有损失语义

\[
U_E(a)=I(E;Z_a\mid D_t)
=H(p_E)-E_ZH(p_{E|Z}).
\]
它优化事件信息，而不是自动优化错误率、材料 MAE 或下游风险。

本次有限概率代码构造了明确反例。三事件先验 \((.49,.49,.02)\)。行动 A 完美识别稀有事件、但不区分两个常见事件；行动 B 对两个常见事件作微弱区分。结果：
- A 的 EIG = .098039 nats，0–1 风险下降 .020000；
- B 的 EIG = .004908 nats，0–1 风险下降 .049000。

因此 EIG 可以偏好 A，而分类风险偏好 B。不是哪一个目标“错误”，而是目标不同。也不能通过这个一般反例断言 EIG 在 Maxwell 中必然更差。

对嵌套事件粒度，
\[
I(E_{\rm fine};Z)
=I(E_{\rm coarse};Z)+I(E_{\rm fine};Z\mid E_{\rm coarse}).
\]
更细分割增加的可以只是 ridge 内任意阈值信息。建议按下游决定物理事件，固定分割版本，报告阈值平移、层级分割及其 action 排名稳定性；不要把十二 bins 称为十二物理 modes。

## 5.3 有用的上下界与不可保证事项

在同一正确模型下，因为测量后可以忽略新数据，
\[
0\le U_L(a)\le R_L(\pi_t).
\]
这是条件期望风险的价值，不保证每个 noisy realization 的真实 MAE 降低。

若 \(L\in[0,L_{\max}]\) 且只依赖任务量 \(T\)，Bayes risk 对 TV 的 Lipschitz 性及 Pinsker 给出
\[
U_L(a)\le L_{\max}\sqrt{I(T;Z_a\mid D_t)/2}.
\]
该式没有逆向的普遍正下界：大量信息可能落在不改变当前最优决策的方向。若 MAE 不仅依赖声明事件 \(E\)，不能把 \(I(E;Z)\) 放进这个界代替 \(I(T;Z)\)。

对于平方误差，条件总方差给出熟知恒等式：
\[
U_{\rm MSE}(a)
=E\|E[T\mid D_t,Z_a]-E[T\mid D_t]\|^2.
\]

另一个必须写入接口的事实是后验质量的 martingale：
\[
\boxed{E_{Z_a|D_t}[p(E_j\mid D_t,Z_a)]=p(E_j\mid D_t).}
\]
因此，在不知道 truth、从当前 posterior predictive 平均时，不应承诺“某事件的期望概率从 .4 上升至 .8”。应给随机更新后概率的 5/50/95% 分位、超过决策阈值的概率，以及期望熵/风险变化。知道 truth 的评价条件与在线 predictive 条件不可混用。

---

# 6. 候选与权重的可靠性分级

| 输入层级 | 合法目标 | 不能声称 |
|---|---|---|
| 只有候选集合 | 集合内 maximin/复合检验、声明设计先验、最坏 regret | 真实 posterior EIG、真实事件质量 |
| 有 proposal score/weights | proposal 下的预期效用；权重扰动敏感性 | proposal score 已校准为后验 |
| 已知 proposal 且合法 full-likelihood 修正 | 对声明 target 的有限重要性近似及数值误差 | 全模式已覆盖、零访问事件被排除 |
| 有独立受控 reference | 声明支持和数值误差内的 Bayesian EIG/任务风险 | 现实模型外概率为零 |
| 怀疑遗漏重要世界 | 扩支持、鲁棒分布集合、独立挑战/保守弃权 | 从已有 bank 推断未知质量的可靠数值 |

对已知 \(q\) 的候选，正确重要性权重形如
\[
w_i\propto p_0(\theta_i)L(D_t|\theta_i)/q(\theta_i).
\]
优化器找到的有限候选通常没有可计算 proposal 密度；不能仅把 residual softmax 改名 full posterior correction。

可采用明确声明的权重集合
\[
\mathcal W_\epsilon=\{w:\operatorname{TV}(w,w_0)\le\epsilon,\ w\in\Delta\}
\]
或多组物理先验，优化最坏未来损失/最坏 regret。此处 \(\epsilon\) 是稳健性假设，除非另有证明，不能称为估计出的 posterior 误差。

若未知质量用 contamination \((1-\epsilon)\pi+\epsilon Q\) 表示，必须声明 \(Q\) 的物理支持和 \(\epsilon\) 来源。完全任意的未知材料类可以令所有设计的 minimax 值退化；此时应承认无法保证，而不是编造 outside-model probability。

覆盖挑战至少包括：预先留出的结构 family；从防御性 prior 生成并做 full likelihood 筛选的候选；少量不由当前二图差异驱动的 sentinel measurements。异常残差只能提示模型/候选不足，不会自动指出遗漏的正确材料。

---

# 7. Maxwell 特有、可证伪的测量规律

## 7.1 跨频共享补偿冲突，而非“频率越高越好”

同一个材料/色散模型、同一个共享 gain 要同时解释多个频率。新频率若改变不同材料之间的频谱比、相位或高阶多极结构，就可能破坏原补偿。需要比较联合复合预测，而非单个频率上 medoid 的最大差。

反面条件：在严格的可分离偶极模型中
\[
F_a(\theta)=c_a\,\alpha_{\rm eff}(\theta),
\]
两个材料若有效极化率相同，则改变仅影响 \(c_a\) 的角度、偏振或低频比例因子不能提供新的无噪声结构辨识。这说明“增加频率”不是自动有效；需要打破这种响应因子化。

近共振可能使分支差异更大，也可能同时放大 gain/几何漂移、离散误差、噪声和 ROM 缺陷。应比较完整风险与稳定裕量，不把共振峰直接当 optimal frequency。

## 7.2 角度与偏振：先做对称性/信息张成检查

对各向同性同心球，全装置共同旋转只是已知坐标变换，没有新增无噪声结构信息。但**固定接收器只改变入射方向**会改变相对散射角，可能产生新函数值；不能把两者混为一谈。

若同一几何已经测得两个独立横向极化通道，其他线性极化组合不增加同一源/接收配置的无噪声张成空间。新硬件测量可以改变噪声、能量分配或动态范围，因此仍可能有成本/SNR 价值；在缓存里数字组合已有读数则不是独立新数据。

所以 E1 径向台架的 polarization 应首先是**冗余负控制**。要证明真正的极化结构辨识优势，需要非球对称/各向异性对象，或初始确实缺少一个极化通道，并另立预算。

## 7.3 接收器、近场和反馈

新的接收映射 \(S_a\) 可以使当前 receiver-dark 的辐射成分可见；不同照明则改变材料注入 \(B_{\chi,s}\)。二者都不能只由裸 current 奇异值判断。

近场可能对局部边界和小结构更敏感，但其位置误差、源耦合、动态范围、模型误差成本也可能更大。不同接收点还需更新伴随缓存；若可从已算全电流快速评价接收器，应把这当实在的复用优势，但不能把未知材料的全电流当免费已知。

材料重要方向需要看
\(S_aL^{-1}\mathrm{diag}(\alpha'\delta\chi)E_s\)，
不是只看 \(S_a u\) 或内部 \(G\) 的能量。

## 7.4 多源组合、互易与重复

固定线性介质状态下 Maxwell 对激励线性。已知完整 source-response bank 时，新 coherent superposition 不增加无噪声输入张成空间，但能在约束的发射功率/驻留时间/接收噪声下优化可测 SNR。不能重复计入数字合成通道的 likelihood。

互易材料与校准几何下，交换源/接收对可能给冗余响应；真正重复测量仍可以降低独立热噪声。重复不能消除持久校准 gauge。例如
\[
y_i=\theta+b+\epsilon_i,\quad b\sim N(0,\tau^2),\quad
\epsilon_i\overset{iid}{\sim}N(0,\sigma^2).
\]
平均观测的方差为 \(\tau^2+\sigma^2/n\)，信息趋于 \(1/\tau^2\)，而不是无界增长。直接校准 \(b\) 可以改变这个上限。

小空腔在一个有效模型中若响应差为 \(O(r^3)\)，独立噪声下达到固定检测能力所需重复次数通常为 \(O(r^{-6})\)。这是该幅度标度的推论，不是普遍 Maxwell 定律；P2 已否定当前数据支持普遍 \(r^5\) 补偿机制。[P2-GATE]

## 7.5 宽频材料模型

窄频基准可以明确假设 \(\epsilon\) 常数。更宽频须共享合适的物理参数，例如采用 \(e^{-i\omega t}\) 约定的被动 Debye family：
\[
\epsilon_r(\omega)=\epsilon_\infty+
\sum_l\frac{\Delta\epsilon_l}{1-i\omega\tau_l}
+\frac{i\sigma_c}{\epsilon_0\omega},
\quad
\Delta\epsilon_l\ge0,\;\tau_l>0,\;\sigma_c\ge0.
\]
这是可选的声明模型，并非所有材料都由少数 Debye poles 描述。要做 Lorentz/模型失配负控制，不能用一个过窄色散 family 制造虚假的唯一性。物理频谱参数共享本身提供约束；把每个频率的 ε 独立自由化会改变问题。[R18]

---

# 8. 三个算法

## B1. Reference-anchored compensation-aware task-risk acquisition

输入：当前联合后验/条件 nuisance、事件内样本、候选行动与成本。

对每个行动：
1. 在同一组联合粒子上计算响应；同频多接收器/照明按实际可复用部分缓存。
2. 从当前 predictive 模拟未来读数，不使用评价 truth。用 common random numbers 降低行动差的 Monte Carlo 方差。
3. 用正确 conditional likelihood 更新全部粒子；事件内分布与 remainder 一起保留。
4. 求新后验材料 MAE 决策或受限候选 medoid，计算 expected risk；同时报告事件错误率和 EIG。
5. 选择预算约束下的风险最优行动，或用明确的损失—成本换算选净价值最大者。对校准+测量做二步 lookahead。
6. 真实读数到来后再更新；若 predictive 冲突、权重退化或新模式挑战成功，扩充支持并重做设计。

对于 \(N\) 粒子、\(M\) 外层 draws、\(A\) 行动、\(q\) 数据维，直接 kernel 评价约 \(O(AMNq)\)。K 个候选决策的损失矩阵预计算后，每轮风险约 \(O(KN)\)。逐像素 L1 中位数可复用排序；复模长 geometric median 需要另计优化误差。低维 reference 可用 cubature 替代纯 MC。

**不需要训练新网络。** C1 的旧单频 C conditioner 不能未经适配用于新的 D/history；合法 reweight 是可先用的更新器，但有限样本仍可能需要 rejuvenation/结构扩展。

## B2. Robust composite-hypothesis acquisition with coverage challenge

没有可靠 posterior 时，使用 A 的候选及其允许补偿邻域。先做联合局部补偿冲突筛选，再对少量 finalist 做受约束非线性补偿与 full-Maxwell 核验。

可以优化声明集合内的最坏任务风险：
\[
a^*=\arg\min_a\inf_\delta\sup_{\pi\in\mathcal P}
E_{\pi,P_a}L(\delta(Z),\theta),
\]
或最坏 regret。必须保持 minimax 的顺序；“先给每个 prior 选自己的 oracle Bayes rule 再取最坏”是更乐观的另一目标。

对 \(d_{jk}\) 的 maximin 是筛选代理，不一般等同于以上风险。若多个分支共享难以打破的对，返回“在所测菜单与预算下保留多解”。同时分配明确的候选覆盖挑战额度，检查 bank 是否遗漏重要结构。

这条路线最适合 A 的早期输出；无需知道真实 event masses，但也不能返回精确 Bayesian 概率改善数字。[R01,R06,R08]

## B3. Multi-fidelity utility interval ranking（条件性）

A17 的 \(D_C=Y_C\Gamma_C^{-1}N_C\) 类双侧因子提示：近似资源应覆盖任务相关的材料注入与候选接收反馈两侧，而不是仅最小化裸 current 能量误差。其 paired GN 支持结论是固定局部模型改变的结构，不是新观测 EIG 的现成定理。

A23/A24 相对反馈在固定频率背景下给
\[
\Delta j=R_b w,\quad
w=(I-HG_b)^{-1}He_b,\quad
\Delta y=D_bw,
\quad G_b=GR_b,\quad D_b=SR_b.
\]
可对多个材料候选复用，但跨频 \(G,\alpha,S,e_b,R_b\) 都可能改变，必须重新验证/刷新。近共振、非正规放大、材料超出训练域时不能只依赖 residual。

资源应先投向可能改变排名的行动。若区间不能分开，增加 full predictions；若不值得支付，则输出“排名未分辨”，而不是选择小数点更漂亮的代理分数。

---

# 9. ROM 误差如何进入真正的 acquisition 排名？

设当前分布 \(\pi\) 固定，真实/近似 conditional observation kernel 为 \(P_a,Q_a\)，并定义
\[
\tau_a=\int\operatorname{TV}(P_a(\cdot|\theta),Q_a(\cdot|\theta))\,\pi(d\theta).
\]
对 \(L\in[0,L_{\max}]\)，新观测后最优风险作为 joint law 上的决策问题，有
\[
\boxed{|U_L(a)-\widehat U_L(a)|\le L_{\max}\tau_a.}
\]
证明：对任意同一决策规则 \(\delta(z)\)，两 joint laws 的期望损失差至多 \(L_{\max}\tau_a\)；分别取 inf，得到双向界。当前风险项相同，直接抵消。

若额外有实际可证明的同支持分布误差
\(\operatorname{TV}(\pi,\widehat\pi)\le\epsilon\)，则保守总效用误差增加 \(2L_{\max}\epsilon\)。**有限原子粒子与连续密度不能默认有小 TV**；此时需要量化积分误差、合适的 Wasserstein+Lipschitz 控制，或仅声明相同离散 bank 内的排名保证。

对同协方差高斯响应，若 white 均值误差为 \(\eta(\theta,a)\)，
\[
\tau_a\le E_\pi[2\Phi(\eta/2)-1]
\le E_\pi\eta/\sqrt{2\pi}.
\]
联合 nuisance 的混合可用 TV 的凸性，但必须控制同一个正确的条件核。固定响应 bank 的经验 holdout 最大误差不能直接升为 uniform \(\eta\)。

若最终误差半径为
\[
\beta_a=\beta_{\rm forward}+\beta_{\rm quadrature}
+\beta_{\rm weights}+\beta_{\rm decision}+\beta_{\rm MC},
\]
且都是同时有效的界，便可在
\[
\widehat U_a-\beta_a>\max_{b\ne a}(\widehat U_b+\beta_b)
\]
时确认排名；否则继续细化或保守弃权。直接选 \(\arg\max\widehat U\) 的 regret 至多 \(\beta_{a^*}+\beta_{\hat a}\)。若只有经验区间，名称应为 empirical ranking confidence，而非 certificate。

离散状态 residual \(r=L\widehat j-Ae^{inc}\) 的输出误差为
\[
F-\widehat F=-SL^{-1}r
\]
（按该 residual 符号定义）。小 \(\|r\|\) 需要 \(\|WSL^{-1}\|\) 或有效 adjoint/resolvent 控制才意味着小 white output error。全局 Neumann 范数界可能过于松，应先测是否能在真实候选上给出非平凡排序区间。

有限事件 EIG 的 entropy 表达有 \(\log J\) 上界，可另用有限 alphabet 条件熵稳定性；不能将任意连续参数 KL/EIG 当成有界损失直接套同一常数。一般 EIG、IPM 与 surrogate stability 已有广泛先例，本式的贡献应是与实在 Maxwell 误差预算结合，而非宣称首次给效用扰动界。[R09,R10,R15,R16]

---

# 10. Ambiguity Graph：可选压缩，不是必需新理论

图节点应是带条件样本的材料事件/分支，边表示在旧数据下存在任务相关的补偿歧义。边可附：材料距离、可行近邻见证、补偿约束、预测分布重叠、目前是否具有全局下界。

图能用于稀疏化候选对、找到未覆盖的结构群、分配 full-verification 次数和选择代表性 batch。但 pairwise overlap 的和不是多类 Bayes risk；一条不可能在预算内解决的边也可能绑架 maximin，使所有行动都显得无用。

**建议先不把图列为主贡献。** 固定同候选、同 full-cost 下比较“直接任务风险”与“图筛选+风险复核”。图若没有节省计算或改善遗漏检测，主动删除。不要只因为网络图看起来有物理解释就保留模块。

---

# 11. 与 Research A 的完整接口

两份机器可读 JSON Schema 在 `schemas/`。它们要求每个数组有格式、shape、dtype、单位与 provenance；没有实际数组时不能仅靠 JSON 验证器宣称物理有效。

## 11.1 MaterialHypothesisPackage

核心字段：

| 组 | 必需内容 |
|---|---|
| 身份 | schema_version、package_id、producer/source_revision、target_id、data_history_id |
| 物理合同 | solver/discretization、材料域/度量、时间约定、色散模型、参考频率 |
| 历史观测 | 数据数组、行动历史、实化顺序、联合噪声/跨历史协方差、reference noise |
| nuisance | 名称、proper prior、跨频/源/接收器 sharing graph、联合粒子或条件 posterior、漂移模型 |
| 材料 | 基础解、K 候选、物理结构说明、full-Maxwell 已有预测、材料/数据距离、物理可行性 |
| 补偿 | 每个候选允许改变的参数、边界、必须保留的结构事件、不能偷偷填洞的保护约束 |
| 概率 | none/design_prior/proposal/importance_corrected/independent_reference；log prior/q/L、测度、数值误差与适用数据 |
| 事件 | 版本、定义、质量、条件粒子/weights、代表候选、是否真正证明为分离 mode |
| remainder | 未展示但已建模质量及其样本/积分表示；outside-model mass 默认 null/unquantified |
| ROM | 使用与否、有效材料/行动域、输出误差性质、刷新规则与成本 |
| 行动与费用 | 可实施菜单、noise、持久 nuisance groups、独立新噪声标志、时间/能量/金额/RHS/分解次数 |
| 完整性 | 未探索结构、coverage audit、资源账本、限制 |

只有候选时 weights.kind=`none`，概率相关值为 null，不妨碍 B2 运行。若 Top-K 未包含全部权重，必须传 remainder；只把显示项归一化为 1 会改变 target。

接口中“每个材料的 full-Maxwell prediction”只指实际计算并有 provenance 的配置；未计算的新行动不能用一个 `full=true` 标志冒充全波结果。

## 11.2 AcquisitionRecommendation

输出包括行动或两步 bundle、推荐状态、概率可靠性、planning horizon、损失定义、当前/未来风险/增益区间、event EIG、更新后事件概率分位数组、预计解决与保留的歧义、校准风险、完整成本、效用误差分项、排名适用范围，以及正确 joint-likelihood 更新配方。

状态至少包括：
`RECOMMEND_MEASURE`、`STOP_TASK_SUFFICIENT`、
`STOP_NOT_WORTH_COST`、`PRESERVE_AMBIGUITY`、
`ABSTAIN_POSTERIOR_UNRELIABLE`、`EXPAND_SUPPORT_FIRST`、
`ABSTAIN_UTILITY_RANK_UNRESOLVED`、`BUDGET_EXHAUSTED`、
`SIMULATION_DESIGN_ONLY`。

“预算用完”“当前算法不能排序”“物理不可辨识”必须是不同状态。没有硬件型号/带宽/极化切换/动态范围/校准能力时，kR 菜单只能是 simulation design，不能称为现实仪器最优方案。

---

# 12. E0：不新增 Maxwell 调用即可启动

**硬上限：0。** 不训练新网络，不访问未授权 GPU，不借用 P2 剩余额度。

E0 的候选新观测只能来自已保存的 1.6/2.4 响应或明确的解析 context 模型。缓存没有的新频率/角度，返回 `NOT_AVAILABLE`，不能为了补图悄悄调用 solver。

## E0-A：合同与恒等式
先核对 P2 accepted reference、gain sharing、complex packing、prior pushforward、所有历史数据 ID。执行本包 11 项概率/代数检查。schema 不能替代数值精度审计。

## E0-B：条件预测风险，而非看已知第二频率的结果
以已存在的单频 observation 为条件，对缓存 2.4 的完整块及合法子块模拟未来读数，比较 raw medoid separation、event EIG、material risk、joint projected Fisher 和 robust composite screening。

另列 repeat-1.6（新独立噪声）与 reference calibration。后两者是模拟设计，不是假装仓库里已经做过对应硬件测量。

历史已经观测到的 2.4 数据只可做暴露案例的回顾性检查，不能作为“未来最有利噪声”的选择依据。**E0 无法在未缓存频率中寻找最优频率。**

## E0-C：代表点与事件内分布
相同当前 event masses，比较 medoid-only、每事件有限条件样本、完整带权 reference。报告 acquisition 排名差及独立 reference 下的 risk regret。这里直接检验 P2 的 conditional CDF 缺口是否会改变下一测量，而非只评价生成器本身。

## E0-D：对称性、校准与 prior stress
分别比较共享 gain 正确模型、每频独立 gain 的明确负控制，以及错误固定 gain；不能合并为公平算法胜负。按同样响应、不同合法 prior 做敏感性，而不是混成单一排名。

角度/极化 retrospective subset 只能通过**移除这些通道后重算 posterior**来评估；已经用于当前 posterior 的旧通道不具有第二份独立证据。

产物：`E0_ACTION_INVENTORY.json`、`E0_CRITERION_COMPARISON.json`、`E0_WEIGHT_AND_PARTITION_SENSITIVITY.json`、`E0_DECISION.md`。已有已找到难对包括中心孔 .08/.25，环形高距离是负控制；最终最难事件对排名须计算预测重叠，不以这张端点表代替。

---

# 13. E1：主动频率和极化的最小付费验证

本报告**建议**独立 E1 radial pilot 上限 4096 material-frequency attempts。它不是已经获得的执行授权，也不是从 P2 的 541 剩余额度自动转入。

冻结示例菜单：
\[
kR\in\{1.2,1.6\text{ repeat},2.0,2.4,2.8,3.2\}.
\]
该菜单来自预先规定的尺寸参数区间，不由 truth 或已有表现选择。仪器带宽未核实时标为 simulation-only；物理测试需把它替换为可实施且仍在冻结前确定的菜单。

| 分配 | attempt 数 |
|---|---:|
| 四个新频率各 625 个 response interpolation 节点 | 2500 |
| 每新频率 64 个独立 full-field audit | 256 |
| 32 个新 truth parent × 四个新频率 | 128 |
| finalist/关键 predictive 点 full verification | 512 |
| 失败、重新积分与修正预留 | 700 |
| 合计 | 4096 |

所有 failed/reserved attempts、实际 factorization、RHS、wall time 分开计账。若参考精度不能在预算内达到，E1 结论为 reference unresolved，不能降低标准、换噪声或把 proxy 排名改称 full target 胜利。旧频率的新 truth 响应若不在已有合法缓存，也须从预留支付，不能遗漏费用。

**第一轮先隔离 criterion。** 所有方法使用同一个候选/参考、同一 full/proxy predictor、同一成本模型；之后才比较 ROM，否则无法判断提升来自目标还是求解精度。

至少保留 random、uniform sweep、preset2.4、raw separation、经典 D-optimal、正确 joint projected Fisher、event EIG、expected material risk、proposed robust compensation、no-new-measurement。Fisher 在可在线取得的基点或 posterior 样本上评价，不用 truth-state Jacobian；参数尺度及 prior 在看到结果前冻结。

径向球已有双极化，先运行张成/旋转冗余负控制。要开展真正矢量极化选择，另申请最小非同心测试，不能把换线性极化基当成新增独立结构维数。

32 fresh parents × 两份噪声不能当成 64 独立物体。先冻结 policy/阈值，再生成新 parents；全部 policy 共用按 parent/action/repeat 编号的噪声 potential outcomes，持久 gain 每个 episode 只抽一次。新的现实噪声相关性和 drift 不得被这个合成 coupling 隐藏。

---

# 14. E2/E3：闭环与候选遗漏

## E2
最大三次新增测量块，允许 stop。流程为：
当前 package → conditional predictive design → 实际/模拟观测 → joint update → 支持检查/必要时 rejuvenation → 继续或交接。

至少对 calibration+frequency 做二步 lookahead。原单频生成器不能接收新历史就假装仍是同一个 calibrated posterior。候选更新若只 reweight，不会自动从零权重支持中创造新 mode；必须保留扩支持路径。

报告每一步和最终：真实全材料 MAE、实/虚分量、空腔/边界/support、Top-K 的参考质量覆盖与真实材料近邻覆盖、事件分类、过度自信、漏候选、Bayes/reference regret、实际噪声导致的错误更新、测量块数、通道数、积分时间、setup/online/full simulation cost。

前瞻 policy 的 reference oracle 只能基于同一历史预测最佳行动；truth-best-action 可作为离线 ceiling，但不能用来训练/挑选已评价策略。

## E3
在相同物理实验、action menu、测量预算和可比计算预算下比较：
- 完整 B1 reference；
- A 的真实候选集合；
- A + 有限且合法的 full-likelihood correction；
- 仅基础 TSOM 图像的局部设计。

若 A 没有 proposal density，第三组需使用明确定义的生成/采样机制或返回“无法构造该权重”，不能对优化端点直接套 importance sampling。纯单点材料 posterior 的材料 EIG 恒为 0；基于 TSOM 的 Fisher 设计可作为点估计局部 baseline，但不能包装成完整 Bayesian EIG。

预先设 mode-deletion 负控制：删除 reference 中重要的一个结构类，B 在线不知道删除对象；评价方比较错误 acquisition、恢复/未恢复、最终 MAE 与虚假自信。少量 sentinel action 能否识别遗漏是实验问题，不是保证。

将失败分解为四项：A 的候选容量不足、权重错误、更新器失效、acquisition 选错。不能把“候选集中不存在真解”全部归罪于主动测量，也不能用更低 entropy 掩盖漏 mode 后的错误自信。

---

# 15. 独立物理、噪声与成本协议

物理层：设计端与评价端分离。径向首轮沿用 P2 vector-Mie n14，但新频率需重新做 n14/n18 或更高阶与独立 admittance recursion；不能将旧 1.6/2.4 的资格外推至 3.2。一般矢量测试需要独立网格/离散或独立 solver；同一 kernel 换一次 nmax 是离散收敛核对，不等于完全独立 Maxwell 实现。

测量层：每个 action 写清发射模式、复通道数量、驻留、噪声谱、reference target 和切换成本。比较等测量块数之外，必须比较等总驻留/能量或现实金额；不同 q 的 action 不应免费获得更多通道和同样小噪声。

噪声层：固定参考定义，不随 candidate/truth 缩放。共享 gain、频率相关 phase、clock、位置校准与 drift 的 sharing graph 显式保存。独立 thermal repeats 与重复使用缓存数据分别标识。

预算层：full material-frequency evaluation、矩阵分解、RHS、ROM setup/refresh、无效尝试、MC 积分、候选生成和下游决策全部记账。一次 Mie expansion 可服务多个接收角；一次同频分解可服务多个 RHS；计费应忠实反映复用而非“每个角都一整次 full solve”。也不能把评估成本排除后宣布 end-to-end 最快。

部署层：报告 cold single-object 和预先声明的 \(N_{\rm deploy}\) 摊销曲线。若 uniform sweep 很便宜而 B 选测计算很贵，应允许“直接扫频”获胜。若已存在独立制造公差/几何扫描且性价比更高，应选择它，但须保证 context 的条件依赖模型正确；制造 prior 不是一条新独立传感器读数。

---

# 16. 停止测量不是认输：Bellman 与诚实状态

将损失和成本用明确换算系数统一后，有限预算决策可写为
\[
V_t(\pi,b)=\min\left\{
R_L(\pi),\
\min_{a:c(a)\le b}\left[
c(a)+E_ZV_{t+1}(\pi_Z,b-c(a))
\right]\right\}.
\]
不能未经声明把秒、焦耳和 MAE 相加。若无法确定换算系数，提供 Pareto frontier：给定 MAE 的最低成本，或给定成本的最低 MAE。

单步 \(U_L/c\) 是方便的贪心 heuristic，不是一般异质成本顺序问题的全局最优性定理。[R13,R14]

停止依据：
- 风险已经低于下游要求：`STOP_TASK_SUFFICIENT`。
- 声明菜单和规划 horizon 的净收益上界非正：`STOP_NOT_WORTH_COST`。
- 预算内有明确硬对/结构对称性：保留多个世界，`PRESERVE_AMBIGUITY`。
- 输入概率或覆盖不可信：先修支持/先验/校准，不要凭错误 entropy 停止。
- 只有数值排名不稳定或预算用完：如实单列。

**所有单步收益很小不证明两步组合无价值。** 参考校准加结构测量可以互补。至少检查关键二步组合；未检查的 horizon 要写入 stopping scope。“当前菜单下无法稳定区分”不等于一般 Maxwell inverse problem 永远不唯一。

---

# 17. Novelty matrix 与论文排序

详细文献层级和来源在来源账本。以下不把新名字当新理论：

| 已覆盖先例 | 已有内容 | B 仍可能补上的证据 |
|---|---|---|
| Atkinson–Fedorov 1975 [R01] | 非线性 rival-model discrimination、顺序设计 | 有界材料/仪器补偿的 Maxwell 分支与真实 full-wave 风险 |
| Fouda–Teixeira 2014 [R02] | UWB Bayesian inverse scattering、后验协方差引导 adaptive sensing、distorted-Born 扩展 | 多分支复合事件、共享 nuisance、漏候选与任务 MAE |
| Sandhu 等 2023 [R03] | 2D 非线性 EM VIE、传感器 EIG、Laplace | 非单峰材料结构、校准补偿、闭环与完整成本 |
| nonlinear GO-OED [R04] | QoI mutual information、MCMC/KDE 计算 | 任务 loss 与材料 feasible decisions、事件内部积分 |
| nuisance OED / REIG [R05,R06] | nuisance 边缘化、prior ambiguity robust design | 具体 Maxwell 跨频共享补偿冲突及可复算反例 |
| risk/control/decision OED [R07,R08,R17] | risk-aware 与直接下游决策目标 | 强散射多世界的结构机制与成像闭环 |
| offline-online GOOED [R09] | 低秩算子、目标 EIG 误差、免 online PDE | candidate-dependent Maxwell 有用的排序半径和真实回本 |
| BAE OED [R10] | nuisance/model-error-aware covariance 设计 | 不是 BAE Gaussian 近似的完整分支风险及精度核验 |
| IPM BOED 2026 [R15] | surrogate/prior 稳定性、非 KL utilities | 特定物理输出误差到材料决策/排序的可用常数 |
| certified EIG 2026 [R16] | EIG–SNR/Fisher 界（本次仅摘要） | 不能先声称首个 certified utility；需全文对照后界定 |
| sequential testing/OED [R13,R14] | 停止、误判成本、lookahead、异质成本 | 实际 Maxwell nuisance gauge、候选遗漏与不确定性交接 |

方向排序：
**第一：共享补偿冲突 + 材料风险闭环。**
它最接近已有物理见证，也最容易被直接反驳/验证。

**第二：候选不完整时的稳健主动测量与拒绝过度自信。**
重要但需明确支持假设；不要追求无条件开放世界保证。

**第三：反馈算子辅助的精度感知 action ranking。**
有价值但只能在 full-cost 获胜后成为贡献，不能再次以恒等式代替速度。

**不建议作为主线：**新网络预测 action、Ambiguity Graph 本身、把 event bins 当 modes、仅展示第二频率提高 true-event probability、仅减少 observation residual。

---

# 18. Go / No-Go：实验前冻结的建议标准

以下是工程与科研相关性门槛建议，不是由定理推出的普遍数字：

**E0 Go：**至少一个固定任务/校准情景中，代表点/raw separation/Fisher 与完整任务风险确有排名差，并能由 event conditional variation 或 shared compensation conflict 解释；reference 和 noise 合同可复现。若完全没有差异，先缩减算法，不增加网络。

**E1 Go：**在相同总成本下，最终 material MAE 相对最强实用 baseline 至少改善 10%，或在匹配 MAE 下 acquisition cost 至少下降 20%；parent-paired interval 排除零，并且 coverage/false confidence 不显著恶化。低维小样本若 interval 太宽，结论为不确定，不以显著性规则伪装失败机理。若 relative MAE baseline 接近零，则预先使用物理绝对改进阈值替代不稳定比值。

**机制 Go：**观察到预测的“补偿冲突大→有效”“保持因子化/对称性→无额外辨识”“校准动作虽 raw material contrast 小但有效”之一，并有对应负控制；否则没有 Maxwell-specific 主贡献。

**ROM Go：**同一被接受的任务质量下，包含 setup、refresh、finalist full checks 的总时间至少节省 20%，或在相同成本出现更好非支配点。失败则关闭 ROM，不必终止 B。

**强 TAP Go：**至少两个物理上不同的歧义机制，含非球对称矢量场景或适当独立实验；有完整成本前沿、正确共享校准、未解决世界、漏 mode 负控制。只有径向两频正例时，保持受限方法/benchmark 结论，不宣称一般成像突破。

**停止扩张：**两轮冻结预算内没有材料风险/成本优势，或者选测成本大于节省的实测成本；此时优先输出多假设或直接均匀采集。不要通过增加噪声、删除难例、扩大 gain prior、调 event partition 或事后更换 primary metric 来制造 active sensing 胜利。

---

# 19. 本次检查结果与直接可用文件

`code/acquisition_core.py` 提供严格的有限 bank 权重更新、事件质量、离散观测 exact utility、高斯 predictive MC 风险、共享 gain proper quadrature、joint nuisance Schur 和排名区间弃权。

`code/verify_theory.py` 实际执行了 11 项检查，全部通过；Maxwell 调用数为 0。包括：精确 current secant、共享 nuisance 先堆叠再投影、补偿冲突 parallel-sum、proper shared-gain 条件似然、EIG 与分类风险反例、后验质量 martingale、100 次有限概率风险扰动界检查、常数似然不更新、排名区间弃权、高斯理论错误率与 predictive MC、单次 true-event probability 可下降。

其中 complex algebra secant 误差约 \(1.86\times10^{-16}\)；高斯二点 MC 错误率 .158753，对应解析 .158655。这些只是实现核对，不是 Maxwell 性能证据。

本包没有偷偷放入一个可自行执行全波 campaign 的 runner。下一步先按 `CODEX_PROMPT_ZH.md` 做零调用 E0；E1 的 cap 在配置中标记 `execution_authorized=false`。

**最终裁决：推进 B 的物理机制和任务风险设计；让 A 的多解 package 成为输入而不是概率完美性的前置要求。任何没有被预算内测量稳定解决的材料差异，都应继续作为材料世界保留，而不是被软件强制消失。**
