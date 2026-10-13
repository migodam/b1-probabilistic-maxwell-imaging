# A32 — Material-First Probabilistic Maxwell Inverse Scattering
## 有限尺度拓扑辨识、材料后验几何与 Maxwell 模型修正的概率作用

研究日期：2026-10-12。状态：理论研究与小规模数值核验；不是已完成的成像算法或发表优先权声明。

## 结论先行

最值得推进的方向不是“生成十张图”，也不是重新训练一个带 Maxwell penalty 的网络，而是：

> **以有限尺度的材料/拓扑事件为决策对象，建立保留完整假设支持的分层材料后验；用全波有限割线描述拓扑辨识，用双侧材料签名描述模型近似修正究竟改变哪些后验概率。**

建议把主问题写成：**哪些拓扑事件受数据支持，哪些只由先验支持；在有限 Maxwell 预算下，能否准确计算这些事件的后验概率，并辨别值得支付的模型修正？**

本轮得到一个值得继续发展的具体连接：A12/A17 的双侧材料结构不仅约束 Gaussian GN 均值步。在冻结、共同残差的似然模型中，它还给出一个适用于任意共同材料先验的**精确低维似然比签名**。模型改变只通过该签名重新加权后验，签名给定后的条件材料分布不变。拓扑概率是否改变，进一步由“拓扑事件与该签名的条件相关性”决定。这个结论不要求高斯先验，也允许固定的被动性截断和多个材料模式；但它不自动覆盖完整非线性跨拓扑响应，不保证签名维数小，也不保证条件采样便宜。

一般的条件重加权、KL 链式法则、Bayesian evidence、Le Cam 界、LIS、multifidelity MCMC 都是已知技术。本研究可能形成的科学贡献，是 **Maxwell 物理对签名的具体确定、有限尺度拓扑的辨识边界、以及模型修正与材料拓扑决策之间可检验的解耦规律**。它们能否组合成强 TAP 论文，仍取决于新的三维矢量实验与计算收益，不能仅凭本报告中的恒等式宣布成功。

---

## 目录

| 章节 | 内容 |
|---:|---|
| 0 | [原始证据核对与研究路线修订](#sec-0) |
| 1 | [Scientific problem formulation：推断对象和概率空间](#sec-1) |
| 2 | [T1/T3：观测等价、材料几何与先验诚实性](#sec-2) |
| 3 | [T2：全波有限割线与有限尺度拓扑辨识](#sec-3) |
| 4 | [从 A17 到非高斯材料后验：双侧似然比签名](#sec-4) |
| 5 | [T4/T6：forward error、evidence 与材料决策误差](#sec-5) |
| 6 | [T5：物理先验能做什么](#sec-6) |
| 7 | [T9：表示应压缩计算，而非删除未知世界](#sec-7) |
| 8 | [文献调查与原创性地图](#sec-8) |
| 9 | [四条实质不同的算法路线](#sec-9) |
| 10 | [T8：预算分配能证明到什么程度](#sec-10) |
| 11 | [从 posterior 到材料决策与统计校准](#sec-11) |
| 12 | [T10：新测量应区分什么](#sec-12) |
| 13 | [本轮实际完成的小规模验证](#sec-13) |
| 14 | [下一阶段的最小决定性实验](#sec-14) |
| 15 | [T1–T10 的最终裁决表](#sec-15) |
| 16 | [Final research decision](#sec-16) |
| 17 | [原始来源、文献与阅读边界](#sec-17) |

---

<a id="sec-0"></a>

## 0. 原始证据核对与研究路线修订

### 0.1 本轮实际读取范围

完整读取 A31 ZIP 中的理论报告，核对其中数学式与原有 22 项检查，并在独立目录重新运行这 22 项小矩阵/解析/微型偶极检查，结果仍为 22/22。阅读 A17 冻结主稿、补充材料的主要推导、物理约定、成本与限制，并检查相关整页图表。A12 的核心推导由这些冻结稿中的完整证明和研究汇报核对，不声称获得了 A12 所有历史代码。A14/A15、A18–A21、A23–A29 读取了原始入口/裁决，其中 A28/A29 进一步核对主要实验报告。A22 独立原始档案没有取得。A30 按固定 commit 读取 E0/E1 裁决、门槛 JSON、父样本 CSV，并重新聚合四项父样本统计。没有重跑 A30 的优化器、Maxwell solver 或 bootstrap 全部流水线。

因此，这是**针对 A32 关键结论的原始证据审计**，不是对 A12–A31 每个仓库逐行代码的完整复现。资料来源和固定版本列于文末 [S1]–[S13]。

### 0.2 A17 可以保留什么，不能推广什么

冻结稿建立了
\[
J_1-J_0=Y\Gamma^{-1}N,
\qquad
s_1-s_0\in H_0^{-1}\operatorname{Range}[J_0^TY,N^T].
\]

这里所有照明共享同一个实材料变量；复电流因子必须先正确提升/实化；材料、注入、残差、正则项和状态冻结。原 GN 支撑结论要求精确无约束端点步；活跃约束改变时不能原样套用。[S1,S2]

十个完整额外对象的对象级 GN 二次差距中位比 0.4928，表示局部材料步相对完整 GN 的差距下降约 50.7%，不是图像误差下降 50.7%。三个完成的非线性配对材料改善分别为 10.84%、0.22%、7.16%；对应墙钟增加到约 2.679、1.743、1.459 倍。额外对象并非全新盲测。[S1,S2]

补充材料明确写明：输出误差半径没有在实验中实例化。在线 FP64 容差是数值显著性检查，不是严格 Maxwell 输出误差证书。其数据使用 clean-data 范数归一化，不是经过实验标定的噪声白化；既有 information-volume 不能直接解释为 mutual information。[S2]

这些结果支持“指定模型改动怎样影响材料”的理论，不支持“任何 current 修正都改善材料后验”的主张。

### 0.3 A30 E0/E1 的最新事实

固定版本：`717d3d26a290248bd2b9aeb55f3d03694b407b84`。真实后端是 24×24、8 源、32 接收器的 **二维 scalar TM**；1152 实材料坐标、512 实观测。不是三维 vector Maxwell。[S3]

下表重新使用公开 `PARENT_STATISTICS.csv` 的 12 个复杂父样本 gain 做等权平均；是原始评分的再聚合，不是新增恢复实验。

| 指标 | 本轮重新聚合 | 原结果的意义 |
|---|---:|---|
| 完整候选凸包 offline oracle | +13.1225% | 凸组合具有容量，不是部署效果 |
| 实际 full-L weighted mean | −3.3154% | 已实施决策没有兑现同 bank 容量 |
| Sparse-64 offline oracle | +35.6896% | 用真值拟合的表示潜力；未实现该在线重建器 |
| TRAIN-PCA-64 offline oracle | +21.9074% | 同样是表示潜力，不是后验质量 |

Sparse-64 原报告父样本区间约 [29.50%,42.81%]；full hull 约 [11.51%,14.65%]。这些都是已暴露 development 样本上的探索性区间。它们没有证明真实拓扑具有足够 posterior support，也没有证明新噪声、新形状下能找到这些表示。[S3]

E1 的 32 条件中有 29 个合法配对，3 个无注册 cavity edit。复杂父样本达到预定“原距离至少 3、补偿后距离不超过 1”的个数是 **0/12**，门槛需要 6/12；最小合法全波补偿距离约 8.84，复杂组最小约 9.18。局部补偿明显减小数据距离，但没有通过该实验的歧义门槛。优化有封顶/未收敛情况，取得的距离是可行解值，不是全局最小距离下界。因此既不能说补偿歧义已找到，也不能说它不存在。[S3]

uniform Neumann/Taylor 界没有给出有用分类区间。`q>=1` 是保守充分条件失败，不是证明 Maxwell 奇异或材料不可区分。E2 未运行。当前 bounded-compensation 和 certificate-guided F 路线的实验裁决为 NO-GO。[S3]

**必须修订 A31：**附件 A31 的 §0 仍把 A30 hull 容量实验写成未完成；该状态已经被以上 E0/E1 取代。A31 的代数证明可以保留，其“耗散认证作为首要部署路线”的排序不能不经新门槛而直接继承。

### 0.4 其他负结果怎样约束 A32

A28：复杂材料六步 F 相对 N3 仅约 +1.01%，相对 L 约 +0.40%，区间跨零；数据残差却下降约 27.33%。指定结构反事实经常可分，不等于整个材料类可唯一恢复。其历史预算至少 807/768 的缺陷应保留，不把程序测试免费化。[S4]

A29：hard candidate oracle 仅约 +2.1468%，但不能据此排除 convex averaging 或连续材料表示；A30 已给出这一反例。[S5]

A18-B/A19：结构或代理精度通过，并没有使整体成本、安全性或实际效益通过。A23：一次二阶反馈的 3D 材料图可行性、质量、端到端成本均失败。A24–A27 的二维结果显示结构先验确有价值，但显式反馈并未稳定优于匹配基线；普通 full-adjoint+同 TV 和所谓 exact-F TV 是同一物理算法，不能制造一个虚假的被击败对照。[S6–S10]

A14 的参考后验与 Flow 未运行；A15 的 33³/41³ 参考积分比较只有 94/384 查询通过，当前两模态解释未达到噪声精度。其后验风险和采集排序是 provisional。**A32 首先要解决可信参考后验，不应再以未充分收敛的 sampler 给另一个 sampler 当真值。**[S11,S12]

---

<a id="sec-1"></a>

## 1. Scientific problem formulation：推断对象和概率空间

### 1.1 材料与实验

用 \(d\) 表示实验设计，包括频率、源、接收器与极化。材料为函数或张量场 \(\chi\)，实验 nuisance 为 \(\gamma\)。在实化后的测量空间中，
\[
Y=F_d(\chi,\gamma)+\eta,\qquad \eta\sim N(0,\Sigma_d).
\]

若噪声 proper complex Gaussian 的复协方差为 \(\sigma_c^2I\)，实部/虚部方差各为 \(\sigma_c^2/2\)；实白化是 \(\sqrt2/\sigma_c\)，不能漏掉该因子。若 covariance 随材料变化，likelihood 必须包含 \(\tfrac12\log\det\Sigma(\chi)\)；以每个候选自己的 signal norm 重新归一化会改变概率模型。

材料参数空间优先定义为具有明确本构约束的函数空间或其一致离散化。对各向同性介电常数，可以用 \(\Re\epsilon,\Im\epsilon\)；各向异性材料用矩阵场；多频推断可用共享的 Lorentz/Drude 系数。无限维问题使用先验测度 \(\mu_0\)，不写不存在的“无限维 Lebesgue 均匀密度”。[L6,L7]

### 1.2 连续材料后验与层级拓扑并不冲突

令
\[
\mathcal X=\bigsqcup_{H\in\mathcal H}\{H\}\times\Theta_H\times\Gamma,
\qquad \chi=\mathcal D_H(\theta).
\]

其中 \(H\) 可以包含对象数、封闭空腔数、连通关系；\(\theta\) 描述边界、分区材料与细尺度场。给定 proper prior
\[
\mu_0(dH,d\theta,d\gamma)
=\pi_H\,\mu_H(d\theta,d\gamma),
\]
后验为
\[
\frac{d\mu^y}{d\mu_0}(x)=\frac{e^{-\Phi_y(x)}}{Z(y)},
\quad
\Phi_y=\tfrac12\|\Sigma^{-1/2}(F_d(x)-y)\|^2.
\tag{1}
\]

材料后验是 \(\mathcal D\) 的 pushforward；拓扑后验是事件的概率。不必把所有拓扑硬塞进同一个光滑 manifold。不同维数 \(\Theta_H\) 可比较，因为各分支有自己的规范化测度和离散质量 \(\pi_H\)。但不能比较未经规范化的惩罚值或跨维 density peak。

多个 chart 表示同一物理材料时，需要不重叠分区、partition-of-unity 权重或明确的 latent prior；否则一个形状因有更多参数化方式而被重复加权。多个同类组件的标签交换需要规范化或 \(K!\) 对称因子。采用一个“other”标签而不给它可采样/可积分的材料分布，不构成真实的开放支持。

### 1.3 六个对象不能混为一谈

点估计是一个行动；局部极小值是优化景观的性质；拓扑模式是材料事件/分支；Bayesian posterior 是指定生成模型下的条件概率测度；不确定性集合可以是可信集或频率学置信集，两者含义不同；观测等价类则由相同测量分布定义。

多个局部极小值不自动对应多个有显著质量的后验模式：窄峰可以有很高 density 而质量很低。一个模式也可能包含多个拓扑。MAP 随 density 的参考坐标改变，事件概率在正确推送测度后不改变。

当主要决策是“有无内部空腔”，且相互竞争的拓扑都能解释数据时，单图 MSE 估计不是足够的终点。后验均值仍是平方损失下的正确 Bayes action，但通常不是某个有代表性的物理世界。

---

### 1.4 拓扑 evidence：为什么最小 residual 不能给概率

对两个拓扑分支，严格的比较对象是
\[
Z_H(y)=\int p(y\mid\theta,\gamma,H)\,\mu_H(d\theta,d\gamma),
\qquad
\frac{P(H_1\mid y)}{P(H_0\mid y)}
=\frac{P(H_1)}{P(H_0)}\frac{Z_{H_1}(y)}{Z_{H_0}(y)}.
\tag{1a}
\]

取最小 residual 只保留积分中的最高峰，丢掉了峰宽、可行体积、prior density、其他局部峰及 nuisance 的积分。当两个候选类都能拟合噪声时，拟合能力更强的分支不应因此自动获得更大 posterior probability。

**最小反例。**观测 y=0、噪声标准差 σ。H0 固定预测0；H1 的预测 θ 在 [−A,A] 上有 proper uniform prior。两类最小 residual 都为0，但
\[
\frac{Z_1}{Z_0}
=\frac{\sqrt{2\pi}\sigma}{2A}
\{2\Phi(A/\sigma)-1\}.
\tag{1b}
\]
取 A=2、σ=0.1，Bayes factor 约0.06267，等分支 prior 下 P(H1|y) 约0.05897。它不是对复杂模型的任意惩罚，而是 H1 大部分 prior 预测并不支持观测。随包脚本直接计算了该例；(1a) 和 Occam-volume 原理属于经典 Bayesian model selection，不作新颖性。[L1]

若分支负 log posterior 在内部非退化极小点 θ̂ 附近可二阶展开，且该峰占据主要质量，则 Laplace 近似为
\[
Z_H\approx p(y\mid\widehat\theta,H)p(\widehat\theta\mid H)
(2\pi)^{d_H/2}|\nabla^2[-\log(p(y\mid\theta,H)p(\theta\mid H))]_{\widehat\theta}|^{-1/2}.
\tag{1c}
\]

这里 θ 可同时包含 nuisance；proper prior 的归一化不能省略。多峰应相加；边界峰、退化 Hessian、拓扑出生点和弱可辨参数可能使 (1c) 不适用。Laplace 不是完整 posterior 的默认替代。

### 1.5 prior sensitivity 不可能完全消除

若 Z1=Z0，则后验 odds 等于 prior odds；弱数据时不存在一个同时“完全不依赖 prior”又能给唯一 Bayesian topology probability 的规则。应报告 Bayes factor、明确 prior odds，以及在一个事先声明的合理 prior 家族下的 probability range，而不是调整 prior 直到出现确定判断。

若只改变分支 prior odds、保持分支条件 prior 不变，后验 odds 与它成正比；若连材料/几何 conditional prior 也改变，ZH 本身也改变。即使两个假设下最优拟合相同，Occam factor 仍可能对 prior 尺度敏感。数据不足应表达为 likelihood 不能主导比较，而不是宣布概率模型失去意义。

---

<a id="sec-2"></a>

## 2. T1/T3：观测等价、材料几何与先验诚实性

### 2.1 精确观测等价

对完整参数 \(x=(\chi,\gamma)\)，定义
\[
x\sim_d x'\quad\Longleftrightarrow\quad
\mathcal L(Y\mid x,d)=\mathcal L(Y\mid x',d).
\tag{2}
\]
相同正定 Gaussian covariance 时，这等价于 \(F_d(x)=F_d(x')\)。它是等价关系；“距离小于 \(\varepsilon\)”一般不满足传递性，不能直接称 quotient。

若先边缘化 nuisance，\(K_\chi=\int K_{\chi,\gamma}\mu_0(d\gamma\mid\chi)\)，则 \(K_\chi=K_{\chi'}\) 同时依赖实验和所用 nuisance prior，不再是只由 Maxwell 决定的几何。

设有限维 \(F\) 在邻域内具有常秩，则其局部等值纤维的切空间为 \(\ker DF\)。秩变化、自交与跨拓扑粘合处不保证形成光滑流形。局部核也不能识别全局碰撞，例如 \(F(x)=x^2\) 的正负分支。

### 2.2 条件纤维定理：坐标不变的 prior–likelihood 分解

**命题 1。**在标准 Borel 参数空间上，若 likelihood 仅依赖 \(z=F(x)\)，则在 regular conditional distributions 意义下，
\[
\mu^y(dx\mid F(x)=z)=\mu_0(dx\mid F(x)=z).
\tag{3}
\]

**证明。**Bayes 密度 \(e^{-\Phi_y(F(x))}\) 在同一纤维上为常数，在该纤维的条件归一化中约掉。等式为几乎处处意义，不要求 \(\{F=z\}\) 有正先验质量。□

这是一般 Bayes 条件化的直接结论，不是新定理。它精确表达：测量重新分配不同可观测响应的质量，但不会在完全相同响应的材料之间创造区分信息。

不能把 (3) 误写成“弱观测材料的边缘后验等于边缘先验”。如果先验把可观测变量 \(u\) 与不可观测变量 \(v\) 相关联，测到 \(u\) 后 \(v\) 的边缘也会改变；新增内容来自先验相关，而不是对 \(v\) 的独立物理观测。随包矩阵例中，\(P(v=1)\) 从 0.5 变为 0.82，但 \(P(v\mid u)\) 不变。

### 2.3 自然材料度量是什么

实白化 Jacobian \(J=\Sigma^{-1/2}D_\chi F\) 给出半正定 pullback metric
\[
g_\chi(h,k)=\langle Jh,Jk\rangle,
\qquad I_\chi=J^TJ.
\tag{4}
\]

材料坐标 \(\chi=\mathcal D(\theta)\) 改变后，\(I_\theta=D\mathcal D^T I_\chi D\mathcal D\)。度量对象具有坐标协变性，但裸矩阵的 eigenvalues/rank threshold 会受单位和坐标尺度影响。应与物理质量度量、声明的先验尺度或明确任务度量构成广义特征问题。A30 的 mass-normalized material coordinates 应在这里保留，而不是改用 current SVD 排序。

Gaussian 线性问题可写
\[
C_y=(C_0^{-1}+J^TJ)^{-1}.
\]
这只在指定 prior/noise/线性模型下成立；高频近共振、多模式和活跃被动性边界不能用一个 Hessian inverse 代表完整后验。[L7,L8]

### 2.4 nuisance 的局部影响

若实白化模型为 \(y=f(\chi)+B\gamma+\eta\)，\(\gamma\sim N(0,\Gamma)\)，\(\eta\sim N(0,I)\)，且 B 固定，则边缘精度
\[
P=(I+B\Gamma B^T)^{-1},\qquad I_{\chi,\mathrm{eff}}=J^TPJ.
\tag{5}
\]
无先验自由 nuisance 的 profile 极限可出现 \(J^T(I-P_B)J\)，但 bounded nuisance 和 informative nuisance prior 不等价于正交投影。若 B 或 covariance 随 \(\chi\) 变化，必须保留 log-determinant 和 covariance derivative 的信息。

噪声、nuisance、结构先验、数值模型误差和非线性混叠通常有交互，不能唯一地拆成几项可加的“百分比”。条件方差分解可以相对于一个明确层级成立，但其分解次序和条件变量必须声明。

---

<a id="sec-3"></a>

## 3. T2：全波有限割线与有限尺度拓扑辨识

### 3.1 不依赖 Born 的精确材料割线

固定几何、源与接收器，材料源方程为
\[
j_i=X_i(e^{inc}+Gj_i),\qquad L_i=I-X_iG.
\]
当两个端点可逆时，直接相减得
\[
\boxed{j_1-j_0=L_1^{-1}(X_1-X_0)E_0^{tot}},
\qquad
\boxed{F_1-F_0=S L_1^{-1}\Delta X E_0^{tot}}.
\tag{6}
\]

证明只需将 \(j_0\) 代入 \(L_1\)：\(L_1(j_1-j_0)=\Delta X E_0^{tot}\)。因此 (6) 包含完整多重散射，不要求 \(\Delta X\) 幅度小。DDA 中 X 是 polarizability，不是裸 \(\chi\)；必须使用实际 \(\Delta a\)，切向极限才使用 \(a'(\chi)\)。[S1,S2,S13]

这比单点 Jacobian 更适合讨论“把某个区域挖空”。shape derivative 描述既有边界的光滑运动；topological derivative 描述小区域出生的渐近效应；有限割线描述两个具体拓扑端点的真实响应。三者不能互相替代。小孔洞可以是有限幅度、微小体积扰动，而不是在所有材料范数中都小的 tangent。

### 3.2 固定两材料的可辨识性

在相同 Gaussian covariance、等先验下，定义
\[
d(\chi_0,\chi_1)=\|\Sigma^{-1/2}(F_1-F_0)\|.
\]
则
\[
D_{KL}(P_0\|P_1)=d^2/2,
\quad P_e^*=\Phi_{N(0,1)}(-d/2),
\quad TV(P_0,P_1)=2\Phi_{N(0,1)}(d/2)-1.
\tag{7}
\]

证明：Gaussian likelihood ratio 是沿 \(\Sigma^{-1}(F_1-F_0)\) 的一维统计量，在两个假设下均值为 \(\pm d^2/2\)、方差为 \(d^2\)。阈值零给出 (7)。这属于经典二点检验，不是 A32 原创。其 full-wave 内容在于以 (6) 而非 Born 响应代入。

对复合拓扑，定义 \(P_H^Y=\int P_\theta^Y\mu_H(d\theta)\)。最优错误仍由这两个**混合测量分布**的 TV 决定，不能把每类 MAP 的距离或某个 \(d_{comp}\) 直接代入 (7)。存在一对很近的跨拓扑材料，足以构造 minimax 困难实例；它不表示两类在给定先验下大多数材料都难区分。

对任意两分支参数的 coupling \(\lambda\)，相同 covariance 下的联合凸性给出
\[
D_{KL}(P_0^Y\|P_1^Y)
\le \int \tfrac12 d(\theta_0,\theta_1)^2\,d\lambda.
\tag{8}
\]
若能用一个覆盖先验质量的 coupling 让右侧小，则是复合类难区分的充分条件。仅最小距离小不够。反过来，仅有限优化得到较大补偿距离也不提供整个类的分离下界；这正是 A30 E1 的证据边界。

### 3.3 “任意小孔洞都可检测”不可成立

**命题 2：有限尺度 no-go。**设 \(\chi_a\) 是从 \(\chi_0\) 的内部移除半径 a 的有限对比材料区域得到的物体；噪声 covariance 固定正定。假设相关全波状态和目标伴随场在该区域及 \(a\downarrow0\) 下具有统一界，使
\[
\|\Sigma^{-1/2}S L_a^{-1}\Delta X_a E_0^{tot}\|\le C a^d.
\tag{9}
\]
这里 d=3 为三维体积缩放，二维无限圆柱截面为 d=2。则
\[
P_e^*(\chi_0,\chi_a)\ge \Phi(-Ca^d/2)\longrightarrow 1/2.
\tag{10}
\]

**证明。**将 (9) 代入 (7)。若 \(\Sigma=\sigma^2I\)，可把噪声尺度显式写成 \(C_0a^d/\sigma\)。对单位目标向量 q，(6) 的伴随形式是局部积分 \(\langle z_{a,q},\Delta X_a E_0^{tot}\rangle\)，有界场与有限对比给出常数乘体积。统一场界是条件，而不是由“被动”二字自动推出。□

因此无最小孔洞尺度的材料类不能有统一可靠拓扑判别。必要的可辨尺度受 \(\sigma\)、频率、孔径、极化、照明场和背景反馈共同决定；不能把 (10) 写成一个普适的 \(\lambda/2\) 分辨率公式。

近共振或弱耗散会改变 C，不能一边假定统一稳定一边宣称覆盖任意共振。共振可能提高特定结构敏感度，也可能同时放大 nuisance、模型误差和积分困难。

### 3.4 拓扑事件必须带物理定义

建议主事件是
\[
H_\rho(\chi)=\mathbf 1\{\text{存在尺度至少 }\rho\text{ 的内部封闭空腔}\},
\tag{11}
\]
并明确“空腔”的材料阈值、最小厚度/直径、与外部是否连通。二维内部洞通常与 \(b_1\) 有关；三维封闭空腔通常与 \(b_2\) 有关，三维 \(b_1\) 描述隧道/handle。不能把两者混写。

如果半径采用连续 prior 且只取 \(a>0\)，那么“任何孔洞存在”的先验概率就是 1；要给无孔非零概率，需要 \(P(K=0)>0\) 的离散原子，或把事件定义成 \(a\ge\rho\)。这不是数值小细节，而是模型是否能回答用户问题的基本条件。

---

<a id="sec-4"></a>

## 4. 从 A17 到非高斯材料后验：双侧似然比签名

本节是本轮最值得继续推进的理论连接。它不是把 current reconstruction 作为推断目标，而是利用一个已经存在的 Maxwell 近似差异，确定其在**材料后验**中的作用。

### 4.1 精确陈述

令 s 是共同的实材料坐标，\(\mu_\star(ds)\) 是任意 proper prior，包括非高斯、混合、固定材料可行域截断。两种冻结的似然模型为
\[
L_i(s)=\exp[-\tfrac12\|r+J_i s\|^2],\quad i=0,1,
\]
并满足
\[
J_1=J_0+Y\Gamma^{-1}N.
\tag{12}
\]
这里 r 已真正白化，\(\Gamma\) 可逆；两模型使用相同 r、材料映射和 prior。假设证据有限且非零。

定义材料签名
\[
U(s)=(a,b)=(Ns,\;Y^TJ_0s),
\]
及
\[
h(a,b)=(Y^Tr+b)^T\Gamma^{-1}a
+\tfrac12 a^T\Gamma^{-T}Y^TY\Gamma^{-1}a.
\tag{13}
\]

**命题 3：双侧后验改变定理。**两后验满足
\[
\boxed{\frac{d\pi_1}{d\pi_0}(s)
=\frac{e^{-h(U(s))}}{\mathbb E_{\pi_0}e^{-h(U)}}.}
\tag{14}
\]
并且
\[
\boxed{\pi_1(ds\mid U=u)=\pi_0(ds\mid U=u).}
\tag{15}
\]

如果 KL 有限，则
\[
D_{KL}(\pi_1\|\pi_0)=D_{KL}(\pi_1^U\|\pi_0^U),
\]
而总变差也满足
\[
TV(\pi_1,\pi_0)=TV(\pi_1^U,\pi_0^U).
\tag{16}
\]

**证明。**设 \(c=\Gamma^{-1}Ns\)。展开共同 residual：
\[
\tfrac12\|r+J_0s+Yc\|^2-\tfrac12\|r+J_0s\|^2
=(Y^Tr+Y^TJ_0s)^Tc+\tfrac12 c^TY^TYc=h(U(s)).
\]
两个 Bayes 密度相除即得 (14)。在 U 的纤维上重加权因子恒定，得到 (15)。对 (14) 的 log 和绝对偏离按 U 条件积分，分别给出 KL 和 TV 等式。□

### 4.2 这比原 GN 支撑多出了什么

原 A17 结论是“Gaussian 二次优化的均值步差位于哪个线性空间”。(14)–(16) 则表明：**整个概率测度的改变由哪些材料统计量决定**，不要求材料先验高斯，不要求 posterior 单模态，不需要计算 MAP。

第一组 \(Ns\) 是材料激发签名；第二组 \(Y^TJ_0s\) 是新增接收响应与已有材料响应的干涉签名。只保留其中一组通常不能得到正确 likelihood ratio。这个双侧性正是 A17 Maxwell 因子在概率问题中的具体内容，而不是任意挑选一个 latent。

固定被动性约束可通过先验支持直接纳入 (14)。即使 truncated posterior 的边界很活跃，条件重加权恒等式仍成立；这与原“无约束 GN 线性支撑”需要区分。对非高斯先验，不能进而声称 posterior mean 的差仍位于同一 Gaussian 线性支撑。

### 4.3 Gaussian 特例与 A12/A17 严格连续

若 \(\pi_0=N(m_0,C_0)\)，令 \(Q=[J_0^TY,N^T]\)，则条件 Gaussian 的均值沿 \(C_0Q\) 变化，从而
\[
m_1-m_0\in\operatorname{Range}(C_0Q),
\]
\[
\operatorname{Range}(C_1-C_0)\subseteq\operatorname{Range}(C_0Q).
\tag{17}
\]
共同 Gaussian prior 下 \(C_0=H_0^{-1}\)，恢复原 paired material support。signed Hessian defect 的 rank 至多 2r，并不意味着它半正定。改变同一观测的近似模型不是新增独立测量。[S1,S2,S13]

该特例中的 Woodbury 与 Gaussian 低秩更新已在 A31 和既有文献中存在，不是 A32 新颖性。更值得研究的是非高斯、拓扑事件及条件误差的部分。

**非高斯反例。**令共同 prior 在 s=(0,0),(1,1),(−1,1) 上等质量，J0=0、J1=[1,0]、r=0。签名仅含第一坐标，但旧均值为 (0,2/3)，新均值为 (0,2e^(−1/2)/(1+2e^(−1/2)))。均值差沿第二坐标，第一坐标的均值始终为0；旧 covariance 的交叉项也为0。因此不能把 Gaussian 的 C0Q 均值支撑推广给任意非高斯 prior。另一方面，给定第一坐标后的第二坐标依旧确定为其平方，所以 (15) 完全成立。

### 4.4 拓扑概率是否改变：需要第二个条件

令 \(H\) 为任意可测材料事件，\(p=\pi_0(H)\in(0,1)\)，
\[
m_H(u)=\pi_0(H\mid U=u),\qquad w(u)=e^{-h(u)}.
\]
由 (14)，
\[
\boxed{\pi_1(H)-\pi_0(H)
=\frac{\operatorname{Cov}_{\pi_0}(m_H(U),w(U))}
{\mathbb E_{\pi_0}w(U)}.}
\tag{18}
\]
定义
\[
\tau_H(U)=\frac{\operatorname{Var}_{\pi_0}(m_H(U))}{p(1-p)}\in[0,1].
\]
在 \(\mathbb E w^2<\infty\) 时，Cauchy–Schwarz 给出
\[
\boxed{|\pi_1(H)-\pi_0(H)|
\le\sqrt{p(1-p)\tau_H(U)}\,
\frac{\operatorname{sd}(w)}{\mathbb Ew}.}
\tag{19}
\]

**解释。**一次物理修正即使显著改变材料后验，只要改变的签名与空腔事件条件独立，空腔概率就完全不变。如果相关性弱，拓扑概率变化也受到抑制。这样就把“forward 修正大但 material decision 变化小”变成一个可证伪的机制，而不是只重复观察。

随包非高斯离散反例中，整个材料后验 TV 改变约 0.360818，但 \(P(H)\) 精确维持 0.3。另一个相关例中，\(P(H)\) 从 0.6225 变成 0.871794，实际变化 0.249294，小于 (19) 的 0.253221。它们核验代数，不是证明真实孔洞通常与 Maxwell 签名无关。

### 4.5 必须写进论文的限制

(14) 描述的是**两种后验之间的改变**，不是说完整后验本身低维。U 给定后的所有材料不确定性仍然存在，不能把它们设为零。

多个照明的真实实化签名维数可能远大于电流列数，接近完整材料维度。A12 的少量 current 列产生近满维材料修复库的成本负结果仍然有效。

如果两个模型不共享 forward anchor，额外偏移 a 会产生 \(s^TJ_0^Ta\) 等项。忽略这些项会错误地改变 topology odds；共同 field offset 一般不能约掉。

若每个 topology branch 使用自己的 \(J_H,U_H\)，分支内条件恒等式并不自动给出分支间证据。必须保留 H、各分支 offset 及其 normalization。把 H 加入 U 会使 H 可由签名直接确定，此时 \(\tau_H=1\)；不能再宣称“拓扑与签名独立”。

条件采样 \(\pi_0(ds\mid U=u)\) 对任意非高斯复杂 prior 可能很难。因此本节没有免费给出高效生成器，也没有证明在线反馈值得支付。

### 4.6 非线性余项与可证明边界

设在同一材料邻域内
\[
F_i(s)=f_\star+J_i s+\rho_i(s),
\quad \|\rho_i(s)\|\le\tfrac12 M_i\|s\|^2.
\]
令 \(\delta(s)=\Phi_1(s)-\Phi_0(s)\)，则
\[
\delta(s)=h(U(s))+\epsilon(s),
\]
\[
|\epsilon(s)|\le
\tfrac12\{M_1\|r+J_1s\|+M_0\|r+J_0s\|\}\|s\|^2
+\tfrac18(M_1^2+M_0^2)\|s\|^4.
\tag{20}
\]

证明是分别展开 residual 与 remainder 的二次项。若 \(\operatorname{osc}\epsilon\le\omega\) 在完整支持上成立，则真实 \(\pi_1\) 与仅用 (14) 的签名后验相差
\[
TV\le\tanh(\omega/4).
\tag{21}
\]
若只在 \(\Omega\) 成立，还必须加上真实/近似后验在 \(\Omega^c\) 的尾部质量；有限样本上误差小不能代替这两个尾部质量上界。

A30 的 uniform bound 失败说明不能假定 (20) 在整个材料 box 足够小。大孔洞、拓扑出生、弱观测造成的宽后验都可能超出局部范围。因此应把 (14) 首先用作**诊断与 proposal/preconditioning 结构**，而不是未经核验地替换完整 likelihood。

### 4.7 与已有理论的重合和真正待解决部分

若似然比只依赖 U，条件分布不变，这是一般统计事实；KL 链式法则也是已知结果。Zahm 等的条件期望后验降维、Spantini 等的低维 transport couplings 已覆盖一般的降维和条件运输思想。[L7,L8]

A32 不应把这些一般事实包装成新数学。可争取的贡献是：(a) Maxwell feedback 的双侧材料因子恰好确定具体 likelihood-ratio 统计量；(b) 它与拓扑事件之间的耦合决定材料概率变化；(c) 在真实矢量 Maxwell 中得到便宜、非空的非线性条件误差控制。前两项本报告已给出推导，第三项没有闭合，必须作为下一阶段的核心研究门槛。

---

<a id="sec-5"></a>

## 5. T4/T6：forward error、evidence 与材料决策误差

### 5.1 从数值模型误差到 likelihood

使用真正白化的量，令 \(e=F-\widetilde F\)、\(\widetilde r=\widetilde F-y\)。
\[
\delta=\Phi_F-\Phi_{\widetilde F}
=\widetilde r^Te+\tfrac12\|e\|^2.
\tag{22}
\]
若 \(\|e\|\le b\)，则
\[
|\delta|\le\|\widetilde r\|b+\tfrac12b^2.
\]
同先验下
\[
\pi_F\propto e^{-\delta}\widetilde\pi.
\tag{23}
\]

若 \(a\le\delta\le b\) 全支持成立，\(Z_F/Z_{\widetilde F}\in[e^{-b},e^{-a}]\)。设 \(\omega=\operatorname{ess\,osc}\delta\)，则
\[
TV(\pi_F,\widetilde\pi)\le\tanh(\omega/4).
\tag{24}
\]
它对参数坐标不变，且是 sharp 的 bounded-tilt 界。A31 已给出该结果；本报告不重复申报新颖性。一般 forward-to-posterior stability 也已有系统理论。[S13,L9]

“\(J_F-J_A\) 很小”不能推出 (24) 很小：取 \(F_A=F+c\)，Jacobian 完全相同，但 likelihood 与 evidence 可以严重错误。

### 5.2 为什么绝对 field error 不是 topology odds 的正确目标

对分支 H，
\[
\frac{Z_{F,H}}{Z_{A,H}}
=\mathbb E_{\pi_A(\cdot\mid H,y)}e^{-\delta}=:M_H.
\]
二分支有
\[
\boxed{\operatorname{logit}p_F(H)-\operatorname{logit}p_A(H)
=\log M_H-\log M_{\bar H}.}
\tag{25}
\]

若 \(\delta\) 是所有材料共有的常数，posterior 不变；绝对 evidence 会乘共同常数，Bayes factor 不变。如果 \(\delta\) 只在每个分支内近似常数，但分支之间不同，拓扑概率仍可强烈改变。

沿桥接 \(\pi_t\propto e^{-t\delta}\pi_A\)，
\[
\frac{d}{dt}\operatorname{logit}\pi_t(H)
=-\{\mathbb E_{\pi_t}[\delta\mid H]-\mathbb E_{\pi_t}[\delta\mid\bar H]\}.
\tag{26}
\]

因此要保护 topology odds，应控制**跨分支的 likelihood error 对比**，而不是只最小化平均 current error。若每分支有统一界 \(a_H\le\delta\le b_H\)，log-odds 修正位于
\[
[a_{\bar H}-b_H,\ b_{\bar H}-a_H].
\]
单个 MAP 点上的误差界不能代替整个分支的 evidence 积分界。

### 5.3 误差到行动

任何材料事件 A 满足 \(|\pi_F(A)-\pi_A(A)|\le TV\)。对值域为 \([0,L_{max}]\) 的损失，使用近似后验的最优行动 \(a_A\) 有
\[
R_F(a_A)-\inf_aR_F(a)\le2L_{max}TV(\pi_F,\pi_A).
\tag{27}
\]
证明将两后验在行动 \(a_A,a_F\) 上的风险差分别加减，利用 bounded-expectation 的 TV 界。

对无约束平方损失，后验均值为 Bayes action，excess risk 精确为两个均值之差的平方；材料支持直径有界时可进一步控制。若行动受集合 A 约束而损失仍为固定正定度量下的平方损失，则只要全局投影存在，后验均值在 A 上的全局最近点仍是 Bayes action，即使 A 非凸；非凸性使投影可能非唯一或难计算。任意阈值化、局部投影或换成拓扑检测损失后，不能沿用这个最优性。

### 5.4 measurement noise 与 prior 改变

若测量均值集的白化直径为 \(D_F\)，y 改为 \(y'\)，势能差的 oscillation 至多
\[
D_F\|\Sigma^{-1/2}(y-y')\|.
\]
于是 (24) 给出所有材料事件概率的稳定界。若均值集无界，需要已知 Bayesian well-posedness 理论中的可积性条件，不能直接使用有限直径常数。[L9]

先验变化若在共同支持上具有有界 log-density-ratio oscillation，可以与 \(\delta\) 合并。若真实拓扑在采用先验中概率零，这种稳定界不可用，任何精确 Maxwell correction 都无法创造已被排除的后验质量。

### 5.5 Maxwell 残差与耗散证书的正确位置

对 \(Lj=b\)，近似解 \(\widetilde j\) 的 residual \(\rho=b-L\widetilde j\) 给出
\[
F-\widetilde F=S L^{-1}\rho.
\tag{28}
\]
已知稳定常数可把它变成输出界；但 \(\|L^{-1}\|\) 或合适的 coercivity 常数不是免费数据。

在适用的非零响应支持上，用阻抗形式 \(Z=X^{-1}-G\)。固定 \(e^{-i\omega t}\)，若实际离散模型满足
\[
D=-\operatorname{Im}_H Z\succeq D_0\succ0,
\]
则
\[
\|D_0^{1/2}Z^{-1}D_0^{1/2}\|\le1.
\]
对 \(Zj=b\)，\(Z^*\psi=S^*q\)，近似 primal/dual 的残差 \(\rho,\tau\) 满足
\[
|q^*S(j-\widetilde j)-\widetilde\psi^*\rho|
\le\|\tau\|_{D_0^{-1}}\|\rho\|_{D_0^{-1}}.
\tag{29}
\]

这是 A31 中已有的 dual-weighted residual/耗散分析，不是新的 Bayesian 技术。harmonic Maxwell 的 certified reduced-basis output bounds 至少在 2010 年已有专门研究。[S13,L10]

合法性要由实际 polarizability、自项、辐射边界和质量度量确认；辐射反作用不能计两次。零 contrast 区域不宜强行使用 \(X^{-1}\)；可退回 L 形式。弱耗散可让 D0 下界趋零，近共振与非正规性使普通范数界松；谱半径小也不足以控制非正规 resolvent。全域 Neumann 界失败时，不能加虚假 loss floor 获得“认证”。

**算法决策：**A32 主流程不应依赖 (29) 一定有用。先保留准确 full-likelihood fallback；只有 certificate 构造与 dual 成本确实小于节省的 full calls，且界不空泛，才启用在线证书。无法认证时返回未认证，而不是将经验 residual 当误差半径。

---

<a id="sec-6"></a>

## 6. T5：物理先验能做什么

### 6.1 材料约束与 forward 恒等式必须分开

采用 \(e^{-i\omega t}\)  convention，局部各向异性被动介质要求损耗矩阵 \(\operatorname{Im}_H\epsilon\succeq0\)；互易常用 \(\epsilon^T=\epsilon\) 等本构对称性，而不是 Hermitian 性。被动不等于互易，被动也不等于无损。

一旦每个候选由正确 Maxwell solver 精确求解，Poynting balance 是 forward 模型已满足的关系。再把同一 balance 加为独立 prior/likelihood，会重复计算物理信息。若额外测到总吸收功率、热响应或能量通量，那是新的观测，可增加 likelihood，而不是“免费能量先验”。

材料本构限制、可制造几何、已知对称性和材料类别是真正的 prior 信息。它们能排除不物理候选，但不能替代有限测量中缺失的结构信息。先验收缩不等于从测量中取得新证据。

### 6.2 因果、色散与多频

可构造共享 Lorentz 型参数化，例如
\[
\epsilon(\omega)=\epsilon_\infty+
\sum_\ell\frac{A_\ell}{\omega_\ell^2-\omega^2-i\gamma_\ell\omega}
+\frac{i\sigma}{\omega},
\]
\(A_\ell\succeq0,\gamma_\ell>0,\sigma\succeq0\)，再给这些物理系数 proper prior。实际单位中的 \(\epsilon_0\) 因子应一致恢复。这比各频率完全独立猜测复介电常数更符合一个声明的材料族。[L11]

Kramers–Kronig 关系关联全频响应；只测有限频带并不能唯一补全带外响应。Lorentz 阶数、带外极点和模型错配仍应纳入不确定性。不能用有限频带数据自动证明采用的色散模型正确。

绝对材料 \(\chi_a,\chi_b\) 均被动时，相对变化 \(h=\chi_a-\chi_b\) 的虚部可以为负。把 \(\operatorname{Im}h\ge0\) 当绝对被动性会错误排除合法材料变化。

### 6.3 physics-native probability measure

一个可解释构造是
\[
\mu_0(dH,d\theta,dv)\propto
\pi_H\,\mathbf1_{\mathcal C_{phys}}\,
\exp[-\lambda_{int}\mathcal R_{interface}-\lambda_{fine}\mathcal R_{fine}]
\,\nu_H(d\theta)\nu(dv).
\tag{30}
\]

\(\mathcal C_{phys}\) 来自被动性、材料上下界与色散；界面/多尺度能量来自统计结构假设，而不伪称 Maxwell 强制要求。若 λ 是待推断 hyperparameter，还需相应 partition function，不能省掉依赖 λ 的 normalization。

物理约束的作用可以精确描述为条件化：若 C 是物理可行集，\(\pi_C=\pi(\cdot\mid C,y)\)。若基于正确信息排除了不可能材料，按正确联合生成模型平均，附加信息不能增大最优 Bayes risk；但对某个固定真实对象、错误物理约束或错误 prior，风险完全可能更差。passivity 本身也不保证 posterior 集中、唯一性或低数值条件数。

### 6.4 离散结构与连续物理统一

建议将 topology prior 与材料本构分开：K 的 prior 给无孔/多孔非零质量；每个拓扑分支再给几何和本构系数；细尺度变量使用弱收缩而非必然为零。可采用重尾尺度混合避免过度压制稀有结构，但重尾并不保证真实结构被合理覆盖。

当先验完全无法表示真实结构时，posterior predictive check 可能发现失配，也可能因为观测本来不可辨识而发现不了。不能保证“错误 prior 总会由 data residual 揭露”。

---

<a id="sec-7"></a>

## 7. T9：表示应压缩计算，而非删除未知世界

### 7.1 支持遗漏的下界

若近似后验 q 只能位于 S，则
\[
TV(\pi,q)\ge\pi(S^c).
\tag{31}
\]
取事件 \(S^c\) 即证。特别地，若目标后验在高维材料空间具有密度，而 q 被限制在严格低维光滑 manifold 上，则 \(\pi(S)=0\)，所以 TV=1。图像 MSE 很小与这个事实并不矛盾。

因此，“固定 64D 表示对 point reconstruction 很好”不能推出它忠实近似 full material posterior。可以改为声明一个 coarse/task pushforward 的误差，或者显式保留细尺度余项，而不是声称整个材料分布的 TV 很小。

### 7.2 安全降维的已知原则

写 \(x=(u,v)\)，保持 prior conditional \(\mu_0(dv\mid u)\)。已知 KL 最优的 reduced likelihood 形式是
\[
\overline L(u)=\mathbb E_{\mu_0}[L(u,v)\mid u],
\qquad q(dx)\propto\overline L(u)\mu_0(dx).
\tag{32}
\]
它保留正确的 u 边缘，v 条件采用 prior；不能用 \(L(u,\mathbb E[v\mid u])\) 替代该积分。Zahm 等已系统研究条件期望形式和可认证降维误差，不作 A32 新定理。[L7]

对模型修正也可令参考测度为旧后验：
\[
g(u)=\mathbb E_{\pi_0}[e^{-\delta(x)}\mid U=u],
\quad q_1(dx)\propto g(U(x))\pi_0(dx).
\tag{33}
\]
则 \(q_1^U=\pi_1^U\)，且
\[
D_{KL}(\pi_1\|q_1)
=\mathbb E_{\pi_1^U}D_{KL}(\pi_1(dx\mid U)\|\pi_0(dx\mid U)).
\tag{34}
\]

当 U 包含 H 时，q1 的 topology probability 与 π1 完全相同，**条件是 g 的积分正确计算**。这不是免费的压缩定理：难点被定位到条件 evidence，而不是消失了。

若每个纤维上的 \(\operatorname{osc}\delta(\cdot\mid u)\le\omega(u)\)，可得
\[
TV(\pi_1,q_1)\le
\mathbb E_{\pi_1^U}\tanh(\omega(U)/4).
\tag{35}
\]
对 topology 还可以只控制条件事件概率，不必认证所有材料细节。这提供了比全 box uniform bound 更具体的研究对象，但条件概率质量与误差界仍需要有效计算。

### 7.3 推荐混合表示，不指定唯一赢家

| 表示 | 优点 | 对 posterior 的主要风险 | A32 用法 |
|---|---|---|---|
| voxel/FEM | 不依赖固定 topology，易本构约束 | 维数、网格依赖、随机细孔 | 参考表示与细尺度余项 |
| region/piecewise constant | 界面与材料可解释 | 类别数/几何限制，复杂物体遗漏 | 最小 topology pilot |
| level-set/signed-distance | 自然允许合并和分裂 | latent prior 诱导的形状分布难读，阈值不连续 | 对照；不默认最优 |
| wavelet/multiresolution | 局部尖锐结构、可分尺度 | 截断会删除微小 topology | 主材料坐标之一，保留随机尾部 |
| PCA/Gaussian chart | 计算便宜、局部预条件有效 | 固定支持、平滑化、漏模态 | branch 内 proposal，不做全局唯一支持 |
| implicit neural field | 灵活连续表示 | 频谱偏好、复杂 normalization、难采样 | 暂不优先训练 |
| variable-dimension mixture | 支持对象/孔洞出生 | RJ 混合、label symmetry、evidence 成本 | 明确 topology posterior 的主框架 |

Bayesian level-set 的函数空间构造早有工作；persistent-homology prior 也已有新研究。但 2026 PHG 工作的对象是二维声学障碍的 log-radial 边界，其 persistence 主要约束径向函数，不等于推断任意三维材料内部空腔。[L5,L6]

拓扑不对材料 L² 误差连续：缩小孔洞可让材料差趋零，但任何“有孔”标签仍不同。要证明 topology-preserving approximation，需 feature-size/reach、阈值裕量和足够强的边界误差控制；L² 收敛本身不够。针对有限 \(\rho\) 的任务，可转为对 persistent feature 的稳定性分析，仍需显式裕量而非一个无条件保证。

**Measurement-adaptive coordinates 的合法性。**给定观测 y 后，以可逆映射 x=T_y(z) 组织计算是允许的，但必须把同一个原始 prior/likelihood 一起变换，并保留所需 Jacobian；这是改变坐标，不是获得一个新的免费先验。数据自适应的非可逆截断、只保留 N3 附近 chart 或删掉某些 topology，则改变了支持，必须显式建模或量化偏差。数据依赖的 proposal q(x|y) 也允许，但 importance/MH correction 必须使用 q 及其支持；不应将 q 当独立 prior 再乘一次相同数据的 likelihood。

---

<a id="sec-8"></a>

## 8. 文献调查与原创性地图

检索截至 2026-10-12，结合原论文全文、公式、图表与出版记录。下面区分“已核对全文关键部分”与“只取得摘要/元数据”，不把没有读到某项结果当成证明其不存在。该检索足以排除若干宽泛新颖性声明，不构成全球文献穷尽证明。

| 工作 | 实际对象与概率机制 | 对 A32 的约束 |
|---|---|---|
| Carpio–Iakunin–Stadler 2020 [L1] | 2D scalar Helmholtz；topological-derivative 初始化先验；Laplace 与 nonlinear MCMC；另做对象数 evidence 比较。全文 §8 使用先验 Monte Carlo evidence，不是联合 RJ sampler | “拓扑先验+Bayes scattering”不是新题。其 star-shaped 对象模型排除嵌套边界，不是内部空腔后验；没有报告校准的对象数概率 |
| Deep Injective Prior 2023 [L2] | 2D TMz；材料图训练、不需要配对场；式(22)–(25)用以 MAP 为固定均值的 diagonal Gaussian latent variational approximation | 是真实 Bayesian target 的受限近似，不应贬为任意候选；但单中心 latent、高维支持限制和调节 β 不能等同于跨拓扑概率校准 |
| DiffusionEMIS 2024 [L3] | 3D 点云与电参数的 conditional diffusion；取得出版社摘要和 DOI，未取得可核对的完整公式 | 不能声称“第一个3D EM generative inverse solver”；本轮无法确认其全部 vector-Maxwell 实现、精确 posterior correction 或 topology calibration |
| Back-Projection Diffusion [L4] | 2D constant-density acoustic Helmholtz；明确学习 conditional posterior；利用回投、旋转/平移结构；式(55)为 conditional score flow | 应承认其后验目标与理论结构；不能把声学结果当3D矢量 Maxwell，也不能把有限离散 BP 的无损性自动继承 |
| Physics-Guided Conditional Diffusion 2025 [L5] | 2D permittivity maps；生成多个候选，逐个 forward 检查，选最低 field MSE，全文 §II-B4 明确 | 这是物理筛选，不是计算各拓扑的 evidence；“多样图+最小残差”不能成为 A32 核心贡献 |
| PHG prior 2026 [L6] | 2D acoustic obstacle，log-radial boundary 上 persistence penalty 的 Gaussian Gibbs tilt；证明后验稳定/逼近 | “persistent homology prior+Bayes inverse scattering”也已有近作；并非任意内部拓扑/三维材料后验 |
| Zahm–Cui–Law–Spantini–Marzouk [L7] | 条件期望 likelihood 与 KL-optimal posterior reduction，具有误差认证 | 安全保留尾部、LIS、条件积分不是新发明；A32 必须说明 Maxwell 专属签名和任务收益 |
| Spantini–Bigoni–Marzouk [L8] | 非高斯低维 couplings、稀疏/可分解 transport | “用低维 transport 表示后验”不新；双侧 Maxwell 因子只是可能的具体物理来源 |
| Cotter–Dashti–Stuart [L9] | forward approximation 到 Hellinger posterior error 的理论 | 一般 posterior stability 不是新结论 |
| Chen–Hesthaven–Maday–Rodríguez 2010 [L10] | harmonic Maxwell certified RB outputs，出版记录/摘要已核对，完整 PDF 获取失败 | 不能宣称首个 Maxwell residual certificate；本轮不对其全部定理常数作优劣判断 |
| Christen–Fox；Peherstorfer–Marzouk [L12,L13] | delayed acceptance、低保真 transport proposal+高保真接受，保留目标后验 | exact correction 是必要基线，不是新颖性 |
| Jakeman 等 2026 [L14] | risk-aware、goal-oriented Bayesian design；已核对摘要与研究定位 | 普通“从信息量改为目标风险的测量设计”不能单独申报原创 |

最终补充检索还发现 2026 年《Uncertainty Quantification in Inverse Scattering Problems》[L19]，出版社摘要提出从低维到高维的渐进 Bayesian formulation。全文获取失败，本轮只确认该宽泛组织思路也已有研究，不对其未读的公式与算法作优劣判断。

补充检索发现 Wang–Liu 2021 TAP 的 hierarchical Bayesian piecewise-homogeneous electromagnetic imaging [L18]，以及 2026 aperture-flexible 3D diffusion 工作的元数据。未拿到完整原文的项目列为待核查，不据其摘要判断具体概率校准的有无。

对 Carpio 2020 的数据使用应公平描述：论文把 detector 数据分成 \(d^{(1)}\) 构造 prior 与 \(d^{(2)}\) 建 likelihood，避免同一观测的直接重复使用；这不自动等同于完整 hierarchical Bayes，但不能错误指控它完全没有 data split。Deep Injective Prior 的 β 是 diversity/未知噪声相关的调节量；其公式有明确 Bayesian 动机，局限在近似族及 calibration evidence，而不是没有概率定义。

**真正仍值得寻找的空白：**有限尺度内部 topology 在真实三维矢量 Maxwell、复材料、nuisance 存在时的 posterior evidence；以及一个可预测“模型修正会不会改变该 topology probability”的 Maxwell-specific 结构，最终减少必要的 full-wave 工作。当前文献审计没有为“首个”提供充分依据，论文应避免 first-ever 表述。

---

<a id="sec-9"></a>

## 9. 四条实质不同的算法路线

### 路线 A：无大型网络的 topology-stratified posterior inference

材料先验采用离散 H、连续几何/复材料 \(\theta_H\)、随机细尺度余项 v。先运行每分支的 tempered SMC 或可靠的多链 MCMC；分支多且变化维数时使用 reversible-jump birth/death/split/merge。SMC 用权重递推估计分支证据，再计算 H 的概率。每个 topology 都必须有正提议概率，不能只围绕 CNN 已生成的结构采样。

RJMCMC 需要匹配维数的可逆映射、proposal density、Jacobian 与正确的分支 prior；从高分图像库随机跳转不是自动合法。分支内可用 pCN 处理 Gaussian reference field，避免随网格细化全面失效；已有函数空间 MCMC 理论应作为基线。[L15]

优点是无大型训练、概率目标明确，最适合建立可信参考。弱点是复杂高维/多模态下的积分代价。必须报告各分支 ESS、证据重复误差与漏模态检查，不能只给全链一个 ESS。

**判定：主线的 reference 与第一版算法。**先证明在新的小型 vector-Maxwell 问题上能算对概率，再谈摊销学习。

### 路线 B：保留拓扑与随机尾部的 observable-material posterior

每个分支/温度下使用白化、nuisance-aware 的 material metric，构造 u；残余 v 保持 prior conditional 或通过 pCN 更新，不把 v 设为零。用 (32) 或 (33) 的条件 likelihood 组织推断；A17 双侧签名可用来引导桥接阶段或分析模型修正，而非强迫所有材料位于签名张成空间。

对光滑有限维密度 \(\pi(x)\propto e^{-V(x)}\)，正定 mobility \(M(x)\) 的 Itô 动力学
\[
dx=[-M\nabla V+\operatorname{div}M]dt+\sqrt{2M}\,dW
\tag{36}
\]
在适当边界/可积条件下保持相对于 Lebesgue 测度的 π。候选 \(M=(C_0^{-1}+J^TPJ)^{-1}\) 来自材料而非电流几何。由 Fokker–Planck 可验证 probability flux 为 \(-M(\pi\nabla V+\nabla\pi)=0\)。漏掉 divergence 项通常改变 stationary law。若采用 Riemannian volume 而非 Lebesgue，需要相应密度变换，不能混用。[L16]

离散积分还需要 MH correction 或可控偏差；边界反射/变换 Jacobian、nonsmooth prior 都要处理。曲率矩阵只用作 proposal/preconditioner，不能误把局部 Gaussian 近似当最终概率。

**判定：最值得与路线 A 组合的小算法增量。**必须与普通 branch-wise pCN/SMC 作同成本比较；如果构造 material metric 比采样省下的成本更多，则只作解释性诊断。

### 路线 C：physics-native generative proposal，而不是生成模型直接宣称 posterior

可以先生成 topology token，再生成 geometry、material coefficients 和 fine-scale field；也可以联合生成。若先采样 H，其权重必须包含 likelihood evidence，不能先按训练集 H 比例抽样然后永不跨分支修正。

自回归 token 优先用区域、孔洞/组件、材料属性和多尺度细节，而不是把 voxel 顺序假装物理因果。组件排序需 canonicalization 或 permutation-aware prior；旋转/平移要与实验坐标和 nuisance 一起处理。其每步 normalized probability 定义的是 proposal/prior 的概率，只有完成 Bayes correction 才能称 posterior。

推荐第一版可计算密度的 conditional flow/proposal q，并用完整 likelihood 做 importance/MH correction：
\[
\alpha(x,x')=\min\left\{1,
\frac{\pi_0(x')L_F(x')q(x\mid x')}
{\pi_0(x)L_F(x)q(x'\mid x)}\right\}.
\tag{37}
\]
这样流的作用是效率，不是替换统计目标。proposal support 必须覆盖目标；MH 的正确 stationary distribution 不等于有限链已混合。

Diffusion/flow-matching 只有在正确 conditional score/velocity 和相应数值条件下才产生目标分布。向 clean reconstruction 的 forward residual 求梯度，通常不等价于噪声层 \(p(y\mid x_t)\) 的精确 score。最终 data consistency 也不证明 mode mass 正确。对 implicit diffusion proposal，若密度不可得，不能直接套普通 importance weights；需要额外的合法 augmented-state 或其他推断构造。

**判定：高风险备选，不启动大型训练。**只有路线 A 建立可靠 posterior reference、路线 B 明确存在可利用的运输结构后，才考虑少量数据的 conditional transport。

### 路线 D：pointwise-exact 的 multifidelity Bayesian inference

廉价模型可以是 single-anchor calibrated-L、固定 ROM、多项式近似或已有 surrogate，不必是在线 current feedback。

标准 delayed acceptance 先用 cheap posterior ratio \(r_A\)，再用 \(r_F/r_A\) 做第二阶段接受，保留 full target；这是已有方法。允许对**本次提案与当前状态**的 log ratio 构造合法区间，若 uniform random 的接受阈值已落在区间外，可不求全精度；否则 refinement/full fallback。

该策略只需被使用的点/接受比上的有效区间，不一定需要 A30 未通过的整个材料 box uniform certificate。不过这只是避开一个过强条件，不是证明本实现会加速。若每次都 fallback，算法仍可正确，但可能更慢。

适应性的 surrogate 在生产采样时宜冻结，或者使用明确的 adaptive-MCMC 收敛条件。只说“每一步都用了 MH”不能自动解决任意历史依赖 adaptation。

**判定：可选工程层，不作为 A32 首个必须成功的物理贡献。**记录包含失败的总成本；始终保留 full-likelihood 与 calibrated-L 强对照。[L12,L13]

---

<a id="sec-10"></a>

## 10. T8：预算分配能证明到什么程度

### 10.1 有限粒子上的目标敏感度

对规范化权重 \(w_i\) 和 topology indicator \(h_i\)，\(p_H=\sum_iw_ih_i\)。若 \(\delta_i\) 是负 log-likelihood correction，
\[
\frac{\partial p_H}{\partial\delta_i}=-w_i(h_i-p_H).
\tag{38}
\]

因此 \(w_i|h_i-p_H|\times\text{error width}/\text{cost}\) 可作为 first-order refinement 优先级。它是启发式/局部敏感度，不是全局最优 allocation。两个修正或两个新测量可以互补，greedy 不必最优。

若有限 bank 的未归一化权重在 \([l_i,u_i]\)，H 的确定区间为
\[
p_H^- =\frac{\sum_{i\in H}l_i}{\sum_{i\in H}l_i+\sum_{i\notin H}u_i},
\quad
p_H^+ =\frac{\sum_{i\in H}u_i}{\sum_{i\in H}u_i+\sum_{i\notin H}l_i}.
\tag{39}
\]
这只认证这个有限 bank；连续模型还需要积分误差、漏模态和尾部误差。

### 10.2 一个可以精确证明的有限 allocation 结果

若各分支有独立无偏估计量，目标估计方差为 \(\sum_H V_H/N_H\)，每样本成本 cH，总预算 \(\sum_H c_HN_H=B\)，连续放松下最优
\[
N_H\propto\sqrt{V_H/c_H}.
\tag{40}
\]
对 \(N_H\) 做 Lagrange 一阶条件即可。这里 VH 已包含所估任务与分支权重，不能拿普通场 variance 代替。离散、相关 MCMC、估计证据 ratio、适应性选择时不直接适用。(40) 是已知 stratified allocation 原理，不是新的 Maxwell 最优性。

### 10.3 不可能免费的“最优预算定理”

只有误差上界而没有误差分布，无法证明一次 full call 会产生正收益：真实误差可能为零。只能认证“最大可能决策损失已经很小，可以跳过”，不能把一个大上界解释为“值得计算的收益很大”。

新增 full Maxwell evaluation 不增加实验本身的信息 \(I(\chi;Y)\)，只降低我们的计算近似误差。计算价值与新测量价值必须分别记账。

---

<a id="sec-11"></a>

## 11. 从 posterior 到材料决策与统计校准

### 11.1 输出应是什么

推荐输出：\(P(H_\rho=1\mid y)\)、不同拓扑分支的质量、每分支的代表性材料样本、空腔 occupancy field、材料参数的条件/边缘可信区间、prior sensitivity 与拒判状态。

孔洞位置场可定义
\[
p_{void}(r)=P(r\text{ 位于满足尺度定义的内部空腔}\mid y).
\]
孔洞中心的**无条件**位置测度质量是 \(P(H=1\mid y)\)；条件在“孔洞存在”后的位置 density 才归一为 1。多孔情况下需要随机有限集/计数测度或标签不变总结，不能把相互排斥的多个孔洞位置平均成一个假孔洞。

连续材料中，单张精确图像的概率通常为零；可给 density、邻域质量、mode/basin mass 或 topology mass，但必须声明是哪一个。正确目标且充分混合的样本频率可以估计事件概率；任意 best-of-N 生成频率不能。

### 11.2 不同任务的 Bayes rule

对 hole detection，误报成本 CFP、漏报成本 CFN：报告有孔当
\[
p_H>\frac{C_{FP}}{C_{FP}+C_{FN}}.
\]
平方材料损失用 posterior mean；绝对损失用相应 median；非凸材料行动集合需要最小化期望损失而不是直接选 MAP。整体 topology 对象与像素误差是不同任务，不应共享一个未声明的评分规则。

多模式 posterior mean 可能把“实心”和“空心”平均成半密度核心；它是平方损失下的合法行动，却不应作为“物体最可能真实形态”展示。应同时返回 topology-conditioned summary。

### 11.3 校准含义

在正确的联合模型下，令 \(p(Y)=P(H=1\mid Y)\)，则
\[
\mathbb E[\mathbf1_H\mid p(Y)]=p(Y).
\tag{41}
\]
这是 conditional expectation tower property，不是每个固定真物体的高概率保证。模型错误、支持遗漏、未混合采样和数值偏差会破坏实际校准。

使用 Brier/log score 等 proper scoring rules，按 topology、feature size、材料损耗、nuisance/OOD 分组报告。SBC 用于检查在声明 prior-generative model 下的推断实现；SBC 通过不能证明真实物体来自该 prior。continuous ranks 与 discrete topology probability 的检查需要分别设计。[L17]

pointwise credible interval 与 simultaneous material band 不同；posterior predictive checks 能发现部分失配，但无法保证识别所有错误 topology。必须单列 reference integration error、Monte Carlo error、forward discretization error、prior misspecification；不能把 parent bootstrap 区间当成包含这些误差。

### 11.4 拒判应区分原因

若对称 0–1 错误损失、拒判成本 \(c_A<1/2\)，则当
\[
\min(p_H,1-p_H)>c_A
\]
应拒绝二选一。若概率的计算区间跨过决策阈值，返回 NUMERICAL_UNRESOLVED；若概率可靠但集中在中间，返回 OBSERVATION_AMBIGUOUS；若不同合理先验给出相反结论，返回 PRIOR_SENSITIVE；若表示缺少合理支持，返回 COVERAGE_UNRESOLVED。

观测不可辨时应报告 Bayes factor 接近 1 与后验接近先验，而不是一律输出 50%——只有先验本身 50/50 时才应为 50%。

---

<a id="sec-12"></a>

## 12. T10：新测量应区分什么

新设计 d 的 topology information 为
\[
I(H;Y_d\mid y),
\]
等价于预测测量混合分布的 generalized Jensen–Shannon divergence。它不同于 \(I(\chi,\gamma;Y_d\mid y)\)；后者可能主要学习对最终任务无关的细节或 calibration。

更直接的任务目标是
\[
\Delta R(d)=R^*(y)-\mathbb E_{Y_d\mid y}[R^*(y,Y_d)].
\tag{42}
\]
正确联合模型下 \(\Delta R\ge0\)，因为可以选择忽略新数据。保证是对未来数据的期望，不是每个已观测 realization 都会降低 posterior risk。一般 goal-oriented/risk-aware design 已有丰富研究及 2026 新作，不作为 A32 的独立新颖性。[L14]

对于固定两材料、独立 Gaussian 通道，
\[
d_{all}^2=\sum_k d_k^2,
\]
可把 (6) 的 full-wave secant 直接用来筛选频率/极化/源。对复合分支，应比较 posterior predictive mixtures，不能只用两个 MAP 的 secant 代替。

未知 gain、shared clock、translation、频率色散会改变联合 covariance 与可区分方向。给定 nuisance 后独立的频率，边缘化 nuisance 后可能相关；不能盲目相加 Fisher information。高频对细节更敏感也可能对误差更敏感，最高频不必有最高任务收益。

理论上可研究以跨拓扑有限割线为主、局部 tangent 为辅的设计代理，并用误差界控制其真实风险排序。若取得的只是普通 D-optimal design 在一个模型上的改名，则没有新的贡献。

---

<a id="sec-13"></a>

## 13. 本轮实际完成的小规模验证

所有新代码见 `a32_checks.py`，结果见 `A32_CHECKS.json`，初版保留于 `A32_CHECKS_INITIAL.json`，逐调用记录见 `A32_RESOURCE_LEDGER.json`。未读取历史 truth arrays、locked ID/OOD；未训练网络；未调用 GPU；未运行 current-selection sweep。

### 13.1 代数与已有结果复核

A31 原 22 项有限维/解析检查在独立目录重跑，22/22 通过；其中只有一个 9 维矢量偶极物理状态、两条 RHS，不能称为新的三维材料反演。结果保留为 `A31_RECHECK.json`。

A32 新增 14 项代数/解析/数值检查通过。关键数值：

| 核验 | 结果 | 正确解释 |
|---|---:|---|
| paired potential identity | 最大误差 2.66e−15 | 冻结 Gaussian likelihood 的精确代数 |
| paired posterior weights | 最大误差 5.55e−17 | 同一有限非高斯 prior 下重加权一致 |
| Gaussian mean/covariance 的支撑外分量 | 2.38e−16 | 仅 Gaussian 特例 |
| 签名独立的 topology | 材料 TV=0.360818，P(H)=0.3 不变 | 存在“后验变化大、该任务概率不变”的反例 |
| 签名相关的 topology | 概率变化 0.249294；界 0.253221 | 任务耦合界在一个小例中非空 |
| A30 父样本再聚合 | 13.1225%、−3.3154%、35.6896%、21.9074% | 与 pinned 原 CSV 一致，不是重跑原实验 |

随包还检查了固定 field offset 在零 Jacobian defect 下改变概率、同峰值 residual 不同 evidence、fibers 上 prior conditional 不变但 marginal 改变、遗漏 Itô divergence 的偏差等反例。

### 13.2 新建二维解析全波圆环例

独立实现无限介质圆柱截面的 TM \(E_z\) 分离变量解。外半径 R=0.35，材料 \(\epsilon_r=2.5+0.08i\)，内部为真空空腔，外部真空。分别匹配两个界面的场与径向导数连续性；每个 Fourier harmonic 解 2×2 或 4×4 线性系统。使用 \(e^{-i\omega t}\)，包含完整多重散射，不是 Born 近似。物理对象仍是 **二维 scalar TM reduction，不是三维 vector-Maxwell 实验**。

真空散射控制约 5.96e−17；16 与 20 阶截断在所测接收场上相同至浮点精度。边界 backward residual 有记录。它们是实现核验，不是 rigorous interval error proof。

固定 k=4、8 个角度接收、已知材料、10% RMS complex noise，比较 solid 与不同空腔固定对。最优等先验错误率直接按 (7) 计算：

| 空腔半径/R | 白化距离 d | 最优固定二点判错率 |
|---:|---:|---:|
| 0.01 | 0.004966 | 49.9010% |
| 0.03 | 0.044503 | 49.1124% |
| 0.10 | 0.479297 | 40.5301% |
| 0.30 | 3.821851 | 2.8007% |
| 0.55 | 12.389921 | 约 2.91e−10 |

这支持有限尺度问题的必要性：物体都被动，solver 全波，极小空腔仍几乎不可能二选一。它不是一个普适频率/孔径分辨率曲线，也不是把 unknown material/geometry marginalize 后的 topology error。

### 13.3 有限材料 prior 的概率检查

建立 18 个明确 prior states：无孔/有孔两类，每类 9 个介电实部值，虚部固定 0.08。这里**声明的是有限 prior**，不是用 18 个粒子冒充连续材料后验。

初始双频例几乎完全可分，不能检验中间概率的校准；还发现 probability=1+roundoff 落出最后 bin 的边界统计问题。保留初版输出，不把它称 calibration PASS。修复纯数值 bin guard，并另设单频、30% 噪声的探索性难度检查；该修改不是独立确认实验。

在单频例的一个有孔真值 realization 上，计算得到 P(hole|y)=0.20233；加入 proper Gaussian gain uncertainty 后约 0.20682。真实有孔不意味着每次噪声 realization 都应给高概率；校准是重复生成层面的性质。2000 次 prior-predictive 数据模拟只复用已算场，没有再解物理；有限模型的 Brier score 约 0.10281。五个概率 bin 的 mean probability 与频率分别约为：0.056/0.041、0.300/0.291、0.502/0.494、0.706/0.715、0.935/0.925。只把它当生成—推断链条的一致性检查，不声明连续后验或真实物理分布已经校准。

两次脚本实际共调用 **74 次解析材料–频率评估**，包括重复控制和初版，低于本轮独立设定的 80 次上限；两次核心计算分别约 0.127 与 0.090 秒。该时间只覆盖脚本核心，不包括文献检索、代码编写和报告。A31 的微型偶极复核另列，不混入新圆环计数。

**本轮没有获得的结果：**任意连续材料后验参考、真实三维复杂拓扑校准、SMC/flow 的同成本胜出、在线 feedback/certificate 加速。这些不是被上述算例证明的结论。

---

<a id="sec-14"></a>

## 14. 下一阶段的最小决定性实验

以下是设计，不是本轮已经执行的研究。所有新对象先定义分布与随机种子，单独登记，不打开此前 locked sets；精确/近似方法共用 prior、表示与数据。

### E0：先让参考后验可信

优先使用**三维同心分层介质球的矢量 Mie 解**，同时包含电/磁多极系数及极化响应，比较实心与有固定最小尺度内部空腔。它是三维矢量物理，但几何受限，不能代表一般复杂目标。

第一批只设两种 H、每分支两个连续材料变量 Re ε 与 Im ε，保留相同 proper material prior；先不同时引入十几个边界自由度。用独立实现/解析极限核对 forward，再用 nested quadrature 或可验证的自适应积分核对 ZH。不要用某一次 MCMC 的样本均值当 reference truth。

冻结门槛可采用：两独立积分/细化方式的 topology probability 差不超过 0.02，log Bayes factor 差不超过 0.1；这些是研究容差，不是保证全局真实误差的定理。若能给严格 quadrature remainder 则另标为认证。控制在一个预先批准的小预算内，例如总计 384 次新增材料–频率评估，包含标签、重复、失败和检查。若预算不足以建立 reference，返回 REFERENCE_UNRESOLVED，禁止转入“我们的后验更准”的比较。

同预算下优先付给少量可靠概率点，而不是增加 parent 数。A15 的失败表明“画出 posterior”不是建立 posterior reference。

### E1：直接检验本轮最有价值的机制

在上述共同 posterior 上比较 ordinary full F、single-anchor calibrated-L 和一个预先固定的近似模型。冻结 U 后，在离线评价中计算真实 \(\delta\)、paired \(h(U)\)、剩余 \(\epsilon\) 和 \(\tau_H\)。该离线分析可以使用 full 信息，但不能当在线特征。

比较四个量：forward error、oscillatory likelihood-ratio error、签名–topology 耦合、实际 \(\Delta P(H)\)。设置共同 offset、可被 gain 吸收的偏差、只影响分支内材料的偏差和跨分支偏差四个干预。必须包含跨拓扑有限变化，不能只看零残差附近。

**否证条件：**若局部签名的非线性余项通常与主项同量级，或 \(\tau_H\) 计算需要接近完整 posterior 成本，签名不能作为主部署加速器。可保留其解释定理，但不要追加网络硬学一个难以验证的 predictor。

### E2：安全表示 versus 强制小 latent

比较固定 Sparse/PCA chart、同 chart+随机多尺度尾部、显式 topology mixture、level-set。使用相同物理 prior 信息，不把某一表示的 oracle truth directions 注入在线。

主指标不是最小 MSE，而是是否遗漏正确拓扑、topology probability 的 reference error、task calibration 和材料 credible coverage。必须有一个超出主表示的负控制。

**否证条件：**压缩降低 MSE 却系统性漏掉一个真实 topology，不能称 probabilistic representation 改进。复杂 geometry 的支持若不足，停止模型精度比较，标为 REPRESENTATION_LIMITED。

### E3：是否真的值得多算 Maxwell

只在 E0–E2 合格后比较 always-F、calibrated-L、标准 delayed acceptance 和 candidate pointwise certificate/decision policy。所有方法目标为相同 full posterior；若使用近似目标，另列偏差。

统计总 wall、setup、material metric、certificate、primal/adjoint、失败、缓存、积分与 mixing 成本。以 probability tolerance 或 Bayes decision regret 对齐质量，不能用相同 current rank 对齐成本。

**否证条件：**相同概率精度下 proposed 全成本不优于标准 DA/always-F；或者 full F 本身相对 calibrated-L 不改变关键材料决策。遇到后一种情况，应采用 calibrated-L，而不是为保留 feedback 继续加大模型。

### E4：主动测量只在确有分支竞争时测试

用冻结 posterior 的两个或多个竞争性 topology predictive distributions，选择少量新频率/极化。对照同预算随机设计、material Fisher 设计、topology MI、task risk 设计。新数据独立生成后验证实际 Bayes risk；在 free gains/dispersion misspecification 中做负控制。

若初始后验已几乎确定，主动测量的增益小属于合理结果。若错误 prior 排除了真结构，所谓设计收益不能解释为消除了完整物理歧义。

---

<a id="sec-15"></a>

## 15. T1–T10 的最终裁决表

| 候选 | 本轮可成立的具体结论 | 假设/失败点 | 原创性与算法用途 |
|---|---|---|---|
| T1 observational equivalence | 同测量分布的等价关系；常秩局部纤维 | nuisance marginal 的等价依赖 prior；近似距离不传递 | 基础定义/已知几何；约束研究目标 |
| T2 topology identifiability | 精确全波割线+固定二点 Bayes error；无最小尺度时 uniform detection no-go | 统一场/伴随稳定；复合类需积分而非单对距离 | 组成工具已知；真实 Maxwell 有效常数/材料机制可能有价值 |
| T3 prior–likelihood decomposition | 精确响应纤维上的 prior conditional 不变 | 弱方向 marginal 可由先验相关更新 | 已知 Bayes 事实；提供 prior honesty 的可检查表述 |
| T4 topology stability | bounded log-likelihood/prior tilt 控制任意事件概率 | 需共同支持/有效尾部；material L²/Wasserstein 小不够 | A31/既有稳定性理论；不能作为唯一新定理 |
| T5 physical posterior | proper physical prior、条件化、色散共享 | 被动不保证 coercivity/唯一性；Poynting 常为重复 | 物理建模与实测约束有价值，普通截断不新 |
| T6 feedback-to-posterior | (22)–(29) 及条件 evidence/odds | 实际误差界与积分需要支付成本 | 普通证书已知；Maxwell-specific 低成本对比界仍待证 |
| T7 material-aware transport | (14)–(19) 双侧签名的精确 posterior change；非线性扩展(20) | 冻结共同模型，维数/条件采样成本，跨分支 normalization | 最值得继续的理论连接；一般条件 tilt 不申报新颖 |
| T8 budget optimality | 特定独立 strata 下(40)；有限 weights 下(38)–(39) | 一般 adaptive/global 最优性不成立；大误差上界不代表收益 | 只给有限条件结果，拒绝过度声明 |
| T9 topology-preserving representation | 支持遗漏下界、保留随机尾部、有限尺度稳定条件 | 固定 manifold 对全后验可能 TV=1；任意小孔不稳定 | 安全表示的具体 Maxwell/task 实证可贡献，通用条件期望已知 |
| T10 information-guided design | fixed-pair secant、predictive mixture MI、risk reduction | 非线性复合类、shared nuisance、未知色散 | 已有 Bayesian design；须有新的 EM 几何/代价优势 |

最小反例已分别给出：\(x^2\) 的全局分支；微孔洞的尺度消失；先验相关导致 marginal 更新；零支持使 posterior 永远排除真结构；零 Jacobian error 的常数 forward offset；非高斯先验破坏 Gaussian mean-support 推论；独立 topology 导致物理修正无决策收益；低维 manifold 的 TV=1；纯 positivity 不提供严格正耗散；单点 residual winner 不等于 evidence；未纠正 variable mobility 改变 invariant law。

---

<a id="sec-16"></a>

## 16. Final research decision

### 16.1 最值得推进的一条主线

**有限尺度 topology-stratified material posterior + Maxwell 模型修正的双侧材料签名。**

第一版推断用无大型网络的 SMC/MCMC/低维积分，完整 likelihood 作为可检查目标；material geometry 用于 proposal/preconditioning，随机尾部保留不确定性；feedback 先用于解释与判定 probability change，而不是要求每一步在线执行。

这比“再做一个 posterior generator”更具体：它把 A17 的解释对象从 GN step 延伸到 likelihood ratio 与材料事件，又直接面对 A30 的 capacity-versus-selection gap 和 A15 的 reference 失败。

### 16.2 高风险、高理论回报的备选

**非线性、topology-conditioned likelihood-ratio transport**：寻找 U=(H,少量 material coordinates)，在跨拓扑 credible regions 上对 \(\delta-h(U)\) 的条件 oscillation 或条件 evidence 误差给出紧界，保留未压缩材料条件分布。

如果成功，可获得比全局 forward norm certificate 更任务相关的可信推断。但它可能因弱耗散、近共振、宽 posterior、复杂 prior conditional 和 signature 构造成本失败。允许最后结论是：该结构适合离线 diagnostics，不适合在线加速。

### 16.3 明确终止或不启动的路线

不再扩展 OPM/current-selection sweep；不把固定 Sparse-64 作为整个 posterior support；不对旧 CNN candidate bank 的 softmax 赋予未经修正的 Bayesian 概率；不把 energy residual 作为独立新证据；不重复 A30 已失败的全 box uniform Neumann certificate 主流程；不在 reference posterior 未建立前启动大型 diffusion/flow training。

这些停止针对当前研究投入，不是断言相关技术在所有问题上不可能有效。

### 16.4 最先需要证明什么

先把命题 3 写成完整 measure-theoretic statement，明确 arbitrary common prior、realification、固定约束、分支 normalization 和非线性 remainder；然后证明一个**在实际 Maxwell credible region 上可计算的任务级误差控制**。若后者只能依赖完整高精度 posterior 或全 resolvent，它没有部署优势，必须直说。

### 16.5 最先需要运行什么

不是新网络，也不是大批复杂图像，而是新建三维 vector-Mie/小型 vector-Maxwell 的实心–空腔**低维可信 evidence reference**。先检验“多个分支是否真的在噪声下竞争”和“我们的概率积分是否可靠”，再看 dual signature 是否解释 full-versus-cheap posterior 的差异。

### 16.6 论文最可能成立的贡献与失败原因

可能成立的贡献组合：明确的有限尺度 material/topology 推断对象；full-wave secant 与辨识难度的物理分析；双侧材料 likelihood-ratio 签名与 topology coupling 定理；保留支持的实际推断算法；在真实三维矢量、复材料、nuisance 条件下验证 calibrated event probabilities 和可计费的计算效益。

可能失败的原因：相关理论只是已知工具直接组合而缺少新的 EM 内容；拓扑竞争在所选实验中根本不出现；连续 reference 无法在预算内收敛；表示遗漏真结构；签名不足以描述有限非线性变化；校准只在过窄 prior predictive 分布成立；或 proposed 成本高于标准 full/DA。

**研究 GO：建立概率对象、完成上述理论连接与小型参考实验。算法 superiority、复杂三维校准、在线认证加速：尚未建立，必须由新证据决定。**

---

<a id="sec-17"></a>

## 17. 原始来源、文献与阅读边界

### 17.1 项目来源

[S1] `PAPER_TAP_FREEZE.pdf`，Conditional Current Exchange for Material Updates in Maxwell Inverse Scattering，上传冻结主稿。重点：§II–VI、式(2)–(7)、表I。

[S2] `PAPER_SUPPLEMENT_FREEZE.pdf`，上传补充材料。重点：§I–IV 物理/实化/证明；§V.3 未实例化证书；§VI 信息指标；§VIII–X 数据范围与非线性；§XII–XIV 成本与物理局限。

[S3] `migodam/a30-maxwell-observable-material-manifold`，commit `717d3d26a290248bd2b9aeb55f3d03694b407b84`。读取 `A30_E0_E1_DECISION.md`、`results/a30/FINAL_DECISION.json`、`E0_GATE.json`、`PARENT_STATISTICS.csv`、capacity summary 的相关部分。再聚合的 CSV 行：full_hull 为50–61；full_L_mean 为114–125；sparse_64 为274–285；single_chart_64 为322–333。以仓库实际 header 计行，下载器的显示行号可能不同。

[S4] `migodam/a28-maxwell-material-information`，commit `c7819fa142bb70f9920835e8e91ba145f15ed8e3`，`E1_DECISION.md`。

[S5] `migodam/a29-material-hypothesis-coverage`，commit `5d5e31fcfc56ac8ef9d64c46b4b7bbdedef96bb4`，`A29_DECISION.md`。

[S6] A18-B、A19、A20/A21 原仓库 README/裁决入口；仅核对已声明结果与范围，不声称复核全部训练或代码。

[S7] `a23-fast-material-imaging`，commit `a47cda424dfb55792e69d952a25f54ce9a2a5613`，README/实验裁决。

[S8] `a24-source-anchored-om-feedback`，commit `ca80ca2a29867e1b03589061dd5c17688c9acbe0`，README。

[S9] A25 commit `59d8cc57c9da99ccb7978febc308193c43d6f3e5`；A26 commit `2603f3bd16c48e48587fedc4363975779def17f5`；各 README。

[S10] A27 commit `325cc8a5f99d686ce1b1570296f94b0df1e90141`，README。

[S11] `a14-task-relevant-posterior-coordinates/README.md`，本轮读取的 main 入口。

[S12] `a15-fullwave-posterior-ambiguity/START_HERE.md`，本轮读取的 main 入口；完整参考积分记录未重新计算。

[S13] 用户上传 `A31_Research_Package(2).zip` 内 `A31_THEORY_REPORT_ZH.md`、`a31_matrix_checks.py`、`A31_MATRIX_CHECKS.json`。报告全文读取，脚本独立目录复核。

[S14] `REPORT_TO_XUDONG_FINAL_ZH.pdf`，上传中文研究汇报；重点核对 A12/Krylov 成本、广义 current-selection 失败、A17 改进的证据层级。

### 17.2 文献标识

[L1] A. Carpio, S. Iakunin, G. Stadler. Bayesian approach to inverse scattering with topological priors. Inverse Problems (2020). DOI `10.1088/1361-6420/abaa30`；arXiv `2003.09318v2`。全文关键方法、模型选择及数值范围已核对。

[L2] A. Khorashadizadeh et al. Deep Injective Prior for Inverse Scattering. IEEE TAP (2023), 8894–8906. DOI `10.1109/TAP.2023.3312818`；arXiv `2301.03092v2`。全文后验式(22)–(25)与实验已核对。

[L3] X. Bi, Y.-J. Chen, L.-L. Li. DiffusionEMIS: Diffusion Model for 3-D Electromagnetic Inverse Scattering. IEEE TGRS (2024). DOI `10.1109/TGRS.2024.3349681`。出版社记录、摘要与元数据已核对；全文未取得。

[L4] Back-Projection Diffusion: Solving the Wideband Inverse Scattering Problem with Diffusion Models. arXiv `2408.02866v4`。模型、架构、posterior flow 的关键全文部分已读；并非全文每项证明逐行复核。

[L5] S. Chehelgami, J. Lovetri, V. Khoshdel. Physics-Guided Conditional Diffusion Networks for Microwave Image Reconstruction. arXiv `2510.25729v1`。全文二维设定与 §II-B4 选择机制已核对。

[L6] Z. Wang et al. A persistent-homology-Gaussian prior for solving infinite-dimensional Bayesian inverse scattering problems. arXiv `2607.15289v1`；以及 M. A. Iglesias, Y. Lu, A. M. Stuart. A Bayesian Level Set Method for Geometric Inverse Problems. arXiv `1504.00313`。前者已核对 log-radial 模型和 posterior 定理假设，未逐行审查所有证明；后者作为既有方法来源。

[L7] O. Zahm, T. Cui, K. Law, A. Spantini, Y. Marzouk. Certified dimension reduction in nonlinear Bayesian inverse problems. arXiv `1807.03712`；Mathematics of Computation，DOI `10.1090/mcom/3737`。条件 expectation/KL 最优性式与 PDF 第6页已直接核对。

[L8] A. Spantini, D. Bigoni, Y. Marzouk. Inference via Low-Dimensional Couplings. JMLR 19(66), 1–71 (2018)，paper identifier `17-747`。低维 transport 理论范围已核对，不声称全部71页证明完成审计。

[L9] S. L. Cotter, M. Dashti, A. M. Stuart. Approximation of Bayesian Inverse Problems for PDEs. arXiv `0909.2126`。用于已知 forward-to-posterior stability 归属。

[L10] Y. Chen, J. S. Hesthaven, Y. Maday, J. Rodríguez. Certified Reduced Basis Methods and Output Bounds for the Harmonic Maxwell’s Equations. SIAM J. Sci. Comput. 32(2), 970–996 (2010). DOI `10.1137/09075250X`。本轮全文下载失败，仅用来确认既有主题，未据其未读公式比较常数。

[L11] M. Cassier, P. Joly. An operator approach to the analysis of electromagnetic wave propagation in dispersive media. Part 1: general results. arXiv `2504.00239`。因果/被动色散的背景文献；本报告的 Lorentz 参数例是明确给定的建模构造。

[L12] J. A. Christen, C. Fox. Markov Chain Monte Carlo Using an Approximation. JCGS 14(4), 795–810 (2005)。既有 delayed acceptance。

[L13] B. Peherstorfer, Y. Marzouk. Multifidelity sampling with transport maps. Advances in Computational Mathematics (2019). DOI `10.1007/s10444-019-09711-y`。

[L14] J. D. Jakeman et al. Risk-Aware Goal-Oriented Bayesian Optimal Experimental Design. arXiv `2609.08023`。已读定位/摘要与PDF首页，不声称重证全部风险定理。

[L15] S. L. Cotter, G. O. Roberts, A. M. Stuart, D. White. MCMC Methods for Functions: Modifying Old Algorithms to Make Them Faster. Statistical Science (2013). DOI `10.1214/13-STS421`；arXiv `1202.0709`。

[L16] M. Girolami, B. Calderhead. Riemann manifold Langevin and Hamiltonian Monte Carlo methods (2011)。本报告 (36) 的 Lebesgue 目标平稳性由 Fokker–Planck 直接推导，不依赖含糊的“natural gradient 自动正确”主张。

[L17] S. Talts et al. Validating Bayesian Inference Algorithms with Simulation-Based Calibration. arXiv `1804.06788v2`。

[L18] F.-F. Wang, Q.-H. Liu. A Hierarchical Bayesian Inversion Method for Electromagnetic Imaging of Inhomogeneous Objects With Piecewise Homogeneities. IEEE TAP (2021), 2903–2912. DOI `10.1109/TAP.2020.3028212`。仅元数据/摘要，全文待核对。

[L19] Uncertainty Quantification in Inverse Scattering Problems. Entropy 28(4), 461 (2026). DOI `10.3390/e28040461`。最终补充检索取得出版社摘要：讨论由低维到高维 Bayesian formulation 的渐进框架；全文页面获取失败，只记录为进一步核查项，不据此比较具体定理或实现。

### 17.3 交付与复现

本报告是独立 A32 理论包，不修改任何上游项目。`A32_RESOURCE_LEDGER.json` 已记录本轮74次解析调用；在原目录再次完整运行会触发80次累计上限。独立复现应复制**脚本本身**到明确的新空目录，生成新的结果与账本，不删除或覆盖本轮历史账本。依赖为 NumPy/SciPy；无 GPU 依赖，无网络/用户文件读取。

本轮失败/边界记录：容器外网下载不可用；部分论文全文无法获取；最初 interactive file-writing 会话不受环境支持，改为常规文件写入；初版双频 posterior 太容易且 calibration bin 缺少浮点边界 guard，原结果保留并在新单频探索性检查中修复。任何这些修复都没有被计为成像或理论 superiority 的证据。

最终交付检查：修复分段写入造成的3个标题换行，补全 topology evidence/Occam 推导、measurement-adaptive coordinates 的测度要求和非高斯均值支撑反例；平方损失下受约束 Bayes action 明确为存在时的全局度量投影。support-exclusion 的预设数值断言改成显式有限分布 TV 计算，并仅重跑代数，新增物理调用0；其他核验数值保持一致。
