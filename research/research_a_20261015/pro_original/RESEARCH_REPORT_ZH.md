# Research A：Maxwell 补偿驱动的多解材料成像
## 从基础重建到可实现、非线性补偿、可交接的材料候选集合

**研究裁决与下一轮执行规范｜2026-10-14**

**状态：理论研究完成；独立八偶极 CPU 代数检查完成；新的 Mie、二维成像、三维成像及训练 campaign 均未运行。**

本文把已有实验证据、标准数学推论、本轮条件推导、待验证算法建议分开。公式成立不等于成像有优势；在离散方程下验证不等于连续介质误差已经受控。文献编号 R 与项目编号 P 的出处及访问层级见文末。配套文件包括执行 Prompt、候选接口 JSON Schema、可运行的小型检查和原始结果。

---

## 0. 最终判断

**值得推进，但不应将主线命名为“TSOM ambiguous-current generator”。最值得测试的是：material-first 的有限结构编辑，经过受约束非线性补偿，形成全波核验的候选集合。TSOM 是可替换的基础重建器；ROM 是需要证明省算的计算部件；神经网络暂不进入主算法。**

推荐主算法：

$$

\boxed{\text{基础材料重建}
\to\text{有限几何/材料编辑}
\to\text{secant 可补偿性诊断}
\to\text{非线性材料/仪器补偿}
\to\text{全波核验}
\to\text{结构覆盖选择}}

$$


结构不同的备选是**不依赖单一基础解的 population search / tempered sampling**，使用全域多尺度材料参数和显式 birth/death 编辑。作为经验性候选搜索时可以 candidate-only；只有正确处理 prior、proposal 与接受概率时才谈 posterior sampling。

最可能有价值的机制不是“存在很多图片”，而是：

1. **Material realizability**：严格限制多照明电流扰动来自同一个合法材料，避免把电流核空间当作材料核空间。
2. **Finite-edit compensation**：直接分析一个有限结构编辑产生的完整散射 secant，并判断允许的材料和 nuisance 是否能吸收它，而不是把冻结 Jacobian 的弱方向无限外推。
3. **Localized full-response reuse**：固定离散网格上的局部材料编辑可通过小型 Schur 核保留全部重复散射。其优势只在修改支撑小、参考算子可复用、验证与刷新成本不吞掉收益时存在。

**本轮不能给出“Best-of-4 将改善多少”的实测数字。** A30 的 35.69% 是特定表示和旧评价指标的离线真值拟合收益，既不是当前算法的实测收益，也不是所有方法的上界。下一轮必须直接测量完整材料 MAE，而不是替换指标后沿用旧百分比。[P4]

**首个停止试验**：在 A30 原二维困难条件下，用相同预算比较有限编辑补偿、随机弱方向补偿和普通 multi-start material GN。若有限编辑仍产生不了更多合法结构族，或 Best-of-4 MAE–成本没有优势，不进入 ROM 优化和网络训练。

---

## 1. 证据审查：哪些结果可以继承

| 证据 | 已经支持的内容 | 不能继承的主张 | 对 Research A 的约束 |
|---|---|---|---|
| A17 closure_final | 给定状态和二次目标的 current exchange 改善 material-step fidelity；十个完整额外对象 gap 中位下降 50.7% | 不是 50.7% 成像 MAE 改善；不是普遍加速或材料多解定理 | 复用反馈结构和模型误差分析，不把 current basis change 当成材料世界变化 |
| A17 nonlinear transfer | 三组完整轨迹复材料误差相对下降 10.84%、0.22%、7.16% | 总时间分别增加到 2.68、1.74、1.46 倍；近接触实部略差 | 新算法必须核算拒绝、修正、验证等额外动作 |
| A23 | 有限相对反馈恒等式；局部二阶修正；source anchor 下特定 P seed 消失 | P 的 forcing 冗余不等于传播消失；低参数维数不等于低成本 | 可以不用 OPM，以真实 finite material response 为起点 |
| A24 | 24×24 TM 平台揭示初始化、成像质量及 ROM 任务误差问题 | Gate2 FAIL，质量门槛未达；没有同质量 speedup；三维 NOT_RUN | ROM 的默认状态应是“尚未证明值得使用” |
| A30 | Sparse-64 有离线容量；局部补偿能缩短距离 | 最小 full-F 距离仍 8.84；复杂对象 9.18；未发现原门槛所需模糊对 | 原失败必须保留；大编辑是新假设，不是重新解释旧失败为成功 |
| Joint P2 | 同心分层球中，中心空腔可由介电常数及实幅度 gain 部分补偿；受限后验保留不同材料族 | 不证明一般三维非唯一性；环形编辑有失败；未支持普遍 $r^3\to r^5$ 消除 | 用作校准、物理机制和接口试验，不作开放三维成像证据 |
| C1-P1/P2 | 低维密度蒸馏能够保留某些材料事件；P2 full/diag 事件 MAE 为 1.3215%/1.3028% | 事件内部条件 CDF 误差约 14.18%/16.18%；没有稳定 state/memory 增量；没有证明材料 MAE 好 | 学习 compensation proposal 可能更经济，但仍须单独证明 |

来源：[P1–P6]。A17、A30、C1 的原“material error”“event MAE”和本报告建议的 volume-normalized material MAE 不应混为同一指标。

### 1.1 Joint P2 的精确解释

四个 latent 为 $z_c,z_s,\epsilon',\epsilon''$，不是四个始终可辨的材料坐标。空腔未激活区间中不同 latent 解码为同一材料，因此必须在**材料场商空间**上去重、比较距离和统计候选数。latent KL 的变化不能直接解释成 material-field KL 的同等变化。

单频、$r=0.25R$ 的固定对：

$$

d_{\rm raw}=1.529652,\quad d_{\epsilon}=1.295221,\quad
d_{\epsilon,g}=0.900170.

$$

加入第二频率并使用共享 gain 后，距离为 2.467944。这里 $d$ 是固定白化 convention 下的**两预测均值距离**，不是单个候选对有噪声观测的残差阈值。[P5]

此外，$0.25R$ 球形空腔体积分数只有 $0.25^3=1.5625\%$。因此“全域材料 MAE 已很低”和“空腔恢复正确”可以同时相反。必须加入局部结构指标。


原 target 还应原样冻结：$z_c,z_s\sim U[-3,1]$，$\epsilon'\in[2.2,2.8]$，$\epsilon''\in[0.04,0.16]$；未激活空腔区间各有先验质量 1/2，实心材料的先验质量为 1/4。gain 为 $[0.7,1.3]$ 上截断的实数 $N(1,0.15^2)$，两频共享。每频 8 个角度、两偏振、16 个复通道，按原 Re/Im 顺序实化；每实通道噪声标准差是固定 solid 参考响应实化 RMS 的 0.12 倍，复数总方差为 $2\sigma^2$。不能对每个候选重设噪声或把系数均匀先验换成半径均匀。[P5]

### 1.2 两项术语修正

**MR-FFT-TSOM 的 MR 是 multiplicative regularization，不是 multiresolution。** 多分辨率 SOM/IMSA-SOM 是另一条文献线；其中二维 PEC 的重建结果不能直接用于可穿透介电材料。[R3,R5]

C1-P1 的“二维”主要指两个连续材料参数。P0 使用矢量偶极模型及强度读数；这不是一般二维 TM 材料图像 posterior，更不是完整三维材料成像。[P6]

---

## 2. 统一目标：可行集合，不是假装已知的 posterior

设固定 acquisition 与 physics identity 为 $a$，材料为 $\chi$，nuisance 为 $\nu$：

$$

y=F_a(\chi,\nu)+n.

$$

将复数据按冻结顺序实化，以真实噪声协方差白化：

$$

f(\chi,\nu)=W\,{\cal R}F_a(\chi,\nu),\qquad
y_w=W\,{\cal R}y,\qquad n_w\sim{\cal N}(0,I_m).

$$


定义

$$

{\cal C}_{\tau}(y)=
\left\{\chi\in{\cal X}_{\rm phys}\cap{\cal X}_{\rm support}:
\exists\nu\in{\cal N},\
\|f(\chi,\nu)-y_w\|_2\le\tau\right\}.

$$

${\cal X}_{\rm support}$ 是事前声明的材料支持，不等同于已经校准的 prior。可额外声明 $R(\chi)\le R_{\max}$，但必须说明它删掉了哪些物理世界。

数据相容、材料可实现、结构显著不同、prior 合理、高 posterior density、高区域 posterior mass、低真值 MAE，是七个不同命题：

- residual 小只约束预测；
- 被动性和本构一致性只约束合法性；
- diversity 衡量候选之间差别；
- prior 是外加知识或建模假设；
- density 依赖参数化，单个连续图像通常无点概率；
- mass 需要对区域积分并保留区域内变化；
- truth MAE 只在独立评分时可知。

**candidate-only 是完整且诚实的交付类型。** 不把 $\exp(-r_k^2/2)$ 对自适应搜索端点归一化后称为 posterior weights：搜索密度、参数体积、prior、重复表示及漏掉的区域均可能未处理。

### 2.1 残差阈值与二元可辨性必须分开

在完全正确的固定 Gaussian 模型下，真实材料和真实 nuisance 的残差平方服从 $\chi^2_m$。可以预注册 $\tau^2=\chi^2_{m,0.95}$ 作为模型内相容性阈值；这不是估计后残差的通用自由度公式，更不是 posterior credible set。

不能机械使用 $m-p$：约束、正则、非线性、nuisance profile 和病态有效维数都会改变解释。存在离散误差时另设 numerical/model discrepancy 预算；不允许发现补偿失败后调大噪声。

两个候选的均值距离

$$

d_{ij}=\|f(\chi_i,\nu_i)-f(\chi_j,\nu_j)\|

$$

回答另一个问题。二者都在 ${\cal C}_\tau$ 内，不代表二者难区分。数据维数大时 $\tau$ 可能很大。

候选对还需注明 nuisance 策略：固定仪器参数、各自受限 profile，或将共同 proper nuisance prior 积分成混合分布。三者不能交换。

---

## 3. 从连续 Maxwell 到离散材料可实现性

### 3.1 连续形式

取 $e^{-i\omega t}$ 约定，已知磁导率、源、边界和观测算子：

$$

{\cal M}_\epsilon E_s
=\left(\nabla\times\mu^{-1}\nabla\times-\omega^2\epsilon\right)E_s=b_s.

$$

在两个材料的散射问题均适定、相关算子域和边界保持一致时，

$$

{\cal M}_{\epsilon_1}(E_{1,s}-E_{0,s})
=\omega^2(\epsilon_1-\epsilon_0)E_{0,s}.
\tag{1}

$$

这是有限差值关系，不是 Born 截断。连续形式要在相应 $H(\mathrm{curl})$ 或 VIE 空间中解释，要求乘法和求解算子有界；不能从离散矩阵的一个条件数推广到所有连续介质。

对归一化 contrast current，写

$$

j_s=\chi E_s,\qquad E_s=E_s^{\rm inc}+Gj_s.

$$

任何候选必须同时满足内部场关系、同一材料本构关系、所有照明共享 $\chi$，以及材料支持、互易性、被动性等条件。

被动性约束作用于**绝对介电常数**：标量为 $\operatorname{Im}\epsilon\ge0$，张量为 $\operatorname{Im}_H\epsilon\succeq0$。有耗参考背景下 contrast 的虚部允许为负。多频需共同色散模型或明确有限带常数近似；逐频随意取值不等于全频因果材料。[P2]

### 3.2 离散 VIE / DDA 的共同形式

令

$$

j_s=A(\chi)(e_s^{\rm inc}+Gj_s),\quad
L_\chi=I-A(\chi)G.

$$

VIE 的 $A$ 可以线性依赖 $\chi$；DDA 的 $A$ 往往是 CM/radiation-reaction 极化率的块对角矩阵。因此 $\Delta A\neq DA[\Delta\chi]$ 一般成立，有限编辑必须计算**真实极化率差**。[P1,P2]

#### 命题 1：标量材料的逐单元、多照明可实现性

给定所有候选电流 $j'_s$，先算 $E'_s=e_s^{\rm inc}+Gj'_s$。将单元 $i$ 的全部三分量及全部照明堆叠为 $u_i,v_i$，分别表示候选电流与总场。

若 $\|v_i\|>0$，存在一个共同标量极化率 $a_i$ 当且仅当

$$

\left(I-\frac{v_iv_i^*}{\|v_i\|^2}\right)u_i=0,
\quad
a_i=\frac{v_i^*u_i}{\|v_i\|^2},
\tag{2}

$$

并且这个 $a_i$ 属于合法材料通过本构映射得到的集合，能够映回允许的 $\chi_i$。

**证明**：合法性要求 $u_i=a_iv_i$。对 $v_i$ 的正交分解直接给出必要性；投影残差为零时，式(2)构造的 $a_i$ 使本构关系成立。逐单元成立，加上已经构造的内部场关系，等价于全部离散 Maxwell 方程。再检查材料域即可。这里“充分”针对所声明离散模型，不是网格收敛证明。

若 $v_i=0$，必须 $u_i=0$；此时当前照明不能确定该单元材料。若总场很小，材料估计对电流误差至少按 $1/\|v_i\|$ 放大。多个照明通过 $\sum_s\|E_{i,s}\|^2$ 缓和问题，但不保证消除所有场节点。

对于未知张量材料，写 $J_i=A_iE_i$，其中矩阵列对应照明。纯线性相容条件为

$$

J_i(I-E_i^\dagger E_i)=0.

$$

还需存在满足对称性、被动性和参数支持的 $A_i$。张量的自由度更多，不能悄悄把标量重建失败改成张量重建成功。

#### 命题 2：合法 material tangent 的精确条件

在参考材料处，合法材料方向 $h$ 产生

$$

L_0\,\delta j_s=B_sh,\qquad
B_sh=DA(\chi_0)[h]E_{0,s}.
\tag{3}

$$

因此 current perturbation 合法的条件是：存在**同一个**实材料向量 $h$，属于当前物理约束的切空间或切锥，使式(3)对全部源成立。

定义堆叠的 ${\cal K}h=[L_0^{-1}B_sh]_s$。则

$$

\text{material-realizable current tangents}=\operatorname{Range}{\cal K},

$$


$$

\text{receiver-dark tangents}=\ker {\cal S},
\qquad
\text{material-null directions}=\ker({\cal S}{\cal K}).
\tag{4}

$$

在材料合法内点、允许全部线性方向时，可实现暗方向是 $\operatorname{Range}{\cal K}\cap\ker{\cal S}$，而不是整个 $\ker{\cal S}$。在约束边界，应将 $\operatorname{Range}{\cal K}$ 改成切锥的像 ${\cal K}(T_{{\cal X}_{\rm phys}}(\chi_0))$，它一般不再是线性子空间。

**反例**：两个单元，接收器只看第一个单元。第二单元的两个照明总场为 1、2。请求第二单元电流增量均为 1，确实接收器不可见，但分别要求材料增量 1、1/2，不存在共同材料。这是逻辑反例，不把这个简化接收器当作真实阵列设计。

### 3.3 一项容易写错的“feedback 变可见”

若 $\Delta j=j_1-j_0$ 已经是两个**完整精确电流**之差且 $S\Delta j=0$，它就仍然不可见，不能声称之后 feedback 又使同一个 $\Delta j$ 可见。

会变可见的是：一个原始注入方向、冻结的一阶方向、或 basis atom $C$，经过物理反馈后形成另一个方向。例如 A17 中 $SC=0$ 不阻止

$$

S\bar C=S(I-R_0L)C\ne0.

$$

两者不是同一个完整电流扰动。[P1]

### 3.4 TSOM 能提供什么

SOM/TSOM 的可观测电流、ambiguous completion、内部传播约束，天然适合作为**提出或预条件材料补偿方向的工具**；不能直接将 arbitrary ambiguous current division 变成多解算法。

可以使用

$$

\min_{h,\{c_s\}}\sum_s
\|L_0U_{\rm amb}c_s-B_sh\|^2
+\lambda\|h\|_{M_\chi}^2,

$$

同时约束 $h$ 的材料合法性和散射变化。消去 $c_s$ 后，本质上又回到材料空间的约束问题。若比直接 $J_\chi h$ 方法更贵，current latent 应删除。

第二重 subspace 扩张主要增强内部状态表示和一致性，不自动产生分支探索。只有加入多起点、拓扑编辑、排斥或延续，才成为多解搜索机制。基底旋转 $U\to UQ$、相位 gauge、重复参数化均不是新材料。

---

## 4. 局部材料解流形：成立条件和失败边界

### 4.1 常秩命题

令 $\theta=(\chi,\nu)\in\mathbb R^d$，$f$ 为 $C^1$。如果在内点 $\theta_0$ 的邻域内 $\operatorname{rank}Df=r$ 恒定，则

$$

f^{-1}(f(\theta_0))

$$

在该邻域是维数 $d-r$ 的嵌入流形，切空间为 $\ker Df(\theta_0)$。这由常秩定理得到，是标准数学，不应宣称新的 Maxwell 非唯一性定理。

需要注意：$d-r$ 可能主要来自 nuisance 或重复 latent，不能直接称为不同材料的维数。材料投影后的维数还需扣除纯 nuisance/gauge 方向。物理边界会把流形截成带边界或分层可行集。

连续材料空间需 Banach/Hilbert 版本的 split-surjectivity、闭值域和补空间等条件，不能仅用一个离散 SVD 证明连续结论。

### 4.2 生成变量与补偿变量分开

写 $\chi=\chi(c,z)$，$z$ 表示编辑或弱结构，$c$ 为可补偿坐标。若选定 $r$ 个局部独立测量约束 $f_s$，且 $\partial_c f_s$ 可逆，则

$$

c=\psi(z),\qquad
D\psi(z)=-(\partial_cf_s)^{-1}\partial_zf_s.
\tag{5}

$$


生成器只需提出 $z$，物理 corrector 解出 $c$。但只保持“强通道”并不自动保持完整数据。只有这些约束局部刻画完整等响应集合，才得到完整等响应；否则仍需核验弱通道 residual。

这也说明 A23 的 recoil 应修改稳定材料坐标以抵消弱结构的二阶响应，不能把强坐标永久固定。[P2]

### 4.3 Nuisance 的处理

白化线性化：

$$

J_\chi h+J_\nu v\approx0.

$$

无约束 nuisance 下，消去 $v$ 得

$$

J_{\rm eff}=(I-P_{\operatorname{Range}J_\nu})J_\chi.
\tag{6}

$$

对 proper Gaussian nuisance 正则，得到的是 Schur 权重

$$

I-J_\nu(J_\nu^TJ_\nu+\Lambda_\nu)^{-1}J_\nu^T,

$$

一般不是正交投影。对有界 gain、已知标定区间，应直接做受约束子问题；式(6)至多提供乐观的线性可补偿下界。

### 4.4 零切向为什么不够

$f(a,b)=(a,b^2)$ 在原点的 $b$ 方向属于 Jacobian nullspace，但 $f=(0,0)$ 的局部解只有原点。这里 rank 不恒定。

同一模型对观测 $f=(0,1)$ 有 $b=\pm1$ 两个不连通分支；在一个分支内做足够小的 continuation 不能发现另一个。局部 Fisher/Hessian 也无法告诉你另一个分支在哪里。

有限噪声下应沿 feasible tube，而不是强制精确插值 noisy data。可采用预测—校正：

$$

\theta_{\rm pred}=\theta_t+\alpha v_t,\quad
\theta_{t+1}=\arg\min_{\theta\in{\cal X}}
\|\theta-\theta_{\rm pred}\|^2_M
\quad\text{s.t. }\|f(\theta)-y_w\|\le\tau.

$$

实际解为局部受约束优化，不保证全局投影或唯一 retraction。步长由 residual 预测误差和合法性约束调整；fold 附近换坐标或用 pseudo-arclength，但 rank drop 和 disconnected components 仍需单独处理。

### 4.5 拓扑改变不是普通边界形变

保持光滑可逆形变的 shape chart 保持拓扑，不能出生新孔洞。密度/level-set 可穿过零水平集临界点改变拓扑，但该点通常不属于同一个光滑 shape chart。数字图像拓扑还依赖阈值、网格和最小特征尺度。

因此采用两级搜索：**离散或非局部 structural edit 跳转**到新族，再**连续 Maxwell compensation**恢复数据相容。孔洞 birth/death、分裂合并、移动、层析改变，应由预注册编辑族覆盖，而不是依靠一个局部 chart 自动完成。

---

## 5. 有限材料变化：真正应复用的 Maxwell 关系

### 5.1 精确相对反馈

令 $A_1=A_0+\Delta A$，固定 $G,S,e^{\rm inc}$，且 $L_0,L_1$ 可逆。由两状态相减：

$$

(L_0-\Delta AG)\Delta j_s=\Delta A E_{0,s}.

$$

记

$$

R_0=L_0^{-1},\quad X_s=R_0\Delta AE_{0,s},
\quad T=R_0\Delta AG,

$$

则

$$

\boxed{\Delta j_s=(I-T)^{-1}X_s.}
\tag{7}

$$

式(7)只需相关逆存在，不要求 $\|T\|<1$。但使用 Neumann 多项式时需要更强条件；若 $\eta=\|T\|<1$，

$$

\|\Delta y_s-SX_s-STX_s\|
\le \|S\|\frac{\eta^2}{1-\eta}\|X_s\|.
\tag{8}

$$

绝对背景强散射但相对编辑弱，可以有小 $\eta$；大孔洞、近接触或共振附近通常不能默认如此。谱半径小也不控制非正规有限阶放大。

若整条材料路径合法且 $F$ 可微，

$$

F(\chi_1)-F(\chi_0)=
\int_0^1J_\chi(\chi_0+t(\chi_1-\chi_0))(\chi_1-\chi_0)\,dt.
\tag{9}

$$

式(9)说明，有限编辑取决于路径上旋转的 tangent，而不是基础解处一张 SVD。

### 5.2 命题 3：局部编辑的有限幅度 Schur 实现

在固定网格上，若

$$

\Delta A=UDV^*,

$$

其中 $r$ 是修改支撑对应的电流秩，定义

$$

Q_0=R_0U,\quad K_0=V^*GQ_0,\quad b_{0,s}=V^*E_{0,s}.

$$

则

$$

\boxed{\Delta j_s=
Q_0(I-DK_0)^{-1}D b_{0,s},\qquad
\Delta y_s=SQ_0(I-DK_0)^{-1}D b_{0,s}.}
\tag{10}

$$


**证明**：式(7)的增量由 $\Delta A$ 注入，故写成 $\Delta j_s=Q_0c_s$。代回差值方程得到
$c_s=D b_{0,s}+DK_0c_s$，即式(10)。由 determinant identity，若 $L_0$ 可逆，核心 $I-DK_0$ 可逆等价于对应 $L_1$ 可逆。

对于 $q$ 个各向同性 DDA 单元，取 $U=V$ 为三分量选择矩阵，$r\le3q$。公式不需要 $D^{-1}$，因此包含将某些单元移除为真空、contrast 变零等编辑。所有照明共享核心，只更换右端项。

**这是 Woodbury/Schur 消元的物理化使用，不是声称首次发明低秩更新。** 潜在新意在可实现结构编辑、补偿位置选择、多候选共用、覆盖与实测成本的组合。

四个重要限制：

- 小材料参数维数不等于小 $r$。一个全域介电常数参数可改变全部对角块，矩阵增量满秩。
- 更新支撑扩散至大量区域，或为补偿加入全域 contrast，精确小核优势会消失。
- 改网格、源、接收器、磁导率、自项约定时，不能继续假定 $G,S$ 不变。
- 获得 $Q_0=R_0U$ 需要真实求解动作；一次参考场求解不等于免费得到了任意 $R_0$ 作用。

### 5.3 命题 4：有限编辑的局部可补偿性夹逼

在编辑后材料 $\chi_e$ 处，令完整白化 residual 为 $r_e$。允许补偿坐标 $c$ 的 Jacobian 为

$$

A_e=[J_\chi C_\chi,\ J_\nu C_\nu],\qquad
u_e=(I-A_eA_e^\dagger)r_e.

$$

在受限球 $\|c\|\le\rho$ 内，假设

$$

f_e(c)-y_w=r_e+A_ec+R(c),
\quad \|R(c)\|\le L_c\|c\|^2/2.
\tag{11}

$$

这里 $c$ 的尺度须冻结，包含材料和 nuisance 单位归一化。

对任何允许 $c$，

$$

\|f_e(c)-y_w\|\ge\|u_e\|-L_c\rho^2/2.
\tag{12}

$$

因此若右侧大于 $\tau$，该**局部补偿域**内无解。

若 $c_*=-A_e^\dagger r_e$ 满足全部约束、且

$$

\|u_e\|+L_c\|c_*\|^2/2\le\tau,
\tag{13}

$$

则它给出一个可行候选。含 forward 数值误差时在(12)减去、(13)加上相应保守误差界。

证明只用投影消去线性补偿项与三角不等式。真实约束会缩小可补偿集合，故受约束最小二乘应优先；无约束投影可以作乐观诊断，不能据此接受超界解。

**若 $L_c$ 和误差半径只有经验估计，这就是 screening indicator，不是 certificate。** 更不能从局部无解宣布所有远方世界不存在。

此命题将 A30 与 Joint P2 放进同一框架：结构编辑 secant 是否落在允许材料/nuisance 补偿响应内，才是机制问题。某个空腔可补偿，不代表某个 ring 可补偿。

### 5.4 A17 在这里的准确位置


$$

\Delta J=Y\Gamma^{-1}N

$$

能解释并修复**给定状态下**的 current-ROM tangent 缺失，帮助选择增广方向、构造低秩 correction、诊断 receiver-side 与 injection-side 的双侧误差。

它不能独自预测一个大编辑后的 $F(\chi_1)$，因为材料改变同时修改 $A,L,E,B$。跨材料有限响应应使用(7)/(10)或重新求解，不能将 A17 的 frozen model-change 定理改名为 material-world equivalence。

---

## 6. 本轮独立小型验证

`validation/verify_algebra.py` 使用八个矢量电偶极、六个平面波照明、六个远场标量通道，固定网格，CM 加 radiation-reaction 极化率。修改包括一个单元变为真空、另一个单元改变 contrast。没有读取上游真值、没有训练或运行成像实验。

| 检查 | 实测 |
|---|---:|
| 精确 relative feedback 相对误差 | $6.26\times10^{-16}$ |
| 局部 Schur 完整电流相对误差 | $6.41\times10^{-16}$ |
| 局部 Schur 输出相对误差 | $1.10\times10^{-15}$ |
| 合法电流多照明材料一致性误差 | $2.40\times10^{-16}$ |
| 构造的暗电流接收响应相对误差 | $2.09\times10^{-16}$ |
| 同一暗电流的有限材料一致性缺陷 | 18.75% |
| 同一暗扰动的共同材料 tangent 缺陷 | 94.75% |
| 双侧残差校正恒等式误差 | $7.69\times10^{-16}$ |
| 一阶响应误差 log-log 斜率 | 1.99989 |
| 二阶响应误差 log-log 斜率 | 3.00055 |
| 单个全域材料参数的矩阵更新秩 | 24/24 |

$\operatorname{cond}(L_0)\approx1.123$，$\|T\|\approx0.03786$。这是适定小矩阵示例，**不支持强散射鲁棒性、空腔成像成功、算法 speedup 或全部材料覆盖**。完整代码、参数、断言、版本和输出已交付。


## 7. 八类架构比较：哪些是真正不同的选择

| 路线 | 真正的搜索变量与机制 | 主要优点 | 主要风险 | 本轮排序 |
|---|---|---|---|---|
| A：TSOM + ambiguous-current completion | 每源电流系数，再投影到共同材料 | 可复用现有 current pipeline | 不合法电流比例可能高；多源本构修复使采样优势消失 | 诊断对照，不作主线 |
| B：material-null continuation | 合法材料切向 + 非线性 corrector | 局部分支探索直接、物理明确 | rank 旋转、边界、断开的族；弱维数未必小 | 主算法的局部扩展 |
| C：结构编辑 + 材料/nuisance 补偿 | 有限编辑变量及可补偿坐标 | 可跨孔洞/连通性；易保留明确结构区别 | 编辑族偏置，优化可能难；需要 full response | **首选主线** |
| D：material-space deflation / multistart | inverse objective/KKT，多盆地搜索 | 不依赖一个 anchor；简单强基线 | 排斥可能只找同族或破坏拟合 | 核心基线、备选部件 |
| E：source/receiver adaptive ROM | 改变候选评估器，不一定改变搜索策略 | 可能摊薄大量近邻评价成本 | setup、组装、刷新、误拒绝；不是新世界生成器 | 条件开启 |
| F：conditional learned compensation | 学习 edit、弱变量或 corrector 初值 | 可以降低学习输出维数 | teacher 成本；低维 latent 排除材料；没有自动概率语义 | deterministic 成功后再做 |
| G：hypothesis distillation + correction | 从高质量候选学习 proposal | 把已有搜索经验摊销到多 query | 候选搜索偏差被蒸馏；无全模式保证 | 最低成本学习分支 |
| H：wavelet / sparse / multiscale completion | 生成高频/弱可观测材料系数 | 不强迫固定拓扑；可接 A30 容量证据 | 固定 support 仍会漏；稀疏度不等于 weak rank | 主算法的通用参数表示 |
| I：局部 exact Green/Schur 更新 | 求解有限局部材料改动的完整离散响应 | 可避免小扰动截断；多源共享核心 | 大支撑失去低秩，预计算不免费 | **优先于通用 ROM 试验** |

B+H 往往是同一 material-space 方法采用不同基底；F 与 G 是 proposal 来源不同，不是两种新物理。E 和 I 是可搭配各搜索器的 forward engine，不应把其加速全算成某一生成器的优势。A 经 material consistency 消元后，可能只是更昂贵的 B。D 中“找到不同逆解的 nonlinear deflation”和线性 Krylov deflation 不是一回事。[R11]

### 7.1 能否避免完整高维 posterior 学习

能。可以输出

$$

z\sim q(z\mid\chi_0,y),\qquad
\chi={\cal P}_{y,{\cal X}}\bigl(\operatorname{Edit}(\chi_0,z)\bigr),

$$

其中 ${\cal P}$ 是物理 corrector，不要求 $q$ 是 posterior。

但少学几个坐标并不自动保留全部材料自由度。若合法局部材料纤维维数为 $d_f$，一个光滑、秩不超过 $d_z<d_f$ 的 decoder 不能覆盖其开放子集；它只能取一个较薄的选择。若 latent 支撑连通、decoder 连续且严格保持可行，输出也连通，不能同时覆盖互不连通的完整分支。混合离散 edit labels、多 chart、全域随机跳转或允许提议暂时不可行再纠正，可避免这项结构限制。

对于本目标，合理 latent 是“编辑类型、位置、尺度、幅度、多尺度残差种子”，而不是假定整个高维补偿空间必然只有 4 或 8 维。

### 7.2 MDN、Flow 与物理优化的分工

MDN 适合低维、多簇的 edit/compensation 初值；有界 variable transform 和混合成分易审计。Flow 适合确实弯曲、连续且维度适中的补偿分布；需要明显多 query 摊销时才有意义。离散孔洞出生/消失应显式建模混合变量，不能指望一个可逆连续 Flow 自行改变支撑拓扑。

强可观测坐标的精确修正、nuisance 约束、被动性和全波一致性，优先交给优化器。网络可预测初值或预条件器，不应替代 acceptance test。

物理投影后的分布通常不再具有原 proposal 的解析密度。除非计算相应 Jacobian、coarea 权重或用合法 MCMC/IS correction，交付 reliability 仍为 candidate-only 或 learned-proposal，而非 corrected posterior。[P2,R6–R9]

---

## 8. 推荐主算法：有限编辑的受约束补偿集合

### 8.1 基础解与支持

用当前最便宜、在共同任务上达到质量标准的 TSOM、SOM、GN、已有 N3 或 ROM 得到 $\chi_0,\nu_0$。保留初始化失败；不能换一个容易的 truth 继续测试。

基础解不是“真实可观测部分”的无误读出。允许它在后续修正中改变。若其 full residual 不合格，先修复或重新初始化；不可把不合格基础解算作合法 $K=1$。

材料表示使用完整合法域上的 multiscale coordinates，允许按预注册规则扩张 wavelet support。第一轮保留 16/32/64 活跃系数，但不永久锁死在一个训练 PCA chart；记录所有已经冻结为零的自由度。

### 8.2 编辑字典

按固定顺序轮询：中心/偏心空腔、空腔移动与缩放、多孔 birth/death、材料岛 split/merge、内部分层、边界低频形变、局部 contrast 对换、高对比块和损耗再分配。位置不依赖 truth，可来自均匀覆盖、当前材料 support 和测量相关的有限 secant 分数。

保留约 20% 的 proposal 预算给非贪心随机/全域编辑，避免 response score 只探索最容易的同一族。这个比例是拟议的设计参数，应在 validation 冻结，不是理论最优。

需要结构保持约束：例如要求新空腔至少达到预注册体积与 persistence，优化时不能直接把它填回原样。必须把这个临时子域公开为搜索约束；其他编辑臂仍允许 hole death，避免整体支持被偷偷缩小。

### 8.3 先可行，后 diversity

对编辑 $e$ 求

$$

\min_{c,\nu}\
\frac12\|f(\chi(e,c),\nu)-y_w\|^2+\lambda R_{\rm soft}(\chi)

$$

并受材料 box、被动性、仪器界、最小特征和编辑保持约束限制。

先到 $\|r\|\le\tau-\varepsilon_{\rm num}$ 的相容区，再沿可行方向增加

$$

D(\chi;{\cal H})=\min_{\psi\in{\cal H}}d_{\rm struct}(\chi,\psi).

$$

不能用一个可任意调权的 residual–diversity 加权和，允许“大多样性补偿很坏的数据拟合”。更适合 lexicographic / constrained Pareto 形式：硬保持合法性和残差，软优化结构覆盖与 prior。

$d_{\rm struct}$ 同时包括 volume material distance、support、边界、孔洞位置/大小和 persistence。尺度与权重只在训练/validation 冻结，不根据 test 真值调参。

### 8.4 Corrector 与 engine

优先 LM/SQP/trust-region + material adjoint，必要时 relinearize。每次更新共用所有照明的材料。实际/预测 reduction 不匹配、core 变病态或白化响应误差超 margin 时缩步、切 full engine 或 refresh。

小支撑编辑先试式(10)的 exact local engine；大支撑使用原 full solver。只有 engine profiling 显示 setup 可摊销，才加入 approximate ROM。对 ROM 预测会被拒绝的点进行预注册随机 full audit，避免“只核验被接受点”掩盖漏掉真解的错误。

### 8.5 多分支 archive

保存全部已核验候选及其编辑父节点、corrector 路径和失败原因。连续路径经检查才能标记为“存在已验证连接”；仅仅两材料场接近或预测接近，不能证明分支连接，也不能证明断开。

使用固定 beam 宽度和每族配额，而不是发现第一条容易路径就耗尽预算。弱方向 continuation 负责扩展已发现分支；全域编辑负责发现远处族。

### 8.6 固定 K 的展示集

从所有合法候选建立嵌套 $K=1,2,4,8$ 集合，primary 为 $K=4$。保留合法基础解，然后按固定 farthest-point/medoid 规则增加代表。该规则只看允许的测量、材料自身描述和预声明 prior，不看 truth。

不足四个物理不同候选时输出真实 $K_{\rm valid}<4$，不复制同一图凑数。比较时显示 requested K、valid K 和 distinct-family K。对不同算法，proposal 数、全波/RHS预算及候选筛选规则同时匹配。

### 8.7 不以输出多张图冒充概率

同一优化 basin 找出的十个端点不能贡献十倍 mass。候选图、conditional samples、proper weights 和 unexplored remainder 分栏保存。仅有 candidates 时，把质量解释限制为“这些材料尚未被本测量排除”。

---

## 9. 结构不同的备选：population material search

建立不围绕一个 $\chi_0$ 的粒子群，从共同材料 support 的多尺度随机初始化和拓扑编辑出发。使用 tempered residual level 或 proper posterior tempering，配合局部材料 GN/MALA/随机游走及全域 birth/death move。

与主线的区别是：它不假设其他重要材料世界可由当前基础解的少量编辑到达。代价是更分散的全波预算、混合困难和重复收敛。

可以有两个明确实现：

- **经验性集合搜索**：population + 不同材料簇保持 + full validation，输出 candidate-only；不要求 posterior 正确。
- **合法 SMC/MCMC**：指定 proper prior、温度、proposal、birth/death 维度变换 Jacobian、接受概率；importance correction 使用真实 proposal density，ESS 与多模式访问分别报告。

对连续解流形，不建议直接给全部 noisy residual 乘 singular deflation 因子。通常应对 inverse objective 的局部驻点/KKT 使用 deflation 或加 material-space exclusion constraints，然后重新核验原始残差。固定材料下的线性 Maxwell 场唯一性，不是这里要 deflate 的“多解”。[R11]

在 P2 四维空间，参考积分、有限网格/低阶代理和传统多起点应是强基线。先击败它们再谈 MDN/Flow；不能用新网络得到更多监督后，宣称架构本身超过普通 sampling。[P5,P6]

---

## 10. ROM 省算审计与可靠性

### 10.1 用途排序

| 用途 | 所需正确性 | 本研究建议 |
|---|---|---|
| 基础重建 | 达到共同图像质量和合法残差 | TSOM/ROM 必须和 full GN 实测比较 |
| proposal | 允许不精确，但需保留探索机会 | 最宽松，可用 Born/局部近似 |
| screening | 错误排除重要候选会破坏 coverage | 需边界 margin、随机拒绝点核验 |
| 补偿优化 | 方向要足以到达真实可行域 | trust region + full fallback |
| 正式 likelihood / posterior | 对全部保留区域有受控 log-likelihood 误差 | 本轮不默认 ROM 合格 |

### 10.2 双侧残差校正

对 $Lj=b$，近似电流 $\hat j$、近似 adjoint $\hat Z$，令

$$

r_p=b-L\hat j,\quad r_d=S^*-L^*\hat Z,\quad
\hat y_c=S\hat j+\hat Z^*r_p.

$$

则

$$

\boxed{y-\hat y_c=r_d^*L^{-1}r_p.}
\tag{14}

$$

所以

$$

\|y-\hat y_c\|\le\|r_d\|\|L^{-1}\|\|r_p\|.

$$

这是精确代数恒等式及条件上界。若没有可靠 $\|L^{-1}\|$ 或 inf-sup 下界，小 primal residual、本构被动性、有限点核验均不形成 uniform certificate。Maxwell 非正规性和近共振尤其需要谨慎。[P1,P2,R10]

若白化 forward 误差不超过 $\varepsilon$，Gaussian likelihood 满足

$$

|\log L-\log\hat L|
\le \|\hat f-y_w\|\varepsilon+\varepsilon^2/2.
\tag{15}

$$

因此“1% 相对场误差”不一定适合 posterior。高维、小噪声条件下，即使适合 proposal 的误差也可能严重影响 likelihood ranking。

### 10.3 成本公式

设全电流维数 $n$，源数 $s$，活跃材料维数 $p$，ROM 秩 $r$，局部修改支撑 $q$，提议数 $P$，合法输出数 $K$：

$$

C_{\rm total}=C_{\rm base}+C_{\rm setup}
+\sum_{i=1}^{P}(C_{{\rm proposal},i}+C_{{\rm corrector},i})
+C_{\rm refresh}+C_{\rm audit}+C_{\rm final}+C_{\rm select}.
\tag{16}

$$

失败、拒绝和中断都进入分母。

稠密直接求解一般是 $O(n^3+n^2s)$；有 FFT 的迭代 VIE 才可讨论每次 matvec 的 $O(n\log n)$，实际迭代数和预条件费用仍须付。不能给未实现 FFT 的 DDA 后端记上 FFT 复杂度。

ROM 求解 $O(r^3+r^2s)$ 不包含算子组装。一般材料变化下，投影组装可能需要 $O(nr^2)$，或存储 $p$ 个仿射块、付 $O(pr^2)$ 在线组合。Material-Jacobian 的全部弱方向搜索本身也不免费；优先 matrix-free JVP/VJP 与有限随机 sketch，并记录 sketch 漏方向风险。

式(10)的 $Q_0$ 要付 $r\le3q$ 个 anchor RHS；多个源共享小核。全局一参数修改仍可能有 $r=n$。缓存完整逆或大 Green 块的内存与 setup 要明示。

适用的摊销条件是

$$

P(C_{\rm full}-C_{\rm reduced})
>
C_{\rm setup}+C_{\rm refresh}+C_{\rm extra\ failures}+C_{\rm audit},
\tag{17}

$$

在相同成功率、验证规则和 output K 下比较。每个候选都需最终 full check 时，那部分成本不能藏起来。

### 10.4 精确局部响应缓存的已有先例

2026 年 Sun 等的 PNGF 工作已经采用预计算 numerical Green function 和 Woodbury 更新进行全波局部逆设计；同年可更新 multiport connection 工作也复用低秩更新。[R17,R18]

这显著收窄本研究的原创性：**局部变化、一次环境预计算、多个候选低秩更新，本身已有电磁先例。** Research A 应借用并比较这些思想，不应把式(10)包装为新的电磁求解范式。

同时，这些局部设计框架与未知材料成像的条件不同：静态环境和可变区域必须可靠；未知材料可能分布在全域。是否仍有小可变区域、怎样搜索而不漏重要世界，是本项目需要实测的部分。

### 10.5 何时立即放弃 ROM

若在预注册工作量 $P,Q,K$ 下总 cold cost 不下降，或只在未实测的 $Q=100$ 摊销投影下获胜；若刷新频繁到接近每个候选一次 full solve；若正确 screening 需要的误差控制比 full evaluation 更贵；若同质量图像尚未达到——直接停 ROM，保留材料补偿搜索及原 full solver。

原始 VIE/DDA 加缓存、直接复用因子、block RHS、局部精确更新，可能比另建一个低秩 current ROM 更自然。


## 11. 以材料 MAE、覆盖和选择代价为主的评价

### 11.1 材料指标

冻结体素体积 $v_i$、评价域 $D$ 和材料动态范围 $\Delta_R,\Delta_I>0$：

$$

E_R(\chi,\chi^*)=\frac{\sum_i v_i|\epsilon'_i-\epsilon_i^{*\prime}|}
{\sum_i v_i},\quad
E_I(\chi,\chi^*)=\frac{\sum_i v_i|\epsilon''_i-\epsilon_i^{*\prime\prime}|}
{\sum_i v_i},

$$


$$

E_{\rm norm}=\frac12(E_R/\Delta_R+E_I/\Delta_I).
\tag{18}

$$

同时报告非归一化 $E_R,E_I$ 及复模长 MAE；不能只报平均值掩盖虚部失败。若只恢复实部，不给虚部虚假的零误差成绩。径向材料采用 $4\pi r^2dr$ 体积积分，不能使用四个 latent 的普通 MAE。

评价域固定，不随候选缩小；同时提供固定对象包围域/结构 ROI 的误差，防止大背景稀释。IoU、boundary distance、空腔体积/位置/数量误差和多阈值拓扑 persistence 单独报告。

### 11.2 三种误差不能互换


$$

E_{\rm Best@K}=\min_{k\le K}E(\chi_k,\chi^*),

$$

仅供模拟评测，禁止用于在线排序。


$$

E_{\rm Selected@K}=E(\chi_{\hat k(y,\mathrm{prior},\mathrm{context})},\chi^*),

$$

selector 只使用实际可得信息。至少保留“最低 full residual”“proper prior + likelihood 的离散 MAP 候选”“存在参考时的候选条件风险最小化”三种清楚命名的选择器。

$E_{\rm Top1}$ 是基础 TSOM/GN/MAP 等独立重建的误差，不等于从候选池事后用 truth 挑出的最优图。MAE 的 Bayes 最优决策通常与 MAP 不同；低 residual 也不是 MAE 最优选择器。


$$

{\rm Coverage@K}(\varepsilon)
=\Pr[\min_k E(\chi_k,\chi^*)\le\varepsilon,\
\text{相关结构误差也合格}],
\tag{19}

$$


$$

{\rm Regret@K}=E_{\rm Selected@K}-E_{\rm Best@K}.

$$

报告 residual、合法候选率、结构 diversity、有效不同材料族数、每新增合法解的累计成本、全波核验失败率，以及先验/噪声/gain/基础解敏感性。

base 合法且保留时，Best-of-K 在该嵌套集合中不会比 base 差。这只是集合包含关系，不是新算法优势；必须和同预算 multi-start 的 Best-of-K 比较。不同 pool 的 selected regret 不能用“oracle 潜力”偷换。

### 11.3 可预注册的“合格 MAE”

目前没有提供具体应用容许的介电常数偏差和最小空腔尺度，不能凭空宣布固定 NMAE=0.05 为工程合格。

建议分成两级：

**研究合格**：按预声明噪声/acquisition 难度层 $h$，在 validation 上取得最强常规基线的材料 MAE 中位数 $b_h$，冻结

$$

\varepsilon_h=0.90\,b_h.
\tag{20}

$$

10% 是本轮建议的最低研究收益，不是物理常数。若数值重建重复产生的材料误差波动 $u_h$ 已接近 $0.10b_h$，先收紧数值误差或宣布该层无法分辨这一收益，不根据 test 调门槛。

**工程合格**：额外给定材料实/虚部允许误差、最小可见结构尺寸、允许边界偏差和 cavity 漏检率。二者同时达到才能称应用合格。若有该阈值，最终 $\varepsilon_h$ 取(20)与应用上限中更严格的一个。

确认试验的拟议 GO：
- 同 full-forward/RHS 预算下，Best-of-4 相对最强多解基线改善至少 10%，按 parent 配对的 95% 区间排除零；
- 在冻结 $\varepsilon_h$ 下，Coverage@4 提升至少 10 个百分点，或以至少 20% 更低总成本达到非劣覆盖和 MAE；
- Selected-MAE 不作不可能的必胜承诺，但必须报告；若恶化，不能以 Best-of-K 掩盖；
- 所有正式交接候选通过规定 full-Maxwell/材料核验；失败尝试全部计费；
- 几何、实部、虚部指标不能出现预注册的严重退化。

这些是拟议门槛，试验前冻结；不表示本研究已达到。小 pilot 只估计效应和失败率，不能靠少量 bootstrap 宣称总体稳定优势。

### 11.4 选择真值的信息下界

#### 命题 5：两个混淆材料的不可避免选择风险

设两材料的 MAE 距离为 $\Delta$，两测量分布为 $P_0,P_1$，等先验，任何基于本次观测的估计器 $\hat\chi(y)$ 满足

$$

\frac12{\mathbb E}_0d(\hat\chi,\chi_0)
+\frac12{\mathbb E}_1d(\hat\chi,\chi_1)
\ge\frac{\Delta}{2}(1-\operatorname{TV}(P_0,P_1)).
\tag{21}

$$

证明：逐 $y$ 下界为
$\tfrac12\min(p_0,p_1)[d(\hat\chi,\chi_0)+d(\hat\chi,\chi_1)]$，用三角不等式，再积分。

若都是相同单位协方差的 Gaussian，均值距离为 $d_y$，最优二元错误率为 $\Phi(-d_y/2)$，故

$$

{\cal R}\ge\Delta\,\Phi(-d_y/2).
\tag{22}

$$

对 P2 的固定补偿端点，$d_y=0.900170$ 给出二元错误率 32.63%；加入第二频率的固定端点距离 2.467944 对应 10.86%。这些是固定两世界、等先验的说明，不是 P2 实际多事件 posterior 的错误率。nuisance 积分后若变为混合分布，应重新计算 TV 或界。

这说明：候选集含两个真实可行世界时 Best-of-2 可以为零，但在线 Selected-MAE 仍严格为正。算法不能从没有的信息中稳定选中真值；Route B 的额外测量才可能缩小该下界。

#### 命题 6：有限 K 的根本容量限制

若存在 $M$ 个两两材料距离超过 $2\varepsilon$、测量分布完全相同的等先验世界，则任意含最多 K 个候选的集合满足

$$

{\rm Coverage@K}(\varepsilon)\le\min(1,K/M).
\tag{23}

$$

证明：每个候选至多覆盖一个世界，而候选集合在各真值下的观测分布相同。

因此不能保证有限 K 包含所有重要解。只采一个局部流形、固定低维 chart、有限 topology 编辑、有限 MCMC、有限 bank、连续低维 generator，都可能漏掉远处世界。模型支持之外的世界更不在任何内部保证内。


### 11.5 “最多能改善多少”：条件集合风险与 posterior quantization

在 P2 这种已有可信低维参考后验的平台，可以比“oracle headroom”更精确地回答这个问题。假设固定观测 $y$ 的参考材料后验为 $\pi(d\chi\mid y)$，合法动作集合 ${\cal C}_\tau(y)$ 非空，且距离具有有限期望，定义

$$

D_K(y)=
\inf_{\substack{H\subset{\cal C}_\tau(y)\\1\le |H|\le K}}
\int\min_{\psi\in H}d(\psi,\chi)\,\pi(d\chi\mid y).
\tag{24}

$$

这是在给定材料距离和合法性约束下的 K 点量化风险，不是新的 Maxwell 定理。它允许算法不学习 posterior，却让已有 posterior reference 提供评价基准。

对任意生成集合 $H(y)$，条件期望 Best-of-K 风险至少为 $D_K(y)$。而任何只根据本次 $y$ 选出一个合法材料的策略，其条件期望损失至少为 $D_1(y)$。因此相对基础解条件风险 $R_0(y)$，理想集合改善至多为 $R_0-D_K$，理想单图选择改善至多为 $R_0-D_1$。二者都是信息与支持相关的条件量，不能从 A30 的旧 35.69% 推出。

对**同一个**候选集合，直接写成

$$

B(H\mid y)=\int\min_k d(\chi_k,\chi)\,\pi(d\chi\mid y),

$$


$$

S(H\mid y)=\min_k\int d(\chi_k,\chi)\,\pi(d\chi\mid y).
\tag{25}

$$

总有 $B(H\mid y)\le S(H\mid y)$。前者允许评测器看到真值后选图；后者是知道正确后验、但看不到真值的最优 candidate selector。差 $S-B$ 是该集合的条件 oracle-selection gap。即使生成器和选择器没有实现错误，这个差也可以严格大于零。

P2 可以利用已有带权材料样本，零新增 Maxwell 地近似计算式(24)–(25)。但有三个限制：

1. 在有限候选 bank 上求 K-medoids 得到的是连续 $D_K$ 的**一个可行上界/最佳已找到值**，不是已认证的全局最优值。没有全局最优证据，就不能将它当成严格的“最多改善多少”数值上限。
2. 参考后验积分、有限候选支持、非凸 K-medoids 优化的误差分开记录。不要为了估算这个辅助量开启大型枚举或新的响应 bank。
3. 使用测试观测的高质量参考后验作 selector 只能标为 **reference-assisted diagnostic**，不进入普通在线 selector 的公平主比较；它不偷读该观测真值，但使用了额外计算和参考信息。

这也解释为什么单纯 constrained diversity maximization 可能浪费 K：最远的点未必能降低期望材料距离。主算法没有可信 posterior 时，可用事先固定的 prior 合理性约束和结构分层覆盖；P2 则可以额外检查其是否接近参考材料距离下的合理 K 点配置。

“多个候选由真值择优”“winner-takes-all/oracle loss”已有 Multiple Choice Learning 先例，不能当作本任务的原创目标。[R19] 本研究改变的是有限 Maxwell 材料集合如何在物理和预算限制下被发现，不是首次定义 Best-of-K。


---

## 12. 最小决定性实验与预算

以下是**下一轮建议，不是本次已经运行的 campaign，也不继承旧项目剩余额度**。所有预算采用“每个材料状态 × 频率尝试”计数；一个状态多源 RHS 另计，已失败尝试也收费。wall 上限是保护开关，不是运行时间预测。

### E0：零新增全波的审计与单元测试

读取并冻结 commits、noise packing、材料支持、prior/gain、旧 bank、teacher、reference、已暴露测试集合。核对 P2 inactive latent 去重及材料场 MAE。对已有响应做 cache identity 检查。

运行本包八偶极恒等式检查，并用仓库自己的最小 backend 做对应 smoke tests；后者也属于新增调用，必须单独记账，不能混进“零新增全波”一栏。

**停止**：发现 whitening、polarizability、自项、source stack、latent prior 或 target identity 不一致，先修复，不发布候选覆盖成绩。已有任务剩余 541 次等数字不构成新授权。

### E1：Joint P2 四维机制与 MAE

物理 support、原四 latent prior、实 gain 的截断正规先验、共享两频 gain、$kR=1.6/2.4$、复数 packing 全部沿用原定义，不偷偷换 radius-uniform prior。

已有 64 个 parent × 2 noise 的 fresh 集现已暴露，标记 regression，不再次称 blind。新的独立 panel 在 generator 和 selector 冻结后，另生成 **32 个独立 parent × 1 noise**；truth 仅由评分进程持有。可复用旧 response bank，但新真值响应、端点和 audit 都计费。

四主臂：有限-edit compensation；全局 multi-start material optimizer；固定小预算的 prior/tempered population；冻结 C1 MDN 加相同 corrector。所有臂可用同一已付费响应代理做 proposal；端点均使用原始 Mie。SMC 若预算太小不足以混合，明确标成 budget-limited sampler，不宣布其统计收敛。

- primary：volume-normalized material Best-of-4 MAE、Coverage@4、Selected-MAE；
- 同时报告 K=1/2/4，K=8 仅预算允许的探索项；
- 每 parent 每臂最多 16 次新增原始 Mie 材料×频率求值，合计 2048；
- 另 256 次共用预算包含新真值、固定前 8 parent 的第二频率负控制及复查；
- **总 cap=2304 次材料×频率尝试，GPU=0，CPU wall 保护上限 1800 秒**；
- 原 teacher、训练和 response bank 的历史成本分栏，不伪装零成本新方法；
- 三个算法 seed 可先在 regression panel 检查稳定性；确认 panel 不能把 seed 当新 parent。

GO：数学/接口正确、真实材料指标有改善信号、核验通过且不靠 latent 重复凑覆盖。这里不要求 TSOM/ROM 加速，也不以四维结果宣称一般三维成像成功。若主法被同预算简单优化全面支配，先停止复杂化，而不是加网络挽救。

### E2：A30 原困难条件的二维反证试验

固定原 24×24 TM、8 源/32 接收器、complex128、原噪声和材料支持。先选原四 family 各第一个预声明 parent 的 weak/ring 两条件，共 8 条件；这是 **exposed regression pilot**，不是独立最终结论。

三主臂：有限结构编辑 + secant 指导补偿；普通 material multi-start；随机 material weak/null 编辑 + 相同 corrector。共同材料域、相同对源和 gain 的知识、相同基础解、same-prior 及相同 forward engine 权限。

- 每条件每臂最多 12 个 full 材料状态，8×3×12=288；
- 共用最多 32 个 base/恢复状态、48 个预注册拒绝点 audit、16 个数值复查；
- **总 cap=384 状态×频率尝试，vector RHS≤9216，GPU=0，CPU wall≤3600 秒**；
- primary K=4；不合法或不足 K 如实报告；
- Jacobian/adjoint/局部更新 setup 全部进入 RHS 与 wall 账本；预算耗尽就停止，不免费构造完整 J。

至少加入两个机制消融：完全相同编辑不做 nonlinear compensation；同一搜索使用 full engine 和 local cached engine。可在固定小列表上执行 engine 消融，避免另开全矩阵实验。相同 engine 也允许给基线使用，防止把通用求解器升级当成生成器优势。

**最短 NO-GO**：四个 parent 均没有新增结构显著、数据相容的世界；或者只是增加微扰图、Selected/Best-MAE 无意义、同预算多起点全面更优。预算内未收敛只构成当前机制 NO-GO，不证明全域无解。

若找到远大于 A30 原局部范围的合法对，完整披露编辑幅度、保护区域、材料支持是否扩大以及原旧规则仍然失败。不能重写 A30 的历史 Gate。

### E3：条件式三维 vector-DDA/VIE 验证

仅在 E2 通过后开启，不默认运行。复用已实现、已核对自项和 polarizability 的 A17/A23 三维 backend。

至少四个非同心 parent：偏心空腔、双孔/分裂、非对称分层、局部高对比损耗块。材料输出必须是三维完整场，不能只有中心球四参数。建议起始至少 216 个三维单元、多个独立偏振照明；实际网格必须由原 backend 的稳定性和离散误差检查决定，而非为了达数字增加病态单元。

三主臂、每 parent 每臂最多 8 个 full state，共 96；最多 32 个 base/独立精度检查：
**pilot cap=128 个状态×频率，vector RHS≤3072，wall≤3600 秒，GPU 默认0。**

primary K=4；目标是证明多照明共同材料补偿机制在非平凡矢量模型中存在并具备预算可行性，四 parent 不足以证明广泛重建优势。至少三个 parent 中出现新的合法结构族，才值得设计确认试验；所有失败保留。

离散一致性不足会造成“伪不可辨”。对预选候选进行更细网格或独立更严容差核验；若误差已接近噪声/成对距离量级，就不能把离散歧义宣称连续物理歧义。

### E4：确认性比较，不自动启动

至少 **24 个、优先32个独立三维 parent**，每 parent 两份噪声，三算法 seeds 作为技术重复；最终样本量依据 pilot 的 parent 方差与最小效应做 power/budget 设计，不由“32”这个数字保证充分性。

比较 SOM/TSOM、multi-start material GN、random-null、简单 sampler、冻结 C1 可兼容基线及主法。C1 无法直接适配高维 target 时标 NOT_APPLICABLE，另行训练成本完全单列，不能把低维 checkpoint 当公平高维基线。

所有方法采用相同 K、相同材料/nuisance/prior support、相同 wall 与 full-physics 预算曲线。第一阶段筛出两三个强对手后再做锁定比较，避免六臂全组合浪费费用。

总预算按

$$

N_{\rm parent}\times N_{\rm noise}\times N_{\rm arms}
\times B_{\rm per\ query}

$$

加固定 setup/audit 明确登记。没有该登记，不进入 E4。主要结论使用 parent 配对区间，noise、seed、source 不视为独立物体。

只有 E4 通过，才开展更大复杂形状、跨分辨率/OOD 和多 query learned proposal 摊销。不会为追求“完整故事”自动开始大型训练。

---

## 13. 反例与失败假设台账

| 待推翻的假设 | 最小反例或诊断 | 决策 |
|---|---|---|
| receiver-dark current 都能材料化 | 多源场为1、2，暗电流增量为1、1；式(2)失败 | 禁止 raw current-null division |
| Jacobian null 意味有限等响应 | $f(a,b)=(a,b^2)$ 原点 | 要求常秩或 full nonlinear correction |
| 一个局部分支足以覆盖全域 | $f(a,b)=(a,b^2)$，目标 $(0,1)$ | 引入 global/edit/population |
| 小参数维数意味着低秩物理更新 | 一个全域 contrast 参数更新全部对角块 | 计入实际 current rank |
| 保持强材料坐标就保持测量 | $f(a,b)=a+b^2$ | 使用强坐标 recoil |
| P seed 为零可以删除全部反馈 | 特定 forcing 为零而传播 Schur 项非零 | 不按 O/P/M 标签删传播 |
| passivity 或小状态 residual 给安全误差界 | 接近奇异/非正规被动离散算子 | 检查 resolvent/inf-sup，不宣称 uniform |
| 几何事件覆盖100%就是材料恢复好 | P2 内部 conditional CDF 仍大；同族有连续变化 | 主指标改完整场 MAE |
| 优化候选 softmax 是 posterior | 两盆地 proposal 到达率不同、体积不同 | candidate-only 或正确校正 |
| diversity 越高 coverage 越好 | 极端稀薄 prior 区域也能远离已有候选 | 可行性+prior合理性+覆盖共同评价 |
| 加第二频率总使每次真值误差下降 | 单次噪声、模型错配、错误共享/独立gain | 区分期望信息与单次估计，不承诺逐例改善 |
| 局部 exact update 没有相关先例 | PNGF、multiport low-rank updates | 把 engine 当既有技术比较 |
| 本次有限搜索没找到另一族证明唯一 | 编辑字典、预算、优化半径不足 | 只报告searched/exhausted/unknown |
| 图像有洞就代表连续介质有洞 | 一个像素阈值翻转、网格伪拓扑 | min feature + persistence + mesh check |

---

## 14. 文献原创性矩阵

检索截止 2026-10-14。该表是证据驱动的相关性审计，不是“未找到即首次”的证明。

| 文献线 | 已有机制 | 与本任务重叠 | 仍需解决 |
|---|---|---|---|
| SOM / TSOM / FFT-TSOM [R1,R2] | 外部电流分解、内部传播约束、粗细子空间 | 基础解 + current ambiguity 已有 | 材料可实现的多分支集合、MAE–覆盖–成本 |
| MR-FFT-TSOM [R3] | 乘性正则与FFT-TSOM | 正则和高效基础成像已有 | 不是多解生成理论 |
| SOM-Net [R4] | unrolled current/material 更新，物理损失 | physics-native current network 已有 | 单个重建网络不自动提供多材料世界 |
| IMSA-SOM [R5] | multiresolution，二维PEC未知形状/数量 | 多尺度表示与搜索已有 | penetrable vector material 与多解区别 |
| e-GLASS [R6] | 从初始解沿 GAN latent 的测量弱方向产生不同解 | 与“基础解→null方向→多图”直接重叠 | 非线性Maxwell、共同材料约束、有限结构编辑 |
| Bhadra等 [R7] | StyleGAN潜空间反复求数据相容的替代图 | empirical multi-solution generation 已有 | 不等于真实posterior或Maxwell材料 |
| Deep Null Space / nonlinear DC [R8,R9] | 线性及非线性数据一致学习、收敛分析 | “非线性null补全”也不能单独称首次 | 多拓扑覆盖与物理成本 |
| Interpolatory / adaptive TR ROM [R10] | output/Jacobian匹配、模型更新、full fallback | receiver-adjoint/刷新框架已有 | 非强制coercive Maxwell 的可信误差与成本 |
| Nonlinear deflation [R11] | 去除已知根寻找其他根 | 多解排斥搜索已有 | noisy feasible tubes、材料流形与完整验证 |
| Bayesian topological priors [R12] | 二维波散射形状后验、MCMC多峰与对象数推断 | 散射多解、拓扑与Bayesian已重叠 | 高维三维矢量材料、多query经济性 |
| NR-SOM+SQP [R13] | 重建 non-radiating subspace 改善成像 | 暗电流物理优化已有 | 未由已审查材料建立多世界覆盖贡献 |
| 2025 conditional null-space score / 2026 DC综述 [R14] | 物理条件、range-null生成与数据一致理论 | 只生成弱部分不是新范式 | Maxwell特定合法性和finite changes |
| ScoreField 2026 preprint [R15] | contrast/current INR + score prior + Lippmann–Schwinger | Maxwell+生成先验+current coupling 已有 | 摘要并未建立本任务的材料多世界集合保证 |
| cloaking / illusion optics [R16] | 用材料结构改变或抵消散射 | “不同物体能有相近/相同外部响应”并非新物理 | 普通受限介质、有限噪声逆问题内的重要解搜索 |
| PNGF / updatable multiport [R17,R18] | 预计算Green、局部全波低秩响应复用 | 式(10)式计算思想已有直接电磁先例 | 同预算可行材料集合与真值覆盖 |
| Multiple Choice Learning [R19] | 多个结构输出、oracle/winner-takes-all损失 | Best-of-K优化和下游择优已有 | 不能由这个损失推导真实后验或Maxwell可行性 |
| Set inversion / bounded-error estimation [R20] | 以有界观测误差定义逆像集合、区间排除 | feasible tube和集合反演不是新概念 | 全波可靠包围算子、高维复杂度及可部署覆盖 |

**最关键的原创性警告：** 不是只需要避开“通用生成模型”；还必须正面比较 e-GLASS、非线性 data-consistent learning、散射多峰形状反演，以及 PNGF 的有限局部全波更新。把它们拼起来可以形成有用算法，但“拼起来”本身不保证强 TAP。

目前未核实到一个同时完成本任务全部条件、且具有相同三维材料覆盖评价的现成系统；这只表示本轮检索没有建立该事实，绝不是不存在或首次的证明。经典2009/2011/2015全文访问仍有缺口，不能宣称已经逐式核对所有版本。详见文末来源访问记录。

---

## 15. Route B：MaterialHypothesisPackage 的科学契约

包内 observation/noise/acquisition/material law/grid/solver tolerance/prior/nuisance 的 identity 必须不可混淆。推荐 manifest + NPZ/Zarr/HDF5 数组及 SHA-256；仅靠文件名不足以共享物理 target。

每个 candidate 必须有完整材料数组、结构描述、全波预测、data fit 与误差类型、补偿参数、合法性、多源一致性、父编辑、cost 和 ROM validity。候选对包含真实 material distance、measurement distance、nuisance 策略、可验证连接状态；graph cluster 只是分组，不自动代表独立 mode。

四种 weight reliability：
- `reference`：与声明reference一致的质量，仍附reference数值限制；
- `corrected`：有合法校正、ESS/误差与漏模式诊断；
- `learned_proposal`：生成分布权重，不冒充posterior；
- `candidate_only`：无概率，权重字段必须为空。

`unexplored_uncertainty` 默认 unknown；candidate-only 不能给“剩余posterior质量0”。若 material support 之外未建模，记录 `outside_support = unquantified`。

**新测量接口需要额外信息。** 当前预测 $F_a(\chi_k)$ 不能代替新 action $a'$ 的预测。包须给出可复现 solver/material identity，以及 action-response cache 或 evaluator 契约；Route B 真正调用新 forward 时另计预算。

当仅有候选、不具有 proper weights，B 可做候选条件的 minimax pair separation、worst-case decision regret 或明确标为 heuristic 的覆盖设计，不能把 uniform candidate weights 的 entropy gain 写成 calibrated Bayesian information gain。已知 pair distance 还应在新增 action 的噪声与仪器模型下重算。

---

## 16. 对十三个核心问题的明确回答

| 问题 | 回答 |
|---|---|
| 1. TSOM暗电流能否成为材料补偿生成器？ | 只能在共同材料可实现性、全部照明和物理约束闭合后；不是直接采样暗电流再相除。 |
| 2. ROM最适合哪里？ | 首先是proposal或受审计screening，其次局部补偿；基础成像需赢同质量成本；full inference最严格，不默认使用。 |
| 3. 局部流形何时存在？ | 可微/常秩、适定forward、合法内点；材料投影还要去掉nuisance/gauge维数。 |
| 4. 怎么找远处分支？ | 显式非局部结构编辑、多chart、population/global moves；局部continuation本身不保证。 |
| 5. 怎么保真且多样？ | 材料合法和full残差作硬约束，再做结构covering；拒绝点随机核验，不能让diversity买掉数据误差。 |
| 6. 能否不学完整posterior？ | 能，先做确定性/随机candidate-only compensation，且本轮应如此。 |
| 7. Flow/MDN用在哪里？ | 低维edit/compensator warm-start和跨query摊销；强数据坐标、物理约束与最终核验留给优化。 |
| 8. 同预算Best-of-K最多改善多少？ | 目前无可信数值预测；式(24)定义信息/支持下的理想集合风险。有限bank最佳值不是连续全局界；实际还受预算限制。 |
| 9. 在线选择能兑现多少oracle？ | 式(21)–(23)给不可辨下界，式(25)把同一集合的Best-of-K与最优在线选择分开；一般不能全部兑现。 |
| 10. 哪些路线不能保证发现全部重要解？ | 所有有限预算局部/有限chart/有限编辑/有限采样方法，以及支持外世界；即使算法正确仍可能漏。 |
| 11. A17/A23有何不可替代增量？ | 没有证据证明“不可替代”；它们提供可复用的反馈解释、双侧模型修正与finite response框架，但标准full solver/Green更新也可实现。 |
| 12. 何时停TSOM/ROM？ | 材料一致性修复太贵、base质量差、刷新频繁、冷启动不省算、覆盖不及material-first时立即删除，不牵连整个多解方向。 |
| 13. 最有希望的论文路线？ | 见下一节；首选受约束有限编辑集合，备选population材料探索，第三是与B结合的可区分性驱动闭环。 |

---

## 17. 论文定位与最后裁决

### 路线 I：Finite-edit compensated material ensembles

核心是：在受限、可实现的三维 Maxwell 材料域内，发现哪些有限结构变化可被真实材料/nuisance 吸收，利用 exact/controlled response 评估，优于等预算多起点的 MAE–覆盖–成本。

最低要求：非平凡条件机制、material realizability、同预算强基线、非同心三维 vector 实验。式(10)已有低秩更新先例，式(12)是标准局部余项论证；两者不能单独承担原创定理。若只有漂亮候选图而没有覆盖成本优势，强 TAP 证据不足。

### 路线 II：Global multiscale material-set exploration

不依赖TSOM、不强制NN，针对单一anchor和固定chart的盲区，研究显式拓扑moves与受约束材料优化如何提高有效解族发现率。必须超过普通population、多起点和预算受限采样，而非换名。

### 路线 III：Compensation-aware imaging and measurement handoff

A负责保留真正难区分的材料世界；B选择能打破补偿的下一次测量。在相同总成像+测量成本下，降低最终 material MAE/误选风险。可以利用式(21)的选择下界和实际新增数据收益形成闭环，但A的候选漏解及未校准权重必须显式传播。

**建议顺序：I 主线，II 强对照/备选，III 在A交接有效后开展。** 不将新网络塞回 A17，不因历史投入强制复用TSOM/ROM，不把极低维P2的成功升级为开放三维成像。

**当前最终状态：GO for bounded falsification; NO CLAIM of a strong-TAP result yet.**
先做 E0–E2；若能越过A30原失败并在真实材料指标上击败简单多起点，再支付三维与摊销训练成本。若不能，就停止这项算法组合，保留物理见证、反例与可复现负结果。

---

## 附录：来源与访问层级

### 项目证据（P）

P1. `migodam/a17-maxwell-certified-exchange`：
`closure_final/START_HERE.md`；
`closure_final/manuscript/PAPER_TAP_FINAL.tex`；
`CORE_THEORY.tex`；
`FINAL_REFERENCES.tex`；
`closure_final/reports/A17_CLOSURE_REPORT_ZH.md`。
https://github.com/migodam/a17-maxwell-certified-exchange
访问：GitHub 原文；核心理论全文，中文报告主要结果段。没有重新运行外部依赖完整回放。

P2. `migodam/a23-fast-material-imaging`：
`A23_EXECUTIVE_DECISION.md`、
`BAYESIAN_FEEDBACK_INVERSION.md`、
`OPM_REDUNDANCY_AND_MINIMALITY.md`、
`PHYSICS_NATIVE_PRIORS.md`、
`THEOREM_AND_COUNTEREXAMPLE_LEDGER.md`。
https://github.com/migodam/a23-fast-material-imaging
访问：上述原文；部分长响应显示有尾部截断，未用不可见段支持结论。

P3. `migodam/a24-source-anchored-om-feedback/README.md`。
https://github.com/migodam/a24-source-anchored-om-feedback
访问：README 的完整Gate、方法、预算和结果摘要；没有宣称完整审查所有实验代码。

P4. `migodam/a30-maxwell-observable-material-manifold/A30_E0_E1_DECISION.md`。
https://github.com/migodam/a30-maxwell-observable-material-manifold/blob/main/A30_E0_E1_DECISION.md
访问：原始主要决策、数字、方法、失败与执行偏差；长文末尾部分未展开。旧具体冻结commit应由执行阶段manifest固定。

P5. Joint P2，branch `joint-p2-20261014`：
https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/joint-p2-20261014/research/joint_p2/RESEARCH_REPORT_ZH.md
访问：报告主要全文分段核对，包括第7节C1数字。没有 materialize 原始response数组或新跑Mie。

P6. C1-P1：
https://github.com/migodam/c1-physics-state-generative-imaging/blob/main/p1/C1_P1_RESEARCH_REPORT_ZH.md
访问：原报告主要正文，方法、参考精度、成绩与失败。P2新结果以P5为准。

### 外部文献（R）

R1. Chen, *Subspace-Based Optimization Method for Solving Inverse-Scattering Problems*, IEEE TGRS 48(1),42–49,2010. DOI:10.1109/TGRS.2009.2025122.
Zhong & Chen, *Twofold subspace-based optimization method for solving inverse scattering problems*, Inverse Problems25,085003,2009. DOI:10.1088/0266-5611/25/8/085003.
https://doi.org/10.1088/0266-5611/25/8/085003
访问：出版记录/摘要、项目对原文的引用和相关SOM-Net表述；经典全文未成功取得。本报告不对未见定理编号、精确截断规则作主张。

R2. Zhong & Chen, *An FFT Twofold Subspace-Based Optimization Method for Solving Electromagnetic Inverse Scattering Problems*, IEEE TAP59(3),914–927,2011.
https://doi.org/10.1109/TAP.2010.2103027
访问：出版信息和摘要层；不是完整原算法复现。

R3. Xu et al., *Multiplicative-Regularized FFT Twofold Subspace-Based Optimization Method for Inverse Scattering Problems*, IEEE TGRS53(2),841–850,2015.
https://doi.org/10.1109/TGRS.2014.2329032
访问：出版记录/元数据。MR缩写按原题名。

R4. Liu et al., *SOM-Net: Unrolling the Subspace-Based Optimization for Solving Inverse Scattering Problems*,2022 preprint.
https://arxiv.org/abs/2209.03567
访问：作者HTML正文的模型、算法和loss部分。全文访问后续重试有缓存失败，不将其训练约束称为full-Maxwell certificate。

R5. Ye et al., *Multiresolution Subspace-Based Optimization Method for the Retrieval of 2-D Perfect Electric Conductors*,2022 preprint / IEEE TMTT publication.
https://arxiv.org/abs/2208.10944
https://doi.org/10.1109/TMTT.2022.3220252
访问：作者摘要/出版信息；避免将preprint年份与期刊卷年混用。

R6. Montanaro, Valsesia, Magli, *Exploring the Solution Space of Linear Inverse Problems with GAN Latent Geometry*,ICIP2022.
https://arxiv.org/abs/2207.00460
https://doi.org/10.1109/ICIP46576.2022.9897211
访问：5页作者PDF解析正文，重点算法与范围。PDF截图服务报错，未独立视觉核验图表；本文不依赖其图像测量或速度数字。

R7. Bhadra, Villa, Anastasio, *Mining the manifolds of deep generative models for multiple data-consistent solutions of ill-posed tomographic imaging problems*,2022.
https://arxiv.org/abs/2202.05311
访问：作者摘要；明确empirical sampling，不声称posterior correctness。

R8. Schwab, Antholzer, Haltmeier, *Deep Null Space Learning for Inverse Problems: Convergence Analysis and Rates*,Inverse Problems35,025008,2019.
https://arxiv.org/abs/1806.06137
访问：作者HTML正文的框架和收敛设定。

R9. Boink, Haltmeier, Holman, Schwab, *Data-consistent neural networks for solving nonlinear inverse problems*,Inverse Problems and Imaging17(1),203–229,2023；2022early access.
https://doi.org/10.3934/ipi.2022037
https://arxiv.org/abs/2003.11253
访问：出版社摘要、方法说明与figure/table captions；没有宣称其适用条件已在本Maxwell问题验证。

R10. de Sturler et al., *Nonlinear Parametric Inversion Using Interpolatory Model Reduction*,SIAM SISC37,B495–B517,2015.
https://arxiv.org/abs/1311.0922
https://doi.org/10.1137/130946320
Kartmann et al., *Adaptive Reduced Basis Trust Region Methods for Parameter Identification Problems*,2024.
https://arxiv.org/abs/2309.07627
https://doi.org/10.1007/s44207-024-00002-z
访问：前者作者摘要，后者HTML正文的算法与适用条件；非Maxwell coercivity结果不可直接移植。

R11. Farrell, Birkisson, Funke, *Deflation Techniques for Finding Distinct Solutions of Nonlinear Partial Differential Equations*,SIAM SISC37,A2026–A2045,2015.
https://arxiv.org/abs/1410.5620
https://doi.org/10.1137/140984798
访问：作者全文的deflation定义、continuation比较与适用限制。

R12. Carpio, Iakunin, Stadler, *Bayesian approach to inverse scattering with topological priors*,Inverse Problems36,105001,2020.
https://arxiv.org/abs/2003.09318
https://doi.org/10.1088/1361-6420/abaa30
访问：30页PDF解析正文，重点MCMC多峰与二维Helmholtz范围；截图失败，不将图像细节当独立证据。

R13. Siampour & Nezhad, *Revealing the Invisible: Imaging Through Non-Radiating Subspace*,JOPR1(4),159–169,2024.
https://doi.org/10.47852/bonviewJOPR42022785
访问：出版社原始摘要及作者机构记录，非全文结论审计。

R14. Zhang et al., *Score-Based Generative Model With Conditional Null-Space Learning for Limited-Angle Tomographic Reconstruction in Medical Imaging*,IEEE TCI,2025.
https://doi.org/10.1109/TCI.2025.3562059
Haltmeier & Hwang, *Data-Consistent Learning of Inverse Problems*,2026preprint.
https://arxiv.org/abs/2601.12831
访问：前者学术索引摘要，后者作者HTML正文。二者不是本任务完成证据。

R15. Guo, Gao, Sun, *ScoreField: Neural Inverse Scattering with Score-Based Generative Priors*,2026-08-03preprint.
https://arxiv.org/abs/2608.02937
访问：作者摘要；PDF未成功取得。仅用作需比较的相近方向，不据摘要宣布其全部实验性质或未做的工作。

R16. Alù & Engheta, *Achieving transparency with plasmonic and metamaterial coatings*,Phys.Rev.E72,016623,2005；另有2006erratum.
https://doi.org/10.1103/PhysRevE.72.016623
Lai et al., *Illusion Optics: The Optical Transformation of an Object into Another Object*,PRL102,253902,2009.
https://doi.org/10.1103/PhysRevLett.102.253902
访问：APS原始摘要与记录。材料实现域不同，不把其理想装置当普通介质逆问题的可行候选。

R17. Sun et al., *Near real-time full-wave inverse design of electromagnetic devices*,Nature Communications17,2372,2026-02-14.
https://doi.org/10.1038/s41467-026-69477-y
访问：出版社HTML正文的current equivalence、precomputation与low-rank update推导；本报告不复用其特定speedup数字作为本算法成绩。

R18. Prod'homme et al., *Updatable Closed-Form Evaluation of Arbitrarily Complex Multiport Network Connections*,Advanced Electronic Materials,2026-05-13.
https://doi.org/10.1002/aelm.202500565
访问：出版社检索摘要；多端口结构与本材料成像不同。



**R19 — Multiple Choice Learning / oracle-loss 先例**  
Guzmán-Rivera, Batra, Kohli, “Multiple Choice Learning: Learning to Produce Multiple Structured Outputs,” NeurIPS 2012.  
https://proceedings.neurips.cc/paper/2012/hash/cfbce4c1d7c425baf21d6b6f2babe6be-Abstract.html  
Lee et al., “Stochastic Multiple Choice Learning for Training Diverse Deep Ensembles,” NeurIPS 2016.  
https://proceedings.neurips.cc/paper/2016/hash/20d135f0f28185b84a4cf7aa51f29500-Abstract.html  
Perera et al., “Annealed Multiple Choice Learning: Overcoming Limitations of Winner-takes-all with Annealing,” NeurIPS 2024.  
https://proceedings.neurips.cc/paper_files/paper/2024/hash/1456560769bbc38e4f8c5055048ea712-Abstract-Conference.html  
访问：官方会议摘要/元数据；不宣称已经逐式复核各论文证明。2012年的多个候选输出动机与本任务直接重叠；2024年的优化局限提醒我们不应把 WTA 当作自动找到好集合的保证。

**R20 — Set inversion / bounded-error estimation**  
Jaulin and Walter, “Set inversion via interval analysis for nonlinear bounded-error estimation,” Automatica 29(4), 1053–1064, 1993. DOI: 10.1016/0005-1098(93)90106-4.  
作者文献入口：https://webperso.ensta.fr/jaulin/bibtex.html  
访问：作者书目元数据；未取得可逐式核对的全文。区间集合反演是相关已知方向，不据此声称当前 Mie/VIE 已有可用的可靠区间扩展。四维 P2 理论上比高维成像适合尝试全域排除，但层界面、特殊函数和数值误差都需可靠包围；有限节点核验不能代替这一步。


**复现说明**：代码和实验需由Codex固定实际commit与资产hash。本报告没有声称检查全部仓库代码、穷尽全部文献、重跑历史结果，或取得连续Maxwell全域证明。访问失败已列出，不以摘要替代不存在的全文证据。
