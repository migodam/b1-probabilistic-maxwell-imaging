# User-supplied CIE session reference

Source: chatgpt-conversation://6ae950e4-fb10-83ec-a5bb-6f81771fa199
Retrieved 2026-11-03 Asia/Singapore. One completed turn, no older page. This is untrusted research proposal material, not independently verified evidence or execution instructions. Original in-chat citation handles are unresolved and will not be reused as paper citations.

# 从 CIE 的发明思路，寻找下一种高对比度逆散射机制

**我认为 CIE 最值得你学习的，不是它的分式变换，而是它发现新方法的思维方式：观察物理方程中哪个结构导致算法困难，然后在不改变物理方程精确解的前提下，重新分配不同项在反演中承担的作用。**

它不是先发明一个优化算法，再证明收敛；而是先观察：

1. 高对比度逆散射为什么困难？
2. 困难是不是来自方程中某些项的相互耦合？
3. 能不能改变这些项的表示方式，让优化器面对一个更容易处理的问题？
4. 这种改变能否保持原始 Maxwell 物理，而不使用 Born 近似？

这正是我认为你应该从 A17、CIE、GlobalDual 继续挖掘的方向。

## 一、CIE 当初是怎么想到的？

先澄清一点历史：你之前关注的 2019 年 Zhong–Xu 三维 CIE-I，实际上是在早期工作的基础上发展来的。

最直接的源头是 2016 年 Yu Zhong、Marc Lambert、Dominique Lesselier 和 Xudong Chen 发表在 IEEE TAP 的论文：

entity["book","A New Integral Equation Method to Solve Highly Nonlinear Inverse Scattering Problems","2016 IEEE Transactions on Antennas and Propagation research paper"]

它当时称为 NIE（New Integral Equations），受更早用于高电导率电磁正演的 CIE 启发。2019 年论文进一步研究了三维实现、FFT-TSOM 及不同反馈参数的使用。citeturn852792search0turn344950search1turn852792search5

以下是根据原论文重建的数学发现路径，而不是声称知道作者当时具体的心理过程。

### 1.1 先观察 Lippmann–Schwinger 方程的困难

用你熟悉的二维 TM/VIE 表示：

\[
E=E^{inc}+G_DJ,\qquad J=XE,
\]

其中 \(X=\operatorname{diag}(\chi)\)。

消去总场：

\[
\boxed{J=X(E^{inc}+G_DJ)}
\]

这是标准 LS 方程。

问题出在这个结构：

\[
\underbrace{J}_{\text{未知电流}}
=
\underbrace{XE^{inc}}_{\text{一次散射}}
+
\underbrace{XG_DJ}_{\text{多次散射反馈}}.
\]

在弱散射情形，多次散射反馈比较容易处理。

但在高对比度情形，电流 \(J\) 不仅依赖本单元材料，还依赖整个区域的材料分布。每个单元都会通过 Green 算子影响其他单元，并产生多轮相互散射。

传统 CSI/SOM 要同时估计材料和电流，因而面对一个强耦合的双线性约束。

**这里的关键观察：困难不只是 \(\chi\) 很大，而是一个较大的材料系数直接乘在全局耦合项 \(G_DJ\) 上。**

### 1.2 CIE 的关键灵感：人为加入一个局部项

观察下面的恒等式：

\[
J=X(E^{inc}+G_DJ).
\]

在等式两边乘以参数 \(\beta\)，再通过增加一个局部电流项进行整理：

\[
\beta J
=
\frac{\beta X}{I+\beta X}
\left(\beta J+E^{inc}+G_DJ\right).
\]

这里分式表示对角矩阵运算，定义

\[
\boxed{R_\beta=\beta X(I+\beta X)^{-1}}.
\]

得到

\[
\boxed{
\beta J
=
R_\beta\left(\beta J+E^{inc}+G_DJ\right)
}
\tag{CIE}
\]

这就是 CIE-I 的核心形式。这里先取非零实数常量 \(\beta\)，空间变化的 \(\beta(r)\) 则有相应推广。

注意它做了两个动作：

- 将原来的材料参数 \(X\) 改写成有界化倾向的参数 \(R_\beta\)；
- 在原来的全局散射场 \(G_DJ\) 之外，加入局部反馈 \(\beta J\)。

而且这是**严格等价变换，没有丢掉高阶散射**。

对于标量材料单元：

\[
R_\beta(\chi)=\frac{\beta\chi}{1+\beta\chi}.
\]

当 \(\chi\geq0,\ \beta>0\) 时：

\[
0\leq R_\beta<1.
\]

即使材料对比度非常高，变换后的系数也不会随 \(\chi\) 无限制增大。

但我要强调，\(R_\beta<1\) 本身不代表整个迭代算子一定收缩；这一点后面会直接影响我们能否设计真正的新算法。

### 1.3 最值得观察的是这个结构变化

原始 LS

\[
J=X(E^{inc}+G_DJ)
\]

材料系数直接耦合

```svg
<svg viewBox="0 0 220 110" width="100%" xmlns="http://www.w3.org/2000/svg"><rect x="8" y="35" width="78" height="40" rx="7" fill="#DDA457" fill-opacity=".24" stroke="#C18B3D" /><text x="47" y="58" text-anchor="middle" font-size="13" fill="currentColor">X</text><path d="M86 55 H124" stroke="currentColor" stroke-width="1.6" /><path d="M117 50 L125 55 L117 60" fill="none" stroke="currentColor" stroke-width="1.6" /><rect x="126" y="35" width="85" height="40" rx="7" fill="#CD6969" fill-opacity=".20" stroke="#C16868" /><text x="168" y="58" text-anchor="middle" font-size="13" fill="currentColor">Gᴅ J</text><path d="M168 75 V96 H47 V75" fill="none" stroke="#C16868" stroke-width="1.5" stroke-dasharray="5 3" /><path d="M42 83 L47 75 L52 83" stroke="#C16868" stroke-width="1.5" fill="none" /></svg>
```

全局多次散射直接反馈材料更新

变换后的 CIE

\[
\beta J=R_\beta(\beta J+E)
\]

引入局部反馈分支

```svg
<svg viewBox="0 0 220 110" width="100%" xmlns="http://www.w3.org/2000/svg"><rect x="7" y="35" width="70" height="40" rx="7" fill="#36968A" fill-opacity=".20" stroke="#36968A" /><text x="42" y="58" text-anchor="middle" font-size="12" fill="currentColor">Rβ</text><path d="M77 55 H111" stroke="currentColor" stroke-width="1.6" /><path d="M105 50 L113 55 L105 60" fill="none" stroke="currentColor" stroke-width="1.6" /><rect x="115" y="35" width="98" height="40" rx="7" fill="#DDA457" fill-opacity=".20" stroke="#C18B3D" /><text x="164" y="52" text-anchor="middle" font-size="12" fill="currentColor">β J + E</text><text x="164" y="65" text-anchor="middle" font-size="9" fill="currentColor">local + global</text><path d="M164 75 V95 H42 V75" fill="none" stroke="#36968A" stroke-width="1.5" stroke-dasharray="5 3" /><path d="M37 83 L42 75 L47 83" stroke="#36968A" stroke-width="1.5" fill="none" /></svg>
```

局部项与全局项共同参与更新

原来只有 \(G_DJ\) 体现未知电流的反馈；现在又有了局部、无须跨空间传播的 \(\beta J\)。

当 \(\beta J\) 在相应的材料更新关系中占主导时，全局未知散射场的相对影响便可能降低。

**CIE 的真正思维突破，是通过一个局部代数项重新分配逆问题中的非线性耦合，而非简单让 Maxwell 散射变弱。**

这为你寻找新机制提供了一个很好的起点。

## 二、把 CIE 再拆一层：一个更重要的恒等式

我觉得这里有一个比 \(R_\beta\) 本身还值得观察的公式。

定义原始 LS 的物理残差：

\[
r_{\mathrm{LS}}(X,J)
=
J-X(E^{inc}+G_DJ).
\]

定义 CIE 残差：

\[
r_{\mathrm{CIE}}
=
\beta J-R_\beta(\beta J+E^{inc}+G_DJ).
\]

将 \(R_\beta=\beta X(I+\beta X)^{-1}\) 代入，整理得到：

\[
\boxed{
r_{\mathrm{CIE}}
=
\underbrace{\beta(I+\beta X)^{-1}}_{P_\beta(X)}
r_{\mathrm{LS}}
}
\tag{1}
\]

**也就是说，CIE 还可以被理解为：对原始物理方程施加一个随材料变化的左预条件器。**

这个观察非常重要。

### 2.1 它改变了什么？

原始 CSI 的状态方程损失可以写成：

\[
\mathcal L_{\mathrm{LS}}
=
\|r_{\mathrm{LS}}\|_2^2.
\]

变换后：

\[
\boxed{
\mathcal L_{\mathrm{CIE}}
=
r_{\mathrm{LS}}^*
P_\beta^*P_\beta
r_{\mathrm{LS}}
}
\tag{2}
\]

原来的物理约束没有改变：

\[
r_{\mathrm{LS}}=0
\iff
r_{\mathrm{CIE}}=0,
\]

前提是变换非奇异。

但是当你尚未满足 Maxwell 方程、还在联合优化 \(X,J\) 时，两个损失的梯度和曲率可以非常不同。

**于是，虽然两者物理解相同，优化器走过的中间状态可以完全不同。**

这是 CIE 能够改善某些逆散射优化过程的一个更精确的数学解释。

### 2.2 高对比度区域究竟发生了什么？

对于实数材料：

\[
P_\beta(\chi)
=
\frac{\beta}{1+\beta\chi}.
\]

例如取 \(\beta=2\)：

| 物理对比度 \(\chi\) | 修改后 \(R_\beta\) | 残差权重 \(P_\beta\) |
|--------|--------|--------|
| 0 | 0 | 2.000 |
| 1 | 0.667 | 0.667 |
| 5 | 0.909 | 0.182 |
| 10 | 0.952 | 0.095 |
| 20 | 0.976 | 0.049 |

当 \(\chi\) 增大时，两个现象同时出现：

- 修改后的材料系数趋向于 1；
- 原始物理残差在 CIE 损失中的权重下降。

这也解释了为什么 CIE 可能改变高对比度区域的优化行为。

但它存在代价：

\[
\boxed{
\chi=\frac{R_\beta}{\beta(1-R_\beta)}
}
\]

\[
\frac{d\chi}{dR_\beta}
=
\frac{1}{\beta(1-R_\beta)^2}.
\]

当 \(R_\beta\to1\) 时，反变换越来越敏感。

因此，不应该简单得出“\(\beta\) 越大越好”的结论。原始 2016 年论文事实上采用了较大 \(\beta\) 到较小 \(\beta\)、同时逐步增加电流子空间维度的策略，并明确指出强正则化和复杂电流表示能力之间的冲突。citeturn303007view1

### 2.3 这还揭示了一个科研上非常重要的限制

假设你直接用真实 Maxwell 正演消去电流：

\[
F(\chi)=G_S(I-XG_D)^{-1}XE^{inc},
\]

然后最小化：

\[
f(\chi)=\frac12\|F(\chi)-y\|^2.
\]

如果仅仅进行一个光滑、可逆的材料重参数化：

\[
\chi=T(\theta),
\]

那么：

\[
\nabla_\theta f
=
DT(\theta)^T\nabla_\chi f.
\]

在驻点处：

\[
\boxed{
H_\theta=DT^T H_\chi DT
}
\tag{3}
\]

这里将复数观测视为实变量问题，且 \(DT\) 可逆。

由 Sylvester 惯性定律，正、负、零曲率方向的数目不会因为这样的坐标变化而改变。

所以：

**单纯将 \(\chi\) 换成 \(R_\beta\)，并不能消灭原始物理数据目标中的错误局部极小值。**

它能改变优化步的几何和有限预算下的吸引域，但不能仅凭这个代数变换就证明目标已经被凸化。

这也解释了你 FBG-R2 中一个很有价值的现象：

你已有的 FBG-R2 新鲜高对比度实验

Native TRF

# 0/12

GlobalDual

# 3/12

固定 β=.2 的反馈 Gram

# 8/12

这证明了你测试的辅助目标几何具有实际效果，而不是证明 CIE 分式坐标本身能够消灭局部极小值。

urlFBG-R2 冻结实验报告https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/research_fbg_r2_adaptive_feedback/RESEARCH_REPORT_ZH.md

你的项目也已经测试过 \(B=\beta I\)、对角/块反馈和切向能量回归。因此，下一步不能再简单宣称“把 CIE 的 \(\beta\) 改成矩阵 \(B\)”就是新方法。

必须找出一个**此前尚未充分利用的物理结构**。

## 三、再换一个角度：CIE 实际上改善了什么除法？

观察最普通的材料恢复关系：

\[
J_i=\chi_i E_i.
\]

如果已经有电流 \(J_i\) 和总场 \(E_i\)，那么：

\[
\chi_i=\frac{J_i}{E_i}.
\]

对于强散射体，总场可能因为干涉而变得很小，导致材料估计极其敏感。

CIE 则使用：

\[
\boxed{
R_{\beta,i}
=
\frac{\beta J_i}{\beta J_i+E_i}
}
\tag{4}
\]

它在分母中人为引入 \(\beta J_i\)。

于是，**当当前估计的总场 \(E_i\) 很小时，CIE 不必立即通过一个很小的 \(E_i\) 计算很大的材料参数。**

不过这只是改善中间的参数估计：如果物理上真的不存在足够信息，CIE 并不会凭空恢复可辨识性。

这启发了一个特别值得研究的问题：

> 我们能否不再只修改“材料除以场”的关系，而是修改整个系统里最容易引发错误搜索的物理耦合结构？

下面是我认为最值得你盯着看的几组公式。

---

## 四、候选机制 A：不要压低对比度，而是吸收真正困难的散射循环

(我最推荐先研究的结构之一)

### 4.1 两个散射单元就能暴露问题

设一个非常简单的双单元模型：

\[
X=
\begin{bmatrix}
\chi_1&0\\
0&\chi_2
\end{bmatrix},
\qquad
G_D=
\begin{bmatrix}
0&g_{12}\\
g_{21}&0
\end{bmatrix}.
\]

这里为突出单元间耦合，暂时忽略自作用项。

那么：

\[
(I-XG_D)^{-1}
=
\frac{1}{1-\chi_1\chi_2g_{12}g_{21}}
\begin{bmatrix}
1&\chi_1g_{12}\\
\chi_2g_{21}&1
\end{bmatrix}.
\]

最应该观察的是分母：

\[
\boxed{
1-\underbrace{\chi_1\chi_2g_{12}g_{21}}_{\text{两单元散射回路增益}}
}
\tag{5}
\]

这说明一个很重要的事实：

**高对比度问题的困难未必主要来自每个 \(\chi_i\) 的数值，而可能来自不同散射单元之间形成的强反馈回路。**

例如即使 \(\chi_1,\chi_2\) 都不是特别大，只要这个回路增益接近 1，就会出现强烈的响应和敏感性。

CIE 的 \(\beta I\) 并没有直接针对 \(g_{12}g_{21}\) 这个结构。

那能不能针对它？

### 4.2 推广 CIE：吸收一部分真正的 Green 算子

回到完整方程：

\[
J=X(E^{inc}+G_DJ).
\]

对于任意适当的线性算子 \(B\)，有：

\[
(I+XB)J=X[E^{inc}+(G_D+B)J].
\]

若 \(I+XB\) 可逆：

\[
\boxed{
J=(I+XB)^{-1}X
[E^{inc}+(G_D+B)J]
}
\tag{6}
\]

标准 CIE 就是 \(B=\beta I\) 的特例。

但如果我们把 Green 算子分解为：

\[
G_D=G_{\mathrm{in}}+G_{\mathrm{cross}},
\]

其中：

- \(G_{\mathrm{in}}\) 表示选定单元组内部的散射；
- \(G_{\mathrm{cross}}\) 表示不同单元组之间的传播；

然后令

\[
B=-G_{\mathrm{in}},
\]

就得到：

\[
\boxed{
J=
\underbrace{(I-XG_{\mathrm{in}})^{-1}X}_{T_{\mathrm{in}}(X)}
(E^{inc}+G_{\mathrm{cross}}J)
}
\tag{7}
\]

这很有意思。

原来需要反复处理的组内多次散射，现在被完整吸收到一个有效算子 \(T_{\mathrm{in}}\) 里面。

而优化器面对的是组间相互作用。

从局部反馈走向物理散射吸收

这与普通对角 \(\beta\) 的区别在于：

- \(\beta\) 是人为添加的局部代数反馈；
- \(G_{\mathrm{in}}\) 来自真实 Maxwell 传播结构；
- \(T_{\mathrm{in}}\) 可以精确保留组内任意阶多次散射。

如果 \(G_{\mathrm{in}}\) 是块对角的，那么各个 \(T_{\mathrm{in}}\) 可以由小规模局部求解得到。

如果它是低秩的，还可以使用 Woodbury 恒等式：

\[
(I-XUKV^*)^{-1}
=
I+XU(I-KV^*XU)^{-1}KV^*.
\]

这样可能把昂贵的全局运算转移到小矩阵中。

**但这里真正需要的新东西，是如何依据可观测散射回路构造这些块，而不是普通的空间分块。** 块消元、T-matrix 和 Woodbury 本身都不是新数学。

还要特别小心：如果直接取 \(B=-G_D\)，尽管右侧不再有显式多次散射，求解 \((I-XG_D)^{-1}\) 仍然是原始困难问题。这只是转移了计算量。

所以你应该研究的是：

\[
\text{哪些散射回路最难，但又最便宜地被精确吸收？}
\]

---

## 五、候选机制 B：从“材料反馈”转向“观测反馈”

这个方向与 GlobalDual、P3L 和你之前的 SOM/TSOM relevance mismatch 最接近。

我认为它有可能比不断优化 \(\beta\) 更有意义。

### 5.1 观察一个精确的误差传播恒等式

定义：

\[
A(X)=I-XG_D.
\]

给定材料 \(X\)，真实物理电流满足：

\[
A(X)J_*=XE^{inc}.
\]

但假设我们在联合反演中有一个尚未满足 Maxwell 方程的电流候选 \(J\)，它的残差是：

\[
r=A(X)J-XE^{inc}.
\]

则必然有：

\[
\boxed{
J-J_*=A(X)^{-1}r
}
\tag{8}
\]

进一步投影到接收端：

\[
\boxed{
G_S(J-J_*)=G_SA(X)^{-1}r
}
\tag{9}
\]

这表示：

**一个物理残差到底有多严重，不只是看它的欧氏范数，而是看它经过真实 Maxwell 系统传播以后，在接收端产生多大的影响。**

CIE 的残差加权主要取决于局部材料 \(X\)。

但式 (9) 告诉我们，也可以从整个测量系统出发设计权重。

定义：

\[
T_{\mathrm{port}}(X)=G_SA(X)^{-1}.
\]

便得到一个候选目标：

\[
\boxed{
\mathcal L_{\mathrm{port}}
=
\|T_{\mathrm{port}}(X)r\|^2
+\eta\|r\|^2
}
\tag{10}
\]

第二项防止算法完全忽略当前观测不到的物理残差。

### 5.2 这个结构可以产生什么机制？

假设有两种电流残差：

\[
r_1,\quad r_2,
\qquad
\|r_1\|=\|r_2\|.
\]

但是：

\[
\|T_{\mathrm{port}}r_1\|
\gg
\|T_{\mathrm{port}}r_2\|.
\]

从普通物理残差看，它们一样严重；从真实测量看，第一种误差更值得优先修正。

这意味着可以构造一种：

**Measurement-aware Physics Residual Geometry**

它不是单纯减小总体物理残差，而是优先纠正会明显破坏数据一致性的内部物理偏差。

不过 \(T_{\mathrm{port}}\) 本身涉及全波逆算子，不能每一步都昂贵地重新构造。可以研究少量接收端口方向、低秩近似或者复用 GlobalDual 中已经计算的量。

它最有价值的延伸是：有些内部电流方向几乎不可见，但会通过后续多次散射影响可见方向，所以不能直接把它们设为零。

例如把方程分成可观测电流 \(J_o\) 和其余电流 \(J_n\)：

\[
\begin{bmatrix}
A_{oo}&A_{on}\\
A_{no}&A_{nn}
\end{bmatrix}
\begin{bmatrix}
J_o\\J_n
\end{bmatrix}
=
\begin{bmatrix}
b_o\\b_n
\end{bmatrix}.
\]

消元得到：

\[
\boxed{
\left(A_{oo}-A_{on}A_{nn}^{-1}A_{no}\right)J_o
=
b_o-A_{on}A_{nn}^{-1}b_n
}
\tag{11}
\]

注意：

\[
\boxed{
\Sigma_{\mathrm{hidden}}
=
A_{on}A_{nn}^{-1}A_{no}
}
\]

这一项代表隐藏电流对子空间内可观测电流的反馈。

这是我觉得很适合你继续想的地方：

> TSOM 为什么要限制电流子空间？如果不可观测电流通过 Schur 补仍然能影响观测，那么我们是不是可以不恢复这些电流，却精确保留它们对观测的反馈？

这和你已经研究的 A17 retained-core Schur、GlobalDual 以及未知墙体的传播等价有自然联系。

但与此前的 Schur-Galerkin 探针相比，新实验必须证明一个明确的新收益，例如保留不可见电流的反馈后，产生了原来拿不到的材料解，或者以更少的全波 RHS 达到相同精度。

---

## 六、候选机制 C：利用 Maxwell 正演的有理函数结构，而不只依赖局部 Jacobian

(我认为最适合进一步做原创探索)

这一方向与 CIE 的分式思路有关，但不再停留于修改材料坐标。

### 6.1 从一个材料更新方向开始

设当前材料为 \(X_0\)，沿一个方向 \(D=\operatorname{diag}(d)\) 更新：

\[
X(t)=X_0+tD.
\]

定义：

\[
A_0=I-X_0G_D,
\qquad
J_0=A_0^{-1}X_0E^{inc},
\]

\[
E_0=E^{inc}+G_DJ_0.
\]

那么更新后的电流：

\[
J(t)=(A_0-tDG_D)^{-1}(X_0+tD)E^{inc}.
\]

经过严格整理，可以得到：

\[
\boxed{
J(t)=J_0+
t(A_0-tDG_D)^{-1}DE_0
}
\tag{12}
\]

继续定义：

\[
M=A_0^{-1}DG_D,
\qquad
u=A_0^{-1}DE_0.
\]

则：

\[
\boxed{
J(t)=J_0+t(I-tM)^{-1}u
}
\tag{13}
\]

于是接收端响应沿这个材料方向的精确变化为：

\[
\boxed{
F(X_0+tD)-F(X_0)
=
tG_S(I-tM)^{-1}u
}
\tag{14}
\]

**这是一个非常值得研究的公式。**

在有限维离散 Maxwell 系统中，沿材料直线的正演响应不是任意非线性函数，而是一个具有明确 resolvent 结构的有理函数。

### 6.2 GN 实际上只看到了第一项

当 Neumann 展开收敛时：

\[
(I-tM)^{-1}
=
I+tM+t^2M^2+\cdots.
\]

所以：

\[
\Delta F(t)
=
tG_Su+
t^2G_SMu+
t^3G_SM^2u+\cdots
\]

第一阶：

\[
DF[d]=G_Su.
\]

第二阶：

\[
D^2F[d,d]=2G_SMu.
\]

更高阶则不断出现：

\[
G_SM^ku.
\]

这些高阶项实际上编码了材料更新方向引起的反复散射反馈。

因此：

**GN 只利用了材料响应有理函数在 \(t=0\) 的切线，而在高对比度下，真正决定能否跨过错误吸引区域的，可能正是后面的分母结构。**

### 6.3 新机制：寻找对接收端真正有影响的极点

式 (14) 的极点满足：

\[
\det(I-tM)=0.
\]

如果 \(M\) 可以对角化，响应可以形式化写成：

\[
\Delta F(t)
=
t\sum_k
\frac{c_k}{1-t\lambda_k},
\]

其中 \(c_k\) 包含源端激励和接收端投影。

非常有意思的是：

某个 \(\lambda_k\) 对应的极点即使很接近真实材料搜索路径，如果它的 \(c_k\) 很小，测量端也可能几乎看不到它。

也就是说，我们真正关心的不是全部内部共振，而是：

\[
\boxed{
\text{Port-visible nonlinear poles}
}
\]

这就把你之前研究的三个东西联系起来了：

- A17-X 的有理材料路径；
- P3L 的高阶可观测性；
- Z-ALG 的弱方向、多解与跨盆地搜索。

### 6.4 怎样把它变成算法？

我设想的是一种 *Pole-aware Rational Search*：

1. 从当前材料产生若干方向 \(d\)，包括普通下降方向和弱可观测方向。
2. 不立刻进行 GN 更新，而是构造 \(M=A_0^{-1}DG_D\) 的矩阵向量作用。
3. 用小规模 Krylov/Arnoldi 子空间建立 \(G_S(I-tM)^{-1}u\) 的有理近似。
4. 在合法材料区间内寻找值得试探的有限步长或非局部候选，而不只在 \(t=0\) 附近线性搜索。
5. 对候选执行真实 Maxwell 正演验收，失败则分析是有理模型截断、近极点不稳定，还是方向本身不合适。

主要计算成本可以是少量复 RHS 求解，而不是每次构造完整材料 Jacobian。

当然，若需要大量方向或模型接近共振，成本也可能迅速增加；非正规算子的伪谱效应也不能只用特征值描述。

这里真正值得研究的新机制不是“使用 Padé”或“使用 Krylov”本身——这些都是成熟技术——而是：

**利用物理响应中接收端可见的极点和散射回路，改变材料更新的有限步探索策略。**

它比再调一个 \(\beta\) 更有可能获得本质不同的搜索行为。

---

## 七、把这些公式与你已有的结果放在一起

你 FBG-R2 有一个很值得反过来研究的发现：

固定标量 \(\beta=0.2\) 表现很好，但通过最小化电流切向能量获得的复杂反馈矩阵并没有更好，甚至表现更差。

这表明“更弱的反馈”不是正确的充分解释。

你已有的条件梯度分解：

\[
\nabla\widehat f_\beta(b\mathbf1)
=
\nabla f_0(b\mathbf1)
-
c\frac{\beta}{1+\beta b}
\operatorname{diag}(QW^*WQ^*)
\]

也给出了一个更深的提示。

右侧修正是一个由观测兼容残差产生、具有特定符号的材料更新驱动项。

所以新的反馈设计不能只问：

\[
\min_B\|(G_D+B)T\|^2.
\]

还应当问：

\[
\text{改变后的更新方向是否能够脱离错误吸引区域？}
\]

而且这种脱离是否能在真实 Maxwell 数据损失和固定计算预算下得到验证。

我会特别注意以下区别：

| 研究思路 | 它改变什么 | 主要风险 |
|--------|--------|--------|
| CIE 分式参数 | 材料表示、联合残差几何 | 反变换敏感，依赖正则 |
| FBG 反馈 Gram | 材料候选的辅助目标几何 | 有效反馈难以预测 |
| 散射回路吸收 | 内部多次散射的计算表示 | 消元可能更贵 |
| Port-aware Schur | 隐藏电流对观测的作用 | 子空间截断失真 |
| 有理极点搜索 | 有限步材料搜索路径 | 共振敏感、近似成本 |

## 八、我建议你真正用来思考新理论的五个问题

**1\. 我能否找到比 \(\chi\) 更本质的“非线性强度”？**

比如是 \(XG_D\) 的谱、特定散射回路增益，还是接收端可见的有理极点？高对比度本身未必是足够好的难度指标。

---

**2\. 我能否提前吸收一部分多次散射？**

不是 Born 近似丢弃它，而是把复杂的高阶反馈精确装进更小的 T-matrix、Schur 补或 resolvent 中。

---

**3\. 哪些内部物理变化是真正值得接收端关心的？**

不是所有大电流、强散射或小奇异值都代表对材料反演有用的信息。应观察端口可见的影响及其高阶变化。

---

**4\. 我能否从局部线性更新转向有限步的物理结构搜索？**

尤其沿 GN 认为很弱、但二阶和高阶响应很强的方向，使用真正的全波有理结构探索新的材料区域。

---

**5\. 这个机制到底改变了物理求解、辅助目标，还是仅仅换了坐标？**

三者需要分别证明与测试。只有这样才能避免把一个等价重参数化误认为新的全局收敛机制。

## 九、如果是我来安排你的下一轮研究

我不会继续让 Codex 大规模搜索 \(\beta\)、拟合 \(B\) 或盲目做更多参数消融。你已经从 FBG-R2 得到足够多证据，说明单纯的切向能量衰减不是好的设计准则。

我会先做两个相对独立、可以互相竞争的理论原型。

**第一条：基于散射回路的局部 T-matrix 吸收。** 先在两单元、四单元模型中完整推导散射回路与 Schur 补之间的关系，寻找比标量 CIE 更有针对性的反馈结构。然后在 24×24 高对比度数字、双岛和环状散射体上检验，比较真实全波调用和收敛结果。

**第二条：基于接收端有理响应的非局部材料探索。** 从式 (14) 出发，研究有限步的极点、零点和高阶响应。先比较 GN 的方向是否真的因为一阶信息失效而错过可恢复材料，再试验低阶有理近似能否发现新的候选路径。

两条都必须与已有的固定 \(\beta=0.2\) FBG、GlobalDual 和 Native Full-J 比较，且最终采用同一真实物理验收与计算账本。成熟 CIE-I/FFT-TSOM 也应补上，避免错误认定新颖性。

**如果只选一个公式让我继续深挖，我会选：**

\[
\boxed{
F(X_0+tD)-F(X_0)
=
tG_S
\left(I-tA_0^{-1}DG_D\right)^{-1}
A_0^{-1}DE_0
}
\]

因为 CIE 已经告诉我们：一个好的分式结构能够显著改变算法表现。

而上面这个公式更进一步——它不是人为构造的分式材料变换，而是**完整离散 Maxwell 方程本身沿材料方向具有的精确有理结构**。

这意味着我们有机会直接从散射的数学结构设计搜索机制，而不是不断在传统 GN 外面套新的启发式优化层。
