# 一接收、2s 实材料的通用 quotient 内核

状态：条件定理与可执行公式；未运行实验，未读取新 private GT，未修改已冻结的四材料实现。结论是现有 determinantal / Eagon–Northcott / multiplication-matrix 工具对本配对系统的应用，不是一般代数定理的优先权主张。

## 1. 结论先行及精确条件

设有 \(2s\) 个独立共享实材料参数、一个复接收通道和 \(s\) 个独立复照明，\(s\ge1\)。精确消去线性测量后，对偶系统为 \(2s\) 条共轭配对双线性方程，未知 \(c,d\in\mathbb C^s\)。将它写成

\[
F(c,d)=A(c)(1,d_1,\ldots,d_s)^T,
\quad A(c)=A_0+\sum_{j=1}^s c_jA_j\in\mathbb C^{2s\times(s+1)}.
\tag{1}
\]

取 **全部** \((s+1)\times(s+1)\) 最大子式，其数量是

\[
K=\binom{2s}{s+1}.
\]

它们的最高次为 \(s+1\)。同次单项式 \(c^\alpha\)、\(|\alpha|=s+1\) 也恰有 \(K\) 个。令 \(T\in\mathbb C^{K\times K}\) 为这些最高次系数构成的方阵。

**主结论：若精确 \(T\) 可逆，则 determinantal 商代数维数恰为**

\[
\boxed{D=\binom{2s}{s},}
\tag{2}
\]

**一个基是全部次数不超过 \(s\) 的单项式。** 配对结构进一步保证：所有 determinantal 项目根都在 affine \(c\) 图内，每个根的 \(A\) 恰为 rank \(s\)，且唯一的 \(d\) 也有限。因而原双线性 incidence scheme 同样为长度 \(D\) 的有限 affine scheme。

若另外 reduced 且选取分离线性形式，则一次 \(K\times K\) solve（\(D\) 个右端）与一次 \(D\times D\) eigensolve 恢复全部 \(D\) 个简单复候选。重数、非共轭根、材料节点、原 resolvent 极点和合法盒必须分别处理。

对 \(s=3\)，这正是 **15 条 quartic、15×15 solve、20 维 quotient、20×20 eigensolve**。它是六材料/三独立照明的条件公式，不能直接沿用四材料实验的实测时间结论。

## 2. 固定 Green 问题如何产生这个矩阵

固定基点或 Schur 消元后，仍采用完整非线性响应

\[
Y-Y_0=R(I-D_zH)^{-1}D_zB,
\quad z\in\mathbb R^{2s},\quad
R\in\mathbb C^{1\times2s},\quad B\in\mathbb C^{2s\times s},
\tag{3}
\]

要求 \(\operatorname{rank}_{\mathbb C}B=s\)。定义接收传递行

\[
w=R(I-D_zH)^{-1}D_z,
\quad v=R+wH,
\quad wB=Y-Y_0,
\quad w_i=z_iv_i.
\tag{4}
\]

取满足线性测量的特解 \(w_0\) 及完整零空间行基 \(N\in\mathbb C^{s\times2s}\)：

\[
w_0B=Y-Y_0,\quad NB=0,\quad \operatorname{rank}N=s,
\quad w=w_0+cN,
\]

\[
v=v_0+cM,\quad v_0=R+w_0H,\quad M=NH.
\tag{5}
\]

当每个 cell 的材料为同一个实数时，对偶约束为

\[
F_i(c,d)=w_i(c)\bigl(\overline{v_{0i}}+d\overline M_{:,i}\bigr)
-v_i(c)\bigl(\overline{w_{0i}}+d\overline N_{:,i}\bigr)=0,
\quad i=1,\ldots,2s.
\tag{6}
\]

真实材料点要求 \(d=\overline c\)；复代数求解中 \(c,d\) 独立，只有固定系数取共轭。矩阵 \(A\) 的第零列与第 \(a\) 列分别是

\[
A_{i0}=w_i\overline{v_{0i}}-v_i\overline{w_{0i}},\qquad
A_{ia}=w_i\overline{M_{ai}}-v_i\overline{N_{ai}},\quad a=1,\ldots,s.
\tag{7}
\]

这些列均对 \(c\) affine。方程不需要 reciprocity 或几何对称，只需要同一实材料被所有照明共享。

## 3. 最大子式数量与一般系数配方

记 \(\mathcal I\) 为所有大小 \(s+1\) 的行子集，\(g_I(c)=\det A_I(c)\)。固定一致行序即可，不需要四材料情形的删除一行符号。

记

\[
\mathcal B=\{\beta\in\mathbb N^s:|\beta|\le s\},\qquad
\mathcal C=\{\alpha\in\mathbb N^s:|\alpha|=s+1\},
\]

\[
B_{\rm mon}=(c^\beta)_{\beta\in\mathcal B},\quad
C_{\rm mon}=(c^\alpha)_{\alpha\in\mathcal C}.
\]

计数给出 \(|\mathcal B|=D\)、\(|\mathcal C|=K\)，所以

\[
\boxed{g=TC_{\rm mon}+QB_{\rm mon},\qquad Q\in\mathbb C^{K\times D}.}
\tag{8}
\]

系数可直接由行列式多线性计算。对每个 \(I\)，令 \(L_j=(A_j)_I\)，约定 \(c_0=1\)。则

\[
g_I(c)=\sum_{j_1,\ldots,j_{s+1}=0}^s
\det[L_{j_1}[:,1],\ldots,L_{j_{s+1}}[:,s+1]]
\prod_{a=1}^{s+1}c_{j_a}.
\tag{9}
\]

按多指标归并即可获得 \(T,Q\)，不需要拟合响应。式 (9) 是精确配方，不是高维时最优的系数装配算法；预计算、稀疏卷积和装配成本都应计入实测账本。

同一切片内固定 \(H,R,N,M\) 时，最高次数部分的第零列对共轭目标 affine，其余列固定。因此 \(T\) 对共轭目标 affine，可缓存模板。这不意味着全部低次系数也只对目标 affine。

## 4. 主定理的完整证明

### 4.1 最高次可逆排除项目无穷远并确定 height

在 \(S=\mathbb C[c_0,c_1,\ldots,c_s]\) 中齐次化

\[
\mathsf A(C)=c_0A_0+\sum_{j=1}^s c_jA_j,
\quad I=I_{s+1}(\mathsf A).
\]

\(I\) 由全部最大子式 \(G_I(C)\) 生成。令 \(c_0=0\)，最高次系数关系变成

\[
(G_I(0,c))_I=T(c^\alpha)_{|\alpha|=s+1}.
\tag{10}
\]

若 \(T\) 可逆，右边只能在 \(c=0\) 时全为零。因此 \(\operatorname{Proj}(S/I)\) 与超平面 \(c_0=0\) 不相交。任何正维项目子簇都与超平面相交，所以项目零集至多零维。

另一方面，一个 \(2s\times(s+1)\) 矩阵最大子式的 determinantal height 上界是

\[
2s-(s+1)+1=s.
\]

齐次正次数理想 proper，故不可能只在原点有 height \(s+1\) 的零集。结合上面的零维项目支持，得到

\[
\boxed{\operatorname{ht}I=s.}
\tag{11}
\]

这里不再另加“泛型矩阵”假设：精确 \(T\) 可逆已经强到足以推出所需 expected height。

### 4.2 Eagon–Northcott 给出 Cohen–Macaulay 性

在 expected height \(s\) 下，最大子式的 Eagon–Northcott complex 是自由分解，因此 \(S/I\) 是维数 1 的 Cohen–Macaulay 环。其齐次自由分解的 Betti 数为

\[
\beta_i=\binom{2s}{s+i}\binom{s+i-1}{i-1},
\quad i=1,\ldots,s,
\tag{12}
\]

相应生成次数为 \(s+i\)。因此

\[
\operatorname{Hilb}_{S/I}(t)=
\frac{1+\sum_{i=1}^s(-1)^i\beta_i t^{s+i}}{(1-t)^{s+1}}.
\tag{13}
\]

更直接的长度推导是：Cohen–Macaulay 性排除嵌入的 irrelevant 分量，而无穷远无点保证 \(c_0\) 不属于任何 associated prime。所以 \(c_0\) 是非零因子。

### 4.3 模掉 c0 后是明确的截断多项式环

由式 (10) 的可逆系数变换，

\[
\frac{S/I}{c_0(S/I)}
\cong\frac{\mathbb C[c_1,\ldots,c_s]}{(c_1,\ldots,c_s)^{s+1}}.
\tag{14}
\]

因此利用 \(c_0\) 的非零因子精确序列，

\[
\boxed{\operatorname{Hilb}_{S/I}(t)=
\frac{\displaystyle\sum_{j=0}^s\binom{s+j-1}{j}t^j}{1-t}.}
\tag{15}
\]

长度即分子在 \(t=1\) 的值：

\[
\sum_{j=0}^s\binom{s+j-1}{j}=\binom{2s}{s}=D.
\tag{16}
\]

\(S/I\) 作为 \(\mathbb C[c_0]\) 模可由 \(c^\beta\)、\(|\beta|\le s\) 自由生成；这是式 (14)、graded Nakayama 和相同 Hilbert series 的直接结果。令 \(c_0=1\) 后，affine 商仍为维数 \(D\)，而这些单项式成为一个基。证毕。

### 4.4 s=3 的可核对特例

此时 Betti 数为 \((15,24,10)\)，自由分解为

\[
0\to S(-6)^{10}\to S(-5)^{24}\to S(-4)^{15}\to S\to S/I\to0,
\]

\[
\operatorname{Hilb}_{S/I}(t)=
\frac{1-15t^4+24t^5-10t^6}{(1-t)^4}
=\frac{1+3t+6t^2+10t^3}{1-t}.
\tag{17}
\]

这给出 \(1+3+6+10=20\) 个基单项式，而不是把 15 条 quartic 当成相互独立的一般四次方程使用总次数 Bézout。

## 5. 共轭交换使投影消元保持完整 scheme

令 \(C=(c_0,c)\)、\(D_h=(d_0,d)\)，并将式 (6) 双齐次化。其系数给出

\[
\boxed{F(C,D_h)=-\overline{F(\overline{D_h},\overline C)}.}
\tag{18}
\]

这是复系数的反线性配对恒等式，不是物理群作用。

若 \(T\) 可逆，已经没有 \(c_0=0\) 的 determinantal 根。若存在 \(d_0=0\) 的 incidence 根，式 (18) 会将其交换为一个 \(c_0=0\) 的 incidence 根，矛盾。

对任一 determinantal 根，\(\operatorname{rank}\mathsf A\le s\)。如果其 rank 至多 \(s-1\)，核的向量维数至少为 2；相应项目核至少是一条直线，必与超平面 \(d_0=0\) 相交，又产生矛盾。因此

\[
\boxed{\operatorname{rank}\mathsf A(C)=s,\quad
\dim\ker\mathsf A(C)=1,\quad d_0\ne0.}
\tag{19}
\]

不仅每点有唯一有限 lift，scheme 也相同：在某个根附近选择尾部 \(s\) 列的非零 \(s\times s\) 子式，局部可逆行列变换将矩阵写成

\[
\begin{pmatrix}0&I_s\\ r&0\end{pmatrix},
\]

其中 \(r\) 有 \(s\) 个分量。最大子式理想局部就是 \(\langle r_1,\ldots,r_s\rangle\)，而 incidence 方程先唯一解出 \(d\)，余下方程也是 \(r=0\)。所以消元投影是 scheme 同构，长度和重数均保留。

结论不需要额外假设每个根的 rank；它由可逆 \(T\) 加配对恒等式推出。若配对恒等式被近似、截断或独立系数扰动破坏，则只能保留 determinantal 商定理，不能自动保留本节的有限唯一 lift。

## 6. 一次 coefficient solve 与 s 个 multiplication matrix

由式 (8) 求解

\[
\boxed{L=-\operatorname{solve}(T,Q),\qquad
c^\alpha\equiv\sum_{\beta\in\mathcal B}L_{\alpha\beta}c^\beta,
\quad |\alpha|=s+1.}
\tag{20}
\]

对每个变量 \(c_j\)，定义 \(D\times D\) 矩阵 \(M_j\)，它的第 \(\beta\) 列为

\[
\boxed{(M_j)_{:,\beta}=
\begin{cases}
e_{\beta+e_j},&|\beta|\le s-1,\\
L_{\beta+e_j,:}^{T},&|\beta|=s.
\end{cases}}
\tag{21}
\]

这里外侧 \(e_\gamma\) 表示单项式基对应的标准向量，内侧 \(e_j\) 表示多指标增量。矩阵对多项式**系数列向量**作用。

精确代数给出

\[
[M_i,M_j]=0,\quad g_I(M_1,\ldots,M_s)=0.
\tag{22}
\]

全部最高次单项式已经降到次数不超过 \(s\)，且商维数正好是 \(D\)，因此这些 normalized reduction 是任意 degree-compatible 顺序下的 Gröbner basis。高阶 S-polynomial 不会引入额外低次关系。

反过来，即使不调用 Eagon–Northcott 代码，精确的两两 commutation 加上述构造也足以证明完整商表示：单位元 cyclic，所有 degree-\(s+1\) reduction 按构造 annihilate 单位元；对易性使它们 annihilate 全空间；reduction 给维数上界 \(D\)，cyclic 表示给维数下界 \(D\)。浮点 commutator 小只是诊断，不是这项精确等式的认证。

## 7. joint eigenvalue 恢复及六材料的明确尺寸

若 scheme reduced，取一般线性形式 \(\ell(c)=\sum_j\eta_jc_j\)，令

\[
M_\ell=\sum_j\eta_jM_j.
\]

不同根被某个 \(\ell\) 混淆的 \(\eta\) 只落在有限多个真超平面上。因此一般 \(\eta\) 分离全部 \(D\) 个点。根 \(c^{(a)}\) 的 evaluation vector \(b^{(a)}=(c^{(a)\beta})_\beta\) 满足

\[
(b^{(a)})^TM_j=c_j^{(a)}(b^{(a)})^T.
\tag{23}
\]

所以取 \(M_\ell^T\) 的右特征向量，归一化常数分量为 1，再读取一次单项式分量。须验证其余分量确为对应单项式，并验所有共同特征方程与原始子式。

\(s=3\) 时可以按总次数顺序使用

\[
1;\quad x,y,z;\quad x^2,xy,xz,y^2,yz,z^2;
\]

\[
x^3,x^2y,x^2z,xy^2,xyz,xz^2,y^3,y^2z,yz^2,z^3.
\]

degree-4 border 有 15 项，故主系数 solve 是 15×15、20 右端；\(M_x,M_y,M_z\) 均为 20×20。然后最多对 20 个候选解

\[
A(c)[:,1:s+1]d=-A(c)[:,0],
\tag{24}
\]

即 6×3 一致满列秩 backsub。

若特征值重复，应区分所选 \(\ell\) 不分离与 nonreduced 真重根；更换 \(\eta\) 只能处理前者。后者仍需局部代数、deflation 或已有同伦/认证方法，不能强读成 20 个不同点。

## 8. 缺失假设、反例与成本边界

**可逆 T 是充分条件，不是所有可逆材料问题都满足的条件。** 纯代数例取

\[
H=0,\quad R=(1,\ldots,1),\quad
B=\begin{pmatrix}I_s\\ iI_s\end{pmatrix}.
\]

响应是 \(Y_j=z_j+iz_{s+j}\)，所以实材料有全局唯一线性反解，resolvent 为单位阵；但双线性系数退化为 affine，最高次子式块 \(T=0\)。这证明 quotient chart 失败不等于物理方程无解或不唯一。此例是代数边界反例，未声称该 \(H,B,R\) 来自本轮真实 Maxwell Green 几何。

**source rank 不能省略。** 若实际 \(\operatorname{rank}B<s\)，零空间维数变大，式 (1) 的平衡尺寸不再成立。SVD 丢弃非零奇异值会改变问题，不能作为本定理的精确化简。

**T 可逆不保证 reduced。** 定理计入代数重数。非共轭候选、材料恢复节点和原 resolvent 极点仍可能存在；合法实材料只能在 \(d=\overline c\)、\(v_i\ne0\)、\(w_i=z_iv_i\) 和原完整 forward 验收后接受。

**浮点条件数仍可能主导成本与可靠性。** 应记录 \(\kappa(T)\)、coefficient solve backward error、所有 commutator、eigen separation/conditioning、lift 与 constitutive 残差和独立物理残差。投影恒等式与共轭交换若只近似成立，也必须计入严格证书的误差；不能只认证最终 rounded polynomial。

**维度增长仍是组合量级。**

| s | 材料数 2s | 系数 solve 尺寸 K | quotient/eigen 尺寸 D |
|---:|---:|---:|---:|
| 1 | 2 | 1 | 2 |
| 2 | 4 | 4 | 6 |
| 3 | 6 | 15 | 20 |
| 4 | 8 | 56 | 70 |
| 5 | 10 | 210 | 252 |
| 6 | 12 | 792 | 924 |

\(D\sim4^s/\sqrt{\pi s}\)，且 \(K/D=s/(s+1)\)。这比总次数或 factorial 路径界紧，但并非大规模材料成像的多项式复杂度算法。在线系数装配、预计算摊销、失败/高精度回退和物理验收均应计费。

## 9. 先例与本文件的归属边界

最高次最大子式张成完整齐次形式空间并非新的一般结构。Ito、Noma、Ohno 的 *Maximal minors of a matrix with linear form entries*（2015，在线 2014）明确研究：齐次线性矩阵在非零点处始终满秩时，其最大子式张成相应次数的全部齐次形式空间。这里的 \(T\) 可逆条件正落在这类结构中。本轮核验到出版社摘要和书目信息，未宣称已取得该文完整证明全文。[出版社页面](https://www.tandfonline.com/doi/full/10.1080/03081087.2014.959516)

Gopalakrishnan、Neiger、Safey El Din 的 *Optimized Gröbner basis algorithms for maximal determinantal ideals and critical point computations*（2024）也明确使用 Eagon–Northcott 分析最大子式的 syzygy 并改进 Gröbner 计算。它的研究对象和算法不能自动等同本文件的物理问题，但足以排除“利用 Eagon–Northcott 求 determinantal 系统本身是新方法”的说法。[作者预印本](https://arxiv.org/abs/2402.07353)

本文件的工作是写清本任务的平衡对偶参数化、共轭配对 lift、固定 monomial quotient 和验收边界。一般代数定理有先例；物理适用性、实测收益与最终新颖性由独立证据和主线程审查决定。
