# 六根对偶内核的 quotient solver：一次 4×4 线性解与 6×6 特征分解

状态：精确代数推导，未运行实验，未读取新 private GT。本文研究如何求解已经导出的四条双线性方程；固定 Green、共享实材料、原 resolvent 可逆的物理合同不变。Hilbert–Burch、商代数与 multiplication matrix 都是标准工具；本文件不作最终新颖性或实测加速判断。

## 1. 从已有双线性系统形成 4×3 affine matrix

沿用 `algebraic_structure.md` 的一接收、四材料、两独立照明记号。测量零空间被固定参数化为

\[
w=w_0+cN,\quad v=v_0+cM,\quad
v_0=R+w_0H,\quad M=NH,
\]

其中 \(c=(x,y)\in\mathbb C^2\)、\(N,M\in\mathbb C^{2\times4}\)。把共轭参数独立化为 \(d=(d_1,d_2)\)，四条方程是

\[
F_i(c,d)=w_i(c)(\overline{v_{0i}}+d\overline M_{:,i})
-v_i(c)(\overline{w_{0i}}+d\overline N_{:,i})=0.
\tag{1}
\]

对每个 \(i=1,\ldots,4\)，定义

\[
\begin{aligned}
A_{i0}(c)&=w_i(c)\overline{v_{0i}}-v_i(c)\overline{w_{0i}},\\
A_{i1}(c)&=w_i(c)\overline{M_{1i}}-v_i(c)\overline{N_{1i}},\\
A_{i2}(c)&=w_i(c)\overline{M_{2i}}-v_i(c)\overline{N_{2i}}.
\end{aligned}
\tag{2}
\]

因此

\[
\boxed{F(c,d)=A(c)(1,d_1,d_2)^T,\qquad A(c)=A_0+xA_1+yA_2.}
\tag{3}
\]

这里 \(A_j\) 为 \(4\times3\) 常矩阵。复共轭只作用于已冻结系数，不作用于未知 \(x,y,d\)。矩阵 \(A\) 在本节是消元矩阵，不是原 Maxwell 状态 resolvent。

令

\[
g_i(x,y)=(-1)^{i+1}\det A_{\widehat i}(x,y),\qquad i=1,\ldots,4,
\tag{4}
\]

其中 \(A_{\widehat i}\) 删除第 \(i\) 行。四个 \(g_i\) 的总次数至多 3，且满足三条线性系数 syzygy

\[
A(c)^Tg(c)=0.
\tag{5}
\]

这是它们属于一个 determinantal ideal 的实质结构；把四个独立随机 cubic 代入下述算法不会自动得到相同结论。

## 2. 构造系数，不需要多项式拟合

采用列向量单项式

\[
B=(1,x,y,x^2,xy,y^2)^T,\qquad
C=(x^3,x^2y,xy^2,y^3)^T.
\]

唯一写成

\[
g=TC+QB,
\quad T\in\mathbb C^{4\times4},\quad Q\in\mathbb C^{4\times6}.
\tag{6}
\]

这里 \(Q\) 只是低次数系数矩阵，勿与状态 current 矩阵混淆。可直接用 3×3 行列式的列多线性得到全部系数：对每个删除行的子矩阵 \(L_j=(A_j)_{\widehat i}\)，

\[
g_i=(-1)^{i+1}\sum_{k_1,k_2,k_3\in\{0,1,2\}}
\det[L_{k_1}[:,1],L_{k_2}[:,2],L_{k_3}[:,3]]
x^{\#\{a:k_a=1\}}y^{\#\{a:k_a=2\}}.
\tag{7}
\]

这只有 27 个固定组合；可进一步按次数合并。四个最高次系数分别来自三个 \(L_1\) 列、两个 \(L_1\) 加一个 \(L_2\)、一个 \(L_1\) 加两个 \(L_2\)、三个 \(L_2\)。无需采样响应或拟合 cubic。

当 \(H,N,M,R\) 在同一切片内固定，而目标只通过 \(w_0\) 变化时，\(A_1,A_2\) 的第 1 列对共轭目标 affine，其余两列与目标无关。行列式对该列线性，所以最高次块 \(T\) 也对共轭目标 affine；可缓存对应模板。这一缓存不改变下面的可逆性要求。

## 3. 主定理：最高 cubic 系数块可逆即给出六维商代数

**定理。** 设 \(A(x,y)\) 为 \(4\times3\) affine-linear 矩阵，\(g_i\) 是它的四个带符号最大子式。若式 (6) 的 \(T\) 可逆，则

\[
\mathcal A=\mathbb C[x,y]/\langle g_1,g_2,g_3,g_4\rangle
\]

作为复向量空间恰有维数 6，且 \(B\) 的六个单项式是一个基。该结论按代数重数计数，不声称六个点互异。

**证明。** 引入齐次坐标 \(C_h=(c_0,c_1,c_2)\)，将矩阵齐次化为

\[
\mathsf A(C_h)=c_0A_0+c_1A_1+c_2A_2,
\]

最大子式 \(G_i\) 是三变量齐次 cubic。限制在 \(c_0=0\) 后，式 (6) 成为

\[
G(0,c_1,c_2)=T(c_1^3,c_1^2c_2,c_1c_2^2,c_2^3)^T.
\]

\(T\) 可逆意味着这些子式在无穷远直线没有共同项目零点。正维项目曲线必与任意项目直线相交，因此该 determinantal 零集至多零维。另一方面，4×3 最大子式的 determinantal height 上界为 2；齐次正次数理想是 proper，且 \(T\ne0\) 保证子式不全为零。因此其 height 恰为 2。

Hilbert–Burch 给出精确分解

\[
0\longrightarrow S(-4)^3
\xrightarrow{\mathsf A}
S(-3)^4\longrightarrow S\longrightarrow S/I\longrightarrow0,
\quad S=\mathbb C[c_0,c_1,c_2].
\tag{8}
\]

于是

\[
\operatorname{Hilb}_{S/I}(t)=
\frac{1-4t^3+3t^4}{(1-t)^3}
=\frac{1+2t+3t^2}{1-t}.
\tag{9}
\]

项目零维 scheme 的长度是 6。该 height-2 perfect ideal 是 Cohen–Macaulay，无嵌入的 irrelevant 分量；无穷远没有点意味着 \(c_0\) 是非零因子。故去齐次化得到长度 6 的 affine 商代数，并且没有在去齐次化后凭空出现低于三次的关系。式 (6) 的可逆 \(T\) 又将每个三次单项式写成 \(B\) 的组合，因此 \(B\) 张成且线性无关。证毕。

这也证明，下述四个三次 reduction 构成任意 degree-compatible 单项式序下的 Gröbner basis：它们的最高单项式恰好是全部四个三次单项式，而实际商维数等于六个 standard monomial 的数量。无需额外运行 Buchberger 才能得到这个结构结论。

## 4. 一次 4×4 solve 构造两个 6×6 multiplication matrix

求解六个右端的线性系统

\[
\boxed{L=-T^{-1}Q,\qquad C\equiv LB\pmod I.}
\tag{10}
\]

实际实现用线性 solve，不显式形成 \(T^{-1}\)。令 \(e_j\) 为六维标准基，按零起始索引 \(e_0,\ldots,e_5\) 对应 \(1,x,y,x^2,xy,y^2\)。以“多项式系数列向量”为约定，定义

\[
\boxed{
M_x=[e_1,e_3,e_4,L_{0,:}^T,L_{1,:}^T,L_{2,:}^T],}
\tag{11}
\]

\[
\boxed{
M_y=[e_2,e_4,e_5,L_{1,:}^T,L_{2,:}^T,L_{3,:}^T].}
\tag{12}
\]

它们分别实现商代数中的乘 \(x\) 和乘 \(y\)。因此精确满足

\[
M_xM_y=M_yM_x,\qquad g_i(M_x,M_y)=0.
\tag{13}
\]

### 4.1 不依赖 Hilbert–Burch 引用的代数自检

对任何四个最高次系数块可逆的 cubic，都可以造出式 (11)–(12)，但一般未必对易。如果它们**精确对易**，则已经足以证明这组 reduction 完整：

1. 单位元 \(e_0\) 自动 cyclic，因为 \(B(M_x,M_y)e_0\) 的六项恰为 \(e_0,\ldots,e_5\)。
2. 四条 cubic reduction 按构造都 annihilate \(e_0\)；对易性使它们 annihilate 每个 \(B_j(M_x,M_y)e_0\)，因此 annihilate 全空间。
3. 任意高次单项式不断使用 cubic reduction 都降到 \(\operatorname{span}B\)，故多项式商维数至多 6。
4. cyclic 的六维表示给出从该商到 \(\mathbb C^6\) 的满射，故维数至少 6。

于是该表示与商代数同构。对精确最大子式，主定理已经保证对易；浮点 commutator 是有用的构造/条件诊断，不能把“小”替代“精确为零”作为严格证书。

## 5. 该 Maxwell 配对结构排除无限或非唯一的 d-lift

对一般 affine matrix，最大子式只是 \(\operatorname{rank}A\le2\) 的条件，可能给出无穷远或多值的 \(d\)。本任务的共轭配对结构更强。

把式 (1) 齐次化成 \(F(C_h,D_h)\)，其中 \(D_h=(d_0,d_1,d_2)\)。系数构造直接给出反线性 involution

\[
\boxed{F(C_h,D_h)=-\overline{F(\overline{D_h},\overline{C_h})}.}
\tag{14}
\]

**引理。** 在式 (14) 成立且 \(T\) 可逆时，每个 determinantal 项目根 \(C_h\) 都有 \(c_0\ne0\)，矩阵 \(\mathsf A(C_h)\) 恰为 rank 2，且它唯一的项目核向量 \(D_h\) 满足 \(d_0\ne0\)。

**证明。** 主定理的开头已经排除全部 \(c_0=0\) 的 determinantal 根，也就排除了全部这样的 incidence 根。若存在 \(d_0=0\) 的 incidence 根，式 (14) 把它换成 \(c_0=0\) 的 incidence 根，矛盾。若某个根处 rank 至多 1，则项目核至少是一条直线，必与 \(d_0=0\) 的项目超平面相交，又得矛盾。因此 rank 恰为 2，项目核唯一且有限。证毕。

非零的 2×2 minor 给出局部有理反解，因此这不只是点集双射，而是 determinantal scheme 与原项目双线性 incidence scheme 的同构。结论：在**配对结构与精确 \(T\) 可逆**这两个条件下，原双线性系统恰有长度 6 的有限 affine scheme，没有隐藏的 c/d 无穷远根，也没有 rank 降低导致的额外 d-fiber。

该结论没有排除重根，也没有自动排除材料恢复节点或原 Maxwell 极点。

### 5.1 为何一般矩阵的反例不能省略配对假设

一般 Hilbert–Burch 矩阵

\[
A(x,y)=\begin{pmatrix}
y&0&0\\
-x&y&0\\
0&-x&y\\
0&0&-x
\end{pmatrix}
\]

的最大子式是 \(\pm(x^3,x^2y,xy^2,y^3)\)，所以 \(T\) 可逆，但 \(A(0,0)=0\)，原双线性方程在那里允许任意 \(d\)。determinantal 商仍是长度 6 的 fat point。此例不满足式 (14)，因此不反驳本任务的加强引理；它说明一般的六维商定理与原双线性系统的等价需要正确结构。

## 6. 从一次 6×6 eigensolve 恢复全部简单根

若六维商代数 reduced，即六个点互异，则存在一个分离线性形式

\[
\ell=x+\eta y,
\qquad M_\ell=M_x+\eta M_y.
\tag{15}
\]

对六个不同根，造成特征值碰撞的 \(\eta\) 至多来自 15 个根对的比值，因此泛型 \(\eta\) 分离全部点。\(M_\ell\) 的六个特征值是 \(\ell\) 在六个根的值。

**左右特征向量约定必须正确。** 因为式 (11)–(12) 对系数列向量作用，根 \((x_j,y_j)\) 的 evaluation vector

\[
b_j=(1,x_j,y_j,x_j^2,x_jy_j,y_j^2)^T
\]

满足

\[
b_j^TM_x=x_jb_j^T,\qquad b_j^TM_y=y_jb_j^T.
\tag{16}
\]

因此，应取 \(M_\ell^T\) 的**右特征向量**，归一化第零坐标为 1，然后读出 \(x_j=(b_j)_1,y_j=(b_j)_2\)。直接对 \(M_\ell\) 的右特征向量取坐标比一般错误。随后检查剩余三项确为 \(x_j^2,x_jy_j,y_j^2\)，以及两个共同特征方程和原始最大子式残差。

重复特征值有两种不同原因：所选 \(\eta\) 没有分离不同根，或者商代数本身 nonreduced。可以更换 \(\eta\) 诊断前者；非零重数不会被任意换 \(\eta\) 消除。重根可能只有一个特征向量，需要局部代数、deflation、严格包围或已有同伦端点处理；不能把六个数值特征值强行当成六个不同根。

## 7. 恢复 d、实材料与原物理验收

对每个 \(c_j=(x_j,y_j)\)，原配对理论保证 rank 2 且有限 lift。可解

\[
A(c_j)[:,1:3]\,d_j=-A(c_j)[:,0]
\tag{17}
\]

这个 \(4\times2\) 一致满列秩系统；数值实现使用 QR/SVD 并保存残差。等价地，取 \(A\) 右零向量 \(r\)，用 \(d=(r_1/r_0,r_2/r_0)\)。若使用 NumPy 式 SVD，右零向量是 `Vh[-1].conj()`，不能把 `Vh[-1]` 直接当列零向量。

实际物理根须满足

\[
d_j=\overline{c_j},\qquad
w_i=z_iv_i,\qquad
z_i=\frac{\operatorname{Re}(w_i\overline v_i)}{|v_i|^2},
\tag{18}
\]

以及 \(v_i\ne0\)、材料合法盒、原 resolvent 可逆和独立完整 forward 的残差要求。特征求解只完成代数候选枚举，不代替这些物理过滤。一个共轭实点若 \(w_i=v_i=0\)，其材料恢复存在连续自由度；它不是普通正则有限材料根。

## 8. 精确且便宜的测量零空间参数化

式 (1) 需要 \(w_0B=T_*\)、\(NB=0\) 精确成立。若希望避免在线 SVD 并保留有理系数解释，选两个源矩阵行 \(P\) 使 \(B_P\in\mathbb C^{2\times2}\) 可逆；其余行记为 \(F\)。按 \((P,F)\) 排序，置自由变量 \(c=w_F\)，则

\[
w_P=(T_*-cB_F)B_P^{-1},\qquad
w_{0,P}=T_*B_P^{-1},\quad w_{0,F}=0,
\]

\[
\boxed{N=[-B_FB_P^{-1}\mid I_2].}
\tag{19}
\]

这只需要缓存一个稳定的 2×2 solve。取 pivot 时考虑条件数；该坐标并不自动优于正交 SVD 坐标，但它在精确算术中严格保留线性测量等式。归档 binary64 算子可被当作精确有理数进行这样的代数计算。

如果使用浮点 SVD 取得 \(N,w_0\)，则 \(NB\approx0\)、\(w_0B\approx T_*\) 只是数值关系。对随后 rounded cubic 的严格根证书不能自动变成原冻结物理方程的严格证书；必须包括投影误差，或使用精确参数化。独立 forward 验收无论采用哪种坐标都保留。

## 9. 数值条件与成本边界

固定维数的主要在线工作是：构造四个 cubic 的系数、一次 4×4 solve（六右端）、一次 6×6 eigen、最多六次 4×2 backsub，以及材料/物理验收。它消除了逐路径 ODE/Newton continuation，但不因此自动证明实际总时间优势。

至少记录以下实际诊断：

- \(\kappa(T)\) 与线性求解相对 backward residual \(\|TL+Q\|/(\|T\|\|L\|+\|Q\|)\)；
- 相对 commutator \(\|M_xM_y-M_yM_x\|/(\|M_x\|\|M_y\|)\)；
- eigenpair residual、特征值最小分离度、特征向量矩阵条件数，以及 evaluation vector 第零分量的归一化稳定性；
- 原始最大子式、双线性方程、线性测量和复材料构成残差；
- 共轭误差、\(\min_i|v_i|\)、原 resolvent 稳定性、合法材料盒与独立 forward 结果。

这些指标的阈值应服从既有协议和误差尺度，不在看到结果后放宽。高 \(\kappa(T)\)、近重根和小 \(|v_i|\) 可造成明显放大；行/变量尺度平衡或更高精度有时有帮助，但“小残差”不能单独保证根位置误差小。

精确 \(T\) 不可逆时，本文的特定基 \(B\) 与直接 eigensolver 不适用；这不是原物理方程无根。应保留失败状态，按授权范围回退既有六路径同伦、另一个项目坐标图或一般 quotient basis 构造。普通 affine 坐标变换不能凭空提高最高次子式张成空间的秩。

`QUOTIENT_SIX_SIMPLE_NUMERICAL` 可以表示六个数值上分离的共同根和完整物理分类；它仍不是严格的精确系数全根证书。严格完整性需要认证 \(T\) 可逆、配对/投影恒等式、根及分离条件，并处理物理恢复域。没有证据时保留 `UNKNOWN`，不把一个 eigen call 的完成解释为这些门已经通过。

## 10. 结论与明确不成立的升级

对于本任务的共轭配对四条双线性方程，**精确最高 cubic 系数块 \(T\) 可逆**给出长度 6 的有限 affine incidence scheme、固定的二次单项式商基和正确的 multiplication matrices。若该 scheme reduced 且选取分离的 \(\ell\)，一次 6×6 特征分解恢复全部六个复候选，不必在线同伦。

该结论不保证候选全部为合法实材料，不保证 resolvent 非奇异，不保证在 branch collision 附近稳定，不保证成本超过强普通基线，也不把一个材料切片扩展为完整六材料域的覆盖。以上均需分别给出证据。
