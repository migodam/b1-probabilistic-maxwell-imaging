# 推导公式集：从固定 Green 到可计算的材料根

维护原则：只保留推导得到、说明假设、可以改变算法的公式。每条标明“标准工具”“本轮推导组合”或“待证”。没有把重新命名的算子作为成果。本文讨论固定几何、独立标量材料的有限离散模型；不声称连续PDE或三维矢量Maxwell定理。

## F01　Born双Green积分，以及它没有包含的多次散射

标量均匀背景Green满足辐射条件。设入射点源在 r_s，接收在 r_r，材料对比度 χ：

\[
 E_s^{\rm sc}(r_r)=k_b^2\int_\Omega g_b(r_r,x)\chi(x)E_s(x)\,dx,
\quad E_s(x)=g_b(x,r_s)+k_b^2\int_\Omega g_b(x,z)\chi(z)E_s(z)\,dz.
\]

用入射场替代第二式的总场得到Born：

\[
 E_{rs}^{\rm B}=k_b^2\int_\Omega g_b(r_r,x)\chi(x)g_b(x,r_s)\,dx,
\qquad Y_B=S D(m)E_0.
\]

离散权重、k_b²和基函数积分吸收到G、S内。G、S、E0由几何/背景固定；材料仅在D中。**标准工具。** 本轮最终候选用完整非线性模型验收，不把Born拟合当物理解。

## F02　远场是幅度与相位的共同渐近式

三维标量 \(g(R\hat r,x)=e^{ik|R\hat r-x|}/(4\pi|R\hat r-x|)\)。对 \(|x|\le a\)，

\[
 |R\hat r-x|=R-\hat r\cdot x+O(a^2/R),\qquad
 g(R\hat r,x)=\frac{e^{ikR}}{4\pi R}e^{-ik\hat r\cdot x}
 [1+O(a/R)+O(ka^2/R)].
\]

因此要同时有 \(a/R\ll1\)、\(ka^2/R\ll1\)；只提出1/r不够。二维本轮模型 \(g=(i/4)H_0^{(1)}(k|r-x|)\)，远场幅度是 \(R^{-1/2}\)，不是 \(R^{-1}\)。矢量Maxwell还带横向投影/并矢结构。**标准渐近工具，本轮不以其替代完整Hankel算子。**

## F03　精确VIE与可复用的Green矩阵

令 j=D(m)E，E=E0+Gj，则

\[
 M(m)j=D(m)E_0,\quad M(m)=I-D(m)G,\quad
 Y=S M(m)^{-1}D(m)E_0.\tag{3}
\]

仅假定M可逆，不要求Born收敛。若额外有谱半径ρ(DG)<1，才能用Neumann级数展开；本轮求解器不依赖此条件。

圆盘常数pulse、半径a、单元中心距离d的实际代码系数为

\[
 G_{ij}=\frac{i\pi ka}{2}J_1(ka)H_0^{(1)}(kd_{ij})\;(i\ne j),\quad
 G_{ii}=\frac{i\pi ka}{2}H_1^{(1)}(ka)-1.
\]

外部S取同一圆盘积分互项；E0点源取(i/4)H0。有限pulse的严格代数不等于连续介质离散误差认证。

## F04　固定采集SVD、投影与TSOM表达的边界

S=UΣV*，保留前r个非零奇异值，则

\[
 j_d=V_r\Sigma_r^{-1}U_r^*Y,\quad
 P_r=V_rV_r^*,\quad Sj_d=U_rU_r^*Y.
\]

保留全部非零奇异值且数据相容时，所有观测相容电流是 j=j_d+N c，其中SN=0。截断后j_d是正则化近似，丢弃通道残差不可称精确零空间。**SVD、确定/模糊电流来自SOM/TSOM已有结构。**

两个子空间可以用投影P、Q计算。顺序作用PQ容易实现，但一般不等于交集投影：只有正交投影可交换等额外条件下PQ才是交集正交投影。若要在第一个子空间内提取第二结构，可对QB或B*QB分解；这不是宣称两个子空间有未经计算的精确交集。有限维交替投影极限可给交集，但速率由夹角决定。本轮不用任意“twofold”命名替代此检验。

构成残差也能直接暴露截断错误。若j=j_d+j_a为真电流，则

\[
 j_d-D(E_0+Gj_d)=-(I-DG)j_a.
\]

弱接收电流不一定是弱内部材料效应。

## F05　充分观测下的精确非线性直接逆及条件数

**命题（标准线性逆+本构恢复）。** S满列秩，Y精确来自式(3)，且第i单元跨照明场不全零，则

\[
 J=S^\dagger Y,\quad E=E_0+GJ,\quad
 m_i=\frac{\operatorname{Re}\sum_s\overline{E_{is}}J_{is}}
 {\sum_s|E_{is}|^2}.\tag{5}
\]

证明：左逆唯一恢复J，然后J_is=m_i E_is；按行内积即可。该表达保留多次散射，不是Born。若所有E_is为零，所有J_is也为零，在非极点邻域改变m_i仍可保留同一J，不能唯一恢复该坐标。

噪声η经h=S†η传到电流，E'=E+Gh。行拟合估计满足精确误差恒等式

\[
 \widehat m_i-m_i=
 \frac{\operatorname{Re}\langle E'_i,h_i-m_i(Gh)_i\rangle}{\|E'_i\|^2},\qquad
 |\widehat m_i-m_i|\le\frac{\|h_i-m_i(Gh)_i\|}{\|E'_i\|}.
\]

因此小奇异值、强内部传播和弱局部照明都会放大误差。满秩不是鲁棒成像保证。若S不满列秩但E0满行秩，也有对偶直接逆：W=Y E0†，由 W(I-GD)=SD 得 W=(S+WG)D，按列同样恢复m。此式不要求G对称。**直接逆不列为原创性主张。**

## F06　固定活跃材料集合的精确有理闭包

U必须是选择独立材料单元的坐标列；一般decoder稠密方向不能冒充这种低秩更新。设D=D0+UΔUᵀ，Δ=diag(δ)，M0可逆，V=UᵀG。

从两状态方程相减得到

\[
 M_0(j-j_0)=U\Delta[U^TE_0+Vj_0+V(j-j_0)].
\]

令 H=M0^{-1}U，K=VH，q=UᵀE0+Vj0。写j−j0=Hz，得

\[
 (I-\Delta K)z=\Delta q,\quad
 Y-Y_0=L(I-\Delta K)^{-1}\Delta q,\quad L=SH.\tag{6}
\]

矩阵行列式引理同时给出

\[
 \det M=\det M_0\det(I-\Delta K).
\]

M0固定可逆时两端奇异性完全等价。Δ可含零，不需要Δ^{-1}。setup支付一次M0分解及照明+活跃列RHS。每候选只求k×k系统；完整状态只在最终验收时新求。**SMW/Schur标准恒等式；本轮价值必须来自其后的根结构和成本。**

## F07　行列式多重仿射、数据仿射系数与饱和

每个m_i只进入M及右端D e_s的第i行，故

\[
 P_{rs}(m;y)=\det\begin{bmatrix}M(m)&D(m)e_s\\S_r&y_{rs}\end{bmatrix}
 =y_{rs}\det M-S_r\operatorname{adj}(M)D e_s
\]

对每个m_i次数≤1。Schur补证明在detM≠0时P=0等价真实观测。

\[
 \det(I-DG)=\sum_{A\subseteq\{1,\ldots,n\}}
 (-1)^{|A|}\det G_{A,A}\prod_{i\in A}m_i.
\]

空子式取1。各观测的分母系数相同、分子系数固定，所以新y只更新系数的仿射组合；2^k系数展开只适合小活跃块。

对实m取ReP、ImP，再把其**实系数多项式**复化。不能在复变量上直接调用数值real/imag，因为那不是全纯多项式。设D+为原复系数行列式、D−为系数逐项共轭的行列式，则合法复化理想使用

\[
 I_{\rm sat}=\langle P_{\rm Re},P_{\rm Im}\rangle:(D_+D_-)^\infty.
\]

实物理解上D−=conj(D+)；复根上一般不同。只过滤D+可能留下另一套分母伪根。**多项式化已有Everett等先例，不是原创。**

## F08　两变量直接消元与四变量24路径

两实未知x,z的两个方程p_l=a_l xz+b_l x+c_l z+d_l。消去z得到

\[
 (a_1x+c_1)(b_2x+d_2)-(a_2x+c_2)(b_1x+d_1)=0,
\]

至多二次，不是四次。恢复z前需单独处理a_l x+c_l=0，保留可能的垂直分支/伪根。

k个独立实材料、k个多重仿射方程在(P¹)^k的孤立复根数计重上界为k!。标准线性乘积起始系统g_i=∏_j(x_j−a_ij)，每个变量被一条方程选中，k!个置换给全部泛型起点。γ同伦H=(1−t)γg+tP在t<1泛型正则；终点仍可能奇异、在无穷远或分母基点。**4材料24路径是标准多齐次工具的应用，不是已认证24个物理解。**

## F09　关键表示：测量先消元的对偶双线性材料系统

一接收，k活跃材料，s独立照明，固定核(6)中q∈C^{k×s}满列秩。定义实际接收传递行

\[
 w=L(I-\Delta K)^{-1}\Delta.
\]

由(I−ΔK)^{-1}Δ=Δ(I−KΔ)^{-1}推得

\[
 wq=t=Y-Y_0,\qquad w=(L+wK)\Delta.\tag{9a}
\]

这是一条接收传递行，**不是每照明分别引入完整电流**。测量线性方程先精确消去：

\[
 w(c)=w_0+c^TN,\quad w_0q=t,\quad Nq=0,\quad c\in\mathbb C^{k-s},\quad
 v(c)=L+w(c)K.
\]

N由固定q的SVD一次得到；不同t只改变w0，N、NK不变。实δ共享于全部照明的必要且在v_i≠0时充分条件是

\[
 f_i(c,d)=w_i(c)\,\overline v_i(d)-v_i(c)\,\overline w_i(d)=0,
 \quad i=1,\ldots,k.\tag{9b}
\]

横线只共轭系数，d为独立变量；实际物理解要求d=conj(c)。每式是两组变量的双仿射，齐次化后bidegree(1,1)。恢复δ_i=w_i/v_i，再检查实性、材料box、两分母、完整Maxwell。

**推导组合，原创性待文献裁决。** 它把测量约束和共享实材料先消元，而非扩大current零空间后再做局部回缩。

## F10　平衡采集的中心二项式根界与六起点

**命题（有条件、已推导）。** k=2s，rankq=s；在v_i及两复化resolvent分母非零域，若目标复化饱和系统零维，则非奇异物理材料根可嵌入(P^s×P^s)中2s个(1,1)方程的孤立根集，故

\[
 \#\{\text{孤立复根，计重}\}\le [h_c^s h_d^s](h_c+h_d)^{2s}
 =\binom{2s}{s}.\tag{10}
\]

更局部地，上界适用于对应于孤立复零点的根；物理实根孤立**本身不够**推出复化孤立。可用满秩实材料Jacobian作为正则充分条件。4材料2照明时上界6，优于4!=24；6材料3照明时20对720。无宣称每个Maxwell系统实际达到此界。

起始系统g_i=A_i(c)B_i(d)，每个A、B为泛型仿射线性式。选s条方程令A=0、其余s条令B=0，分别解两个s×s线性系统，总计binom(2s,s)起点。标准γ同伦在此表示上执行。

若物理解v_i=0，则w_i=δ_i v_i=0；保持w不变并独立改变δ_i，式(9a)仍成立，且可逆性在小邻域保持，所以形成真实材料连续族。这样不能是独立坐标内的孤立物理解。此论证不替代复零维前提，也不排除box边界或受限decoder的特殊情况。

## F11　固定核在多个材料切片之间的Schur复用

把(6)的活跃单元分成固定B和待求A，固定δ_B。分块状态z_B=Δ_B(q_B+K_BA z_A+K_BB z_B)。若T_B=I−Δ_B K_BB可逆，

\[
 z_B=T_B^{-1}\Delta_B q_B+T_B^{-1}\Delta_B K_{BA}z_A,
\]
\[
 K_{\rm eff}=K_{AA}+K_{AB}T_B^{-1}\Delta_BK_{BA},\quad
 q_{\rm eff}=q_A+K_{AB}T_B^{-1}\Delta_Bq_B,
\]
\[
 L_{\rm eff}=L_A+L_BT_B^{-1}\Delta_BK_{BA},\quad
 Y_{\rm base}=Y_0+L_BT_B^{-1}\Delta_Bq_B.
\]

于是Y=Y_base+L_eff(I−Δ_A K_eff)^{-1}Δ_A q_eff。这避免每个切片重做大型M0分解；若T_B奇异需换参考/分块，不能略去。**标准Schur推导，后续实现状态单列。**

## F12　单坐标的精确有理最小化（已推导，未作为主要实验）

固定其它材料，只改一个实增量δ，由(6)有Y(δ)=Y0+aδ/(1−κδ)，a=Lq。令r=Y0−Yobs、b=a−κr，则

\[
 \Phi(\delta)=\frac{\|r+b\delta\|^2}{|1-\kappa\delta|^2}
 =\frac{A\delta^2+B\delta+C}{D\delta^2+E\delta+1},
\]

A=||b||²，B=2Re〈r,b〉，C=||r||²，D=|κ|²，E=−2Reκ。其导数分子严格简化为

\[
 (AE-BD)\delta^2+2(A-CD)\delta+(B-CE).
\]

因此一个有界坐标的全局线搜索只需二次实根、区间端点和极点排除，不需多次完整forward。多坐标迭代仍可能陷于局部极小；本轮不以这条标准有理最小化创造全局枚举主张。

## F13　Jacobian、计费与覆盖声明

完整非线性微分给出

\[
 \partial_{m_i}Y_{rs}=(SM^{-1})_{ri}E_{is}.
\]

一次forward需s个source RHS，一次共享LU下接收转置伴随需r个RHS；也可选择正向列推导，成本应按实现记录。重复验证、准备及失败路径都计价。

\[
 T(N)=T_{\rm setup}+N T_{\rm target},\quad
 N_* =\left\lceil\frac{T_{\rm setup,A}-T_{\rm setup,B}}
 {T_{\rm target,B}-T_{\rm target,A}}\right\rceil
\]

只在分母正且成本模型对应同任务时定义有意义交叉点；若A同时setup更低可从N=1受益。不能仅用RHS比值替代wall。

全部路径数值结束 ≠ 认证全部根；切片上的全部候选 ≠ 全材料域；多个切片材料点 ≠ 多个连通分支；参数闭环返回同根 ≠ 非平凡monodromy。需要证书时必须进一步认证零维、无穷远、重根、分母基点及路径完整性。

## 依据与归因

- [Everett1996，JCP](https://doi.org/10.1006/jcph.1996.0070)：小Helmholtz多项式逆问题与homotopy先例，当前可访问证据级别见audit/novelty.md。
- [Chen/Zhong/Agarwal，2010，SOM/TSOM综述](https://archive.intlpress.com/site/pub/files/_fulltext/journals/maa/2010/0017/0004/MAA-2010-0017-0004-a006.pdf)：固定SVD与确定/模糊电流已有先例。
- [Verschelde，Homotopy Methods教程](https://homepages.math.uic.edu/~jan/tutorial.pdf)：多齐次Bézout、linear-product start、γ技巧、参数同伦与witness均为标准数学工具。

F09–F10的特定共享实材料/接收对偶组合是本轮候选数学机制；“第一次”“普适更快”“全部材料域遍历”均未建立。

## F14　从四个三次式直接得到6×6商代数：本轮主要算法

令F09的c=(x,y)，将四个双线性方程整理成

\[
 A(x,y)\begin{bmatrix}1\\d_1\\d_2\end{bmatrix}=0,
 \qquad A_{ib}=C_{i0b}+xC_{i1b}+yC_{i2b},\quad A\in\mathbb C^{4\times3}.
\]

存在非零齐次d当且仅当四个3×3 maximal minors全部为零。各minor为三次多项式。按

\[
 B=(1,x,y,x^2,xy,y^2)^T,\qquad H=(x^3,x^2y,xy^2,y^3)^T
\]

排列系数，四式共同写成TH+RB=0。**T可逆时，直接计算Q=−T^{-1}R即可得到H=QB。**

### 命题及证明机制

假设系数精确、四式确实来自同一4×3线性矩阵的maximal minors、detT≠0。把x,y齐次化到P²，T可逆使无穷远上四个minor张成全部二元三次式，因此不存在无穷远公共零点。maximal-minor理想高度为2；标准Hilbert–Burch分解给出

\[
 0\to S(-4)^3\to S(-3)^4\to S\to S/I\to0,
 \qquad \operatorname{Hilb}_{S/I}(t)=\frac{1+2t+3t^2}{1-t}.
\]

于是仿射商环维数严格为6，B是一组基；这里长度计重，不等于六个不同点。

列为“乘以变量后在B下的系数”，构造

\[
 M_x=[e_1,e_3,e_4,Q_0^T,Q_1^T,Q_2^T],\qquad
 M_y=[e_2,e_4,e_5,Q_1^T,Q_2^T,Q_3^T].
\]

下标从0计；例如e1表示x的坐标。精确代数中MxMy=MyMx。取分离线性式ℓ=x+ηy，求(Mx+ηMy)^T的右特征向量u；若根简单、η确实分离，则u与(1,x,y,x²,xy,y²)^T成比例。按u0归一化或联合特征值恢复x,y，再回代A求d。**这是标准商环乘法矩阵求根原理在特定材料消元上的具体实现。**

### Maxwell共轭结构为何排除一个额外的提升漏洞

一般4×3线性矩阵即使T可逆，也可能在某投影根rankA≤1，产生无穷多d或仅d0=0的提升。不能把一般情况说成六个双线性解。

本问题却有额外的、由实材料产生的恒等式

\[
 f_i(c,d)=-\overline{f_i(\bar d,\bar c)}.
\]

因此任一d在无穷远的射影零点会经共轭交换成为c在无穷远的零点，与T可逆矛盾。若rankA≤1，其至少二维核必与d0=0相交，也矛盾。故在这一**保持精确共轭配对**的结构下，T可逆还保证每个投影根唯一提升到有限d。仍须另行处理重根、材料节点、resolvent极点与实物理box；它们没有被该引理消除。

### 计算后果和停止条件

新测量的在线工作变为：更新w0及9系数→4个小行列式的系数→4×4线性求解→6×6特征分解→有限回代→合法材料过滤→完整Maxwell验收。无需沿同伦参数走六条路径。

T近奇异、重根/特征值碰撞、不稳定回代或多项式残差未达预定数值资格时，代码转入已计费的六路径fallback或保留UNRESOLVED；不删去失败候选。浮点交换子小不是严格证书。源码：`src/quotient_solver.py`；精确推导和反例：`theory/quotient_solver.md`、`audit/quotient_review.md`。

**归因边界：** Hilbert–Burch、结式、商环特征值法、双线性求根均为成熟数学；本轮候选贡献在特定Maxwell实材料对应、固定核可复用系数、共轭提升条件、实际不同材料恢复与总成本证据，优先权仍须看直接先例。

## F15　平衡采集的行列式商公式族：6、20、70、……

F14的4材料机制推广为2s材料、s独立源、1复接收。A(c)为2s×(s+1)仿射矩阵，全部最大子式g有K=binom(2s,s+1)个。令H为s变量的全部s+1次单项式，B为全部次数≤s的单项式：

\[
 g=T H+R B,\quad T\in\mathbb C^{K\times K},\quad
 D=\dim B=\binom{2s}{s},\qquad H=-T^{-1}RB.
\]

T可逆时，Eagon–Northcott给出CM商；模掉齐次变量后的Hilbert级数是截断多项式环，所以B确为商基而非任意坐标命名。每个乘c_j算子只需上述最高次归约：

\[
 M_j[:,\alpha]=
 \begin{cases}
 e_{\alpha+e_j},&|\alpha|<s,\\
 (-T^{-1}R)[\alpha+e_j,:]^T,&|\alpha|=s.
 \end{cases}
\]

在简单根和分离条件下，\(\big(\sum_j\eta_jM_j\big)^T\) 的共同特征向量给全部c。反Hermitian兼容系数的共轭交换再保证有限唯一d提升，之后才过滤物理合法性。

**算法后果：** s=2用4×4线性解+6×6特征解；s=3用15×15线性解+20×20特征解。本轮两者均已运行；s=4的70维仅为理论维数，**NOT_RUN**。D仍按中央二项式增长，不能包装成多项式复杂度的大规模全域算法。

新目标的最高次块T对共轭观测仿射，因为每个最大子式只有首列的最高次系数依赖该观测；R及乘法矩阵仍变化。当前实现按目标生成系数并完整计价，未把一次性模板优化当成已实现加速。

**标准代数工具的物理专门化；原创性待决。** 严格假设、非必要性反例与证明见[数学理论第10节](MATHEMATICAL_THEORY.md)。
