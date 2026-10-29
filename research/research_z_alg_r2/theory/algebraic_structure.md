# 固定 Green 材料切片的计算代数结构

状态：数学推导；未在本文声称新的物理对称、全域覆盖或计算优势。只读取历史 E1 的公开物理定义；没有读取新的 private GT。所有结论保留固定 Green、每 cell 一个共享实标量材料以及可逆 resolvent 的假设。

## 1. 从真实状态方程到可计算的小块

历史 `research/research_z_alg_asymmetric_validation/src/physics.py` 使用

\[
J=D_m(E^{\rm inc}+GJ),\qquad
Y=SJ,\qquad A(m)=I-D_mG,
\tag{1}
\]

其中 \(m\in\mathbb R^n\)、\(D_m=\operatorname{diag}m\)、\(J,E^{\rm inc}\in\mathbb C^{n\times p}\)、\(S\in\mathbb C^{r\times n}\)。同一个 \(m_i\) 用于全部照明；不允许逐源材料。固定带标签的 \(G,S,E^{\rm inc}\) 不随候选材料改变。\(A(m)\) 可逆时

\[
J=A(m)^{-1}D_mE^{\rm inc},\qquad
E=E^{\rm inc}+GJ=(I-GD_m)^{-1}E^{\rm inc}.
\tag{2}
\]

### 1.1 任意基点上的精确小块增量

选取要变化的 \(k\) 个 cell，以列选择矩阵 \(P\in\mathbb R^{n\times k}\) 表示。令 \(m=m^0+Pz\)，基点 \(A_0=I-D_{m^0}G\) 可逆。先计算 \(J_0,E_0,Y_0\)，并定义

\[
U=A_0^{-1}P,\quad
H=P^TGU,\quad B=P^TE_0,\quad R=SU.
\tag{3}
\]

直接相减状态方程得到

\[
A_0(J-J_0)=PD_zP^T(E^{\rm inc}+GJ).
\]

置 \(Q=D_z[B+P^TG(J-J_0)]\)，则

\[
J-J_0=UQ,\qquad Q=D_z(B+HQ),
\]

所以

\[
\boxed{Y-Y_0=R(I-D_zH)^{-1}D_zB.}
\tag{4}
\]

这是精确块消元，不是 Born 截断。行列式引理同时给出

\[
\det A(m)=\det A_0\det(I-D_zH).
\tag{5}
\]

故基点有效时，原 resolvent 的可逆性与这个 \(k\times k\) resolvent 完全等价。式 (4) 的低秩求逆本身是标准公式；本轮可执行的新路线在后续的低度有限根方程，而不在重命名该恒等式。

### 1.2 固定其余材料的 Schur 形式

也可以不引入增量基点。将 cell 分成固定的 \(A\) 与未知的 \(C\)，若 \(I-D_{m_A}G_{AA}\) 可逆，定义

\[
M_A=(I-D_{m_A}G_{AA})^{-1}D_{m_A}.
\]

由 \(J_A=M_A(E^{\rm inc}_A+G_{AC}J_C)\) 得

\[
\begin{aligned}
H&=G_{CC}+G_{CA}M_AG_{AC},\\
B&=E^{\rm inc}_C+G_{CA}M_AE^{\rm inc}_A,\\
R&=S_C+S_AM_AG_{AC},\\
Y_0&=S_AM_AE^{\rm inc}_A.
\end{aligned}
\tag{6}
\]

此时式 (4) 中 \(z=m_C\)。Schur 子块可能奇异而完整 \(A(m)\) 可逆；这种失败属于消元坐标失效，不能判定物理材料非法。需换固定块，或用式 (3) 的可逆完整基点。

## 2. 分母与分子的多重仿射系数

以下写 \(\Delta(z)=\det(I-D_zH)\)，并以 \(z_I=\prod_{i\in I}z_i\)、\(z_\varnothing=1\) 表示单项式。行列式对每一行线性，且第 \(i\) 行仅含一个 \(z_i\)，所以

\[
\boxed{\Delta(z)=\sum_{I\subseteq[k]}(-1)^{|I|}\det H_{II}\,z_I.}
\tag{7}
\]

每个变量的次数不超过 1；总次数最多 \(k\)。空集的行列式约定为 1。

对接收行 \(R_\rho\) 和照明列 \(B_\sigma\)，定义

\[
N_{\rho\sigma}(z)=R_\rho\operatorname{adj}(I-D_zH)D_zB_\sigma.
\tag{8}
\]

为明确分子也多重仿射，考察

\[
\det(I-D_zH+tD_zB_\sigma R_\rho)
=\Delta(z)+tN_{\rho\sigma}(z).
\]

对每个主子式取 \(t\) 的一次系数，即得

\[
\boxed{N_{\rho\sigma}(z)=
\sum_{\varnothing\ne I\subseteq[k]}
(-1)^{|I|-1}R_{\rho,I}\operatorname{adj}(H_{II})B_{I,\sigma}\,z_I.}
\tag{9}
\]

该式在 \(H_{II}\) 奇异时仍成立，因为使用 adjugate，未除以其行列式。对单元素集，\(\operatorname{adj}([h])=[1]\)。测量目标 \(T=Y_*-Y_0\) 给出

\[
f_{\rho\sigma}(z)=T_{\rho\sigma}\Delta(z)-N_{\rho\sigma}(z)=0.
\tag{10}
\]

实材料时把每式分成实部和虚部，仍为实系数多重仿射多项式。清分母只在 \(\Delta\ne0\) 上等价；极点基点必须饱和并在候选阶段重新排除。

### 2.1 四变量材料切片的 24 路径界

E1 有六个实材料、两个复照明输出。固定两个事前选择的 cell 材料后，剩余四个变量对应四条实多重仿射方程。齐次化到 \((\mathbb P^1)^4\)，每个方程的多次数至多 \((1,1,1,1)\)。其交数为

\[
[h_1h_2h_3h_4](h_1+h_2+h_3+h_4)^4=4!=24,
\quad h_i^2=0.
\tag{11}
\]

因此，零维复交集最多有 24 个根，按重数计并包括无穷远；存在高维分量时，同一 Bézout 界用于孤立复根，但不能把高维分量称为已枚举。稠密总次数界 \(4^4=256\) 在这里太松。

可取起始方程

\[
g_i(z)=\prod_{j=1}^4(a_{ij}z_j+b_{ij}),\qquad i=1,\ldots,4.
\]

泛型系数下，一个根给每行恰选一个零因子，且每列恰使用一次；24 个排列给出 24 个显式起点。泛型 gamma 同伦保持多次数。最终极点、无穷远、重复根与未跟踪路径均须独立记账。这个 24 界不是 24 个合法实材料根的承诺。

## 3. 两材料块：精确消元一般是二次，不是四次

令 \(z=(x,y)\)，\(H=(h_{ij})\)。直接求逆给

\[
\Delta=1-h_{11}x-h_{22}y+(h_{11}h_{22}-h_{12}h_{21})xy,
\tag{12}
\]

以及每个通道的分子

\[
N_{\rho\sigma}=R_{\rho1}B_{1\sigma}x+R_{\rho2}B_{2\sigma}y
+[R_{\rho1}(h_{12}B_{2\sigma}-h_{22}B_{1\sigma})
+R_{\rho2}(h_{21}B_{1\sigma}-h_{11}B_{2\sigma})]xy.
\tag{13}
\]

任取两个独立实约束，写成

\[
f_i=(a_i+b_ix)+(c_i+d_ix)y=0,\qquad i=1,2.
\]

消去 \(y\) 的 resultant 为

\[
\boxed{r(x)=(a_1+b_1x)(c_2+d_2x)-(a_2+b_2x)(c_1+d_1x),}
\tag{14}
\]

次数最多 2。若某个 \(c_i+d_ix\ne0\)，则
\(y=-(a_i+b_ix)/(c_i+d_ix)\)，再验全部通道及 \(\Delta\ne0\)。若两个 \(y\) 系数同时为零，不能直接丢掉该 \(x\)：须检查两常数项，决定无解、整条 \(y\) 线或额外方程决定的解。若 resultant 恒为零，则该方程对不够独立或共享分量，应更换方程或作分量处理。

对两个泛型双线性方程，\(\mathbb P^1\times\mathbb P^1\) 的交数是 2。总次数 Bézout 的 4 只是松界；把精确复振幅观测消元强写成四次没有依据。改用强度、平方残差或其它非等价统计量会改变次数与问题，不能借此声称本模型有四次内核。

## 4. 接收传递满列秩时直接恢复

式 (4) 中若 \(R\in\mathbb C^{r\times k}\) 满列秩，测量直接确定

\[
Q=R^\dagger T,\qquad (I-RR^\dagger)T=0.
\]

令 \(V=B+HQ\)。每个 cell 的全部源须满足 \(Q_{i,:}=z_iV_{i,:}\)。若 \(\|V_{i,:}\|>0\)，唯一实候选是

\[
z_i=\frac{\operatorname{Re}\langle V_{i,:},Q_{i,:}\rangle}{\|V_{i,:}\|^2},
\tag{15}
\]

并严格检查复向量残差 \(Q_{i,:}-z_iV_{i,:}\)。这里内积定义为 \(\langle V,Q\rangle=\sum_s\overline V_sQ_s\)。若 \(V_{i,:}=0\)，必须有 \(Q_{i,:}=0\)，不能除零。该满秩恢复是标准线性代数，不是新颖性主张。接收矩阵只有满行秩时此论证不成立。

## 5. 一接收器的对偶内核：六路径而非二十四路径

这一步利用多照明的全部共享材料约束。它不需要 \(G=G^T\)、几何对称或物理群。

当 \(r=1\) 时，定义一个接收传递行

\[
w=R(I-D_zH)^{-1}D_z\in\mathbb C^{1\times k}.
\tag{16}
\]

由 \((I-D_zH)^{-1}D_z=D_z(I-HD_z)^{-1}\) 得

\[
\boxed{T=wB,\qquad w=(R+wH)D_z.}
\tag{17}
\]

令 \(v=R+wH\)，则每个材料满足 \(w_i=z_iv_i\)。反向也成立：构成方程可重排为 \(w(I-HD_z)=RD_z\)。在 \(I-HD_z\) 可逆时，该方程唯一给出式 (16)，继而恢复全部照明输出。

### 5.1 先消掉全部线性测量

设 \(\operatorname{rank}_{\mathbb C}B=s\)。若 \(k=2s\)，对一致目标 \(wB=T\)，取特解 \(w_0\) 和行零空间基 \(N\in\mathbb C^{s\times 2s}\)，满足

\[
w_0B=T,\qquad NB=0,\qquad \operatorname{rank}N=s.
\]

SVD 可用于稳定求取这两个量，但应对 \(B^T\) 取零空间；这里的转置来自代数方程，不能擅自换成共轭转置。一次固定投影之后

\[
w=w_0+cN,\quad c\in\mathbb C^{1\times s},\qquad
v=v_0+cM,\quad v_0=R+w_0H,\ M=NH.
\tag{18}
\]

共享实材料的条件变成 \(\operatorname{Im}(w_i\overline v_i)=0\)。在复代数化中把 \(\overline c\) 政名为独立变量 \(d\)，而只共轭固定系数：

\[
\begin{aligned}
\widetilde w&=\overline{w_0}+d\overline N,\\
\widetilde v&=\overline{v_0}+d\overline M,\\
\boxed{f_i(c,d)&=w_i\widetilde v_i-\widetilde w_i v_i=0,\qquad i=1,\ldots,2s.}
\end{aligned}
\tag{19}
\]

每式在 \(c\) 组和 \(d\) 组各为一次；实际实材料根还必须满足 \(d=\overline c\)。当 \(v_i\ne0\) 时，式 (19) 的实固定点等价于

\[
\boxed{z_i=\frac{\operatorname{Re}(w_i\overline v_i)}{|v_i|^2}\in\mathbb R,
\qquad w_i-z_iv_i=0.}
\tag{20}
\]

复数商 \(w_i/v_i\) 可用于代数恢复；实际验收用式 (20) 加复残差，避免把小虚部静默截断。

### 5.2 通用 \(\binom{2s}{s}\) 界与六个显式起点

将式 (19) 齐次化到 \(\mathbb P^s_c\times\mathbb P^s_d\)。每条方程属于 \((1,1)\) 类。其零维交数为

\[
\boxed{[h_c^sh_d^s](h_c+h_d)^{2s}=\binom{2s}{s}.}
\tag{21}
\]

因此，在非节点、可逆 resolvent 域内，只要相应根是孤立复根，其个数按重数最多 \(\binom{2s}{s}\)。零维交集可以直接使用此界；实 Jacobian 满秩的物理解也必为孤立复根。仅仅“在实域孤立”不足以代替这一假设，例如 \(x^2+y^2=0\) 的实原点孤立但复零集是曲线。

对 E1 的四材料切片，\(k=4,s=2\)，所以界为 **6**。材料多重仿射的 24 界仍正确，但没有使用一接收器、两个独立照明和实材料的这组更强结构。

显式起始系统可取

\[
g_i(c,d)=(a_{i0}+a_i^Tc)(b_{i0}+b_i^Td),\qquad i=1,\ldots,2s.
\tag{22}
\]

泛型系数下选择一个大小为 \(s\) 的指标集 \(I\)：令 \(i\in I\) 的 \(c\) 因子为零，令 \(i\notin I\) 的 \(d\) 因子为零。两组各解一个 \(s\times s\) 线性系统。全部 \(\binom{2s}{s}\) 个选择就是全部起点。保持双次数的 gamma 同伦为

\[
H_i(c,d,t)=(1-t)\gamma g_i(c,d)+t f_i(c,d).
\tag{23}
\]

路径中一般 \(d\ne\overline c\)；只在目标端筛选物理实固定点。泛型复 gamma 与起始系数用于避开内部非横截问题，不是对目标端奇异性或追踪数值误差的证书。

### 5.3 实现所需的九个系数

对 \(s=2\)，按常数、\(c_1,c_2\)、\(d_1,d_2\)、\(c_ad_b\) 编码。每个 cell \(i\) 的系数为

\[
\begin{aligned}
f_{i,00}&=w_{0i}\overline{v_{0i}}-\overline{w_{0i}}v_{0i},\\
f_{i,c_a}&=N_{ai}\overline{v_{0i}}-\overline{w_{0i}}M_{ai},\\
f_{i,d_b}&=w_{0i}\overline{M_{bi}}-\overline{N_{bi}}v_{0i},\\
f_{i,c_ad_b}&=N_{ai}\overline{M_{bi}}-\overline{N_{bi}}M_{ai}.
\end{aligned}
\tag{24}
\]

这些系数完全由冻结算子、切片和观测构成，不使用未知真材料。求解器应保存全部起点、路径状态、候选端点、共轭误差、节点阈值、复构成残差、材料合法性、resolvent 最小奇异值和独立 forward 残差。

## 6. 两个关键退化结论

### 6.1 对偶节点意味着真实连续材料自由度

假设一个合法可逆物理解有 \(v_i=0\)。式 (17) 给 \(w_i=0\)。仅改变 \(z_i\) 时，同一个 \(w\) 仍满足所有 \(w_j=z_jv_j\)，因为第 \(i\) 个等式仍为 \(0=z_i0\)，其它等式未改变。可逆性是开条件，因此附近一段实 \(z_i\) 都有相同输出 \(wB\)。

所以，只要该 cell 的合法区间有非零宽度，这个点不是孤立物理解；即便在普通盒边界，至少保留一个方向的局部材料区间。排除 \(v_i=0\) 不会丢失可逆、正则的孤立实根，但会丢失真正的连续纤维，必须单列，不能叫“无解”。

### 6.2 对偶 Jacobian 与材料 Jacobian 正则性等价

在一个物理解，定义 \(\phi_i(c)=\operatorname{Im}(w_i\overline v_i)\)。因为 \(w_i=z_iv_i\) 且 \(z_i\in\mathbb R\)，对测量零空间扰动 \(\delta w=\delta cN\) 有

\[
\delta\phi_i=
\operatorname{Im}[(\delta w_i-z_i\delta v_i)\overline v_i],
\qquad \delta v=\delta wH.
\tag{25}
\]

当所有 \(v_i\ne0\) 时，\(\delta\phi=0\) 当且仅当存在实 \(\delta z\) 满足

\[
\delta w(I-HD_z)=(\delta z)\odot v.
\tag{26}
\]

因此

\[
\delta w=[(\delta z)\odot v](I-HD_z)^{-1},\qquad
\delta Y=\delta wB.
\tag{27}
\]

而 \(\delta wB=0\) 正是已消去的测量约束。两边的核通过式 (26) 双射，故在 \(k=2s\) 下，\(2s\times2s\) 的实对偶 Jacobian 可逆，当且仅当同一材料切片的实响应 Jacobian 可逆。这是局部正则性的等价，不保证条件数相同；小 \(|v_i|\) 会恶化反向恢复。

## 7. 一般接收数及不能推广的部分

原问题也可写传递矩阵 \(W=S(I-D_mG)^{-1}D_m\)，满足

\[
Y=WE^{\rm inc},\qquad W=(S+WG)D_m.
\tag{28}
\]

每列的同一个实 \(m_i\) 需要同时适用于全部接收行。这是源/接收方向的代数转置处理，不要求传播 reciprocity。若源数多于接收数，\(W\) 的测量零空间坐标可能少于逐源 current 坐标；但上述 \(\binom{2s}{s}\) 界只对“一接收行、\(k=2s\)、\(\operatorname{rank}B=s\)”直接成立。不能把多接收的列共享约束忽略后套用六根结论。

如果 \(k>2s\)，对偶坐标有 \(2(k-s)\) 个实自由度、只有 \(k\) 个实材料约束，预期材料纤维维数为 \(k-2s\)。如果 \(B\) 降秩，必须按实际秩重新建模；固定声称六路径会漏掉高维结构。数值 SVD 的截断若丢弃非零奇异值，已经改变观测等式，不能叫精确换基。

## 8. 可直接执行的求解路线及失败口径

1. 冻结 \(G,S,E^{\rm inc}\)、任意事前 cell 切片、合法盒和目标；从式 (3) 或 (6) 建 \(H,B,R,Y_0\)。
2. 先验证所用消元块有效，以及实际 \(\operatorname{rank}B=2\)。一接收四变量时按式 (18)、(24) 建双线性系数。
3. 解析生成六个起点，逐条跟踪式 (23)，保留成功、失败、无穷远、合并和奇异端点。
4. 对全部目标根检查 \(d=\overline c\)；用式 (20) 恢复实材料，验全部复构成等式、合法盒与非奇异 resolvent。
5. 用原完整 forward 独立验收；不能用已强制满足的 \(wB=T\) 给自己作唯一物理证据。
6. 数值完成六路径可记 `NUMERICAL_ALL_PATHS_COMPLETED`。`COMPLETE_CERTIFIED` 还需要有效根/路径验证及目标系统维数、奇异性、饱和域与无穷远的完整处理。任一未解路径都保留 `UNKNOWN`，不得用剩余五根的干净残差补足。

本文件给出有限根内核及其上界；它不证明 E1 的实际代数次数恰为 6，也不证明六个端点全为实数、合法材料或互异，更不证明成本优于普通优化器。
