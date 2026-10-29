# R3 全材料替代结构：观测先消元与共同实材料约束

状态：`THEORY_DERIVED / AUTHORIZED_FOR_ONE_DEV_PIVOT / IMPLEMENTATION_IN_PROGRESS_AT_HANDOFF`。主线程已明确授权一次全部接收global-dual/variable-projection DEV比较，具体实现交给R3 solver。此授权不代表该臂已实现、通过或进入holdout。本理论工作新增 Maxwell RHS=0，数值实验=0，未读取私有 GT，未改变冻结R2。最终取舍由主线程完成。

## 1. 建议与实际合同

初始 DEV 合同由主线程明确为 **N=12 个独立实材料全部未知，3 个源、3 个复接收，共18个实观测，材料 box [.5,30]**。不是 N12/6源/1接收的平衡问题。N16 的采集未在本文件中补设；下文一般式由实际源数与接收数决定维数。

建议保留原物理目标和所有观测，以 **全部接收的观测先消元 + 共同实材料 variable projection** 产生全材料候选，随后交给原目标的 F12 坐标最小化或 Full-J 校正与独立 Maxwell 验收。备选的价值是每步允许所有材料共同变化、没有冻结块外材料、没有丢弃观测；它不是全局求解保证，也没有已证明的速度优势。

这与局部 Q6/Q20 的区别是表示与耦合范围：Q6/Q20 对一个固定补集的4/6材料切片可给全部小代数候选；这里保留全部未知材料和接收，仅将满足观测的辅助变量线性消元。它不免除非凸搜索。

## 2. 固定有理 family 不要求冻结真实材料

统一采用 R2 F06 的记号：

\[
\Delta=\operatorname{diag}(\delta),\qquad
Y(\delta)=Y_0+L(I-\Delta K)^{-1}\Delta q
=Y_0+L\Delta(I-K\Delta)^{-1}q.
\]

其中 K∈C^(N×N)、q∈C^(N×s)、L∈C^(r×N)，全部由已知几何和一个固定参考材料决定。参考 m0 可为一个公开辅助状态，δ=m−m0；N个坐标全部允许变化。若直接取零参考，M0=I，K=G、q=Einc、L=S、Y0=0，即可覆盖整个允许材料 box；零参考不必是允许的最终材料，也不是已知真实背景。保留非零参考时需要支付其固定 factorization 与列求解。

假设 rank(q)=s≤N，并在物理候选处要求

\[
R_\delta=I-K\Delta\quad\text{可逆}.
\]

不要求 Born/Neumann 级数收敛，不用材料倒数，不删除零增量。由于这里所有 N 个材料均未知，F06 是固定系数复用，**不是 N→k 的状态降维**。

## 3. 将全部观测精确放入一个仿射空间

定义真实接收传递矩阵

\[
W_\delta=L\Delta R_\delta^{-1}\in\mathbb C^{r\times N},\qquad T=Y_{obs}-Y_0.
\]

真实数据一致性和共同材料方程为

\[
Wq=T,\qquad W=(L+WK)\Delta.\tag{A1}
\]

第二式由 W Rδ=LΔ 得到，所有接收行必须使用同一个 Δ。这是 R2 F09 的多接收扩展，不是逐接收生成不同材料。

取一次固定 QR/SVD 得到 Z∈C^((N−s)×N)，使 Zq=0、ZZ*=I；取 W0=Tq†。因为 q 满列秩，W0q=T，所有精确满足观测的 W 都唯一表示为

\[
W(C)=W_0+CZ,\qquad C\in\mathbb C^{r\times(N-s)}.\tag{A2}
\]

Z、ZK 可按采集缓存；新数据只更新 W0 与 V0=L+W0K。注意这里左零空间是 q 的**转置**约束 Zq=0，若使用复 SVD/库函数必须检查其实际返回的转置与共轭约定，不能将 Zq=0误写成 Z*q=0。

当前 N12/s3/r3：C为3×9复矩阵，共54实坐标。它不是一个12维坐标压缩。若仅保留一个接收，则C有18实坐标，而12个材料只接受6实观测，在正常秩条件下可有6维局部材料纤维；它不能被称为全18观测的成像结果。

## 4. 共同实材料的精确可行性

令 V(C)=L+W(C)K。对于每个材料格 i，记 u_i=W[:,i]、v_i=V[:,i]∈C^r。要求存在一个共同实数 δ_i∈[a_i,b_i]，使

\[
u_i=\delta_i v_i.\tag{A3}
\]

若 h_i=||v_i||²>0，则唯一最小二乘材料为

\[
d_i(C)=\frac{\operatorname{Re}(v_i^*u_i)}{h_i},\qquad
\delta_i(C)=\operatorname{clip}(d_i(C),a_i,b_i).\tag{A4}
\]

没有 box 时，严格实材料相容等价于

\[
||u_i||^2||v_i||^2-[\operatorname{Re}(v_i^*u_i)]^2=0.
\]

这是 C 上的四次式。若希望保留二次约束，可改用所有接收分量的复交叉式

\[
u_{ia}v_{ib}-u_{ib}v_{ia}=0\;(a<b),\quad
\operatorname{Im}(v_i^*u_i)=0,
\]

加上二次 box 条件

\[
\operatorname{Re}(v_i^*u_i)-a_i||v_i||^2\ge0,\quad
b_i||v_i||^2-\operatorname{Re}(v_i^*u_i)\ge0.
\]

在 v_i≠0 时这些条件充分必要。**v_i=0 必须另行处理：u_i=0时该状态不确定此材料；u_i≠0时不相容。** 不能用分母 epsilon-floor 把节点伪装成唯一可恢复材料。A17X `row_fit` 的显式 undetermined 处理比简单截断分母更适合本处。

一接收时交叉式消失，只有 N 个实二次方程 Im(w_i conj(v_i))=0；只有再满足 N=2s，才成为 N个实变量/N个方程的平衡情形。这是最初理论备选的特殊情况，不用于当前 DEV 的维数声明。

**精确对应命题。** 假设 q满列秩、Rδ可逆。物理材料 δ使Y(δ)=Yobs，当且仅当存在C使(A3)在全部格点成立。证明：物理解提供 Wδ及(A2)；反之(A3)给 W Rδ=LΔ，所以 W=Wδ，再由Wq=T得到原完整观测。节点只能使某些材料恢复不唯一，不能免除(A3)。这是有限离散模型的精确对应，不是连续 Maxwell 的唯一性定理。

## 5. 一个可以直接执行的全材料候选生成器

定义非负兼容性目标

\[
\mathcal E(C,\delta)=||W(C)-(L+W(C)K)\Delta||_F^2.
\tag{A5}
\]

它处处满足全部观测约束 Wq=T，但只有 E=0且物理可逆时，W才是实际材料的传递矩阵。中间的 W不是可接受物理解。

固定δ时，记

\[
A(\delta)=Z(I-K\Delta)\in\mathbb C^{(N-s)\times N},
\quad B_0(\delta)=W_0-(L+W_0K)\Delta.
\]

求

\[
C_*(\delta)=\arg\min_C||B_0+CA||_F^2.\tag{A6}
\]

计算上对 A^T（不是未经说明的共轭置换）做一次复 QR/SVD，将r个接收作为共享分解的多个右端，求 A^T C^T≈−B0^T。若Rδ可逆，A满行秩，解唯一。使用QR/SVD，不形成AA*正规方程。当前一次分解尺寸为12×9，3个右端；每次材料变化都要重新分解，不能只报9阶三角回代。

固定C时，全部δ_i可按(A4)同时更新。v_i=0时保留先前合法δ_i并标记未定；若u_i非零，本轮残差不能消失。两个步骤都精确最小化同一个目标，因此交替目标值单调不增；只有目标值收敛，没有已证明的全局最优、根完整性或材料序列收敛。

也可将C完全消去，以

\[
\psi(\delta)=\mathcal E(C_*(\delta),\delta)
\]

作为一个N变量、box约束的 variable-projection 目标。记E=W−VΔ在C*处的残差；最小二乘正交性消去dC*项，得到精确梯度

\[
\frac{\partial\psi}{\partial\delta_i}
=-2\operatorname{Re}\sum_{a=1}^r\overline{E_{ai}}V_{ai}.\tag{A7}
\]

可选普通 box L-BFGS 或交替最小化作为实现，二者都属于标准非凸优化。任何微分资格都只在保持内层满秩的域中成立。

建议接口只需要公开内容：

```text
candidate_generator(K, q, L, Y0, Yobs, m0, box,
                    public_material_starts, frozen_caps)
  fixed setup: rank(q), q-pseudoinverse, Z, ZK
  per observation: T, W0, V0
  per public material start:
    minimize compatibility by shared QR + all-material updates
    retain stop reason, QR/rank diagnostics and final material
    refine against original all-channel objective if budget permits
    run original independent Maxwell validator
  return every retained candidate and every failure, with costs
```

使用与强基线相同的公开材料起点，由(A6)初始化C，不从GT或已评分目标选初值。不在此文件中发明次数、阈值、预算或评分门；采用主线程冻结合同。构造重复点、失败点、触及cap和rank失败都要保留。

## 6. 消元后到底优化什么：精确残差度量与反例

这是本备选最重要的限制。记

\[
Q_\delta=R_\delta^{-1}q,\qquad
R_d=T-W_\delta q=Y_{obs}-Y(\delta).
\]

对满足测量的任何W，兼容性残差E满足

\[
E=(W-W_\delta)R_\delta,
\qquad E Q_\delta=R_d.\tag{A8}
\]

反之满足EQδ=Rd的E均可反推一个W=Wδ+ERδ^-1，使Wq=T，因此(A6)等价于在EQδ=Rd下求最小Frobenius范数E。Qδ满列秩，最小范数解与消元目标为

\[
E_*=R_d(Q_\delta^*Q_\delta)^{-1}Q_\delta^*,
\boxed{\psi(\delta)=\operatorname{tr}\{R_d(Q_\delta^*Q_\delta)^{-1}R_d^*\}.}\tag{A9}
\]

所以ψ是**随材料变化的观测残差度量**，不是原Frobenius噪声似然。它与原目标具有相同的精确零点（在可逆域），但一般具有不同的非零极值与优化轨迹。直接从(A8)有

\[
||R_d||_F\le ||E||_F\,||Q_\delta||_2.
\]

没有稳定性上界时，小兼容性残差不能保证小数据残差。用Qδ构造严格等价的物理权重又需要求变化的resolvent，不能将其视为免费修复。

**解析反例，仅用于否定通用停止规则，不冒称实际 Hankel 几何。** 取N2/s1/r1，K=diag(1,0)、q=(1,1)^T、L=(0,1)、Δ=diag(1−τ,1)，0<τ<.5，T=2。则Rδ=diag(τ,1)可逆，物理预测恒为1，真实误差恒为1；Qδ=(τ^-1,1)^T，而

\[
\psi=\frac{\tau^2}{1+\tau^2}\to0.
\]

材料仍在[.5,30]内，却可沿近极点状态把兼容性目标任意压低。因而该生成器必须保留原全观测残差、原 resolvent 条件/极点门和独立求解验收。不能仅靠ψ或ALS停滞宣告找到成像结果。

含噪时真实观测通常不严格位于材料响应流形；把Wq=T写成硬约束只定义了一个候选生成方式，并未建立正确噪声 likelihood。最终资格始终相对于原观测与原噪声合同。

## 7. F12 是必须面对的廉价强对照

固定其它所有材料，改变一个实坐标t。对当前精确状态，由F06得到

\[
Y(t)=Y_c+\frac{a t}{1-\kappa t},\quad
a=(SM^{-1}e_i)E_{i,:},\quad
\kappa=G_{i,:}M^{-1}e_i,
\]

其中多接收、多源的a按Frobenius范数合并。令r=Yc−Yobs、b=a−κr，则纯数据目标为

\[
\phi(t)=\frac{A t^2+B t+C}{D t^2+E t+1},
\]

其驻点只需解

\[
(AE-BD)t^2+2(A-CD)t+(B-CE)=0.
\]

其中A=||b||²、B=2Re〈r,b〉、C=||r||²、D=|κ|²、E=−2Reκ。比较所有合法实根与区间端点，并处理极点，即可精确最小化该坐标上的有理数据目标。固定观测白化可吸收进r与a；不能擅自把新增非二次regularizer也称为相同二次根公式。若整个坐标区间保持可逆，候选清单完整；若额外规定σmin等稳定性硬门，需要处理这些门形成的可行边界，仅拒绝不合格驻点并不保证得到受该门约束的最小值。

维护当前H=M^-1时，接受t后的Sherman–Morrison更新为

\[
H^+=H+\frac{t}{1-\kappa t}(He_i)(G_{i,:}H).
\]

它可减少重复大LU，但每次更新有O(N²)工作和数值漂移；必须记录更新与必要重新分解，终点独立求解。循环F12只保证每步原目标不增，仍可能困于坐标极小。它是廉价强基线/校正器，不能包装成新的全局根枚举。

## 8. 为什么局部代数块不等于全材料结构

F11的Schur消元在冻结补集B、相应TB可逆时精确；此处“精确”指该**条件切片上的响应**。当当前块外材料错误时，4/6个块变量一般不足以满足全部18实观测。若Q6/Q20只用选定接收/源的平衡子系统，就必须承认其余观测尚未满足，完整校正及其费用不可省略。

全材料备选的优化变量和(A5)同时涉及全部材料及全部接收，它可能避开一个错误补集造成的条件障碍；这只是结构上的可能收益。ALS/varpro自身有非凸驻点、近极点诱导错误度量及内层QR成本，未证明优于Q6/Q20或普通Full-J。

N12/6源/1接收的F15商维数是924，N16/8源/1接收是12870，说明直接平衡推广会很快扩大；这些数字只是**另外的假设采集**，不是当前3源3接收模型的根数上界。当前合同不能直接套F15中央二项式证明，也不能通过将多接收数据随意打包成“更多独立源”获得同一结构。

## 9. 复用、拒绝与先例边界

直接复用：F06全材料有理族、F09测量先消元、R9 F-OP-008/010固定算子与正确多源/接收求导、F12原目标单坐标精确最小化。F11/F14/F15保留为条件块proposal与成本参照，不把其覆盖范围扩大。

明确拒绝：统一固定对角化处理全部非均匀材料；删除receiver-dark current即声称完整物理降维；使用单接收根作为全部观测解；将有限兼容性残差或profile当作物理likelihood/posterior；为了展示该公式另跑小型恒等式实验；未经证据声称新理论或更快。

先例定位为 **CSI式本构残差优化与标准variable projection的接收对偶、数据先消元版本**。CSI的明确历史来源为 [van den Berg 与 Kleinman 1997, A contrast source inversion method](https://doi.org/10.1088/0266-5611/13/6/013)；可分离非线性最小二乘的变量消元来源为 [Golub 与 Pereyra 1973](https://doi.org/10.1137/0710036)。本次仅核对出版社元数据/可见摘要，没有逐式证明与历史CSI完全相同或不同。因此没有独立新颖性结论；“观测先消元的特定组合是否已出现”仍须主线程的先例审查。

本文件已经给出可实施接口、精确对应证明、真实目标关系与致命停止规则反例，达到理论交付停止点。主线程已授权一次DEV pivot并要求向solver传递全矩阵公式与梯度；该接口已交付。运行结论、继续/停止与holdout冻结由主线程依据新收据裁决。
