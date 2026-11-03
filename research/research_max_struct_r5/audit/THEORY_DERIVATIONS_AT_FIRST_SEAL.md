# R5 实际使用的结构推导

## 1. 本构切向与材料相关传播：根处公式的遗漏项

以下是有限二维 scalar-TM VIE；同形线性代数可用于合法矢量离散，但三维本构参数和离散自项必须另行说明。固定辅助电流 J，Q=Einc+GJ，r=J−diag(m)Q。材料为实数，源 metric Ω_s 固定。把复矩阵乘 Ω_s 的平方根后实化，以欧氏内积表示全部公式。

记 C_mat h=diag(h)Q，P(m) 为可逆空间残差传播，R=P(m)r。则

```text
dR[h] = −P C_mat h + (dP[h])r.
g_state = −C_matᵀ Pᵀ P r + a_r,
<a_r,h> = <P r,(dP[h])r>.
```

这只是乘法微分，但给出一个不能省略的实验区分。r=0 时 a_r=0；非根处 a_r 对 r 二次，对真实优化轨迹未必小。冻结 metric 的读出公式不能自动成为非零缺陷、带 prior/box 的完整材料更新公式。

令 T_r h=(dP[h])r，则局部真实线性化是 R−(P C_mat−T_r)h；冻结传播线性化为 R−P C_mat h。二者连 Gauss–Newton 切向都不同。a_r 的上界依赖 ||r||²、P 和 dP 的范数；没有可用的尺度比较时，不用“残差小”代替测量。根处一致并不等于离根路径接近。

## 2. one-hop 的正反相分解与可判别的反例

b_i=β−G_ii，d_i=1+b_i m_i，a_i=m_i/d_i，u_i=β r_i/d_i，v=diag(a)Goff u；相位控制 s=±1，R_s=u+s v。

```text
du_i[h] = −β (Q_i+b_i J_i) h_i / d_i²,
dv[h] = diag(h/d²)Goff u + diag(a)Goff du[h].
f_s = 1/2||u||² + s Re<u,v> + 1/2||v||²,
g_s = g_self + s g_cross + g_quad.
g_cross·h = Re<du[h],v> + Re<u,dv[h]>.
```

内积均含同一个固定 Ω_s；prior 不属于上面三项，必须单列。正反相控制保持 self 与 quad，改变 cross；但 cross 的**值**与 cross 的**导数方向**不是同一量。cross<0 不能推出 −g_cross 是原物理 Phi 的下降方向，更不能推出材料误差下降。其直接可反驳预测是：即使两个检查点都有负 cross，相对物理梯度的方向关系仍可相反。

更一般地，只依赖 (G,J,Q,m,Ω_s) 的局部 residual 指标没有接收数据方向。在固定辅助状态下改变合法目标产生的 Y，原物理梯度会变，而 cross 完全不变。在内部合法 m 邻域取目标 m±t h，数据为对应完整正演，则小 t 时数据梯度的主项符号相反。任何不使用接收数据的 cross 符号规则都不能普遍判定同一非零材料步的物理收益。这是局部反例结构，不证明任一统计预测量无用，也不需重复机器精度恒等式实验。

## 3. 冻结 metric 的材料读出：保留方向而非只看能量

C_mat 满列秩，M=PᵀP 正定，Π_mat=C_mat C_mat†。固定线性问题的解满足

```text
D_M−C_mat†=(C_matᵀ M C_mat)⁻¹ C_matᵀ M(I−Π_mat).
```

此 R4 已证公式直接复用。对当前 b=r，输运后的不可解释分量驱动
δ_transport=(C_matᵀ M C_mat)⁻¹ C_matᵀ M(I−Π_mat)r。
可观测量需要 δ_transport 与实际 dm、负原物理梯度的夹角及投影，而非单独 ||δ_transport||。只有逐cell共享源权重时，局部标量权重在无正则自由读出中抵消；空间耦合通常不抵消。相同正定 metric 可在不同 b 上有益或有害。

当前四cell块允许每块4×4实 Gram 解，不必构造全域(2NsN)²投影矩阵。若场切向弱，报告其 Gram 条件和最弱值；不要靠任意截断把弱方向改成证明的精确核。真实材料内层含 prior、box 和非线性 P，另用原程序作相同状态反事实，不能把此局部LS当已执行算法。

## 4. 原物理方向与曲率：继承完整双边传播

物理状态 A=I−XG，j*=A⁻¹XEinc，E*=Einc+Gj*，O=S A⁻¹。材料导数为 dF[h]=O diag(h)E*。Gram 是 Re[(O*O) ∘ (conj(E*)E*ᵀ)]/||Y||²；这是 Atlas/P3L 已有双边结构，不是本轮新信息。

沿固定 h，δF(t)=t O(I−t diag(h)G A⁻¹)⁻¹diag(h)E*。原目标曲率含 ||dF[h]||² 与 residual-weighted second response，不能把 GN 曲率称精确 Hessian。使用已有 tangent/HVP 和同一保留LU测必要方向，不为完整谱构造新的近似物理。

这些式子能预测局部 Phi 的线性/二次变化；要检验有限材料步，仍算真实完整 forward。低 Phi 不是材料恢复证书。所有新解释是否有经验预测力，待轨迹和反事实结果填写，当前不预宣称成功。

## 5. 可检验的有限缺陷压力：它不是成功判据

普通 CIE 的 P_i=β/(1+βm_i)，m_i>0。在固定原始本构 r 的乘积微分中，权重导数那一项单独给出

```text
(a_r)_i = −β³ sum_s w_s |r_is|²/(1+βm_i)³ ≤ 0.
```

因此负的该梯度分量指向更大 contrast；它不读取接收目标 Y，也不判断这一升高是否正确。对 self 形式 b_i=β−G_ii、d_i=1+b_i m_i，有

```text
(a_r,self)_i = −β² sum_s w_s |r_is|²
               (Re b_i + m_i |b_i|²)/|d_i|⁴.
```

若 Re b_i≥0，这一分量也非正。证明直接微分 β²|r_i|²/(2|d_i|²)，此处只取 P 的导数项；**总梯度还包含 r 随 m 的变化、prior，并受 box 和电流更新影响**。不能把分量的负号直接命名为算法发散、材料偏差或大contrast必错。

此式可检验两件事：根处解释是否遗漏主要有限缺陷贡献，以及这种上推压力是否在成功/失败轨迹中都存在。如果两类都有，则它只能解释共同的动力学因素，不能作为分类成功的预测规则。已知 CIE/NIE 文献的大β饱和/读出病态是相关先例，但不自动等于上述具体梯度分解；参见文献账本，暂不主张首次。

## 6. Phase 值与轨迹因果是两个问题

首批168个已冻结检查点的全域cross均为正；因此本轮旧成功不能由“负cross降低辅助能量”解释。正交叉罚项也可能通过其梯度改变路径。冻结metric的transport与实际步夹角变化很大，并且metric导数项常不可忽略，尚无可靠成败判别。

为区别作用位置，固定四种因子组合：(普通current,普通material)、(onehop current,onehop material)、(onehop current,普通material)、(普通current,onehop material)。前两臂复用R4/本轮重放，后两臂先只在旧开发对象运行；不调β或预算，接同一原物理corrector。若某单独位置复制恢复，仅说明该实例另一位置不必要；如果均不能，则支持需要配合的实例解释，但仍不能把单例称普遍预测规律。后续是否扩展由实际结果决定。
