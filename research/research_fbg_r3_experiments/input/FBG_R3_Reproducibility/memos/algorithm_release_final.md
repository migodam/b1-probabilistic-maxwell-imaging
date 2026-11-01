# 可执行候选：受控辅助电流释放

本文件规定一个可以直接实现的小型机制算法。它没有被证明优于成熟 CIE-I，也没有运行新的成像基准。核心区别是：**把共享材料不相容的 current normal 当作有预算的探索自由度，用其对材料步的影响选择释放动作，并用完整物理重新验收。** Schur、variable projection、罚方法和 current lifting 都作为已有工具使用。

## 1. 六元组与运行状态

|项目|明确规定|
|---|---|
|优化变量|完整实材料 \(m\in[m_{\min},m_{\max}]^n\)，以及所有源的辅助电流 \(j\in\mathbb C^{n\times p}\) 在当前仿射空间中的系数。没有深度材料 latent。|
|物理约束|各单元及所有源共享同一个实材料。最终候选必须完整求解 \(L(m)j_{\rm phys}=X(m)E\)。中间态允许 state defect，但受独立预算控制。|
|目标|在每个固定阶段使用 \(\Psi_{\theta,\rho}(m,j)=\frac12\|\mathcal S j-y\|^2+\frac\rho2\|P_\theta(m)e(m,j)\|^2+\mathcal R(m)\)，其中 \(e=j-X(E+Gj)\)。探索阶段可以暂时固定一个非零新电流坐标，明确支付离开物理流形的代价。|
|更新|正确实化的 joint Gauss–Newton/LM；对 current 系数剖面化。新增模式先按真实材料位移计算候选，再用独立物理任务模型筛选。单 mode 采用坐标不变的 ridge 软释放。|
|子空间|共同 receiver-visible backbone 与投影 Fourier 模式。新增候选可来自 Fourier、随机、原 A17 字典、普通 Krylov/primal-dual 残差。所有方法使用相同候选获取权利及成本账本。|
|全局化|少量有限 normal seed、有限 rollout、branch 去重及完整 Maxwell corrector。保留最佳已验收物理 incumbent。允许辅助路径中原物理目标暂时升高。|

\(\mathcal S\) 表示接收、实化与正确噪声白化的组合；不是新增物理照明。

本文件采用标量对比度 \(X(m)=\operatorname{diag}(m)\)，\(\rho>0\)，且阶段内 \(P_\theta\) 可微并可逆；对 \(P_\beta\) 需要在允许材料 box 上保持 \(1+\beta m_i\) 远离零。若实际后端使用 DDA 极化率 \(a(m)\)，须将材料注入改为 \(DX[h]Q\)，并保留 \(a'(m)\)。

运行状态包括：物理 incumbent、多个已验收 branch、各 branch 的物理 anchor、当前 current basis、固定阶段参数 \((\theta,\rho)\)、剩余预算以及拒绝原因。

### 1.1 明确 current affine space

第一版建议在物理 anchor \(m_*\) 求得 \(j_*=j_{\rm phys}(m_*)\)，使用
\[
j=j_*+\mathcal U a.
\]
\(\mathcal U\) 可由接收可见方向 \(U_{\rm vis}\) 与 Fourier-APCS 方向共同构成。可靠接收模式受数据项直接约束；Fourier-APCS 控制域内 current 复杂度。

这不是逐字复刻“固定 DPCS + APCS”的原生 TSOM。若要直接在原生 TSOM 仿射空间 \(j_D+U_Aa\) 内应用，则必须确认物理 anchor 属于该空间；否则增加已计费的 anchor modes，或承认当前态非物理并使用下一节的一般线性化。不得把一个不属于表示空间的物理 anchor 当成已满足 state constraint。

成熟 CIE-I/FFT-TSOM 按原论文完整算法另外运行，保留其真正的 DPCS/APCS、阶段日程和双 β 能力。

## 2. 物理锚点与非物理迭代必须使用不同的线性化

写
\[
Q=E+Gj,\qquad L=I-XG,\qquad e=j-XQ.
\]
对实材料增量 \(h\) 和复 current 增量 \(\delta j\)，在 \(X=\operatorname{diag}(m)\) 时，
\[
De[h,\delta j]=L\,\delta j-\operatorname{diag}(h)Q.
\]
对加权 state residual：
\[
\boxed{
D(P_\theta e)[h,\delta j]
=P_\theta L\,\delta j
-P_\theta\operatorname{diag}(h)Q
+(DP_\theta[h])e.
}
\]
最后一项不能在非物理迭代点删除。例如取固定阶段的
\[
P_\beta(m)=\operatorname{diag}\{(1+\beta m_i)^{-1}\},
\]
则
\[
(DP_\beta[h])e
=-\operatorname{diag}\left\{\frac{\beta h_i}{(1+\beta m_i)^2}\right\}e.
\]

这里选择 \(P_\beta\) 是已声明的状态度量；若与原 CIE 的 \(\beta/(1+\beta m)\) 残差比较，须连同全局归一化和 \(\rho\) 对齐。不同 β 阶段的目标改变必须明确记录，不能跨阶段复用旧 BFGS secant 并声称仍在优化一个固定目标。

固定 current 基 \(\mathcal U\) 后，正确实化的联合残差线性模型为
\[
b+B h+A_Uc,
\]
其中
\[
b=
\begin{bmatrix}
\mathcal Sj-y\\
\sqrt\rho\,\mathcal R(P_\theta e)
\end{bmatrix},
\]
\[
A_Uc=
\begin{bmatrix}
\mathcal S(\mathcal Uc)\\
\sqrt\rho\,\mathcal R(P_\theta L\mathcal Uc)
\end{bmatrix},
\]
\[
Bh=
\begin{bmatrix}
0\\
\sqrt\rho\,\mathcal R\{(DP_\theta[h])e
-P_\theta\operatorname{diag}(h)Q\}
\end{bmatrix}.
\]
\(\mathcal R\) 表示一致的实化。一个复空间 current atom 在 \(p\) 个源下通常产生 \(2p\) 个实系数，不得按一个实自由度计价或比较。

### 2.1 在物理 anchor

若 \(e=0\)，则 \(DP_\theta[h]e=0\)，可以使用简化 material injection。此处计算的物理 tangent/normal 是严格的。

### 2.2 在非物理 rollout

一般 \(e\ne0\)，必须使用完整 \(B\)。当前 \(Q=E+Gj\) 只是辅助 current 产生的场，并非 \(m\) 的完整物理总场。

若此时仍用 \(\Pi_Q^\perp L\delta j\)，只能称当前线性 state-compatibility 的诊断量，不能直接称物理流形的法向投影。第一版在物理 anchor 建立 normal 字典；rollout 中刷新正确的 \(b,B,A_U\)，下一次完整物理投影后才重建几何 anchor。

例如非物理点的 real-CSI readout \(\widehat m\) 还满足
\[
D\widehat m[V]
=\frac{\operatorname{Re}\langle q_i,[(I-\operatorname{diag}(\widehat m)G)V]_i\rangle
+\operatorname{Re}\langle (GV)_i,(e_{\rm CSI})_i\rangle}{\|q_i\|^2}
\quad\text{（第 \(i\) 行）},
\]
不能无条件沿用物理点中删除第二项的公式。

## 3. Normal 指标是什么，怎样避免只按电流能量选模

物理 anchor 每行 \(q_i\ne0\) 时，定义
\[
\Pi_{q_i}^{\perp,\mathbb R}v
=v-q_i\frac{\operatorname{Re}\langle q_i,v\rangle}{\|q_i\|^2},
\qquad
N_*(V)=\Pi_{Q_*}^{\perp,\mathbb R}L_*V.
\]
这是共享实材料不相容的一阶缺陷。物理 tangent \(V=T_h=L_*^{-1}\operatorname{diag}(h)Q_*\) 满足 \(N_*(V)=0\)。

normal 分数也必须消除 retained-space 的坐标任意性。令 \(\mathsf L_*\) 表示 \(V\mapsto L_*V\) 在所有源上的实化算子，\(M_j=\mathsf L_*^T\mathsf L_*\) 是正定 current 度量，即 \(\|V\|_{M_j}=\|L_*V\|_F\)。先取候选的商空间代表
\[
\mathcal C_q=(I-\Pi_{\mathcal U}^{M_j})\mathcal C,
\qquad
W_C=\mathcal C_q(\mathcal C_q^TM_j\mathcal C_q)^{-1/2},
\]
其中只保留非零秩 \(r_C\)，\(\Pi_{\mathcal U}^{M_j}\) 是 retained current range 在该度量下的正交投影。定义
\[
\nu_C=\frac1{r_C}\|N_*(W_C)\|_F^2.
\]
在使用一致的逐行源内积时 \(0\le\nu_C\le1\)。整块白化后 \(W_C^TM_jW_C=I\)，所以对 \(\mathcal C\mapsto\mathcal CT+\mathcal UK\) 的可逆重参数化，该分数不变；仅逐列归一化不够。它区分模式的角色，不将 \(\nu_C\) 最大者自动视为最好。总场节点没有合法逐行除法，应单独保留状态/材料变量并标记该行，而不是加一个未经解释的 epsilon。

**命题 B：分组 trace 可能完全没有选择能力。** 假设一个空间 atom 对全部 \(p\) 源各有独立复系数，使用相同空间算子和等权源 Frobenius 度量；retained quotient 仍按源可分；各相关行 \(q_i\ne0\)。则上述整组均值恰为 \(1-1/(2p)\)。

**证明：** 取经过商投影后的空间 atom \(c\)，令 \(a=L_*c/\|L_*c\|\)，所以 \(\sum_i|a_i|^2=1\)。每行 \(p\) 个复源的完整系数基有 \(2p\) 个实方向。共享实材料只占 \(\operatorname{span}_{\mathbb R}\{q_i\}\) 的一个实方向，因此 \(\operatorname{tr}\Pi_{q_i}^{\perp,\mathbb R}=2p-1\)。乘以复数 \(a_i\) 在该实系数空间中表示为 \(|a_i|O_i\)，其中 \(O_i\) 正交。整组 normal 能量为
\[
\sum_i|a_i|^2
\operatorname{tr}\!\left(O_i^T\Pi_{q_i}^{\perp,\mathbb R}O_i\right)
=(2p-1)\sum_i|a_i|^2=2p-1.
\]
整组总能量为 \(2p\)，故 \(\nu_C=(2p-1)/(2p)\)。多个按该空间度量正交的 atom 逐组相加，仍得到同一均值。证毕。

所以此类 Fourier 组之间不能用 \(\nu_C\) 排序。它仅区分自由 current 组与共享材料耦合的 physical tangent 组。缺少源全基、源间系数受约束、零场行或不同度量时，不可直接套用这个常数公式。

实际算法必须先求候选产生的联合动作。令
\[
V_C=\delta j_\gamma-\delta j_{\rm old},
\qquad
\nu_C^{\rm act}=
\frac{\|N_*(V_C)\|_F^2}{\|L_*V_C\|_F^2},
\]
其中 \(V_C=0\) 单列为 no-op。完整 current 动作而非裸候选系数也消除了 retained-space dressing 的任意性。finite seed 则使用 seed 及旧 current 优化共同产生的真实 current 位移。对这一实际预测动作 \(V_C\)，还可记录 real-CIE 有限偏置的二阶预测
\[
[b_C]_i=
\frac{\beta\|N_*(V_C)_i\|^2}
{\|q_i\|^2(1+\beta m_{*,i})}.
\]
若当前独立物理模型给材料梯度 \(\widehat g_F\)，则
\[
a_C=-\widehat g_F^T b_C
\]
描述这一偏置与当前任务下降方向是否一致。正值只是一阶任务对齐诊断，不能当成 finite basin 证明。若 \(\widehat g_F\) 是 full adjoint，必须支付每源 adjoint RHS；若来自共同 reduced anchor，所有对照共享该 anchor 并记录它的获取成本。

**实际选模用原物理任务收益排序，normal量作为角色标签和预算约束。** 不把 current energy、\(\nu_C\) 或 \(a_C\) 单独变成接受规则。

## 4. 模式释放的可计算材料后果

设旧 current range 的 projector 为 \(P_U=A_UA_U^\dagger\)，材料正则/阻尼局部模型为
\[
\tfrac12h^T\Lambda h+\ell^Th,\quad \Lambda\succ0.
\]
旧剖面模型：
\[
q_U(h)=\tfrac12\|(I-P_U)(b+Bh)\|^2
+\tfrac12h^T\Lambda h+\ell^Th.
\]
令
\[
H=B^T(I-P_U)B+\Lambda,\quad
g=B^T(I-P_U)b+\ell,\quad
p=-H^{-1}g.
\]
新增候选 \(\mathcal C\) 对应联合残差矩阵 \(A_C\)。对
\[
Z_C=(I-P_U)A_C
\]
做 rank-revealing 薄 QR，仅保留非零列空间：
\[
Z_C=Q_CR_C,\qquad Q_C^TQ_C=I.
\]
定义
\[
D_C=Q_C^TB,\qquad v_C=Q_C^Tb.
\]
硬释放给
\[
H_C=H-D_C^TD_C,\qquad
g_C=g-D_C^Tv_C,
\]
\[
p_C-p=H_C^{-1}D_C^T(v_C+D_Cp).
\]
计算的是冻结联合 GN/二次模型，不是任意非线性 profile Hessian。

固定 box 活动面上可用相同公式。若预测步触发新的活动约束，应求实际约束 QP，并按实际执行位移重新打分，不把无约束 \(p_C\) 的收益转移给 clip 后的步。

## 5. 单 mode 不能靠缩放 basis 逐渐释放

### 命题 A：A17 dressing 本身不产生新的增广空间

若 retained current 空间为 \(V\)，A17 dressed candidate 满足
\[
\bar C=C+VK,
\]
则
\[
\operatorname{span}[V,\bar C]=\operatorname{span}[V,C].
\]
在相同冻结联合残差导数下，
\[
A_{\bar C}=A_C+A_VK,\qquad
(I-P_{A_V})A_{\bar C}=(I-P_{A_V})A_C.
\]
因此硬释放产生完全相同的 profile Hessian、梯度和精确材料步。对原非线性 joint objective，只要 current/材料可行集合和正则保持相同，两个 current 坐标也表示同一个搜索集合。

**证明：** 新旧候选仅相差 retained-space 列组合；对应的系数变换可逆。投影消去 retained range 后该差为零。证毕。

A17 可能通过**选择不同的 raw candidate、同秩替换、有限任务评分或更便宜的更新**产生价值。不能把同一个 \(C\) 换成 \(\bar C\) 后就声称释放了新的物理方向。

同样，将一个非零 mode 乘以 \(t\ne0\) 不改变 projector。列归一化使某个固定 coefficient ridge 改变有效惩罚时，所得差异属于未运输的正则，并非新的 Maxwell 效应。

### 5.1 使用 QR 残差坐标做真正连续的软释放

把当前线性模型中的新增 current 增量系数记为 \(a_C\)，定义
\[
u_C=R_Ca_C.
\]
在这组坐标中，新模式对剖面联合残差的贡献是 \(Q_Cu_C\)。加入
\[
\frac\gamma2\|u_C\|^2,\quad \gamma\ge0
\]
作为明确的数值释放控制。它等价于原系数上的 \(\Gamma_C=\gamma Z_C^TZ_C\)，而不是任意 basis 中固定 \(\gamma I\)。

消去 \(u_C\) 得
\[
u_C^*(h)=-\frac{v_C+D_Ch}{1+\gamma}.
\]
令 \(\alpha=(1+\gamma)^{-1}\)，则
\[
\boxed{q_\gamma(h)=q_U(h)-\frac\alpha2\|v_C+D_Ch\|^2,}
\]
\[
\boxed{
H_\gamma=H-\alpha D_C^TD_C,\quad
g_\gamma=g-\alpha D_C^Tv_C.
}
\]
所以新增自由度仍产生 PSD downdate，但可以通过 \(\gamma\) 连续控制：
\[
p_\gamma-p
=H_\gamma^{-1}\alpha D_C^T(v_C+D_Cp).
\]
若
\[
\eta_0=\|D_CH^{-1/2}\|^2,\quad
\eta_\gamma=\alpha\eta_0<1,
\]
则
\[
\|p_\gamma-p\|_H
\le
\frac{\alpha\sqrt{\eta_0}}{1-\alpha\eta_0}
\|v_C+D_Cp\|.
\]
给定暂时的步稳定目标 \(0<\eta_{\rm target}<1\)，可以直接取
\[
\boxed{\gamma=\max\{0,\eta_0/\eta_{\rm target}-1\}.}
\]
这不是“最优 γ”理论，只是一个明确可复现的软释放尺度。第一版将此 ridge 用作每次 joint GN 的局部步正则；原非线性阶段目标仍是 \(\Psi_{\theta,\rho}\)。若改为惩罚累计 current 系数，则必须另写出相应固定阶段目标及其梯度，不能继续照搬上述零中心增量公式。候选未释放时只优化旧 current，释放后每次重新构造正确线性模型，再用该局部 ridge 控制新增坐标增量。

QR 基可以改变正交旋转而不改变惩罚。对任意可逆新系数变换 \(C\mapsto CT+VK\)，只要同时运输 \(R_C\)，该 penalty 与 profile 更新不变。另一种允许的做法是按真实新增 current correction 的物理范数定义 penalty，并使用其 Gram 运输；不能在不同归一化系数中固定同一裸 ridge 后宣称公平。

### 5.2 大 \(\eta\) 不等于模式有害

单单元 \(G=0,S=\epsilon,L=B_m=1\) 中，
\[
H_{\rm old}=\rho+\lambda,\quad
H_{\rm full}=\frac{\rho\epsilon^2}{\rho+\epsilon^2}+\lambda.
\]
\(\epsilon\to0\) 时，旧曲率大多是冻结 current 产生的假刚性。释放后曲率下降反而使模型更诚实。

因此 \(\eta_{\rm target}\) 只控制一个局部 trial。不能永久拒绝全部大 \(\eta\) 模式。可以软释放后刷新线性化、修改 \(\rho/\beta\)、扩大允许步或进入 finite exploration。对真实信息仍另看完整 \(\mathcal J_F\)；大的 \(\Lambda\) 也不能冒充信息。

## 6. Normal-aware 选模必须付出真实物理任务评价成本

每个 candidate action 保存：
\[
\mathcal A_C=(p_\gamma,\delta j_\gamma,\nu_C,\nu_C^{\rm act},b_C,a_C,\eta_\gamma,
\text{predicted gain},\text{incremental cost}).
\]

### 6.1 廉价预筛

使用一个共同、独立于候选 current basis 的 reduced physical anchor \(\widehat{\mathcal J}_F\)。令 \(s_0=(m-m_*)+p\) 是旧 trial 相对物理 anchor 的总材料位移；选择动作恰在 anchor 时 \(s_0=p\)。写该 anchor 的原物理二次模型为
\[
q_F^*(s)=\tfrac12\|r_*+\mathcal J_Fs\|^2
+\ell_F^Ts+\tfrac12s^T\Lambda_Fs.
\]
\(\ell_F,\Lambda_F\) 来自共同原物理 prior 的局部模型；joint-step 的 LM 阻尼不计入此物理收益。对不光滑 prior 使用同一精确 prior 差值。廉价版本给
\[
\widehat{\mathcal G}_C
=
-\widehat g_F(s_0)^T(p_\gamma-p)
-\tfrac12(p_\gamma-p)^T\widehat H_F(p_\gamma-p).
\]
按单位实际增量成本的 \(\widehat{\mathcal G}_C\) 排序，同时记录 normal activation 与材料偏置对齐。该 reduced anchor 可复用，但所有方法支付同样 setup；不得把 Full-J 免费送给新 selector。

如果 anchor 没有可靠误差界，该分数只是 proposal，不是证书。

### 6.2 最终少量 JVP 核验

在物理 anchor 的线性模型中，缓存
\[
u=r_*+\mathcal J_Fs_0.
\]
对 finalist 位移 \(\Delta p=p_\gamma-p\)，计算
\[
z=\mathcal J_F\Delta p.
\]
原完整物理二次模型收益为
\[
\boxed{
\mathcal G_C
=-\operatorname{Re}\langle u,z\rangle
-(\ell_F+\Lambda_Fs_0)^T\Delta p
-\tfrac12\|z\|^2-\tfrac12\Delta p^T\Lambda_F\Delta p.
}
\]
这里 residual 与 JVP 已一致白化。每个 JVP 需要**全部 \(p\) 个源的 tangent RHS**，利用 anchor LU 只省 factorization，不会把 RHS 变成零。缓存 \(u\) 的 \(\mathcal J_Fs_0\) 同样支付一次 \(p\)-RHS block。无 source-separable covariance 时必须按完整数据白化。若已远离 anchor，这仍仅是旧 anchor 的物理二次预测；须重新付费重锚，或将其降为不保证精度的提案分数。第一版仅在已验证物理 anchor 选择下一轮释放，非物理 rollout 内不再次选新组。

A17 所保留的 residual 与 curvature 两侧作用在这个 gain 中都存在。单看 normal energy 或反馈场衰减无法替代它。

### 6.3 两个预算槽，避免把全局搜索变成局部接受

**局部槽：**用已核验 \(\mathcal G_C\) 选择最好的正收益动作；再用真实非线性物理目标检查最终 material trial。

**探索槽：**每个 episode 最多保留一个具有非零 normal activation、预定有限材料移动尺度、合法 box 和可控 defect 的候选组。它可以有 \(\mathcal G_C\le0\)。在同一预定移动尺度内仍按物理任务预测排序，不能通过任意缩小步长把 no-op 包装成成功。

探索槽的首要约束是可实现性风险和预算；实际动作的 normal activation 决定它是否确实测试了一个不同机制，不是越大越好。记录 \(\mathcal G_C\)、\(\nu_C\)、\(\nu_C^{\rm act}\)、\(a_C\) 的独立数值，不把它们合成事后调出来的万能评分。

## 7. 必须显式启动有限探索；只开放模式不会自动跨根

若 anchor 已是 exact-data 解、prior 梯度为零，则 \(b=0\)、\(p=0\)，无论新增多少自由 current，局部 joint 最小二乘仍给零步。**开放 normal mode 本身不生成第二个解。**

第一版使用一个明确的 finite normal seed：

1. 对已选候选的 QR 坐标，若 \(v_C+D_Cp\ne0\)，取其归一化方向；否则取小矩阵 \(D_CH^{-1}D_C^T\) 的最大特征向量。它选择会影响材料的 current 坐标，而不单按 current energy。
2. 分别测试 \(u_C=+\tau w\) 和 \(u_C=-\tau w\)。两支有相同幅度和支付上限。
3. 对固定 \(u_C\)，联合局部 material step 为
   \[
   h(u_C)=p-H^{-1}D_C^Tu_C.
   \]
   用预注册目标 \(\|h(u_C)-p\|_{M_\chi}=\delta_m\) 选 \(\tau\)，并要求 current/defect 预算均满足。若一阶材料作用为零，可将它保留为一个专门的二阶机制对照；不要把无材料移动算成有效 branch。
4. 固定这一新 current 坐标最多若干 inner iterations，优化材料和旧 current。固定的是 seed 时 \(a_C=R_{C,*}^{\dagger}u_C\) 所表示的物理 current 分量，保持 seed 时的映射 \(R_{C,*}\)；不能重算 QR 后继续固定同一数值 \(u_C\)，从而暗中改变 seed。这样明确允许走过一个有限 state-defect 障碍。
5. 解除固定坐标，用前述每步坐标不变的局部 ridge 和逐步提高的物理可行性要求运行有限 rollout。在固定阶段内以真实 \(\Psi_{\theta,\rho}\) 决定一般 inner step 是否接受；seed 插入及明示的阶段切换另行计费，不要求每步原物理目标下降。
6. 用所得材料重新求解完整 Maxwell，并在该物理材料上执行公共 corrector。若回到原 anchor，记 duplicate；若不满足物理条件，记 fail；两者全部计费。

固定 current 坐标是探索约束，不是物理先验，也不是新增 illumination。它可以暂时增加原物理数据目标和联合目标。若拒绝所有此类离开，exact root 的多分支探索在定义上就不可能发生。

这套有限 seed 也不能保证找到另一支。对照必须包含相同幅度、相同预算的随机/Fourier current seed 和普通材料多起点。

## 8. 可执行伪代码

~~~text
输入:
  E, G, S, Y, noise metric, material box and common prior
  physical incumbents/branches, current backbone U
  fixed candidate pool rule, chart family Pθ
  budget B, fixed shortlist and exploration limits

for each permitted physical anchor m*:
    j* <- FULL_MAXWELL(m*)                     # 计 LU 与 p 个 RHS
    set j = j*, affine base = j*
    build/refresh current backbone and common proposal reference
    account for all dictionary/reference acquisition

    for each permitted release episode:
        restore a verified physical anchor before selecting a new group
        form EXACT b, B, AU at current (m,j)
        # 非物理点中 B 含 DP[h]e
        solve old constrained joint GN/LM block
        on a fixed active face obtain H, g, p

        for each candidate group C in the fixed pool:
            Z <- (I - P_AU) AC
            rank-revealing QR; drop exact duplicate residual directions
            D <- QZ^T B; v <- QZ^T b
            eta0 <- ||D H^(-1/2)||²
            gamma <- max(0, eta0 / eta_target - 1)
            compute softly released step pγ and actual current action
            compute quotient-normal signature and finite-bias diagnostic
            score predicted PHYSICAL task gain using common reference
            record all actions and candidate costs

        verify at most the declared finalists with FULL JVPs
        local slot <- highest valid positive task-gain action, if any
        exploration slot <- best permitted finite-scale normal candidate
                            even if its local physical gain is nonpositive

        execute local trial under its declared nonlinear acceptance rule
        for signs ± in the permitted exploration seed pair:
            seed a finite current coordinate at declared material scale
            hold seed in its original physical current coordinates
            optimize for a small fixed number of iterations
            release via coordinate-invariant local increment ridge
            run only the allowed joint rollout, with exact derivatives
            if basis/chart changes: reset local derivative history
            if budget/defect/material constraints fail: stop this branch

            m_trial <- current material
            j_phys <- FULL_MAXWELL(m_trial)
            run common physical corrector within remaining budget
            evaluate original Maxwell loss and physical admissibility
            update best physical incumbent only after this validation
            store distinct qualified material candidate or record failure

        merge only declared material/response near-duplicates
        keep at most Kmax physical branches with diversity
        stop if episode budget exhausted or no permitted action remains

输出:
  best verified physical incumbent
  distinct verified material candidates with qualification labels
  full failed/duplicate/rejected costs and mechanism diagnostics
~~~

## 9. 默认小型机制预算及停止规则

以下只是下一轮 E0/E1 的可执行初始配置，不表示已经运行，亦不是 576 参数基准的最终预算：

- 每个 anchor 固定候选池 8 组。各组使用相同数量的实 current 自由度。
- 每个 episode 同时新增至多 1 组；下一次物理重锚后恢复共同 backbone，再重新选组。这样保持一个明确的 active-rank 上限，避免把累计扩容误报为 selector 收益。固定秩 replacement 是可另测分支，须另计剔除模式后的重解成本。
- 每个 episode 最多核验 2 个 JVP finalist；缓存旧步的 JVP 最多一次。没有额外免费 Full-J。
- 每个 anchor 最多 2 个 release episodes。
- 每个 episode 最多一对 \(+\tau/-\tau\) finite seeds。
- 每个 seed 最多 3 次固定坐标 inner updates、9 次解除后的 joint updates。
- 每个 seed 最多 2 个完整物理 checkpoint，包含末端 corrector 的评估；若线搜索需要更多 full evaluations，必须从同一预算扣除。
- 存活物理 branch 数上限 \(K_{\max}=3\)。 incumbent 永不因探索失败被覆盖。
- 预注册总 full-state、LU、复 RHS、Green-vector action 和 wall-clock 上限；任一耗尽立即停止该轮。不同矩阵大小和 source 数要先校准这些 primitive costs，不能将“一次 block action”都当同样费用。
- 使用同一个 final physical corrector、相同停止阈值与相同最终选择规则。Best-of-K oracle 仅用于离线报告，不能冒充可部署选择器。

state-defect 预算必须与数据域传播风险同时解释。若没有 \(\|\mathcal S L^{-1}\|\) 的可靠界，\(\|e\|\) 小只是一项诊断；full checkpoints 才是物理判断。CC-CSI 的交叉项 \(R-Se\) 使用背景 receiver map，不能当作完整 \(\mathcal S L(m)^{-1}e\) 误差证书。

## 10. 最小决定性实验：让算法与近邻真正可区分

### E0a：表示等价与软释放

在同一个合成复 VIE 上比较：

1. raw C；
2. A17 dressed \(C+VK\)；
3. \(CT+VK\) 的非正交重参数化。

硬释放必须给相同 profile update。使用 QR 残差坐标 ridge 后也应相同。裸 coefficient ridge 若不同，应标为改变正则，而非 A17 收益。

已完成的独立脚本 supplement_checks.py：

|量|结果|
|---|---:|
|dressing 后 projected Z 差|2.38×10^-16|
|硬释放 Hessian 差|1.23×10^-16|
|硬释放材料步差|2.47×10^-17|
|坐标不变软释放 Hessian 差|6.33×10^-17|
|坐标不变软释放材料步差|7.37×10^-18|
|裸 coefficient ridge 重参数化后的 Hessian 差|0.05335|
|裸 coefficient ridge 重参数化后的材料步差|0.02085|
|硬释放 η|0.47657|
|按目标 η=0.05 选 γ|8.53138|
|实际软释放 η|0.050000000000000024|
|current 商空间 normal trace 的重参数化差|0（本次浮点结果）|
|不取商空间时 normal trace 的 dressing 差|0.00581|

n=7、p=2、r=3 的随机复对称 VIE，仅验证代数，不是成像基准。本例按源可分的单空间 atom 整组 normal trace 为 0.75，恰等于 \(1-1/(2p)\)，同时验证了该组均值不能承担选模任务。

### E0b：normal 机制与完整物理增益

用小型真实二维 TM 后端或精确两单元/少单元 VIE，先固定同一材料、观测、噪声、current backbone、候选池及成本。

比较：

- 随机 current 选模；
- Fourier/current energy 选模；
- receiver singular value 选模；
- 原 A17 physical tangent-step relevance；
- 新方法：正确 joint-profile task relevance，normal预算及软释放；
- 可支付的小模型 Full-J/完整 current 参考。

所有选模使用同一个候选池的版本是主消融；另做各自原生候选生成的完整成本比较，不能混为一个结果。

评价：

- 预测材料步差与实际步差；
- hard/soft \(\eta\)、当前 residual 和原 physical quadratic gain；
- normal activation 与 CIE 二阶/四阶偏置；
- full Maxwell 后的 data residual 和材料误差；
- 每个原子的产生、未入选评分、失败和验证费用。

若新方法只获得更大 normal energy、较低辅助目标，却在真实物理任务上不占优，应停止该 selector。

### E0c：跨盆地不能由局部下降伪装

同预算比较：

- 不释放 normal，仅物理 Full-J corrector；
- 开放 normal 但没有 finite seed；
- 固定 \(+\tau/-\tau\) Fourier seed；
- 固定 \(+\tau/-\tau\) 随机 seed；
- 以上任务相关 finite normal seed；
- 相同材料距离的普通 material 多起点；
- 小支持 exact-response/Z-ALG 枚举提案。

所有臂使用相同 rollout 次数与末端物理corrector。在 exact root 上，“只开放 normal”应给零局部步；它是必要负控制。报告新的物理解数量、重复率、资格失败和有限预算成功率。已知另一根只能作为离线 evaluator，不能进入 seed 构造。

## 11. 新贡献与停止条件

本算法不会因为写出了 profile Schur complement 就成为新算法。必须证明的是：

1. 共享实材料的 normal 定义和有限 CIE 偏置，能指导选择比随机、Fourier能量或单纯物理tangent评分更有价值的辅助方向；
2. 这种选择不是 dressed/raw 同空间的坐标效应，也不是额外 ridge 的不公平差异；
3. 有限探索产生的收益超过同预算普通多起点；
4. 完整原物理和 material 合法性验收后仍有优势。

若没有这些增量，保留成熟 CIE/FFT-TSOM 或普通 joint-current 方法。若 normal 选择有增量但 A17 的候选/评分没有独立作用，论文应诚实以共享材料法向与受控释放为中心，而不强行宣称 A17 的算法继承关系。

