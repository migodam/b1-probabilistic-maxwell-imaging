# 06 · 由算子结构推出的优化与裁决

本章回答“固定什么才真的减少计算”。优先级由第02章的精确结构决定：先消去固定背景，复用多源 factorization 和接收端口，再决定是否需要更低维近似。冻结四维弱方向只是搜索策略；如果仍然每步计算相同的完整 Jacobian，它并未自动减少物理线性求解。以下候选均区分恒等式、近似与尚未证实的速度收益；没有新增物理实验。

## F-OPT-001 · 固定宿主消元是所有搜索方法共享的强基线

由 F-OP-008，原非均匀 Full66 支持上的完整预测为

$$L_Dh=DE_R,\qquad L_D=I-DC,\qquad F=F_h+Oh.$$

固定项是 $C\in\mathbb C^{256\times256}$、$E_R\in\mathbb C^{256\times8}$、$O\in\mathbb C^{32\times256}$、$F_h$。每个材料只更新对角 $D$，解256阶完整系统。它保留原 decoder、噪声模型和 likelihood，所以 fixed-B、dynamic-C、随机探索、普通 Full-J corrector 以及概率积分都应共享它；不能只把这一优化赋给欲证明优越的一臂。

相比576阶稠密 LU，其 leading factorization 规模由 $N^3$ 变成 $q^3$；这不是实际速度比，因为 setup、BLAS、端口收缩、I/O 和状态数量均有成本。已发生的16个暴露 DEV 状态、3次计时重复的核比较，只支持该窗口的效率判断。它不是新的代数身份，也不是新对象上的完整算法加速。

**直接决策：** 后续方法强基线应为精确256 ROI solver。只有原有瓶颈转移后仍存在的成本，才值得进一步压缩。固定算子低秩近似与真实 residual 校正见 F-ERR-002—006。

## F-OPT-002 · 多源、多参数与接收端口：三种导数成本不同

由 F-OP-010，$\delta F=T D_{\delta d}E$，其中 $T=OL_D^{-1}$，本节 $E=E_{ROI}$。同一材料的所有源和所有导数共享一个 $L_D$ factorization。

1. **完整66列 Jacobian：** 8 forward RHS 得 $E$，32 transpose RHS 得 $T$，再与 $M=\partial d/\partial x$ 收缩。传播费用为40个 reduced RHS，而非 $8\times66$ 个 tangent RHS。
2. **仅一个实标量目标的梯度：** 设 $\Phi=\|F-Y\|_F^2/\sigma^2$，解

$$L_D^*\Lambda=O^*(F-Y)/\sigma^2,\qquad c_i=\sum_s E_{is}\overline{\Lambda_{is}}.$$

则 $d\Phi=2\operatorname{Re}\sum_i c_i\,\delta d_i$，所以

$$\nabla_x\Phi=2\operatorname{Re}(M^Tc).$$

只需8个 adjoint RHS；若含标准 Gaussian prior，梯度再加 $x$。这里 $M^T$ 而非 $M^*$ 与定义的 $c$ 配对；也可采用全共轭的等价写法，但不能混配。
3. **矩阵自由 GN 作用：** 一个 $v\mapsto\mathcal J^T\mathcal Jv$ 通常需8 tangent RHS 加8 adjoint RHS，在已有 forward/factorization 之后约16个 reduced RHS。对66维当前问题，若每步迭代次数超过少数次，可能比一次32 RHS构造全 Jacobian更贵；需要计入实际迭代、预条件与收缩成本。

因此，“矩阵自由”“只更新少数方向”不天然更便宜。应先写清目标到底需要梯度、全 Jacobian，还是某几个作用，再选择算子。不同材料状态通常需要不同 factorization；不同源在同状态则不应重复 factorization。

## F-OPT-003 · 电流残差消元与 variable projection 的精确边界

原材料目标可以等价写成带硬物理约束的联合问题

$$\min_{m,J_c}\ \|WSJ_c-WY\|_F^2+R(m),\qquad L_mJ_c=D_mE^{inc}.$$

当 $L_m$ 可逆时，约束唯一决定 $J_c=L_m^{-1}D_mE^{inc}$。消去 $J_c$ 正好返回原 forward 目标；它没有凭空免除求解。

若改为有限罚参数 $\beta>0$，固定 $m$ 的内层是线性最小二乘：

$$\min_{J_c}\ \left\|\underbrace{\begin{bmatrix}WS\\\sqrt\beta L_m\end{bmatrix}}_{A_\beta(m)}J_c-\underbrace{\begin{bmatrix}WY\\\sqrt\beta D_mE^{inc}\end{bmatrix}}_{b_\beta(m)}\right\|_F^2.$$

在 $L_m$ 可逆时 $A_\beta$ 满列秩，故

$$J_*(m)=A_\beta^\dagger b_\beta,\qquad \psi_\beta(m)=\|(I-A_\beta A_\beta^\dagger)b_\beta\|_F^2.$$

这是标准 variable projection。令 $r_\beta=A_\beta J_*-b_\beta$，由 $A_\beta^*r_\beta=0$，微分时含 $dJ_*$ 的项消去：

$$d\psi_\beta=2\operatorname{Re}\operatorname{tr}\{r_\beta^*[(dA_\beta)J_*-db_\beta]\}.$$

这给出准确的 reduced 梯度，并解释为何联合变量的内层可消元。但有限 $\beta$ 允许不满足 Maxwell 的电流；它改变了 likelihood。$\beta\to\infty$ 在适当条件下逼近硬约束，同时可能恶化条件数。用它产生初值后再做原物理验证，可以作为优化流程；不能把有限罚目标直接称为原 posterior。此处未主张它比256阶精确求解更快。

实现上优先 QR/SVD 或合适的增广求解，而不是无条件形成 $A_\beta^*A_\beta$；后者平方条件数。若原材料用66 latent 表示，增加每源256复 current 会引入4096实状态变量，变量数上未必占优。

## F-OPT-004 · 能否只解接收可见的少量端口

为说明障碍，令固定接收的完整正交基为 $[V_o,V_n]$，$V_o\in\mathbb C^{q\times r}$，$V_n\in\mathbb C^{q\times(q-r)}$，$OV_n=0$。写 $h=V_oa+V_nb$，原约束为

$$L_DV_oa+L_DV_nb=DE_R.$$

因为 $L_D$ 可逆，$B_D=L_DV_n$ 满列秩。取其左零空间正交基 $U_D\in\mathbb C^{q\times r}$，即 $U_D^*B_D=0$。乘上 $U_D^*$ 消去 $b$：
$$\underbrace{U_D^*L_DV_o}_{T_D^{port}}a=U_D^*DE_R.$$

$T_D^{port}$ 可逆。若 $T_D^{port}a=0$，则 $L_DV_oa\in\operatorname{range}(L_DV_n)$，由 $L_D$ 可逆得 $V_oa\in\operatorname{range}V_n$，正交性迫使 $a=0$。求得 $a$ 后

$$b=B_D^\dagger(DE_R-L_DV_oa),\qquad F=F_h+OV_oa.$$

这是精确端口消元，但 $U_D$ 随材料变化；每次为 $q\times(q-r)$ 矩阵求 QR/零空间可能花费与原 LU 同阶的工作。把最终 $r\times r$ 方程很小当作整体计算很小，是漏记消元成本。固定 $CV_o,CV_n$ 可以缓存，$D$ 的行缩放也很便宜，但不能因此免费得到变化的 $U_D$。

另一种块 Schur 方法需要 $V_n^*L_DV_n$ 可逆，这是比 $L_D$ 可逆更强的假设：矩阵 $\begin{bmatrix}0&1\\1&0\end{bmatrix}$ 可逆，而任一单坐标对角块均为零。不可未经检查使用该块逆。

真正免费的固定闭合需要额外结构，例如 F-OP-006 的 $P_oDCP_n=0$；原非均匀材料族没有给出这种闭合。**结论：接收 rank 小是输出维数事实，不能单独证明可同样便宜地求解非线性状态。**

## F-OPT-005 · 精确可达子空间与近似端口压缩

均匀幅度 $D=vI$ 时，固定可达子空间

$$\mathcal R=\operatorname{span}\{E_R,CE_R,\ldots,C^{q-1}E_R\}$$

包含所有源列，并由 Cayley–Hamilton 定理对 $C$ 不变。若 $Q_c$ 为其正交基，$C Q_c=Q_cC_c$、$E_R=Q_cB_c$，则

$$h=Q_c\,v(I-vC_c)^{-1}B_c,\qquad F=F_h+OQ_c\,v(I-vC_c)^{-1}B_c.$$

这是精确可达部分消元；若可达 rank 仍为 $q$，没有降维收益。进一步删除不可观测状态需要正确的不变子空间商或最小实现条件，不能仅把 $\ker O$ 删掉。有限阶 Krylov 截断通常只是矩匹配近似，除非已达到完整不变闭包。

非均匀材料则需一个对所有允许 $D C$ 共同不变、且包含所有 $DE_R$ 的固定空间。受限 decoder 也许有特殊结构，但“latent仅66维”本身并不证明这种不变性；对角材料变化可激发许多独立空间方向。其一般障碍见 F-ERR-001。

**实际优化路线：** 若没有精确小闭包，只把固定压缩作为预条件/初值，计算真实 residual，再依 F-ERR-003—006 校正或回退。判据应是接收或任务误差以及总费用，不是单独的 $C$ 奇异值尾。当前256阶系统已经很小，不预设低秩路线必胜。

## F-OPT-006 · 全列稳定基的 GN 求解，不必先寻找弱方向

在当前实 latent $x$，白化 residual 为 $r$、Jacobian 为 $\mathcal J$。标准 Gaussian prior 下 GN 二次模型是

$$q(p)=\tfrac12\|r+\mathcal Jp\|^2+\tfrac12\|x+p\|^2=\tfrac12\|Ap+b\|^2,$$
$$A=\begin{bmatrix}\mathcal J\\I\end{bmatrix},\qquad b=\begin{bmatrix}r\\x\end{bmatrix}.$$

取任意满列秩变化基 $T\in\mathbb R^{66\times k}$，令 $p=T\xi$，正确问题是

$$\min_\xi\left\|\begin{bmatrix}\mathcal JT\\T\end{bmatrix}\xi+\begin{bmatrix}r\\x\end{bmatrix}\right\|^2.$$

只有 $T^TT=I$ 时，步长的 prior 二次项才能当作单位 metric；若 $T$ 是材料归一化或 generalized eigenvector 基，通常不能把下块 $T$ 换成 $I$。即便正交，prior 线性项也必须保留 $T^Tx$。

用薄 QR $AT=UR$，解 $R\xi=-U^Tb$；数值 rank 存疑时用 pivoted QR/SVD。这样避免显式形成正规方程导致条件数平方。若 $k=66$ 且 $T$ 可逆，这是同一个 GN 步的重参数化；若 $T$ 完整正交，奇异值和二范数条件数甚至完全不变，不能宣称旋转本身改善条件数。若 $k<66$，才是限制搜索空间，需明确遗漏下降量。

**直接决策：** corrector 可以直接用完整稳定基求解；只有 proposal、诊断、预条件或明确的低维近似需要额外弱谱。不要为了一个普通 GN 步每次额外计算无用的 generalized eigensystem。

## F-OPT-007 · 仿射约束、局部消元与最小刷新判据

若要求 $Cp=c$，先检查可行性，取一个特解 $p_0$ 和 $\ker C$ 的正交基 $Z$。所有可行步为 $p=p_0+Z\xi$，所以只需求解

$$\min_\xi\|AZ\xi+(Ap_0+b)\|^2.$$

固定弱**绝对坐标** $Q_w^T(x+p-x_0)=w_*$ 对应 $Q_w^Tp=w_*-Q_w^T(x-x_0)$，通常是仿射约束；不能总写成齐次 $Q_w^Tp=0$。同时满足两组约束应叠成一个线性系统。连续施加两个投影，只有相应正交投影对易且仿射偏置相容时才可当作交集投影。非线性半径与 retraction 的细节另见第03章。

令 $H=A^TA\succ0$、$g=A^Tb$、无约束最优步 $p_*=-H^{-1}g$。对固定可行**线性**子空间 $\mathcal S$，约束最优步 $p_S$ 是 $p_*$ 的 $H$-正交投影，且

$$q(p_S)-q(p_*)=\tfrac12\|p_S-p_*\|_H^2.$$

证明：完成平方，$q(p)=q(p_*)+\tfrac12\|p-p_*\|_H^2$。仿射可行集也有相同最近点解释。于是刷新不应只看弱空间“旋转角”：若当前固定空间已经含有有效下降步，旋转很大也可能无关；旋转很小而刚好遗漏重要下降，也可能有影响。

若现成 full Jacobian 已给出 $H,g$，可在不新增物理调用的条件下估计固定限制损失。**最小有意义的刷新触发**是：当前约束仍可行，且可释放的二次下降超过已声明的刷新成本/目标容差，同时真实步的接受比提示局部模型仍可信。它是决策准则，不是全局最优刷新定理；强非线性、硬事件或 retraction 失真可使局部收益预测失败。若已经算了全 Jacobian，刷新只增加小矩阵 CPU，不能把它虚记成新的 Maxwell RHS。

## F-OPT-008 · 消去强坐标：profile 与概率积分不相同

完整正交 $[Q_s,Q_w]$ 下写 $p=Q_ss+Q_ww$，令 $A_s=AQ_s$、$A_w=AQ_w$。固定 $w$，对强变量作 QR：$A_s=U_sR_s$。则

$$s_*(w)=-R_s^{-1}U_s^T(b+A_ww),$$
$$\min_s\|A_ss+A_ww+b\|^2=\|(I-U_sU_s^T)(A_ww+b)\|^2.$$

这给出稳定的强坐标消元与弱坐标二次模型。等价 Schur Hessian 是 $H_{ww}-H_{ws}H_{ss}^{-1}H_{sw}$，但计算上可用 QR 避免显式平方条件数。

若把同一严格二次函数作为 Gaussian 积分，则

$$\int_{\mathbb R^{d_s}}e^{-q(s,w)}ds=(2\pi)^{d_s/2}\det(H_{ss})^{-1/2}e^{-q(s_*(w),w)}.$$

该式假设无限强变量域、正定强块和准确二次形式。非线性模型中 $H_{ss}$ 可随 $w$ 变化，边界、多个分支和非二次余项也会改变积分。只输出 corrector 最优点是 profile，不是 marginalization；把 determinant 当常数或省掉分支会改变概率。对第07章“可能性云”的候选构造，profile 可以是合法优化手段，但不因此获得 posterior 权重。

## F-OPT-009 · 与材料兼容的电流投影

在合法材料 $x$，所有可实现的一阶 current 变化来自

$$\delta h_s=L_D^{-1}\operatorname{diag}(E_{:,s})M\,\delta x.$$

按实部/虚部堆叠得到 $\mathcal T\in\mathbb R^{2qn_s\times66}$，即 $\delta h_{real}=\mathcal T\delta x$。任意 current proposal $u$ 若要一阶材料兼容，应投到 $\operatorname{range}\mathcal T$，而不是直接接受接收暗空间中的任意向量。Euclidean 正交投影为 $\mathcal T\mathcal T^\dagger$；若 current metric 为 $W_c\succ0$，对应

$$P_{\mathcal T}^{W_c}=\mathcal T(\mathcal T^TW_c\mathcal T)^\dagger\mathcal T^TW_c.$$

该式给出加权投影，但实现宜对白化后的 $W_c^{1/2}\mathcal T$ 做 QR/SVD。它一般不是普通 Euclidean 对称投影。

要同时一阶材料合法且数据不变，应先解 $\mathcal J\delta x=0$，再映射 $\mathcal T\delta x$；其空间是 $\operatorname{range}\mathcal T\cap\ker A_{obs}$，其中 $\mathcal J=A_{obs}\mathcal T$。材料投影与接收暗投影通常不对易，因此两次串联不等于该交集。若 $\mathcal J$ 在实 latent 空间满列秩，就没有非零 latent 一阶精确数据零方向；小奇异值只给近似弱方向。

构造完整 $\mathcal T$ 本身可能比已需的512×66观测 Jacobian贵得多。因此该式主要用于合法性解释；不建议为了“做一个投影”无条件构造巨大 current tangent。已知 decoder 的 latent proposal 加真实 corrector，通常更直接地维持原材料支持。

## F-OPT-010 · 嵌套缓存与错误的低秩推断

可按依赖顺序缓存：固定几何 $K,S,E^{inc}$ → 固定宿主 $C,E_R,O,F_h$ → 当前 $D,L_D$ factorization → 当前 $h,E,T$ → decoder 导数与 GN 小矩阵。该层次明确指出缓存何时失效。

材料增量若只改变 $k$ 个格点，$\Delta D\,C$ 的秩至多 $k$，可以考虑对应 Woodbury 更新；但 latent维数66不意味着 $\operatorname{rank}(\Delta D C)\le66$。一个 latent方向可能同时改变全部256个格点，使对角 $\Delta D$ 满秩。低维参数流形与每次矩阵更新低秩是不同性质。

完整 current 的重建 $J_c=J_h+V_hh$ 只在物理诊断或输出确实需要时做；接收预测与 ROI material Jacobian 不必每次重建576格场。若合同要求某项完整诊断，其成本仍必须保留。各源、各receiver和各 latent 参数应共享现成 factorization；重复 setup、重复求源和遗漏回退均会扭曲成本。

## F-OPT-011 · 由数学直接裁决的优先级

|候选|代数能确认的事实|仍需证据的事实|当前建议|
|---|---|---|---|
|完整固定ROI消元|原模型完全等价，全部非均匀支持保留|实际全流程费用|作为各臂公共强基线|
|8源共享因子、32receiver端口导数|完整66列导数无需528 tangent RHS|收缩和I/O的实际比例|直接采用正确结构|
|完整稳定基 QR corrector|保留原GN目标，避免正规方程条件数平方|是否改进具体求解器稳定性|普通稳定实现，不作为新理论|
|均匀幅度固定Schur|精确三角复用|特定任务是否值得复用|只用于真正均匀的条件模型|
|固定receiver SVD/TSOM表示|固定投影可组织current；截断有明确改变|整体速度与合法材料覆盖|不把表示优势冒称完整物理降维|
|current variable projection|内层可消元；有限罚目标不同|是否帮助难优化与总费用|低优先级，先与精确ROI强基线比较|
|低秩C、端口Krylov、网络初值|可接真实residual与fallback|误差界是否有用、是否摊销|见第05章，不先开rank sweep|
|冻结少量weak坐标|定义一种受限搜索；可能有可行性障碍|候选覆盖和收敛率|不宣称减少完整Jacobian成本|

## F-OPT-012 · 只保留能改变论文结论的实验

以下物理核端到端对照和额外严格证书实验是设计建议，不是新授权。可能性云E5则已另行批准并冻结，当前实施以 [最终E5协议](../config/possibility_volume_protocol.json) 与 [最终统计协议](../config/possibility_statistics_frozen.json) 为准；不能把它仍称为DRAFT或待批准。

**优化器/物理核主张：** 同一原始材料支持、同一任务、同一候选规则与接受条件，对比完整576系统和精确256系统的端到端费用；计入固定 setup、所有实际状态、导数、fallback和必要诊断。若只是相同16个 DEV 材料状态计时，只能报告核层结果，不能扩张成完整重建或新对象泛化。

**可能性体积成像主张：** 应比较相同候选数量/费用下的 GT邻域覆盖、数据相容性和集合大小；其具体无后验前提设计见第07章。这是当前主方向，不被严格后验资格阻塞。

**严格概率或证书主张：** 仅在明确提出该额外主张时，才要求原prior全支持、数值参考资格或严格积分包络，并与使用同一精确solver的强普通积分对照。该要求不应自动套在有限合法候选云上。

不为标准投影身份、Schur恒等式或QR换基追加小型公式实验；已批准的32个parent来自明确的新E5合同，不是DEV结果自动扩展，也不授权进一步扩容；不把新增名字当作先例不存在。最终科学取舍由主线程依据目标与成本决定。

## 来源与资格

以上为标准背景消元、伴随求导、variable projection、正交投影、QR最小二乘、Schur补与条件Gaussian积分的显式推导。SOM/TSOM背景见第02章及 `PRIOR_ART_AND_CLAIMS.md`。与概率降维相关的已读原始作者论文为 [Cui等2014 LIS](https://arxiv.org/html/1403.4680v2) 和 [Zahm等 certified reduction](https://arxiv.org/html/1807.03712v4)；二者不支持把当前有限候选云或profile直接叫成准确posterior。

本章执行新增物理RHS为0，没有运行代数数值测试。近似算法只是推导，不因写出公式而成为已验证实现。
