# 02 · 固定传播算子、投影与合法材料响应

本章先从波动方程说明哪些量固定，再区分三件经常混淆的事：固定接收投影、材料允许的电流变化、固定宿主的精确消元。接收端看不见某个电流，并不意味着该电流能由共同材料产生；固定 Green 算子也不意味着完整非线性响应固定。以下是标准散射与线性代数在本合同中的推导，不提出新的物理定律。符号与尺寸采用 F-BASE-001—004。

## F-OP-001 · 连续 Helmholtz 与体积分方程

采用时间因子 $e^{-i\omega t}$。在均匀背景、标量介电对比、固定磁导率的模型中，令 $k_b>0$、$\chi(\mathbf x)=k(\mathbf x)^2/k_b^2-1$，有

$$[\Delta+k_b^2(1+\chi)]u=-s,\qquad (\Delta+k_b^2)u^{inc}=-s.$$

出射 Green 函数满足 $(\Delta+k_b^2)g_b=-\delta$。因此

$$u(\mathbf x)=u^{inc}(\mathbf x)+k_b^2\int_\Omega g_b(\mathbf x,\mathbf x')\chi(\mathbf x')u(\mathbf x')\,d\mathbf x'.$$

证明只需对右侧作用 $\Delta+k_b^2$：积分项变成 $-k_b^2\chi u$，移项即得原方程；出射条件选择相应唯一分支。定义对比电流密度 $j_c=\chi u$ 后，非线性仅通过构成关系进入；传播 $j_c\mapsto u-u^{inc}$ 是固定线性算子。

二维 scalar-TM 与三维 scalar Helmholtz 的核分别为

$$g_2(R)=\frac{i}{4}H_0^{(1)}(k_bR),\qquad g_3(R)=\frac{e^{ik_bR}}{4\pi R}.$$

三维**矢量** Maxwell 一般需要 dyadic Green 张量及极化分量，不能把这里的标量 $1/R$ 核直接当作完整矢量模型。吸收背景可以使用相应复波数，但下面的远场幅相估计写在实 $k_b>0$ 条件下；复波数还包含距离相关衰减。

**资格：** `CONDITIONAL_THEOREM`，标准标量散射表示；它没有替离散积分误差、边界模型或连续 Maxwell 精度提供证书。

## F-OP-002 · 二维与三维远场：振幅近似和相位近似分开

令接收点为 $R\widehat{\mathbf r}$、目标坐标满足 $|\mathbf x'|\le a$。几何展开为

$$|R\widehat{\mathbf r}-\mathbf x'|=R-\widehat{\mathbf r}\cdot\mathbf x'+\frac{|\mathbf x'|^2-(\widehat{\mathbf r}\cdot\mathbf x')^2}{2R}+O(a^3/R^2).$$

因此，在 $R\gg a$ 且 $k_ba^2/R\ll1$ 时，三维有

$$g_3(R\widehat{\mathbf r},\mathbf x')\simeq\underbrace{\frac{1}{4\pi R}}_{\text{振幅}}\underbrace{e^{ik_bR}e^{-ik_b\widehat{\mathbf r}\cdot\mathbf x'}}_{\text{相位}}.$$

二维还需 $k_bR\gg1$ 才能使用 Hankel 渐近式，得到

$$g_2(R\widehat{\mathbf r},\mathbf x')\simeq\underbrace{\frac{1}{\sqrt{8\pi k_bR}}}_{\text{振幅}}\underbrace{e^{i\pi/4}e^{ik_bR}e^{-ik_b\widehat{\mathbf r}\cdot\mathbf x'}}_{\text{相位}}.$$

振幅在目标内近似常数的相对误差规模是 $a/R$；丢弃二次路径差造成的相位误差规模是 $k_ba^2/R$。二维另有 Hankel 渐近误差规模 $(k_bR)^{-1}$。这些条件不同：距离足以让振幅变化小，仍可能不足以让相位线性化。近场二维核也不能简化为精确的 $R^{-1/2}$ 振幅。

若入射为平面波 $u^{inc}=a_s e^{ik_b\widehat{\mathbf s}\cdot\mathbf x'}$，Born 远场的空间积分为

$$\int\chi(\mathbf x')e^{-ik_b(\widehat{\mathbf r}-\widehat{\mathbf s})\cdot\mathbf x'}d\mathbf x'.$$

这解释了固定收发几何对应的 Fourier 采样；它同时使用弱散射、远场接收和平面波照明。有限距离点源需保留自身 Green 振幅与相位，不能自动套用单一 Fourier 频率。当前离散核保留原合同，不因上述解释而替换成远场核。

## F-OP-003 · 离散 Born 双边结构及其非线性边界

对体积分作固定单元离散，令 $m_i$ 为单元对比材料。把 $k_b^2$、体积权重和自单元积分处理纳入内部 $K$ 与接收 $S$，而不再重复乘入 $m$。以所有源为列：

$$E^{inc}\in\mathbb C^{N\times n_s},\quad K\in\mathbb C^{N\times N},\quad S\in\mathbb C^{n_r\times N},\quad D_m=\operatorname{diag}(m).$$

零对比背景下，第一次 Born 替换 $E$ 为 $E^{inc}$，所以

$$F_B(m)=S D_m E^{inc}.$$

若显式写源传播 $G_{source}$ 与源振幅矩阵 $A_s$，则 $E^{inc}=G_{source}A_s$，故 $F_B=G_{receiver}D_mG_{source}A_s$，其中 $G_{receiver}=S$。两侧固定、中央对角材料变化，正是固定双边结构。Born 的来源是

$$J_c=D_m(E^{inc}+KJ_c)=D_mE^{inc}+D_mKD_mE^{inc}+\cdots.$$

Neumann 展开需要相应收敛条件，例如 $\|D_mK\|<1$ 是充分条件；方程本身可在该条件不成立而 $I-D_mK$ 仍可逆时求解。截掉高阶项是物理近似，不是单纯换基。非零宿主不能直接套用零背景 Born；F-OP-008 给出宿主上的精确响应，其一阶项才是 distorted Born。

## F-OP-004 · 固定接收 SVD 到底能确定什么

先考虑同一固定、可逐源使用的白化 $W\in\mathbb C^{n_r\times n_r}$。当前合同是相同方差的独立复噪声，可取 $W=I/\sigma$；若噪声在源间非可分离相关，应对白化后的整体向量算子做 SVD，不能未经说明沿用一个共享接收矩阵。

令 $A=WS=U_r\Sigma_rV_r^*$ 为保留全部非零奇异值的薄 SVD，$r=\operatorname{rank}A$，$V_0$ 是其正交补。任意电流有唯一分解

$$J_c=V_r a+V_0b,\qquad a\in\mathbb C^{r\times n_s},\quad b\in\mathbb C^{(N-r)\times n_s}.$$

正交分解观测残差得到

$$\|AJ_c-WY\|_F^2=\|\Sigma_ra-U_r^*WY\|_F^2+\|(I-U_rU_r^*)WY\|_F^2.$$

第二项与材料无关，所以在**固定模型的归一化 posterior**中可以约掉相应常数；若比较模型 evidence，必须保留其值。第一项仍依赖 $a$，含噪推断不能把 $a$ 无条件固定为 $\Sigma_r^{-1}U_r^*WY$，否则把有限噪声似然改成了精确数据约束。

固定投影为 $P_o=V_rV_r^*$、$P_n=I-P_o$，且 $SP_n=0$。若只取前 $\ell<r$ 个奇异值，则所谓“模糊”补空间还含弱辐射方向；它不是严格接收零空间。观测截断误差见 F-ERR-007。

**固定项：** $W,S,U,\Sigma,V$ 在采集与噪声固定时固定。**变化项：** 合法材料实际激发的 $a,b$、总场、resolvent 与 latent Jacobian。即使 $S$ 固定，也不能据此冻结材料 Jacobian。

## F-OP-005 · SOM/TSOM 的投影不是两个空间的交集

作者综述 [Chen、Zhong、Agarwal 2010](https://doi.org/10.4310/MAA.2010.v17.n4.a6) §3.1 式28—34描述接收 SVD 电流分解，§3.3 式46—47及其后续讨论描述内部作用筛选与 Fourier 表示。这些是已有方法构件。以下用统一符号解释其几何含义，不把不同篇章的算法目标函数混为一式。

设 $P_w$ 为去除已选接收强方向的正交投影，$V_D$ 是内部传播所选的方向矩阵。一个典型二次筛选构造是

$$Z=P_wV_D.$$

其列空间是内部方向**投影到**接收模糊空间后的像 $\operatorname{range}Z$，一般不是 $\operatorname{range}(V_D)\cap\operatorname{range}(P_w)$。最小反例：在 $\mathbb R^2$ 中，$P_w$ 投到 $e_1$，而 $V_D=(e_1+e_2)/\sqrt2$。则 $\operatorname{range}(P_wV_D)=\operatorname{span}(e_1)$，两原空间交集却只有零。

投影后的列可以相关或消失，应对 $Z$ 作一次 rank-revealing QR/SVD。真正的正交投影是

$$P_Z=ZZ^\dagger=Z(Z^*Z)^\dagger Z^*,$$

而一般不是 $P_wV_DV_D^*P_w$。后一个矩阵虽为半正定，但通常不幂等；把它当投影会重复缩放方向。

若两个正交投影 $P,Q$ 对易，则 $PQ=QP$ 为交集的正交投影；反之，$PQ$ 若为正交投影，由 $(PQ)^*=QP$ 可知必须对易。一般 $PQ$ 是串联线性作用，不是交集投影。需要**真正交集**时可设

$$C_\cap=\begin{bmatrix}I-P\\I-Q\end{bmatrix},\qquad P_\cap=I-C_\cap^\dagger C_\cap.$$

因为 $\ker C_\cap=\operatorname{range}P\cap\operatorname{range}Q$，该式才给交集。上述都是标准正交几何，不应冠以新的子空间定理。

## F-OP-006 · 固定投影为何通常不与非均匀材料交换

对任意固定矩阵 $P_o$ 和对角材料 $D_m$，逐元素直接得到

$$[D_m,P_o]_{ij}=(m_i-m_j)(P_o)_{ij}.$$

标量材料 $D_m=vI$ 与任何投影对易；一般非均匀材料只有在投影不混合不同材料值的格点时才对易。接收 SVD 通常混合许多空间格点，因此不能把 $P_oD_m$ 静默换成 $D_mP_o$。

进一步，将精确 ROI 方程的电流按某个固定 $P_o,P_n$ 分块，观测子块独立于暗子块需要

$$P_oD C P_n=0,\qquad D=\operatorname{diag}(m_{ROI}-m_h).$$

仅有 $OP_n=0$ 不足以满足该条件。暗电流经内部传播和材料相乘后可能回到可观测方向，因而一般有 $OL_D^{-1}P_n\ne0$。F-ERR-001讨论 $[D,C]$ 对固定对角化的限制；这里强调的是接收不可见性并不在非线性构成方程下自动封闭。

## F-OP-007 · 多源共享材料是电流合法性的硬条件

给定任意候选 $J_c$，先由固定传播计算 $E=E^{inc}+KJ_c$。它能来自一幅共同标量材料，必须在每个格点满足

$$J_{c,is}=m_iE_{is}\quad\text{对所有源 }s.$$

若该行 $\|E_{i,:}\|^2>0$，无约束复数最小二乘的唯一最佳材料为

$$\widehat m_i=\frac{\sum_s\overline{E_{is}}J_{c,is}}{\sum_s|E_{is}|^2},\qquad R_{i,:}=J_{c,i,:}-\widehat m_iE_{i,:}.$$

共同复材料存在当且仅当所有这些行的 $R_{i,:}=0$，且在 $E_{i,:}=0$ 的行有 $J_{c,i,:}=0$。后一情形材料未被该状态识别；交叉小行列式 $E_{is}J_{it}-E_{it}J_{is}=0$ 在 $E=0,J\ne0$ 时也成立，因此只检验小行列式会漏掉非法电流。

无约束逐行合法还不够：实际材料必须处于原 decoder 像、满足宿主/ROI 支持及允许的实虚部范围。若要求实材料，最佳无约束系数改为分子实部除分母；若要求被动材料，在本时间约定下介电对比虚部非负。但逐格截断或投影不能替代耦合 decoder 约束，更不能自动产生原 latent prior 的样本。

源共享的 $m_i$ 是关键。为每个源各算 $m_{is}$，即使都能复现数据，也通常对应多幅不同材料。SOM 的电流分解须与状态/构成关系一起使用；其已有先例不能被简化成“从接收零空间自由采样就是合法材料”。

## F-OP-008 · 固定宿主与完整 ROI 的精确消元

令固定背景材料为 $m_h\in\mathbb C^N$、$D_h=\operatorname{diag}(m_h)$，ROI 选择嵌入为 $P\in\mathbb R^{N\times q}$，$P^TP=I_q$。允许原 ROI 内任意非均匀增量 $d\in\mathbb C^q$：

$$D_m=D_h+PDP^T,\qquad D=\operatorname{diag}(d).$$

定义固定背景 field/current 系统与背景解

$$A_h=I-KD_h,\quad L_h=I-D_hK,\quad E_h=A_h^{-1}E^{inc},\quad J_h=D_hE_h,\quad F_h=SJ_h.$$

假设背景系统可逆。预计算

$$B=A_h^{-1}KP\in\mathbb C^{N\times q},\quad V_h=L_h^{-1}P\in\mathbb C^{N\times q},$$
$$C=P^TB\in\mathbb C^{q\times q},\quad E_R=P^TE_h\in\mathbb C^{q\times n_s},\quad O=SV_h\in\mathbb C^{n_r\times q}.$$

两种记号完全相同，因为

$$L_h^{-1}=I+D_hA_h^{-1}K,\qquad A_h^{-1}K=KL_h^{-1},$$
$$V_h=P+D_hB,\qquad C=P^TKV_h,\qquad O=S(P+D_hB).$$

第一恒等式可由左乘 $L_h$ 并用 $(I-KD_h)A_h^{-1}=I$ 验证；第二来自 $A_hK=KL_h$。这同时说明不能丢掉观测端口中的宿主再散射 $SD_hB$。当前代码固定标量宿主 `host`，因此其 `Ah=I-host*G`、`B=Ah^-1 G[:,ROI]`、`C=B[ROI]`、`O=S[:,ROI]+host*S@B` 与上述一般式严格一致。

令 $h=DP^TE$ 为 ROI 新增对比产生的电流。原方程给出

$$E=E_h+Bh,\quad J_c=J_h+V_hh,\quad P^TE=E_R+Ch,$$

所以所有材料状态只需求解

$$\boxed{L_Dh=DE_R,\quad L_D=I_q-DC,\quad F=F_h+Oh.}$$

等价 field 形式是 $(I_q-CD)E_{ROI}=E_R$、$h=DE_{ROI}$。这两个矩阵顺序不同，不能混写。它们可逆性等价，且在背景可逆时与完整材料系统可逆性等价：由矩阵行列式引理，$\det L_m=\det L_h\det(I_q-DC)$。

推导没有除以 $d_i$，因此 $d_i=0$ 或部分材料恰等于宿主都合法。保留全 $q=256$ ROI 的消元适用于原 Full66 decoder 的非均匀材料；它不是均匀幅度假设，也没有削减材料支持、改变 prior 或忽略多重散射。

**成本：** 稠密实现一次背景 $N$ 阶 LU，$q+n_s=264$ 个完整 RHS，以及构造 $C,O$；之后每状态一个 $q$ 阶 LU、$n_s=8$ 个 reduced RHS 与观测乘法。固定存储可包括 $B$ 或 $V_h$，但只求接收预测时不必重建完整 $E,J_c$。这些固定项只在背景、频率、采集与 ROI 不变时复用。改变这些条件需要重新 setup。

**资格：** `EXACT_DISCRETE_IDENTITY`，标准背景消元。低秩截断 $C$、截断 $O$ 或减小 ROI 才引入新的近似；其校正/误差控制见第05章。矩阵浮点残差合格不等于连续 PDE 认证。

## F-OP-009 · 均匀幅度是精确 ROI 响应的特殊情形

当 $D=vI$ 时，一次 $C=QTQ^*$ 的复 Schur 给出

$$F(v)=F_h+(OQ)v(I-vT)^{-1}(Q^*E_R).$$

这保留全部 $q$ 维，参数变化只进入三角系统。$v=0$ 时自然返回宿主。它可用于 R6 的均匀复幅度条件模型；非均匀 Full66 中 $Q^*DQ$ 通常稠密，不能再声称每状态仍只解同样的三角系统。详见 F-ERR-001。

在满足收敛条件的邻域内，$h=DE_R+DCDE_R+\cdots$，所以宿主 distorted Born 为 $F_h+ODE_R$。这里的有效源 $E_R$ 与有效观测 $O$ 都已包含背景散射；把 $O$ 换回裸 $SP$ 会漏掉宿主反馈。

## F-OP-010 · 原非线性导数也沿固定端口传播

对 F-OP-008 的完整非均匀方程微分。因 $\delta L_D=-D_{\delta d}C$，

$$L_D\,\delta h=D_{\delta d}(E_R+Ch)=D_{\delta d}E_{ROI},$$
$$\boxed{\delta F=O L_D^{-1}D_{\delta d}E_{ROI}.}$$

令实 latent $x\in\mathbb R^{66}$，$M(x)=\partial d/\partial x\in\mathbb C^{q\times66}$；定义当前端口 $T(x)=OL_D^{-1}\in\mathbb C^{n_r\times q}$。逐源复导数为

$$\frac{\partial F_{:,s}}{\partial x}=T\operatorname{diag}(E_{ROI,:,s})M.$$

因此全部 latent 列只需一次当前 forward 场和一次端口计算；并不需要对66列逐列做8源 tangent solve。若按源堆叠复观测，材料导数写成 $\operatorname{vec}(\delta F)=(E_{ROI}^T\odot T)\delta d$，其中 $\odot$ 是逐列 Kronecker 积，尺寸为 $(n_rn_s)\times q$。实化与 $\sqrt2/\sigma$ 因子严格沿 F-BASE-002，最终 $\mathcal J\in\mathbb R^{512\times66}$。

求全部 $T$ 可解 $L_D^T X=O^T$，再取 $T=X^T$，共32个 reduced RHS。这里是**转置**求解；若程序用共轭伴随模式，应解 $L_D^*X=O^*$ 再取 $T=X^*$。把 $O^T$ 与共轭伴随 solve 混配会产生错误。物理互易的复对称性也不等于 Hermitian 性。

所以每状态 full Jacobian 的传播成本是一次 $q$ 阶 factorization、8 forward RHS、32 transpose/adjoint RHS，加 decoder 与张量收缩。固定某个四维弱框架本身不会省掉这40个 RHS，除非算法真正改变了需要的导数信息。对应优化判据见 F-OPT-002、006。

## 来源、边界与后续使用

连续 VIE、Born、SVD、电流构成关系与投影均为标准理论。可读作者综述定位见 F-OP-005；[TSOM原篇](https://doi.org/10.1088/0266-5611/25/8/085003)、[FFT-TSOM原篇](https://doi.org/10.1109/TAP.2010.2103027)在本轮核验身份，实质方法支持主要来自已读作者综述，不能冒称已逐式核对这两篇原文。全文与检索错误的证据等级见 `PRIOR_ART_AND_CLAIMS.md`。

精确 ROI 和导数公式由上文独立代数推导并对应 `src/stable_roi_operator.py`。本章没有新增物理调用、公式小测试或实验。其结论决定可复用的算子与必须保留的约束；实际速度、材料覆盖和概率资格分别需要相应证据，不能从等式本身推断。
