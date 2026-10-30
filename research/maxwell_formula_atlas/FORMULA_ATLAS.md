# Maxwell 成像公式总集

这是持续维护的推导主干：固定传播算子、完整非线性材料、投影式SOM/TSOM、合法候选体积，以及按需增强的概率任务积分。公式可被解析证明时不再安排机器精度小实验。先明确哪些量固定、哪些量变化，再决定需要哪一项论文级实验。

当前版本覆盖实际已读、已使用的P1—R12主链和直接相关的F/F3、A17-X、Z-ALG公式；没有声称穷尽未提供的全部历史推导。每章保留源文件、成立条件与资格。完整posterior认证不作为材料可能性体积与风险成像的前置门。

# 01 · 统一符号、依赖与公式资格

本公式集的用途是把“哪里固定、哪里变化、哪里值得计算”直接从方程中读出来。它不是一个把所有历史等式重新编号的综述。当前版本覆盖 P1—R9 主计算链以及直接相关的 F/F3、A17-X、SOM/TSOM、R5/R6 概率与风险公式；没有读到的历史原始推导不冒称已纳入。后续方法必须引用公式编号、注明改变哪一项假设，然后才决定是否需要实验。

## F-BASE-001 · 物理离散空间与矩阵尺寸

采用时间约定 $e^{-i\omega t}$，被动材料的介电虚部非负。原冻结二维 scalar-TM 离散合同为：

$$N=576,\quad q=256,\quad n_s=8,\quad n_r=32,\quad d_x=66.$$

设 $K\in\mathbb C^{N\times N}$ 是包含 $k_b^2$ 与单元积分权重的内部 Green 算子，$S\in\mathbb C^{n_r\times N}$ 为接收算子，$E^{inc}\in\mathbb C^{N\times n_s}$ 为入射场；$m\in\mathbb C^N$ 为实际对比材料，$D_m=\operatorname{diag}(m)$。$J_c\in\mathbb C^{N\times n_s}$ 表示诱导电流，$\mathcal J$ 表示实化 latent Jacobian，两者不混用。

$$E=E^{inc}+KJ_c,\qquad J_c=D_mE,\qquad F(m)=SJ_c.$$

由代入得到两种完全等价的求解形式：

$$A_mE=E^{inc},\quad A_m=I-KD_m;\qquad L_mJ_c=D_mE^{inc},\quad L_m=I-D_mK.$$

这里不可交换 $KD_m$ 与 $D_mK$。二者可逆性等价，因为 $\det(I-AB)=\det(I-BA)$；“field form”与“current form”的逆矩阵本身不相同。原代码中的 `G` 对应本章 $K$，避免与材料 decoder $\mathcal G$ 同名。

**资格：** 冻结离散模型恒等式，非连续 PDE 误差证书。来源：R4/P1 `core.py`、A30 frozen VIE、R9 `src/stable_roi_operator.py`。

## F-BASE-002 · 原始观测、实化和白化

每个 $F_{rs}$ 是一个复观测；$\eta_{rs}\sim\mathcal{CN}(0,\sigma^2)$，所以每个实部和虚部方差均为 $\sigma^2/2$。固定

$$\sigma=0.00119211798514519.$$

令 $\operatorname{pack}$ 按 source-major 顺序逐源排列实部32项、虚部32项。它是实线性等距映射：

$$\|\operatorname{pack}Z\|_2^2=\|Z\|_F^2.$$

物理实化 $f=\sqrt2\operatorname{pack}F$，数据 $d=\sqrt2\operatorname{pack}Y$，白化 residual

$$r(x)=\frac{f(x)-d}{\sigma},\quad \Phi(x)=\tfrac12\|r(x)\|^2=\frac{\|F(\mathcal G(x))-Y\|_F^2}{\sigma^2}.$$

含噪观测残差、两个无噪预测之间的距离与白化导数必须使用同一 packing。$\Delta\Phi\le2$ 是当前似然竞争判据；它不是 posterior credible mass。

## F-BASE-003 · 材料损失与空间度量

材料作为实 Hilbert 空间，内积 $\langle a,b\rangle_M=\operatorname{Re}(a^*b)/576$。于是

$$\|m_1-m_2\|_M=\frac{\|m_1-m_2\|_2}{24},\qquad R(\hat m)=\mathbb E\|\hat m(Y)-M\|_M^2.$$

ROI-only RMS 的分母为16；除非明确改写，R7/R9半径0.01/0.02使用全网格分母24。固定宿主之外没有材料差异，并不允许静默改分母。复数 MSE 同时包含实部与虚部平方，不能只比较实部图。

## F-BASE-004 · 哪些算子真正固定

在背景、频率、网格、源接收位置及归一化固定时，$K,S,E^{inc}$ 固定；ROI嵌入矩阵 $P\in\mathbb R^{576\times256}$ 固定，$P^TP=I$。材料 $m$、$D_m$、总场 $E$、电流 $J_c$、resolvent $L_m^{-1}$、decoder导数以及 $\mathcal J(x)$ 一般变化。

由此有三个不同的复用层：

1. **固定几何层：** $K,S,E^{inc}$、receiver SVD、Fourier/Wavelet基。
2. **固定背景消元层：** 给定宿主后预计算有效ROI传播与观测端口；不同材料仍保留完整非线性。
3. **变化状态层：** 当前材料的小系统、场与 adjoint；可用近似预条件，但不能无误差地冻结全部内容。

SOM current 子空间、材料表示基与当前 $\mathcal J$ 的弱方向分别属于不同空间。它们之间需要显式映射，维数相同也不构成等价。

## F-BASE-005 · 公式依赖图

```text
Helmholtz Green + 单元积分
  ├─ Born: 固定双边传播 × 材料（仅一阶近似）
  └─ Full VIE: 固定 K/S + 变化 D_m
       ├─ 固定宿主精确ROI消元 → 256维完整非线性
       │    ├─ 均匀幅度 vI → 固定Schur三角化
       │    └─ 非均匀 D → 低秩预条件 + 实际residual校正
       ├─ receiver SVD → current投影 → 多源共同材料约束
       └─ decoder支持 → 材料有效切向 → 普通非线性corrector
原prior + 上述完整likelihood
  ├─ 固定Gaussian坐标 → 条件积分（不是profile）
  ├─ 商空间/混合秩 → 原像质量/coarea/原子
  └─ 材料任务 → 条件矩、任务误差、风险与计算停止
```

## 阅读规则与维护方式

每个条目包含等式、推导、成立条件及其计算含义。纯代数正确性由证明审核；不为“再看到一次机器精度相等”追加物理调用。实现纠错应针对实际错误。只有成本排序、真实材料覆盖、统计风险和误差界是否足够紧等无法由代数决定的问题，才进入路径裁决或论文实验。

公式资格分为 `EXACT_DISCRETE_IDENTITY`、`CONDITIONAL_THEOREM`、`CONTROLLED_APPROXIMATION`、`EMPIRICAL_ONLY`。没有统一的“全部PASS”标签。标准恒等式可用于重要优化，但不因此成为本研究的新定理。


---

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


---

# 03 原材料先验、有效几何与实际校正器

本章只整理源码与推导，没有新 Maxwell 求解、数值 QA、训练或 fresh parents。2026-10-29 formula-first steering下，历史 E2/E3/REV1 的 fresh confirmation 保持 NOT_RUN；后来用户另行批准了独立 E5 possibility-volume 合同，见最终 config/possibility_volume_protocol.json 与 config/possibility_statistics_frozen.json。此处不从合同状态推断 E5 结果。稳定编号 F-MAT-001 起；函数、行号和历史身份见 SOURCE_FORMULA_INDEX.csv。原 R9 E0–E3 与后来单独冻结的 REV1 分开记录，后者不覆盖前者的负面证据。

物理上先分清三件事：潜变量变化可能完全不改变材料；非零材料变化可能只有很小接收响应；局部弱方向走到有限距离后也可能不再弱。因此 decoder gauge、材料灵敏度和有限校正器不是同一条零空间结论。

## 1. 坐标、decoder 与先验

### F-MAT-001 原 P1 维度与参考测度

\[
q=(z,u_R,u_I)\in\mathbb R^{64}\times[0,1]^2,\quad
z\sim N(0,I_{64}),\quad u_R,u_I\overset{\rm iid}{\sim}U[0,1].
\]

材料 \(\chi\in\mathbb C^{576}\)，24×24 网格，8个照明、每照明32个复接收通道；\(F(\chi)\in\mathbb C^{32\times8}\)。先验密度相对于 \(dz\,du_Rdu_I\)，不是材料 Lebesgue 体积。原 P1 不是1152个独立实材料像素。

### F-MAT-002 probit Gaussian66

用 \(\Phi_N,\varphi_N\) 表示标准正态 CDF/PDF，与似然代价 \(\Phi_y\) 区分：

\[
x=(z,\xi_R,\xi_I),\quad u_a=\Phi_N(\xi_a),\quad
\varphi_N(t)=(2\pi)^{-1/2}e^{-t^2/2},\quad
du_a=\varphi_N(\xi_a)d\xi_a.
\]

内点精确变换使 \(x\sim N(0,I_{66})\)，没有额外添加物性先验。原生 \(u_a=0,1\) 对应无限 probit 坐标。R4 仅初始化时将原生边界移入 \([10^{-4},1-10^{-4}]\)；这不是普遍的边界等价证明。

### F-MAT-003 坐标 MAP 与 GN Hessian

\[
U_q=\Phi_y+\tfrac12\|z\|^2,\qquad
U_x=\Phi_y+\tfrac12\|x\|^2.
\]

正确 Jacobian 描述同一后验测度，但两个参考测度下的密度众数通常不同。R6 native 阶段使用 \(U_q\)，probit 阶段使用 \(U_x\)；R7/R9 的实际校正选择用 \(U_x\)。probit-coordinate mode 不能称为 native MAP。

\[
\nabla^2U_x=J_x^TJ_x+I+\sum_{i=1}^{512}r_i\nabla^2r_i,\qquad
H_{\rm GN}=J_x^TJ_x+I.
\]

GN 矩阵不是自动等于精确后验 Hessian，数值驻点也不是 posterior calibration。

### F-MAT-004 原相关场

令 \(X_i=p_{i1}+1/2,\;Y_i=p_{i2}+1/2\)，63个非恒定模态为 \(0\le k,l\le7,\;(k,l)\ne(0,0)\)：

\[
a_{kl}=\{1+(\pi\ell)^2(k^2+l^2)\}^{-1},\quad \ell=.08,
\]
\[
b_{kl}(i)=c_kc_l\cos(\pi kX_i)\cos(\pi lY_i),\quad
c_0=1,\;c_j=\sqrt2\ (j>0),
\]
\[
A=[\,\sigma_0\mathbf1,\ \sigma_f[b_{kl}]\operatorname{diag}(a)/\|a\|\,]
\in\mathbb R^{576\times64},\quad f=\mu\mathbf1+Az.
\]

实际 \(\sigma_0=.65,\sigma_f=1,\mu=0\)。这是固定 DCT 相关先验，没有逐对象峰值归一化。

### F-MAT-005 原 C2 饱和 decoder

\[
t(s)=\operatorname{clip}((s+1)/2,0,1),\quad h(s)=10t^3-15t^4+6t^5,
\]
\[
m_i=1_{\|p_i\|\le.38},\quad a_i=m_i h(f_i/\tau),\quad \tau=.5,
\]
\[
\delta=(.8+u_R)+i(.03+.15u_I),\qquad
\boxed{\chi_i=\chi_h+\delta a_i},\quad \chi_h=.1+.04i.
\]

\(a=0\) 是 host，不是空气；\(a=1\) 是纯 inclusion。\(h\) 在 \(s=\pm1\) 为 C2，但三阶导数一般跳变，跨 seam 的统一 C3 展开没有由此成立。

### F-MAT-006 一至三阶饱和导数

在 \(|s|<1\)：

\[
h'=\frac{15}{16}(1-s^2)^2,\quad
h''=\frac{15}{4}(s^3-s),\quad
h'''=\frac{15}{4}(3s^2-1).
\]

plateau 内导数为零。代码三阶 jet 使用严格内区指示，其 seam 取值约定不建立三阶可微性。

### F-MAT-007 decoder Jacobian 与 probit 链式因子

\[
B=\operatorname{diag}(m h'(f/\tau)/\tau)A\in\mathbb R^{576\times64},
\]
\[
D_q\chi=[\,\delta B,\ \Delta_Ra,\ i\Delta_Ia\,]\in\mathbb C^{576\times66},
\quad \Delta_R=1,\ \Delta_I=.15,
\]
\[
D_x\chi=D_q\chi\operatorname{diag}(1,\ldots,1,\varphi_N(\xi_R),\varphi_N(\xi_I)).
\]

链式因子只乘最后两列。有限实数下精确 \(\varphi_N>0\)；尾部浮点下溢或有效秩截断不能作为精确新 gauge 的证明。

### F-MAT-008 decoder jet 的无阶乘系数

代码采用 \(x(h)=x_0+hv+h^2x_2+h^3x_3\)，而非直接以二、三阶导数作系数。令 \(s_j=A(x_j)_z/\tau,\;x_1=v\)：

\[
a_1=mh's_1,\quad a_2=m(h's_2+\tfrac12h''s_1^2),
\]
\[
a_3=m(h's_3+h''s_1s_2+\tfrac16h'''s_1^3).
\]

任一尾物性坐标 \(\xi\) 的区间宽度为 \(\Delta\)，\(p=\Delta\varphi_N(\xi)\)：

\[
d_1=pv_\xi,\quad d_2=p((x_2)_\xi-\tfrac12\xi v_\xi^2),
\]
\[
d_3=p((x_3)_\xi-\xi v_\xi(x_2)_\xi+(\xi^2-1)v_\xi^3/6).
\]

组成复数 \(\delta_j\) 后，

\[
\chi_1=\delta_0a_1+\delta_1a_0,\quad
\chi_2=\delta_0a_2+\delta_1a_1+\delta_2a_0,
\]
\[
\chi_3=\delta_0a_3+\delta_1a_2+\delta_2a_1+\delta_3a_0.
\]

这部分 decoder 曲率与物性非线性不能全部归为 Maxwell 多次散射。

## 2. 实化、目标和有限端点条件

### F-MAT-009 两种 packing

\[
\mathcal R_m(v)=\frac1{24}\binom{\Re v}{\Im v}\in\mathbb R^{1152},
\qquad\|\mathcal R_m(v)\|=\|v\|/24.
\]

材料实化把全部实部与全部虚部分块。接收数据 pack 则是 source-major：每个照明先放32个实部、再放32个虚部，然后进入下一照明。因此前256个白化实分量是前4个照明，不是所有照明的实部块。

### F-MAT-010 白化与 likelihood

\[
\mathcal W(Y)=\frac{\sqrt2}{\sigma}\operatorname{pack}(Y),\quad
\sigma=.00119211798514519,
\]
\[
r(x)=\mathcal W(F(\chi(x))-y)\in\mathbb R^{512},\quad
\Phi_y(x)=\tfrac12\|r(x)\|^2=\|F(\chi(x))-y\|_F^2/\sigma^2,
\]
\[
U_x(x)=\Phi_y(x)+\tfrac12\|x\|^2.
\]

sigma 是固定绝对仪器噪声，不按每个 parent 信号幅度重缩放。原 P1 只从历史 F3 acquisition 复制 sigma，未复制 P0 decoder 或固定每对象 SNR。

### F-MAT-011 原 R7/R9 有限门槛

令 \(\Delta m=\mathcal R_m(\chi(x)-\chi(x_0))\)，\(u_0=M_0v_0/\|M_0v_0\|\)，目标半径 \(\rho>0\)：

\[
|\|\Delta m\|-\rho|\le.02\rho,\qquad u_0^T\Delta m\ge.5\rho,
\]
\[
\Delta\Phi_y\le2,\quad \Phi_y(x)\le\tfrac12\chi^2_{512,.95},\quad
\Delta P=\tfrac12(\|x\|^2-\|x_0\|^2)\le10,
\]

加原生 prior.legal(q)。R9 qualified 是这些条件的合取；DeltaPrior≤10 是 heuristic，不是先验硬支撑的定理。

pair 向量 \(d=\mathcal W(F(x)-F(x_0))\) 满足

\[
\boxed{\Delta\Phi_y=r_0^Td+\tfrac12\|d\|^2}.
\]

故 pair≤1、DeltaPhi≤2、绝对 noise ball 不可互相替代。

## 3. 有效材料秩、gauge 和谱

### F-MAT-012 有效材料截断

\[
M(x)=\mathcal R_m(D_x\chi)\in\mathbb R^{1152\times66}
=U\Sigma V^T,
\]
\[
\epsilon_M=\max(10^{-8}\sigma_1(M),10^{-10}),\qquad
r_M=\#\{\sigma_j(M)>\epsilon_M\}.
\]

这是 R7 material_jets.geometry 和 R9 的实际阈值。较早 decoder_geometry.rank 的默认相对阈值是 \(10^{-10}\)，F 为 \(10^{-7}\)，不能混作同一个数学秩。被截断的小非零奇异值不等于精确 gauge。

### F-MAT-013 plateau 同材料 gauge

设 \(\mathcal T=\{i:m_i=1,|f_i|<\tau\}\)。固定物性且 \(A_{\mathcal T}\dot z=0\) 时 transition 场精确不变。只要其余 ROI 点不跨 plateau seam，材料也精确不变。一个保守对称区间为

\[
|t|\le\min_{i\notin\mathcal T,\ m_i=1,\ (A\dot z)_i\ne0}
\frac{|f_i|-\tau}{|(A\dot z)_i|}.
\]

margin 为0时无正区间保证。这是同材料坐标冗余，不能计入非零 material RMS 的散射多解。

### F-MAT-014 物性径向补偿与 rank

对 \(\delta\ne0\)、有限 probit 内点、\(a\ne0\)，

\[
\operatorname{rank}_{\mathbb R}D_x\chi
=\operatorname{rank}B+1+1_{\{a\notin\operatorname{range}B\}}.
\]

证明：\(\delta B\) 贡献 \(r_B\) 维；两物性列形成 \(\{ca:c\in\mathbb C\}\) 的2维实空间；与前者在 \(a\in\operatorname{range}B\) 时交1维，否则交0维。全 host 的 \(a=0\) 应直接用原矩阵，不能套非零公式。

若试图 \(\dot\delta=c\delta\)，一阶抵消要求 \(B\dot z=-ca\)。不在 range 时只能给残差；即便一阶可消，也没有普遍的有限 gauge，因为 quintic 不满足一般齐次律 \(\alpha h(f/\alpha)=h(f)\)。

### F-MAT-015 材料单位 lift 与 prior 度量

\[
T=V_r\Sigma_r^{-1}\in\mathbb R^{66\times r},\quad
Q_m=MT=U_r,\quad Q_m^TQ_m=I_r,\quad
K=J_xT\in\mathbb R^{512\times r}.
\]

\(T\) 一般不正交，其坐标中的 Gaussian prior 精度是 \(T^TT\)，不是 I。若 r<66，仅在 \(\operatorname{range}V_r\) 采样还删掉原坐标支撑；材料截断几何不是完整后验替代。

### F-MAT-016 两种广义弱谱

纯响应材料单位谱为 K 的奇异值：

\[
\sigma_{\min}^2(K)
=\min_{\substack{v\in\operatorname{range}V_r\\v^TM^TMv=1}}
v^TJ_x^TJ_xv.
\]

R7/R9 实际初始方向却取

\[
A_M=K^TK+T^TT,\quad A_Mw_j=\lambda_jw_j,\quad v_0=Tw_1,
\]
\[
\lambda_1=\min_{\substack{v\in\operatorname{range}V_r\\v^TM^TMv=1}}
v^T(J_x^TJ_x+I)v.
\]

prior-aware eigenvalues 不等于 generalized_joint_singular 的平方。这个选择最小化局部 GN 二次曲率，不包含 \(r^TJv\) 的当前梯度交叉项，也不保证有限 likelihood 最优。符号按最大绝对 native 分量为正固定，独立于 GT。

### F-MAT-017 Gaussian 正交 frame 与 anchor 平移

下列按有效材料秩 r≥4 写。代码实际取 k=min(4,r) 个最小本征向量；若 r<4，以下4与62应分别替换为 k 与66-k。rank=0时这些非零方向公式不适用。取 A_M 最小4个本征向量 W_w：

\[
Q_{m,w}=MTW_w\in\mathbb R^{1152\times4},\quad
Q_w=\operatorname{qr}(TW_w)\in\mathbb R^{66\times4},
\]
\[
Q_s=\operatorname{null}(Q_w^T)\in\mathbb R^{66\times62},\qquad
Q=[Q_w,Q_s]\in O(66).
\]

材料 frame 与 Gaussian frame 是不同对象。固定完整正交 Q 时，\(x=Qz'\) 保留 \(z'\sim N(0,I)\)。但 anchor-centred \(x=x_0+Q\eta\) 对应

\[
\eta\sim N(-Q^Tx_0,I),\quad
\tfrac12\|x_0+Q\eta\|^2=\tfrac12\|Q^Tx_0+\eta\|^2.
\]

不能把平移后的 eta 又设成中心标准正态。随状态刷新的 Q(x) 是移动校正 frame，不自动构成单位 Jacobian 的全局密度变换。

### F-MAT-018 source-group 谱和共同补偿

R7 对前4照明 \(K_A=K_{1:256,:}\)、后4照明 \(K_B=K_{257:512,:}\)、以及 joint K 分别 SVD；各组最弱4个右奇异向量通过 Q_m 提升为材料 frame，再算主角余弦。

\[
\|Kv\|^2=\|K_Av\|^2+\|K_Bv\|^2.
\]

不同观测块不存在用相反符号相互抵消 residual 的自由。各组必须共用材料方向与补偿系数。对共同补偿 \(z_*=-A_c^+b\)，若各组单独最优为 z_g，

\[
\min_z\sum_g\|A_gz+b_g\|^2
=\sum_g\|A_gz_g+b_g\|^2+
\sum_g(z_g-z_*)^TA_g^TA_g(z_g-z_*).
\]

第二项是共享材料的组间冲突；为不同 source 单独设置材料不是同一问题。

### F-MAT-019 projector、gap 与资格

对同秩正交 frame \(Q_a,Q_b\)，令 \(c_j=\sigma_j(Q_a^TQ_b)=\cos\theta_j\)：

\[
\|Q_aQ_a^T-Q_bQ_b^T\|_2
=\sin\theta_{\max}=\sqrt{1-c_{\min}^2},
\]
\[
\|Q_aQ_a^T-Q_bQ_b^T\|_F^2=2\sum_j(1-c_j^2).
\]

R9 最终 projector_drift 用第一式；不需要组装1152阶 projector。不同 rank 时不能静默套同秩公式。

简并谱簇内的单个向量可任意旋转，应比较整个簇。按 prior-aware 本征值选弱4维时，边界 gap 是 \(\lambda_5(A_M)-\lambda_4(A_M)\)，不能误用 K 的响应奇异 gap；按纯响应 K 选簇时才用对应奇异值分离。稳定性推论还需要同一度量、固定截断和算子扰动界。当前程序未计算 outward gap 证书或沿有限路径的统一分离界，因此没有已认证的稳定子空间保证。

作为未在当前程序认证的条件化理论：若已经把两个 Hermitian 矩阵放到同一个正交坐标空间，\(H_1=H_0+\mathcal E\)，目标为最小k个本征值的簇，边界 gap \(\delta=\lambda_{k+1}(H_0)-\lambda_k(H_0)>0\)，且 \(\varepsilon=\|\mathcal E\|_2<\delta/2\)，标准谱投影扰动界给出一个保守形式

\[
\|P_1-P_0\|_2\le 2\varepsilon/\delta.
\]

本式只是条件清单，不是已执行数值结论。尤其当前 \(A_M=T^T(J^TJ+I)T\) 的T也随材料变化：首先必须固定rank并对齐r维frame；即使本征向量W的子空间稳定，Gaussian lift \(\operatorname{range}(TW)\) 还需要T的变化界和lift的非退化条件。材料frame \(\operatorname{range}(MTW)\) 又包含M的变化。因此“一个 A_M gap 很大”本身不能直接证明实际 Gaussian 或材料 weak projector 沿有限路径稳定。

## 4. 实际 GN、有限 RMS retraction 与 revision

### F-MAT-020 R7 native/Frozen-J 子问题

令 E_k 为 M_k 保留的右奇异向量：

\[
Z_k=E_k\operatorname{null}(u_0^TM_kE_k),\quad
a_*=-\binom{J_kZ_k}{Z_k}^{+}\binom{r_k}{x_k},\quad d_k=Z_ka_*.
\]

native A 用当前 J_k；Frozen-J 用 J_0，但依然计算真实非线性 F、当前 decoder M 和真实 prior。Frozen-J 不等于固定探索 frame，也不是只算一次 forward。

\[
d_k\leftarrow\gamma d_k,\quad
\gamma=\min\{1,\Delta_k/\|d_k\|,(\rho/2)/\|M_kd_k\|\}.
\]

这是最小二乘后统一缩步，不是精确求解两个球约束的 trust-region 子问题。

### F-MAT-021 R7 原 radial retraction

对 proposal p=x_k+d_k：

\[
x^+=x_0+h(p-x_0),\quad
\|\mathcal R_m(\chi(x^+)-\chi(x_0))\|=\rho,\quad
\|x^+-x_0\|\le4.
\]

代码从小正 h 起以1.35倍扩展，首个符号变化后用 Brent 求根，失败保留。decoder 距离未证明沿射线全局单调；此操作不是最近点投影、测地线或全局最优匹配。

### F-MAT-022 实际 retraction 后的 trust model

用实际增量 e=x^+-x_k 重查 \(\|e\|\le\Delta_k\)、\(\|M_ke\|\le\rho/2\) 和有限 signed progress；必要时原 d 减半、重做 retraction，最多9轮。

\[
\operatorname{pred}(e)=\tfrac12\left[
\left\|\binom{r_k}{x_k}\right\|^2-
\left\|\binom{r_k}{x_k}+\binom{J_ke}{e}\right\|^2\right].
\]

只有 pred>0 且步长不退化才支付下个真实状态。下次实际改善比为

\[
\varrho_k=\{U_x(x_k)-U_x(x^+)\}/\operatorname{pred}(e).
\]

小于.25则 trust 减半，大于.75则乘1.5并截到1。最多3个真实试点，保留最低 U_x 的可行点；不保证收敛 profile。拒绝/失败成本不能删掉。

### F-MAT-023 原 R9 E2/E3 固定弱系数

原 B 用固定 anchor Q_w，C 用当前刷新 Q_w,k。除了 signed row，又添加 radial tangent 行：

\[
n_k=\Delta m_k/\|\Delta m_k\|,\quad
Z_k=E_k\operatorname{null}
\begin{bmatrix}Q_{w,k}^TE_k\\u_0^TM_kE_k\\n_k^TM_kE_k\end{bmatrix}.
\]

再在 \(P_s=I-Q_{w,k}Q_{w,k}^T\) 中用 decoder Newton 恢复 RMS。令 \(g=\|\Delta m\|-\rho\)：

\[
a=P_sM^T\Delta m/\|\Delta m\|,\qquad
\delta x=-g\,a/(a^Ta).
\]

实际另有步长上限、回溯和30轮上限。这固定弱系数，不只是固定 basis。

**零旋转也可能阻断的解析反例：** M=I，\(x-x_0=\alpha e_1+s,\;s\perp e_1\)。若固定 \(\alpha=\rho\)，又要求 \(\|x-x_0\|=\rho\)，则 \(\rho^2+\|s\|^2=\rho^2\Rightarrow s=0\)。因此原停止不能直接归因于弱子空间旋转。这个反例也不自动证明每个真实 nonlinear 试点的具体失败原因。

### F-MAT-024 REV1 固定弱射线、幅度自由

REV1 删去 radial tangent 行，恢复原 R7 radial root：

\[
Z_k=E_k\operatorname{null}
\begin{bmatrix}Q_{w,k}^TE_k\\u_0^TM_kE_k\end{bmatrix},\quad
x^+=x_0+h(x_k+d_k-x_0).
\]

B 固定 \(Q_{w,k}=Q_{w,0}\)。若初始位移为 \(\alpha_0v_0\)，\(v_0\in\operatorname{range}Q_{w,0}\)，由 \(Q_{w,0}^Td_k=0\) 归纳：

\[
x_k-x_0=\alpha_kv_0+Q_{s,0}s_k,\quad
Q_{w,0}^T(x^+-x_0)=hQ_{w,0}^T(x_k-x_0).
\]

固定的是弱射线方向，幅度随 radial root 合法改变。固定 basis、固定 weak coefficient、固定 weak ray 是三件事。

这不是在 \((\alpha,s)\in\mathbb R^{1+62}\) 的63个变量上共同做 GN。代码先做 Qs 内受 signed-row 限制的 GN，再做一维 decoder root。满材料秩且约束独立时，GN 校正空间至多61维；原 E2多一行后至多60维。不能称完整63维 GN、精确 constrained profile 或其收敛解。

### F-MAT-025 动态 C 究竟保持什么

C 每次需要当前 J_k 时重建 Q_w,k；若 \(Q_{w,k}^TQ_{w,0}=U\Sigma V^T\)，作 \(Q_{w,k}\leftarrow Q_{w,k}UV^T\) 是 Procrustes 对齐，不改变子空间/补空间。

对当前一步 \(w_k=Q_{w,k}^T(x_k-x_0)\in\mathbb R^4\)：

\[
Q_{w,k}^T(x^+-x_0)=hw_k.
\]

w_k 一般有4个非零分量；代码没有把其中3个置零。下次换 frame 后也没有显式平行运输或投回 anchor ray。准确描述是当前弱补空间刷新，不是跟踪唯一动态弱方向或覆盖整个动态4维弱空间。

### F-MAT-026 有限弱射线缺陷

\[
w_0=Q_{w,0}^Tv_0,\quad w=Q_{w,0}^T(x-x_0),\quad
c=w^Tw_0/(w_0^Tw_0),
\]
\[
\epsilon_{\rm ray}=\|w-cw_0\|,\quad
\epsilon_{\rm rel}=\epsilon_{\rm ray}/\max(\|w\|,10^{-30}).
\]

它检查方向，允许幅度 c 改变。还需单独看 c 的符号、有限 signed progress 和 RMS，不能以一个投影残差替代原门槛。

## 5. 其他实际 arms 与 F/F3/A17X 接缝

### F-MAT-027 global DEV frame 与随机对照

\[
C_G=\sum_{j=1}^{4}Q_{w,j}Q_{w,j}^T,\qquad
Q_{w,G}=\operatorname{eigvec}_{\rm top4}(C_G).
\]

global frame 仅用原4个 development anchors；当前 anchor 只用于材料速度归一化。随机对照在 Gaussian native 坐标取 \(\eta\sim N(0,I)\)，投到指定 material strong 约束的 native nullspace 后材料单位归一化；不是在 whitened 材料坐标各向同性取样再通过病态 T 提升。

### F-MAT-028 Fourier/Haar 目标回到原 latent 支撑

16列目标 \(B_m\in\mathbb R^{1152\times16}\) 先乘 ROI mask 与 anchor inclusion increment：

\[
L=M^+B_m,\quad \widehat B_m=ML,\quad
\epsilon_B=\|\widehat B_m-B_m\|_F/\|B_m\|_F.
\]

对 \(\widehat B_m=U_B\Sigma_BV_B^T\) 截断，\(T_B=LV_{B,r}\Sigma_{B,r}^{-1}\)，再对 \(JT_B\) 和 \(T_B^TT_B\) 的 prior-aware 二次型选方向。端点仍由原 decoder 生成。“Wavelet”实际为4×4 Haar scaling/piecewise-constant目标，不是完整多尺度变换。受限子空间与投影误差意味着它们不是完整66维问题的无损换基。

### F-MAT-029 SOM dark-current 的共享材料 lift

若 \(S=U\Sigma V^*\)，所选候选为

\[
\dot j=(I-VV^*)j_0,\quad
b=(I-\operatorname{diag}(\chi_0)G)\dot j.
\]

所有照明共用 \(\dot\chi\) 应满足 \(b_{i,p}=\dot\chi_iE_{i,p}\)。逐格 LS 候选：

\[
\dot\chi_i=
\frac{\sum_p\overline{E_{i,p}}b_{i,p}}{\sum_p|E_{i,p}|^2}.
\]

检查共享材料残差及 \(\mathcal R_m(\dot\chi)\) 对 range(M) 的 lift 残差，两个相对阈值均为 \(10^{-6}\)。R9 分母用 \(10^{-30}\) 下限，而 A17X row_fit 对零场显式 undetermined；epsilon-clipping 不证明数学唯一性。MATERIAL_LIFT_FAIL 只否定所选候选，不否定全部 SOM null directions。

### F-MAT-030 旧 F/F3 P0 quotient 与当前 P1 不同

F 的8个均匀参数对应 max-normalized Gaussian 固定 family；其 material frame 采用 \(10^{-7}\) relative cutoff 和纯响应最弱奇异方向。F3 积分掉共同正权重缩放 gauge 得7维 quotient。原权重独立均匀于 [.4,.9]，a=log(9/4)：

\[
\rho_w=\log(w_0/w_1),\quad t=(\rho_w+a)/(2a),
\]
\[
p_{\rho_w}(\rho_w)=
\{.81e^{-|\rho_w|}-.16e^{|\rho_w|}\}/.5,\quad
p_t(t)=2a\,p_{\rho_w}(2at-a).
\]

这不是给7维 chart 均匀 prior。边界材料合法而连续先验点密度为0。固定最大格 j 的局部层内，F3 max-normalization 导数为

\[
D_q\chi=1.2V\{(D_qz)/z_j-z(D_qz)_j/z_j^2\}.
\]

最大格改变需分层；该 P0 gauge 不能搬给不做峰值归一化的当前 P1。7维 latent 也不等于物理 \(\Delta X=\operatorname{diag}(\Delta\chi)\) 的 rank=7。

### F-MAT-031 A17X 的物理接缝与边界

\[
L=I-\operatorname{diag}(\chi)G,\quad
\dot j=L^{-1}\operatorname{diag}(\dot\chi)E,\quad
\dot F=S\dot j.
\]

必须再通过 \(D_x\chi\) 拉回原 prior，才能得到 R7/R9 的 \(J_x\in\mathbb R^{512\times66}\)。ambient current nullspace 不自动给合法材料自由度。

声明 \(\Delta X=UDU^*\) 的端口模型用 \(Q_H(D)=(I-DH)^{-1}D\)。A17X current-space exchange 为

\[
R_V=V(V^*LV)^{-1}V^*,\quad
\bar C=(I-R_VL)C,\quad \Gamma_C=C^*L\bar C,
\]
\[
\Delta P=\bar C\Gamma_C^{-1}C^*(U-LR_VU).
\]

这些需要 core/Schur 可逆条件，是 current-space 公式，不是材料 decoder、prior-preserving换基或 finite-RMS retraction。R9 B/C 没有使用该 exchange 核，不能称已验证 A17 exchange。

### F-MAT-032 条件化局部无精确 fiber 界

历史 F3/A17X 兼容性文件给出：在合法固定度量 chart，若 \(f=\mathcal W\circ F\circ\chi\) 为 C1，\(\sigma_{\min}(Df(x_0))=\beta>0\)，段上 Df 为 L-Lipschitz，则

\[
\|f(x_0+h)-f(x_0)\|\ge\beta\|h\|-\tfrac12L\|h\|^2.
\]

这在 \(0<\|h\|<2\beta/L\) 排除该 chart 的非零 exact-isodata endpoint，不排除有限噪声宽度、远处模式或别的 chart/prior。R9 没有统一 outward beta/L 界；数值 full-column-rank SVD 不能升级为已认证应用。

## 6. 从公式得出的简化与未成立的声称

- 材料度量 T 与 Gaussian 正交 Q 分开；TᵀT 与 anchor 平移不能漏。
- 固定基不要求固定幅度；原 E2固定弱系数可能阻断校正，REV1固定弱射线是不同算法。
- 已执行的是 Qs GN + scalar radius root，未求解完整63变量 GN。
- 动态 C 只刷新 complement；没有唯一弱方向的显式运输或4维弱空间覆盖。
- pred≤0 或初始点停止会降低成本，不能自动叫速度优势；拒绝、失败、endpoint-J 和路径诊断均需保留成本。
- 历史 E2/E3/REV1 结果未达到20%成本优势门槛，其原 fresh confirmation 没有开展；后续 E5 是用户另行批准的可能性覆盖任务，不能回填原门槛或覆盖原失败。
- 真实 finite-path spectral gap、rank变化、seam、可微retraction和候选完整性仍未获得统一保证。五节点检查不是连续路径证明，候选集合不是 posterior。

本章源码核对与解析推导新增 Maxwell RHS=0；未读取大型原始数组。


---

# 04 概率、原材料先验与任务积分的公式集

本页是可持续增补的推导目录，编号 `F-PROB-xxx` 保持稳定。它整理原 P1/R5/R6 与 R8P/R9 实际使用的表达式、条件、源码位置和当前资格；不新增计算、物理实验或公式级微型测试。下文“精确”默认指给定离散矩阵及明确 prior 下的数学恒等式，不自动包含浮点舍入、Hankel生成、离散化或连续Maxwell误差。Full66始终为 `FULL66_UNQUALIFIED`。

记 \(\varphi,\Phi_N\) 为标准正态密度/CDF；用 \(\mathcal V\) 或 \(v\) 表示复材料幅度，避免与CDF混名。\(G_{\rm dec}\) 为材料decoder，\(K\) 为离散Green算子，\(F\) 为完整8源32接收forward。所有任务积分都针对同一个固定观测 \(y\)，除非另立模型合同。

## F-PROB-001｜原圆对称复噪声、白化与likelihood

**合同。** 256个复测量的噪声独立满足 \(\eta_k\sim\mathcal{CN}(0,\sigma_c^2)\)，即 \(E|\eta_k|^2=\sigma_c^2\)，实、虚部方差各为 \(\sigma_c^2/2\)。原值为 \(\sigma_c=0.00119211798514519\)。令 \(\operatorname{pack}\) 逐源排列Re32、Im32，则

\[
f(z)=\frac{\sqrt2}{\sigma_c}\operatorname{pack}F(G_{\rm dec}(z)),\quad y_w=\frac{\sqrt2}{\sigma_c}\operatorname{pack}y,
\quad \mathcal P(z)=\frac12\|f(z)-y_w\|^2=\frac{\|F(G_{\rm dec}(z))-y\|_F^2}{\sigma_c^2},
\quad L(z)=e^{-\mathcal P(z)}.
\]

**推导。** 每个复密度为 \((\pi\sigma_c^2)^{-1}e^{-|\eta|^2/\sigma_c^2}\)。乘256个测量，负log的材料相关项恰为上式。白化后的512实分量协方差为单位阵。删去常数 \((\pi\sigma_c^2)^{-256}\) 不改变同模型posterior；跨噪声/采集/model evidence比较必须恢复常数。Z-ALG“每实/虚σ=.001”是另一个合同，不能代入本式。

**源码。** [R4 `Target.phi/white`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport/src/p3r4/core.py:236)；[R6 `Target.white/phi`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/shared/physics.py:143)；[R9 `pack`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_transfer.py:44)。

## F-PROB-002｜原Gaussian64×Uniform2与Gaussianized66

原native坐标 \(q=(z_{0:63},u_R,u_I)\)，\(z_{0:63}\sim N(0,I_{64})\)，\(u_R,u_I\stackrel{\rm iid}{\sim}U(0,1)\)。令 \(x_{64}=\Phi_N^{-1}(u_R),x_{65}=\Phi_N^{-1}(u_I)\)，则 **完整原先验** \(x\sim N(0,I_{66})\)：因为 \(du_j=\varphi(x_j)dx_j\)，两项Jacobian正好补成两个Gaussian密度。

原decoder为

\[
s_i=\frac{b_i+A_i x_{0:63}}{\theta},\quad Q_i=H(s_i)1_{i\in ROI},\quad
D(x)=.8+\Phi_N(x_{64})+i\{.03+.15\Phi_N(x_{65})\},\quad
m_i=.1+.04i+D(x)Q_i.
\]

这里 \(\theta=.5\)，offset列 \(A_{i0}=.65\)。\(H\) 是原C² quintic saturation，不能换成sigmoid/clip后仍声称原prior相同。native Uniform边界对应无穷probit坐标；仅初始化时的向内clip不是精确概率变换。

**源码。** [R4 `to_probit/from_probit/probit_derivative`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport/src/p3r4/core.py:294)；[原 `MaterialPrior.decode/derivative`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport/src/p1_material/prior.py:50)；[`saturation`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport/src/p1_material/prior.py:10)。

## F-PROB-003｜材料推前、evidence与同材料条件核

\[
\mu_0=(G_{\rm dec})_\#\pi_0,\quad
Z_y=\int L_y(m)\mu_0(dm)>0,\quad
\pi_M(dm\mid y)=Z_y^{-1}L_y(m)\mu_0(dm).
\]

对任何有界可测材料任务 \(\tau\)，

\[
E[\tau(M)\mid y]=\frac{\int\tau(G_{\rm dec}(x))L_y(G_{\rm dec}(x))\pi_0(dx)}{\int L_y(G_{\rm dec}(x))\pi_0(dx)}.
\]

若 \(\kappa_0(dx\mid m)\) 是原prior对decoder的disintegration，由 \(L_y\) 仅依赖 \(m\) 可得 \(\pi(dx\mid m,y)=\kappa_0(dx\mid m)\)：数据不会在**精确同材料**gauge内重新赋权。相反，不同材料即使同forward，仍按各自材料prior质量分配；finite-noise近邻关系不能当作精确等价关系。

**源码/理论。** [R8P measure章节](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/02_measure.tex:5)；R9实现为[ `summary`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:12)，它计算有限求积的分子/分母，未实现Full66 disintegration。

## F-PROB-004｜平方损失的Bayes动作与固定数据风险

对实任务 \(T=\tau(M)\)，若 \(E[T^2\mid y]<\infty\)，令 \(\mu_y=E[T\mid y]\)，则

\[
R_y(a)=E[(a-T)^2\mid y]=\operatorname{Var}(T\mid y)+(a-\mu_y)^2,
\quad a^*=\mu_y.
\]

**推导。** 展开 \(a-T=(a-\mu_y)+(\mu_y-T)\)，条件交叉项因 \(E[T-\mu_y\mid y]=0\) 消失。复材料向量实化后，对半正定损失矩阵 \(W\)，同理

\[
R_y(a)=\operatorname{tr}(W\Sigma_y)+\|a-\mu_y\|_W^2.
\]

所以更精确计算动作只减少第二项，固定同一个 \(y\) 的posterior variance不因计算轮数下降。若均值区间半宽≤ε，中点动作的excess risk≤ε²。

**源码/理论。** [R8P square loss与停止](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/06_algorithm_cost.tex:35)；R6已实现的有限VAL恒等式是另一个对象，见F-PROB-021与[ `track_b/report.py`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/track_b/report.py:24)。

## F-PROB-005｜有限特征Bayes动作与信息损失

若 \(S=S(Y)\)，总体平方风险最优动作分别是 \(a_Y=E[T\mid Y]\)、\(a_S=E[T\mid S]\)。条件期望正交性给

\[
E\|T-a_S\|^2=E\|T-a_Y\|^2+E\|a_Y-a_S\|^2.
\]

进一步对任何只依赖S的学习动作 \(\hat a(S)\)，

\[
E\|T-\hat a(S)\|^2=E\|T-a_Y\|^2+E\|a_Y-a_S\|^2+E\|a_S-\hat a(S)\|^2.
\]

**条件与边界。** S必须是Y的可测函数、二阶矩有限。有限已暴露VAL上的teacher/student损失，不等于这里未知的population conditional mean。Gaussian四通道可作为任务有用特征，但特征有用不能反推其生成Gaussian是校准posterior。

**来源。** [R5任务理论说明](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/THEORY_AND_NUMERICAL_CHECKS.md:27)；这是解释学习目标的恒等式，R5/R6未计算未知 \(a_Y,a_S\) 的真实总体分解。

## F-PROB-006｜固定正交Full66坐标与anchor平移

设 \(Q=[Q_s,Q_w]\in O(66)\) **固定**。若 \(x=Q_s s+Q_w w\)，原 \(x\sim N(0,I)\) 等价于独立 \(s\sim N(0,I),w\sim N(0,I)\)。若为了局部几何改成

\[
x=x_*+Q_s s+Q_w w,
\]

则

\[
(s,w)=Q^T(x-x_*)\sim N(-Q^Tx_*,I),\qquad
\rho(s,w)\propto\exp\{-\tfrac12\|x_*+Q_s s+Q_w w\|^2\}.
\]

**推导。** 正交变换保留范数和Lebesgue体积；平移改变均值。忽略anchor-induced线性prior项会丢失prior tilt。固定 \(w=w_0\) 的slice是正则条件模型，不能代替对w的prior conditional积分，也不等于Full66材料边缘。Gaussian盒截断时需另保留截断归一化与domain；旋转后的原盒一般不再是独立轴盒。

**来源。** [R8P正交坐标原则](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/02_measure.tex:54)；[R5 anchor中心响应坐标约定](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/THEORY_AND_NUMERICAL_CHECKS.md:5)。固定正交与移动校正frame的区别亦见同目录03章；本式未被R9用来生成Full66 posterior。

## F-PROB-007｜非正交lift与非线性chart的prior

若满秩线性chart \(x=x_*+T\xi\) 且T为方阵可逆，则

\[
\rho_\xi(\xi)=\varphi_{66}(x_*+T\xi)|\det T|,
\quad \text{precision}=T^TT,\quad E\xi=-T^{-1}x_*.
\]

薄矩阵T只描述子空间slice；可与补空间组成完整变换后再条件化。材料lift \(T=V\Sigma^{-1}\) 通常不正交，不能赋予iid标准Gaussian系数。对可逆非线性 \(x=\Psi(s,w)\)，必须使用 \(\rho=\varphi_{66}(\Psi)|\det D\Psi|\)。状态依赖的移动frame \(Q(x)\) 不自动定义一个全域单射变换，不能逐点写det=1就略去全局映射。

**来源。** [R8P chart Jacobian](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/02_measure.tex:58)；实际实现资格只在明确固定chart内成立；一般全域非线性chart在R9 `NOT_RUN`。

## F-PROB-008｜条件积分不是profile：截面体积、Hessian与分支

固定弱坐标w，令含原prior/Jacobian的联合密度为ρ，则

\[
Z(w)=\int e^{-\mathcal P(s,w)}\rho(s,w)ds,
\quad N_\tau(w)=\int\tau(G_{\rm dec}(s,w))e^{-\mathcal P(s,w)}\rho(s,w)ds.
\]

只有对全部w/源层求和再取 \(N/Z\) 才是任务。\(\mathcal P_{\rm prof}(w)=\min_s\mathcal P(s,w)\) 缺少截面宽度、prior tilt和其他极小点。若s维数r，且每个孤立内部峰 \(s_k(w)\) 的likelihood Hessian \(H_k\succ0\)，Laplace主项是

\[
Z(w)\approx\sum_k e^{-\mathcal P(s_k,w)}\rho(s_k,w)(2\pi)^{r/2}|H_k(w)|^{-1/2}.
\]

若对 \(\mathcal P-\log\rho\) 做MAP展开，应改用该势能的Hessian并避免重复乘ρ。边界截断峰、seam、奇异Hessian需另处理；要成为误差证书还须控制局部余项、峰外evidence及遗漏峰，不能仅比较局部GH阶数。

**来源/实现。** [R8P条件积分](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/02_measure.tex:66)；[R9 fresh mode/GN](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:18)与[ `gh`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:27)使用真实likelihood校正，未把Laplace当oracle；[1D局部Laplace](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:45)只标局部近似。

## F-PROB-009｜响应坐标的体积因子

为避免把原噪声再除一次，使用已白化f。固定w及固定输出正交frameU，设 \(h=U^T(f(s,w)-y_w)\)。若 \(s\mapsto h\) 在第k分支上单射可微，且 \(B_k=U^TD_s f\) 可逆，则

\[
ds=|\det B_k|^{-1}dh,
\quad Z(w)=\sum_k\int e^{-\frac12\|h\|^2-\frac12\|P_{U^\perp}(f-y_w)\|^2}\,
\frac{\rho(s_k(h,w),w)}{|\det B_k|}dh.
\]

若使用未白化输出r及 \(h=U^T(r-y)/\sigma\)，则对应 \(\sigma^r/|\det(U^TD_s r)|\)。这里r为实噪声尺度所在合同，不能把σc与每实分量σ混用。多个逆分支必须相加；chart重叠需原域不交分割/partition of unity，不能重复同一prior质量。有限求根成功不是全域单射或覆盖证明。

**来源。** [R5响应坐标原式](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/THEORY_AND_NUMERICAL_CHECKS.md:5)。R9主线使用原source/材料幅度积分，没有新建非线性响应atlas。

## F-PROB-010｜材料quotient/coarea与混合秩参考测度

在一个d维原latent源层 \(S\) 上，若 \(G_{\rm dec}\) 为常秩r并满足coarea适用条件，设 \(J_rG\) 为非零奇异值之积，则相对于材料像的 \(r\) 维Hausdorff测度，

\[
\frac{d\mu_{0,S}}{d\mathcal H^r}(m)
=\int_{S\cap G^{-1}(m)}\frac{p_0(x)}{J_rG(x)}d\mathcal H^{d-r}(x).
\]

**推导。** 对coarea公式选择被积函数 \(p_0/J_rG\)，先沿原像纤维积分，再沿材料像积分。材料图的几何体积 \(\sqrt{\det(Dm^TDm)}\) 本身不是原材料prior；还需原latent密度与所有原像。

不同秩的源层分别推送并以**测度**相加，不能把不同维density作为同单位普通函数相加。最稳妥的通用表示为 \(\mu_0=\sum_jG_\#(1_{S_j}\pi_0)\)，原 \(S_j\) 不重叠。材料像重叠时各原prior贡献相加；同一个原prior区域被多个chart重复积分才是错误。

**来源。** [R8P coarea及混合层](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/02_measure.tex:62)、[同节源层公式](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/02_measure.tex:122)；[R6 `source_nodes/quotient_nodes`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/track_a/stratified.py:71)。

## F-PROB-011｜host原子、pure层与Full66的区别

设 \(H(s)=0\) 于 \(s\le-1\)、\(H(s)=1\) 于 \(s\ge1\)。原Full66中

\[
p_h=\Pr_0\{s_i\le-1\;\forall i\in ROI\},\quad
\mu_0\supset p_h\delta_{m_h},\quad
\Pr(M=m_h\mid y)=p_h L(m_h)/Z_y.
\]

这是一个高维Gaussian多面体概率，本轮没有计算。全部ROI纯相时若两个属性仍活动，材料像是二维连续层而不是单个原子；只有属性也固定的条件slice才把pure平台变为单材料原子。Full66不能套用uniform slice的单变量PH。finite-noise下若forward有限，则L严格正；浮点下溢0不能证明任一正prior平台的posterior质量为0。

**来源/实现。** [原 `MaterialPrior.decode`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport/src/p1_material/prior.py:50)；R6 uniform平台质量在[ `stratified.py`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/track_a/stratified.py:18)；R9下溢修复在[1D atom bounds](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:17)。

## F-PROB-012｜R6复幅度V=QD的正确联合prior与Q⁻² fibre

R6 S3明确采用三个活动Gaussianized坐标独立截断[-4,4]，其余场坐标0。记

\[
Z_0=\Phi_N(4)-\Phi_N(-4),\ A=1/1.3,\ U_- =\Phi_N(-4),\ U_+=\Phi_N(4),
\quad P_H=\frac{\Phi_N(-A)-\Phi_N(-4)}{Z_0},
\]

\[
Q(z)=H(1.3z),\quad R_D\in[R_-,R_+]=[.8+U_-,.8+U_+],\quad
I_D\in[I_-,I_+]=[.03+.15U_-,.03+.15U_+].
\]

由于 \(du=\varphi(x)dx\)，截断probit属性推成这些区间的**归一化Uniform**。在transition \(-A<z<A\)，联合prior为 \(\varphi(z)dz\,dR_DdI_D/(Z_0\Delta R\Delta I)\)。固定Q的变换 \((R_D,I_D)\mapsto(R,I)=Q(R_D,I_D)\) 的实Jacobian是 **Q²**。所以

\[
\rho_{\rm tr}(R,I)=\frac1{\Delta R\Delta I}
\int_{z_-}^{z_+}\frac{\varphi(z)}{Z_0Q(z)^2}dz,
\]

\[
q_- =\max(R/R_+,I/I_+),\quad q_+=\min(R/R_-,I/I_-,1),\quad
z_\pm=\frac{2I^{-1}_{q_\pm}(3,3)-1}{1.3}.
\]

当 \(0<q_-<q_+\) 时上式成立，否则transition密度0。pure层密度为 \(\rho_{\rm pure}=P_H1_{[R_-,R_+]\times[I_-,I_+]}/(\Delta R\Delta I)\)，host质量 \(P_H\delta_0\)。总材料prior为

\[
\mu_V=P_H\delta_0+(\rho_{\rm tr}+\rho_{\rm pure})\,dR\,dI.
\]

**源码。** [R6 `fibre`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/track_a/stratified.py:31)；[R9 `fibre/prior`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_transfer.py:50)。\(Q^{-1}\)只对应一维实乘积例，不能替代这里的 \(Q^{-2}\)。

## F-PROB-013｜V固定后的条件依赖、源层权重与copula

在transition条件于 \(V=R+iI\)，

\[
p(z\mid V,\mathrm{tr})\propto\varphi(z)Q(z)^{-2}1_{z_-<z<z_+},
\quad D=V/Q(z).
\]

因此Q与D不再独立。总source标签的条件概率为 \(\Pr(\mathrm{tr}\mid V)=\rho_{\rm tr}/(\rho_{\rm tr}+\rho_{\rm pure})\)，pure同理；这些同材料贡献不应被MERGE删掉。事件 \(Q>.5\) 的transition任务密度只把积分区间改为 \((\max(z_-,0),z_+)\)，再加pure密度。原source坐标直接积分时没有Q⁻²，因为尚未换成材料Lebesgue坐标。

历史实标量控制 \(Q,D\stackrel{\rm iid}{\sim}U[1,2]\)、条件QD=2，则

\[
p(q\mid QD=2)=\frac1{q\log2},\quad D=2/q,\quad
E(QD)=2,\quad \operatorname{Cov}(Q,D)=2-(1/\log2)^2.
\]

若独立生成正确的两个边缘，反而得到 \(E(QD)=(1/\log2)^2\ne2\)。这是已有条件依赖反例，不是本轮新测试。

**源码。** [R6 fibre phase numerator](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/track_a/stratified.py:43)；[R9 `source_integral`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:50)；[既有 `qd_reproduction`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/experiments/r8p_benchmarks.py:166)。

## F-PROB-014｜条件prior平均与真实conditional likelihood的差

固定s，以原conditional prior记 \(\bar L(s)=E_0[L\mid s]\)。直接代入归一化定义可得

\[
E[\tau\mid s,y]-E_0[\tau\mid s]
=\frac{\operatorname{Cov}_0(\tau,L\mid s)}{\bar L(s)}.
\]

若任务跨度为Δτ，且整个条件域有 \(0<\ell\le L\le u\)，则sharp bounded-tilt界给

\[
|E_L\tau-E_0\tau|\le\Delta\tau\frac{\sqrt u-\sqrt\ell}{\sqrt u+\sqrt\ell}
=\Delta\tau\tanh\{\tfrac14\log(u/\ell)\}.
\]

条件kernel不变需要**likelihood对该变量恒定**，不只是任务恒定。任务与弱变量无关时，弱方向likelihood平均仍可能改变强方向的边缘权重，不能固定成一个代表。

**来源。** [R8P条件均值与tilt界](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/02_measure.tex:98)。R9未以单点残差估计替代整个条件域 \(\ell,u\)。

## F-PROB-015｜Gaussian checkpoint返回原坐标后的均值/协方差

若网络坐标 \(\xi\sim N(m_\xi,\operatorname{diag}v_\xi)\)，历史row-vector约定为 \(x=(\xi\odot s+b)F\)，则

\[
m_x=(m_\xi\odot s+b)F,\quad C_x=F^T\operatorname{diag}(v_\xi\odot s^2)F.
\]

diagonal网络输出不意味着原材料坐标独立。对场 \(S_i=(b_i+A_ix_{0:63})/\theta\)，

\[
E S_i=(b_i+A_i m_z)/\theta,\quad
\operatorname{Cov}(S_i,S_j)=A_iC_{zz}A_j^T/\theta^2,\quad
\operatorname{Cov}(S_i,x_{64:65})=A_iC_{z,u}/\theta.
\]

只有在检查到相应交叉块为0后，Gaussian独立性才允许把属性/场期望相乘。这个Gaussian是近似q；其pushforward矩计算准确不等于q是真posterior。

**源码。** [R5 `raw_gaussian`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:73)；[ `diagonal_image`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:81)显式拒绝非零场–属性cross；[ `gaussian_image`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:92)。

## F-PROB-016｜Gaussian probit属性的一阶、二阶与相关矩

若 \((U,V)\) 联合Gaussian，均值μU、μV，方差vU、vV，covariance cUV，则

\[
E\Phi_N(U)=\Phi_N\!\left(\frac{\mu_U}{\sqrt{1+v_U}}\right),
\quad E\Phi_N(U)^2=\Phi_2(t_U,t_U;v_U/(1+v_U)),
\]

\[
E[\Phi_N(U)\Phi_N(V)]=\Phi_2\!\left(t_U,t_V;\frac{c_{UV}}{\sqrt{(1+v_U)(1+v_V)}}\right),\quad t_U=\mu_U/\sqrt{1+v_U}.
\]

**推导。** 引入与(U,V)独立的标准正态E1,E2，把CDF乘积写为事件 \(E_1\le U,E_2\le V\)，再标准化Gaussian差值。对 \(D=a+b\Phi_N(U)\)，\(ED=a+bE\Phi_N(U)\)，\(ED^2=a^2+2abE\Phi_N(U)+b^2E\Phi_N(U)^2\)。

**源码。** [R5 `bvn_same/bvn`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:17)；[ `correlated_cell`条件probit矩](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:67)。二元CDF数值积分本身仍需精度资格。

## F-PROB-017｜相关Gaussian通过quintic decoder的材料矩

对单像素联合Gaussian \((S,U_R,U_I)\)，记 \(v_S>0\)、\(c=\operatorname{Cov}(U,S)\)，则

\[
U\mid S=s\sim N\left(\mu_U+\frac c{v_S}(s-\mu_S),\ C_{UU}-\frac{cc^T}{v_S}\right).
\]

由F-PROB-016得到 \(d_a(s)=E[D_a\mid S=s]\)、\(d_{ab}(s)=E[D_aD_b\mid S=s]\)，并有

\[
\mu_{\Delta m,a}=\int H(s)d_a(s)\varphi_{\mu_S,v_S}(s)ds,
\quad E[\Delta m_a\Delta m_b]=\int H(s)^2d_{ab}(s)\varphi_{\mu_S,v_S}(s)ds.
\]

减去 \(\mu_a\mu_b\) 得实虚2×2 covariance。vS=0时PSD要求c=0，取S确定值极限即可。只有场–属性独立时才简化为 \(E[HD_a]=E[H]E[D_a]\)。

四通道只含Re/Im均值和各自variance；复Hermitian variance为VarR+VarI，pseudo-covariance为VarR−VarI+2iCov(R,I)，二者不同。跨像素相关还需 \((S_i,S_j,U)\) 的联合积分；逐像素四通道不是完整空间posterior sampler。

**源码。** [R5 `correlated_cell`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:58)；[ `gaussian_image`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:92)。

## F-PROB-018｜quintic的Gaussian截断矩与平台

\[
H(s)=\begin{cases}0&s\le-1,\\
\frac12+\frac{15}{16}s-\frac58s^3+\frac3{16}s^5&-1<s<1,\\1&s\ge1.\end{cases}
\]

对S=μ+σZ，设a=(-1−μ)/σ、b=(1−μ)/σ，\(I_k(a,b)=\int_a^b z^k\varphi(z)dz\)。分部积分给

\[
I_0=\Phi_N(b)-\Phi_N(a),\quad I_1=\varphi(a)-\varphi(b),\quad
I_k=a^{k-1}\varphi(a)-b^{k-1}\varphi(b)+(k-1)I_{k-2}.
\]

将 \(H(\mu+\sigma z)=\sum_{k=0}^5c_kz^k\) 展开，则 \(EH(S)=\Phi_N(-b)+\sum c_kI_k\)；H²用10阶矩，平台tail同样只加一次 \(\Phi_N(-b)\)。σ=0取H(μ)。极窄/极尾区间有相消，应使用稳定CDF差或seam-split积分，不能从符号精确推断浮点实现自动严格。

另有 \(H(-1+\epsilon)=\frac54\epsilon^3-\frac{15}{16}\epsilon^4+\frac3{16}\epsilon^5\)，说明decoder接缝的三次起始；它不证明Full66存在物理精确解纤维。

**源码。** [R5 `truncated_phase_moments/phase_moments`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:35)；[R9 1D Gaussian moments](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:48)。

## F-PROB-019｜材料阈值与256-ROI硬事件的prior矩

规定

\[
\tau_{\rm mean}=\frac1{256}\sum_{i\in ROI}\Re m_i,\quad
E=\{\#\{i\in ROI:\Re m_i\ge.75\}\ge64\}.
\]

固定合法正实属性 \(D_R\)，单像素阈值转成 \(H(s_i)\ge q_*=(.75-h_R)/D_R\)。若0<q*<1，\(s_i\ge 2I_{q_*}^{-1}(3,3)-1\)。沿原Gaussian坐标t，若 \(s_i=a_it+b_i\)，每像素至多一个阈值根，ai<0时不等号反向。排序至多256个根后，每段计数恒定，prior质量用真实conditional Gaussian CDF差相加；ai=0或q*越界需单独判定，原子阈值使用≥。

这给 \(\alpha_j=\pi_0(C_j\cap E)\)，而不是posterior事件概率。属性仍活动时外层还要对属性conditional prior积分；不能把逐像素Bernoulli边缘当独立概率相乘。平滑均值很接近不保证硬事件接近。

**来源/实现。** [R8P event/prior slicing](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/03_certificates.tex:67)；[R9 1D event masses](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:76)；[R9 uniform summary](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:12)。

## F-PROB-020｜已使用RRR：predictable residual不是GT covariance PCA

将任务损失权重吸入 \(U=W^{1/2}(M-m_0(Y))\)，用只依赖Y的居中特征ψ，设 \(B=E[\psi\psi^T]\)、\(C=E[U\psi^T]\)。对于rank≤k的修正A，考虑

\[
\min_{\operatorname{rank}A\le k}E\|U-A\psi\|^2+\lambda\|A\|_F^2.
\]

令 \(B_\lambda=B+\lambda I\succ0\)，配方得 \(\|AB_\lambda^{1/2}-CB_\lambda^{-1/2}\|_F^2+\mathrm{const}\)。截断SVD因此给

\[
A_k=[CB_\lambda^{-1/2}]_kB_\lambda^{-1/2}.
\]

λ=0时训练/总体最优线性目标的收益为前k个奇异值平方和；λ>0该谱收益对应惩罚目标，不可直接当未惩罚验证风险改善。代码使用 \(\lambda=\mathrm{ridge}\cdot\operatorname{tr}(B)/p\)，并保留输入/输出均值作为intercept。GT residual PCA只看 \(E[UU^T]\)，没有可从Y预测的交叉矩C，不能替代RRR。历史嵌套OOF没有收益，这个负结论只约束冻结线性类，不否定一般非线性表示。

**源码。** [R5 `covariance_fit/prediction`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/rrr_controls.py:15)；[nested RRR](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/nested_rrr.py:51)。

## F-PROB-021｜feature/action蒸馏及有限风险恒等式

R5 teacher四通道为 \(T_q(Y)\)，student provider \(g_\theta(Y)\)。归一化feature损失是 \(L_F=E\|S_T^{-1}(g_\theta-T_q)\|^2\)。若冻结head为 \(h(Y,T)\)，记student动作 \(a_S=h(Y,g_\theta)\)、teacher动作 \(a_T=h(Y,T_q)\)，则R6使用

\[
L_{GT}=E\|a_S-M\|^2/s_A,\quad
L_{F+A}=L_F+E\|a_S-a_T\|^2/s_A,\quad
L_{GT+A}=(E\|a_S-M\|^2+E\|a_S-a_T\|^2)/s_A,
\]

其中sA是冻结训练目标component variance标度。对任何同一评价分布，记d=aS−aT、e=aT−M，则

\[
R_S-R_T=E\|d\|^2+2E\langle e,d\rangle,
\quad |R_S-R_T|\le D^2+2\sqrt{R_T}D,\quad D^2=E\|d\|^2.
\]

**推导。** 展开 \(\|e+d\|^2\)，再Cauchy–Schwarz。若head对T具有已证明统一Lipschitz常数L，可再用 \(D\le L\|g-T_q\|_{L^2}\)；本项目没有提供紧全局深网L。teacher不是Bayes真值，故交叉项通常不为0；有限VAL展开只是empirical identity，不能升级为真实posterior风险分解。

**源码。** [R5 width128 feature蒸馏](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/distill_compact_w128.py:19)、[loss](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/distill_compact_w128.py:38)；[R6 action losses](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/track_b/run.py:47)；[R6 empirical expansion](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/track_b/report.py:24)。

## F-PROB-022｜整个cell的响应范数界转likelihood界

若对整个cell C有 \(\|f(x)-f(a)\|\le e_C\)，令 \(r_a=f(a)-y_w\)，三角不等式给

\[
\max(0,\|r_a\|-e_C)\le\|f(x)-y_w\|\le\|r_a\|+e_C,
\]

\[
\ell_C=e^{-(\|r_a\|+e_C)^2/2}\le L(x)\le
u_C=e^{-\max(0,\|r_a\|-e_C)^2/2}.
\]

**关键条件。** eC必须覆盖整个cell、全部测量和正确白化；锚点残差、少量节点误差或中心条件数不能替代它。若 \(\mathcal P\) 用公共shift P*缩放，ℓ、u、Z全部乘 \(e^{P_*}\)，任务比值不变；跨记录恢复evidence时 \(\log Z=\log Z_{\rm scaled}-P_*\)。

**源码/来源。** [R8P cell residual bounds](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/04_maxwell.tex:5)；[R9 `global_audit`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:86)、[1D `cell`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:54)。

## F-PROB-023｜evidence分母、遗漏质量与数值代理的区别

若不重叠cells覆盖全部原prior，\(p_j=\pi_0(C_j)\)，且 \(\ell_j\le L\le u_j\)，则

\[
Z_L=\sum p_j\ell_j\le Z\le\sum p_ju_j=Z_U.
\]

遗漏域U的evidence上界 \(B_U\ge\int_U Ld\pi_0\)，已处理域D的正下界 \(Z_{D,L}>0\)，给

\[
\pi(U\mid y)\le\frac{B_U}{Z_{D,L}+B_U}\le B_U/Z_{D,L}.
\]

只有prior尾质量δ而L≤1时只能给BU≤δ；高SNR下Z极小，δ小未必表示posterior tail小。若U为所有被排除cell的并，整体prior质量≤1，且每个cell满足 \(\mathcal P\ge P_*+\gamma\)，则 \(B_{U,\rm scaled}\le e^{-\gamma}\)；这不需要给每个cell各分配质量1再相加。

另一量 \(Z_L/Z_U\) 是保证纳入的子概率质量，\(1-Z_L/Z_U\) 是未分配预算，不是已发现新模式的质量。R9 denominator为 \(\min Z_{\rm numerical}-2(\max Z_{\rm numerical}-\min Z_{\rm numerical})\) 的有限细化代理，**不是证明的Z下界**；因此相应遗漏posterior值只标numerical proxy，strict字段null。

**来源/源码。** [R8P evidence与遗漏界](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/03_certificates.tex:4)；[R9 whole-support audit](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:86)；[v3 qualification postprocess](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_qualification_v3.py:10)。

## F-PROB-024｜正负任务矩与归一化中心化包络

对有界任务τ，定义原prior正负矩

\[
A_j^+(c)=\int_{C_j}(\tau-c)_+d\pi_0,\quad
A_j^-(c)=\int_{C_j}(c-\tau)_+d\pi_0.
\]

由ℓ≤L≤u分别控制正负部分，得

\[
H_L(c)=\sum_j[\ell_jA_j^+(c)-u_jA_j^-(c)]
\le\int(\tau-c)Ld\pi_0
\le H_U(c)=\sum_j[u_jA_j^+(c)-\ell_jA_j^-(c)].
\]

若ZL>0，\(H_L,H_U\)各自零点cL、cU满足 \(c_L\le E_\pi\tau\le c_U\)。

**推导与sharpness。** 每个可行likelihood的中心化积分是 \(Z_L^{\rm actual}(E\tau-c)\)，对c严格递减，斜率≤−ZL。为了最小化积分，在τ>c处取ℓ、τ<c处取u；最大化相反，故上下零点是给定包络信息类的可达端点。该“最紧”并不保证极值likelihood能由Maxwell产生。

若A矩也有区间，安全下函数取 \(\ell A_{+,L}-uA_{-,U}\)，上函数取 \(uA_{+,U}-\ell A_{-,L}\)。只在符号被严格确定时排除根区间，不能对未控制误差的moment点估计声称严格bisection。

**源码/来源。** [R8P theorem](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/03_certificates.tex:25)；[R9 `probability_envelopes.run`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_envelopes.py:13)的prior矩是GL数值点估计；[1D `H`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:67)使用quintic/Gaussian矩，仍无定向舍入资格。

## F-PROB-025｜硬事件保持共享分母的sharp区间

令 \(\alpha_j=\pi_0(C_j\cap E)\)，\(p_j-\alpha_j\) 是同cell补事件质量。则

\[
P_L=\frac{\sum\ell_j\alpha_j}{\sum\ell_j\alpha_j+\sum u_j(p_j-\alpha_j)},\quad
P_U=\frac{\sum u_j\alpha_j}{\sum u_j\alpha_j+\sum\ell_j(p_j-\alpha_j)}.
\]

**推导。** 概率为事件evidence/(事件+补事件evidence)，分别让事件端取下/上、补事件取上/下即得。若α只有区间且pj精确，下端取αL、上端取αU。区间分子与分母互相包含，不应把N与Z当独立误差任意相减。任务在cell内恒为0/1，只省了任务边界分辨率，层间evidence仍决定概率。

**源码。** [R8P event公式](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/03_certificates.tex:67)；[R9 1D EL/EU/CL/CU](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:76)；[二维numerical envelope](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_envelopes.py:13)。

## F-PROB-026｜物理近似误差与任务中心化，不能漏掉求积项

令 \(\widetilde L\) 为近似物理likelihood，\(\hat N,\hat Z\) 为其数值积分，\(\hat\mu=\hat N/\hat Z\)。在 \(\hat Z>0,Z\ge Z_L>0\) 下，有恒等式

\[
\mu-\hat\mu=\frac{\int(\tau-\hat\mu)L\,d\pi_0}{Z}.
\]

若中心化求积误差满足

\[
\left|\int_D(\tau-\hat\mu)\widetilde Ld\pi_0-(\hat N-\hat\mu\hat Z)\right|\le\epsilon_{Q,c},
\]

则由括号内最后一项为0得到

\[
|\mu-\hat\mu|\le\frac{\epsilon_{Q,c}+\int_D|\tau-\hat\mu||L-\widetilde L|d\pi_0+\int_{D^c}|\tau-\hat\mu|Ld\pi_0}{Z_L}.
\]

这分别包含求积、物理、遗漏域三项。若每cell \(|\log L-\log\widetilde L|\le\delta_j\)，则 \(|L-\widetilde L|\le\widetilde L(e^{\delta_j}-1)\)。若只给统一δ且积分精确，bounded likelihood tilt还给 \(\mathrm{TV}(\pi,\widetilde\pi)\le\tanh(\delta/2)\)；但单点forward误差不能提供统一δ。

**来源。** [R8P对旧task-centered物理误差的审计](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/01_audit.tex:39)、[cell/task envelope实际统一处理](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/03_certificates.tex:25)。R9不填strict总误差证书，原因正是这些项未全部外包。

## F-PROB-027｜任务停止、事件动作与新观测价值

若可信任务区间 \([\mu_L,\mu_U]\) 半宽≤ε，中点均值误差≤ε，平方loss excess risk≤ε²。事件误报/漏报代价为cFP、cFN时，动作1阈值为 \(p\ge c_{FP}/(c_{FP}+c_{FN})\)；整个概率区间在阈值同侧即可停止该动作，但概率值可能仍不精确。

对一般有限动作集合，有可信风险区间时，若 \(R_U(a)-\min_bR_L(b)\le\epsilon\)，该动作regret≤ε。真正新增数据Yd才有

\[
\operatorname{Var}(T\mid y)-E[\operatorname{Var}(T\mid y,Y_d)\mid y]
=\operatorname{Var}(E[T\mid y,Y_d]\mid y)\ge0.
\]

**推导。** 前两式分别来自F-PROB-004和线性事件风险；最后式为条件total variance。重新打包同8源并不新增观测。计算完成只缩小数值动作误差，不会降低固定y的真实posterior variance。

**来源。** [R8P stopping/action](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/06_algorithm_cost.tex:29)、[new data value](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/06_algorithm_cost.tex:79)；[1D numerical停止轨迹](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:79)。

## F-PROB-028｜固定ROI均匀复幅度的精确传递

设field方程 \(A_hE_h=E^{inc}\)，\(A_h=I-KD_h\)，更新 \(D_{\chi(v)}=D_h+vPP^T\)，P选256个ROI。预计算

\[
B=A_h^{-1}KP,\quad C=P^TB,\quad E_R=P^TE_h,\quad
O=S(P+D_hB),\quad F_h=SD_hE_h.
\]

场差满足 \(E-E_h=BvP^TE\)，故 \(P^TE=(I-vC)^{-1}E_R\)，而current差为 \(D_h(E-E_h)+vPP^TE\)，得到

\[
F(v)=F_h+O\,v(I-vC)^{-1}E_R.
\]

它保留全部离散多次散射，不是Taylor/Born截断。必须包含原非零host \(.1+.04i\)；shape固定但latent活动维数小不等于材料更新rank1。改变ROI形状或非均匀幅度时，标量v与C的Schur复用不再直接适用。

**源码。** [R9 `Physics.setup`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_transfer.py:116)、[ `Physics.reduced/direct`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_transfer.py:143)；[R8P exact amplitude theorem](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/04_maxwell.tex:89)。完整物理推导另见02章，此处只固定概率节点使用的接口。

## F-PROB-029｜Schur全cell majorant、导数与每概率节点费用

令 \(C=QTQ^*\)，\(O'=OQ,E'=Q^*E_R\)。则 \(F(v)=F_h+O'v(I-vT)^{-1}E'\)，

\[
F'(v)=O'(I-vT)^{-2}E',\qquad
F(v)-F(a)=(v-a)O(I-vC)^{-1}(I-aC)^{-1}E_R.
\]

对中心a、半径r的圆盘，\(d_i=|1-aT_{ii}|-r|T_{ii}|>0\)，vmax=|a|+r。若b是右端绝对值，逆序递推

\[
x_i=\frac{|b_i|+v_{\max}\sum_{j>i}|T_{ij}|x_j}{d_i}
\]

由逐行三角不等式包住整个圆盘的 \(|(I-vT)^{-1}b|\)。对|E′|应用两次，再乘|O′|并取Frobenius norm，给 \(e_C=\sqrt2\,r\||O'|x\|_F/\sigma_c\)。d≤0需细分/保留未知域；非正规性可能使界很松。普通Schur残差、Q非酉、基列误差与系数rounding均须另传到输出，才能升级strict。

费用按数学对象分开：setup=1个576 LU+8 host源RHS+256基列RHS+1个256 Schur；每v节点=8个256三角RHS+读出；majorant还需两次256逆序绝对值递推。RHS意义不同，reduced工作不能记“0物理费用”。

**源码。** [R9 `majorant_derivative`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_transfer.py:159)、[ `reduced`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_transfer.py:143)、[ `Book.reserve/span`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_transfer.py:85)；[R8P majorant/cost](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/report_source/04_maxwell.tex:123)。

## F-PROB-030｜原P1一维正交条件模型：Gaussian尾部变平台原子

本轮追加条件模型固定 \(x_{1:63}=0,x_{64:65}=0\)，native属性均为.5，只保留原 \(x_0=z\sim N(0,1)\) 未截断。于是

\[
D=1.3+.105i,\quad v(z)=D H(1.3z),\quad A=1/1.3,
\]

\[
\mu_v=\Phi_N(-A)\delta_0+\Phi_N(-A)\delta_D+
(v)_\#(1_{-A<z<A}\varphi(z)dz).
\]

任意任务的evidence分子（τ=1给分母）为

\[
N_\tau=\Phi_N(-A)\{\tau(0)L(0)+\tau(D)L(D)\}
+\int_{-A}^{A}\tau(v(z))L(v(z))\varphi(z)dz.
\]

令t=ΦN(z)，有限积分亦等于 \(\int_{\Phi_N(-A)}^{\Phi_N(A)}\tau(v(\Phi_N^{-1}t))L(v(\Phi_N^{-1}t))dt\)，构成不同坐标的普通数值交叉。均值任务为.1+1.3H(1.3z)，ROI硬事件在此均匀slice中恰为z≥0；它是很容易的阈值配置，不能代表一般空间事件。固定观测为已暴露R6 REGULAR，不是原Full66边缘。

**源码。** [R9 1D contract与integrals](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:6)、[CDF coordinate](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:23)、[exact polynomial moments](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_original_1d.py:48)。Gaussian无限tail质量被显式纳入，未以[-4,4]改写原prior。

## F-PROB-031｜任务–响应covariance、观测敏感性与源层分解

用 \(r=\sqrt2\operatorname{pack}F\)、\(y_r=\sqrt2\operatorname{pack}y\)，原势能为 \(\|r-y_r\|^2/(2\sigma_c^2)\)。允许微分穿过积分时，商法则给

\[
\nabla_{y_r}E[\tau\mid y]=\sigma_c^{-2}\operatorname{Cov}_{\pi_y}(\tau,r).
\]

对白化观测y_w，等价 \(\nabla_{y_w}E\tau=\operatorname{Cov}(\tau,f)\)。这说明“弱响应变化”与任务相关性也有关；一个单点奇异值并不直接给posterior任务误差。

源层标签J可用total covariance分解

\[
\operatorname{Cov}(\tau,r)=E[\operatorname{Cov}(\tau,r\mid J)]
+\operatorname{Cov}(E[\tau\mid J],E[r\mid J]).
\]

R6因此分别保存within、between、total任务响应covariance。pure/transition材料像有重叠，所以between-source项不能改名between-physical-mode项。

**源码。** [R6 `summarize`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/track_a/stratified.py:125)，其中fw=√2packF，再除σc²得到上述未白化r观测的敏感性口径。

## F-PROB-032｜数值资格、严格证书与公平费用分解

严格任务证书需要F-PROB-023的**真分母下界**，加完整prior/strata/尾域、整个cell的响应/likelihood界、可靠prior任务矩及算术外包。finite refinement差、QUADPACK error estimate、浮点Schur/majorant和正数denominator proxy均不自动满足这些条件。

R6旧门为均值max absolute≤1e-5 **且** ≤.01参考SD；源层/冻结事件绝对≤.002。reference代理使用再紧五倍的2e-6、.002SD、.0004。本轮数值达门不改变strict字段null。近零SD必须同时报告尺度退化。

费用定义

\[
C_{total}=C_{setup}+C_{online\ physics}+C_{prior/tasks}+C_{profile}+C_{enclosure}+C_{reference/qualification}+C_{rejected/replay}.
\]

二因子分别是“相同概率规则下direct/同一reduced物理”和“相同物理下不同概率规则”。同GH节点、prior与任务的direct/reduced对照只检验物理复用；不能归于新概率理论。setup+积分不等于完整cold start，共同reference/验证也不能不记。已存在结果只生成费用CSV，不重跑任何计时。

**源码/证据。** [R9 canonical v3](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/results/probability/final_numerical_comparison_v3.json)、[aggregate costs](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/results/probability/aggregate_costs.json)、[ledger](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/audit/probability/ledger.jsonl)。详见本页配套 `04_PROBABILITY_COST_TWO_FACTOR.csv`；未执行格保持 `NOT_RUN`，无完整cold wall处留空。

## F-PROB-033｜GN适应Gaussian只作求积坐标：GH importance校正

R9先用真实forward残差定位 \(\hat v\in\mathbb R^2\)，取GN矩阵 \(H_{GN}=J^TJ\succ0\)，\(C=H_{GN}^{-1}=LL^T\)。定义辅助密度 \(q(v)=N(\hat v,C)\)，则原材料积分恒等地改写为

\[
Z_{\rm tr+pure}=E_{q}\!\left[\frac{L_y(V)\{\rho_{\rm tr}(V)+\rho_{\rm pure}(V)\}}{q(V)}\right].
\]

若 \((g_i,w_i)\) 为 \(\int e^{-g^2}h(g)dg\) 的Gauss–Hermite规则，节点 \(v_{ij}=\hat v+\sqrt2L(g_i,g_j)^T\)，缩放权重为

\[
a_{ij}=\frac{w_iw_j}{\pi}\,e^{P_*-\mathcal P(v_{ij})}
\frac{\rho_{\rm tr}(v_{ij})+\rho_{\rm pure}(v_{ij})}{q(v_{ij})},
\quad \log q=-\log(2\pi)-\tfrac12\log\det C-\tfrac12\|L^{-1}(v-\hat v)\|^2.
\]

材料mean/variance、source概率与事件取同一组归一化 \(a_{ij}/\sum a\)，host原子另加。**推导。** 先乘除q，再标准化Gaussian；GH的 \(\sqrt2\) 与 \(1/\pi\) 来自二维标准Gaussian。这里GN不是精确势能Hessian也无妨，因为节点上仍用真实L/q修正；但有限阶GH仍有求积误差与域外遗漏，不能由坐标良好推断reference合格。直接把q或profile当posterior会删除真实likelihood非线性与prior fibre。

**源码。** [R9 `locate`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:18)、[ `gh`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:27)；原R6 [ `local_material_reference.main`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/track_a/local_material_reference.py:22)。

## F-PROB-034｜已知交换对称quotient的分支条件概率

既有Z-ALG二单元例有材料 \((x_1,x_2)\in[.5,4]^2\)、非对称prior \(\rho(x_1,x_2)=x_1/C\)，forward在交换下不变。对有序根 \(a>b\)，两个原材料分支是(a,b)、(b,a)，故

\[
\rho_Q(a,b)=\rho(a,b)+\rho(b,a)=(a+b)/C,
\quad \Pr(x_1=a\mid\{a,b\},y)=a/(a+b).
\]

**推导。** 每个交换分支Jacobian绝对值为1，两个分支L相同，所以Bayes条件概率中的L抵消，只剩原prior比。于是

\[
E[x_1\mid\{a,b\},y]=\frac{a^2+b^2}{a+b},\qquad
\Pr(x_1>x_2\mid\{a,b\},y)=a/(a+b).
\]

若保留正确总densityρQ却错误给两分支各1/2，evidence完全相同，事件却恒为1/2。固定无序pair{1,3}上的分支权重为1/4、3/4；有限噪声时还须对整个有序quotient积分，其事件概率不是固定pair的3/4。这个反例说明evidence正确不能单独资格材料任务。

**源码。** [既有 `integrate_ordered_quotient`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/inputs/Research_P3LR8P_Research_Kit/inputs/zalg/Research_Z_ALG_Research_Kit/zalgo_probability/e4_two_cell_posterior.py:148)；固定pair与有限噪声结果引用既有记录，没有本次重算。此对称例不提供一般非对称Maxwell图册。

## F-PROB-035｜有限求积权重、population variance与K32的ddof0

对给定非负有限求积 \(a_n\)，\(p_n=a_n/\sum a_n\)，将材料向量实化为\(m_n\)，数值输出为

\[
\hat\mu=\sum_np_nm_n,\quad
\hat\Sigma=\sum_np_n(m_n-\hat\mu)(m_n-\hat\mu)^T,
\quad \hat P(E)=\sum_np_n1_E(m_n).
\]

这描述有限积分近似的矩，不是一个连续posterior的TV证书：若真实材料posterior在该层非原子，有限节点测度与该连续测度互相奇异，在该层通常TV=1。合适的资格是连续积分误差、任务误差/包络或另外构造的连续近似，而不是把节点权重相近称为posterior TV小。

另对q的一个实材料分量M的K个iid样本，历史K32使用

\[
\hat v_{K,0}=K^{-1}\sum_{k=1}^K(M_k-\bar M)^2,
\quad E\hat v_{K,0}=\frac{K-1}{K}\operatorname{Var}_q(M).
\]

**推导。** 展开 \(\sum(M_k-\bar M)^2=\sum(M_k-\mu)^2-K(\bar M-\mu)^2\)，取期望。确定性decoder积分给population q variance，没有凭空乘31/32；有限K统计差异必须随接口口径保留。二者仍都只是近似q的矩，不认证真实posterior。

**源码。** [R9 `summary`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/probability_integrate.py:12)；[R6 `summarize`](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell/track_a/stratified.py:125)；[R5 Gaussian矩与ddof边界](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/THEORY_AND_NUMERICAL_CHECKS.md:21)；[R5 `material_moments.py` population定义](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/material_moments.py:1)。

## 最小局限账

1. **响应残差界与遗漏prior质量是不同输入。** F-PROB-022给每cell likelihood界；还需原prior质量/任务矩及F-PROB-023的正分母。大量小残差候选或宽域“没有找到新峰”均不能替代这些量。
2. R9是普通浮点全支持审查与数值交叉，不是strict certificate。1D upper likelihood floor1e-300是保留下溢不确定性的数值处理，不能认证Hankel/Schur系数或极小概率的有效位数。
3. 原R6先验截断模型与原1D未截断条件模型分别有自己的normalizer/atoms；二者都不提供Full66 posterior。Gaussian decoder矩、RRR、蒸馏风险收益也不提供该资格。
4. 全部historical负结果、first-run下溢报告、qualification修复、未运行direct矩阵格保留。本文新增完整RHS、reduced likelihood调用、LU、benchmark调用均为0。

## 配套费用CSV的读法

`04_PROBABILITY_COST_TWO_FACTOR.csv`只读取既有canonical v3、1D v2与receipt，未调用physics模块。56行区分：共同一次setup、在线同节点direct/reduced积分、独立reference、两轮全支持审查、数值任务包络、下溢修复重放与未执行格。`aggregation_role=COMPONENT`是组件；`NONADDITIVE_VIEW/CUMULATIVE_TRAJECTORY/JOB_SUBTOTAL/CAMPAIGN_TOTAL`等带NONADDITIVE的行与组件重叠，禁止直接把全CSV相加。

`denominator_rigorous=FALSE`贯穿数值comparison；分母列写清有限细化代理或普通quad值，不能读取为严格Z下界。`full_cold_start_wall=NOT_RECORDED`保留完整冷启动未测的事实；既有setup+integral组件只在自己的列出现。direct source积分与direct task envelope均明确`NOT_RUN`。总账仍是3624完整complex RHS、421个完整LU、1个Schur、1452408个reduced RHS；这次公式与CSV整理新增费用全部为0个forward/likelihood/benchmark调用。


---

# 05 · 从固定算子直接推到可校正的压缩算法

本章继续推进“固定载体”的数学后果。目标不是把材料固定在某个弱方向，而是让固定算子承担离线计算，让变化材料只进入必要的小系统，并让误差控制面向实际观测/任务。以下推导不需要再跑一个恒等式数值测试。尚未实施的算法与已运行的精确ROI消元分开。

## F-ERR-001 · 一次对角化为何不能消除所有非线性求解

精确固定宿主消元后写为

$$L_Dh=b_D,\quad L_D=I-DC,\quad b_D=DE_R,\quad F=F_h+Oh,$$

其中 $C,E_R,O,F_h$ 固定，$D=\operatorname{diag}(m_{ROI}-m_h)$ 变化。对于任意两格，

$$[D,C]_{ij}=(d_i-d_j)C_{ij}.$$

因此，在物理 Green 有跨格耦合时，一般非均匀材料不与 $C$ 对易。若一个固定子空间要对**所有对角材料矩阵**都不变，它必须由一组坐标轴张成：证明是让 $D$ 依次取每个坐标投影，子空间内向量的每个非零坐标分量均须留在子空间。若还要求对连接这些格点的稠密 $C$ 不变，除空集/全空间外通常没有这样的公共子空间。

**限制范围：** 上述是对“全部对角材料族”的结构结果，不是证明 P1受限decoder像中不存在特殊小不变子空间。原P1未给出这样一个共同不变性，因此不能无条件沿用均匀幅度的对角化。

均匀材料 $D=vI$ 则不同。对 $C=QTQ^*$ 做一次复Schur后，

$$F(v)=F_h+(OQ)v(I-vT)^{-1}(Q^*E_R).$$

每次仅解三角系统，保留全部256维，不是rank-one更新。此公式已用于本轮R6条件模型；一般Full66材料只能复用固定 $C$，不能复用同一个三角resolvent。

## F-ERR-002 · 固定低秩 Green 近似与小系统

取固定离线近似 $\widetilde C=UV^*$，$U,V\in\mathbb C^{q\times k}$。可以来自截断SVD、物理端口Krylov或其它明确的拟合；此处不预设哪种最好。

由 Woodbury，若 $I-V^*DU$ 可逆，

$$\widetilde h=(I-DUV^*)^{-1}DE_R
=DE_R+DU(I-V^*DU)^{-1}V^*DE_R.$$

这是真实近似算子的精确解。固定 $U,V$；只有 $k\times k$ 系统依赖材料。其复杂度包含 $O(qk^2+k^3+qkn_s)$，不能只报 $k^3$。低秩近似支持仍是全256格对角材料；它没有改变 decoder 或 prior，但改变了forward，故未经校正不能作为原likelihood。

## F-ERR-003 · 原 Maxwell 方程的精确缺陷

定义实际原算子 residual：

$$r=b_D-L_D\widetilde h=D(C-\widetilde C)\widetilde h.$$

直接代回得

$$h-\widetilde h=L_D^{-1}r,\qquad F-\widetilde F=OL_D^{-1}r.$$

因此，低秩近似的价值不由 $\|C-\widetilde C\|$ 单独决定，而由 **材料加权、当前场激发、最终接收端口**共同决定。左端 $D$ 与右端 $\widetilde h$ 均不可丢掉。计算真实 $r$ 要付出一次 $C\widetilde h$ 的费用；这仍可能远低于重新LU，但不能记为免费。

## F-ERR-004 · 可收敛的缺陷校正与精确 fallback

把 $\widetilde L=I-D\widetilde C$ 当预条件器，

$$h_{j+1}=h_j+\widetilde L^{-1}(b_D-L_Dh_j),\qquad
 e_{j+1}=T_De_j,\quad T_D=\widetilde L^{-1}D(C-\widetilde C).$$

若有可计算上界 $\|T_D\|\le\gamma<1$，则

$$\|e_j\|\le\frac{\|h_{j+1}-h_j\|}{1-\gamma}.$$

这是由 $e_j=(I-T_D)^{-1}(h_{j+1}-h_j)$ 的Neumann级数得到。非正规系统中只看 $\rho(T_D)<1$ 不能把该范数界替换成同样常数；无法得到有用 $\gamma$ 时，可以使用预条件Krylov及真实residual，或回退完整256LU。**本章没有假设这种迭代一定更快。** 原P1规模小，直接256LU很强，只有整个求解链成本比较才能决定是否值得实施压缩。

## F-ERR-005 · 接收任务的双边误差表示

对线性观测任务 $\ell(F)=\operatorname{Re}\operatorname{tr}(W^*F)$，令

$$L_D^*Z=O^*W.$$

则严格恒等式为

$$\ell(F)-\ell(\widetilde F)=\operatorname{Re}\operatorname{tr}(Z^*r).$$

若使用近似dual $\widetilde Z$，并有 $\|L_D^{-1}\|\le\beta$，令 $r_d=O^*W-L_D^*\widetilde Z$，则

$$\left|\ell(F)-\ell(\widetilde F)-\operatorname{Re}\operatorname{tr}(\widetilde Z^*r)\right|
\le\beta\|r_d\|_F\|r\|_F.$$

证明来自 $Z-\widetilde Z=L_D^{-*}r_d$ 和Cauchy–Schwarz。dual也随材料变化；固定离线dual只可作近似或预条件。要节省成本，应只为少量实际任务求dual，不为每个材料格点建立巨大免费证书。没有稳定性上界时，该式仍是误差表示，但不是可认证误差区间。

## F-ERR-006 · 从物理误差传到 likelihood，而非追求机器精度

设白化预测误差 $\|f_w-\widetilde f_w\|\le\varepsilon$，近似 residual $\widetilde r=\widetilde f_w-d_w$。由三角不等式，

$$\tfrac12\max(0,\|\widetilde r\|-\varepsilon)^2\le\Phi\le
\tfrac12(\|\widetilde r\|+\varepsilon)^2.$$

也有

$$|\Phi-\widetilde\Phi|\le\|\widetilde r\|\varepsilon+\tfrac12\varepsilon^2.$$

似然为 $L=e^{-\Phi}$，直接反转两端指数即可得 $L$ 包络。与接受阈值或任务积分区间比较后，才能决定是否继续解物理方程。注意 $\varepsilon$ 已包含 $\sqrt2/\sigma$；原始复预测误差小不等于高SNR下likelihood误差小。

## F-ERR-007 · 固定 receiver 投影遗漏的信息

设 $S=U\Sigma V^*$，保留前 $p$ 个左奇异方向，$P_U=U_pU_p^*$。

$$\|(I-P_U)F\|_F=\|\Sigma_{tail}V_{tail}^*J_c\|_F
\le\sigma_{p+1}(S)\|J_c\|_F.$$

这是**输出截断**界，不是材料误差界。令复 residual $R=(F-Y)/\sigma\in\mathbb C^{32\times8}$，原512维实白化 residual 为 $r=\sqrt2\operatorname{pack}R$，则

$$\|r\|^2=2\|P_UR\|_F^2+2\|(I-P_U)R\|_F^2.$$

丢掉第二项就改变了likelihood，除非它与材料无关；即使尾部预测很小，尾部数据噪声与变化预测的交叉项也须处理。对于严格receiver-null current，数据的确不变，但 $J_c=D_m(E^{inc}+KJ_c)$ 仍是所有源共享的硬构成关系，不能省略。

## F-ERR-008 · 投影式 TSOM 的正确压缩对象

固定接收弱投影 $P_w$ 后，真正影响内部场的是 $KP_w$。对其作rank-$k$近似或选固定Fourier作用基，控制的是

$$\|K P_w-KP_w Q_kQ_k^*\|,$$

而不是计算两个子空间的几何交集。若先选 $V_k$ 再投影，$P_wV_k$ 通常不正交；应一次QR/SVD得到 $Q$，保留其实际rank。串联投影 $P_wP_b$ 一般不是正交投影，因为 $(P_wP_b)^*=P_bP_w$；只有两者对易时才是交集投影。使用投影复合是合法算法构造，但须明确改变了哪些方向及其metric。

## F-ERR-009 · NN 的数学位置：预测状态，不替换物理目标

如果未来用网络预测 $\widehat h_\theta(D,E_R)$ 或低维系数，其资格与上述低秩近似相同：

$$r_\theta=b_D-L_D\widehat h_\theta,\qquad
h=\widehat h_\theta+L_D^{-1}r_\theta.$$

可以先用便宜预测，再用真实residual决定继续Krylov或完整回退；也可以由 F-ERR-005/006 控制任务。网络训练误差不是 $\beta$ 或全域置信界，不能把经验误差直接代入严格证书。这样网络承担难以解析表示的状态映射，原likelihood、prior和材料支持保持明确。当前未训练此网络。

## F-ERR-010 · 计算是否值得做的代数门

设一次固定setup成本 $C_0$，原每状态成本 $c_f$，替代路径含fallback的平均成本 $c_r$。对同一固定采集的 $n$ 个实际状态：

$$C_{new}=C_0+nc_r,\qquad C_{old}=nc_f.$$

盈亏平衡 $n>C_0/(c_f-c_r)$，只在 $c_f>c_r$ 时成立。达到20%费用门要求

$$n\ge\frac{C_0}{0.8c_f-c_r},\qquad c_r<0.8c_f.$$

材料优化产生的状态数也可能变化，所以完整系统需要实际 $n_{new},n_{old}$；不能用同state核函数加速比冒称整个重建加速。R9精确ROI核已有匹配对照。后续E5按[正式冻结协议](../config/possibility_volume_protocol.json)在32个新对象运行完整anchor和候选云；全部方法共用ROI后端，因此它不提供同对象direct576对ROI256的整系统随机化速度比较。最终执行和覆盖结论以E5封存统计为准。

## 由推导直接得到的决策

- 优先保留精确固定宿主消元；它不改变材料/概率目标，适用于完整原ROI支持。
- 均匀幅度可继续固定Schur；一般非均匀材料不可强行沿用标量可交换性。
- 固定SOM/TSOM投影用于组织current或预条件，不把非法电流当材料样本。
- 暂不新增压缩rank sweep或NN训练：只有在256LU确为剩余瓶颈，且上述缺陷校正有实际余量时，才做一项完整端到端裁决实验。
- 若做该实验，强基线须为精确256ROI solver，而不是旧576系统；固定与动态探索均共享此优化。比较完整anchor+探索+任务成本，保留失败与fallback。

**来源与新颖性：** 本章为Woodbury、residual correction、dual-weighted residual与正交投影性质在本合同下的显式推导。它们是标准工具；尚无独立Maxwell新定理主张。源码对应精确分支 `src/stable_roi_operator.py`，压缩/NN分支仅推导未实验。


---

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


---

# 07 · 材料可能性体积：候选、覆盖、规模与成本

本章的主目标是：从一份观测产生多幅原材料模型内、与数据相容的候选，使真实材料在这些候选的邻域内出现的频率足够高。准确 posterior 是可选的更强目标，不是开始验证可能性体积成像的前置门。有限候选集合可以有真实实用价值，同时必须避免三种虚增：用越来越大的集合换覆盖却不报告体积，用 GT 选择候选却称可部署结果，以及把启发式权重称为后验概率。

下文给出完整定义与代数先行的候选方法。当前E5实施身份为已批准的 [最终生成/运行协议](../config/possibility_volume_protocol.json) `R9_POSSIBILITY_VOLUME_V1` 与 [最终统计协议](../config/possibility_statistics_frozen.json) `R9_POSSIBILITY_VOLUME_STATISTICS_V1`；`DRAFT`文件只保留设计历史，不作为当前实施依据。本文不读取E5结果或GT，运行/封存状态由对应正式receipt给出。以下标注 `NOT_IMPLEMENTED` 的备选不属于E5，也不修改已冻结的R7/R9历史门。

## F-VOL-001 · 有限候选云及其资格

对观测 $y$，算法输出至多 $K$ 个不同材料

$$\mathcal C_K(y)=\{m_1(y),\ldots,m_{K_y}(y)\},\qquad 0\le K_y\le K.$$

每个 $m_j=\mathcal G(x_j)$ 必须处于原 decoder 像、保留固定宿主/ROI 支持，并用原非线性物理解验证。设数据相容集合为 $\mathcal A_\tau(y)$，实际合格云为

$$\mathcal C_K^{qual}(y)=\mathcal C_K(y)\cap\mathcal A_\tau(y).$$

不能用同一个 anchor 的重复拷贝补齐 $K$，也不能把失败候选从分母静默移走。应记录提案数、完整求解数、合格数、去重后数、终止原因与成本。候选云并不天然是一张分布；在高维连续材料空间中，有限点集合的体积为零，所谓“体积”需要明确邻域厚度或另一个集合构造。

## F-VOL-002 · 覆盖的事件与 ε 管状集合

采用原全576格复材料 RMS 距离

$$d_M(m,m')=\|m-m'\|_2/24.$$

对真实材料 $m^*$，定义候选距离与覆盖指示

$$D_K(y,m^*)=\min_{m\in\mathcal C_K^{qual}(y)}d_M(m,m^*),\qquad I_{K,\varepsilon}=\mathbf1\{D_K\le\varepsilon\}.$$

空集时 $D_K=+\infty$。这等价于 $m^*$ 属于

$$\mathcal V_{K,\varepsilon}(y)=\bigcup_{m\in\mathcal C_K^{qual}(y)}B_M(m,\varepsilon),$$

即以候选为中心的材料邻域并集。实验中的经验覆盖为 $\widehat C_{K,\varepsilon}=n^{-1}\sum_{i=1}^nI_{K,\varepsilon}^{(i)}$。它必须注明 parent来源、噪声分布和 $\varepsilon$；不能从有限样本升级为全材料族覆盖定理。

固定前缀顺序时，$D_{K+1}\le D_K$，所以更多候选不会降低 best-of-K覆盖。这是集合单调性，不是算法学到更多物理的证据。必须同时报告 $K=1,2,4,8,16$ 等预定前缀、实际计算预算和至少一组预定 $\varepsilon$。顺序应按半径/方向交错，不能看 GT 后重新排序。

## F-VOL-003 · 覆盖与体积必须成对报告

候选云直径

$$\operatorname{diam}_M\mathcal C=\max_{a,b\in\mathcal C}d_M(a,b)$$

以及相对 anchor 的最大半径、最近邻间距、去重后数量，描述了云是否坍缩或无限扩散。一个更大的云更容易碰到 GT，但也更不确定。二维材料图的逐格最小/最大包络可作可视化；其笛卡尔积可能包含大量不满足 decoder 或物理的组合，不能称那些组合都是合法材料。

仅有有限点的 covariance 或凸包，也不自动给出真实体积：若点数 $K$ 小于材料局部维数，covariance秩至多 $K-1$；在非线性材料支持中，点间线性插值还可能非法。可报告 covariance有效秩作为离散分散度，但不把其行列式当作 Full66 posterior体积。

推荐联合结果是：同一资格下的 $C_{K,\varepsilon}$、$D_K$分布、云直径/实际半径、合格率/去重数量、总费用。若在相同 $K$ 与费用下覆盖改善只是来自更大半径，应明确展示这一取舍。

## F-VOL-004 · GT覆盖容量与可部署决策不同

离线 oracle 最佳候选为 $j_*(y,m^*)=\arg\min_jd_M(m_j,m^*)$。其损失衡量云的表示容量，但需要 GT，不可部署。部署时的选择规则 $\widehat j(y)$ 只能依赖观测、候选、预先固定的目标和预算，不能读取 $m^*$。

应分别报告

$$R_{oracle,K}=\mathbb E\min_jd_M(m_j,M)^2,\qquad R_{deploy}=\mathbb E d_M(m_{\widehat j(Y)},M)^2.$$

显然 $R_{oracle,K}\le R_{deploy}$，差距是选择难度，不是可以免费兑现的成像收益。若输出候选加权平均，材料 decoder 像非凸时该平均可能不是可实现材料；可将其称决策估计，并另行验证或投回原支持，不能自动把它列为物理候选。

主实验可先证明有限云的 coverage/capacity，而不假装已经解决部署选择。若论文标题或产品主张涉及最终单幅重建，必须增加不看 GT 的选择器结果。

## F-VOL-005 · 数据资格：绝对噪声一致性与历史相对竞争门

原噪声为256个独立复观测，每个复方差 $\sigma^2$。对正确物理模型生成的真实材料 $M$，F-BASE-002 的 residual 满足

$$r(M)\sim\mathcal N(0,I_{512}),\qquad \|r(M)\|^2\sim\chi^2_{512}.$$

因此可预定 $\alpha$，用

$$\mathcal A_\tau(y)=\{\mathcal G(x):\|r(x;y)\|^2\le\tau\},\qquad \tau=\chi^2_{512}(1-\alpha).$$

在精确模型和该噪声分布下，$\Pr\{M\in\mathcal A_\tau(Y)\}=1-\alpha$。这是数据相容**集合**的重复噪声资格，不是云已找到其中GT邻域的保证，更不是posterior credible mass。模型误差、噪声估计误差、相关噪声和数据预处理都会改变这一资格。

历史 $\Delta\Phi\le2$ 相对优化 anchor 是另一指标。拟合参数会降低 anchor residual；真实 GT 的 loss 可以合理地高于优化最优值，尤其在许多可拟合方向存在时。因此不能未经分析把该门同时当作“真实材料必然应通过”的噪声资格。旧结果与旧门保持原样；最终E5协议已经采用绝对噪声门，并将旧门通过比例列为辅助诊断。

一个可选的解析保守阈值由标准Gaussian平方和尾界给出：对 $t>0$，$\Pr\{\chi^2_n\ge n+2\sqrt{nt}+2t\}\le e^{-t}$。可取 $n=512,t=\log(1/\alpha)$；精确分位数通常更紧。此界不需要把真实 GT 用于阈值选择。

## F-VOL-006 · 固定弱方向负责提案，完整 corrector 负责合法修正

在 anchor $x_a$，计算完整 whitened $\mathcal J_a$ 与材料导数 $M_a$。令

$$G_a=\operatorname{Re}(M_a^*M_a)/576,\qquad H_a=\mathcal J_a^T\mathcal J_a+I.$$

若按 material单位变化寻找低数据敏感方向，可在可识别材料切向上分析 $\mathcal J_a^T\mathcal J_a v=\mu G_av$；若希望同时抑制 prior代价，则分析 $H_av=\nu G_av$。二者目标不同，必须预先选择。$G_a$ 可能奇异；不能把 decoder零方向当作高价值弱材料方向，须在实际材料非零的切向中求解或显式正则化并说明影响。

取所选 $r_w=8$ 方向的 latent张成空间，做 Euclidean QR 得正交 $Q_w\in\mathbb R^{66\times8}$。这保留空间但不保留逐列 material归一化；实际提案幅度由真实 $d_M(\mathcal G(x_a+t v),\mathcal G(x_a))$ 决定。沿初始射线搜索所需幅度只涉及 decoder，不调用 Maxwell；若给定bracketing与latent cap下未找到目标半径，记录搜索失败，不把它提升为非线性材料全局不可行证明，也不修改材料支持。

**固定B：** 所有候选从同一 anchor frame 生成。**动态C：** 只有在预先声明的 continuation状态/当前候选处刷新 frame，并明确下一次提案如何从该状态产生。若所有候选都直接从同一个 anchor独立发射，而且 corrector完全相同，所谓动态frame并未影响任何动作，B/C应数学上相同；不能靠多算一个分解制造“动态臂”。

**当前实现的准确限制见F-VOL-012：** B/C使用完整66列物理Jacobian，但corrector步限制在 $Q_s$ 加一条weak直线的59维空间，再施加signed-progress切向行；它不是 unrestricted Full66 corrector。弱幅度允许为正或负，最终协议明确历史 `*_RAY` 名称应解释为weak line；signed material progress不等于weak系数为正。这避开旧E2“固定weak幅度加精确半径”的特定几何障碍。R随机控制采用完整66维可行切向，属于更强的corrector几何类。合法材料由decoder保证，数据资格由原forward验证。

## F-VOL-007 · 为何普通MAP修正可能让云坍缩，以及一个平滑修正

如果每个提案最终都无约束地最小化同一个 $\Phi(x)+\|x\|^2/2$，所有轨迹可能返回同一局部MAP。多方向、多半径不能单独阻止坍缩。**当前实现已经使用硬半径与正signed progress防止返回anchor，见F-VOL-012；本节tether状态为 `NOT_IMPLEMENTED`，不能写成当前方法。** 若未来硬半径导致低可行率，可另行考虑到达数据资格立即停止并限制提案附近的trust region，或增加平滑tether。后者为

$$\Psi(x;x_p)=\Phi(x)+\tfrac12\|x\|^2+\frac\lambda2\|\mathcal G(x)-\mathcal G(x_p)\|_M^2.$$

它仍可用 ordinary Full-J增广GN求解，仅多出 decoder导数块；没有新增 Maxwell参数方向。材料残差用实/虚拼接后除以24，不能混入观测的 $\sqrt2/\sigma$。

在 anchor驻点附近，设 $p=x-x_a$、提案位移 $\delta=x_p-x_a$，局部模型为

$$\tfrac12p^TH_ap+\tfrac\lambda2(p-\delta)^TG_a(p-\delta).$$

其最优位移为

$$p_\lambda=(H_a+\lambda G_a)^{-1}\lambda G_a\delta.$$

在 $G_a$ 正定的有效切向及 $H_av=\nu G_av$ 的方向上，位移保留比例为 $\lambda/(\nu+\lambda)$。取 $\lambda=\nu$ 得局部一半保留；取 $\lambda$ 不小于所选弱空间最大 $\nu$，可使该空间每个广义分量局部保留至少一半。这给 tether尺度一个可解释的起点，而不是凭经验为每个GT调参。

这只是局部二次解释，不保证非线性最终半径，也不保证数据门可通过。实际半径、tether目标、物理数据loss、prior代价须分别记录；若所有候选仍坍缩，判该策略失败。$\Psi$ 是候选生成目标，不是新定义的 posterior。

## F-VOL-008 · 全秩混合支持与其有限保证

纯 $Q_w$ 低维提案对其张成空间外没有 latent Lebesgue支持。**本节状态为 `NOT_IMPLEMENTED`；当前B/C没有混合全秩分支，R是独立对照臂，不能借R的方向支持给B/C作保证。** 后续可在 Gaussianized prior坐标采用

$$\delta\sim(1-\eta)\,\mathcal N(0,Q_w\Lambda Q_w^T)+\eta\,\mathcal N(0,\Sigma_0),\quad 0<\eta\le1,\quad\Sigma_0\succ0.$$

它虽含低维分量，但整体对任意非空 latent开集给正概率，因为全秩分量密度处处为正。另一方案是单一满秩 Gaussian协方差 $Q_w\Lambda Q_w^T+\eta I$，此时每个方向都有小的全空间扰动。$\Lambda$ 应在prior尺度上定义；material归一化的非正交列不自动维持标准Gaussian prior。

这里保证的是**提案支持**，不是有限 $K$ 的高覆盖，也不是正确prior抽样。对某个GT邻域 $B$，若一次合格保留提案进入它的概率为 $p_B>0$，独立 $K$ 次至少命中一次的概率才是 $1-(1-p_B)^K$；有限 $p_B$ 可能极小。动态、antithetic或去重的提案并不独立，不能直接套用此式。

此外，归一到有限几个固定半径后，绝对连续支持会变成若干壳面；固定半径列表不能再声称完整latent支持。确定性corrector也可能把满支持提案压到一个点。若需要保留“合格开集仍可被发现”的结构保证，可用概率 $\eta$ 独立抽原完整prior或满支持局部proposal，保留已经通过数据门的原提案，不强制修正；不合格者再进入corrector或按预算记失败。这一分支仍可能接受率很低，不是免费覆盖。

有限 $K=16$ 下可预定少数全秩探索槽位来防止完全零探索，但它不能保证每个邻域被看到。候选物理支持、proposal支持和posterior权重是三个不同问题。

## F-VOL-009 · 去重、信赖区域与费用控制

在预定材料阈值 $\delta_{dup}$ 下，若新候选与已接受候选最小距离小于阈值，则标为重复；可保留数据loss较低者，同时保留原重复记录。候选是否入云只依赖观测与候选，不依赖GT。去重前后数量、重复率和追加尝试必须计入费用；不能不断重试直到凑够16个而称各臂预算相同。

每个候选限制相同最大 forward/full-J次数、line-search次数与总费用；通过数据门后可提前停。若设置提案附近的 trust region，需报告它在latent还是material metric中定义。固定提案半径不是最终半径；correction、边界和prior均可改变后者。

总费用应包含 anchor、固定operator setup、初始/刷新frame、所有失败提案、solver、corrector和验证。固定setup可跨同采集复用，但按每观测与跨观测摊销分别报告。GT离线评估不计入部署推断费用，但也不能从评估中导出一个未付成本的部署选择器。

## F-VOL-010 · 已冻结E5与设计取舍

当前实施以本章开头链接的最终协议为准：32个fresh parent、B/C/R三臂、每臂16个非anchor槽位、四个固定半径。旧8个parent已暴露，不能冒称本轮fresh总体；“新噪声但同材料”与“新材料”也须分别注明。本节保留设计理由，不再把已批准协议写成待冻结建议。

**实际对照：** 共享anchor A、固定弱偏置云 B、动态弱偏置云 C、全空间随机云 R，均采用精确ROI solver、相同噪声资格与半径/费用合同。B/C的corrector几何为弱补空间加一条weak line，R为完整有效66维，不能称三臂corrector空间完全相同。C更新当前允许步空间，具体见F-VOL-012；其作用并非只重算一个不影响算法的frame。平滑tether未实施。

冻结的 $K=16$ 表示非anchor方向槽位，不是保证得到16个不同合格材料；最多再加入一个满足噪声/合法性门的共享anchor。当前采用8个随机弱方向的正负配对，每个半径安排2对，共4个半径、16个提案；连续Gaussian系数在8维弱空间中几乎必然线性独立，但有限角度/半径覆盖仍很稀疏。这只是有限方向设计，不具有F-VOL-008的完整支持保证。若未来需要全秩混合探索，应另行预先冻结槽位安排，不在当前结果出来后额外补点。

半径 $0.02/0.05/0.10/0.15$ 已冻结为广域探索安排，不是已证明处处可达的材料尺度。给定decoder搜索未找到bracket时直接消耗该槽位并记失败，不能读GT后改半径或补点。coverage主容差已冻结为 $\varepsilon=0.075$，辅助容差为0.03/0.05/0.10/0.15；$\varepsilon$ 与提案半径是不同对象。

**parent数量：** 冻结的32个新parent可支持有意义的配对效果描述，但不足以用点估计宣称接近保证的高覆盖。即使32/32全部覆盖，一侧95%二项下界仅为 $0.05^{1/32}\approx0.911$；若要求95%覆盖的一侧95%下界，在零失败时也至少需要59个独立parent。这个算式决定结论强度，不改变已批准32个parent的合同，也不授权扩展59个parent。

如果实际执行因硬预算提前终止，必须保留未完成状态，并遵守最终协议对全部32个parent封存后评分的门；不能把少数完成parent悄然改成原确认性分析。32个parent中的配对差异、失败与coverage区间，比堆更多次计时重复有信息。冻结K前缀为0/1/4/8/16，复用同一批候选，不另跑候选数sweep。

**最终合同的关键项：** parent生成与暴露状态、绝对噪声门、旧门辅助诊断、K与固定前缀、RMS分母、coverage容差、候选与corrector规则、去重阈值 $10^{-6}$、全费用上限、封存后离线GT评估均已写入最终协议。资格云为空时保留parent失败，不作complete-case删样；任务动作可明确标记anchor fallback，但fallback不恢复资格云coverage。

## F-VOL-011 · 近似权重可以有用途，但不能改名为posterior

对已合格候选，可以给出均匀权重、基于loss的软权重或明确的质量评分，用于浏览和排序。若提案密度、corrector映射的多对一质量、拒绝规则和有限覆盖未被处理，则

$$w_j\propto e^{-\Phi(x_j)-\|x_j\|^2/2}$$

通常不是原posterior的离散近似权重：候选出现次数本身已受到提案与优化吸引域影响。去重也会改变出现频率。把这些权重标为“候选评分”即可，不应显示为可信事件概率。

若后续确实需要准确posterior，才进入原prior条件积分、重要性校正、MCMC或已证实的传输/测度处理；相关资格见第04章及 `PRIOR_ART_AND_CLAIMS.md`。当前主方向首先检验：有限物理合格云能否以可接受的大小和成本覆盖更多真实材料。

## F-VOL-012 · 当前候选引擎的静态证明核对

本节依据 `src/possibility_volume.py`、已批准的 [最终E5协议](../config/possibility_volume_protocol.json) 与 [冻结统计协议](../config/possibility_statistics_frozen.json) 静态阅读，不运行求解或读取GT。源身份按最终协议的 `source_files` 记录；本节优先于前文 `NOT_IMPLEMENTED` 备选，DRAFT不再作为当前方法依据。

**提案与模型。** 每parent先从零latent求一个prior权重1的anchor。令 $M_k$ 为完整材料导数实化后除以24的矩阵，$M_k=U\Sigma V^T$。在保留的非零材料切向上取 $T=V\Sigma^{-1}$，对

$$T^T(\mathcal J_k^T\mathcal J_k+I)T$$

求最小8个特征方向并映回latent，得到material归一化 `weaklift`；其Euclidean QR给 $Q_{w,k}$。因此当前弱方向明确考虑data加prior，权重都是1。数值秩阈值定义的是实际实现切向，不是数学证明decoder具有精确零方向。

四个半径 $\rho\in\{0.02,0.05,0.10,0.15\}$ 各有两条独立Gaussian弱方向及其正负号；B/C共享原提案，R使用独立66维Gaussian方向。初始方向 $v_0$ 令 $\|M_av_0\|=1$。decoder-only标量搜索得到 $x=x_a+t v_0$，要求实际材料RMS为 $\rho$，并限制 $\|x-x_a\|\le4$。不同半径独立从anchor开始，不沿半径串行continuation。

**B的弱直线。** 固定 $Q_w=Q_{w,a}$，取 $Q_s$ 为其正交补。每步允许空间

$$\mathcal S_B=\operatorname{range}[Q_s,v_0],\qquad \dim\mathcal S_B=59.$$

由初始 $x-x_a\in\operatorname{span}(v_0)$，若 $p\in\mathcal S_B$，则 $Q_w^T(x+p-x_a)$ 仍平行于 $Q_w^Tv_0$。径向retraction $x^+=x_a+h(x+p-x_a)$ 只乘一个实标量，故继续保留该弱直线，系数允许改变并可为负。它没有锁住弱投影长度，也没有证明保持正向射线。最终协议为保留历史ID而使用 `FIXED_Q8_RAY`，数学含义以 `weak_line_sign_policy` 为准。

**C的逐步弱直线。** 在当前状态重算 $Q_{w,k}$，令

$$a_k=Q_{w,k}Q_{w,k}^T(x_k-x_a),\quad \mathcal S_{C,k}=\operatorname{range}[Q_{s,k},a_k].$$

当 $a_k\ne0$ 时，若 $p\in\mathcal S_{C,k}$，则对某标量 $b$，

$$Q_{w,k}^T(x_k+p-x_a)=(1+b)Q_{w,k}^T(x_k-x_a).$$

径向retraction再乘 $h$，所以保留**当前frame下**的直线，$h(1+b)$不被要求为正。下次frame变化可以改变其坐标；没有跨frame运输，也不将其余7个当前弱坐标清零。若 $a_k$ 过小，代码用 $Q_{w,k}Q_{w,k}^Tv_0$ 作已记录fallback；仍为零则该候选几何失败。C虽不改变初始提案，但会改变corrector的允许空间，因此不是名义上的动态臂。

**signed-row GN与真正的接受规则。** 定义固定材料方向 $u=M_av_0$，$\|u\|=1$。在上述allowed空间内再去除数值材料零方向，并施加

$$u^TM_kp=0.$$

这是signed material progress的一阶保持条件，不是材料RMS半径的一阶切向条件。随后用完整数据Jacobian求解

$$\min_{p\in\mathcal S_k,\ u^TM_kp=0}\ \tfrac12\|r_k+\mathcal J_kp\|^2+\tfrac12\|x_k+p\|^2,$$

限制latent步长与材料线性步长，再作decoder径向retraction。真正评估的步为 $e=x_{trial}-x_k$；代码用该 $e$ 而非retraction前的 $p$ 计算GN预测下降，这避免错用切向模型。只有预测MAP下降为正才付forward；实际接受要求data $\Phi$ 降低，允许总MAP目标上升，并用总MAP实际/预测比调trust。这是明确的混合启发式：不能称MAP单调下降算法，也不能称纯data最小化。

**最终几何与数据资格。** 对每个真实visited state计算 $\Delta m=\operatorname{pack}_{material}(m-m_a)/24$，需满足

$$\big|\|\Delta m\|-\rho\big|\le0.02\rho,\qquad u^T\Delta m\ge0.5\rho,$$

以及decoder合法和 $\Phi\le282.8737921622988=\chi^2_{512}(0.95)/2$。正半径和signed progress排除返回anchor，故当前方法没有无约束MAP坍缩到anchor的问题；不同射线仍可能落到相同/相近材料，须保留去重结果。半径加signed门也可能降低合格率，这是要检验的结果，不是必须救回的实现错误。

每方向在实际visited且几何合法状态中选择最小dataPhi点；若最小者也不满足噪声门，则该slot失败。存储的 `PRIMARY_CLOUD.npz` 是allselected加资格flags，不是已过滤的posterior样本；统计再按flags、合法性和噪声门构造合格云。原始点和所有visited噪声合法点保留为secondary变量K云，不与最多16个非anchor点加合格anchor的主云混算。最终统计已要求不合格anchor从qualified云过滤，同时单独保留无条件anchor baseline与显式action fallback。

**不作过强支持主张。** B/C受初始弱空间、corrector allowed空间、半径壳面与latent半径4限制；R也受后两项限制。三者均不是全prior支持抽样。`original_prior_support=True`应解释为每个输出位于原模型支持内，绝不是算法覆盖整个prior支持。

**已明确命名的搜索失败。** `match_radius`逐渐增大标量并在末端检查距离。一般非线性decoder的沿射线材料距离不保证单调；端点低于 $\rho$ 不能排除两个采样点间曾超过 $\rho$。最终代码使用 `RADIUS_BRACKET_NOT_FOUND_WITHIN_LATENT_CAP`，只表示当前bracketing/cap下未找到目标半径，不能报告为数学不可行证书。失败保留，不追加未授权密集搜索。

**费用。** 32parent×3臂×16提案=1536提案，每提案最多6次forward和4次完整J，费用上界为 $8\times6+32\times4=176$ reduced RHS。加每parent最多64个anchor状态以及32次数据生成：

$$1536\times176+32\times64\times40+32\times8=352512,$$
$$1536\times6+32\times64+32=11296\quad\text{个 reduced LU}.$$

它们位于最终批准的400000 RHS256与12000 LU256上限内；7200秒与内存上限仍可能先终止，不能由RHS算式保证完成。已有264 RHS576 setup作为已支付复用成本，完整研究历史2184 RHS576与当前新增RHS256分别列账，不混称同规模物理费用。

**静态结论。** 该设计作为有界、无GT在线选择的可能性云生成器，在以上声明下没有发现致命代数矛盾。它检验固定/动态受限corrector与更强随机corrector的整体策略差异；不是只隔离proposal方向的因果试验，也不认证posterior、全support或总体高覆盖。是否启动及最终科学判断仍由主线程负责。

## F-VOL-013 · 噪声资格与算法漏覆盖的分解

令 $A_\alpha(Y)$ 为F-VOL-005中原噪声模型下的材料相容集合，$\Pr\{M\in A_\alpha(Y)\}\ge1-\alpha$。把算法本身的漏覆盖量定义为**联合概率**

$$\beta_{K,\varepsilon}=\Pr\{M\in A_\alpha(Y),\ d_M(M,\mathcal C_K^{qual}(Y))>\varepsilon\}.$$

事件按GT是否噪声合格分割，立即得到

$$\Pr\{d_M(M,\mathcal C_K^{qual}(Y))\le\varepsilon\}\ge1-\alpha-\beta_{K,\varepsilon}.$$

该下界还允许 $M\notin A_\alpha(Y)$ 时恰好被某个合格候选的 $\varepsilon$ 邻域覆盖；这些额外命中没有被下界计入。$\beta$不是条件漏失率，不能再乘一次 $1-\alpha$。若选择定义条件漏失 $\widetilde\beta=\Pr\{D_K>\varepsilon\mid M\in A_\alpha(Y)\}$，才有相应 $(1-\alpha)(1-\widetilde\beta)$ 型下界。

这个分解把两类任务分开：噪声阈值决定物理模型下的资格失配 $\alpha$，候选算法决定剩余的搜索/表示漏失 $\beta$。当前方法没有证明 $\mathcal C_K$ 是 $A_\alpha$ 的全局 $\varepsilon$-net，所以 $\beta$必须通过独立parent实测并给出不确定区间；不能因为K变大或有8个弱方向就假定它接近零。该分解不要求posterior准确，也不提供免费的coverage保证。

如未来额外裁去prior概率至多 $\delta$ 的材料区域，可用union bound得到 $1-\alpha-\delta-\beta$ 型下界，其中 $\beta$须相应定义在剩余区域内。当前没有以已证明prior质量为依据的硬裁剪；latent距anchor≤4也不是一个已计算质量的prior置信区域，不能自行填入某个 $\delta$。

## F-VOL-014 · 集合行动风险与覆盖曲线

若云始终非空，定义集合行动损失 $L_{set}(M,\mathcal C)=\min_{m\in\mathcal C}d_M(M,m)^2$。由非负随机变量的尾积分恒等式，

$$R_{set,K}=\mathbb E D_K^2=\int_0^\infty\Pr\{D_K^2>t\}\,dt
=2\int_0^\infty\varepsilon\,[1-C_{K,\varepsilon}]\,d\varepsilon.$$

所以一整条coverage曲线对应一个集合容量风险；单个 $\varepsilon$ 的高覆盖不足以控制尾部损失。若另有全局 $D_K\le R_{max}$，可得 $R_{set,K}\le\varepsilon^2+R_{max}^2\Pr(D_K>\varepsilon)$，但这里的有界性必须真实成立。

若qualified云可能为空，按F-VOL-002定义其距离为无穷，则该风险可以为无穷。最终E5统计保留全部32个parent，空云coverage=false，最近材料误差记未定义/无穷右删失；禁止complete-case主对比，也没有发明有限失败罚值。任务动作使用明确标记的anchor fallback，与集合coverage分开。未来合同可另选有界任务罚损，但它是 `NOT_IMPLEMENTED` 备选。原始anchor始终可用的未过滤云可另报有限容量风险，但它与qualified云不是同一对象。

$R_{set,K}$不是单幅Bayes决策风险：它允许损失计算在看到GT后从集合中取最好者。部署选择器的风险仍是F-VOL-004的 $R_{deploy}$。同时报告集合直径/规模，才能解释降低集合风险付出了多少不确定性。以上是标准集合覆盖与尾积分关系，不作新概率定理主张。

## F-VOL-015 · 任务加权覆盖—分散度界及完整推导

**假设。** 固定一次观测与非空有限合格云 $\{m_j\}_{j=1}^{K_y}$，真实材料为 $M$。任务 $T:\mathcal M\to\mathbb R^p$、固定权重矩阵 $W\succeq0$，写 $z(m)=W^{1/2}T(m)$。假设在所讨论材料域有任务Lipschitz条件 $\|z(m)-z(m')\|\le L_Td_M(m,m')$。候选权重 $w_j\ge0$、$\sum_jw_j=1$，定义

$$\bar z=\sum_jw_jz(m_j),\qquad V_w=\sum_jw_j\|z(m_j)-\bar z\|^2.$$

令 $G_\varepsilon=\{j:d_M(m_j,M)\le\varepsilon\}$，并假设其总权重 $a=\sum_{j\in G_\varepsilon}w_j>0$。则

$$\boxed{\|\bar z-z(M)\|\le L_T\varepsilon+\sqrt{\frac{1-a}{a}V_w}.}$$

**推导。** 当 $0<a<1$ 时，定义两组条件加权均值

$$\bar z_G=a^{-1}\sum_{j\in G_\varepsilon}w_jz(m_j),\qquad
\bar z_B=(1-a)^{-1}\sum_{j\notin G_\varepsilon}w_jz(m_j).$$

由Lipschitz条件、非负权重及三角不等式，$\|\bar z_G-z(M)\|\le L_T\varepsilon$。整体均值为 $\bar z=a\bar z_G+(1-a)\bar z_B$。以各组均值展开平方，组内交叉项为零，得到完整方差分解

$$V_w=\sum_{j\in G_\varepsilon}w_j\|z_j-\bar z_G\|^2
+\sum_{j\notin G_\varepsilon}w_j\|z_j-\bar z_B\|^2
+a(1-a)\|\bar z_G-\bar z_B\|^2.$$

前两项非负，而 $\bar z_G-\bar z=(1-a)(\bar z_G-\bar z_B)$，所以

$$\|\bar z_G-\bar z\|^2\le\frac{1-a}{a}V_w.$$

再对 $\bar z-z(M)=(\bar z-\bar z_G)+(\bar z_G-z(M))$ 用三角不等式，即得结论。$a=1$ 时所有正权重候选都位于GT的ε邻域，直接有 $\|\bar z-z(M)\|\le L_T\varepsilon$；$a=0$ 时该界不适用，即使云中存在权重为零的近GT候选也不能使用。

均匀权重且至少一个近GT候选时，$a\ge1/K_y$，有保守界 $L_T\varepsilon+\sqrt{(K_y-1)V_w}$。另一个恒等式

$$V_w=\tfrac12\sum_{i,j}w_iw_j\|z_i-z_j\|^2$$

把分散度与加权任务pair距离直接联系起来。两式共同说明：更多候选可以改善集合覆盖，却不自动改善加权动作；近GT候选被极小权重压制或其他候选过于分散时，集合容量可能无法兑现为单一决策收益。

**资格与部署限制。** 该结论是标准组间方差分解与三角不等式的条件结果，不作为原创数学定理。$V_w$可在部署时由候选计算，但 $G_\varepsilon$ 与 $a$依赖未知GT，当前没有它们的在线下界。因此本式用于理论解释和封存后的分析，**不是在线误差证书**，不授权读GT、更改冻结权重或追加统计门。硬阈值任务未必Lipschitz，不能无条件套用；空云须按最终E5失败/fallback规则单独处理。任务平均 $\sum_jw_jT(m_j)$ 也不等于一般非线性任务在材料平均上的值 $T(\sum_jw_jm_j)$。更一般的任务margin与总体风险推导见 `PAPER_CORE_AND_ORIGINALITY.md` 的F-PAPER-002—004。

## 结论与资格

本章给出的方向是一个可验证的有限集合成像目标，不宣称新概率定理。当前机制是固定/动态弱偏置提案、保持弱直线而允许signed幅度变化的合法修正、随机全空间方向对照、固定候选/费用预算，以及GT不可见的去重和资格规则。平滑tether与全秩混合支持为 `NOT_IMPLEMENTED` 备选。它与严格posterior任务并行，而不是等待后者认证才开始。

本章维护没有执行新物理调用或数值测试；E5的parent数、K、半径与对照已经按最终协议冻结，不能再当作可随结果调整的建议。覆盖、范围和费用的实际结论必须取自正式封存后的E5证据；本文不据配置或代数推导宣称实验已经完成或成功。


---

# 07 可能材料云：新目标、算法与证据流

这是新 possibility-volume namespace 的源码流程说明，写作时没有读取新 truth 或运行新物理。它不改变03章保留的原 E2/E3与REV1。最终协议已由 root批准，分别是 config/possibility_volume_protocol.json 与 config/possibility_statistics_frozen.json。生成、在线运行、封存及单独评分仍由 root负责；本章不根据批准状态宣称实验已完成。

## 目标变化必须明示

原探索主要使用相对锚点 DeltaPhi≤2与DeltaPrior≤10评估局部候选。新任务的主目标是实际材料可能性覆盖：候选满足原绝对95%noise ball、原decoder合法性、指定有限材料半径和signed progress。旧相对门槛保留为secondary fields，不参与新primary reject。

\[
\Phi_y(x)=\|F(\chi(x))-y\|_F^2/\sigma^2
\le\tfrac12\chi^2_{512,.95},
\quad |\|\mathcal R_m(\chi(x)-\chi(x_0))\|-\rho|\le.02\rho,
\]
\[
u_0^T\mathcal R_m(\chi(x)-\chi(x_0))\ge.5\rho.
\]

校正器的建议步仍使用原MAP二次模型
\[
m_k(d)=\tfrac12\|r_k+J_kd\|^2+\tfrac12\|x_k+d\|^2.
\]
但每个proposal的最终选择是已访问、合法且有限半径/progress合格状态中的最小真实dataPhi。模型目标和选择目标不同，必须同时记录。原raw proposal始终保留；若没有更小dataPhi的合格访问点，不强制用修正点替代raw。最小dataPhi仍不在noise ball时，点保存但primary_qualified=false。

## 固定8维弱空间、59维允许family和更强随机对照

用03章相同的material-effective阈值、prior-aware广义谱，取最弱8维，构造Gaussian正交frame \(Q_w\in\mathbb R^{66\times8}\)，其补空间 \(Q_s\in\mathbb R^{66\times58}\)。初始方向由8个标准Gaussian系数作用于材料单位weak lift \(TW_8\)，B/C完全共用。

固定方法允许
\[
x-x_0=\alpha v_0+Q_s s,\qquad
(\alpha,s)\in\mathbb R^{1+58}.
\]
其GN允许空间由 \([Q_s,v_0/\|v_0\|]\) 张成，再限制到材料有效空间并施加原signed-progress tangent row；通常维数至多58。随后使用原R7 radial decoder root匹配有限材料半径。

这与REV1“只在Qs中做GN，再由radial root调整幅度”不同：新runner的GN本身包含所选weak line的幅度方向。也不能称为完整66维GN。Gaussian Qs与材料方向并非材料度量正交，所以signed material progress不自动保证native weak coefficient为正；冻结合同允许signed amplitude，历史名称中的RAY应理解为weak line，实际系数/符号另存。

动态方法每次有已付费当前J时刷新rank8 weak frame，取当前位移在该frame中的weak projection作为本步weak line；不声称跨frame平行运输。随机方法用独立nativeGaussian66 antithetic方向与原ordinary full-effective66 signed-row校正器，是更强几何class，不能人为削弱成59维后再声称公平。

## K16不是每半径16

每个arm/parent总共16个非anchor proposal：两个antithetic pairs乘四个半径 \(.02,.05,.10,.15\)，每个pair有+/-两个方向。顺序为pair0的全部半径+/-，再pair1的全部半径+/-。因此prefix K=1/4/8/16是确定且公开的。

每个半径独立从同一anchor出发，没有从较小半径的已选点warm-start。latent距离上限固定4，不自动放宽；返回RADIUS_BRACKET_NOT_FOUND_WITHIN_LATENT_CAP只表示给定射线搜索和cap未找到bracket，不证明非单调decoder距离下的全局不可行性。

共享anchor加16个slot至多17个点，可能重复、可能不合格、可能slot unavailable。PRIMARY_CLOUD.npz存放所有selected states及资格flags，并非“仅合格点文件”，更不是posterior samples。anchor也只有原decoder合法且在绝对noise ball内才进入qualified cloud。

## 真物理与成本分开写

固定host exact ROI elimination对所有原P1非均匀ROI材料有效：
\[
L_R=I_{256}-\operatorname{diag}(\delta)C,\quad
j_R=L_R^{-1}\operatorname{diag}(\delta)E_R^h,\quad
F=F_h+Oj_R.
\]
\[
E_R=E_R^h+Cj_R,\qquad
J_x=\mathcal W\{O L_R^{-1}\operatorname{diag}(E_R)D_x\chi_R\}.
\]

一次forward实际支付8个RHS256和一次256LU；full66 J用同一LU的32个reduced adjoint RHS，不是528个forward directional solves。只保留一个live LU。这个精确块消元是共同的普通physics backend，不是分配给fixed Q的独占科学优势。

共享SETUP缓存已经支付264个full576 RHS；此前完整stable ROI timing study共2184个full576 RHS。两者分别注明历史身份，不在新campaign再次收费。当前协议没有新full576审计或fallback。

上界是1536个proposal，各最多6次forward和4个J-bearing states：
\[
1536(6\cdot8+4\cdot32)=270336\ {\rm RHS256}.
\]
加32个cap64 anchor的81920 RHS256及生成256 RHS256，总352512 RHS256；最坏11296次LU256。硬cap为400000 RHS256、12000 LU256、7200 controlled seconds、8GB。不同维度RHS不相互冒充。拒绝、失败、cache恢复、源几何/SVD、decoder radius search与I/O都保留成本；sourcegeometry/SVD通常不增加RHS但增加CPU/wall。

## 数据流和seal边界

1. 独立generator使用20261030_R9_POSSIBILITY_VOLUME_V1 SeedSequence命名空间，prior与noise分角色。所有truth、generator预测及generator求解attempt arrays写private output；online文件只有observation/sigma/parent。
2. online anchor从zeros66开始，原R7 direct-probit TRF设置不变；online代码不导入generator，不打开truth。
3. 每个parent记录48行（3arm×16slot），包括不可用与失败slot；保存raw、每个实际visited trial、接受标志、最终selected数组、J是否真正计算、timing window/arm order和所有计费。
4. 每个parent的seal同时包含parents/<id>和jobs/<id>，因为VISITED.x_path指向后者；绝不加入private generator目录。
5. 只有32个parent全部完成并产生PREDICTIONS_SEALED.json后，下游统计/离线truth评分才能读相应数组。运行器不自动把候选云变成概率参考或posterior结论。

源冻结记录BLAKE2b与文件stat，包含实际engine/generator、原P1源码/target、共享SETUP缓存及root批准的statistics frozen config。没有SHA要求，也没有额外微型数值测试来替代最终32-parent语义验证。


---

# 08 · 径向—角向覆盖、随机化目标与任务方向

本章把 `POST_E5_THEORY_DECISION.md` 中F-POST-001—003正式纳入持续维护的公式集。保持原编号，不把标准推导重新命名为新定理。F-POST-001是材料空间中的精确几何关系；引用E5的解释属于 `POSTHOC_EXPLANATORY_ONLY`。F-POST-002与003是已知线性代数/随机优化构件的推导，其作为新候选生成机制的状态均为 **`NOT_IMPLEMENTED / NOT_RUN`**。

E5实际实施与评分仍以 [最终运行协议](../config/possibility_volume_protocol.json) 和 [冻结统计协议](../config/possibility_statistics_frozen.json) 为准。添加公式不改变候选、预算、主统计门或已有结果。本章新增物理RHS、训练、数值测试均为0。

## F-POST-001 · 材料邻域覆盖的精确径向—角向分解

### 空间、假设与推导

原材料 $m\in\mathbb C^{576}$，采用全网格复材料RMS距离。令

$$z(m)=\frac1{24}\begin{bmatrix}\operatorname{Re}m\\\operatorname{Im}m\end{bmatrix}\in\mathbb R^{1152},\qquad d_M(m,m')=\|z(m)-z(m')\|_2.$$

固定anchor $m_a$、真实材料 $m_*$ 和一个候选 $m_j$。定义 $g=z(m_*)-z(m_a)$、$u=z(m_j)-z(m_a)$，以及 $R=\|g\|$、$\rho=\|u\|$。若 $R\rho>0$，令 $\cos\theta=g^Tu/(R\rho)$。展开平方即得

$$\boxed{d_M(m_*,m_j)^2=\|g-u\|^2=R^2+\rho^2-2R\rho\cos\theta.}$$

这在材料嵌入空间中精确成立，不要求decoder线性、弱散射、固定Jacobian或GT方向各向同性。若 $R=0$，误差就是 $\rho$；若 $\rho=0$，候选就是anchor，误差就是 $R$，两种情形均无需定义角度。

当 $R\rho>0$ 时，覆盖 $d_M(m_*,m_j)\le\varepsilon$ 等价于

$$\cos\theta\ge c(R,\rho,\varepsilon):=\frac{R^2+\rho^2-\varepsilon^2}{2R\rho}.$$

由于 $\cos\theta\le1$，必要条件为 $|R-\rho|\le\varepsilon$。当 $c>1$ 时任意方向都不可能命中；当 $c\le-1$ 时任意方向都可命中；$c=1$时只有严格同向的退化角帽。**径向可达只是必要条件，不能代替角向覆盖。**

与anchor比较，候选改善要求

$$d_M(m_*,m_j)^2<R^2\quad\Longleftrightarrow\quad \rho<2R\cos\theta\qquad(\rho>0).$$

因此沿一个固定正方向，误差关于长度为二次函数 $R^2+\rho^2-2R\rho\cos\theta$，其在 $\rho\ge0$ 上的最优长度为

$$\rho_* =\max(0,R\cos\theta).$$

若允许沿同一直线取signed长度，则最优signed系数是 $R\cos\theta$；E5的weak-line系数允许signed变化，但其材料方向也随decoder/corrector变化，不能把这个一维oracle直接当作当前算法。$R$与$\theta$都依赖未知GT，以上最优长度是解释工具，不是部署规则。

### 线性子空间与有限K的条件解释

仅在另行指定的固定线性材料空间 $\mathcal S$ 内，分解 $g=g_\parallel+g_\perp$，其中 $g_\parallel=P_\mathcal Sg$。设候选增量为 $u=\rho v$、$v\in\mathcal S$、$\|v\|=1$。由正交性，

$$\|g-\rho v\|^2=\|g_\perp\|^2+\|g_\parallel\|^2+\rho^2-2\rho\|g_\parallel\|\cos\theta_\parallel.$$

若 $\|g_\perp\|>\varepsilon$，无限增加该空间内的点也不能覆盖GT。若 $\rho\|g_\parallel\|>0$，覆盖对应一个角帽，其阈值为

$$c=\frac{\|g_\perp\|^2+\|g_\parallel\|^2+\rho^2-\varepsilon^2}{2\rho\|g_\parallel\|}.$$

若进一步假设 $d=\dim\mathcal S\ge2$，提案方向独立且均匀分布于 $S^{d-1}$，那么对固定GT方向，$t=\cos\theta_\parallel$ 的密度为

$$f_d(t)=\frac{\Gamma(d/2)}{\sqrt\pi\Gamma((d-1)/2)}(1-t^2)^{(d-3)/2},\qquad -1<t<1.$$

这是单位球面按第一坐标切片得到的标准密度；其归一化由Beta积分给出。于是对 $-1<c<1$，单次命中概率与独立K次至少命中一次的概率分别为

$$p_d(c)=\int_c^1 f_d(t)dt,\qquad P_K=1-[1-p_d(c)]^K\le Kp_d(c).$$

这个概率解释假设的是**提案方向均匀**，不需要GT方向随机或均匀。不同半径对应不同 $p_k$ 时，独立性下为 $1-\prod_k(1-p_k)$；没有独立性时最多保留适用的union bound，不能继续使用乘积。antithetic配对、动态frame、非线性decoder、corrector以及资格筛选都可能破坏均匀性或独立性。$d=1$只有正负两个方向，需按离散分布单独处理。

**不把此条件模型冒充E5概率。** 不能直接代入d=8或66推断本轮命中率，也不能以“高维”作为全部失败原因。真实候选不保证留在初始线性空间；当前corrector也会改变方向。

### E5的posthoc含义与边界

既有结果显示：B的新增主ε=.075命中来自.05半径，.10/.15没有新增主命中；保存同一批路径全部visited噪声合格点后仍为6/32。较大半径不是自动收益，与上述精确几何关系一致；没有据此声称已经逐候选计算了GT角度。

实际候选仿射空间的GT投影是另一个事后容量诊断。若 $U$的列为实际候选相对anchor的增量，则

$$\|(I-UU^\dagger)g\|\le\min_j\|g-Ue_j\|.$$

这是因为每个实际候选均位于该仿射空间。B有19/32个parent的这一oracle剩余RMS仍高于.075，另7个parent的仿射oracle可进入.075但有限候选未命中；前者揭示实际方向张成的遗漏，后者揭示点/系数放置的遗漏。数值秩依据现有诊断阈值；仿射投影可能使用负系数、违反decoder或数据门，并且使用GT，不能作为新成像结果。事实来源：[posthoc汇总](../results/possibility_volume/posthoc_geometry/SUMMARY.json)、[parent诊断](../results/possibility_volume/posthoc_geometry/PARENT_DIAGNOSTICS.csv)、[封存后理论裁决](../POST_E5_THEORY_DECISION.md)。

**资格：** 精确材料几何身份 + 明示假设的球面概率解释；E5连接为 `POSTHOC_EXPLANATORY_ONLY`。平方展开、正交投影与球面角帽公式均为标准数学，不主张原创定理。

## F-POST-002 · 随机化data/prior目标产生完整局部方向

### 原目标与局部模型

采用Gaussianized实latent $x\in\mathbb R^{66}$，原prior为标准Gaussian。白化数据residual $r(x)\in\mathbb R^{512}$，原目标为

$$\Psi(x)=\tfrac12\|r(x)\|^2+\tfrac12\|x\|^2.$$

在一个固定anchor $x_a$，记 $r_a=r(x_a)$、$J=\partial r/\partial x|_{x_a}\in\mathbb R^{512\times66}$。固定这些量并作局部线性化 $r(x_a+p)\approx r_a+Jp$。令

$$H=J^TJ+I_{66}\succ0,\qquad g=J^Tr_a+x_a.$$

取彼此独立的 $\eta\sim\mathcal N(0,I_{512})$ 与 $\xi\sim\mathcal N(0,I_{66})$，分别随机化白化data目标与prior中心。局部随机目标为

$$q_{\eta,\xi}(p)=\tfrac12\|r_a+Jp-\eta\|^2+\tfrac12\|x_a+p-\xi\|^2.$$

对p求导并置零，得到唯一最优位移

$$\boxed{p=H^{-1}(J^T\eta+\xi-g).}$$

均值为 $-H^{-1}g$；只有anchor恰为该局部目标驻点、$g=0$时才是零均值。这个项不能因anchor被称作MAP而未经检查省略。

独立性给出

$$\operatorname{Cov}(J^T\eta+\xi)=J^TJ+I=H,$$

因而

$$\boxed{\operatorname{Cov}(p)=H^{-1}HH^{-1}=H^{-1}.}$$

它在全部66个latent方向上正定，并按数据曲率缩放，而不预先把weak8之外的方向置零。只有随机化起点、却把同一严格凸目标优化到底，仍会返回同一点；这里改变的是优化目标的随机线性项。

### 一般线性Gaussian问题中的精确分布资格

若完整forward本身线性，$y=Ax+\epsilon_{obs}$、$\epsilon_{obs}\sim\mathcal N(0,\Gamma)$，prior为 $\mathcal N(m,C_0)$，且 $\Gamma,C_0\succ0$，独立抽 $\eta_y\sim\mathcal N(0,\Gamma)$、$\xi_x\sim\mathcal N(0,C_0)$，解

$$x_* =\arg\min_x\left[\tfrac12\|Ax-y-\eta_y\|_{\Gamma^{-1}}^2+\tfrac12\|x-m-\xi_x\|_{C_0^{-1}}^2\right].$$

令 $H_0=A^T\Gamma^{-1}A+C_0^{-1}$。由正规方程，

$$x_*=H_0^{-1}\{A^T\Gamma^{-1}(y+\eta_y)+C_0^{-1}(m+\xi_x)\}.$$

它是Gaussian，均值为 $H_0^{-1}(A^T\Gamma^{-1}y+C_0^{-1}m)$，协方差为 $H_0^{-1}$，与该线性Gaussian posterior一致。这项精确资格来自完整线性模型和精确优化，不自动延伸到非线性Maxwell、截断域或受限corrector。

### 实现位置、成本与失败条件

固定anchor局部采样可以复用H的一次factorization；但若沿每个随机目标继续非线性优化，r、J、H都随状态变化，需要支付真实物理与导数费用。固定host/ROI算子仍可共享，各候选的变化材料系统不能免费冻结。使用近似解、局部极小点、硬几何约束、拒绝或去重，均会改变输出分布。

作为可能性云方法时，候选最终必须对**原观测y**检查数据资格；不能用其随机化后的观测目标放行。它可能扩展当前局部探索，但不保证数据合格、高GT覆盖或正确posterior。随机化目标、其有限非线性优化与过滤分布的关系，需要独立分析与实验。

**当前状态：`NOT_IMPLEMENTED / NOT_RUN`。** E5未运行此机制，不把局部Gaussian公式当作已发生的Full66采样。提出它是为了给下一次机制裁决一个成熟强基线，而不是替当前weak8失败补造结果，也不授权追加实验。

### 先例与原创性边界

优化产生随机候选的直接先例为 [Bardsley等，RTO，2014](https://doi.org/10.1137/140964023) 与 [Wang、Bui-Thanh、Ghattas，rMAP，2018](https://doi.org/10.1137/16M1060625)。rMAP作者明确说明特定随机data/prior中心下与RML的关系；RTO作者明确区分线性直接posterior样本与非线性需密度/接受校正的情形。本轮新核验范围是出版社/作者摘要，不冒称逐式审核其全部非线性算法。上述局部和线性Gaussian公式为独立标准推导，不把它等同于完整非线性RTO实现，也不主张原创。证据层级与检索429记录见 `PAPER_CORE_AND_ORIGINALITY.md`。

## F-POST-003 · 在曲率预算下选择任务敏感方向

### 假设与推导

令任务 $T(m)\in\mathbb R^s$，$W\succeq0$为预先确定的任务单位/重要性metric，$m=\mathcal G(x)$。假设复材料decoder作为实映射以及任务在anchor附近可微。定义

$$B=\left.\frac{\partial[W^{1/2}T(\mathcal G(x))]}{\partial x}\right|_{x_a}\in\mathbb R^{s\times66},\qquad H=J^TJ+\lambda I\succ0,\quad\lambda>0.$$

当λ=1且J按原噪声白化时，H是当前Gaussian prior对应的GN曲率；其他λ是另行声明的正则尺度，不能静默当作原prior。给定 $b\ge0$，考虑局部方向问题

$$\max_{p^THp\le b^2}\ \|Bp\|^2.$$

令 $u=H^{1/2}p$。约束变成 $\|u\|\le b$，目标变成

$$u^T\underbrace{H^{-1/2}B^TBH^{-1/2}}_{C_T\succeq0}u.$$

由对称矩阵Rayleigh商，最优值为

$$\boxed{b^2\lambda_{max}(C_T).}$$

若最大特征向量 $v_T$ 取单位范数，则一个最优位移为 $p_*=bH^{-1/2}v_T$。将 $C_Tv_T=\mu v_T$ 换回原坐标，得到等价广义特征问题

$$\boxed{B^TBp=\mu Hp.}$$

当B为零时最优值为零，不能因此声称任务在有限非线性区域恒定；当最大特征值重复时，最大特征子空间内方向不唯一，不应给单向量附加无根据的物理身份。

### 固定弱空间中的受限版本

若仅允许 $p=Qa$，其中 $Q\in\mathbb R^{66\times k}$满列秩，则 $Q^THQ\succ0$，问题变为

$$\max_{a^TQ^THQa\le b^2}a^TQ^TB^TBQa,$$

对应 $Q^TB^TBQ\,a=\mu Q^THQ\,a$。Q不必Euclidean正交，但必须保留真实 $Q^THQ$，不能因为用了材料归一化基便把该metric替换为单位阵。限制在weak8内只比较其中任务方向，不能解决其空间外的覆盖遗漏。

### “预算”与实际物理资格不是一回事

$p^THp$是明确的局部二次曲率预算。若anchor梯度 $g\ne0$，原二次目标变化为 $g^Tp+\tfrac12p^THp$，而不是仅有二次项；即使g=0，有限步的非线性余项仍需处理。因此该约束不能自动解释为真实likelihood变化上限，也不等于材料半径约束或绝对噪声资格。

任务重要性与数据不可辨识性是不同指标：某方向数据很弱但任务几乎不变，另一个方向可能在相同曲率代价下产生明显任务变化。该广义谱给出这种区别的局部量化。W必须由任务尺度/容差事先确定，不得看到GT后调整。

H、B通常随状态变化；若任务仅是材料decoder的函数，构造B主要需要任务/decoder导数，但构造当前J仍要付物理成本。hard threshold在阈值处不可微，不可直接套用；soft occupancy可以求导，但局部导数不保证硬事件正确。任务数小使 $\operatorname{rank}(B^TB)\le s$，这不意味着剩余材料方向可删除，它们仍可能影响GT邻域、其他任务或非线性行为。

最大task spread也不保证围绕GT展开；它可能仅把错误候选拉得更远。因此该方向准则必须与真实物理资格、集合范围、有限候选预算及独立GT/task覆盖一起裁决。F-VOL-015说明为什么task dispersion与near-GT权重还会影响最终加权动作。

### 资格与已知背景

**当前状态：`NOT_IMPLEMENTED / NOT_RUN`。** E5计算了任务表现，没有实施这里的任务导向proposal或槽位分配。Rayleigh商/广义特征值是标准数学；目标导向降维也不是本轮新概念。已核验的 [Cui等2014 LIS](https://arxiv.org/html/1403.4680v2) 与 [Zahm等 certified reduction](https://arxiv.org/html/1807.03712v4)构成prior/likelihood降维背景，但不能把它们说成已经逐式提出本章具体任务目标，也不能因只核验这几篇就声称本构造没有先例。

若未来把任务重复性、原物理误差界与这个局部方向准则组合成独立算法，需要补做最接近的目标导向降维/实验设计先例核验，并与成熟随机目标方法和当前精确ROI强基线比较。这里仅保留推导与成立条件，不进行新颖性认证，不启动新的物理实验。

## 使用与状态规则

F-POST-001说明现有材料点为何可能“大而不覆盖”；F-POST-002说明完整局部方向可以通过已知随机目标机制被激发；F-POST-003说明任务方向与单纯数据弱方向不同。三个公式分别属于几何解释、未运行的随机目标机制和未运行的任务导向机制，不应合并成一个已经验证的新算法。

实际E5结果只能由封存后的证据支持；局部Gaussian、球面独立均匀提案、任务可微性与曲率预算等条件必须随公式一起引用。新增本章不改变任何历史门或已冻结配置。

---

# 09 · 双边传递、形状—幅度分离与按任务控制的物理计算

本章沿着“固定的两端—变化的中间材料”展开，不把标准恒等式当成新定理。R10 的实测问题是：这些公式能否减少完整非线性 Maxwell 的重复费用，或者生成更有用的材料世界？公式成立与算法胜出分别报告。

统一尺寸：ROI 阶数 $q=256$，接收端 $n_r=32$，源端 $n_s=8$。固定宿主消元后的 $C\in\mathbb C^{q\times q}$、$O\in\mathbb C^{n_r\times q}$、$E_R\in\mathbb C^{q\times n_s}$、$F_h$ 均固定；$D$ 是相对宿主的材料增量对角矩阵。时间因子、白化和原先验沿用 F-BASE-001—003。本章只对冻结离散系统作结论。

## F-BI-001 · 双边非线性传递

由 $h=D(E_R+Ch)$ 消去电流得到

$$\mathcal T(D)=(I-DC)^{-1}D=D(I-CD)^{-1},\qquad F(D)=F_h+O\mathcal T(D)E_R.$$

推导使用 $(I-AB)^{-1}A=A(I-BA)^{-1}$，只要求对应 resolvent 可逆，不要求 $D$ 可逆。接收端和输入端的维数小，不意味着中间状态在所有材料下都可以无损缩小。**资格：标准离散 T-operator 恒等式。** 实现：R9 精确 ROI solver；R10 `src/transfer_campaign.py`。

## F-BI-002 · 有限材料差的双边条件

$$\mathcal T(D_2)-\mathcal T(D_1)=(I-D_2C)^{-1}(D_2-D_1)(I-CD_1)^{-1}.$$

在等式左乘 $I-D_2C$、右乘 $I-CD_1$ 后，两侧均化为 $D_2-D_1$，给出直接证明。因此等观测要求

$$O(I-D_2C)^{-1}\Delta D(I-CD_1)^{-1}E_R=0.$$

这是有限非线性响应等价的必要且充分条件，但未知 $D_2$ 仍在左侧 resolvent 中；它没有自动变成固定线性核约束。材料必须仍属于原 decoder 支持。**计算含义：** 允许构造有限响应和误差分析；不能据此宣布找到 exact fiber。

## F-BI-003 · 固定形状的精确加权实现

原 P1 材料增量 $D=\delta A(z)$，其中 $A=\operatorname{diag}(a_i)$，$0\le a_i\le1$，$\delta\in[.8,1.8]+i[.03,.18]$。令

$$P=A^{1/2},\quad B=PCP,\quad L=OP,\quad R=PE_R.$$

则

$$F(z,\delta)-F_h=\delta L(I-\delta B)^{-1}R.$$

证明：$\mathcal T(\delta P^2)=\delta P(I-\delta PCP)^{-1}P$，可直接左右乘原方程验证，完全不需要 $P^{-1}$。所以零占据率 cell 不造成除零。**资格：精确离散恒等式，保留所有多次散射。** 形状改变后 $B,L,R$ 也改变；仅复幅度查询可以共用同一分解。

## F-BI-004 · 传递矩、可达与可观测状态

在 $|\delta|\rho(B)<1$ 的收敛域内，

$$F-F_h=\sum_{j\ge0}\delta^{j+1}LB^jR.$$

系数 $LB^jR$ 是输入输出矩。有限维实现的可达空间为 $\mathcal R=\operatorname{span}\{R,BR,\ldots,B^{q-1}R\}$；不可观测空间为 $\mathcal N=\cap_{j=0}^{q-1}\ker(LB^j)$。最小实现阶数是 $\dim\mathcal R-\dim(\mathcal R\cap\mathcal N)$，也等于充分长块 Hankel 矩阵的秩。这里的交集是理论定义，算法可用 QR/SVD 的投影实现，不要求显式集合求交。

解析延拓可将同一有理响应用于级数收敛域外的非极点位置；不能把级数截断在该域外当成有效近似。R10 的 64 阶结果只是被测试的压缩阶数，未认证最小阶数。**来源：标准有限维实现与 Krylov 理论；最近材料模态先例见 R10 查新表。**

## F-BI-005 · 固定宿主的双边材料 Gram

宿主处一阶变化为 $dF=O\,dD\,E_R$，第 $i$ 个材料列 $J_{(r,s),i}=O_{ri}(E_R)_{is}$，故

$$J^*J=(O^*O)\circ\overline{E_RE_R^*}.$$

证明逐项求和：$\sum_{r,s}\overline{O_{ri}E_{is}}O_{rj}E_{js}=(O^*O)_{ij}\sum_s\overline{E_{is}}E_{js}$。这比仅 receiver SVD 同时保留了输入照明，但仍是宿主一阶信息。若 $D_x=\partial d/\partial x\in\mathbb C^{q\times66}$，真实高斯 latent 上的白化 Gram 为

$$H_x=\frac2{\sigma^2}\operatorname{Re}\{D_x^*(J^*J)D_x\}.$$

**条件：** 参数是实坐标，材料按实 Hilbert 度量；复 Gram 不得直接当成实 Hessian。R10 的谱更贴近 Full-J 并未转化为更好的 GT 覆盖。

## F-BI-006 · 完整非线性灵敏度的同型结构

令 $L_D=O(I-DC)^{-1}$，$R_D=(I-CD)^{-1}E_R$。由 F-BI-002 的极限有

$$dF=L_D\,dD\,R_D,\qquad J_D^*J_D=(L_D^*L_D)\circ\overline{R_DR_D^*}.$$

固定两端宿主 Gram 与当前 Full-J 的差别恰在两侧内部反馈。它解释何时 Born 双边坐标失配，但不允许免费得到状态依赖的两侧场。**资格：标准 distorted-Born 导数结构。**

## F-BI-007 · 共同状态压缩的条件性障碍

考虑允许所有独立对角 $D$ 的材料族。若输入矩阵 $E_R$ 每行非零，则 $\operatorname{span}_{D}\{DE_R\}=\mathbb C^q$：用只保留第 $i$ 个对角元的矩阵即可产生第 $i$ 个坐标方向。因而不存在一个小于 $q$ 的固定状态子空间包含所有这些注入。

这只是充分反例：它限制所有材料通用的精确状态压缩，不否定固定形状的压缩，也不直接适用于受限的 P1 Gaussian64 decoder 像。更小的输入输出商仍要另审可观测性，不能由注入满张成单独排除。

## F-BI-008 · 交换子不是端口误差证书

$$[A,C]_{ij}=(a_i-a_j)C_{ij},\qquad \|[A,C]\|_F^2=\sum_{ij}|a_i-a_j|^2|C_{ij}|^2.$$

它衡量交换性；$[A,C]=0$ 仅说明 $C$ 不混合 $A$ 的不同特征值子空间，不保证共同对角化或低秩（例如 $A=I$、$C$ 为 Jordan 块）。大交换子可能主要落在端口不敏感方向。实际观测误差还需要 $L,R$、resolvent 稳定性和被舍弃状态的响应。**资格：恒等式加反推限制，不是新的普适误差界。**

## F-BI-009 · 固定形状 Schur 与摊销

复杂非 Hermitian $B=QTQ^*$，$Q$ 酉、$T$ 上三角。令 $\widehat L=LQ,\widehat R=Q^*R$，则

$$F-F_h=\delta\widehat L(I-\delta T)^{-1}\widehat R.$$

每节点只作三角求解。若 setup 为 $s$、直接节点费用 $c_d$、复用节点费用 $c_r<c_d$，同一计时边界下的盈亏点为 $M>s/(c_d-c_r)$。不可拿 warm 三角解与 cold 完整求解直接作全系统加速比。R10 精确256阶 Schur 的实测盈亏点21—58节点；压缩后的盈亏点另测。

## F-BI-010 · 复幅度解析导数

在 resolvent 非奇异处，因 $B$ 与其 resolvent 交换，

$$\frac{\partial F}{\partial\delta}=L(I-\delta B)^{-2}R,\qquad \frac{\partial^2F}{\partial\delta^2}=2LB(I-\delta B)^{-3}R.$$

这来自 $\partial[\delta(I-\delta B)^{-1}]=(I-\delta B)^{-2}$。对 $\delta=\alpha+i\beta$，实 Jacobian 两列为 $F_\delta,iF_\delta$。likelihood 包含复共轭绝对值，因而不是全纯函数；优化和概率积分仍在两个实属性维度上进行。

## F-BI-011 · 正确属性条件质量与 profile 的区别

原 Uniform2 经线性映射后密度为 $1/.15$，所以

$$\bar L(z)=\frac1{.15}\int_{.8}^{1.8}\int_{.03}^{.18}\exp[-\Phi(z,\alpha+i\beta)]\,d\beta\,d\alpha.$$

它不等于 $\exp[-\min_\delta\Phi]$。如果可交换导数和积分，

$$\nabla_z\log\bar L=-\mathbb E[\nabla_z\Phi\mid z,y].$$

保留 Gaussian64 先验后才得到形状 posterior。一个由 $y$ 选出的固定形状条件积分，不是 Full66 posterior。R10 对3个固定形状建立数值参考；没有为全部 $z$ 积分。

## F-BI-012 · 条件化能减少什么不确定性

对材料任务 $\varphi$，

$$\mathbb E[\varphi\mid y]=\mathbb E_{z\mid y}\mathbb E[\varphi\mid z,y],$$
$$\operatorname{Var}(\varphi\mid y)=\mathbb E_{z\mid y}\operatorname{Var}(\varphi\mid z,y)+\operatorname{Var}_{z\mid y}\mathbb E[\varphi\mid z,y].$$

精确处理属性只解决第一层条件计算；缺失形状分支仍影响外层质量与均值。Rao–Blackwell 化降低给定正确外层采样器的积分方差，不保证随机优化生成的形状具有正确 posterior 分布。

## F-BI-013 · 随机目标改变的条件协方差

在同一局部线性化，$H=J^TJ+I$，随机数据目标 $\eta\sim N(0,I)$、随机 prior 目标 $\xi\sim N(0,\Sigma_\xi)$ 独立，则 GN 更新中的随机部分满足

$$\operatorname{Cov}(\Delta x)=H^{-1}(J^TJ+\Sigma_\xi)H^{-1}.$$

完整随机目标 $\Sigma_\xi=I$ 时化为 $H^{-1}$。去掉两属性目标后 $\Sigma_\xi=I-P_{amp}$，差为 $-H^{-1}P_{amp}H^{-1}$，虽至多秩2却通过耦合影响全部坐标。修改目标协方差不等于保留原 posterior 的采样器。R10 的 bilateral trace-normalized covariance 是提案启发式，非新增物理信息。

## F-BI-014 · 非正规矩阵的损耗稳定性界

令 $K_\delta=\delta^{-1}I-B$、$\operatorname{Im}_H B=(B-B^*)/(2i)$，设 $\beta_B=\lambda_{min}(\operatorname{Im}_H B)$，则

$$-\operatorname{Im}_H K_\delta\succeq b_\delta I,\quad b_\delta=\frac{\operatorname{Im}\delta}{|\delta|^2}+\beta_B.$$

若 $b_\delta>0$，由 $|v^*K_\delta v|\ge b_\delta\|v\|^2$ 和 Cauchy–Schwarz 得

$$\sigma_{min}(K_\delta)\ge b_\delta,\qquad \|K_\delta^{-1}\|_2\le b_\delta^{-1}.$$

不要求 $B$ 正规，也没有假设 $\operatorname{Im}_HB\succeq0$。R10 的最小特征值实际略负，材料损耗仍使 $b_\delta$ 正。float64 特征值仅给数值资格；没有 outward rounding 就不能称严格机器证书。

## F-BI-015 · 全状态 residual 到 likelihood 误差

对近似解 $\tilde u$，$r=R-K_\delta\tilde u$，真实输出误差为 $LK_\delta^{-1}r$，所以

$$e:=\|W_\sigma(F-\tilde F)\|\le\frac{\sqrt2\|L\|_2\|r\|_F}{\sigma b_\delta}.$$

代码可用较保守但便宜的 $\|L\|_F$ 代替 $\|L\|_2$。若 $\tilde r_y=W_\sigma(\tilde F-y)$，

$$|\Phi-\tilde\Phi|\le\|\tilde r_y\|e+\tfrac12e^2.$$

这是从同一物理 residual 到任务判断的桥梁。R10 冻结门为 $e\le.01$ 且势误差上界 $\le.1$；不合格才 exact fallback。门只覆盖被评价的节点，不凭有限节点推断全先验支持。

## F-BI-016 · 可达 Krylov 压缩与不能省略的检查

取正交 $V_r$ 张成 $\operatorname{span}\{R,BR,\ldots,B^{p-1}R\}$ 的有效列，$H_r=V_r^*BV_r$，

$$\tilde F-F_h=LV_r(\delta^{-1}I-H_r)^{-1}V_r^*R.$$

若前 $p$ 个 Krylov 块完全包含在空间中，Galerkin 传递匹配相应低阶输入矩。这个局部展开性质不能独自保证所有合法复幅度的准确度，故实际算法仍计算完整 $r$，用 F-BI-014—015 决定是否复用或回退。压缩可减少内部解阶数，不减少原材料自由度；一次固定形状的 basis 不跨任意形状通用。

## F-BI-017 · 势误差到概率误差的条件

若同一可归一化先验支持上，$\operatorname{osc}(\Phi-\tilde\Phi)\le\Delta$，则标准 likelihood ratio 界给

$$\|\pi-\tilde\pi\|_{TV}\le\tanh(\Delta/4).$$

若只有 $|\Phi-\tilde\Phi|\le\varepsilon$，可用 $\Delta\le2\varepsilon$。有界实任务误差至多 $(\sup\varphi-\inf\varphi)TV$。必须同时处理积分截断、求积和遗漏支持；R10 节点级误差门尚不满足这个全支持前提。因此没有把节点通过率当作 posterior 认证。

## F-BI-018 · 被动性、互易性与有限观测的边界

对 $e^{-i\omega t}$ 标量 TM、实 $\mu$，$H=(\partial_yu,-\partial_xu)/(i\omega\mu)$，

$$\langle S\rangle=\frac{\operatorname{Im}(\bar u\nabla u)}{2\omega\mu},\quad \nabla\cdot\langle S\rangle=-\frac\omega2\operatorname{Im}\varepsilon|u|^2-\frac12\operatorname{Re}(\bar uJ).$$

若离散 $C^T=C$、$D$ 对角，$\mathcal T(D)^T=\mathcal T(D)$；这不要求矩形带标签 $OTE_R$ 自身对称。32个有限 receiver 不是完整边界功率积分。原 prior 已要求被动性，不能重复把它作为独立新观测。能量型稳定性用于误差控制有计算价值，但没有在本轮新增材料信息。

## 从公式到本轮决策

保留：F-BI-003/009/014—016 组合的固定形状响应复用，检查 setup 摊销与完整 residual。停止扩大：仅凭 F-BI-005 的宿主双边弱方向主张 GT 覆盖改善。保留为独立数值任务：F-BI-011 的属性积分。未启动：Full66 SDE、通用最小状态、跨所有形状的精确公共压缩。

完整推导、先例与实测来源：[R10 理论](../p3l_r10_bilateral_structure/BILATERAL_OPERATOR_THEORY.md)、[能量与误差](../p3l_r10_bilateral_structure/PHYSICS_ENERGY_AND_RECIPROCITY.md)、[10形状物理复用](../p3l_r10_bilateral_structure/REDUCED_RESPONSE_FRESH10.md)、[原创性审查](../p3l_r10_bilateral_structure/PRIOR_ART_AND_NOVELTY.md)。


---

# 10 · TSOM 如何作用于双边 Maxwell 传递

2026-10-30，回应用户关于 $O(I-DC)^{-1}DE_R$ 的追问。本章为解析推导与算法定位，新增 Maxwell RHS=0；不改写 R10 已冻结实验。结论是：**可以使用 SOM/TSOM，材料与 Green 不对易并不构成否决理由。** 要区分绕开内层完整正演的电流辅助反演，与降低每次完整正演成本的近似算子求解。

实质文献依据为作者公开的 [Xu、Zhong、Song、Chen、Ran，MR-FFT-TSOM，DOI 10.1109/TGRS.2014.2329032](https://renchengsong.github.io/papers/TGRS2329032.pdf)，§II，式(3)—(11)。该文明确重述接收 SVD、内部传播基的投影及数据/状态双残差。以下符号吸收了原文电流中的频率常数；原论文算法、这里的噪声修正、variable projection 和低秩正演必须分开归属。没有声称已逐式读到2009原篇。

## F-BI-019 · 从双边逆算子还原电流方程

令 $N=256$，$n_s=8$，$h\in\mathbb C^{N\times n_s}$ 为固定宿主消元后的 ROI 对比电流，$y=Y-F_h$。定义 $K_D=I-DC$。当 $K_D$ 可逆时：

$$h=K_D^{-1}DE_R\quad\Longleftrightarrow\quad h=D(E_R+Ch),\qquad F-F_h=Oh.$$

这是精确等价式。原完整正演只需解8个右端，并不显式构造256×256逆矩阵或整个 $\mathcal T(D)$。SOM/TSOM 的入口是将 $h$ 作为辅助未知量，以固定投影组织它，同时满足上述状态方程与数据方程。它不要求 $D$ 与 $C$ 对易，也不要求真实材料解集线性。

## F-BI-020 · 两次筛选是投影像，不是交集

设 $S=WO$，其中 $W$ 为已固定的逐源噪声白化；$S=U\Sigma V^*$。保留噪声下可靠的前 $\ell$ 个方向：

$$h_d=V_\ell\Sigma_\ell^{-1}U_\ell^*Wy,\qquad P_\ell=V_\ell V_\ell^*.$$

令 $V_{C,r}$ 为固定内部算子 $C$ 的所选右奇异向量，或使用明确声明的 Fourier 基 $V_{F,r}$。TSOM式表示为

$$Z=(I-P_\ell)V_{C,r},\qquad Q=\operatorname{orth}(Z),\qquad h=h_d+Q\alpha.$$

这里 $\operatorname{orth}$ 删除相关/零列，保留相同列空间。$\operatorname{range}Z$ 是投影后的像，一般不等于两空间交集；$ZZ^\dagger$ 才是该像的正交投影。截断接收强方向后，$OQ$ 通常非零，不能省略数据残差。原作者的 twofold 内部基/Fourier构造支持这个位置；它不是从材料 Jacobian 的 weak8 方向构造电流。

**固定缓存：** $Q,CQ,OQ,V_\ell,CV_\ell$；同一观测的 $h_d,Ch_d$。宿主消元后的 $C$ 一般不再具有完整均匀网格的卷积形式，Fourier作为基仍合法，但不能未经推导声称 $C$ 乘法自动具有FFT卷积加速。

## F-BI-021 · 共享材料与含噪数据的联合目标

原 P1 材料保持 $D=D(x)=\operatorname{diag}(\delta(x)a(x))$，$x\in\mathbb R^{66}$。不得逐源恢复不同材料，或用自由像素代替原 decoder。为不把含噪接收强系数硬锁在 $h_d$ 上，可加入其修正：

$$R=[V_\ell,Q],\qquad h=h_d+R\beta,\qquad \beta\in\mathbb C^{q\times n_s}.$$

采用明确的候选生成目标（归一化权重需单独冻结）：

$$\mathcal J_\lambda(x,\beta)=\|S(h_d+R\beta)-Wy\|_F^2+\lambda\|K_{D(x)}(h_d+R\beta)-D(x)E_R\|_F^2+\mathcal R(x),\quad\lambda>0.$$

对全部源共享同一个 $D(x)$；$\mathcal R$ 可为原先验惩罚或事先固定的随机化目标。没有 $N$ 阶逆出现在目标计算中：$K_DR=R-D(CR)$，因此固定的 $CR$ 可预计算，材料更新主要是对角行缩放。

有限 $\lambda$ 允许状态残差，所以这是优化辅助目标，不等于原 Maxwell likelihood。作为材料可能性生成器完全可以使用；候选需由原完整正演验证。不要求先解决精确 posterior 才能研究这个生成器。若把接收强系数固定不动，则应将其标为传统截断数据近似，而非原含噪模型的精确等价。

## F-BI-022 · 电流系数的小规模消元

对固定 $D$，将上一目标的两项堆叠：

$$T_\lambda(D)=\begin{bmatrix}SR\\\sqrt\lambda K_DR\end{bmatrix},\qquad B_\lambda(D)=\begin{bmatrix}Wy-Sh_d\\\sqrt\lambda(DE_R-K_Dh_d)\end{bmatrix}.$$

若 $R$ 满列秩且 $K_D$ 可逆，则 $T_\lambda$ 满列秩：其下块已经满列秩。于是

$$\beta_*(D)=T_\lambda(D)^\dagger B_\lambda(D),\qquad \min_\beta\mathcal J_\lambda=\|(I-T_\lambda T_\lambda^\dagger)B_\lambda\|_F^2+\mathcal R(x).$$

采用薄 QR/SVD 求解，不应为了书写方便在实现中显式形成 $(T^*T)^{-1}$。若 $q\ll N$，主要内层费用约为 $O((N+n_r)q^2+q^3)$，多个源共享同一QR；不再需要每一步 $N\times N$ LU。若 $q$ 很大或需要大量恢复迭代，收益可能消失。

令 $r=T_\lambda\beta_*-B_\lambda$，正交性 $T_\lambda^*r=0$ 给出

$$d\|r\|_F^2=2\operatorname{Re}\operatorname{tr}\{r^*[(dT_\lambda)\beta_*-dB_\lambda]\}.$$

因此约化目标求导不用显式求 $d\beta_*$，但它仍是上述有限罚目标的导数，不是原 Maxwell posterior 的梯度。这是标准 variable projection；是否比精确256阶Full-J更好需要真正的同质量成本比较。

## F-BI-023 · 若目的是直接压缩中间逆：低秩 Green 与 Woodbury

另一条路线固定近似 $C_r=AB^*$，$A,B\in\mathbb C^{N\times r}$（例如截断SVD因子）。令 $g=DE_R$。在 $I-B^*DA$ 可逆时：

$$\widetilde h=(I-DC_r)^{-1}g=g+DA(I-B^*DA)^{-1}B^*g.$$

这是标准 Woodbury 恒等式，对近似算子 $C_r$ 精确。对每个新材料仍需更新 $B^*DA$ 并解 $r$ 阶系统，但不用重新分解 $N$ 阶矩阵。若 $C=C_r$ 是实际精确低秩，则也是原离散物理的精确表达；一般截断情况下只是近似。

这与 TSOM 有共同的固定传播压缩哲学，但算法对象不同：TSOM主要减少反演中未知电流系数；此式直接压缩正演求解器。不能用“把 $C$ 分解一次”推出“所有 $D$ 下的逆都已分解一次”。同样不要求 $D$ 与 $C$ 对易。

## F-BI-024 · 完整物理缺陷给出投影的失效位置

对任意近似电流 $\widetilde h$ 定义真实缺陷

$$r_D=DE_R-K_D\widetilde h.$$

当 $K_D$ 可逆，精确误差恒为

$$h-\widetilde h=K_D^{-1}r_D,\qquad F-\widetilde F=OK_D^{-1}r_D.$$

在 F-BI-023 的低秩算子求解中，进一步有

$$r_D=D(C-C_r)\widetilde h.$$

这说明被删传播的实际影响同时受材料、内部场和resolvent放大控制；只看固定Green的奇异值尾不够。若有稳定性界 $\|K_D^{-1}\|\le b_D^{-1}$，则 $\|W(F-\widetilde F)\|_F\le\|WO\|_2\|r_D\|_F/b_D$；无可靠界时用完整校正/验证。

对纯TSOM，缺陷可能来自电流子空间不足、有限罚优化未收敛或材料与电流失配，不能全部叫“流形不稳定”。误差超门时先考虑扩大固定电流基或精确校正，不必默认更新材料弱方向。

## 本次判断与未执行边界

1. **成立：** 双边式天然允许TSOM式固定电流表示；完整非线性仍在状态方程中，没有使用Born替代。
2. **成立：** 固定基缓存+小系数QR有实际避开内层完整LU的路径；这比只说“可能做SVD”更具体。
3. **尚未实测：** 原66维decoder内的上述TSOM/variable-projection候选生成器是否改善材料覆盖及总费用。本次没有执行，不将R10双边Gram随机目标的失败当作它的失败。
4. **已有但不同：** R10固定形状的rank64传递压缩与条件幅度积分得到有范围的正面计算结果；它不是这里完整的TSOM反演实验。
5. **下一项有辨别力的实验：** 共享同一精确ROI验证器和材料decoder，对照普通Full-J随机目标与固定TSOM电流系数+材料优化，按合格材料覆盖、遗漏/拒绝及完整费用比较。先不扩大弱8射线，不把准确posterior作为启动门。

本章数学结论属于已知状态分裂、投影、最小二乘消元和Woodbury的具体组合。潜在论文增量只能来自原材料支持下更好的合法可能性生成或总费用；不能将这些标准公式本身写成新定理。

---

# 11 · 精确活动支撑、总场投影与全形状损耗界

R11-S3，2026-10-30。接续 F-BI-019—024。以下为离散 Maxwell 的解析推导、标准投影与条件矩阵界；不以机器精度重放证明其新颖性。计算收益与数字图像资格见 [R11 报告](../p3l_r11_stable_structure_digits/RESEARCH_REPORT_ZH.md)。R11 的全阶 LU 自设上限失守，实验总状态为 `FAIL_CAP_EXCEEDED`；这不改写以下恒等式，也不允许把数值报告升级为预算合规 PASS。

统一尺寸：$n=256$，$n_r=32$，$n_s=8$。$C\in\mathbb C^{n\times n}$、$O\in\mathbb C^{n_r\times n}$、$E_R\in\mathbb C^{n\times n_s}$、$F_h\in\mathbb C^{n_r\times n_s}$ 固定；$a\in[0,1]^n$、$A=\operatorname{diag}(a)$、$D=\delta A$ 改变。$\delta$ 在原被动物性矩形内，$\sigma$ 和 source-major 白化沿用 F-BASE-002。$^*$ 表示共轭转置；解析导数转置求解与 Hermitian 内积不可混同。

## F-SS-001 · 严格零支撑的精确电流消元

设 $\mathcal S=\{i:a_i>0\}$，$\mathcal I$ 为其补集。由 $h=D(E_R+Ch)$ 逐行得 $h_\mathcal I=0$，故

$$M_S=I_S-\delta A_SC_{SS},\quad h_S=M_S^{-1}\delta A_SE_{R,S},\quad F=F_h+O_Sh_S.$$

按 $S,I$ 排序，$I-DC$ 是 $\begin{bmatrix}M_S&-D_SC_{SI}\\0&I_I\end{bmatrix}$，因此它与 $M_S$ 的可逆性等价。只删**精确零**；小正数删减是另一项近似。全部源共享一次分解。空支持时 $F=F_h$。资格：标准精确块消元；实现 `agent_a/solver.py:Solver.evaluate`。

## F-SS-002 · 完整端口灵敏度的活动/非活动分块

定义 $T=O(I-DC)^{-1}$ 和 $Z=O_SM_S^{-1}$。对上述块三角矩阵求逆直接得

$$T_{:,S}=Z,\qquad T_{:,I}=O_I+ZD_SC_{SI},\qquad e=E_R+C_{:,S}h_S.$$

只需活动系统的 $n_s$ 个正演和 $n_r$ 个转置右端；不需要另做全 $n$ 阶 LU。即使材料有零格点，仍要保留全部 $n$ 个材料导数列。资格：标准逆矩阵分块；不能把活动电流状态数当材料参数数。

## F-SS-003 · 支撑出生的非零材料导数

由 $(I-DC)\,dh=(dD)e$，得到

$$dF=T(dD)e,\qquad \frac{\partial F_{rs}}{\partial a_i}=\delta T_{ri}e_{is},\qquad \frac{\partial F}{\partial\delta}=T A e.$$

当 $a_i=0$ 时，$h_i=0$，但 $\partial F/\partial a_i$ 一般非零。因此原生占比优化可以用这些列激活新支持。对原 P1 latent 则须再乘 $Da/Dz$；平台内该因子可为零，不能由本式推出 latent 梯度必能出生新笔画。实虚物性导数分别为 $F_\delta$ 与 $iF_\delta$，Gaussianized 属性还须乘正态 CDF 导数。资格：精确一阶导数，不是新测量信息。

## F-SS-004 · 不随支持改变的总场注入

从 $e=E_R+Ch$、$h=\delta Ae$ 得

$$M_e e=E_R,\quad M_e=I-\delta CA,\quad F=F_h+\delta OA e.$$

$E_R$ 始终固定，避免旧加权状态 $u=\delta\sqrt A e$ 的注入 $\sqrt A E_R$ 随形状改变。这里不包含对低维共同不变空间的任何保证。资格：精确状态变换；零 $a_i$ 无需除法。

## F-SS-005 · 固定总场 Galerkin 与真实缺陷

取固定正交基 $V\in\mathbb C^{n\times r}$，例如 $\operatorname{orth}[E_R,CA_0E_R,\ldots]$ 或 $\operatorname{orth}[E_R,CE_R,\ldots]$。每个当前材料解

$$(I_r-\delta V^*CAV)t=V^*E_R,\quad\widetilde e=Vt,\quad\widetilde F=F_h+\delta OAVt,$$

$$r_e=E_R-(I-\delta CA)\widetilde e.$$

固定基不等于固定 $V^*CAV$。每次形状变化仍需更新 $CAV$、$OAV$ 和投影；其成本必须计入。非正规 Galerkin 的稳定性不由基正交自动保证。资格：标准近似；实现 `agent_b/fixed_total_field.py`。

## F-SS-006 · 全部占比形状的 Hermitian 虚部下界

令 $H_C=(C-C^*)/(2i)$、$P=\sqrt A$、$B=PCP$，有 $\operatorname{Im}_H B=PH_CP$。记 $\lambda=\lambda_{\min}(H_C)$，$\beta_*=\min(0,\lambda)$。对任意 $v$，

$$v^*PH_CP v\ge\lambda\|Pv\|^2\ge\beta_*\|v\|^2,$$

最后一步在 $\lambda<0$ 时用 $\|Pv\|\le\|v\|$，在 $\lambda\ge0$ 时直接用非负性。因此 $\lambda_{\min}(\operatorname{Im}_H B)\ge\beta_*$。资格：条件标准矩阵序不等式；不能把负 $\lambda$ 随意抹为零。

## F-SS-007 · 被动物性矩形上的统一强制性常数

对 $\delta=x+iy$，本合同 $0<y_{min}\le y\le y_{max}<x_{min}\le x\le x_{max}$。$y/(x^2+y^2)$ 随 $x$ 递减、随 $y$ 递增，故

$$b_*={y_{min}\over x_{max}^2+y_{min}^2}+\beta_*.$$

当 $b_*>0$，令 $K_\delta=\delta^{-1}I-B$，则 $-\operatorname{Im}_H K_\delta\succeq b_*I$。Cauchy–Schwarz 给

$$b_*\|v\|^2\le|v^*K_\delta v|\le\|v\|\,\|K_\delta v\|\Longrightarrow\|K_\delta^{-1}\|_2\le b_*^{-1}.$$

若 $b_*\le0$，本充分界失效，不能改用经验正数。本合同有限浮点结果 $\lambda=-4.5739299047\times10^{-5}$，$b_*=0.009210948658>0$；不是向外舍入的连续 PDE 证书。

## F-SS-008 · 总场缺陷到端口误差的精确矩阵顺序

由 push-through 恒等式，$P K_\delta^{-1}P=\delta A(I-\delta CA)^{-1}$。因此

$$F-\widetilde F=OP K_\delta^{-1}P r_e.$$

这个式子在 $A$ 奇异时仍成立，不需要 $P^{-1}$。直接展开 Neumann 多项式并由可逆区域延拓，或用 $(I-AB)^{-1}A=A(I-BA)^{-1}$，均可得到。资格：离散精确恒等式，不需要 $C,A$ 对易。

## F-SS-009 · 全形状白化端口误差界

将前两式取范数并用 source-major 等距实化，有

$$\|W(F-\widetilde F)\|_2\le\epsilon_W={\sqrt2\over\sigma b_*}\|OP\|_2\|P r_e\|_F,$$

$$\|OP\|_2\le\|O\|_2.$$

当前 $\|OP\|_2$ 可用 $32\times32$ 行 Gram 的最大本征值求得；全统一 $\|O\|_2$ 更便宜但更保守。资格：在上述耗散假设下的矩阵范数界；实现是 float64 数值条件界。

## F-SS-010 · likelihood 势误差与回退门

$\Phi(F;y)=\frac12\|W(F-y)\|^2$。若 $\|W(F-\widetilde F)\|\le\epsilon_W$，则

$$|\Phi(F;y)-\Phi(\widetilde F;y)|\le\|W(\widetilde F-y)\|\epsilon_W+\tfrac12\epsilon_W^2=\epsilon_\Phi.$$

接受门为 $\epsilon_W\le0.01$ 且 $\epsilon_\Phi\le0.1$；不通过有限增阶后精确回退。实际配对误差只能做事后评估，不能反过来替代在线界。这是误差控制，不是 posterior 认证。

## F-SS-011 · 近似物理导数不自动继承响应误差门

对精确 $T$、$e$，$J_i=\delta T_{:,i}e_{i,:}$。若用近似 $\widetilde T,\widetilde e$，则

$$J_i-\widetilde J_i=\delta[(T-\widetilde T)_{:,i}e_{i,:}+\widetilde T_{:,i}(e-\widetilde e)_{i,:}].$$

只有响应残差门并不控制全部 Jacobian。R11 的 B/C 在需要 Full-J 时使用精确一致的 F/J；只对 forward trial 尝试降阶。因此不能宣称降阶同时省了32个导数右端。资格：精确误差分解；受控降阶导数尚未实现。

## F-SS-012 · 三类成本与图像收益不可互换

对于相同观测、参数化、正则、初值和迭代配置，报告

$$C_{system}=C_{load}+C_{basis}+\sum C_{trial}+\sum C_J+C_{gate}+C_{fallback}+C_{validation}+C_{optimizer}+C_{save}.$$

活动支撑减少的是物理状态阶数，不改变可用观测；固定总场与耗散界也不增加材料信息。故“相同材料结果、较低 $C_{system}$”是求解器收益；“同预算更低材料风险/更好拓扑”需另有实验。R11 已有 72 份真实反演，但 B/C 全部 forward trial 回退，图像没有新增改善；支撑的正面核结果不能掩盖全阶 LU 预算失败。资格：评价定义及可审计计算分解，不是物理定理。

### 本轮原生反演目标与IRLS的具体含义

N模型令 $x=(a,\xi)$，$a\in[0,1]^{256}$，$\delta(\xi)=0.8+\Phi_N(\xi_1)+i[0.03+0.15\Phi_N(\xi_2)]$。$L$为ROI相邻水平/竖直差分，冻结目标为

$$\mathcal J_N=\Phi(F;y)+0.4\sum_j[\sqrt{(La)_j^2+0.05^2}-0.05]+0.01\|a\|^2+0.05\mathbf1^Ta+\tfrac12\|\xi\|^2.$$

$w_j=\sqrt{(La)_j^2+0.05^2}$，TV梯度为 $0.4L^T[(La)/w]$，代码采用正定IRLS majorizer $0.4L^T\operatorname{diag}(1/w)L$，而非TV精确Hessian $0.4L^T\operatorname{diag}(0.05^2/w^3)L$。阻尼GN解线性步后对$a$裁到$[0,1]$并用原目标回溯；未读GT。该新原生正则不等于原P1 Gaussian64 prior。

P模型仍为原Gaussianized66目标 $\mathcal J_P=\Phi(F(G(x));y)+\tfrac12\|x\|^2$。这一点估计目标不是posterior均值，也不等于native Uniform2坐标下的点目标。N与P不能当成纯后端对照。

### 新N模型的尺度gauge与材料任务

当 $a$ 为自由占比而 $\delta$ 未知，对任意使两端均合法的 $c>0$，有

$$a'=a/c,\qquad\delta'=c\delta,\qquad\delta'\operatorname{diag}(a')=\delta A.$$

这是同一材料的参数化gauge，不是新的Maxwell材料世界。占比阈值事件 $\{a_i>0.5\}$ 一般不沿该gauge不变；复材料误差或固定 $\chi_i$ 阈值才是材料任务。这一结论仅针对本轮新N自由占比参数化，不能从中推断原Gaussian64 quintic image具有普遍全局重缩放gauge。

---

# 12 · 双边材料模态、原先验几何与固定点谱

R12-BF，2026-10-30。一个材料格点怎样影响观测，要同时问两件事：源能否照亮它，以及它激起的电流能否传到接收端。两端共同决定数据灵敏度；材料 decoder 和原先验再决定这种改变能否以合理代价实际实现。把同一个 Jacobian 反复分解不会增加第三种物理信息。

本章接续 F-OP-010、F-MAT 系列及 F-BI-001—024，补齐 R12 实际使用的原 Gaussianized66 拉回、材料局部商、固定弱坐标和完整固定点导数。以下九条是**标准微分与线性代数在冻结离散模型上的具体推导**，不主张新定理。公式成立、代码实现、有限实验支持和方法新颖性是四个独立层次。本次维护新增 Maxwell RHS=0、LU=0；没有为恒等式安排数值重放。

统一尺寸：ROI 材料格点 $q=256$，接收端 $n_r=32$，源 $n_s=8$，原先验实坐标 $p=66$；完整材料网格 $N=576$。固定 $C\in\mathbb C^{q\times q}$、$O\in\mathbb C^{n_r\times q}$、$E_R\in\mathbb C^{q\times n_s}$、宿主响应 $F_h$ 与噪声 $\sigma=0.00119211798514519$。$d\in\mathbb C^q$ 是相对宿主的复材料，$D=\operatorname{diag}(d)$；$D_x=\partial d/\partial x\in\mathbb C^{q\times p}$ 是 decoder 导数，二者不混用。$^*$ 为共轭转置，$^T$ 为实转置或注明的普通转置。所有优化坐标均为实数。

## F-BF-001 · 当前总场与接收传播共同决定精确材料导数

电流和数据满足

$$h=D(E_R+Ch),\qquad F=F_h+Oh.$$

令 $M=I-DC$、$N_c=I-CD$，假设当前状态可逆；由 $\det(I-AB)=\det(I-BA)$，两种可逆性等价。定义

$$R_D=E_R+Ch=N_c^{-1}E_R\in\mathbb C^{q\times n_s},\qquad L_D=OM^{-1}\in\mathbb C^{n_r\times q}.$$

对电流方程求微分，注意材料与场都改变：

$$dh=(dD)(E_R+Ch)+DC\,dh,$$
$$\boxed{dF=L_D(dD)R_D.}$$

因此同一格点 $i$ 的材料改变由接收列 $L_{D,:,i}$ 与照明行 $R_{D,i,:}$ 相乘传到全部端口。无需 $C,D$ 对易，也无需 $D$ 可逆。这里 $dD$ 表示材料矩阵的微分；若使用占比 $a_i$ 且 $d_i=\delta a_i$，其列导数必须再乘 $\delta$。已含该因子的占比 Jacobian 不能再次乘它。

**假设与资格：** 冻结有限维状态方程、同一材料下全部源共享 $D$、当前 resolvent 可逆；这是精确一阶导数，不是有限步等响应曲线定理。宿主 Born 仅取 $D=0$ 时的 $L_0=O,R_0=E_R$，一般不等于当前导数。

**计算含义：** 完整当前 $L_D,R_D$ 必须由已付费状态取得；拿到缓存后进行 Gram、SVD 或坐标变换不新增 Maxwell RHS，取得当前场本身仍有真实费用。来源：[R12 双边推导第1节](../p3l_r12_bilateral_fixedpoint/BILATERAL_MODE_DERIVATION.md)、[物理适配器](../p3l_r12_bilateral_fixedpoint/shared/physics.py)；与 F-SS-002/003 的活动支撑导数一致。

## F-BF-002 · 对角材料的端口对算子、共轭 Hadamard Gram 与截断边界

按逐源、每源逐接收端排列复数据。材料 Jacobian $J_m\in\mathbb C^{(n_rn_s)\times q}$ 的第 $i$ 列为

$$j_i=R_{D,i,:}^{T}\otimes L_{D,:,i},\qquad J_m=R_D^T\odot L_D,$$

其中 $\odot$ 指逐列 Kronecker 积。逐项求内积得到

$$\boxed{G_m=J_m^*J_m=(L_D^*L_D)\circ\overline{R_DR_D^*}.}$$

证明是 $\langle j_i,j_j\rangle=(\sum_r\overline L_{ri}L_{rj})(\sum_s\overline R_{is}R_{js})$；第二因子必须为 $\overline{RR^*}=\overline R R^T$。两个因子均为 Hermitian 半正定，Schur 乘积也半正定。

进一步取薄 SVD，秩分别为 $r_L,r_R$：

$$L=U_L\Sigma_LV_L^*,\qquad R=U_R\Sigma_RV_R^*,$$
$$\boxed{K_{(a,b),i}=\sigma_{L,a}\sigma_{R,b}\overline{(V_L)_{ia}}(U_R)_{ib}.}$$

因为

$$L\operatorname{diag}(v)R=U_L[\Sigma_LV_L^*\operatorname{diag}(v)U_R\Sigma_R]V_R^*,$$

左右薄酉因子保持 Frobenius 范数，故 $\|Kv\|_2=\|J_mv\|_2$ 对任意复材料微分 $v$ 成立，完整 $K^*K=G_m$。这是对**独立接收与照明物理端口**的分解；最多 $r_Lr_R\le256$ 个复端口对，并非两个材料弱子空间的集合交集。

**假设与资格：** 使用同一当前 $L,R$ 和全部端口对时为精确表示恒等式。保留最大 $\sigma_{L,a}\sigma_{R,b}$ 的 top32 对是明确的截断启发式；这些乘积并未考虑对角材料限制及 $D_x$，不保证在原 66 维先验上实现最优秩32近似。单边 $L^*L$ 或 $\overline{RR^*}$ 均丢掉另一端信息。

**计算含义：** 完整 pair 与 Full-J 不是两种独立数据来源。top32 复行实化后秩至多64，可在66维材料像上制造截断零方向；source-only 至多16个实行，纯数据 weak16 的选择可能高度非唯一。来源：[双边推导第2—3节](../p3l_r12_bilateral_fixedpoint/BILATERAL_MODE_DERIVATION.md)、[modes.py 的 bilateral_gram / pair_modes](../p3l_r12_bilateral_fixedpoint/agent_a/modes.py)。新增的是本模型的显式对应，不是 Khatri–Rao、Schur 乘积或 SVD 原理。

## F-BF-003 · 原 Gaussianized66 的复材料导数与实白化拉回

原先验坐标为 $x=(z,\xi)\in\mathbb R^{64}\times\mathbb R^2$，$u_j=\Phi_N(\xi_j)$；$\Phi_N,\varphi_N$ 分别为标准正态 CDF、PDF。设冻结 quintic saturation 为 $S$，$s=(\mu+A_{ROI}z)/\tau$，则

$$a=S(s),\qquad \delta=0.8+\Phi_N(\xi_1)+i[0.03+0.15\Phi_N(\xi_2)],\qquad d=\delta a.$$

链式法则给出 $D_x\in\mathbb C^{256\times66}$：

$$D_{x,:,1:64}=\delta\operatorname{diag}(S'(s)/\tau)A_{ROI},$$
$$D_{x,:,65}=a\varphi_N(\xi_1),\qquad D_{x,:,66}=i\,0.15a\varphi_N(\xi_2).$$

令 $\operatorname{pack}$ 按每源32实部、32虚部排列，$W_{real}Z=(\sqrt2/\sigma)\operatorname{pack}Z$；对 Jacobian 逐列应用该实线性映射。对任意实方向 $v$，

$$\|J_wv\|^2={2\over\sigma^2}\|J_mD_xv\|^2,
\qquad J_w=W_{real}(J_mD_x)\in\mathbb R^{512\times66},$$
$$\boxed{H_{data}=J_w^TJ_w={2\over\sigma^2}\operatorname{Re}(D_x^*G_mD_x).}$$

证明用 $v$ 为实数及 $\|\operatorname{pack}Z\|^2=\|Z\|_F^2$，再由实二次型唯一确定对称矩阵。不能直接把复矩阵 $D_x^*G_mD_x$ 当实参数 Hessian，也不能漏掉两个 CDF 导数。$H_{data}$ 是 Gauss–Newton 数据曲率；非零残差下它不等于精确目标 Hessian。

**假设与资格：** 冻结原 decoder、固定噪声白化及实参数约定下的精确一阶拉回。平台内 $S'$ 可为零，因此局部材料秩不一定等于66。native 自由占比258模型与此先验支持不同。

**计算含义：** 完整当前 pair、完整当前 Gram 和真实 Full-J 拉回给出同一局部数据几何；不能把不同的归一化或重复 $\delta$ 因子写成方法收益。来源：[双边推导第4节](../p3l_r12_bilateral_fixedpoint/BILATERAL_MODE_DERIVATION.md)、[physics.py 的 Problem.decode / white_jac](../p3l_r12_bilateral_fixedpoint/shared/physics.py)、[modes.py 的 pullback_gram](../p3l_r12_bilateral_fixedpoint/agent_a/modes.py)。

## F-BF-004 · 材料有效局部商与原先验增广因子

ROI 外材料不变，但材料 RMS 仍按完整576格计算。定义实化材料微分

$$B={1\over24}\begin{bmatrix}\operatorname{Re}D_x\\\operatorname{Im}D_x\end{bmatrix}\in\mathbb R^{512\times66},\qquad G_x=B^TB={1\over576}\operatorname{Re}(D_x^*D_x).$$

取 $B=U\Sigma V^T$，保留明确阈值以上的奇异值，令 $r$ 为保留数：

$$T=V_+\Sigma_+^{-1}\in\mathbb R^{66\times r},\quad Q=BT=U_+,\quad T^TG_xT=I_r.$$

这是把单位坐标变化变为单位材料 RMS 微分的局部 lift。精确零奇异值对应局部 decoder 不改变材料的方向；数值阈值以下则只称不可分辨方向。R12 阈值为相对 $10^{-10}$ 与浮点尺寸阈值的较大者，不能由此证明全局 gauge 或全局商流形。

对于原 Gaussianized66 先验，局部目标的 GN 曲率为 $H_{data}+I$，故直接构造

$$\boxed{A_{aug}=\begin{bmatrix}J_wT\\T\end{bmatrix},\quad A_{aug}^TA_{aug}=T^T(H_{data}+I)T.}$$

在该显式实因子的右奇异向量中选最弱 $k$ 列 $V_{weak}$，令 $W=TV_{weak}$，则 $W^TG_xW=I_k$。纯数据版本仅分解 $J_wT$；不同消融以自己声明的数据因子替换上块，统一保留原先验下块 $T$。

**假设与资格：** 当地材料导数及固定原 Gaussian 参考测度下的标准白化/增广最小二乘。加入 $I$ 是保留原先验曲率，没有新增观测。该标准分解不能证明 weak16 的有限材料步可达，也不能证明 $0.02W$ 位于可信线性区；$\|W\|$ 可能很大。

**计算含义：** 直接 SVD 显式因子，避免先形成正规矩阵再经病态 $T$ 拉回时弱谱精度恶化。原欧氏阻尼 $\mu I$ 在材料坐标上成为 $\mu T^TT$，不能把 $J_wT$ 的奇异值直接代入带 $\mu I_r$ 的固定点公式。来源：[双边推导第5、7节](../p3l_r12_bilateral_fixedpoint/BILATERAL_MODE_DERIVATION.md)、[modes.py 的 material_canonical_basis / proposal_basis](../p3l_r12_bilateral_fixedpoint/agent_a/modes.py)。

## F-BF-005 · 同一个 Jacobian 的嵌套 SVD 不增加筛选信息

设实白化 Jacobian $J=U\Sigma V^T$，$V_I$ 是它自身一个奇异索引集 $I$ 的右奇异向量，则

$$JV_I=U_I\Sigma_I.$$

所以再对 $JV_I$ 做 SVD，其右奇异向量只是坐标轴；相同奇异值的块内可任选正交基，映回仍为原子空间。若另取它自身左奇异向量 $U_K$，

$$U_K^TJV_I=\Sigma_{K,I},\qquad P_LJP_R=U_K\Sigma_{K,I}V_I^T.$$

这是原对角奇异结构的子块。对严格不重叠的索引集，子空间正交；有交叠时交集也由原索引确定。重复分解、强弱空间求交或等值块重命名不会产生新的物理信息。

**假设与资格：** 两次使用同一状态、同一度量、同一 $J$ 的奇异向量。若第二次加入独立材料度量、原先验、不同状态的 Jacobian 或独立接收/照明算子，问题已经改变；这种变化须明确写出。本条不否定 F-BF-002 的双物理端口分解。

**计算含义：** 不为该纯冗余构造分配 Maxwell 预算，也不把它列为发现新材料方向的独立方法。来源：[R12 嵌套 SVD 解析说明](../p3l_r12_bilateral_fixedpoint/NESTED_SVD_REDUNDANCY_AUDIT.md)。这是既有 SVD 的直接推论，本轮新增的是对特定提案的排除条件。

## F-BF-006 · 阻尼数据映射在非根处的完整导数

令 $r:\mathbb R^p\to\mathbb R^{512}$ 为冻结白化残差，$r$ 二次连续可微，$J=Dr$，$\mu>0$ 为固定标量。用户给定的数据映射是

$$A(x)=J(x)^TJ(x)+\mu I,\quad g(x)=J(x)^Tr(x),\qquad \mathcal T_\mu(x)=x-A(x)^{-1}g(x).$$

对方向 $h$，记 $\dot J=DJ[h]$。由 $D(A^{-1})[h]=-A^{-1}(DA[h])A^{-1}$ 及

$$DA[h]=\dot J^TJ+J^T\dot J,\qquad Dg[h]=J^TJh+\dot J^Tr,$$

得到

$$\boxed{D\mathcal T_\mu[h]=h-A^{-1}(J^TJh+\dot J^Tr)+A^{-1}(\dot J^TJ+J^T\dot J)A^{-1}J^Tr.}$$

在响应根 $r(x_*)=0$，两项残差贡献同时消失，因而

$$D\mathcal T_\mu(x_*)=\mu(J^TJ+\mu I)^{-1},\qquad \lambda_i={\mu\over\sigma_i(J)^2+\mu}.$$

这只适用于与阻尼相同的欧氏参数度量。近单位特征值可以来自小的正奇异值，也可由增大 $\mu$ 人为造成；精确单位特征值只指一阶数据零方向。若 $J^Tr=0$ 但 $r\ne0$，最后一项消失，仍有

$$D\mathcal T_\mu=I-A^{-1}(J^TJ+S_r),\qquad S_r=\sum_a r_a\nabla^2r_a.$$

因此正残差驻点不能当作精确响应根。若阻尼本身随 $x$ 改变，还须在 $DA$ 中加入其导数；本式未作该假设。

**假设与资格：** 有限离散、固定阻尼、指定坐标下的精确微分。局部一阶盲方向不保证能积成有限精确根曲线；需另有正则秩/隐函数条件。满列秩根局部孤立，但不排除相隔较远的第二个根。

**计算含义：** 计算非根处固定点谱必须提供完整 $D\mathcal T$，不能仅由 $J$ 猜测。按 F-BF-004 选定的局部 lift 可报告 $Q^TB(D\mathcal T)T$ 的压缩谱；在没有零方向不变性等条件时，它只是所选材料表示的压缩，不能冒称全局商映射导数。来源：[固定点谱说明第1、3—4节](../p3l_r12_bilateral_fixedpoint/FIXED_POINT_SPECTRUM_AUDIT.md)、[engine.py 的 effective_spectrum](../p3l_r12_bilateral_fixedpoint/agent_b/engine.py)。已有两格控制分别保留近单位但孤立根与非零残差曲率证据，本次不重复运行。

## F-BF-007 · 数据阻尼、Gaussianized MAP 与 RML 的三个不同目标

F-BF-006 的 $\mu I$ 只在步矩阵中出现，梯度仍是 $J^Tr$。即使 $\mu=1$，缺少 $+x$ 梯度便不能把它称为标准 Gaussian MAP。

原 Gaussianized66 MAP 目标为

$$\mathcal L_{MAP}(x)=\tfrac12\|r(x)\|^2+\tfrac12\|x\|^2,\qquad \nabla\mathcal L_{MAP}=J^Tr+x.$$

固定随机目标 $\eta\in\mathbb R^{512},\xi\in\mathbb R^p$ 的 RML 目标及其 GN 映射为

$$\mathcal L_{\eta,\xi}(x)=\tfrac12\|r(x)-\eta\|^2+\tfrac12\|x-\xi\|^2,$$
$$\mathcal T_{RML}(x)=x-(J^TJ+I)^{-1}\{J^T(r-\eta)+(x-\xi)\}.$$

令 $A=J^TJ+I$、$g_{RML}=\nabla\mathcal L_{\eta,\xi}$，直接微分得

$$D\mathcal T_{RML}[h]=h-A^{-1}\nabla^2\mathcal L_{\eta,\xi}\,h+A^{-1}DA[h]A^{-1}g_{RML}.$$

在该目标驻点，$\nabla^2\mathcal L=A+S_{r-\eta}$，所以

$$\boxed{D\mathcal T_{RML}=-A^{-1}S_{r-\eta}.}$$

这与数据根处的 $\mu A^{-1}$ 不是同一谱公式。原 native $(z,u)$ 参考测度中的 Uniform2 经 $u=\Phi_N(\xi)$ 后变成 Gaussian2；概率先验可由此正确推送，但相对于不同坐标 Lebesgue 测度的 MAP 点一般不保持不变。不得把两个点目标静默混用。

**假设与资格：** 原 Gaussianized 坐标、固定随机目标、未加主动约束的光滑映射。RML 是随机优化候选机制；对非线性正演，有限 RML 端点不会自动成为精确后验样本。主动约束和回溯引擎还会改变实际迭代映射。

**计算含义：** 主臂应分别声明提案几何、校正目标与候选资格门。先验增广可用于选择方向，但不能据此把数据阻尼映射称为 MAP，或把合格有限云称为后验体积。来源：[固定点谱说明第2节](../p3l_r12_bilateral_fixedpoint/FIXED_POINT_SPECTRUM_AUDIT.md)、[engine.py 的 _objective / transverse_correct](../p3l_r12_bilateral_fixedpoint/agent_b/engine.py)、[开发运行脚本](../p3l_r12_bilateral_fixedpoint/shared/run_development.py)。

## F-BF-008 · 固定弱坐标的横向正确器与人为单位模态

锚点 $x_0$ 处固定 $G_0$ 和材料正交弱基 $W_0$，定义

$$C_w=W_0^TG_0\in\mathbb R^{k\times p},\quad C_wW_0=I_k,\qquad c(x)=C_w(x-x_0).$$

预测点 $x_{pred}=x_0+W_0t$ 满足 $c(x_{pred})=t$。若以 $Z$ 张成 $\ker C_w$，则任意有限校正 $x=x_{pred}+Z\alpha$ 都严格保持相同仿射坐标：

$$\boxed{C_w(x-x_0)=t.}$$

这是有限坐标恒等式，不只是每一步的一阶近似。若正确器还限制在当前有效材料切空间 $\operatorname{range}T(x)$，使用

$$Z(x)=T(x)\operatorname{null}(C_wT(x)),$$

同样有 $C_wZ(x)=0$；即使切空间变化，原 $C_w$ 必须在该次正确器中保持不变。仅使用 $\operatorname{null}C_w$ 则还包括局部 decoder 零方向。原生占比边界须把向外活动坐标作为额外约束或合法线段处理，事后裁剪一般会破坏固定坐标。

若把此规则延拓为全空间映射并满足 $C_w\mathcal T(x)=C_wx$，求导即得

$$C_wD\mathcal T=C_w.$$

这些是算法主动保持的单位**左**模态。判断横向收缩应在其实际自由域 $\ker C_w$ 及合法切空间中进行，不能把人为冻结的 $k$ 个坐标算作 Maxwell 连续根流形维数。

**假设与资格：** 固定锚点仿射泛函、校正步精确落在其核内。$c(x)$ 只是一组声明的有限 latent 坐标，不保证等于非线性材料投影，不保证目标半径可达，也不保证正确器能找到噪声合格点。

**计算含义：** H1 固定 frame、H2 分层刷新和 HA 自适应刷新改变的是锚点/提案策略；每次校正仍须保持自己声明的固定 $C_w$。边界处没有可行预测或横向步时应报告未到达机制试验。来源：[双边推导第5节](../p3l_r12_bilateral_fixedpoint/BILATERAL_MODE_DERIVATION.md)、[modes.py 的 coordinate_preserving_complement](../p3l_r12_bilateral_fixedpoint/agent_a/modes.py)、[engine.py 的 weak_frame / transverse_correct](../p3l_r12_bilateral_fixedpoint/agent_b/engine.py)。

## F-BF-009 · 有限 secant、投影变化和覆盖只提供经验诊断

对两个已实际求值的点 $x_a,x_b$，令 $\Delta x=x_b-x_a$、$\Delta r=r(x_b)-r(x_a)$。冻结起点 Jacobian 的 secant 缺陷是

$$e_{a\to b}=\Delta r-J_a\Delta x,\qquad \delta_{a\to b}={\|e_{a\to b}\|\over\max(\|\Delta r\|,1)}.$$

若整条线段位于光滑合法域，微积分基本定理给

$$e_{a\to b}=\int_0^1[J(x_a+t\Delta x)-J_a]\Delta x\,dt.$$

只有额外获得全线段 Jacobian 的 Lipschitz 常数 $L_J$，才能推出 $\|e\|\le\tfrac12L_J\|\Delta x\|^2$。有限端点比较本身不提供这样的常数。终点刷新 $J_b$ 的缺陷 $\Delta r-J_b\Delta x$ 可衡量这个 secant 的经验改善，不能证明全邻域的二阶控制或刷新必有材料收益。

把不同点的材料弱空间放入同一完整材料实空间，令 $Q_a,Q_b$ 的列正交；则

$$\rho_{a,b}={\|(I-Q_aQ_a^T)Q_b\|_F\over\sqrt{k}}$$

衡量新 $k$ 维空间在旧空间之外的比例。它是子空间变化，不是新增观测量；data-only weak8 的诊断也不是 prior-aware weak16 的提案轨迹。

有限合格云 $\mathcal C_{p,h}$ 的冻结经验覆盖为

$$\widehat{coverage}_h(\epsilon)={1\over n_{DEV}}\sum_p\mathbf1\!\left\{\min_{m\in\mathcal C_{p,h}}{\|m-m_p^{GT}\|\over24}\le\epsilon\right\}.$$

空云按失败计；最近 GT 成员是离线容量诊断，不能用来选择在线输出。此处的 $\epsilon$、噪声门、去重尺度和先封存预测后读 GT 的顺序须由既有协议给定，不能看到结果后移动。

**假设与资格：** secant 和投影差来自有限已存节点；只有在连续线段已满足光滑假设时，积分恒等式才适用。经验覆盖只针对声明的暴露 DEV 对象，不能升级为总体高概率保证、连续根族或后验覆盖。均匀或 $\exp(-\Phi)$ 有限云行动仍须注明有限云语义。

**计算含义：** 现有 R12 诊断保存128条 secant、96个材料投影对，且采样包含被拒 trial 与分支切换，不能称为连续接受路径。H2/HA 刷新次数更多但是否改善行动，必须看已有真实有限材料结果与总成本；局部方向转动本身不能代替收益。来源：[在线诊断已有收据](../p3l_r12_bilateral_fixedpoint/agent_c/p1_diagnostics/ONLINE_DIAGNOSTIC_RECEIPT.json)、[FINITE_SECANTS.json](../p3l_r12_bilateral_fixedpoint/agent_c/p1_diagnostics/FINITE_SECANTS.json)、[冻结评分汇总](../p3l_r12_bilateral_fixedpoint/agent_c/p1_score/SUMMARY.json)。本次没有重算这些诊断。

### 既有 R12 结果对本章的约束

截至本章维护时，P1 的4个暴露 DEV 对象、8个臂在冻结 $\epsilon=0.075$ 下均为覆盖0/4。普通 H0 RML 的平均离线最近材料 MSE 为0.0119253，H1/H2/HA 完整16维分别为0.0130177、0.0133596、0.0131985；该基线的云容量指标更好。所有32组的在线选择仍为原锚点，平均所选 MSE 同为0.0148154。这里支持“现有开发结果没有显示分层/自适应刷新带来在线行动改善”，不支持“所有双边方法都不可能有效”。数字来自已有[评分汇总](../p3l_r12_bilateral_fixedpoint/agent_c/p1_score/SUMMARY.json)与[候选门收据](../p3l_r12_bilateral_fixedpoint/agent_c/P1_ADVERSARIAL_QA.json)，本次没有重新评分。

原生 N1 的既有45组运行停留在初始化边界，未到达非锚点非线性机制比较；刷新是同一锚点的代数 frame 重算。其状态保留为 `INITIALIZATION_BOUNDARY_FAILURE_NOT_PHYSICAL_TEST`，不能作为新原生材料方法的正面或负面性能结论。后续 N2 的实现、封存和真值顺序由主线程及对应工作目录独立管理，本章不提前赋予结果。来源：[N1 已有边界收据](../p3l_r12_bilateral_fixedpoint/agent_c/N1_BOUNDARY_AUDIT.json)。

### 已知结构、待比较先例与仍未建立的贡献

F-BF-001—009 使用 resolvent 微分、Khatri–Rao/Hadamard 结构、SVD、最小二乘增广、链式法则、隐式约束和有限经验评分等既有工具。本章贡献是把这些工具在同一个原先验/白化/材料 RMS 合同中写清，并排除不成立的谱与新颖性推论；它不声称创造了这些数学工具。

已有 ScholarQA 记录核实的最近主题先例为 Zhong、Chen 的 [Twofold subspace-based optimization method for solving inverse scattering problems (2009)](https://doi.org/10.1088/0266-5611/25/8/085003)、[FFT Twofold Subspace-Based Optimization (2011)](https://doi.org/10.1109/TAP.2010.2103027)，以及 Xu 等的 [Multiplicative-Regularized FFT Twofold Subspace-Based Optimization (2015)](https://doi.org/10.1109/TGRS.2014.2329032)。这里只复用已完成的 **Tier C 出版身份** 核实；没有可用逐段全文证据证明它们与本轮具体算法相同或不同。429 和未覆盖的 distorted-Born、张量分解、广义奇异值与先验感知流形优化文献缺口均保留，不能宣称“首次”。来源：[既有最近工作身份说明](../p3l_r12_bilateral_fixedpoint/agent_a/CLOSEST_PRIOR_WORK_METADATA.md)、[已有五条身份核实记录](../p3l_r12_bilateral_fixedpoint/agent_a/SCHOLARQA_VERIFIED.json)。

状态：上述公式为 `DERIVED_FINITE_MODEL`；实现仅在所引源码的声明范围内；P1 改善结论为 `NOT_DEMONSTRATED_ON_EXPOSED_DEV`；原 P1 连续根流形、完整材料后验与总体覆盖为 `NOT_ESTABLISHED`；方法新颖性为 `UNRESOLVED_WITH_METADATA_ONLY_PRECEDENTS`。最终科学判断由主线程作出。
