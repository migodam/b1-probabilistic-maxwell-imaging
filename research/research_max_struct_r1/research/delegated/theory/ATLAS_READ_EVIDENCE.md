# Scoped source reading evidence

All paths relative to repository root. Excerpts are local audit material; no historical edits or reruns.


## research/maxwell_formula_atlas/01_CONVENTIONS_AND_DEPENDENCIES.md


### Lines 5-57

5: ## F-BASE-001 · 物理离散空间与矩阵尺寸
6: 
7: 采用时间约定 $e^{-i\omega t}$，被动材料的介电虚部非负。原冻结二维 scalar-TM 离散合同为：
8: 
9: $$N=576,\quad q=256,\quad n_s=8,\quad n_r=32,\quad d_x=66.$$
10: 
11: 设 $K\in\mathbb C^{N\times N}$ 是包含 $k_b^2$ 与单元积分权重的内部 Green 算子，$S\in\mathbb C^{n_r\times N}$ 为接收算子，$E^{inc}\in\mathbb C^{N\times n_s}$ 为入射场；$m\in\mathbb C^N$ 为实际对比材料，$D_m=\operatorname{diag}(m)$。$J_c\in\mathbb C^{N\times n_s}$ 表示诱导电流，$\mathcal J$ 表示实化 latent Jacobian，两者不混用。
12: 
13: $$E=E^{inc}+KJ_c,\qquad J_c=D_mE,\qquad F(m)=SJ_c.$$
14: 
15: 由代入得到两种完全等价的求解形式：
16: 
17: $$A_mE=E^{inc},\quad A_m=I-KD_m;\qquad L_mJ_c=D_mE^{inc},\quad L_m=I-D_mK.$$
18: 
19: 这里不可交换 $KD_m$ 与 $D_mK$。二者可逆性等价，因为 $\det(I-AB)=\det(I-BA)$；“field form”与“current form”的逆矩阵本身不相同。原代码中的 `G` 对应本章 $K$，避免与材料 decoder $\mathcal G$ 同名。
20: 
21: **资格：** 冻结离散模型恒等式，非连续 PDE 误差证书。来源：R4/P1 `core.py`、A30 frozen VIE、R9 `src/stable_roi_operator.py`。
22: 
23: ## F-BASE-002 · 原始观测、实化和白化
24: 
25: 每个 $F_{rs}$ 是一个复观测；$\eta_{rs}\sim\mathcal{CN}(0,\sigma^2)$，所以每个实部和虚部方差均为 $\sigma^2/2$。固定
26: 
27: $$\sigma=0.00119211798514519.$$
28: 
29: 令 $\operatorname{pack}$ 按 source-major 顺序逐源排列实部32项、虚部32项。它是实线性等距映射：
30: 
31: $$\|\operatorname{pack}Z\|_2^2=\|Z\|_F^2.$$
32: 
33: 物理实化 $f=\sqrt2\operatorname{pack}F$，数据 $d=\sqrt2\operatorname{pack}Y$，白化 residual
34: 
35: $$r(x)=\frac{f(x)-d}{\sigma},\quad \Phi(x)=\tfrac12\|r(x)\|^2=\frac{\|F(\mathcal G(x))-Y\|_F^2}{\sigma^2}.$$
36: 
37: 含噪观测残差、两个无噪预测之间的距离与白化导数必须使用同一 packing。$\Delta\Phi\le2$ 是当前似然竞争判据；它不是 posterior credible mass。
38: 
39: ## F-BASE-003 · 材料损失与空间度量
40: 
41: 材料作为实 Hilbert 空间，内积 $\langle a,b\rangle_M=\operatorname{Re}(a^*b)/576$。于是
42: 
43: $$\|m_1-m_2\|_M=\frac{\|m_1-m_2\|_2}{24},\qquad R(\hat m)=\mathbb E\|\hat m(Y)-M\|_M^2.$$
44: 
45: ROI-only RMS 的分母为16；除非明确改写，R7/R9半径0.01/0.02使用全网格分母24。固定宿主之外没有材料差异，并不允许静默改分母。复数 MSE 同时包含实部与虚部平方，不能只比较实部图。
46: 
47: ## F-BASE-004 · 哪些算子真正固定
48: 
49: 在背景、频率、网格、源接收位置及归一化固定时，$K,S,E^{inc}$ 固定；ROI嵌入矩阵 $P\in\mathbb R^{576\times256}$ 固定，$P^TP=I$。材料 $m$、$D_m$、总场 $E$、电流 $J_c$、resolvent $L_m^{-1}$、decoder导数以及 $\mathcal J(x)$ 一般变化。
50: 
51: 由此有三个不同的复用层：
52: 
53: 1. **固定几何层：** $K,S,E^{inc}$、receiver SVD、Fourier/Wavelet基。
54: 2. **固定背景消元层：** 给定宿主后预计算有效ROI传播与观测端口；不同材料仍保留完整非线性。
55: 3. **变化状态层：** 当前材料的小系统、场与 adjoint；可用近似预条件，但不能无误差地冻结全部内容。
56: 
57: SOM current 子空间、材料表示基与当前 $\mathcal J$ 的弱方向分别属于不同空间。它们之间需要显式映射，维数相同也不构成等价。


### Lines 76-80

76: ## 阅读规则与维护方式
77: 
78: 每个条目包含等式、推导、成立条件及其计算含义。纯代数正确性由证明审核；不为“再看到一次机器精度相等”追加物理调用。实现纠错应针对实际错误。只有成本排序、真实材料覆盖、统计风险和误差界是否足够紧等无法由代数决定的问题，才进入路径裁决或论文实验。
79: 
80: 公式资格分为 `EXACT_DISCRETE_IDENTITY`、`CONDITIONAL_THEOREM`、`CONTROLLED_APPROXIMATION`、`EMPIRICAL_ONLY`。没有统一的“全部PASS”标签。标准恒等式可用于重要优化，但不因此成为本研究的新定理。


## research/maxwell_formula_atlas/02_FIXED_OPERATORS_AND_PROJECTIONS.md


### Lines 5-61

5: ## F-OP-001 · 连续 Helmholtz 与体积分方程
6: 
7: 采用时间因子 $e^{-i\omega t}$。在均匀背景、标量介电对比、固定磁导率的模型中，令 $k_b>0$、$\chi(\mathbf x)=k(\mathbf x)^2/k_b^2-1$，有
8: 
9: $$[\Delta+k_b^2(1+\chi)]u=-s,\qquad (\Delta+k_b^2)u^{inc}=-s.$$
10: 
11: 出射 Green 函数满足 $(\Delta+k_b^2)g_b=-\delta$。因此
12: 
13: $$u(\mathbf x)=u^{inc}(\mathbf x)+k_b^2\int_\Omega g_b(\mathbf x,\mathbf x')\chi(\mathbf x')u(\mathbf x')\,d\mathbf x'.$$
14: 
15: 证明只需对右侧作用 $\Delta+k_b^2$：积分项变成 $-k_b^2\chi u$，移项即得原方程；出射条件选择相应唯一分支。定义对比电流密度 $j_c=\chi u$ 后，非线性仅通过构成关系进入；传播 $j_c\mapsto u-u^{inc}$ 是固定线性算子。
16: 
17: 二维 scalar-TM 与三维 scalar Helmholtz 的核分别为
18: 
19: $$g_2(R)=\frac{i}{4}H_0^{(1)}(k_bR),\qquad g_3(R)=\frac{e^{ik_bR}}{4\pi R}.$$
20: 
21: 三维**矢量** Maxwell 一般需要 dyadic Green 张量及极化分量，不能把这里的标量 $1/R$ 核直接当作完整矢量模型。吸收背景可以使用相应复波数，但下面的远场幅相估计写在实 $k_b>0$ 条件下；复波数还包含距离相关衰减。
22: 
23: **资格：** `CONDITIONAL_THEOREM`，标准标量散射表示；它没有替离散积分误差、边界模型或连续 Maxwell 精度提供证书。
24: 
25: ## F-OP-002 · 二维与三维远场：振幅近似和相位近似分开
26: 
27: 令接收点为 $R\widehat{\mathbf r}$、目标坐标满足 $|\mathbf x'|\le a$。几何展开为
28: 
29: $$|R\widehat{\mathbf r}-\mathbf x'|=R-\widehat{\mathbf r}\cdot\mathbf x'+\frac{|\mathbf x'|^2-(\widehat{\mathbf r}\cdot\mathbf x')^2}{2R}+O(a^3/R^2).$$
30: 
31: 因此，在 $R\gg a$ 且 $k_ba^2/R\ll1$ 时，三维有
32: 
33: $$g_3(R\widehat{\mathbf r},\mathbf x')\simeq\underbrace{\frac{1}{4\pi R}}_{\text{振幅}}\underbrace{e^{ik_bR}e^{-ik_b\widehat{\mathbf r}\cdot\mathbf x'}}_{\text{相位}}.$$
34: 
35: 二维还需 $k_bR\gg1$ 才能使用 Hankel 渐近式，得到
36: 
37: $$g_2(R\widehat{\mathbf r},\mathbf x')\simeq\underbrace{\frac{1}{\sqrt{8\pi k_bR}}}_{\text{振幅}}\underbrace{e^{i\pi/4}e^{ik_bR}e^{-ik_b\widehat{\mathbf r}\cdot\mathbf x'}}_{\text{相位}}.$$
38: 
39: 振幅在目标内近似常数的相对误差规模是 $a/R$；丢弃二次路径差造成的相位误差规模是 $k_ba^2/R$。二维另有 Hankel 渐近误差规模 $(k_bR)^{-1}$。这些条件不同：距离足以让振幅变化小，仍可能不足以让相位线性化。近场二维核也不能简化为精确的 $R^{-1/2}$ 振幅。
40: 
41: 若入射为平面波 $u^{inc}=a_s e^{ik_b\widehat{\mathbf s}\cdot\mathbf x'}$，Born 远场的空间积分为
42: 
43: $$\int\chi(\mathbf x')e^{-ik_b(\widehat{\mathbf r}-\widehat{\mathbf s})\cdot\mathbf x'}d\mathbf x'.$$
44: 
45: 这解释了固定收发几何对应的 Fourier 采样；它同时使用弱散射、远场接收和平面波照明。有限距离点源需保留自身 Green 振幅与相位，不能自动套用单一 Fourier 频率。当前离散核保留原合同，不因上述解释而替换成远场核。
46: 
47: ## F-OP-003 · 离散 Born 双边结构及其非线性边界
48: 
49: 对体积分作固定单元离散，令 $m_i$ 为单元对比材料。把 $k_b^2$、体积权重和自单元积分处理纳入内部 $K$ 与接收 $S$，而不再重复乘入 $m$。以所有源为列：
50: 
51: $$E^{inc}\in\mathbb C^{N\times n_s},\quad K\in\mathbb C^{N\times N},\quad S\in\mathbb C^{n_r\times N},\quad D_m=\operatorname{diag}(m).$$
52: 
53: 零对比背景下，第一次 Born 替换 $E$ 为 $E^{inc}$，所以
54: 
55: $$F_B(m)=S D_m E^{inc}.$$
56: 
57: 若显式写源传播 $G_{source}$ 与源振幅矩阵 $A_s$，则 $E^{inc}=G_{source}A_s$，故 $F_B=G_{receiver}D_mG_{source}A_s$，其中 $G_{receiver}=S$。两侧固定、中央对角材料变化，正是固定双边结构。Born 的来源是
58: 
59: $$J_c=D_m(E^{inc}+KJ_c)=D_mE^{inc}+D_mKD_mE^{inc}+\cdots.$$
60: 
61: Neumann 展开需要相应收敛条件，例如 $\|D_mK\|<1$ 是充分条件；方程本身可在该条件不成立而 $I-D_mK$ 仍可逆时求解。截掉高阶项是物理近似，不是单纯换基。非零宿主不能直接套用零背景 Born；F-OP-008 给出宿主上的精确响应，其一阶项才是 distorted Born。


### Lines 63-169

63: ## F-OP-004 · 固定接收 SVD 到底能确定什么
64: 
65: 先考虑同一固定、可逐源使用的白化 $W\in\mathbb C^{n_r\times n_r}$。当前合同是相同方差的独立复噪声，可取 $W=I/\sigma$；若噪声在源间非可分离相关，应对白化后的整体向量算子做 SVD，不能未经说明沿用一个共享接收矩阵。
66: 
67: 令 $A=WS=U_r\Sigma_rV_r^*$ 为保留全部非零奇异值的薄 SVD，$r=\operatorname{rank}A$，$V_0$ 是其正交补。任意电流有唯一分解
68: 
69: $$J_c=V_r a+V_0b,\qquad a\in\mathbb C^{r\times n_s},\quad b\in\mathbb C^{(N-r)\times n_s}.$$
70: 
71: 正交分解观测残差得到
72: 
73: $$\|AJ_c-WY\|_F^2=\|\Sigma_ra-U_r^*WY\|_F^2+\|(I-U_rU_r^*)WY\|_F^2.$$
74: 
75: 第二项与材料无关，所以在**固定模型的归一化 posterior**中可以约掉相应常数；若比较模型 evidence，必须保留其值。第一项仍依赖 $a$，含噪推断不能把 $a$ 无条件固定为 $\Sigma_r^{-1}U_r^*WY$，否则把有限噪声似然改成了精确数据约束。
76: 
77: 固定投影为 $P_o=V_rV_r^*$、$P_n=I-P_o$，且 $SP_n=0$。若只取前 $\ell<r$ 个奇异值，则所谓“模糊”补空间还含弱辐射方向；它不是严格接收零空间。观测截断误差见 F-ERR-007。
78: 
79: **固定项：** $W,S,U,\Sigma,V$ 在采集与噪声固定时固定。**变化项：** 合法材料实际激发的 $a,b$、总场、resolvent 与 latent Jacobian。即使 $S$ 固定，也不能据此冻结材料 Jacobian。
80: 
81: ## F-OP-005 · SOM/TSOM 的投影不是两个空间的交集
82: 
83: 作者综述 [Chen、Zhong、Agarwal 2010](https://doi.org/10.4310/MAA.2010.v17.n4.a6) §3.1 式28—34描述接收 SVD 电流分解，§3.3 式46—47及其后续讨论描述内部作用筛选与 Fourier 表示。这些是已有方法构件。以下用统一符号解释其几何含义，不把不同篇章的算法目标函数混为一式。
84: 
85: 设 $P_w$ 为去除已选接收强方向的正交投影，$V_D$ 是内部传播所选的方向矩阵。一个典型二次筛选构造是
86: 
87: $$Z=P_wV_D.$$
88: 
89: 其列空间是内部方向**投影到**接收模糊空间后的像 $\operatorname{range}Z$，一般不是 $\operatorname{range}(V_D)\cap\operatorname{range}(P_w)$。最小反例：在 $\mathbb R^2$ 中，$P_w$ 投到 $e_1$，而 $V_D=(e_1+e_2)/\sqrt2$。则 $\operatorname{range}(P_wV_D)=\operatorname{span}(e_1)$，两原空间交集却只有零。
90: 
91: 投影后的列可以相关或消失，应对 $Z$ 作一次 rank-revealing QR/SVD。真正的正交投影是
92: 
93: $$P_Z=ZZ^\dagger=Z(Z^*Z)^\dagger Z^*,$$
94: 
95: 而一般不是 $P_wV_DV_D^*P_w$。后一个矩阵虽为半正定，但通常不幂等；把它当投影会重复缩放方向。
96: 
97: 若两个正交投影 $P,Q$ 对易，则 $PQ=QP$ 为交集的正交投影；反之，$PQ$ 若为正交投影，由 $(PQ)^*=QP$ 可知必须对易。一般 $PQ$ 是串联线性作用，不是交集投影。需要**真正交集**时可设
98: 
99: $$C_\cap=\begin{bmatrix}I-P\\I-Q\end{bmatrix},\qquad P_\cap=I-C_\cap^\dagger C_\cap.$$
100: 
101: 因为 $\ker C_\cap=\operatorname{range}P\cap\operatorname{range}Q$，该式才给交集。上述都是标准正交几何，不应冠以新的子空间定理。
102: 
103: ## F-OP-006 · 固定投影为何通常不与非均匀材料交换
104: 
105: 对任意固定矩阵 $P_o$ 和对角材料 $D_m$，逐元素直接得到
106: 
107: $$[D_m,P_o]_{ij}=(m_i-m_j)(P_o)_{ij}.$$
108: 
109: 标量材料 $D_m=vI$ 与任何投影对易；一般非均匀材料只有在投影不混合不同材料值的格点时才对易。接收 SVD 通常混合许多空间格点，因此不能把 $P_oD_m$ 静默换成 $D_mP_o$。
110: 
111: 进一步，将精确 ROI 方程的电流按某个固定 $P_o,P_n$ 分块，观测子块独立于暗子块需要
112: 
113: $$P_oD C P_n=0,\qquad D=\operatorname{diag}(m_{ROI}-m_h).$$
114: 
115: 仅有 $OP_n=0$ 不足以满足该条件。暗电流经内部传播和材料相乘后可能回到可观测方向，因而一般有 $OL_D^{-1}P_n\ne0$。F-ERR-001讨论 $[D,C]$ 对固定对角化的限制；这里强调的是接收不可见性并不在非线性构成方程下自动封闭。
116: 
117: ## F-OP-007 · 多源共享材料是电流合法性的硬条件
118: 
119: 给定任意候选 $J_c$，先由固定传播计算 $E=E^{inc}+KJ_c$。它能来自一幅共同标量材料，必须在每个格点满足
120: 
121: $$J_{c,is}=m_iE_{is}\quad\text{对所有源 }s.$$
122: 
123: 若该行 $\|E_{i,:}\|^2>0$，无约束复数最小二乘的唯一最佳材料为
124: 
125: $$\widehat m_i=\frac{\sum_s\overline{E_{is}}J_{c,is}}{\sum_s|E_{is}|^2},\qquad R_{i,:}=J_{c,i,:}-\widehat m_iE_{i,:}.$$
126: 
127: 共同复材料存在当且仅当所有这些行的 $R_{i,:}=0$，且在 $E_{i,:}=0$ 的行有 $J_{c,i,:}=0$。后一情形材料未被该状态识别；交叉小行列式 $E_{is}J_{it}-E_{it}J_{is}=0$ 在 $E=0,J\ne0$ 时也成立，因此只检验小行列式会漏掉非法电流。
128: 
129: 无约束逐行合法还不够：实际材料必须处于原 decoder 像、满足宿主/ROI 支持及允许的实虚部范围。若要求实材料，最佳无约束系数改为分子实部除分母；若要求被动材料，在本时间约定下介电对比虚部非负。但逐格截断或投影不能替代耦合 decoder 约束，更不能自动产生原 latent prior 的样本。
130: 
131: 源共享的 $m_i$ 是关键。为每个源各算 $m_{is}$，即使都能复现数据，也通常对应多幅不同材料。SOM 的电流分解须与状态/构成关系一起使用；其已有先例不能被简化成“从接收零空间自由采样就是合法材料”。
132: 
133: ## F-OP-008 · 固定宿主与完整 ROI 的精确消元
134: 
135: 令固定背景材料为 $m_h\in\mathbb C^N$、$D_h=\operatorname{diag}(m_h)$，ROI 选择嵌入为 $P\in\mathbb R^{N\times q}$，$P^TP=I_q$。允许原 ROI 内任意非均匀增量 $d\in\mathbb C^q$：
136: 
137: $$D_m=D_h+PDP^T,\qquad D=\operatorname{diag}(d).$$
138: 
139: 定义固定背景 field/current 系统与背景解
140: 
141: $$A_h=I-KD_h,\quad L_h=I-D_hK,\quad E_h=A_h^{-1}E^{inc},\quad J_h=D_hE_h,\quad F_h=SJ_h.$$
142: 
143: 假设背景系统可逆。预计算
144: 
145: $$B=A_h^{-1}KP\in\mathbb C^{N\times q},\quad V_h=L_h^{-1}P\in\mathbb C^{N\times q},$$
146: $$C=P^TB\in\mathbb C^{q\times q},\quad E_R=P^TE_h\in\mathbb C^{q\times n_s},\quad O=SV_h\in\mathbb C^{n_r\times q}.$$
147: 
148: 两种记号完全相同，因为
149: 
150: $$L_h^{-1}=I+D_hA_h^{-1}K,\qquad A_h^{-1}K=KL_h^{-1},$$
151: $$V_h=P+D_hB,\qquad C=P^TKV_h,\qquad O=S(P+D_hB).$$
152: 
153: 第一恒等式可由左乘 $L_h$ 并用 $(I-KD_h)A_h^{-1}=I$ 验证；第二来自 $A_hK=KL_h$。这同时说明不能丢掉观测端口中的宿主再散射 $SD_hB$。当前代码固定标量宿主 `host`，因此其 `Ah=I-host*G`、`B=Ah^-1 G[:,ROI]`、`C=B[ROI]`、`O=S[:,ROI]+host*S@B` 与上述一般式严格一致。
154: 
155: 令 $h=DP^TE$ 为 ROI 新增对比产生的电流。原方程给出
156: 
157: $$E=E_h+Bh,\quad J_c=J_h+V_hh,\quad P^TE=E_R+Ch,$$
158: 
159: 所以所有材料状态只需求解
160: 
161: $$\boxed{L_Dh=DE_R,\quad L_D=I_q-DC,\quad F=F_h+Oh.}$$
162: 
163: 等价 field 形式是 $(I_q-CD)E_{ROI}=E_R$、$h=DE_{ROI}$。这两个矩阵顺序不同，不能混写。它们可逆性等价，且在背景可逆时与完整材料系统可逆性等价：由矩阵行列式引理，$\det L_m=\det L_h\det(I_q-DC)$。
164: 
165: 推导没有除以 $d_i$，因此 $d_i=0$ 或部分材料恰等于宿主都合法。保留全 $q=256$ ROI 的消元适用于原 Full66 decoder 的非均匀材料；它不是均匀幅度假设，也没有削减材料支持、改变 prior 或忽略多重散射。
166: 
167: **成本：** 稠密实现一次背景 $N$ 阶 LU，$q+n_s=264$ 个完整 RHS，以及构造 $C,O$；之后每状态一个 $q$ 阶 LU、$n_s=8$ 个 reduced RHS 与观测乘法。固定存储可包括 $B$ 或 $V_h$，但只求接收预测时不必重建完整 $E,J_c$。这些固定项只在背景、频率、采集与 ROI 不变时复用。改变这些条件需要重新 setup。
168: 
169: **资格：** `EXACT_DISCRETE_IDENTITY`，标准背景消元。低秩截断 $C$、截断 $O$ 或减小 ROI 才引入新的近似；其校正/误差控制见第05章。矩阵浮点残差合格不等于连续 PDE 认证。


## research/maxwell_formula_atlas/03_MATERIAL_GEOMETRY_AND_CORRECTION.md


### Lines 9-44

9: ### F-MAT-001 原 P1 维度与参考测度
10: 
11: \[
12: q=(z,u_R,u_I)\in\mathbb R^{64}\times[0,1]^2,\quad
13: z\sim N(0,I_{64}),\quad u_R,u_I\overset{\rm iid}{\sim}U[0,1].
14: \]
15: 
16: 材料 \(\chi\in\mathbb C^{576}\)，24×24 网格，8个照明、每照明32个复接收通道；\(F(\chi)\in\mathbb C^{32\times8}\)。先验密度相对于 \(dz\,du_Rdu_I\)，不是材料 Lebesgue 体积。原 P1 不是1152个独立实材料像素。
17: 
18: ### F-MAT-002 probit Gaussian66
19: 
20: 用 \(\Phi_N,\varphi_N\) 表示标准正态 CDF/PDF，与似然代价 \(\Phi_y\) 区分：
21: 
22: \[
23: x=(z,\xi_R,\xi_I),\quad u_a=\Phi_N(\xi_a),\quad
24: \varphi_N(t)=(2\pi)^{-1/2}e^{-t^2/2},\quad
25: du_a=\varphi_N(\xi_a)d\xi_a.
26: \]
27: 
28: 内点精确变换使 \(x\sim N(0,I_{66})\)，没有额外添加物性先验。原生 \(u_a=0,1\) 对应无限 probit 坐标。R4 仅初始化时将原生边界移入 \([10^{-4},1-10^{-4}]\)；这不是普遍的边界等价证明。
29: 
30: ### F-MAT-003 坐标 MAP 与 GN Hessian
31: 
32: \[
33: U_q=\Phi_y+\tfrac12\|z\|^2,\qquad
34: U_x=\Phi_y+\tfrac12\|x\|^2.
35: \]
36: 
37: 正确 Jacobian 描述同一后验测度，但两个参考测度下的密度众数通常不同。R6 native 阶段使用 \(U_q\)，probit 阶段使用 \(U_x\)；R7/R9 的实际校正选择用 \(U_x\)。probit-coordinate mode 不能称为 native MAP。
38: 
39: \[
40: \nabla^2U_x=J_x^TJ_x+I+\sum_{i=1}^{512}r_i\nabla^2r_i,\qquad
41: H_{\rm GN}=J_x^TJ_x+I.
42: \]
43: 
44: GN 矩阵不是自动等于精确后验 Hessian，数值驻点也不是 posterior calibration。


### Lines 46-77

46: ### F-MAT-004 原相关场
47: 
48: 令 \(X_i=p_{i1}+1/2,\;Y_i=p_{i2}+1/2\)，63个非恒定模态为 \(0\le k,l\le7,\;(k,l)\ne(0,0)\)：
49: 
50: \[
51: a_{kl}=\{1+(\pi\ell)^2(k^2+l^2)\}^{-1},\quad \ell=.08,
52: \]
53: \[
54: b_{kl}(i)=c_kc_l\cos(\pi kX_i)\cos(\pi lY_i),\quad
55: c_0=1,\;c_j=\sqrt2\ (j>0),
56: \]
57: \[
58: A=[\,\sigma_0\mathbf1,\ \sigma_f[b_{kl}]\operatorname{diag}(a)/\|a\|\,]
59: \in\mathbb R^{576\times64},\quad f=\mu\mathbf1+Az.
60: \]
61: 
62: 实际 \(\sigma_0=.65,\sigma_f=1,\mu=0\)。这是固定 DCT 相关先验，没有逐对象峰值归一化。
63: 
64: ### F-MAT-005 原 C2 饱和 decoder
65: 
66: \[
67: t(s)=\operatorname{clip}((s+1)/2,0,1),\quad h(s)=10t^3-15t^4+6t^5,
68: \]
69: \[
70: m_i=1_{\|p_i\|\le.38},\quad a_i=m_i h(f_i/\tau),\quad \tau=.5,
71: \]
72: \[
73: \delta=(.8+u_R)+i(.03+.15u_I),\qquad
74: \boxed{\chi_i=\chi_h+\delta a_i},\quad \chi_h=.1+.04i.
75: \]
76: 
77: \(a=0\) 是 host，不是空气；\(a=1\) 是纯 inclusion。\(h\) 在 \(s=\pm1\) 为 C2，但三阶导数一般跳变，跨 seam 的统一 C3 展开没有由此成立。


### Lines 187-255

187: ## 3. 有效材料秩、gauge 和谱
188: 
189: ### F-MAT-012 有效材料截断
190: 
191: \[
192: M(x)=\mathcal R_m(D_x\chi)\in\mathbb R^{1152\times66}
193: =U\Sigma V^T,
194: \]
195: \[
196: \epsilon_M=\max(10^{-8}\sigma_1(M),10^{-10}),\qquad
197: r_M=\#\{\sigma_j(M)>\epsilon_M\}.
198: \]
199: 
200: 这是 R7 material_jets.geometry 和 R9 的实际阈值。较早 decoder_geometry.rank 的默认相对阈值是 \(10^{-10}\)，F 为 \(10^{-7}\)，不能混作同一个数学秩。被截断的小非零奇异值不等于精确 gauge。
201: 
202: ### F-MAT-013 plateau 同材料 gauge
203: 
204: 设 \(\mathcal T=\{i:m_i=1,|f_i|<\tau\}\)。固定物性且 \(A_{\mathcal T}\dot z=0\) 时 transition 场精确不变。只要其余 ROI 点不跨 plateau seam，材料也精确不变。一个保守对称区间为
205: 
206: \[
207: |t|\le\min_{i\notin\mathcal T,\ m_i=1,\ (A\dot z)_i\ne0}
208: \frac{|f_i|-\tau}{|(A\dot z)_i|}.
209: \]
210: 
211: margin 为0时无正区间保证。这是同材料坐标冗余，不能计入非零 material RMS 的散射多解。
212: 
213: ### F-MAT-014 物性径向补偿与 rank
214: 
215: 对 \(\delta\ne0\)、有限 probit 内点、\(a\ne0\)，
216: 
217: \[
218: \operatorname{rank}_{\mathbb R}D_x\chi
219: =\operatorname{rank}B+1+1_{\{a\notin\operatorname{range}B\}}.
220: \]
221: 
222: 证明：\(\delta B\) 贡献 \(r_B\) 维；两物性列形成 \(\{ca:c\in\mathbb C\}\) 的2维实空间；与前者在 \(a\in\operatorname{range}B\) 时交1维，否则交0维。全 host 的 \(a=0\) 应直接用原矩阵，不能套非零公式。
223: 
224: 若试图 \(\dot\delta=c\delta\)，一阶抵消要求 \(B\dot z=-ca\)。不在 range 时只能给残差；即便一阶可消，也没有普遍的有限 gauge，因为 quintic 不满足一般齐次律 \(\alpha h(f/\alpha)=h(f)\)。
225: 
226: ### F-MAT-015 材料单位 lift 与 prior 度量
227: 
228: \[
229: T=V_r\Sigma_r^{-1}\in\mathbb R^{66\times r},\quad
230: Q_m=MT=U_r,\quad Q_m^TQ_m=I_r,\quad
231: K=J_xT\in\mathbb R^{512\times r}.
232: \]
233: 
234: \(T\) 一般不正交，其坐标中的 Gaussian prior 精度是 \(T^TT\)，不是 I。若 r<66，仅在 \(\operatorname{range}V_r\) 采样还删掉原坐标支撑；材料截断几何不是完整后验替代。
235: 
236: ### F-MAT-016 两种广义弱谱
237: 
238: 纯响应材料单位谱为 K 的奇异值：
239: 
240: \[
241: \sigma_{\min}^2(K)
242: =\min_{\substack{v\in\operatorname{range}V_r\\v^TM^TMv=1}}
243: v^TJ_x^TJ_xv.
244: \]
245: 
246: R7/R9 实际初始方向却取
247: 
248: \[
249: A_M=K^TK+T^TT,\quad A_Mw_j=\lambda_jw_j,\quad v_0=Tw_1,
250: \]
251: \[
252: \lambda_1=\min_{\substack{v\in\operatorname{range}V_r\\v^TM^TMv=1}}
253: v^T(J_x^TJ_x+I)v.
254: \]
255: 


## research/maxwell_formula_atlas/05_STATIC_OPERATOR_ERROR_CONTROL.md


### Lines 25-104

25: ## F-ERR-002 · 固定低秩 Green 近似与小系统
26: 
27: 取固定离线近似 $\widetilde C=UV^*$，$U,V\in\mathbb C^{q\times k}$。可以来自截断SVD、物理端口Krylov或其它明确的拟合；此处不预设哪种最好。
28: 
29: 由 Woodbury，若 $I-V^*DU$ 可逆，
30: 
31: $$\widetilde h=(I-DUV^*)^{-1}DE_R
32: =DE_R+DU(I-V^*DU)^{-1}V^*DE_R.$$
33: 
34: 这是真实近似算子的精确解。固定 $U,V$；只有 $k\times k$ 系统依赖材料。其复杂度包含 $O(qk^2+k^3+qkn_s)$，不能只报 $k^3$。低秩近似支持仍是全256格对角材料；它没有改变 decoder 或 prior，但改变了forward，故未经校正不能作为原likelihood。
35: 
36: ## F-ERR-003 · 原 Maxwell 方程的精确缺陷
37: 
38: 定义实际原算子 residual：
39: 
40: $$r=b_D-L_D\widetilde h=D(C-\widetilde C)\widetilde h.$$
41: 
42: 直接代回得
43: 
44: $$h-\widetilde h=L_D^{-1}r,\qquad F-\widetilde F=OL_D^{-1}r.$$
45: 
46: 因此，低秩近似的价值不由 $\|C-\widetilde C\|$ 单独决定，而由 **材料加权、当前场激发、最终接收端口**共同决定。左端 $D$ 与右端 $\widetilde h$ 均不可丢掉。计算真实 $r$ 要付出一次 $C\widetilde h$ 的费用；这仍可能远低于重新LU，但不能记为免费。
47: 
48: ## F-ERR-004 · 可收敛的缺陷校正与精确 fallback
49: 
50: 把 $\widetilde L=I-D\widetilde C$ 当预条件器，
51: 
52: $$h_{j+1}=h_j+\widetilde L^{-1}(b_D-L_Dh_j),\qquad
53:  e_{j+1}=T_De_j,\quad T_D=\widetilde L^{-1}D(C-\widetilde C).$$
54: 
55: 若有可计算上界 $\|T_D\|\le\gamma<1$，则
56: 
57: $$\|e_j\|\le\frac{\|h_{j+1}-h_j\|}{1-\gamma}.$$
58: 
59: 这是由 $e_j=(I-T_D)^{-1}(h_{j+1}-h_j)$ 的Neumann级数得到。非正规系统中只看 $\rho(T_D)<1$ 不能把该范数界替换成同样常数；无法得到有用 $\gamma$ 时，可以使用预条件Krylov及真实residual，或回退完整256LU。**本章没有假设这种迭代一定更快。** 原P1规模小，直接256LU很强，只有整个求解链成本比较才能决定是否值得实施压缩。
60: 
61: ## F-ERR-005 · 接收任务的双边误差表示
62: 
63: 对线性观测任务 $\ell(F)=\operatorname{Re}\operatorname{tr}(W^*F)$，令
64: 
65: $$L_D^*Z=O^*W.$$
66: 
67: 则严格恒等式为
68: 
69: $$\ell(F)-\ell(\widetilde F)=\operatorname{Re}\operatorname{tr}(Z^*r).$$
70: 
71: 若使用近似dual $\widetilde Z$，并有 $\|L_D^{-1}\|\le\beta$，令 $r_d=O^*W-L_D^*\widetilde Z$，则
72: 
73: $$\left|\ell(F)-\ell(\widetilde F)-\operatorname{Re}\operatorname{tr}(\widetilde Z^*r)\right|
74: \le\beta\|r_d\|_F\|r\|_F.$$
75: 
76: 证明来自 $Z-\widetilde Z=L_D^{-*}r_d$ 和Cauchy–Schwarz。dual也随材料变化；固定离线dual只可作近似或预条件。要节省成本，应只为少量实际任务求dual，不为每个材料格点建立巨大免费证书。没有稳定性上界时，该式仍是误差表示，但不是可认证误差区间。
77: 
78: ## F-ERR-006 · 从物理误差传到 likelihood，而非追求机器精度
79: 
80: 设白化预测误差 $\|f_w-\widetilde f_w\|\le\varepsilon$，近似 residual $\widetilde r=\widetilde f_w-d_w$。由三角不等式，
81: 
82: $$\tfrac12\max(0,\|\widetilde r\|-\varepsilon)^2\le\Phi\le
83: \tfrac12(\|\widetilde r\|+\varepsilon)^2.$$
84: 
85: 也有
86: 
87: $$|\Phi-\widetilde\Phi|\le\|\widetilde r\|\varepsilon+\tfrac12\varepsilon^2.$$
88: 
89: 似然为 $L=e^{-\Phi}$，直接反转两端指数即可得 $L$ 包络。与接受阈值或任务积分区间比较后，才能决定是否继续解物理方程。注意 $\varepsilon$ 已包含 $\sqrt2/\sigma$；原始复预测误差小不等于高SNR下likelihood误差小。
90: 
91: ## F-ERR-007 · 固定 receiver 投影遗漏的信息
92: 
93: 设 $S=U\Sigma V^*$，保留前 $p$ 个左奇异方向，$P_U=U_pU_p^*$。
94: 
95: $$\|(I-P_U)F\|_F=\|\Sigma_{tail}V_{tail}^*J_c\|_F
96: \le\sigma_{p+1}(S)\|J_c\|_F.$$
97: 
98: 这是**输出截断**界，不是材料误差界。令复 residual $R=(F-Y)/\sigma\in\mathbb C^{32\times8}$，原512维实白化 residual 为 $r=\sqrt2\operatorname{pack}R$，则
99: 
100: $$\|r\|^2=2\|P_UR\|_F^2+2\|(I-P_U)R\|_F^2.$$
101: 
102: 丢掉第二项就改变了likelihood，除非它与材料无关；即使尾部预测很小，尾部数据噪声与变化预测的交叉项也须处理。对于严格receiver-null current，数据的确不变，但 $J_c=D_m(E^{inc}+KJ_c)$ 仍是所有源共享的硬构成关系，不能省略。
103: 
104: ## F-ERR-008 · 投影式 TSOM 的正确压缩对象


### Lines 112-119

112: ## F-ERR-009 · NN 的数学位置：预测状态，不替换物理目标
113: 
114: 如果未来用网络预测 $\widehat h_\theta(D,E_R)$ 或低维系数，其资格与上述低秩近似相同：
115: 
116: $$r_\theta=b_D-L_D\widehat h_\theta,\qquad
117: h=\widehat h_\theta+L_D^{-1}r_\theta.$$
118: 
119: 可以先用便宜预测，再用真实residual决定继续Krylov或完整回退；也可以由 F-ERR-005/006 控制任务。网络训练误差不是 $\beta$ 或全域置信界，不能把经验误差直接代入严格证书。这样网络承担难以解析表示的状态映射，原likelihood、prior和材料支持保持明确。当前未训练此网络。


## research/maxwell_formula_atlas/06_DERIVED_OPTIMIZATIONS.md


### Lines 5-33

5: ## F-OPT-001 · 固定宿主消元是所有搜索方法共享的强基线
6: 
7: 由 F-OP-008，原非均匀 Full66 支持上的完整预测为
8: 
9: $$L_Dh=DE_R,\qquad L_D=I-DC,\qquad F=F_h+Oh.$$
10: 
11: 固定项是 $C\in\mathbb C^{256\times256}$、$E_R\in\mathbb C^{256\times8}$、$O\in\mathbb C^{32\times256}$、$F_h$。每个材料只更新对角 $D$，解256阶完整系统。它保留原 decoder、噪声模型和 likelihood，所以 fixed-B、dynamic-C、随机探索、普通 Full-J corrector 以及概率积分都应共享它；不能只把这一优化赋给欲证明优越的一臂。
12: 
13: 相比576阶稠密 LU，其 leading factorization 规模由 $N^3$ 变成 $q^3$；这不是实际速度比，因为 setup、BLAS、端口收缩、I/O 和状态数量均有成本。已发生的16个暴露 DEV 状态、3次计时重复的核比较，只支持该窗口的效率判断。它不是新的代数身份，也不是新对象上的完整算法加速。
14: 
15: **直接决策：** 后续方法强基线应为精确256 ROI solver。只有原有瓶颈转移后仍存在的成本，才值得进一步压缩。固定算子低秩近似与真实 residual 校正见 F-ERR-002—006。
16: 
17: ## F-OPT-002 · 多源、多参数与接收端口：三种导数成本不同
18: 
19: 由 F-OP-010，$\delta F=T D_{\delta d}E$，其中 $T=OL_D^{-1}$，本节 $E=E_{ROI}$。同一材料的所有源和所有导数共享一个 $L_D$ factorization。
20: 
21: 1. **完整66列 Jacobian：** 8 forward RHS 得 $E$，32 transpose RHS 得 $T$，再与 $M=\partial d/\partial x$ 收缩。传播费用为40个 reduced RHS，而非 $8\times66$ 个 tangent RHS。
22: 2. **仅一个实标量目标的梯度：** 设 $\Phi=\|F-Y\|_F^2/\sigma^2$，解
23: 
24: $$L_D^*\Lambda=O^*(F-Y)/\sigma^2,\qquad c_i=\sum_s E_{is}\overline{\Lambda_{is}}.$$
25: 
26: 则 $d\Phi=2\operatorname{Re}\sum_i c_i\,\delta d_i$，所以
27: 
28: $$\nabla_x\Phi=2\operatorname{Re}(M^Tc).$$
29: 
30: 只需8个 adjoint RHS；若含标准 Gaussian prior，梯度再加 $x$。这里 $M^T$ 而非 $M^*$ 与定义的 $c$ 配对；也可采用全共轭的等价写法，但不能混配。
31: 3. **矩阵自由 GN 作用：** 一个 $v\mapsto\mathcal J^T\mathcal Jv$ 通常需8 tangent RHS 加8 adjoint RHS，在已有 forward/factorization 之后约16个 reduced RHS。对66维当前问题，若每步迭代次数超过少数次，可能比一次32 RHS构造全 Jacobian更贵；需要计入实际迭代、预条件与收缩成本。
32: 
33: 因此，“矩阵自由”“只更新少数方向”不天然更便宜。应先写清目标到底需要梯度、全 Jacobian，还是某几个作用，再选择算子。不同材料状态通常需要不同 factorization；不同源在同状态则不应重复 factorization。


### Lines 35-92

35: ## F-OPT-003 · 电流残差消元与 variable projection 的精确边界
36: 
37: 原材料目标可以等价写成带硬物理约束的联合问题
38: 
39: $$\min_{m,J_c}\ \|WSJ_c-WY\|_F^2+R(m),\qquad L_mJ_c=D_mE^{inc}.$$
40: 
41: 当 $L_m$ 可逆时，约束唯一决定 $J_c=L_m^{-1}D_mE^{inc}$。消去 $J_c$ 正好返回原 forward 目标；它没有凭空免除求解。
42: 
43: 若改为有限罚参数 $\beta>0$，固定 $m$ 的内层是线性最小二乘：
44: 
45: $$\min_{J_c}\ \left\|\underbrace{\begin{bmatrix}WS\\\sqrt\beta L_m\end{bmatrix}}_{A_\beta(m)}J_c-\underbrace{\begin{bmatrix}WY\\\sqrt\beta D_mE^{inc}\end{bmatrix}}_{b_\beta(m)}\right\|_F^2.$$
46: 
47: 在 $L_m$ 可逆时 $A_\beta$ 满列秩，故
48: 
49: $$J_*(m)=A_\beta^\dagger b_\beta,\qquad \psi_\beta(m)=\|(I-A_\beta A_\beta^\dagger)b_\beta\|_F^2.$$
50: 
51: 这是标准 variable projection。令 $r_\beta=A_\beta J_*-b_\beta$，由 $A_\beta^*r_\beta=0$，微分时含 $dJ_*$ 的项消去：
52: 
53: $$d\psi_\beta=2\operatorname{Re}\operatorname{tr}\{r_\beta^*[(dA_\beta)J_*-db_\beta]\}.$$
54: 
55: 这给出准确的 reduced 梯度，并解释为何联合变量的内层可消元。但有限 $\beta$ 允许不满足 Maxwell 的电流；它改变了 likelihood。$\beta\to\infty$ 在适当条件下逼近硬约束，同时可能恶化条件数。用它产生初值后再做原物理验证，可以作为优化流程；不能把有限罚目标直接称为原 posterior。此处未主张它比256阶精确求解更快。
56: 
57: 实现上优先 QR/SVD 或合适的增广求解，而不是无条件形成 $A_\beta^*A_\beta$；后者平方条件数。若原材料用66 latent 表示，增加每源256复 current 会引入4096实状态变量，变量数上未必占优。
58: 
59: ## F-OPT-004 · 能否只解接收可见的少量端口
60: 
61: 为说明障碍，令固定接收的完整正交基为 $[V_o,V_n]$，$V_o\in\mathbb C^{q\times r}$，$V_n\in\mathbb C^{q\times(q-r)}$，$OV_n=0$。写 $h=V_oa+V_nb$，原约束为
62: 
63: $$L_DV_oa+L_DV_nb=DE_R.$$
64: 
65: 因为 $L_D$ 可逆，$B_D=L_DV_n$ 满列秩。取其左零空间正交基 $U_D\in\mathbb C^{q\times r}$，即 $U_D^*B_D=0$。乘上 $U_D^*$ 消去 $b$：
66: $$\underbrace{U_D^*L_DV_o}_{T_D^{port}}a=U_D^*DE_R.$$
67: 
68: $T_D^{port}$ 可逆。若 $T_D^{port}a=0$，则 $L_DV_oa\in\operatorname{range}(L_DV_n)$，由 $L_D$ 可逆得 $V_oa\in\operatorname{range}V_n$，正交性迫使 $a=0$。求得 $a$ 后
69: 
70: $$b=B_D^\dagger(DE_R-L_DV_oa),\qquad F=F_h+OV_oa.$$
71: 
72: 这是精确端口消元，但 $U_D$ 随材料变化；每次为 $q\times(q-r)$ 矩阵求 QR/零空间可能花费与原 LU 同阶的工作。把最终 $r\times r$ 方程很小当作整体计算很小，是漏记消元成本。固定 $CV_o,CV_n$ 可以缓存，$D$ 的行缩放也很便宜，但不能因此免费得到变化的 $U_D$。
73: 
74: 另一种块 Schur 方法需要 $V_n^*L_DV_n$ 可逆，这是比 $L_D$ 可逆更强的假设：矩阵 $\begin{bmatrix}0&1\\1&0\end{bmatrix}$ 可逆，而任一单坐标对角块均为零。不可未经检查使用该块逆。
75: 
76: 真正免费的固定闭合需要额外结构，例如 F-OP-006 的 $P_oDCP_n=0$；原非均匀材料族没有给出这种闭合。**结论：接收 rank 小是输出维数事实，不能单独证明可同样便宜地求解非线性状态。**
77: 
78: ## F-OPT-005 · 精确可达子空间与近似端口压缩
79: 
80: 均匀幅度 $D=vI$ 时，固定可达子空间
81: 
82: $$\mathcal R=\operatorname{span}\{E_R,CE_R,\ldots,C^{q-1}E_R\}$$
83: 
84: 包含所有源列，并由 Cayley–Hamilton 定理对 $C$ 不变。若 $Q_c$ 为其正交基，$C Q_c=Q_cC_c$、$E_R=Q_cB_c$，则
85: 
86: $$h=Q_c\,v(I-vC_c)^{-1}B_c,\qquad F=F_h+OQ_c\,v(I-vC_c)^{-1}B_c.$$
87: 
88: 这是精确可达部分消元；若可达 rank 仍为 $q$，没有降维收益。进一步删除不可观测状态需要正确的不变子空间商或最小实现条件，不能仅把 $\ker O$ 删掉。有限阶 Krylov 截断通常只是矩匹配近似，除非已达到完整不变闭包。
89: 
90: 非均匀材料则需一个对所有允许 $D C$ 共同不变、且包含所有 $DE_R$ 的固定空间。受限 decoder 也许有特殊结构，但“latent仅66维”本身并不证明这种不变性；对角材料变化可激发许多独立空间方向。其一般障碍见 F-ERR-001。
91: 
92: **实际优化路线：** 若没有精确小闭包，只把固定压缩作为预条件/初值，计算真实 residual，再依 F-ERR-003—006 校正或回退。判据应是接收或任务误差以及总费用，不是单独的 $C$ 奇异值尾。当前256阶系统已经很小，不预设低秩路线必胜。


### Lines 142-176

142: ## F-OPT-009 · 与材料兼容的电流投影
143: 
144: 在合法材料 $x$，所有可实现的一阶 current 变化来自
145: 
146: $$\delta h_s=L_D^{-1}\operatorname{diag}(E_{:,s})M\,\delta x.$$
147: 
148: 按实部/虚部堆叠得到 $\mathcal T\in\mathbb R^{2qn_s\times66}$，即 $\delta h_{real}=\mathcal T\delta x$。任意 current proposal $u$ 若要一阶材料兼容，应投到 $\operatorname{range}\mathcal T$，而不是直接接受接收暗空间中的任意向量。Euclidean 正交投影为 $\mathcal T\mathcal T^\dagger$；若 current metric 为 $W_c\succ0$，对应
149: 
150: $$P_{\mathcal T}^{W_c}=\mathcal T(\mathcal T^TW_c\mathcal T)^\dagger\mathcal T^TW_c.$$
151: 
152: 该式给出加权投影，但实现宜对白化后的 $W_c^{1/2}\mathcal T$ 做 QR/SVD。它一般不是普通 Euclidean 对称投影。
153: 
154: 要同时一阶材料合法且数据不变，应先解 $\mathcal J\delta x=0$，再映射 $\mathcal T\delta x$；其空间是 $\operatorname{range}\mathcal T\cap\ker A_{obs}$，其中 $\mathcal J=A_{obs}\mathcal T$。材料投影与接收暗投影通常不对易，因此两次串联不等于该交集。若 $\mathcal J$ 在实 latent 空间满列秩，就没有非零 latent 一阶精确数据零方向；小奇异值只给近似弱方向。
155: 
156: 构造完整 $\mathcal T$ 本身可能比已需的512×66观测 Jacobian贵得多。因此该式主要用于合法性解释；不建议为了“做一个投影”无条件构造巨大 current tangent。已知 decoder 的 latent proposal 加真实 corrector，通常更直接地维持原材料支持。
157: 
158: ## F-OPT-010 · 嵌套缓存与错误的低秩推断
159: 
160: 可按依赖顺序缓存：固定几何 $K,S,E^{inc}$ → 固定宿主 $C,E_R,O,F_h$ → 当前 $D,L_D$ factorization → 当前 $h,E,T$ → decoder 导数与 GN 小矩阵。该层次明确指出缓存何时失效。
161: 
162: 材料增量若只改变 $k$ 个格点，$\Delta D\,C$ 的秩至多 $k$，可以考虑对应 Woodbury 更新；但 latent维数66不意味着 $\operatorname{rank}(\Delta D C)\le66$。一个 latent方向可能同时改变全部256个格点，使对角 $\Delta D$ 满秩。低维参数流形与每次矩阵更新低秩是不同性质。
163: 
164: 完整 current 的重建 $J_c=J_h+V_hh$ 只在物理诊断或输出确实需要时做；接收预测与 ROI material Jacobian 不必每次重建576格场。若合同要求某项完整诊断，其成本仍必须保留。各源、各receiver和各 latent 参数应共享现成 factorization；重复 setup、重复求源和遗漏回退均会扭曲成本。
165: 
166: ## F-OPT-011 · 由数学直接裁决的优先级
167: 
168: |候选|代数能确认的事实|仍需证据的事实|当前建议|
169: |---|---|---|---|
170: |完整固定ROI消元|原模型完全等价，全部非均匀支持保留|实际全流程费用|作为各臂公共强基线|
171: |8源共享因子、32receiver端口导数|完整66列导数无需528 tangent RHS|收缩和I/O的实际比例|直接采用正确结构|
172: |完整稳定基 QR corrector|保留原GN目标，避免正规方程条件数平方|是否改进具体求解器稳定性|普通稳定实现，不作为新理论|
173: |均匀幅度固定Schur|精确三角复用|特定任务是否值得复用|只用于真正均匀的条件模型|
174: |固定receiver SVD/TSOM表示|固定投影可组织current；截断有明确改变|整体速度与合法材料覆盖|不把表示优势冒称完整物理降维|
175: |current variable projection|内层可消元；有限罚目标不同|是否帮助难优化与总费用|低优先级，先与精确ROI强基线比较|
176: |低秩C、端口Krylov、网络初值|可接真实residual与fallback|误差界是否有用、是否摊销|见第05章，不先开rank sweep|


## research/maxwell_formula_atlas/09_BILATERAL_SHAPE_AMPLITUDE.md


### Lines 5-77

5: 统一尺寸：ROI 阶数 $q=256$，接收端 $n_r=32$，源端 $n_s=8$。固定宿主消元后的 $C\in\mathbb C^{q\times q}$、$O\in\mathbb C^{n_r\times q}$、$E_R\in\mathbb C^{q\times n_s}$、$F_h$ 均固定；$D$ 是相对宿主的材料增量对角矩阵。时间因子、白化和原先验沿用 F-BASE-001—003。本章只对冻结离散系统作结论。
6: 
7: ## F-BI-001 · 双边非线性传递
8: 
9: 由 $h=D(E_R+Ch)$ 消去电流得到
10: 
11: $$\mathcal T(D)=(I-DC)^{-1}D=D(I-CD)^{-1},\qquad F(D)=F_h+O\mathcal T(D)E_R.$$
12: 
13: 推导使用 $(I-AB)^{-1}A=A(I-BA)^{-1}$，只要求对应 resolvent 可逆，不要求 $D$ 可逆。接收端和输入端的维数小，不意味着中间状态在所有材料下都可以无损缩小。**资格：标准离散 T-operator 恒等式。** 实现：R9 精确 ROI solver；R10 `src/transfer_campaign.py`。
14: 
15: ## F-BI-002 · 有限材料差的双边条件
16: 
17: $$\mathcal T(D_2)-\mathcal T(D_1)=(I-D_2C)^{-1}(D_2-D_1)(I-CD_1)^{-1}.$$
18: 
19: 在等式左乘 $I-D_2C$、右乘 $I-CD_1$ 后，两侧均化为 $D_2-D_1$，给出直接证明。因此等观测要求
20: 
21: $$O(I-D_2C)^{-1}\Delta D(I-CD_1)^{-1}E_R=0.$$
22: 
23: 这是有限非线性响应等价的必要且充分条件，但未知 $D_2$ 仍在左侧 resolvent 中；它没有自动变成固定线性核约束。材料必须仍属于原 decoder 支持。**计算含义：** 允许构造有限响应和误差分析；不能据此宣布找到 exact fiber。
24: 
25: ## F-BI-003 · 固定形状的精确加权实现
26: 
27: 原 P1 材料增量 $D=\delta A(z)$，其中 $A=\operatorname{diag}(a_i)$，$0\le a_i\le1$，$\delta\in[.8,1.8]+i[.03,.18]$。令
28: 
29: $$P=A^{1/2},\quad B=PCP,\quad L=OP,\quad R=PE_R.$$
30: 
31: 则
32: 
33: $$F(z,\delta)-F_h=\delta L(I-\delta B)^{-1}R.$$
34: 
35: 证明：$\mathcal T(\delta P^2)=\delta P(I-\delta PCP)^{-1}P$，可直接左右乘原方程验证，完全不需要 $P^{-1}$。所以零占据率 cell 不造成除零。**资格：精确离散恒等式，保留所有多次散射。** 形状改变后 $B,L,R$ 也改变；仅复幅度查询可以共用同一分解。
36: 
37: ## F-BI-004 · 传递矩、可达与可观测状态
38: 
39: 在 $|\delta|\rho(B)<1$ 的收敛域内，
40: 
41: $$F-F_h=\sum_{j\ge0}\delta^{j+1}LB^jR.$$
42: 
43: 系数 $LB^jR$ 是输入输出矩。有限维实现的可达空间为 $\mathcal R=\operatorname{span}\{R,BR,\ldots,B^{q-1}R\}$；不可观测空间为 $\mathcal N=\cap_{j=0}^{q-1}\ker(LB^j)$。最小实现阶数是 $\dim\mathcal R-\dim(\mathcal R\cap\mathcal N)$，也等于充分长块 Hankel 矩阵的秩。这里的交集是理论定义，算法可用 QR/SVD 的投影实现，不要求显式集合求交。
44: 
45: 解析延拓可将同一有理响应用于级数收敛域外的非极点位置；不能把级数截断在该域外当成有效近似。R10 的 64 阶结果只是被测试的压缩阶数，未认证最小阶数。**来源：标准有限维实现与 Krylov 理论；最近材料模态先例见 R10 查新表。**
46: 
47: ## F-BI-005 · 固定宿主的双边材料 Gram
48: 
49: 宿主处一阶变化为 $dF=O\,dD\,E_R$，第 $i$ 个材料列 $J_{(r,s),i}=O_{ri}(E_R)_{is}$，故
50: 
51: $$J^*J=(O^*O)\circ\overline{E_RE_R^*}.$$
52: 
53: 证明逐项求和：$\sum_{r,s}\overline{O_{ri}E_{is}}O_{rj}E_{js}=(O^*O)_{ij}\sum_s\overline{E_{is}}E_{js}$。这比仅 receiver SVD 同时保留了输入照明，但仍是宿主一阶信息。若 $D_x=\partial d/\partial x\in\mathbb C^{q\times66}$，真实高斯 latent 上的白化 Gram 为
54: 
55: $$H_x=\frac2{\sigma^2}\operatorname{Re}\{D_x^*(J^*J)D_x\}.$$
56: 
57: **条件：** 参数是实坐标，材料按实 Hilbert 度量；复 Gram 不得直接当成实 Hessian。R10 的谱更贴近 Full-J 并未转化为更好的 GT 覆盖。
58: 
59: ## F-BI-006 · 完整非线性灵敏度的同型结构
60: 
61: 令 $L_D=O(I-DC)^{-1}$，$R_D=(I-CD)^{-1}E_R$。由 F-BI-002 的极限有
62: 
63: $$dF=L_D\,dD\,R_D,\qquad J_D^*J_D=(L_D^*L_D)\circ\overline{R_DR_D^*}.$$
64: 
65: 固定两端宿主 Gram 与当前 Full-J 的差别恰在两侧内部反馈。它解释何时 Born 双边坐标失配，但不允许免费得到状态依赖的两侧场。**资格：标准 distorted-Born 导数结构。**
66: 
67: ## F-BI-007 · 共同状态压缩的条件性障碍
68: 
69: 考虑允许所有独立对角 $D$ 的材料族。若输入矩阵 $E_R$ 每行非零，则 $\operatorname{span}_{D}\{DE_R\}=\mathbb C^q$：用只保留第 $i$ 个对角元的矩阵即可产生第 $i$ 个坐标方向。因而不存在一个小于 $q$ 的固定状态子空间包含所有这些注入。
70: 
71: 这只是充分反例：它限制所有材料通用的精确状态压缩，不否定固定形状的压缩，也不直接适用于受限的 P1 Gaussian64 decoder 像。更小的输入输出商仍要另审可观测性，不能由注入满张成单独排除。
72: 
73: ## F-BI-008 · 交换子不是端口误差证书
74: 
75: $$[A,C]_{ij}=(a_i-a_j)C_{ij},\qquad \|[A,C]\|_F^2=\sum_{ij}|a_i-a_j|^2|C_{ij}|^2.$$
76: 
77: 它衡量交换性；$[A,C]=0$ 仅说明 $C$ 不混合 $A$ 的不同特征值子空间，不保证共同对角化或低秩（例如 $A=I$、$C$ 为 Jordan 块）。大交换子可能主要落在端口不敏感方向。实际观测误差还需要 $L,R$、resolvent 稳定性和被舍弃状态的响应。**资格：恒等式加反推限制，不是新的普适误差界。**


### Lines 124-146

124: ## F-BI-014 · 非正规矩阵的损耗稳定性界
125: 
126: 令 $K_\delta=\delta^{-1}I-B$、$\operatorname{Im}_H B=(B-B^*)/(2i)$，设 $\beta_B=\lambda_{min}(\operatorname{Im}_H B)$，则
127: 
128: $$-\operatorname{Im}_H K_\delta\succeq b_\delta I,\quad b_\delta=\frac{\operatorname{Im}\delta}{|\delta|^2}+\beta_B.$$
129: 
130: 若 $b_\delta>0$，由 $|v^*K_\delta v|\ge b_\delta\|v\|^2$ 和 Cauchy–Schwarz 得
131: 
132: $$\sigma_{min}(K_\delta)\ge b_\delta,\qquad \|K_\delta^{-1}\|_2\le b_\delta^{-1}.$$
133: 
134: 不要求 $B$ 正规，也没有假设 $\operatorname{Im}_HB\succeq0$。R10 的最小特征值实际略负，材料损耗仍使 $b_\delta$ 正。float64 特征值仅给数值资格；没有 outward rounding 就不能称严格机器证书。
135: 
136: ## F-BI-015 · 全状态 residual 到 likelihood 误差
137: 
138: 对近似解 $\tilde u$，$r=R-K_\delta\tilde u$，真实输出误差为 $LK_\delta^{-1}r$，所以
139: 
140: $$e:=\|W_\sigma(F-\tilde F)\|\le\frac{\sqrt2\|L\|_2\|r\|_F}{\sigma b_\delta}.$$
141: 
142: 代码可用较保守但便宜的 $\|L\|_F$ 代替 $\|L\|_2$。若 $\tilde r_y=W_\sigma(\tilde F-y)$，
143: 
144: $$|\Phi-\tilde\Phi|\le\|\tilde r_y\|e+\tfrac12e^2.$$
145: 
146: 这是从同一物理 residual 到任务判断的桥梁。R10 冻结门为 $e\le.01$ 且势误差上界 $\le.1$；不合格才 exact fallback。门只覆盖被评价的节点，不凭有限节点推断全先验支持。


## research/maxwell_formula_atlas/10_TSOM_ON_BILATERAL_TRANSFER.md


### Lines 5-87

5: 实质文献依据为作者公开的 [Xu、Zhong、Song、Chen、Ran，MR-FFT-TSOM，DOI 10.1109/TGRS.2014.2329032](https://renchengsong.github.io/papers/TGRS2329032.pdf)，§II，式(3)—(11)。该文明确重述接收 SVD、内部传播基的投影及数据/状态双残差。以下符号吸收了原文电流中的频率常数；原论文算法、这里的噪声修正、variable projection 和低秩正演必须分开归属。没有声称已逐式读到2009原篇。
6: 
7: ## F-BI-019 · 从双边逆算子还原电流方程
8: 
9: 令 $N=256$，$n_s=8$，$h\in\mathbb C^{N\times n_s}$ 为固定宿主消元后的 ROI 对比电流，$y=Y-F_h$。定义 $K_D=I-DC$。当 $K_D$ 可逆时：
10: 
11: $$h=K_D^{-1}DE_R\quad\Longleftrightarrow\quad h=D(E_R+Ch),\qquad F-F_h=Oh.$$
12: 
13: 这是精确等价式。原完整正演只需解8个右端，并不显式构造256×256逆矩阵或整个 $\mathcal T(D)$。SOM/TSOM 的入口是将 $h$ 作为辅助未知量，以固定投影组织它，同时满足上述状态方程与数据方程。它不要求 $D$ 与 $C$ 对易，也不要求真实材料解集线性。
14: 
15: ## F-BI-020 · 两次筛选是投影像，不是交集
16: 
17: 设 $S=WO$，其中 $W$ 为已固定的逐源噪声白化；$S=U\Sigma V^*$。保留噪声下可靠的前 $\ell$ 个方向：
18: 
19: $$h_d=V_\ell\Sigma_\ell^{-1}U_\ell^*Wy,\qquad P_\ell=V_\ell V_\ell^*.$$
20: 
21: 令 $V_{C,r}$ 为固定内部算子 $C$ 的所选右奇异向量，或使用明确声明的 Fourier 基 $V_{F,r}$。TSOM式表示为
22: 
23: $$Z=(I-P_\ell)V_{C,r},\qquad Q=\operatorname{orth}(Z),\qquad h=h_d+Q\alpha.$$
24: 
25: 这里 $\operatorname{orth}$ 删除相关/零列，保留相同列空间。$\operatorname{range}Z$ 是投影后的像，一般不等于两空间交集；$ZZ^\dagger$ 才是该像的正交投影。截断接收强方向后，$OQ$ 通常非零，不能省略数据残差。原作者的 twofold 内部基/Fourier构造支持这个位置；它不是从材料 Jacobian 的 weak8 方向构造电流。
26: 
27: **固定缓存：** $Q,CQ,OQ,V_\ell,CV_\ell$；同一观测的 $h_d,Ch_d$。宿主消元后的 $C$ 一般不再具有完整均匀网格的卷积形式，Fourier作为基仍合法，但不能未经推导声称 $C$ 乘法自动具有FFT卷积加速。
28: 
29: ## F-BI-021 · 共享材料与含噪数据的联合目标
30: 
31: 原 P1 材料保持 $D=D(x)=\operatorname{diag}(\delta(x)a(x))$，$x\in\mathbb R^{66}$。不得逐源恢复不同材料，或用自由像素代替原 decoder。为不把含噪接收强系数硬锁在 $h_d$ 上，可加入其修正：
32: 
33: $$R=[V_\ell,Q],\qquad h=h_d+R\beta,\qquad \beta\in\mathbb C^{q\times n_s}.$$
34: 
35: 采用明确的候选生成目标（归一化权重需单独冻结）：
36: 
37: $$\mathcal J_\lambda(x,\beta)=\|S(h_d+R\beta)-Wy\|_F^2+\lambda\|K_{D(x)}(h_d+R\beta)-D(x)E_R\|_F^2+\mathcal R(x),\quad\lambda>0.$$
38: 
39: 对全部源共享同一个 $D(x)$；$\mathcal R$ 可为原先验惩罚或事先固定的随机化目标。没有 $N$ 阶逆出现在目标计算中：$K_DR=R-D(CR)$，因此固定的 $CR$ 可预计算，材料更新主要是对角行缩放。
40: 
41: 有限 $\lambda$ 允许状态残差，所以这是优化辅助目标，不等于原 Maxwell likelihood。作为材料可能性生成器完全可以使用；候选需由原完整正演验证。不要求先解决精确 posterior 才能研究这个生成器。若把接收强系数固定不动，则应将其标为传统截断数据近似，而非原含噪模型的精确等价。
42: 
43: ## F-BI-022 · 电流系数的小规模消元
44: 
45: 对固定 $D$，将上一目标的两项堆叠：
46: 
47: $$T_\lambda(D)=\begin{bmatrix}SR\\\sqrt\lambda K_DR\end{bmatrix},\qquad B_\lambda(D)=\begin{bmatrix}Wy-Sh_d\\\sqrt\lambda(DE_R-K_Dh_d)\end{bmatrix}.$$
48: 
49: 若 $R$ 满列秩且 $K_D$ 可逆，则 $T_\lambda$ 满列秩：其下块已经满列秩。于是
50: 
51: $$\beta_*(D)=T_\lambda(D)^\dagger B_\lambda(D),\qquad \min_\beta\mathcal J_\lambda=\|(I-T_\lambda T_\lambda^\dagger)B_\lambda\|_F^2+\mathcal R(x).$$
52: 
53: 采用薄 QR/SVD 求解，不应为了书写方便在实现中显式形成 $(T^*T)^{-1}$。若 $q\ll N$，主要内层费用约为 $O((N+n_r)q^2+q^3)$，多个源共享同一QR；不再需要每一步 $N\times N$ LU。若 $q$ 很大或需要大量恢复迭代，收益可能消失。
54: 
55: 令 $r=T_\lambda\beta_*-B_\lambda$，正交性 $T_\lambda^*r=0$ 给出
56: 
57: $$d\|r\|_F^2=2\operatorname{Re}\operatorname{tr}\{r^*[(dT_\lambda)\beta_*-dB_\lambda]\}.$$
58: 
59: 因此约化目标求导不用显式求 $d\beta_*$，但它仍是上述有限罚目标的导数，不是原 Maxwell posterior 的梯度。这是标准 variable projection；是否比精确256阶Full-J更好需要真正的同质量成本比较。
60: 
61: ## F-BI-023 · 若目的是直接压缩中间逆：低秩 Green 与 Woodbury
62: 
63: 另一条路线固定近似 $C_r=AB^*$，$A,B\in\mathbb C^{N\times r}$（例如截断SVD因子）。令 $g=DE_R$。在 $I-B^*DA$ 可逆时：
64: 
65: $$\widetilde h=(I-DC_r)^{-1}g=g+DA(I-B^*DA)^{-1}B^*g.$$
66: 
67: 这是标准 Woodbury 恒等式，对近似算子 $C_r$ 精确。对每个新材料仍需更新 $B^*DA$ 并解 $r$ 阶系统，但不用重新分解 $N$ 阶矩阵。若 $C=C_r$ 是实际精确低秩，则也是原离散物理的精确表达；一般截断情况下只是近似。
68: 
69: 这与 TSOM 有共同的固定传播压缩哲学，但算法对象不同：TSOM主要减少反演中未知电流系数；此式直接压缩正演求解器。不能用“把 $C$ 分解一次”推出“所有 $D$ 下的逆都已分解一次”。同样不要求 $D$ 与 $C$ 对易。
70: 
71: ## F-BI-024 · 完整物理缺陷给出投影的失效位置
72: 
73: 对任意近似电流 $\widetilde h$ 定义真实缺陷
74: 
75: $$r_D=DE_R-K_D\widetilde h.$$
76: 
77: 当 $K_D$ 可逆，精确误差恒为
78: 
79: $$h-\widetilde h=K_D^{-1}r_D,\qquad F-\widetilde F=OK_D^{-1}r_D.$$
80: 
81: 在 F-BI-023 的低秩算子求解中，进一步有
82: 
83: $$r_D=D(C-C_r)\widetilde h.$$
84: 
85: 这说明被删传播的实际影响同时受材料、内部场和resolvent放大控制；只看固定Green的奇异值尾不够。若有稳定性界 $\|K_D^{-1}\|\le b_D^{-1}$，则 $\|W(F-\widetilde F)\|_F\le\|WO\|_2\|r_D\|_F/b_D$；无可靠界时用完整校正/验证。
86: 
87: 对纯TSOM，缺陷可能来自电流子空间不足、有限罚优化未收敛或材料与电流失配，不能全部叫“流形不稳定”。误差超门时先考虑扩大固定电流基或精确校正，不必默认更新材料弱方向。


## research/maxwell_formula_atlas/11_ACTIVE_SUPPORT_TOTAL_FIELD_AND_GLOBAL_BOUNDS.md


### Lines 5-109

5: 统一尺寸：$n=256$，$n_r=32$，$n_s=8$。$C\in\mathbb C^{n\times n}$、$O\in\mathbb C^{n_r\times n}$、$E_R\in\mathbb C^{n\times n_s}$、$F_h\in\mathbb C^{n_r\times n_s}$ 固定；$a\in[0,1]^n$、$A=\operatorname{diag}(a)$、$D=\delta A$ 改变。$\delta$ 在原被动物性矩形内，$\sigma$ 和 source-major 白化沿用 F-BASE-002。$^*$ 表示共轭转置；解析导数转置求解与 Hermitian 内积不可混同。
6: 
7: ## F-SS-001 · 严格零支撑的精确电流消元
8: 
9: 设 $\mathcal S=\{i:a_i>0\}$，$\mathcal I$ 为其补集。由 $h=D(E_R+Ch)$ 逐行得 $h_\mathcal I=0$，故
10: 
11: $$M_S=I_S-\delta A_SC_{SS},\quad h_S=M_S^{-1}\delta A_SE_{R,S},\quad F=F_h+O_Sh_S.$$
12: 
13: 按 $S,I$ 排序，$I-DC$ 是 $\begin{bmatrix}M_S&-D_SC_{SI}\\0&I_I\end{bmatrix}$，因此它与 $M_S$ 的可逆性等价。只删**精确零**；小正数删减是另一项近似。全部源共享一次分解。空支持时 $F=F_h$。资格：标准精确块消元；实现 `agent_a/solver.py:Solver.evaluate`。
14: 
15: ## F-SS-002 · 完整端口灵敏度的活动/非活动分块
16: 
17: 定义 $T=O(I-DC)^{-1}$ 和 $Z=O_SM_S^{-1}$。对上述块三角矩阵求逆直接得
18: 
19: $$T_{:,S}=Z,\qquad T_{:,I}=O_I+ZD_SC_{SI},\qquad e=E_R+C_{:,S}h_S.$$
20: 
21: 只需活动系统的 $n_s$ 个正演和 $n_r$ 个转置右端；不需要另做全 $n$ 阶 LU。即使材料有零格点，仍要保留全部 $n$ 个材料导数列。资格：标准逆矩阵分块；不能把活动电流状态数当材料参数数。
22: 
23: ## F-SS-003 · 支撑出生的非零材料导数
24: 
25: 由 $(I-DC)\,dh=(dD)e$，得到
26: 
27: $$dF=T(dD)e,\qquad \frac{\partial F_{rs}}{\partial a_i}=\delta T_{ri}e_{is},\qquad \frac{\partial F}{\partial\delta}=T A e.$$
28: 
29: 当 $a_i=0$ 时，$h_i=0$，但 $\partial F/\partial a_i$ 一般非零。因此原生占比优化可以用这些列激活新支持。对原 P1 latent 则须再乘 $Da/Dz$；平台内该因子可为零，不能由本式推出 latent 梯度必能出生新笔画。实虚物性导数分别为 $F_\delta$ 与 $iF_\delta$，Gaussianized 属性还须乘正态 CDF 导数。资格：精确一阶导数，不是新测量信息。
30: 
31: ## F-SS-004 · 不随支持改变的总场注入
32: 
33: 从 $e=E_R+Ch$、$h=\delta Ae$ 得
34: 
35: $$M_e e=E_R,\quad M_e=I-\delta CA,\quad F=F_h+\delta OA e.$$
36: 
37: $E_R$ 始终固定，避免旧加权状态 $u=\delta\sqrt A e$ 的注入 $\sqrt A E_R$ 随形状改变。这里不包含对低维共同不变空间的任何保证。资格：精确状态变换；零 $a_i$ 无需除法。
38: 
39: ## F-SS-005 · 固定总场 Galerkin 与真实缺陷
40: 
41: 取固定正交基 $V\in\mathbb C^{n\times r}$，例如 $\operatorname{orth}[E_R,CA_0E_R,\ldots]$ 或 $\operatorname{orth}[E_R,CE_R,\ldots]$。每个当前材料解
42: 
43: $$(I_r-\delta V^*CAV)t=V^*E_R,\quad\widetilde e=Vt,\quad\widetilde F=F_h+\delta OAVt,$$
44: 
45: $$r_e=E_R-(I-\delta CA)\widetilde e.$$
46: 
47: 固定基不等于固定 $V^*CAV$。每次形状变化仍需更新 $CAV$、$OAV$ 和投影；其成本必须计入。非正规 Galerkin 的稳定性不由基正交自动保证。资格：标准近似；实现 `agent_b/fixed_total_field.py`。
48: 
49: ## F-SS-006 · 全部占比形状的 Hermitian 虚部下界
50: 
51: 令 $H_C=(C-C^*)/(2i)$、$P=\sqrt A$、$B=PCP$，有 $\operatorname{Im}_H B=PH_CP$。记 $\lambda=\lambda_{\min}(H_C)$，$\beta_*=\min(0,\lambda)$。对任意 $v$，
52: 
53: $$v^*PH_CP v\ge\lambda\|Pv\|^2\ge\beta_*\|v\|^2,$$
54: 
55: 最后一步在 $\lambda<0$ 时用 $\|Pv\|\le\|v\|$，在 $\lambda\ge0$ 时直接用非负性。因此 $\lambda_{\min}(\operatorname{Im}_H B)\ge\beta_*$。资格：条件标准矩阵序不等式；不能把负 $\lambda$ 随意抹为零。
56: 
57: ## F-SS-007 · 被动物性矩形上的统一强制性常数
58: 
59: 对 $\delta=x+iy$，本合同 $0<y_{min}\le y\le y_{max}<x_{min}\le x\le x_{max}$。$y/(x^2+y^2)$ 随 $x$ 递减、随 $y$ 递增，故
60: 
61: $$b_*={y_{min}\over x_{max}^2+y_{min}^2}+\beta_*.$$
62: 
63: 当 $b_*>0$，令 $K_\delta=\delta^{-1}I-B$，则 $-\operatorname{Im}_H K_\delta\succeq b_*I$。Cauchy–Schwarz 给
64: 
65: $$b_*\|v\|^2\le|v^*K_\delta v|\le\|v\|\,\|K_\delta v\|\Longrightarrow\|K_\delta^{-1}\|_2\le b_*^{-1}.$$
66: 
67: 若 $b_*\le0$，本充分界失效，不能改用经验正数。本合同有限浮点结果 $\lambda=-4.5739299047\times10^{-5}$，$b_*=0.009210948658>0$；不是向外舍入的连续 PDE 证书。
68: 
69: ## F-SS-008 · 总场缺陷到端口误差的精确矩阵顺序
70: 
71: 由 push-through 恒等式，$P K_\delta^{-1}P=\delta A(I-\delta CA)^{-1}$。因此
72: 
73: $$F-\widetilde F=OP K_\delta^{-1}P r_e.$$
74: 
75: 这个式子在 $A$ 奇异时仍成立，不需要 $P^{-1}$。直接展开 Neumann 多项式并由可逆区域延拓，或用 $(I-AB)^{-1}A=A(I-BA)^{-1}$，均可得到。资格：离散精确恒等式，不需要 $C,A$ 对易。
76: 
77: ## F-SS-009 · 全形状白化端口误差界
78: 
79: 将前两式取范数并用 source-major 等距实化，有
80: 
81: $$\|W(F-\widetilde F)\|_2\le\epsilon_W={\sqrt2\over\sigma b_*}\|OP\|_2\|P r_e\|_F,$$
82: 
83: $$\|OP\|_2\le\|O\|_2.$$
84: 
85: 当前 $\|OP\|_2$ 可用 $32\times32$ 行 Gram 的最大本征值求得；全统一 $\|O\|_2$ 更便宜但更保守。资格：在上述耗散假设下的矩阵范数界；实现是 float64 数值条件界。
86: 
87: ## F-SS-010 · likelihood 势误差与回退门
88: 
89: $\Phi(F;y)=\frac12\|W(F-y)\|^2$。若 $\|W(F-\widetilde F)\|\le\epsilon_W$，则
90: 
91: $$|\Phi(F;y)-\Phi(\widetilde F;y)|\le\|W(\widetilde F-y)\|\epsilon_W+\tfrac12\epsilon_W^2=\epsilon_\Phi.$$
92: 
93: 接受门为 $\epsilon_W\le0.01$ 且 $\epsilon_\Phi\le0.1$；不通过有限增阶后精确回退。实际配对误差只能做事后评估，不能反过来替代在线界。这是误差控制，不是 posterior 认证。
94: 
95: ## F-SS-011 · 近似物理导数不自动继承响应误差门
96: 
97: 对精确 $T$、$e$，$J_i=\delta T_{:,i}e_{i,:}$。若用近似 $\widetilde T,\widetilde e$，则
98: 
99: $$J_i-\widetilde J_i=\delta[(T-\widetilde T)_{:,i}e_{i,:}+\widetilde T_{:,i}(e-\widetilde e)_{i,:}].$$
100: 
101: 只有响应残差门并不控制全部 Jacobian。R11 的 B/C 在需要 Full-J 时使用精确一致的 F/J；只对 forward trial 尝试降阶。因此不能宣称降阶同时省了32个导数右端。资格：精确误差分解；受控降阶导数尚未实现。
102: 
103: ## F-SS-012 · 三类成本与图像收益不可互换
104: 
105: 对于相同观测、参数化、正则、初值和迭代配置，报告
106: 
107: $$C_{system}=C_{load}+C_{basis}+\sum C_{trial}+\sum C_J+C_{gate}+C_{fallback}+C_{validation}+C_{optimizer}+C_{save}.$$
108: 
109: 活动支撑减少的是物理状态阶数，不改变可用观测；固定总场与耗散界也不增加材料信息。故“相同材料结果、较低 $C_{system}$”是求解器收益；“同预算更低材料风险/更好拓扑”需另有实验。R11 已有 72 份真实反演，但 B/C 全部 forward trial 回退，图像没有新增改善；支撑的正面核结果不能掩盖全阶 LU 预算失败。资格：评价定义及可审计计算分解，不是物理定理。


## research/maxwell_formula_atlas/12_BILATERAL_MATERIAL_MODES_AND_FIXED_POINTS.md


### Lines 5-116

5: 本章接续 F-OP-010、F-MAT 系列及 F-BI-001—024，补齐 R12 实际使用的原 Gaussianized66 拉回、材料局部商、固定弱坐标和完整固定点导数。以下九条是**标准微分与线性代数在冻结离散模型上的具体推导**，不主张新定理。公式成立、代码实现、有限实验支持和方法新颖性是四个独立层次。本次维护新增 Maxwell RHS=0、LU=0；没有为恒等式安排数值重放。
6: 
7: 统一尺寸：ROI 材料格点 $q=256$，接收端 $n_r=32$，源 $n_s=8$，原先验实坐标 $p=66$；完整材料网格 $N=576$。固定 $C\in\mathbb C^{q\times q}$、$O\in\mathbb C^{n_r\times q}$、$E_R\in\mathbb C^{q\times n_s}$、宿主响应 $F_h$ 与噪声 $\sigma=0.00119211798514519$。$d\in\mathbb C^q$ 是相对宿主的复材料，$D=\operatorname{diag}(d)$；$D_x=\partial d/\partial x\in\mathbb C^{q\times p}$ 是 decoder 导数，二者不混用。$^*$ 为共轭转置，$^T$ 为实转置或注明的普通转置。所有优化坐标均为实数。
8: 
9: ## F-BF-001 · 当前总场与接收传播共同决定精确材料导数
10: 
11: 电流和数据满足
12: 
13: $$h=D(E_R+Ch),\qquad F=F_h+Oh.$$
14: 
15: 令 $M=I-DC$、$N_c=I-CD$，假设当前状态可逆；由 $\det(I-AB)=\det(I-BA)$，两种可逆性等价。定义
16: 
17: $$R_D=E_R+Ch=N_c^{-1}E_R\in\mathbb C^{q\times n_s},\qquad L_D=OM^{-1}\in\mathbb C^{n_r\times q}.$$
18: 
19: 对电流方程求微分，注意材料与场都改变：
20: 
21: $$dh=(dD)(E_R+Ch)+DC\,dh,$$
22: $$\boxed{dF=L_D(dD)R_D.}$$
23: 
24: 因此同一格点 $i$ 的材料改变由接收列 $L_{D,:,i}$ 与照明行 $R_{D,i,:}$ 相乘传到全部端口。无需 $C,D$ 对易，也无需 $D$ 可逆。这里 $dD$ 表示材料矩阵的微分；若使用占比 $a_i$ 且 $d_i=\delta a_i$，其列导数必须再乘 $\delta$。已含该因子的占比 Jacobian 不能再次乘它。
25: 
26: **假设与资格：** 冻结有限维状态方程、同一材料下全部源共享 $D$、当前 resolvent 可逆；这是精确一阶导数，不是有限步等响应曲线定理。宿主 Born 仅取 $D=0$ 时的 $L_0=O,R_0=E_R$，一般不等于当前导数。
27: 
28: **计算含义：** 完整当前 $L_D,R_D$ 必须由已付费状态取得；拿到缓存后进行 Gram、SVD 或坐标变换不新增 Maxwell RHS，取得当前场本身仍有真实费用。来源：[R12 双边推导第1节](../p3l_r12_bilateral_fixedpoint/BILATERAL_MODE_DERIVATION.md)、[物理适配器](../p3l_r12_bilateral_fixedpoint/shared/physics.py)；与 F-SS-002/003 的活动支撑导数一致。
29: 
30: ## F-BF-002 · 对角材料的端口对算子、共轭 Hadamard Gram 与截断边界
31: 
32: 按逐源、每源逐接收端排列复数据。材料 Jacobian $J_m\in\mathbb C^{(n_rn_s)\times q}$ 的第 $i$ 列为
33: 
34: $$j_i=R_{D,i,:}^{T}\otimes L_{D,:,i},\qquad J_m=R_D^T\odot L_D,$$
35: 
36: 其中 $\odot$ 指逐列 Kronecker 积。逐项求内积得到
37: 
38: $$\boxed{G_m=J_m^*J_m=(L_D^*L_D)\circ\overline{R_DR_D^*}.}$$
39: 
40: 证明是 $\langle j_i,j_j\rangle=(\sum_r\overline L_{ri}L_{rj})(\sum_s\overline R_{is}R_{js})$；第二因子必须为 $\overline{RR^*}=\overline R R^T$。两个因子均为 Hermitian 半正定，Schur 乘积也半正定。
41: 
42: 进一步取薄 SVD，秩分别为 $r_L,r_R$：
43: 
44: $$L=U_L\Sigma_LV_L^*,\qquad R=U_R\Sigma_RV_R^*,$$
45: $$\boxed{K_{(a,b),i}=\sigma_{L,a}\sigma_{R,b}\overline{(V_L)_{ia}}(U_R)_{ib}.}$$
46: 
47: 因为
48: 
49: $$L\operatorname{diag}(v)R=U_L[\Sigma_LV_L^*\operatorname{diag}(v)U_R\Sigma_R]V_R^*,$$
50: 
51: 左右薄酉因子保持 Frobenius 范数，故 $\|Kv\|_2=\|J_mv\|_2$ 对任意复材料微分 $v$ 成立，完整 $K^*K=G_m$。这是对**独立接收与照明物理端口**的分解；最多 $r_Lr_R\le256$ 个复端口对，并非两个材料弱子空间的集合交集。
52: 
53: **假设与资格：** 使用同一当前 $L,R$ 和全部端口对时为精确表示恒等式。保留最大 $\sigma_{L,a}\sigma_{R,b}$ 的 top32 对是明确的截断启发式；这些乘积并未考虑对角材料限制及 $D_x$，不保证在原 66 维先验上实现最优秩32近似。单边 $L^*L$ 或 $\overline{RR^*}$ 均丢掉另一端信息。
54: 
55: **计算含义：** 完整 pair 与 Full-J 不是两种独立数据来源。top32 复行实化后秩至多64，可在66维材料像上制造截断零方向；source-only 至多16个实行，纯数据 weak16 的选择可能高度非唯一。来源：[双边推导第2—3节](../p3l_r12_bilateral_fixedpoint/BILATERAL_MODE_DERIVATION.md)、[modes.py 的 bilateral_gram / pair_modes](../p3l_r12_bilateral_fixedpoint/agent_a/modes.py)。新增的是本模型的显式对应，不是 Khatri–Rao、Schur 乘积或 SVD 原理。
56: 
57: ## F-BF-003 · 原 Gaussianized66 的复材料导数与实白化拉回
58: 
59: 原先验坐标为 $x=(z,\xi)\in\mathbb R^{64}\times\mathbb R^2$，$u_j=\Phi_N(\xi_j)$；$\Phi_N,\varphi_N$ 分别为标准正态 CDF、PDF。设冻结 quintic saturation 为 $S$，$s=(\mu+A_{ROI}z)/\tau$，则
60: 
61: $$a=S(s),\qquad \delta=0.8+\Phi_N(\xi_1)+i[0.03+0.15\Phi_N(\xi_2)],\qquad d=\delta a.$$
62: 
63: 链式法则给出 $D_x\in\mathbb C^{256\times66}$：
64: 
65: $$D_{x,:,1:64}=\delta\operatorname{diag}(S'(s)/\tau)A_{ROI},$$
66: $$D_{x,:,65}=a\varphi_N(\xi_1),\qquad D_{x,:,66}=i\,0.15a\varphi_N(\xi_2).$$
67: 
68: 令 $\operatorname{pack}$ 按每源32实部、32虚部排列，$W_{real}Z=(\sqrt2/\sigma)\operatorname{pack}Z$；对 Jacobian 逐列应用该实线性映射。对任意实方向 $v$，
69: 
70: $$\|J_wv\|^2={2\over\sigma^2}\|J_mD_xv\|^2,
71: \qquad J_w=W_{real}(J_mD_x)\in\mathbb R^{512\times66},$$
72: $$\boxed{H_{data}=J_w^TJ_w={2\over\sigma^2}\operatorname{Re}(D_x^*G_mD_x).}$$
73: 
74: 证明用 $v$ 为实数及 $\|\operatorname{pack}Z\|^2=\|Z\|_F^2$，再由实二次型唯一确定对称矩阵。不能直接把复矩阵 $D_x^*G_mD_x$ 当实参数 Hessian，也不能漏掉两个 CDF 导数。$H_{data}$ 是 Gauss–Newton 数据曲率；非零残差下它不等于精确目标 Hessian。
75: 
76: **假设与资格：** 冻结原 decoder、固定噪声白化及实参数约定下的精确一阶拉回。平台内 $S'$ 可为零，因此局部材料秩不一定等于66。native 自由占比258模型与此先验支持不同。
77: 
78: **计算含义：** 完整当前 pair、完整当前 Gram 和真实 Full-J 拉回给出同一局部数据几何；不能把不同的归一化或重复 $\delta$ 因子写成方法收益。来源：[双边推导第4节](../p3l_r12_bilateral_fixedpoint/BILATERAL_MODE_DERIVATION.md)、[physics.py 的 Problem.decode / white_jac](../p3l_r12_bilateral_fixedpoint/shared/physics.py)、[modes.py 的 pullback_gram](../p3l_r12_bilateral_fixedpoint/agent_a/modes.py)。
79: 
80: ## F-BF-004 · 材料有效局部商与原先验增广因子
81: 
82: ROI 外材料不变，但材料 RMS 仍按完整576格计算。定义实化材料微分
83: 
84: $$B={1\over24}\begin{bmatrix}\operatorname{Re}D_x\\\operatorname{Im}D_x\end{bmatrix}\in\mathbb R^{512\times66},\qquad G_x=B^TB={1\over576}\operatorname{Re}(D_x^*D_x).$$
85: 
86: 取 $B=U\Sigma V^T$，保留明确阈值以上的奇异值，令 $r$ 为保留数：
87: 
88: $$T=V_+\Sigma_+^{-1}\in\mathbb R^{66\times r},\quad Q=BT=U_+,\quad T^TG_xT=I_r.$$
89: 
90: 这是把单位坐标变化变为单位材料 RMS 微分的局部 lift。精确零奇异值对应局部 decoder 不改变材料的方向；数值阈值以下则只称不可分辨方向。R12 阈值为相对 $10^{-10}$ 与浮点尺寸阈值的较大者，不能由此证明全局 gauge 或全局商流形。
91: 
92: 对于原 Gaussianized66 先验，局部目标的 GN 曲率为 $H_{data}+I$，故直接构造
93: 
94: $$\boxed{A_{aug}=\begin{bmatrix}J_wT\\T\end{bmatrix},\quad A_{aug}^TA_{aug}=T^T(H_{data}+I)T.}$$
95: 
96: 在该显式实因子的右奇异向量中选最弱 $k$ 列 $V_{weak}$，令 $W=TV_{weak}$，则 $W^TG_xW=I_k$。纯数据版本仅分解 $J_wT$；不同消融以自己声明的数据因子替换上块，统一保留原先验下块 $T$。
97: 
98: **假设与资格：** 当地材料导数及固定原 Gaussian 参考测度下的标准白化/增广最小二乘。加入 $I$ 是保留原先验曲率，没有新增观测。该标准分解不能证明 weak16 的有限材料步可达，也不能证明 $0.02W$ 位于可信线性区；$\|W\|$ 可能很大。
99: 
100: **计算含义：** 直接 SVD 显式因子，避免先形成正规矩阵再经病态 $T$ 拉回时弱谱精度恶化。原欧氏阻尼 $\mu I$ 在材料坐标上成为 $\mu T^TT$，不能把 $J_wT$ 的奇异值直接代入带 $\mu I_r$ 的固定点公式。来源：[双边推导第5、7节](../p3l_r12_bilateral_fixedpoint/BILATERAL_MODE_DERIVATION.md)、[modes.py 的 material_canonical_basis / proposal_basis](../p3l_r12_bilateral_fixedpoint/agent_a/modes.py)。
101: 
102: ## F-BF-005 · 同一个 Jacobian 的嵌套 SVD 不增加筛选信息
103: 
104: 设实白化 Jacobian $J=U\Sigma V^T$，$V_I$ 是它自身一个奇异索引集 $I$ 的右奇异向量，则
105: 
106: $$JV_I=U_I\Sigma_I.$$
107: 
108: 所以再对 $JV_I$ 做 SVD，其右奇异向量只是坐标轴；相同奇异值的块内可任选正交基，映回仍为原子空间。若另取它自身左奇异向量 $U_K$，
109: 
110: $$U_K^TJV_I=\Sigma_{K,I},\qquad P_LJP_R=U_K\Sigma_{K,I}V_I^T.$$
111: 
112: 这是原对角奇异结构的子块。对严格不重叠的索引集，子空间正交；有交叠时交集也由原索引确定。重复分解、强弱空间求交或等值块重命名不会产生新的物理信息。
113: 
114: **假设与资格：** 两次使用同一状态、同一度量、同一 $J$ 的奇异向量。若第二次加入独立材料度量、原先验、不同状态的 Jacobian 或独立接收/照明算子，问题已经改变；这种变化须明确写出。本条不否定 F-BF-002 的双物理端口分解。
115: 
116: **计算含义：** 不为该纯冗余构造分配 Maxwell 预算，也不把它列为发现新材料方向的独立方法。来源：[R12 嵌套 SVD 解析说明](../p3l_r12_bilateral_fixedpoint/NESTED_SVD_REDUNDANCY_AUDIT.md)。这是既有 SVD 的直接推论，本轮新增的是对特定提案的排除条件。


### Lines 199-239

199: ## F-BF-009 · 有限 secant、投影变化和覆盖只提供经验诊断
200: 
201: 对两个已实际求值的点 $x_a,x_b$，令 $\Delta x=x_b-x_a$、$\Delta r=r(x_b)-r(x_a)$。冻结起点 Jacobian 的 secant 缺陷是
202: 
203: $$e_{a\to b}=\Delta r-J_a\Delta x,\qquad \delta_{a\to b}={\|e_{a\to b}\|\over\max(\|\Delta r\|,1)}.$$
204: 
205: 若整条线段位于光滑合法域，微积分基本定理给
206: 
207: $$e_{a\to b}=\int_0^1[J(x_a+t\Delta x)-J_a]\Delta x\,dt.$$
208: 
209: 只有额外获得全线段 Jacobian 的 Lipschitz 常数 $L_J$，才能推出 $\|e\|\le\tfrac12L_J\|\Delta x\|^2$。有限端点比较本身不提供这样的常数。终点刷新 $J_b$ 的缺陷 $\Delta r-J_b\Delta x$ 可衡量这个 secant 的经验改善，不能证明全邻域的二阶控制或刷新必有材料收益。
210: 
211: 把不同点的材料弱空间放入同一完整材料实空间，令 $Q_a,Q_b$ 的列正交；则
212: 
213: $$\rho_{a,b}={\|(I-Q_aQ_a^T)Q_b\|_F\over\sqrt{k}}$$
214: 
215: 衡量新 $k$ 维空间在旧空间之外的比例。它是子空间变化，不是新增观测量；data-only weak8 的诊断也不是 prior-aware weak16 的提案轨迹。
216: 
217: 有限合格云 $\mathcal C_{p,h}$ 的冻结经验覆盖为
218: 
219: $$\widehat{coverage}_h(\epsilon)={1\over n_{DEV}}\sum_p\mathbf1\!\left\{\min_{m\in\mathcal C_{p,h}}{\|m-m_p^{GT}\|\over24}\le\epsilon\right\}.$$
220: 
221: 空云按失败计；最近 GT 成员是离线容量诊断，不能用来选择在线输出。此处的 $\epsilon$、噪声门、去重尺度和先封存预测后读 GT 的顺序须由既有协议给定，不能看到结果后移动。
222: 
223: **假设与资格：** secant 和投影差来自有限已存节点；只有在连续线段已满足光滑假设时，积分恒等式才适用。经验覆盖只针对声明的暴露 DEV 对象，不能升级为总体高概率保证、连续根族或后验覆盖。均匀或 $\exp(-\Phi)$ 有限云行动仍须注明有限云语义。
224: 
225: **计算含义：** 现有 R12 诊断保存128条 secant、96个材料投影对，且采样包含被拒 trial 与分支切换，不能称为连续接受路径。H2/HA 刷新次数更多但是否改善行动，必须看已有真实有限材料结果与总成本；局部方向转动本身不能代替收益。来源：[在线诊断已有收据](../p3l_r12_bilateral_fixedpoint/agent_c/p1_diagnostics/ONLINE_DIAGNOSTIC_RECEIPT.json)、[FINITE_SECANTS.json](../p3l_r12_bilateral_fixedpoint/agent_c/p1_diagnostics/FINITE_SECANTS.json)、[冻结评分汇总](../p3l_r12_bilateral_fixedpoint/agent_c/p1_score/SUMMARY.json)。本次没有重算这些诊断。
226: 
227: ### 既有 R12 结果对本章的约束
228: 
229: 截至本章维护时，P1 的4个暴露 DEV 对象、8个臂在冻结 $\epsilon=0.075$ 下均为覆盖0/4。普通 H0 RML 的平均离线最近材料 MSE 为0.0119253，H1/H2/HA 完整16维分别为0.0130177、0.0133596、0.0131985；该基线的云容量指标更好。所有32组的在线选择仍为原锚点，平均所选 MSE 同为0.0148154。这里支持“现有开发结果没有显示分层/自适应刷新带来在线行动改善”，不支持“所有双边方法都不可能有效”。数字来自已有[评分汇总](../p3l_r12_bilateral_fixedpoint/agent_c/p1_score/SUMMARY.json)与[候选门收据](../p3l_r12_bilateral_fixedpoint/agent_c/P1_ADVERSARIAL_QA.json)，本次没有重新评分。
230: 
231: 原生 N1 的既有45组运行停留在初始化边界，未到达非锚点非线性机制比较；刷新是同一锚点的代数 frame 重算。其状态保留为 `INITIALIZATION_BOUNDARY_FAILURE_NOT_PHYSICAL_TEST`，不能作为新原生材料方法的正面或负面性能结论。后续 N2 的实现、封存和真值顺序由主线程及对应工作目录独立管理，本章不提前赋予结果。来源：[N1 已有边界收据](../p3l_r12_bilateral_fixedpoint/agent_c/N1_BOUNDARY_AUDIT.json)。
232: 
233: ### 已知结构、待比较先例与仍未建立的贡献
234: 
235: F-BF-001—009 使用 resolvent 微分、Khatri–Rao/Hadamard 结构、SVD、最小二乘增广、链式法则、隐式约束和有限经验评分等既有工具。本章贡献是把这些工具在同一个原先验/白化/材料 RMS 合同中写清，并排除不成立的谱与新颖性推论；它不声称创造了这些数学工具。
236: 
237: 已有 ScholarQA 记录核实的最近主题先例为 Zhong、Chen 的 [Twofold subspace-based optimization method for solving inverse scattering problems (2009)](https://doi.org/10.1088/0266-5611/25/8/085003)、[FFT Twofold Subspace-Based Optimization (2011)](https://doi.org/10.1109/TAP.2010.2103027)，以及 Xu 等的 [Multiplicative-Regularized FFT Twofold Subspace-Based Optimization (2015)](https://doi.org/10.1109/TGRS.2014.2329032)。这里只复用已完成的 **Tier C 出版身份** 核实；没有可用逐段全文证据证明它们与本轮具体算法相同或不同。429 和未覆盖的 distorted-Born、张量分解、广义奇异值与先验感知流形优化文献缺口均保留，不能宣称“首次”。来源：[既有最近工作身份说明](../p3l_r12_bilateral_fixedpoint/agent_a/CLOSEST_PRIOR_WORK_METADATA.md)、[已有五条身份核实记录](../p3l_r12_bilateral_fixedpoint/agent_a/SCHOLARQA_VERIFIED.json)。
238: 
239: 状态：上述公式为 `DERIVED_FINITE_MODEL`；实现仅在所引源码的声明范围内；P1 改善结论为 `NOT_DEMONSTRATED_ON_EXPOSED_DEV`；原 P1 连续根流形、完整材料后验与总体覆盖为 `NOT_ESTABLISHED`；方法新颖性为 `UNRESOLVED_WITH_METADATA_ONLY_PRECEDENTS`。最终科学判断由主线程作出。
