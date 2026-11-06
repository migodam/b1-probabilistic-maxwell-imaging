# GN–Green R2：可搜索材料空间的数学定义

本文件由主线程推导。固定二维 scalar-TM disk-pulse VIE、576个共享实材料参数、源为列；全部反演与有限步验收使用完整非线性物理。下列均为条件线性代数结论，不声称新的全局恢复定理。

## 1. 从暗电流空间拉回共享实材料

令 $L=I-D_mG$ 可逆，$Q=E+Gj(m)$，$W=L^{-1}$。真实材料切向是

$$\mathcal T_mh=W D_hQ,\quad h\in\mathbb R^n,\qquad J_mh=\operatorname{pack}_{\mathbb R}(S\mathcal T_mh).$$

沿用原固定谱阈值，定义 $\mathcal D=V_{\mathrm{in},W}\cap\ker S$，$U$ 为其正交补基。由于内部强、弱空间互补，

$$U=\operatorname{orth}[V_{\mathrm{in},S},V_{\mathrm{ob},\mathrm{nonzero}}].$$

历史配置中 $U$ 有67个复方向，$\mathcal D$ 有509个复方向。这些维数属于单个电流列，尚未施加共享实材料约束。

定义实线性算子 $\mathcal A_mh=\operatorname{pack}_{\mathbb R}(U^*\mathcal T_mh)$，则

$$\mathcal H_m=\ker\mathcal A_m\subseteq\ker J_m.$$

证明：$\mathcal A_mh=0$ 意味着每个源的电流切向均在 $\mathcal D\subseteq\ker S$，故所有接收切向为零。因此若材料 Jacobian 满列秩，则 $\mathcal H_m=\{0\}$。这不是说不存在有用的弱响应方向；它说明必须将数值核、噪声容差下的近似方向和有限材料步区分开。

实验以实际SVD的相对阈值 $10^{-12}$ 报告数值秩，不用 $\mathcal A_m^T\mathcal A_m$ 的小特征值冒充高精度零空间，也不称其为符号认证。多源子集固定为1/3/6/12个嵌套源；新增行不能扩大精确核，但相对阈值或数据归一化后的有效维数不一定单调。

若 $q=\dim\mathcal D^\perp$、源数为 $p$，实化的泄漏算子至多有 $2pq$ 行，因此秩—零度定理给出 $\dim\ker\mathcal A_m\ge\max(0,n-2pq)$。本例 $q=67$，1源和3源的下界分别是442、174；6源已允许数值满列秩，但此计数不保证满秩。实测443/189/0/0属于当前采集的结果，不能作为任意Maxwell几何的源数定理。

## 2. 高效构造实材料的电流能量与泄漏度量

记 $C=QQ^*$。对实材料 $h$，有

$$\|\mathcal T_mh\|_F^2=h^TH_Th,\quad H_T=\Re[(W^*W)\odot C^T],$$

$$\|U^*\mathcal T_mh\|_F^2=h^TH_Ch,\quad H_C=\Re[((U^*W)^*(U^*W))\odot C^T].$$

这里是普通转置 $C^T$，不能漏掉复数共轭。逐项看，$H_{ij}$ 包含 $\langle W_{:,i},W_{:,j}\rangle\sum_s\overline{Q_{is}}Q_{js}$。两个Gram均在实材料域对称半正定，且 $0\preceq H_C\preceq H_T$。

当每个材料行的总场向量非零时，$H_T$ 正定：$W D_hQ=0$ 由 $W$ 可逆推出 $D_hQ=0$，继而 $h=0$。因此可求

$$H_Cv_i=\mu_i H_Tv_i,\qquad 0\le\mu_i\le1.$$

取 $\mu_i\le\zeta$ 的张成空间，则其中每个非零方向均满足电流泄漏能量比不超过 $\zeta$。主门为1%，敏感性门为10%。这保证的是相对切向电流分布，不保证绝对数据漂移小。后者还依赖材料步长和 $J_m$。

广义特征向量在 $H_T$ 度量下正交，计算Euclidean材料误差覆盖时必须对所选列另做Euclidean正交化；实现使用列Gram的Cholesky白化，和QR得到相同投影。按特征值事前排序，不能用GT重新排序。

显式 $W$ 需要一个现有LU上的576个复RHS；不能把它说成免费的固定Green投影。比对时计入这项准备费用。无需存储所有源的巨大材料—电流三阶张量。

## 3. 与材料尺度绑定的有效自由度

令 $K=J_m/\|Y\|_F$ 的奇异值为 $s_i$。沿单位Euclidean方向移动材料RMS为 $r$，该轴的线性相对数据变化为 $r\sqrt n\,s_i$。

$$d_{r,\varepsilon}=\#\{i:r\sqrt n\,s_i\le\varepsilon\},\qquad
d_{\mathrm{soft}}=\sum_i\frac{1}{1+(r\sqrt n\,s_i/\varepsilon)^2}.$$

$d_{r,\varepsilon}$ 是指定局部线性模型、半径与容差的轴数；$d_{soft}$ 是软权重之和，不是整数流形维数。沿每个轴的线性允许RMS半径为 $\varepsilon/(\sqrt n s_i)$；零轴仍受材料box及有限非线性限制。clean中的容差是假定可接受的数据变化，不代表观测含有这种噪声。

GT误差覆盖定义为 $\|P_{\mathcal U_r}(m_*-\widehat m)\|^2/\|m_*-\widehat m\|^2$。这是离线表示上限，不是在线方法已经知道正确系数。原始复Green投影、共享实材料空间与真正Full-J方向分别报告。

## 4. 同一结构提案的两种物理过滤

Green方法将同一实Haar提案投到上述1%泄漏材料空间；若空间为空，记录空空间，不以测试GT选择放松阈值。

照明感知方法解

$$\min_h\tfrac12\|Kh\|^2+\tfrac\lambda2\|h-p\|^2,\quad
h=\lambda(K^TK+\lambda I)^{-1}p,$$

$$\lambda=(\varepsilon/(r\sqrt n))^2.$$

它是标准软谱滤波，权重为 $\lambda/(s_i^2+\lambda)$。过滤后不再把提案放大回原RMS，否则会取消抑制可观测方向的效果。满足box时使用单一射线缩放，保持方向；全部最终材料再由原始完整Maxwell验证。

过滤器只决定允许怎样移动。预测正确移动符号和形状仍来自结构先验或以后训练的NN；暗空间本身不提供缺失的GT标签。

## 5. 有限非线性响应与局部横向收缩

若 $L(m+h)$ 亦可逆，完整响应满足

$$F(m+h)-F(m)=S L(m+h)^{-1}D_hQ(m).$$

所以 $J_mh$ 小不保证有限 $h$ 的数据变化小。实验记录实际漂移、线性漂移及两者的余项，不能仅凭切向分类接受候选。

在内点精确数据解 $m_*$、无正则条件下，阻尼GN映射

$$T_y(m)=m-(J_m^TJ_m+\lambda I)^{-1}J_m^T[f(m)-y]$$

具有

$$DT_y(m_*)=\lambda(J_*^TJ_*+\lambda I)^{-1}.$$

强可观测方向收缩较快；近零方向保留更多位移；精确零方向导数为1。该结论是局部的，且box活动面变化、非零残差、正则梯度以及实际TRF全局化都会改变行为。本实验不把这个导数当作SciPy TRF的精确迭代导数。

实际corrector测量同一批材料提案在固定anchor强/弱投影中的校正前后位移和观测响应。弱位移被保留不代表它朝GT；GT距离仅在评分阶段评价。

## 6. 归因边界

- 既有GN没有使用28维TSOM截断；误差投影不能证明TSOM造成其材料偏差。
- 只读历史GT作离线诊断，在线提案与排序只取公开观测和合法材料。
- 普通继续TRF、正则减弱和数据项对照用于识别优化未完成与正则偏差。
- 不同噪声、方法端点属于同一材料的相关重复，不扩大独立样本数。
- 本轮不训练NN、不增加真实测量信息，也不提出概率质量或全局等响应纤维覆盖结论。

## 7. 材料box让可走方向成为锥，而不是对称子空间

在下界活动点必须有 $h_i\ge0$，在上界活动点必须有 $h_i\le0$。因此真正的一阶可行集合是 $\mathcal U\cap T_{\rm box}(m)$。即使 $\mathcal U$ 维数很大，任意正负射线也可能有一端或两端被不同活动像素阻断。

本轮预定的射线缩放是 $a=\min(1,\min_{h_i>0}(u_i-m_i)/h_i,\min_{h_i<0}(l_i-m_i)/h_i)$。若某个下界像素对应 $h_i<0$，则 $a=0$。这是该条射线不可行，不能推出整个子空间与box的交集只有原点。部分有效GT投影射线仍可以移动；报告必须将它们与非GT的Haar提案分别计数。

下一步可考虑把box纳入软滤波二次问题联合求解，而不是事后逐像素裁剪或把全方向缩到零。本轮没有执行这个新算法，不把建议计为结果。
