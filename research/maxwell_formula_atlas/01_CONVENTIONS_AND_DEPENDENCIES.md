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
