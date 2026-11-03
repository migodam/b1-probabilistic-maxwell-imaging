# R6：A-material 与 FFT-TSOM 的组合边界

## 结论和身份

本轮没有提出新的 Maxwell 方程，也没有把 Fourier 表示改名为新算法。A-material 是 R4 的局部本构缺陷传播、经 R5 作用位置分离后的固定算法：普通 global-PR CIE 更新辅助电流；仅材料读出使用一次真实 Green 邻域传播。FFT-TSOM 则是这个辅助电流空间的已有表示与实现框架，两者作用位置不同，因而可组合。是否更好恢复材料必须由本轮封存后的同数据比较回答。

以下均是有限 disk-pulse scalar-TM 矩阵结论；不外推连续三维 Maxwell。约定 exp(-iωt)，背景归一化 εb=μb=1，m=εr−1 为实正 contrast；G 包括原自由空间 Green 与原有限脉冲自项。源接收均在未知域外。冻结物理为

```text
Q = Einc + G J,
r(m,J) = J − diag(m) Q,
Jphys(m) = [I − diag(m)G]⁻¹ diag(m)Einc,
F(m) = S Jphys(m).
```

最终 F 使用完整 LU 求解，包含该离散模型的全部多次散射；一次传播只作用于辅助缺陷，不是一次散射近似正演。

## 1. A 改的是材料块的离流形几何

令 Goff 为预定不重叠2×2 cell块的非对角 Green，相应自项为 Gii。固定辅助 J、Q 后，定义

```text
b = β − diag(G),       d = 1 + b m,
a = m/d,               u = β (J − mQ)/d,
L = diag(a) Goff,      RA = (I+L)u,
fA(m | J,Q) = ρ/2 tr(RA Ωs RA*) + Rprior(m).
```

乘除均按cell进行，Ωs 是固定源metric。R4已有整个合法box上的无极点/一次传播可逆充分界，本轮复用其适用网格和材料范围，不重复验证既有恒等式。相同零集合不代表相同离流形目标、材料更新或搜索盆地；加上prior后更不能从零集合相同推出同一最优材料。

精确实材料导数为

```text
du_i/dm_i = −β (Q_i + b_i J_i)/d_i²,
da_i/dm_i = 1/d_i²,
dRA = (I+L)du + diag(dm/d²) Goff u.
```

因此不能把一次传播当成冻结的常数metric再忽略第二项。实现调用 R4 已有精确 value/gradient；R5 的两项缓存只省去完全相同 m 的重复计算，不改变梯度。材料读出用同一box、同一公开prior、最多60内迭代；每次标量CIE解只是已收费的seed，不以GT选择。

R5已经证明并诊断：残差零点附近，权重只能通过材料切向正交补影响材料读出；离流形还有上述 dL 项。R6检验其实际适用域与代价，不把这个已知加权最小二乘关系再次记作原创。

## 2. 为什么能放到 FFT-TSOM 上

取 S 的保留右奇异空间 Vd，Qd=I−VdVd*，DPCS current 为 JD。原始Fourier坐标系数 A 下，

```text
J = JD + U A,              U = Qd Fmask,
U A = Qd IFFTortho(mask·A),
U* z = mask* FFTortho(Qd z).
```

固定 m 的普通CIE current残差为

```text
p = βcurrent/(1+βcurrent m),
H A = p·[UA − m·G(UA)],
H* z = U*{conj(p)·z − G*[m·conj(p)·z]}.
```

材料读出只需要 J 和 Q，不需要显式保存 U、GU、H 或 SU。因此替换 U/G 的action与adjoint后，A-material公式及其prior完全不变。CPU实现用每轴至少2n−1的零填充嵌入实现有限Toeplitz G，保留原自项；不是把24×24 ROI当周期物理域。接收 S、低秩DPCS投影、最终物理LU和corrector仍为dense。完整技术合同见 [FFT_TSOM_INTEGRATION.md](FFT_TSOM_INTEGRATION.md)。

这是已有FFT-TSOM框架上的可运行组合，属于算法实现与成本验证，不是FFT或新的电磁表示的发明。历史显式U还用SVD清理数值零作用；native FFT保持同一精确算术action，但浮点非线性长轨迹不承诺逐位相同。终点改变不能仅凭更好的某张图归因为FFT提供了新信息。

## 3. 表示容量与成本不能混算

mf3/6/10对应36/144/400个native Fourier系数；mf12对应576，但扣除DPCS后APCS秩在当前rank17接收配置中为559。每个源的current容量是559复自由度，材料始终576独立实自由度。576个native系数有17个近零作用，不是多了17个可观测方向。

即使Fourier mask全满，Vd*J仍被固定为Vd*JD；有噪声时该DPCS约束不等于真实current在相同子空间的投影。因此“full mode”不能称辅助current完全无约束，也不能用它证明算法已经达到无条件能力上限。

固定邻域尺寸时，A每个材料value/gradient的额外工作为O(Ns)，但每次读出可能调用多次。current更新的 FFT 主要为O(sN log N)，另有O(sNr_d)的DPCS及O(sNNr)的接收作用；dense物理corrector的分解仍为O(N³)、多RHS回代O(N²s)。不能把只加速辅助块说成整个反演都变成FFT复杂度。

## 4. 公平比较及可辨识性的边界

四臂都用历史原物理目标 Phi=0.5||F(m)−Y||²/||Y||²+Rprior(m) 作在线选择。这里W=I/||Y||是原基线的归一化，不冒称按σ写出的完整Gaussian负对数似然；公开σ另用于噪声discrepancy评分。改用σ且不同比例缩放prior会改变问题，本轮没有这样改基线。

CIE和Dual-CIE保持成熟sourcewise PR；Global-Dual是分离global PR影响的必要控制。A与Global-Dual差在材料残差；A与Dual-CIE还差current PR的聚合方式，不能把这两个差别混成一次传播的净效果。所有方法同mode日程、同corrector上限，实际LU/RHS、内层读出和wall另算；名义迭代数相等不意味着算量相等。

无论A或FFT是否获得更低Phi，都不能推出材料更准确、全局唯一或更大可辨识空间。本轮仅有一个固定初始化、一个新的离散合成噪声等级，不能据此证明跨初始化、损耗、孔径、连续网格、三维或实测泛化。

## 5. 与NIE的关系及新颖性判断

NIE2016与CIE2019已有非线性重写、FFT-TSOM与交替读出的明确先例。R6不把self吸收等同原NIE整个算法，也不把缺少本地成熟源码当成NIE不存在。详见 [NIE_BASELINE_AUDIT.md](NIE_BASELINE_AUDIT.md)。本轮主比较采用已验证成熟CIE/Dual-CIE，最直接回答A作为该框架插件是否值得使用。更强的独立NIE复现是另一个实现/验证任务，不能用一个弱替代实现补足表格。

状态：公式与组合结论 **DERIVED（有限模型、给定可逆性条件）**；FFT数值接线 **EMPIRICAL**；材料优势与成本竞争力由结果报告单列，不能从推导直接宣布。
