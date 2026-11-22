# 原生 96² 网格的分块光学证书

**#03 半环、10% 噪声通过原生网格复核。** 对 **optfresh_03__pw180_10x20__n10**，排除「全域 \(0\le\chi\le0.5\)」的残差下界为 **0.622783219965**。冻结噪声半径加 \(10^{-6}\|Y\|_F\) 后的上界为 **0.594048781081**，严格裕量至少 **0.028734438885**。该结论针对生成观测所用的离散 FFT Green 算子和端口，涉及全部 9216 个材料格点。

## 分块乘子与频域下界

保留 32² 扫描选出的 \(4\times4\) 光学乘子，每块在 96² 网格上对应一个 \(24\times24\) 方区。\(\mu\) 和每块的 \(\lambda\) 均除以 9，对应每个粗格覆盖的细格数。材料上限 \(u\) 保持不变。#03 的两档 \(\lambda/\mu\) 约为 0.22 和 0.58。将其对角矩阵写成 \(tI+dM\)，其中 \(t\) 为两档中点、\(d>0\) 为半差，\(M\) 为对角符号矩阵，逐格取 \(-1,+1\)。实对角权重与实材料乘法交换，因此按 \(M\) 加权的虚部矩阵等式仍成立。

设 \(T\) 为从 96² 零填充到 192² 后进行单位化二维 DFT 的等距映射。生成算子为 \(G=T^*\operatorname{diag}(g)T\)，\(S\) 为接收算子，\(B\) 为入射场矩阵。记 \(\operatorname{Herm}(G)=(G+G^*)/2\)、\(\operatorname{Im}_H(G)=(G-G^*)/(2i)\)。矩阵 \(A=K/\mu\) 为

\[
A=I+S^*S/\mu-u\operatorname{Herm}(G)+t\operatorname{Im}_H(G)
  +\tfrac{id}{2}(G^*M-MG).
\]

空间项中的 \(G\) 可替换为 \(G-cI\)，因为实常数 \(c\) 的贡献相消。令 \(z=Tx,\ w=TMx\)，则 \(\|z\|_2=\|w\|_2=\|x\|_2\)。记 \(a_q=1-u\operatorname{Re}(g_q)+t\operatorname{Im}(g_q)\)。以下 \(2\times2\) 块作用于 \((z_q,w_q)\)，其二次型求和等于 \(x^*(A-\alpha I-S^*S/\mu)x\)。

\[
C_q=
\begin{pmatrix}
a_q-\alpha-\beta d^2 & \tfrac{id}{2}\overline{(g_q-c)}\\
-\tfrac{id}{2}(g_q-c)&\beta d^2
\end{pmatrix}.
\]

两项 \(\beta d^2\|z\|^2\) 与 \(\beta d^2\|w\|^2\) 相消。取 \(\beta=0.875\)、\(c=0.234375\)、\(\alpha=0.6875\)，对方形频域集合 \(Q=\{-10,\dots,10\}^2\) 外的每一项，用 112 位 Arb 验证其 Schur 补

\[
a_q-\alpha-\beta d^2-\frac{|g_q-c|^2}{4\beta}>0.
\]

最小严格裕量为 \(8.7750\times10^{-4}\)。保留 \(Q\) 内的 441 个频点，令 \(U=T^*|_Q\) 和 \(F=[U,MU,S^*]\)。把保留块装入 \(C\)，并加上接收对角块 \(I/\mu\)，得到

\[
A\succeq A_0=\alpha I+FCF^*.
\]

为了证明 \(A\succeq\gamma I\)，取 \(\gamma=1/256\)。由同一分块矩阵的两种 Schur 消元，
\(n_-(A_0-\gamma I)=n_+(C^{-1}+F^*F/(\alpha-\gamma))-n_+(C)\)。
用有限几何级数和区间 DFT 封住 \(F^*F\)，再用位移 Cholesky 残差和负子空间二次型严格复核 902 维矩阵的惯量。其正、负方向数为 649、253，正方向数等于 \(n_+(C)\)。该矩阵与 \(C\) 都已证非奇异，因此 \(A_0-\gamma I\) 为正定。

## 完整算子的驻点残差

频域下界用于证明谱底 \(\gamma\)。数据界使用保留所有频点的原生算子。记 \(L'=L/\mu=S^*Y/\mu+[uI+i(tI+dM)]B/2\)。取候选解 \(X\)，令 \(R=AX-L'\)。有

\[
H=Y^*Y-\mu L'^*A^{-1}L'
 =Y^*Y+\mu(X^*AX-L'^*X-X^*L')-\mu R^*A^{-1}R.
\]

沿候选 \(H\) 的前两个照明方向，用未截断的区间 FFT 重算两列驻点残差。由 \(A\succeq\gamma I\)，修正量上界为 \(\mu\|R\|_F^2/\gamma\)，在这两个方向所张的试探子空间上得到 \(3.473\times10^{-17}\)。对一维、二维试探子空间，分别用投影分子的 Gershgorin 下界及基 Gram 上界作 Rayleigh 复核。所证前两正特征值下界为 0.31884891284595 和 0.06901002622428。取两者之和开平方即得报告的残差下界。

## 配对结果与适用范围

| 原始观测与检验上限 | 冻结半径 | 原生分块下界 | 原生同预算标量下界 |
|---|---:|---:|---:|
| #03 半环，10%，\(u=0.5\) | 0.594042806767 | **≥0.622783219965（区间），\(\nu=0\)** | 0.544721500760（数值），\(\nu=0\) |
| #04 全环，2%，\(u=0.75\) | 0.095726084788 | 0.056036870114（数值），\(\nu=4\) | 0（数值），\(\nu=4\) |

标量对照使用粗网格同搜索预算中选定的冻结参数，按相同比例迁移，走相同的原生数值流程。#03 的分块和标量分别使用 1622、1636 次 Green FFT 卷积。#03 证实冻结分块剖面可跨过该半径，冻结标量在数值对照中低于半径。任意全局标量的最优性未由这组比较判定。#04 的分块数值界距半径尚差 0.039689214674，因此该剖面未进入区间复核。这个负结果适用于冻结剖面的迁移。

[冻结粗网格候选及证书](../results/core/gradient_certificates.json) · [两例原生数值回放](../results/theory/fine_spatial_followup.json) · [#03 严格复核记录](../results/theory/fine_spatial_verified.json) · [复算源码](spatial_followup.py)

在项目根目录使用具有 NumPy、SciPy、threadpoolctl 和 python-flint 的 Python 环境依次运行下列命令。

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B research/z_alg_optical_inertia/fine/spatial_followup.py
python -B research/z_alg_optical_inertia/fine/spatial_followup.py --verify
~~~
