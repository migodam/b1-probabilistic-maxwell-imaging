# 原生端口压缩 Jensen 补证

试验在固定的 96² 全环标量 TM 模型上使用上限 0.75。两个端点是 Paper 1 已证明正定的 16 区权重 0 和 2。整个试验只在同一条权重线段上比较。验证用 12 件原验证材料的 2% 视图，独立批次用 32 件材料的三个配对噪声视图，历史压力测试另用 28 件旧材料的 84 个视图。历史批次包含提供端点权重的材料。

## 如何获得 30×30 区间矩阵

令 $T_j=[S^*,Z_jB]$、$A_j=K_j/\alpha$、$\ell_j=T_j/\alpha$。两端点的冻结物理尺度均为 $\alpha=0.00035$，归一化正谱下界均为 $g=1/1024$。保留权重记录里的二进制数，再用 Arb 精确乘以尺度。中点在 Arb 内对归一化权重取算术半和，因而由正定性的凸性继承同一个归一化下界。

原代码给出 30 个驻点右端的数值提案 $X_j$，再通过原生 FFT 复核 $R_j^{\mathrm{def}}=A_jX_j-\ell_j$ 及二次型 $M_j=X_j^*A_jX_j-X_j^*\ell_j-\ell_j^*X_j$。实际响应满足

$$
Q_j = T_j^*K_j^{-1}T_j
= -\alpha M_j+\alpha (R_j^{\mathrm{def}})^*A_j^{-1}R_j^{\mathrm{def}},
\qquad
0\preceq Q_j+\alpha M_j\preceq
\alpha (R_j^{\mathrm{def}})^*R_j^{\mathrm{def}}/g .
$$

用 PSD 余项的对角上界及其 Cauchy–Schwarz 关系包住每个复条目。`results/port_enclosures.json` 保存实部与虚部的 Arb 中点、半径，每个数以精确二进制尾数和指数表示。两端点共支付 60 个原生右端，独立中点对照额外支付 30 个。没有求解任何新材料的散射正问题。在线评估只访问这三份压缩矩阵和观测。

## 端口下界与整段上界

对新观测取 $\Xi=[Y;I]$、$H_j=Y^*Y-\Xi^*Q_j\Xi$、$R_j=Q_j[1:20,1:20]$、$E=(Q_1-Q_0)[1:20,:]\Xi$。以 $t$ 为右端点权重，令 $A_t=tR_0+(1-t)R_1$。对任意数值提出的端口解 $X$，用严格区间评价

$$
H(t)\succeq
(1-t)H_0+tH_1+t(1-t)\bigl(E^*X+X^*E-X^*A_tX\bigr).
$$

投影项限制在 $S^*$ 生成的端口方向，它的优化值由 $A_tX=E$ 给出。下面的材料下界可以直接用有限精度 $X$ 的变分值，无需求逆的公差假设。源矩阵正谱之和是对所有收缩投影的迹的最大值。程序认证数值子空间 Gram 的谱上界，缩放成收缩投影，再用区间计算迹的下端。中点能量诊断另以 $A_t$ 的正谱间隔包住端口驻点缺陷，给出端口吸收比例的双侧区间。

对整段的上界，$H(t)$ 的矩阵凹性给出两段仿射上包络。从左到右，三个顶点分别是 $2H_{1/2}-H_1$、$H_{1/2}$、$2H_{1/2}-H_0$。正谱迹对 Loewner 序单调且凸，因此取三个顶点的严格正谱迹上端的最大值，就覆盖所有实数 $t\in[0,1]$。这项界沿用精确中点矩阵，无需新增原生右端。本批数据的外推上界较松，没有据此判定整段未过半径的观测。

在线路径在固定的七个内部二进制网格点作数值提案，严格核验最强内部提案及两个端点。中点对照、整段上端和全网格区间诊断另计时。完整回执中，`projected_midpoint.distance_upper` 上界的是该端口试探产生的下界矩阵，`native_midpoint.distance_upper` 上界的才是原生中点的完整证书。每个已核验的网格下端均与整段上端比较。

## 复现入口

从工作区根执行。使用此前区间复核的 NumPy、SciPy、python-flint 和 threadpoolctl 环境。省略 `--rebuild` 时重用压缩矩阵。

~~~bash
PYTHONDONTWRITEBYTECODE=1 .artifacts/z_alg_diag_python/bin/python -B \
  research/z_alg_neural_dual_pilot/followup/jensen_port.py --validation
PYTHONDONTWRITEBYTECODE=1 .artifacts/z_alg_diag_python/bin/python -B \
  research/z_alg_neural_dual_pilot/followup/jensen_port.py \
  --observations research/z_alg_neural_dual_pilot/data/followup_observations.npz
PYTHONDONTWRITEBYTECODE=1 .artifacts/z_alg_diag_python/bin/python -B \
  research/z_alg_neural_dual_pilot/followup/jensen_port.py --historical
~~~

三个批次分别保存在 `results/validation_jensen.json`、`results/fresh_jensen.json` 和 `results/historical_jensen.json`。对历史批次的两个端点，168 次压缩区间均与冻结的原生端点区间相交，过线判断全部一致。核对记录见 `results/archived_endpoint_crosscheck.json`。
