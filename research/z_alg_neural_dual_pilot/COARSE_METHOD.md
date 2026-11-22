# 32² 谱梯度的推导与核查

保留 32² 全部 1,024 个状态。接收和入射分别用冻结的 \(S\in\mathbb C^{20\times1024}\)、\(B\in\mathbb C^{1024\times10}\)。\(D_b\) 标记第 \(b\) 个空间块。两种块网格分别为 4×4 和 8×8，每个块包含 64 和 16 个粗网格像素，按先 \(y\) 后 \(x\) 排列。

令 \(Z=(uD_\mu+iD_\lambda)/2\)，状态矩阵 \(K=S^*S+D_\mu-ZG-G^*Z^*\)，右端项 \(L=S^*Y+ZB\)。固定一组乘子，求 \(J_0=K^{-1}L\)、\(E_0=B+GJ_0\) 和 \(H=Y^*Y-L^*K^{-1}L\)。逐次完整 Hermitian \(LDL^*\) 分解给出数值惯性和驻点解。接近数值秩边界时沿用旧验证器的密集谱分解。\(H\) 的特征值从大到小排列，丢弃前 \(\nu=i_-(K)\) 项，将其余正特征值相加，所得 \(d^2\) 保持观测场平方单位。

由对 \(H\) 求微分得到
\[
\mathrm dH=J_0^*(\mathrm dK)J_0-2\operatorname{Re}_H\bigl(J_0^*\mathrm dL\bigr).
\]
因此在惯性与正方向不变、截断边界有谱间隔的区域，
\[
\begin{aligned}
\partial_{\mu_b}d^2&=\operatorname{tr}P\left[J_0^*D_bJ_0-u\operatorname{Re}_H(J_0^*D_b E_0)\right],\\
\partial_{\lambda_b}d^2&=\operatorname{tr}P\operatorname{Im}_H(J_0^*D_b E_0).
\end{aligned}
\]
这里 \(P\) 投影到保留的正特征向量。零特征值、丢弃和保留两群之间的重根、状态惯性跳变处，这个经典导数不保证存在。

训练时用温度 \(\tau=0.01\|Y\|_F^2/10\) 的 softplus 近似保留谱上的正部，导数权重为 \(\operatorname{sigmoid}(h_i/\tau)\)。前 \(\nu\) 个模式仍按实际惯性截去。如果 \(\nu\ge10\)，硬下界为零，代理值改为 \((\|Y\|_F^2/10)\tanh(\kappa_{10}/10^{-3})\)，其中 \(\kappa_{10}\) 是 \(K\) 的第十小特征值。对简单特征值用其单位特征向量 \(v\) 求梯度，分别为 \(\|D_bv\|^2-u\,\mathrm{Re}(v^*D_bGv)\) 和 \(\mathrm{Im}(v^*D_bGv)\)，再乘代理链式因子。相同 \(K\) 的代理特征对在同一进程中缓存，收据分别记录计算和缓存命中。这两种代理分数只用于生成候选。

纯梯度和网络共用非负 \(\mu\) 与带符号 \(\lambda\) 参数化。参数化零点给 \(\mu=0.001\)、\(\lambda=0.0005\)。纯梯度从此出发。精修从网络的实际输出出发，最初评价和最多三次候选检验都计费。用参数化梯度最大分量归一化，初始步长 0.5。接受后乘 1.35（上限 2），拒绝后减半。在有正值时用硬目标的梯度。硬值为零时可按代理上升，在所有候选硬值均为零时也用代理择优。正硬值相同则保留较早候选。纯梯度由同一最长 5 秒运行保存两个快照，分别从在 1 秒、5 秒内已算完的候选里择优。未赶上预算的计算单列并仍计入实际开销。

[核查回执](results/coarse_gradient_check.json)使用固定种子的复高斯场，没有读取留出观测或材料。4 个工况包括正定、保留了正谱的不定状态，以及 10 个源方向被完全耗尽。与独立 1,024 阶特征分解比较硬分数，并对 \(\mu\)、\(\lambda\)、联合方向及参数化方向逐一做中心差分。16 个方向的最大相对误差为 \(3.67\times10^{-6}\)。本机这批核查中，普通一次评估的分组中位墙钟为 45–56 ms，未缓存且需要特征值代理的评估约 181 ms。数字和微分均是粗网格上的数值证据，做严格排除时要以原生 96² 的向外区间为准。

在项目根目录安装 NumPy、SciPy 和 threadpoolctl 后，可用下列命令重复梯度与预算检查。固定对照和 SRAP 补充生成另外使用 python-flint。固定对照的主盲测与补充输出目录分别为 `results/` 和 `results/supplement/`。

```bash
PYTHONDONTWRITEBYTECODE=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  python research/z_alg_neural_dual_pilot/src/check_coarse_gradients.py
PYTHONDONTWRITEBYTECODE=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  python research/z_alg_neural_dual_pilot/src/fixed_baselines.py
PYTHONDONTWRITEBYTECODE=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  python research/z_alg_neural_dual_pilot/src/generate_srap_bands.py
PYTHONDONTWRITEBYTECODE=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  python research/z_alg_neural_dual_pilot/src/fixed_baselines.py --supplement
```

固定对照使用已经封存的常权和两个 16 块乘子。相同 96² 物理算子和乘子允许先把共同的 30 个右端项做 Arb 向外驻点检查。每份观测随后只用这个共同结果计算源谱区间。共同准备和每条观测的重用耗时分别记录。补充生成仅在原模型冻结后执行。按原始随机种子复现 SRAP 的复高斯方向并按存档的 10% 扰动模长缩放，然后沿同一方向从存档读数构造 1% 与 2% 视图。再独立用同一原生几何和材料生成器证明新误差球包含精确前向响应。补充包只保存有噪场、上舍入半径、对象编号和档位。
