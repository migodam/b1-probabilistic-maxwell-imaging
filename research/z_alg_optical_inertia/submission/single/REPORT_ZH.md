# 单源加入空间约束后，四宫格仍有差距，十六宫格已经追上

在开槽块 01 和错列沟槽 05 的原生 96² 全环 2% 观测上，整体上限均为 0.75。对任意单位复相干照明，四宫格单源 SDP 的残差最优值严格为零。改成十六宫格后，两件样品都存在一个相干照明及一组空间乘子，能严格排除该材料上限。相同十六宫格乘子在联合矩阵中各只保留一个正照明方向，单源证书已经达到这两组固定乘子的联合界。

## 同一份观测上的结果

所有数均为复接收残差的范数。两份数据均使用原生 9,216 格、20 接收口、10 照明的同一个有限算子。冻结的误差半径包含前向生成补量。对单位照明组合 \(a\)，同一噪声球给出 \(\|Ea\|_2\le\|E\|_F\le\delta\)。相干方向可从全部 10 照明测量及入射场选择，因此这里比较的是在相同 200 个复读数下的认证方法。原始某一束照明的单次采集无法自行完成这种数据依赖的方向选择。

| 证据 | 01 开槽块 | 05 错列沟槽 |
|---|---:|---:|
| 同一个扩大噪声半径 | 0.090967634889 | 0.088025053207 |
| 四宫格单源对任意相干组合的精确 SDP 最优值 | 0 | 0 |
| 联合四宫格或十六宫格可用的常权重曲线，\(\nu=4\) | ≥0.160698781023 | ≥0.110607872802 |
| 十六宫格单源相干证书，完整 32 个分区乘子 | ≥0.189674078449691 | ≥0.124521510770564 |
| 按同一组十六宫格乘子做联合，正谱合成 | [0.189674078449691, 0.189674078449705] | [0.124521510770564, 0.124521510770574] |

联合常权重证书在两种分区上都可取每块相同的 \(\mu=0.00035,\lambda/\mu=0.3\)，故分区约束完全可用。01 的常权重实部对照另外给出 ≥0.168329834663。十六宫格的两份数据分别用一组权重。对各自固定的权重，原始十束照明的二次对偶值上界分别为 −0.7450 和 −0.2621，取非负下界后均为零。严格越线来自相干合成。两份数据的最优十六宫格 SDP 值及最优联合权重仍待求解。 [投稿图 PDF](partition_comparison.pdf) · [网页图 PNG](partition_comparison.png) · [原联合证书](../../results/theory/fine_verified.json) · [同一误差球](../../results/theory/fine_feasibility.json) · [四宫格](local_moment.json) · [十六宫格](coherent_general_verified.json) · [匹配联合回放](matched_joint_replay.json)

## 四宫格同时施加八条约束

每块 \(b\) 的扇区量为 \(q_{r,b}=\|P_bj\|^2-u\operatorname{Re}j^*P_b(Ba+Gj)\)，光学量为 \(q_{o,b}=\operatorname{Im}j^*P_b(Ba+Gj)\)。单源 SDP 用均值 \(x\) 和协方差 \(C\succeq0\) 表示一次与二次矩，同时要求四条 \(q_{r,b}\le0\) 和四条 \(q_{o,b}=0\)。目标为 \(\|Sx-Ya\|^2+\operatorname{tr}(SCS^*)\)。原二次问题对应零协方差。若存在真实材料，其电流在这个松弛中必然可行。

先由公开读数给出均值矩阵 \(X\)，再用准确的接收 Gram 解强制 \(SX=Y\)。四个暗电流 \(v_\ell\) 也经过同样的准确投影，满足 \(Sv_\ell=0\)。以 \(v_\ell\) 在各块的齐次光学量组成 \(4\times4\) 实矩阵 \(\tau\)，再解

\[
\sum_{\ell=1}^4\tau_{b\ell} W_\ell=-Q_{o,b}(X).
\]

Arb 验证了 \(\tau\) 的可逆性、每个 \(W_\ell\succ0\) 和四块加权实部的严格负定性。给定任何照明 \(a\)，取 \(x=Xa\) 与 \(C(a)=\sum_\ell(a^*W_\ell a)v_\ell v_\ell^*\)，即可同时满足全部局部条件并使目标严格等于零。提升矩阵的秩至多为 5。四宫格的最大实部扇区谱上界为 \(-5.27\times10^8\) 与 \(-8.71\times10^7\)。协方差系数的规模很大，该可行性属于给定的局部矩 SDP，其物理秩一问题仍然开放。 [区间复核](local_moment.json) · [冻结提出量](local_moment.npz)

## 十六宫格候选和证明

固定空间权重 \(w_b\ge0,t_b\in\mathbb R\) 后，设 \(Z=\operatorname{diag}((u w+i t)/2)\)、\(A=S^*S/\mu+\operatorname{diag}w-ZG-G^*Z^*\) 与 \(L=S^*Y/\mu+ZB\)。若 \(A\succ0\)，单源给定方向的残差平方下界是 \(a^*Ha\)，其中 \(H=Y^*Y-\mu L^*A^{-1}L\)。对固定权重取最大特征向量就可选择最佳相干组合。联合界是 \(H\) 的正特征值和的平方根。分区优化采用每块 36 个 Legendre 张量模态作为搜索子空间。独立的原生数值回放选出相干方向，区间算术严格证明原生算子的正性、相干残差和相同乘子下的全源谱。

严格验证保留两组傅里叶场。一组是电流场，另一组是经分区复权调整的电流场。取 \(Z=z_0 I+cM\) 后，将 \(\gamma I+\beta MM^*\) 加到傅里叶二次型，并从准确的分区对角项扣除。在每个遗漏模态上检查

\[
\gamma-2\operatorname{Re}(z_0g)-c^2|g|^2/\beta>0.
\]

保留 \([-13,13]^2\) 两组模态和 20 个接收列，形成 1,478 阶区间惯性系统。两次均严格得 \(A\succeq I/1024\)。随后在原生算子上计算驻点缺陷，并以 \(\mu\|Ax-La\|^2/\operatorname{gap}(A)\) 控制消元损失。对联合回放，计算 \(H\) 的一个 Loewner 上界，分别证明九维子空间严格负定。此前已知 \(H\) 有正特征值，因此两案均恰有一个正源方向，联合与相干单源的固定乘子界相等。 [求乘子](coherent_general_search.json) · [冻结乘子](coherent_general_profiles.json) · [原生数值回放](coherent_general_native.json) · [原生区间地板与残差](coherent_general_verified.json) · [匹配联合正谱](matched_joint_replay.json)

## 实际花费与证书前沿

时间均为同一台 Mac 的实际墙钟。共享阶段按两件合计报告。原标量联合验证的两线程与十六宫格的四线程 Arb 设置在各自回执中列出。

| 方法及阶段 | 实际墙钟 | 物理与计算计数 |
|---|---:|---|
| 四宫格单源零最优值构造，共享预处理 | 1.69 s | 80 个浮点 FFT 卷积列、4 个区间列、1 次 20 阶接收 Gram 求解 |
| 四宫格单源逐件验证，复用预处理 | 2.54 s / 2.60 s | 每件 10 个浮点 FFT 列、10 个区间列、1 次 20 阶接收 Gram 求解 |
| 原生常权重联合证书，冻结 \(t=0.3\) 的共享证明 | 6.61 s | 645 阶区间惯性核，复用几何；两件共享同一因子，逐案后处理时间未测 |
| 十六宫格单源对称 84 组原生诊断 | 66.29 s | 137,632 个浮点 FFT 卷积列；这批扫描在 05 上未越线 |
| 十六宫格单源一般权重候选 | 16.21 s | 576 个浮点 FFT 列建 576 阶搜索子空间，466 次子空间最小特征值及 Cholesky 评估 |
| 十六宫格单源四组权重原生回放 | 6.28 s | 13,202 个浮点 FFT 卷积列，包含用全部公开读数选方向 |
| 十六宫格单源两件区间验证，含共享准备 | 69.44 s | 34.75 s / 34.42 s 逐件，分别验证 1,478 阶区间惯性核；对应 154 / 152 个浮点 FFT 列和每件 2 个区间 FFT 列 |
| 十六宫格固定权重的全源谱回放，复用两组正性证据 | 6.09 s | 每件 20 个区间 FFT 列，01/05 各 1,534/1,526 个浮点 FFT 列 |

常权重联合在两件上以一份共享的 645 阶证明越线。一般十六宫格单源提高了已证实的残差界，付出分区搜索和两次较大的区间证明。直接复算十六宫格冻结证明花 69.44 s。按可复算候选路径把 16.21 s 搜索和 6.28 s 回放也计入为 91.93 s，再计探索性的 84 组对称前缀共 158.22 s。原 64.44 s 联合标量扫描覆盖 2,832 条记录及多个噪声、阈值和几何，不能按这两件拆分。四宫格的共享项加逐件项给出 4.23 s / 4.29 s 冷启动，合计实测 6.95 s。全程真实材料正演重建求解次数为零。涉及认证二次型的 Krylov 迭代已列在卷积成本中。各权重族均只比较已支付的证书，早期试算的费用没有按样品精确摊销。 [四宫格冷启动](local_moment.json) · [对称前缀](coherent_dual.json) · [一般权重回放](coherent_general_native.json) · [原联合扫描](../../results/theory/fine_scan.json)

## 复制冻结结果

在本目录运行下列命令，环境安装 NumPy、SciPy、python-flint 和 threadpoolctl。命令读取同一目录内已有的冻结乘子、原生相干方向和四宫格提出量。发布包还需带上该研究的 synthetic/public/fine_geometry.npz、两份 optfresh_{01,05}__pw360_10x20__n02.npz、results/theory/fine_feasibility.json、core/certificate.py 和 core/inertia.py。证明只使用公开测量、公开算子和保存的误差球。

    python -B local_moment.py --proposal local_moment.npz --output local_moment_replay.json
    python -B general_fourier_verify.py
    python -B matched_joint_replay.py

已安装 Matplotlib 时，运行下列命令即可从封存回执重画英文投稿图。

    python -B plot_partition_comparison.py

需要重跑选方向和搜索时，可运行下面的命令。依赖本机数值优化器的步骤输出到另外的文件，便于同冻结提议核对。

    python -B coherent_dual.py --profiles coherent_general_profiles.json --output native_replay.json
    python -B coherent_profile_search.py --output search_replay.json --profiles candidate_replay.json
    python -B coherent_dual.py --output symmetric_replay.json
