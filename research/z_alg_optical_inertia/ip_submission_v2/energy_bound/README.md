# 原生网格上的相干电流上界

对保存的 96×96 二维标量 TM 算子、全环 10 发 20 收、材料上限 \(u=0.75\)，取 \(b=Ba\)、\(y=Ya\)、\(\|a\|_2=1\)。下表的上界覆盖所有满足本构关系 \(j=\operatorname{diag}(m)(b+Gj)\)、\(0\le m_i\le u\) 且 \(\|Sj-y\|_2\le r\) 的相干电流，包括由整个测量矩阵的误差球投影得到的每个相干电流。\(r\) 沿用原 2% 观测的扩大半径。

上界和真值分数均按 \(\|j\|_2^2=\sum_{i=1}^{9216}|j_i|^2\) 报告。\(j_i\) 是归档算子的单元对比度电流 \(m_i e_i\)，沿用归档入射场振幅。相干向量以十束入射列的系数欧氏范数归一化。数据项是 20 个接收端口复场的平方和。若把电流平方和换为对 \(2\,\mathrm m\times2\,\mathrm m\) 区域的按面积 \(L^2\) 平方，需统一乘 \((2/96)^2\,\mathrm m^2\)。SDP 的对应量是 \(\operatorname{tr}(C)+\|x\|_2^2\)。

| 原始 2% 观测 | 扩大半径 \(r\) | 数据球内严格上界 | 复用原正性证明时的上界 | 真材料离线相干峰值区间 |
|---|---:|---:|---:|---:|
| `optfresh_01__pw360_10x20__n02` | 0.09096763488916171 | 4,982,632.805676025 | 53,240,809.03608654 | [66,192.78950620674, 66,192.78950848749] |
| `optfresh_05__pw360_10x20__n02` | 0.0880250532073907 | 1,035,000.28174651 | 11,880,925.739471074 | [31,647.962400386616, 31,647.96240156425] |

01 开槽块和 05 错列沟槽来自两件独立材料。较早的噪声曲线也逐元素确认两个 `pair2pct` 中心与各自 `n02` 完全相同。 [数据冻结](../../synthetic/FREEZE.md) · [公开输入登记](../../synthetic/public/manifest.json) · [原生算子](../../synthetic/public/fine_geometry.npz) · [半径回执](../../results/theory/fine_feasibility.json) · [配对观测审计](../../submission/curve/observation_audit.json)

## 证明

固定十六区乘子 \(\mu w_i\ge0,\mu t_i\in\mathbb R\)，其中 \(\mu>0\)。记
\[
 Z=\operatorname{diag}((uw_i+\mathrm i t_i)/2),\quad
 A=K/\mu=S^*S/\mu+\operatorname{diag}(w_i)-ZG-G^*Z^*,\quad
 L=S^*Y/\mu+ZB .
\]
把扇区不等式和光学等式加到数据残差平方，得到
\[
 \|Sj-Ya\|_2^2 \ \ge\
 a^*Y^*Ya-2\mu\operatorname{Re}(j^*La)+\mu j^*Aj .
\]
若 \(\delta I\preceq A\)，选 \(0<\sigma<\delta\)，令 \(A_\sigma=A-\sigma I\)、\(H_\sigma=Y^*Y-\mu L^*A_\sigma^{-1}L\)。配方后
\[
 \|Sj-Ya\|_2^2 \ge \mu\sigma\|j\|_2^2+a^*H_\sigma a .
\]
因此 \(H_\sigma\succeq-cI\) 和噪声球一起给出 \(\|j\|_2^2\le (r^2+c)/(\mu\sigma)\)。这里的 \(A\) 含 \(S^*S/\mu\)，证书依赖残差半径。证明使用十六区的物理约束和正性下界，涵盖每个相干方向。

计算中以十个驻点试列组成 \(X\)，令 \(D=A_\sigma X-L\)。恒等式
\[
 H_\sigma=Y^*Y+\mu(X^*A_\sigma X-X^*L-L^*X)-\mu D^*A_\sigma^{-1}D
\]
允许用 \(\mu\|D\|_F^2/(\delta-\sigma)\) 界住最后一项。Arb 在保存的二进制数组上重算 FFT、缺陷和十阶矩阵。浮点特征向量用来提出基，再用区间 Gershgorin 界验证基和合同变换。扩大正性界时复用原十六区傅里叶下包络并重做紧凑惯性证明，遗漏模式也逐个用区间验证。 [原十六区验证器](../../submission/single/general_fourier_verify.py)

真值段单独打开 `synthetic/evaluation/optfresh_{01,05}.npz`。用真材料的活跃网格重解十束原生 VIE，再用验证过的逆算子上界和区间原生 FFT 残差围住 \(\max_{\|a\|=1}\|J_{\rm true}a\|^2\)。真材料峰值分别为 3.631929583738103 和 3.4083102321253347，均超过 0.75。真值电流给出参考尺度，竞争者电流上界由公开观测和算子确定。

## 材料箱内能达到的电流平方

另取三组与观测无关的合法材料和固定单位相干照明，原生 VIE 的电流平方严格落在下列区间。常值材料取 \(m_i=0.75\)。两份非均匀材料以归档 \(B,G\) 建议，并在验证前固定。优化目标是当前照明的电流平方，材料变量保持在 \([0,0.75]\)，随后用该材料对十束照明的响应提出相干方向。 [显式材料与照明](receipts/coarse_12x12.npz) · [三案区间验证](receipts/physical_current.json)

| 材料 | 固定相干方向的电流平方区间 |
|---|---:|
| 常值 0.75 | [7,288.608223053841, 7,288.610430242617] |
| 原生像素优化构造 | [8,705.792643632329, 8,705.795202097323] |
| 12×12 块状优化构造 | [8,886.624046063722, 8,886.626033220236] |

这给任何适用于整个材料箱、所有单位相干照明的统一标量电流上限设下严格下限 8,886.624046063722。它高于四宫格零残差见证的上界 8,499.987374387498 和 6,999.999930290155。因此这类统一上限的取值都保留两个见证。此推论的范围是与观测无关、对全部材料和照明同时有效的标量上限。 [四宫格见证回执](../energy_sdp/verified_alt025.json)

合法例的存在性和误差控制用同一个去掉接收 Gram 的正性证明。令 \(A_0 \succeq \delta_0 I>0\) 为经辐射界转移后的状态型，\(Z_0\) 为相应的乘子。齐次 VIE 的任何解满足 \(j^*A_0j\le0\)，故齐次解为零，有限维 VIE 可逆。设 \(x=Mf\) 为试电流，\(h=x-j\)。在活跃材料格上有 \((M^{-1}-G)h=f-Gx-Ba\)，再用 \(2\operatorname{Re}(h^*Z_0(M^{-1}-G)h)\ge\delta_0\|h\|^2\) 和 Cauchy–Schwarz 围住电流误差。区间程序直接由固定试总场 \(f\) 形成 \(x\)，只在活跃格计算缺陷，无需除以接近零的 \(m_i\)。 [去接收项的证明回执](receipts/radiation_transfer.json) · [合法材料验证器](source/verify_physical_current.py)

## 复算与费用

从项目根运行，环境需 NumPy、SciPy、python-flint、threadpoolctl。记录用 Apple M4、Python 3.13.13、NumPy 2.5.5、SciPy 1.18.1 和 python-flint 0.9.0。

```sh
python -B research/z_alg_optical_inertia/ip_submission_v2/energy_bound/source/verify_energy_bound.py \
  --inherited-gap --output research/z_alg_optical_inertia/ip_submission_v2/energy_bound/receipts/energy_bound_inherited_gap.json
python -B research/z_alg_optical_inertia/ip_submission_v2/energy_bound/source/verify_energy_bound.py
python -B research/z_alg_optical_inertia/ip_submission_v2/energy_bound/source/score_truth_energy.py
python -B research/z_alg_optical_inertia/ip_submission_v2/energy_bound/source/verify_physical_current.py \
  --proposal research/z_alg_optical_inertia/ip_submission_v2/energy_bound/receipts/coarse_12x12.npz \
  --proposal research/z_alg_optical_inertia/ip_submission_v2/energy_bound/receipts/pixel_profile.npz
```

复用原十六区的 \(A\succeq I/1024\) 并取 \(\sigma=1/2048\)，新增认证花 7.03 秒墙钟，含 01/05 的 1,534/1,536 个浮点 FFT 卷积列和各 20 个区间 FFT 卷积列，真实材料正演次数为零。 [原十六区回执](../../submission/single/coherent_general_verified.json) · [复用正性界的电流回执](receipts/energy_bound_inherited_gap.json)

较强的数据球内界复算换用 33×33 傅里叶保留方块，两案均严格验证 \(\delta=12/1024\) 并取 \(\sigma=9/1024\)。新增认证含两次 2,198 阶紧凑惯性证明，合计墙钟 283.57 秒、进程 CPU 509.99 秒，01/05 的紧凑惯性证明各需 138.73/138.29 秒，驻点界另需 3.25/3.16 秒及 1,580/1,632 个浮点 FFT 卷积列、各 20 个区间 FFT 卷积列。该数值未计各试探网格的探索，原 27×27 傅里叶保留方块试图验证 \(1/64\) 的正性界未通过，随后改用当前配置。 [紧正性界的电流回执](receipts/energy_bound.json)

离线真值复算花 3.12 秒墙钟，每件一次活跃状态 LU、十列真材料正演及十列 Arb FFT 检查，另含各自活跃状态矩阵的近似逆验证。 [真值区间回执](receipts/truth_energy_offline.json)

三件显式材料的在线无关验证合计花 1.19 秒墙钟，含常值材料的十列数值选向和三件各一列的原生数值求解及区间 FFT。单独重放最终 12×12 提案的搜索花 4.95 秒、59 次目标与梯度评估及 10,466 个浮点 FFT 卷积列，在第 40 次 L-BFGS-B 迭代达到预算。搜索前还试过其他初始条件及像素优化。复现证明直接读取冻结材料和照明。 [合法材料验证回执](receipts/physical_current.json) · [12×12 搜索代码](source/search_physical_current.py) · [12×12 搜索重放](receipts/physical_search_replay.json)
