# 有限介质板的 FFT-TSOM 与全波材料反演

在 24 mm × 24 mm 介质板上恢复 576 个独立的相对介电偏差。输入是带缺陷物体和同一宿主空板的复散射场之差。两种方法都从均匀板出发，可以逐格改变材料。宽频和角度都能改善恢复，作用随采样布局和方法变化。前向网格太粗时，降低数据误差会产生明显的材料伪影。

## 同一个有限板模型

无损板的相对介电常数为 2.25。每格未知实值 q 代表相对已知板的介电偏差，范围为 [−1.25, 1.5]。8 个线源和 16 个接收器分别位于半径 45 mm 和 50 mm 的圆上。原始试验有 3 个开发对象、8 个独立非空对象和 1 个空板。另有 9 个新随机对象用于配对盲测。

在每个频率上，用体积分方程消去整个已知有限板，得到宿主入射场、域内传播算子和接收映射。扰动电流和差分读数满足

```math
D_f=\operatorname{diag}(a_f q),\qquad
J_f=D_f(E_{h,f}+G_{h,f}J_f),\qquad Y_f=S_{h,f}J_f.
```

无损组 a_f 为 1。Debye 组 a_f = (ε_h(f)−1)/1.25，宿主和扰动均随频率变化。所有频率仍共享同一实值材料图。原观测在 48×48 网格生成，评分材料是其面积平均的 24×24 图。每个等预算设计用总能量 1、8 次独立发射和 64 个复数读数。2% 噪声相对 16 GHz 开发圆孔的差分场 RMS 0.0022717125 定义，然后按每发能量重标。 [采集模型](src/ndt.py) [采样与归一化](src/designs.py) [配对数据](HOLDOUT_DATA.md)

## FFT-TSOM 的电流空间

对每个频率选中的接收矩阵做 SVD。确定电流由读数直接反演。余下电流先在低空间频率的二维傅里叶基中展开，再投影到确定电流子空间的正交补。

```math
S_h=U\Sigma V^H,\qquad
J_\mathrm{det}=V_L\Sigma_L^{-1}U_L^H Y,\qquad
J=J_\mathrm{det}+(I-V_LV_L^H)\operatorname{IFFT}_2(\beta\odot M).
```

作者的 `mymusk.m` 保留未移位二维 FFT 数组的四角，每角为 M_F×M_F。每轴频率为 [−M_F, M_F−1]。作者的 `TwoD_Inversion_MOM_FTSOM.m` 实现了上述接收投影、逐阶段扩展掩膜，以及交替更新电流与材料。陈书式 (6.71)–(6.72) 和共同作者的 MR-FFT-TSOM 论文式 (5)–(11) 同时支持该构造。2011 年原论文全文未取得。 [原论文 DOI](https://doi.org/10.1109/TAP.2010.2103027) [作者教材附录](https://www.ece.nus.edu.sg/stfpage/elechenx/Book.htm) [作者 MATLAB 实现](https://www.ece.nus.edu.sg/stfpage/elechenx/My_book/Chapter10_Fig4_6.zip) [MR-FFT-TSOM 作者稿](https://renchengsong.github.io/papers/TGRS2329032.pdf)

本轮 M_F 依次为 3、6、9，各阶段最多 10、14、60 个材料步。事先规定保留奇异值至少为首奇异值 13%、读数投影的发射平均 RMS 超过噪声标准差 2.75 倍的最后一个方向之前的全部方向。给定材料时，以 Cholesky 解各频电流的线性最小二乘子问题，再用包络梯度和 L-BFGS-B 更新全部盒内材料。有限宿主的 G_h 没有循环平移结构，程序预计算投影 IFFT 字典经过 G_h、S_h 的响应。

目标包含加权数据残差、加权本构残差和材料正则。数据权重为 w = 1/(s₀√N)，其中 N 为复数读数总数、s₀ 为 16 GHz 开发圆孔的基准 RMS。本构权重为 0.32 w σ₁(S_h)。在每个输出材料上重新解完整波动方程，并独立记全波数据残差。以 30 GHz 单频为例，8 个接收方向全部由确定电流拟合，辅助电流的数据残差中位数只有约 7.0×10⁻¹⁴ 倍噪声标准差，输出材料的全波数据残差中位数为 **4.75 倍**。 [实现](src/inversion.py)

## 全波投影 Gauss–Newton

全波 GN 每步重新求解 A_f J_f = D_f E_h,f，其中 A_f = I − D_f G_h,f。对第 j 格材料的解析导数为

```math
\partial_{q_j}Y_f
=a_f\,(S_{h,f}A_f^{-1})_{:,j}\,
(E_{h,f}+G_{h,f}J_f)_{j,:}.
```

复残差和雅可比分别叠放实部和虚部，再解盒约束 GN 子问题，用真实全波目标做 Armijo 回溯，最多 24 步。两种方法的同权对照使用平滑 TV 与近均匀板的平滑 L1，平均权重分别为 0.0055 和 0.0035，平滑幅度为 0.04。强 GN 另用原三开发对象，在共同倍率 1、2、4、8、16、32、64 中按平均材料 RMSE 选权重。多数高频臂选中 64 倍，落在已测范围的上界。盲测原单频、均分四频、28/30 GHz 定位设计分别冻结 16、8、64 倍。

对非零材料做中心差分，常数宿主的全波雅可比相对误差为 7.97×10⁻¹¹，Debye 宿主为 2.25×10⁻¹⁰，Debye FFT-TSOM 的变投影方向梯度为 1.36×10⁻¹⁰。时间列包含已知宿主冷建和单对象反演，独立记录校准和评分用时。 [反演入口](run_inversion.py) [开发校准](results/inversion_metrics_gn_calibrated.json) [盲测缺少的两臂校准](results/inversion_metrics_gn_calibration_holdout.json)

## 原仿真结果

表中为原八个独立非空对象的中位数，每格依次为材料 RMSE、符号 IoU、冷建与反演秒数。以 |q| ≥ 0.20 定义支持，IoU 的重叠还要求正负号相同。

| 设计 | 复读数 | FFT-TSOM | GN 同权 | GN 开发校准 |
|---|---:|---:|---:|---:|
| 16 GHz · 8 发 × 8 收 | 64 | .087 / .525 / .349 | .081 / .536 / .288 | 见盲测 |
| 2 / 30 GHz · 4 发 × 8 收 | 64 | .139 / .238 / .737 | .065 / .775 / .419 | .063 / .755 / .448 |
| 4 频 · 2 发 × 8 收 | 64 | .106 / .430 / 1.386 | .072 / .674 / 1.001 | 见盲测 |
| 8 频 · 1 发 × 8 收 | 64 | .092 / .568 / 2.000 | .066 / .707 / 1.124 | — |
| 15.68 / 16.32 GHz · 4 发 × 8 收 | 64 | .092 / .537 / .479 | .085 / .596 / .342 | — |
| 30 GHz · 8 发 × 8 收 | 64 | .108 / .374 / .240 | .088 / .432 / .185 | .063 / .708 / .262 |
| 4 频 · 8 发 × 16 收 | 512 | .060 / .667 / .859 | .131 / .217 / .955 | .036 / .856 / .926 |
| 4 频 · 2 发 × 2 收 | 16 | .133 / .162 / 1.025 | .117 / .467 / .572 | — |
| 4 频 · 2 发 × 16 收 | 128 | .105 / .433 / 1.202 | .105 / .255 / .935 | .061 / .811 / .936 |
| 2 / 30 GHz · 单侧 | 64 | .138 / .191 / .477 | .052 / .768 / .297 | .051 / .770 / .593 |
| 28 / 30 GHz · 定位设计 | 64 | .093 / .536 / .797 | .078 / .528 / .560 | .056 / .864 / .483 |

四频全阵列的 GN 校准后，数据相对残差中位数从 .126 升到 .141，材料误差从 .131 降到 .036，空板标记从 3 格降到 0 格。其余同权 GN 设计和全部 FFT-TSOM 设计的空板标记均为 0。两峰判据要求最强的两个负对比局部极大分别在两个真实圆心 1.5 mm 内，近双孔还要求小于半间距。30 GHz、四频全阵列和 28/30 GHz 均使两种方法分开原组 3、4.4、6.2 mm 的三个双孔。 [主矩阵逐例](results/inversion_metrics.json) [28/30 GHz 逐例](results/inversion_metrics_task_design.json) [增强正则前后图 PDF](figures/reconstruction_dense_regularization.pdf)

30 GHz 无噪圆孔将观测从 24 网格换成 48 网格生成，GN 材料 RMSE 从 .010 升到 .063，数据相对残差仍约为 .001。近双孔和细槽的材料 RMSE 分别从 .036 升到 .086、从 .043 升到 .082。FFT-TSOM 在无噪同网格三例仍有 .076、.090、.086 的材料 RMSE，表明即使消除了噪声和网格差，受限电流表示仍会留下误差。 [同网格无噪](results/inversion_metrics_30ghz_same_clean.json) [跨网格无噪](results/inversion_metrics_30ghz_cross_clean.json) [同网格带噪](results/inversion_metrics_30ghz_same_noise.json)

## 新随机对象配对盲测

三臂的总能量、发射次数、复读数数相同，两种物理模型共享相同连续材料和噪声实现。读取材料评分前，先逐模型写出 135 条预测。表中为八个非空对象的中位数，每格依次为材料 RMSE、符号 IoU、对比度保留。保留率等于真支持上估计与真实材料的点积除以真实材料平方和，值 1 表示投影幅度匹配真值。

| 输入 / 使用的宿主 / 方法 | 16 GHz | 2、12、20、30 GHz | 28 / 30 GHz |
|---|---:|---:|---:|
| 常数 / 常数 / FFT-TSOM | .072 / .515 / .554 | .091 / .406 / .546 | .085 / .529 / .656 |
| 常数 / 常数 / 校准 GN | .059 / .667 / .656 | .052 / .686 / .709 | .051 / .820 / .749 |
| Debye / Debye / FFT-TSOM | .066 / .549 / .606 | .074 / .582 / .556 | .074 / .626 / .618 |
| Debye / Debye / 校准 GN | .061 / .667 / .643 | .054 / .682 / .664 | .047 / .826 / .728 |
| Debye / 常数 / 校准 GN | .080 / .617 / .359 | .083 / .609 / .338 | .090 / .483 / .305 |

28/30 GHz 在看读数前选定。本批只有此设计在两种宿主、两种算法下都分出 2.885 mm 双孔，三臂都分出 5.858 mm 双孔。误用常数宿主时，28/30 GHz 仍分出近双峰，对比度保留只有 .305。正确 Debye 三臂的校准 GN 冷建加反演中位数约 0.24、0.72、0.36 秒。每个模型与算法组合的八个非空对象均有标记，唯一配对空板均无标记，这批样本不足以估计实物误报率。 [常数盲测](results/inversion_metrics_holdout_constant.json) [正确 Debye](results/inversion_metrics_holdout_debye.json) [错误常数](results/inversion_metrics_holdout_debye_wrong_host.json) [Debye 材料图](figures/reconstruction_projected_gn_holdout_debye.png)

## 复现

```sh
python run_inversion.py
python run_inversion.py --design-json results/task_design.json --arms loc_nominal --suffix _task_design
python run_inversion.py --design-json results/task_design.json --arms f4_t8_r16 f4_t2_r16 high_t8_r8 loc_nominal f2_t4_r8 oneside_f2_t4_r8 --methods projected_gn --calibrate-gn --suffix _gn_calibrated
python run_inversion.py --arms f1_t8_r8 f4_t2_r8 --methods projected_gn --calibrate-gn --development-only --suffix _gn_calibration_holdout
```

盲测先给 `--data` 传各组 `*_observations.npz`，并给 `--design-json results/task_design.json --arms f1_t8_r8 f4_t2_r8 loc_nominal --gn-penalty-json results/inversion_metrics_gn_calibrated.json results/inversion_metrics_gn_calibration_holdout.json`。正确 Debye 加 `--host-model debye`，错误常数只运行 `--methods projected_gn`。每组指定独立的 `--suffix`。全部预测完成后，分别用 `--score-blind --data` 传带真值的对应文件，并给 `--geometry results/holdout.json` 和该组后缀。 [复现代码](run_inversion.py)
