# Z-ALG R14　给下一轮 GPT Pro 的研究交接

**最值得复核的论文主张，是在同一已知物理模型中同时回答两个问题。现有端口和能量无论怎样选频，都可能漏掉多大的材料差？一条实际新回波到达后，又能保证一个具体材料估计准到什么程度？** R14 已经有整段频率的反例，也有独立噪声对照和全材料输出误差界。旧数据直接主动选频没有胜过固定宽带。

[可读网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r14/) · [完整研究包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/research-z-alg-r14-20261119/Z_ALG_R14_2026-11-19.zip) · [固定版本源码与结果](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/research-z-alg-r14-20261119/research/research_z_alg_r14)

## 已知的计算结果

1. **一条带噪旧读数确实可以对应两片断开的材料域。** 在 $[.5,30]^6$ 上，R10 低维余子式排除比 R13 的全枢轴扰动界更有效。62 份历史旧数据中，有 16 份找到两份严格相容材料，同时在包含所有可行材料的外包上证明一道完整的空投影带。其余记录的连通性未定。R12 的 36 份是六套条件几何的重复，R10 自然、有限孔径各有 18、8 个目标。[旧域比较](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/research-z-alg-r14-20261119/research/research_z_alg_r14/PIVOT_RESULTS.md) · [断开的逐案见证](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r14-20261119/research/research_z_alg_r14/results/disconnected_old.json)
2. **测量前，没有选频优势。** 即使在材料盒内先计算两频分式之差、保留分子抵消和共同增益，再比较所有区域对，62/62 仍回退固定 ±2%。18/18 个自然目标对固定双频都有候选材料相距超过 1 RMS、又可能共用未来读数的严格见证。部分低噪案例的上下界仍相距很远。真正选频失败和证明力不足在不同记录中占比不同。[事前与实际截面对照 CSV](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r14-20261119/research/research_z_alg_r14/results/case_summary.csv)
3. **有限带宽的能量限制可以对一切反馈策略成立。** 在自然目标 `fresh6-16` 的固定旧读数下，已知增益 1 时有两种材料相距至少 7.9733032 RMS。对连续相对带宽 ±2%、任意三端口相干馈电、有限次采集、逐路径最多一轮扫频能量，两个假设可能得到相同的完整测量史。任何估计器的有界误差最坏 RMS 至少 3.98665。若未来噪声独立圆复高斯、两假设等先验，给定这条旧读数之后的最优二元判错率至少 1.014%。允许两个固定合法增益时另有一对相差 11.568 RMS。552 段计算包括 Hankel 及功率 Gram 的数值误差。已知增益反例对任意方向的小接收器平移至少保持到 $2.00\times10^{-6}\lambda_0$ 左右。空间扰动界很小，尚无网格收敛或真实硬件误差界。[频带绝对上界与风险](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r14-20261119/research/research_z_alg_r14/results/continuous_summary.json) · [位置扰动原始数值](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r14-20261119/research/research_z_alg_r14/results/receiver_shift_check.json)
4. **收到回波后，全域缩盒能帮助决定停机，也能修复一个数值错解。** 每频率记 $F_j=N_j/q_j$，用 $N_j-y_jq_j=-\zeta N_j-q_je_j$，各块共享复增益 $1+\zeta$。叠加实测读数后，每个残差对各独立材料仍至多一次。首批半能量后 R12 33/36 能证整个集合直径 ≤0.5，两批后 36/36；平均用 0.5417 轮新增能量。包括旧测量，总能量节省 22.9%。自然与有限孔径均 0/26。重新生成八套几何材料，首批低噪 7/8、高噪 0/8。新低噪 `val6-06` 的八起点误差 15.71，同一回波用区域初始化、限制同一 TRF 后为 0.000904，可靠误差半径 0.0728 以内。旧 R10 盒直接加入实测数据，在 36 份条件数据上可为实际输出全部给出半径 ≤0.5，CPU 中位约 0.90 秒。固定半轮的普通八起点在这 36 条也全都达到实际 RMS <0.1，耗时仅几十毫秒。[完整同预算恢复对照](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/research-z-alg-r14-20261119/research/research_z_alg_r14/EXPERIMENT_MEASUREMENTS.md)

## 请优先解决的两个问题

**空间接收自由度能否解除已证的频带歧义？** 对上述已知增益材料对，令 $d_r(\delta)=S_r(\delta)[J_1(\delta)-J_2(\delta)]$。同发射端口和同能量下，所有可实现新接收位置的差别都在接收行 $S_r$。请研究真实有限孔径上的空间上、下界，给出足以解除这对歧义且对小定位误差稳定的接收位置或接收组合。应比较固定几何选择、等成本增加带宽和接收通道；一对材料能分开后，仍需检查整盒剩余材料。镜像双材料有一个固定单侧馈电就能成功的严格对照，能检验收益是否来自测量几何本身。[镜像与错极小实验](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/research-z-alg-r14-20261119/research/research_z_alg_r14/GEOMETRY_RESULTS.md)

**能否把给定输出的误差半径算得既可靠又便宜？** 六维盒内，每个 $N_j-y_jq_j$ 用 64 个角值即可精确表达。现在使用少数对偶方向、按轴缩盒，再求输出到所有保留盒的最远角距离。请判断哪些传播或传感器结构能降低成本，是否有物理上自然的块低秩条件，使证书在有限材料加载或多个已知缺陷参数中可扩展。比较对象应有成熟的通用集合反演和廉价固定宽带求解，指标取相同输出误差半径、相同相容性以及全部成本。六个大区域各自包含许多网格单元时，材料参数未必只改变一行，需重新确定参数次数。

## 原创性审查要点

请把已经存在的集合反演、Bernstein 投影、参数电路多线性、反馈二元检验与本系统具体可计算的界逐项比较。R14 没有发现新的通用矩阵恒等式，也没有建立全面优于 CIE、A+ 或成熟宽带 GN 的恢复算法。最接近 Inverse Problems 的问题是有限设备和能量下可证明的材料误差下限，以及观测到达后可证明的误差上限，二者能否在实际可校准的物理模型上接近。连续介质、改变接收位置和跨网格模型的适用条件，需要新证据。用当前源重复一般选频、扩大无约束分支树或在 576 自由材料上直接枚举角点，现有证据均未显示优势。[数学推导及一手文献](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/research-z-alg-r14-20261119/research/research_z_alg_r14/MATHEMATICAL_THEORY.md) · [实验条件和复现入口](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/research-z-alg-r14-20261119/research/research_z_alg_r14/EXPERIMENTAL_RESULTS.md)
