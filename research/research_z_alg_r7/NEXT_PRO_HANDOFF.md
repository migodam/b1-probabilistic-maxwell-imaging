# 给下一轮 GPT Pro 的交接：哪一种材料量可以用有限照明识别？

R7 研究一个已知区域内由多个单元共享的实材料参数。在区外材料未知时，能否用少量固定的外部照明测出这个参数？研究完成了理论推导、全波数值实验、第二种收发几何的确认及宿主标定误差分析。全部观测由物理模型数值生成。建议依次阅读中文报告、推导第 3、8、9、10 节及全波对照。

[R7 网页报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r7/) · [ChatGPT 交接网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r7/handoff.html) · [发布版本及下载](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-z-alg-r7-20261118)

## 可直接读取的原始文件

- [当前交接全文 Markdown](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r7-20261118/research/research_z_alg_r7/NEXT_PRO_HANDOFF.md)
- [中文主报告](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r7-20261118/research/research_z_alg_r7/RESEARCH_REPORT_ZH.md) · [完整推导及定理](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r7-20261118/research/research_z_alg_r7/THEORY_DERIVATIONS.md) · [实验、计算费用和复现入口](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r7-20261118/research/research_z_alg_r7/EXPERIMENTAL_RESULTS.md) · [文献对照](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r7-20261118/research/research_z_alg_r7/LITERATURE.md)
- [两种区域的全波对照](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r7-20261118/research/research_z_alg_r7/results/RESULTS_PENCIL.md) · [逐样本结果 CSV](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r7-20261118/research/research_z_alg_r7/results/pencil_individual.csv) · [复读数和发射波形 NPZ](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r7-20261118/research/research_z_alg_r7/results/pencil_observations.npz) · [全波可行材料和残差 JSON](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r7-20261118/research/research_z_alg_r7/results/pencil_feasible_followup.json)
- [宿主失配 JSON](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r7-20261118/research/research_z_alg_r7/results/geometry_mismatch.json) · [空载标定和加载观测 JSON](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r7-20261118/research/research_z_alg_r7/results/geometry_calibrated.json) · [计入标定误差的逐样本界 JSON](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r7-20261118/research/research_z_alg_r7/results/host_error_bound.json)
- [连续介质数值实验](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r7-20261118/research/research_z_alg_r7/results/RESULTS_CONTINUUM.md) · [环形补偿及独立积分网格结果 JSON](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r7-20261118/research/research_z_alg_r7/results/continuum_mimic.json)

## 可用于论文的三项结果

1. **设计的对象是整个连续候选材料区间。** 在 $m$ 个活动坐标中，已知支撑的 $r$ 个坐标共享一个实材料 $x$，其余可以独立变化。一只外部接收器给出 $W=VD_z$ 的投影。对 $L=\operatorname{ran}B$，求 $(I-qH)Ep(q)\in L$ 的向量多项式。实际发射其系数场所张空间的一组基，由同一批复读数形成 $g(q)=(x-q)VEp(q)$。$r>m-\dim L$ 给出次数上界。若一般材料值的聚焦方向由一组模板张成且对应系数列整体独立，则任一在整个区间都能线性消除任意复干扰的固定波形组，曝光数至少为最低系数曲线次数加一。不要把这个下界转述成全部实材料反演器的下界。
2. **对照能复制一半结论。** 两个 $n=36,m=12,d=10$ 的全波几何中，四单元方形实际使用三束，两个独立局部聚焦方向之间需要协调；三单元 L 形实际使用四束，各锚点的局部设计已完全复制。确认方形同总能量 2、同联合半径 $2\times10^{-4}$ 的同源 TRF 误差中位数为 0.0183，局部锚点 0.1004，随机 2.505，全部十束 0.00919。一个 L 形样本的同一数据能接受目标 1.64787 和 6.65485。原始数据包含完整实材料见证，并保留求解失败与残差。[关键歧义](results/pencil_feasible_followup.json)
3. **连续材料改变结论。** 均匀二维全空间、任意有限个不同外部单极线源、独立相干波形、单独一个外部近场复接收器。在任意非空干扰开集，双侧实材料扰动到全部复通道的导数实满秩。“实解析延拓＋两个出射核的双倍径向相位＋Rellich”给出简短证明。用 $2s$ 个干扰材料函数补偿不相交的预设小目标，非线性 IFT 给出精确局部同数据族。对 +0.10 圆盘目标的 12 束有限幅度实例，960 点拟合后在 2,640 点独立积分得到未补偿信号的 $6.84\times 10^{-4}$。这是跨网格近似值；理论没有为 +0.10 给出显式充分半径。[材料图及保留的粗网格结果](results/continuum_mimic.json)

## 最需要 Pro 挑战的命题

**一，系数空间给出的 shot 最小性。** 请检查第 8 节的左逆/Vandermonde 论证是否需要比“在一般 $q$ 张成核＋所有系数列满秩”更强的假设，特别是秩跳跃、馈电零空间、可见性消失、不同模板次数与实数 $q$。更有价值的推广是系数列线性相关时，给出所有候选值都有非零聚焦方向的最小物理子空间维数，或给出这一推广的反例。请对照最小基干扰解耦、虚拟源实验和多源压缩；将滤波器的动态阶数直接等同于照明数会产生错误。

**二，标定不确定性到底损失多少分辨率。** 已知宿主材料统一提高 1% 后，真值覆盖由 32/32 降至 14/32。每种布局增加一次空载标定，再按估计的宿主倍数重新设计波形，样本覆盖恢复至 32/32。第 10 节用同一次标定读数、无耗散实材料的 $T$ 算子界和精确二点式，推导出对所有允许的活动材料都成立的读数误差上界。计入这一项后，目标区间的中位宽度由 0.1745 增至 1.665。观测到的系统误差投影最大仅约为其上界的 3%。请优先研究保留发射方向和接收方向相关性的更紧界，例如受限负载 $T$ 算子双线性支持函数。目标是用同一批数据降低外包宽度并保留模型内包含性。扩成多维漂移时，要说明一束标定能约束哪些参数方向。[标定记录](results/geometry_calibrated.json) · [统一误差界](results/host_error_bound.json)

**三，连续近场自动实满秩究竟推进了什么。** 请优先比对远场有限方向不可见材料构造和点电极 EIT 的局部灵敏度积独立性，它们已经提供了非线性校正与局部 Gram 对偶的证明框架。可争取的增量是近场点源自动满秩所用的出射相位，以及与实补偿代价的联系。应分析当接收/发射距离增大、源几何变坏，或干扰只许非负变化时，Gram 最小奇异值、补偿幅度和定理的有效半径如何变化。独立坐标的幅值约束与连续干扰允许的 $L^\infty$ 幅度也应分别陈述。[远场先例](https://arxiv.org/pdf/1410.8382) · [EIT 先例](https://arxiv.org/pdf/1412.6768)

## 论文应保留的边界

论文主问题可以定为“已知加载支撑、未知扰动和受限物理通道下的材料读出”。正文先陈述目标量和正反可辨条件，再展示方形胜例、L 形等价和双分量反例，最后用宿主标定及连续环形补偿检验物理外推。作为强 TAP 稿件，真实馈电标定、互耦和重复加载仍缺实测。作为 Inverse Problems 或 SIAM JIS 稿件，最紧要的是把上述三个命题与最近先例的差异写实，并说明结论对实际材料精度的影响。

已有的 Q20 根计数、fold、单点两束比值、一般 Fisher/多源叠加、动态最小基和局部隐身都应逐项归因。R7 没有证明未知支持成像、全 box 噪声材料集重建、任意实材料估计的普适采样下界、3D 或 PEC 的成果。强基线已复现的结果应继续用作正对照。[文献与原实验设计](LITERATURE.md) · [运行入口和原始结果](EXPERIMENTAL_RESULTS.md)
