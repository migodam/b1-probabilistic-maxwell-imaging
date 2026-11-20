# Z-ALG R16｜给下一轮 ChatGPT 的研究交接

我们在感应电流域中降阶，每轮在选定的 Fourier 或物理切线方向上更新电流，接着让全部 576 个实体素共同解释所有发射。主实验使用 MAX-STRUCT 原存档，并直接调用原生 FFT-CIE 做同预算消融。当前最值得复核的主张是，有限电流空间对一个已选材料方向的表示能力有可达到的上限，达到它对实际恢复的帮助仍受当前下降方向、传播和优化目标限制。

[中文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r16/) · [可复制交接页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r16/handoff.html) · [含原始输入、代码、逐案结果与图的 ZIP](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/research-z-alg-r16-20261120/Z_ALG_R16_2026-11-20.zip)

## 已经证明和已经测到的

1. 记 \(L_m=I-D_mG\)，\(E=E_0+GJ,\ T=D_b E\)。对于电流基 \(Q\)，\(\min_a\|L_mQa-T\|^2=x^\top(I-K)x\)，\(x_i=\|E_i\|b_i\)。所有格均被照亮且同频各源共享电流基时，\(\operatorname{tr}K=\mathrm{rank}(Q)\)。这是实共享材料本构与复电流投影的直接计算，作用域是固定线性化点。多频按各格的场能量加权。[固定版本推导](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r16-20261120/research/research_z_alg_r16/MATHEMATICAL_THEORY.md)
2. 固定强接收分量 \(V^*X=0\)，令 \(Z=L_m^{-*}V\)。全体暗电流的最优方向是 \(X_*=L_m^{-1}(I-P_Z)T\)，最佳误差为 \(\|P_ZT\|\)。在旧暗基上补入 \(X_*\) 的列，至多新增发射数那么多个方向即可达到这个上限。56 份带噪历史观测中，固定 48 暗方向在均匀背景处可表示实际下降方向能量中位 33.11%，无限制增加暗基的上限 33.48%。体素拟合后 55 个非零下降方向上分别为 33.95%、78.96%。在所有非零测试点上，按 \(L_m\) 度量选模与上限的差最大为 \(1.23\times10^{-15}\)。[实验报告及逐案入口](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r16-20261120/research/research_z_alg_r16/EXPERIMENTAL_RESULTS.md)
3. 保留 576 个独立材料，96 案从同样旧观测比较。使用 MAX-STRUCT 的原物理目标和优化器尺度、至多 180 次完整分解，直接拟合像素、普通 48 模交替、候选读数传递、按度量取最优暗电流的材料达标数依次为 8、11、9、9。前两者 CPU 中位数均约 1.38 秒。8 个跨网格对象均未达标。不能从“对一个方向完全表示”推出同预算材料恢复改进。[180 次逐案](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r16-20261120/research/research_z_alg_r16/results/current_cases_archive180.csv)
4. 在真正的 R10 原生 GlobalDual-CIE 上，同 24 案、同 200 步，只缩减末段 Fourier 坐标至四分之一，CPU 节省 3.2%，接收误差中位由 2.84% 变为 14.32%。用首阶段梯度换走六个已有 Fourier 方向使 21/24 案材料 RMSE 上升。保留原方向再加 6 个，同维 Fourier 在 18/24 案比梯度方向的材料 RMSE 更低。原目标下的环、开缝壳和包覆包体坏端点均未被本轮有限电流提案救回。[原生对照](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r16-20261120/research/research_z_alg_r16/NATIVE_CURRENT.md)

5. R15 既有板上 8 份缺陷各在无损和已知色散下取数，另有空板对照。固定 8 次发射、8 个接收、64 个复读数及总能量。使用双频 (2,30) GHz/4 发射角，完整体素 GN 的无损/色散 RMSE 中位 .0568/.0594，优于单频 16 GHz/8 角和 4 频/2 角。每频 49 个低维暗电流、所有频率共同读回全图，两种宽带安排的中位支撑 IoU 都为零；平衡频率权重仍遗漏所有超过阈值的缺陷。保持每频总阶数49并允许强接收模一起更新，32次非空配置比较的完整物理残差均更高。多频加权迹随场分配在 48.982 至 49.149 间。此例支持区分设备读数的有效信息与某个当前空间利用它的能力。[宽带实算](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/research-z-alg-r16-20261120/research/research_z_alg_r16/BROADBAND_CURRENT.md)

## 真正剩下的数学问题

- 对**整个候选材料变化集**，能否从数据给出有物理依据的权重，使有限电流更新在相关的多个方向上都保持近似精度？茫然把随机方向加入当前空间，或对单个梯度做到 100%，都没有解决“当前数据应支持哪一组有限材料变化”。可从试探材料的可计算两点差 \(D_{m'-m}E_{m'}\) 与接收可分离性联合入手。应与加同样数量 Fourier 模式、普通物理线搜索和原生 CIE 同算力比较。
- 什么幅度内，式 (6) 的本构方向上限能控制**有限步**的全波数据误差？式 (7) 已把传递误差分解为材料步幅、本构误差和两个传播因子。若其可见部分能在不重复高成本正演的条件下得到紧界，才有可能用它决定何时采用电流更新、何时直接更新体素。现有逐候选反例表明，沿用辅助电流总场到晚期，实际会比普通一阶近似更差。
- 宽带下的 \(w_{fi}(P_f)_{ii}\) 量化同一材料怎样在若干频率的电流基里被表示。观测有效秩还要乘上各频率的接收传播。下一步的物理问题是：在等能量、等曝光的设计下，何种公开的散射结构能同时改善这两个量并改善恢复？R15 已有材料对排序反转，无法仅按频率数或发射数给出统一排序。[上轮宽带比较](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r15/)

历史 Q6/Q20 根数、R3 冻结补集与 R10 单接收照明零空间均已与此处区分。数据和源码支持继续判断这些命题的适用域。目前的证据支持一个方向诊断和普通两步更新的小幅收益，尚未建立超越成熟 CIE 的通用成像算法。
