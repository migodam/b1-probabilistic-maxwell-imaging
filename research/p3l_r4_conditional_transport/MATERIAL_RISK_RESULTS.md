# 材料风险：已执行的估计器与独立评价

主表和区间见 [完整报告](RESEARCH_REPORT_ZH.md) 第4节及 [全部45方法原始表](evidence_summary/ALL45_METHOD_ACTIONS.csv)。主要动作是逐组件 clipping 后图像；raw 同时报告，不用 clipping 证明 decoder 合法。

## 实际模型与公平归因

R2 冻结 baseline 为5,842参数局部 CNN。新首轮直接 CNN 与它的函数类别相同，只需将 R2 residual_scale 吸收到最后一层权重和 bias；数值复播误差约1.21e−7。因此重新训练的收益记为强基线改善，不能叫新架构。

strong_spectral 使用1,024个付费 TRAIN physics contexts、256个 validation；decoder/random控制共享数据、容量和训练预算。静态谱相对随机空间包络的两批改善为1.20%/0.83%，不是动态解流形贡献。

最终主要 fusion head 有三层16通道3×3卷积、GELU及2通道head；输入8个原重建/物理通道，加4个全局材料均值实/虚与方差实/虚通道，共6,418参数。TRAIN-only normalizer，batch64、lr5e−4、weight decay1e−3、4,000 updates，三个固定seed。原基础图加残差输出。另需冻结provider约0.68M参数、K32及原物理反演费用。

provider 比较包括显式 Spline、同容量 diagonal Gaussian、直接全场矩网络，以及同图同值但固定按 parent 打乱空间对应的控制。主要模型按原R2VAL256选择，不按新DEV或新GT重选。两批均保留此选择。

## 独立主结果

Spline fusion 的 MSE 为0.00680623/0.00737782；R2为0.00783202/0.00832320；strong CNN为0.00719170/0.00775824。Gaussian fusion为0.00682170/0.00739262。Spline与Gaussian没有显著风险差，主要成像增益不能称为Flow独有贡献。

直接生成均值很差，但其全局信息能帮助局部模型。这个现象通过空间打乱控制与普通provider控制检验，属于有限信息表示与监督任务风险，不要求输出已校准的posterior。

## 风险与物理输出分离

新主要seed1 clipped图仅242/512、125/256满足noise-ball；它们是风险动作，尚未证明在原decoder image。合法代表在固定64对象上保留部分风险收益，却损失大部分原图精度；所有优化预算停止和fallback都保留。

不能用GT挑生成样本或候选当在线结果；本轮没有把Best-of-K作为主指标。统计单位为parent，三个seed是训练重复而非新增样本。所有原始每parent loss及canonical batch预测已由独立进程复核，loss最大差5.55e−17，bootstrap最大差2.15e−16，失败0。
