# R10 真矢量3D执行交付

实际18条：v02/v09/v04/v08各3臂 matched12，以及v02/v09各3臂 cross18→12。36条是未执行草案。所有18原预测已封存；工人未读fresh GT，root独立评分。

源CIE与A+既有current beta/聚合差别；GlobalDual与A+ current完全匹配，仅材料机制不同。全源5184×60复场，材料共享1728实未知量；每stage实际APCS/Fourier自由DOF另见representation_rank回执。Phi3为逐源相对残差平方均值，不与2D原Phi或iid白化ML混同。

最终平台：v04三臂Mac CPU retained-dense/FFT；v08三臂XINAN CUDA retained-dense/FFT；v02/v09的matched及cross均Mac。全部源多RHS与adjoint共享各自同材料LU；辅助仍native FFT/PR。原18公共终态见PRIMARY_PUBLIC_TERMINALS.csv，表内费用为whole worker的nested视图，不能再加总到canonical。

旧Windows CPU v04、旧CUDA FFT/GMRES v08完整/partial与已付尾部保留RESULTS/vector_engineering_attempts。原18独立费用账为RESULTS/vector/ACCOUNTING_FINAL.json：Mac4095.307625、XINAN5962.796875 CPU秒已知下界；CUDA-event、prefix、nestedsetup及STOP重复GetProcessTimes均不重加。仍有未单计手工/初始helper成本null。root评分/作图及后续v04追加诊断独立收费。

模型限于有限CM+RR DDA，高介电连续误差不以旧Mie失败掩盖，两个fine18独立生成只刻画模型误差。0.6λ、12³是有三分量耦合的缩放案例，不是2019 3λ/30³→60³原尺寸复现。
