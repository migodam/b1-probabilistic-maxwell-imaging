# R6：A路径扩大比较与FFT-TSOM接入

这是用户明确要求的扩展研究，保留R5证据，不能以其局部负结果代替广泛算法比较。算法/预算在读取新对象GT前冻结。主比较完成之前不根据某个case的评分更改方法。

## 数据和方法

15个已公开开发对象：FBG-R3的数字1–9及其余3目标，加R5旧弯条和两个R4失败弯条。另由独立生成端按固定设计生成6个新材料对象，GT与label留sealed_truth，求解只接公共观测清单。共21个材料、84个主任务。

共同576个独立实材料，12源/48收，原Maxwell VIE、box、host与公开prior不改。历史噪声按各原合同；6个新对象统一complex CN .1% RMS，公开σ/噪声维数。新观测使用同一离散，明确是合成同离散确认，不冒充连续PDE或实测。

四臂：mature sourcewise CIE、mature sourcewise Dual-CIE、Global-Dual控制、A-material（普通global current + one-hop材料读出）。CIE是每个阶段current/material用同一个β；Dual-CIE是后两阶段current β=1、material β=.2。CIE单β指每次更新的同β，不是本轮搜索某个固定β。

## 充分运行与公平性

所有臂同mf3/6/10各500步，再在同一调用中自然扩到mf12共1500步，总3000辅助步。Fourier原系数按原代码零填充，PR仅在原阶段边界重置；不以更换材料初始化代替延续。mf12覆盖完整Fourier网格，APCS保留其原DPCS正交补语义。保存1500步后的阶段状态，但最终不据GT选早停。

每臂从原Phi最低的初始/阶段候选出发，接同一个original-Maxwell corrector最多3000次接受迭代，原ftol/gtol不变。记录100/500/1000/2000/3000的已付费接受点、全部已付费目标/梯度历史。没有45秒墙限。预算与mode增长在本次比较中同时给全部方法，不能宣称未收敛任务已达到算法能力上限。

同日程/同容量不等于等实际费用。比较最终材料与原Phi，同时用实际LU、complex RHS、接受步和wall/CPU轨迹显示费用；source归一化和prior始终同合同。旧R5结果仅作为历史参照，不假装是这次较长日程的结果。

## FFT工程臂

A的材料读出不依赖current用显式Fourier矩阵或FFT实现，因此结构上可接到FFT-TSOM。实际matrix-free实现必须同时替换U、H、G及其adjoint，不能只把IDFT列生成换成FFT便声称FFT复杂度。零填充Toeplitz卷积保留自由空间离散Green及原自项，禁止周期物理替代。

在旧b20d及一个预先选定历史数字（public清单第一项）上，完整运行FFT Dual-CIE/A与同日程dense参照。这是固定两例的计算接线与成本比较，非GT挑选快或成功的案例。局部浮点接线一致不保证非线性长轨迹逐位一致；终点实际差异保留。最终物理验收和corrector仍用相同dense full-wave LU，FFT只加速辅助current块。

NIE只有在查到实际相容实现/约定后才作为额外实测算法。self absorption不自动叫NIE，不以更弱替代品规避CIE/Dual-CIE主比较。

## 评分、成本与执行

主任务全部完成并封存后，由独立生成/评分端读GT。RMSE分母为公开contrast amplitude，IoU沿用原Otsu合同并保留绝对材料指标；新噪声discrepancy从公开复噪声取得。所有已发布对象仍标开发，6新材料单列；不把噪声样本或source列当独立材料。

运行在本地CPU，四个有限worker、每个BLAS单线程。远端只读盘点发现已有Python进程，不占用或终止别的任务；本轮无需移植A的细粒度材料readout到CUDA才能完成主比较。保留失败、导入/setup、观察/绘图、每次LU/RHS与总队列历时。源码快照而非SHA验证。
