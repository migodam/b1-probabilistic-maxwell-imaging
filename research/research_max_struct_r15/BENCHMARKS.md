# R15 数据、观测与评分合同

实际输入为原R14已公开24对象，完整ID为r14_new_0000至0023。R15在算法冻结后重跑；材料对象已经曝光，因此只能称后续配对实验。

- 二维scalar-TM、时间exp(-i omega t)、真空归一化lambda0=1，k=2pi，背景epsilon_b=mu_b=1；m=epsilon_r-1为实无损contrast。
- 代码中J=mQ是以epsilon_b归一化的对比极化/current变量；物理对比电流为-i omega epsilon_b J，不与外加源电流混记。内部G已经包含原disk-pulse积分与k²，S为其外部观测算子。
- ROI有24x24=576个独立实材料变量，cell中心坐标范围[-.46,.46]^2，圆盘pulse半径.018；box[.02,20]，m_init=.1。保留原smooth/TV正则，全部full forward与corrector用同一模型。
- E0为实际公开入射场矩阵；16对象12照明，8对象6照明；均48个接收通道，16对象全接收孔径、8对象半孔径。sources/receivers实际坐标在public包，不用非法极化或隐藏shape label。
- 0000/0001/0004/0005观测由36x36生成、在24x24反演。跨网格失配不被悄悄计入噪声球delta；noise-only CN99资格不能自动覆盖该失配。
- 原复杂圆对称高斯噪声以复幅度RMS sigma定义，relative RMS为1%或2%；若实虚各独立方差sigma²/2，则delta=sigma sqrt(chi2_(2d,.99)/2)。主目标D=.5||F(m)-Y||²/||Y||²与原R(m)相加，和原有限候选selector保持完全一致。

公开metadata逐例见audit/PUBLIC_ACQUISITION_SUMMARY.json。它只读取实际public标量与尺寸，没有从GT反推合同。S/G谱分区取自同public包；全/半孔径的实与复四空间分别缓存，不混用。

## 选择与评分

主候选是已经完成、合法、在预算内的原Phi最小者；secondary为同一已付费候选池中原CN99球内最小R，无合格者明确缺失并保留primary回退。所有预测封存后，才单独授权读取原R14 sealed_truth/independent24的material_coarse。

连续RMSE使用全部576cell；幅度NRMSE沿原评分器以max(coarse truth)-host为分母，boxrange NRMSE和relativeL2另列；支持IoU仅作预先存在的固定chi>1描述。本轮保留原事后联合统计NRMSE<=.15且IoU>=.70，不用于选择、停止或改参。主要结构/部分结构/幅值偏差的图像解释单独报告，不能倒改旧统计。

raw、corrected、primary、secondary与完整deployment回退分开。Phi下降、noise球可行与材料误差下降不是一回事。跨网格m_coarse只是逆网格参考，不是精细PDE的完全真值。

## 本轮未扩展的证据

没有新随机盲材料、没有新的3D vector或Fresnel实测反演。用户这次提出的twofold/周期DR机制先在同批数据上与原G19、S19和成熟基线做严格区分。历史R13/R14三维和外部数据保留原范围，不挪用为新算法的成功证据。
