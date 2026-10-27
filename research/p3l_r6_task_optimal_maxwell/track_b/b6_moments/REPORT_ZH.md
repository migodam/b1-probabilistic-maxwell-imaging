# B6：真实缓存响应下的小型条件矩机制试验

**结果：完成真实回归试验，但没有获得可用的后验协方差或导数一致性。** 居中残差显著缩小raw c−ab的减法/回归误差放大，仍比同一个a模型的方向差分大很多；没有通过放宽目标、换VAL赢家或追加训练来救结果。

## 输入与资产资格

只读历史TRAIN32768、VAL2048。R4 joint archive没有noiseless字段，但可追溯到原R2 complex128 observation与每个TRAIN/VAL parent的offline_only/prediction、noise、q。实际读了全部34816份allowlist缓存：F+noise等于保存的observation逐值完全一致，cached q与packed source q逐值一致。直接使用保存的Maxwell prediction，不是从量化joint y32重建，更没有重新求解物理。源complex128减法与保存noise的最大roundoff为4.60e-17；乘sigma^-2约3.23e-11，远小于当前回归不一致。

三个目标是ROI材料实部均值、虚部均值、ROI内满足Re(chi)−Re(host)>.02的单元比例。第三项是明确材料阈值事件的空间比例，**不是decoder latent fraction**。原始平方损失是对应条件均值的proper监督目标；验证标签方差不等于posterior真值误差。

响应坐标统一为sqrt(2)×source-major pack(F)，噪声坐标方差为sigma²。固定3组各4个单位投影：canonical receivers、结构化source/receiver求和、固定seed正交随机。仅研究这3×4个函数对，不声称恢复Full66协方差。

## 冻结的小模型与交叉拟合

每个模型127特征加截距，共128特征，ridge=.01，未用VAL选超参数。raw_y_rbf使用所有512测量坐标的TRAIN标准化与127个固定TRAIN landmark径向特征；random_projection使用127个固定测量投影；decoder_only把标准化y的随机66维code经原C2/probit decoder变成材料，再作127个固定材料投影。后者是**未校准pseudo-latent的decoder-shaped随机特征控制**，没有anchor/J/Maxwell residual，不能冒充历史强anchor的decoder-only控制或后验。

TRAIN按parent索引奇偶分两半。每个方向先在另一半训练a0，再在本半训练a,b,c,bres,cres；a0映射在该方向固定，VAL标签不进入拟合。两个方向的a预测及各自cov估计平均。每条交叉拟合路径包含16384个moment训练对象，不能把两折称为独立随机seed复现。

|输入控制|mean Re normalized MSE|mean Im|材料阈值比例|
|---|---:|---:|---:|
|raw_y_rbf|0.014276|0.232618|0.045526|
|random_projection|0.002894|0.468539|0.112722|
|decoder_only|0.134359|0.557471|0.200161|

normalized MSE的分母仅用TRAIN标签方差，表中均为未裁剪VAL2048结果；不同任务的赢家不同，没有统一decoder优势。

## 条件矩与导数检验

拟合目标为a=E[T|y]、b=E[g|y]、c=E[Tg|y]，以及r=g−project(y)、bres=E[r|y]、cres=E[(T−a0(y))r|y]。若这些都是同一个真实posterior下的精确条件量，则Cov(T,g|y)=c−ab=cres−(a−a0)bres，沿单位观测投影的Da=Cov/sigma²。**这个恒等式只属于精确posterior；估计器之间不自动满足。**

在预先固定的VAL前128对象上，直接对同一个a预测器作12个观测方向的中央差分，步长sigma×1e-3与sigma×1e-4。两步差分相对差约1.6e-9至5.0e-9，说明差分本身稳定；它不提供真实posterior导数。

三个输入×三投影组的raw相对一致性gap范围1744.8–37952.5；居中后仍为5.28–22.54。该gap表示当前有限回归的条件矩与a导数未相互一致，不能单独证明a或cov的总体bias，也不能把居中相对改善称为posterior校准通过。sigma^-2=703657.8会放大很小的矩回归误差。

|输入|a normalized MSE|b|c|bres|cres|
|---|---:|---:|---:|---:|---:|
|raw_y_rbf|0.09988|0.05878|0.08344|1.01630|1.03927|
|random_projection|0.19515|0.00496|0.07268|0.94229|1.00254|
|decoder_only|0.29810|0.38912|0.42107|1.01875|1.02762|

此表是每条合法交叉拟合路径在独立VAL上的proper moment MSE，再平均两折；a与上表ensemble risk口径不同。bres/cres约1意味着残差矩的conditional信号未被这些小回归器有效提取。raw c−ab的中位相消比约11.7–20.8、p95约195–315；这是有限回归误差与相消在小噪声下被放大的实际证据，不是机器精度导致的单一失败。ridge正则化矩阵条件数最大5685.6，没有触发数值ill-conditioning停止线。

## 成本与停止边界

实际12个线性fit，总CPU 21.304s，wall 32.449s，单线程，peakRSS 1.385GB；新增物理RHS=0、MPS=0。只做固定pilot，没有再调feature数、ridge、投影或数据。无旧TEST/D truth读取，没有新Full MCMC或Full66 covariance声明。

实现：track_b/b6_moments.py；本报告与图只读已存结果，未增加fit。PLAN_FROZEN.json、ASSET_RECOVERY.json、FIXED_PROJECTIONS_AND_FOLDS.npz、JOB.json、RESULTS.json、DIRECTIONAL_CONSISTENCY.csv、MOMENT_PROPER_MSE.csv、COVARIANCE_AND_FD.npz保留完整证据。图源为TASK_PROPER_MSE.csv和DIRECTIONAL_CONSISTENCY.csv。

![真实小型机制与一致性失败](B6_TASK_AND_CONSISTENCY.png)
