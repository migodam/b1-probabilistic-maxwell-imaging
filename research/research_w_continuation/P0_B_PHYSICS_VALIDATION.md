# 新候选数值复核

每个已完成context按在线cal失配最小、stable method ID打破并列选择一个端点，未用truth选择。原P3像素网格与全部三角形再细分后的P2做相同材料/源/噪声对比。

|seed / 带宽 /状态|细化白化差|
|---|---:|
|2026102201 / narrow / cal|0.004319|
|2026102201 / narrow / task|0.198189|
|2026102201 / wide / cal|0.000586|
|2026102201 / wide / task|0.237604|
|2026102202 / wide / cal|0.003690|
|2026102202 / wide / task|0.270484|


cal一致不等于task一致。未贴合圆形目标的积分产生额外几何误差，不能把原参考世界的贴合网格资格扩展给所有新材料/任务。目标补充复核在既有validation预算，使用同一目标、相同连续总材料，只修正圆周界面处理。

补充贴合目标界面的实际结果（同一像素世界、同一个probe8，P3，圆周256→512段）：

|seed / 带宽|完整目标响应细化白化差|最细目标数据距离|最大RHS相对残差|
|---|---:|---:|---:|
|2026102201 / narrow|0.000644318|26.817818|1.256e-10|
|2026102201 / wide|0.000922609|47.754853|2.783e-10|
|2026102202 / wide|0.000915547|51.558657|1.670e-10|

三项均低于原0.1。为保持全部未知材料像素界面，圆周与像素网格交点会产生细瘦三角形；没有平滑物理几何、放宽1e-9求解残差或加jitter。网格最小角及aspect诊断公开，P3 DOF为35899/38203，复核峰值约1.065GB。这只是三端点的圆周细化一致性，不是全部12候选合格，也不是连续PDE上界。见[原始补充结果](results/p0b/fitted-task-endpoints/RESULTS.json)及[MESH_DIAGNOSTICS](results/p0b/fitted-task-endpoints/MESH_DIAGNOSTICS.json)。

[原复核结果](results/p0b/finer-endpoints/RESULTS.json)；实际NPZ复场/观测见Release的[数组清单](ARRAY_MANIFEST.json)。没有任何candidate达到cal<=1，没有最终GO端点微扰/全部24probe或joint复核；这些项目为NOT_RUN，不写成已成功。
