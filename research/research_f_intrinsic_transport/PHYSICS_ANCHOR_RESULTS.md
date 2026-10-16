# 独立 physics-only anchor 实验

本轮真正运行了非神经网络、非GT初始化的原先验内 full-physics 反演。与旧一步 physics_R_base 不同，每次更新均求解原24×24 scalar-TM/VIE，构造完整接收伴随Jacobian，最小化原噪声白化的 Gaussian 负log likelihood。

| 条件 | 初始白化残差 | 最终残差 | 最优性误差 | 停止原因 |
|---|---:|---:|---:|---|
| p250073_ring_a1_snr20 | 91.360919 | 22.239555 | 0.045273396 | The maximum number of function evaluations is exceeded. |
| p250073_weak_a1_snr20 | 71.360219 | 21.512112 | 0.0020009887 | `ftol` termination condition is satisfied. |


冻结 tau=23.7854490041；两例均达到数据相容，并保持 P0_EXACT。冻结的 bound-scaled first-order optimality 门槛为0.001；两例均未达到，因此 GATE2=false。weak虽触发ftol，也不能用目标变化小替代实际驻点检查；ring在12次evaluation额度耗尽停止。没有根据结果改门槛或增加迭代。

## 已知条件与初始化

只选原validation元数据顺序中最早两个 smooth/ID、幅度1.2、SNR20条件：rows67/71，同parent250073的ring/weak背景。八个待估参数是两Gaussian的中心、宽度、权重；全局平移、旋转、scale、axis-ratio被当作已知条件，来自公开生成法的元数据回放。这比未知完整物体更容易，不能与N3的监督信息量当作相同。

初始Gaussian中心(-.052,0)/(.052,0)，宽度.0725/.0725，权重.65/.65；不依赖测量、N3、GT或GT误差。在线输入NPZ只有observation/sigma；JSON移除了八个真实未知参数与原seed。源码不读取offline_only；这是审查过的路径隔离，并非操作系统权限隔离。

## 如何判断可观测信息

残差进入noise球只说明数据相容。还报告原uniform坐标中的梯度、bound-scaled最优性、材料正交frame的响应谱及真实材料RMS的可行局部下降量。原参数GN步长受gauge与条件数影响，不能直接当作未恢复信息。

两最终点的七维材料响应均数值满列秩，最小每单位材料RMS响应分别265.43与52.82；这没有建立exact fiber。两Gaussian权重同时接近原支撑的.4/.9比值边界，共同尺度gauge不能消除这个真实比例边界。完整分析见evidence/INDEPENDENT_MATH_REVIEW.md。

驻点未达标不等于材料非法，也不等于不存在局部噪声不确定性。目标是真正Gaussian likelihood；95%noise球不构成posterior硬支持。没有用GT决定在线停止。

所有24个state和24份接收Jacobian均已持久化，含材料、电流、场、prediction、LU/pivots；额外真值误差只由离线图表脚本读取。
