# 合法局部路径与材料chart曲率

两个原P0 TRAIN smooth材料只用于数学校准，不是新反演或GT发现。原8个独立uniform生成参数经Gaussian CDF坐标化，不扩support；共同Gaussian权重缩放在峰值归一化中抵消，所以材料rank7。

| TRAIN anchor，prior步长.75 | 方法 | 材料RMS位移 | 白化prediction位移 | 观测残差 | 显著新结构 |
|---|---|---:|---:|---:|---|
| 250033_weak | frozen_weak | 0.00134209 | 1.765752 | 22.313779 | NO |
| 250033_weak | second_order_curvature | 0.00095589 | 0.700696 | 22.298618 | NO |
| 250033_weak | updated_J_corrector | 0.00079099 | 0.036521 | 22.314053 | NO |
| 250034_ring | frozen_weak | 0.01019053 | 3.725476 | 23.009726 | NO |
| 250034_ring | second_order_curvature | 0.00941457 | 3.217260 | 22.977520 | NO |
| 250034_ring | updated_J_corrector | 0.00874272 | 2.765352 | 22.895485 | NO |


冻结弱方向取每单位材料RMS响应最小的非gauge方向，再归一成prior坐标单位。初始rank7奇异值全非零；两个最小值12.849和297.076。有限noise球内小宽度与精确等数据fiber不同。

二阶恒等式针对F(g(u))：D²(F∘g)[v,v]=D²F[Dg v,Dg v]+DF D²g[v,v]。第一项用真实R0反馈，第二项为decoder曲率；decoder二阶用两步长有限差分，非解析exact curvature。

| TRAIN chart | Maxwell二阶项范数 | decoder曲率项范数 | decoder/Maxwell | t=1冻结chart响应 | t=1实际响应 |
|---|---:|---:|---:|---:|---:|
| 250033_weak | 9.5209e-05 | 3.65231 | 38361.0 | 0.00642902 | 3.73365 |
| 250034_ring | 0.177677 | 7.52726 | 42.4 | 4.53431 | 5.3759 |


更多tangent→finite材料偏离和向量误差见POSTHOC_CONSISTENCY.json。第一例t=1材料secant/tangent cosine仅.5343，原弱切向已显著转动；Maxwell对实际材料增量的一阶误差却仅.446%。有限弱方向失效主要是decoder弯曲，不是Maxwell突然近奇异。

方法均不读取GT。updated-J每步刷新后保持弱方向定号，再作一次有界法向corrector；非收敛水平集追踪器。第二阶曲线修正和updated-J的材料位移不同；不能在相同prior t下宣称相同材料距离，更不能忽略其额外Jacobian/RHS。26条记录皆noise-compatible但0显著新结构，updated-J累计材料路径长度.000791/.009167，端点净位移.000791/.008743。

独立轨迹重放验证了执行代码定号缺陷没有改变本次六次sign，12个endpoint maxabs重放误差0。真实tangent角变化约.65°/1.06°与15.38°/6.62°。rank与active max cell在采样点不变，不是连续路径全区间的严格rank/noise证明。路径中点未穷尽检查，不称certified path。

known TRAIN自己的GT就是起点，离开anchor只会增大“到GT”的距离；这不构成成像风险改进。N3与physics R基线没有合法P0局部path启动条件，不运行其E3。
