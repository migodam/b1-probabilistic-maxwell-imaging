# E3：回访能纠错，但新调度没有净收益

本轮实际执行24个含噪案例、12个独立prior世界、每世界2个噪声副本，三次嵌套观测，共11种更新规则、792行结果。所有方法看到相同数据和同一个2950原子有限联合先验。已校准仪器、已知目标内禀材料、静态目标是本实验条件；这不是连续任意拓扑墙体posterior。

## 新数据、地图修正和目标变化

早期只有一个复端口，后续依次加入另一发射源的同端口和两个基准源剩余接收器。三份likelihood具有唯一ID；刷新旧factor只替换版本，不新增证据。普通增量联合更新与从原prior重算batch reference的最大权重差为浮点舍入量级。

完整24案例中：早期MAP墙组错误18例；后期MAP墙组正确21例；15例发生错误墙组到正确墙组的修正。早期支撑错误14例，后期支撑MAP全部正确；20例后期目标平方误差比早期下降，11例同时具有墙组纠错和目标误差下降。其中7例的真实墙非空。所有案例都留在`results/pilot_summary_tables/E3_ALL24_OFFLINE_MECHANISM_CASES.csv`，不以这些离线标记决定在线更新。

主cohort有12个不同原子真值，parent2/4为空墙，因此主统计同时包含10个非空墙和2个合法空墙。没有因结果不利或墙为空而重抽样。不能把这套先验描述成12堵复杂有孔墙。

例如非空墙parent5/repeat0：墙组MAP从116修正为真实95；目标平方误差从1.11523降至0.715323 λ0²，硬冻结终点为1.24456 λ0²。绝对目标误差仍大，地图修正不等于高精度成像。parent6/repeat0同样修正墙组，但普通联合终点平方误差0.250219，反而高于硬冻结的0.020646；单次真实误差不保证随正确Bayes更新单调下降。全cohort中硬冻结终点的实际位置误差更小有9/24例。

总体位置误差：硬冻结0.530718 λ0，普通增量联合更新0.346138 λ0。平方误差均值分别0.533227与0.354537 λ0²。这里比较的是同一批新观测下是否允许墙边际修正；不能把全部收益归因于新调度。普通smoother已经具备这种回访能力。

## 概率精度、风险和费用分开比较

下表为最终状态；CDF/联合TV先在parent内平均两个副本。费用是每条三次观测轨迹的归因cached CPU，包含belief summary、实际likelihood及必要reference依赖，不是新进程完整成像延迟。全部方法另外共享9288.023秒CPU、11800 RHS、2950 FE的参考表准备。

| 规则 | 48个后续决策中重开/刷新数 | 最终边际概率门通过 | 任务风险门通过 | 目标CDF差均值 | 联合TV均值 | 归因CPU ms |
|---|---:|---:|---:|---:|---:|---:|
| 普通incremental joint | 48 | 24/24 | 24/24 | ≈0 | ≈0 | 10.663 |
| 硬冻结墙 | 0 | 7/24 | 12/24 | .356879 | .726732 | 9.831 |
| 普通两轮墙–pose坐标MAP | 48条更新记录 | 9/24 | 14/24 | .214743 | .323338 | 8.797 |
| batch MAP点 | 48 | 10/24 | 15/24 | .144859 | .269602 | 8.521 |
| 周期全刷新 | 24 | 0/24 | 6/24 | .279387 | .683839 | 13.772 |
| 一直soft freeze | 0 | 0/24 | 1/24 | .416381 | .864597 | 13.064 |
| spectral-only | 48 | 24/24 | 24/24 | 0 | 0 | 14.711 |
| innovation-only | 2 | 1/24 | 2/24 | .396582 | .823227 | 14.542 |
| task-weighted，当前reference阈值 | 27 | 24/24 | 24/24 | 0 | 0 | 14.409 |
| qualified reopen，含完整reference审核 | 48 | 24/24 | 24/24 | 0 | 0 | 15.808 |
| online task proxy，旧belief阈值 | 25 | 22/24 | 22/24 | .038014 | .060083 | 14.112 |

原有边际门是目标坐标CDF、事件和重要**支撑模式**质量误差≤.02；风险门为reference Bayes risk的5%加`.005²`。它们不控制所有联合分支，故另报联合TV及墙–目标协方差误差，没有事后添加联合TV通过门。

`TASK_WEIGHTED`必须先得到当前完整reference风险才能设阈值，计算已计费；`QUALIFIED_REOPEN`还调用完整reference审核。它们不是便宜的在线证书。`ONLINE_TASK_WEIGHTED_PROXY`在观测生成前登记，只读自身旧belief，其2个失败案例保留。2,000次parent paired bootstrap仅描述12簇不确定性，详见`results/pilot_summary_tables/E3_PAIRED_BOOTSTRAP.csv`；e−16舍入差不能当显著优势。

所有方法实际新增Maxwell RHS均为0，因为采用同一预计算表。没有任何方法减少共享物理准备，合格的规则cached CPU也高于普通联合更新。更少名义REOPEN次数并未变成更少真实费用。当前实现不支持新的调度加速主张；这不证明未来稀疏、惰性刷新Maxwell算法不可能有效。

## 保留数值失败及辅助负控制

第一次回访程序在高斯moment surrogate的精度矩阵对称性审核处退出，未保存在线posterior、未新增PDE。修复只让不合格surrogate保守REOPEN，不放宽原门槛、不加jitter、不修改precision。正式结果中18次对称性审核拒绝、8次奇异/高度集中的covariance拒绝均保留原因。局部D6–D8只控制声明的小高斯系统，不能认证真实多峰posterior。

辅助负控制实际执行：已知空墙条件下25个pose原子，conditional更新与batch最大差3.33e−15，墙无不确定性时不需回访；继承的一维板–PEC精确补偿在6次独立重复观测后仍保持两分支各1/2，最大偏差3.38e−14。前者是已知墙先验，后者是一维解析反例，均不替代新二维未知墙体的精确非恢复世界对。

## 判决与比较名称

**受控回访机制：EXPERIMENTALLY_SUPPORTED。当前新调度独立优势：FAIL。完整E3合同：INCONCLUSIVE。** 尚缺新二维精确非恢复负控制；有限表的全局factor refresh不等于重新组装FE系统。`ITERATIVE_LAYER_POINT`实际是合法原子上的两轮墙–pose坐标MAP，不能称成熟连续layer-stripping；`BATCH_JOINT_MAP`是有限类穷举MAP，不能称新实现的连续batch FWI。真实连续FWI在独立A0中执行且未获得稳定锚点。

主张限于本实现和本有限prior。保留普通可回访联合推断，不以击败hard freezing证明原创调度，也不据此自动启动更宽材料类、NN或三维。

![误差与归因CPU](figures/pilot_actual/E3_error_vs_CPU.png)

来源：`results/e2e3/inference/`、`results/SUMMARY_FOR_PARENT.json`、`results/pilot_summary_tables/`、`results/e2e3/negative_controls/`和原始job receipt。
