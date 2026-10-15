# 方法比较与未测部分

|方法|实际运行|可支持的结论|
|---|---|---|
|同N3全材料GN／Haar64资格|两臂8case各一次有效更新|均获得8/8数据相容；Haar较少材料偏离|
|相同有限编辑，无补偿|40个实际raw端点|0/40相容|
|普通受约束full GN|三轮40条轨迹|32条相容且保住结构；成像MAE无改善|
|额外正材料相保护|相同16raw与步数|10/16 vs15/16；额外near-air减少，精度无收益|
|minimum-change合格base|同4parent固定8编辑|7/8；baseline质量较好但候选MAE仍大|
|constrained score排序|同paid端点轨迹离线重放|有完整取得score成本；不是新online加速试验|
|A17 paired／同成本Krylov|NOT_RUN|未建立独立机制增量；完整普通Jacobian已有feedback|
|QD/CEM／全搜索multistart|NOT_RUN，本任务属于G分工|不能宣称T胜过其搜索效率|
|新3D／GPU NN|NOT_RUN／不授权|无一般3D外推|

RANKING_COST_REPLAY.csv中U排序先支付每个raw的1state/8RHS与32-adjoint初始评分，再支付被修正端点的剩余成本；不会把已知成功顺序当免费oracle。每个编辑实际约128RHS，其中首次Jacobian可由score与corrector共享，另有两次refresh。score的额外小矩阵优化wall未独立分离；所有pipeline wall包含在总job费用。

OFFLINE_BEST_OF_K.json按基础图再固定edit次序，K=1/2/4/8；不足K标实际数量，不复制补齐。所有case含base的Best-of-up-to-K均无MAE收益。完整MAE、support、ROI、Re/Im、边界和几何结果见offline/MATERIAL_METRICS.csv。truth只作离线评分。
